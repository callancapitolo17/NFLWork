"""Empirically verify each scraper's recorded game_time TZ against Odds API.

For each (home_team, away_team) appearing in both the scraper's mlb_odds and
the Odds API response, computes:

    Δ = (odds_api_commence_time_UTC) − (scraper_recorded_dt_treated_as_UTC)

The modal Δ across games reveals the actual source TZ of the scraper.

Sign convention: Δ is HOW MANY HOURS THE SCRAPER'S CLOCK IS BEHIND UTC.
  - Δ ≈ 0    → scraper writes UTC
  - Δ ≈ +4   → scraper writes EDT (4h behind UTC)
  - Δ ≈ +5   → scraper writes EST
  - Δ ≈ +7   → scraper writes PDT (7h behind UTC)
  - Δ ≈ +8   → scraper writes PST

Doubleheader / multi-day disambiguation: a single (home, away) can match
multiple Odds API rows (e.g. a 4-game weekend series). We pick the API row
whose commence_time is closest in calendar-time to the scraper's recorded
local time, ignoring TZ, so that a 7pm-local game on 5/22 maps to the 5/22
evening API row even if their UTC offsets differ.

This is a one-shot diagnostic. Run it, eyeball the output, fold the modal Δ
into Phase 2's migration logic.

Findings get recorded in tools/TZ_AUDIT_FINDINGS.md.

NOTE: Scraper DuckDBs live in the main NFLWork directory, not the worktree.
Per CLAUDE.md, we never symlink DuckDBs (WAL writes corrupt the original on
worktree removal), so we open them read-only at their canonical absolute path.
"""
import duckdb
import json
import os
from collections import Counter
from datetime import datetime, timezone, timedelta
from pathlib import Path
from urllib.request import urlopen, Request

# Where scrapers actually store their DBs — outside the worktree.
NFLWORK_ROOT = Path("/Users/callancapitolo/NFLWork")

# Relative path (inside NFLWORK_ROOT) → scraper label
SCRAPERS = {
    "wagerzon_odds/wagerzon.duckdb":   "wagerzon",
    "hoop88_odds/hoop88.duckdb":       "hoop88",
    "bookmaker_odds/bookmaker.duckdb": "bookmaker",
    "bet105_odds/bet105.duckdb":       "bet105",
}


def _load_api_key():
    """Pull ODDS_API_KEY from env first, fall back to ~/.Renviron."""
    if os.environ.get("ODDS_API_KEY"):
        return os.environ["ODDS_API_KEY"].strip().strip('"')
    renviron = Path("~/.Renviron").expanduser()
    if not renviron.exists():
        raise RuntimeError("ODDS_API_KEY not in env and ~/.Renviron not found")
    for line in renviron.read_text().splitlines():
        if line.strip().startswith("ODDS_API_KEY"):
            # e.g. ODDS_API_KEY=abc123  or  ODDS_API_KEY="abc123"
            _, _, val = line.partition("=")
            return val.strip().strip('"').split()[0]
    raise RuntimeError("ODDS_API_KEY entry not found in ~/.Renviron")


ODDS_API_KEY = _load_api_key()


def fetch_odds_api_games():
    url = (
        "https://api.the-odds-api.com/v4/sports/baseball_mlb/odds"
        f"?apiKey={ODDS_API_KEY}&regions=us&markets=h2h&oddsFormat=american"
    )
    return json.loads(urlopen(Request(url)).read())


def parse_scraper_dt(date_str, time_str):
    """Parse scraper-written (game_date, game_time) as naive UTC.

    Handles both 'MM/DD' (year-less; what WZ/Hoop88/BKM/Bet105 actually write)
    and 'YYYY-MM-DD'. For year-less dates, picks the year that produces the
    closest match to 'today' (handles Dec→Jan rollover).
    """
    time_part = time_str.strip()
    date_part = date_str.strip()

    today_utc = datetime.now(timezone.utc).date()
    if "/" in date_part and len(date_part.split("/")) == 2:
        # MM/DD — infer year
        mm, dd = date_part.split("/")
        # Try current year first, then ±1, pick closest to today
        candidates = []
        for yr in (today_utc.year - 1, today_utc.year, today_utc.year + 1):
            try:
                d = datetime(yr, int(mm), int(dd))
            except ValueError:
                continue
            candidates.append(d)
        if not candidates:
            raise ValueError(f"unparseable date {date_part!r}")
        best = min(candidates, key=lambda d: abs((d.date() - today_utc).days))
        base_date = best.date()
    else:
        # YYYY-MM-DD
        base_date = datetime.strptime(date_part, "%Y-%m-%d").date()

    hh, mm = time_part.split(":")
    return datetime(
        base_date.year, base_date.month, base_date.day,
        int(hh), int(mm), tzinfo=timezone.utc,
    )


def tz_guess(delta_h):
    """Given Δ = (api_utc − scraper_treated_as_utc), guess source TZ.

    Δ here measures how many hours scraper's clock lags UTC.
    """
    mapping = {
        0.0: "UTC",
        5.0: "America/New_York (EST)",
        4.0: "America/New_York (EDT)",
        8.0: "America/Los_Angeles (PST)",
        7.0: "America/Los_Angeles (PDT)",
        6.0: "America/Chicago (CDT) or America/Denver (MDT)",
    }
    return mapping.get(round(delta_h, 1), f"unknown offset {delta_h:+.1f}h")


def _best_api_match(scraper_dt_naive_utc, api_candidates):
    """Pick the API game whose commence_time is closest to scraper's wall-clock.

    api_candidates is a list of (commence_dt_utc,) values for the same
    (home, away). We compare them to scraper_dt_naive_utc and return the
    one with smallest |abs diff|, plus that diff in hours.
    """
    best = None
    best_diff_h = None
    for api_dt in api_candidates:
        diff_h = (api_dt - scraper_dt_naive_utc).total_seconds() / 3600.0
        if best_diff_h is None or abs(diff_h) < abs(best_diff_h):
            best = api_dt
            best_diff_h = diff_h
    return best, best_diff_h


def audit():
    api_games_raw = fetch_odds_api_games()
    # Key (home, away) → list of commence_times (handles doubleheaders / series)
    api_games = {}
    for g in api_games_raw:
        key = (g["home_team"], g["away_team"])
        ct = datetime.fromisoformat(g["commence_time"].replace("Z", "+00:00"))
        api_games.setdefault(key, []).append(ct)

    print(f"\nOdds API returned {len(api_games_raw)} MLB rows across "
          f"{len(api_games)} distinct matchups")
    print(f"Now UTC: {datetime.now(timezone.utc).isoformat()}")
    print("Sample API entries:")
    for i, (k, v_list) in enumerate(api_games.items()):
        if i >= 3:
            break
        print(f"  {k[1]} @ {k[0]}: {[dt.isoformat() for dt in v_list]}")

    for rel_path, label in SCRAPERS.items():
        full_path = NFLWORK_ROOT / rel_path
        if not full_path.exists():
            print(f"\n[{label}] missing DB at {full_path} — skip")
            continue

        con = duckdb.connect(str(full_path), read_only=True)
        try:
            rows = con.execute(
                "SELECT DISTINCT home_team, away_team, game_date, game_time "
                "FROM mlb_odds WHERE game_time IS NOT NULL"
            ).fetchall()
        finally:
            con.close()

        deltas = []
        unmatched = []
        unparseable = []
        for ht, at, gd, gt in rows:
            candidates = api_games.get((ht, at))
            if candidates is None:
                unmatched.append((ht, at, gd, gt))
                continue
            try:
                scraper_dt = parse_scraper_dt(gd, gt)
            except (ValueError, AttributeError) as e:
                unparseable.append((ht, at, gd, gt, str(e)))
                continue
            api_dt, delta_h = _best_api_match(scraper_dt, candidates)
            # delta_h = api − scraper. POSITIVE = scraper is behind UTC by N hours.
            deltas.append((ht, at, gd, gt, round(delta_h, 2), api_dt))

        print(f"\n=== {label} ({len(deltas)} matched / {len(rows)} distinct rows in DB) ===")
        if unparseable:
            print(f"  unparseable date/time on {len(unparseable)} rows; first: {unparseable[0]}")
        if unmatched:
            print(f"  {len(unmatched)} DB rows had no Odds API match")
            for u in unmatched[:5]:
                print(f"    unmatched: {u[1]} @ {u[0]}  ({u[2]} {u[3]})")
        if not deltas:
            print("  no matches — cannot infer TZ")
            continue

        for ht, at, gd, gt, d, api_dt in deltas:
            print(f"  {at} @ {ht}  scraper={gd} {gt}  api={api_dt.isoformat()}  Δ = {d:+.2f}h")
        delta_values = [d for _, _, _, _, d, _ in deltas]
        # Round to nearest 0.5h to absorb minor differences (some books round
        # game times to the 5min, while Odds API gives exact first-pitch).
        rounded = [round(d * 2) / 2 for d in delta_values]
        counter = Counter(rounded)
        modal, freq = counter.most_common(1)[0]
        print(f"  → MODAL Δ: {modal:+.1f}h  ({freq}/{len(deltas)} games)  ⇒  {tz_guess(modal)}")
        if freq < len(deltas):
            print(f"    other rounded-to-0.5h Δ values: {dict(counter)}")


if __name__ == "__main__":
    audit()
