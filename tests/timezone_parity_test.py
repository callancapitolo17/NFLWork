"""Regression gate — every scraper's game_start_time matches Odds API commence_time.

Cross-references each MLB scraper's stored game_start_time (TIMESTAMPTZ, UTC) against
the Odds API's commence_time for the same (home_team, away_team) pair within a 60s
tolerance. Run after any scraper-touching change.

Exit code 0 = all scrapers within tolerance (or no overlap to test).
Exit code 1 = at least one (scraper × game) pair outside tolerance.

Compares ALL 8 MLB-relevant scrapers:
  - DK, FD, BFA, Kalshi (Group A — always wrote ISO UTC internally)
  - Wagerzon (EDT → UTC), Hoop88 (PDT → UTC, surprise from audit),
    Bookmaker (PDT → UTC, the headline fix), Bet105 (UTC, no conversion)
"""
import os
import sys
import json
from datetime import datetime, timezone, timedelta
from pathlib import Path
from urllib.request import urlopen, Request

import duckdb

# Tolerance: 60s covers Odds API's ~minute-level precision + intra-cycle jitter.
TOLERANCE_S = 60

# Path → scraper label. Mirror tests/timezone_audit.py SCRAPERS dict for parity.
REPO_ROOT = Path(__file__).resolve().parent.parent
SCRAPERS = {
    "dk_odds/dk.duckdb":             "dk",
    "fd_odds/fd.duckdb":             "fd",
    "bfa_odds/bfa.duckdb":           "bfa",
    "kalshi_odds/kalshi.duckdb":     "kalshi",
    "wagerzon_odds/wagerzon.duckdb": "wagerzon",
    "hoop88_odds/hoop88.duckdb":     "hoop88",
    "bookmaker_odds/bookmaker.duckdb": "bookmaker",
    "bet105_odds/bet105.duckdb":     "bet105",
}


def _load_api_key():
    """ODDS_API_KEY env first, then ~/.Renviron fallback.

    ~/.Renviron may use either 'KEY=value' or 'KEY = value' (R is whitespace-lenient),
    so we match the prefix loosely and partition on '='.
    """
    key = os.environ.get("ODDS_API_KEY")
    if key:
        return key.strip().strip('"').strip("'")
    renviron = Path.home() / ".Renviron"
    if renviron.exists():
        for line in renviron.read_text().splitlines():
            if line.strip().startswith("ODDS_API_KEY"):
                _, _, val = line.partition("=")
                return val.strip().strip('"').strip("'").split()[0]
    raise RuntimeError("ODDS_API_KEY not found in env or ~/.Renviron")


def fetch_odds_api_games():
    """Fetch next 24h of MLB games from Odds API, keyed by (home, away)."""
    key = _load_api_key()
    url = (f"https://api.the-odds-api.com/v4/sports/baseball_mlb/odds"
           f"?apiKey={key}&regions=us&markets=h2h&oddsFormat=american")
    games = json.loads(urlopen(Request(url)).read())
    # Doubleheader-safe: keep a list per team pair, pick closest at compare time.
    api = {}
    for g in games:
        ht, at = g["home_team"], g["away_team"]
        commence = datetime.fromisoformat(g["commence_time"].replace("Z", "+00:00"))
        api.setdefault((ht, at), []).append(commence)
    return api


def closest_api_match(api_list, scraper_dt):
    """Pick the api commence_time closest to scraper_dt (doubleheader-safe)."""
    if not api_list:
        return None
    return min(api_list, key=lambda d: abs((d - scraper_dt).total_seconds()))


def check_scraper(db_path, label, api_games):
    """Returns list of (label, ht, at, delta_s) failures (empty = all good)."""
    full = REPO_ROOT / db_path
    if not full.exists():
        print(f"[{label}] DB missing — skipping (run scraper first)", flush=True)
        return []
    con = duckdb.connect(str(full), read_only=True)
    try:
        # If the DB exists but has no mlb_odds table yet, the scraper hasn't been
        # run post-migration. Treat as a soft skip (not a failure).
        tables = [r[0] for r in con.execute(
            "SELECT table_name FROM information_schema.tables "
            "WHERE table_schema='main'"
        ).fetchall()]
        if "mlb_odds" not in tables:
            print(f"[{label}] mlb_odds table missing — skipping (run scraper first)",
                  flush=True)
            return []
        # Verify the column exists; if not, this is a stale-schema DB — surface loudly.
        cols = [r[0] for r in con.execute("DESCRIBE mlb_odds").fetchall()]
        if "game_start_time" not in cols:
            print(f"[{label}] FAIL — schema does not have game_start_time column "
                  f"(found {cols}); migration not applied?", flush=True)
            return [(label, None, None, None)]
        rows = con.execute(
            "SELECT DISTINCT home_team, away_team, game_start_time "
            "FROM mlb_odds WHERE game_start_time IS NOT NULL"
        ).fetchall()
    finally:
        con.close()

    if not rows:
        print(f"[{label}] no rows with game_start_time — skipping (empty scrape?)", flush=True)
        return []

    failures = []
    matched = 0
    for ht, at, gst in rows:
        # game_start_time should be tz-aware UTC; defend against naive.
        if gst.tzinfo is None:
            gst = gst.replace(tzinfo=timezone.utc)
        api_list = api_games.get((ht, at))
        if api_list is None:
            continue   # game not in Odds API (e.g. Kalshi event without an Odds API match)
        match = closest_api_match(api_list, gst)
        delta_s = abs((gst - match).total_seconds())
        if delta_s > TOLERANCE_S:
            failures.append((label, ht, at, delta_s))
        else:
            matched += 1
    print(f"[{label}] {matched} game(s) within {TOLERANCE_S}s tolerance, "
          f"{len(failures)} failure(s)", flush=True)
    return failures


def main():
    print(f"Running TZ parity test (tolerance: {TOLERANCE_S}s) at {datetime.now(timezone.utc)}")
    api_games = fetch_odds_api_games()
    print(f"Odds API: {sum(len(v) for v in api_games.values())} game commence_time(s) "
          f"across {len(api_games)} team-pair(s)", flush=True)
    all_failures = []
    for db_path, label in SCRAPERS.items():
        all_failures.extend(check_scraper(db_path, label, api_games))
    if all_failures:
        print("\n=== FAILURES ===", file=sys.stderr)
        for label, ht, at, delta_s in all_failures:
            if delta_s is None:
                print(f"  {label}: schema missing game_start_time", file=sys.stderr)
            else:
                print(f"  {label}: {ht} vs {at} — Δ {delta_s:.1f}s "
                      f"(>{TOLERANCE_S}s tolerance)", file=sys.stderr)
        sys.exit(1)
    print("\nPASS — all scrapers within tolerance.", flush=True)
    sys.exit(0)


if __name__ == "__main__":
    main()
