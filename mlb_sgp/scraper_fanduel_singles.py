"""FanDuel single-leg odds scraper.

Walks every MLB event today, fetches all runners + market metadata via
FanDuelClient, transforms to the offshore mlb_odds schema (wagerzon-style),
and writes to fd_odds/fd.duckdb.

Mirror of scraper_draftkings_singles.py (Tasks 8+9) but with two FD-specific
quirks the DK parser doesn't need to handle:

  1. FD ALT-spread runners encode the line in the NAME, not in `handicap`.
     `Runner("Cardinals +3.5", line=None, ...)` — we have to regex the name.
     (DK puts the line in both name AND handicap, so DK's parser uses handicap.)

  2. FD ALT-total runners use the format `"Over (8.5)"` — name has the line in
     parens, `handicap` is 0/missing. Again the parser regex-parses the name.

For MAIN markets FD uses `handicap` (line=-1.5 on the Runner), so those flow
through identically to DK.

Other FD specifics:
  - classify_market is a keyword classifier (mirrors DK's approach) that
    rejects FD-specific junk via exclusion keywords + event-team-name checks
    and auto-covers F7/F3 periods and any new FD markets going forward.
  - No team-name canonicalization map — Phase 1 confirmed FD team names
    already match the canonical Odds-API forms (e.g. "St. Louis Cardinals",
    "Athletics", "New York Yankees").
"""
from __future__ import annotations
import argparse
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import duckdb

from fd_client import FanDuelClient, Event, Market, Runner


# Regex helpers for FD's alt-market name formats.
# - Spread: 'St. Louis Cardinals +3.5' or 'Athletics -1.5'
# - Total:  'Over (8.5)' or 'Under (7.5)'
_FD_ALT_SPREAD_RE = re.compile(r"^(?P<team>.+?)\s+(?P<line>[+-]\d+(?:\.\d+)?)\s*$")
_FD_ALT_TOTAL_RE  = re.compile(r"^(?P<side>Over|Under)\s*\((?P<line>\d+(?:\.\d+)?)\)\s*$", re.IGNORECASE)


# Single-inning markets ("7th Inning Total Runs", "1st Inning Run Line") are
# out of scope — only cumulative "First N Innings" periods map to bet cards.
# Plural "Innings" (First 5 Innings) deliberately does NOT match this.
_SINGLE_INNING_RE = re.compile(r"\b\d+(st|nd|rd|th)\s+inning\b", re.IGNORECASE)

# Substrings that disqualify a market. Ported from DraftKings'
# classify_market exclusion list, plus FD-specific junk that keyword-collides
# with real game lines (FD posts ~150 markets/event):
#   - "parlay"  -> "First 5 Innings Run Line / Total Runs Parlay", "Line / Total Parlay N"
#   - "listed"  -> "Moneyline Away/Home/Both Listed" (pitcher-conditional; keep plain "Moneyline")
#   - "bands"   -> "Total Runs (Bands)"
#   - "tri-bet", "specials" -> FD novelty markets
_FD_EXCLUDE_KEYWORDS = (
    "team total", "player", "prop", "futures", "to record", "to score",
    "to hit", "first to", "race to", "correct score", "winning margin",
    "total bases", "rbis", "hits o/u", "strikeouts thrown", "odd/even",
    "score last", "bat bottom", "highest scoring", "most innings",
    "last run", "both teams to score",
    "parlay", "listed", "bands", "tri-bet", "specials",
)


def classify_market(
    name: str,
    home_team: str | None = None,
    away_team: str | None = None,
) -> tuple[str, str] | None:
    """Map an FD market name to (period, market_type), or None to skip.

    Keyword classifier (mirrors scraper_draftkings_singles.classify_market) so
    FD picks up every FG/F5/F7/F3 main + alt line it posts and auto-covers new
    ones. period is one of "FG", "F3", "F5", "F7".

    home_team / away_team are optional; when provided, any market name that
    contains a team name is treated as a per-team market (team total) and
    excluded. They default to None so existing single-arg callers/tests work.
    Game lines ("Run Line", "Total Runs", "Moneyline", "Alternate ...") never
    contain a team name, so this never false-excludes them.
    """
    n = name.lower()

    # Team totals: "<Team> Total Runs", "<Team> Alt. Total Runs".
    for team in (home_team, away_team):
        if team and team.lower() in n:
            return None

    if any(k in n for k in _FD_EXCLUDE_KEYWORDS):
        return None
    if _SINGLE_INNING_RE.search(n):
        return None

    # Period detection. Default FG; F3/F5/F7 if matched explicitly.
    if "first 7 innings" in n or "1st 7 innings" in n:
        period = "F7"
    elif "first 5 innings" in n or "1st 5 innings" in n:
        period = "F5"
    elif "first 3 innings" in n or "1st 3 innings" in n:
        period = "F3"
    else:
        period = "FG"

    # Market-type detection. Check "alternate" before main. Match both
    # "moneyline" (FD's FG name) and "money line" (FD's F-period name, e.g.
    # "First 5 Innings Money Line") so F-period MLs are not silently dropped.
    if "alternate" in n and "run line" in n:
        return (period, "alternate_spreads")
    if "alternate" in n and "total" in n:
        return (period, "alternate_totals")
    if "run line" in n:
        return (period, "main")
    if "moneyline" in n or "money line" in n:
        return (period, "main")
    if "total" in n:
        return (period, "main")
    return None


def parse_runners_to_wide_rows(
    event: Event,
    runners: list[Runner],
    market_meta: dict[str, tuple[str, str]],   # market_id -> (period, market_type)
    fetch_time: datetime,
) -> list[dict[str, Any]]:
    """Group runners by (period, market_type, line); emit wide rows.

    Logic mirrors parse_selections_to_wide_rows in scraper_draftkings_singles
    exactly — Runner and Selection have identical relevant attributes
    (name, line, american_odds) so the body is reusable verbatim.
    """
    buckets: dict[tuple[str, str, float | None], dict[str, Any]] = {}

    def _row_skeleton(period: str, market_type: str) -> dict:
        return {
            "fetch_time": fetch_time,
            "sport_key": "baseball_mlb",
            "game_id": event.event_id,
            "game_date": fetch_time.strftime("%Y-%m-%d"),
            "game_time": event.start_time,
            "away_team": event.away_team,
            "home_team": event.home_team,
            "market": market_type,
            "period": period,
            "away_spread": None, "away_spread_price": None,
            "home_spread": None, "home_spread_price": None,
            "total": None, "over_price": None, "under_price": None,
            "away_ml": None, "home_ml": None,
        }

    for r in runners:
        meta = market_meta.get(r.market_id)
        if meta is None:
            continue
        period, market_type = meta

        # FD-specific line extraction: for alt-spreads and alt-totals, FD often
        # encodes the line in the runner NAME instead of on `handicap`. We
        # resolve `effective_line` here so the rest of the loop can treat both
        # main and alt rows uniformly. (For main markets, FD always sets
        # handicap on the Runner, so r.line is non-None and no name-parse runs.)
        effective_line: float | None = r.line
        if r.line is None and market_type == "alternate_spreads":
            m = _FD_ALT_SPREAD_RE.match(r.name.strip())
            if m:
                try:
                    effective_line = float(m.group("line"))
                except ValueError:
                    effective_line = None
        elif r.line is None and market_type == "alternate_totals":
            m = _FD_ALT_TOTAL_RE.match(r.name.strip())
            if m:
                try:
                    effective_line = float(m.group("line"))
                except ValueError:
                    effective_line = None

        # Bucket key: main rows coalesce by period only; alt rows split by line.
        # For alt-spreads, both sides of a paired line (e.g. Yankees -2.5 and
        # Red Sox +2.5) share the same row, BUT the opposite-direction pair
        # (Yankees +2.5 and Red Sox -2.5) must be a DIFFERENT row. We bucket
        # by the home-team-signed line so both sides of a paired line collapse
        # to one key while the opposite-direction pair stays separate.
        if market_type == "main":
            bucket_line: float | None = None
        elif market_type == "alternate_spreads" and effective_line is not None:
            # Determine signed home-team line from the runner name. The home
            # team's line is its own when the runner is home; otherwise it's
            # the negation of the away-team runner's line.
            if event.home_team in r.name:
                bucket_line = effective_line
            elif event.away_team in r.name:
                bucket_line = -effective_line
            else:
                # Runner name doesn't match either team -- shouldn't happen,
                # but degrade gracefully: bucket by absolute line (legacy
                # behaviour) so this row at least doesn't crash.
                bucket_line = abs(effective_line)
        else:
            bucket_line = effective_line
        key = (period, market_type, bucket_line)
        if key not in buckets:
            buckets[key] = _row_skeleton(period, market_type)
        row = buckets[key]

        name_lower = r.name.lower()

        # Totals detection. Both FD main ("Over") and alt ("Over (8.5)") names
        # start with "over"/"under"; effective_line carries the correct line
        # in both cases.
        if name_lower.startswith("over"):
            row["total"] = effective_line
            row["over_price"] = r.american_odds
        elif name_lower.startswith("under"):
            row["total"] = effective_line
            row["under_price"] = r.american_odds
        # Spread detection: an effective line is present AND name contains a
        # team name. Note we check effective_line (not r.line) so alt-spread
        # runners with line embedded in the name still route here.
        elif effective_line is not None:
            if event.home_team in r.name:
                row["home_spread"] = effective_line
                row["home_spread_price"] = r.american_odds
            elif event.away_team in r.name:
                row["away_spread"] = effective_line
                row["away_spread_price"] = r.american_odds
            else:
                # Spread-shaped runner but doesn't match either team — skip
                continue
        # Moneyline detection: no line at all, just team name. This branch only
        # fires for genuine ML markets ("Moneyline" / "First 5 Innings Money
        # Line") because alt-spread runners with line=None+team-name in the
        # name had their effective_line populated by the regex above.
        else:
            if event.home_team in r.name:
                row["home_ml"] = r.american_odds
            elif event.away_team in r.name:
                row["away_ml"] = r.american_odds

    return list(buckets.values())


# FD's event-page returns different market slices per tab; no single tab has
# everything. The default tab carries F7/F3 + per-inning lines; the SGP tab
# carries FG-alts + all F5. Coverage = union of both. Verified 2026-05-23.
FD_TABS = ("", "same-game-parlay-")


def fetch_merged_markets_and_runners(
    client: FanDuelClient,
    event_id: str,
    tabs: tuple[str, ...] = FD_TABS,
) -> tuple[list[Market], list[Runner]]:
    """Fetch each tab once, union markets (dedup by market_id) and runners
    (dedup by runner_id). Returns (list[Market], list[Runner])."""
    markets_by_id = {}
    runners_by_id = {}
    for tab in tabs:
        markets, runners = client.fetch_event_page(event_id, tab)
        for m in markets:
            markets_by_id.setdefault(m.market_id, m)
        for r in runners:
            runners_by_id.setdefault(r.runner_id, r)
    return list(markets_by_id.values()), list(runners_by_id.values())


def scrape_singles(verbose: bool = False) -> int:
    """Scrape all MLB events from FD and atomically write singles to DuckDB.

    Per-game isolation: a single event's API failure does NOT tank the scrape.
    Returns the total number of rows written.
    """
    client = FanDuelClient(verbose=verbose)
    events = client.list_events()
    print(f"[fd_singles] {len(events)} events to scrape", flush=True)

    # Use timezone-aware UTC so DuckDB's TIMESTAMPTZ column receives an
    # explicit UTC instant. A naive datetime (datetime.utcnow()) would be
    # interpreted as local time and silently shifted by the host TZ offset
    # on insert — see TZ_AUDIT_FINDINGS for the class of bug this avoids.
    fetch_time = datetime.now(timezone.utc)
    all_rows: list[dict] = []
    failed: list[str] = []

    for event in events:
        try:
            markets, runners = fetch_merged_markets_and_runners(
                client, event.event_id)

            market_meta: dict[str, tuple[str, str]] = {}
            for m in markets:
                classified = classify_market(
                    m.name, event.home_team, event.away_team)
                if classified is not None:
                    market_meta[m.market_id] = classified

            rows = parse_runners_to_wide_rows(event, runners, market_meta, fetch_time)
            all_rows.extend(rows)
            if verbose:
                print(
                    f"  [{event.event_id}] {event.away_team} @ {event.home_team}: "
                    f"{len(rows)} rows ({len(market_meta)} in-scope markets, "
                    f"{len(runners)} runners total)",
                    flush=True,
                )
        except Exception as e:
            print(f"  [{event.event_id}] FAILED: {e}", flush=True)
            failed.append(event.event_id)
            continue

    write_to_duckdb(all_rows)
    print(
        f"[fd_singles] wrote {len(all_rows)} rows "
        f"({len(failed)} events failed)",
        flush=True,
    )
    return len(all_rows)


def write_to_duckdb(rows: list[dict]) -> None:
    """Atomic write: CREATE IF NOT EXISTS -> BEGIN -> DELETE -> INSERT -> COMMIT.

    The fd_odds/mlb_odds table is fully rewritten each cycle (offshore-style
    scraper convention — MLB.R is the consumer and snapshots into mlb.duckdb
    as needed).
    """
    db_path = Path(__file__).resolve().parent.parent / "fd_odds" / "fd.duckdb"
    db_path.parent.mkdir(exist_ok=True)

    con = duckdb.connect(str(db_path))
    try:
        con.execute(
            """
            CREATE TABLE IF NOT EXISTS mlb_odds (
                fetch_time        TIMESTAMP,
                sport_key         VARCHAR,
                game_id           VARCHAR,
                game_date         VARCHAR,
                game_time         VARCHAR,
                away_team         VARCHAR,
                home_team         VARCHAR,
                market            VARCHAR,
                period            VARCHAR,
                away_spread       FLOAT,
                away_spread_price INTEGER,
                home_spread       FLOAT,
                home_spread_price INTEGER,
                total             FLOAT,
                over_price        INTEGER,
                under_price       INTEGER,
                away_ml           INTEGER,
                home_ml           INTEGER
            )
            """
        )
        con.execute("BEGIN TRANSACTION")
        try:
            con.execute("DELETE FROM mlb_odds")
            if rows:
                cols = [
                    "fetch_time", "sport_key", "game_id", "game_date", "game_time",
                    "away_team", "home_team", "market", "period",
                    "away_spread", "away_spread_price",
                    "home_spread", "home_spread_price",
                    "total", "over_price", "under_price",
                    "away_ml", "home_ml",
                ]
                tuples = [tuple(r.get(c) for c in cols) for r in rows]
                placeholders = ", ".join(["?"] * len(cols))
                con.executemany(
                    f"INSERT INTO mlb_odds VALUES ({placeholders})", tuples
                )
            con.execute("COMMIT")
        except Exception:
            con.execute("ROLLBACK")
            raise
    finally:
        con.close()


def main() -> None:
    p = argparse.ArgumentParser(description="FanDuel MLB singles scraper")
    p.add_argument("sport", nargs="?", default="mlb",
                   help="Sport key (from run.py orchestrator). Only 'mlb' runs the scrape.")
    p.add_argument("--verbose", action="store_true", help="Per-event row logging")
    args = p.parse_args()
    if args.sport != "mlb":
        print(f"[fd_singles] sport={args.sport!r} not supported, exiting", flush=True)
        return
    scrape_singles(verbose=args.verbose)


if __name__ == "__main__":
    main()
