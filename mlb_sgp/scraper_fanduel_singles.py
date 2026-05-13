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
  - classify_market is a strict whitelist against `_FD_MARKET_WHITELIST` —
    FD posts ~150 markets per event and naive keyword matching catches dozens
    of bogus "Line / Total Parlay" and team-total markets.
  - No team-name canonicalization map — Phase 1 confirmed FD team names
    already match the canonical Odds-API forms (e.g. "St. Louis Cardinals",
    "Athletics", "New York Yankees").
"""
from __future__ import annotations
import argparse
import re
from datetime import datetime
from pathlib import Path
from typing import Any

import duckdb

from fd_client import FanDuelClient, Event, Runner


# Regex helpers for FD's alt-market name formats.
# - Spread: 'St. Louis Cardinals +3.5' or 'Athletics -1.5'
# - Total:  'Over (8.5)' or 'Under (7.5)'
_FD_ALT_SPREAD_RE = re.compile(r"^(?P<team>.+?)\s+(?P<line>[+-]\d+(?:\.\d+)?)\s*$")
_FD_ALT_TOTAL_RE  = re.compile(r"^(?P<side>Over|Under)\s*\((?P<line>\d+(?:\.\d+)?)\)\s*$", re.IGNORECASE)


# FD market name whitelist — these are the ONLY market names we accept.
# Whitelist (vs DK's keyword blacklist) because FD posts ~150 markets per
# event and a huge fraction match the obvious keywords ("Run Line", "Total
# Runs", "Moneyline") but are NOT the markets we want:
#   - "Run Line / Total Runs Parlay" (2-leg parlay, not a single bet)
#   - "Money Line / Total Runs Parlay"
#   - "Line / Total Parlay 1..12"
#   - "Moneyline Away Listed" / "Moneyline Home Listed" / "Moneyline Both Listed"
#     (pitcher-conditional ML variants — ignore; the plain "Moneyline" is the
#     'action regardless of pitcher' version we use everywhere else)
#   - "Total Runs (Bands)" (3+ way bands market, not over/under)
#   - "<Team> Total Runs", "<Team> Alt. Total Runs" (team totals — out of scope)
# Trying to blacklist all these is fragile; the canonical names are stable
# and FD has never added new core line markets.
_FD_MARKET_WHITELIST: dict[str, tuple[str, str]] = {
    # FG main
    "Run Line":                              ("FG", "main"),
    "Total Runs":                            ("FG", "main"),
    "Moneyline":                             ("FG", "main"),
    # FG alt
    "Alternate Run Lines":                   ("FG", "alternate_spreads"),
    "Alternate Total Runs":                  ("FG", "alternate_totals"),
    # F5 main
    "First 5 Innings Run Line":              ("F5", "main"),
    "First 5 Innings Total Runs":            ("F5", "main"),
    "First 5 Innings Money Line":            ("F5", "main"),
    # F5 alt
    "First 5 Innings Alternate Run Lines":   ("F5", "alternate_spreads"),
    "First 5 Innings Alternate Total Runs":  ("F5", "alternate_totals"),
}


def classify_market(name: str) -> tuple[str, str] | None:
    """Map FD market name to (period, market_type), or None to skip.

    Strict whitelist against `_FD_MARKET_WHITELIST`. Returns the same shape
    as DK's classify_market — (period, market_type) — so the parser is
    agnostic to which book the rows came from.
    """
    return _FD_MARKET_WHITELIST.get(name)


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
        # For alt-spreads, both sides (e.g. Yankees -2.5 and Red Sox +2.5) share
        # the same row, so bucket by absolute line. For alt-totals, Over/Under
        # already share the same line value.
        if market_type == "main":
            bucket_line: float | None = None
        elif market_type == "alternate_spreads" and effective_line is not None:
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


def scrape_singles(verbose: bool = False) -> int:
    """Scrape all MLB events from FD and atomically write singles to DuckDB.

    Per-game isolation: a single event's API failure does NOT tank the scrape.
    Returns the total number of rows written.
    """
    client = FanDuelClient(verbose=verbose)
    events = client.list_events()
    print(f"[fd_singles] {len(events)} events to scrape", flush=True)

    fetch_time = datetime.utcnow()
    all_rows: list[dict] = []
    failed: list[str] = []

    for event in events:
        try:
            markets = client.fetch_event_markets(event.event_id)
            runners = client.fetch_event_runners(event.event_id)

            market_meta: dict[str, tuple[str, str]] = {}
            for m in markets:
                classified = classify_market(m.name)
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
    p.add_argument("--verbose", action="store_true", help="Per-event row logging")
    args = p.parse_args()
    scrape_singles(verbose=args.verbose)


if __name__ == "__main__":
    main()
