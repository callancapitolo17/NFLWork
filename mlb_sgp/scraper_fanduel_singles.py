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
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import duckdb

from fd_client import FanDuelClient, Event, Runner


# Regex helpers for FD's alt-market name formats.
# - Spread: 'St. Louis Cardinals +3.5' or 'Athletics -1.5'
# - Total:  'Over (8.5)' or 'Under (7.5)'
_FD_ALT_SPREAD_RE = re.compile(r"^(?P<team>.+?)\s+(?P<line>[+-]\d+(?:\.\d+)?)\s*$")
_FD_ALT_TOTAL_RE  = re.compile(r"^(?P<side>Over|Under)\s*\((?P<line>\d+(?:\.\d+)?)\)\s*$", re.IGNORECASE)


# Unicode-minus variants FD has been observed using interchangeably with the
# ASCII '-' on alt-spread runners. The alt-spread regex above only matches
# `[+-]`, so any name with U+2212 (true minus), U+2013 (en-dash), or U+2014
# (em-dash) would silently fail to parse and the row would degrade. The
# normalizer below rewrites them to ASCII '-' and collapses whitespace runs
# so the regex sees a stable shape.
_UNICODE_MINUS = "−–—"  # U+2212, en-dash, em-dash


def _normalize_alt_name(raw: str) -> str:
    """Normalize unicode minus to ASCII minus, collapse whitespace."""
    out = raw
    for ch in _UNICODE_MINUS:
        out = out.replace(ch, "-")
    return " ".join(out.split())


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
    # F3 main — BLIND ADD: 2026-05-22 live probe showed FD posts ONLY
    # "First 3 Innings Result" (3-way Home/Away/Tie ML, not the 2-way ML
    # we want) and no F3 spread / total / alt markets. Names below follow
    # FD's F5 convention so they activate the day FD adds parity; today
    # they match nothing and produce zero F3 rows. Mirror of `fd_no_f7_
    # line_markets` situation, just at the F3 boundary.
    "First 3 Innings Run Line":              ("F3", "main"),
    "First 3 Innings Total Runs":            ("F3", "main"),
    "First 3 Innings Money Line":            ("F3", "main"),
    # F3 alt — same caveat as above
    "First 3 Innings Alternate Run Lines":   ("F3", "alternate_spreads"),
    "First 3 Innings Alternate Total Runs":  ("F3", "alternate_totals"),
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
            # event.start_time is FD's ISO 8601 UTC string (open_date field,
            # e.g. "2026-05-23T01:41:00.000Z"). DuckDB coerces ISO 8601
            # strings to TIMESTAMPTZ on INSERT — the stored value is the
            # absolute UTC instant; display tz is just rendering. Mirrors
            # the DK scraper pattern (scraper_draftkings_singles.py).
            "game_start_time": event.start_time,
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
        #
        # We pre-normalize the name (unicode minus -> ASCII minus, whitespace
        # collapse) so a stray U+2212 / en-dash / em-dash doesn't silently
        # bypass the regex and degrade the row to a ML-shaped "no line" path.
        effective_line: float | None = r.line
        if r.line is None and market_type == "alternate_spreads":
            m = _FD_ALT_SPREAD_RE.match(_normalize_alt_name(r.name))
            if m:
                try:
                    effective_line = float(m.group("line"))
                except ValueError:
                    effective_line = None
        elif r.line is None and market_type == "alternate_totals":
            m = _FD_ALT_TOTAL_RE.match(_normalize_alt_name(r.name))
            if m:
                try:
                    effective_line = float(m.group("line"))
                except ValueError:
                    effective_line = None

        # Skip-on-parse-fail guard for alt rows. If we couldn't recover a line
        # for an alt-spread or alt-total runner, we MUST NOT fall through to
        # the ML branch below — that would stamp a price onto home_ml/away_ml
        # on an alt-classified row, silently corrupting the bucket. Warn so
        # silent FD format drift is visible in logs.
        if effective_line is None and market_type in (
            "alternate_spreads", "alternate_totals"
        ):
            print(
                f"[fd_singles] WARN: alt-line parse failed for "
                f"runner={r.name!r} market_type={market_type}",
                flush=True,
            )
            continue

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

    # Finalization pass — mirrors the DK scraper:
    #  (a) Paired-side guard: a spread row with only one side resolved is
    #      dropped. Downstream get_fd_odds reads home_spread as canonical
    #      and derives away from -home_spread; a half-row would either
    #      silently drop both sides or pair a real line with a NULL price.
    #  (b) Sign-symmetry guard: the two sides must sum to ~0. An asymmetric
    #      pair (e.g. home=-1.5 away=+2.5) signals parser misgrouping or FD
    #      posting mismatched alt lines into the same bucket — log + skip.
    out: list[dict[str, Any]] = []
    for row in buckets.values():
        is_spread_row = row["market"] in ("main", "alternate_spreads") and (
            row.get("home_spread") is not None or row.get("away_spread") is not None
        )
        if is_spread_row:
            if row.get("home_spread") is None or row.get("away_spread") is None:
                # Half-row — drop silently. Asymmetric posting is common
                # enough on F5 alts that warning would be noise.
                continue
            if abs(row["home_spread"] + row["away_spread"]) > 1e-9:
                print(
                    f"[fd_singles] WARN: asymmetric spread "
                    f"home={row['home_spread']} away={row['away_spread']} "
                    f"on event={event.event_id}; skipping",
                    flush=True,
                )
                continue
        out.append(row)
    return out


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
    """Atomic write: stage to TEMP table, then CREATE OR REPLACE the live one.

    The fd_odds/mlb_odds table is fully rewritten each cycle (no history kept
    inside this DB — offshore-style scraper convention; MLB.R is the consumer
    and snapshots into mlb.duckdb as needed).

    Empty-scrape guard: if no rows were produced this cycle, we DO NOT touch
    the existing table. An empty scrape is almost always a transient failure
    (network blip, FD quiet period); blowing away the prior snapshot would
    leave downstream MLB.R consumers with no FD pills until the next cycle.
    Leaving the old snapshot in place is the safer default — staleness is
    visible via fetch_time, but absence is invisible. Mirrors the DK pattern.
    """
    db_path = Path(__file__).resolve().parent.parent / "fd_odds" / "fd.duckdb"
    db_path.parent.mkdir(exist_ok=True)

    con = duckdb.connect(str(db_path))
    try:
        con.execute(
            """
            CREATE TABLE IF NOT EXISTS mlb_odds (
                fetch_time        TIMESTAMPTZ,
                sport_key         VARCHAR,
                game_id           VARCHAR,
                game_start_time   TIMESTAMPTZ,
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
        if rows:
            cols = [
                "fetch_time", "sport_key", "game_id", "game_start_time",
                "away_team", "home_team", "market", "period",
                "away_spread", "away_spread_price",
                "home_spread", "home_spread_price",
                "total", "over_price", "under_price",
                "away_ml", "home_ml",
            ]
            # Stage into a TEMP table cloned from the live schema, then
            # atomically swap. The CREATE OR REPLACE TABLE ... AS SELECT step
            # is the closest DuckDB equivalent of an atomic rename — readers
            # see either the entire old snapshot or the entire new one, never
            # a half-written state.
            con.execute(
                "CREATE OR REPLACE TEMP TABLE mlb_odds_new AS "
                "SELECT * FROM mlb_odds LIMIT 0"
            )
            tuples = [tuple(r.get(c) for c in cols) for r in rows]
            placeholders = ", ".join(["?"] * len(cols))
            con.executemany(
                f"INSERT INTO mlb_odds_new VALUES ({placeholders})", tuples
            )
            con.execute(
                "CREATE OR REPLACE TABLE mlb_odds AS SELECT * FROM mlb_odds_new"
            )
            con.execute("DROP TABLE IF EXISTS mlb_odds_new")
        else:
            print(
                "[fd_singles] empty scrape — leaving prior snapshot in place",
                flush=True,
            )
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
