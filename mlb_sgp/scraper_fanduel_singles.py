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
  - The scrape fetches BOTH FD event-page tabs (default + same-game-parlay-)
    and merges them — no single tab returns every market (F7/F3 lines live
    only on the default tab; FG-alts + all F5 live only on the SGP tab).
  - classify_market is a keyword classifier (mirrors DK) with FD-specific junk
    exclusions ("Line / Total Parlay", "Total Runs (Bands)", "... Listed",
    per-inning, 3-way "Result") and event-team team-total filtering — FD posts
    ~150 markets per event, so naive matching catches dozens of bogus markets.
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
#   - "result" -> 3-way "First N Innings Result" (Home/Tie/Away ML) — we only
#     take 2-way lines; load-bearing so a rename like "Result Total" can't slip
#     through the keyword fall-through. No 2-way line market name contains it.
_FD_EXCLUDE_KEYWORDS = (
    "team total", "player", "prop", "futures", "to record", "to score",
    "to hit", "first to", "race to", "correct score", "winning margin",
    "total bases", "rbis", "hits o/u", "strikeouts thrown", "odd/even",
    "score last", "bat bottom", "highest scoring", "most innings",
    "last run", "both teams to score",
    "parlay", "listed", "bands", "tri-bet", "specials", "result",
)


def classify_market(
    name: str,
    home_team: str | None = None,
    away_team: str | None = None,
) -> tuple[str, str] | None:
    """Map an FD market name to (period, market_type), or None to skip.

    Keyword classifier (mirrors scraper_draftkings_singles.classify_market) so
    FD picks up every FG/F5/F7/F3 main + alt line it posts and auto-covers new
    ones — replacing the old exact-name whitelist that silently missed markets
    FD posts under names not in the list. period is one of "FG", "F3", "F5", "F7".

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
    markets_by_id: dict[str, Market] = {}
    runners_by_id: dict[str, Runner] = {}
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
        # Migrate naive-TIMESTAMP schema to TIMESTAMPTZ if needed.
        # DuckDB does not support ALTER COLUMN TYPE between these, so drop+create.
        existing = con.execute(
            "SELECT column_name, data_type FROM information_schema.columns "
            "WHERE table_name = 'mlb_odds' AND column_name = 'fetch_time'"
        ).fetchone()
        if existing is not None and "WITH TIME ZONE" not in (existing[1] or "").upper():
            print(f"[fd] Migrating mlb_odds.fetch_time TIMESTAMP -> TIMESTAMPTZ "
                  f"(existing snapshot will be re-populated this run)")
            con.execute("DROP TABLE mlb_odds")

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
