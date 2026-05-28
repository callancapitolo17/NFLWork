"""DraftKings single-leg odds scraper.

Walks every MLB event today, fetches all selections + market metadata via
DraftKingsClient, transforms to the offshore mlb_odds schema (wagerzon-style),
and writes to dk_odds/dk.duckdb.

Task 8 implements the parser. Task 9 wires the scrape loop, classify_market
function, DK team-name canonicalization, and DuckDB write.
"""
from __future__ import annotations
import argparse
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import duckdb

from dk_client import DraftKingsClient, Event, Selection


# DK team-name canonicalization. DK uses abbreviated-city prefixes
# (e.g. "CLE Guardians") that don't match canonical Odds-API names used by
# mlb_consensus_temp. Built from a live DK scrape on 2026-05-12.
# "Athletics" is identical on both sides (DK does not prefix it).
DK_TEAM_MAP: dict[str, str] = {
    "ARI Diamondbacks": "Arizona Diamondbacks",
    "ATL Braves": "Atlanta Braves",
    "Athletics": "Athletics",
    "BAL Orioles": "Baltimore Orioles",
    "BOS Red Sox": "Boston Red Sox",
    "CHI Cubs": "Chicago Cubs",
    "CHI White Sox": "Chicago White Sox",
    "CIN Reds": "Cincinnati Reds",
    "CLE Guardians": "Cleveland Guardians",
    "COL Rockies": "Colorado Rockies",
    "DET Tigers": "Detroit Tigers",
    "HOU Astros": "Houston Astros",
    "KC Royals": "Kansas City Royals",
    "LA Angels": "Los Angeles Angels",
    "LA Dodgers": "Los Angeles Dodgers",
    "MIA Marlins": "Miami Marlins",
    "MIL Brewers": "Milwaukee Brewers",
    "MIN Twins": "Minnesota Twins",
    "NY Mets": "New York Mets",
    "NY Yankees": "New York Yankees",
    "PHI Phillies": "Philadelphia Phillies",
    "PIT Pirates": "Pittsburgh Pirates",
    "SD Padres": "San Diego Padres",
    "SEA Mariners": "Seattle Mariners",
    "SF Giants": "San Francisco Giants",
    "STL Cardinals": "St. Louis Cardinals",
    "TB Rays": "Tampa Bay Rays",
    "TEX Rangers": "Texas Rangers",
    "TOR Blue Jays": "Toronto Blue Jays",
    "WAS Nationals": "Washington Nationals",
}


def canonicalize_team(dk_name: str) -> str:
    """Translate DK team name to canonical (Odds-API style). Falls back to original on miss."""
    return DK_TEAM_MAP.get(dk_name, dk_name)


# Match single-inning markets DK posts (4th Inning, 5th Inning, ... 9th Inning,
# plus "(3 Way)" variants). These are NOT in scope — only "1st N Innings"
# cumulative periods are. We never accept "Nth Inning" (singular ordinal suffix).
_SINGLE_INNING_RE = re.compile(r"\b\d+(st|nd|rd|th)\s+inning\b", re.IGNORECASE)


def classify_market(name: str) -> tuple[str, str] | None:
    """Map DK market name to (period, market_type), or None to skip.

    Returns:
      - (period, "main") for main spread/total/ML markets
      - (period, "alternate_spreads") for alt spread markets
      - (period, "alternate_totals") for alt total markets
      - None for markets out of scope (props, team totals, single-inning, futures)

    period is one of "FG", "F3", "F5", "F7".
    """
    n = name.lower()

    # Out-of-scope categories first.
    if "team total" in n:
        return None
    # Some DK per-team markets don't carry the "team total" substring
    # (e.g. "Alternate HOU Astros Total Runs", "CHI Cubs Home Runs"). These
    # used to slip through and collide with the real game-totals market in
    # parse_selections_to_wide_rows's (period, market_type, line) bucket,
    # silently overwriting game-total odds with one team's team-total odds.
    # Filter any market name containing a DK team prefix as team-specific.
    if any(team.lower() in n for team in DK_TEAM_MAP):
        return None
    if any(k in n for k in (
        "player", "prop", "futures", "to record", "to score",
        "to hit", "first to", "race to", "correct score", "winning margin",
        "total bases", "rbis", "hits o/u", "strikeouts thrown",
        "odd/even", "score last", "bat bottom", "highest scoring",
        "most innings", "last run", "both teams to score",
    )):
        return None
    # Single-inning markets (e.g. "Run Line - 5th Inning", "Total Runs - 6th
    # Inning", "7th Inning (3 Way)") — exclude. Note "1st N Innings" stays
    # because it's plural "Innings", not singular "Inning".
    if _SINGLE_INNING_RE.search(n):
        return None

    # Period detection. Default FG; F3/F5/F7 if matched explicitly.
    if "1st 3 innings" in n or "first 3 innings" in n:
        period = "F3"
    elif "1st 5 innings" in n or "first 5 innings" in n:
        period = "F5"
    elif "1st 7 innings" in n or "first 7 innings" in n:
        period = "F7"
    else:
        period = "FG"

    # Market-type detection. Order matters: check "alternate" before main
    # because "Total Alternate" contains both "alternate" and "total".
    if "alternate" in n and "run line" in n:
        return (period, "alternate_spreads")
    if "alternate" in n and "total" in n:
        return (period, "alternate_totals")
    if "run line" in n:
        return (period, "main")
    if "moneyline" in n:
        return (period, "main")
    if "total" in n:
        return (period, "main")
    # Bare period name "1st 3 Innings" / "1st 5 Innings" / "1st 7 Innings" is
    # DK's 2-way period winner (no tie on first-N-innings). Its selections have
    # no line, so they flow through the moneyline branch as home_ml/away_ml.
    if n in ("1st 3 innings", "1st 5 innings", "1st 7 innings",
             "first 3 innings", "first 5 innings", "first 7 innings"):
        return (period, "main")
    return None


def parse_selections_to_wide_rows(
    event: Event,
    selections: list[Selection],
    market_meta: dict[str, tuple[str, str]],   # market_id -> (period, market_type)
    fetch_time: datetime,
) -> list[dict[str, Any]]:
    """Group selections by (period, market_type, line); emit wide rows.

    DK's F3/F5/F7 run-line and total markets bundle the main line and all alt
    lines into ONE market (e.g. "Total Runs - 1st 7 Innings" carries Over/
    Under at 5.5, 6.5, AND 7.5 as 6 selections in one market id). Without
    detecting that, the parser would coalesce every selection into the same
    (period, "main", None) bucket and lose all but the last line. The
    pre-pass below counts distinct lines per market and any "main"-classified
    market carrying >1 distinct line is reclassified per selection as
    alternate_* so each line becomes its own DB row. FG markets are unaffected
    (DK splits FG main from FG alt into separate markets).
    """
    # Pre-pass: count distinct totals lines and distinct spread |lines| per
    # market_id. Used below to detect bundled multi-line "main" markets.
    totals_lines_by_market: dict[str, set[float]] = {}
    spread_lines_by_market: dict[str, set[float]] = {}
    for sel in selections:
        if sel.market_id not in market_meta:
            continue
        if sel.line is None:
            continue
        if sel.name.lower().startswith(("over", "under")):
            totals_lines_by_market.setdefault(sel.market_id, set()).add(sel.line)
        else:
            spread_lines_by_market.setdefault(sel.market_id, set()).add(abs(sel.line))

    buckets: dict[tuple[str, str, float | None], dict[str, Any]] = {}

    def _row_skeleton(period: str, market_type: str) -> dict:
        return {
            "fetch_time": fetch_time,
            "sport_key": "baseball_mlb",
            "game_id": event.event_id,
            # event.start_time is an ISO 8601 UTC string straight from DK's
            # startEventDate field (see dk_client.Event docstring). DuckDB
            # coerces ISO 8601 strings to TIMESTAMPTZ on INSERT — the value
            # is the absolute UTC instant; display tz is just rendering.
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

    for sel in selections:
        meta = market_meta.get(sel.market_id)
        if meta is None:
            continue
        period, market_type = meta
        name_lower = sel.name.lower()

        # If this is a "main"-classified market whose selections bundle multiple
        # distinct lines (DK F3/F5/F7 quirk), reclassify per-selection as alt.
        effective_market_type = market_type
        if market_type == "main":
            if name_lower.startswith(("over", "under")):
                if len(totals_lines_by_market.get(sel.market_id, ())) > 1:
                    effective_market_type = "alternate_totals"
            elif sel.line is not None:
                if len(spread_lines_by_market.get(sel.market_id, ())) > 1:
                    effective_market_type = "alternate_spreads"

        # Bucket key: main rows coalesce by period only; alt-spread rows split
        # by the SIGNED home-equivalent line. DK posts the two opposite
        # directions at one magnitude as separate two-way markets — e.g.
        # ARI -2.5/SFG +2.5 AND ARI +2.5/SFG -2.5. Bucketing by abs(line)
        # collapsed both into one row, so last-write-wins silently dropped a
        # whole direction (the favorite-laying-runs ladder vanished, leaving
        # only a lone deep line whose mirror DK didn't post). Keying on the
        # signed home line keeps each side of one market together while
        # separating the two directions. Alt-totals: Over/Under already share
        # one line value, so key on it directly.
        if effective_market_type == "main":
            bucket_line: float | None = None
        elif effective_market_type == "alternate_spreads" and sel.line is not None:
            sel_name = sel.name.strip()
            if (sel_name.startswith(event.home_team + " ")
                    or sel_name == event.home_team):
                bucket_line = sel.line           # home line IS the key
            elif (sel_name.startswith(event.away_team + " ")
                    or sel_name == event.away_team):
                bucket_line = -sel.line          # away line -> implied home line
            else:
                # Unrecognized team on a spread selection — skip rather than
                # create a malformed half-row (the finalization guard would
                # drop it anyway).
                continue
        else:
            bucket_line = sel.line
        key = (period, effective_market_type, bucket_line)
        if key not in buckets:
            buckets[key] = _row_skeleton(period, effective_market_type)
        row = buckets[key]

        # Totals detection: name starts with "over" or "under"
        if name_lower.startswith("over"):
            row["total"] = sel.line
            row["over_price"] = sel.american_odds
        elif name_lower.startswith("under"):
            row["total"] = sel.line
            row["under_price"] = sel.american_odds
        # Spread detection: line present AND name matches a team name. We
        # use exact-prefix matching (with a separator space or full equality)
        # rather than substring containment — DK names like "LA Angels" can
        # appear inside "Los Angeles Angels" headers and substring tests
        # silently mis-bucket. Prefix-with-space anchors to the team-name
        # token, leaving the spread value to follow ("LA Angels +1.5").
        elif sel.line is not None:
            sel_name = sel.name.strip()
            if (sel_name.startswith(event.home_team + " ")
                    or sel_name == event.home_team):
                row["home_spread"] = sel.line
                row["home_spread_price"] = sel.american_odds
            elif (sel_name.startswith(event.away_team + " ")
                    or sel_name == event.away_team):
                row["away_spread"] = sel.line
                row["away_spread_price"] = sel.american_odds
            else:
                # Unrecognized — skip, do not silently mis-bucket.
                continue
        # Moneyline detection: no line, just team name. Same prefix anchor.
        else:
            sel_name = sel.name.strip()
            if (sel_name.startswith(event.home_team + " ")
                    or sel_name == event.home_team):
                row["home_ml"] = sel.american_odds
            elif (sel_name.startswith(event.away_team + " ")
                    or sel_name == event.away_team):
                row["away_ml"] = sel.american_odds

    # Finalization pass over spread rows:
    #  (a) Paired-side guard — if a spread row has only one side resolved,
    #      drop it. Tools.R::get_dk_odds reads home_spread as canonical and
    #      derives away from -home_spread; a half-row would either silently
    #      drop both sides downstream or pair a real line with a NULL price.
    #  (b) Sign-symmetry guard — both sides must sum to ~0. An asymmetric
    #      pair (e.g. home=-1.5 away=+2.5) indicates a parser misgrouping
    #      or DK posting mismatched alt lines into the same market bucket.
    out: list[dict[str, Any]] = []
    for row in buckets.values():
        is_spread_row = row["market"] in ("main", "alternate_spreads") and (
            row.get("home_spread") is not None or row.get("away_spread") is not None
        )
        if is_spread_row:
            if row.get("home_spread") is None or row.get("away_spread") is None:
                # Half-row — drop silently. Asymmetric posting is common
                # enough on F5/F7 alts that warning would be noise.
                continue
            if abs(row["home_spread"] + row["away_spread"]) > 1e-9:
                print(
                    f"[dk_singles] WARN: asymmetric spread "
                    f"home={row['home_spread']} away={row['away_spread']} "
                    f"on event={event.event_id}; skipping",
                    flush=True,
                )
                continue
        out.append(row)
    return out


def scrape_singles(verbose: bool = False) -> int:
    """Scrape all MLB events from DK and atomically write singles to DuckDB.

    Per-game isolation: a single event's API failure does NOT tank the scrape.
    Returns the total number of rows written.
    """
    client = DraftKingsClient(verbose=verbose)
    events = client.list_events()
    print(f"[dk_singles] {len(events)} events to scrape", flush=True)

    # Use timezone-aware UTC so DuckDB's TIMESTAMPTZ column receives an
    # explicit UTC instant. A naive datetime would be interpreted as local
    # time and silently shifted on insert (see TZ_AUDIT_FINDINGS for the
    # class of bug this avoids).
    fetch_time = datetime.now(timezone.utc)
    all_rows: list[dict] = []
    failed: list[str] = []
    unmapped_teams: set[str] = set()

    # Late import to avoid circular dep at module import time.
    from scraper_draftkings_sgp import DK_SGP_PARLAYS_URL

    for event in events:
        # Track DK names that aren't in the map so we can surface them in logs.
        for raw in (event.home_team, event.away_team):
            if raw not in DK_TEAM_MAP:
                unmapped_teams.add(raw)

        # Canonicalize team names BEFORE constructing the parser's Event. The
        # parser uses event.home_team / event.away_team for substring matching
        # against selection names — DK's selection names use the SAME
        # abbreviated DK form (e.g. "CLE Guardians +1.5"), so the parser
        # substring-match must run against DK names, NOT canonical names.
        # Solution: parser sees DK names; we swap to canonical when stamping
        # rows. To do that without changing the parser, we pass the DK event
        # to the parser and then re-stamp home_team / away_team on each row.
        try:
            markets = client.fetch_event_markets(event.event_id)
            selections = client.fetch_event_selections(event.event_id)

            market_meta: dict[str, tuple[str, str]] = {}
            for m in markets:
                classified = classify_market(m.name)
                if classified is not None:
                    market_meta[m.market_id] = classified

            # fetch_event_markets only covers subcats 4519 + 15628. The parlays
            # endpoint exposes more markets (e.g. F7 totals/alts) that we still
            # want. Pull market names directly from the parlays payload and
            # enrich market_meta for any IDs not already classified.
            try:
                r = client.session.get(
                    f"{DK_SGP_PARLAYS_URL}/{event.event_id}", timeout=60
                )
                mkts = ((r.json() or {}).get("data") or {}).get("markets") or []
                for m in mkts:
                    mid = str(m.get("id", ""))
                    if mid and mid not in market_meta:
                        cls = classify_market(m.get("name", "") or "")
                        if cls is not None:
                            market_meta[mid] = cls
            except Exception as e:
                if verbose:
                    print(f"  [{event.event_id}] parlays-enrich failed: {e}", flush=True)

            rows = parse_selections_to_wide_rows(event, selections, market_meta, fetch_time)
            # Re-stamp DK team names to canonical Odds-API names.
            canonical_home = canonicalize_team(event.home_team)
            canonical_away = canonicalize_team(event.away_team)
            for row in rows:
                row["home_team"] = canonical_home
                row["away_team"] = canonical_away
            all_rows.extend(rows)
            if verbose:
                print(f"  [{event.event_id}] {event.away_team} @ {event.home_team}: {len(rows)} rows", flush=True)
        except Exception as e:
            print(f"  [{event.event_id}] FAILED: {e}", flush=True)
            failed.append(event.event_id)
            continue

    write_to_duckdb(all_rows)
    print(
        f"[dk_singles] wrote {len(all_rows)} rows "
        f"({len(failed)} events failed)",
        flush=True,
    )
    if unmapped_teams:
        print(
            f"[dk_singles] WARNING — DK team names missing from DK_TEAM_MAP: "
            f"{sorted(unmapped_teams)}",
            flush=True,
        )
    return len(all_rows)


def write_to_duckdb(rows: list[dict]) -> None:
    """Atomic write: stage to TEMP table, then CREATE OR REPLACE the live one.

    The dk_odds/mlb_odds table is fully rewritten each cycle (no history kept
    inside this DB — that's the offshore-style scraper convention; MLB.R is
    the consumer and snapshots into mlb.duckdb as needed).

    Empty-scrape guard: if no rows were produced this cycle, we DO NOT touch
    the existing table. An empty scrape is almost always a transient failure
    (network blip, DK quiet period); blowing away the prior snapshot would
    leave downstream MLB.R consumers with no DK pills until the next cycle.
    Leaving the old snapshot in place is the safer default — staleness is
    visible via fetch_time, but absence is invisible.
    """
    db_path = Path(__file__).resolve().parent.parent / "dk_odds" / "dk.duckdb"
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
            print(f"[dk] Migrating mlb_odds.fetch_time TIMESTAMP -> TIMESTAMPTZ "
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
            # see either the entire old snapshot or the entire new one,
            # never a half-written state.
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
                "[dk_singles] empty scrape — leaving prior snapshot in place",
                flush=True,
            )
    finally:
        con.close()


def main() -> None:
    p = argparse.ArgumentParser(description="DraftKings MLB singles scraper")
    p.add_argument("sport", nargs="?", default="mlb",
                   help="Sport key (from run.py orchestrator). Only 'mlb' runs the scrape.")
    p.add_argument("--verbose", action="store_true", help="Per-event row logging")
    args = p.parse_args()
    if args.sport != "mlb":
        print(f"[dk_singles] sport={args.sport!r} not supported, exiting", flush=True)
        return
    scrape_singles(verbose=args.verbose)


if __name__ == "__main__":
    main()
