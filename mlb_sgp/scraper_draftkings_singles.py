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
from datetime import datetime
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

    period is one of "FG", "F5", "F7".
    """
    n = name.lower()

    # Out-of-scope categories first.
    if "team total" in n:
        return None
    if any(k in n for k in (
        "player", "prop", "futures", "to record", "to score",
        "to hit", "first to", "race to", "correct score", "winning margin",
        "total bases", "rbis", "hits o/u", "strikeouts thrown",
        "odd/even", "score last", "bat bottom", "highest scoring",
        "most innings", "last run", "both teams to score",
    )):
        return None
    # Per plan: skip F3 markets entirely (DK posts them, but they're not in scope).
    if "1st 3 innings" in n or "first 3 innings" in n:
        return None
    # Single-inning markets (e.g. "Run Line - 5th Inning", "Total Runs - 6th
    # Inning", "7th Inning (3 Way)") — exclude. Note "1st 5 Innings" stays
    # because it's plural "Innings", not singular "Inning".
    if _SINGLE_INNING_RE.search(n):
        return None

    # Period detection. Default FG; F5/F7 if matched explicitly.
    if "1st 5 innings" in n or "first 5 innings" in n:
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
    return None


def parse_selections_to_wide_rows(
    event: Event,
    selections: list[Selection],
    market_meta: dict[str, tuple[str, str]],   # market_id -> (period, market_type)
    fetch_time: datetime,
) -> list[dict[str, Any]]:
    """Group selections by (period, market_type, line); emit wide rows."""
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

    for sel in selections:
        meta = market_meta.get(sel.market_id)
        if meta is None:
            continue
        period, market_type = meta

        # Bucket key: main rows coalesce by period only; alt rows split by line.
        # For alt-spreads, both sides (e.g. Yankees -2.5 and Red Sox +2.5) share
        # the same row, so bucket by absolute line. For alt-totals, Over/Under
        # already share the same line value.
        if market_type == "main":
            bucket_line: float | None = None
        elif market_type == "alternate_spreads" and sel.line is not None:
            bucket_line = abs(sel.line)
        else:
            bucket_line = sel.line
        key = (period, market_type, bucket_line)
        if key not in buckets:
            buckets[key] = _row_skeleton(period, market_type)
        row = buckets[key]

        name_lower = sel.name.lower()

        # Totals detection: name starts with "over" or "under"
        if name_lower.startswith("over"):
            row["total"] = sel.line
            row["over_price"] = sel.american_odds
        elif name_lower.startswith("under"):
            row["total"] = sel.line
            row["under_price"] = sel.american_odds
        # Spread detection: line present AND name contains a team name
        elif sel.line is not None:
            if event.home_team in sel.name:
                row["home_spread"] = sel.line
                row["home_spread_price"] = sel.american_odds
            elif event.away_team in sel.name:
                row["away_spread"] = sel.line
                row["away_spread_price"] = sel.american_odds
            else:
                # Spread-shaped selection but doesn't match either team — skip
                continue
        # Moneyline detection: no line, just team name
        else:
            if event.home_team in sel.name:
                row["home_ml"] = sel.american_odds
            elif event.away_team in sel.name:
                row["away_ml"] = sel.american_odds

    return list(buckets.values())


def scrape_singles(verbose: bool = False) -> int:
    """Scrape all MLB events from DK and atomically write singles to DuckDB.

    Per-game isolation: a single event's API failure does NOT tank the scrape.
    Returns the total number of rows written.
    """
    client = DraftKingsClient(verbose=verbose)
    events = client.list_events()
    print(f"[dk_singles] {len(events)} events to scrape", flush=True)

    fetch_time = datetime.utcnow()
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
    """Atomic write: CREATE IF NOT EXISTS -> BEGIN -> DELETE -> INSERT -> COMMIT.

    The dk_odds/mlb_odds table is fully rewritten each cycle (no history kept
    inside this DB — that's the offshore-style scraper convention; MLB.R is
    the consumer and snapshots into mlb.duckdb as needed).
    """
    db_path = Path(__file__).resolve().parent.parent / "dk_odds" / "dk.duckdb"
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
    p = argparse.ArgumentParser(description="DraftKings MLB singles scraper")
    p.add_argument("--verbose", action="store_true", help="Per-event row logging")
    args = p.parse_args()
    scrape_singles(verbose=args.verbose)


if __name__ == "__main__":
    main()
