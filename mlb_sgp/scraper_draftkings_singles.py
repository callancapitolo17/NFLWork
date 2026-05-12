"""DraftKings single-leg odds scraper.

Walks every MLB event today, fetches all selections + market metadata via
DraftKingsClient, transforms to the offshore mlb_odds schema (wagerzon-style),
and writes to dk_odds/dk.duckdb.

Task 8 implements the parser. Task 9 wires the scrape loop, classify_market
function, DK team-name canonicalization, and DuckDB write.
"""
from __future__ import annotations
from datetime import datetime
from typing import Any

from dk_client import Event, Selection


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
