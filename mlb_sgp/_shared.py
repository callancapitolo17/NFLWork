"""Shared types and utilities for the MLB SGP library.

Module-private (`_` prefix) but stable API consumed by per-book modules
(draftkings.py, fanduel.py, prophetx.py, novig.py) and by callers
(dashboard scraper shims, kalshi_mlb_rfq/sgp_runner.py).
"""
from __future__ import annotations
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Literal

Period = Literal["FG", "F5"]


@dataclass(frozen=True)
class TargetLine:
    """One (game, period, spread, total) tuple the bot/dashboard wants priced."""
    game_id: str
    home_team: str
    away_team: str
    commence_time: datetime
    period: Period
    spread: float   # signed, home-perspective (negative = home favored)
    total: float


@dataclass(frozen=True)
class PricedRow:
    """One book's SGP price for a (game, period, spread, total, combo) tuple."""
    game_id: str
    combo: str
    period: Period
    spread_line: float
    total_line: float
    bookmaker: str
    source: str
    sgp_decimal: float
    sgp_american: int
    fetch_time: datetime


def decimal_to_american(dec: float) -> int:
    """Convert decimal odds to American format. Favorites are negative."""
    if dec >= 2.0:
        return int(round((dec - 1.0) * 100))
    return int(round(-100.0 / (dec - 1.0)))


def american_to_decimal(am: int) -> float:
    """Inverse of decimal_to_american."""
    if am > 0:
        return 1.0 + am / 100.0
    return 1.0 + 100.0 / abs(am)


def _utc_bucket(ts: datetime | str) -> str:
    """Extract a UTC "YYYY-MM-DDTHH" bucket string from a timestamp.

    Used as a match key when correlating events across data sources at
    date+hour granularity (avoids spurious matches across doubleheaders).
    Accepts both datetime objects and ISO-8601 strings.

    Empty/None input returns "" to match existing scraper behavior
    (callers pass `lines.get("commence_time", "")` and rely on "" to
    trigger their team-only fallback matcher).
    """
    if not ts:
        return ""
    if isinstance(ts, str):
        # Normalize Z suffix to +00:00 for fromisoformat
        normalized = ts.replace("Z", "+00:00") if ts.endswith("Z") else ts
        ts = datetime.fromisoformat(normalized)
    if ts.tzinfo is None:
        ts = ts.replace(tzinfo=timezone.utc)
    else:
        ts = ts.astimezone(timezone.utc)
    return ts.strftime("%Y-%m-%dT%H")
