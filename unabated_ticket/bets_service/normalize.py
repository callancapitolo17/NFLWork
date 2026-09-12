"""Shared record helpers for every bet source: price, money and Eastern-time
conventions of the normalised bet record (contract in
docs/2026-09-11-issue-114-bet-history-plan.md).

Pure functions, no I/O. They mirror the helpers in extension/bets.js exactly,
including JavaScript's Math.round (half toward +infinity, not Python's
banker's rounding), so a Python-normalised record is byte-equivalent to the
one the node tests pin.
"""
import math
from datetime import datetime, timezone
from zoneinfo import ZoneInfo

EASTERN = ZoneInfo("America/New_York")



def js_round(value: float) -> int:
    """JavaScript Math.round: halves round toward +infinity."""
    return int(math.floor(value + 0.5))


def to_number(value: object) -> float:
    """Kalshi sends numbers as strings ("250.00"); unparseable -> 0 (bets.js toNumber)."""
    if isinstance(value, bool) or value is None:
        return 0.0
    try:
        number = float(value)
    except (TypeError, ValueError):
        return 0.0
    return number if math.isfinite(number) else 0.0


def cents_to_american(cents: float | None) -> int | None:
    """Contract price in cents -> American odds; 50c is -100, outside (0, 100) is None."""
    if cents is None or not (0 < cents < 100):
        return None
    if cents >= 50:
        return -js_round(cents / (100 - cents) * 100)
    return js_round((100 - cents) / cents * 100)


def round_cents(dollars: float) -> float:
    return js_round(dollars * 100) / 100


def parse_iso_ms(iso: str) -> int | None:
    """ISO timestamp -> epoch milliseconds, truncated like JavaScript Date.parse."""
    try:
        parsed = datetime.fromisoformat(iso.replace("Z", "+00:00"))
    except (TypeError, ValueError):
        return None
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return math.floor(parsed.timestamp() * 1000)


def eastern_wall_clock_to_iso(year: int, month: int, day: int, hour: int, minute: int) -> str:
    """Eastern wall-clock time -> ISO UTC with millisecond zeros ("...:00.000Z"),
    the string JavaScript's Date.toISOString produces."""
    local = datetime(year, month, day, hour, minute, tzinfo=EASTERN)
    return local.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.000Z")


def json_clean(value: object) -> object:
    """Integral floats -> ints (400.0 -> 400) so the JSON matches what JavaScript
    prints; applied recursively to a record before it is stored or served."""
    if isinstance(value, float) and value.is_integer():
        return int(value)
    if isinstance(value, dict):
        return {key: json_clean(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_clean(item) for item in value]
    return value


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
