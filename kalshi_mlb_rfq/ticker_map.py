"""Map internal game/leg descriptors to Kalshi ticker strings.

Per recon (2026-04-27):
  KXMLBGAME   single ticker per (event, team) — winner
  KXMLBSPREAD multiple tickers per event: -{TEAM}{N}, "team wins by over (N-1).5"
  KXMLBTOTAL  multiple tickers per event: -{N}, "total over (N-1).5"
  KXMLBRFI    single ticker per event, no suffix — run in 1st inning

Event suffix format: YYMMMDDHHMM{AwayCode}{HomeCode}, all in ET.
"""

from datetime import datetime
from zoneinfo import ZoneInfo

ET = ZoneInfo("America/New_York")
_MONTHS = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN",
           "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]


def format_event_suffix(commence_time: datetime, away_code: str, home_code: str) -> str:
    """commence_time may be UTC or any tz-aware datetime; converted to ET for the suffix."""
    if commence_time.tzinfo is None:
        raise ValueError("commence_time must be timezone-aware")
    et = commence_time.astimezone(ET)
    return (
        f"{et.year % 100:02d}"
        f"{_MONTHS[et.month - 1]}"
        f"{et.day:02d}"
        f"{et.hour:02d}"
        f"{et.minute:02d}"
        f"{away_code}{home_code}"
    )


def spread_ticker(event_suffix: str, team_code: str, line: float) -> str:
    """Spread ticker. line is the team's spread (negative = favorite).

    -1.5 → suffix N=2 (team wins by over 1.5)
    -2.5 → suffix N=3
    -3.5 → suffix N=4
    """
    if line >= 0:
        raise ValueError(f"spread_ticker expects negative line for favorite; got {line}")
    n = int(round(-line + 0.5))
    return f"KXMLBSPREAD-{event_suffix}-{team_code}{n}"


def total_ticker(event_suffix: str, line: float) -> str:
    """Total ticker. line is the .5-suffixed total (e.g. 7.5).

    7.5 → -8
    8.5 → -9
    10.5 → -11
    """
    n = int(round(line + 0.5))
    return f"KXMLBTOTAL-{event_suffix}-{n}"


def game_ticker(event_suffix: str, team_code: str) -> str:
    return f"KXMLBGAME-{event_suffix}-{team_code}"


def rfi_ticker(event_suffix: str) -> str:
    return f"KXMLBRFI-{event_suffix}"
