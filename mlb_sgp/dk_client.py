"""DraftKings API client — extracted from scraper_draftkings_sgp.py.

Owns the curl_cffi Chrome-TLS session and exposes the three operations
both the SGP scraper and the singles scraper need:
  - list_events()           — all MLB events today
  - fetch_event_markets()   — market metadata (FG vs F5 vs alt)
  - fetch_event_selections()— all selections with prices, in one call
"""
from __future__ import annotations
from dataclasses import dataclass


@dataclass
class Event:
    event_id: str
    home_team: str
    away_team: str
    start_time: str  # ISO UTC string


@dataclass
class Market:
    market_id: str
    name: str
    subcategory: str  # "game_lines" | "alt_lines" | "innings" etc.


@dataclass
class Selection:
    selection_id: str
    market_id: str
    name: str            # e.g. "Yankees -1.5"
    line: float | None   # spread/total value; None for moneylines
    american_odds: int


class DraftKingsClient:
    def __init__(self, verbose: bool = False) -> None:
        from scraper_draftkings_sgp import init_session
        self.session = init_session()
        self.verbose = verbose

    def list_events(self) -> list[Event]:
        raise NotImplementedError("Task 3")

    def fetch_event_markets(self, event_id: str) -> list[Market]:
        raise NotImplementedError("Task 4")

    def fetch_event_selections(self, event_id: str) -> list[Selection]:
        raise NotImplementedError("Task 4")
