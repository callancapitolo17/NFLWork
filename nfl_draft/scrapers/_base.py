"""Shared scraper data structures."""

from dataclasses import dataclass
from datetime import datetime


@dataclass
class OddsRow:
    """A single row of raw scraper output (one binary outcome from one venue)."""
    book: str
    book_label: str
    book_subject: str
    american_odds: int
    fetched_at: datetime
    market_group: str = ""


@dataclass
class TradeRow:
    """A single Kalshi trade event (for the trade tape)."""
    trade_id: str
    ticker: str
    side: str
    price_cents: int
    count: int
    traded_at: datetime
    fetched_at: datetime
