"""Shared scraper data structures."""

from dataclasses import dataclass
from datetime import datetime
from typing import Optional


@dataclass
class OddsRow:
    """A single row of raw scraper output (one binary outcome from one venue).

    The `implied_prob` and `devig_prob` fields are optional overrides. When
    set, `quarantine.write_or_quarantine` uses them verbatim for the
    `draft_odds.implied_prob` / `devig_prob` columns instead of deriving both
    from `american_to_implied(american_odds)`. Used by Kalshi so that the
    take-price (buy / yes_ask) and fair-value (mid) can be persisted without
    rounding loss through the American-odds round-trip.
    """
    book: str
    book_label: str
    book_subject: str
    american_odds: int
    fetched_at: datetime
    market_group: str = ""
    implied_prob: Optional[float] = None
    devig_prob: Optional[float] = None


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
