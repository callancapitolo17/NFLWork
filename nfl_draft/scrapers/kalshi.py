"""Kalshi adapter: calls legacy kalshi_draft/fetcher.py for raw market data,
returns List[OddsRow] for the new nfl_draft pipeline.

The legacy fetcher's side-effect writes (kalshi_odds, draft_series,
market_info, positions, resting_orders) happen as a byproduct of the calls
we make here - this keeps the Kalshi dashboard + MM bot fed from the same
code path that the new portal uses.
"""

import sys
from datetime import datetime
from pathlib import Path
from typing import List

# Ensure both the repo root AND kalshi_draft/ itself are importable - the
# legacy fetcher uses bare `from auth import ...` and `from db import ...`.
_REPO_ROOT = Path(__file__).resolve().parent.parent.parent
_KALSHI_DIR = _REPO_ROOT / "kalshi_draft"
for p in (_REPO_ROOT, _KALSHI_DIR):
    sp = str(p)
    if sp not in sys.path:
        sys.path.insert(0, sp)

from kalshi_draft import fetcher as legacy_fetcher  # noqa: E402
from kalshi_draft.auth import public_request  # noqa: E402

from nfl_draft.scrapers._base import OddsRow, TradeRow


def _extract_yes_bid_cents(market: dict):
    """Return Kalshi yes_bid as an integer number of cents, or None.

    Kalshi's v2 response includes both legacy (int cents) and new
    (`_dollars` string) keys depending on the account's `response_price_units`.
    We accept either and normalize to cents.
    """
    raw_cents = market.get("yes_bid")
    if raw_cents not in (None, 0):
        try:
            return int(raw_cents)
        except (TypeError, ValueError):
            pass
    raw_dollars = market.get("yes_bid_dollars")
    if raw_dollars is not None:
        try:
            return int(round(float(raw_dollars) * 100))
        except (TypeError, ValueError):
            return None
    return raw_cents  # may be 0 or None


def parse_markets_response(raw_response: dict, series_ticker: str) -> List[OddsRow]:
    """Convert a raw Kalshi /markets response into a list of OddsRow.

    One row per market with a usable yes_bid. Skips markets with zero/empty
    bids (no live market price) and markets missing a candidate name.
    """
    rows: List[OddsRow] = []
    if not isinstance(raw_response, dict):
        return rows
    now = datetime.now()

    for market in raw_response.get("markets", []) or []:
        yes_bid_cents = _extract_yes_bid_cents(market)
        if yes_bid_cents is None or yes_bid_cents <= 0:
            continue
        p = yes_bid_cents / 100.0
        if p <= 0 or p >= 1:
            continue
        # Convert implied probability to American odds
        if p > 0.5:
            american = int(round(-100 * p / (1 - p)))
        else:
            american = int(round(100 * (1 - p) / p))

        candidate = (
            market.get("yes_sub_title")
            or market.get("subtitle")
            or (market.get("custom_strike") or {}).get("Person")
            or (market.get("custom_strike") or {}).get("Team")
            or ""
        )
        candidate = (candidate or "").strip()
        if not candidate:
            continue

        rows.append(OddsRow(
            book="kalshi",
            book_label=series_ticker,
            book_subject=candidate,
            american_odds=american,
            fetched_at=now,
        ))
    return rows


def fetch_draft_odds() -> List[OddsRow]:
    """Discover all NFL Draft series, fetch their open markets, return rows.

    Side effect: calling the legacy fetcher writes kalshi_odds + draft_series
    + market_info into nfl_draft.duckdb (the legacy fetcher has been
    repointed to the new DB in Task 3-9).
    """
    rows: List[OddsRow] = []
    series_list = legacy_fetcher.discover_draft_series()
    for series in series_list:
        ticker = series["series_ticker"]
        # Legacy side-effect write (fills kalshi_odds via the run pipeline).
        # discover + fetch_markets_for_series together are what the legacy
        # dashboard depends on - we mirror that call pattern so nothing
        # downstream regresses.
        legacy_fetcher.fetch_markets_for_series(ticker)
        # Re-parse once more to produce normalized OddsRow output for the
        # new portal. Re-fetching is cheap (and now throttled via auth.py).
        raw = public_request(
            f"/markets?series_ticker={ticker}&status=open&limit=100"
        )
        if raw:
            rows.extend(parse_markets_response(raw, ticker))
    return rows


def fetch_trades() -> List[TradeRow]:
    """Poll Kalshi /markets/trades for each NFL draft series since last cursor.

    Stub - implemented in Task 21 (trade-tape poller).
    """
    raise NotImplementedError("Implemented in Task 21 (trade-tape poller)")
