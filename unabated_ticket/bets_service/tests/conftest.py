"""Shared fixture loaders. The Kalshi fixture is the one the node tests use
(unabated_ticket/tests/fixtures/bets/kalshi_fixture.json): 19 real anonymised
fills, 8 positions, 12 market and 12 event payloads keyed by ticker / event_ticker."""
import json
from pathlib import Path

import pytest

FIXTURE_PATH = Path(__file__).parents[2] / "tests" / "fixtures" / "bets" / "kalshi_fixture.json"
FETCHED_AT = "2026-09-11T21:00:00Z"
# The Novig fixture (#116, REST rebuild 2026-09-22): live Portfolio cards from
# GET /nbx/v1/portfolio/{active,settled}, trimmed (see its "provenance").
NOVIG_FIXTURE_PATH = Path(__file__).parents[2] / "tests" / "fixtures" / "bets" / "novig_portfolio.json"
NOVIG_FETCHED_AT = "2026-09-22T21:00:00Z"


@pytest.fixture
def kalshi_fixture() -> dict:
    return json.loads(FIXTURE_PATH.read_text())


def fill(ticker: str, side: str, count_fp: float, yes_price: float, created_time: str,
         action: str = "buy", trade_id: str | None = None) -> dict:
    """A fill on the wire grammar (bets.test.js `fill`), with an optional trade_id."""
    record = {
        "ticker": ticker, "market_ticker": ticker, "side": side, "action": action,
        "count_fp": f"{count_fp}", "yes_price_dollars": f"{yes_price:.4f}",
        "no_price_dollars": f"{1 - yes_price:.4f}", "created_time": created_time,
        "is_taker": True, "fee_cost": "0.000000", "exchange_index": 0, "book_side": "bid",
        "outcome_side": side,
    }
    if trade_id:
        record["trade_id"] = trade_id
    return record


def position(ticker: str, position_fp: float, total_traded: float,
             last_updated_ts: str = "2026-09-11T15:05:04.400227Z") -> dict:
    return {
        "ticker": ticker, "position_fp": f"{position_fp:.2f}", "total_traded_dollars": f"{total_traded:.6f}",
        "market_exposure_dollars": "0.000000", "realized_pnl_dollars": "0.000000",
        "fees_paid_dollars": "0.000000", "exchange_index": 0, "last_updated_ts": last_updated_ts,
    }


def market(ticker: str, event_ticker: str, strike_type: str, floor_strike: float | None,
           title: str, result: str = "") -> dict:
    return {
        "ticker": ticker, "event_ticker": event_ticker, "title": title, "yes_sub_title": title,
        "strike_type": strike_type, "floor_strike": floor_strike, "custom_strike": None,
        "status": "finalized" if result else "active", "result": result,
        "expected_expiration_time": "2026-09-13T21:00:00Z", "close_time": "2026-09-15T17:00:00Z",
        "rules_primary": "",
    }


def event(event_ticker: str, series_ticker: str, title: str, sub_title: str) -> dict:
    return {"event_ticker": event_ticker, "series_ticker": series_ticker, "title": title,
            "sub_title": sub_title, "mutually_exclusive": series_ticker.endswith("GAME")}


@pytest.fixture
def novig_fixture() -> dict:
    """{provenance, active: [card], settled: [card]} — the two Portfolio lists."""
    return json.loads(NOVIG_FIXTURE_PATH.read_text())
