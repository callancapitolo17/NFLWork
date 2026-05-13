"""Tests for mlb_sgp/prophetx.py::price_sgps() orchestrator.

The orchestrator is a thin composition layer over the well-tested
helpers in scraper_prophetx_sgp.py and the ProphetXClient. These tests
focus on the orchestration contract (input/output shape, period
filtering, source-label tagging, sanity-filter behavior) — the
underlying helper logic is covered by test_prophetx_client.py and
test_sgp_regression.py.
"""
from datetime import datetime, timezone
from unittest.mock import MagicMock

from mlb_sgp._shared import TargetLine
from mlb_sgp.prophetx import (
    BOOK_NAME,
    SANITY_MULT_RATIO,
    SOURCE_LABEL,
    SOURCE_LABEL_FALLBACK,
    _passes_sanity_mult_ratio,
    price_sgps,
)


def test_price_sgps_empty_returns_empty():
    """Empty input -> empty output, no API calls."""
    assert price_sgps([], periods=("FG",)) == []


def test_price_sgps_period_filter_drops_excluded():
    """If periods=('FG',), F5 TargetLines are filtered out before any HTTP.

    The MagicMock client proves the function never reaches the client's
    list_events()/fetch_event_markets()/submit_parlay_rfq() methods once
    the period filter empties the targets list.
    """
    targets = [
        TargetLine(
            game_id="g1", home_team="NYY", away_team="BOS",
            commence_time=datetime(2026, 5, 13, 23, 0, tzinfo=timezone.utc),
            period="F5", spread=-0.5, total=4.5,
        ),
    ]
    mock_client = MagicMock()
    out = price_sgps(targets, periods=("FG",), client=mock_client)
    assert out == []
    mock_client.list_events.assert_not_called()
    mock_client.fetch_event_markets.assert_not_called()
    mock_client.submit_parlay_rfq.assert_not_called()


def test_sanity_mult_ratio_passes_normal():
    """Normal correlated parlay (parlay/naive ~= 1.0x) passes the filter."""
    # Two -110 legs (decimal 1.91). Naive = 1.91 * 1.91 ~= 3.65. A
    # normal correlated parlay of decimal 2.50 falls well below the
    # 1.5x naive ceiling, so it should pass.
    assert _passes_sanity_mult_ratio(
        parlay_decimal=2.50, leg1_dec=1.91, leg2_dec=1.91
    ) is True


def test_sanity_mult_ratio_blocks_f5_over_anomaly():
    """F5-Over bug: parlay decimal is 5-7x naive. Block it."""
    # Naive independent: 1.5 * 2.0 = 3.0. Bug-priced parlay at 25.0 is
    # ~8.3x naive, far above the 1.5x ceiling -> drop.
    assert _passes_sanity_mult_ratio(
        parlay_decimal=25.0, leg1_dec=1.5, leg2_dec=2.0
    ) is False


def test_sanity_mult_ratio_handles_zero():
    """Defensive: zero or negative naive product -> drop the combo.

    Returning False here means a missing single-leg odds (-> 0.0 from
    _safe_leg_decimal) causes the combo to be silently dropped, which
    is the safe behavior: we can't sanity-check what we can't measure.
    """
    assert _passes_sanity_mult_ratio(
        parlay_decimal=2.0, leg1_dec=0, leg2_dec=2.0
    ) is False


def test_source_labels_exposed():
    """Module exports the source-label / book-name constants used by the
    Task 15 shim and the dashboard reader. Lockstep with the legacy
    scraper output keeps the dashboard contract byte-identical.
    """
    assert SOURCE_LABEL == "prophetx_direct"
    assert SOURCE_LABEL_FALLBACK == "prophetx_interpolated"
    assert BOOK_NAME == "prophetx"
    assert SANITY_MULT_RATIO == 1.5
