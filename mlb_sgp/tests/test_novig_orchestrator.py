"""Tests for mlb_sgp/novig.py::price_sgps() orchestrator.

The orchestrator is a thin composition layer over the well-tested
helpers in scraper_novig_sgp.py and the NovigClient. These tests
focus on the orchestration contract (input/output shape, period
filtering, source-label tagging, sanity-filter behavior) — the
underlying helper logic is covered by test_novig_client.py and
test_sgp_regression.py.
"""
from datetime import datetime, timezone
from unittest.mock import MagicMock

from mlb_sgp._shared import TargetLine
from mlb_sgp.novig import (
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
    list_events()/fetch_event_legs()/submit_parlay() methods once the
    period filter empties the targets list.
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
    mock_client.fetch_event_legs.assert_not_called()
    mock_client.submit_parlay.assert_not_called()


def test_sanity_mult_ratio_passes_normal():
    """Normal correlated parlay (parlay/naive ~= 1.0x) passes the filter.

    Both legs implied prob 0.524 (decimal ~1.91). Naive product decimal
    = 1 / (0.524 * 0.524) = ~3.64. A normal correlated parlay of
    decimal 2.50 falls well below the 1.5x naive ceiling -> pass.
    """
    assert _passes_sanity_mult_ratio(
        parlay_decimal=2.50, sp_available=0.524, to_available=0.524
    ) is True


def test_sanity_mult_ratio_blocks_anomaly():
    """Anomalous mispricing: parlay decimal multiples larger than naive.

    Both legs implied prob 0.5 (decimal 2.0). Naive product decimal
    = 4.0. Bug-priced parlay at 25.0 is ~6.25x naive, far above the
    1.5x ceiling -> drop.
    """
    assert _passes_sanity_mult_ratio(
        parlay_decimal=25.0, sp_available=0.5, to_available=0.5
    ) is False


def test_sanity_mult_ratio_missing_available_keeps_combo():
    """Defensive: if either leg's implied prob is missing, keep the combo.

    Matches the legacy scraper's behavior — the sanity check only runs
    when both leg implied probs are available. Without them we can't
    compute the naive product, so the safe move is to keep the combo
    (the alternative would be silently dropping every priced row, which
    would mask data issues rather than surface them).
    """
    assert _passes_sanity_mult_ratio(
        parlay_decimal=2.0, sp_available=None, to_available=0.5
    ) is True
    assert _passes_sanity_mult_ratio(
        parlay_decimal=2.0, sp_available=0.5, to_available=None
    ) is True


def test_source_labels_exposed():
    """Module exports the source-label / book-name constants used by the
    Task 16 shim and the dashboard reader. Lockstep with the legacy
    scraper output keeps the dashboard contract byte-identical.
    """
    assert SOURCE_LABEL == "novig_direct"
    assert SOURCE_LABEL_FALLBACK == "novig_interpolated"
    assert BOOK_NAME == "novig"
    assert SANITY_MULT_RATIO == 1.5
