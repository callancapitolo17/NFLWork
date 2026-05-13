"""Tests for mlb_sgp/draftkings.py::price_sgps() orchestrator.

The orchestrator is a thin composition layer over the well-tested
helper functions in scraper_draftkings_sgp.py. These tests focus on
the orchestration contract (input/output shape, period filtering,
source-label tagging) — the underlying helper logic is covered by
test_sgp_regression.py and test_integration_integer_line.py.
"""
from datetime import datetime, timezone
from unittest.mock import MagicMock

from mlb_sgp._shared import TargetLine
from mlb_sgp.draftkings import (
    BOOK_NAME,
    SOURCE_LABEL,
    SOURCE_LABEL_FALLBACK,
    price_sgps,
)


def test_price_sgps_empty_targets_returns_empty():
    """Empty input -> empty output, no API calls."""
    out = price_sgps([], periods=("FG",))
    assert out == []


def test_price_sgps_period_filter_drops_excluded():
    """If periods=('FG',), F5 TargetLines are filtered out.

    Even with no real client, the function must return [] without
    crashing because all targets get dropped *before* the client is
    touched. We pass a MagicMock client to prove the function never
    bothers to look at it once the period filter has emptied targets.
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
    # The client's session attribute must never have been read — proves
    # the early-exit short-circuits before any HTTP setup.
    mock_client.session.assert_not_called()


def test_source_labels_exposed():
    """Module exports the two source labels for downstream consumers.

    Task 13 (scraper_draftkings_sgp shim) and the dashboard reader both
    filter rows by ``source IN (..._direct)`` — keeping these constants
    in lockstep with the legacy scraper output protects that contract.
    """
    assert SOURCE_LABEL == "draftkings_direct"
    assert SOURCE_LABEL_FALLBACK == "draftkings_interpolated"
    assert BOOK_NAME == "draftkings"
