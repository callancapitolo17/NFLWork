"""Tests for mlb_sgp/fanduel.py::price_sgps() orchestrator.

The orchestrator is a thin composition layer over the well-tested
helper functions in scraper_fanduel_sgp.py. These tests focus on the
orchestration contract (input/output shape, period filtering,
source-label tagging) — the underlying helper logic is covered by
test_sgp_regression.py and test_integration_integer_line.py.
"""
from datetime import datetime, timezone
from unittest.mock import MagicMock

from mlb_sgp._shared import TargetLine
from mlb_sgp.fanduel import (
    BOOK_NAME,
    SOURCE_LABEL,
    SOURCE_LABEL_FALLBACK,
    price_sgps,
)


def test_price_sgps_empty_returns_empty():
    """Empty input -> empty output, no API calls."""
    assert price_sgps([], periods=("FG",)) == []


def test_price_sgps_period_filter_drops_excluded():
    """If periods=('FG',), F5 TargetLines are filtered out.

    Even with no real client, the function must return [] without
    crashing because all targets get dropped *before* the client is
    touched. The MagicMock client proves the function never bothers to
    issue any HTTP calls once the period filter empties targets.
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
    # Stronger assertion than the DK test: neither GET nor POST should
    # have been issued. FD's scraper uses both verbs (POST for the scan
    # endpoint, GET for event-page), so we check each independently.
    mock_client.session.get.assert_not_called()
    mock_client.session.post.assert_not_called()


def test_source_labels_exposed():
    """Module exports the two source labels for downstream consumers.

    Task 14 (scraper_fanduel_sgp shim) and the dashboard reader both
    filter rows by ``source IN (..._direct)`` — keeping these constants
    in lockstep with the legacy scraper output protects that contract.
    """
    assert SOURCE_LABEL == "fanduel_direct"
    assert SOURCE_LABEL_FALLBACK == "fanduel_interpolated"
    assert BOOK_NAME == "fanduel"
