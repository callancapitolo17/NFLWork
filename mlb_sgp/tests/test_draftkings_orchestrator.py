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


def test_filter_a_extracts_signed_home_spreads_dk():
    """Filter A regression: ``_extract_offered_lines_dk`` must return the
    home-perspective signed spread set + flat total set.

    DK keys spreads as (sign, abs_line, participant). Sign 'N' is the
    negative-line side; participant '1' is home, '3' is away. The helper
    must:
      - flip sign to negative when sign == 'N' on participant '1',
      - keep sign positive when sign == 'P' on participant '1',
      - ignore participant '3' (away mirrors home; double-counting).

    If this sign convention drifts, Filter A silently drops valid
    targets (the dropped target's spread won't appear in the offered
    set, so price_sgps skips it). This test pins the convention.
    """
    from mlb_sgp.draftkings import _extract_offered_lines_dk

    # Fixture mirrors what scraper_draftkings_sgp.fetch_selection_ids
    # returns post-parse: tuple keys, not JSON-encoded strings.
    sel_ids_all = {
        "fg": {
            "spreads": {
                # Home favored: home gets ('N', 1.5, '1') -> -1.5
                ("N", 1.5, "1"): "id_h_main",
                # Away mirror at +1.5 on participant '3' — ignored.
                ("P", 1.5, "3"): "id_a_main",
                # Alt home -2.5
                ("N", 2.5, "1"): "id_h_alt",
                ("P", 2.5, "3"): "id_a_alt",
                # Home underdog +1.5 (sign 'P' on home)
                ("P", 1.5, "1"): "id_h_dog",
                ("N", 1.5, "3"): "id_a_dog",
            },
            "totals": {
                ("O", 8.5): "id_o_main",
                ("U", 8.5): "id_u_main",
                ("O", 9.0): "id_o_int",
                ("U", 9.0): "id_u_int",
                ("O", 9.5): "id_o_alt",
                ("U", 9.5): "id_u_alt",
            },
        },
        "f5": {"spreads": {}, "totals": {}},
    }

    offered = _extract_offered_lines_dk(sel_ids_all, "fg")

    # Home-perspective signed lines: -1.5 (fav), -2.5 (fav alt), +1.5 (dog).
    assert offered["spreads"] == {-1.5, -2.5, 1.5}, (
        f"home-perspective signed spreads wrong: {offered['spreads']}"
    )
    # Totals flatten (over/under collapse).
    assert offered["totals"] == {8.5, 9.0, 9.5}

    # Missing period -> empty sets (defensive).
    empty = _extract_offered_lines_dk(sel_ids_all, "f5")
    assert empty == {"spreads": set(), "totals": set()}
