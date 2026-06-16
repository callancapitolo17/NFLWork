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


def test_filter_a_extracts_signed_home_spreads_px():
    """Filter A regression: ``_extract_offered_lines_px`` must return the
    home-perspective signed spread set + flat total set.

    PX stores ``line`` on each selection from the *outcome's* own
    perspective. If home is favored at -1.5, the home selection carries
    ``line = -1.5`` and the away selection carries ``line = +1.5``. To
    get one home-perspective signed line per market, the helper must
    restrict to ``competitorId == home_id``. If we accidentally
    iterated away selections we'd get the mirrored sign, which would
    silently let through targets that don't match what PX actually
    offers from the home side.

    Totals are filtered to selection names starting with "over"/"under"
    so auxiliary buckets (if any) don't pollute the line set.
    """
    from mlb_sgp.prophetx import _extract_offered_lines_px

    # Build a minimal PX market list mirroring the shape that
    # ProphetXClient.fetch_event_markets emits (after `_iter_selections`
    # walks marketLines[].outcomes for the flat-list shape).
    markets = [
        {
            "name": "Run Line",
            "marketLines": [
                {
                    "outcomes": [
                        # Main -1.5 — home perspective.
                        {"competitorId": "home123", "line": -1.5,
                         "name": "Yankees -1.5"},
                        {"competitorId": "away456", "line": 1.5,
                         "name": "Red Sox +1.5"},
                        # Alt: home -2.5
                        {"competitorId": "home123", "line": -2.5,
                         "name": "Yankees -2.5"},
                        {"competitorId": "away456", "line": 2.5,
                         "name": "Red Sox +2.5"},
                    ],
                },
            ],
        },
        {
            "name": "Total Runs",
            "marketLines": [
                {
                    "outcomes": [
                        {"name": "Over 8.5", "line": 8.5},
                        {"name": "Under 8.5", "line": 8.5},
                        {"name": "Over 9.0", "line": 9.0},
                        {"name": "Under 9.0", "line": 9.0},
                        # Auxiliary garbage selection — must be filtered out.
                        {"name": "Push", "line": 9.0},
                    ],
                },
            ],
        },
    ]

    # Pass find_market + the FG market names directly so the helper
    # doesn't need to import scraper_prophetx_sgp (which would force
    # mlb_sgp/ onto sys.path and shadow the project-root db module).
    def find_market(mkts, name):
        for m in mkts:
            if m.get("name") == name:
                return m
        return None

    offered = _extract_offered_lines_px(
        markets, "fg", home_id="home123",
        find_market=find_market,
        spread_market_name="Run Line",
        total_market_name="Total Runs",
    )

    # Only home-side lines (negative signs preserved from outcome view).
    assert offered["spreads"] == {-1.5, -2.5}, (
        f"home-perspective signed spreads wrong: {offered['spreads']}"
    )
    # "Push" filtered out; over/under both collapse to {8.5, 9.0}.
    assert offered["totals"] == {8.5, 9.0}

    # No matching F5 markets in the fixture -> empty sets.
    empty = _extract_offered_lines_px(
        markets, "f5", home_id="home123",
        find_market=find_market,
        spread_market_name="1st-5th Inning Spread",
        total_market_name="1st-5th Inning Total Runs",
    )
    assert empty == {"spreads": set(), "totals": set()}


def test_px_parallelism_resolution(monkeypatch):
    """PX pool width: explicit arg > env > module default (4)."""
    from mlb_sgp import prophetx
    monkeypatch.setenv("MLB_SGP_PX_PARALLELISM", "5")
    assert prophetx._resolve_parallelism(None) == 5
    assert prophetx._resolve_parallelism(3) == 3
    monkeypatch.delenv("MLB_SGP_PX_PARALLELISM")
    assert prophetx._resolve_parallelism(None) == 4


def test_px_price_sgps_accepts_parallelism_kwarg():
    """Contract: the kwarg exists and the empty-input early-exit honors it."""
    from mlb_sgp import prophetx
    assert prophetx.price_sgps([], parallelism=4) == []
