"""Regression tests locking 2-leg equivalence between router.combo_fair and
the original _book_fairs + statistics.median path, plus basic cross-game
and Phase-1 scope checks.
"""
import statistics

import pandas as pd
import pytest

from kalshi_common import legset
from kalshi_common.leg_types import combo_descriptor, SPREAD_TOTAL_FAMILY
from kalshi_mlb_mm import router

# ---------------------------------------------------------------------------
# Shared fixture helpers
# ---------------------------------------------------------------------------

GAME_ID = "gTEXLAA-integration"  # DB-side game id (passed via resolve_game)
EVT = "KXMLBGAME-25JUN271905TEXLAA"  # Kalshi event ticker
EVT2 = "KXMLBGAME-25JUN271905NYYBOS"  # Second game for cross-game tests

# Valid 2-leg spread×total legs that parse through legset.parse_legs:
#   LAA=home (3-letter), TEX=away; LAA2 → home spread -1.5; 9 → total 8.5
LEGS_SPREAD_TOTAL = [
    {"market_ticker": "KXMLBSPREAD-25JUN271905TEXLAA-LAA2",
     "event_ticker": EVT, "side": "yes"},  # Home Spread -1.5 YES
    {"market_ticker": "KXMLBTOTAL-25JUN271905TEXLAA-9",
     "event_ticker": EVT, "side": "yes"},  # Over 8.5 YES
]

# Balanced 4-cell grid: each cell gets the same decimal so vig is near-zero.
ST_CELLS = {
    "Home Spread + Over": 4.2,
    "Home Spread + Under": 4.2,
    "Away Spread + Over": 3.8,
    "Away Spread + Under": 3.8,
}


def _grid_rows(book, game_id, spread_line, total_line, cells):
    return [
        {"game_id": game_id, "combo": c, "period": "FG", "bookmaker": book,
         "sgp_decimal": d, "fetch_time": None,
         "spread_line": spread_line, "total_line": total_line}
        for c, d in cells.items()
    ]


# ---------------------------------------------------------------------------
# Core regression: router.combo_fair == median(_book_fairs) for 2-leg grid
# ---------------------------------------------------------------------------

def test_router_combo_fair_equals_book_fairs_median(monkeypatch):
    """For a 2-leg spread×total combo on one game, router.combo_fair must
    equal statistics.median(main._book_fairs(game_id, desc).values()) to
    within 1e-9.  This locks the two pricing paths to identical arithmetic.
    """
    from kalshi_mlb_mm import main
    import kalshi_mlb_mm.config as cfg

    # Two books with a balanced grid (no outlier — both within BAND).
    df = pd.DataFrame(
        _grid_rows("dk", GAME_ID, -1.5, 8.5, ST_CELLS)
        + _grid_rows("fd", GAME_ID, -1.5, 8.5, ST_CELLS)
    )

    # Inject into the module global that _book_fairs reads.
    monkeypatch.setattr(main, "_SGP_ODDS", df)

    legs = LEGS_SPREAD_TOTAL
    desc = combo_descriptor(legs)
    assert desc is not None, "Fixture legs must produce a valid ComboDescriptor"

    # Old path: _book_fairs → consensus filter → median
    book_fairs = main._book_fairs(GAME_ID, desc)
    assert book_fairs, "Fixture must yield at least 1 agreeing book"
    old_fair = statistics.median(book_fairs.values())

    # New path: router.combo_fair with identity resolve_game
    new_fair = router.combo_fair(legs, df, lambda gl: GAME_ID,
                                 cfg.MIN_AGREEING_BOOKS, cfg.BOOK_CONSENSUS_BAND)
    assert new_fair is not None, "router.combo_fair must price this combo"

    assert abs(new_fair - old_fair) < 1e-9, (
        f"Pricing paths diverge: old={old_fair:.12f} new={new_fair:.12f} "
        f"delta={abs(new_fair - old_fair):.2e}")


def test_router_combo_fair_equals_book_fairs_with_outlier(monkeypatch):
    """Outlier rejection: one rogue book should be filtered by both paths
    identically, yielding equal medians to 1e-9.
    """
    from kalshi_mlb_mm import main
    import kalshi_mlb_mm.config as cfg

    # Three books: dk + fd agree near 0.25; px is a clear outlier (0.60).
    OUTLIER_CELLS_PX = {k: 1 / 0.40 for k in ST_CELLS}  # ~0.40 each → devigged ~0.40+
    df = pd.DataFrame(
        _grid_rows("dk", GAME_ID, -1.5, 8.5, ST_CELLS)
        + _grid_rows("fd", GAME_ID, -1.5, 8.5, ST_CELLS)
        + _grid_rows("px", GAME_ID, -1.5, 8.5, OUTLIER_CELLS_PX)
    )
    monkeypatch.setattr(main, "_SGP_ODDS", df)

    legs = LEGS_SPREAD_TOTAL
    desc = combo_descriptor(legs)

    book_fairs = main._book_fairs(GAME_ID, desc)
    old_fair = statistics.median(book_fairs.values())

    new_fair = router.combo_fair(legs, df, lambda gl: GAME_ID,
                                 cfg.MIN_AGREEING_BOOKS, cfg.BOOK_CONSENSUS_BAND)
    assert new_fair is not None
    assert abs(new_fair - old_fair) < 1e-9, (
        f"Outlier-filtered paths diverge: old={old_fair:.12f} new={new_fair:.12f}")


# ---------------------------------------------------------------------------
# Cross-game 2-game combo prices through _priceable_in_phase1 + router
# ---------------------------------------------------------------------------

def test_cross_game_combo_prices_via_priceable_in_phase1():
    """A 2-game combo (one ML leg per game) should be priceable in Phase 1
    (each game is a single-leg marginal) and router.combo_fair multiplies
    the per-game consensus fairs.
    """
    # ML family rows (spread_line=None) for two different games.
    ML_CELLS = {
        "Home ML + Over": 3.6, "Home ML + Under": 4.0,
        "Away ML + Over": 4.4, "Away ML + Under": 4.0,
    }

    GAME1 = "g-game1"
    GAME2 = "g-game2"

    df = pd.DataFrame(
        _grid_rows("dk", GAME1, None, 8.5, ML_CELLS)
        + _grid_rows("fd", GAME1, None, 8.5, ML_CELLS)
        + _grid_rows("dk", GAME2, None, 8.5, ML_CELLS)
        + _grid_rows("fd", GAME2, None, 8.5, ML_CELLS)
    )

    # Cross-game legs: one home ML from each game.
    legs_dicts = [
        {"market_ticker": f"KXMLBGAME-25JUN271905TEXLAA-LAA",
         "event_ticker": EVT, "side": "yes"},   # Home ML game1
        {"market_ticker": f"KXMLBGAME-25JUN271905NYYBOS-BOS",
         "event_ticker": EVT2, "side": "yes"},  # Home ML game2
    ]

    # resolve_game maps event_ticker (canon.game_id) → DB game_id
    resolve = lambda gl: GAME1 if gl[0].game_id == EVT else GAME2

    # Phase-1 scope check: each game has 1 leg → "single"
    from kalshi_mlb_mm.main import _priceable_in_phase1
    canon = legset.parse_legs(legs_dicts)
    assert canon is not None, "Legs must parse"
    assert _priceable_in_phase1(canon), "Cross-game single legs must be priceable in Phase 1"

    out = router.combo_fair(legs_dicts, df, resolve, min_books=2, band=0.02)
    assert out is not None, "router.combo_fair must price 2-game single-leg combo"

    # Verify it's the product of the two per-game fairs.
    parts = legset.partition_by_game(canon)
    game_legs_per_game = list(parts.values())
    f1 = router.subcombo_fair(GAME1, game_legs_per_game[0], df, 2, 0.02)
    f2 = router.subcombo_fair(GAME2, game_legs_per_game[1], df, 2, 0.02)
    assert f1 is not None and f2 is not None
    assert abs(out - f1 * f2) < 1e-9, (
        f"Cross-game product mismatch: combo={out:.10f} f1*f2={f1*f2:.10f}")


# ---------------------------------------------------------------------------
# _priceable_in_phase1 scope gate
# ---------------------------------------------------------------------------

def test_priceable_in_phase1_accepts_grid_combos():
    from kalshi_mlb_mm.main import _priceable_in_phase1
    # Single-game spread×total
    canon = legset.parse_legs(LEGS_SPREAD_TOTAL)
    assert _priceable_in_phase1(canon)


def test_priceable_in_phase1_rejects_on_demand():
    """3+ legs same game → on_demand → not priceable in Phase 1."""
    from kalshi_mlb_mm.main import _priceable_in_phase1
    # 3-leg combo: spread + total + ml = on_demand
    three_legs = [
        {"market_ticker": "KXMLBSPREAD-25JUN271905TEXLAA-LAA2",
         "event_ticker": EVT, "side": "yes"},
        {"market_ticker": "KXMLBTOTAL-25JUN271905TEXLAA-9",
         "event_ticker": EVT, "side": "yes"},
        {"market_ticker": "KXMLBGAME-25JUN271905TEXLAA-LAA",
         "event_ticker": EVT, "side": "yes"},
    ]
    canon = legset.parse_legs(three_legs)
    assert canon is not None
    assert not _priceable_in_phase1(canon)


def test_priceable_in_phase1_rejects_untypeable():
    """Untypeable leg → parse_legs returns None → _priceable_in_phase1 not called."""
    from kalshi_mlb_mm.main import _priceable_in_phase1
    bad = [{"market_ticker": "KXMLBPLAYER-foo", "event_ticker": EVT, "side": "yes"}]
    canon = legset.parse_legs(bad)
    assert canon is None   # parse_legs guards before we can call _priceable_in_phase1
