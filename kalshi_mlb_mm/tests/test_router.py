import pandas as pd
import pytest
from kalshi_common import legset
from kalshi_mlb_mm import router

EVT = "KXMLBGAME-25JUN271905NYYBOS"

def _grid_rows(book, game_id, spread_line, total_line, family_decimals):
    # family_decimals: dict of combo-name -> decimal odds (4 cells)
    return [
        {"game_id": game_id, "combo": c, "period": "FG", "bookmaker": book,
         "sgp_decimal": d, "fetch_time": None,
         "spread_line": spread_line, "total_line": total_line}
        for c, d in family_decimals.items()
    ]

# A balanced 4-cell grid (~4.0 decimal each => ~0.25 raw, sums ~1.0 -> tiny vig)
ST_CELLS = {"Home Spread + Over": 4.2, "Home Spread + Under": 4.2,
            "Away Spread + Over": 3.8, "Away Spread + Under": 3.8}

def test_consensus_returns_none_below_min():
    assert router.consensus({"dk": 0.30}, min_books=2, band=0.02) is None

def test_consensus_median_of_agreeing():
    bf = {"dk": 0.30, "fd": 0.31, "px": 0.50}  # px is an outlier
    out = router.consensus(bf, min_books=2, band=0.02)
    assert out == pytest.approx(0.305)  # median of dk+fd (px rejected)

def test_grid_spec_spread_total():
    legs = legset.parse_legs([
        {"market_ticker": "KXMLBSPREAD-25JUN271905NYYBOS-BOS2",
         "event_ticker": EVT, "side": "yes"},
        {"market_ticker": "KXMLBTOTAL-25JUN271905NYYBOS-9",
         "event_ticker": EVT, "side": "yes"}])
    family, spread_line, total_line, target = router.grid_spec(legs)
    assert spread_line == -1.5 and total_line == 8.5
    assert target == "Home Spread + Over"

def test_grid_cell_fairs_devigs_per_book():
    rows = (_grid_rows("dk", EVT, -1.5, 8.5, ST_CELLS)
            + _grid_rows("fd", EVT, -1.5, 8.5, ST_CELLS))
    df = pd.DataFrame(rows)
    from kalshi_common.leg_types import SPREAD_TOTAL_FAMILY
    out = router.grid_cell_fairs(EVT, SPREAD_TOTAL_FAMILY, -1.5, 8.5,
                                 "Home Spread + Over", df)
    assert set(out) == {"dk", "fd"}
    assert 0.20 < out["dk"] < 0.30   # devigged ~0.25-ish

ML_CELLS = {"Home ML + Over": 4.0, "Home ML + Under": 4.0,
            "Away ML + Over": 4.0, "Away ML + Under": 4.0}

def test_grid_spec_ml_total():
    legs = legset.parse_legs([
        {"market_ticker": "KXMLBGAME-25JUN271905NYYBOS-BOS",
         "event_ticker": EVT, "side": "yes"},          # home ML
        {"market_ticker": "KXMLBTOTAL-25JUN271905NYYBOS-9",
         "event_ticker": EVT, "side": "yes"}])         # over 8.5
    from kalshi_common.leg_types import ML_TOTAL_FAMILY
    family, spread_line, total_line, target = router.grid_spec(legs)
    assert family == ML_TOTAL_FAMILY
    assert spread_line is None and total_line == 8.5
    assert target == "Home ML + Over"

def test_grid_cell_fairs_ml_total_null_spread():
    # ml-total rows carry spread_line = None; isna() match must find them
    rows = (_grid_rows("dk", EVT, None, 8.5, ML_CELLS)
            + _grid_rows("fd", EVT, None, 8.5, ML_CELLS))
    df = pd.DataFrame(rows)
    from kalshi_common.leg_types import ML_TOTAL_FAMILY
    out = router.grid_cell_fairs(EVT, ML_TOTAL_FAMILY, None, 8.5,
                                 "Home ML + Over", df)
    assert set(out) == {"dk", "fd"}
    assert 0.0 < out["dk"] < 1.0

def test_single_marginal_total_over():
    # Over fair = Home Spread+Over + Away Spread+Over, devigged
    rows = (_grid_rows("dk", EVT, -1.5, 8.5, ST_CELLS)
            + _grid_rows("fd", EVT, -1.5, 8.5, ST_CELLS))
    df = pd.DataFrame(rows)
    leg = legset.CanonicalLeg(EVT, "total", 8.5, "over")
    out = router.single_marginal_fairs(EVT, leg, df)
    assert set(out) == {"dk", "fd"}
    assert 0.40 < out["dk"] < 0.60     # ~0.5 for this balanced grid

def test_single_marginal_spread_home():
    # Home Spread fair = Home Spread+Over + Home Spread+Under, devigged
    rows = (_grid_rows("dk", EVT, -1.5, 8.5, ST_CELLS)
            + _grid_rows("fd", EVT, -1.5, 8.5, ST_CELLS))
    df = pd.DataFrame(rows)
    leg = legset.CanonicalLeg(EVT, "spread", -1.5, "home")
    out = router.single_marginal_fairs(EVT, leg, df)
    assert set(out) == {"dk", "fd"}
    assert 0.0 < out["dk"] < 1.0

def test_single_marginal_ml_home():
    rows = _grid_rows("dk", EVT, None, 8.5, ML_CELLS)  # ml family => spread NULL
    df = pd.DataFrame(rows)
    leg = legset.CanonicalLeg(EVT, "ml", None, "home")
    out = router.single_marginal_fairs(EVT, leg, df)
    # Home ML = Home ML+Over + Home ML+Under, devigged; = 0.5 (balanced 4.0 grid)
    assert abs(out["dk"] - 0.5) < 1e-6

def test_single_marginal_missing_grid_returns_empty():
    df = pd.DataFrame(_grid_rows("dk", EVT, -1.5, 8.5, ST_CELLS))
    leg = legset.CanonicalLeg(EVT, "spread", -2.5, "home")  # no -2.5 grid
    assert router.single_marginal_fairs(EVT, leg, df) == {}


def test_subcombo_grid_spread_total():
    legs = legset.parse_legs([
        {"market_ticker": "KXMLBSPREAD-25JUN271905NYYBOS-BOS2",
         "event_ticker": EVT, "side": "yes"},
        {"market_ticker": "KXMLBTOTAL-25JUN271905NYYBOS-9",
         "event_ticker": EVT, "side": "yes"}])
    df = pd.DataFrame(_grid_rows("dk", EVT, -1.5, 8.5, ST_CELLS)
                      + _grid_rows("fd", EVT, -1.5, 8.5, ST_CELLS))
    out = router.subcombo_fair(EVT, legs, df, min_books=2, band=0.02)
    assert out is not None and 0.20 < out < 0.30


def test_subcombo_on_demand_returns_none_in_phase1():
    legs = legset.parse_legs([
        {"market_ticker": "KXMLBSPREAD-25JUN271905NYYBOS-BOS2",
         "event_ticker": EVT, "side": "yes"},
        {"market_ticker": "KXMLBTOTAL-25JUN271905NYYBOS-9",
         "event_ticker": EVT, "side": "yes"},
        {"market_ticker": "KXMLBGAME-25JUN271905NYYBOS-BOS",
         "event_ticker": EVT, "side": "yes"}])
    df = pd.DataFrame(_grid_rows("dk", EVT, -1.5, 8.5, ST_CELLS))
    assert router.subcombo_fair(EVT, legs, df, 2, 0.02) is None


def test_combo_fair_cross_game_multiplies():
    EVT2 = "KXMLBGAME-25JUN271905LADSF"
    legs_dicts = [
        {"market_ticker": "KXMLBGAME-25JUN271905NYYBOS-BOS",
         "event_ticker": EVT, "side": "yes"},                 # game1 home ML
        {"market_ticker": "KXMLBGAME-25JUN271905LADSF-LAD",
         "event_ticker": EVT2, "side": "yes"}]                # game2 away ML
    df = pd.DataFrame(
        _grid_rows("dk", EVT, None, 8.5, ML_CELLS)
        + _grid_rows("fd", EVT, None, 8.5, ML_CELLS)
        + _grid_rows("dk", EVT2, None, 8.5, ML_CELLS)
        + _grid_rows("fd", EVT2, None, 8.5, ML_CELLS))  # EVT2 = LADSF
    resolve = lambda game_legs: game_legs[0].game_id  # identity for the test
    out = router.combo_fair(legs_dicts, df, resolve, 2, 0.02)
    # both games priced & multiplied -> strictly less than either single fair
    g1 = router.subcombo_fair(EVT, [l for l in legset.parse_legs(legs_dicts)
                                    if l.game_id == EVT], df, 2, 0.02)
    assert out is not None and abs(out - g1 * g1) < 1e-9  # symmetric grids


def test_combo_fair_skips_when_a_game_unresolved():
    legs_dicts = [{"market_ticker": "KXMLBGAME-25JUN271905NYYBOS-BOS",
                   "event_ticker": EVT, "side": "yes"}]
    df = pd.DataFrame(_grid_rows("dk", EVT, None, 8.5, ML_CELLS)
                      + _grid_rows("fd", EVT, None, 8.5, ML_CELLS))
    assert router.combo_fair(legs_dicts, df, lambda gl: None, 2, 0.02) is None


def test_combo_fair_skips_untypeable():
    legs_dicts = [{"market_ticker": "KXMLBPLAYER-foo",
                   "event_ticker": EVT, "side": "yes"}]
    assert router.combo_fair(legs_dicts, pd.DataFrame(), lambda gl: EVT, 2, 0.02) is None
