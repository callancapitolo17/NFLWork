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

ML_CELLS = {"Home ML + Over": 3.6, "Home ML + Under": 4.0,
            "Away ML + Over": 4.4, "Away ML + Under": 4.0}

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
