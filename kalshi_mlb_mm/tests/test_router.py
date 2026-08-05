import pandas as pd
import pytest
from kalshi_common import legset
from kalshi_mlb_mm import router
from kalshi_mlb_mm.tests.conftest import leg


EVT = "25JUN271905NYYBOS"

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
    cons, reason = router.consensus({"dk": 0.30}, min_books=2, sigma_z_max=0.07)
    assert cons is None and reason == "too_few_books"

def test_consensus_declines_on_dispersion():
    # Issue #20: a loud dissenter is no longer outvoted — the whole combo
    # declines (the outlier may be the informed book).
    bf = {"dk": 0.30, "fd": 0.31, "px": 0.50}
    cons, reason = router.consensus(bf, min_books=2, sigma_z_max=0.07)
    assert cons is None and reason == "consensus_dispersion"

def test_consensus_median_of_all_books():
    bf = {"dk": 0.30, "fd": 0.31, "px": 0.32}   # tight cluster -> quote
    cons, reason = router.consensus(bf, min_books=2, sigma_z_max=0.07)
    assert reason == "ok"
    assert cons.fair == pytest.approx(0.31)      # median of ALL books
    assert cons.n_books == 3 and cons.sigma_pts == pytest.approx(0.01)

def test_grid_spec_spread_total():
    legs = legset.parse_legs([
        leg("KXMLBSPREAD-25JUN271905NYYBOS-BOS2", "yes"),
        leg("KXMLBTOTAL-25JUN271905NYYBOS-9", "yes")])
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

def test_grid_spec_away_favourite_positive_line():
    # Issue #70: NYY is the away team; its margin ticker is away -1.5, which is
    # +1.5 home-perspective. Before the fix this collapsed onto spread_line
    # -1.5, whose "Away Spread" cell is away +1.5 — a different, easier event.
    legs = legset.parse_legs([
        leg("KXMLBSPREAD-25JUN271905NYYBOS-NYY2", "yes"),
        leg("KXMLBTOTAL-25JUN271905NYYBOS-9", "yes")])
    family, spread_line, total_line, target = router.grid_spec(legs)
    assert spread_line == 1.5 and total_line == 8.5
    assert target == "Away Spread + Over"


# Two DISTINGUISHABLE grids for the same game: the home-favourite (-1.5) grid
# devigs every cell to 0.25; the away-favourite (+1.5) grid devigs its Away
# cells to 0.30 and its Home cells to 0.20. A test asserting ~0.30/0.60 can
# therefore only pass if the +1.5 grid was selected.
NEG_GRID_CELLS = {"Home Spread + Over": 4.0, "Home Spread + Under": 4.0,
                  "Away Spread + Over": 4.0, "Away Spread + Under": 4.0}
POS_GRID_CELLS = {"Home Spread + Over": 5.0, "Home Spread + Under": 5.0,
                  "Away Spread + Over": 10.0 / 3.0, "Away Spread + Under": 10.0 / 3.0}


def _both_sign_grids_df():
    return pd.DataFrame(
        _grid_rows("dk", EVT, -1.5, 8.5, NEG_GRID_CELLS)
        + _grid_rows("fd", EVT, -1.5, 8.5, NEG_GRID_CELLS)
        + _grid_rows("dk", EVT, 1.5, 8.5, POS_GRID_CELLS)
        + _grid_rows("fd", EVT, 1.5, 8.5, POS_GRID_CELLS))


def test_away_favourite_combo_prices_from_positive_grid():
    # away -1.5 + Over 8.5 must price from the +1.5 grid's "Away Spread + Over"
    # (~0.30), NOT the -1.5 grid's away +1.5 cell (~0.25).
    legs = legset.parse_legs([
        leg("KXMLBSPREAD-25JUN271905NYYBOS-NYY2", "yes"),
        leg("KXMLBTOTAL-25JUN271905NYYBOS-9", "yes")])
    cons, reason = router.subcombo_consensus(EVT, legs, _both_sign_grids_df(),
                                             min_books=2, sigma_z_max=0.07)
    assert reason == "ok"
    assert cons.fair == pytest.approx(0.30, abs=0.02)


def test_away_favourite_single_marginalizes_positive_grid():
    # The SINGLE route (most common in the corpus): a lone away -1.5 leg must
    # marginalize the +1.5 grid's Away cells (~0.60), not the -1.5 grid's
    # away +1.5 cells (~0.50).
    canon = legset.parse_legs([leg("KXMLBSPREAD-25JUN271905NYYBOS-NYY2", "yes")])
    out = router.single_marginal_fairs(EVT, canon[0], _both_sign_grids_df())
    assert set(out) == {"dk", "fd"}
    assert out["dk"] == pytest.approx(0.60, abs=0.02)


def test_home_underdog_leg_prices_from_away_favourite_grid():
    # NYY2 NO == home +1.5: lives on the AWAY team's margin market, so it must
    # look up the +1.5 grid's Home cells (~0.40), not the -1.5 grid's Home
    # cells (~0.50).
    canon = legset.parse_legs([leg("KXMLBSPREAD-25JUN271905NYYBOS-NYY2", "no")])
    out = router.single_marginal_fairs(EVT, canon[0], _both_sign_grids_df())
    assert set(out) == {"dk", "fd"}
    assert out["dk"] == pytest.approx(0.40, abs=0.02)


def test_grid_spec_ml_total():
    legs = legset.parse_legs([
        leg("KXMLBGAME-25JUN271905NYYBOS-BOS", "yes"),          # home ML
        leg("KXMLBTOTAL-25JUN271905NYYBOS-9", "yes")])         # over 8.5
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
    # Home ML = Home ML+Over + Home ML+Under, devigged; asymmetric grid yields ~0.525
    assert 0.50 < out["dk"] < 0.55

def test_single_marginal_missing_grid_returns_empty():
    df = pd.DataFrame(_grid_rows("dk", EVT, -1.5, 8.5, ST_CELLS))
    leg = legset.CanonicalLeg(EVT, "spread", -2.5, "home")  # no -2.5 grid
    assert router.single_marginal_fairs(EVT, leg, df) == {}


def test_subcombo_grid_spread_total():
    legs = legset.parse_legs([
        leg("KXMLBSPREAD-25JUN271905NYYBOS-BOS2", "yes"),
        leg("KXMLBTOTAL-25JUN271905NYYBOS-9", "yes")])
    df = pd.DataFrame(_grid_rows("dk", EVT, -1.5, 8.5, ST_CELLS)
                      + _grid_rows("fd", EVT, -1.5, 8.5, ST_CELLS))
    cons, reason = router.subcombo_consensus(EVT, legs, df, min_books=2,
                                             sigma_z_max=0.07)
    assert reason == "ok" and 0.20 < cons.fair < 0.30


def test_subcombo_on_demand_returns_none_in_phase1():
    legs = legset.parse_legs([
        leg("KXMLBSPREAD-25JUN271905NYYBOS-BOS2", "yes"),
        leg("KXMLBTOTAL-25JUN271905NYYBOS-9", "yes"),
        leg("KXMLBGAME-25JUN271905NYYBOS-BOS", "yes")])
    df = pd.DataFrame(_grid_rows("dk", EVT, -1.5, 8.5, ST_CELLS))
    cons, reason = router.subcombo_consensus(EVT, legs, df, 2, 0.07)
    assert cons is None and reason == "unpriceable"


def test_combo_fair_cross_game_multiplies():
    EVT2 = "25JUN271905LADSF"
    legs_dicts = [
        leg("KXMLBGAME-25JUN271905NYYBOS-BOS", "yes"),                 # game1 home ML
        leg("KXMLBGAME-25JUN271905LADSF-LAD", "yes")]                # game2 away ML
    df = pd.DataFrame(
        _grid_rows("dk", EVT, None, 8.5, ML_CELLS)
        + _grid_rows("fd", EVT, None, 8.5, ML_CELLS)
        + _grid_rows("dk", EVT2, None, 8.5, ML_CELLS)
        + _grid_rows("fd", EVT2, None, 8.5, ML_CELLS))  # EVT2 = LADSF
    resolve = lambda game_legs: game_legs[0].game_id  # identity for the test
    out = router.combo_fair(legs_dicts, df, resolve, 2, 0.07)
    # both games priced & multiplied -> product of the two per-game fairs
    canon = legset.parse_legs(legs_dicts)
    g_home, _ = router.subcombo_consensus(
        EVT, [l for l in canon if l.game_id == EVT], df, 2, 0.07)
    g_away, _ = router.subcombo_consensus(
        EVT2, [l for l in canon if l.game_id == EVT2], df, 2, 0.07)
    assert out is not None
    assert g_home is not None and g_away is not None
    assert abs(out - g_home.fair * g_away.fair) < 1e-9


def test_combo_fair_detail_reports_n_games_and_sigma():
    # Same 2-game setup as above: detail must carry n_games=2 and a
    # non-negative combo sigma (identical books here -> sigma 0).
    EVT2 = "25JUN271905LADSF"
    legs_dicts = [
        leg("KXMLBGAME-25JUN271905NYYBOS-BOS", "yes"),
        leg("KXMLBGAME-25JUN271905LADSF-LAD", "yes")]
    df = pd.DataFrame(
        _grid_rows("dk", EVT, None, 8.5, ML_CELLS)
        + _grid_rows("fd", EVT, None, 8.5, ML_CELLS)
        + _grid_rows("dk", EVT2, None, 8.5, ML_CELLS)
        + _grid_rows("fd", EVT2, None, 8.5, ML_CELLS))
    detail, reason = router.combo_fair_detail(
        legs_dicts, df, lambda gl: gl[0].game_id, 2, 0.07)
    assert reason == "ok"
    assert detail.n_games == 2
    assert detail.sigma_pts == pytest.approx(0.0, abs=1e-12)


def test_combo_fair_skips_when_a_game_unresolved():
    legs_dicts = [leg("KXMLBGAME-25JUN271905NYYBOS-BOS", "yes")]
    df = pd.DataFrame(_grid_rows("dk", EVT, None, 8.5, ML_CELLS)
                      + _grid_rows("fd", EVT, None, 8.5, ML_CELLS))
    detail, reason = router.combo_fair_detail(legs_dicts, df, lambda gl: None,
                                              2, 0.07)
    assert detail is None and reason == "unresolved_game"


def test_combo_fair_skips_untypeable():
    legs_dicts = [leg("KXMLBPLAYER-foo", "yes")]
    detail, reason = router.combo_fair_detail(legs_dicts, pd.DataFrame(),
                                              lambda gl: EVT, 2, 0.07)
    assert detail is None and reason == "unparseable"


# FIX A: _priceable_in_phase1 must reject lone single-leg RFQs (< 2 legs).
def test_priceable_in_phase1_rejects_lone_single():
    from kalshi_mlb_mm import main
    one_leg = legset.parse_legs([leg("KXMLBGAME-25JUN271905NYYBOS-BOS", "yes")])
    assert one_leg is not None and len(one_leg) == 1
    assert main._priceable_in_phase1(one_leg) is False


def test_priceable_in_phase1_accepts_two_leg_combo():
    from kalshi_mlb_mm import main
    two_legs = legset.parse_legs([
        leg("KXMLBSPREAD-25JUN271905NYYBOS-BOS2", "yes"),
        leg("KXMLBTOTAL-25JUN271905NYYBOS-9", "yes"),
    ])
    assert two_legs is not None and len(two_legs) == 2
    assert main._priceable_in_phase1(two_legs) is True
