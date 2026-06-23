"""Phase 2: general cross-game / N-leg combo pricer."""
import pandas as pd

from kalshi_mlb_mm import combo_pricer as cp


# ---- helpers to build synthetic book frames ---------------------------------

def _singles_df(rows):
    cols = ["game_id", "market", "line", "outcome", "decimal", "bookmaker"]
    return pd.DataFrame(rows, columns=cols)


def _ml_rows(game_id, p_home, books):
    """Fair-coin-ish: emit decimals that devig to ~p_home across `books` books."""
    out = []
    dh = 1.0 / p_home
    da = 1.0 / (1.0 - p_home)
    for b in books:
        out.append([game_id, "moneyline", None, "home", dh, b])
        out.append([game_id, "moneyline", None, "away", da, b])
    return out


def _resolver(mapping):
    return lambda suffix: mapping.get(suffix)


_BOOKS = ["dk", "fd", "bm", "br"]


def _ml_leg(suffix, team_code, home_code, side="yes"):
    return {"event_ticker": f"KXMLBGAME-{suffix}",
            "market_ticker": f"KXMLBGAME-{suffix}-{team_code}", "side": side,
            "_home_code": home_code}


# Use real-ish suffixes that _parse_event_suffix can decode (date + away+home).
SUF_A = "26MAY232205TEXLAA"   # away TEX, home LAA
SUF_B = "26MAY232205NYYBOS"   # away NYY, home BOS


def test_single_moneyline_home():
    singles = _singles_df(_ml_rows("gA", 0.60, _BOOKS))
    legs = [{"event_ticker": f"KXMLBGAME-{SUF_A}",
             "market_ticker": f"KXMLBGAME-{SUF_A}-LAA", "side": "yes"}]
    res = cp.combo_fair(legs, _resolver({SUF_A: "gA"}),
                        sgp_df=pd.DataFrame(), singles_df=singles,
                        min_agreeing=3, band=0.05)
    assert res is not None
    fair, agree = res
    assert abs(fair - 0.60) < 1e-6
    assert agree == 4


def test_two_game_moneyline_parlay_multiplies():
    singles = _singles_df(_ml_rows("gA", 0.60, _BOOKS) + _ml_rows("gB", 0.50, _BOOKS))
    legs = [
        {"event_ticker": f"KXMLBGAME-{SUF_A}", "market_ticker": f"KXMLBGAME-{SUF_A}-LAA", "side": "yes"},
        {"event_ticker": f"KXMLBGAME-{SUF_B}", "market_ticker": f"KXMLBGAME-{SUF_B}-BOS", "side": "yes"},
    ]
    res = cp.combo_fair(legs, _resolver({SUF_A: "gA", SUF_B: "gB"}),
                        sgp_df=pd.DataFrame(), singles_df=singles,
                        min_agreeing=3, band=0.05)
    assert res is not None
    fair, agree = res
    assert abs(fair - 0.30) < 1e-6   # 0.60 * 0.50


def test_away_moneyline_no_side_complement():
    singles = _singles_df(_ml_rows("gA", 0.60, _BOOKS))
    # away team TEX, side no -> P(TEX does NOT win) = P(home wins) = 0.60
    legs = [{"event_ticker": f"KXMLBGAME-{SUF_A}",
             "market_ticker": f"KXMLBGAME-{SUF_A}-TEX", "side": "no"}]
    fair, _ = cp.combo_fair(legs, _resolver({SUF_A: "gA"}),
                            pd.DataFrame(), singles, min_agreeing=3, band=0.05)
    assert abs(fair - 0.60) < 1e-6


def test_skips_when_game_unresolved():
    singles = _singles_df(_ml_rows("gA", 0.60, _BOOKS))
    legs = [{"event_ticker": f"KXMLBGAME-{SUF_A}",
             "market_ticker": f"KXMLBGAME-{SUF_A}-LAA", "side": "yes"}]
    res = cp.combo_fair(legs, _resolver({}), pd.DataFrame(), singles,
                        min_agreeing=3, band=0.05)
    assert res is None


def test_skips_when_too_few_books():
    singles = _singles_df(_ml_rows("gA", 0.60, ["dk", "fd"]))  # only 2 books
    legs = [{"event_ticker": f"KXMLBGAME-{SUF_A}",
             "market_ticker": f"KXMLBGAME-{SUF_A}-LAA", "side": "yes"}]
    res = cp.combo_fair(legs, _resolver({SUF_A: "gA"}), pd.DataFrame(), singles,
                        min_agreeing=3, band=0.05)
    assert res is None


def test_total_single_leg():
    rows = []
    for b in _BOOKS:
        rows.append(["gA", "total", 8.5, "over", 1.0 / 0.45, b])
        rows.append(["gA", "total", 8.5, "under", 1.0 / 0.55, b])
    singles = _singles_df(rows)
    # TotalLeg over 8.5 -> line_n = 9 (n-0.5 = 8.5). side yes = over.
    legs = [{"event_ticker": f"KXMLBTOTAL-{SUF_A}",
             "market_ticker": f"KXMLBTOTAL-{SUF_A}-9", "side": "yes"}]
    fair, _ = cp.combo_fair(legs, _resolver({SUF_A: "gA"}), pd.DataFrame(),
                            singles, min_agreeing=3, band=0.05)
    assert abs(fair - 0.45) < 1e-6


def test_spread_total_same_game_uses_sgp_grid():
    # Build a 4-cell SGP grid for one game/line across 4 books.
    grid = []
    # devigged target ~0.25 for "Home Spread + Over"
    for b in _BOOKS:
        grid.append(["gA", "Home Spread + Over", "FG", b, 1.0 / 0.25, -1.5, 8.5])
        grid.append(["gA", "Home Spread + Under", "FG", b, 1.0 / 0.25, -1.5, 8.5])
        grid.append(["gA", "Away Spread + Over", "FG", b, 1.0 / 0.25, -1.5, 8.5])
        grid.append(["gA", "Away Spread + Under", "FG", b, 1.0 / 0.25, -1.5, 8.5])
    sgp = pd.DataFrame(grid, columns=["game_id", "combo", "period", "bookmaker",
                                      "sgp_decimal", "spread_line", "total_line"])
    # legs: home spread -1.5 (N=2) yes + over 8.5 (N=9) yes -> "Home Spread + Over"
    legs = [
        {"event_ticker": f"KXMLBSPREAD-{SUF_A}", "market_ticker": f"KXMLBSPREAD-{SUF_A}-LAA2", "side": "yes"},
        {"event_ticker": f"KXMLBTOTAL-{SUF_A}", "market_ticker": f"KXMLBTOTAL-{SUF_A}-9", "side": "yes"},
    ]
    res = cp.combo_fair(legs, _resolver({SUF_A: "gA"}), sgp, _singles_df([]),
                        min_agreeing=3, band=0.05)
    assert res is not None
    fair, _ = res
    assert abs(fair - 0.25) < 1e-6


def _spread_rows(game_id, home_line, p_home_cover, books):
    rows = []
    dh = 1.0 / p_home_cover
    da = 1.0 / (1.0 - p_home_cover)
    for b in books:
        rows.append([game_id, "spread", home_line, "home", dh, b])
        rows.append([game_id, "spread", home_line, "away", da, b])
    return rows


def test_home_spread_single_main_line():
    # Home -1.5 (line_n=2), home favored: singles home-perspective line = -1.5.
    singles = _singles_df(_spread_rows("gA", -1.5, 0.55, _BOOKS))
    legs = [{"event_ticker": f"KXMLBSPREAD-{SUF_A}",
             "market_ticker": f"KXMLBSPREAD-{SUF_A}-LAA2", "side": "yes"}]
    fair, _ = cp.combo_fair(legs, _resolver({SUF_A: "gA"}), pd.DataFrame(),
                            singles, min_agreeing=3, band=0.05)
    assert abs(fair - 0.55) < 1e-6


def test_away_spread_single_uses_positive_home_line():
    # Away -1.5 (line_n=2) == home +1.5: singles home-perspective line = +1.5.
    # p_home_cover(+1.5)=0.40 -> P(away -1.5)=0.60.
    singles = _singles_df(_spread_rows("gA", 1.5, 0.40, _BOOKS))
    legs = [{"event_ticker": f"KXMLBSPREAD-{SUF_A}",
             "market_ticker": f"KXMLBSPREAD-{SUF_A}-TEX2", "side": "yes"}]
    fair, _ = cp.combo_fair(legs, _resolver({SUF_A: "gA"}), pd.DataFrame(),
                            singles, min_agreeing=3, band=0.05)
    assert abs(fair - 0.60) < 1e-6


def test_away_spread_alt_line_drops_not_misprices():
    # Away -1.5 leg but the singles only carry the home-favorite main line -1.5.
    # The away leg looks up +1.5 -> no rows -> combo dropped (NOT mis-priced
    # against the -1.5 main rows).
    singles = _singles_df(_spread_rows("gA", -1.5, 0.55, _BOOKS))
    legs = [{"event_ticker": f"KXMLBSPREAD-{SUF_A}",
             "market_ticker": f"KXMLBSPREAD-{SUF_A}-TEX2", "side": "yes"}]
    res = cp.combo_fair(legs, _resolver({SUF_A: "gA"}), pd.DataFrame(),
                        singles, min_agreeing=3, band=0.05)
    assert res is None


def _props_df(rows):
    cols = ["game_id", "market_type", "player", "line", "outcome", "decimal", "bookmaker"]
    return pd.DataFrame(rows, columns=cols)


def _hr_rows(game_id, player, p_over, books, line=0.5):
    out = []
    do = 1.0 / p_over
    du = 1.0 / (1.0 - p_over)
    for b in books:
        out.append([game_id, "home_runs", player, line, "over", do, b])
        out.append([game_id, "home_runs", player, line, "under", du, b])
    return out


def _prop_resolver(mapping):
    # mapping: market_ticker -> dict(market_type, player, line_n)
    return lambda mt: mapping.get(mt)


def test_prop_home_run_yes():
    props = _props_df(_hr_rows("gA", "cal raleigh", 0.30, _BOOKS))
    mt = f"KXMLBHR-{SUF_A}-SEACRALEIGH29-1"   # 1+ HR => Over 0.5
    legs = [{"event_ticker": f"KXMLBHR-{SUF_A}",
             "market_ticker": mt, "side": "yes"}]
    res = cp.combo_fair(legs, _resolver({SUF_A: "gA"}), pd.DataFrame(), _singles_df([]),
                        min_agreeing=3, band=0.05, props_df=props,
                        resolve_prop=_prop_resolver({mt: dict(market_type="home_runs",
                                                              player="cal raleigh", line_n=1)}))
    assert res is not None
    fair, _ = res
    assert abs(fair - 0.30) < 1e-6


def test_prop_no_side_complement():
    props = _props_df(_hr_rows("gA", "cal raleigh", 0.30, _BOOKS))
    mt = f"KXMLBHR-{SUF_A}-SEACRALEIGH29-1"
    legs = [{"event_ticker": f"KXMLBHR-{SUF_A}", "market_ticker": mt, "side": "no"}]
    fair, _ = cp.combo_fair(legs, _resolver({SUF_A: "gA"}), pd.DataFrame(), _singles_df([]),
                            min_agreeing=3, band=0.05, props_df=props,
                            resolve_prop=_prop_resolver({mt: dict(market_type="home_runs",
                                                                 player="cal raleigh", line_n=1)}))
    assert abs(fair - 0.70) < 1e-6   # P(0 HR) = 1 - P(1+)


def test_prop_dropped_when_unresolved():
    props = _props_df(_hr_rows("gA", "cal raleigh", 0.30, _BOOKS))
    mt = f"KXMLBHR-{SUF_A}-SEACRALEIGH29-1"
    legs = [{"event_ticker": f"KXMLBHR-{SUF_A}", "market_ticker": mt, "side": "yes"}]
    # resolve_prop returns None (e.g. GET /markets failed) -> combo dropped.
    res = cp.combo_fair(legs, _resolver({SUF_A: "gA"}), pd.DataFrame(), _singles_df([]),
                        min_agreeing=3, band=0.05, props_df=props,
                        resolve_prop=_prop_resolver({}))
    assert res is None


def test_prop_dropped_when_player_not_in_book_data():
    props = _props_df(_hr_rows("gA", "someone else", 0.30, _BOOKS))
    mt = f"KXMLBHR-{SUF_A}-SEACRALEIGH29-1"
    legs = [{"event_ticker": f"KXMLBHR-{SUF_A}", "market_ticker": mt, "side": "yes"}]
    res = cp.combo_fair(legs, _resolver({SUF_A: "gA"}), pd.DataFrame(), _singles_df([]),
                        min_agreeing=3, band=0.05, props_df=props,
                        resolve_prop=_prop_resolver({mt: dict(market_type="home_runs",
                                                              player="cal raleigh", line_n=1)}))
    assert res is None


def test_prop_plus_moneyline_cross_game_parlay():
    # game A: moneyline (singles); game B: HR prop. Independent -> product.
    singles = _singles_df(_ml_rows("gA", 0.60, _BOOKS))
    props = _props_df(_hr_rows("gB", "cal raleigh", 0.30, _BOOKS))
    mtp = f"KXMLBHR-{SUF_B}-BOSCRALEIGH29-1"
    legs = [
        {"event_ticker": f"KXMLBGAME-{SUF_A}", "market_ticker": f"KXMLBGAME-{SUF_A}-LAA", "side": "yes"},
        {"event_ticker": f"KXMLBHR-{SUF_B}", "market_ticker": mtp, "side": "yes"},
    ]
    fair, _ = cp.combo_fair(legs, _resolver({SUF_A: "gA", SUF_B: "gB"}), pd.DataFrame(),
                            singles, min_agreeing=3, band=0.05, props_df=props,
                            resolve_prop=_prop_resolver({mtp: dict(market_type="home_runs",
                                                                  player="cal raleigh", line_n=1)}))
    assert abs(fair - 0.18) < 1e-6   # 0.60 * 0.30


def test_prop_plus_same_game_leg_dropped():
    # HR prop + moneyline in the SAME game -> correlated, no joint price -> drop.
    props = _props_df(_hr_rows("gA", "cal raleigh", 0.30, _BOOKS))
    singles = _singles_df(_ml_rows("gA", 0.60, _BOOKS))
    mt = f"KXMLBHR-{SUF_A}-LAACRALEIGH29-1"
    legs = [
        {"event_ticker": f"KXMLBGAME-{SUF_A}", "market_ticker": f"KXMLBGAME-{SUF_A}-LAA", "side": "yes"},
        {"event_ticker": f"KXMLBHR-{SUF_A}", "market_ticker": mt, "side": "yes"},
    ]
    res = cp.combo_fair(legs, _resolver({SUF_A: "gA"}), pd.DataFrame(), singles,
                        min_agreeing=3, band=0.05, props_df=props,
                        resolve_prop=_prop_resolver({mt: dict(market_type="home_runs",
                                                             player="cal raleigh", line_n=1)}))
    assert res is None


def test_skips_threeway_same_game():
    # 3 legs same game that aren't a clean priceable shape -> skip whole combo.
    legs = [
        {"event_ticker": f"KXMLBGAME-{SUF_A}", "market_ticker": f"KXMLBGAME-{SUF_A}-LAA", "side": "yes"},
        {"event_ticker": f"KXMLBSPREAD-{SUF_A}", "market_ticker": f"KXMLBSPREAD-{SUF_A}-LAA2", "side": "yes"},
        {"event_ticker": f"KXMLBTOTAL-{SUF_A}", "market_ticker": f"KXMLBTOTAL-{SUF_A}-9", "side": "yes"},
    ]
    res = cp.combo_fair(legs, _resolver({SUF_A: "gA"}), pd.DataFrame(),
                        _singles_df([]), min_agreeing=3, band=0.05)
    assert res is None
