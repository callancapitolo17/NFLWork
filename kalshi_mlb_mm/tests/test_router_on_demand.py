"""Phase 2: router on_demand seam (grid path must stay byte-identical)."""
import pandas as pd

from kalshi_common import legset
from kalshi_mlb_mm import router

EVT = "KXMLBGAME-25JUN271905NYYBOS"


def _legs3():
    return [legset.CanonicalLeg(EVT, "spread", -1.5, "home"),
            legset.CanonicalLeg(EVT, "total", 8.5, "over"),
            legset.CanonicalLeg(EVT, "ml", None, "home")]


def test_on_demand_route_prices_via_injected_lookup():
    legs = _legs3()
    h = legset.leg_set_hash(legs)
    fairs = {"draftkings": 0.140, "fanduel": 0.142, "betmgm": 0.145}
    lookup = lambda hash_: fairs if hash_ == h else None
    cons, reason = router.subcombo_consensus("g1", legs, sgp_df=None,
                                             min_books=2, sigma_z_max=0.07,
                                             on_demand_fairs=lookup)
    assert reason == "ok"
    assert abs(cons.fair - 0.142) < 1e-12    # median of ALL books


def test_on_demand_route_none_without_lookup():
    # Default arg: Phase 1 behavior (on_demand -> unpriceable)
    cons, reason = router.subcombo_consensus("g1", _legs3(), sgp_df=None,
                                             min_books=2, sigma_z_max=0.07)
    assert cons is None and reason == "unpriceable"


def test_on_demand_route_none_when_lookup_misses():
    lookup = lambda hash_: None
    cons, reason = router.subcombo_consensus("g1", _legs3(), sgp_df=None,
                                             min_books=2, sigma_z_max=0.07,
                                             on_demand_fairs=lookup)
    assert cons is None and reason == "too_few_books"


def test_on_demand_route_consensus_gate_applies():
    legs = _legs3()
    lookup = lambda hash_: {"draftkings": 0.14}       # 1 book < min_books=2
    cons, reason = router.subcombo_consensus("g1", legs, sgp_df=None,
                                             min_books=2, sigma_z_max=0.07,
                                             on_demand_fairs=lookup)
    assert cons is None and reason == "too_few_books"


def test_combo_fair_threads_lookup_through_multi_game():
    # game A: 2-leg grid (from sgp_df), game B: 3-leg on_demand (from lookup)
    evt_b = "KXMLBGAME-25JUN272005SDLAD"
    legs_a = [{"market_ticker": f"KXMLBSPREAD-25JUN271905NYYBOS-BOS2",
               "event_ticker": EVT, "side": "yes"},
              {"market_ticker": f"KXMLBTOTAL-25JUN271905NYYBOS-9",
               "event_ticker": EVT, "side": "yes"}]
    legs_b = [{"market_ticker": "KXMLBSPREAD-25JUN272005SDLAD-LAD2",
               "event_ticker": evt_b, "side": "yes"},
              {"market_ticker": "KXMLBTOTAL-25JUN272005SDLAD-8",
               "event_ticker": evt_b, "side": "yes"},
              {"market_ticker": "KXMLBGAME-25JUN272005SDLAD-LAD",
               "event_ticker": evt_b, "side": "yes"}]
    # grid rows for game A: full 4-cell grid, two books
    rows = []
    for book in ("draftkings", "fanduel"):
        for combo, dec in (("Home Spread + Over", 3.10), ("Home Spread + Under", 4.20),
                           ("Away Spread + Over", 3.60), ("Away Spread + Under", 3.05)):
            rows.append(dict(game_id="gA", combo=combo, period="FG",
                             bookmaker=book, sgp_decimal=dec,
                             spread_line=-1.5, total_line=8.5))
    sgp_df = pd.DataFrame(rows)
    canon_b = legset.parse_legs(legs_b)
    hash_b = legset.leg_set_hash(canon_b)
    lookup = lambda h: {"draftkings": 0.20, "fanduel": 0.21} if h == hash_b else None
    resolve = lambda gl: "gA" if gl[0].game_id == EVT else "gB"
    got = router.combo_fair(legs_a + legs_b, sgp_df, resolve, 2, 0.07,
                            on_demand_fairs=lookup)
    cons_a, _ = router.subcombo_consensus("gA", legset.parse_legs(legs_a),
                                          sgp_df, 2, 0.07)
    assert got is not None and cons_a is not None
    assert abs(got - cons_a.fair * 0.205) < 1e-12


def test_consensus_detail_matches_consensus_and_names_books():
    fairs = {"draftkings": 0.140, "fanduel": 0.142, "betmgm": 0.145}
    det = router.consensus_detail(fairs, min_books=2, sigma_z_max=0.07)
    cons, _ = router.consensus(fairs, min_books=2, sigma_z_max=0.07)
    assert det is not None and cons is not None
    fair, books = det
    assert abs(fair - cons.fair) < 1e-12
    assert sorted(books) == ["betmgm", "draftkings", "fanduel"]  # ALL books
    # dispersion decline and count decline both -> None
    outlier = {"draftkings": 0.140, "fanduel": 0.142, "novig": 0.30}
    assert router.consensus_detail(outlier, min_books=2, sigma_z_max=0.07) is None
    assert router.consensus_detail({"dk": 0.5}, min_books=2, sigma_z_max=0.07) is None
