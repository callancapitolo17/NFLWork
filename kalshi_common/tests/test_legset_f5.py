"""Issue #84 regression tests: F5 leg parsing — period-aware CanonicalLeg,
on_demand routing, FG-only grid guard, F5-winner unpriceable guard.

Ticker grammar verified against the live Kalshi API (2026-08-11):

- ``KXMLBF5TOTAL-26AUG122210KCLAD-7`` has floor_strike 6.5 ("Over 6.5 runs in
  the first 5 innings") => integer suffix N maps to line N - 0.5, same as FG.
- ``KXMLBF5SPREAD-26AUG122210KCLAD-LAD3`` is "Los Angeles D -2.5 first 5
  innings" => same ``-{TEAM}{N}`` grammar and sign semantics as KXMLBSPREAD.
- ``KXMLBF5-26AUG122210KCLAD-{KC,LAD,TIE}`` — THREE markets per event: the
  F5 winner family is 3-way (team markets resolve NO on a tie). This is the
  #86 tie trap; F5-winner legs parse but stay unpriceable until #86.
"""
from kalshi_common import legset
from kalshi_common.leg_types import combo_descriptor

GAME = "26AUG122210KCLAD"            # away=KC, home=LAD (live 2026-08-11 slate)
EVT = GAME


def _fg_spread(team, n, side):
    return {"market_ticker": f"KXMLBSPREAD-{GAME}-{team}{n}",
            "event_ticker": f"KXMLBSPREAD-{GAME}", "side": side}


def _fg_total(n, side):
    return {"market_ticker": f"KXMLBTOTAL-{GAME}-{n}",
            "event_ticker": f"KXMLBTOTAL-{GAME}", "side": side}


def _fg_ml(team, side):
    return {"market_ticker": f"KXMLBGAME-{GAME}-{team}",
            "event_ticker": f"KXMLBGAME-{GAME}", "side": side}


def _f5_spread(team, n, side):
    return {"market_ticker": f"KXMLBF5SPREAD-{GAME}-{team}{n}",
            "event_ticker": f"KXMLBF5SPREAD-{GAME}", "side": side}


def _f5_total(n, side):
    return {"market_ticker": f"KXMLBF5TOTAL-{GAME}-{n}",
            "event_ticker": f"KXMLBF5TOTAL-{GAME}", "side": side}


def _f5_winner(team, side):
    return {"market_ticker": f"KXMLBF5-{GAME}-{team}",
            "event_ticker": f"KXMLBF5-{GAME}", "side": side}


# --------------------------------------------------------------------- #
# parse_leg: all three F5 series, both sides                             #
# --------------------------------------------------------------------- #

def test_parse_f5_spread_home_minus_1_5():
    leg = legset.parse_leg(_f5_spread("LAD", 2, "yes"))   # home -1.5 F5
    assert leg == legset.CanonicalLeg(EVT, "spread", -1.5, "home", "F5")


def test_parse_f5_spread_home_no_flips_to_away():
    leg = legset.parse_leg(_f5_spread("LAD", 2, "no"))
    assert leg == legset.CanonicalLeg(EVT, "spread", -1.5, "away", "F5")


def test_parse_f5_spread_away_team_yes():
    # KC is away: away margin ticker -> positive home-perspective line (#70)
    leg = legset.parse_leg(_f5_spread("KC", 2, "yes"))
    assert leg == legset.CanonicalLeg(EVT, "spread", 1.5, "away", "F5")


def test_parse_f5_spread_away_team_no():
    leg = legset.parse_leg(_f5_spread("KC", 2, "no"))
    assert leg == legset.CanonicalLeg(EVT, "spread", 1.5, "home", "F5")


def test_parse_f5_total_over_6_5():
    # Live-verified: suffix 7 <=> floor_strike 6.5 (N - 0.5 convention)
    leg = legset.parse_leg(_f5_total(7, "yes"))
    assert leg == legset.CanonicalLeg(EVT, "total", 6.5, "over", "F5")


def test_parse_f5_total_under_6_5():
    leg = legset.parse_leg(_f5_total(7, "no"))
    assert leg == legset.CanonicalLeg(EVT, "total", 6.5, "under", "F5")


def test_parse_f5_winner_yes_is_half_run_line_favorite():
    """#86 rework: 'team wins F5' IS 'team -0.5 on the F5 run line' (win by
    >=1; a tie loses both). Team-winner legs re-encode as spread legs at
    +-0.5, which books price two-sided with NO push — the tie mass is
    already inside the +0.5 side, so no conversion is ever needed."""
    assert legset.parse_leg(_f5_winner("LAD", "yes")) == \
        legset.CanonicalLeg(EVT, "spread", -0.5, "home", "F5")
    assert legset.parse_leg(_f5_winner("KC", "yes")) == \
        legset.CanonicalLeg(EVT, "spread", 0.5, "away", "F5")


def test_parse_f5_winner_no_is_other_team_plus_half():
    """NO on a team market pays on 'other team wins OR tie' — which is
    exactly the other team +0.5 (this is what the ml-complement encoding
    got wrong: it dropped the tie mass a NO holder wins)."""
    assert legset.parse_leg(_f5_winner("LAD", "no")) == \
        legset.CanonicalLeg(EVT, "spread", -0.5, "away", "F5")
    assert legset.parse_leg(_f5_winner("KC", "no")) == \
        legset.CanonicalLeg(EVT, "spread", 0.5, "home", "F5")


def test_parse_f5_winner_tie_market_is_not_a_team_side():
    """KXMLBF5-...-TIE exists live (3-way family). It must never mislabel as
    an away ML — that would be a silently wrong leg in telemetry AND a
    mispriced side if it ever reached pricing."""
    tie_yes = legset.parse_leg(_f5_winner("TIE", "yes"))
    tie_no = legset.parse_leg(_f5_winner("TIE", "no"))
    assert tie_yes is not None and tie_no is not None
    assert tie_yes.side not in ("home", "away")
    assert tie_no.side not in ("home", "away")
    assert tie_yes.side != tie_no.side
    assert tie_yes.period == "F5" and tie_yes.market_type == "ml"


def test_fg_legs_default_period_fg():
    leg = legset.parse_leg(_fg_spread("LAD", 2, "yes"))
    assert leg.period == "FG"
    # defaulted trailing field: 4-positional construction stays valid
    assert leg == legset.CanonicalLeg(EVT, "spread", -1.5, "home")


# --------------------------------------------------------------------- #
# hash distinctness + duplicate-market guard                             #
# --------------------------------------------------------------------- #

def test_hash_distinguishes_fg_and_f5_same_line():
    """FG spread -1.5 and F5 spread -1.5 agree on (market_type, line, side) —
    only period separates them. Their hashes must differ."""
    fg = legset.parse_legs([_fg_spread("LAD", 2, "yes")])
    f5 = legset.parse_legs([_f5_spread("LAD", 2, "yes")])
    assert legset.leg_set_hash(fg) != legset.leg_set_hash(f5)


def test_dup_guard_allows_fg_and_f5_totals_same_game():
    legs = legset.parse_legs([_fg_total(9, "yes"), _f5_total(5, "yes")])
    assert legs is not None
    assert legset.classify_subcombo(legs) == "on_demand"


def test_dup_guard_allows_fg_and_f5_spreads_at_same_line():
    """The sharpened guard key (market_type, period, line): FG -1.5 with
    F5 -1.5 is NOT a contradiction — different periods, both can win."""
    legs = legset.parse_legs([_fg_spread("LAD", 2, "yes"),
                              _f5_spread("LAD", 2, "yes")])
    assert legs is not None
    assert legset.classify_subcombo(legs) == "on_demand"


def test_dup_guard_rejects_two_f5_totals_same_line():
    over = legset.parse_leg(_f5_total(5, "yes"))
    under = legset.parse_leg(_f5_total(5, "no"))
    assert legset.classify_subcombo([over, under]) == "unpriceable"
    assert legset.classify_subcombo([over, over]) == "unpriceable"


# --------------------------------------------------------------------- #
# classify_subcombo routing                                              #
# --------------------------------------------------------------------- #

def test_f5_spread_total_same_game_routes_on_demand():
    legs = legset.parse_legs([_f5_spread("LAD", 2, "yes"), _f5_total(7, "yes")])
    assert legs is not None
    by_game = legset.partition_by_game(legs)
    assert len(by_game) == 1
    assert legset.classify_subcombo(next(iter(by_game.values()))) == "on_demand"


def test_mixed_fg_f5_same_game_never_grids():
    """Grids are FG-only: an FG spread + F5 total 2-leg set must not classify
    grid_spread_total (the grid tables are full-game surfaces)."""
    legs = legset.parse_legs([_fg_spread("LAD", 2, "yes"), _f5_total(7, "yes")])
    assert legset.classify_subcombo(legs) == "on_demand"
    legs = legset.parse_legs([_fg_ml("LAD", "yes"), _f5_total(7, "yes")])
    assert legset.classify_subcombo(legs) == "on_demand"


def test_fg_grids_unchanged():
    legs = legset.parse_legs([_fg_spread("LAD", 2, "yes"), _fg_total(9, "yes")])
    assert legset.classify_subcombo(legs) == "grid_spread_total"
    legs = legset.parse_legs([_fg_ml("LAD", "yes"), _fg_total(9, "yes")])
    assert legset.classify_subcombo(legs) == "grid_ml_total"


GAME2 = "26AUG121940CINCWS"          # away=CIN, home=CWS (live slate)


def test_cross_game_f5_single_plus_fg_single_routes_per_game_singles():
    f5_leg = {"market_ticker": f"KXMLBF5SPREAD-{GAME2}-CWS2",
              "event_ticker": f"KXMLBF5SPREAD-{GAME2}", "side": "yes"}
    legs = legset.parse_legs([f5_leg, _fg_ml("LAD", "yes")])
    assert legs is not None
    by_game = legset.partition_by_game(legs)
    assert set(by_game) == {GAME, GAME2}
    assert all(legset.classify_subcombo(gl) == "single"
               for gl in by_game.values())


# --------------------------------------------------------------------- #
# F5 winner: team legs price as +-0.5 run lines (#86); TIE legs decline  #
# --------------------------------------------------------------------- #

def test_f5_winner_combo_routes_on_demand_same_game():
    """#86: F5-winner TEAM legs are priceable — re-encoded as +-0.5 spread
    legs, the multi-leg set routes on_demand like other F5 shapes."""
    legs = legset.parse_legs([_f5_winner("LAD", "yes"), _f5_total(7, "yes")])
    assert legs is not None
    assert legset.classify_subcombo(legs) == "on_demand"


def test_f5_winner_lone_partition_leg_routes_single():
    """A cross-game combo puts the F5 winner alone in its game's partition —
    it prices as a marginalized single (live routing fetches the 1-leg set
    on-demand as a +-0.5 run-line single)."""
    leg = legset.parse_leg(_f5_winner("LAD", "yes"))
    assert legset.classify_subcombo([leg]) == "single"


def test_f5_winner_tie_leg_is_unpriceable():
    """No book prices a bare F5-tie leg in an SGP — TIE-side legs stay
    fail-closed (they're the one leg the +-0.5 re-encoding can't express)."""
    tie = legset.parse_leg(_f5_winner("TIE", "yes"))
    assert legset.classify_subcombo([tie, legset.parse_leg(_f5_total(7, "yes"))]) \
        == "unpriceable"
    assert legset.classify_subcombo([tie]) == "unpriceable"


def test_f5_winner_home_plus_away_dies_on_contradiction_guard():
    """Home F5 winner + away F5 winner is a joint-probability-zero pair.
    Pre-rework both keyed as ("ml","F5",None) and the dup guard caught it;
    re-encoded they sit at DIFFERENT keys ((-0.5,home) / (+0.5,away)), so
    the spread contradiction guard must catch it instead — a book that
    product-priced it would hand RFQ creators a worthless-YES pick-off."""
    legs = [legset.parse_leg(_f5_winner("LAD", "yes")),
            legset.parse_leg(_f5_winner("KC", "yes"))]
    assert legset.classify_subcombo(legs) == "unpriceable"


def test_f5_winner_duplicates_equivalent_f5spread_ticker():
    """KXMLBF5SPREAD at N=1 encodes the SAME +-0.5 market — combining it
    with the winner ticker is the same leg twice (dup guard)."""
    legs = [legset.parse_leg(_f5_winner("LAD", "yes")),
            legset.parse_leg(_f5_spread("LAD", 1, "yes"))]
    assert legset.classify_subcombo(legs) == "unpriceable"


def test_spread_contradiction_guard_catches_opposed_favorites_fg_too():
    """The guard is general: FG home -1.5 + FG away -1.5 (both favorites,
    margin >=2 both ways) is impossible and previously slipped past the
    same-key dup guard."""
    legs = [legset.parse_leg(_fg_spread("LAD", 2, "yes")),
            legset.parse_leg(_fg_spread("KC", 2, "yes"))]
    assert legset.classify_subcombo(legs) == "unpriceable"


def test_nested_same_side_spreads_still_allowed():
    """Home -0.5 with home -2.5 nests (both can win) — must NOT trip the
    contradiction guard."""
    legs = [legset.parse_leg(_f5_winner("LAD", "yes")),
            legset.parse_leg(_f5_spread("LAD", 3, "yes"))]
    assert legset.classify_subcombo(legs) == "on_demand"


def test_total_contradiction_guard_catches_impossible_band():
    """Over 8.5 + Under 4.5 (same period) is empty; Over 4.5 + Under 8.5
    is a legitimate band and stays priceable."""
    impossible = [legset.parse_leg(_fg_total(9, "yes")),
                  legset.parse_leg(_fg_total(5, "no"))]
    assert legset.classify_subcombo(impossible) == "unpriceable"
    band = [legset.parse_leg(_fg_total(5, "yes")),
            legset.parse_leg(_fg_total(9, "no"))]
    assert legset.classify_subcombo(band) == "on_demand"


def test_cross_period_bounds_do_not_interact():
    """F5 and FG margins are different variables — F5 home -0.5 with FG
    away -1.5 is a legitimate correlated combo, not a contradiction."""
    legs = [legset.parse_leg(_f5_winner("LAD", "yes")),
            legset.parse_leg(_fg_spread("KC", 2, "yes"))]
    assert legset.classify_subcombo(legs) == "on_demand"


# --------------------------------------------------------------------- #
# flip_leg / enumerate_partition preserve period                         #
# --------------------------------------------------------------------- #

def test_flip_leg_round_trips_f5_leg():
    leg = legset.parse_leg(_f5_spread("LAD", 2, "yes"))
    flipped = legset.flip_leg(leg)
    assert flipped.period == "F5" and flipped.side == "away"
    assert legset.flip_leg(flipped) == leg


def test_enumerate_partition_preserves_period_in_every_cell():
    legs = [legset.parse_leg(_f5_spread("LAD", 2, "yes")),
            legset.parse_leg(_f5_total(7, "yes"))]
    for cell in legset.enumerate_partition(legs):
        assert all(l.period == "F5" for l in cell)


# --------------------------------------------------------------------- #
# taker regression: combo_descriptor stays FG-grid-only                  #
# --------------------------------------------------------------------- #

def test_combo_descriptor_returns_none_for_f5_legs():
    assert combo_descriptor([_f5_spread("LAD", 2, "yes"),
                             _f5_total(7, "yes")]) is None
    assert combo_descriptor([_fg_spread("LAD", 2, "yes"),
                             _f5_total(7, "yes")]) is None
    # FG pair still resolves (regression)
    assert combo_descriptor([_fg_spread("LAD", 2, "yes"),
                             _fg_total(9, "yes")]) is not None
