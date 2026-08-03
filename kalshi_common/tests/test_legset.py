from kalshi_common import legset

EVT = "KXMLBGAME-25JUN271905NYYBOS"  # away=NYY, home=BOS

def _spread(team, n, side):   # team code, ticker N (2 => ±1.5), yes/no
    return {"market_ticker": f"KXMLBSPREAD-25JUN271905NYYBOS-{team}{n}",
            "event_ticker": EVT, "side": side}

def _total(n, side):
    return {"market_ticker": f"KXMLBTOTAL-25JUN271905NYYBOS-{n}",
            "event_ticker": EVT, "side": side}

def _ml(team, side):
    return {"market_ticker": f"KXMLBGAME-25JUN271905NYYBOS-{team}",
            "event_ticker": EVT, "side": side}

def test_parse_spread_home_minus_1_5():
    leg = legset.parse_leg(_spread("BOS", 2, "yes"))   # home -1.5, YES
    assert leg == legset.CanonicalLeg(EVT, "spread", -1.5, "home")

def test_parse_total_over_8_5():
    leg = legset.parse_leg(_total(9, "yes"))            # Over 8.5
    assert leg == legset.CanonicalLeg(EVT, "total", 8.5, "over")

def test_parse_ml_home():
    leg = legset.parse_leg(_ml("BOS", "yes"))           # home ML
    assert leg == legset.CanonicalLeg(EVT, "ml", None, "home")

def test_parse_spread_no_side_flips_to_away():
    # home -1.5 NO == away covers +1.5
    leg = legset.parse_leg(_spread("BOS", 2, "no"))
    assert leg == legset.CanonicalLeg(EVT, "spread", -1.5, "away")

def test_unparseable_leg_returns_none():
    bad = {"market_ticker": "KXMLBPLAYER-foo", "event_ticker": EVT, "side": "yes"}
    assert legset.parse_leg(bad) is None

def test_parse_legs_all_or_nothing():
    good = [_spread("BOS", 2, "yes"), _total(9, "yes")]
    assert legset.parse_legs(good) is not None
    assert len(legset.parse_legs(good)) == 2
    bad = good + [{"market_ticker": "KXMLBPLAYER-foo", "event_ticker": EVT, "side": "yes"}]
    assert legset.parse_legs(bad) is None

def test_parse_leg_missing_side_returns_none():
    """Fail-safe: malformed leg dict missing side key -> None"""
    bad = {"market_ticker": "KXMLBSPREAD-25JUN271905NYYBOS-BOS2", "event_ticker": EVT}
    assert legset.parse_leg(bad) is None

def test_parse_total_under_8_5():
    leg = legset.parse_leg(_total(9, "no"))            # Under 8.5
    assert leg == legset.CanonicalLeg(EVT, "total", 8.5, "under")

def test_parse_ml_no_flips_to_away():
    leg = legset.parse_leg(_ml("BOS", "no"))           # home NO == away ML
    assert leg == legset.CanonicalLeg(EVT, "ml", None, "away")

def test_parse_spread_away_team_yes():
    # NYY is away team in EVT; -1.5 YES on away = away covers
    leg = legset.parse_leg(_spread("NYY", 2, "yes"))
    assert leg == legset.CanonicalLeg(EVT, "spread", -1.5, "away")

def test_hash_is_order_independent():
    a = legset.parse_legs([_spread("BOS", 2, "yes"), _total(9, "yes")])
    b = legset.parse_legs([_total(9, "yes"), _spread("BOS", 2, "yes")])
    assert legset.leg_set_hash(a) == legset.leg_set_hash(b)

def test_hash_distinguishes_different_sets():
    a = legset.parse_legs([_spread("BOS", 2, "yes"), _total(9, "yes")])
    c = legset.parse_legs([_spread("BOS", 2, "yes"), _total(9, "no")])  # under
    assert legset.leg_set_hash(a) != legset.leg_set_hash(c)

def test_canonical_legs_sorts_with_none_line():
    legs = legset.parse_legs([_total(9, "yes"), _ml("BOS", "yes"),
                              _spread("BOS", 2, "yes")])
    ordered = legset.canonical_legs(legs)
    assert [l.market_type for l in ordered] == sorted(l.market_type for l in legs)

EVT2 = "KXMLBGAME-25JUN271905LADSFG"  # a different game

def _spread2(team, n, side):
    return {"market_ticker": f"KXMLBSPREAD-25JUN271905LADSFG-{team}{n}",
            "event_ticker": EVT2, "side": side}

def test_partition_groups_two_games():
    legs = legset.parse_legs([_spread("BOS", 2, "yes"), _spread2("SF", 2, "yes")])
    parts = legset.partition_by_game(legs)
    assert set(parts) == {EVT, EVT2}
    assert len(parts[EVT]) == 1 and len(parts[EVT2]) == 1

def test_classify_single():
    legs = legset.parse_legs([_ml("BOS", "yes")])
    assert legset.classify_subcombo(legs) == "single"

def test_classify_grid_spread_total():
    legs = legset.parse_legs([_spread("BOS", 2, "yes"), _total(9, "yes")])
    assert legset.classify_subcombo(legs) == "grid_spread_total"

def test_classify_grid_ml_total():
    legs = legset.parse_legs([_ml("BOS", "yes"), _total(9, "yes")])
    assert legset.classify_subcombo(legs) == "grid_ml_total"

def test_classify_three_leg_is_unpriceable():
    """FLIPPED by issue #66 (was: on_demand). An MLB same-game 3-leg can only
    draw on spread/total/ml, so it always carries the spread+ml pair, and a
    spread leg IMPLIES an ml leg — the 2^N partition is not a partition."""
    legs = legset.parse_legs([_spread("BOS", 2, "yes"), _total(9, "yes"),
                              _ml("BOS", "yes")])
    assert legset.classify_subcombo(legs) == "unpriceable"

def test_classify_spread_ml_is_unpriceable():
    """FLIPPED by issue #66 (was: on_demand)."""
    legs = legset.parse_legs([_spread("BOS", 2, "yes"), _ml("BOS", "yes")])
    assert legset.classify_subcombo(legs) == "unpriceable"

def test_classify_two_totals_is_on_demand():
    legs = legset.parse_legs([_total(9, "yes"), _total(11, "no")])
    assert legset.classify_subcombo(legs) == "on_demand"

def test_classify_empty_is_unpriceable():
    assert legset.classify_subcombo([]) == "unpriceable"


# ---------------------------------------------------------------- #
# Phase 2: enumerate_partition (2^N side-combination cells)        #
# ---------------------------------------------------------------- #

def test_flip_leg_spread_home_away():
    leg = legset.CanonicalLeg(EVT, "spread", -1.5, "home")
    assert legset.flip_leg(leg) == legset.CanonicalLeg(EVT, "spread", -1.5, "away")
    assert legset.flip_leg(legset.flip_leg(leg)) == leg

def test_flip_leg_total_over_under():
    leg = legset.CanonicalLeg(EVT, "total", 8.5, "over")
    assert legset.flip_leg(leg) == legset.CanonicalLeg(EVT, "total", 8.5, "under")

def test_flip_leg_ml_line_stays_none():
    leg = legset.CanonicalLeg(EVT, "ml", None, "away")
    flipped = legset.flip_leg(leg)
    assert flipped.side == "home" and flipped.line is None

def test_enumerate_partition_count_and_target_first():
    legs = [legset.CanonicalLeg(EVT, "spread", -1.5, "home"),
            legset.CanonicalLeg(EVT, "total", 8.5, "over"),
            legset.CanonicalLeg(EVT, "ml", None, "home")]
    cells = legset.enumerate_partition(legs)
    assert len(cells) == 8
    assert cells[0] == legs                     # cell 0 == the target

def test_enumerate_partition_bitmask_ordering():
    # bit j of cell index i flips leg j (LSB = leg 0)
    legs = [legset.CanonicalLeg(EVT, "spread", -1.5, "home"),
            legset.CanonicalLeg(EVT, "total", 8.5, "over")]
    cells = legset.enumerate_partition(legs)
    assert cells[1][0].side == "away" and cells[1][1].side == "over"   # i=1: flip leg0
    assert cells[2][0].side == "home" and cells[2][1].side == "under"  # i=2: flip leg1
    assert cells[3][0].side == "away" and cells[3][1].side == "under"  # i=3: flip both

def test_enumerate_partition_deterministic_and_lines_unchanged():
    legs = [legset.CanonicalLeg(EVT, "total", 7.5, "over"),
            legset.CanonicalLeg(EVT, "total", 8.5, "over")]
    a = legset.enumerate_partition(legs)
    b = legset.enumerate_partition(legs)
    assert a == b
    for cell in a:
        assert [l.line for l in cell] == [7.5, 8.5]


def test_duplicate_market_same_game_is_unpriceable():
    """Adversarial review #3: a creator-crafted combo repeating one market
    (e.g. Over 8.5 AND Under 8.5 — joint probability exactly 0, or the same
    leg twice) must never reach a pricing route: Route B would happily
    transfer-price a contradictory pair a multiplying book returns."""
    over = legset.CanonicalLeg(EVT, "total", 8.5, "over")
    under = legset.CanonicalLeg(EVT, "total", 8.5, "under")
    assert legset.classify_subcombo([over, under]) == "unpriceable"
    assert legset.classify_subcombo([over, over]) == "unpriceable"
    ml_h = legset.CanonicalLeg(EVT, "ml", None, "home")
    ml_a = legset.CanonicalLeg(EVT, "ml", None, "away")
    assert legset.classify_subcombo([ml_h, ml_a]) == "unpriceable"
    sp = legset.CanonicalLeg(EVT, "spread", -1.5, "home")
    assert legset.classify_subcombo([sp, over, ml_h, ml_a]) == "unpriceable"
    # nested totals at DIFFERENT lines stay allowed (Fréchet handles them)
    over75 = legset.CanonicalLeg(EVT, "total", 7.5, "over")
    assert legset.classify_subcombo([over, over75]) == "on_demand"


# ------------------------------------------------------------------ #
# Issue #66: spread+ML implication guard                              #
# ------------------------------------------------------------------ #
#
# Working the margin m = home - away at a 1.5 spread, the four side
# combinations of (spread, ml) within ONE game are:
#
#   spread home + ml home   m>1.5 & m>0    spread => ml     implication
#   spread away + ml away   m<1.5 & m<0    ml => spread     implication
#   spread home + ml away   m>1.5 & m<0    empty            contradiction
#   spread away + ml home   m<1.5 & m>0    0<m<1.5          legitimate overlap
#
# Three of four make the 2^N enumeration a non-partition, which is Route A's
# precondition; the fourth is legitimate but no book prices it (verified live
# across all six on 2026-08-03), so refusing all four costs no coverage.

# A SECOND, well-formed event (away=LAD, home=SF). EVT2 above cannot carry a
# moneyline leg: its suffix parses to no valid home code, so _moneyline_side
# returns None and parse_leg drops it.
EVT3 = "KXMLBGAME-25JUN271905LADSF"

def _ml3(team, side):
    return {"market_ticker": f"KXMLBGAME-25JUN271905LADSF-{team}",
            "event_ticker": EVT3, "side": side}

_SPREAD_SIDES = ("home", "away")
_ML_SIDES = ("home", "away")

def _spread_leg_for(side, n=2, use_away_ticker=False):
    """A spread leg whose CANONICAL side is `side`, minted from either team's
    ticker. Both mintings must be refused — the guard may not depend on which
    team's market the RFQ creator happened to pick."""
    if use_away_ticker:                       # NYY = away in EVT
        return _spread("NYY", n, "yes" if side == "away" else "no")
    return _spread("BOS", n, "yes" if side == "home" else "no")

def _ml_leg_for(side, use_away_ticker=False):
    if use_away_ticker:
        return _ml("NYY", "yes" if side == "away" else "no")
    return _ml("BOS", "yes" if side == "home" else "no")


def test_spread_plus_ml_same_game_is_unpriceable_all_four_sides():
    for s_side in _SPREAD_SIDES:
        for m_side in _ML_SIDES:
            legs = legset.parse_legs([_spread_leg_for(s_side),
                                      _ml_leg_for(m_side)])
            assert [l.side for l in legs] == [s_side, m_side]
            assert legset.classify_subcombo(legs) == "unpriceable", \
                f"spread {s_side} + ml {m_side} must be refused"


def test_spread_plus_ml_refused_from_either_teams_ticker():
    """The away-team mirror: same four canonical shapes, minted off NYY's
    spread/moneyline tickers instead of BOS's."""
    for s_side in _SPREAD_SIDES:
        for m_side in _ML_SIDES:
            legs = legset.parse_legs([
                _spread_leg_for(s_side, use_away_ticker=True),
                _ml_leg_for(m_side, use_away_ticker=True)])
            assert [l.side for l in legs] == [s_side, m_side]
            assert legset.classify_subcombo(legs) == "unpriceable"


def test_three_leg_spread_total_ml_is_unpriceable_all_sides():
    """The shape from issue #66: Novig returns the total+ml price under a
    3-leg label; ProphetX/Caesars devig it through Route B, which does not
    know the spread implies the ML."""
    for s_side in _SPREAD_SIDES:
        for m_side in _ML_SIDES:
            for ou in ("yes", "no"):
                legs = legset.parse_legs([_spread_leg_for(s_side),
                                          _total(9, ou),
                                          _ml_leg_for(m_side)])
                assert legset.classify_subcombo(legs) == "unpriceable"


def test_spread_plus_ml_refusal_is_line_independent():
    """Both spread magnitudes (N=2 => 1.5, N=3 => 2.5) and a two-spread set."""
    for n in (2, 3, 4):
        legs = legset.parse_legs([_spread("BOS", n, "yes"), _ml("BOS", "yes")])
        assert legset.classify_subcombo(legs) == "unpriceable"
    nested = legset.parse_legs([_spread("BOS", 2, "yes"), _spread("BOS", 3, "yes"),
                                _ml("BOS", "yes")])
    assert legset.classify_subcombo(nested) == "unpriceable"


def test_spread_plus_ml_refusal_is_order_independent():
    a = legset.parse_legs([_spread("BOS", 2, "yes"), _total(9, "yes"),
                           _ml("BOS", "yes")])
    b = legset.parse_legs([_ml("BOS", "yes"), _spread("BOS", 2, "yes"),
                           _total(9, "yes")])
    c = legset.parse_legs([_total(9, "yes"), _ml("BOS", "yes"),
                           _spread("BOS", 2, "yes")])
    assert all(legset.classify_subcombo(x) == "unpriceable" for x in (a, b, c))


def test_spread_plus_ml_in_DIFFERENT_games_is_not_refused():
    """The guard is per-game. A spread in game A and an ML in game B are
    independent events — refusing them would destroy real cross-game coverage.
    partition_by_game runs first, so each game classifies as a lone single."""
    legs = legset.parse_legs([_spread("BOS", 2, "yes"), _ml3("SF", "yes")])
    by_game = legset.partition_by_game(legs)
    assert set(by_game) == {EVT, EVT3}
    assert all(legset.classify_subcombo(gl) == "single" for gl in by_game.values())


def test_cross_game_spread_ml_pairs_stay_priceable_all_four_sides():
    """All four side combinations, split across two games, stay priceable."""
    for s_side in _SPREAD_SIDES:
        for m_side in _ML_SIDES:
            legs = legset.parse_legs([
                _spread_leg_for(s_side),
                _ml3("SF" if m_side == "home" else "LAD", "yes")])
            routes = [legset.classify_subcombo(gl)
                      for gl in legset.partition_by_game(legs).values()]
            assert routes == ["single", "single"]


# ---- classify_subcombo_detail: the named reason ---- #

def test_detail_reason_spread_ml_implication():
    legs = legset.parse_legs([_spread("BOS", 2, "yes"), _ml("BOS", "yes")])
    assert legset.classify_subcombo_detail(legs) == \
        ("unpriceable", "spread_ml_implication")


def test_detail_reason_duplicate_market():
    over = legset.CanonicalLeg(EVT, "total", 8.5, "over")
    under = legset.CanonicalLeg(EVT, "total", 8.5, "under")
    assert legset.classify_subcombo_detail([over, under]) == \
        ("unpriceable", "duplicate_market")


def test_detail_reason_empty_leg_set():
    assert legset.classify_subcombo_detail([]) == ("unpriceable", "empty_leg_set")


def test_detail_reason_ok_for_priceable_routes():
    single = legset.parse_legs([_ml("BOS", "yes")])
    grid_st = legset.parse_legs([_spread("BOS", 2, "yes"), _total(9, "yes")])
    grid_mt = legset.parse_legs([_ml("BOS", "yes"), _total(9, "yes")])
    two_tot = legset.parse_legs([_total(9, "yes"), _total(11, "no")])
    assert legset.classify_subcombo_detail(single) == ("single", "ok")
    assert legset.classify_subcombo_detail(grid_st) == ("grid_spread_total", "ok")
    assert legset.classify_subcombo_detail(grid_mt) == ("grid_ml_total", "ok")
    assert legset.classify_subcombo_detail(two_tot) == ("on_demand", "ok")


def test_detail_and_classify_can_never_drift():
    """classify_subcombo MUST be the first element of the detail tuple for
    every shape — two functions answering 'what route' differently is the
    drift this pair is designed to make impossible."""
    shapes = [
        [], legset.parse_legs([_ml("BOS", "yes")]),
        legset.parse_legs([_spread("BOS", 2, "yes"), _total(9, "yes")]),
        legset.parse_legs([_ml("BOS", "yes"), _total(9, "yes")]),
        legset.parse_legs([_total(9, "yes"), _total(11, "no")]),
        legset.parse_legs([_spread("BOS", 2, "yes"), _ml("BOS", "yes")]),
        legset.parse_legs([_spread("BOS", 2, "yes"), _total(9, "yes"),
                           _ml("BOS", "yes")]),
        [legset.CanonicalLeg(EVT, "total", 8.5, "over"),
         legset.CanonicalLeg(EVT, "total", 8.5, "under")],
    ]
    for shape in shapes:
        assert legset.classify_subcombo(shape) == \
            legset.classify_subcombo_detail(shape)[0]


# ---- non-regression: the routes that must NOT move ---- #

def test_issue_66_does_not_move_existing_routes():
    assert legset.classify_subcombo(
        legset.parse_legs([_ml("BOS", "yes")])) == "single"
    assert legset.classify_subcombo(
        legset.parse_legs([_spread("BOS", 2, "yes")])) == "single"
    assert legset.classify_subcombo(
        legset.parse_legs([_spread("BOS", 2, "yes"), _total(9, "yes")])) \
        == "grid_spread_total"
    assert legset.classify_subcombo(
        legset.parse_legs([_spread("NYY", 2, "yes"), _total(9, "no")])) \
        == "grid_spread_total"
    assert legset.classify_subcombo(
        legset.parse_legs([_ml("BOS", "yes"), _total(9, "yes")])) \
        == "grid_ml_total"
    assert legset.classify_subcombo(
        legset.parse_legs([_ml("NYY", "yes"), _total(9, "no")])) \
        == "grid_ml_total"
    # totals-only stays on_demand — out of scope for #66 (see PR body)
    assert legset.classify_subcombo(
        legset.parse_legs([_total(9, "yes"), _total(11, "no")])) == "on_demand"
    # two spreads at different lines, no ML: still on_demand
    assert legset.classify_subcombo(
        legset.parse_legs([_spread("BOS", 2, "yes"), _spread("BOS", 3, "yes")])) \
        == "on_demand"
