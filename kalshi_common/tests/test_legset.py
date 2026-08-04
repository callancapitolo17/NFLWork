from kalshi_common import legset

# Kalshi puts every market family in its OWN event, so the three helpers below
# deliberately emit three DIFFERENT event_tickers for one game. An earlier
# version of this file shared one event_ticker across all three families —
# input Kalshi never emits — which is what hid issue #71 (the same-game grids
# were unreachable in production while these tests passed).
GAME = "25JUN271905NYYBOS"          # away=NYY, home=BOS
EVT = GAME                           # the canonical game_id: family-independent

def _spread(team, n, side):   # team code, ticker N (2 => ±1.5), yes/no
    return {"market_ticker": f"KXMLBSPREAD-{GAME}-{team}{n}",
            "event_ticker": f"KXMLBSPREAD-{GAME}", "side": side}

def _total(n, side):
    return {"market_ticker": f"KXMLBTOTAL-{GAME}-{n}",
            "event_ticker": f"KXMLBTOTAL-{GAME}", "side": side}

def _ml(team, side):
    return {"market_ticker": f"KXMLBGAME-{GAME}-{team}",
            "event_ticker": f"KXMLBGAME-{GAME}", "side": side}

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

GAME2 = "25JUN271905LADSFG"          # a different game
EVT2 = GAME2

def _spread2(team, n, side):
    return {"market_ticker": f"KXMLBSPREAD-{GAME2}-{team}{n}",
            "event_ticker": f"KXMLBSPREAD-{GAME2}", "side": side}

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

def test_classify_three_leg_is_on_demand():
    legs = legset.parse_legs([_spread("BOS", 2, "yes"), _total(9, "yes"),
                              _ml("BOS", "yes")])
    assert legset.classify_subcombo(legs) == "on_demand"

def test_classify_spread_ml_is_on_demand():
    legs = legset.parse_legs([_spread("BOS", 2, "yes"), _ml("BOS", "yes")])
    assert legset.classify_subcombo(legs) == "on_demand"

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


# --------------------------------------------------------------------- #
# Issue #71: game_id must be family-independent                          #
# --------------------------------------------------------------------- #

def test_same_game_legs_share_one_partition_across_market_families():
    """THE #71 regression. Kalshi gives each market family its own event, so
    a spread leg and a total leg on one game arrive with DIFFERENT
    event_tickers. Keying game_id on the event_ticker put them in separate
    partitions, which made grid_spread_total unreachable in production and
    silently priced 86 real same-game RFQs as independent games multiplied
    together (the zero-correlation error)."""
    legs = legset.parse_legs([_spread("BOS", 2, "yes"), _total(9, "yes")])
    assert legs is not None
    assert legs[0].game_id == legs[1].game_id == GAME, (
        "spread and total on one game must canonicalize to one game_id")
    by_game = legset.partition_by_game(legs)
    assert len(by_game) == 1, f"expected 1 partition, got {list(by_game)}"
    assert legset.classify_subcombo(next(iter(by_game.values()))) == \
        "grid_spread_total"


def test_all_three_families_partition_together():
    legs = legset.parse_legs([_spread("BOS", 2, "yes"), _total(9, "yes"),
                              _ml("BOS", "yes")])
    by_game = legset.partition_by_game(legs)
    assert len(by_game) == 1
    assert legset.classify_subcombo(next(iter(by_game.values()))) == "on_demand"


def test_different_games_still_partition_apart():
    """The fix must not over-merge: two genuinely different games stay split
    even when they share a market family."""
    other = {"market_ticker": "KXMLBGAME-25JUN271905LADSF-SF",
             "event_ticker": "KXMLBGAME-25JUN271905LADSF", "side": "yes"}
    legs = legset.parse_legs([_ml("BOS", "yes"), other])
    by_game = legset.partition_by_game(legs)
    assert len(by_game) == 2
    assert all(legset.classify_subcombo(v) == "single"
               for v in by_game.values())


def test_leg_set_hash_still_separates_the_two_families():
    """game_id no longer distinguishes a spread leg from a total leg on the
    same game — market_type must still do it, or the on-demand store would
    collide two different combos onto one hash."""
    spread_only = legset.parse_legs([_spread("BOS", 2, "yes")])
    total_only = legset.parse_legs([_total(9, "yes")])
    assert (legset.leg_set_hash(spread_only)
            != legset.leg_set_hash(total_only))
