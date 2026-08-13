"""Issue #87 regression tests: KXMLBRFI (run in 1st inning) leg parsing.

Ticker grammar + settlement verified against the live Kalshi API
(2026-08-13):

- ``KXMLBRFI-26AUG161610MILLAD`` — ONE binary market per event; the market
  ticker IS the event ticker (no trailing team/line suffix, unlike every
  other MLB family).
- rules_primary: "If either team scores a run in the first inning ...
  resolves to Yes", with floor_strike 1 / strike_type greater_or_equal.
  YES is therefore EXACTLY "1st-inning total runs >= 1" = Over 0.5 — an
  identity, not an approximation, so books' 1st-inning 0.5 totals price
  both sides with no push and no conversion (the #86 re-encoding pattern).

RFI legs re-encode at parse as total legs: line 0.5, period "I1",
side yes -> over / no -> under.
"""
import pytest

from kalshi_common import legset
from kalshi_common.leg_types import combo_descriptor

GAME = "26AUG132210MILLAD"            # away=MIL, home=LAD (live 2026-08-13 slate)
EVT = GAME


def _rfi(side):
    # Market ticker == event ticker (live grammar: no suffix).
    return {"market_ticker": f"KXMLBRFI-{GAME}",
            "event_ticker": f"KXMLBRFI-{GAME}", "side": side}


def _fg_total(n, side):
    return {"market_ticker": f"KXMLBTOTAL-{GAME}-{n}",
            "event_ticker": f"KXMLBTOTAL-{GAME}", "side": side}


def _fg_ml(team, side):
    return {"market_ticker": f"KXMLBGAME-{GAME}-{team}",
            "event_ticker": f"KXMLBGAME-{GAME}", "side": side}


def _f5_total(n, side):
    return {"market_ticker": f"KXMLBF5TOTAL-{GAME}-{n}",
            "event_ticker": f"KXMLBF5TOTAL-{GAME}", "side": side}


# --------------------------------------------------------------------- #
# parse_leg: both sides re-encode as I1 totals at 0.5                    #
# --------------------------------------------------------------------- #

def test_parse_rfi_yes_is_first_inning_over_0_5():
    leg = legset.parse_leg(_rfi("yes"))
    assert leg == legset.CanonicalLeg(EVT, "total", 0.5, "over", "I1")


def test_parse_rfi_no_is_first_inning_under_0_5():
    leg = legset.parse_leg(_rfi("no"))
    assert leg == legset.CanonicalLeg(EVT, "total", 0.5, "under", "I1")


def test_parse_rfi_missing_or_junk_side_returns_none():
    """Fail-safe like every other family: an RFQ leg with no usable side
    must not default to a side (a silently wrong side is a mispriced quote)."""
    bad = {"market_ticker": f"KXMLBRFI-{GAME}",
           "event_ticker": f"KXMLBRFI-{GAME}"}
    assert legset.parse_leg(bad) is None
    assert legset.parse_leg(_rfi("maybe")) is None


def test_parse_rfi_game_id_partitions_with_other_families():
    """game_id is the family-independent code (#71): an RFI leg and an FG
    total on the same game must share one partition."""
    legs = legset.parse_legs([_rfi("yes"), _fg_total(9, "yes")])
    assert legs is not None
    by_game = legset.partition_by_game(legs)
    assert set(by_game) == {GAME}
    assert len(by_game[GAME]) == 2


# --------------------------------------------------------------------- #
# hash distinctness + duplicate-market guard                             #
# --------------------------------------------------------------------- #

def test_hash_distinguishes_rfi_from_fg_and_f5_totals():
    """Period is in the hash payload: I1 total 0.5 must hash apart from a
    hypothetical FG/F5 total at the same line."""
    rfi = legset.parse_legs([_rfi("yes")])
    fg = [legset.CanonicalLeg(EVT, "total", 0.5, "over", "FG")]
    f5 = [legset.CanonicalLeg(EVT, "total", 0.5, "over", "F5")]
    assert legset.leg_set_hash(rfi) != legset.leg_set_hash(fg)
    assert legset.leg_set_hash(rfi) != legset.leg_set_hash(f5)


def test_dup_guard_rejects_rfi_pair_same_game():
    """YES+NO same game is a contradictory pair (joint prob 0); YES+YES is
    the same leg twice. Both collide on ("total","I1",0.5) -> unpriceable."""
    over = legset.parse_leg(_rfi("yes"))
    under = legset.parse_leg(_rfi("no"))
    assert legset.classify_subcombo([over, under]) == "unpriceable"
    assert legset.classify_subcombo([over, over]) == "unpriceable"


def test_rfi_bounds_do_not_interact_with_fg_or_f5_totals():
    """I1 and FG/F5 run counts are different variables: NRFI + FG Over 8.5
    is a legitimate correlated combo, not a contradiction."""
    legs = [legset.parse_leg(_rfi("no")), legset.parse_leg(_fg_total(9, "yes"))]
    assert legset.classify_subcombo(legs) == "on_demand"
    legs = [legset.parse_leg(_rfi("no")), legset.parse_leg(_f5_total(5, "yes"))]
    assert legset.classify_subcombo(legs) == "on_demand"


# --------------------------------------------------------------------- #
# classify_subcombo routing                                              #
# --------------------------------------------------------------------- #

def test_rfi_plus_ml_same_game_routes_on_demand():
    """Grids are FG-only surfaces: an RFI-containing 2-leg set must never
    classify grid_ml_total / grid_spread_total."""
    legs = legset.parse_legs([_rfi("yes"), _fg_ml("LAD", "yes")])
    assert legset.classify_subcombo(legs) == "on_demand"


def test_rfi_plus_fg_total_same_game_routes_on_demand():
    legs = legset.parse_legs([_rfi("yes"), _fg_total(9, "yes")])
    assert legset.classify_subcombo(legs) == "on_demand"


def test_lone_rfi_partition_leg_routes_single():
    """A cross-game combo puts the RFI leg alone in its game's partition —
    live routing prices it as a 1-leg on-demand set."""
    leg = legset.parse_leg(_rfi("yes"))
    assert legset.classify_subcombo([leg]) == "single"


# --------------------------------------------------------------------- #
# flip_leg / enumerate_partition                                         #
# --------------------------------------------------------------------- #

def test_flip_leg_round_trips_rfi():
    leg = legset.parse_leg(_rfi("yes"))
    flipped = legset.flip_leg(leg)
    assert flipped.side == "under" and flipped.period == "I1"
    assert flipped.line == 0.5
    assert legset.flip_leg(flipped) == leg


def test_flip_leg_unknown_side_raises_named_error():
    """#87 insurance: flip_leg on a side outside the vocabulary must fail
    LOUDLY with the side and leg named — not a bare KeyError — so the next
    vocab addition can't silently half-land."""
    junk = legset.CanonicalLeg(EVT, "total", 0.5, "maybe", "I1")
    with pytest.raises(ValueError, match="maybe"):
        legset.flip_leg(junk)


def test_enumerate_partition_preserves_i1_period():
    legs = [legset.parse_leg(_rfi("yes")),
            legset.parse_leg(_fg_total(9, "yes"))]
    for cell in legset.enumerate_partition(legs):
        assert [l.period for l in cell] == [l.period for l in legs]
        assert all(l.line in (0.5, 8.5) for l in cell)


# --------------------------------------------------------------------- #
# taker regression: combo_descriptor stays FG-grid-only                  #
# --------------------------------------------------------------------- #

def test_combo_descriptor_returns_none_for_rfi_legs():
    assert combo_descriptor([_rfi("yes"), _fg_total(9, "yes")]) is None
