"""BetMGM Phase 2 on-demand pricing: resolve_legs + price_selection_set.

resolve_legs maps CanonicalLegs (signed home-perspective spread lines,
home/away/over/under sides) onto the parse_markets(...)["FG"] structure,
returning ResolvedLegs whose ref = (market_id, option_id) and whose
single/opposite decimals come straight from the stored leg tuples.
price_selection_set wraps client.price_picks(fixture_id, refs) into a
bare decimal-or-None. Both are append-only and must never raise.
"""
import sys
from pathlib import Path
from unittest.mock import MagicMock

sys.path.insert(0, str(Path(__file__).parent.parent))

from kalshi_common.legset import CanonicalLeg  # noqa: E402
from mlb_sgp._shared import ResolvedLeg  # noqa: E402


# ------------------------------------------------------------------ #
# Fixtures: raw optionMarkets -> parse_markets structure             #
# ------------------------------------------------------------------ #

def _om(name, opts, mid, bb=True):
    return {"name": {"value": name}, "id": mid,
            "isBetBuilder": bb, "options": opts}


def _opt(name, oid, odds, attr=None):
    return {"name": {"value": name}, "id": oid, "attr": attr,
            "price": {"odds": odds}}


HOME = "San Francisco Giants"
AWAY = "Athletics"


def _raw_markets():
    return [
        _om("Money Line",
            [_opt("Giants", 11, 2.5), _opt("Athletics", 12, 1.53)], mid=100),
        # Home favored: home -1.5 (signed home-perspective key = -1.5)
        _om("Run Line Spread",
            [_opt("San Francisco Giants -1.5", 13, 2.4, "-1,5"),
             _opt("Athletics +1.5", 14, 1.6, "+1,5")], mid=200),
        # Away favored alt line: home +1.5 (signed key = +1.5)
        _om("Run Line Spread",
            [_opt("San Francisco Giants +1.5", 17, 1.55, "+1,5"),
             _opt("Athletics -1.5", 18, 2.45, "-1,5")], mid=201),
        _om("Totals",
            [_opt("Over 8.5", 15, 1.9), _opt("Under 8.5", 16, 1.9)], mid=300),
    ]


def _structure():
    from betmgm import parse_markets
    return parse_markets(_raw_markets(), HOME, AWAY)["FG"]


def _leg(market_type, line, side):
    return CanonicalLeg(game_id="e1", market_type=market_type,
                        line=line, side=side)


# ------------------------------------------------------------------ #
# resolve_legs — happy paths                                         #
# ------------------------------------------------------------------ #

def test_spread_home_signed_lookup():
    from betmgm import resolve_legs
    out = resolve_legs(_structure(), [_leg("spread", -1.5, "home")], HOME, AWAY)
    assert out is not None and len(out) == 1
    rl = out[0]
    assert isinstance(rl, ResolvedLeg)
    assert rl.ref == (200, 13)
    assert rl.opposite_ref == (200, 14)
    assert rl.single_decimal == 2.4
    assert rl.opposite_decimal == 1.6


def test_spread_away_side_same_bucket():
    from betmgm import resolve_legs
    # Away +1.5 lives in the home-signed -1.5 bucket.
    out = resolve_legs(_structure(), [_leg("spread", -1.5, "away")], HOME, AWAY)
    assert out is not None
    rl = out[0]
    assert rl.ref == (200, 14)
    assert rl.opposite_ref == (200, 13)
    assert rl.single_decimal == 1.6
    assert rl.opposite_decimal == 2.4


def test_spread_positive_home_line():
    from betmgm import resolve_legs
    # Away-favorite grid: home +1.5 keyed at +1.5, away side = away -1.5.
    out = resolve_legs(_structure(), [_leg("spread", 1.5, "away")], HOME, AWAY)
    assert out is not None
    assert out[0].ref == (201, 18)
    assert out[0].opposite_ref == (201, 17)
    assert out[0].single_decimal == 2.45


def test_total_over_under():
    from betmgm import resolve_legs
    out = resolve_legs(
        _structure(),
        [_leg("total", 8.5, "over"), _leg("total", 8.5, "under")],
        HOME, AWAY)
    assert out is not None and len(out) == 2
    over, under = out
    assert over.ref == (300, 15) and over.opposite_ref == (300, 16)
    assert under.ref == (300, 16) and under.opposite_ref == (300, 15)
    assert over.single_decimal == 1.9 and over.opposite_decimal == 1.9


def test_moneyline_both_sides():
    from betmgm import resolve_legs
    out = resolve_legs(
        _structure(),
        [_leg("ml", None, "home"), _leg("ml", None, "away")],
        HOME, AWAY)
    assert out is not None and len(out) == 2
    home, away = out
    assert home.ref == (100, 11) and home.opposite_ref == (100, 12)
    assert home.single_decimal == 2.5 and home.opposite_decimal == 1.53
    assert away.ref == (100, 12) and away.opposite_ref == (100, 11)


def test_multi_leg_mixed_set():
    from betmgm import resolve_legs
    legs = [_leg("spread", -1.5, "home"), _leg("total", 8.5, "over"),
            _leg("ml", None, "away")]
    out = resolve_legs(_structure(), legs, HOME, AWAY)
    assert out is not None and len(out) == 3
    assert [rl.ref for rl in out] == [(200, 13), (300, 15), (100, 12)]


# ------------------------------------------------------------------ #
# resolve_legs — misses and fail-safety                              #
# ------------------------------------------------------------------ #

def test_chosen_line_miss_returns_none():
    from betmgm import resolve_legs
    assert resolve_legs(_structure(), [_leg("spread", -2.5, "home")],
                        HOME, AWAY) is None
    assert resolve_legs(_structure(), [_leg("total", 9.5, "over")],
                        HOME, AWAY) is None


def test_missing_moneyline_returns_none():
    from betmgm import resolve_legs
    st = _structure()
    st = {**st, "moneyline": None}
    assert resolve_legs(st, [_leg("ml", None, "home")], HOME, AWAY) is None


def test_one_leg_miss_fails_whole_set():
    from betmgm import resolve_legs
    legs = [_leg("spread", -1.5, "home"), _leg("total", 9.5, "over")]
    assert resolve_legs(_structure(), legs, HOME, AWAY) is None


def test_opposite_side_missing_yields_opposite_ref_none():
    from betmgm import resolve_legs
    # parse_markets only emits two-sided buckets, so hand-build a one-sided
    # structure (Route B shape: chosen side present, opposite absent).
    st = {
        "spreads": {-1.5: {"home": (200, 13, 2.4)}},
        "totals": {8.5: {"over": (300, 15, 1.9)}},
        "moneyline": {"home": (100, 11, 2.5)},
    }
    out = resolve_legs(
        st,
        [_leg("spread", -1.5, "home"), _leg("total", 8.5, "over"),
         _leg("ml", None, "home")],
        HOME, AWAY)
    assert out is not None and len(out) == 3
    for rl in out:
        assert rl.opposite_ref is None
        assert rl.opposite_decimal is None
        assert rl.single_decimal is not None


def test_unknown_market_type_returns_none():
    from betmgm import resolve_legs
    assert resolve_legs(_structure(), [_leg("first_inning", None, "home")],
                        HOME, AWAY) is None


def test_garbage_structure_never_raises():
    from betmgm import resolve_legs
    leg = _leg("spread", -1.5, "home")
    assert resolve_legs(None, [leg], HOME, AWAY) is None
    assert resolve_legs({}, [leg], HOME, AWAY) is None
    assert resolve_legs({"spreads": "oops"}, [leg], HOME, AWAY) is None
    assert resolve_legs({"spreads": {-1.5: {"home": "notatuple"}}},
                        [leg], HOME, AWAY) is None


def test_empty_legs_resolves_to_empty_list():
    from betmgm import resolve_legs
    assert resolve_legs(_structure(), [], HOME, AWAY) == []


# ------------------------------------------------------------------ #
# price_selection_set                                                #
# ------------------------------------------------------------------ #

def test_price_selection_set_happy_path():
    from betmgm import price_selection_set
    client = MagicMock()
    client.price_picks.return_value = {"decimal": 3.5, "american": 250}
    refs = [(200, 13), (300, 15)]
    dec = price_selection_set(client, refs, "fx1")
    assert dec == 3.5
    client.price_picks.assert_called_once_with("fx1", [(200, 13), (300, 15)])


def test_price_selection_set_single_ref():
    from betmgm import price_selection_set
    client = MagicMock()
    client.price_picks.return_value = {"decimal": 1.9}
    assert price_selection_set(client, [(300, 15)], "fx1") == 1.9


def test_price_selection_set_declined_returns_none():
    from betmgm import price_selection_set
    client = MagicMock()
    client.price_picks.return_value = None
    assert price_selection_set(client, [(200, 13)], "fx1") is None


def test_price_selection_set_bad_decimal_returns_none():
    from betmgm import price_selection_set
    client = MagicMock()
    client.price_picks.return_value = {"decimal": 1.0}
    assert price_selection_set(client, [(200, 13)], "fx1") is None
    client.price_picks.return_value = {"decimal": None}
    assert price_selection_set(client, [(200, 13)], "fx1") is None


def test_price_selection_set_empty_refs_returns_none():
    from betmgm import price_selection_set
    client = MagicMock()
    assert price_selection_set(client, [], "fx1") is None
    client.price_picks.assert_not_called()


def test_price_selection_set_never_raises():
    from betmgm import price_selection_set
    client = MagicMock()
    client.price_picks.side_effect = RuntimeError("boom")
    assert price_selection_set(client, [(200, 13)], "fx1") is None
