"""Caesars Phase 2 on-demand resolution + single-combo pricing.

`resolve_legs(structure, legs, home_team, away_team)` maps CanonicalLegs to
the ready-made /bets/details leg dicts that `parse_markets` already builds
(the ref IS the leg dict). The critical Caesars property under test: the AWAY
spread leg dict must carry the NEGATED home line (`caesars.py::_leg` /
`parse_markets` store `line=-ml` for the away selection) — a wrong sign makes
Caesars silently decline the combo.

`price_selection_set(client, refs)` is one `client.price_combo(refs)` call
returning a decimal > 1.0 or None. Both functions are fail-safe (never raise).

No network: structure fixtures are built through the real `parse_markets`
on a hand-built event tree (pattern from test_caesars_moneyline.py).
"""
import sys
from pathlib import Path
from unittest.mock import MagicMock

sys.path.insert(0, str(Path(__file__).parent.parent))

from kalshi_common.legset import CanonicalLeg  # noqa: E402
from mlb_sgp._shared import ResolvedLeg        # noqa: E402


# ------------------------------------------------------------------ #
# Fixtures (mirrors test_caesars_moneyline.py)                       #
# ------------------------------------------------------------------ #

def _sel(t, sid, d):
    return {"id": sid, "type": t, "name": t, "price": {"d": d}}


def _mkt(name, line, sels, mid):
    return {"id": mid, "name": name, "line": line, "selections": sels}


def _event():
    return {"id": "e1", "competitionId": "c1", "keyMarketGroups": [{"markets": [
        _mkt("Moneyline", None,
             [_sel("home", "mh", 1.8), _sel("away", "ma", 2.1)], "m_ml"),
        _mkt("Run Line", -1.5,
             [_sel("home", "sh", 2.4), _sel("away", "sa", 1.6)], "m_rl"),
        _mkt("Total Runs", 8.5,
             [_sel("over", "ov", 1.9), _sel("under", "un", 1.95)], "m_tot"),
    ]}]}


def _structure():
    """The per-period dict resolve_legs consumes: parse_markets(event)["FG"]."""
    from caesars import parse_markets
    return parse_markets(_event())["FG"]


def _leg(market_type, line, side):
    return CanonicalLeg("EV1", market_type, line, side)


HOME, AWAY = "Detroit Tigers", "Houston Astros"


# ------------------------------------------------------------------ #
# resolve_legs — spreads (incl. the away-line negation property)     #
# ------------------------------------------------------------------ #

def test_resolve_spread_home():
    from caesars import resolve_legs
    out = resolve_legs(_structure(), [_leg("spread", -1.5, "home")], HOME, AWAY)
    assert out is not None and len(out) == 1
    rl = out[0]
    assert isinstance(rl, ResolvedLeg)
    assert rl.ref["selectionId"] == "sh"
    assert rl.ref["line"] == -1.5           # home leg carries the market line
    assert rl.opposite_ref["selectionId"] == "sa"
    assert rl.opposite_ref["line"] == 1.5   # away leg carries the NEGATED line
    assert rl.single_decimal == 2.4
    assert rl.opposite_decimal == 1.6


def test_resolve_spread_away_negation():
    """Away chosen side: the returned ref must carry the negated home line —
    Caesars silently declines the combo if the away selection's `line` field
    keeps the home sign."""
    from caesars import resolve_legs
    out = resolve_legs(_structure(), [_leg("spread", -1.5, "away")], HOME, AWAY)
    assert out is not None
    rl = out[0]
    assert rl.ref["selectionId"] == "sa"
    assert rl.ref["line"] == 1.5            # negated: home -1.5 => away +1.5
    assert rl.opposite_ref["selectionId"] == "sh"
    assert rl.opposite_ref["line"] == -1.5
    assert rl.single_decimal == 1.6
    assert rl.opposite_decimal == 2.4


# ------------------------------------------------------------------ #
# resolve_legs — totals / moneyline                                  #
# ------------------------------------------------------------------ #

def test_resolve_total_over_under():
    from caesars import resolve_legs
    out = resolve_legs(_structure(), [_leg("total", 8.5, "over")], HOME, AWAY)
    rl = out[0]
    assert rl.ref["selectionId"] == "ov"
    assert rl.opposite_ref["selectionId"] == "un"
    assert rl.single_decimal == 1.9
    assert rl.opposite_decimal == 1.95

    out = resolve_legs(_structure(), [_leg("total", 8.5, "under")], HOME, AWAY)
    rl = out[0]
    assert rl.ref["selectionId"] == "un"
    assert rl.opposite_ref["selectionId"] == "ov"
    assert rl.single_decimal == 1.95
    assert rl.opposite_decimal == 1.9


def test_resolve_moneyline():
    from caesars import resolve_legs
    out = resolve_legs(_structure(), [_leg("ml", None, "home")], HOME, AWAY)
    rl = out[0]
    assert rl.ref["selectionId"] == "mh"
    assert rl.opposite_ref["selectionId"] == "ma"
    assert rl.single_decimal == 1.8
    assert rl.opposite_decimal == 2.1

    out = resolve_legs(_structure(), [_leg("ml", None, "away")], HOME, AWAY)
    rl = out[0]
    assert rl.ref["selectionId"] == "ma"
    assert rl.opposite_ref["selectionId"] == "mh"


def test_resolve_multi_leg_order_preserved():
    from caesars import resolve_legs
    legs = [_leg("spread", -1.5, "away"), _leg("total", 8.5, "over"),
            _leg("ml", None, "home")]
    out = resolve_legs(_structure(), legs, HOME, AWAY)
    assert [rl.ref["selectionId"] for rl in out] == ["sa", "ov", "mh"]


# ------------------------------------------------------------------ #
# resolve_legs — one-sided / misses / fail-safe                      #
# ------------------------------------------------------------------ #

def test_resolve_one_sided_total_gives_none_opposite():
    from caesars import resolve_legs
    st = _structure()
    del st["totals"][8.5]["under"]
    out = resolve_legs(st, [_leg("total", 8.5, "over")], HOME, AWAY)
    rl = out[0]
    assert rl.ref["selectionId"] == "ov"
    assert rl.opposite_ref is None
    assert rl.single_decimal == 1.9
    assert rl.opposite_decimal is None


def test_resolve_chosen_miss_returns_none():
    from caesars import resolve_legs
    st = _structure()
    # Line not offered at the book.
    assert resolve_legs(st, [_leg("spread", -2.5, "home")], HOME, AWAY) is None
    assert resolve_legs(st, [_leg("total", 9.5, "over")], HOME, AWAY) is None
    # Chosen SIDE missing from an offered line.
    del st["spreads"][-1.5]["home"]
    assert resolve_legs(st, [_leg("spread", -1.5, "home")], HOME, AWAY) is None
    # Moneyline absent entirely.
    st["moneyline"] = None
    assert resolve_legs(st, [_leg("ml", None, "home")], HOME, AWAY) is None
    # One bad leg poisons the whole set (all-or-nothing).
    assert resolve_legs(_structure(),
                        [_leg("total", 8.5, "over"), _leg("spread", -2.5, "home")],
                        HOME, AWAY) is None


def test_resolve_unknown_market_type_returns_none():
    from caesars import resolve_legs
    assert resolve_legs(_structure(), [_leg("prop", 1.5, "over")], HOME, AWAY) is None


def test_resolve_never_raises_on_garbage():
    from caesars import resolve_legs
    assert resolve_legs(None, [_leg("spread", -1.5, "home")], HOME, AWAY) is None
    assert resolve_legs({}, [_leg("total", 8.5, "over")], HOME, AWAY) is None
    assert resolve_legs({"spreads": "not-a-dict", "totals": {}, "moneyline": None},
                        [_leg("spread", -1.5, "home")], HOME, AWAY) is None


def test_resolve_missing_price_gives_none_decimal():
    """A leg dict without a usable price still resolves (ref intact); the
    decimal is just None — _price_d semantics."""
    from caesars import resolve_legs
    st = _structure()
    st["totals"][8.5]["over"]["price"] = {}
    out = resolve_legs(st, [_leg("total", 8.5, "over")], HOME, AWAY)
    rl = out[0]
    assert rl.ref["selectionId"] == "ov"
    assert rl.single_decimal is None
    assert rl.opposite_decimal == 1.95


# ------------------------------------------------------------------ #
# price_selection_set                                                #
# ------------------------------------------------------------------ #

def test_price_selection_set_happy():
    from caesars import price_selection_set
    client = MagicMock()
    client.price_combo.return_value = {"decimal": 3.5, "american": 250}
    st = _structure()
    refs = [st["spreads"][-1.5]["away"], st["totals"][8.5]["over"]]
    assert price_selection_set(client, refs) == 3.5
    client.price_combo.assert_called_once_with(refs)


def test_price_selection_set_none_paths():
    from caesars import price_selection_set
    client = MagicMock()
    client.price_combo.return_value = None            # book declined
    assert price_selection_set(client, [{"selectionId": "x"}]) is None

    client.price_combo.return_value = {"decimal": 1.0}  # not > 1.0
    assert price_selection_set(client, [{"selectionId": "x"}]) is None

    client.price_combo.side_effect = RuntimeError("boom")  # never raises
    assert price_selection_set(client, [{"selectionId": "x"}]) is None

    fresh = MagicMock()
    assert price_selection_set(fresh, []) is None      # empty refs
    fresh.price_combo.assert_not_called()
