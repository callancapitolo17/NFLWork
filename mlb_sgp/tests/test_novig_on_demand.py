"""Phase 2 on-demand pricing: Novig `resolve_legs` + `price_selection_set`.

Structure fixture convention (documented on `resolve_legs`):

    {
      "home_spread": {home_line: leg},   # keyed by HOME team's signed line
      "away_spread": {away_line: leg},   # keyed by AWAY team's line (= -home)
      "over":        {total_line: leg},
      "under":       {total_line: leg},
      "home_ml":     leg | None,
      "away_ml":     leg | None,
    }

where each leg is the `_find_outcome_in_spread` / `_find_outcome_in_total`
shape from scraper_novig_sgp.py: {"id": <outcome uuid>, "available": <implied
prob or None>}.

CanonicalLeg lines are SIGNED home-perspective (negative = home favored), so
an away-side lookup must mirror the sign into the away_spread bucket.

No network anywhere — price_selection_set is exercised with a fake client.
"""
from kalshi_common.legset import CanonicalLeg
from mlb_sgp._shared import ResolvedLeg
from mlb_sgp.novig import price_selection_set, resolve_legs


GAME = "KXMLBGAME-26JUL09NYYBOS"
HOME = "Boston Red Sox"
AWAY = "New York Yankees"


def _leg(uuid, available):
    return {"id": uuid, "available": available}


def _structure():
    """Home favored -1.5 (so away bucket keys +1.5), total 8.5, both MLs."""
    return {
        "home_spread": {-1.5: _leg("hs-15", 0.60)},
        "away_spread": {1.5: _leg("as+15", 0.44)},
        "over": {8.5: _leg("ov85", 0.52)},
        "under": {8.5: _leg("un85", 0.50)},
        "home_ml": _leg("hml", 0.62),
        "away_ml": _leg("aml", 0.42),
    }


# ---------------------------------------------------------------------------
# resolve_legs — spread sign mirroring
# ---------------------------------------------------------------------------

def test_spread_home_side_resolves_with_opposite():
    out = resolve_legs(
        _structure(),
        [CanonicalLeg(GAME, "spread", -1.5, "home")],
        HOME, AWAY,
    )
    assert out is not None and len(out) == 1
    rl = out[0]
    assert isinstance(rl, ResolvedLeg)
    assert rl.ref == "hs-15"
    assert rl.opposite_ref == "as+15"
    assert rl.single_decimal == 1.0 / 0.60
    assert rl.opposite_decimal == 1.0 / 0.44


def test_spread_away_side_mirrors_sign():
    """Canonical line -1.5 + side='away' means away +1.5: lookup must hit
    away_spread[+1.5] and the opposite must be home_spread[-1.5]."""
    out = resolve_legs(
        _structure(),
        [CanonicalLeg(GAME, "spread", -1.5, "away")],
        HOME, AWAY,
    )
    assert out is not None
    rl = out[0]
    assert rl.ref == "as+15"
    assert rl.opposite_ref == "hs-15"
    assert rl.single_decimal == 1.0 / 0.44
    assert rl.opposite_decimal == 1.0 / 0.60


def test_spread_positive_home_line_away_favored():
    """Home +2.5 (away favored): home bucket keyed +2.5, away keyed -2.5."""
    structure = {
        "home_spread": {2.5: _leg("hs+25", 0.55)},
        "away_spread": {-2.5: _leg("as-25", 0.49)},
        "over": {}, "under": {},
        "home_ml": None, "away_ml": None,
    }
    out = resolve_legs(
        structure, [CanonicalLeg(GAME, "spread", 2.5, "home")], HOME, AWAY)
    assert out[0].ref == "hs+25" and out[0].opposite_ref == "as-25"

    out = resolve_legs(
        structure, [CanonicalLeg(GAME, "spread", 2.5, "away")], HOME, AWAY)
    assert out[0].ref == "as-25" and out[0].opposite_ref == "hs+25"


# ---------------------------------------------------------------------------
# resolve_legs — totals and moneyline
# ---------------------------------------------------------------------------

def test_total_over_and_under_same_line():
    out = resolve_legs(
        _structure(),
        [CanonicalLeg(GAME, "total", 8.5, "over"),
         CanonicalLeg(GAME, "total", 8.5, "under")],
        HOME, AWAY,
    )
    assert out is not None and len(out) == 2
    over, under = out
    assert over.ref == "ov85" and over.opposite_ref == "un85"
    assert under.ref == "un85" and under.opposite_ref == "ov85"
    assert over.single_decimal == 1.0 / 0.52
    assert over.opposite_decimal == 1.0 / 0.50


def test_moneyline_home_and_away():
    out = resolve_legs(
        _structure(),
        [CanonicalLeg(GAME, "ml", None, "home"),
         CanonicalLeg(GAME, "ml", None, "away")],
        HOME, AWAY,
    )
    assert out is not None and len(out) == 2
    hml, aml = out
    assert hml.ref == "hml" and hml.opposite_ref == "aml"
    assert aml.ref == "aml" and aml.opposite_ref == "hml"
    assert hml.single_decimal == 1.0 / 0.62


def test_multi_leg_order_preserved():
    out = resolve_legs(
        _structure(),
        [CanonicalLeg(GAME, "total", 8.5, "under"),
         CanonicalLeg(GAME, "spread", -1.5, "home"),
         CanonicalLeg(GAME, "ml", None, "away")],
        HOME, AWAY,
    )
    assert [rl.ref for rl in out] == ["un85", "hs-15", "aml"]


# ---------------------------------------------------------------------------
# resolve_legs — misses
# ---------------------------------------------------------------------------

def test_chosen_side_miss_returns_none():
    """Line not offered at the book -> whole resolution fails (None)."""
    out = resolve_legs(
        _structure(),
        [CanonicalLeg(GAME, "spread", -2.5, "home")],  # only -1.5 offered
        HOME, AWAY,
    )
    assert out is None


def test_chosen_miss_in_any_leg_fails_all():
    out = resolve_legs(
        _structure(),
        [CanonicalLeg(GAME, "total", 8.5, "over"),
         CanonicalLeg(GAME, "total", 9.5, "under")],  # 9.5 not offered
        HOME, AWAY,
    )
    assert out is None


def test_no_integer_line_fallback():
    """Integer total 8.0 with 7.5/8.5 neighbors offered must NOT fall back
    to the adjacent-half-point derivation — exact line match only."""
    structure = {
        "home_spread": {}, "away_spread": {},
        "over": {7.5: _leg("ov75", 0.6), 8.5: _leg("ov85", 0.5)},
        "under": {7.5: _leg("un75", 0.42), 8.5: _leg("un85", 0.52)},
        "home_ml": None, "away_ml": None,
    }
    out = resolve_legs(
        structure, [CanonicalLeg(GAME, "total", 8.0, "over")], HOME, AWAY)
    assert out is None


def test_opposite_miss_yields_none_opposite_ref():
    """One-sided line: chosen resolves, opposite_ref/decimal are None."""
    structure = {
        "home_spread": {-1.5: _leg("hs-15", 0.60)},
        "away_spread": {},                            # away side not offered
        "over": {}, "under": {},
        "home_ml": _leg("hml", 0.62), "away_ml": None,  # ML one-sided too
    }
    out = resolve_legs(
        structure,
        [CanonicalLeg(GAME, "spread", -1.5, "home"),
         CanonicalLeg(GAME, "ml", None, "home")],
        HOME, AWAY,
    )
    assert out is not None and len(out) == 2
    sp, ml = out
    assert sp.ref == "hs-15" and sp.opposite_ref is None
    assert sp.opposite_decimal is None
    assert ml.ref == "hml" and ml.opposite_ref is None


def test_chosen_ml_missing_returns_none():
    structure = dict(_structure(), home_ml=None)
    out = resolve_legs(
        structure, [CanonicalLeg(GAME, "ml", None, "home")], HOME, AWAY)
    assert out is None


# ---------------------------------------------------------------------------
# resolve_legs — available -> decimal edge values
# ---------------------------------------------------------------------------

def test_available_edge_values():
    """available=0, 1, None, negative -> single_decimal None; 0.5 -> 2.0."""
    structure = {
        "home_spread": {-1.5: _leg("a", 0.0)},
        "away_spread": {1.5: _leg("b", 1.0)},
        "over": {8.5: _leg("c", None)},
        "under": {8.5: _leg("d", -0.2)},
        "home_ml": _leg("e", 0.5),
        "away_ml": _leg("f", 0.25),
    }
    out = resolve_legs(
        structure,
        [CanonicalLeg(GAME, "spread", -1.5, "home"),   # avail 0 / opp 1.0
         CanonicalLeg(GAME, "total", 8.5, "over"),     # avail None / opp -0.2
         CanonicalLeg(GAME, "ml", None, "home")],      # avail 0.5 / opp 0.25
        HOME, AWAY,
    )
    assert out is not None
    sp, tot, ml = out
    assert sp.single_decimal is None          # available == 0
    assert sp.opposite_decimal is None        # available == 1
    assert tot.single_decimal is None         # available is None
    assert tot.opposite_decimal is None       # available < 0
    assert ml.single_decimal == 2.0           # 1 / 0.5
    assert ml.opposite_decimal == 4.0         # 1 / 0.25


def test_available_non_numeric_is_none_not_raise():
    structure = {
        "home_spread": {-1.5: _leg("a", "garbage")},
        "away_spread": {1.5: _leg("b", 0.5)},
        "over": {}, "under": {}, "home_ml": None, "away_ml": None,
    }
    out = resolve_legs(
        structure, [CanonicalLeg(GAME, "spread", -1.5, "home")], HOME, AWAY)
    assert out is not None
    assert out[0].single_decimal is None
    assert out[0].opposite_decimal == 2.0


# ---------------------------------------------------------------------------
# resolve_legs — fail-safety
# ---------------------------------------------------------------------------

def test_unknown_market_type_returns_none():
    out = resolve_legs(
        _structure(), [CanonicalLeg(GAME, "prop", 1.5, "home")], HOME, AWAY)
    assert out is None


def test_unknown_side_returns_none():
    out = resolve_legs(
        _structure(), [CanonicalLeg(GAME, "total", 8.5, "home")], HOME, AWAY)
    assert out is None


def test_garbage_structure_never_raises():
    for structure in (None, {}, {"home_spread": "not-a-dict"}, 42):
        out = resolve_legs(
            structure, [CanonicalLeg(GAME, "spread", -1.5, "home")],
            HOME, AWAY)
        assert out is None


def test_empty_legs_returns_empty_list():
    assert resolve_legs(_structure(), [], HOME, AWAY) == []


def test_leg_with_missing_id_is_a_miss():
    """A chosen outcome without a usable UUID can't be priced -> None."""
    structure = dict(_structure())
    structure["home_spread"] = {-1.5: {"id": None, "available": 0.6}}
    out = resolve_legs(
        structure, [CanonicalLeg(GAME, "spread", -1.5, "home")], HOME, AWAY)
    assert out is None


def test_opposite_with_missing_id_treated_one_sided():
    structure = dict(_structure())
    structure["away_spread"] = {1.5: {"id": None, "available": 0.44}}
    out = resolve_legs(
        structure, [CanonicalLeg(GAME, "spread", -1.5, "home")], HOME, AWAY)
    assert out is not None
    assert out[0].ref == "hs-15" and out[0].opposite_ref is None


# ---------------------------------------------------------------------------
# price_selection_set
# ---------------------------------------------------------------------------

class _FakeClient:
    def __init__(self, response):
        self.response = response
        self.calls = []

    def submit_parlay(self, outcome_ids, stake=1.0):
        self.calls.append(list(outcome_ids))
        if isinstance(self.response, Exception):
            raise self.response
        return self.response


def test_price_selection_set_happy_path():
    client = _FakeClient({"decimal": 2.85, "american": 185, "status": "OPEN"})
    out = price_selection_set(client, ["u1", "u2", "u3"])
    assert out == 2.85
    assert client.calls == [["u1", "u2", "u3"]]


def test_price_selection_set_single_ref():
    client = _FakeClient({"decimal": 1.91})
    assert price_selection_set(client, ["u1"]) == 1.91


def test_price_selection_set_empty_dict_returns_none():
    client = _FakeClient({})
    assert price_selection_set(client, ["u1", "u2"]) is None


def test_price_selection_set_decimal_not_above_one_returns_none():
    assert price_selection_set(_FakeClient({"decimal": 1.0}), ["u1"]) is None
    assert price_selection_set(_FakeClient({"decimal": 0.9}), ["u1"]) is None


def test_price_selection_set_never_raises():
    assert price_selection_set(
        _FakeClient(RuntimeError("boom")), ["u1"]) is None
    assert price_selection_set(
        _FakeClient({"decimal": "garbage"}), ["u1"]) is None
    assert price_selection_set(_FakeClient(None), ["u1"]) is None


def test_price_selection_set_empty_refs_returns_none():
    client = _FakeClient({"decimal": 2.0})
    assert price_selection_set(client, []) is None
    assert client.calls == []
