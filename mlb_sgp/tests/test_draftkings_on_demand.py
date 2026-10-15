"""DraftKings on-demand pricing (Phase 2): resolve_legs + price_selection_set.

Mirrors the inline-fixture style of test_dk_moneyline.py — hand-built
structure dicts, no network, fake duck-typed session for the wire call.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from kalshi_common.legset import CanonicalLeg  # noqa: E402


HOME, AWAY = "New York Yankees", "Boston Red Sox"


def _structure():
    """Two-period structure (#85) with the FG bucket populated. Home
    favored -1.5.

    DK keys spreads as (sign, abs_line, participant): sign "N" = the
    negative-line side, participant "1" = home, "3" = away. With home
    favored at -1.5, home is ("N", 1.5, "1") and away ("P", 1.5, "3").
    F5 coverage lives in test_f5_on_demand_resolvers.py.
    """
    return {
        "FG": {
            "spreads": {
                ("N", 1.5, "1"): ["0HC100N150_1", "0HC101N150_1"],
                ("P", 1.5, "3"): ["0HC100P150_3"],
            },
            "totals": {
                ("O", 8.5): ["0OU200O850_1"],
                ("U", 8.5): ["0OU200U850_1"],
            },
            "moneyline": {"1": ["0ML300_1"], "3": ["0ML300_3"]},
            "canonical": {"100", "200", "300"},
        },
        "F5": None,
    }


def _structure_home_dog():
    """Home underdog +1.5: home is the P side, away the N side."""
    return {
        "FG": {
            "spreads": {
                ("P", 1.5, "1"): ["0HC100P150_1"],
                ("N", 1.5, "3"): ["0HC100N150_3"],
            },
            "totals": {},
            "moneyline": {"1": [], "3": []},
            "canonical": set(),
        },
        "F5": None,
    }


# ---------------------------------------------------------------------------
# resolve_legs — spreads
# ---------------------------------------------------------------------------
def test_resolve_home_spread_hit_with_opposite():
    from mlb_sgp.draftkings import resolve_legs
    legs = [CanonicalLeg("g1", "spread", -1.5, "home")]
    out = resolve_legs(_structure(), legs, HOME, AWAY)
    assert out is not None and len(out) == 1
    rl = out[0]
    assert rl.ref == "0HC100N150_1"           # FIRST sel id of chosen side
    assert rl.opposite_ref == "0HC100P150_3"  # away at mirrored line
    assert rl.single_decimal is None          # no straight_odds in structure
    assert rl.opposite_decimal is None


def test_resolve_carries_straight_odds_when_the_structure_has_them():
    """1-leg fairs (kalshi_rfi) and Route B devig DK's straight prices: the
    SGP builder's 1-leg price is not the straight price."""
    from mlb_sgp.draftkings import resolve_legs
    structure = _structure()
    structure["straight_odds"] = {"0HC100N150_1": 1.95, "0HC100P150_3": 1.87}
    rl = resolve_legs(structure, [CanonicalLeg("g1", "spread", -1.5, "home")],
                      HOME, AWAY)[0]
    assert (rl.single_decimal, rl.opposite_decimal) == (1.95, 1.87)


def test_resolve_leaves_a_suspended_side_without_odds():
    """A selection missing from straight_odds (suspended) gets None, never a
    borrowed number."""
    from mlb_sgp.draftkings import resolve_legs
    structure = _structure()
    structure["straight_odds"] = {"0HC100N150_1": 1.95}
    rl = resolve_legs(structure, [CanonicalLeg("g1", "spread", -1.5, "home")],
                      HOME, AWAY)[0]
    assert rl.single_decimal == 1.95
    assert rl.opposite_decimal is None


def test_resolve_away_spread_mirrored_sign():
    from mlb_sgp.draftkings import resolve_legs
    # Same signed home line -1.5, but the AWAY side: sign flips to P,
    # participant flips to "3"; opposite is the home N/"1" entry.
    legs = [CanonicalLeg("g1", "spread", -1.5, "away")]
    out = resolve_legs(_structure(), legs, HOME, AWAY)
    assert out is not None
    assert out[0].ref == "0HC100P150_3"
    assert out[0].opposite_ref == "0HC100N150_1"


def test_resolve_home_dog_positive_line():
    from mlb_sgp.draftkings import resolve_legs
    # Home +1.5 (positive signed line): home is "P"/"1", away "N"/"3".
    out = resolve_legs(
        _structure_home_dog(), [CanonicalLeg("g1", "spread", 1.5, "home")],
        HOME, AWAY)
    assert out[0].ref == "0HC100P150_1"
    assert out[0].opposite_ref == "0HC100N150_3"
    out = resolve_legs(
        _structure_home_dog(), [CanonicalLeg("g1", "spread", 1.5, "away")],
        HOME, AWAY)
    assert out[0].ref == "0HC100N150_3"
    assert out[0].opposite_ref == "0HC100P150_1"


# ---------------------------------------------------------------------------
# resolve_legs — totals
# ---------------------------------------------------------------------------
def test_resolve_total_over_and_under():
    from mlb_sgp.draftkings import resolve_legs
    out = resolve_legs(
        _structure(), [CanonicalLeg("g1", "total", 8.5, "over")], HOME, AWAY)
    assert out[0].ref == "0OU200O850_1"
    assert out[0].opposite_ref == "0OU200U850_1"
    out = resolve_legs(
        _structure(), [CanonicalLeg("g1", "total", 8.5, "under")], HOME, AWAY)
    assert out[0].ref == "0OU200U850_1"
    assert out[0].opposite_ref == "0OU200O850_1"


# ---------------------------------------------------------------------------
# resolve_legs — moneyline
# ---------------------------------------------------------------------------
def test_resolve_ml_both_participants():
    from mlb_sgp.draftkings import resolve_legs
    out = resolve_legs(
        _structure(), [CanonicalLeg("g1", "ml", None, "home")], HOME, AWAY)
    assert out[0].ref == "0ML300_1"           # participant "1" = home
    assert out[0].opposite_ref == "0ML300_3"
    out = resolve_legs(
        _structure(), [CanonicalLeg("g1", "ml", None, "away")], HOME, AWAY)
    assert out[0].ref == "0ML300_3"
    assert out[0].opposite_ref == "0ML300_1"


# ---------------------------------------------------------------------------
# resolve_legs — multi-leg + failure routing
# ---------------------------------------------------------------------------
def test_resolve_multi_leg_order_preserved():
    from mlb_sgp.draftkings import resolve_legs
    legs = [
        CanonicalLeg("g1", "spread", -1.5, "home"),
        CanonicalLeg("g1", "total", 8.5, "over"),
        CanonicalLeg("g1", "ml", None, "away"),
    ]
    out = resolve_legs(_structure(), legs, HOME, AWAY)
    assert [rl.ref for rl in out] == [
        "0HC100N150_1", "0OU200O850_1", "0ML300_3"]


def test_resolve_chosen_side_miss_fails_whole_set():
    from mlb_sgp.draftkings import resolve_legs
    # Line not offered at all → whole resolution fails.
    legs = [CanonicalLeg("g1", "spread", -2.5, "home")]
    assert resolve_legs(_structure(), legs, HOME, AWAY) is None
    # One good leg + one bad leg still fails as a whole.
    legs = [CanonicalLeg("g1", "total", 8.5, "over"),
            CanonicalLeg("g1", "total", 9.5, "over")]
    assert resolve_legs(_structure(), legs, HOME, AWAY) is None
    # Chosen side present-but-empty list also fails.
    s = _structure()
    s["FG"]["totals"][("O", 8.5)] = []
    legs = [CanonicalLeg("g1", "total", 8.5, "over")]
    assert resolve_legs(s, legs, HOME, AWAY) is None


def test_resolve_opposite_side_miss_yields_none_opposite():
    from mlb_sgp.draftkings import resolve_legs
    s = _structure()
    del s["FG"]["spreads"][("P", 1.5, "3")]   # book one-sides the spread
    out = resolve_legs(
        s, [CanonicalLeg("g1", "spread", -1.5, "home")], HOME, AWAY)
    assert out is not None
    assert out[0].ref == "0HC100N150_1"
    assert out[0].opposite_ref is None        # routes book to Route B
    # Same for a one-sided total.
    s = _structure()
    s["FG"]["totals"][("U", 8.5)] = []
    out = resolve_legs(
        s, [CanonicalLeg("g1", "total", 8.5, "over")], HOME, AWAY)
    assert out[0].opposite_ref is None


def test_resolve_malformed_structure_returns_none():
    from mlb_sgp.draftkings import resolve_legs
    legs = [CanonicalLeg("g1", "spread", -1.5, "home")]
    assert resolve_legs(None, legs, HOME, AWAY) is None
    assert resolve_legs({}, legs, HOME, AWAY) is None
    assert resolve_legs({"spreads": ["not", "a", "dict"]}, legs, HOME, AWAY) is None
    assert resolve_legs("garbage", legs, HOME, AWAY) is None
    # Unknown market_type / side never raise either.
    assert resolve_legs(_structure(),
                        [CanonicalLeg("g1", "runline", -1.5, "home")],
                        HOME, AWAY) is None
    assert resolve_legs(_structure(),
                        [CanonicalLeg("g1", "total", 8.5, "home")],
                        HOME, AWAY) is None


# ---------------------------------------------------------------------------
# price_selection_set — fake session wire tests (DK's SGP-builder GET)
# ---------------------------------------------------------------------------
class FakeResp:
    def __init__(self, status_code, payload=None, text=""):
        self.status_code = status_code
        self._payload = payload or {}
        self.text = text
        self.headers = {}

    def json(self):
        return self._payload


class FakeSession:
    def __init__(self, resp):
        self._resp = resp
        self.calls = []

    def get(self, url, params=None, headers=None, **kwargs):
        self.calls.append({"url": url, "params": params, "headers": headers})
        if isinstance(self._resp, BaseException):
            raise self._resp
        return self._resp


def _client(resp):
    from mlb_sgp.dk_client import DraftKingsClient
    client = DraftKingsClient.__new__(DraftKingsClient)
    client.session = FakeSession(resp)
    return client


def _builder_answer(base, candidates):
    """DK's SGP-builder body: base taken whole, candidate prices listed."""
    return {"data": {
        "trueOdds": 2.0, "selectionsMapped": list(base),
        "selectionsNotMapped": [],
        "compatibleMarkets": [{"id": "m", "selections": [
            {"id": sid, "trueOdds": odds} for sid, odds in candidates.items()]}],
    }}


def test_price_selection_set_happy_path():
    from mlb_sgp.dk_client import DK_SGP_PRICE_URL
    from mlb_sgp.draftkings import price_selection_set
    refs = ["0HC100N150_1", "0ML300_1", "0OU200O850_1"]
    client = _client(FakeResp(200, _builder_answer(
        refs[:2], {"0OU200O850_1": 6.42, "0OU200U850_3": 5.1})))
    assert price_selection_set(client, refs, "evt9") == 6.42
    # Base = every leg but the last; candidate = the last leg's market.
    (call,) = client.session.calls
    assert call["url"] == DK_SGP_PRICE_URL
    assert call["headers"] == {"X-SportId": "7"}
    assert call["params"] == {"eventId": "evt9",
                              "selections": "0HC100N150_1,0ML300_1",
                              "marketCandidates": "3_200",
                              "oddsStyle": "american"}


def test_price_selection_set_single_ref_declines_without_a_call():
    """The builder prices nothing on an empty base, and its 1-leg price is
    not DK's straight price — a lone leg prices from straight_odds."""
    from mlb_sgp.draftkings import price_selection_set
    client = _client(FakeResp(200, {}))
    assert price_selection_set(client, ["0ML300_1"], "evt9") is None
    assert client.session.calls == []


def test_price_selection_set_400_returns_none():
    """400 InvalidSelection: a line DK pulled since the structure was read."""
    from mlb_sgp.draftkings import price_selection_set
    client = _client(FakeResp(400, {"code": "InvalidSelection"}))
    assert price_selection_set(client, ["0HC100N150_1", "0OU200O850_1"],
                               "evt9") is None


def test_price_selection_set_dropped_base_leg_returns_none():
    """DK drops a base leg it cannot combine and prices the rest — that
    smaller-base price must never come back as this set's price."""
    from mlb_sgp.draftkings import price_selection_set
    body = _builder_answer(["0HC100N150_1"], {"0OU200O850_1": 3.5})
    body["data"]["selectionsNotMapped"] = [
        {"selectionId": "0ML300_1", "removalReasonCode": "IncompatibleSelections"}]
    client = _client(FakeResp(200, body))
    assert price_selection_set(
        client, ["0HC100N150_1", "0ML300_1", "0OU200O850_1"], "evt9") is None


def test_price_selection_set_candidate_absent_returns_none():
    """A candidate missing from compatibleMarkets is one DK declines."""
    from mlb_sgp.draftkings import price_selection_set
    client = _client(FakeResp(200, _builder_answer(
        ["0HC100N150_1"], {"0OU200U850_3": 3.76})))
    assert price_selection_set(client, ["0HC100N150_1", "0OU200O850_1"],
                               "evt9") is None


def test_price_selection_set_wire_exception_returns_none():
    from mlb_sgp.draftkings import price_selection_set
    client = _client(RuntimeError("network down"))
    assert price_selection_set(client, ["0HC100N150_1", "0OU200O850_1"],
                               "evt9") is None
