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
    """FG bucket of fetch_selection_ids output. Home favored -1.5.

    DK keys spreads as (sign, abs_line, participant): sign "N" = the
    negative-line side, participant "1" = home, "3" = away. With home
    favored at -1.5, home is ("N", 1.5, "1") and away ("P", 1.5, "3").
    """
    return {
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
    }


def _structure_home_dog():
    """Home underdog +1.5: home is the P side, away the N side."""
    return {
        "spreads": {
            ("P", 1.5, "1"): ["0HC100P150_1"],
            ("N", 1.5, "3"): ["0HC100N150_3"],
        },
        "totals": {},
        "moneyline": {"1": [], "3": []},
        "canonical": set(),
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
    assert rl.single_decimal is None          # DK structure carries no odds
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
    s["totals"][("O", 8.5)] = []
    legs = [CanonicalLeg("g1", "total", 8.5, "over")]
    assert resolve_legs(s, legs, HOME, AWAY) is None


def test_resolve_opposite_side_miss_yields_none_opposite():
    from mlb_sgp.draftkings import resolve_legs
    s = _structure()
    del s["spreads"][("P", 1.5, "3")]         # book one-sides the spread
    out = resolve_legs(
        s, [CanonicalLeg("g1", "spread", -1.5, "home")], HOME, AWAY)
    assert out is not None
    assert out[0].ref == "0HC100N150_1"
    assert out[0].opposite_ref is None        # routes book to Route B
    # Same for a one-sided total.
    s = _structure()
    s["totals"][("U", 8.5)] = []
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
# price_selection_set — fake session wire tests
# ---------------------------------------------------------------------------
class FakeResp:
    def __init__(self, status_code, payload=None):
        self.status_code = status_code
        self._payload = payload or {}

    def json(self):
        return self._payload


class FakeSession:
    def __init__(self, resp):
        self._resp = resp
        self.last_url = None
        self.last_json = None

    def post(self, url, json=None, **kwargs):
        self.last_url = url
        self.last_json = json
        return self._resp


class FakeClient:
    def __init__(self, resp):
        self.session = FakeSession(resp)


def test_price_selection_set_happy_path():
    from mlb_sgp.draftkings import price_selection_set
    refs = ["0HC100N150_1", "0OU200O850_1", "0ML300_1"]
    client = FakeClient(FakeResp(200, {
        "combinabilityRestrictions": [],
        "bets": [{"trueOdds": 6.42, "displayOdds": "+542",
                  "selectionsMapped": [{}, {}, {}]}],
    }))
    assert price_selection_set(client, refs) == 6.42
    # Body follows the calculateBets contract: one entry per ref, all in
    # yourBetGroup 0 (one SGP), american odds style.
    body = client.session.last_json
    assert body["selectionsForYourBet"] == [
        {"id": sid, "yourBetGroup": 0} for sid in refs]
    assert body["oddsStyle"] == "american"


def test_price_selection_set_single_ref():
    from mlb_sgp.draftkings import price_selection_set
    client = FakeClient(FakeResp(200, {
        "bets": [{"trueOdds": 1.87, "selectionsMapped": [{}]}],
    }))
    assert price_selection_set(client, ["0ML300_1"]) == 1.87


def test_price_selection_set_422_returns_none():
    from mlb_sgp.draftkings import price_selection_set
    client = FakeClient(FakeResp(422, {"statusCode": 400, "description": "no"}))
    assert price_selection_set(client, ["a", "b"]) is None


def test_price_selection_set_combinability_restrictions_returns_none():
    from mlb_sgp.draftkings import price_selection_set
    client = FakeClient(FakeResp(200, {
        "combinabilityRestrictions": [{"reason": "NonCombinable"}],
        "bets": [{"trueOdds": 3.0, "selectionsMapped": [{}, {}]}],
    }))
    assert price_selection_set(client, ["a", "b"]) is None


def test_price_selection_set_no_bets_or_exception_returns_none():
    from mlb_sgp.draftkings import price_selection_set
    client = FakeClient(FakeResp(200, {"bets": []}))
    assert price_selection_set(client, ["a", "b"]) is None

    class BoomSession:
        def post(self, *a, **k):
            raise RuntimeError("network down")

    class BoomClient:
        session = BoomSession()

    assert price_selection_set(BoomClient(), ["a", "b"]) is None
