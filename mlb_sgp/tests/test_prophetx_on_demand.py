"""ProphetX on-demand pricing (Phase 2): resolve_legs + price_selection_set.

Mirrors the inline-fixture style of test_prophetx_orchestrator.py /
test_draftkings_on_demand.py — hand-built PX market-tree dicts, no network,
fake duck-typed client for the wire call.

PX line-sign convention under test: each selection carries `line` from the
OUTCOME's own perspective, so a home -1.5 favorite appears as line=-1.5 on
the home selection and line=+1.5 on the away selection. CanonicalLeg lines
are signed home-perspective, so:
  side "home" -> pick home competitor at  leg.line
  side "away" -> pick away competitor at -leg.line
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from kalshi_common.legset import CanonicalLeg  # noqa: E402
from mlb_sgp._shared import american_to_decimal  # noqa: E402


HOME, AWAY = "New York Yankees", "Boston Red Sox"
HOME_ID, AWAY_ID = 111, 222
EVENT_ID = "ev900"


def _markets():
    """Raw-dict market list mirroring the markets_cache shape price_sgps
    builds from ProphetXClient.fetch_event_markets (id / name /
    marketLines / outcomes / selections)."""
    return [
        {
            "id": 10, "name": "Run Line",
            "marketLines": [
                {"outcomes": [
                    # Main line: home favored -1.5.
                    {"id": 1001, "lineID": "L-h-15", "line": -1.5,
                     "competitorId": HOME_ID, "name": "Yankees -1.5",
                     "odds": -115},
                    {"id": 1002, "lineID": "L-a+15", "line": 1.5,
                     "competitorId": AWAY_ID, "name": "Red Sox +1.5",
                     "odds": -105},
                    # Alt line the other way: home +1.5 dog.
                    {"id": 1003, "lineID": "L-h+15", "line": 1.5,
                     "competitorId": HOME_ID, "name": "Yankees +1.5",
                     "odds": -260},
                    {"id": 1004, "lineID": "L-a-15", "line": -1.5,
                     "competitorId": AWAY_ID, "name": "Red Sox -1.5",
                     "odds": 210},
                ]},
            ],
            "outcomes": [], "selections": [],
        },
        {
            "id": 20, "name": "Total Runs",
            "marketLines": [
                {"outcomes": [
                    {"id": 2001, "lineID": "L-o85", "line": 8.5,
                     "name": "Over 8.5", "odds": -110},
                    {"id": 2002, "lineID": "L-u85", "line": 8.5,
                     "name": "Under 8.5", "odds": -110},
                ]},
            ],
            "outcomes": [], "selections": [],
        },
        {
            "id": 30, "name": "Moneyline",
            "marketLines": [],
            "outcomes": [
                {"competitorId": HOME_ID}, {"competitorId": AWAY_ID},
            ],
            # ML order book: one best-first ladder per outcome (top-level
            # `selections`), each entry carrying odds + id + lineID.
            "selections": [
                [{"id": 3001, "lineID": "L-mlh", "line": 0.0,
                  "competitorId": HOME_ID, "odds": -140},
                 {"id": 3009, "lineID": "L-mlh2", "line": 0.0,
                  "competitorId": HOME_ID, "odds": -145}],
                [{"id": 3002, "lineID": "L-mla", "line": 0.0,
                  "competitorId": AWAY_ID, "odds": 120}],
            ],
        },
    ]


def _structure(markets=None):
    return {
        "event_id": EVENT_ID,
        "markets": _markets() if markets is None else markets,
        "home_id": HOME_ID,
        "away_id": AWAY_ID,
    }


def _leg_ids(rl):
    """(outcome_id, line_id) of a ResolvedLeg's chosen SelectionLeg ref."""
    return rl.ref.outcome_id, rl.ref.line_id


# ---------------------------------------------------------------------------
# resolve_legs — spreads (sign mirroring)
# ---------------------------------------------------------------------------
def test_resolve_home_spread_hit_with_opposite():
    from mlb_sgp.prophetx import resolve_legs
    legs = [CanonicalLeg("g1", "spread", -1.5, "home")]
    out = resolve_legs(_structure(), legs, HOME, AWAY)
    assert out is not None and len(out) == 1
    rl = out[0]
    # Chosen = home competitor at leg.line (-1.5).
    assert _leg_ids(rl) == ("1001", "L-h-15")
    assert rl.ref.sport_event_id == EVENT_ID
    assert rl.ref.market_id == "10"
    assert rl.ref.line == -1.5
    # Opposite = away competitor at the MIRRORED line (+1.5).
    assert (rl.opposite_ref.outcome_id, rl.opposite_ref.line_id) == \
        ("1002", "L-a+15")
    assert rl.opposite_ref.line == 1.5
    # Single odds captured american -> decimal.
    assert rl.single_decimal == american_to_decimal(-115)
    assert rl.opposite_decimal == american_to_decimal(-105)


def test_resolve_away_spread_mirrored_sign():
    from mlb_sgp.prophetx import resolve_legs
    # Same signed home line -1.5, AWAY side: away outcome carries +1.5.
    legs = [CanonicalLeg("g1", "spread", -1.5, "away")]
    out = resolve_legs(_structure(), legs, HOME, AWAY)
    assert out is not None
    assert _leg_ids(out[0]) == ("1002", "L-a+15")
    assert out[0].ref.line == 1.5
    assert out[0].opposite_ref.outcome_id == "1001"
    assert out[0].single_decimal == american_to_decimal(-105)
    assert out[0].opposite_decimal == american_to_decimal(-115)


def test_resolve_home_dog_positive_line():
    from mlb_sgp.prophetx import resolve_legs
    # Home +1.5 (positive signed home line): home outcome carries +1.5,
    # away outcome carries -1.5.
    out = resolve_legs(
        _structure(), [CanonicalLeg("g1", "spread", 1.5, "home")], HOME, AWAY)
    assert _leg_ids(out[0]) == ("1003", "L-h+15")
    assert out[0].opposite_ref.outcome_id == "1004"
    out = resolve_legs(
        _structure(), [CanonicalLeg("g1", "spread", 1.5, "away")], HOME, AWAY)
    assert _leg_ids(out[0]) == ("1004", "L-a-15")
    assert out[0].opposite_ref.outcome_id == "1003"
    assert out[0].single_decimal == american_to_decimal(210)


# ---------------------------------------------------------------------------
# resolve_legs — totals
# ---------------------------------------------------------------------------
def test_resolve_total_over_and_under():
    from mlb_sgp.prophetx import resolve_legs
    out = resolve_legs(
        _structure(), [CanonicalLeg("g1", "total", 8.5, "over")], HOME, AWAY)
    assert _leg_ids(out[0]) == ("2001", "L-o85")
    assert out[0].ref.market_id == "20"
    assert out[0].opposite_ref.outcome_id == "2002"   # under, SAME line
    assert out[0].single_decimal == american_to_decimal(-110)
    out = resolve_legs(
        _structure(), [CanonicalLeg("g1", "total", 8.5, "under")], HOME, AWAY)
    assert _leg_ids(out[0]) == ("2002", "L-u85")
    assert out[0].opposite_ref.outcome_id == "2001"


# ---------------------------------------------------------------------------
# resolve_legs — moneyline (best-first order book ladder)
# ---------------------------------------------------------------------------
def test_resolve_ml_best_ladder_entry_both_sides():
    from mlb_sgp.prophetx import resolve_legs
    out = resolve_legs(
        _structure(), [CanonicalLeg("g1", "ml", None, "home")], HOME, AWAY)
    # Best (FIRST) ladder entry, not the deeper 3009 one.
    assert _leg_ids(out[0]) == ("3001", "L-mlh")
    assert out[0].ref.market_id == "30"
    assert out[0].opposite_ref.outcome_id == "3002"
    assert out[0].single_decimal == american_to_decimal(-140)
    assert out[0].opposite_decimal == american_to_decimal(120)
    out = resolve_legs(
        _structure(), [CanonicalLeg("g1", "ml", None, "away")], HOME, AWAY)
    assert _leg_ids(out[0]) == ("3002", "L-mla")
    assert out[0].opposite_ref.outcome_id == "3001"


# ---------------------------------------------------------------------------
# resolve_legs — Market dataclass input (fetch_event_markets return shape)
# ---------------------------------------------------------------------------
def test_resolve_accepts_market_dataclasses():
    from mlb_sgp.prophetx import resolve_legs
    from mlb_sgp.prophetx_client import _parse_event_markets
    raw = {"data": {"markets": [
        {"id": m["id"], "name": m["name"], "marketLines": m["marketLines"],
         "outcomes": m["outcomes"], "selections": m["selections"]}
        for m in _markets()
    ]}}
    market_objs = _parse_event_markets(raw)   # list[Market], as the client returns
    out = resolve_legs(
        _structure(markets=market_objs),
        [CanonicalLeg("g1", "spread", -1.5, "home"),
         CanonicalLeg("g1", "ml", None, "away")],
        HOME, AWAY)
    assert out is not None
    assert _leg_ids(out[0]) == ("1001", "L-h-15")
    assert _leg_ids(out[1]) == ("3002", "L-mla")


# ---------------------------------------------------------------------------
# resolve_legs — multi-leg + failure routing
# ---------------------------------------------------------------------------
def test_resolve_multi_leg_order_preserved():
    from mlb_sgp.prophetx import resolve_legs
    legs = [
        CanonicalLeg("g1", "spread", -1.5, "home"),
        CanonicalLeg("g1", "total", 8.5, "over"),
        CanonicalLeg("g1", "ml", None, "away"),
    ]
    out = resolve_legs(_structure(), legs, HOME, AWAY)
    assert [rl.ref.outcome_id for rl in out] == ["1001", "2001", "3002"]


def test_resolve_chosen_side_miss_fails_whole_set():
    from mlb_sgp.prophetx import resolve_legs
    # Line not offered at all -> whole resolution fails.
    assert resolve_legs(
        _structure(), [CanonicalLeg("g1", "spread", -2.5, "home")],
        HOME, AWAY) is None
    # One good leg + one bad leg still fails as a whole.
    legs = [CanonicalLeg("g1", "total", 8.5, "over"),
            CanonicalLeg("g1", "total", 9.5, "over")]
    assert resolve_legs(_structure(), legs, HOME, AWAY) is None
    # Market missing entirely -> None.
    no_ml = _structure([m for m in _markets() if m["name"] != "Moneyline"])
    assert resolve_legs(
        no_ml, [CanonicalLeg("g1", "ml", None, "home")], HOME, AWAY) is None
    # ML side with an empty order book -> None (can't price it).
    mkts = _markets()
    mkts[2]["selections"] = [mkts[2]["selections"][1]]   # away ladder only
    assert resolve_legs(
        _structure(mkts), [CanonicalLeg("g1", "ml", None, "home")],
        HOME, AWAY) is None


def test_resolve_opposite_side_miss_yields_none_opposite():
    from mlb_sgp.prophetx import resolve_legs
    # Book one-sides the spread: drop the away +1.5 outcome.
    mkts = _markets()
    mkts[0]["marketLines"][0]["outcomes"] = [
        s for s in mkts[0]["marketLines"][0]["outcomes"]
        if s["id"] != 1002]
    out = resolve_legs(
        _structure(mkts), [CanonicalLeg("g1", "spread", -1.5, "home")],
        HOME, AWAY)
    assert out is not None
    assert out[0].ref.outcome_id == "1001"
    assert out[0].opposite_ref is None        # routes book to Route B
    assert out[0].opposite_decimal is None
    # Same for a one-sided total (no Under 8.5).
    mkts = _markets()
    mkts[1]["marketLines"][0]["outcomes"] = [
        s for s in mkts[1]["marketLines"][0]["outcomes"]
        if s["id"] != 2002]
    out = resolve_legs(
        _structure(mkts), [CanonicalLeg("g1", "total", 8.5, "over")],
        HOME, AWAY)
    assert out[0].opposite_ref is None
    # ML with only the chosen side's ladder.
    mkts = _markets()
    mkts[2]["selections"] = [mkts[2]["selections"][0]]   # home ladder only
    out = resolve_legs(
        _structure(mkts), [CanonicalLeg("g1", "ml", None, "home")],
        HOME, AWAY)
    assert out[0].ref.outcome_id == "3001"
    assert out[0].opposite_ref is None


def test_resolve_odds_absent_or_zero_gives_none_decimal():
    from mlb_sgp.prophetx import resolve_legs
    mkts = _markets()
    outs = mkts[0]["marketLines"][0]["outcomes"]
    del outs[0]["odds"]          # home -1.5: odds absent
    outs[1]["odds"] = 0          # away +1.5: odds zero
    out = resolve_legs(
        _structure(mkts), [CanonicalLeg("g1", "spread", -1.5, "home")],
        HOME, AWAY)
    assert out is not None                    # legs still resolve...
    assert out[0].ref.outcome_id == "1001"
    assert out[0].single_decimal is None      # ...but odds are None-safe
    assert out[0].opposite_decimal is None


def test_resolve_malformed_structure_returns_none():
    from mlb_sgp.prophetx import resolve_legs
    legs = [CanonicalLeg("g1", "spread", -1.5, "home")]
    assert resolve_legs(None, legs, HOME, AWAY) is None
    assert resolve_legs({}, legs, HOME, AWAY) is None
    assert resolve_legs("garbage", legs, HOME, AWAY) is None
    assert resolve_legs({"event_id": EVENT_ID, "markets": "nope",
                         "home_id": HOME_ID, "away_id": AWAY_ID},
                        legs, HOME, AWAY) is None
    # Unknown market_type / side never raise either.
    assert resolve_legs(_structure(),
                        [CanonicalLeg("g1", "runline", -1.5, "home")],
                        HOME, AWAY) is None
    assert resolve_legs(_structure(),
                        [CanonicalLeg("g1", "total", 8.5, "home")],
                        HOME, AWAY) is None
    assert resolve_legs(_structure(),
                        [CanonicalLeg("g1", "ml", None, "over")],
                        HOME, AWAY) is None


# ---------------------------------------------------------------------------
# price_selection_set — fake client wire tests
# ---------------------------------------------------------------------------
class FakeClient:
    """Duck-typed ProphetXClient: records the refs passed to
    submit_parlay_rfq and returns a canned (offer, used_fallback)."""

    def __init__(self, offer, used_fallback=False, raise_exc=None):
        self._offer = offer
        self._used_fallback = used_fallback
        self._raise = raise_exc
        self.last_legs = None

    def submit_parlay_rfq(self, legs):
        if self._raise is not None:
            raise self._raise
        self.last_legs = legs
        return self._offer, self._used_fallback


def _refs():
    from mlb_sgp.prophetx_client import SelectionLeg
    return [
        SelectionLeg(EVENT_ID, "10", "1001", "L-h-15", -1.5),
        SelectionLeg(EVENT_ID, "20", "2001", "L-o85", 8.5),
    ]


def test_price_selection_set_happy_path():
    from mlb_sgp.prophetx import price_selection_set
    client = FakeClient({"odds": 542, "stake": 200})
    refs = _refs()
    # offer["odds"] is AMERICAN parlay odds -> decimal, same as price_sgps.
    assert price_selection_set(client, refs) == american_to_decimal(542)
    assert client.last_legs == refs           # refs passed through unmodified


def test_price_selection_set_single_ref():
    from mlb_sgp.prophetx import price_selection_set
    client = FakeClient({"odds": -140, "stake": 500})
    assert price_selection_set(client, _refs()[:1]) == american_to_decimal(-140)


def test_price_selection_set_fallback_offer_still_counts():
    from mlb_sgp.prophetx import price_selection_set
    # used_fallback=True (thin market, teaser tier) is still a price.
    client = FakeClient({"odds": 615, "stake": 50}, used_fallback=True)
    assert price_selection_set(client, _refs()) == american_to_decimal(615)


def test_price_selection_set_none_offer_or_bad_odds_returns_none():
    from mlb_sgp.prophetx import price_selection_set
    assert price_selection_set(FakeClient(None), _refs()) is None
    assert price_selection_set(FakeClient({"stake": 200}), _refs()) is None
    assert price_selection_set(
        FakeClient({"odds": "garbage", "stake": 200}), _refs()) is None


def test_price_selection_set_never_raises():
    from mlb_sgp.prophetx import price_selection_set
    client = FakeClient(None, raise_exc=RuntimeError("network down"))
    assert price_selection_set(client, _refs()) is None
    assert price_selection_set(FakeClient({"odds": 542}), []) is None
