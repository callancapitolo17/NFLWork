"""DraftKings on-demand pricing through the SGP-builder GET (2026-10-04).

calculateBets — the betslip's price call — has been Akamai-gated to a real
Chrome since ~2026-08-20 (issue #102). DK's SGP-builder widget prices the
same parlays through an ungated GET that takes a base plus ONE candidate
market and answers for every outcome of that market. This file pins the
pieces that make that safe to run from the maker:

- ``dk_market_id_of``: the candidate market is derived from the last leg's
  selection id, and an id we do not recognise must fail, never guess.
- ``parse_added_leg_odds``: DK silently drops a base leg it cannot combine
  and prices the rest — that must come back as a decline.
- ``PriceLane``: one call at a time, >= 0.5 s apart (a concurrent burst is
  what drew Akamai's rule onto calculateBets), and a sibling partition cell
  reuses the answer instead of a second call.
- ``fetch_selection_ids`` carries DK's straight odds, because the builder's
  1-leg price is not the straight price.

No network: every test injects a fake session, clock or fetch.
"""
from __future__ import annotations

import sys
import threading
from pathlib import Path

import pytest

_MLB_SGP_DIR = Path(__file__).resolve().parents[1]
if str(_MLB_SGP_DIR) not in sys.path:
    sys.path.insert(0, str(_MLB_SGP_DIR))

from mlb_sgp import draftkings                                   # noqa: E402
from mlb_sgp._shared import FetchCounters                        # noqa: E402
from mlb_sgp.dk_client import (DraftKingsClient, PriceLane,      # noqa: E402
                               parse_added_leg_odds)


# --------------------------------------------------------------------------- #
# dk_market_id_of                                                              #
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("selection_id, market_id", [
    ("0ML86559201_1", "1_86559201"),         # Moneyline
    ("0HC86559201N150_1", "2_86559201"),     # Run Line
    ("0OU86559201O850_1", "3_86559201"),     # Total
    ("0HC86560034P250_3", "2_86560034"),     # Run Line Alternate
    ("0OU86559661U450_3", "3_86559661"),     # Total Runs - 1st 5 Innings
    ("0OU86560036O50_1", "3_86560036"),      # Runs - 1st Inning
])
def test_market_id_matches_the_payload(selection_id, market_id):
    """Pairs read off DK's parlays payload for ATL@LAD, 2026-10-03."""
    assert draftkings.dk_market_id_of(selection_id) == market_id


@pytest.mark.parametrize("selection_id", [
    "0QA375313183#2306244131_13L84240Q12Q21",    # a prop-style market
    "0HC_1", "", "garbage",
])
def test_unknown_selection_shape_raises(selection_id):
    with pytest.raises(ValueError, match="0ML/0HC/0OU"):
        draftkings.dk_market_id_of(selection_id)


def test_unknown_selection_shape_is_a_parse_failure_not_transport():
    """A guessed market would price the wrong leg; the failure must point at
    the parser, not at DK's endpoint."""
    client = DraftKingsClient.__new__(DraftKingsClient)
    client.session = None                    # never reached
    counters = FetchCounters("draftkings", "on_demand")
    assert draftkings.price_selection_set(
        client, ["0HC1N150_1", "0QA9#1"], "evt1", counters=counters) is None
    snap = counters.snapshot()
    assert (snap.parse_failures, snap.transport_errors) == (1, 0)


# --------------------------------------------------------------------------- #
# parse_added_leg_odds                                                         #
# --------------------------------------------------------------------------- #

def _body(mapped, not_mapped=(), candidates=None):
    return {"data": {
        "selectionsMapped": list(mapped),
        "selectionsNotMapped": list(not_mapped),
        "compatibleMarkets": [{"id": "3_1", "selections": [
            {"id": sid, "trueOdds": odds}
            for sid, odds in (candidates or {}).items()]}],
    }}


def test_whole_base_returns_every_candidate_price():
    body = _body(["a"], candidates={"over": 3.5, "under": 3.76})
    assert parse_added_leg_odds(body, frozenset({"a"})) == {
        "over": 3.5, "under": 3.76}


def test_dropped_base_leg_prices_nothing():
    """Live 2026-10-03: both sides of one run line in the base -> DK drops one
    (IncompatibleSelections) and prices the other as if it were the base."""
    body = _body(["a"], not_mapped=[{"selectionId": "b",
                                     "removalReasonCode": "IncompatibleSelections"}],
                 candidates={"over": 3.5})
    assert parse_added_leg_odds(body, frozenset({"a", "b"})) == {}


def test_mapped_set_must_equal_the_base():
    body = _body(["a"], candidates={"over": 3.5})
    assert parse_added_leg_odds(body, frozenset({"a", "b"})) == {}


def test_unpriced_candidate_is_left_out():
    body = _body(["a"], candidates={"over": None, "under": 3.76})
    assert parse_added_leg_odds(body, frozenset({"a"})) == {"under": 3.76}


# --------------------------------------------------------------------------- #
# PriceLane                                                                    #
# --------------------------------------------------------------------------- #

class _FakeClock:
    def __init__(self):
        self.now = 1000.0
        self.sleeps: list[float] = []

    def __call__(self) -> float:
        return self.now

    def sleep(self, seconds: float) -> None:
        self.sleeps.append(seconds)
        self.now += seconds


def _lane(clock, *, max_wait_sec=1.0):
    return PriceLane(min_gap_sec=0.5, max_wait_sec=max_wait_sec,
                     answer_ttl_sec=2.0, clock=clock, sleep=clock.sleep)


def test_same_key_inside_the_ttl_reuses_the_answer():
    clock, fetches = _FakeClock(), []
    lane = _lane(clock)
    fetch = lambda: fetches.append(1) or {"x": 2.0}
    assert lane.answer("k", fetch) == {"x": 2.0}
    clock.now += 1.9
    assert lane.answer("k", fetch) == {"x": 2.0}
    assert len(fetches) == 1


def test_same_key_after_the_ttl_fetches_again():
    clock, fetches = _FakeClock(), []
    lane = _lane(clock)
    fetch = lambda: fetches.append(1) or {"x": 2.0}
    lane.answer("k", fetch)
    clock.now += 2.0
    lane.answer("k", fetch)
    assert len(fetches) == 2


def test_a_decline_is_remembered_too():
    """``{}`` is DK's answer, not a miss: the sibling cell must not re-ask."""
    clock, fetches = _FakeClock(), []
    lane = _lane(clock)
    fetch = lambda: fetches.append(1) or {}
    assert lane.answer("k", fetch) == {}
    assert lane.answer("k", fetch) == {}
    assert len(fetches) == 1


def test_calls_are_spaced_by_the_minimum_gap():
    clock = _FakeClock()
    lane = _lane(clock)
    lane.answer("k1", lambda: {})
    clock.now += 0.2
    lane.answer("k2", lambda: {})
    assert clock.sleeps == [pytest.approx(0.3)]


def test_first_call_does_not_wait():
    clock = _FakeClock()
    _lane(clock).answer("k1", lambda: {})
    assert clock.sleeps == []


def test_a_failed_fetch_is_not_remembered():
    clock = _FakeClock()
    lane = _lane(clock)

    def boom():
        raise RuntimeError("wire")

    with pytest.raises(RuntimeError):
        lane.answer("k", boom)
    assert lane.answer("k", lambda: {"x": 2.0}) == {"x": 2.0}


def _hold_lane(lane):
    """Start a fetch that holds the lane until the returned event is set."""
    entered, release = threading.Event(), threading.Event()

    def slow_fetch():
        entered.set()
        release.wait(5)
        return {}

    holder = threading.Thread(target=lane.answer, args=("held", slow_fetch))
    holder.start()
    assert entered.wait(5)
    return holder, release


def test_busy_lane_returns_none_after_the_wait_budget():
    lane = PriceLane(min_gap_sec=0.0, max_wait_sec=0.05, answer_ttl_sec=2.0)
    holder, release = _hold_lane(lane)
    try:
        assert lane.answer("other", lambda: {"x": 2.0}) is None
    finally:
        release.set()
        holder.join(5)


def test_remembered_answer_does_not_queue_behind_a_busy_lane():
    """A sibling cell's answer is already in hand; it must not wait behind
    another RFQ's wire call."""
    lane = PriceLane(min_gap_sec=0.0, max_wait_sec=0.05, answer_ttl_sec=2.0)
    lane.answer("k", lambda: {"x": 2.0})
    holder, release = _hold_lane(lane)
    try:
        assert lane.answer("k", lambda: {"x": 9.9}) == {"x": 2.0}
    finally:
        release.set()
        holder.join(5)


# --------------------------------------------------------------------------- #
# Client: the two cells that differ only in the last leg share one GET         #
# --------------------------------------------------------------------------- #

class _Resp:
    def __init__(self, body):
        self.status_code = 200
        self._body = body
        self.text = ""
        self.headers = {}

    def json(self):
        return self._body


class _CountingSession:
    def __init__(self, body):
        self._body = body
        self.gets = 0

    def get(self, url, **kwargs):
        self.gets += 1
        return _Resp(self._body)


def test_sibling_cells_share_one_wire_call():
    base, over, under = "0HC100N150_1", "0OU200O850_1", "0OU200U850_3"
    client = DraftKingsClient.__new__(DraftKingsClient)
    client.session = _CountingSession(
        _body([base], candidates={over: 3.5, under: 3.76}))
    assert draftkings.price_selection_set(client, [base, over], "evt1") == 3.5
    assert draftkings.price_selection_set(client, [base, under], "evt1") == 3.76
    assert client.session.gets == 1


def test_base_order_does_not_split_the_shared_answer():
    """SGPService builds cells by flipping legs in place; the same base in a
    different order is still the same question to DK."""
    a, b, over = "0HC100N150_1", "0ML300_1", "0OU200O850_1"
    client = DraftKingsClient.__new__(DraftKingsClient)
    client.session = _CountingSession(_body([a, b], candidates={over: 6.0}))
    draftkings.price_selection_set(client, [a, b, over], "evt1")
    draftkings.price_selection_set(client, [b, a, over], "evt1")
    assert client.session.gets == 1


def test_a_lazily_built_client_shares_one_lane():
    client = DraftKingsClient.__new__(DraftKingsClient)
    assert client.price_lane is client.price_lane


# --------------------------------------------------------------------------- #
# fetch_selection_ids carries DK's straight odds                               #
# --------------------------------------------------------------------------- #

class _ParlaysResp:
    def __init__(self, text, payload):
        self.status_code = 200
        self.text = text
        self._payload = payload
        self.headers = {}

    def json(self):
        if isinstance(self._payload, Exception):
            raise self._payload
        return self._payload


class _ParlaysSession:
    def __init__(self, resp):
        self._resp = resp

    def get(self, url, **kwargs):
        return self._resp


def test_straight_odds_cover_live_selections_only():
    from scraper_draftkings_sgp import fetch_selection_ids
    payload = {"data": {"markets": [
        {"id": "3_86560036", "name": "Runs - 1st Inning", "selections": [
            {"id": "0OU86560036O50_1", "trueOdds": 2.25},
            {"id": "0OU86560036U50_3", "trueOdds": 1.625},
            {"id": "0OU86560036O150_1", "trueOdds": 3.5, "isDisabled": True},
        ]},
    ]}}
    text = "0OU86560036O50_1 0OU86560036U50_3 0OU86560036O150_1"
    out = fetch_selection_ids(_ParlaysSession(_ParlaysResp(text, payload)), "e1")
    assert out["straight_odds"] == {"0OU86560036O50_1": 2.25,
                                    "0OU86560036U50_3": 1.625}


def test_unreadable_payload_degrades_to_no_straight_odds(caplog):
    from scraper_draftkings_sgp import fetch_selection_ids
    resp = _ParlaysResp("0OU1O850_1", ValueError("not json"))
    out = fetch_selection_ids(_ParlaysSession(resp), "e1")
    assert out["straight_odds"] == {}
    assert "no straight odds" in caplog.text
