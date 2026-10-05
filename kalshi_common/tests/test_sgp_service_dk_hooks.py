"""DraftKings through the REAL SGPService on-demand hooks (2026-10-04).

The book-level tests (mlb_sgp/tests/test_dk_sgp_builder_price.py) pin the
GET, the lane and the parser. These pin what only shows up with the service
driving them: the hook hands DK its event id, a 2-leg partition costs two
wire calls even though SGPService fires cells 1..3 on concurrent threads,
and a lone leg (kalshi_rfi's 1st-inning total) prices off DK's straight odds
with no SGP call at all.

No network: DK's structure fetchers are faked and the client's session
answers from a table.
"""
from __future__ import annotations

from datetime import datetime, timezone

import pytest

from kalshi_common import fair_value
from kalshi_common.legset import CanonicalLeg
from kalshi_common.sgp_service import SGPService
from mlb_sgp._shared import GameRef

HOME_M15, AWAY_P15 = "0HC100N150_1", "0HC100P150_3"
OVER, UNDER = "0OU200O850_1", "0OU200U850_3"
I1_OVER, I1_UNDER = "0OU300O50_1", "0OU300U50_3"

# Cells as DK priced ATL@LAD home -1.5 / total 8.5 on 2026-10-03.
BUILDER_ANSWERS = {
    frozenset({HOME_M15}): {OVER: 3.5, UNDER: 3.76},
    frozenset({AWAY_P15}): {OVER: 3.74, UNDER: 3.08},
}

SELECTION_IDS = {
    "fg": {"spreads": {("N", 1.5, "1"): [HOME_M15], ("P", 1.5, "3"): [AWAY_P15]},
           "totals": {("O", 8.5): [OVER], ("U", 8.5): [UNDER]},
           "moneyline": {"1": [], "3": []}, "canonical": {"100", "200"}},
    "f5": None,
    "i1": {"spreads": {}, "totals": {("O", 0.5): [I1_OVER], ("U", 0.5): [I1_UNDER]},
           "moneyline": {"1": [], "3": []}, "canonical": {"300"}},
    "straight_odds": {HOME_M15: 1.99, AWAY_P15: 1.84, OVER: 1.98, UNDER: 1.85,
                      I1_OVER: 2.25, I1_UNDER: 1.625},
}


class _Resp:
    status_code = 200
    text = ""
    headers: dict = {}

    def __init__(self, body):
        self._body = body

    def json(self):
        return self._body


class _BuilderSession:
    """Answers DK's SGP-builder GET from BUILDER_ANSWERS; records each call."""

    def __init__(self):
        self.calls = []

    def get(self, url, params=None, headers=None, **kwargs):
        self.calls.append(params)
        base = frozenset(params["selections"].split(","))
        priced = BUILDER_ANSWERS.get(base, {})
        return _Resp({"data": {
            "selectionsMapped": sorted(base), "selectionsNotMapped": [],
            "compatibleMarkets": [{"id": params["marketCandidates"], "selections": [
                {"id": sid, "trueOdds": odds} for sid, odds in priced.items()]}],
        }})


GAME = GameRef(game_id="26OCT032005ATLLAD", home_team="Los Angeles Dodgers",
               away_team="Atlanta Braves",
               commence_time=datetime(2026, 10, 3, 20, 5, tzinfo=timezone.utc))


@pytest.fixture
def dk_service(monkeypatch):
    """An SGPService whose DK client is real but whose wire is faked."""
    import scraper_draftkings_sgp as legacy
    from mlb_sgp.dk_client import DraftKingsClient

    def fake_fetchers(st, profile=None, *, miss_cb=None):
        return {
            "fetch_dk_events": lambda session: [{"dk_event_id": "evt7"}],
            "fetch_main_market_nums": lambda session, eid: {},
            "fetch_selection_ids":
                lambda session, eid, nums, verbose=False: SELECTION_IDS,
        }

    monkeypatch.setattr(SGPService, "_dk_fetchers", staticmethod(fake_fetchers))
    monkeypatch.setattr(legacy, "match_events", lambda events, lines: events)

    def build(**kwargs):
        svc = SGPService(books=("draftkings",), health_db_path=None, **kwargs)
        client = DraftKingsClient.__new__(DraftKingsClient)
        client.session = _BuilderSession()
        svc._state["draftkings"].client = client
        return svc, client.session
    return build


def test_two_leg_partition_costs_two_wire_calls(dk_service):
    svc, session = dk_service()
    legs = [CanonicalLeg(GAME.game_id, "spread", -1.5, "home"),
            CanonicalLeg(GAME.game_id, "total", 8.5, "over")]
    result = svc.price_on_demand("draftkings", GAME, legs)
    assert result is not None and result.route == "partition"
    assert result.n_cells_priced == 4
    assert result.fair == pytest.approx(
        fair_value.devig_partition([3.5, 3.74, 3.76, 3.08], 2))
    # Cells 1..3 ran on concurrent threads; the lane still sent only one
    # call per distinct base, each with the event id the hook passed.
    assert len(session.calls) == 2
    assert {c["eventId"] for c in session.calls} == {"evt7"}
    assert {c["marketCandidates"] for c in session.calls} == {"3_200"}


def test_lone_leg_prices_off_straight_odds_without_an_sgp_call(dk_service):
    """kalshi_rfi's 1st-inning total: the builder's 1-leg price is not DK's
    straight price, so the fair must come from the straight two-way."""
    svc, session = dk_service(structure_ttl_sec=0.0,
                              single_leg_structure_fair=True)
    leg = CanonicalLeg(GAME.game_id, "total", 0.5, "over", "I1")
    result = svc.price_on_demand("draftkings", GAME, [leg])
    assert result is not None and result.route == "single_two_way"
    assert result.fair == pytest.approx(
        fair_value.devig_partition([2.25, 1.625], 1))
    assert session.calls == []
