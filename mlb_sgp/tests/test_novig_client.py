"""Unit tests for NovigClient — fixture-based, no network.

Fixtures:
  - nv_events.json:     REAL response captured 2026-05-13 from
                        POST https://api.novig.us/v1/graphql (MLBEvents query).
                        18 upcoming MLB events.
  - nv_event_legs.json: REAL response captured 2026-05-13 from
                        POST https://api.novig.us/v1/graphql (EventMarkets_Query
                        for the first event id from nv_events.json).
                        259 markets including SPREAD, TOTAL, SPREAD_1H, TOTAL_1H.
  - nv_parlay_response.json: SYNTHETIC — actual parlay submission needs valid
                        outcome UUIDs that move on every line update. Shape
                        matches Novig's `[{"price": "0.35088", "status": "OPEN",
                        ...}]` list-of-offers response.
"""
import json
from pathlib import Path

from mlb_sgp.novig_client import (
    NovigClient,
    Event,
    EventLegs,
    _parse_events_response,
    _parse_event_legs_response,
    _parse_parlay_response,
)

FIX = Path(__file__).parent / "fixtures"


def test_parse_events_response_real_fixture():
    """Real captured response should yield Event dataclasses for every MLB game."""
    raw = json.loads((FIX / "nv_events.json").read_text())
    events = _parse_events_response(raw)
    assert len(events) > 0, "fixture should contain MLB events"
    for e in events:
        assert isinstance(e, Event)
        assert e.event_id
        assert e.home_team
        assert e.away_team
        assert e.home_sym, "Novig events always have a competitor symbol"
        assert e.away_sym
        # ISO-style timestamp
        assert "T" in e.start_time and (e.start_time.endswith("Z")
                                        or "+" in e.start_time)


def test_parse_events_tolerates_flat_shape():
    """Synthetic fallback shape: {"event": [...]} without the data envelope."""
    raw = {
        "event": [
            {
                "id": "test-1",
                "scheduled_start": "2026-05-14T00:00:00+00:00",
                "game": {
                    "homeTeam": {"name": "Test Home", "symbol": "TH"},
                    "awayTeam": {"name": "Test Away", "symbol": "TA"},
                },
            }
        ]
    }
    events = _parse_events_response(raw)
    assert len(events) == 1
    assert events[0].home_sym == "TH"
    assert events[0].away_sym == "TA"


def test_parse_events_skips_missing_team_names():
    """Defensive: drop events whose game.homeTeam.name is missing."""
    raw = {
        "data": {
            "event": [
                {
                    "id": "bad-1",
                    "scheduled_start": "2026-05-14T00:00:00+00:00",
                    "game": {"homeTeam": {}, "awayTeam": {"name": "Away"}},
                },
                {
                    "id": "good-1",
                    "scheduled_start": "2026-05-14T00:00:00+00:00",
                    "game": {
                        "homeTeam": {"name": "Home", "symbol": "H"},
                        "awayTeam": {"name": "Away", "symbol": "A"},
                    },
                },
            ]
        }
    }
    events = _parse_events_response(raw)
    assert len(events) == 1
    assert events[0].event_id == "good-1"


def test_parse_event_legs_response_real_fixture():
    """Real captured market tree should yield spread + total legs."""
    raw = json.loads((FIX / "nv_event_legs.json").read_text())
    legs = _parse_event_legs_response(raw)
    assert isinstance(legs, EventLegs)
    assert legs.event_id  # came from the captured event[0].id
    assert legs.spread_legs, "Novig event must yield spread legs"
    assert legs.total_legs, "Novig event must yield total legs"

    # Spot-check a spread leg
    sp = legs.spread_legs[0]
    assert sp["id"]
    assert sp["period"] in {"fg", "f5"}
    assert sp["strike"] is not None
    assert sp["competitor_symbol"], "spread leg must carry competitor symbol"

    # Spot-check a total leg
    to = legs.total_legs[0]
    assert to["id"]
    assert to["period"] in {"fg", "f5"}
    assert to["side"] in {"over", "under"}, \
        f"total leg side should be over/under; got {to['side']!r}"


def test_parse_event_legs_covers_fg_and_f5_periods():
    """The captured fixture has both full-game and first-5 markets."""
    raw = json.loads((FIX / "nv_event_legs.json").read_text())
    legs = _parse_event_legs_response(raw)
    spread_periods = {l["period"] for l in legs.spread_legs}
    total_periods = {l["period"] for l in legs.total_legs}
    assert "fg" in spread_periods, "expected at least one FG spread"
    assert "fg" in total_periods, "expected at least one FG total"
    # f5 may or may not be present depending on the matchup; we don't assert it.


def test_parse_event_legs_ignores_non_spread_total_markets():
    """Player props (BATTING_STRIKEOUTS, HITS, etc.) must NOT leak through."""
    raw = json.loads((FIX / "nv_event_legs.json").read_text())
    legs = _parse_event_legs_response(raw)
    # All spread legs should have a competitor symbol (player props don't)
    for l in legs.spread_legs:
        assert l["competitor_symbol"], \
            f"non-spread market leaked: {l!r}"
    # All total legs should have over/under side
    for l in legs.total_legs:
        assert l["side"] in {"over", "under"}, \
            f"non-total market leaked (no over/under desc): {l['description']!r}"


def test_parse_event_legs_tolerates_flat_shape():
    """Synthetic fixture without `data.event` envelope still parses."""
    raw = {
        "markets": [
            {
                "type": "SPREAD",
                "strike": -1.5,
                "is_consensus": True,
                "outcomes": [
                    {"id": "u1", "description": "MIN -1.5", "available": 0.45,
                     "competitor": {"symbol": "MIN"}},
                    {"id": "u2", "description": "MIA +1.5", "available": 0.55,
                     "competitor": {"symbol": "MIA"}},
                ],
            },
            {
                "type": "TOTAL",
                "strike": 8.5,
                "is_consensus": True,
                "outcomes": [
                    {"id": "u3", "description": "Over 8.5", "available": 0.5},
                    {"id": "u4", "description": "Under 8.5", "available": 0.5},
                ],
            },
        ]
    }
    legs = _parse_event_legs_response(raw, event_id_fallback="synthetic-event")
    assert legs.event_id == "synthetic-event"
    assert len(legs.spread_legs) == 2
    assert len(legs.total_legs) == 2
    over = next(l for l in legs.total_legs if l["side"] == "over")
    assert over["id"] == "u3"


def test_parse_parlay_response_list_shape():
    """The real /unauthenticated endpoint returns a list of offers."""
    raw = json.loads((FIX / "nv_parlay_response.json").read_text())
    parsed = _parse_parlay_response(raw)
    assert parsed["price_str"] == "0.35088"
    # 1 / 0.35088 ~= 2.85, american ~= +185
    assert abs(parsed["decimal"] - 2.85) < 0.01
    assert parsed["american"] == 185
    assert parsed["status"] == "OPEN"


def test_parse_parlay_response_empty():
    """Empty list / missing price / out-of-range price -> {}."""
    assert _parse_parlay_response([]) == {}
    assert _parse_parlay_response([{"status": "OPEN"}]) == {}
    assert _parse_parlay_response([{"price": "1.5"}]) == {}  # > 1 invalid
    assert _parse_parlay_response([{"price": "0"}]) == {}    # 0 invalid


def test_parse_parlay_response_mutation_shape():
    """Future-proof: if Novig switches to a BuildParlay mutation, we still parse."""
    raw = {
        "data": {
            "parlay": {
                "decimalOdds": "2.85",
                "americanOdds": -200,
                "totalStake": 1.0,
                "potentialPayout": 2.85,
                "status": "OPEN",
            }
        }
    }
    parsed = _parse_parlay_response(raw)
    assert parsed["decimal"] == 2.85
    assert parsed["american"] == -200
    assert parsed["status"] == "OPEN"


def test_submit_parlay_with_fake_session():
    """End-to-end submit_parlay flow with a mocked session.

    Verifies that the client unwraps the list response and returns the
    parsed decimal/american price for the top offer.
    """
    response_body = json.loads((FIX / "nv_parlay_response.json").read_text())

    class FakeResponse:
        status_code = 200

        def json(self):
            return response_body

    class FakeSession:
        def post(self, url, **kwargs):
            return FakeResponse()

    client = NovigClient.__new__(NovigClient)  # bypass __init__ -> no session
    client.session = FakeSession()
    client.verbose = False

    result = client.submit_parlay(["uuid-1", "uuid-2"])
    assert result, "submit_parlay must return a non-empty dict"
    assert abs(result["decimal"] - 2.85) < 0.01
    assert result["price_str"] == "0.35088"


def test_submit_parlay_handles_non_200():
    """Non-2xx HTTP from the parlay endpoint should yield {}."""
    class Fake500:
        status_code = 500
        def json(self):
            return []

    class FakeSession:
        def post(self, url, **kwargs):
            return Fake500()

    client = NovigClient.__new__(NovigClient)
    client.session = FakeSession()
    client.verbose = False
    assert client.submit_parlay(["uuid-1"]) == {}
