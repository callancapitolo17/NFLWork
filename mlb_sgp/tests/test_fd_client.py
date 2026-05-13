"""Unit tests for FanDuelClient — fixture-based, no network."""
import json
from pathlib import Path

import pytest

from mlb_sgp.fd_client import FanDuelClient, Event, Runner

FIXTURES = Path(__file__).parent / "fixtures"


def test_event_dataclass_has_required_fields():
    e = Event(event_id="abc", home_team="Yankees", away_team="Red Sox",
              start_time="2026-05-12T22:00:00Z")
    assert e.event_id == "abc"
    assert e.home_team == "Yankees"
    assert e.away_team == "Red Sox"
    assert e.start_time == "2026-05-12T22:00:00Z"


def test_runner_dataclass_has_required_fields():
    r = Runner(runner_id="42", market_id="m1", name="Over 9.5",
               line=9.5, american_odds=-110)
    assert r.runner_id == "42"
    assert r.market_id == "m1"
    assert r.name == "Over 9.5"
    assert r.line == 9.5
    assert r.american_odds == -110


def test_list_events_parses_fixture():
    """Verifies list_events parses a captured FD scan response.

    Fixture top-level keys: facets, results, attachments. Real event payloads
    live under attachments.events keyed by eventId. fetch_fd_events uses a
    recursive walk() to find dicts with eventId + openDate + name and rebuilds
    the away/home from the 'Away (P Pitcher) @ Home (P Pitcher)' name format.
    """
    with open(FIXTURES / "fd_events.json") as f:
        fixture = json.load(f)

    class FakeResponse:
        status_code = 200

        def json(self):
            return fixture

        def raise_for_status(self):
            pass

    class FakeSession:
        # fetch_fd_events uses session.post (it's a search endpoint)
        def post(self, url, **kwargs):
            return FakeResponse()

        def get(self, url, **kwargs):
            return FakeResponse()

    client = FanDuelClient.__new__(FanDuelClient)  # bypass __init__
    client.session = FakeSession()
    client.verbose = False

    events = client.list_events()
    assert len(events) > 0
    assert all(isinstance(e, Event) for e in events)
    assert all(e.event_id and e.home_team and e.away_team for e in events)
    # Spot-check the canonical-name passthrough — FD has 'St. Louis Cardinals'
    # and 'Athletics' for the first event in the fixture
    names = {(e.away_team, e.home_team) for e in events}
    assert ("St. Louis Cardinals", "Athletics") in names


def test_fetch_event_runners_parses_fixture():
    """Verifies fetch_event_runners walks an FD event-page payload to a flat
    list of Runner dataclasses.

    fd_runners.json is debug output from the legacy scraper, NOT a raw
    event-page response. The raw market structures live under raw_markets_sample.
    We synthesize a minimal event-page payload by nesting those markets, then
    verify the walk + flatten produces Runner objects with int american_odds
    parsed from winRunnerOdds.americanDisplayOdds.americanOdds.
    """
    with open(FIXTURES / "fd_runners.json") as f:
        debug = json.load(f)
    raw_markets = debug.get("raw_markets_sample") or []
    assert len(raw_markets) > 0, "fixture missing raw_markets_sample"

    # FD event-page payload: any structure containing the markets list will
    # work because the client uses a recursive walk(). Wrap them in a typical
    # 'eventPage.markets' shape.
    fake_payload = {"eventPage": {"markets": raw_markets}}

    class FakeResponse:
        status_code = 200
        text = ""

        def __init__(self, body):
            self._body = body
            self.text = json.dumps(body)

        def json(self):
            return self._body

        def raise_for_status(self):
            pass

    class FakeSession:
        def get(self, url, **kwargs):
            return FakeResponse(fake_payload)

    client = FanDuelClient.__new__(FanDuelClient)
    client.session = FakeSession()
    client.verbose = False

    runners = client.fetch_event_runners(event_id="dummy")

    assert len(runners) > 0, "fixture should yield many runners"
    assert all(isinstance(r, Runner) for r in runners)
    # FD's americanOdds is natively int — no string parsing on the way in
    assert all(isinstance(r.american_odds, int) for r in runners)
    # Every runner must have a non-empty id and market id
    assert all(r.runner_id for r in runners)
    assert all(r.market_id for r in runners)
    # FD prices can be positive or negative; the player-prop fixture is all
    # positive (long-tail player props), so we only require ints here
    assert any(r.american_odds > 0 for r in runners)
