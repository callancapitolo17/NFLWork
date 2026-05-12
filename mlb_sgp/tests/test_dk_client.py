"""Unit tests for DraftKingsClient — fixture-based, no network."""
import json
from pathlib import Path
import pytest
from mlb_sgp.dk_client import DraftKingsClient, Event, Market, Selection

FIXTURES = Path(__file__).parent / "fixtures"


def test_event_dataclass_has_required_fields():
    e = Event(event_id="123", home_team="Yankees", away_team="Red Sox",
              start_time="2026-05-12T22:00:00Z")
    assert e.event_id == "123"
    assert e.home_team == "Yankees"
    assert e.away_team == "Red Sox"
    assert e.start_time == "2026-05-12T22:00:00Z"


def test_selection_dataclass_has_required_fields():
    s = Selection(selection_id="0HC1234N150_1", market_id="1234",
                  name="Yankees -1.5", line=-1.5, american_odds=120)
    assert s.line == -1.5
    assert s.american_odds == 120


def test_market_dataclass_has_required_fields():
    m = Market(market_id="1234", name="Run Line", subcategory="game_lines")
    assert m.subcategory == "game_lines"


def test_list_events_parses_fixture():
    """Verifies list_events parses a captured DK league API response."""
    fixture_path = FIXTURES / "dk_league_response.json"
    with open(fixture_path) as f:
        fake_response_body = json.load(f)

    class FakeResponse:
        status_code = 200

        def json(self):
            return fake_response_body

        def raise_for_status(self):
            pass  # 200 OK — nothing to raise

    class FakeSession:
        def get(self, url, **kwargs):
            return FakeResponse()

    client = DraftKingsClient.__new__(DraftKingsClient)  # bypass __init__
    client.session = FakeSession()
    client.verbose = False

    events = client.list_events()
    assert len(events) > 0
    assert all(isinstance(e, Event) for e in events)
    assert all(e.event_id and e.home_team and e.away_team for e in events)
