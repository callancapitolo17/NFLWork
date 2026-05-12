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


def test_fetch_event_markets_parses_fixture():
    """Verifies fetch_event_markets parses captured DK subcategory responses.

    Fixture shape (dk_event_markets.json):
      { "event_id": str,
        "subcat_market_ids": { "4519": [[mid, name], ...],
                               "15628": [[mid, name], ...] },
        "parlays_payload": {...} }

    The two subcat IDs map to game_lines (4519) and alt_lines (15628).
    """
    with open(FIXTURES / "dk_event_markets.json") as f:
        fixture = json.load(f)

    class FakeResponse:
        status_code = 200

        def __init__(self, body):
            self._body = body

        def json(self):
            return self._body

        def raise_for_status(self):
            pass

    class FakeSession:
        """Returns the matching subcat slice based on the marketsQuery filter."""
        def __init__(self, fixture):
            self.fixture = fixture

        def get(self, url, **kwargs):
            # _fetch_subcat_markets passes the subcategory id via the
            # marketsQuery param. Inspect kwargs["params"] to figure out
            # which subcat is being requested.
            params = kwargs.get("params", {}) or {}
            q = params.get("marketsQuery", "")
            for subcat_id, mlist in self.fixture.get("subcat_market_ids", {}).items():
                if f"subCategoryId eq '{subcat_id}'" in q:
                    return FakeResponse({
                        "markets": [{"id": m[0], "name": m[1]} for m in mlist]
                    })
            return FakeResponse({"markets": []})

    client = DraftKingsClient.__new__(DraftKingsClient)
    client.session = FakeSession(fixture)
    client.verbose = False

    markets = client.fetch_event_markets(event_id=fixture["event_id"])
    # game_lines subcat has 3 markets (ML, RL, Total) per fixture
    assert len(markets) > 0
    assert all(isinstance(m, Market) for m in markets)
    assert all(m.market_id and m.name for m in markets)
    # Verify subcategory labels are applied
    subcats = {m.subcategory for m in markets}
    assert "game_lines" in subcats
    # alt_lines may be empty in some captures; subcat key must still be valid
    assert subcats.issubset({"game_lines", "alt_lines"})


def test_fetch_event_selections_returns_priced_rows():
    """Verifies selections include american_odds parsed correctly from displayOdds.american.

    The DK parlays endpoint returns a payload with shape:
      { "data": { "markets": [
          { "id": str, "name": str, "selections": [
              { "id": str, "marketId": str, "displayOdds": {"american": "+260"},
                "outcomeType": "Over"/"Under"/..., "points": 1.5,
                "name": str, "isDisabled": bool, "status": "Active" },
              ... ]},
          ... ]}}

    Crucially, DK uses the Unicode minus sign U+2212 ("−") for negative odds,
    not the ASCII hyphen "-". The parser must handle both.
    """
    with open(FIXTURES / "dk_event_markets.json") as f:
        fixture = json.load(f)

    parlays_body = fixture["parlays_payload"]

    class FakeResponse:
        status_code = 200
        text = ""

        def __init__(self, body):
            self._body = body
            import json as _json
            self.text = _json.dumps(body)

        def json(self):
            return self._body

        def raise_for_status(self):
            pass

    class FakeSession:
        def __init__(self, body):
            self._body = body

        def get(self, url, **kwargs):
            return FakeResponse(self._body)

    client = DraftKingsClient.__new__(DraftKingsClient)
    client.session = FakeSession(parlays_body)
    client.verbose = False

    selections = client.fetch_event_selections(event_id=fixture["event_id"])

    assert len(selections) > 0, "fixture has 149 markets, should produce many selections"
    assert all(isinstance(s, Selection) for s in selections)
    # Every selection must have an integer american_odds (parsed from string)
    assert all(isinstance(s.american_odds, int) for s in selections)
    # Non-empty selection_id
    assert all(s.selection_id for s in selections)
    # Non-empty market_id
    assert all(s.market_id for s in selections)
    # At least one negative-priced selection — must come from a Unicode-minus
    # ("−370") string in the fixture; proves the parser handles U+2212.
    assert any(s.american_odds < 0 for s in selections)
    # And at least one positive (e.g. "+260")
    assert any(s.american_odds > 0 for s in selections)
