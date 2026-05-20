"""Unit tests for ProphetXClient — fixture-based, no network.

Fixtures are real responses captured 2026-05-13 from the live ProphetX API:
  - px_events.json:        GET /trade/public/api/v1/tournaments
  - px_event_markets.json: GET /trade/public/api/v2/events/10077925/markets
  - px_parlay_response.json: hand-crafted minimal RFQ response (3 offer tiers)

The parlay-response fixture is synthetic because RFQ submission requires
auth + valid lineIDs which expire on every line move. The offer-ladder
shape it tests (teaser tier at $50, then real liquidity at $200+) is
documented in the ProphetX SGP memory note.
"""
import json
from pathlib import Path

from mlb_sgp.prophetx_client import (
    ProphetXClient,
    Event,
    Market,
    SelectionLeg,
    _parse_events_response,
    _parse_event_markets,
    _pick_offer,
)

FIX = Path(__file__).parent / "fixtures"


def test_parse_events_response():
    """Real captured fixture should yield Event dataclasses for every MLB game."""
    raw = json.loads((FIX / "px_events.json").read_text())
    events = _parse_events_response(raw)
    assert len(events) > 0, "fixture should contain MLB events"
    for e in events:
        assert isinstance(e, Event)
        assert e.event_id
        assert e.home_team
        assert e.away_team
    # Spot-check a known game from the capture
    names = {(e.away_team, e.home_team) for e in events}
    assert any("Reds" in h for (_, h) in names), \
        f"expected a Reds home game in fixture; got {names}"


def test_parse_events_filters_non_mlb_sports():
    """Only Baseball events under MLB-named tournaments should pass through.

    The fixture has Tennis, Soccer, Golf, NBA, NHL etc.; none should appear
    in the parsed list. (We also assert no soccer team names slip in via
    "Major League Soccer" — important: the tournament filter must NOT match
    the bare substring "Major League".)
    """
    raw = json.loads((FIX / "px_events.json").read_text())
    events = _parse_events_response(raw)
    team_names = {e.home_team for e in events} | {e.away_team for e in events}
    # MLS team names that would slip in if the filter was too loose
    mls_giveaways = {"Inter Miami CF", "New York City FC", "Los Angeles FC",
                     "Charlotte FC", "FC Dallas"}
    leaked = team_names & mls_giveaways
    assert not leaked, f"MLS teams leaked into MLB events: {leaked}"


def test_parse_event_markets():
    """Real captured market tree should yield Run Line + Total Runs markets."""
    raw = json.loads((FIX / "px_event_markets.json").read_text())
    markets = _parse_event_markets(raw)
    assert len(markets) > 0
    names = {m.name for m in markets}
    # ProphetX MLB FG markets always include spread + total
    assert any("Spread" in n or "Run Line" in n for n in names), \
        f"expected a spread/Run Line market; got {sorted(names)[:5]}"
    assert any("Total" in n or "Run Total" in n for n in names), \
        f"expected a total market; got {sorted(names)[:5]}"
    # Each Market must carry through raw marketLines so downstream resolvers
    # can read line/lineID per outcome
    run_line = next((m for m in markets if m.name == "Run Line"), None)
    assert run_line is not None
    assert len(run_line.selections) > 0, "Run Line should have at least one marketLine"


def test_pick_offer_filters_teasers():
    """offers[0] is often a $50 teaser; we want the first offer with stake >= 150.

    This is the policy memorialized in the ProphetX SGP memory note —
    interleaved low-stake teaser tiers are hidden by PX's UI and shouldn't
    drive our priced parlay decimal.
    """
    offers = [
        {"odds": 5.0, "stake": 50, "estimatedPrices": []},
        {"odds": 2.85, "stake": 200, "estimatedPrices": []},
        {"odds": 2.85, "stake": 500, "estimatedPrices": []},
    ]
    picked = _pick_offer(offers, min_stake=150)
    assert picked is not None
    assert picked["stake"] == 200
    assert picked["odds"] == 2.85


def test_pick_offer_falls_back_to_first_when_all_thin():
    """If nothing clears min_stake, return offers[0] — better thin price than none."""
    offers = [{"odds": 5.0, "stake": 50, "estimatedPrices": []}]
    picked = _pick_offer(offers, min_stake=150)
    assert picked is not None
    assert picked["stake"] == 50, "thin market — return first offer"


def test_pick_offer_empty_returns_none():
    assert _pick_offer([], min_stake=150) is None


def test_submit_parlay_rfq_with_fixture_response():
    """End-to-end RFQ flow with a mocked session.

    Verifies that submit_parlay_rfq:
      - posts the marketLines payload (we don't assert on it but the call shouldn't raise)
      - unwraps the data.offers envelope
      - picks the first $200 tier (not the $50 teaser) when min_offer_stake=150
      - returns used_fallback=False since we cleared the threshold
    """
    response_body = json.loads((FIX / "px_parlay_response.json").read_text())

    class FakeResponse:
        status_code = 200
        text = json.dumps(response_body)

        def json(self):
            return response_body

    class FakeSession:
        def post(self, url, **kwargs):
            return FakeResponse()

    client = ProphetXClient.__new__(ProphetXClient)  # bypass __init__
    client.session = FakeSession()
    client.verbose = False

    legs = [
        SelectionLeg(sport_event_id="10077925", market_id="256",
                     outcome_id="1714", line_id="abc123", line=-1.5),
        SelectionLeg(sport_event_id="10077925", market_id="257",
                     outcome_id="1800", line_id="def456", line=8.5),
    ]
    picked, used_fallback = client.submit_parlay_rfq(legs, min_offer_stake=150)
    assert picked is not None
    assert picked["stake"] == 200
    assert picked["odds"] == 2.85
    assert used_fallback is False
