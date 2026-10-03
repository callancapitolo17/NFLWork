"""Unit tests for NovigClient — fixture-based, no network.

Fixtures:
  - nv_trading_page_mlb.json: REAL response captured 2026-10-02 from
                        GET https://api.novig.us/nbx/v1/trading/MLB/page.
                        4 pregame games on 2026-10-03 plus the Featured
                        Parlays / Series / Futures sections.
  - nv_event_legs.json: REAL response captured 2026-05-13 from
                        POST https://api.novig.us/v1/graphql (EventMarkets_Query).
                        259 markets including SPREAD, TOTAL, SPREAD_1H, TOTAL_1H.
  - nv_parlay_response.json: SYNTHETIC — actual parlay submission needs valid
                        outcome UUIDs that move on every line update. Shape
                        matches Novig's `[{"price": "0.35088", "status": "OPEN",
                        ...}]` list-of-offers response.
"""
import json
from datetime import datetime, timezone
from pathlib import Path

import pytest

from mlb_sgp._shared import BookTransportError
from mlb_sgp.novig_client import (
    NovigClient,
    Event,
    EventLegs,
    _parse_events_response,
    _parse_event_legs_response,
    _parse_parlay_response,
)

FIX = Path(__file__).parent / "fixtures"

# Two hours before the fixture's first game (2026-10-03T17:00Z).
FIXTURE_NOW = datetime(2026, 10, 3, 15, 0, tzinfo=timezone.utc)
WINDOW_HOURS = 48


def _games_page(*cards):
    """A trading page whose Games section holds ``cards``."""
    return {"sections": [{"title": "Games",
                          "content": {"type": "components",
                                      "components": list(cards)}}]}


def _card(event_id, start="2026-10-03T20:00:00.000Z", status="OPEN_PREGAME",
          card_type="game_event_card", home=("Home Team", "HOM", "HOM"),
          away=("Away Team", "AWY", "AWY")):
    return {"type": card_type, "eventId": event_id, "scheduledStart": start,
            "eventStatus": status,
            "homeTeam": {"name": home[0], "symbol": home[1], "shortName": home[2]},
            "awayTeam": {"name": away[0], "symbol": away[1], "shortName": away[2]}}


def test_parse_events_response_real_fixture():
    """The captured MLB page yields one Event per game card."""
    raw = json.loads((FIX / "nv_trading_page_mlb.json").read_text())
    events = _parse_events_response(raw, now=FIXTURE_NOW,
                                    window_hours=WINDOW_HOURS)
    assert len(events) == 4
    for e in events:
        assert isinstance(e, Event)
        assert e.event_id
        assert e.home_team
        assert e.away_team
        assert e.home_sym and e.away_sym
        assert e.start_time.startswith("2026-10-0") and e.start_time.endswith("Z")
    cle = next(e for e in events if e.home_team == "Cleveland Guardians")
    assert cle.event_id == "01a0f758-1b00-79d2-a2ad-323e9c1dbc11"
    assert cle.away_team == "Chicago White Sox"
    assert cle.start_time == "2026-10-03T17:00:00.000Z"


def test_parse_events_uses_symbol_not_short_name():
    """Market-tree outcomes carry `symbol`; the White Sox are CHI there and
    CWS only in shortName, so shortName would orphan every White Sox leg."""
    raw = json.loads((FIX / "nv_trading_page_mlb.json").read_text())
    events = _parse_events_response(raw, now=FIXTURE_NOW,
                                    window_hours=WINDOW_HOURS)
    cle = next(e for e in events if e.home_team == "Cleveland Guardians")
    assert (cle.home_sym, cle.away_sym) == ("CLE", "CHI")


def test_parse_events_keeps_only_pregame_games_inside_the_window():
    """Same filter the retired GraphQL WHERE clause applied server-side."""
    raw = _games_page(
        _card("keep"),
        _card("live", status="OPEN_LIVE"),
        _card("started", start="2026-10-03T14:00:00.000Z"),
        _card("too-far", start="2026-10-06T20:00:00.000Z"),
        _card("series", card_type="future_event_card"),
    )
    events = _parse_events_response(raw, now=FIXTURE_NOW,
                                    window_hours=WINDOW_HOURS)
    assert [e.event_id for e in events] == ["keep"]


def test_parse_events_reads_game_cards_in_any_section():
    """Cards are keyed on their type: renaming the "Games" section (display
    text) must not turn a full slate into an apparent off-day."""
    raw = {"sections": [{"title": "Upcoming Matchups",
                         "content": {"type": "components",
                                     "components": [_card("renamed")]}}]}
    events = _parse_events_response(raw, now=FIXTURE_NOW,
                                    window_hours=WINDOW_HOURS)
    assert [e.event_id for e in events] == ["renamed"]


def test_parse_events_dedupes_a_game_listed_in_two_sections():
    """Two Events for one game make the on-demand matcher decline it as
    ambiguous, so a game shown twice must come back once."""
    raw = _games_page(_card("dup"))
    raw["sections"].append({"title": "Featured",
                            "content": {"type": "components",
                                        "components": [_card("dup")]}})
    events = _parse_events_response(raw, now=FIXTURE_NOW,
                                    window_hours=WINDOW_HOURS)
    assert [e.event_id for e in events] == ["dup"]


def test_parse_events_off_season_page_is_empty():
    """An off-season league answers {"sections": []} — valid, no games."""
    assert _parse_events_response({"sections": []}, now=FIXTURE_NOW,
                                  window_hours=WINDOW_HOURS) == []


def test_parse_events_raises_on_unrecognised_page():
    """No `sections` list = a page we no longer understand: the book is
    down, not empty (the issue #33 contract)."""
    with pytest.raises(BookTransportError) as exc:
        _parse_events_response({"data": {"event": []}}, now=FIXTURE_NOW,
                               window_hours=WINDOW_HOURS)
    assert exc.value.stage == "events"


def test_parse_events_skips_missing_team_names():
    """Defensive: drop cards whose homeTeam.name is missing."""
    bad = _card("bad-1", home=(None, "H", "H"))
    raw = _games_page(bad, _card("good-1"))
    events = _parse_events_response(raw, now=FIXTURE_NOW,
                                    window_hours=WINDOW_HOURS)
    assert [e.event_id for e in events] == ["good-1"]


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


def test_parse_parlay_response_keeps_each_legs_vendor():
    """Novig relays every parlay from a vendor book and names it per leg
    (shape captured 2026-10-02, trimmed)."""
    raw = [{"price": "0.36000", "status": "Unfilled",
            "legs": [{"price": "0.57100", "vendor": "FANDUEL",
                      "outcomeId": "o-1"},
                     {"price": "0.60900", "vendor": "FANDUEL",
                      "outcomeId": "o-2"}]}]
    parsed = _parse_parlay_response(raw)
    assert parsed["vendors"] == ["FANDUEL", "FANDUEL"]
    assert parsed["decimal"] == pytest.approx(2.7778)


def test_parse_parlay_response_leg_without_vendor_reads_none():
    raw = [{"price": "0.36000", "legs": [{"vendor": "BETMGM"}, {"price": "0.5"}]}]
    assert _parse_parlay_response(raw)["vendors"] == ["BETMGM", None]


def test_parse_parlay_response_without_legs_names_no_vendors():
    raw = json.loads((FIX / "nv_parlay_response.json").read_text())
    assert _parse_parlay_response(raw)["vendors"] == []


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
