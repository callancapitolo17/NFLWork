"""Unit tests for NovigClient — fixture-based, no network.

Fixtures:
  - nv_trading_page_mlb.json: REAL response captured 2026-10-02 from
                        GET https://api.novig.us/nbx/v1/trading/MLB/page.
                        4 pregame games on 2026-10-03 plus the Featured
                        Parlays / Series / Futures sections.
  - nv_parlay_response.json: SYNTHETIC — actual parlay submission needs valid
                        outcome UUIDs that move on every line update. Shape
                        matches Novig's `[{"price": "0.35088", "status": "OPEN",
                        ...}]` list-of-offers response.
"""
import json
import logging
from datetime import datetime, timezone
from pathlib import Path

import pytest

from mlb_sgp._shared import BookTransportError
from mlb_sgp.novig_client import (
    EVENT_WINDOW_HOURS,
    NovigClient,
    Event,
    _parse_events_response,
    _parse_parlay_response,
)

FIX = Path(__file__).parent / "fixtures"

# Two hours before the fixture's first game (2026-10-03T17:00Z).
FIXTURE_NOW = datetime(2026, 10, 3, 15, 0, tzinfo=timezone.utc)
WINDOW_HOURS = EVENT_WINDOW_HOURS


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
        assert e.start_time.startswith("2026-10-0")
        assert e.start_time.endswith("+00:00")
    cle = next(e for e in events if e.home_team == "Cleveland Guardians")
    assert cle.event_id == "01a0f758-1b00-79d2-a2ad-323e9c1dbc11"
    assert cle.away_team == "Chicago White Sox"
    assert cle.start_time == "2026-10-03T17:00:00+00:00"


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


def test_parse_events_reads_game_cards_nested_in_subsections():
    """The page already nests Futures as content.type="subsections"; Games
    regrouped by day the same way must still parse."""
    raw = {"sections": [{"title": "Games", "content": {
        "type": "subsections",
        "subSections": [{"title": "Saturday", "components": [_card("sat")]},
                        {"title": "Sunday", "components": [_card(
                            "sun", start="2026-10-04T20:00:00.000Z")]}]}}]}
    events = _parse_events_response(raw, now=FIXTURE_NOW,
                                    window_hours=WINDOW_HOURS)
    assert [e.event_id for e in events] == ["sat", "sun"]


def test_parse_events_normalises_start_time_to_utc():
    """match_events buckets on start_time[:13] as a UTC hour, so an offset
    timestamp must come out in UTC."""
    raw = _games_page(_card("offset", start="2026-10-03T13:00:00.000-04:00"))
    events = _parse_events_response(raw, now=FIXTURE_NOW,
                                    window_hours=WINDOW_HOURS)
    assert events[0].start_time == "2026-10-03T17:00:00+00:00"


def test_parse_events_reads_a_naive_start_as_utc():
    raw = _games_page(_card("naive", start="2026-10-03T20:00:00"))
    events = _parse_events_response(raw, now=FIXTURE_NOW,
                                    window_hours=WINDOW_HOURS)
    assert events[0].start_time == "2026-10-03T20:00:00+00:00"


def test_parse_events_warns_and_skips_an_unparseable_start(caplog):
    raw = _games_page(_card("epoch", start=1791046800000), _card("ok"))
    with caplog.at_level(logging.WARNING, logger="mlb_sgp.novig_client"):
        events = _parse_events_response(raw, now=FIXTURE_NOW,
                                        window_hours=WINDOW_HOURS)
    assert [e.event_id for e in events] == ["ok"]
    assert "unparseable scheduledStart" in caplog.text


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
