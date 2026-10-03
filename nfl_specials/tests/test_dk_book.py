"""DraftKingsBook pricing over DK's SGP-builder endpoint — no network."""
import pytest

from nfl_specials import dk_book
from nfl_specials.books import BookGame
from nfl_specials.dk_book import DraftKingsBook, parse_added_leg_odds

GAME = BookGame("evt1", home="SEA", away="LAC", game_start_time="2099-01-01T00:00:00Z")
MARKETS = [
    {"id": "m_q1", "name": "1st Quarter (3 Way)", "tags": ["SGP"],
     "selections": [{"id": "q1_sea"}, {"id": "q1_tie"}, {"id": "q1_lac"}]},
    {"id": "m_ml", "name": "Moneyline", "tags": ["SGP"],
     "selections": [{"id": "ml_sea"}, {"id": "ml_lac"}]},
]


def _price_answer(base, odds_by_selection, dropped=()):
    return {"data": {"trueOdds": 2.0,
                     "selectionsMapped": [s for s in base if s not in dropped],
                     "selectionsNotMapped": [{"selectionId": s, "removalReasonCode": "IncompatibleSelections"}
                                             for s in dropped],
                     "compatibleMarkets": [{"id": "m_ml", "selections": [
                         {"id": sid, "trueOdds": odds} for sid, odds in odds_by_selection.items()]}]}}


class FakeResponse:
    def __init__(self, status_code, body):
        self.status_code, self._body, self.text = status_code, body, str(body)

    def json(self):
        return self._body


class FakeSession:
    def __init__(self, price_answer, price_status=200):
        self.price_answer, self.price_status = price_answer, price_status
        self.price_calls = []

    def get(self, url, params=None, headers=None, timeout=None):
        if url == dk_book.DK_SGP_PRICE_URL:
            self.price_calls.append((params, headers))
            return FakeResponse(self.price_status, self.price_answer)
        if url == dk_book.DK_LEAGUE_URL:
            return FakeResponse(200, {"events": []})
        if "parlays/v1/sgp/events" in url:
            return FakeResponse(200, {"data": {"markets": MARKETS}})
        return FakeResponse(200, {})                     # sportsbook warm-up page


def _book(monkeypatch, session):
    monkeypatch.setattr(dk_book.cffi_requests, "Session", lambda impersonate: session)
    monkeypatch.setattr(dk_book, "MIN_SECONDS_BETWEEN_PRICE_CALLS", 0.0)
    return DraftKingsBook()


def test_cells_sharing_a_base_cost_one_request(monkeypatch):
    session = FakeSession(_price_answer(["q1_sea"], {"ml_sea": 2.4, "ml_lac": 15.0}))
    book = _book(monkeypatch, session)
    assert book.price(GAME, ("q1_sea", "ml_sea")) == 2.4
    assert book.price(GAME, ("q1_sea", "ml_lac")) == 15.0
    assert len(session.price_calls) == 1
    params, headers = session.price_calls[0]
    assert params["selections"] == "q1_sea" and params["marketCandidates"] == "m_ml"
    assert headers == {"X-SportId": dk_book.DK_FOOTBALL_SPORT_ID}


def test_a_base_dk_shrank_prices_nothing():
    answer = _price_answer(["q1_sea", "q1_lac"], {"ml_sea": 3.3}, dropped=("q1_lac",))
    assert parse_added_leg_odds(answer, frozenset({"q1_sea", "q1_lac"})) == {}


def test_a_candidate_dk_will_not_add_is_a_decline(monkeypatch):
    book = _book(monkeypatch, FakeSession(_price_answer(["q1_sea"], {"ml_sea": 2.4})))
    assert book.price(GAME, ("q1_sea", "ml_lac")) is None


def test_http_failure_raises(monkeypatch):
    book = _book(monkeypatch, FakeSession({"message": "nope"}, price_status=404))
    with pytest.raises(RuntimeError, match="HTTP 404"):
        book.price(GAME, ("q1_sea", "ml_sea"))
