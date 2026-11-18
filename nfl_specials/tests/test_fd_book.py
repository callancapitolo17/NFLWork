"""FanDuelBook: stopping on a block — no network."""
import pytest

from nfl_specials import fd_book
from nfl_specials.books import BookBlocked, BookGame, RequestPacer
from nfl_specials.fd_book import FanDuelBook

GAME = BookGame("e1", home="SEA", away="LAC", game_start_time="2099-01-01T00:00:00Z")
EVENT_PAGE = {"attachments": {"markets": {"m1": {
    "marketId": "m1", "marketName": "1st Half Total", "sgmMarket": True, "marketStatus": "OPEN",
    "runners": [{"selectionId": 11, "runnerName": "Over", "handicap": 23.5},
                {"selectionId": 12, "runnerName": "Under", "handicap": 23.5}]}}}}


class FakeResponse:
    def __init__(self, status_code, body):
        self.status_code, self._body, self.text = status_code, body, str(body)

    def json(self):
        return self._body


class FakeSession:
    def __init__(self, event_status=200, price_status=200):
        self.event_status, self.price_status = event_status, price_status
        self.sent = []

    def post(self, url, headers=None, data=None, timeout=None):
        self.sent.append(url)
        if url == fd_book.FD_SCAN_URL:
            return FakeResponse(200, {})
        return FakeResponse(self.price_status, {})

    def get(self, url, headers=None, timeout=None):
        self.sent.append(url)
        return FakeResponse(self.event_status, EVENT_PAGE)


def _book(monkeypatch, session):
    monkeypatch.setattr(fd_book.cffi_requests, "Session", lambda impersonate: session)
    return FanDuelBook(pacer=RequestPacer("FanDuel", 0.0))


def test_an_implybets_403_stops_fanduel_for_the_refresh(monkeypatch):
    session = FakeSession(price_status=403)
    book = _book(monkeypatch, session)
    with pytest.raises(RuntimeError, match="HTTP 403"):
        book.price(GAME, (("m1", 11), ("m2", 21)))
    sent = len(session.sent)
    with pytest.raises(BookBlocked, match="implyBets HTTP 403"):
        book.price(GAME, (("m1", 12), ("m2", 21)))
    with pytest.raises(BookBlocked):
        book.leg_market(BookGame("e2", "SEA", "LAC", GAME.game_start_time), "home", None)
    assert len(session.sent) == sent


def test_an_event_page_403_stops_fanduel_too(monkeypatch):
    session = FakeSession(event_status=403)
    book = _book(monkeypatch, session)
    with pytest.raises(RuntimeError, match="event page"):
        book.leg_market(GAME, "home", None)
    sent = len(session.sent)
    with pytest.raises(BookBlocked, match="event page HTTP 403"):
        book.price(GAME, (("m1", 11), ("m2", 21)))
    assert len(session.sent) == sent


def test_a_blocked_scan_stops_the_next_fanduel_open(monkeypatch):
    """The fecta board opens FanDuel first; the 1H board's two opens must not scan into the block."""
    class BlockedScan(FakeSession):
        def post(self, url, headers=None, data=None, timeout=None):
            self.sent.append(url)
            return FakeResponse(403, {})
    session = BlockedScan()
    monkeypatch.setattr(fd_book.cffi_requests, "Session", lambda impersonate: session)
    pacer = RequestPacer("FanDuel", 0.0)
    with pytest.raises(RuntimeError, match="event scan returned HTTP 403"):
        FanDuelBook(pacer=pacer)
    with pytest.raises(BookBlocked, match="event scan HTTP 403"):
        FanDuelBook(league="CFB", pacer=pacer)
    assert len(session.sent) == 1
