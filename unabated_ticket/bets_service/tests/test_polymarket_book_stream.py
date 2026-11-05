"""Polymarket US's live books over the markets websocket
(polymarket_book_stream.py), driven frame by frame with a fake connection and
a fake clock: on-demand subscribes with Place's first, whole-book frames as
the recon of 2026-11-05 saw them, renewals of a still market, refusals, idle
markets dropped, the background cap, pings and silence; plus POST
/polymarket_books.json over a real HTTP server."""
import json
import threading
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer

import pytest

from unabated_ticket.bets_service import polymarket_book_stream as pbs
from unabated_ticket.bets_service import service
from unabated_ticket.bets_service.store import BetsStore

DUKE = "tsc-cfb-duke-ncst-2026-11-07-total-48pt5"
OKST = "tsc-cfb-okst-kanst-2026-11-07-total-58pt5"


class Clock:
    def __init__(self):
        self.now = 1000.0

    def __call__(self) -> float:
        return self.now


class FakeClient:
    key_id = "0ad58243-test"

    def __init__(self):
        self.signed = []

    def signed_headers(self, method, path):
        self.signed.append((method, path))
        return {"X-PM-Access-Key": self.key_id, "X-PM-Timestamp": "1", "X-PM-Signature": "sig"}


class FakeConnection:
    def __init__(self):
        self.sent = []
        self.pings = 0
        self.incoming = []
        self.closed = False

    def send(self, text):
        self.sent.append(json.loads(text))

    def recv(self, timeout_sec):
        if not self.incoming:
            return None
        item = self.incoming.pop(0)
        if isinstance(item, Exception):
            raise item
        return item if isinstance(item, str) else json.dumps(item)

    def ping(self):
        self.pings += 1

    def close(self):
        self.closed = True


def market_data(slug=DUKE, bid="0.6200", offer="0.6300", state="MARKET_STATE_OPEN"):
    """The recon's whole-book frame body (trimmed to one level a side)."""
    return {"marketSlug": slug, "bids": [{"px": {"value": bid, "currency": "USD"}, "qty": "3515.1700"}],
            "offers": [{"px": {"value": offer, "currency": "USD"}, "qty": "3857.0000"}], "state": state,
            "stats": {"lastTradePx": {"value": "0.6200", "currency": "USD"}}, "transactTime": "2026-11-05T22:53:38Z"}


def book_frame(request_id, slug=DUKE, **fields):
    return {"requestId": request_id, "subscriptionType": "SUBSCRIPTION_TYPE_MARKET_DATA", "marketData": market_data(slug, **fields)}


def handle(stream, *messages):
    with stream._cond:
        for message in messages:
            stream._handle(message)


@pytest.fixture
def made():
    clock = Clock()
    client = FakeClient()
    connection = FakeConnection()
    stream = pbs.PolymarketBookStream(client, "wss://api.polymarket.us/v1/ws/markets",
                                      connect=lambda url, headers: connection, clock=clock)
    stream._open()
    return stream, connection, clock, client


def subscribes(frames):
    return [(frame["subscribe"]["requestId"], frame["subscribe"]["marketSlugs"]) for frame in frames if "subscribe" in frame]


def subscribe(stream, slug=DUKE, **book):
    """Ask for `slug`, send its group's subscribe, answer it with a whole book; returns the requestId."""
    stream.books_for([slug], 0.0)
    (request_id, _slugs), = subscribes(stream._outgoing_frames())
    handle(stream, book_frame(request_id, slug, **book))
    return request_id


def test_the_upgrade_is_signed_like_a_get_of_the_stream_path(made):
    _stream, _connection, _clock, client = made
    assert client.signed == [("GET", "/v1/ws/markets")]
    assert pbs.stream_url("https://api.polymarket.us") == "wss://api.polymarket.us/v1/ws/markets"


def test_a_new_market_is_subscribed_and_served_as_the_rest_books_body(made):
    stream, _connection, clock, _client = made
    first = stream.books_for([DUKE], 0.0)
    assert "waiting for" in first["books"][DUKE]["error"]
    frames = stream._outgoing_frames()
    assert frames == [{"subscribe": {"requestId": "pm-1", "subscriptionType": "SUBSCRIPTION_TYPE_MARKET_DATA",
                                     "marketSlugs": [DUKE]}}]
    handle(stream, book_frame("pm-1"))
    clock.now += 2
    view = stream.books_for([DUKE], 0.0)["books"][DUKE]
    assert view == {"marketData": market_data(), "ageMs": 2000}
    assert stream._outgoing_frames() == [], "subscribed once, not again every loop"


def test_markets_share_at_most_ten_subscriptions(made):
    """Polymarket US refuses an 11th subscription on a connection (live, 2026-11-05)."""
    stream, _connection, _clock, _client = made
    slugs = [f"slug-{index:03d}" for index in range(35)]
    stream.books_for(slugs, 0.0)
    sent = subscribes(stream._outgoing_frames())
    assert len(sent) == pbs.MAX_SUBSCRIPTIONS
    assert sorted(slug for _request_id, group in sent for slug in group) == slugs
    assert max(len(group) for _request_id, group in sent) - min(len(group) for _request_id, group in sent) <= 1


def test_a_new_market_resubscribes_only_the_group_it_joins(made):
    stream, _connection, _clock, _client = made
    stream.books_for([DUKE], 0.0)
    stream._outgoing_frames()
    handle(stream, book_frame("pm-1"))
    stream.books_for([OKST], 0.0)
    frames = stream._outgoing_frames()
    assert subscribes(frames) == [("pm-2", [OKST])], "the emptiest group takes it; DUKE's group is left alone"
    assert "marketData" in stream.books_for([DUKE], 0.0)["books"][DUKE]


def test_a_full_group_is_resubscribed_with_its_new_member(made, monkeypatch):
    stream, _connection, _clock, _client = made
    monkeypatch.setattr(pbs, "MAX_SUBSCRIPTIONS", 1)
    request_id = subscribe(stream)
    stream.books_for([OKST], 0.0)
    frames = stream._outgoing_frames()
    assert frames == [{"unsubscribe": {"requestId": request_id}},
                      {"subscribe": {"requestId": "pm-2", "subscriptionType": "SUBSCRIPTION_TYPE_MARKET_DATA",
                                     "marketSlugs": [DUKE, OKST]}}]
    assert "marketData" in stream.books_for([DUKE], 0.0)["books"][DUKE], "the held book is served meanwhile"
    handle(stream, book_frame(request_id, OKST))
    assert "waiting for" in stream.books_for([OKST], 0.0)["books"][OKST]["error"], "the old id's frames are ignored"
    handle(stream, book_frame("pm-2", OKST))
    assert "marketData" in stream.books_for([OKST], 0.0)["books"][OKST]


def test_when_every_group_is_full_a_new_market_says_so(made, monkeypatch):
    stream, _connection, _clock, _client = made
    monkeypatch.setattr(pbs, "MAX_SUBSCRIPTIONS", 1)
    monkeypatch.setattr(pbs, "SLUGS_PER_SUBSCRIPTION", 1)
    subscribe(stream)
    stream.books_for([OKST], 0.0)
    assert subscribes(stream._outgoing_frames()) == []
    assert "subscriptions of 1 markets are full" in stream.books_for([OKST], 0.0)["books"][OKST]["error"]


def test_every_new_whole_book_replaces_the_last(made):
    stream, _connection, _clock, _client = made
    request_id = subscribe(stream)
    handle(stream, book_frame(request_id, offer="0.6400"))
    assert stream.books_for([DUKE], 0.0)["books"][DUKE]["marketData"]["offers"][0]["px"]["value"] == "0.6400"


def test_place_is_placed_ahead_of_the_live_pass(made, monkeypatch):
    stream, _connection, _clock, _client = made
    monkeypatch.setattr(pbs, "MAX_SUBSCRIPTIONS", 1)
    monkeypatch.setattr(pbs, "SLUGS_PER_SUBSCRIPTION", 2)
    stream.books_for(["z-live-a", "a-live-b"], 0.0, background=True)
    stream.books_for([OKST], 0.0, priority=True)
    assert subscribes(stream._outgoing_frames()) == [("pm-1", sorted([OKST, "a-live-b"]))]


def test_an_unanswered_subscribe_is_sent_again_under_a_new_id_after_the_wait(made):
    stream, _connection, clock, _client = made
    stream.books_for([DUKE], 0.0)
    assert subscribes(stream._outgoing_frames()) == [("pm-1", [DUKE])]
    clock.now += pbs.ANSWER_WAIT_SEC + 1
    frames = stream._outgoing_frames()
    assert frames == [{"unsubscribe": {"requestId": "pm-1"}},
                      {"subscribe": {"requestId": "pm-2", "subscriptionType": "SUBSCRIPTION_TYPE_MARKET_DATA",
                                     "marketSlugs": [DUKE]}}]
    handle(stream, book_frame("pm-1"))
    assert "waiting for" in stream.books_for([DUKE], 0.0)["books"][DUKE]["error"], "the dropped id's late book is ignored"
    handle(stream, book_frame("pm-2"))
    assert "marketData" in stream.books_for([DUKE], 0.0)["books"][DUKE]


def test_a_market_its_answered_group_sent_no_book_for_reads_as_an_error(made, monkeypatch):
    stream, _connection, clock, _client = made
    monkeypatch.setattr(pbs, "MAX_SUBSCRIPTIONS", 1)
    stream.books_for([DUKE, OKST], 0.0)
    stream._outgoing_frames()
    handle(stream, book_frame("pm-1"))
    clock.now += pbs.ANSWER_WAIT_SEC + 1
    assert "no book of" in stream.books_for([OKST], 0.0)["books"][OKST]["error"]
    assert stream._outgoing_frames() == [], "an answered group is not resent for one silent market"


def test_a_group_is_renewed_for_fresh_books_and_its_last_books_served_meanwhile(made):
    stream, _connection, clock, _client = made
    subscribe(stream)
    clock.now += pbs.RESUBSCRIBE_SEC + 1
    stream.books_for([DUKE], 0.0)
    frames = stream._outgoing_frames()
    assert frames[0] == {"unsubscribe": {"requestId": "pm-1"}} and subscribes(frames) == [("pm-2", [DUKE])]
    assert "marketData" in stream.books_for([DUKE], 0.0)["books"][DUKE], "the renewal does not blank the book"
    clock.now += pbs.BOOK_MAX_AGE_SEC - pbs.RESUBSCRIBE_SEC
    assert "s old" in stream.books_for([DUKE], 0.0)["books"][DUKE]["error"], "an unanswered renewal never serves an old book"
    handle(stream, book_frame("pm-2"))
    assert stream.books_for([DUKE], 0.0)["books"][DUKE]["ageMs"] == 0


def test_a_refused_group_refuses_its_markets_and_asks_again_a_minute_later(made):
    stream, _connection, clock, _client = made
    stream.books_for([DUKE], 0.0)
    stream._outgoing_frames()
    handle(stream, {"requestId": "pm-1", "error": "max subscriptions per connection reached",
                    "subscriptionType": "SUBSCRIPTION_TYPE_MARKET_DATA"})
    assert "max subscriptions" in stream.books_for([DUKE], 0.0)["books"][DUKE]["error"]
    assert stream._outgoing_frames() == []
    clock.now += pbs.REFUSED_RETRY_SEC + 1
    assert subscribes(stream._outgoing_frames()) == [("pm-2", [DUKE])]


def test_a_book_frame_naming_a_market_its_subscription_does_not_name_is_not_served(made):
    stream, _connection, _clock, _client = made
    stream.books_for([DUKE], 0.0)
    stream._outgoing_frames()
    handle(stream, book_frame("pm-1", slug=OKST))
    stream.books_for([OKST], 0.0)
    assert "waiting for" in stream.books_for([OKST], 0.0)["books"][OKST]["error"]


def test_frames_without_a_whole_book_change_nothing(made):
    stream, _connection, _clock, _client = made
    request_id = subscribe(stream)
    handle(stream, {"heartbeat": {}}, {"requestId": request_id, "marketDataLite": {"marketSlug": DUKE}},
           {"requestId": "pm-999", "marketData": market_data(offer="0.9900")})
    assert stream.books_for([DUKE], 0.0)["books"][DUKE]["marketData"] == market_data()


def test_a_market_nobody_asks_about_leaves_its_group(made):
    stream, _connection, clock, _client = made
    subscribe(stream)
    clock.now += pbs.WATCH_IDLE_SEC + 1
    assert stream._outgoing_frames() == [{"unsubscribe": {"requestId": "pm-1"}}]
    status = stream.status()
    assert (status["wantedMarkets"], status["booksHeld"], status["subscriptions"]) == (0, 0, 0)


def test_the_live_pass_is_capped_but_place_never_is(made, monkeypatch):
    stream, _connection, _clock, _client = made
    monkeypatch.setattr(pbs, "MAX_BACKGROUND_MARKETS", 1)
    stream.books_for(["live-a"], 0.0, background=True)
    capped = stream.books_for(["live-b"], 0.0, background=True)["books"]["live-b"]
    assert "already watches 1" in capped["error"]
    assert "waiting for" in stream.books_for(["live-b"], 0.0, priority=True)["books"]["live-b"]["error"]


def test_no_book_is_served_while_disconnected_and_a_reconnect_starts_empty(made):
    stream, _connection, _clock, _client = made
    subscribe(stream)
    stream._drop("ConnectionError: gone")
    assert "not connected (ConnectionError: gone)" in stream.books_for([DUKE], 0.0)["books"][DUKE]["error"]
    stream._open()
    assert "waiting for" in stream.books_for([DUKE], 0.0)["books"][DUKE]["error"]
    assert subscribes(stream._outgoing_frames()) == [("pm-2", [DUKE])], "the wanted markets are subscribed again"


def test_the_pump_pings_and_drops_a_line_that_answers_nothing(made):
    stream, connection, clock, _client = made
    stop = threading.Event()

    def recv(timeout_sec):
        clock.now += 5
        return None

    connection.recv = recv
    with pytest.raises(ConnectionError, match="no frame or pong"):
        stream._pump(stop)
    assert connection.pings == pbs.SILENT_RECONNECT_SEC // pbs.PING_SEC


def test_a_pong_keeps_the_line_alive(made):
    stream, connection, clock, _client = made
    stop = threading.Event()
    calls = []

    def recv(timeout_sec):
        clock.now += 5
        calls.append(1)
        if len(calls) > 30:
            stop.set()
        return ""

    connection.recv = recv
    stream._pump(stop)
    assert stream.status()["connected"] is True


def test_a_waiting_ask_returns_as_soon_as_the_book_lands(made):
    stream, _connection, _clock, _client = made
    stream.books_for([DUKE], 0.0)
    stream._outgoing_frames()
    timer = threading.Timer(0.05, lambda: handle_and_notify(stream, book_frame("pm-1")))
    timer.start()
    # A real clock for the wait: the fake one never moves, so a missed notify would hang to the cap.
    stream._clock = __import__("time").monotonic
    answer = stream.books_for([DUKE], 2.0)
    timer.join()
    assert "marketData" in answer["books"][DUKE]


def handle_and_notify(stream, message):
    with stream._cond:
        stream._handle(message)
        stream._cond.notify_all()


@pytest.fixture
def bets_store(tmp_path):
    store = BetsStore(tmp_path / "bets.duckdb", 7)
    yield store
    store.close()


def _post(url, payload, host=None):
    headers = {"Content-Type": "application/json"}
    if host:
        headers["Host"] = host
    request = urllib.request.Request(url, data=json.dumps(payload).encode(), method="POST", headers=headers)
    try:
        with urllib.request.urlopen(request, timeout=5) as response:
            return response.status, json.loads(response.read())
    except urllib.error.HTTPError as error:
        return error.code, json.loads(error.read())


def _serve(handler):
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    return server, f"http://127.0.0.1:{server.server_address[1]}"


def test_post_polymarket_books_answers_the_streams_books(bets_store, made):
    stream, _connection, _clock, _client = made
    subscribe(stream)
    server, base = _serve(service.make_handler(bets_store, 0.0, ["kalshi"], polymarket_book_stream=stream,
                                               extra_hosts=("mlb-stack.test",)))
    try:
        status, body = _post(f"{base}/polymarket_books.json", {"slugs": [DUKE, DUKE], "waitMs": 0, "priority": True})
        bad = _post(f"{base}/polymarket_books.json", {"slugs": [""]})
        foreign = _post(f"{base}/polymarket_books.json", {"slugs": []}, host="mlb-stack.test")
        with urllib.request.urlopen(f"{base}/health", timeout=5) as response:
            health = json.loads(response.read())
    finally:
        server.shutdown()
        server.server_close()
    assert status == 200 and body["books"] == {DUKE: {"marketData": market_data(), "ageMs": 0}}
    assert body["stream"]["connected"] is True
    assert bad[0] == 400 and "slugs[0]" in bad[1]["error"]
    assert foreign[0] == 403, "the runner's book route is loopback-only"
    assert health["polymarketStream"]["booksHeld"] == 1


def test_post_polymarket_books_is_503_without_the_key(bets_store):
    server, base = _serve(service.make_handler(bets_store, 0.0, ["kalshi"]))
    try:
        status, body = _post(f"{base}/polymarket_books.json", {"slugs": [], "waitMs": 0})
    finally:
        server.shutdown()
        server.server_close()
    assert status == 503 and "POLYMARKET_US_KEY_ID" in body["error"]


def test_the_stream_is_off_without_the_key(monkeypatch):
    monkeypatch.setattr(service.config, "POLYMARKET_US_KEY_ID", "")
    assert service.open_polymarket_book_stream() is None
