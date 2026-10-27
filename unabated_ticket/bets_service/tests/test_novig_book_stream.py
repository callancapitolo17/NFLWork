"""Novig's live books over the websocket (novig_book_stream.py), driven frame
by frame with a fake connection and a fake clock: on-demand event subscribes
paced by the stream throttle, snapshots and deltas, gaps healed by a snapshot,
markets opening under a watched event, go-live, refusals, idle events dropped,
silence and reconnects; plus POST /novig_books.json over a real HTTP server."""
import json
import threading
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer

import pytest

from unabated_ticket.bets_service import novig_book_stream as nbs
from unabated_ticket.bets_service import service
from unabated_ticket.bets_service.store import BetsStore

EVENT = "ev-1"
MARKET = "mk-1"
DUQ = "out-duq"
WSU = "out-wsu"


class Clock:
    def __init__(self):
        self.now = 1000.0

    def __call__(self) -> float:
        return self.now


class FakeClient:
    key_id = "4365ed85-test"

    def __init__(self):
        self.signed = []

    def signed_headers(self, method, path, raw_query="", data=b""):
        self.signed.append((method, path, raw_query, data))
        return {"Novig-Key-Id": self.key_id, "Novig-Timestamp": "1", "Novig-Signature": "sig"}


class FakeConnection:
    def __init__(self, incoming=()):
        self.sent = []
        self.incoming = list(incoming)
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

    def close(self):
        self.closed = True


def snapshot_reply(nonce, orders=None, status="OPEN", seq=48120, market=MARKET, event=EVENT):
    return {"nonce": nonce, "subscribed": {"markets": {}, "events": {event: "book"}, "private": []},
            "snapshot": {market: {"eventId": event,
                                  "book": {"seq": seq, "orders": orders if orders is not None else {
                                      DUQ: [{"order": "o1", "price": "0.445", "qty": 110}]}},
                                  "lifecycle": {"seq": 3, "status": status}}}}


def book_delta(seq, deltas, market=MARKET, event=EVENT, lifecycle=None):
    body = {"eventId": event, "book": {"seq": seq, "deltas": deltas}}
    if lifecycle is not None:
        body["lifecycle"] = lifecycle
    return {"ts": 1, "delta": {market: body}}


@pytest.fixture
def made():
    clock = Clock()
    client = FakeClient()
    connection = FakeConnection()
    stream = nbs.NovigBookStream(client, "wss://novig.test/v3/ws", connect=lambda url, headers: connection, clock=clock)
    stream._open()
    return stream, connection, clock, client


def ask(stream, wait_sec=0.0, market=MARKET, event=EVENT):
    return stream.books_for([{"marketId": market, "eventId": event}], wait_sec)["books"][market]


def subscribed(stream):
    """Send what is due and answer the first subscribe with the default snapshot."""
    frames = stream._outgoing_frames()
    with stream._cond:
        stream._handle(snapshot_reply(frames[0]["nonce"]))
    return frames


def test_the_upgrade_is_signed_as_a_get_of_the_stream_path(made):
    _stream, _connection, _clock, client = made
    assert client.signed == [("GET", "/v3/ws", "", b"")]
    assert nbs.stream_url("https://api.novig.com") == "wss://api.novig.com/v3/ws"


def test_an_event_is_subscribed_on_first_ask_and_its_snapshot_is_the_book(made):
    stream, _connection, _clock, _client = made
    assert ask(stream) == {"error": f"waiting for Novig's stream snapshot of event {EVENT}"}
    frames = stream._outgoing_frames()
    assert frames == [{"nonce": 1, "subscribe": {"events": {EVENT: "book"}}}]
    assert stream._outgoing_frames() == [], "an event in flight is not subscribed twice"
    with stream._cond:
        stream._handle(snapshot_reply(1))
    assert ask(stream) == {"eventId": EVENT, "status": "OPEN", "started": False, "seq": 48120,
                           "orders": {DUQ: [{"order": "o1", "price": "0.445", "qty": 110}]}}
    assert ask(stream, market="mk-other") == {"error": f"Novig's stream lists no market mk-other under event {EVENT}"}


def test_deltas_add_and_remove_orders_in_seq_order_and_a_stale_seq_is_ignored(made):
    stream, _connection, _clock, _client = made
    ask(stream)
    subscribed(stream)
    with stream._cond:
        stream._handle(book_delta(48121, [{"kind": "add", "order": "o2", "outcome": WSU, "price": "0.550", "qty": 40},
                                          {"kind": "remove", "order": "o1", "reason": "fill"}]))
        stream._handle(book_delta(48121, [{"kind": "remove", "order": "o2", "reason": "cancel"}]))
    view = ask(stream)
    assert view["seq"] == 48121
    assert view["orders"] == {WSU: [{"order": "o2", "price": "0.550", "qty": 40}]}


def test_a_gap_serves_no_book_until_a_fresh_snapshot_heals_it(made):
    stream, _connection, _clock, _client = made
    ask(stream)
    subscribed(stream)
    with stream._cond:
        stream._handle(book_delta(48123, [{"kind": "remove", "order": "o1", "reason": "fill"}]))
    assert "resyncing after a gap" in ask(stream)["error"]
    frames = stream._outgoing_frames()
    assert frames == [{"nonce": 2, "snapshot": {"markets": {MARKET: "book"}}}]
    assert stream._outgoing_frames() == [], "one healing snapshot in flight at a time"
    with stream._cond:
        stream._handle({"nonce": 2, **{k: v for k, v in snapshot_reply(2, orders={}, seq=48125).items() if k != "subscribed"}})
        stream._handle(book_delta(48124, [{"kind": "add", "order": "late", "outcome": DUQ, "price": "0.4", "qty": 1}]))
    view = ask(stream)
    assert (view["seq"], view["orders"]) == (48125, {}), "a delta the snapshot covers is dropped"


def test_a_market_opening_under_a_watched_event_starts_at_seq_zero(made):
    stream, _connection, _clock, _client = made
    ask(stream)
    subscribed(stream)
    with stream._cond:
        stream._handle(book_delta(1, [{"kind": "add", "order": "n1", "outcome": DUQ, "price": "0.300", "qty": 5}],
                                  market="mk-new", lifecycle={"seq": 1, "deltas": ["OPEN"]}))
        stream._handle(book_delta(4, [], market="mk-mid"))
    assert ask(stream, market="mk-new")["orders"] == {DUQ: [{"order": "n1", "price": "0.300", "qty": 5}]}
    assert ask(stream, market="mk-new")["status"] == "OPEN"
    assert "resyncing" in ask(stream, market="mk-mid")["error"], "a first delta past seq 1 is a gap"
    with stream._cond:
        stream._handle(book_delta(1, [], market="mk-unwatched", event="ev-unwatched"))
    assert "mk-unwatched" not in stream._books


def test_go_live_marks_the_market_started_and_a_close_its_status(made):
    stream, _connection, _clock, _client = made
    ask(stream)
    subscribed(stream)
    with stream._cond:
        stream._handle(book_delta(48121, [{"kind": "remove", "order": "o1", "reason": "cancel"}],
                                  lifecycle={"seq": 4, "deltas": ["GOLIVE"]}))
    assert ask(stream)["started"] is True
    with stream._cond:
        stream._handle({"delta": {MARKET: {"eventId": EVENT, "lifecycle": {"seq": 5, "deltas": ["CLOSE"]}}}})
    assert ask(stream)["status"] == "CLOSED"


def test_subscribes_are_paced_by_the_stream_throttle(made):
    stream, _connection, clock, _client = made
    events = [f"ev-{index:02d}" for index in range(40)]
    stream.books_for([{"marketId": f"mk-{event}", "eventId": event} for event in events], 0.0)
    first = stream._outgoing_frames()
    # 512 tokens less the 32 the upgrade cost = 30 events at 16 each.
    assert [len(frame["subscribe"]["events"]) for frame in first] == [30]
    assert stream._outgoing_frames() == []
    clock.now += 8  # 32 tokens: two more events
    assert [len(frame["subscribe"]["events"]) for frame in stream._outgoing_frames()] == [2]


def test_a_refused_subscribe_answers_its_code_and_is_retried_after_a_minute(made):
    stream, _connection, clock, _client = made
    ask(stream)
    frames = stream._outgoing_frames()
    with stream._cond:
        stream._handle({"nonce": frames[0]["nonce"], "code": "SUBSCRIPTION_LIMIT_EXCEEDED", "message": "over 2048"})
    assert ask(stream) == {"error": f"Novig refused the stream subscription to event {EVENT}: SUBSCRIPTION_LIMIT_EXCEEDED over 2048"}
    assert stream._outgoing_frames() == []
    clock.now += nbs.REFUSED_RETRY_SEC + 1
    ask(stream)
    assert [frame["subscribe"] for frame in stream._outgoing_frames()] == [{"events": {EVENT: "book"}}]


def test_a_throttled_frame_empties_the_bucket_and_its_subscribe_goes_back_in_the_queue(made):
    stream, _connection, clock, _client = made
    ask(stream)
    stream._outgoing_frames()
    with stream._cond:
        stream._handle({"code": "RATE_LIMIT_EXCEEDED", "message": "slow down"})
    assert stream._bucket.available() == 0
    clock.now += nbs.SUBSCRIBE_ANSWER_SEC + 1
    ask(stream)
    assert [frame["subscribe"] for frame in stream._outgoing_frames()] == [{"events": {EVENT: "book"}}]


def test_an_event_not_asked_about_for_ten_minutes_is_unsubscribed_and_its_books_dropped(made):
    stream, _connection, clock, _client = made
    ask(stream)
    subscribed(stream)
    clock.now += nbs.WATCH_IDLE_SEC + 1
    assert stream._outgoing_frames() == [{"nonce": 2, "unsubscribe": [f"event:{EVENT}"]}]
    assert stream._books == {} and EVENT not in stream._wanted


def test_a_dropped_connection_serves_no_book_and_reconnects_with_a_fresh_subscribe(made):
    stream, connection, _clock, _client = made
    ask(stream)
    subscribed(stream)
    stream._drop("ConnectionError: Novig closed the stream (code 1008, 'SLOW_CONSUMER')")
    assert connection.closed
    assert ask(stream) == {"error": "the Novig stream is not connected (ConnectionError: Novig closed the stream "
                                    "(code 1008, 'SLOW_CONSUMER'))"}
    stream._open()
    assert ask(stream)["error"].startswith("waiting for"), "the wanted event survives the reconnect"
    assert [frame["subscribe"] for frame in stream._outgoing_frames()] == [{"events": {EVENT: "book"}}]


def test_forty_seconds_of_silence_is_a_dead_connection(made):
    stream, _connection, clock, _client = made
    clock.now += nbs.SILENT_RECONNECT_SEC + 1
    with pytest.raises(ConnectionError, match="no frame from Novig"):
        stream._pump(threading.Event())


def test_run_reconnects_after_a_failed_connect_and_stops_on_request():
    clock = Clock()
    stop = threading.Event()
    connection = FakeConnection([ConnectionError("Novig closed the stream (code 1008, 'GEOLOCATION_NOT_FOUND')")])
    attempts = []

    def connect(url, headers):
        attempts.append(url)
        if len(attempts) == 1:
            raise OSError("connection refused")
        if len(attempts) == 3:
            stop.set()
        return connection

    stream = nbs.NovigBookStream(FakeClient(), "wss://novig.test/v3/ws", connect=connect, clock=clock)
    original_wait = stop.wait
    stop.wait = lambda timeout=None: original_wait(0)
    stream.run(stop)
    assert len(attempts) == 3
    assert stream.status()["connected"] is False


def test_books_for_waits_for_a_snapshot_arriving_on_the_stream_thread():
    real = nbs.NovigBookStream(FakeClient(), "wss://novig.test/v3/ws", connect=lambda url, headers: FakeConnection())
    real._open()
    results = []
    asker = threading.Thread(target=lambda: results.append(ask(real, wait_sec=2.0)))
    asker.start()
    while not real._wanted:
        pass
    frames = real._outgoing_frames()
    with real._cond:
        real._handle(snapshot_reply(frames[0]["nonce"]))
        real._cond.notify_all()
    asker.join(timeout=3)
    assert results and results[0]["seq"] == 48120


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


def test_post_novig_books_answers_the_streams_books(bets_store, made):
    stream, _connection, _clock, _client = made
    ask(stream)
    subscribed(stream)
    server, base = _serve(service.make_handler(bets_store, 0.0, ["kalshi"], novig_book_stream=stream,
                                               extra_hosts=("mlb-stack.test",)))
    try:
        status, body = _post(f"{base}/novig_books.json", {"markets": [{"marketId": MARKET, "eventId": EVENT}], "waitMs": 0})
        bad = _post(f"{base}/novig_books.json", {"markets": [{"marketId": MARKET}]})
        foreign = _post(f"{base}/novig_books.json", {"markets": []}, host="mlb-stack.test")
    finally:
        server.shutdown()
        server.server_close()
    assert status == 200 and body["books"][MARKET]["seq"] == 48120 and body["stream"]["connected"] is True
    assert bad[0] == 400 and "markets[0]" in bad[1]["error"]
    assert foreign[0] == 403, "the runner's book route is loopback-only"


def test_post_novig_books_is_503_without_the_bot_key(bets_store):
    server, base = _serve(service.make_handler(bets_store, 0.0, ["kalshi"]))
    try:
        status, body = _post(f"{base}/novig_books.json", {"markets": [], "waitMs": 0})
    finally:
        server.shutdown()
        server.server_close()
    assert status == 503 and "NOVIG_TAKER_KEY_ID" in body["error"]
