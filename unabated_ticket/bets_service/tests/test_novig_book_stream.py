"""Novig's live books over the websocket (novig_book_stream.py), driven frame
by frame with a fake connection and a fake clock: on-demand market subscribes
paced by the stream throttle with Place's first, snapshots and deltas, gaps
buffered and replayed on a snapshot, go-live, refusals, idle markets dropped,
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
O1 = {"order": "o1", "price": "0.445", "qty": 110}


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


def market_snapshot(orders=None, seq=48120, status="OPEN", lifecycle_seq=3):
    return {"eventId": EVENT, "book": {"seq": seq, "orders": {DUQ: [O1]} if orders is None else orders},
            "lifecycle": {"seq": lifecycle_seq, "status": status}}


def answer(nonce, markets=None):
    """The answer to a subscribe or a snapshot: `markets` {id: snapshot body}, MARKET's default one if None."""
    return {"nonce": nonce, "snapshot": {MARKET: market_snapshot()} if markets is None else markets}


def book_delta(seq, deltas, market=MARKET, lifecycle=None):
    body = {"eventId": EVENT, "book": {"seq": seq, "deltas": deltas}}
    if lifecycle is not None:
        body["lifecycle"] = lifecycle
    return {"ts": 1, "delta": {market: body}}


def add(order, outcome=DUQ, price="0.400", qty=5):
    return {"kind": "add", "order": order, "outcome": outcome, "price": price, "qty": qty}


def handle(stream, *messages):
    with stream._cond:
        for message in messages:
            stream._handle(message)


@pytest.fixture
def made():
    clock = Clock()
    client = FakeClient()
    connection = FakeConnection()
    stream = nbs.NovigBookStream(client, "wss://novig.test/v3/ws", connect=lambda url, headers: connection, clock=clock)
    stream._open()
    return stream, connection, clock, client


def ask(stream, market=MARKET, wait_sec=0.0, priority=False):
    return stream.books_for([market], wait_sec, priority)["books"][market]


def subscribe(stream):
    """Ask for MARKET, send what is due and answer it with the default snapshot."""
    ask(stream)
    frames = stream._outgoing_frames()
    handle(stream, answer(frames[0]["nonce"]))
    return frames


def test_the_upgrade_is_signed_as_a_get_of_the_stream_path(made):
    _stream, _connection, _clock, client = made
    assert client.signed == [("GET", "/v3/ws", "", b"")]
    assert nbs.stream_url("https://api.novig.com") == "wss://api.novig.com/v3/ws"
    assert nbs.close_reason(b"\x03\xf0SLOW_CONSUMER") == "code 1008, 'SLOW_CONSUMER'"


def test_a_market_is_subscribed_on_first_ask_and_its_snapshot_is_the_book(made):
    stream, _connection, _clock, _client = made
    assert ask(stream) == {"error": f"waiting for Novig's stream snapshot of market {MARKET}"}
    assert stream._outgoing_frames() == [{"nonce": 1, "subscribe": {"markets": {MARKET: "book"}}}]
    assert stream._outgoing_frames() == [], "a market in flight is not subscribed twice"
    handle(stream, answer(1))
    assert ask(stream) == {"status": "OPEN", "started": False, "seq": 48120, "orders": {DUQ: [O1]}}
    assert stream.books_for([], 1.0) == {"stream": stream.status(), "books": {}}


def test_a_market_the_subscribe_answer_leaves_out_is_refused_and_asked_again_a_minute_later(made):
    stream, _connection, clock, _client = made
    ask(stream, "mk-gone")
    handle(stream, answer(stream._outgoing_frames()[0]["nonce"], markets={}))
    assert "no such market" in ask(stream, "mk-gone")["error"]
    assert stream._outgoing_frames() == []
    clock.now += nbs.REFUSED_RETRY_SEC + 1
    ask(stream, "mk-gone")
    assert [frame["subscribe"] for frame in stream._outgoing_frames()] == [{"markets": {"mk-gone": "book"}}]


def test_deltas_add_and_remove_orders_in_seq_order_and_a_stale_seq_is_ignored(made):
    stream, _connection, _clock, _client = made
    subscribe(stream)
    handle(stream, book_delta(48121, [add("o2", WSU, "0.550", 40), {"kind": "remove", "order": "o1", "reason": "fill"}]),
           book_delta(48121, [{"kind": "remove", "order": "o2", "reason": "cancel"}]))
    assert ask(stream) == {"status": "OPEN", "started": False, "seq": 48121,
                           "orders": {WSU: [{"order": "o2", "price": "0.550", "qty": 40}]}}


def test_a_partial_fill_keeps_the_remainder(made):
    stream, _connection, _clock, _client = made
    subscribe(stream)
    handle(stream, book_delta(48121, [{"kind": "remove", "order": "o1", "reason": "fill"}, add("o1", DUQ, "0.445", 60)]))
    assert ask(stream)["orders"] == {DUQ: [{"order": "o1", "price": "0.445", "qty": 60}]}


def test_a_gap_buffers_the_deltas_after_it_and_replays_those_past_the_snapshot(made):
    stream, _connection, _clock, _client = made
    subscribe(stream)
    handle(stream, book_delta(48122, [add("a")]), book_delta(48123, [add("b")]))
    assert "resyncing after a gap" in ask(stream)["error"]
    frames = stream._outgoing_frames()
    assert frames == [{"nonce": 2, "snapshot": {"markets": {MARKET: "book"}}}]
    assert stream._outgoing_frames() == [], "one healing snapshot in flight at a time"
    handle(stream, book_delta(48124, [add("c")]))
    # Taken at 48122 (it holds a): b and c replay on it.
    handle(stream, answer(2, {MARKET: market_snapshot(orders={DUQ: [{"order": "a", "price": "0.400", "qty": 5}]}, seq=48122)}))
    view = ask(stream)
    assert view["seq"] == 48124
    assert sorted(order["order"] for order in view["orders"][DUQ]) == ["a", "b", "c"]


def test_a_snapshot_short_of_the_buffered_deltas_is_a_new_gap(made):
    stream, _connection, _clock, _client = made
    subscribe(stream)
    handle(stream, book_delta(48125, [add("late")]))
    handle(stream, answer(stream._outgoing_frames()[0]["nonce"], {MARKET: market_snapshot(seq=48122)}))
    assert "resyncing" in ask(stream)["error"], "48123-48124 never arrived: not served"


def test_a_refused_gap_snapshot_waits_a_minute_and_never_starves_new_subscribes(made):
    stream, _connection, clock, _client = made
    subscribe(stream)
    handle(stream, book_delta(48125, []))
    nonce = stream._outgoing_frames()[0]["nonce"]
    handle(stream, {"nonce": nonce, "code": "MARKET_NOT_FOUND", "message": "gone"})
    assert ask(stream) == {"error": f"Novig refused the stream subscription to market {MARKET}: MARKET_NOT_FOUND gone"}
    ask(stream, "mk-2")
    assert [list(frame) for frame in stream._outgoing_frames()] == [["nonce", "subscribe"]]
    clock.now += nbs.REFUSED_RETRY_SEC + 1
    assert [frame["snapshot"] for frame in stream._outgoing_frames() if "snapshot" in frame] == [{"markets": {MARKET: "book"}}]


def test_go_live_marks_the_market_started_and_a_close_its_status(made):
    stream, _connection, _clock, _client = made
    subscribe(stream)
    handle(stream, book_delta(48121, [{"kind": "remove", "order": "o1", "reason": "cancel"}], lifecycle={"seq": 4, "deltas": ["GOLIVE"]}))
    assert ask(stream)["started"] is True
    handle(stream, {"delta": {MARKET: {"eventId": EVENT, "lifecycle": {"seq": 5, "deltas": ["CLOSE"]}}}},
           {"delta": {MARKET: {"eventId": EVENT, "lifecycle": {"seq": 4, "deltas": ["OPEN"]}}}})
    assert ask(stream)["status"] == "CLOSED", "a late, older status change is ignored"


def test_a_lifecycle_gap_resyncs_the_market(made):
    stream, _connection, _clock, _client = made
    subscribe(stream)
    handle(stream, {"delta": {MARKET: {"eventId": EVENT, "lifecycle": {"seq": 6, "deltas": ["CLOSE"]}}}})
    assert "resyncing" in ask(stream)["error"]
    assert [frame.get("snapshot") for frame in stream._outgoing_frames()] == [{"markets": {MARKET: "book"}}]


def test_an_unreadable_part_of_a_frame_resyncs_only_its_market(made):
    stream, _connection, _clock, _client = made
    subscribe(stream)
    handle(stream, book_delta(48121, [{"kind": "replace", "order": "o1"}]))
    assert "resyncing" in ask(stream)["error"]
    with pytest.raises(TypeError):
        handle(stream, ["not", "an", "object"])


def test_deltas_of_a_market_not_subscribed_build_nothing(made):
    stream, _connection, _clock, _client = made
    subscribe(stream)
    handle(stream, book_delta(1, [add("x")], market="mk-other"),
           {"nonce": 99, "snapshot": {"mk-other": market_snapshot()}})
    assert "mk-other" not in stream._books


def test_subscribes_are_paced_by_the_stream_throttle_with_place_first(made):
    stream, _connection, clock, _client = made
    stream.books_for([f"mk-{index:02d}" for index in range(40)], 0.0)
    stream.books_for(["mk-place"], 0.0, priority=True)
    first = stream._outgoing_frames()
    # 512 tokens less the 32 the upgrade cost = 30 markets at 16 each, Place's first.
    assert [len(frame["subscribe"]["markets"]) for frame in first] == [30]
    assert next(iter(first[0]["subscribe"]["markets"])) == "mk-place"
    assert stream._outgoing_frames() == []
    clock.now += 8  # 32 tokens: two more markets
    assert [len(frame["subscribe"]["markets"]) for frame in stream._outgoing_frames()] == [2]


def test_the_live_price_pass_subscribes_behind_the_bot_until_the_bot_asks_too(made):
    stream, _connection, _clock, _client = made
    stream.books_for(["mk-aa-live", "mk-ab-both"], 0.0, background=True)
    stream.books_for(["mk-zz-bot"], 0.0)
    stream.books_for(["mk-ab-both"], 0.0)  # the bot asks for a market the live pass watches: it moves up
    stream.books_for(["mk-ab-both"], 0.0, background=True)  # a later background ask does not move it back
    first = stream._outgoing_frames()
    assert list(first[0]["subscribe"]["markets"]) == ["mk-ab-both", "mk-zz-bot", "mk-aa-live"]


def test_the_live_price_pass_watches_at_most_its_cap_and_the_bot_still_gets_in(made, monkeypatch):
    stream, _connection, _clock, _client = made
    monkeypatch.setattr(nbs, "MAX_BACKGROUND_MARKETS", 2)
    answer = stream.books_for(["mk-a", "mk-b", "mk-c"], 0.0, background=True)
    assert answer["books"]["mk-c"] == {"error": "the live price pass already watches 2 Novig markets"}
    assert "mk-c" not in stream._wanted
    stream.books_for(["mk-bot"], 0.0)
    assert "mk-bot" in stream._wanted and "mk-bot" not in stream._background


def test_a_refused_subscribe_answers_its_code(made):
    stream, _connection, _clock, _client = made
    ask(stream)
    nonce = stream._outgoing_frames()[0]["nonce"]
    handle(stream, {"nonce": nonce, "code": "SUBSCRIPTION_LIMIT_EXCEEDED", "message": "over 2048"})
    assert ask(stream) == {"error": f"Novig refused the stream subscription to market {MARKET}: SUBSCRIPTION_LIMIT_EXCEEDED over 2048"}


def test_a_throttled_frame_empties_the_bucket_and_its_subscribe_is_sent_again_later(made):
    stream, _connection, clock, _client = made
    ask(stream)
    stream._outgoing_frames()
    handle(stream, {"code": "RATE_LIMIT_EXCEEDED", "message": "slow down"})
    assert stream._bucket.available() == 0
    clock.now += nbs.ANSWER_WAIT_SEC + 1
    assert [frame["subscribe"] for frame in stream._outgoing_frames()] == [{"markets": {MARKET: "book"}}]


def test_a_market_not_asked_about_for_ten_minutes_is_unsubscribed_and_its_book_dropped(made):
    stream, _connection, clock, _client = made
    subscribe(stream)
    clock.now += nbs.WATCH_IDLE_SEC + 1
    handle(stream, {"heartbeat": {}})  # the line is alive
    assert stream._outgoing_frames() == [{"nonce": 2, "unsubscribe": [f"market:{MARKET}"]}]
    assert stream._books == {} and MARKET not in stream._wanted


def test_a_subscribe_answered_after_its_wait_is_adopted_while_wanted_else_unsubscribed(made):
    stream, _connection, clock, _client = made
    ask(stream)
    nonce = stream._outgoing_frames()[0]["nonce"]
    clock.now += nbs.ANSWER_WAIT_SEC + 1
    ask(stream)
    resent = stream._outgoing_frames()[0]["nonce"]
    handle(stream, answer(nonce))
    assert ask(stream)["seq"] == 48120, "the late answer is the book"
    # Idle with an answer still to come: forgotten, and the answer that lands is unsubscribed.
    ask(stream, "mk-late")
    late = stream._outgoing_frames()[0]["nonce"]
    clock.now += nbs.WATCH_IDLE_SEC + 1
    stream._outgoing_frames()
    assert "mk-late" not in stream._wanted
    handle(stream, answer(late, {"mk-late": market_snapshot()}), answer(resent))
    assert stream._outgoing_frames()[-1]["unsubscribe"] == [f"market:{MARKET}", "market:mk-late"], "both went idle"


def test_a_silent_line_serves_no_book_and_forty_seconds_drops_it(made):
    stream, _connection, clock, _client = made
    subscribe(stream)
    clock.now += nbs.STALE_VIEW_SEC + 1
    assert ask(stream)["error"].startswith("the Novig stream has sent nothing for")
    clock.now += nbs.SILENT_RECONNECT_SEC
    with pytest.raises(ConnectionError, match="no frame from Novig"):
        stream._pump(threading.Event())


def test_a_dropped_connection_serves_no_book_and_reconnects_with_a_fresh_subscribe(made):
    stream, connection, _clock, _client = made
    subscribe(stream)
    stream._drop("ConnectionError: Novig closed the stream (code 1008, 'SLOW_CONSUMER')")
    assert connection.closed
    assert ask(stream) == {"error": "the Novig stream is not connected (ConnectionError: Novig closed the stream "
                                    "(code 1008, 'SLOW_CONSUMER'))"}
    stream._open()
    assert ask(stream)["error"].startswith("waiting for"), "the wanted market survives the reconnect"
    assert [frame["subscribe"] for frame in stream._outgoing_frames()] == [{"markets": {MARKET: "book"}}]


def test_run_reconnects_after_failures_and_stops_on_request():
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

    stream = nbs.NovigBookStream(FakeClient(), "wss://novig.test/v3/ws", connect=connect, clock=Clock())
    original_wait = stop.wait
    stop.wait = lambda timeout=None: original_wait(0)
    stream.run(stop)
    assert len(attempts) == 3
    assert stream.status()["connected"] is False


def test_books_for_waits_for_a_snapshot_arriving_on_the_stream_thread():
    stream = nbs.NovigBookStream(FakeClient(), "wss://novig.test/v3/ws", connect=lambda url, headers: FakeConnection())
    stream._open()
    results = []
    asker = threading.Thread(target=lambda: results.append(ask(stream, wait_sec=2.0)))
    asker.start()
    while not stream._wanted:
        pass
    frames = stream._outgoing_frames()
    with stream._cond:
        stream._handle(answer(frames[0]["nonce"]))
        stream._cond.notify_all()
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


def test_a_novig_books_ask_is_priority_or_background_never_both():
    assert service.validate_novig_books_request({"marketIds": ["mk-a"], "background": True}) == (["mk-a"], 0.0, False, True)
    assert service.validate_novig_books_request({"marketIds": ["mk-a"], "priority": True}) == (["mk-a"], 0.0, True, False)
    assert service.validate_novig_books_request({"marketIds": ["mk-a"], "priority": True, "background": True}) == "priority and background cannot both be true"
    assert "background must be" in service.validate_novig_books_request({"marketIds": ["mk-a"], "background": "yes"})


def test_post_novig_books_answers_the_streams_books(bets_store, made):
    stream, _connection, _clock, _client = made
    subscribe(stream)
    server, base = _serve(service.make_handler(bets_store, 0.0, ["kalshi"], novig_book_stream=stream,
                                               extra_hosts=("mlb-stack.test",)))
    try:
        status, body = _post(f"{base}/novig_books.json", {"marketIds": [MARKET, MARKET], "waitMs": 0, "priority": True})
        bad = _post(f"{base}/novig_books.json", {"marketIds": [""]})
        foreign = _post(f"{base}/novig_books.json", {"marketIds": []}, host="mlb-stack.test")
    finally:
        server.shutdown()
        server.server_close()
    assert status == 200 and body["books"] == {MARKET: {"status": "OPEN", "started": False, "seq": 48120, "orders": {DUQ: [O1]}}}
    assert body["stream"]["connected"] is True
    assert bad[0] == 400 and "marketIds[0]" in bad[1]["error"]
    assert foreign[0] == 403, "the runner's book route is loopback-only"


def test_post_novig_books_is_503_without_the_bot_key(bets_store):
    server, base = _serve(service.make_handler(bets_store, 0.0, ["kalshi"]))
    try:
        status, body = _post(f"{base}/novig_books.json", {"marketIds": [], "waitMs": 0})
    finally:
        server.shutdown()
        server.server_close()
    assert status == 503 and "NOVIG_TAKER_KEY_ID" in body["error"]


def test_an_answer_that_leaves_a_book_out_of_sync_is_asked_again_a_minute_later_not_every_loop(made):
    stream, _connection, clock, _client = made
    for name, body in (("mk-lifecycle-only", {"eventId": EVENT, "lifecycle": {"seq": 1, "status": "OPEN"}}),
                       ("mk-unreadable", {"eventId": EVENT, "book": {"seq": 1, "orders": {DUQ: [{"order": "x"}]}}})):
        ask(stream, name)
        handle(stream, answer(stream._outgoing_frames()[0]["nonce"], {name: body}))
        assert "out of sync" in ask(stream, name)["error"]
        assert stream._outgoing_frames() == [], f"{name}: no snapshot re-sent at once"
    clock.now += nbs.REFUSED_RETRY_SEC + 1
    stream.books_for(["mk-lifecycle-only", "mk-unreadable"], 0.0)
    assert sorted(market for frame in stream._outgoing_frames() for market in frame["snapshot"]["markets"]) == [
        "mk-lifecycle-only", "mk-unreadable"]


def test_a_gap_snapshot_answer_never_subscribes_a_market_dropped_meanwhile(made):
    stream, _connection, clock, _client = made
    subscribe(stream)
    handle(stream, book_delta(48125, []))
    healing = stream._outgoing_frames()[0]["nonce"]
    clock.now += nbs.WATCH_IDLE_SEC + 1
    handle(stream, {"heartbeat": {}})
    assert stream._outgoing_frames()[0]["unsubscribe"] == [f"market:{MARKET}"]
    ask(stream)
    handle(stream, answer(healing))
    assert MARKET not in stream._books and MARKET not in stream._subscribed
    assert [frame.get("subscribe") for frame in stream._outgoing_frames()] == [{"markets": {MARKET: "book"}}]


def test_a_go_live_survives_an_unreadable_book_part_in_its_frame(made):
    stream, _connection, _clock, _client = made
    subscribe(stream)
    handle(stream, book_delta(48121, [{"kind": "replace"}], lifecycle={"seq": 4, "deltas": ["GOLIVE"]}))
    assert stream._books[MARKET].started is True


def test_a_batch_refused_as_a_whole_goes_again_one_market_a_frame(made):
    stream, _connection, _clock, _client = made
    stream.books_for(["mk-a", "mk-b"], 0.0)
    nonce = stream._outgoing_frames()[0]["nonce"]
    handle(stream, {"nonce": nonce, "code": "MARKET_NOT_FOUND", "message": "mk-b"})
    assert [frame["subscribe"] for frame in stream._outgoing_frames()] == [{"markets": {"mk-a": "book"}}, {"markets": {"mk-b": "book"}}]


def test_ids_named_order_id_and_outcome_id_as_in_the_live_frames_are_read_too(made):
    stream, _connection, _clock, _client = made
    ask(stream)
    rest_named = {"orderId": "o1", "price": "0.445", "qty": 110}
    handle(stream, answer(stream._outgoing_frames()[0]["nonce"], markets={MARKET: market_snapshot(orders={DUQ: [rest_named]})}))
    assert ask(stream)["orders"] == {DUQ: [O1]}
    handle(stream, book_delta(48121, [{"kind": "add", "orderId": "o2", "outcomeId": WSU, "price": "0.550", "qty": 40},
                                      {"kind": "remove", "orderId": "o1", "reason": "fill"}]))
    assert ask(stream)["orders"] == {WSU: [{"order": "o2", "price": "0.550", "qty": 40}]}


def test_a_remove_naming_no_order_resyncs_the_market_and_logs_the_part(made, caplog):
    stream, _connection, _clock, _client = made
    subscribe(stream)
    handle(stream, book_delta(48121, [{"kind": "remove", "reason": "fill"}]))
    assert "resyncing" in ask(stream)["error"], "a remove that names no order must not leave the filled order showing"
    warning = next(record.getMessage() for record in caplog.records if "unreadable frame part" in record.getMessage())
    assert "expected an order id under order or orderId, got keys ['kind', 'reason']" in warning
    assert '"reason": "fill"' in warning, "the part itself is logged so a new shape can be read off the log"


def test_an_update_delta_sets_what_is_left_of_an_order_and_zero_takes_it_off(made):
    stream, _connection, _clock, _client = made
    subscribe(stream)
    handle(stream, book_delta(48121, [{"kind": "update", "orderId": "o1", "remaining": 57}]))
    assert ask(stream)["orders"] == {DUQ: [{"order": "o1", "price": "0.445", "qty": 57}]}
    handle(stream, book_delta(48122, [{"kind": "update", "orderId": "o1", "remaining": 0}]))
    assert ask(stream)["orders"] == {}


def test_an_update_for_an_order_the_book_never_saw_resyncs_the_market(made):
    stream, _connection, _clock, _client = made
    subscribe(stream)
    handle(stream, book_delta(48121, [{"kind": "update", "orderId": "o-unknown", "remaining": 5}]))
    assert "resyncing" in ask(stream)["error"]


# ---- the private orders channel ---------------------------------------------------------

ORDER = "7c9e6679-0000-4000-8000-000000000001"
CLICK = "0b3f5f20-0000-4000-8000-000000000001"


def fill_delta(seq, qty=40, remaining=70, fill_id="1e6a2f94-1", price="0.670", order=ORDER):
    return {"ts": 2, "delta": {}, "orders": {"seq": seq, "deltas": [
        {"kind": "fill", "fillId": fill_id, "orderId": order, "clientId": CLICK, "outcomeId": DUQ, "price": price,
         "qty": qty, "originalQty": 110, "remaining": remaining, "taker": True}]}}


def orders_answer(nonce, seq=907, open_orders=(), with_ack=True):
    frame = {"ts": 1, "nonce": nonce, "snapshot": {}, "orders": {"seq": seq, "open": list(open_orders)}}
    if with_ack:
        frame["subscribed"] = {"markets": {}, "events": {}, "private": ["orders"]}
    return frame


@pytest.fixture
def orders_made():
    clock = Clock()
    client = FakeClient()
    connection = FakeConnection()
    stream = nbs.NovigBookStream(client, "wss://novig.test/v3/ws", connect=lambda url, headers: connection, clock=clock,
                                 orders_channel=True, wall_clock_ms=lambda: 1762600000000)
    stream._open()
    return stream, connection, clock


def subscribe_orders(stream, seq=907):
    [frame] = stream._outgoing_frames()
    assert frame == {"nonce": 1, "subscribe": {"private": ["orders"]}}
    handle(stream, orders_answer(1, seq))
    return frame


def test_the_orders_channel_is_subscribed_first_thing_for_one_token_and_the_snapshot_syncs_it(orders_made):
    stream, _connection, _clock = orders_made
    assert not stream.orders_live()
    ask(stream)  # a market asked for at the same time goes behind the private subscribe
    frames = stream._outgoing_frames()
    assert frames[0] == {"nonce": 1, "subscribe": {"private": ["orders"]}}
    assert "subscribe" in frames[1] and "markets" in frames[1]["subscribe"]
    assert stream._bucket.available() == pytest.approx(nbs.STREAM_CAPACITY - nbs.UPGRADE_COST - 1 - nbs.BOOK_PAIR_COST)
    handle(stream, orders_answer(1, 907, open_orders=[{"orderId": ORDER, "clientId": CLICK, "marketId": MARKET, "outcomeId": DUQ,
                                                       "price": "0.665", "qty": 110, "tif": "GTC"}]))
    assert stream.orders_live()
    status = stream.status()["orders"]
    assert (status["subscribed"], status["synced"], status["seq"], status["epoch"]) == (True, True, 907, 1)
    resting, version = stream.wait_for_order(ORDER, 0, 0.0)
    assert resting["final"] == "OPEN" and resting["remaining"] == 110 and resting["filled"] == 0 and version == 1
    assert stream.fill_events_after(0) == {"events": [], "version": 1, "live": True, "epoch": 1}


def test_a_stream_without_the_channel_sends_no_private_subscribe(made):
    stream, _connection, _clock, _client = made
    assert stream._outgoing_frames() == []
    assert "orders" not in stream.status()
    assert not stream.orders_live()
    assert stream.wait_for_order(ORDER, 0, 0.0) == (None, 0)


def test_fills_then_a_cancel_or_a_full_fill_make_the_orders_state_and_fill_events(orders_made):
    stream, _connection, _clock = orders_made
    subscribe_orders(stream)
    handle(stream, fill_delta(908, qty=40, remaining=70))
    state, version = stream.wait_for_order(ORDER, 0, 0.0)
    assert state["final"] is None and state["filled"] == 40 and state["cost"] == pytest.approx(0.268) and state["remaining"] == 70
    assert stream.wait_for_order(ORDER, version, 0.0) == (None, version), "nothing newer yet"
    # An IOC's unfilled rest is cancelled: fill then cancel, one batch each.
    handle(stream, {"ts": 3, "delta": {}, "orders": {"seq": 909, "deltas": [
        {"kind": "cancel", "orderId": ORDER, "clientId": CLICK, "outcomeId": DUQ, "originalQty": 110, "remaining": 70}]}})
    state, version2 = stream.wait_for_order(ORDER, version, 0.0)
    assert state["final"] == "CANCELED" and state["filled"] == 40 and version2 > version
    # Another order filled whole: remaining 0 is FILLED whether or not Novig sends more.
    other = "7c9e6679-0000-4000-8000-000000000002"
    handle(stream, fill_delta(910, qty=110, remaining=0, fill_id="1e6a2f94-2", order=other))
    state, _ = stream.wait_for_order(other, 0, 0.0)
    assert state["final"] == "FILLED" and state["filled"] == 110 and state["cost"] == pytest.approx(0.737)
    drained = stream.fill_events_after(0)
    assert [(event["fillId"], event["orderId"], event["qty"], event["ts"]) for event in drained["events"]] == [
        ("1e6a2f94-1", ORDER, 40, 1762600000000), ("1e6a2f94-2", other, 110, 1762600000000)]
    assert drained["version"] == stream.status()["orders"]["version"]
    assert stream.fill_events_after(drained["version"])["events"] == []
    # A reject: an IOC that filled nothing.
    third = "7c9e6679-0000-4000-8000-000000000003"
    handle(stream, {"ts": 4, "delta": {}, "orders": {"seq": 911, "deltas": [
        {"kind": "reject", "orderId": third, "clientId": CLICK, "outcomeId": DUQ, "originalQty": 50}]}})
    state, _ = stream.wait_for_order(third, 0, 0.0)
    assert state["final"] == "REJECTED" and state["filled"] == 0 and state["remaining"] == 50


def test_a_repeated_fill_id_is_counted_once(orders_made):
    stream, _connection, _clock = orders_made
    subscribe_orders(stream)
    handle(stream, fill_delta(908), {"ts": 2, "delta": {}, "orders": {"seq": 909, "deltas": [fill_delta(0)["orders"]["deltas"][0]]}})
    state, _ = stream.wait_for_order(ORDER, 0, 0.0)
    assert state["filled"] == 40 and len(stream.fill_events_after(0)["events"]) == 1


def test_a_skipped_orders_seq_or_a_heartbeat_ahead_asks_for_a_snapshot_and_bumps_the_epoch(orders_made):
    stream, _connection, clock = orders_made
    subscribe_orders(stream)
    handle(stream, fill_delta(910))  # 908 and 909 never came
    assert not stream.orders_live()
    assert stream.fill_events_after(0)["live"] is False
    assert stream.wait_for_order(ORDER, 0, 5.0) == (None, 0), "a wait on a channel out of sync returns at once"
    [frame] = stream._outgoing_frames()
    assert frame == {"nonce": 2, "snapshot": {"private": ["orders"]}}
    assert stream._outgoing_frames() == [], "one recovery in flight at a time"
    # The snapshot answer (the kit's mock shape: no `subscribed`): back in sync at its seq, the buffered 910 replayed.
    handle(stream, orders_answer(2, 909, with_ack=False))
    assert stream.orders_live() and stream.status()["orders"] == pytest.approx(
        {"live": True, "subscribed": True, "synced": True, "seq": 910, "epoch": 2, "version": 1, "error": None})
    assert stream.wait_for_order(ORDER, 0, 0.0)[0]["filled"] == 40
    # A heartbeat naming a seq past ours: a dropped frame.
    handle(stream, {"heartbeat": {"orders": 915, "ts": 1757894415123}})
    assert not stream.orders_live()
    assert stream._outgoing_frames() == [{"nonce": 3, "snapshot": {"private": ["orders"]}}]
    handle(stream, {"heartbeat": {"orders": 915, "ts": 1757894415124}})  # another while healing changes nothing
    handle(stream, orders_answer(3, 915, with_ack=False))
    assert stream.status()["orders"]["epoch"] == 3 and stream.orders_live()
    # A snapshot answer that never comes is asked for again after the wait.
    handle(stream, fill_delta(920, fill_id="late"))
    assert stream._outgoing_frames()[0]["nonce"] == 4
    clock.now += nbs.ANSWER_WAIT_SEC + 1
    assert stream._outgoing_frames()[0]["nonce"] == 5


def test_a_refused_private_subscribe_is_tried_again_a_minute_later_and_named_in_health(orders_made):
    stream, _connection, clock = orders_made
    stream._outgoing_frames()
    handle(stream, {"nonce": 1, "code": "SCOPE_INSUFFICIENT", "message": "api key scope is insufficient"})
    assert stream._outgoing_frames() == []
    assert stream.status()["orders"]["error"] == "SCOPE_INSUFFICIENT api key scope is insufficient"
    clock.now += nbs.REFUSED_RETRY_SEC + 1
    assert stream._outgoing_frames() == [{"nonce": 2, "subscribe": {"private": ["orders"]}}]
    # An answer that acks nothing private and carries no orders object is a refusal too.
    handle(stream, {"ts": 1, "nonce": 2, "subscribed": {"markets": {}, "events": {}, "private": []}, "snapshot": {}})
    assert not stream.orders_live() and "named no orders channel" in stream.status()["orders"]["error"]


def test_an_unreadable_orders_delta_resyncs_the_channel_and_logs_the_part(orders_made, caplog):
    stream, _connection, _clock = orders_made
    subscribe_orders(stream)
    handle(stream, {"ts": 2, "delta": {}, "orders": {"seq": 908, "deltas": [{"kind": "fill", "orderId": ORDER, "qty": "forty", "price": "0.67"}]}})
    assert not stream.orders_live()
    assert any("orders channel: unreadable frame part" in record.message and "forty" in record.message for record in caplog.records)
    assert stream._outgoing_frames() == [{"nonce": 2, "snapshot": {"private": ["orders"]}}]


def test_a_reconnect_subscribes_the_channel_again_and_a_wait_returns_at_once_while_down(orders_made):
    stream, connection, clock = orders_made
    subscribe_orders(stream)
    handle(stream, fill_delta(908))
    stream._drop("ConnectionError: closed")
    assert not stream.orders_live()
    assert stream.wait_for_order(ORDER, 0, 5.0)[0]["filled"] == 40, "a state already pushed is still answered"
    assert stream.wait_for_order("never-seen", 0, 5.0) == (None, 0)
    assert stream.fill_events_after(0)["live"] is False and len(stream.fill_events_after(0)["events"]) == 1
    stream._open()
    assert stream._outgoing_frames() == [{"nonce": 1, "subscribe": {"private": ["orders"]}}]
    handle(stream, orders_answer(1, 3))
    assert stream.status()["orders"]["epoch"] == 2, "a new epoch: the source reads REST once"
    assert stream.orders_live()


def test_a_silent_line_takes_the_channel_off_live(orders_made):
    stream, _connection, clock = orders_made
    subscribe_orders(stream)
    clock.now += nbs.STALE_VIEW_SEC + 1
    assert not stream.orders_live()


def test_wait_for_order_wakes_on_the_frame_from_the_stream_thread():
    clock = Clock()
    connection = FakeConnection()
    stream = nbs.NovigBookStream(FakeClient(), "wss://novig.test/v3/ws", connect=lambda url, headers: connection, clock=clock,
                                 orders_channel=True)
    stream._open()
    subscribe_orders(stream)
    got = []

    def waiter():
        got.append(stream.wait_for_order(ORDER, 0, 5.0))

    thread = threading.Thread(target=waiter)
    thread.start()
    import time as _time
    _time.sleep(0.05)
    handle(stream, fill_delta(908, qty=110, remaining=0))
    thread.join(timeout=5)
    assert got and got[0][0]["final"] == "FILLED"


def test_a_heartbeat_below_our_seq_is_a_novig_reset_and_asks_for_a_snapshot(orders_made):
    """Novig's seq resets on its restart; on a socket that stayed up every new
    delta (seq 1, 2, ...) would read as stale, so the heartbeat decides."""
    stream, _connection, _clock = orders_made
    subscribe_orders(stream)
    handle(stream, fill_delta(908))
    handle(stream, fill_delta(2, fill_id="after-reset"))  # ignored as stale, as a replay would be
    assert stream.orders_live() and len(stream.fill_events_after(0)["events"]) == 1
    handle(stream, {"heartbeat": {"orders": 2, "ts": 1757894415123}})
    assert not stream.orders_live()
    assert stream._outgoing_frames() == [{"nonce": 2, "snapshot": {"private": ["orders"]}}]
    handle(stream, orders_answer(2, 2, with_ack=False))
    handle(stream, fill_delta(3, fill_id="after-reset-2"))
    assert stream.orders_live() and stream.status()["orders"]["seq"] == 3
    assert [event["fillId"] for event in stream.fill_events_after(0)["events"]] == ["1e6a2f94-1", "after-reset-2"]


def test_a_fill_with_no_fill_id_is_unreadable_and_resyncs_the_channel(orders_made):
    """A made-up id would never match REST's copy, so the subaccount source would count the fill twice."""
    stream, _connection, _clock = orders_made
    subscribe_orders(stream)
    handle(stream, fill_delta(908))
    frame = fill_delta(909, fill_id=None)
    del frame["orders"]["deltas"][0]["fillId"]
    handle(stream, frame)
    assert not stream.orders_live()
    assert len(stream.fill_events_after(0)["events"]) == 1, "nothing counted off the unreadable frame"
    assert stream._outgoing_frames() == [{"nonce": 2, "snapshot": {"private": ["orders"]}}]
