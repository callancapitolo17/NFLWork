"""Kalshi's live books over the orderbook_delta websocket (kalshi_book_stream.py),
driven frame by frame with a fake connection and a fake clock, in the frame
shapes Kalshi's docs give (read 2026-11-10; the first deploy 2026-11-11 saw
them work): one subscribe a connection with a Ticket's/Place's markets first,
snapshots served as the REST /orderbook's body, deltas applied, seq gaps and
orphan deltas reconnected,
crossed books withheld, adds and deletes on a live subscription, refusals, the
background cap, silence; plus the handshake signature and POST
/kalshi_books.json over a real HTTP server."""
import base64
import json
import threading
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer

import pytest

from unabated_ticket.bets_service import kalshi_book_stream as kbs
from unabated_ticket.bets_service import service
from unabated_ticket.bets_service.store import BetsStore
from unabated_ticket.bets_service.taker_orders import KalshiOrderClient

WSU = "KXNCAAFSPREAD-26SEP19DUQWSU-WSU36"
DUQ = "KXNCAAFSPREAD-26SEP19DUQWSU-DUQ7"


class Clock:
    def __init__(self):
        self.now = 1000.0

    def __call__(self) -> float:
        return self.now


class FakeSigner:
    key_id = "d518d266-test"

    def __init__(self):
        self.signed = []

    def signed_ws_headers(self, ws_path):
        self.signed.append(ws_path)
        return {"KALSHI-ACCESS-KEY": self.key_id, "KALSHI-ACCESS-TIMESTAMP": "1", "KALSHI-ACCESS-SIGNATURE": "sig"}


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


URL = "wss://api.elections.kalshi.com/trade-api/ws/v2"


@pytest.fixture
def made():
    clock = Clock()
    signer = FakeSigner()
    connection = FakeConnection()
    stream = kbs.KalshiBookStream(signer, URL, connect=lambda url, headers: connection, clock=clock)
    stream._open()
    return stream, connection, clock, signer


def handle(stream, *messages):
    with stream._cond:
        for message in messages:
            stream._handle(message)


def tick(stream, clock, sec=0.3):
    """The pump's next send after `sec` (past the housekeeping and command gaps)."""
    clock.now += sec
    return stream._outgoing_frames()


def snapshot(sid, seq, ticker=WSU, yes=(("0.5500", "100.00"),), no=(("0.4300", "40.00"),)):
    msg = {"market_ticker": ticker, "market_id": "m-1", "clock": 7}
    if yes:
        msg["yes_dollars_fp"] = [list(level) for level in yes]
    if no:
        msg["no_dollars_fp"] = [list(level) for level in no]
    return {"type": "orderbook_snapshot", "sid": sid, "seq": seq, "msg": msg}


def delta(sid, seq, price, change, side="no", ticker=WSU):
    return {"type": "orderbook_delta", "sid": sid, "seq": seq,
            "msg": {"market_ticker": ticker, "market_id": "m-1", "clock": 8, "price_dollars": price,
                    "delta_fp": change, "side": side, "ts_ms": 1}}


def subscribed(command_id, sid):
    return {"id": command_id, "type": "subscribed", "msg": {"channel": "orderbook_delta", "sid": sid}}


def subscribe(stream, clock, ticker=WSU, sid=11, **book):
    """Ask for `ticker`, send its subscribe, answer it and snapshot the market; returns the subscribe frame."""
    stream.books_for([ticker], 0.0)
    (frame,) = tick(stream, clock)
    handle(stream, subscribed(frame["id"], sid), snapshot(sid, 1, ticker, **book))
    return frame


def reconnect(stream, why):
    """The pump's planned reconnect: the stream asked for one (`why` in its reason), ends it and opens anew."""
    assert why in (stream._reconnect_reason or ""), stream._reconnect_reason
    stream._drop_for_reconnect(stream._reconnect_reason)
    stream._open()


def add(stream, clock, ticker=DUQ):
    """Ask for `ticker` with WSU's subscription answered: the add it sends."""
    stream.books_for([ticker], 0.0)
    (frame,) = tick(stream, clock)
    assert frame["params"]["action"] == "add_markets"
    return frame


def view(stream, ticker=WSU):
    return stream.books_for([ticker], 0.0)["books"][ticker]


def test_the_upgrade_is_signed_over_the_websocket_path(made):
    _stream, _connection, _clock, signer = made
    assert signer.signed == ["/trade-api/ws/v2"]


def test_the_read_key_signs_the_handshake_as_kalshi_documents(tmp_path):
    """timestamp + "GET" + "/trade-api/ws/v2", RSA-PSS SHA-256, base64 (docs quick_start_websockets)."""
    from cryptography.hazmat.primitives import hashes, serialization
    from cryptography.hazmat.primitives.asymmetric import padding, rsa
    key = rsa.generate_private_key(public_exponent=65537, key_size=2048)
    key_path = tmp_path / "read.pem"
    key_path.write_bytes(key.private_bytes(serialization.Encoding.PEM, serialization.PrivateFormat.PKCS8,
                                           serialization.NoEncryption()))
    client = KalshiOrderClient("d518d266-test", str(key_path), "https://kalshi.test/trade-api/v2")
    headers = client.signed_ws_headers("/trade-api/ws/v2")
    assert headers["KALSHI-ACCESS-KEY"] == "d518d266-test" and client.key_id == "d518d266-test"
    message = f"{headers['KALSHI-ACCESS-TIMESTAMP']}GET/trade-api/ws/v2".encode()
    key.public_key().verify(base64.b64decode(headers["KALSHI-ACCESS-SIGNATURE"]), message,
                            padding.PSS(mgf=padding.MGF1(hashes.SHA256()), salt_length=padding.PSS.DIGEST_LENGTH),
                            hashes.SHA256())


def test_a_new_market_is_subscribed_and_its_snapshot_served_as_the_rest_orderbook(made):
    stream, _connection, clock, _signer = made
    assert "waiting for" in view(stream)["error"]
    assert stream._outgoing_frames() == [{"id": 1, "cmd": "subscribe",
                                          "params": {"channels": ["orderbook_delta"], "market_tickers": [WSU]}}]
    handle(stream, subscribed(1, 11))
    assert "waiting for" in view(stream)["error"], "no book before the snapshot"
    handle(stream, snapshot(11, 1, yes=(("0.5500", "100.00"), ("0.5600", "5.00"))))
    clock.now += 2
    assert view(stream) == {"orderbook_fp": {"yes_dollars": [["0.5600", "5.00"], ["0.5500", "100.00"]],
                                             "no_dollars": [["0.4300", "40.00"]]}, "ageMs": 2000}
    assert tick(stream, clock) == [], "subscribed once, not again every loop"


def test_an_empty_side_is_absent_from_the_snapshot_and_served_empty(made):
    stream, _connection, clock, _signer = made
    subscribe(stream, clock, no=())
    assert view(stream)["orderbook_fp"]["no_dollars"] == []


def test_a_snapshot_whose_sides_are_renamed_stops_the_connection_rather_than_serve_an_empty_book(made):
    stream, _connection, clock, _signer = made
    stream.books_for([WSU], 0.0)
    (frame,) = tick(stream, clock)
    renamed = {"type": "orderbook_snapshot", "sid": 11, "seq": 1,
               "msg": {"market_ticker": WSU, "yes": [[55, 100]], "no": [[43, 40]]}}
    with pytest.raises(TypeError, match="expected Kalshi snapshot sides"):
        handle(stream, subscribed(frame["id"], 11), renamed)


def test_deltas_change_add_and_remove_levels(made):
    stream, _connection, clock, _signer = made
    subscribe(stream, clock)
    clock.now += 5
    handle(stream, delta(11, 2, "0.4300", "-15.00"), delta(11, 3, "0.4400", "7.00"),
           delta(11, 4, "0.5500", "-100.00", side="yes"))
    assert view(stream) == {"orderbook_fp": {"yes_dollars": [], "no_dollars": [["0.4400", "7.00"], ["0.4300", "25.00"]]},
                            "ageMs": 0}


def test_a_seq_gap_withholds_the_books_and_reconnects_after_a_doubling_gap(made):
    """Kalshi ignores a second subscribe on a connection, so fresh snapshots come off a fresh one."""
    stream, connection, clock, _signer = made
    subscribe(stream, clock)
    handle(stream, delta(11, 3, "0.4300", "-15.00"))
    assert "missed a message" in view(stream)["error"]
    assert stream.status()["lastResync"] == "seq 3 after 1"
    assert tick(stream, clock, 0.5) == [] and stream._reconnect_reason is None, "the first reconnect waits RESYNC_MIN_GAP_SEC"
    assert tick(stream, clock, 0.6) == [], "no unsubscribe and no second subscribe on this connection"
    reconnect(stream, "missed a message (seq 3 after 1)")
    assert "waiting for" in view(stream)["error"], "no book from the old connection"
    (frame,) = tick(stream, clock)
    assert frame == {"id": 1, "cmd": "subscribe", "params": {"channels": ["orderbook_delta"], "market_tickers": [WSU]}}
    handle(stream, subscribed(1, 12), snapshot(12, 1, no=(("0.4200", "9.00"),)))
    assert view(stream)["orderbook_fp"]["no_dollars"] == [["0.4200", "9.00"]]
    # A second gap soon after waits twice as long.
    handle(stream, delta(12, 5, "0.4200", "1.00"))
    assert tick(stream, clock, 1.5) == [] and stream._reconnect_reason is None
    tick(stream, clock, 0.6)
    reconnect(stream, "missed a message")
    assert stream.status()["resyncs"] == 2 and stream.status()["plannedReconnects"] == 2
    assert stream.status()["disconnects"] == 0 and stream.status()["error"] is None


def test_a_delta_before_its_snapshot_withholds_the_books(made):
    stream, _connection, clock, _signer = made
    subscribe(stream, clock)
    add(stream, clock)
    handle(stream, {"id": 2, "sid": 11, "seq": 2, "type": "ok", "msg": {"market_tickers": [WSU, DUQ]}})
    handle(stream, delta(11, 3, "0.4300", "1.00", ticker=DUQ))
    assert "missed a message" in view(stream)["error"]


def test_a_level_driven_below_zero_withholds_the_books(made):
    stream, _connection, clock, _signer = made
    subscribe(stream, clock)
    handle(stream, delta(11, 2, "0.4300", "-41.00"))
    assert "missed a message" in view(stream)["error"]


def test_a_crossed_book_is_never_served(made):
    stream, _connection, clock, _signer = made
    subscribe(stream, clock)
    handle(stream, delta(11, 2, "0.4500", "3.00", side="yes"))
    assert "orderbook_fp" in view(stream), "0.55 + 0.43 is not crossed"
    handle(stream, delta(11, 3, "0.5800", "3.00"))
    assert "crossed" in view(stream)["error"]
    handle(stream, delta(11, 4, "0.5800", "-3.00"))
    assert "orderbook_fp" in view(stream)


def test_a_new_market_is_added_to_the_answered_subscription(made):
    stream, _connection, clock, _signer = made
    subscribe(stream, clock)
    frame = add(stream, clock)
    assert frame == {"id": 2, "cmd": "update_subscription",
                     "params": {"sid": 11, "market_tickers": [DUQ], "action": "add_markets"}}
    # The docs' ok carries the sid's next seq: counted, so the next snapshot is no gap.
    handle(stream, {"id": 2, "sid": 11, "seq": 2, "type": "ok", "msg": {"market_tickers": [WSU, DUQ]}},
           snapshot(11, 3, DUQ))
    assert "orderbook_fp" in view(stream, DUQ) and "orderbook_fp" in view(stream, WSU)


def test_a_market_added_to_an_old_subscription_is_waited_for_from_its_add(made):
    """The snapshot wait runs from the add, not from when the subscription was first answered."""
    stream, _connection, clock, _signer = made
    subscribe(stream, clock)
    clock.now += 30
    add(stream, clock)
    assert "waiting for" in view(stream, DUQ)["error"]
    assert stream.books_for([DUQ], 0.0)["books"][DUQ] and stream._view_of(DUQ).get("pending") is True
    clock.now += kbs.ANSWER_WAIT_SEC + 1
    assert "sent no snapshot of" in view(stream, DUQ)["error"]


def test_one_subscription_grows_by_adds_one_command_at_a_time(made, monkeypatch):
    """Kalshi ignores a second orderbook_delta subscribe on a connection (VM, 2026-11-11), so the
    markets past the first command's ride add_markets on the one subscription."""
    stream, _connection, clock, _signer = made
    monkeypatch.setattr(kbs, "MARKETS_PER_COMMAND", 2)
    stream.books_for(["A", "B", "C", "D", "E"], 0.0)
    (first,) = tick(stream, clock)
    assert first["cmd"] == "subscribe" and first["params"]["market_tickers"] == ["A", "B"]
    assert tick(stream, clock) == [], "no second subscribe while the first is unanswered"
    handle(stream, subscribed(first["id"], 11))
    sent = [tick(stream, clock) for _ in range(3)]
    assert [[(frame["cmd"], frame["params"].get("action"), frame["params"]["market_tickers"]) for frame in frames]
            for frames in sent] == [[("update_subscription", "add_markets", ["C", "D"])],
                                    [("update_subscription", "add_markets", ["E"])], []]
    assert stream.status()["subscriptions"] == 1 and stream.status()["subscribedMarkets"] == 5
    assert stream._outgoing_frames() == [], "nothing within the housekeeping gap"


def test_a_renewal_reconnects_then_subscribes_one_command_and_adds_the_rest(made, monkeypatch):
    stream, _connection, clock, _signer = made
    monkeypatch.setattr(kbs, "MARKETS_PER_COMMAND", 1)
    monkeypatch.setattr(kbs, "RENEW_SEC", 30)
    subscribe(stream, clock)
    add(stream, clock)
    assert tick(stream, clock, kbs.RENEW_SEC) == [], "no unsubscribe or subscribe on the old connection"
    reconnect(stream, "renewal")
    (first,) = tick(stream, clock)
    assert first["cmd"] == "subscribe" and first["params"]["market_tickers"] == [DUQ], "one command's worth"
    assert tick(stream, clock) == []
    handle(stream, subscribed(first["id"], 12))
    (again,) = tick(stream, clock)
    assert again["params"] == {"sid": 12, "market_tickers": [WSU], "action": "add_markets"}
    assert tick(stream, clock) == []
    assert stream.status()["subscriptions"] == 1 and stream.status()["subscribedMarkets"] == 2


def test_an_ok_whose_seq_is_not_next_is_not_counted(made):
    stream, _connection, clock, _signer = made
    subscribe(stream, clock)
    handle(stream, {"id": 9, "sid": 11, "seq": 50, "type": "ok"}, delta(11, 2, "0.4300", "1.00"))
    assert view(stream)["orderbook_fp"]["no_dollars"] == [["0.4300", "41.00"]]


def test_an_idle_market_is_deleted_and_a_late_frame_for_it_is_ignored(made, monkeypatch):
    stream, _connection, clock, _signer = made
    monkeypatch.setattr(kbs, "RENEW_SEC", 3 * kbs.WATCH_IDLE_SEC)
    subscribe(stream, clock)
    add(stream, clock)
    handle(stream, {"id": 2, "sid": 11, "seq": 2, "type": "ok"}, snapshot(11, 3, DUQ))
    clock.now += kbs.WATCH_IDLE_SEC - 1
    stream.books_for([WSU], 0.0)
    assert tick(stream, clock, 2) == [{"id": 3, "cmd": "update_subscription",
                                       "params": {"sid": 11, "market_tickers": [DUQ], "action": "delete_markets"}}]
    handle(stream, delta(11, 4, "0.4300", "1.00", ticker=DUQ), delta(11, 5, "0.4300", "1.00"))
    assert view(stream)["orderbook_fp"]["no_dollars"] == [["0.4300", "41.00"]], "seq still counted, no resync"
    clock.now += kbs.WATCH_IDLE_SEC
    assert tick(stream, clock) == [{"id": 4, "cmd": "unsubscribe", "params": {"sids": [11]}}]
    assert stream.status()["subscriptions"] == 0


def test_a_market_re_asked_before_its_delete_goes_back_to_its_own_subscription(made, monkeypatch):
    stream, _connection, clock, _signer = made
    monkeypatch.setattr(kbs, "RENEW_SEC", 3 * kbs.WATCH_IDLE_SEC)
    subscribe(stream, clock)
    clock.now += kbs.WATCH_IDLE_SEC + 1
    with stream._cond:
        stream._forget_idle(clock.now)
    stream.books_for([WSU], 0.0)
    assert tick(stream, clock) == [], "no delete and no second subscription"
    assert "orderbook_fp" in view(stream)


def test_a_market_limit_on_an_add_caps_the_subscription_at_its_size_and_frees_room_as_markets_leave(made, monkeypatch):
    stream, _connection, clock, _signer = made
    monkeypatch.setattr(kbs, "RENEW_SEC", 3 * kbs.WATCH_IDLE_SEC)
    subscribe(stream, clock)
    frame = add(stream, clock)
    handle(stream, {"id": frame["id"], "type": "error", "msg": {"code": 26, "msg": "per-subscription market limit exceeded"}})
    assert stream.status()["marketsPerSubscription"] == 1
    assert "orderbook_fp" in view(stream, WSU), "the subscription's other books stay served"
    assert tick(stream, clock) == [] and stream._reconnect_reason is None, "no second subscribe, no reconnect"
    assert "takes no more markets (at most 1" in view(stream, DUQ)["error"]
    # WSU goes idle and leaves; DUQ, still asked for, takes its seat on the same subscription.
    clock.now += kbs.WATCH_IDLE_SEC - 1
    stream.books_for([DUQ], 0.0)
    clock.now += 2
    frames = [tick(stream, clock), tick(stream, clock)]
    assert [(f["params"]["action"], f["params"]["market_tickers"]) for fs in frames for f in fs] == [
        ("delete_markets", [WSU]), ("add_markets", [DUQ])]


def test_a_market_limit_on_a_subscribe_halves_it_and_keeps_the_best_placed_markets(made, monkeypatch):
    stream, _connection, clock, _signer = made
    monkeypatch.setattr(kbs, "MARKETS_PER_COMMAND", 4)
    monkeypatch.setattr(kbs, "MIN_MARKETS_PER_SUBSCRIPTION", 1)
    stream.books_for(["A", "B", "C", "D", "E", "F"], 0.0, background=True)
    stream.books_for(["Z-PLACE"], 0.0, priority=True)
    (frame,) = tick(stream, clock)
    assert frame["params"]["market_tickers"] == ["Z-PLACE", "A", "B", "C"]
    handle(stream, {"id": frame["id"], "type": "error", "msg": {"code": 26, "msg": "per-subscription market limit exceeded"}})
    assert stream.status()["marketsPerSubscription"] == 2
    assert tick(stream, clock) == []
    tick(stream, clock, kbs.RESYNC_MIN_GAP_SEC + 0.1)
    reconnect(stream, "a new subscription needs a fresh connection")
    (again,) = tick(stream, clock)
    assert again["cmd"] == "subscribe" and again["params"]["market_tickers"] == ["Z-PLACE", "A"]
    handle(stream, subscribed(again["id"], 11))
    assert tick(stream, clock) == [], "the rest wait for room"
    assert stream.status()["marketsPerSubscription"] == 2, "the cap outlives the connection"
    assert "takes no more markets" in view(stream, "B")["error"]


def test_a_market_limit_on_a_small_subscribe_refuses_its_markets_instead_of_halving_for_ever(made):
    stream, _connection, clock, _signer = made
    stream.books_for(["A", "B"], 0.0)
    (frame,) = tick(stream, clock)
    handle(stream, {"id": frame["id"], "type": "error", "msg": {"code": 26, "msg": "per-subscription market limit exceeded"}})
    assert stream.status()["marketsPerSubscription"] == kbs.MARKETS_PER_SUBSCRIPTION
    assert "code 26" in view(stream, "A")["error"]
    assert tick(stream, clock) == [], "refused for a minute, not sent again at once"


def test_an_error_for_a_command_of_an_older_subscription_is_ignored(made):
    stream, _connection, clock, _signer = made
    subscribe(stream, clock)
    old_add = add(stream, clock)
    with stream._cond:
        stream._groups[0]["subscribeId"] = 99
    handle(stream, {"id": old_add["id"], "type": "error", "msg": {"code": 7, "msg": "Unknown subscription ID"}})
    assert "orderbook_fp" in view(stream, WSU) and "refused" not in view(stream, DUQ)["error"]


def test_a_reconnect_further_off_than_any_wait_is_not_waited_for(made):
    stream, _connection, clock, _signer = made
    subscribe(stream, clock)
    handle(stream, delta(11, 13, "0.4300", "1.00"))
    with stream._cond:
        stream._groups[0]["resyncAfter"] = clock.now + 2
    answer = stream._view_of(WSU, 1.5)
    assert "missed a message" in answer["error"] and "pending" not in answer
    assert stream._view_of(WSU, 2.5).get("pending") is True, "due within the wait: worth waiting for"


def test_a_rate_limit_pauses_commands_and_asks_the_markets_again(made):
    stream, _connection, clock, _signer = made
    subscribe(stream, clock)
    frame = add(stream, clock)
    handle(stream, {"id": frame["id"], "type": "error", "msg": {"code": 27, "msg": "exceeded its command rate limit"}})
    assert "waiting for" in view(stream, DUQ)["error"], "not refused"
    assert tick(stream, clock, 1) == []
    assert [f["params"]["action"] for f in tick(stream, clock, kbs.RATE_LIMIT_PAUSE_SEC)] == ["add_markets"]


def test_a_rate_limited_subscribe_is_sent_again_on_a_fresh_connection_after_the_pause(made):
    stream, _connection, clock, _signer = made
    stream.books_for([WSU], 0.0)
    (frame,) = tick(stream, clock)
    handle(stream, {"id": frame["id"], "type": "error", "msg": {"code": 27, "msg": "exceeded its command rate limit"}})
    assert tick(stream, clock, 1) == [] and stream._reconnect_reason is None, "paused"
    tick(stream, clock, kbs.RATE_LIMIT_PAUSE_SEC)
    tick(stream, clock, kbs.RESYNC_MIN_GAP_SEC + 0.1)
    reconnect(stream, "a new subscription needs a fresh connection")
    assert [f["cmd"] for f in tick(stream, clock)] == ["subscribe"]


def test_a_rate_limit_pause_outlives_a_reconnect(made):
    stream, _connection, clock, _signer = made
    stream.books_for([WSU], 0.0)
    (frame,) = tick(stream, clock)
    handle(stream, {"id": frame["id"], "type": "error", "msg": {"code": 27, "msg": "exceeded its command rate limit"}})
    stream._drop("ConnectionError: closed")
    stream._open()
    assert tick(stream, clock, 1) == [], "still paused on the new connection"
    assert [f["cmd"] for f in tick(stream, clock, kbs.RATE_LIMIT_PAUSE_SEC)] == ["subscribe"]


def test_a_refused_subscribe_refuses_its_markets_then_asks_again_after_a_minute(made):
    stream, _connection, clock, _signer = made
    stream.books_for([WSU], 0.0)
    tick(stream, clock)
    handle(stream, {"id": 1, "type": "error", "msg": {"code": 8, "msg": "Unknown channel name"}})
    assert view(stream)["error"] == "Kalshi refused the stream subscription to " + WSU + ": code 8: Unknown channel name"
    assert tick(stream, clock) == [] and stream._reconnect_reason is None
    tick(stream, clock, kbs.REFUSED_RETRY_SEC + 1)
    tick(stream, clock, kbs.RESYNC_MIN_GAP_SEC + 0.1)
    reconnect(stream, "a new subscription needs a fresh connection")
    assert [frame["cmd"] for frame in tick(stream, clock)] == ["subscribe"]


def test_a_subscribe_answered_to_an_older_one_of_ours_is_unsubscribed(made):
    stream, _connection, clock, _signer = made
    stream.books_for([WSU], 0.0)
    (frame,) = tick(stream, clock)
    with stream._cond:
        stream._groups[0]["subscribeId"] = 99
        stream._commands[99] = stream._commands[frame["id"]]
    handle(stream, subscribed(frame["id"], 5))
    assert tick(stream, clock) == [{"id": 2, "cmd": "unsubscribe", "params": {"sids": [5]}}], "sent by the pump, not under the lock"


def test_an_unanswered_subscribe_drops_the_connection(made):
    stream, _connection, clock, _signer = made
    stream.books_for([WSU], 0.0)
    tick(stream, clock)
    assert tick(stream, clock, kbs.ANSWER_WAIT_SEC - 1) == [], "never sent again on this connection"
    with pytest.raises(ConnectionError, match="did not answer the subscribe of 1 market"):
        tick(stream, clock, 2)


def test_the_connection_is_renewed_for_fresh_snapshots(made):
    stream, _connection, clock, _signer = made
    subscribe(stream, clock)
    stream.books_for([WSU], 0.0)
    assert tick(stream, clock, kbs.RENEW_SEC) == []
    assert stream._reconnect_reason == "renewal"


def test_markets_wanted_again_after_the_subscription_emptied_come_on_a_fresh_connection(made, monkeypatch):
    stream, _connection, clock, _signer = made
    monkeypatch.setattr(kbs, "RENEW_SEC", 3 * kbs.WATCH_IDLE_SEC)
    subscribe(stream, clock)
    clock.now += kbs.WATCH_IDLE_SEC + 1
    assert tick(stream, clock) == [{"id": 2, "cmd": "unsubscribe", "params": {"sids": [11]}}]
    stream.books_for([DUQ], 0.0)
    assert tick(stream, clock) == [], "no second subscribe on this connection"
    tick(stream, clock, kbs.RESYNC_MIN_GAP_SEC + 0.1)
    reconnect(stream, "a new subscription needs a fresh connection")
    assert tick(stream, clock)[0]["params"]["market_tickers"] == [DUQ]


def test_a_planned_reconnect_ends_the_pump_and_run_opens_a_new_connection(monkeypatch):
    clock = Clock()
    connections = []

    def connect(url, headers):
        connections.append(FakeConnection())
        return connections[-1]

    stream = kbs.KalshiBookStream(FakeSigner(), URL, connect=connect, clock=clock)
    stop = threading.Event()
    waits = []

    def fake_wait(sec):
        waits.append(sec)
        if len(connections) >= 2:
            stop.set()
        return stop.is_set()

    monkeypatch.setattr(stop, "wait", fake_wait)
    with stream._cond:
        stream._reconnect_reason = "renewal"
    original_open = stream._open

    def open_and_ask():
        original_open()
        with stream._cond:
            stream._reconnect_reason = "renewal"

    monkeypatch.setattr(stream, "_open", open_and_ask)
    stream.run(stop)
    assert len(connections) == 2 and connections[0].closed
    assert waits[0] == kbs.RECONNECT_MIN_SEC
    status = stream.status()
    assert status["plannedReconnects"] >= 1 and status["lastPlannedReconnect"] == "renewal"
    assert status["disconnects"] == 0


def test_a_reread_or_place_market_is_placed_ahead_of_the_live_pass(made, monkeypatch):
    stream, _connection, clock, _signer = made
    monkeypatch.setattr(kbs, "MARKETS_PER_SUBSCRIPTION", 2)
    stream.books_for(["Z-LIVE-A", "A-LIVE-B"], 0.0, background=True)
    stream.books_for([DUQ], 0.0, priority=True)
    (frame,) = tick(stream, clock)
    assert frame["params"]["market_tickers"] == [DUQ, "A-LIVE-B"], "the Place market first"
    assert "takes no more markets (at most 2" in view(stream, "Z-LIVE-A")["error"]


def test_the_live_pass_is_capped_below_the_subscriptions_room(made, monkeypatch):
    stream, _connection, _clock, _signer = made
    monkeypatch.setattr(kbs, "MAX_BACKGROUND_MARKETS", 1)
    books = stream.books_for([WSU, DUQ], 0.0, background=True)["books"]
    assert "waiting for" in books[WSU]["error"]
    assert "already watches 1 Kalshi markets" in books[DUQ]["error"]


def test_a_dropped_connection_serves_no_book_and_a_new_one_subscribes_afresh(made):
    stream, connection, clock, _signer = made
    subscribe(stream, clock)
    stream._drop("ConnectionError: Kalshi closed the stream")
    assert view(stream)["error"] == "the Kalshi stream is not connected (ConnectionError: Kalshi closed the stream)"
    assert stream.status()["lastFrameAgeMs"] is None
    stream._open()
    status = stream.status()
    assert status["error"] is None, "a new connection clears the error"
    assert status["disconnects"] == 1 and status["lastDisconnect"] == "ConnectionError: Kalshi closed the stream"
    assert status["lastDisconnectAgoSec"] == 0
    assert "waiting for" in view(stream)["error"], "no book from the old connection"
    assert tick(stream, clock) == [{"id": 1, "cmd": "subscribe",
                                    "params": {"channels": ["orderbook_delta"], "market_tickers": [WSU]}}]


def test_a_ticket_ask_waits_for_a_new_markets_snapshot(made):
    stream, _connection, clock, _signer = made
    stream.books_for([WSU], 0.0)
    (frame,) = tick(stream, clock)

    def answer():
        with stream._cond:
            stream._handle(subscribed(frame["id"], 11))
            stream._handle(snapshot(11, 1))
            stream._cond.notify_all()

    timer = threading.Timer(0.05, answer)
    timer.start()
    # A real clock for the wait: the fake one never moves, so a missed notify would hang to the cap.
    stream._clock = __import__("time").monotonic
    got = stream.books_for([WSU], 2.0, priority=True)
    timer.join()
    assert "orderbook_fp" in got["books"][WSU]


def test_silence_reconnects(made):
    stream, connection, clock, _signer = made
    stop = threading.Event()
    real_recv = connection.recv

    def silent_recv(timeout_sec):
        clock.now += kbs.SILENT_RECONNECT_SEC + 1
        return real_recv(timeout_sec)

    connection.recv = silent_recv
    with pytest.raises(ConnectionError, match="no frame, ping or pong"):
        stream._pump(stop)


def test_run_reconnects_after_an_unreadable_frame_and_stops_cleanly(monkeypatch):
    clock = Clock()
    connections = []

    def connect(url, headers):
        connection = FakeConnection()
        connection.incoming = ["not json"]
        connections.append(connection)
        return connection

    stream = kbs.KalshiBookStream(FakeSigner(), URL, connect=connect, clock=clock)
    stop = threading.Event()
    waits = []

    def wait(sec):
        waits.append(sec)
        return len(waits) >= 2

    monkeypatch.setattr(stop, "wait", wait)
    stream.run(stop)
    assert len(connections) == 2 and all(connection.closed for connection in connections)
    assert waits == [kbs.RECONNECT_MIN_SEC, kbs.RECONNECT_MIN_SEC * 2], "the backoff doubles"
    assert "JSONDecodeError" in stream.status()["error"]


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


def test_post_kalshi_books_answers_the_streams_books(bets_store, made):
    stream, _connection, clock, _signer = made
    subscribe(stream, clock)
    server, base = _serve(service.make_handler(bets_store, 0.0, ["kalshi"], kalshi_book_stream=stream,
                                               extra_hosts=("mlb-stack.test",)))
    try:
        status, body = _post(f"{base}/kalshi_books.json", {"tickers": [WSU, WSU], "waitMs": 0, "priority": True})
        bad = _post(f"{base}/kalshi_books.json", {"tickers": [""]})
        foreign = _post(f"{base}/kalshi_books.json", {"tickers": []}, host="mlb-stack.test")
        with urllib.request.urlopen(f"{base}/health", timeout=5) as response:
            health = json.loads(response.read())
    finally:
        server.shutdown()
        server.server_close()
    assert status == 200 and body["books"] == {WSU: {"orderbook_fp": {"yes_dollars": [["0.5500", "100.00"]],
                                                                      "no_dollars": [["0.4300", "40.00"]]}, "ageMs": 0}}
    assert body["stream"]["connected"] is True
    assert bad[0] == 400 and "tickers[0]" in bad[1]["error"]
    assert foreign[0] == 403, "the runner's book route is loopback-only"
    assert health["kalshiStream"]["booksHeld"] == 1


def test_post_kalshi_books_is_503_without_the_read_key(bets_store):
    server, base = _serve(service.make_handler(bets_store, 0.0, ["kalshi"]))
    try:
        status, body = _post(f"{base}/kalshi_books.json", {"tickers": [], "waitMs": 0})
    finally:
        server.shutdown()
        server.server_close()
    assert status == 503 and "KALSHI_API_KEY_ID" in body["error"]


def test_the_stream_is_off_without_the_read_key_or_against_a_different_kalshi(monkeypatch):
    monkeypatch.setattr(service.config, "KALSHI_API_KEY_ID", "")
    assert service.open_kalshi_book_stream() is None
    monkeypatch.setattr(service.config, "KALSHI_API_KEY_ID", "d518d266-test")
    monkeypatch.setattr(service.config, "KALSHI_PRIVATE_KEY_PATH", "/data/kalshi.pem")
    monkeypatch.setattr(service.config, "KALSHI_BASE_URL", "https://demo-api.kalshi.co/trade-api/v2")
    assert service.open_kalshi_book_stream(client_factory=lambda *_args: object()) is None
    monkeypatch.setattr(service.config, "KALSHI_BASE_URL", "https://api.elections.kalshi.com/trade-api/v2")
    made_stream = service.open_kalshi_book_stream(client_factory=lambda *_args: FakeSigner())
    assert isinstance(made_stream, kbs.KalshiBookStream)


def test_the_disconnect_count_skips_failed_connects_and_stops_and_ages(made):
    stream, _connection, clock, _signer = made
    status = stream.status()
    assert (status["disconnects"], status["lastDisconnect"], status["lastDisconnectAgoSec"]) == (0, None, None)
    stream._drop("ConnectionError: one")
    stream._drop("OSError: could not connect")
    clock.now += 7
    stream._open()
    stream._drop("stopped")
    stream._open()
    stream._drop("ConnectionError: two")
    clock.now += 3
    status = stream.status()
    assert status["disconnects"] == 2 and status["lastDisconnect"] == "ConnectionError: two"
    assert status["lastDisconnectAgoSec"] == 3


def test_a_place_market_displaces_a_live_pass_market_at_the_cap(made, monkeypatch):
    stream, _connection, clock, _signer = made
    monkeypatch.setattr(kbs, "MARKETS_PER_SUBSCRIPTION", 1)
    subscribe(stream, clock, ticker="A-LIVE")
    with stream._cond:
        stream._background.add("A-LIVE")
    stream.books_for(["Z-PLACE"], 0.0, priority=True)
    frames = [tick(stream, clock), tick(stream, clock)]
    assert [(f["params"]["action"], f["params"]["market_tickers"]) for fs in frames for f in fs] == [
        ("delete_markets", ["A-LIVE"]), ("add_markets", ["Z-PLACE"])]


def test_the_shipped_sizes_hold_the_live_pass_and_send_fifty_a_command():
    assert kbs.MARKETS_PER_SUBSCRIPTION > kbs.MAX_BACKGROUND_MARKETS
    assert kbs.MARKETS_PER_COMMAND == 50


def test_a_subscribe_kalshi_keeps_refusing_reconnects_after_a_doubling_gap(made):
    stream, _connection, clock, _signer = made
    refused = {"type": "error", "msg": {"code": 8, "msg": "Unknown channel name"}}
    gaps = []
    for round_number in range(3):
        stream.books_for([f"M{round_number}"], 0.0)
        (frame,) = tick(stream, clock)
        handle(stream, {**refused, "id": frame["id"]})
        stream.books_for([f"N{round_number}"], 0.0)
        waited = 0.0
        while stream._reconnect_reason is None:
            tick(stream, clock, 0.25)
            waited += 0.25
        gaps.append(waited)
        reconnect(stream, "a new subscription needs a fresh connection")
    assert gaps[1] > gaps[0] * 1.5 and gaps[2] > gaps[1] * 1.5, gaps


def test_a_refused_delete_reconnects_and_a_rate_limited_one_pauses(made, monkeypatch):
    stream, _connection, clock, _signer = made
    monkeypatch.setattr(kbs, "RENEW_SEC", 3 * kbs.WATCH_IDLE_SEC)
    subscribe(stream, clock)
    add(stream, clock)
    clock.now += kbs.WATCH_IDLE_SEC - 1
    stream.books_for([WSU], 0.0)
    (delete,) = tick(stream, clock, 2)
    assert delete["params"]["action"] == "delete_markets"
    handle(stream, {"id": delete["id"], "type": "error", "msg": {"code": 27, "msg": "exceeded its command rate limit"}})
    assert "a refused delete" in stream._reconnect_reason
    assert stream._commands_paused_until > clock.now


def test_a_place_ask_at_a_full_cap_waits_for_its_seat(made, monkeypatch):
    stream, _connection, clock, _signer = made
    monkeypatch.setattr(kbs, "MARKETS_PER_SUBSCRIPTION", 1)
    stream.books_for(["A-LIVE"], 0.0, background=True)
    tick(stream, clock)
    with stream._cond:
        stream._wanted["Z-PLACE"] = clock.now
        stream._priority.add("Z-PLACE")
        assert stream._view_of("Z-PLACE").get("pending") is True
        stream._wanted["B-LIVE"] = clock.now
        stream._background.add("B-LIVE")
        assert "takes no more" in stream._view_of("B-LIVE")["error"], "a live pass market is not waited for"
