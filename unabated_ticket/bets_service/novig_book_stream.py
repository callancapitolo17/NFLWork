"""Novig's live order books over ONE websocket (GET /v3/ws, `book` channel),
signed with the bot's Novig trading key and used for reads only. It replaces
the runner's per-market polling of Novig's public book
(GET /v3/public/catalog/markets/<id>/book). That polling shared one unkeyed
per-IP throttle with the Place button's quote, so Place got HTTP 429 on
2026-10-27. The runner asks POST /novig_books.json (service.py) for the books
it needs, and the stream subscribes to those markets' events on demand.

Novig facts this rests on (docs.novig.com streaming/connection, streaming/book,
streaming/lifecycle and throttling pages, as copied 2026-09-28 into
github.com/Semionsr/novig-v3-kit fixtures/docs; unverified live until the
first run):
  * The upgrade is signed like any GET (path /v3/ws, empty query and body). A
    trading key may open it. A read accepts a stale phone geolocation, so the
    stream runs while Cal is outside the US.
  * Verbs are JSON frames {"nonce": n, "<verb>": payload}, with nonces
    increasing from 1 on each connection. `subscribe` {"events": {id: "book"}}
    answers with {"nonce", "subscribed", "snapshot": {marketId: {eventId,
    book: {seq, orders: {outcomeId: [{order, price, qty}]}}, lifecycle: {seq,
    status}}}}. An event subscription also covers markets that open later:
    such a market arrives with no snapshot, and its first book delta is seq 1.
  * Deltas: {"delta": {marketId: {eventId, book: {seq, deltas: [{kind: "add",
    order, outcome, price, qty} | {kind: "remove", order, reason}]},
    lifecycle: {seq, deltas: ["GOLIVE", ...]}}}}. seq runs per market per
    channel. A skipped seq is a gap, healed with `snapshot` {"markets": {id:
    "book"}}. A batch is applied whole.
  * An error is a frame with `code` and `message` (plus `nonce` unless the
    frame was throttled or failed to parse).
  * The `stream` throttle: 512 tokens, refilling 4/s. The upgrade costs 32, a
    book pair 16 (one event is one pair), and an unsubscribe 1 per subject. A
    frame costing more than 512 passes only on a full bucket and empties it,
    so no frame here costs more than 512. One connection watches at most
    2048 markets (SUBSCRIPTION_LIMIT_EXCEEDED).
  * Novig pings every 15 s. A connection that misses two pongs is dropped
    without a close frame.

A book is served only while the connection is up, has been heard from within
SILENT_RECONNECT_SEC, and the market's seq ran without a gap since its
snapshot. Anything else answers an `error` naming why, never a stale book.

Side effects: one outbound websocket to Novig on the configured key (reads
only: no verb here can place, cancel or move money). Nothing on disk; books
are held in memory only.
"""
import json
import logging
import threading
import time
from collections.abc import Callable
from typing import Protocol

from unabated_ticket.bets_service.novig_orders import NovigOrderClient

log = logging.getLogger(__name__)

STREAM_PATH = "/v3/ws"
BOOK_CHANNEL = "book"
STREAM_CAPACITY = 512
STREAM_REFILL_PER_SEC = 4.0
UPGRADE_COST = 32
BOOK_PAIR_COST = 16
UNSUBSCRIBE_COST_PER_SUBJECT = 1
# 32 events x 16 tokens = 512, the most one frame may cost.
MAX_EVENTS_PER_SUBSCRIBE = STREAM_CAPACITY // BOOK_PAIR_COST
# Novig pings every 15 s: 40 s of silence is a dead connection.
SILENT_RECONNECT_SEC = 40
RECV_POLL_SEC = 0.25
CONNECT_TIMEOUT_SEC = 10
RECONNECT_MIN_SEC = 1
RECONNECT_MAX_SEC = 30
STABLE_CONNECTION_SEC = 60
# An event the runner has not asked about for this long is unsubscribed (games end; the 2048-market cap).
WATCH_IDLE_SEC = 10 * 60
# A subscribe with no answer this long after it was sent goes back in the queue.
SUBSCRIBE_ANSWER_SEC = 10
# An event Novig refused to subscribe is retried no sooner than this.
REFUSED_RETRY_SEC = 60
# The most a caller may wait for a snapshot (the runner's own timeouts sit above it).
MAX_WAIT_SEC = 3.0

STATUS_OPEN = "OPEN"
# Lifecycle transitions and the market status each leaves (streaming/lifecycle).
STATUS_AFTER_TRANSITION = {"OPEN": "OPEN", "CLOSE": "CLOSED", "GRADE": "SETTLED"}
# Either one means the game is under way: resting orders were voided and fees may apply.
STARTED_TRANSITIONS = ("START", "GOLIVE")


class Connection(Protocol):
    """One open websocket. recv returns a text message, "" for a control frame
    (ping, pong: the line is alive, no data), or None when nothing arrived
    within the timeout. Raises OSError (or ConnectionError) once the
    connection is gone, a close frame included."""

    def send(self, text: str) -> None: ...

    def recv(self, timeout_sec: float) -> str | None: ...

    def close(self) -> None: ...


Connect = Callable[[str, dict], Connection]


class _WebsocketClientConnection:
    """websocket-client underneath. Its frame reader keeps a partial frame
    across a timeout, so short recv timeouts on one thread are safe."""

    def __init__(self, ws):
        self._ws = ws

    def send(self, text: str) -> None:
        self._ws.send(text)

    def recv(self, timeout_sec: float) -> str | None:
        import websocket
        self._ws.settimeout(timeout_sec)
        try:
            opcode, frame = self._ws.recv_data_frame(True)
        except websocket.WebSocketTimeoutException:
            return None
        except websocket.WebSocketException as error:
            raise ConnectionError(f"{type(error).__name__}: {error}") from error
        if opcode == websocket.ABNF.OPCODE_CLOSE:
            data = frame.data or b""
            code = int.from_bytes(data[:2], "big") if len(data) >= 2 else None
            raise ConnectionError(f"Novig closed the stream (code {code}, {data[2:].decode(errors='replace')!r})")
        if opcode in (websocket.ABNF.OPCODE_TEXT, websocket.ABNF.OPCODE_BINARY):
            return frame.data.decode(errors="replace") if isinstance(frame.data, bytes) else str(frame.data)
        return ""

    def close(self) -> None:
        import websocket
        try:
            self._ws.close()
        except (OSError, websocket.WebSocketException) as error:  # closing a dead socket has nothing to report
            log.debug("novig stream: close: %s", error)


def websocket_client_connect(url: str, headers: dict) -> Connection:
    """The real connection (websocket-client): answers Novig's pings itself."""
    import websocket
    try:
        ws = websocket.create_connection(url, header=[f"{name}: {value}" for name, value in headers.items()],
                                         timeout=CONNECT_TIMEOUT_SEC)
    except websocket.WebSocketException as error:
        raise ConnectionError(f"{type(error).__name__}: {error}") from error
    return _WebsocketClientConnection(ws)


def stream_url(api_url: str) -> str:
    """wss://host/v3/ws from the REST base (https://host)."""
    base = api_url.rstrip("/")
    if base.startswith("https://"):
        return "wss://" + base[len("https://"):] + STREAM_PATH
    if base.startswith("http://"):
        return "ws://" + base[len("http://"):] + STREAM_PATH
    raise ValueError(f"expected an http(s) Novig API URL, got {api_url!r}")


class TokenBucket:
    """The client-side model of Novig's `stream` throttle (throttling § Pace your requests)."""

    def __init__(self, capacity: float, refill_per_sec: float, clock: Callable[[], float]):
        self.capacity = capacity
        self.refill_per_sec = refill_per_sec
        self._clock = clock
        self._tokens = capacity
        self._at = clock()

    def available(self) -> float:
        now = self._clock()
        self._tokens = min(self.capacity, self._tokens + (now - self._at) * self.refill_per_sec)
        self._at = now
        return self._tokens

    def take(self, cost: float) -> bool:
        if self.available() < cost:
            return False
        self._tokens -= cost
        return True

    def empty(self) -> None:
        self.available()
        self._tokens = 0.0


class _MarketBook:
    """One market's resting orders, as Novig's snapshot and deltas leave them."""

    def __init__(self, event_id: str | None, seq: int):
        self.event_id = event_id
        self.seq = seq
        self.synced = True
        self.status: str | None = None
        self.lifecycle_seq: int | None = None
        self.started = False
        self.orders: dict[str, dict] = {}

    def load(self, seq: int, orders_by_outcome: dict) -> None:
        self.seq = seq
        self.synced = True
        self.orders = {}
        for outcome_id, ladder in (orders_by_outcome or {}).items():
            for order in ladder or []:
                self.orders[str(order["order"])] = {"outcome": outcome_id, "price": order["price"], "qty": order["qty"]}

    def apply(self, deltas: list) -> None:
        for delta in deltas or []:
            if delta.get("kind") == "add":
                self.orders[str(delta["order"])] = {"outcome": delta["outcome"], "price": delta["price"], "qty": delta["qty"]}
            elif delta.get("kind") == "remove":
                self.orders.pop(str(delta.get("order")), None)
            else:
                raise ValueError(f"expected a book delta of kind add or remove, got {delta!r}")

    def ladders(self) -> dict:
        """{outcomeId: [{order, price, qty}]}: the REST book's `orders` shape (taker.novigQuoteOf)."""
        by_outcome: dict[str, list] = {}
        for order_id, order in self.orders.items():
            by_outcome.setdefault(order["outcome"], []).append({"order": order_id, "price": order["price"], "qty": order["qty"]})
        return by_outcome


class NovigBookStream:
    """The books of the events the runner watches, kept live on one websocket.
    run(stop) is the thread's loop; books_for() is called from HTTP threads."""

    def __init__(self, client: NovigOrderClient, url: str, connect: Connect = websocket_client_connect,
                 clock: Callable[[], float] = time.monotonic):
        self._client = client
        self._url = url
        self._connect = connect
        self._clock = clock
        self._cond = threading.Condition()
        self._bucket = TokenBucket(STREAM_CAPACITY, STREAM_REFILL_PER_SEC, clock)
        # Kept across connections: which events the runner wants, and when it last asked.
        self._wanted: dict[str, float] = {}
        self._refused: dict[str, tuple[float, str]] = {}
        self._reset_connection_state()
        self._connection: Connection | None = None
        self._connected_at: float | None = None
        self._last_frame_at: float | None = None
        self._last_error: str | None = None

    def _reset_connection_state(self) -> None:
        self._nonce = 0
        self._subscribed: set[str] = set()
        self._in_flight: dict[int, tuple[list[str], float]] = {}
        self._snapshotted: set[str] = set()
        self._books: dict[str, _MarketBook] = {}
        self._resync: set[str] = set()
        # nonce -> (market id, sent at) of a snapshot healing a gap.
        self._resync_in_flight: dict[int, tuple[str, float]] = {}

    # ---- what the HTTP route calls ----

    def books_for(self, markets: list[dict], wait_sec: float) -> dict:
        """{stream: status, books: {marketId: view}} for [{marketId, eventId}].
        Marks each event wanted (subscribing it if new) and waits up to
        wait_sec (capped at MAX_WAIT_SEC) for every market to be answered. A
        view is {eventId, status, started, seq, orders} or {error}."""
        wait_sec = max(0.0, min(wait_sec, MAX_WAIT_SEC))
        with self._cond:
            now = self._clock()
            for market in markets:
                self._wanted[market["eventId"]] = now
            self._cond.notify_all()
            deadline = now + wait_sec
            while True:
                views = {market["marketId"]: self._view_of(market["marketId"], market["eventId"]) for market in markets}
                remaining = deadline - self._clock()
                if remaining <= 0 or not any(view.get("pending") for view in views.values()):
                    break
                self._cond.wait(remaining)
            for view in views.values():
                view.pop("pending", None)
            return {"stream": self._status_locked(), "books": views}

    def status(self) -> dict:
        with self._cond:
            return self._status_locked()

    def _status_locked(self) -> dict:
        now = self._clock()
        return {
            "connected": self._connection is not None,
            "connectedForSec": None if self._connected_at is None else round(now - self._connected_at),
            "lastFrameAgeMs": None if self._last_frame_at is None else int((now - self._last_frame_at) * 1000),
            "error": self._last_error,
            "watchedEvents": len(self._subscribed),
            "markets": len(self._books),
        }

    def _view_of(self, market_id: str, event_id: str) -> dict:
        if self._connection is None:
            return {"error": f"the Novig stream is not connected ({self._last_error or 'connecting'})"}
        refused = self._refused.get(event_id)
        if refused is not None:
            return {"error": f"Novig refused the stream subscription to event {event_id}: {refused[1]}"}
        book = self._books.get(market_id)
        if book is not None and book.synced:
            return {"eventId": book.event_id, "status": book.status, "started": book.started, "seq": book.seq,
                    "orders": book.ladders()}
        if book is not None:
            return {"pending": True, "error": f"Novig market {market_id}'s stream book is resyncing after a gap"}
        if event_id in self._snapshotted:
            return {"error": f"Novig's stream lists no market {market_id} under event {event_id}"}
        return {"pending": True, "error": f"waiting for Novig's stream snapshot of event {event_id}"}

    # ---- the stream thread ----

    def run(self, stop: threading.Event) -> None:
        """Connect, keep the wanted events subscribed, apply every frame;
        reconnect with backoff whenever the line drops or falls silent."""
        backoff = RECONNECT_MIN_SEC
        while not stop.is_set():
            opened_at = None
            try:
                self._open()
                opened_at = self._clock()
                self._pump(stop)
            except (OSError, ConnectionError, ValueError, KeyError, TypeError) as error:
                self._drop(f"{type(error).__name__}: {error}")
            # Anything else too: a dead thread would leave the bot blind until a restart.
            except Exception as error:
                log.exception("novig stream: unexpected error")
                self._drop(f"unexpected {type(error).__name__}: {error}")
            # Only a connection that held a while resets the backoff: one Novig closes at once (a 451) is not hammered.
            if opened_at is not None and self._clock() - opened_at > STABLE_CONNECTION_SEC:
                backoff = RECONNECT_MIN_SEC
            if stop.wait(backoff):
                break
            backoff = min(backoff * 2, RECONNECT_MAX_SEC)
        self._drop("stopped")

    def _open(self) -> None:
        with self._cond:
            self._bucket.take(UPGRADE_COST)
        headers = self._client.signed_headers("GET", STREAM_PATH)
        connection = self._connect(self._url, headers)
        with self._cond:
            self._reset_connection_state()
            self._connection = connection
            self._connected_at = self._last_frame_at = self._clock()
            self._last_error = None
            self._cond.notify_all()
        log.info("novig stream: connected to %s on key %s", self._url, self._client.key_id[:8])

    def _drop(self, reason: str) -> None:
        with self._cond:
            connection = self._connection
            self._connection = None
            self._connected_at = None
            self._reset_connection_state()
            if reason != "stopped":
                self._last_error = reason
            self._cond.notify_all()
        if connection is not None:
            connection.close()
            log.warning("novig stream: disconnected: %s", reason)
        elif reason != "stopped":
            log.warning("novig stream: could not connect: %s", reason)

    def _pump(self, stop: threading.Event) -> None:
        while not stop.is_set():
            for frame in self._outgoing_frames():
                self._connection.send(json.dumps(frame, separators=(",", ":")))
            text = self._connection.recv(RECV_POLL_SEC)
            now = self._clock()
            if text is None:
                if now - self._last_frame_at > SILENT_RECONNECT_SEC:
                    raise ConnectionError(f"no frame from Novig for {SILENT_RECONNECT_SEC} s")
                continue
            with self._cond:
                self._last_frame_at = now
                if text:
                    self._handle(json.loads(text))
                self._cond.notify_all()

    def _next_nonce(self) -> int:
        self._nonce += 1
        return self._nonce

    def _outgoing_frames(self) -> list[dict]:
        """The verbs due now, each within the bucket: unsubscribes of idle
        events, snapshots healing gaps, subscribes of newly wanted events."""
        frames: list[dict] = []
        with self._cond:
            now = self._clock()
            for in_flight in (self._in_flight, self._resync_in_flight):
                for nonce, (_subject, sent_at) in list(in_flight.items()):
                    if now - sent_at > SUBSCRIBE_ANSWER_SEC:
                        del in_flight[nonce]
            for event_id, (refused_at, _reason) in list(self._refused.items()):
                if now - refused_at > REFUSED_RETRY_SEC:
                    del self._refused[event_id]
            idle = sorted(event_id for event_id, asked_at in self._wanted.items() if now - asked_at > WATCH_IDLE_SEC)
            to_drop = [event_id for event_id in idle if event_id in self._subscribed]
            if to_drop and self._bucket.take(len(to_drop) * UNSUBSCRIBE_COST_PER_SUBJECT):
                frames.append({"nonce": self._next_nonce(), "unsubscribe": [f"event:{event_id}" for event_id in to_drop]})
                self._forget_events(set(to_drop))
            for event_id in idle:
                if event_id not in self._subscribed:
                    del self._wanted[event_id]
            resyncing = {market_id for market_id, _sent_at in self._resync_in_flight.values()}
            for market_id in sorted(self._resync - resyncing):
                if not self._bucket.take(BOOK_PAIR_COST):
                    break
                nonce = self._next_nonce()
                frames.append({"nonce": nonce, "snapshot": {"markets": {market_id: BOOK_CHANNEL}}})
                self._resync_in_flight[nonce] = (market_id, now)
            in_flight = {event_id for events, _sent_at in self._in_flight.values() for event_id in events}
            due = sorted(event_id for event_id in self._wanted
                         if event_id not in self._subscribed and event_id not in in_flight and event_id not in self._refused)
            while due:
                affordable = int(self._bucket.available() // BOOK_PAIR_COST)
                batch = due[:min(affordable, MAX_EVENTS_PER_SUBSCRIBE)]
                if not batch or not self._bucket.take(len(batch) * BOOK_PAIR_COST):
                    break
                nonce = self._next_nonce()
                frames.append({"nonce": nonce, "subscribe": {"events": {event_id: BOOK_CHANNEL for event_id in batch}}})
                self._in_flight[nonce] = (batch, now)
                due = due[len(batch):]
        return frames

    def _forget_events(self, event_ids: set[str]) -> None:
        self._subscribed -= event_ids
        self._snapshotted -= event_ids
        for market_id in [market_id for market_id, book in self._books.items() if book.event_id in event_ids]:
            del self._books[market_id]
            self._resync.discard(market_id)

    # ---- frames from Novig (called with the lock held) ----

    def _handle(self, message: object) -> None:
        if not isinstance(message, dict):
            raise TypeError(f"expected a JSON object frame from Novig, got {str(message)[:200]}")
        nonce = message.get("nonce")
        if isinstance(message.get("code"), str):
            self._handle_error(nonce, message["code"], str(message.get("message", "")))
            return
        sent = self._in_flight.pop(nonce, None) if isinstance(nonce, int) else None
        if sent is not None:
            self._subscribed.update(sent[0])
            self._snapshotted.update(sent[0])
        if isinstance(nonce, int):
            self._resync_in_flight.pop(nonce, None)
        subscribed = message.get("subscribed")
        if isinstance(subscribed, dict) and isinstance(subscribed.get("events"), dict):
            self._subscribed.update(str(event_id) for event_id in subscribed["events"])
        for market_id, body in (message.get("snapshot") or {}).items():
            self._apply_snapshot(str(market_id), body)
        for market_id, body in (message.get("delta") or {}).items():
            self._apply_delta(str(market_id), body)

    def _handle_error(self, nonce: object, code: str, text: str) -> None:
        if isinstance(nonce, int):
            # A refused gap snapshot is sent again on the next turn of the loop.
            self._resync_in_flight.pop(nonce, None)
        sent = self._in_flight.pop(nonce, None) if isinstance(nonce, int) else None
        if sent is not None:
            now = self._clock()
            for event_id in sent[0]:
                self._refused[event_id] = (now, f"{code} {text}".strip())
            log.error("novig stream: subscribe of %d event(s) refused: %s %s", len(sent[0]), code, text)
            return
        if code == "RATE_LIMIT_EXCEEDED" or nonce is None:
            # A throttled frame carries no nonce: the subscribe it was goes back in the queue after SUBSCRIBE_ANSWER_SEC.
            self._bucket.empty()
        log.warning("novig stream: Novig answered %s %s (nonce %s)", code, text, nonce)

    def _apply_snapshot(self, market_id: str, body: dict) -> None:
        event_id = body.get("eventId")
        if event_id is not None:
            self._snapshotted.add(str(event_id))
        book_body = body.get("book")
        book = self._books.get(market_id)
        if isinstance(book_body, dict):
            if book is None:
                book = self._books[market_id] = _MarketBook(event_id, int(book_body["seq"]))
            book.event_id = event_id or book.event_id
            book.load(int(book_body["seq"]), book_body.get("orders") or {})
            self._resync.discard(market_id)
        lifecycle = body.get("lifecycle")
        if book is not None and isinstance(lifecycle, dict):
            book.lifecycle_seq = int(lifecycle["seq"])
            book.status = lifecycle.get("status")

    def _apply_delta(self, market_id: str, body: dict) -> None:
        book = self._books.get(market_id)
        book_body = body.get("book")
        if book is None and str(body.get("eventId")) not in self._subscribed:
            # A late delta of an event just unsubscribed.
            return
        if book is None:
            # A market opened under a watched event: no snapshot, its book starts at seq 0.
            if not isinstance(book_body, dict) or int(book_body["seq"]) != 1:
                self._books[market_id] = book = _MarketBook(body.get("eventId"), 0)
                book.synced = False
                self._resync.add(market_id)
                return
            book = self._books[market_id] = _MarketBook(body.get("eventId"), 0)
            book.lifecycle_seq = 0
        if isinstance(book_body, dict) and book.synced:
            seq = int(book_body["seq"])
            if seq <= book.seq:
                pass
            elif seq != book.seq + 1:
                log.info("novig stream: gap on market %s's book (had seq %d, got %d); resyncing", market_id, book.seq, seq)
                book.synced = False
                self._resync.add(market_id)
            else:
                book.apply(book_body.get("deltas") or [])
                book.seq = seq
        lifecycle = body.get("lifecycle")
        if isinstance(lifecycle, dict):
            self._apply_lifecycle(market_id, book, lifecycle)

    def _apply_lifecycle(self, market_id: str, book: _MarketBook, lifecycle: dict) -> None:
        transitions = [str(transition) for transition in lifecycle.get("deltas") or []]
        # A start or go-live counts even across a gap: it never un-happens.
        if any(transition in STARTED_TRANSITIONS for transition in transitions):
            book.started = True
        seq = int(lifecycle["seq"])
        if book.lifecycle_seq is not None and seq > book.lifecycle_seq + 1:
            book.synced = False
            self._resync.add(market_id)
        for transition in transitions:
            if transition in STATUS_AFTER_TRANSITION:
                book.status = STATUS_AFTER_TRANSITION[transition]
        book.lifecycle_seq = seq
