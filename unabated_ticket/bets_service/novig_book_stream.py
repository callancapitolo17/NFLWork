"""Novig's live order books over ONE websocket (GET /v3/ws, `book` channel),
signed with the bot's Novig trading key and used for reads only. It replaces
the runner's per-market polling of Novig's public book
(GET /v3/public/catalog/markets/<id>/book). That polling shared one unkeyed
per-IP throttle with the Place button's quote, so Place got HTTP 429 on
2026-10-27. The runner asks POST /novig_books.json (service.py) for the books
it needs, and the stream subscribes to those markets on demand.

Novig facts this rests on (docs.novig.com streaming/connection, streaming/book,
streaming/lifecycle and throttling pages, as copied 2026-09-28 into
github.com/Semionsr/novig-v3-kit fixtures/docs; unverified live until the
first run):
  * The upgrade is signed like any GET (path /v3/ws, empty query and body). A
    trading key may open it. A read accepts a stale phone geolocation, so the
    stream runs while Cal is outside the US.
  * Verbs are JSON frames {"nonce": n, "<verb>": payload}, with nonces
    increasing from 1 on each connection. `subscribe` {"markets": {id: "book"}}
    answers {"nonce", "subscribed", "snapshot": {marketId: {eventId, book:
    {seq, orders: {outcomeId: [{order, price, qty}]}}, lifecycle: {seq,
    status}}}}. A market it cannot resolve still costs tokens; here it is one
    the answer's snapshot leaves out.
  * Deltas: {"delta": {marketId: {eventId, book: {seq, deltas: [{kind: "add",
    order, outcome, price, qty} | {kind: "remove", order, reason}]},
    lifecycle: {seq, deltas: ["GOLIVE", ...]}}}}. seq runs per market per
    channel. A skipped seq is a gap: the deltas after it are buffered, a
    `snapshot` {"markets": {id: "book"}} is asked for, and the buffered deltas
    above the snapshot's seq are replayed on it. A batch is applied whole.
  * An error is a frame with `code` and `message` (plus `nonce` unless the
    frame was throttled or failed to parse).
  * The `stream` throttle: 512 tokens, refilling 4/s. The upgrade costs 32, a
    book subscription 16 per market, an unsubscribe 1 per subject. A frame
    costing more than 512 passes only on a full bucket and empties it, so no
    frame here costs more than 512. Subscribing markets, not events, keeps
    the cost 16 per market whichever way the docs' "an event counts as the
    markets it contains" is read, and spends nothing on markets nobody prices.
  * Novig pings every 15 s. A connection that misses two pongs is dropped
    without a close frame. GOLIVE (or START) arrives as a lifecycle delta;
    a snapshot's status stays OPEN in play, so `started` is best effort: the
    runner's own event OPEN_PREGAME read before every order and quote is the
    guard against a started game.

A book is served only while the connection is up and was heard from within
STALE_VIEW_SEC, and the market's seq ran without a gap since its snapshot.
Anything else answers an `error` naming why, never a stale book. The Place
button's markets go ahead of the bot's in the subscribe queue.

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
# Novig pings every 15 s: a book is not served after this long without a frame,
# and the connection is dropped after SILENT_RECONNECT_SEC.
STALE_VIEW_SEC = 25
SILENT_RECONNECT_SEC = 40
RECV_POLL_SEC = 0.25
CONNECT_TIMEOUT_SEC = 10
RECONNECT_MIN_SEC = 1
RECONNECT_MAX_SEC = 30
# Only a connection that held this long resets the backoff, so one Novig closes at once (a 451) is not hammered.
STABLE_CONNECTION_SEC = 60
# The same connect failure is logged again at most this often (Novig unreachable would log every 30 s).
REPEAT_FAILURE_LOG_SEC = 10 * 60
# A market the runner has not asked about for this long is unsubscribed (games end; the 2048-market cap).
WATCH_IDLE_SEC = 10 * 60
# A subscribe or snapshot with no answer this long after it was sent is sent again.
ANSWER_WAIT_SEC = 10
# A market Novig refused (or left out of a subscribe answer) is asked for again no sooner than this.
REFUSED_RETRY_SEC = 60
# Deltas buffered behind a gap, newest kept.
MAX_BUFFERED_DELTAS = 500
# The most a caller may wait for a snapshot (the runner's own timeouts sit above it).
MAX_WAIT_SEC = 3.0

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


# websocket-client is imported where it is used, so the tests' fake connections never need it.
class _WebsocketClientConnection:
    """websocket-client underneath. Its frame reader keeps a partial frame
    across a timeout, so short recv timeouts on one thread are safe."""

    def __init__(self, ws):
        self._ws = ws

    def send(self, text: str) -> None:
        import websocket
        self._ws.settimeout(CONNECT_TIMEOUT_SEC)
        try:
            self._ws.send(text)
        except websocket.WebSocketException as error:
            raise ConnectionError(f"{type(error).__name__}: {error}") from error

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
            raise ConnectionError(f"Novig closed the stream ({close_reason(frame.data)})")
        if opcode in (websocket.ABNF.OPCODE_TEXT, websocket.ABNF.OPCODE_BINARY):
            return frame.data.decode(errors="replace") if isinstance(frame.data, bytes) else str(frame.data)
        return ""

    def close(self) -> None:
        import websocket
        try:
            self._ws.close()
        except (OSError, websocket.WebSocketException) as error:  # closing a dead socket has nothing to report
            log.debug("novig stream: close: %s", error)


def close_reason(data: bytes | None) -> str:
    """'code 1008, SLOW_CONSUMER' from a close frame's payload (2-byte code, then the reason)."""
    data = data or b""
    if len(data) < 2:
        return "no close code"
    return f"code {int.from_bytes(data[:2], 'big')}, {data[2:].decode(errors='replace')!r}"


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
        """Spend `cost` when the bucket holds it; False (nothing spent) when not."""
        if self.available() < cost:
            return False
        self._tokens -= cost
        return True

    def spend(self, cost: float) -> None:
        """Spend `cost` whatever the bucket holds (it may go below zero): a cost already incurred."""
        self.available()
        self._tokens -= cost

    def empty(self) -> None:
        self.available()
        self._tokens = min(self._tokens, 0.0)


class _MarketBook:
    """One market's resting orders, as Novig's snapshot and deltas leave them.
    Out of sync after a gap until a snapshot loads it again."""

    def __init__(self):
        self.seq: int | None = None
        self.synced = False
        self.status: str | None = None
        self.lifecycle_seq: int | None = None
        self.started = False
        self.orders: dict[str, dict] = {}
        # (seq, deltas) received while out of sync, replayed on the next snapshot.
        self.buffered: list[tuple[int, list]] = []

    def load(self, seq: int, orders_by_outcome: dict) -> None:
        self.orders = {}
        for outcome_id, ladder in (orders_by_outcome or {}).items():
            for order in ladder or []:
                self.orders[str(order["order"])] = {"outcome": outcome_id, "price": order["price"], "qty": order["qty"]}
        self.seq = seq
        self.synced = True
        buffered, self.buffered = sorted(self.buffered, key=lambda item: item[0]), []
        for buffered_seq, deltas in buffered:
            if buffered_seq > self.seq:
                self.step(buffered_seq, deltas)

    def step(self, seq: int, deltas: list) -> bool:
        """Apply one batch in sequence. A stale seq is ignored; a skipped one
        puts the book out of sync and buffers the batch. True on a new gap."""
        if not self.synced:
            # The oldest go first: a snapshot past them never needs them, and one short of them leaves a gap the replay finds.
            self.buffered = [*self.buffered, (seq, deltas)][-MAX_BUFFERED_DELTAS:]
            return False
        if seq <= self.seq:
            return False
        if seq != self.seq + 1:
            self.synced = False
            self.buffered = [(seq, deltas)]
            return True
        self._apply(deltas)
        self.seq = seq
        return False

    def _apply(self, deltas: list) -> None:
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
    """The books of the markets the runner prices, kept live on one websocket.
    run(stop) is the thread's loop; books_for() is called from HTTP threads."""

    def __init__(self, client: NovigOrderClient, url: str, connect: Connect = websocket_client_connect,
                 clock: Callable[[], float] = time.monotonic):
        self._client = client
        self._url = url
        self._connect = connect
        self._clock = clock
        self._cond = threading.Condition()
        self._bucket = TokenBucket(STREAM_CAPACITY, STREAM_REFILL_PER_SEC, clock)
        # Kept across connections: the markets asked for (when last asked), Place's among them, and refusals.
        self._wanted: dict[str, float] = {}
        self._priority: set[str] = set()
        self._refused: dict[str, tuple[float, str]] = {}
        self._reset_connection_state()
        self._connection: Connection | None = None
        self._connected_at: float | None = None
        self._last_frame_at: float | None = None
        self._last_error: str | None = None
        self._last_logged_failure: tuple[str, float] | None = None

    def _reset_connection_state(self) -> None:
        self._nonce = 0
        self._subscribed: set[str] = set()
        # Markets Novig answered a subscribe for after we stopped waiting (or wanting): adopted, or unsubscribed.
        self._stray: set[str] = set()
        self._books: dict[str, _MarketBook] = {}
        # nonce -> (market ids, sent at) of a subscribe, and of a snapshot healing a gap.
        self._subscribes_in_flight: dict[int, tuple[list[str], float]] = {}
        self._snapshots_in_flight: dict[int, tuple[list[str], float]] = {}

    # ---- what the HTTP route calls ----

    def books_for(self, market_ids: list[str], wait_sec: float, priority: bool = False) -> dict:
        """{stream: status, books: {marketId: view}}. Marks each market wanted
        (subscribing it if new; `priority` puts it ahead of the rest, for the
        Place button) and waits up to wait_sec (capped at MAX_WAIT_SEC) for
        every market to be answered. A view is {status, started, seq, orders}
        or {error}; `seq` is for reading the logs, not for the runner."""
        wait_sec = max(0.0, min(wait_sec, MAX_WAIT_SEC))
        with self._cond:
            now = self._clock()
            for market_id in market_ids:
                self._wanted[market_id] = now
                if priority:
                    self._priority.add(market_id)
            self._cond.notify_all()
            deadline = now + wait_sec
            while True:
                views = {market_id: self._view_of(market_id) for market_id in market_ids}
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
            "subscribedMarkets": len(self._subscribed),
            "wantedMarkets": len(self._wanted),
        }

    def _view_of(self, market_id: str) -> dict:
        if self._connection is None:
            return {"error": f"the Novig stream is not connected ({self._last_error or 'connecting'})"}
        silent_sec = self._clock() - self._last_frame_at
        if silent_sec > STALE_VIEW_SEC:
            return {"error": f"the Novig stream has sent nothing for {round(silent_sec)} s"}
        refused = self._refused.get(market_id)
        if refused is not None:
            return {"error": f"Novig refused the stream subscription to market {market_id}: {refused[1]}"}
        book = self._books.get(market_id)
        if book is not None and book.synced:
            return {"status": book.status, "started": book.started, "seq": book.seq, "orders": book.ladders()}
        if book is not None:
            return {"pending": True, "error": f"Novig market {market_id}'s stream book is resyncing after a gap"}
        return {"pending": True, "error": f"waiting for Novig's stream snapshot of market {market_id}"}

    # ---- the stream thread ----

    def run(self, stop: threading.Event) -> None:
        """Connect, keep the wanted markets subscribed, apply every frame;
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
            if opened_at is not None and self._clock() - opened_at > STABLE_CONNECTION_SEC:
                backoff = RECONNECT_MIN_SEC
            if stop.wait(backoff):
                break
            backoff = min(backoff * 2, RECONNECT_MAX_SEC)
        self._drop("stopped")

    def _open(self) -> None:
        with self._cond:
            self._bucket.spend(UPGRADE_COST)
        headers = self._client.signed_headers("GET", STREAM_PATH)
        connection = self._connect(self._url, headers)
        with self._cond:
            self._reset_connection_state()
            self._connection = connection
            self._connected_at = self._last_frame_at = self._clock()
            self._last_error = None
            self._last_logged_failure = None
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
            return
        if reason == "stopped":
            return
        now = self._clock()
        last = self._last_logged_failure
        if last is None or last[0] != reason or now - last[1] > REPEAT_FAILURE_LOG_SEC:
            log.warning("novig stream: could not connect: %s", reason)
            self._last_logged_failure = (reason, now)

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
        markets, then subscribes (Place's first) and the snapshots healing
        gaps, oldest need first, at most 512 tokens a frame."""
        frames: list[dict] = []
        with self._cond:
            now = self._clock()
            self._expire(now)
            frames.extend(self._unsubscribe_idle(now))
            in_flight = {market_id for markets, _sent_at in self._subscribes_in_flight.values() for market_id in markets}
            due = [market_id for market_id in self._wanted
                   if market_id not in self._subscribed and market_id not in in_flight and market_id not in self._refused]
            due.sort(key=lambda market_id: (market_id not in self._priority, market_id))
            frames.extend(self._paced_frames("subscribe", due, self._subscribes_in_flight, now))
            healing = {market_id for markets, _sent_at in self._snapshots_in_flight.values() for market_id in markets}
            gaps = sorted(market_id for market_id, book in self._books.items()
                          if not book.synced and market_id not in healing and market_id not in self._refused)
            frames.extend(self._paced_frames("snapshot", gaps, self._snapshots_in_flight, now))
        return frames

    def _expire(self, now: float) -> None:
        for in_flight in (self._subscribes_in_flight, self._snapshots_in_flight):
            for nonce, (_markets, sent_at) in list(in_flight.items()):
                if now - sent_at > ANSWER_WAIT_SEC:
                    del in_flight[nonce]
        for market_id, (refused_at, _reason) in list(self._refused.items()):
            if now - refused_at > REFUSED_RETRY_SEC:
                del self._refused[market_id]

    def _unsubscribe_idle(self, now: float) -> list[dict]:
        idle = sorted(market_id for market_id, asked_at in self._wanted.items() if now - asked_at > WATCH_IDLE_SEC)
        in_flight = {market_id for markets, _sent_at in self._subscribes_in_flight.values() for market_id in markets}
        to_drop = sorted({market_id for market_id in idle if market_id in self._subscribed} | self._stray)
        frames = []
        if to_drop and self._bucket.take(len(to_drop) * UNSUBSCRIBE_COST_PER_SUBJECT):
            frames.append({"nonce": self._next_nonce(), "unsubscribe": [f"market:{market_id}" for market_id in to_drop]})
            for market_id in to_drop:
                self._subscribed.discard(market_id)
                self._books.pop(market_id, None)
            self._stray.clear()
        # A market still subscribed (the bucket was short) or with a subscribe in flight stays wanted until dropped.
        for market_id in idle:
            if market_id not in self._subscribed and market_id not in in_flight:
                del self._wanted[market_id]
                self._priority.discard(market_id)
        return frames

    def _paced_frames(self, verb: str, market_ids: list[str], in_flight: dict, now: float) -> list[dict]:
        frames = []
        while market_ids:
            # The bucket holds at most 512 tokens, so a frame costs at most 512 (32 markets).
            batch = market_ids[:int(self._bucket.available() // BOOK_PAIR_COST)]
            if not batch or not self._bucket.take(len(batch) * BOOK_PAIR_COST):
                break
            nonce = self._next_nonce()
            frames.append({"nonce": nonce, verb: {"markets": {market_id: BOOK_CHANNEL for market_id in batch}}})
            in_flight[nonce] = (batch, now)
            market_ids = market_ids[len(batch):]
        return frames

    # ---- frames from Novig (called with the lock held) ----

    def _handle(self, message: object) -> None:
        if not isinstance(message, dict):
            raise TypeError(f"expected a JSON object frame from Novig, got {str(message)[:200]}")
        nonce = message.get("nonce")
        if isinstance(message.get("code"), str):
            self._handle_error(nonce, message["code"], str(message.get("message", "")))
            return
        snapshot = message.get("snapshot") or {}
        subscribed = self._subscribes_in_flight.pop(nonce, None) if isinstance(nonce, int) else None
        healed = self._snapshots_in_flight.pop(nonce, None) if isinstance(nonce, int) else None
        if subscribed is not None:
            now = self._clock()
            for market_id in subscribed[0]:
                if market_id in snapshot:
                    self._subscribed.add(market_id)
                    self._books[market_id] = _MarketBook()
                else:
                    self._refused[market_id] = (now, "the subscribe answer has no snapshot of it (no such market?)")
        for market_id, body in snapshot.items():
            market_id = str(market_id)
            if market_id not in self._books:
                # An answer that came after its wait ran out: kept if still wanted, else unsubscribed.
                if market_id not in self._wanted:
                    self._stray.add(market_id)
                    continue
                self._subscribed.add(market_id)
                self._books[market_id] = _MarketBook()
            self._guarded(market_id, self._apply_snapshot, body)
        if healed is not None:
            for market_id in healed[0]:
                book = self._books.get(market_id)
                if book is not None and not book.synced and market_id not in snapshot:
                    self._refused[market_id] = (self._clock(), "the gap snapshot's answer has no snapshot of it")
        for market_id, body in (message.get("delta") or {}).items():
            market_id = str(market_id)
            if market_id in self._books:
                self._guarded(market_id, self._apply_delta, body)

    def _guarded(self, market_id: str, apply: Callable[[_MarketBook, dict], None], body: object) -> None:
        """One market's part of a frame; a part it cannot read puts only that market out of sync."""
        book = self._books[market_id]
        try:
            if not isinstance(body, dict):
                raise TypeError(f"expected an object, got {str(body)[:200]}")
            apply(book, body)
        except (KeyError, ValueError, TypeError) as error:
            log.warning("novig stream: market %s: unreadable frame part (%s: %s); resyncing it",
                        market_id, type(error).__name__, error)
            book.synced = False
            book.buffered = []

    def _handle_error(self, nonce: object, code: str, text: str) -> None:
        reason = f"{code} {text}".strip()
        for in_flight, what in ((self._subscribes_in_flight, "subscribe"), (self._snapshots_in_flight, "gap snapshot")):
            sent = in_flight.pop(nonce, None) if isinstance(nonce, int) else None
            if sent is None:
                continue
            now = self._clock()
            for market_id in sent[0]:
                self._refused[market_id] = (now, reason)
            log.error("novig stream: %s of %d market(s) refused: %s", what, len(sent[0]), reason)
            return
        if code == "RATE_LIMIT_EXCEEDED" or nonce is None:
            # A throttled frame carries no nonce: what it asked for is sent again after ANSWER_WAIT_SEC.
            self._bucket.empty()
        log.warning("novig stream: Novig answered %s (nonce %s)", reason, nonce)

    def _apply_snapshot(self, book: _MarketBook, body: dict) -> None:
        lifecycle = body.get("lifecycle")
        if isinstance(lifecycle, dict):
            book.lifecycle_seq = int(lifecycle["seq"])
            book.status = lifecycle.get("status")
        book_body = body.get("book")
        if isinstance(book_body, dict):
            book.load(int(book_body["seq"]), book_body.get("orders") or {})

    def _apply_delta(self, book: _MarketBook, body: dict) -> None:
        book_body = body.get("book")
        if isinstance(book_body, dict) and book.seq is not None:
            book.step(int(book_body["seq"]), book_body.get("deltas") or [])
        lifecycle = body.get("lifecycle")
        if isinstance(lifecycle, dict):
            self._apply_lifecycle(book, lifecycle)

    def _apply_lifecycle(self, book: _MarketBook, lifecycle: dict) -> None:
        transitions = [str(transition) for transition in lifecycle.get("deltas") or []]
        # A start or go-live counts even across a gap: it never un-happens.
        if any(transition in STARTED_TRANSITIONS for transition in transitions):
            book.started = True
        seq = int(lifecycle["seq"])
        if book.lifecycle_seq is not None and seq <= book.lifecycle_seq:
            return
        if book.lifecycle_seq is not None and seq > book.lifecycle_seq + 1:
            # A missed status change: the next snapshot carries the status.
            book.synced = False
        for transition in transitions:
            if transition in STATUS_AFTER_TRANSITION:
                book.status = STATUS_AFTER_TRANSITION[transition]
        book.lifecycle_seq = seq
