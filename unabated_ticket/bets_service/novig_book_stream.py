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

The same socket can carry the key's PRIVATE `orders` channel (`orders_channel`;
docs.novig.com/api/streaming/private, read 2026-11-09, exact frames below):
  * `subscribe` {"private": ["orders"]} costs 1 token; the answer carries
    `subscribed.private` (a list of channel names) and a TOP-LEVEL `orders`
    object, never under `snapshot`: {"seq", "open": [{orderId, clientId,
    marketId, outcomeId, price, qty, tif}]} (open orders only, so a finished
    IOC is never in it). `seq` is per subaccount per channel and resets on a
    Novig restart; 0 means no event yet.
  * Deltas: a frame with a top-level `orders` {"seq", "deltas": [...]} (its
    `delta` is {}), one seq per batch; an order's `fill` entries come before
    its `open` / `cancel` / `reject`. Entries: `open` {orderId, clientId?,
    marketId, outcomeId, price, qty (resting), originalQty, tif, expiresAt?};
    `fill` {fillId, orderId, clientId?, outcomeId, price, qty (this fill),
    originalQty, remaining, taker} — NO marketId and NO fee; `remaining` 0 is
    FILLED; `cancel` {orderId, clientId?, outcomeId, originalQty, remaining,
    reason? (GO_LIVE | MARKET_CLOSED | SETTLED | NEUTRALIZED)}; `reject`
    {orderId, clientId?, outcomeId, originalQty}. An IOC that fills nothing
    arrives as `reject`; one that part-fills as `fill` then `cancel`.
  * A `heartbeat` {"orders": <seq>, "ts"} every 15 s names the last seq on
    each subscribed private channel: one above ours is a dropped frame (Novig
    drops private frames on a full queue). Recovery: `snapshot` {"private":
    ["orders"]}, answered like the subscribe (no example in the docs; the
    kit's mock answers {"snapshot": {}, "orders": {seq, open}}).
  * The channel covers only the subaccount of the key that signed the
    upgrade, so the Place sub-account needs a socket of its own (service.py
    opens one on the Place key, markets never asked of it).
  Fills are pushed to the subaccount's bet source (fill_events_after) so its
  10 s REST re-read of two weeks of fills stops (2026-11-07: those reads, two
  pages of 500 at 8 + 1 per 50 rows out of a 4/s `history` bucket, were what
  Novig 429'd), and to the Place desk (wait_for_order), which no longer polls
  the order's status. A snapshot load bumps `ordersEpoch`: the source reads
  REST once whenever it sees a new epoch, since a fill during the gap or the
  disconnect is in no snapshot.

A book is served only while the connection is up and was heard from within
STALE_VIEW_SEC, and the market's seq ran without a gap since its snapshot.
Anything else answers an `error` naming why, never a stale book. The Place
button's markets go ahead of the bot's in the subscribe queue.

Side effects: one outbound websocket to Novig on the configured key (reads
only: no verb here can place, cancel or move money). Nothing on disk; books
and order states are held in memory only.
"""
import json
import logging
import threading
import time
from collections.abc import Callable
from typing import Protocol

from unabated_ticket.bets_service.novig_orders import QTY_PER_DOLLAR, NovigOrderClient
from unabated_ticket.bets_service.novig_orders import STATUS_CANCELED as STATUS_CANCELED_ORDER
from unabated_ticket.bets_service.novig_orders import STATUS_FILLED as STATUS_FILLED_ORDER
from unabated_ticket.bets_service.novig_orders import STATUS_OPEN as STATUS_OPEN_ORDER
from unabated_ticket.bets_service.novig_orders import STATUS_REJECTED as STATUS_REJECTED_ORDER

log = logging.getLogger(__name__)

STREAM_PATH = "/v3/ws"
BOOK_CHANNEL = "book"
STREAM_CAPACITY = 512
STREAM_REFILL_PER_SEC = 4.0
UPGRADE_COST = 32
BOOK_PAIR_COST = 16
UNSUBSCRIBE_COST_PER_SUBJECT = 1
# The private `orders` channel: 1 token to subscribe, 1 for a recovery snapshot.
ORDERS_CHANNEL = "orders"
PRIVATE_CHANNEL_COST = 1
# A pushed order state is kept this long after its last event (Place waits at most
# PLACE_FINAL_READ_WINDOW_SEC; the source drains fill events every poll).
ORDER_STATE_KEEP_SEC = 10 * 60
# Fill events kept for the source to drain (a poll is 10 s; this is many hours of fills).
MAX_FILL_EVENTS = 5000
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
# Novig caps a connection at 2048 subscribed markets. The runner's live price
# pass (background asks) may hold at most this many, so the bot's and Place's
# markets always have room; a background market past it answers an error.
MAX_BACKGROUND_MARKETS = 1536

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

    def __init__(self, ws, peer: str = "Novig"):
        self._ws = ws
        self._peer = peer

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
            raise ConnectionError(f"{self._peer} closed the stream ({close_reason(frame.data)})")
        if opcode in (websocket.ABNF.OPCODE_TEXT, websocket.ABNF.OPCODE_BINARY):
            return frame.data.decode(errors="replace") if isinstance(frame.data, bytes) else str(frame.data)
        return ""

    def ping(self) -> None:
        """A websocket-level ping; the pong comes back through recv as "" (the Polymarket US
        stream sends nothing while its books are still, so a ping is its only sign of life)."""
        import websocket
        self._ws.settimeout(CONNECT_TIMEOUT_SEC)
        try:
            self._ws.ping()
        except websocket.WebSocketException as error:
            raise ConnectionError(f"{type(error).__name__}: {error}") from error

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


def websocket_client_connect(url: str, headers: dict, peer: str = "Novig") -> Connection:
    """The real connection (websocket-client): answers websocket-level pings itself.
    `peer` names the server in a close error (the ProphetX order stream shares this)."""
    import websocket
    try:
        ws = websocket.create_connection(url, header=[f"{name}: {value}" for name, value in headers.items()],
                                         timeout=CONNECT_TIMEOUT_SEC)
    except websocket.WebSocketException as error:
        raise ConnectionError(f"{type(error).__name__}: {error}") from error
    return _WebsocketClientConnection(ws, peer)


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


# Novig's streaming docs name an order's id `order` and an add's outcome `outcome`; the live
# frames (2026-10-28) carry `orderId` and `outcomeId`, the REST book's names. Either is read.
ORDER_ID_KEYS = ("order", "orderId")
OUTCOME_ID_KEYS = ("outcome", "outcomeId")
LOGGED_PART_CHARS = 400


def _id_under(entry: dict, keys: tuple[str, ...], what: str) -> str:
    if not isinstance(entry, dict):
        raise TypeError(f"expected an order object, got {str(entry)[:200]}")
    for key in keys:
        if entry.get(key) is not None:
            return str(entry[key])
    raise KeyError(f"expected {what} under {' or '.join(keys)}, got keys {sorted(entry)}")


def order_id_of(entry: dict) -> str:
    """An order's id in a snapshot ladder entry or a book delta. A missing id raises:
    a remove that cannot name its order would leave a filled order showing as liquidity."""
    return _id_under(entry, ORDER_ID_KEYS, "an order id")


def outcome_id_of(delta: dict) -> str:
    """The outcome an add delta rests on (snapshot ladders are keyed by outcome instead)."""
    return _id_under(delta, OUTCOME_ID_KEYS, "an outcome id")


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
        # The last load's replay found a seq missing: a new gap, not a bad answer.
        self.replay_gap = False

    def load(self, seq: int, orders_by_outcome: dict) -> None:
        self.orders = {}
        for outcome_id, ladder in (orders_by_outcome or {}).items():
            for order in ladder or []:
                self.orders[order_id_of(order)] = {"outcome": outcome_id, "price": order["price"], "qty": order["qty"]}
        self.seq = seq
        self.synced = True
        self.replay_gap = False
        buffered, self.buffered = sorted(self.buffered, key=lambda item: item[0]), []
        for buffered_seq, deltas in buffered:
            if buffered_seq > self.seq and self.step(buffered_seq, deltas):
                self.replay_gap = True

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
                self.orders[order_id_of(delta)] = {"outcome": outcome_id_of(delta), "price": delta["price"], "qty": delta["qty"]}
            elif delta.get("kind") == "remove":
                self.orders.pop(order_id_of(delta), None)
            elif delta.get("kind") == "update":
                self._update_remaining(order_id_of(delta), delta["remaining"])
            else:
                raise ValueError(f"expected a book delta of kind add, remove or update, got {delta!r}")

    def _update_remaining(self, order_id: str, remaining: object) -> None:
        """A partial fill: Novig's `update` carries only the order and what is left of it (live
        frames, 2026-10-28; not in the streaming docs). An order this book never saw means the
        book is wrong, so it resyncs rather than guess."""
        if order_id not in self.orders:
            raise KeyError(f"update for order {order_id} this book does not hold")
        if not isinstance(remaining, int) or isinstance(remaining, bool) or remaining < 0:
            raise ValueError(f"expected a whole remaining qty >= 0, got {remaining!r}")
        if remaining == 0:
            del self.orders[order_id]
        else:
            self.orders[order_id]["qty"] = remaining

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
                 clock: Callable[[], float] = time.monotonic, orders_channel: bool = False,
                 wall_clock_ms: Callable[[], int] = lambda: int(time.time() * 1000)):
        """`orders_channel`: also subscribe the key's private `orders` channel
        (the subaccount's own fills, cancels and rejects, pushed)."""
        self._client = client
        self._url = url
        self._connect = connect
        self._clock = clock
        self._wall_clock_ms = wall_clock_ms
        self._orders_channel = orders_channel
        self._cond = threading.Condition()
        self._bucket = TokenBucket(STREAM_CAPACITY, STREAM_REFILL_PER_SEC, clock)
        # Kept across connections: every order the private channel told of, by Novig order id, and the
        # fill events the source drains, each tagged with a version that only grows.
        self._order_states: dict[str, dict] = {}
        self._fill_events: list[tuple[int, dict]] = []
        self._orders_version = 0
        # Bumped on every private snapshot load (connect, gap): a new epoch means frames may have been missed.
        self._orders_epoch = 0
        self._orders_error: str | None = None
        # Kept across connections: the markets asked for (when last asked), Place's among them, and refusals.
        self._wanted: dict[str, float] = {}
        self._priority: set[str] = set()
        # Markets only the runner's live price pass asks for: subscribed after every other.
        self._background: set[str] = set()
        self._refused: dict[str, tuple[float, str]] = {}
        # Markets of a batch Novig refused as a whole: sent one a frame, so one bad market refuses only itself.
        self._solo: set[str] = set()
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
        # Every subscribe nonce sent on this connection: only a subscribe's answer subscribes a market.
        self._subscribe_nonces: set[int] = set()
        # The private orders channel on this connection: subscribed, in sync (its snapshot loaded and
        # no seq skipped), its last seq, the nonce + time of the subscribe or snapshot awaiting an
        # answer, deltas buffered behind a gap, and when a refused subscribe may be sent again.
        self._orders_subscribed = False
        self._orders_synced = False
        self._orders_seq: int | None = None
        self._orders_in_flight: tuple[int, float] | None = None
        self._orders_buffered: list[tuple[int, list]] = []
        self._orders_refused_at: float | None = None

    # ---- what the HTTP route calls ----

    def books_for(self, market_ids: list[str], wait_sec: float, priority: bool = False, background: bool = False) -> dict:
        """{stream: status, books: {marketId: view}}. Marks each market wanted
        (subscribing it if new; `priority` puts it ahead of the rest, for the
        Place button; `background` behind the rest, for the runner's live price
        pass, until the bot or Place asks for it too) and waits up to wait_sec (capped at MAX_WAIT_SEC) for
        every market to be answered. A view is {status, started, seq, orders}
        or {error}; `seq` is for reading the logs, not for the runner."""
        wait_sec = max(0.0, min(wait_sec, MAX_WAIT_SEC))
        with self._cond:
            now = self._clock()
            capped: set[str] = set()
            for market_id in market_ids:
                if background and market_id not in self._wanted:
                    if len(self._background) >= MAX_BACKGROUND_MARKETS:
                        capped.add(market_id)
                        continue
                    self._background.add(market_id)
                elif not background:
                    self._background.discard(market_id)
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
            for market_id in capped:
                views[market_id] = {"error": f"the live price pass already watches {MAX_BACKGROUND_MARKETS} Novig markets"}
            return {"stream": self._status_locked(), "books": views}

    def status(self) -> dict:
        with self._cond:
            return self._status_locked()

    def _status_locked(self) -> dict:
        now = self._clock()
        status = {
            "connected": self._connection is not None,
            "connectedForSec": None if self._connected_at is None else round(now - self._connected_at),
            "lastFrameAgeMs": None if self._last_frame_at is None else int((now - self._last_frame_at) * 1000),
            "error": self._last_error,
            "subscribedMarkets": len(self._subscribed),
            "wantedMarkets": len(self._wanted),
        }
        if self._orders_channel:
            status["orders"] = {"live": self._orders_live_locked(), "subscribed": self._orders_subscribed,
                                "synced": self._orders_synced, "seq": self._orders_seq, "epoch": self._orders_epoch,
                                "fillEvents": self._orders_version, "error": self._orders_error}
        return status

    # ---- what the subaccount source and the Place desk read (orders channel) ----

    def orders_live(self) -> bool:
        """Is the private orders channel delivering now: connected, heard from
        within STALE_VIEW_SEC, subscribed and in sync (no seq skipped since its snapshot)?"""
        with self._cond:
            return self._orders_live_locked()

    def _orders_live_locked(self) -> bool:
        return (self._connection is not None and self._orders_subscribed and self._orders_synced
                and self._last_frame_at is not None and self._clock() - self._last_frame_at <= STALE_VIEW_SEC)

    def fill_events_after(self, version: int) -> dict:
        """{events, version, live, epoch}: the fill events pushed since
        `version` (oldest first, each {fillId, orderId, clientId, outcomeId,
        price, qty, originalQty, remaining, taker, ts}), the version to pass
        next time, whether the channel is live now, and the snapshot epoch. A
        caller that sees a new epoch, or `live` False, reads Novig's REST fills
        once instead: a fill during a gap or a disconnect is in no snapshot."""
        with self._cond:
            events = [event for event_version, event in self._fill_events if event_version > version]
            return {"events": events, "version": self._orders_version, "live": self._orders_live_locked(),
                    "epoch": self._orders_epoch}

    def wait_for_order(self, order_id: str, after_version: int, timeout_sec: float) -> tuple[dict | None, int]:
        """The order's pushed state once it is newer than `after_version` (or
        None at the timeout, or at once when the channel is not live), and the
        version to pass next time. A state is {final, filled, cost, remaining,
        fills}: `final` FILLED | CANCELED | REJECTED once Novig is done with it,
        OPEN while it rests, None while only fills have arrived."""
        deadline = self._clock() + timeout_sec
        with self._cond:
            while True:
                held = self._order_states.get(order_id)
                if held is not None and held["version"] > after_version:
                    return {key: held[key] for key in ("final", "filled", "cost", "remaining", "fills")}, held["version"]
                remaining = deadline - self._clock()
                if remaining <= 0 or not self._orders_live_locked():
                    return None, after_version
                self._cond.wait(min(remaining, RECV_POLL_SEC))

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
        log.info("novig stream: connected to %s on key %s%s", self._url, self._client.key_id[:8],
                 " (private orders channel wanted)" if self._orders_channel else "")

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
            frames.extend(self._orders_channel_frames(now))
            frames.extend(self._unsubscribe_idle(now))
            in_flight = {market_id for markets, _sent_at in self._subscribes_in_flight.values() for market_id in markets}
            due = [market_id for market_id in self._wanted
                   if market_id not in self._subscribed and market_id not in in_flight and market_id not in self._refused]
            due.sort(key=lambda market_id: (market_id not in self._priority, market_id in self._background, market_id))
            for market_id in [market_id for market_id in due if market_id in self._solo]:
                frames.extend(self._paced_frames("subscribe", [market_id], self._subscribes_in_flight, now))
            frames.extend(self._paced_frames("subscribe", [market_id for market_id in due if market_id not in self._solo],
                                             self._subscribes_in_flight, now))
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
        if self._orders_in_flight is not None and now - self._orders_in_flight[1] > ANSWER_WAIT_SEC:
            self._orders_in_flight = None
        if self._orders_refused_at is not None and now - self._orders_refused_at > REFUSED_RETRY_SEC:
            self._orders_refused_at = None

    def _orders_channel_frames(self, now: float) -> list[dict]:
        """The private orders channel's one verb due now, if any: its subscribe
        (first thing on a connection, ahead of every market), or the snapshot
        healing a gap. One at a time, each awaited ANSWER_WAIT_SEC."""
        if not self._orders_channel or self._orders_in_flight is not None or self._orders_refused_at is not None:
            return []
        if not self._orders_subscribed:
            verb = "subscribe"
        elif not self._orders_synced:
            verb = "snapshot"
        else:
            return []
        if not self._bucket.take(PRIVATE_CHANNEL_COST):
            return []
        nonce = self._next_nonce()
        self._orders_in_flight = (nonce, now)
        return [{"nonce": nonce, verb: {"private": [ORDERS_CHANNEL]}}]

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
                self._background.discard(market_id)
                self._solo.discard(market_id)
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
            if verb == "subscribe":
                self._subscribe_nonces.add(nonce)
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
        if self._orders_channel:
            self._handle_orders_part(message, nonce)
        snapshot = {str(market_id): body for market_id, body in (message.get("snapshot") or {}).items()}
        subscribed = self._subscribes_in_flight.pop(nonce, None) if isinstance(nonce, int) else None
        healed = self._snapshots_in_flight.pop(nonce, None) if isinstance(nonce, int) else None
        now = self._clock()
        if subscribed is not None:
            for market_id in subscribed[0]:
                if market_id not in snapshot:
                    self._refused[market_id] = (now, "the subscribe answer has no snapshot of it (no such market?)")
        answers_a_subscribe = isinstance(nonce, int) and nonce in self._subscribe_nonces
        for market_id, body in snapshot.items():
            if market_id not in self._books:
                # Only a subscribe's answer subscribes, a late one included: kept if still wanted, else unsubscribed.
                if not answers_a_subscribe:
                    continue
                if market_id not in self._wanted:
                    self._stray.add(market_id)
                    continue
                self._subscribed.add(market_id)
                self._books[market_id] = _MarketBook()
            self._guarded(market_id, self._apply_snapshot, body)
        # An answer that leaves a book out of sync (no book in it, or one unreadable) is asked again a minute
        # later; a book its replay found short is a new gap, healed at once.
        for market_id in (subscribed or ([], 0))[0] + (healed or ([], 0))[0]:
            book = self._books.get(market_id)
            if book is not None and not book.synced and not book.replay_gap:
                self._refused[market_id] = (now, "its snapshot answer left the book out of sync")
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
            log.warning("novig stream: market %s: unreadable frame part (%s: %s); resyncing it; part: %s",
                        market_id, type(error).__name__, error, json.dumps(body, default=str)[:LOGGED_PART_CHARS])
            book.synced = False
            book.buffered = []

    # ---- the private orders channel (called with the lock held) ----

    def _handle_orders_part(self, message: dict, nonce: object) -> None:
        """The private parts of one frame: the subscribe or snapshot answer's
        `orders` {seq, open}, a delta's `orders` {seq, deltas}, or a heartbeat."""
        answers_ours = self._orders_in_flight is not None and nonce == self._orders_in_flight[0]
        if answers_ours:
            self._orders_in_flight = None
        subscribed_private = (message.get("subscribed") or {}).get("private") if isinstance(message.get("subscribed"), dict) else None
        if answers_ours and isinstance(subscribed_private, list) and ORDERS_CHANNEL in subscribed_private:
            self._orders_subscribed = True
        orders = message.get("orders")
        heartbeat = message.get("heartbeat")
        try:
            if isinstance(orders, dict) and "open" in orders:
                if answers_ours and not self._orders_subscribed:
                    # A snapshot answer to our subscribe without the ack list still means subscribed.
                    self._orders_subscribed = True
                self._load_orders_snapshot(orders)
            elif isinstance(orders, dict) and "deltas" in orders:
                self._step_orders(int(orders["seq"]), orders.get("deltas") or [])
            elif answers_ours and not self._orders_subscribed:
                self._orders_refused_at = self._clock()
                self._orders_error = "the private subscribe answer named no orders channel"
                log.error("novig stream: %s: %s", self._orders_error, json.dumps(message, default=str)[:LOGGED_PART_CHARS])
            if isinstance(heartbeat, dict) and isinstance(heartbeat.get(ORDERS_CHANNEL), int):
                self._check_orders_heartbeat(heartbeat[ORDERS_CHANNEL])
        except (KeyError, ValueError, TypeError) as error:
            log.warning("novig stream: orders channel: unreadable frame part (%s: %s); resyncing it; part: %s",
                        type(error).__name__, error, json.dumps(orders if orders is not None else message, default=str)[:LOGGED_PART_CHARS])
            self._orders_synced = False
            self._orders_buffered = []

    def _load_orders_snapshot(self, orders: dict) -> None:
        seq = int(orders["seq"])
        now = self._clock()
        for open_order in orders.get("open") or []:
            if not isinstance(open_order, dict):
                raise TypeError(f"expected an open order object, got {str(open_order)[:200]}")
            state = self._order_state(order_id_of(open_order), open_order, now)
            state["final"] = STATUS_OPEN_ORDER
            state["remaining"] = open_order.get("qty")
            self._bump_order(state, now)
        self._orders_seq = seq
        self._orders_synced = True
        self._orders_epoch += 1
        self._orders_error = None
        buffered, self._orders_buffered = sorted(self._orders_buffered, key=lambda item: item[0]), []
        for buffered_seq, deltas in buffered:
            if buffered_seq > self._orders_seq:
                self._step_orders(buffered_seq, deltas)
        log.info("novig stream: orders channel in sync at seq %d (%d open order(s), epoch %d) on key %s",
                 seq, len(orders.get("open") or []), self._orders_epoch, self._client.key_id[:8])

    def _step_orders(self, seq: int, deltas: list) -> None:
        """One private batch in sequence; a stale seq is ignored; a skipped one
        puts the channel out of sync (a snapshot is asked for) and buffers it."""
        if not self._orders_synced:
            self._orders_buffered = [*self._orders_buffered, (seq, deltas)][-MAX_BUFFERED_DELTAS:]
            return
        if seq <= self._orders_seq:
            return
        if seq != self._orders_seq + 1:
            log.warning("novig stream: orders channel skipped seq %d -> %d; asking for a snapshot", self._orders_seq, seq)
            self._orders_synced = False
            self._orders_buffered = [(seq, deltas)]
            return
        now = self._clock()
        for delta in deltas:
            self._apply_order_delta(delta, now)
        self._orders_seq = seq
        self._prune_order_states(now)

    def _check_orders_heartbeat(self, heartbeat_seq: int) -> None:
        if self._orders_synced and self._orders_seq is not None and heartbeat_seq > self._orders_seq:
            log.warning("novig stream: orders heartbeat at seq %d, ours %d: a frame was dropped; asking for a snapshot",
                        heartbeat_seq, self._orders_seq)
            self._orders_synced = False
            self._orders_buffered = []

    def _order_state(self, order_id: str, entry: dict, now: float) -> dict:
        state = self._order_states.get(order_id)
        if state is None:
            state = {"orderId": order_id, "clientId": None, "outcomeId": None, "marketId": None, "final": None,
                     "filled": 0, "cost": 0.0, "remaining": None, "fills": [], "version": 0, "at": now}
            self._order_states[order_id] = state
        for key in ("clientId", "outcomeId", "marketId"):
            if isinstance(entry.get(key), str) and entry[key]:
                state[key] = entry[key]
        return state

    def _bump_order(self, state: dict, now: float) -> None:
        self._orders_version += 1
        state["version"] = self._orders_version
        state["at"] = now

    def _apply_order_delta(self, delta: dict, now: float) -> None:
        if not isinstance(delta, dict):
            raise TypeError(f"expected an orders delta object, got {str(delta)[:200]}")
        kind = delta.get("kind")
        state = self._order_state(order_id_of(delta), delta, now)
        if kind == "fill":
            qty = delta["qty"]
            price = float(delta["price"])
            if isinstance(qty, bool) or not isinstance(qty, int) or qty < 1 or not 0 < price < 1:
                raise ValueError(f"fill {delta.get('fillId')!r}: qty {qty!r}, price {delta.get('price')!r}")
            fill_id = str(delta.get("fillId") or f"{state['orderId']}:{len(state['fills']) + 1}")
            if any(fill["fillId"] == fill_id for fill in state["fills"]):
                return  # a snapshot replay or a repeated frame: the fill is already counted
            fill = {"fillId": fill_id, "orderId": state["orderId"], "clientId": state["clientId"],
                    "outcomeId": state["outcomeId"], "price": delta["price"], "qty": qty,
                    "originalQty": delta.get("originalQty"), "remaining": delta.get("remaining"),
                    "taker": delta.get("taker"), "ts": self._wall_clock_ms()}
            state["fills"].append(fill)
            state["filled"] += qty
            state["cost"] = round(state["cost"] + qty * price / QTY_PER_DOLLAR, 5)
            state["remaining"] = delta.get("remaining")
            if delta.get("remaining") == 0:
                state["final"] = STATUS_FILLED_ORDER
            self._bump_order(state, now)
            self._fill_events = [*self._fill_events, (self._orders_version, fill)][-MAX_FILL_EVENTS:]
        elif kind == "cancel":
            state["final"] = STATUS_CANCELED_ORDER
            state["remaining"] = delta.get("remaining")
            self._bump_order(state, now)
        elif kind == "reject":
            state["final"] = STATUS_REJECTED_ORDER
            state["remaining"] = delta.get("originalQty")
            self._bump_order(state, now)
        elif kind == "open":
            state["final"] = STATUS_OPEN_ORDER
            state["remaining"] = delta.get("qty")
            self._bump_order(state, now)
        else:
            raise ValueError(f"expected an orders delta of kind open, fill, cancel or reject, got {delta!r}")

    def _prune_order_states(self, now: float) -> None:
        for order_id in [order_id for order_id, state in self._order_states.items() if now - state["at"] > ORDER_STATE_KEEP_SEC]:
            del self._order_states[order_id]

    def _handle_error(self, nonce: object, code: str, text: str) -> None:
        reason = f"{code} {text}".strip()
        if self._orders_in_flight is not None and nonce == self._orders_in_flight[0]:
            self._orders_in_flight = None
            self._orders_refused_at = self._clock()
            self._orders_error = reason
            log.error("novig stream: the private orders %s was refused: %s",
                      "snapshot" if self._orders_subscribed else "subscribe", reason)
            return
        for in_flight, what in ((self._subscribes_in_flight, "subscribe"), (self._snapshots_in_flight, "gap snapshot")):
            sent = in_flight.pop(nonce, None) if isinstance(nonce, int) else None
            if sent is None:
                continue
            now = self._clock()
            if what == "subscribe" and len(sent[0]) > 1:
                # Which market Novig refused is unknown: each goes again in a frame of its own.
                self._solo.update(sent[0])
            else:
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
        # The lifecycle first: a go-live arrives with its removes, and an unreadable remove must not lose it.
        lifecycle = body.get("lifecycle")
        if isinstance(lifecycle, dict):
            self._apply_lifecycle(book, lifecycle)
        book_body = body.get("book")
        if isinstance(book_body, dict) and book.seq is not None:
            book.step(int(book_body["seq"]), book_body.get("deltas") or [])

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
