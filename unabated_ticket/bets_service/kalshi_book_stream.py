"""Kalshi's live order books over ONE websocket (Kalshi's `orderbook_delta`
channel), signed with the bets service's READ-ONLY Kalshi key
(KALSHI_API_KEY_ID) and used for reads only. The runner asks POST
/kalshi_books.json (service.py) for the books its live price pass, a line's
reread and a Place quote need, and the stream subscribes to those markets on
demand. The automated taker reads REST on its own pass, never this.

Kalshi facts this rests on (docs.kalshi.com/websockets/websocket-connection,
/orderbook-updates and /getting_started/quick_start_websockets, read
2026-11-10; NOT yet seen live, since the cloud container cannot reach Kalshi):
  * wss://api.elections.kalshi.com/trade-api/ws/v2 (the host the MLB maker
    streams from live; the docs also list external-api-ws.kalshi.com). The
    upgrade carries KALSHI-ACCESS-KEY / -TIMESTAMP (ms) / -SIGNATURE, the REST
    scheme signed over timestamp + "GET" + "/trade-api/ws/v2".
  * Commands are {"id": n, "cmd": ..., "params": {...}}, ids unique per
    connection. {"cmd": "subscribe", "params": {"channels":
    ["orderbook_delta"], "market_tickers": [..]}} answers {"id", "type":
    "subscribed", "msg": {"channel", "sid"}}. {"cmd": "update_subscription",
    "params": {"sid", "market_tickers", "action": "add_markets" |
    "delete_markets"}} answers {"id", "sid", "seq", "type": "ok"}. {"cmd":
    "unsubscribe", "params": {"sids": [..]}}. An error is {"id", "type":
    "error", "msg": {"code", "msg"}}; code 26 is "per-subscription market
    limit exceeded" (the limit itself is unpublished), 27 a command rate limit.
  * Each market first gets {"type": "orderbook_snapshot", "sid", "seq", "msg":
    {market_ticker, yes_dollars_fp?, no_dollars_fp?}}: every resting BID on
    each side as [price dollars, contracts] strings, a side absent when empty
    (the REST /orderbook's orderbook_fp.yes_dollars / no_dollars). Then each
    change is {"type": "orderbook_delta", "sid", "seq", "msg": {market_ticker,
    price_dollars, delta_fp, side}}: that level's contracts change by delta_fp.
  * `seq` counts the messages of one subscription; a gap means a message was
    missed. The docs show an `ok` answer carrying the sid's seq too, so an
    `ok` that is the next number is counted and one that is not is ignored.

A book is served while the connection is up and its subscription has missed
no message since the market's snapshot; `ageMs` says how long since it last
changed. A still market gets no frames and its book stays exact, but a stream
that quietly stopped sending would look the same, so every subscription is
renewed for fresh snapshots each RENEW_SEC (the runner also checks a few stream
books against REST each pass). A gap, a delta before its snapshot or a level
below zero resubscribes the whole subscription, no sooner than a backoff that
doubles while it keeps happening, and serves none of its books until the new
snapshots land. A crossed book (best yes bid + best no bid at or over $1) is
never served: one exchange action sends several deltas. A snapshot with
neither documented side key but a yes/no key of another name stops the
connection loudly rather than serving an empty book. Anything else answers an
`error` naming why, never an old book; the runner then reads REST.

Commands go out at most one every COMMAND_GAP_SEC across the connection (a
rate-limit refusal pauses them, and the markets are asked again, not refused).

Side effects: one outbound websocket to Kalshi on the read-only key (reads
only: nothing here can place, cancel or move money). Nothing on disk; books
are held in memory only.
"""
import json
import logging
import threading
import time
from collections.abc import Callable
from typing import Protocol

from unabated_ticket.bets_service.novig_book_stream import websocket_client_connect

log = logging.getLogger(__name__)

STREAM_PATH = "/trade-api/ws/v2"
ORDERBOOK_CHANNEL = "orderbook_delta"
PING_SEC = 15
# No frame of any kind (a pong or Kalshi's own ping included) this long: the line is dead.
SILENT_RECONNECT_SEC = 45
RECV_POLL_SEC = 0.25
RECONNECT_MIN_SEC = 1
RECONNECT_MAX_SEC = 30
# Only a connection that held this long resets the backoff.
STABLE_CONNECTION_SEC = 60
REPEAT_FAILURE_LOG_SEC = 10 * 60
# A subscribe with no `subscribed` answer this long is sent again; a market
# sent in a subscribe or add with no snapshot this long reads as an error.
ANSWER_WAIT_SEC = 10
# This many subscribes in a row with no answer: the connection is useless, so it is dropped (and the error shown).
MAX_UNANSWERED_SUBSCRIBES = 3
# One command (subscribe, add, delete, resubscribe) at most this often per connection.
COMMAND_GAP_SEC = 0.2
# After Kalshi's command rate-limit refusal (error 27), no command for this long.
RATE_LIMIT_PAUSE_SEC = 5
# Every subscription is resubscribed this often, for fresh snapshots of every market.
RENEW_SEC = 10 * 60
# A missed message resubscribes no sooner than this after the last one, doubling
# while they keep coming (up to the max) and starting over after a quiet spell.
RESYNC_MIN_GAP_SEC = 1
RESYNC_MAX_GAP_SEC = 60
RESYNC_RESET_SEC = 10 * 60
# A market Kalshi refused is asked for again no sooner than this.
REFUSED_RETRY_SEC = 60
# A market nobody asked about for this long is dropped (games end).
WATCH_IDLE_SEC = 10 * 60
# Kalshi answers only the first orderbook_delta subscribe on a connection and
# ignores the rest (VM log, 2026-11-11), so every market rides ONE subscription,
# grown with update_subscription add_markets, as Kalshi documents.
MAX_SUBSCRIPTIONS = 1
# Kalshi publishes no per-subscription market limit (error 26 names one), so the
# subscription holds at most this many; a subscribe Kalshi answers 26 halves it
# for the connection, an add answered 26 closes the subscription at its size.
MARKETS_PER_SUBSCRIPTION = 1000
# A subscribe or add carries at most this many markets; the rest follow in adds.
MARKETS_PER_COMMAND = 50
# A subscribe this small that Kalshi still answers 26 refuses its markets instead of halving again.
MIN_MARKETS_PER_SUBSCRIPTION = 5
# The most a caller may wait for a new market's first book (the runner's own timeouts sit above it).
MAX_WAIT_SEC = 3.0
# The runner's live price pass may hold at most this many markets (of the
# subscriptions' 1,000); past it a background market answers an error, so a
# reread's or Place's market always has room.
MAX_BACKGROUND_MARKETS = 900
# Kalshi's error codes: "per-subscription market limit exceeded" and "exceeded its command rate limit".
MARKET_LIMIT_ERROR_CODE = 26
COMMAND_RATE_ERROR_CODE = 27
# Prices are held in ten-thousandths of a dollar (Kalshi's sub-cent ticks are 0.001).
PRICE_UNITS_PER_DOLLAR = 10000
# A level whose contracts fall to this or under is gone (Kalshi sends 2-decimal counts).
EMPTY_LEVEL_QTY = 1e-9
# A level driven below this many contracts means a delta was misapplied or missed.
NEGATIVE_LEVEL_TOLERANCE = -0.005
SNAPSHOT_SIDE_KEYS = ("yes_dollars_fp", "no_dollars_fp")


class Signer(Protocol):
    """taker_orders.KalshiOrderClient on the read key: the handshake's signed headers."""

    key_id: str

    def signed_ws_headers(self, ws_path: str) -> dict: ...


class Connection(Protocol):
    """One open websocket (novig_book_stream's websocket-client wrapper).
    recv returns a text message, "" for a control frame (a ping or pong), or
    None when nothing arrived within the timeout; raises ConnectionError once it is gone."""

    def send(self, text: str) -> None: ...

    def recv(self, timeout_sec: float) -> str | None: ...

    def ping(self) -> None: ...

    def close(self) -> None: ...


Connect = Callable[[str, dict], Connection]


def kalshi_connect(url: str, headers: dict) -> Connection:
    return websocket_client_connect(url, headers, peer="Kalshi")


def price_units_of(price_dollars: object) -> int:
    """0.4500 (a string or number of dollars) -> 4500. Raises ValueError outside (0, 1)."""
    units = round(float(price_dollars) * PRICE_UNITS_PER_DOLLAR)
    if not 0 < units < PRICE_UNITS_PER_DOLLAR:
        raise ValueError(f"expected a Kalshi price strictly between 0 and 1 dollar, got {price_dollars!r}")
    return units


def levels_of(raw_levels: object, what: str) -> dict[int, float]:
    """A snapshot side's [[price, contracts], ..] -> {price units: contracts}; None (an empty side) -> {}."""
    if raw_levels is None:
        return {}
    if not isinstance(raw_levels, list):
        raise TypeError(f"expected a list of [price, contracts] for {what}, got {str(raw_levels)[:200]}")
    levels: dict[int, float] = {}
    for level in raw_levels:
        if not isinstance(level, list) or len(level) != 2:
            raise TypeError(f"expected a [price, contracts] pair in {what}, got {str(level)[:200]}")
        quantity = float(level[1])
        if quantity > EMPTY_LEVEL_QTY:
            levels[price_units_of(level[0])] = quantity
    return levels


def check_snapshot_keys(msg: dict) -> None:
    """Raise when a snapshot carries neither documented side key but some other
    yes/no key: Kalshi renamed them, and reading it would serve an empty book."""
    if any(key in msg for key in SNAPSHOT_SIDE_KEYS):
        return
    other = [key for key in msg if "yes" in key or "no" in key]
    if other:
        raise TypeError(f"expected Kalshi snapshot sides {SNAPSHOT_SIDE_KEYS}, got keys {sorted(msg)}")


def rest_shape_of(book: dict) -> dict:
    """A held book as the REST /orderbook's body ({orderbook_fp: {yes_dollars,
    no_dollars}}, best bid first), so the runner's taker.parseKalshiOrderbook reads it unchanged."""
    def side(levels: dict[int, float]) -> list[list[str]]:
        return [[f"{units / PRICE_UNITS_PER_DOLLAR:.4f}", f"{quantity:.2f}"]
                for units, quantity in sorted(levels.items(), reverse=True)]
    return {"orderbook_fp": {"yes_dollars": side(book["yes"]), "no_dollars": side(book["no"])}}


def crossed(book: dict) -> bool:
    """Best yes bid + best no bid at or over $1: two bids that would have matched."""
    if not book["yes"] or not book["no"]:
        return False
    return max(book["yes"]) + max(book["no"]) >= PRICE_UNITS_PER_DOLLAR


class KalshiBookStream:
    """The books of the Kalshi markets the runner prices, kept live on one
    websocket. run(stop) is the thread's loop; books_for() is called from HTTP threads."""

    def __init__(self, signer: Signer, url: str, connect: Connect = kalshi_connect,
                 clock: Callable[[], float] = time.monotonic):
        self._signer = signer
        self._url = url
        self._connect = connect
        self._clock = clock
        self._cond = threading.Condition()
        # Kept across connections: ticker -> when last asked for, a reread's/Place's, the live pass's alone, refusals.
        self._wanted: dict[str, float] = {}
        self._priority: set[str] = set()
        self._background: set[str] = set()
        self._refused: dict[str, tuple[float, str]] = {}
        # Kept across connections, since a cause that misses messages outlives one:
        # group index -> (resyncs in a row, when the last one was).
        self._resync_streaks: dict[int, tuple[int, float]] = {}
        self._resyncs = 0
        self._last_resync: str | None = None
        self._reset_connection_state()
        self._connection: Connection | None = None
        self._connected_at: float | None = None
        self._last_frame_at: float | None = None
        self._last_ping_at: float | None = None
        self._last_error: str | None = None
        # Kept across connections, since _open clears the error: how many times a connection dropped, and why last.
        self._disconnects = 0
        self._last_disconnect: tuple[str, float] | None = None
        self._last_logged_failure: tuple[str, float] | None = None

    def _reset_connection_state(self) -> None:
        # Command ids are per connection.
        self._command_seq = 0
        self._last_command_at = -float("inf")
        self._commands_paused_until = -float("inf")
        self._unanswered_subscribes = 0
        # A smaller markets-per-subscription cap Kalshi taught with error 26 on a subscribe (relearned per connection).
        self._learned_market_limit: int | None = None
        # Frames the frame handler wants sent (it runs under the lock; the pump sends).
        self._queued: list[dict] = []
        # group index -> its subscription on this connection: {subscribeId, sid, tickers (on Kalshi's
        # side, pending adds included), sentAt, answeredAt, lastSeq, resync, resyncAfter, full}.
        self._groups: dict[int, dict] = {}
        # command id -> (group index, "subscribe" | "add_markets" | "delete_markets", tickers,
        # the subscribeId of the subscription it was sent on).
        self._commands: dict[int, tuple[int, str, tuple[str, ...], int]] = {}
        self._group_of_sid: dict[int, int] = {}
        # ticker -> group index it is (or is to be) subscribed in.
        self._group_of_ticker: dict[str, int] = {}
        # ticker -> when it was last sent in a subscribe or add (its snapshot's wait runs from it).
        self._sent_at: dict[str, float] = {}
        # ticker -> {"yes": {price units: contracts}, "no": {...}, "changedAt"}, this connection's books only.
        self._books: dict[str, dict] = {}
        self._housekeeping_due = True
        self._last_housekeeping_at = -float("inf")

    def _markets_per_subscription(self) -> int:
        if self._learned_market_limit is None:
            return MARKETS_PER_SUBSCRIPTION
        return min(MARKETS_PER_SUBSCRIPTION, self._learned_market_limit)

    # ---- what the HTTP route calls ----

    def books_for(self, tickers: list[str], wait_sec: float, priority: bool = False, background: bool = False) -> dict:
        """{stream: status, books: {ticker: view}}. Marks each market wanted
        (subscribing it if new; `priority` ahead of the rest, for a reread or
        Place; `background` behind the rest and capped, for the live price
        pass, until a priority ask names it) and waits up to wait_sec (capped
        at MAX_WAIT_SEC) for every new market's snapshot. A view is
        {orderbook_fp, ageMs} (the REST /orderbook's body plus how long since
        the book last changed) or {error}."""
        wait_sec = max(0.0, min(wait_sec, MAX_WAIT_SEC))
        with self._cond:
            now = self._clock()
            capped: set[str] = set()
            for ticker in tickers:
                if background and ticker not in self._wanted:
                    if len(self._background) >= MAX_BACKGROUND_MARKETS:
                        capped.add(ticker)
                        continue
                    self._background.add(ticker)
                elif not background:
                    self._background.discard(ticker)
                if ticker not in self._wanted:
                    self._housekeeping_due = True
                self._wanted[ticker] = now
                if priority:
                    self._priority.add(ticker)
            self._cond.notify_all()
            deadline = now + wait_sec
            while True:
                remaining = deadline - self._clock()
                views = {ticker: self._view_of(ticker, remaining) for ticker in tickers if ticker not in capped}
                if remaining <= 0 or not any(view.get("pending") for view in views.values()):
                    break
                self._cond.wait(remaining)
            for view in views.values():
                view.pop("pending", None)
            for ticker in capped:
                views[ticker] = {"error": f"the live price pass already watches {MAX_BACKGROUND_MARKETS} Kalshi markets"}
            return {"stream": self._status_locked(), "books": views}

    def status(self) -> dict:
        with self._cond:
            return self._status_locked()

    def _status_locked(self) -> dict:
        now = self._clock()
        return {
            "url": self._url,
            "connected": self._connection is not None,
            "connectedForSec": None if self._connected_at is None else round(now - self._connected_at),
            "lastFrameAgeMs": None if self._connection is None or self._last_frame_at is None
            else int((now - self._last_frame_at) * 1000),
            "error": self._last_error,
            "disconnects": self._disconnects,
            "lastDisconnect": None if self._last_disconnect is None else self._last_disconnect[0],
            "lastDisconnectAgoSec": None if self._last_disconnect is None
            else round(now - self._last_disconnect[1]),
            "subscriptions": len(self._groups),
            "answeredSubscriptions": sum(1 for group in self._groups.values() if group["sid"] is not None),
            "subscribedMarkets": sum(len(group["tickers"]) for group in self._groups.values()),
            "booksHeld": len(self._books),
            "wantedMarkets": len(self._wanted),
            "refusedMarkets": len(self._refused),
            "marketsPerSubscription": self._markets_per_subscription(),
            "resyncs": self._resyncs,
            "lastResync": self._last_resync,
        }

    def _view_of(self, ticker: str, wait_left_sec: float = 0.0) -> dict:
        if self._connection is None:
            return {"error": f"the Kalshi stream is not connected ({self._last_error or 'connecting'})"}
        refused = self._refused.get(ticker)
        if refused is not None:
            return {"error": f"Kalshi refused the stream subscription to {ticker}: {refused[1]}"}
        group = self._groups.get(self._group_of_ticker.get(ticker, -1))
        book = self._books.get(ticker)
        if book is not None and group is not None and not group["resync"]:
            if crossed(book):
                return {"error": f"the stream's book of {ticker} is crossed (mid-update)"}
            view = rest_shape_of(book)
            view["ageMs"] = int((self._clock() - book["changedAt"]) * 1000)
            return view
        if group is not None and group["resync"]:
            why = f"the stream is resubscribing {ticker} after a missed message"
            # A resubscribe further off than the caller waits is no reason to wait.
            if group["resyncAfter"] - self._clock() > wait_left_sec:
                return {"error": why}
            return {"pending": True, "error": why}
        sent_at = self._sent_at.get(ticker)
        if group is not None and group["sid"] is not None and ticker in group["tickers"] and sent_at is not None \
                and self._clock() - sent_at > ANSWER_WAIT_SEC:
            return {"error": f"Kalshi's stream sent no snapshot of {ticker} within {ANSWER_WAIT_SEC} s of subscribing it"}
        if ticker in self._wanted and ticker not in self._group_of_ticker and self._no_room():
            return {"error": "the stream's subscription takes no more markets "
                             f"(at most {self._markets_per_subscription()}, or Kalshi's limit)"}
        return {"pending": True, "error": f"waiting for Kalshi's stream snapshot of {ticker}"}

    def _room_of(self, index: int, sizes: dict[int, int]) -> bool:
        group = self._groups.get(index)
        return sizes[index] < self._markets_per_subscription() and not (group is not None and group["full"])

    def _no_room(self) -> bool:
        sizes = self._assigned_sizes()
        return not any(self._room_of(index, sizes) for index in range(MAX_SUBSCRIPTIONS))

    def _assigned_sizes(self) -> dict[int, int]:
        sizes = {index: 0 for index in range(MAX_SUBSCRIPTIONS)}
        for index in self._group_of_ticker.values():
            sizes[index] += 1
        return sizes

    # ---- the stream thread ----

    def run(self, stop: threading.Event) -> None:
        """Connect, keep the wanted markets subscribed, apply every snapshot
        and delta; reconnect with backoff whenever the line drops or stops answering."""
        backoff = RECONNECT_MIN_SEC
        while not stop.is_set():
            opened_at = None
            try:
                self._open()
                opened_at = self._clock()
                self._pump(stop)
            except (OSError, ConnectionError, ValueError, KeyError, TypeError) as error:
                self._drop(f"{type(error).__name__}: {error}")
            # Anything else too: a dead thread would leave every Kalshi read on REST until a restart.
            except Exception as error:
                log.exception("kalshi stream: unexpected error")
                self._drop(f"unexpected {type(error).__name__}: {error}")
            if opened_at is not None and self._clock() - opened_at > STABLE_CONNECTION_SEC:
                backoff = RECONNECT_MIN_SEC
            if stop.wait(backoff):
                break
            backoff = min(backoff * 2, RECONNECT_MAX_SEC)
        self._drop("stopped")

    def _open(self) -> None:
        headers = self._signer.signed_ws_headers(STREAM_PATH)
        connection = self._connect(self._url, headers)
        with self._cond:
            self._reset_connection_state()
            self._connection = connection
            self._connected_at = self._last_frame_at = self._last_ping_at = self._clock()
            self._last_error = None
            self._cond.notify_all()
        log.info("kalshi stream: connected to %s on read key %s", self._url, self._signer.key_id[:8])

    def _drop(self, reason: str) -> None:
        with self._cond:
            connection = self._connection
            self._connection = None
            self._connected_at = None
            self._reset_connection_state()
            if reason != "stopped":
                self._last_error = reason
                if connection is not None:
                    self._disconnects += 1
                    self._last_disconnect = (reason, self._clock())
            self._cond.notify_all()
        if connection is not None:
            connection.close()
        if reason == "stopped":
            return
        # The same failure on every reconnect is logged once per REPEAT_FAILURE_LOG_SEC.
        now = self._clock()
        last = self._last_logged_failure
        if last is None or last[0] != reason or now - last[1] > REPEAT_FAILURE_LOG_SEC:
            log.warning("kalshi stream: %s: %s", "disconnected" if connection is not None else "could not connect", reason)
            self._last_logged_failure = (reason, now)

    def _pump(self, stop: threading.Event) -> None:
        while not stop.is_set():
            for frame in self._outgoing_frames():
                self._connection.send(json.dumps(frame, separators=(",", ":")))
            now = self._clock()
            if now - self._last_ping_at >= PING_SEC:
                self._connection.ping()
                self._last_ping_at = now
            text = self._connection.recv(RECV_POLL_SEC)
            now = self._clock()
            if text is None:
                if now - self._last_frame_at > SILENT_RECONNECT_SEC:
                    raise ConnectionError(f"no frame, ping or pong from Kalshi for {SILENT_RECONNECT_SEC} s")
                continue
            with self._cond:
                self._last_frame_at = now
                if text:
                    self._handle(json.loads(text))
                self._cond.notify_all()

    def _next_id(self) -> int:
        self._command_seq += 1
        return self._command_seq

    def _command(self, index: int, kind: str, tickers: tuple[str, ...], cmd: str, params: dict,
                 subscribe_id: int | None = None) -> dict:
        command_id = self._next_id()
        self._commands[command_id] = (index, kind, tickers, command_id if subscribe_id is None else subscribe_id)
        return {"id": command_id, "cmd": cmd, "params": params}

    def _unsubscribe_frame(self, sid: int) -> dict:
        return {"id": self._next_id(), "cmd": "unsubscribe", "params": {"sids": [sid]}}

    def _outgoing_frames(self) -> list[dict]:
        """The frames to send now: any the frame handler queued, then, with the
        housekeeping due (a new ask, or RECV_POLL_SEC since the last), drop
        idle and refused markets and place new ones (a reread's or Place's
        first, the live pass's last), then at most ONE command if
        COMMAND_GAP_SEC has passed and no rate-limit pause is on: the first
        subscription (by index) that needs a subscribe, a resubscribe (missed
        message, unanswered subscribe, renewal) or an add or delete."""
        with self._cond:
            frames, self._queued = self._queued, []
            now = self._clock()
            if not self._housekeeping_due and now - self._last_housekeeping_at < RECV_POLL_SEC:
                return frames
            self._housekeeping_due = False
            self._last_housekeeping_at = now
            for ticker, (refused_at, _reason) in list(self._refused.items()):
                if now - refused_at > REFUSED_RETRY_SEC:
                    del self._refused[ticker]
            self._forget_idle(now)
            for ticker in [ticker for ticker in self._group_of_ticker if ticker in self._refused]:
                del self._group_of_ticker[ticker]
            self._assign_new_markets()
            if now - self._last_command_at < COMMAND_GAP_SEC or now < self._commands_paused_until:
                return frames
            for index in range(MAX_SUBSCRIPTIONS):
                command_frames = self._group_frames(index, now)
                if command_frames:
                    self._last_command_at = now
                    return frames + command_frames
            return frames

    def _forget_idle(self, now: float) -> None:
        for ticker in [ticker for ticker, asked_at in self._wanted.items() if now - asked_at > WATCH_IDLE_SEC]:
            del self._wanted[ticker]
            self._priority.discard(ticker)
            self._background.discard(ticker)
            self._group_of_ticker.pop(ticker, None)

    def _placement_order(self, ticker: str) -> tuple[bool, bool, str]:
        """A reread's or Place's market first, the live pass's last."""
        return (ticker not in self._priority, ticker in self._background, ticker)

    def _assign_new_markets(self) -> None:
        new = [ticker for ticker in self._wanted if ticker not in self._group_of_ticker and ticker not in self._refused]
        new.sort(key=self._placement_order)
        # A market still on a subscription (its delete not sent yet) goes back to that one, never a second.
        still_on = {ticker: index for index, group in self._groups.items() for ticker in group["tickers"]}
        sizes = self._assigned_sizes()
        for ticker in new:
            if ticker in still_on:
                self._group_of_ticker[ticker] = still_on[ticker]
                sizes[still_on[ticker]] += 1
                continue
            room = [index for index in range(MAX_SUBSCRIPTIONS) if self._room_of(index, sizes)]
            if not room:
                continue
            index = min(room, key=lambda candidate: (self._groups.get(candidate, {}).get("sid") is None, candidate))
            self._group_of_ticker[ticker] = index
            sizes[index] += 1

    def _subscribe_frame(self, index: int, tickers: list[str], now: float, full: bool = False) -> dict:
        """The subscribe of the first MARKETS_PER_COMMAND markets; the rest are added once Kalshi answers."""
        tickers = tickers[:MARKETS_PER_COMMAND]
        frame = self._command(index, "subscribe", tuple(tickers), "subscribe",
                              {"channels": [ORDERBOOK_CHANNEL], "market_tickers": tickers})
        self._groups[index] = {"subscribeId": frame["id"], "sid": None, "tickers": set(tickers), "sentAt": now,
                               "answeredAt": None, "lastSeq": None, "resync": False, "resyncAfter": now,
                               "full": full}
        for ticker in tickers:
            self._sent_at[ticker] = now
        return frame

    def _end_group(self, index: int) -> list[dict]:
        """Forget a subscription and its books; the unsubscribe to send when Kalshi had answered it."""
        group = self._groups.pop(index)
        for ticker in group["tickers"]:
            self._books.pop(ticker, None)
        if group["sid"] is None:
            return []
        self._group_of_sid.pop(group["sid"], None)
        return [self._unsubscribe_frame(group["sid"])]

    def _group_frames(self, index: int, now: float) -> list[dict]:
        assigned = sorted((ticker for ticker, of in self._group_of_ticker.items() if of == index),
                          key=self._placement_order)
        group = self._groups.get(index)
        if group is None:
            return [self._subscribe_frame(index, assigned, now)] if assigned else []
        if not assigned:
            return self._end_group(index)
        if group["resync"]:
            if now < group["resyncAfter"]:
                return []
            return self._resubscribe(index, assigned, now)
        if group["sid"] is None:
            if now - group["sentAt"] <= ANSWER_WAIT_SEC:
                return []
            self._unanswered_subscribes += 1
            if self._unanswered_subscribes >= MAX_UNANSWERED_SUBSCRIBES:
                raise ConnectionError(f"Kalshi answered none of {self._unanswered_subscribes} subscribes in a row "
                                      f"within {ANSWER_WAIT_SEC} s")
            log.warning("kalshi stream: subscribe %s of %d market(s) unanswered after %d s, sending it again",
                        group["subscribeId"], len(group["tickers"]), ANSWER_WAIT_SEC)
            return self._resubscribe(index, assigned, now)
        if now - group["sentAt"] >= RENEW_SEC:
            return self._resubscribe(index, assigned, now)
        to_add = [ticker for ticker in assigned if ticker not in group["tickers"]][:MARKETS_PER_COMMAND]
        if to_add:
            group["tickers"].update(to_add)
            for ticker in to_add:
                self._sent_at[ticker] = now
            return [self._command(index, "add_markets", tuple(to_add), "update_subscription",
                                  {"sid": group["sid"], "market_tickers": to_add, "action": "add_markets"},
                                  group["subscribeId"])]
        to_delete = sorted(group["tickers"].difference(assigned))
        if to_delete:
            group["tickers"].difference_update(to_delete)
            for ticker in to_delete:
                self._books.pop(ticker, None)
            return [self._command(index, "delete_markets", tuple(to_delete), "update_subscription",
                                  {"sid": group["sid"], "market_tickers": to_delete, "action": "delete_markets"},
                                  group["subscribeId"])]
        return []

    def _resubscribe(self, index: int, assigned: list[str], now: float) -> list[dict]:
        full = self._groups[index]["full"]
        frames = self._end_group(index)
        frames.append(self._subscribe_frame(index, assigned, now, full))
        return frames

    # ---- frames from Kalshi (called with the lock held) ----

    def _handle(self, message: object) -> None:
        if not isinstance(message, dict):
            raise TypeError(f"expected a JSON object frame from Kalshi, got {str(message)[:200]}")
        kind = message.get("type")
        if kind == "subscribed":
            self._handle_subscribed(message)
        elif kind == "ok":
            self._handle_ok(message)
        elif kind == "error":
            self._handle_error(message)
        elif kind in ("orderbook_snapshot", "orderbook_delta"):
            self._handle_book(kind, message)
        # Anything else (an unsubscribed answer, another channel) carries no book.

    def _handle_subscribed(self, message: dict) -> None:
        command = self._commands.pop(message.get("id"), None)
        msg = message.get("msg") or {}
        sid = msg.get("sid") if isinstance(msg, dict) else None
        if command is None or command[1] != "subscribe" or not isinstance(sid, int):
            log.warning("kalshi stream: a subscribed answer to no subscribe of ours: %s", json.dumps(message)[:300])
            return
        self._unanswered_subscribes = 0
        group = self._groups.get(command[0])
        if group is None or group["subscribeId"] != message.get("id"):
            # An older subscribe of a group already sent again: end it on Kalshi's side too.
            self._queued.append(self._unsubscribe_frame(sid))
            self._housekeeping_due = True
            return
        group["sid"] = sid
        group["answeredAt"] = self._clock()
        self._group_of_sid[sid] = command[0]
        self._housekeeping_due = True

    def _handle_ok(self, message: dict) -> None:
        self._commands.pop(message.get("id"), None)
        group = self._groups.get(self._group_of_sid.get(message.get("sid"), -1))
        seq = message.get("seq")
        # Counted only when it is the subscription's next number (the docs do not say ok shares the data's count).
        if group is not None and isinstance(seq, int) and group["lastSeq"] is not None and seq == group["lastSeq"] + 1:
            group["lastSeq"] = seq

    def _handle_error(self, message: dict) -> None:
        """A refused command. A rate-limit refusal pauses commands and asks the
        markets again; a market-limit refusal shrinks the subscriptions (a
        subscribe) or closes that one to new markets (an add), the markets
        placed again; any other refusal of a subscribe or add refuses its
        markets for REFUSED_RETRY_SEC."""
        msg = message.get("msg") if isinstance(message.get("msg"), dict) else {}
        code = msg.get("code")
        reason = f"code {code}: {str(msg.get('msg'))[:200]}"
        command = self._commands.pop(message.get("id"), None)
        if command is None:
            log.warning("kalshi stream: Kalshi answered an error: %s (id %s)", reason, message.get("id"))
            return
        index, kind, tickers, subscribe_id = command
        group = self._groups.get(index)
        self._housekeeping_due = True
        if kind == "delete_markets":
            log.error("kalshi stream: delete of %d market(s) refused: %s", len(tickers), reason)
            return
        if group is None or group["subscribeId"] != subscribe_id:
            # Sent on a subscription since sent again (renewed, resynced): the new one answers for its markets.
            return
        if code == COMMAND_RATE_ERROR_CODE:
            log.warning("kalshi stream: %s of %d market(s) rate-limited (%s), pausing commands %d s",
                        kind, len(tickers), reason, RATE_LIMIT_PAUSE_SEC)
            self._commands_paused_until = self._clock() + RATE_LIMIT_PAUSE_SEC
            self._unsend(index, kind, tickers)
            return
        if code == MARKET_LIMIT_ERROR_CODE and (kind != "subscribe" or len(tickers) > MIN_MARKETS_PER_SUBSCRIPTION):
            if kind == "subscribe":
                self._learned_market_limit = max(1, len(tickers) // 2)
                for ticker in tickers:
                    self._group_of_ticker.pop(ticker, None)
            else:
                group["full"] = True
                for ticker in tickers:
                    if self._group_of_ticker.get(ticker) == index:
                        del self._group_of_ticker[ticker]
            log.error("kalshi stream: %s of %d market(s) over Kalshi's per-subscription limit (%s); %s",
                      kind, len(tickers), reason, f"the subscription now holds {self._markets_per_subscription()}"
                      if kind == "subscribe" else "the subscription takes no more")
            self._unsend(index, kind, tickers)
            return
        log.error("kalshi stream: %s of %d market(s) refused: %s", kind, len(tickers), reason)
        now = self._clock()
        for ticker in tickers:
            self._refused[ticker] = (now, reason)
            if self._group_of_ticker.get(ticker) == index:
                del self._group_of_ticker[ticker]
        self._unsend(index, kind, tickers)

    def _unsend(self, index: int, kind: str, tickers: tuple[str, ...]) -> None:
        """Undo a refused subscribe (the subscription is forgotten) or add (its markets leave it)."""
        if kind == "subscribe":
            self._groups.pop(index, None)
        else:
            self._groups[index]["tickers"].difference_update(tickers)
        for ticker in tickers:
            self._books.pop(ticker, None)
            self._sent_at.pop(ticker, None)

    def _handle_book(self, kind: str, message: dict) -> None:
        index = self._group_of_sid.get(message.get("sid"))
        group = self._groups.get(index) if index is not None else None
        if group is None or group["resync"]:
            # A subscription already ended or being resubscribed: its frames are not ours to apply.
            return
        seq = message.get("seq")
        if not isinstance(seq, int):
            raise TypeError(f"expected an integer seq on a Kalshi {kind}, got {str(message)[:200]}")
        if group["lastSeq"] is not None and seq != group["lastSeq"] + 1:
            self._resync(index, f"seq {seq} after {group['lastSeq']}")
            return
        group["lastSeq"] = seq
        msg = message.get("msg")
        if not isinstance(msg, dict):
            raise TypeError(f"expected a msg object on a Kalshi {kind}, got {str(message)[:200]}")
        ticker = msg.get("market_ticker")
        if ticker not in group["tickers"]:
            # A market just deleted from the subscription.
            return
        if kind == "orderbook_snapshot":
            check_snapshot_keys(msg)
            self._books[ticker] = {"yes": levels_of(msg.get("yes_dollars_fp"), f"{ticker} yes"),
                                   "no": levels_of(msg.get("no_dollars_fp"), f"{ticker} no"),
                                   "changedAt": self._clock()}
            return
        book = self._books.get(ticker)
        if book is None:
            self._resync(index, f"a delta for {ticker} before its snapshot")
            return
        side = msg.get("side")
        if side not in ("yes", "no"):
            raise ValueError(f"expected side yes or no on a Kalshi delta for {ticker}, got {side!r}")
        units = price_units_of(msg.get("price_dollars"))
        quantity = book[side].get(units, 0.0) + float(msg.get("delta_fp"))
        if quantity < NEGATIVE_LEVEL_TOLERANCE:
            self._resync(index, f"{ticker} {side} at {msg.get('price_dollars')} went to {quantity:.2f} contracts")
            return
        if quantity > EMPTY_LEVEL_QTY:
            book[side][units] = quantity
        else:
            book[side].pop(units, None)
        book["changedAt"] = self._clock()

    def _resync(self, index: int, why: str) -> None:
        """Serve none of the subscription's books until it is resubscribed,
        after a gap that doubles while missed messages keep coming."""
        now = self._clock()
        streak, last_at = self._resync_streaks.get(index, (0, -float("inf")))
        streak = 0 if now - last_at > RESYNC_RESET_SEC else streak
        gap = min(RESYNC_MAX_GAP_SEC, RESYNC_MIN_GAP_SEC * 2 ** streak)
        self._resync_streaks[index] = (streak + 1, now)
        group = self._groups[index]
        group["resync"] = True
        group["resyncAfter"] = now + gap
        self._resyncs += 1
        self._last_resync = why
        self._housekeeping_due = True
        log.warning("kalshi stream: subscription %s missed a message (%s), resubscribing its %d market(s) in %d s",
                    group["sid"], why, len(group["tickers"]), gap)
