"""Polymarket US's live order books over ONE websocket
(wss://api.polymarket.us/v1/ws/markets), signed with the POLYMARKET_US_* key
and used for reads only. It replaces the runner's REST polling of the public
gateway's GET /v1/markets/<slug>/book: since the live price pass (#252) that
polling sent up to 60 reads every 30 s from the VM's address, the gateway
answered most of them HTTP 429 (it allows about 5 then a 10 s wait), and the
Place button's own book read from the same address was refused with it
(2026-11-05). The runner asks POST /polymarket_books.json (service.py) for the
books it needs, and the stream subscribes to those markets on demand.

Polymarket US facts this rests on (the official SDK polymarket-us 2.4.1,
websocket/base.py + markets.py, and a live read-only recon on the VM,
2026-11-05 23:17Z, polymarket_stream_recon.py):
  * The upgrade is signed like a REST GET of the path /v1/ws/markets
    (X-PM-Access-Key, X-PM-Timestamp, X-PM-Signature).
  * {"subscribe": {"requestId", "subscriptionType":
    "SUBSCRIPTION_TYPE_MARKET_DATA", "marketSlugs": [..]}} answers, per slug,
    at once, {"requestId", "subscriptionType", "marketData": {marketSlug,
    bids, offers, state, stats, transactTime}}: the WHOLE book, the REST
    /book's marketData shape (bids/offers of the YES contract, each {px:
    {value}, qty}). Every later change pushes the whole book again (live:
    3 frames in 38 s on a busy market, none on a still one). Several
    subscriptions with their own requestIds share one connection, and one
    may be added at any time. {"unsubscribe": {"requestId"}} ends one.
  * An error frame is {"requestId"?, "error", "subscriptionType"?}.
  * The SDK names a heartbeat frame; none came in 60 s live, so silence is
    no sign of a dead line. The stream pings every PING_SEC and drops a
    connection that answers nothing (no pong, no frame) for SILENT_RECONNECT_SEC.
  * One connection holds at most 10 subscriptions: the 11th answers
    "max subscriptions per connection reached" (live, 2026-11-05 23:36Z,
    the first deploy, which subscribed one market per requestId). A
    subscription names a list of slugs, so the markets are packed into at
    most MAX_SUBSCRIPTIONS (9) groups of SLUGS_PER_SUBSCRIPTION. Changing a
    group's markets means unsubscribing it and subscribing the new list under
    a new requestId, which answers every market's whole book again. No limit
    on slugs per subscription is published.

A book is served while the connection is up and answering, and its market's
last whole book came on this connection within BOOK_MAX_AGE_SEC. A still
market gets no frames, so each group is renewed (unsubscribe, subscribe the
same list under a new requestId) once it was sent RESUBSCRIBE_SEC ago.
Anything else answers an `error` naming why, never an old book. The Place
button's markets are placed in a group ahead of the rest, and the live
pass's are capped below the groups' room so Place always finds a place.

Side effects: one outbound websocket to Polymarket US on the configured key
(reads only: no verb here can place, cancel or move money). Nothing on disk;
books are held in memory only.
"""
import json
import logging
import threading
import time
from collections.abc import Callable
from typing import Protocol

from unabated_ticket.bets_service.novig_book_stream import websocket_client_connect
from unabated_ticket.bets_service.polymarket_orders import PolymarketOrderClient

log = logging.getLogger(__name__)

STREAM_PATH = "/v1/ws/markets"
MARKET_DATA_SUBSCRIPTION = "SUBSCRIPTION_TYPE_MARKET_DATA"
PING_SEC = 15
# No frame of any kind (a pong included) this long: the line is dead.
SILENT_RECONNECT_SEC = 45
RECV_POLL_SEC = 0.25
RECONNECT_MIN_SEC = 1
RECONNECT_MAX_SEC = 30
# Only a connection that held this long resets the backoff.
STABLE_CONNECTION_SEC = 60
REPEAT_FAILURE_LOG_SEC = 10 * 60
# A still market's subscription is renewed for a fresh whole book this often...
RESUBSCRIBE_SEC = 120
# ...and its book is not served once older than this (a renewal that went unanswered).
BOOK_MAX_AGE_SEC = 180
# A group with no book this long after it was sent is sent again; a market
# in an answered group with no book of its own this long reads as an error.
ANSWER_WAIT_SEC = 10
# A market Polymarket US refused is asked for again no sooner than this.
REFUSED_RETRY_SEC = 60
# A market nobody asked about for this long is unsubscribed (games end).
WATCH_IDLE_SEC = 10 * 60
# Polymarket US's cap on one connection's subscriptions. One is kept free, so
# a resubscribe never needs an 11th if Polymarket US frees the old one late.
SUBSCRIPTIONS_PER_CONNECTION = 10
MAX_SUBSCRIPTIONS = SUBSCRIPTIONS_PER_CONNECTION - 1
# The markets one subscription names.
SLUGS_PER_SUBSCRIPTION = 50
# The most a caller may wait for a new market's first book (the runner's own timeouts sit above it).
MAX_WAIT_SEC = 3.0
# The runner's live price pass may hold at most this many markets (of the
# groups' 450); past it a background market answers an error, so Place's
# market always has room.
MAX_BACKGROUND_MARKETS = 400


class Connection(Protocol):
    """One open websocket (novig_book_stream's websocket-client wrapper).
    recv returns a text message, "" for a control frame (a pong), or None when
    nothing arrived within the timeout; raises ConnectionError once it is gone."""

    def send(self, text: str) -> None: ...

    def recv(self, timeout_sec: float) -> str | None: ...

    def ping(self) -> None: ...

    def close(self) -> None: ...


Connect = Callable[[str, dict], Connection]


def polymarket_connect(url: str, headers: dict) -> Connection:
    return websocket_client_connect(url, headers, peer="Polymarket US")


def stream_url(api_url: str) -> str:
    """wss://host/v1/ws/markets from the REST base (https://host)."""
    base = api_url.rstrip("/")
    if base.startswith("https://"):
        return "wss://" + base[len("https://"):] + STREAM_PATH
    if base.startswith("http://"):
        return "ws://" + base[len("http://"):] + STREAM_PATH
    raise ValueError(f"expected an http(s) Polymarket US API URL, got {api_url!r}")


class PolymarketBookStream:
    """The books of the markets the runner prices, kept live on one websocket.
    run(stop) is the thread's loop; books_for() is called from HTTP threads."""

    def __init__(self, client: PolymarketOrderClient, url: str, connect: Connect = polymarket_connect,
                 clock: Callable[[], float] = time.monotonic):
        self._client = client
        self._url = url
        self._connect = connect
        self._clock = clock
        self._cond = threading.Condition()
        # Kept across connections: slug -> when last asked for, Place's, the live pass's alone, refusals.
        self._wanted: dict[str, float] = {}
        self._priority: set[str] = set()
        self._background: set[str] = set()
        self._refused: dict[str, tuple[float, str]] = {}
        self._request_seq = 0
        self._reset_connection_state()
        self._connection: Connection | None = None
        self._connected_at: float | None = None
        self._last_frame_at: float | None = None
        self._last_ping_at: float | None = None
        self._last_error: str | None = None
        self._last_logged_failure: tuple[str, float] | None = None

    def _reset_connection_state(self) -> None:
        # group index -> its live subscription on this connection: {requestId, slugs, sentAt, answeredAt}.
        self._groups: dict[int, dict] = {}
        # requestId -> group index, for this connection's live subscriptions.
        self._group_of_request: dict[str, int] = {}
        # slug -> group index it is (or is to be) subscribed in.
        self._group_of_slug: dict[str, int] = {}
        # slug -> (marketData, received at), this connection's books only.
        self._books: dict[str, tuple[dict, float]] = {}

    # ---- what the HTTP route calls ----

    def books_for(self, slugs: list[str], wait_sec: float, priority: bool = False, background: bool = False) -> dict:
        """{stream: status, books: {slug: view}}. Marks each market wanted
        (subscribing it if new; `priority` ahead of the rest, for the Place
        button; `background` behind the rest and capped, for the runner's live
        price pass, until Place asks for it too) and waits up to wait_sec
        (capped at MAX_WAIT_SEC) for every new market's first book. A view is
        {marketData, ageMs}: the REST /book's body plus how old the book is,
        or {error}."""
        wait_sec = max(0.0, min(wait_sec, MAX_WAIT_SEC))
        with self._cond:
            now = self._clock()
            capped: set[str] = set()
            for slug in slugs:
                if background and slug not in self._wanted:
                    if len(self._background) >= MAX_BACKGROUND_MARKETS:
                        capped.add(slug)
                        continue
                    self._background.add(slug)
                elif not background:
                    self._background.discard(slug)
                self._wanted[slug] = now
                if priority:
                    self._priority.add(slug)
            self._cond.notify_all()
            deadline = now + wait_sec
            while True:
                views = {slug: self._view_of(slug) for slug in slugs if slug not in capped}
                remaining = deadline - self._clock()
                if remaining <= 0 or not any(view.get("pending") for view in views.values()):
                    break
                self._cond.wait(remaining)
            for view in views.values():
                view.pop("pending", None)
            for slug in capped:
                views[slug] = {"error": f"the live price pass already watches {MAX_BACKGROUND_MARKETS} Polymarket US markets"}
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
            "subscriptions": len(self._groups),
            "subscribedMarkets": sum(len(group["slugs"]) for group in self._groups.values()),
            "booksHeld": len(self._books),
            "wantedMarkets": len(self._wanted),
        }

    def _view_of(self, slug: str) -> dict:
        if self._connection is None:
            return {"error": f"the Polymarket US stream is not connected ({self._last_error or 'connecting'})"}
        refused = self._refused.get(slug)
        if refused is not None:
            return {"error": f"Polymarket US refused the stream subscription to {slug}: {refused[1]}"}
        held = self._books.get(slug)
        if held is not None:
            market_data, received_at = held
            age_sec = self._clock() - received_at
            if age_sec > BOOK_MAX_AGE_SEC:
                return {"error": f"the stream's book of {slug} is {round(age_sec)} s old (its renewal went unanswered)"}
            return {"marketData": market_data, "ageMs": int(age_sec * 1000)}
        group = self._groups.get(self._group_of_slug.get(slug, -1))
        if group is not None and slug in group["slugs"] and group["answeredAt"] is not None \
                and self._clock() - group["answeredAt"] > ANSWER_WAIT_SEC:
            return {"error": f"Polymarket US's stream answered its subscription with no book of {slug}"}
        if slug in self._wanted and slug not in self._group_of_slug and self._groups_full():
            return {"error": f"the stream's {MAX_SUBSCRIPTIONS} subscriptions of {SLUGS_PER_SUBSCRIPTION} markets are full"}
        return {"pending": True, "error": f"waiting for Polymarket US's stream book of {slug}"}

    def _groups_full(self) -> bool:
        return len(self._group_of_slug) >= MAX_SUBSCRIPTIONS * SLUGS_PER_SUBSCRIPTION

    # ---- the stream thread ----

    def run(self, stop: threading.Event) -> None:
        """Connect, keep the wanted markets subscribed, take every book;
        reconnect with backoff whenever the line drops or stops answering."""
        backoff = RECONNECT_MIN_SEC
        while not stop.is_set():
            opened_at = None
            try:
                self._open()
                opened_at = self._clock()
                self._pump(stop)
            except (OSError, ConnectionError, ValueError, KeyError, TypeError) as error:
                self._drop(f"{type(error).__name__}: {error}")
            # Anything else too: a dead thread would leave the live prices and Place on REST until a restart.
            except Exception as error:
                log.exception("polymarket stream: unexpected error")
                self._drop(f"unexpected {type(error).__name__}: {error}")
            if opened_at is not None and self._clock() - opened_at > STABLE_CONNECTION_SEC:
                backoff = RECONNECT_MIN_SEC
            if stop.wait(backoff):
                break
            backoff = min(backoff * 2, RECONNECT_MAX_SEC)
        self._drop("stopped")

    def _open(self) -> None:
        headers = self._client.signed_headers("GET", STREAM_PATH)
        connection = self._connect(self._url, headers)
        with self._cond:
            self._reset_connection_state()
            self._connection = connection
            self._connected_at = self._last_frame_at = self._last_ping_at = self._clock()
            self._last_error = None
            self._last_logged_failure = None
            self._cond.notify_all()
        log.info("polymarket stream: connected to %s on key %s", self._url, self._client.key_id[:8])

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
            log.warning("polymarket stream: disconnected: %s", reason)
            return
        if reason == "stopped":
            return
        now = self._clock()
        last = self._last_logged_failure
        if last is None or last[0] != reason or now - last[1] > REPEAT_FAILURE_LOG_SEC:
            log.warning("polymarket stream: could not connect: %s", reason)
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
                    raise ConnectionError(f"no frame or pong from Polymarket US for {SILENT_RECONNECT_SEC} s")
                continue
            with self._cond:
                self._last_frame_at = now
                if text:
                    self._handle(json.loads(text))
                self._cond.notify_all()

    def _next_request_id(self) -> str:
        self._request_seq += 1
        return f"pm-{self._request_seq}"

    def _outgoing_frames(self) -> list[dict]:
        """Drop idle and refused markets from their groups, put new markets in
        a group with room (Place's first, the live pass's last), then for every
        group whose markets changed, whose subscribe went unanswered or whose
        renewal is due: unsubscribe its old requestId and subscribe its list
        under a new one."""
        with self._cond:
            now = self._clock()
            for slug, (refused_at, _reason) in list(self._refused.items()):
                if now - refused_at > REFUSED_RETRY_SEC:
                    del self._refused[slug]
            self._forget_idle(now)
            for slug in [slug for slug in self._group_of_slug if slug in self._refused]:
                del self._group_of_slug[slug]
                self._books.pop(slug, None)
            self._assign_new_markets()
            return self._group_frames(now)

    def _forget_idle(self, now: float) -> None:
        for slug in [slug for slug, asked_at in self._wanted.items() if now - asked_at > WATCH_IDLE_SEC]:
            del self._wanted[slug]
            self._priority.discard(slug)
            self._background.discard(slug)
            self._group_of_slug.pop(slug, None)
            self._books.pop(slug, None)

    def _assign_new_markets(self) -> None:
        new = [slug for slug in self._wanted if slug not in self._group_of_slug and slug not in self._refused]
        new.sort(key=lambda slug: (slug not in self._priority, slug in self._background, slug))
        sizes = {index: 0 for index in range(MAX_SUBSCRIPTIONS)}
        for index in self._group_of_slug.values():
            sizes[index] += 1
        for slug in new:
            room = [index for index in range(MAX_SUBSCRIPTIONS) if sizes[index] < SLUGS_PER_SUBSCRIPTION]
            if not room:
                return
            # The emptiest group, so one new market resends as few books as it can.
            index = min(room, key=lambda candidate: (sizes[candidate], candidate))
            self._group_of_slug[slug] = index
            sizes[index] += 1

    def _group_frames(self, now: float) -> list[dict]:
        frames: list[dict] = []
        wanted_slugs: dict[int, list[str]] = {}
        for slug, index in self._group_of_slug.items():
            wanted_slugs.setdefault(index, []).append(slug)
        for index in range(MAX_SUBSCRIPTIONS):
            slugs = sorted(wanted_slugs.get(index, []))
            group = self._groups.get(index)
            if group is not None:
                unchanged = group["slugs"] == slugs
                unanswered = group["answeredAt"] is None and now - group["sentAt"] > ANSWER_WAIT_SEC
                renewal_due = now - group["sentAt"] > RESUBSCRIBE_SEC
                if unchanged and not unanswered and not renewal_due:
                    continue
                frames.append({"unsubscribe": {"requestId": group["requestId"]}})
                del self._group_of_request[group["requestId"]]
                del self._groups[index]
            if not slugs:
                continue
            request_id = self._next_request_id()
            self._groups[index] = {"requestId": request_id, "slugs": slugs, "sentAt": now, "answeredAt": None}
            self._group_of_request[request_id] = index
            frames.append({"subscribe": {"requestId": request_id, "subscriptionType": MARKET_DATA_SUBSCRIPTION,
                                         "marketSlugs": slugs}})
        return frames

    # ---- frames from Polymarket US (called with the lock held) ----

    def _handle(self, message: object) -> None:
        if not isinstance(message, dict):
            raise TypeError(f"expected a JSON object frame from Polymarket US, got {str(message)[:200]}")
        if "heartbeat" in message:
            return
        if "error" in message:
            self._handle_error(message)
            return
        market_data = message.get("marketData")
        if market_data is None:
            # Anything else (a lite price, a trade, an unsubscribe answer) carries no whole book.
            return
        request_id = message.get("requestId")
        index = self._group_of_request.get(request_id) if isinstance(request_id, str) else None
        if index is None:
            # An old requestId's last frame after its group was resubscribed: the new one answers on its own.
            return
        group = self._groups[index]
        group["answeredAt"] = group["answeredAt"] or self._clock()
        slug = market_data.get("marketSlug") if isinstance(market_data, dict) else None
        if slug not in group["slugs"]:
            log.warning("polymarket stream: a book frame for %s, which subscription %s does not name (%s)",
                        slug, request_id, json.dumps(market_data, default=str)[:300])
            return
        self._books[slug] = (market_data, self._clock())

    def _handle_error(self, message: dict) -> None:
        """A refused group: its markets are refused for REFUSED_RETRY_SEC (the
        reason names the group's size, so a cap on slugs per subscription shows)."""
        reason = str(message.get("error"))[:300]
        request_id = message.get("requestId")
        index = self._group_of_request.pop(request_id, None) if isinstance(request_id, str) else None
        if index is None:
            log.warning("polymarket stream: Polymarket US answered an error: %s (requestId %s)", reason, request_id)
            return
        group = self._groups.pop(index)
        now = self._clock()
        for slug in group["slugs"]:
            self._refused[slug] = (now, reason)
            self._group_of_slug.pop(slug, None)
            self._books.pop(slug, None)
        log.error("polymarket stream: subscription %s of %d market(s) refused: %s", request_id, len(group["slugs"]), reason)
