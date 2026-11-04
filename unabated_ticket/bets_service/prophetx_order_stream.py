"""ProphetX's live order updates over ONE websocket on the account's own key
pair (the ProphetX source's client): ProphetX pushes an `orders` event on the
account's private channel when one of its orders is activated, fills or is
cancelled. The Place desk (prophetx_place.py) waits on it after a send, so a
fill is reported the moment ProphetX matches it instead of on the next poll of
the order history (Cal 2026-11-04, option A; on 2026-11-03 a match took up to
~15 s and the order history, read at 1 request/s shared with the bet source,
was seen late).

ProphetX facts this rests on (docs.prophetx.co "Integrate with Websockets",
"WebSocket Events", "Various Websocket Use Cases", read 2026-11-04; not a live
capture, so every `orders` frame is logged until its shape is confirmed):
  * The transport speaks the Pusher wire protocol (ProphetX's "Talaria").
    GET /partner/websocket/connection-config -> {app_id, cluster, key,
    service, ws_host}; connect to ws_host (else the Pusher cluster host).
    The URL form wss://<host>/app/<key>?protocol=7 is Pusher's standard one,
    not stated in ProphetX's docs.
  * pusher:connection_established carries socket_id (its `data` is JSON
    encoded inside JSON). POST /partner/v4/mm/websocket {socket_id, service,
    subscriptions} with header X-CFTC-Terminology: true answers
    {authenticated: {auth, user_data}, authorized_channel: [{channel_name,
    auth}]}; then pusher:signin with `authenticated`, and pusher:subscribe to
    the private user channel (private-service=...-user=<account id>) with its
    own `auth`. pusher:ping is answered with pusher:pong.
  * Events: {timestamp, change_type, payload: base64 JSON, op}. `orders` on
    the private channel decodes to {"info": {id, external_id, status,
    matching_status, filled_quantity, total_filled_quantity, open_quantity,
    price, fill_price, fee, ...}}. `filled_quantity` is read as the latest
    fill and `total_filled_quantity` as the order's running total (the docs
    do not define either), so the total is what Place counts.
  * health_check arrives every 5 s; a socket silent for over 30 s is dead.

What it holds: the newest pushed state of each order seen in the last
ORDER_KEEP_SEC, in memory only. Side effects: one outbound websocket and,
per connect, two REST calls on the shared client (config + register); reads
only, nothing here can place, cancel or move money. Nothing on disk.
"""
import base64
import functools
import json
import logging
import threading
import time
from collections.abc import Callable

from unabated_ticket.bets_service.novig_book_stream import Connect, Connection, websocket_client_connect
from unabated_ticket.bets_service.prophetx_client import ProphetXClient

log = logging.getLogger(__name__)

PUSHER_PROTOCOL = 7
PUSHER_CLIENT_QUERY = f"protocol={PUSHER_PROTOCOL}&client=nflwork-bets&version=1.0"
ORDERS_EVENT = "orders"
HEALTH_EVENT = "health_check"
PRIVATE_USER_CHANNEL_MARK = "user="
# ProphetX sends health_check every 5 s; the docs say reconnect after 30 s of silence.
SILENT_RECONNECT_SEC = 30
HANDSHAKE_TIMEOUT_SEC = 15
RECV_POLL_SEC = 0.5
RECONNECT_MIN_SEC = 2
RECONNECT_MAX_SEC = 120
STABLE_CONNECTION_SEC = 60
# A pushed order state is kept this long after its last update (Place waits at most ~20 s).
ORDER_KEEP_SEC = 600
LOGGED_FRAME_CHARS = 600


def pusher_url(config: dict) -> str:
    """The socket URL from GET /websocket/connection-config."""
    host = str(config.get("ws_host") or "").strip()
    if not host:
        cluster = str(config.get("cluster") or "").strip()
        if not cluster:
            raise ValueError(f"connection-config names neither ws_host nor cluster: {str(config)[:200]}")
        host = f"ws-{cluster}.pusher.com"
    for scheme in ("wss://", "https://"):
        if host.startswith(scheme):
            host = host[len(scheme):]
    return f"wss://{host.rstrip('/')}/app/{config['key']}?{PUSHER_CLIENT_QUERY}"


def parse_frame(text: str) -> tuple[str, str | None, object]:
    """(event, channel, data) of one Pusher frame; `data` JSON-decoded when it is a JSON string."""
    frame = json.loads(text)
    if not isinstance(frame, dict) or not isinstance(frame.get("event"), str):
        raise ValueError(f"not a Pusher frame: {text[:200]}")
    data = frame.get("data")
    if isinstance(data, str):
        try:
            data = json.loads(data)
        except json.JSONDecodeError:
            pass
    return frame["event"], frame.get("channel"), data


def order_of_event(data: object) -> dict | None:
    """The order an `orders` event carries, as Place reads orders: `order_id`
    set from ProphetX's `id`, and `filled_quantity` the running total. None
    when the event holds no order id (logged by the caller)."""
    payload = data.get("payload") if isinstance(data, dict) else None
    if isinstance(payload, str):
        payload = json.loads(base64.b64decode(payload))
    if isinstance(payload, dict) and isinstance(payload.get("info"), dict):
        payload = payload["info"]
    if not isinstance(payload, dict):
        return None
    order_id = payload.get("id") or payload.get("order_id")
    if not order_id:
        return None
    order = dict(payload)
    order["order_id"] = str(order_id)
    if order.get("total_filled_quantity") is not None:
        order["filled_quantity"] = order["total_filled_quantity"]
    return order


class ProphetXOrderStream:
    """The account's pushed order updates. Thread-safe: run() on its own
    thread, Place calls wait_for_order() from request threads."""

    def __init__(self, client: ProphetXClient, connect: Connect = functools.partial(websocket_client_connect, peer="ProphetX"),
                 clock: Callable[[], float] = time.monotonic):
        self._client = client
        self._connect = connect
        self._clock = clock
        self._cond = threading.Condition()
        self._connection: Connection | None = None
        self._live = False
        self._connected_at: float | None = None
        self._last_frame_at: float | None = None
        self._last_error: str | None = None
        self._orders: dict[str, tuple[float, int, dict]] = {}
        self._version = 0
        self._orders_events = 0

    # ---- what Place reads ----

    def is_live(self) -> bool:
        """Signed in, subscribed, and heard from within SILENT_RECONNECT_SEC."""
        with self._cond:
            return self._live_locked()

    def _live_locked(self) -> bool:
        return (self._live and self._last_frame_at is not None
                and self._clock() - self._last_frame_at <= SILENT_RECONNECT_SEC)

    def wait_for_order(self, order_id: str, after_version: int, timeout_sec: float) -> tuple[dict | None, int]:
        """The order's newest pushed state once it is newer than `after_version`
        (or None at the timeout), and its version to pass next time."""
        deadline = self._clock() + timeout_sec
        with self._cond:
            while True:
                held = self._orders.get(order_id)
                if held is not None and held[1] > after_version:
                    return dict(held[2]), held[1]
                remaining = deadline - self._clock()
                if remaining <= 0 or not self._live_locked():
                    return None, after_version
                self._cond.wait(min(remaining, RECV_POLL_SEC))

    def status(self) -> dict:
        with self._cond:
            now = self._clock()
            return {"live": self._live_locked(),
                    "connectedForSec": None if self._connected_at is None else round(now - self._connected_at),
                    "lastFrameAgeMs": None if self._last_frame_at is None else int((now - self._last_frame_at) * 1000),
                    "ordersEvents": self._orders_events, "error": self._last_error}

    # ---- the stream thread ----

    def run(self, stop: threading.Event) -> None:
        """Connect, sign in, follow the private channel; reconnect with
        backoff whenever the socket drops or falls silent."""
        backoff = RECONNECT_MIN_SEC
        while not stop.is_set():
            opened_at = None
            try:
                self._open()
                opened_at = self._clock()
                self._pump(stop)
            except Exception as error:  # a dead thread would leave Place on the slow path until a restart
                self._drop(f"{type(error).__name__}: {error}")
            if opened_at is not None and self._clock() - opened_at > STABLE_CONNECTION_SEC:
                backoff = RECONNECT_MIN_SEC
            if stop.wait(backoff):
                break
            backoff = min(backoff * 2, RECONNECT_MAX_SEC)
        self._drop("stopped")

    def _open(self) -> None:
        config = self._client.websocket_config()
        url = pusher_url(config)
        connection = self._connect(url, {})
        with self._cond:
            self._connection = connection
            self._connected_at = self._last_frame_at = self._clock()
        established = self._await_event(connection, "pusher:connection_established")
        socket_id = established.get("socket_id") if isinstance(established, dict) else None
        if not socket_id:
            raise ValueError(f"connection_established without socket_id: {str(established)[:200]}")
        registered = self._client.register_websocket(str(socket_id), config.get("service"))
        authenticated = registered.get("authenticated") if isinstance(registered, dict) else None
        if not isinstance(authenticated, dict) or not authenticated.get("auth"):
            raise ValueError(f"websocket register: no authenticated block in {str(registered)[:300]}")
        connection.send(json.dumps({"event": "pusher:signin", "data": {
            "auth": authenticated["auth"], "user_data": authenticated.get("user_data")}}))
        self._await_event(connection, "pusher:signin_success")
        channel = _private_user_channel(registered.get("authorized_channel"))
        connection.send(json.dumps({"event": "pusher:subscribe", "data": {
            "channel": channel["channel_name"], "auth": channel.get("auth")}}))
        self._await_event(connection, "pusher_internal:subscription_succeeded")
        with self._cond:
            self._live = True
            self._last_error = None
            self._cond.notify_all()
        log.info("prophetx stream: live on %s (%s)", url.split("?")[0], channel["channel_name"])

    def _await_event(self, connection: Connection, wanted: str) -> object:
        """Read frames until `wanted` (answering pings), within HANDSHAKE_TIMEOUT_SEC."""
        deadline = self._clock() + HANDSHAKE_TIMEOUT_SEC
        while self._clock() < deadline:
            text = connection.recv(RECV_POLL_SEC)
            if not text:
                continue
            self._heard()
            event, _channel, data = parse_frame(text)
            if event == wanted:
                return data
            if event == "pusher:error" or event.endswith("subscription_error"):
                raise ConnectionError(f"ProphetX's socket answered {event}: {str(data)[:300]}")
            self._handle(connection, event, data)
        raise TimeoutError(f"no {wanted} within {HANDSHAKE_TIMEOUT_SEC} s")

    def _pump(self, stop: threading.Event) -> None:
        connection = self._connection
        while not stop.is_set():
            text = connection.recv(RECV_POLL_SEC)
            if text:
                self._heard()
                event, _channel, data = parse_frame(text)
                self._handle(connection, event, data)
            with self._cond:
                silent_sec = self._clock() - self._last_frame_at
            if silent_sec > SILENT_RECONNECT_SEC:
                raise TimeoutError(f"ProphetX's socket sent nothing for {round(silent_sec)} s")

    def _heard(self) -> None:
        with self._cond:
            self._last_frame_at = self._clock()

    def _handle(self, connection: Connection, event: str, data: object) -> None:
        if event == "pusher:ping":
            connection.send(json.dumps({"event": "pusher:pong", "data": {}}))
            return
        if event != ORDERS_EVENT:
            return
        # Logged whole until the shape is confirmed live: the docs describe no example.
        log.info("prophetx stream: orders event %s", json.dumps(data)[:LOGGED_FRAME_CHARS])
        try:
            order = order_of_event(data)
        except (ValueError, TypeError) as error:
            log.warning("prophetx stream: unreadable orders event (%s): %s", error, str(data)[:LOGGED_FRAME_CHARS])
            return
        if order is None:
            log.warning("prophetx stream: orders event names no order: %s", str(data)[:LOGGED_FRAME_CHARS])
            return
        now = self._clock()
        with self._cond:
            self._orders_events += 1
            self._version += 1
            self._orders[order["order_id"]] = (now, self._version, order)
            for order_id in [key for key, held in self._orders.items() if now - held[0] > ORDER_KEEP_SEC]:
                del self._orders[order_id]
            self._cond.notify_all()

    def _drop(self, reason: str) -> None:
        with self._cond:
            connection = self._connection
            was_live = self._live
            repeated = reason == self._last_error
            self._connection = None
            self._live = False
            self._connected_at = None
            if reason != "stopped":
                self._last_error = reason
            self._cond.notify_all()
        if connection is not None:
            connection.close()
        if reason == "stopped":
            return
        if was_live:
            log.warning("prophetx stream: down: %s", reason)
        elif not repeated:
            # The same connect failure every backoff is logged once; Place reads the order history meanwhile.
            log.error("prophetx stream: cannot connect (logged once until it changes): %s", reason)


def _private_user_channel(channels: object) -> dict:
    """The account's private channel among the registration's authorized channels."""
    for channel in channels if isinstance(channels, list) else []:
        name = channel.get("channel_name") if isinstance(channel, dict) else None
        if isinstance(name, str) and name.startswith("private-") and PRIVATE_USER_CHANNEL_MARK in name:
            return channel
    raise ValueError(f"websocket register: no private user channel among {str(channels)[:300]}")
