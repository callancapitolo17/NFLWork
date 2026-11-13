"""Kalshi RFQ feed reader for the dry-run quoter. Reads only; sends nothing
but the subscribe command.

Why its own reader (2026-11-13): the feed is ~9,900 frames a second (half
rfq_created, half rfq_deleted). The websocket-client library the MLB maker's
reader uses spends ~76 us a frame in pure-Python frame handling, kept up with
~4,200 frames a second at 93% of a core, and Kalshi cut the socket every
15-20 s ("buffer overflow"). The `websockets` library read the whole feed at
42% of a core with no disconnects (60 s live runs).

One daemon thread: connect -> subscribe -> recv loop. Each frame goes through
door(): a byte check drops ~95% of frames before any JSON decode, then the
full scope rule and the size floor run, so only RFQs the quoter could price
reach the queue. Everything else is counted, never stored. A dead or silent
socket reconnects with backoff; RFQs created while it is down are not
recovered (dry run: the minute counts show `reconnects`).
"""
import json
import logging
import queue
import random
import re
import threading
import time
from collections import Counter
from typing import Callable, Protocol

from kalshi_parlay_quoter import scope

log = logging.getLogger("kalshi_parlay_quoter.feed")

COMMUNICATIONS_CHANNEL = "communications"
SUBSCRIBE_CMD_ID = 1
SUBSCRIBE_ACK_TIMEOUT_SEC = 10.0
RECV_TIMEOUT_SEC = 5.0
RECONNECT_BASE_SEC = 1.0
RECONNECT_MAX_SEC = 30.0

# A leg ticker of one of scope.SERIES as it appears in a frame, e.g.
# "KXNFLSPREAD-. A frame without one cannot be in scope.
_SERIES_IN_FRAME = re.compile(
    b'"(?:' + b"|".join(re.escape(s.encode()) for s in sorted(scope.SERIES)) + b")-")
_CREATED = b'"rfq_created"'
_DELETED = b'"rfq_deleted"'

KEPT = "kept"


def below_floor(rfq: dict, min_cost_usd: float) -> bool:
    """A dollar-sized RFQ under the floor, or a contract-sized one with fewer
    contracts than the floor's dollars (a contract pays at most $1, so its
    cost is under the floor whatever the price)."""
    if min_cost_usd <= 0:
        return False
    try:
        cost = float(rfq.get("target_cost_dollars") or 0)
        contracts = float(rfq.get("contracts_fp") or 0)
    except (TypeError, ValueError):
        return False
    if cost > 0:
        return cost < min_cost_usd
    return 0 < contracts < min_cost_usd


def door(raw: bytes, min_cost_usd: float) -> tuple[str, dict | None]:
    """(count key, rfq). rfq is the frame's msg when the key is KEPT, the
    decoded envelope when the key is frames_other (acks and errors), else None.

    Keys: frames_deleted, frames_other, malformed, door_other_league (no leg
    in a series we price; legless frames land here too), door_no_legs,
    door_unsupported_leg, door_bad_side, door_single_leg, door_same_game,
    door_below_floor, kept."""
    if _CREATED not in raw:
        if _DELETED in raw:
            return "frames_deleted", None
        try:
            envelope = json.loads(raw)
        except ValueError:
            return "malformed", None
        return "frames_other", envelope if isinstance(envelope, dict) else None
    if not _SERIES_IN_FRAME.search(raw):
        return "door_other_league", None
    try:
        envelope = json.loads(raw)
    except ValueError:
        return "malformed", None
    rfq = envelope.get("msg") if isinstance(envelope, dict) else None
    if (not isinstance(rfq, dict) or envelope.get("type") != "rfq_created"
            or not isinstance(rfq.get("id"), str) or not rfq["id"]):
        return "malformed", None
    legs = rfq.get("mve_selected_legs")
    if not isinstance(legs, list) or not all(isinstance(leg, dict) for leg in legs):
        return "door_no_legs", None
    verdict = scope.evaluate(legs)
    if not verdict.in_scope:
        return f"door_{verdict.reason}", None
    if below_floor(rfq, min_cost_usd):
        return "door_below_floor", None
    return KEPT, rfq


class Connection(Protocol):
    """What the reader needs of a websocket. recv raises TimeoutError when
    nothing arrived in time and any other exception when the socket is gone."""
    def send(self, text: str) -> None: ...
    def recv(self, timeout: float) -> bytes: ...
    def close(self) -> None: ...


class _WebsocketsConnection:
    def __init__(self, url: str, headers: dict[str, str]):
        from websockets.sync.client import connect   # lazy: tests need no network lib
        # compression off: Kalshi has always sent this client plain frames,
        # and inflating ~10k frames a second would cost what the swap saved.
        self._ws = connect(url, additional_headers=headers, compression=None,
                           max_size=2 ** 22, open_timeout=10,
                           ping_interval=20, ping_timeout=20)

    def send(self, text: str) -> None:
        self._ws.send(text)

    def recv(self, timeout: float) -> bytes:
        return self._ws.recv(timeout=timeout, decode=False)

    def close(self) -> None:
        self._ws.close()


def kalshi_connect(url: str, headers_fn: Callable[[], dict[str, str]]) -> Callable[[], Connection]:
    """A connect() for RFQFeed: a fresh signed handshake per connection."""
    return lambda: _WebsocketsConnection(url, headers_fn())


class RFQFeed:
    """Kept RFQs land in `rfqs` (queue.Queue of the frame's msg dict). When the
    queue is full the RFQ is counted `queue_dropped` and dropped, so a stalled
    consumer can never stall the socket."""

    def __init__(self, connect: Callable[[], Connection], *, min_cost_usd: float,
                 queue_max: int, idle_timeout_sec: float,
                 clock: Callable[[], float] = time.monotonic):
        self._connect = connect
        self._min_cost_usd = min_cost_usd
        self._idle_timeout_sec = idle_timeout_sec
        self._clock = clock
        self.rfqs: queue.Queue = queue.Queue(maxsize=queue_max)
        self._counts: Counter = Counter()
        self._lock = threading.Lock()
        self._stop = threading.Event()
        self._paused = threading.Event()
        self._thread: threading.Thread | None = None
        self.connected = False

    def start(self) -> None:
        if self._thread is None:
            self._thread = threading.Thread(target=self._run, name="rfq-feed", daemon=True)
            self._thread.start()

    def stop(self) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=RECV_TIMEOUT_SEC + 2)

    def set_paused(self, paused: bool) -> None:
        """Paused, the reader closes its socket and reads nothing (the stop
        file must stop the whole feed, not only the pricing)."""
        if paused:
            self._paused.set()
        else:
            self._paused.clear()

    def take_counts(self) -> dict[str, int]:
        """Counts since the last call."""
        with self._lock:
            counts, self._counts = self._counts, Counter()
        return dict(counts)

    def _count(self, key: str) -> None:
        with self._lock:
            self._counts[key] += 1

    def handle_frame(self, raw: bytes) -> dict | None:
        """Count one frame and queue it when kept. Returns the decoded
        envelope of a non-RFQ frame (ack, error), else None."""
        key, payload = door(raw, self._min_cost_usd)
        if key == KEPT:
            try:
                self.rfqs.put_nowait(payload)
            except queue.Full:
                key = "queue_dropped"
        self._count(key)
        return payload if key == "frames_other" else None

    def _read_until_dead(self, conn: Connection) -> None:
        conn.send(json.dumps({"id": SUBSCRIBE_CMD_ID, "cmd": "subscribe",
                              "params": {"channels": [COMMUNICATIONS_CHANNEL]}}))
        subscribed = False
        opened_at = last_frame_at = self._clock()
        while not self._stop.is_set() and not self._paused.is_set():
            now = self._clock()
            if not subscribed and now - opened_at > SUBSCRIBE_ACK_TIMEOUT_SEC:
                raise ConnectionError(f"subscribe unacknowledged after {SUBSCRIBE_ACK_TIMEOUT_SEC:.0f}s")
            try:
                raw = conn.recv(RECV_TIMEOUT_SEC)
            except TimeoutError:
                if self._clock() - last_frame_at > self._idle_timeout_sec:
                    raise ConnectionError(
                        f"no frames in {self._idle_timeout_sec:.0f}s, socket presumed dead") from None
                continue
            last_frame_at = self._clock()
            if isinstance(raw, str):
                raw = raw.encode()
            envelope = self.handle_frame(raw)
            if envelope is None:
                continue
            if envelope.get("type") == "error":
                raise ConnectionError(f"server error frame: {envelope.get('msg')}")
            if envelope.get("type") == "subscribed" and not subscribed:
                subscribed = True
                self.connected = True
                log.info("[feed] connected, %s subscription acked", COMMUNICATIONS_CHANNEL)

    def _run(self) -> None:
        backoff = RECONNECT_BASE_SEC
        while not self._stop.is_set():
            if self._paused.is_set():
                self._stop.wait(1.0)
                continue
            conn = None
            try:
                conn = self._connect()
                self._read_until_dead(conn)
                backoff = RECONNECT_BASE_SEC
            except Exception as e:                      # any failure reconnects
                if self._stop.is_set():
                    break
                if self.connected:
                    backoff = RECONNECT_BASE_SEC        # it had been working
                self._count("reconnects")
                log.error("[feed] connection lost or failed: %s: %s", type(e).__name__, e)
            finally:
                self.connected = False
                if conn is not None:
                    try:
                        conn.close()
                    except Exception as e:
                        log.warning("[feed] close failed: %s", e)
            if self._stop.is_set() or self._paused.is_set():
                continue
            self._stop.wait(min(backoff, RECONNECT_MAX_SEC) + random.uniform(0, backoff * 0.25))
            backoff = min(backoff * 2, RECONNECT_MAX_SEC)
