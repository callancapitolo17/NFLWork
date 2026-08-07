"""WebSocket RFQ discovery (issue #56): RFQSource adapter over Kalshi's
`communications` channel.

The discovery tick's feed model is LEVEL-TRIGGERED: `poll()` must return ALL
open RFQs every tick, because that re-entry drives on-demand re-fetch,
quorum refinement, and pulls — an RFQ that vanishes from poll() output stops
being fetched and refined by design. A WebSocket stream is EDGE-TRIGGERED
(rfq_created / rfq_deleted events), so this adapter is a stateful mirror:

  - a daemon reader thread owns the connection (connect → subscribe → recv
    loop → reconnect with backoff+jitter) and mutates the open-RFQ set under
    a lock; the raw `msg` payload is stored as-is, so poll() output is
    byte-compatible with RestRFQSource (id, market_ticker, creator_id,
    contracts_fp as STRING, target_cost_dollars — all verbatim from Kalshi);
  - poll() runs on the tick thread, never blocks on the network, and snapshots
    the mirror. A last-frame watchdog inside poll() (so it fires even if the
    reader thread is wedged) trips into a LOUD automatic REST fallback
    (ERROR log + notify + `rfq_source_transition` research event) and flips
    back the same way on recovery — never silently deaf;
  - on every (re)connect the reader performs ONE REST gap-fill AFTER
    subscribing: wholesale reconcile of the mirror (additions AND removals,
    so a missed rfq_deleted cannot leave a zombie being re-fetched forever).
    Events buffered on the socket during the gap-fill re-apply on top, so no
    create/delete between subscribe and gap-fill is lost;
  - get_market() stays pure REST (delegated to the injected rest source).

Measurement: every WS rfq_created emits `rfq_seen_latency` (created_ts from
the payload vs. receipt wall clock) to the research firehose.

Side effects: no DuckDB writes of its own — research events go through the
buffered firehose, notify goes through kalshi_mlb_mm.notify.
"""
import json
import logging
import random
import threading
import time
from datetime import datetime, timezone
from typing import Callable, Protocol

from kalshi_common import auth_client
from kalshi_mlb_mm import config, notify, research
from kalshi_mlb_mm.rfq_source import RFQSource

log = logging.getLogger("kalshi_mlb_mm.ws_rfq_source")

COMMUNICATIONS_CHANNEL = "communications"
SUBSCRIBE_CMD_ID = 1


class TransportTimeout(Exception):
    """recv() saw nothing within its idle timeout — the connection may still
    be alive (liveness is judged by the watchdog, not the transport)."""


class WSTransport(Protocol):
    """Minimal blocking transport the reader thread drives. recv() returns a
    text frame, None for non-message activity (server pings), raises
    TransportTimeout on idle, and raises any other exception when dead."""

    def connect(self) -> None: ...
    def send(self, text: str) -> None: ...
    def recv(self) -> str | None: ...
    def close(self) -> None: ...


class KalshiWSTransport:
    """Real transport over the `websocket-client` package.

    Imported lazily so REST mode and the test suite never need the dependency
    installed. Auth follows the REST pattern: signed headers from
    kalshi_common.auth_client (configure()-injected credentials, no env reads
    here). Server pings (~10s cadence) are answered by recv_data and surfaced
    as None so the caller can count them as liveness.
    """

    def __init__(self, *, url: str, recv_timeout_sec: float = 5.0,
                 handshake_timeout_sec: float = 10.0):
        self._url = url
        self._recv_timeout_sec = recv_timeout_sec
        self._handshake_timeout_sec = handshake_timeout_sec
        self._ws = None

    def connect(self) -> None:
        import websocket  # websocket-client — lazy, see class docstring
        headers = auth_client.ws_auth_headers()
        ws = websocket.WebSocket()
        ws.connect(self._url,
                   header=[f"{k}: {v}" for k, v in headers.items()],
                   timeout=self._handshake_timeout_sec)
        ws.settimeout(self._recv_timeout_sec)
        self._ws = ws

    def send(self, text: str) -> None:
        self._ws.send(text)

    def recv(self) -> str | None:
        import websocket
        try:
            # control_frame=True: pings reach us (recv_data already pongs
            # them) so heartbeat frames count as liveness upstream.
            opcode, data = self._ws.recv_data(control_frame=True)
        except websocket.WebSocketTimeoutException:
            raise TransportTimeout()
        if opcode == 1:  # ABNF.OPCODE_TEXT
            return data.decode() if isinstance(data, bytes) else data
        return None  # ping/pong/binary — activity, not a message

    def close(self) -> None:
        if self._ws is not None:
            self._ws.close()


def _parse_created_ts(value) -> datetime | None:
    """Kalshi timestamps: ISO-8601 string ("...Z") per the AsyncAPI spec;
    tolerate epoch seconds in case the shape drifts."""
    if isinstance(value, (int, float)):
        return datetime.fromtimestamp(float(value), tz=timezone.utc)
    if isinstance(value, str):
        return datetime.fromisoformat(value.replace("Z", "+00:00"))
    return None


class WebSocketRFQSource:
    """RFQSource over the WS `communications` channel with loud REST fallback.

    Thread model: `start()` launches one daemon reader thread; `poll()` and
    `get_market()` are called from the main tick thread only. `_lock` guards
    the mirror + liveness fields shared between the two.
    """

    def __init__(self, *, rest_source: RFQSource,
                 transport_factory: Callable[[], WSTransport] | None = None,
                 heartbeat_timeout_sec: float | None = None,
                 reconnect_base_sec: float | None = None,
                 reconnect_max_sec: float | None = None,
                 connect_alert_streak: int | None = None):
        self._rest = rest_source
        self._transport_factory = transport_factory or (
            lambda: KalshiWSTransport(url=config.KALSHI_WS_URL))
        self._heartbeat_timeout_sec = (
            config.WS_HEARTBEAT_TIMEOUT_SEC if heartbeat_timeout_sec is None
            else heartbeat_timeout_sec)
        self._reconnect_base_sec = (
            config.WS_RECONNECT_BASE_SEC if reconnect_base_sec is None
            else reconnect_base_sec)
        self._reconnect_max_sec = (
            config.WS_RECONNECT_MAX_SEC if reconnect_max_sec is None
            else reconnect_max_sec)
        self._connect_alert_streak = (
            config.WS_CONNECT_ALERT_STREAK if connect_alert_streak is None
            else connect_alert_streak)
        self._lock = threading.Lock()
        self._open_rfqs: dict[str, dict] = {}
        self._ws_ready = False          # connect + subscribe + gap-fill done
        self._last_msg_at = 0.0         # any WS frame, pings included
        self._connect_failures = 0
        self._outage_notified = False   # one notify pair per outage
        # poll()-side serving mode, tick thread only. "rest_startup" = WS has
        # never served (quiet REST until the feed proves healthy);
        # "ws" = serving the mirror; "rest_fallback" = tripped, loud.
        self._mode = "rest_startup"
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self.malformed_count = 0
        self.gap_fill_count = 0

    # -- public RFQSource interface ------------------------------------------

    @property
    def ws_ready(self) -> bool:
        with self._lock:
            return self._ws_ready

    def start(self) -> None:
        if self._thread is not None:
            return
        self._thread = threading.Thread(target=self._run,
                                        name="ws-rfq-reader", daemon=True)
        self._thread.start()

    def stop(self) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=2.0)

    def poll(self) -> list[dict]:
        now = time.time()
        with self._lock:
            ws_ready = self._ws_ready
            last_msg_age = now - self._last_msg_at
        healthy = ws_ready and last_msg_age <= self._heartbeat_timeout_sec
        if healthy and self._mode != "ws":
            self._mode = "ws"
            log.info("[ws_rfq] WS feed healthy — discovery served from WS mirror")
            research.emit("rfq_source_transition",
                          payload=dict(to="ws", last_msg_age_sec=last_msg_age))
            with self._lock:
                was_notified = self._outage_notified
                self._outage_notified = False
            if was_notified:
                notify.resume("ws_rfq_source", detail="WS RFQ feed recovered")
        elif not healthy and (self._mode == "ws"
                              or (self._mode == "rest_startup" and ws_ready)):
            # Never silently deaf: the trip is evaluated HERE on the tick
            # thread so a wedged reader thread still degrades loudly. The
            # rest_startup+ws_ready arm catches a socket that connected and
            # then went silent BEFORE its first healthy poll — without it
            # that outage would serve REST forever without a peep.
            self._mode = "rest_fallback"
            log.error("[ws_rfq] WS feed stale/dead (ready=%s, last_msg_age=%.1fs)"
                      " — automatic REST fallback active", ws_ready, last_msg_age)
            research.emit("rfq_source_transition",
                          payload=dict(to="rest_fallback", ws_ready=ws_ready,
                                       last_msg_age_sec=last_msg_age))
            self._notify_outage_once(
                f"WS RFQ feed stale (last frame {last_msg_age:.0f}s ago); "
                "discovery served by REST fallback")
        if self._mode != "ws":
            return self._rest.poll()
        with self._lock:
            return list(self._open_rfqs.values())

    def get_market(self, market_ticker: str) -> dict | None:
        return self._rest.get_market(market_ticker)

    # -- reader thread -------------------------------------------------------

    def _run(self) -> None:
        backoff = self._reconnect_base_sec
        while not self._stop.is_set():
            transport = None
            try:
                transport = self._transport_factory()
                transport.connect()
                transport.send(json.dumps(
                    {"id": SUBSCRIBE_CMD_ID, "cmd": "subscribe",
                     "params": {"channels": [COMMUNICATIONS_CHANNEL]}}))
                # Gap-fill AFTER subscribe: events created in the meantime
                # buffer on the socket and re-apply on top of the snapshot.
                self._gap_fill("connect")
                with self._lock:
                    self._ws_ready = True
                    self._last_msg_at = time.time()
                    self._connect_failures = 0
                backoff = self._reconnect_base_sec
                log.info("[ws_rfq] connected, subscribed to %s, gap-filled",
                         COMMUNICATIONS_CHANNEL)
                while not self._stop.is_set():
                    try:
                        raw = transport.recv()
                    except TransportTimeout:
                        continue  # idle — watchdog in poll() judges liveness
                    with self._lock:
                        self._last_msg_at = time.time()
                    if raw:
                        self._handle_message(raw)
            except Exception as e:
                if self._stop.is_set():
                    break
                with self._lock:
                    self._ws_ready = False
                    self._connect_failures += 1
                    failures = self._connect_failures
                log.error("[ws_rfq] connection lost/failed (streak=%d): %s",
                          failures, e)
                if failures == self._connect_alert_streak:
                    self._notify_outage_once(
                        f"{failures} consecutive WS connect failures; "
                        "discovery served by REST fallback")
            finally:
                if transport is not None:
                    try:
                        transport.close()
                    except Exception:
                        pass
            delay = (min(backoff, self._reconnect_max_sec)
                     + random.uniform(0, backoff * 0.25))
            backoff = min(backoff * 2, self._reconnect_max_sec)
            self._stop.wait(delay)

    def _notify_outage_once(self, detail: str) -> None:
        with self._lock:
            if self._outage_notified:
                return
            self._outage_notified = True
        notify.halt("ws_rfq_source", detail=detail)

    def _gap_fill(self, reason: str) -> None:
        """ONE REST poll reconciling the mirror wholesale — additions AND
        removals. REST is authoritative for what is open at this instant; a
        REST hiccup (RestRFQSource returns [] on error) can transiently empty
        the mirror, which loses quoting opportunities until the next event or
        reconnect but can never resurrect a dead RFQ — the safe direction."""
        rest_rfqs = self._rest.poll()
        snapshot = {rfq["id"]: rfq for rfq in rest_rfqs if rfq.get("id")}
        with self._lock:
            had = len(self._open_rfqs)
            self._open_rfqs.clear()
            self._open_rfqs.update(snapshot)
        self.gap_fill_count += 1
        if had and not snapshot:
            log.warning("[ws_rfq] gap-fill (%s) emptied a %d-entry mirror — "
                        "REST reported no open RFQs (or errored)", reason, had)
        else:
            log.info("[ws_rfq] gap-fill (%s): %d open RFQs from REST",
                     reason, len(snapshot))
        research.emit("rfq_ws_gap_fill",
                      payload=dict(reason=reason, open_count=len(snapshot),
                                   prior_count=had))

    def _handle_message(self, raw: str) -> None:
        try:
            envelope = json.loads(raw)
            if not isinstance(envelope, dict):
                raise ValueError(f"expected object envelope, got {type(envelope).__name__}")
            msg_type = envelope.get("type")
            if msg_type not in ("rfq_created", "rfq_deleted"):
                return  # subscribed ack, quote_* lifecycle, errors — not ours
            msg = envelope.get("msg")
            if not isinstance(msg, dict):
                raise ValueError(f"{msg_type} without object msg payload")
            rfq_id = msg.get("id")
            if not rfq_id or not isinstance(rfq_id, str):
                raise ValueError(f"{msg_type} missing rfq id")
            if msg_type == "rfq_created":
                with self._lock:
                    self._open_rfqs[rfq_id] = msg
                self._emit_seen_latency(msg)
            else:
                with self._lock:
                    self._open_rfqs.pop(rfq_id, None)
        except Exception as e:
            self.malformed_count += 1
            log.warning("[ws_rfq] malformed payload skipped (count=%d): %s",
                        self.malformed_count, e)

    def _emit_seen_latency(self, msg: dict) -> None:
        """Issue #56 measurement: created→seen latency per RFQ, so the WS
        discovery win over the 2s REST poll is quantified in the firehose."""
        seen = datetime.now(timezone.utc)
        created_raw = msg.get("created_ts")
        latency_ms = None
        try:
            created = _parse_created_ts(created_raw)
            if created is not None:
                latency_ms = (seen - created).total_seconds() * 1000.0
        except (ValueError, TypeError, OSError):
            pass  # unparseable created_ts — latency stays None, event still lands
        research.emit("rfq_seen_latency", rfq_id=msg.get("id"),
                      ticker=msg.get("market_ticker"),
                      payload=dict(created_ts=created_raw,
                                   seen_ts=seen.isoformat(),
                                   latency_ms=latency_ms, source="ws"))
