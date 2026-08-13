"""Regression tests for issue #56: WebSocket RFQ discovery.

WebSocketRFQSource must mirror Kalshi's `communications` channel into a
stateful open-RFQ set so `poll()` keeps the level-triggered contract the
discovery tick depends on (all open RFQs, every tick). These tests drive the
adapter through a scripted in-memory transport (the "stubbed WS server") —
no real sockets, no websocket-client dependency.
"""
import inspect
import json
import queue
import threading
import time

import pytest

from kalshi_mlb_mm import config, notify, research
from kalshi_mlb_mm import ws_rfq_source
from kalshi_mlb_mm.rfq_source import RestRFQSource
from kalshi_mlb_mm.ws_rfq_source import TransportTimeout, WebSocketRFQSource


# --- test doubles -----------------------------------------------------------

def _rfq(rid: str, ticker: str = "KXMVE-TEST-1", **extra) -> dict:
    base = {"id": rid, "creator_id": "", "market_ticker": ticker,
            "event_ticker": ticker.rsplit("-", 1)[0],
            # $500: above the option-B door floor so generic feed tests pass the door
            "contracts_fp": "0.00", "target_cost_dollars": "500.00",
            "created_ts": "2026-08-07T18:00:00Z"}
    base.update(extra)
    return base


class FakeTransport:
    """Scripted transport: recv() blocks on an inbox queue with a short
    timeout so the reader thread stays responsive to stop(). Acks the
    subscription on connect (like the real server) unless auto_ack=False.
    Set connect_exc to make connect() fail (e.g. a DNS resolver error)."""

    def __init__(self, auto_ack: bool = True):
        self.inbox: queue.Queue = queue.Queue()
        self.sent: list[str] = []
        self.connect_calls = 0
        self.closed = False
        self.auto_ack = auto_ack
        self.connect_exc: Exception | None = None
        # While set, recv() blocks — simulates a WEDGED reader thread (stuck
        # inside the transport), which is the one failure the reader's own
        # staleness teardown cannot see; only the poll() watchdog can.
        self.wedged = threading.Event()

    def connect(self) -> None:
        self.connect_calls += 1
        if self.connect_exc is not None:
            raise self.connect_exc
        if self.auto_ack:
            self.push({"id": 1, "type": "subscribed",
                       "msg": {"channel": "communications", "sid": 1}})

    def send(self, text: str) -> None:
        self.sent.append(text)

    def recv(self) -> str | None:
        while self.wedged.is_set():
            time.sleep(0.005)
        try:
            item = self.inbox.get(timeout=0.02)
        except queue.Empty:
            raise TransportTimeout()
        if isinstance(item, Exception):
            raise item
        return item

    def close(self) -> None:
        self.closed = True

    # helpers for tests
    def push(self, obj) -> None:
        self.inbox.put(json.dumps(obj) if isinstance(obj, dict) else obj)

    def drop(self) -> None:
        self.inbox.put(ConnectionError("connection dropped"))


class FakeTransportFactory:
    """Hands out one FakeTransport per connect cycle; remembers them all."""

    def __init__(self, count: int = 4, auto_ack: bool = True):
        self.transports = [FakeTransport(auto_ack=auto_ack) for _ in range(count)]
        self.created = 0

    def __call__(self):
        t = self.transports[self.created]
        self.created += 1
        return t


class FakeRestSource:
    """RFQSource stand-in for the gap-fill / fallback path."""

    def __init__(self):
        self.poll_results: list[list[dict]] = [[]]
        self.poll_calls = 0
        self.get_market_calls = []

    def poll(self) -> list[dict]:
        self.poll_calls += 1
        idx = min(self.poll_calls - 1, len(self.poll_results) - 1)
        return [dict(r) for r in self.poll_results[idx]]

    def get_market(self, market_ticker: str) -> dict | None:
        self.get_market_calls.append(market_ticker)
        return {"ticker": market_ticker}


def wait_until(predicate, timeout: float = 2.0) -> bool:
    deadline = time.time() + timeout
    while time.time() < deadline:
        if predicate():
            return True
        time.sleep(0.005)
    return predicate()


@pytest.fixture
def capture(monkeypatch):
    """Record notify + research calls; never hit webhooks or the firehose."""
    events: list[tuple] = []
    monkeypatch.setattr(notify, "halt",
                        lambda reason, detail="": events.append(("halt", reason, detail)))
    monkeypatch.setattr(notify, "resume",
                        lambda reason, detail="": events.append(("resume", reason, detail)))
    monkeypatch.setattr(research, "emit",
                        lambda event_type, **kw: events.append(("emit", event_type, kw)))
    return events


def make_source(factory, rest, **kw):
    defaults = dict(rest_source=rest, transport_factory=factory,
                    heartbeat_timeout_sec=5.0,
                    reconnect_base_sec=0.01, reconnect_max_sec=0.05,
                    connect_alert_streak=3)
    defaults.update(kw)
    return WebSocketRFQSource(**defaults)


@pytest.fixture
def running():
    """Track sources so every test stops its reader thread."""
    sources = []
    yield sources
    for s in sources:
        s.stop()


# --- feed semantics ---------------------------------------------------------

def test_create_event_appears_in_poll(capture, running):
    factory, rest = FakeTransportFactory(), FakeRestSource()
    src = make_source(factory, rest)
    running.append(src)
    src.start()
    assert wait_until(lambda: src.ws_ready)
    factory.transports[0].push({"type": "rfq_created", "sid": 1, "seq": 1,
                                "msg": _rfq("r1")})
    assert wait_until(lambda: any(r["id"] == "r1" for r in src.poll()))
    out = src.poll()
    assert len(out) == 1 and out[0]["market_ticker"] == "KXMVE-TEST-1"


def test_expire_event_removes_from_poll(capture, running):
    factory, rest = FakeTransportFactory(), FakeRestSource()
    src = make_source(factory, rest)
    running.append(src)
    src.start()
    assert wait_until(lambda: src.ws_ready)
    t = factory.transports[0]
    t.push({"type": "rfq_created", "sid": 1, "msg": _rfq("r1")})
    assert wait_until(lambda: len(src.poll()) == 1)
    t.push({"type": "rfq_deleted", "sid": 1, "msg": _rfq("r1")})
    assert wait_until(lambda: src.poll() == [])


def test_contracts_fp_stays_string_and_parses(capture, running):
    """Byte-compat: the WS dict must round-trip through the tick's size
    helper exactly like a REST dict (contracts_fp is a fixed-point STRING)."""
    from kalshi_mlb_mm.main import _rfq_requested_contracts
    factory, rest = FakeTransportFactory(), FakeRestSource()
    src = make_source(factory, rest)
    running.append(src)
    src.start()
    assert wait_until(lambda: src.ws_ready)
    factory.transports[0].push({"type": "rfq_created", "sid": 1,
                                "msg": _rfq("r1", contracts_fp="21.00")})
    assert wait_until(lambda: len(src.poll()) == 1)
    rfq = src.poll()[0]
    assert rfq["contracts_fp"] == "21.00"
    assert isinstance(rfq["contracts_fp"], str)
    assert _rfq_requested_contracts(rfq) == 21.0


def test_non_rfq_message_types_ignored(capture, running):
    factory, rest = FakeTransportFactory(), FakeRestSource()
    src = make_source(factory, rest)
    running.append(src)
    src.start()
    assert wait_until(lambda: src.ws_ready)
    t = factory.transports[0]
    t.push({"id": 1, "type": "subscribed", "msg": {"channel": "communications", "sid": 1}})
    t.push({"type": "quote_created", "sid": 1, "msg": {"id": "q1"}})
    t.push({"type": "rfq_created", "sid": 1, "msg": _rfq("r1")})
    assert wait_until(lambda: len(src.poll()) == 1)
    assert src.malformed_count == 0


# --- gap-fill ---------------------------------------------------------------

def test_gap_fill_on_connect_seeds_mirror(capture, running):
    factory, rest = FakeTransportFactory(), FakeRestSource()
    rest.poll_results = [[_rfq("pre1"), _rfq("pre2", ticker="KXMVE-TEST-2")]]
    src = make_source(factory, rest)
    running.append(src)
    src.start()
    assert wait_until(lambda: src.ws_ready)
    ids = {r["id"] for r in src.poll()}
    assert ids == {"pre1", "pre2"}
    # WS mode must not touch REST on subsequent polls (only the 1 gap-fill).
    assert rest.poll_calls == 1
    src.poll()
    assert rest.poll_calls == 1


def test_reconnect_gap_fill_reconciles_additions_and_removals(capture, running):
    """Dropped connection → reconnect + exactly ONE REST gap-fill that both
    ADDS RFQs created during the gap and REMOVES ones that expired (a missed
    rfq_deleted must not leave a zombie in the mirror forever)."""
    factory, rest = FakeTransportFactory(), FakeRestSource()
    # 1st connect: empty book. 2nd connect: r1 vanished during the gap
    # (missed delete) and r2 appeared (missed create).
    rest.poll_results = [[], [_rfq("r2")]]
    src = make_source(factory, rest)
    running.append(src)
    src.start()
    assert wait_until(lambda: src.ws_ready)
    factory.transports[0].push({"type": "rfq_created", "sid": 1, "msg": _rfq("r1")})
    assert wait_until(lambda: len(src.poll()) == 1)
    factory.transports[0].drop()
    assert wait_until(lambda: factory.created == 2 and src.ws_ready)
    assert rest.poll_calls == 2
    ids = {r["id"] for r in src.poll()}
    assert ids == {"r2"}, "gap-fill must reconcile removals AND additions"
    # subscribe re-sent on the new connection
    assert any("subscribe" in s for s in factory.transports[1].sent)


# --- subscription ack gate --------------------------------------------------

def test_unacked_subscribe_forces_reconnect(capture, running):
    """A live socket whose subscribe is silently dropped must NOT count as
    ready — server pings would keep the heartbeat watchdog quiet forever."""
    factory, rest = FakeTransportFactory(auto_ack=False), FakeRestSource()
    # 2nd transport acks normally, so recovery is also proven.
    factory.transports[1].auto_ack = True
    src = make_source(factory, rest, subscribe_ack_timeout_sec=0.05)
    running.append(src)
    src.start()
    assert wait_until(lambda: factory.created >= 2)
    assert wait_until(lambda: src.ws_ready)
    assert rest.poll_calls == 1  # gap-fill only ran for the ACKED connection


def test_error_frame_forces_reconnect(capture, running):
    factory, rest = FakeTransportFactory(), FakeRestSource()
    src = make_source(factory, rest)
    running.append(src)
    src.start()
    assert wait_until(lambda: src.ws_ready)
    factory.transports[0].push({"id": 6, "type": "error",
                                "msg": {"code": 9, "msg": "subscription dropped"}})
    assert wait_until(lambda: factory.created == 2 and src.ws_ready)
    assert rest.poll_calls == 2  # one gap-fill per (re)connect
    assert src.malformed_count == 0  # an error frame is fatal, not malformed


def test_connect_failure_streak_alerts_once(capture, running):
    class DeadFactory:
        def __init__(self):
            self.calls = 0

        def __call__(self):
            self.calls += 1
            raise ConnectionError("refused")

    factory, rest = DeadFactory(), FakeRestSource()
    src = make_source(factory, rest, connect_alert_streak=3)
    running.append(src)
    src.start()
    assert wait_until(lambda: factory.calls >= 5)
    halts = [e for e in capture if e[0] == "halt"]
    assert len(halts) == 1 and halts[0][1] == "ws_rfq_source"
    assert src.poll() == []  # REST fallback keeps serving (rest is empty)
    assert rest.poll_calls >= 1


# --- watchdog / fallback ----------------------------------------------------

def test_silent_stale_ws_falls_back_to_rest_loudly(capture, running):
    """A WEDGED reader (stuck inside the transport — cannot run its own
    staleness teardown) must still degrade loudly: the poll()-side watchdog
    trips into REST fallback. A merely-dead socket is the reader's own job
    (test_silent_dead_socket_forces_reconnect)."""
    factory, rest = FakeTransportFactory(), FakeRestSource()
    rest.poll_results = [[], [_rfq("from_rest")]]
    src = make_source(factory, rest, heartbeat_timeout_sec=0.05)
    running.append(src)
    src.start()
    assert wait_until(lambda: src.ws_ready)
    assert src.poll() == []  # healthy WS, empty book
    factory.transports[0].wedged.set()
    time.sleep(0.15)  # exceed heartbeat timeout with zero frames
    out = src.poll()
    assert [r["id"] for r in out] == ["from_rest"], "stale WS must serve REST"
    halts = [e for e in capture if e[0] == "halt"]
    assert len(halts) == 1 and halts[0][1] == "ws_rfq_source"
    transitions = [e for e in capture
                   if e[0] == "emit" and e[1] == "rfq_source_transition"]
    assert transitions and transitions[-1][2]["payload"]["to"] == "rest_fallback"
    factory.transports[0].wedged.clear()  # let stop() join the reader promptly


def test_recovery_flips_back_to_ws_and_notifies(capture, running):
    factory, rest = FakeTransportFactory(), FakeRestSource()
    rest.poll_results = [[_rfq("seed")]]
    src = make_source(factory, rest, heartbeat_timeout_sec=0.05)
    running.append(src)
    src.start()
    assert wait_until(lambda: src.ws_ready)
    t = factory.transports[0]
    t.wedged.set()  # wedged reader: only the poll() watchdog can see this
    time.sleep(0.15)
    src.poll()  # trips watchdog → REST fallback
    assert any(e[0] == "halt" for e in capture)
    # Any WS frame refreshes liveness (server pings count as activity too).
    t.push({"type": "subscribed", "msg": {"sid": 1}})
    t.wedged.clear()
    assert wait_until(
        lambda: src.poll() and src.poll()[0]["id"] == "seed" and
        any(e[0] == "resume" for e in capture))
    transitions = [e for e in capture
                   if e[0] == "emit" and e[1] == "rfq_source_transition"]
    assert transitions[-1][2]["payload"]["to"] == "ws"


def test_fallback_notifies_once_per_outage(capture, running):
    factory, rest = FakeTransportFactory(), FakeRestSource()
    src = make_source(factory, rest, heartbeat_timeout_sec=0.05)
    running.append(src)
    src.start()
    assert wait_until(lambda: src.ws_ready)
    factory.transports[0].wedged.set()
    time.sleep(0.15)
    for _ in range(5):
        src.poll()
    assert len([e for e in capture if e[0] == "halt"]) == 1
    factory.transports[0].wedged.clear()  # let stop() join the reader promptly


# --- silent dead socket (2026-08-12 incident) --------------------------------

def test_silent_dead_socket_forces_reconnect(capture, running):
    """Live incident 2026-08-12 ~18:59: a machine-wide network blip killed
    the TCP connection WITHOUT any close/error frame reaching us. recv()
    then raises only TransportTimeout forever, which the reader treated as
    benign idle — subscribed stayed True, ws_ready stayed True, and the
    reader never reconnected for 13+ hours (the poll() watchdog only flips
    the serving mode; it cannot tear down the socket). Kalshi pings every
    ~10s, so a subscribed connection with no frames for a full heartbeat
    timeout is dead: the reader must tear down and reconnect on its own."""
    factory, rest = FakeTransportFactory(), FakeRestSource()
    src = make_source(factory, rest, heartbeat_timeout_sec=0.1)
    running.append(src)
    src.start()
    assert wait_until(lambda: src.ws_ready)
    # No frames ever arrive on transport 0 after the ack — silence, not an
    # exception. The reader must notice staleness and reconnect by itself.
    assert wait_until(lambda: factory.created >= 2), \
        "reader must reconnect off a silently-dead socket"
    assert factory.transports[0].closed
    assert wait_until(lambda: src.ws_ready)  # transport 1 acked + gap-filled
    # Transport 1 is also silent, so further teardowns may land before this
    # line on a slow machine — pin the invariant, not a fixed count. Reading
    # poll_calls BEFORE created keeps the bound race-free (gap-fills only
    # happen after their transport is created, so poll_calls <= created).
    polls = rest.poll_calls
    assert 2 <= polls <= factory.created  # one gap-fill per (re)connect


def test_frames_within_heartbeat_do_not_reconnect(capture, running):
    """The staleness teardown must not fire on a healthy-but-quiet feed:
    any frame (server pings surface as recv() None) refreshes liveness."""
    factory, rest = FakeTransportFactory(), FakeRestSource()
    src = make_source(factory, rest, heartbeat_timeout_sec=0.3)
    running.append(src)
    src.start()
    assert wait_until(lambda: src.ws_ready)
    for _ in range(6):  # 0.6s of sparse traffic, each gap < heartbeat
        time.sleep(0.1)
        factory.transports[0].push({"type": "quote_created", "sid": 1,
                                    "msg": {"id": "q1"}})
    time.sleep(0.1)
    assert factory.created == 1, "no reconnect while frames keep arriving"
    assert src.ws_ready


def test_dns_failure_during_reconnect_retries_until_recovery(capture, running):
    """Same incident, hypothesized path (disproven but pinned here): a DNS
    outage during reconnect raises OSError from connect(). Every such
    failure must be caught, logged, and retried with backoff until the
    resolver recovers — never allowed to kill the reader thread."""
    factory, rest = FakeTransportFactory(count=5), FakeRestSource()
    dns_error = OSError("nodename nor servname provided, or not known")
    factory.transports[1].connect_exc = dns_error
    factory.transports[2].connect_exc = dns_error
    factory.transports[3].connect_exc = dns_error
    src = make_source(factory, rest)
    running.append(src)
    src.start()
    assert wait_until(lambda: src.ws_ready)
    factory.transports[0].drop()
    # Three DNS-failed attempts, then transport 4 connects and acks.
    assert wait_until(lambda: factory.created == 5 and src.ws_ready)
    assert factory.transports[4].connect_calls == 1
    # The alert streak (3) was hit exactly once during the outage, and the
    # first healthy poll() after recovery sends the matching resume.
    halts = [e for e in capture if e[0] == "halt"]
    assert len(halts) == 1 and halts[0][1] == "ws_rfq_source"
    src.poll()
    assert any(e[0] == "resume" for e in capture)


# --- malformed payloads -----------------------------------------------------

def test_malformed_payload_counted_skipped_no_crash(capture, running):
    factory, rest = FakeTransportFactory(), FakeRestSource()
    src = make_source(factory, rest)
    running.append(src)
    src.start()
    assert wait_until(lambda: src.ws_ready)
    t = factory.transports[0]
    t.push("this is not json")
    t.push({"type": "rfq_created", "sid": 1, "msg": {"no_id_field": True}})
    t.push({"type": "rfq_created", "sid": 1, "msg": _rfq("good")})
    assert wait_until(lambda: any(r["id"] == "good" for r in src.poll()))
    assert src.malformed_count == 2
    assert len(src.poll()) == 1  # reader thread survived both


# --- latency measurement ----------------------------------------------------

def test_seen_latency_emitted_per_rfq(capture, running):
    factory, rest = FakeTransportFactory(), FakeRestSource()
    src = make_source(factory, rest)
    running.append(src)
    src.start()
    assert wait_until(lambda: src.ws_ready)
    factory.transports[0].push({"type": "rfq_created", "sid": 1, "msg": _rfq("r1")})
    assert wait_until(lambda: any(
        e[0] == "emit" and e[1] == "rfq_seen_latency" for e in capture))
    ev = [e for e in capture if e[0] == "emit" and e[1] == "rfq_seen_latency"][-1]
    payload = ev[2]["payload"]
    assert ev[2]["rfq_id"] == "r1"
    assert payload["created_ts"] == "2026-08-07T18:00:00Z"
    assert payload["latency_ms"] is not None and payload["latency_ms"] > 0


# --- source construction (main.build_rfq_source) ----------------------------

def test_build_rfq_source_is_ws_only_with_rest_fallback_injected(monkeypatch):
    """WS-only by user decision 2026-08-07: no mode switch exists. The REST
    source survives only as the injected automatic-fallback/gap-fill path."""
    from kalshi_mlb_mm import main
    started = []

    class StubWS:
        def __init__(self, *, rest_source):
            self.rest_source = rest_source

        def start(self):
            started.append(self)

    monkeypatch.setattr(ws_rfq_source, "WebSocketRFQSource", StubWS)
    src = main.build_rfq_source()
    assert isinstance(src, StubWS)
    assert started == [src]
    assert isinstance(src.rest_source, RestRFQSource)
    assert not hasattr(config, "RFQ_SOURCE")  # the knob must stay deleted


# --- interface discipline ---------------------------------------------------

def test_new_params_are_keyword_only():
    from kalshi_common import auth_client
    init_params = inspect.signature(WebSocketRFQSource.__init__).parameters
    for name, p in init_params.items():
        if name == "self":
            continue
        assert p.kind == inspect.Parameter.KEYWORD_ONLY, f"{name} must be keyword-only"
    ws_hdr_params = inspect.signature(auth_client.ws_auth_headers).parameters
    for name, p in ws_hdr_params.items():
        assert p.kind == inspect.Parameter.KEYWORD_ONLY, f"{name} must be keyword-only"


def test_get_market_delegates_to_rest(capture, running):
    factory, rest = FakeTransportFactory(), FakeRestSource()
    src = make_source(factory, rest)
    running.append(src)
    m = src.get_market("KXMVE-TEST-1")
    assert m == {"ticker": "KXMVE-TEST-1"}
    assert rest.get_market_calls == ["KXMVE-TEST-1"]


# --- ingestion filter + mirror TTL (2026-08-10) ------------------------------

def _mlb_legs():
    return [{"event_ticker": "KXMLBGAME-26AUG10NYMATL",
             "market_ticker": "KXMLBGAME-26AUG10NYMATL-ATL", "side": "yes"},
            {"event_ticker": "KXMLBTOTAL-26AUG10NYMATL",
             "market_ticker": "KXMLBTOTAL-26AUG10NYMATL-T8.5", "side": "yes"}]


def _nba_legs():
    return [{"event_ticker": "KXNBAGAME-X",
             "market_ticker": "KXNBAGAME-X-LAL", "side": "yes"}]


def test_ingestion_filter_drops_non_mlb_keeps_mlb_and_legless(capture, running):
    """The communications channel is the whole exchange's firehose
    (~5-10k creates/min at peak) and expirations never arrive as
    rfq_deleted — everything not plausibly MLB must die at the door.
    Legless frames pass (fail-open: the tick's fallback classifies)."""
    factory, rest = FakeTransportFactory(), FakeRestSource()
    src = make_source(factory, rest)
    running.append(src)
    src.start()
    assert wait_until(lambda: src.ws_ready)
    t = factory.transports[0]
    t.push({"type": "rfq_created", "sid": 1,
            "msg": _rfq("mlb1", mve_selected_legs=_mlb_legs())})
    t.push({"type": "rfq_created", "sid": 1,
            "msg": _rfq("nba1", mve_selected_legs=_nba_legs())})
    t.push({"type": "rfq_created", "sid": 1, "msg": _rfq("legless1")})
    assert wait_until(lambda: len(src.poll()) == 2)
    ids = {r["id"] for r in src.poll()}
    assert ids == {"mlb1", "legless1"}
    assert src.ingestion_dropped == 1


def test_mirror_ttl_evicts_expired_entries(capture, running):
    """Kalshi sends no rfq_deleted on expiry, so without a TTL the mirror
    grows forever (48k entries in 5 min observed live). Entries past the
    TTL disappear from poll()."""
    factory, rest = FakeTransportFactory(), FakeRestSource()
    src = make_source(factory, rest, mirror_ttl_sec=0.2)
    running.append(src)
    src.start()
    assert wait_until(lambda: src.ws_ready)
    factory.transports[0].push({"type": "rfq_created", "sid": 1,
                                "msg": _rfq("r1", mve_selected_legs=_mlb_legs())})
    assert wait_until(lambda: len(src.poll()) == 1)
    time.sleep(0.3)
    assert src.poll() == [], "TTL-expired entry must leave the mirror"


def test_gap_fill_applies_ingestion_filter(capture, running):
    """REST gap-fill snapshots go through the same door: a reconcile must
    not resurrect the non-MLB firehose into the mirror."""
    rest = FakeRestSource()
    rest.poll_results = [[_rfq("mlb1", mve_selected_legs=_mlb_legs()),
                          _rfq("nba1", mve_selected_legs=_nba_legs())]]
    factory = FakeTransportFactory()
    src = make_source(factory, rest)
    running.append(src)
    src.start()
    assert wait_until(lambda: src.ws_ready)
    assert wait_until(lambda: {r["id"] for r in src.poll()} == {"mlb1"})


def test_door_floor_and_leg_count_gates(capture, running):
    """Option B (2026-08-11): with no fetch ceiling, the door's dollar floor
    and leg-count cap are the only limits on book traffic. Below-floor and
    4+-leg RFQs die at the door; contracts-denominated RFQs (target_cost 0)
    pass the floor and meet the tick's size gate instead."""
    factory, rest = FakeTransportFactory(), FakeRestSource()
    src = make_source(factory, rest, door_max_legs=3, door_min_cost_usd=250.0)
    running.append(src)
    src.start()
    assert wait_until(lambda: src.ws_ready)
    t = factory.transports[0]
    four_legs = _mlb_legs() + [
        {"event_ticker": "KXMLBGAME-A", "market_ticker": "KXMLBGAME-A-AAA",
         "side": "yes"},
        {"event_ticker": "KXMLBGAME-B", "market_ticker": "KXMLBGAME-B-BBB",
         "side": "yes"}]
    t.push({"type": "rfq_created", "sid": 1,
            "msg": _rfq("big", mve_selected_legs=_mlb_legs(),
                        target_cost_dollars="300.00")})
    t.push({"type": "rfq_created", "sid": 1,
            "msg": _rfq("small", mve_selected_legs=_mlb_legs(),
                        target_cost_dollars="25.00")})
    t.push({"type": "rfq_created", "sid": 1,
            "msg": _rfq("basket", mve_selected_legs=four_legs,
                        target_cost_dollars="900.00")})
    t.push({"type": "rfq_created", "sid": 1,
            "msg": _rfq("contracts", mve_selected_legs=_mlb_legs(),
                        target_cost_dollars="0.00", contracts_fp="40.00")})
    assert wait_until(lambda: {r["id"] for r in src.poll()} == {"big", "contracts"})
    assert src.ingestion_drops["below_floor"] == 1
    assert src.ingestion_drops["too_many_legs"] == 1
    assert src.ingestion_kept == 2
