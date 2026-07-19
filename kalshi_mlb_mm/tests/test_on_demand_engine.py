"""Phase 2 OnDemandEngine: queue/worker/result store, poll-driven feed."""
import threading
import time

import pytest

from kalshi_common import legset
from mlb_sgp._shared import GameRef, OnDemandBookResult
from kalshi_mlb_mm.on_demand import OnDemandEngine, QUOTE_FRESH_SEC

EVT = "KXMLBGAME-25JUN271905NYYBOS"
GAME = GameRef(game_id="g1", home_team="Boston Red Sox",
               away_team="New York Yankees", commence_time=None)


def _legs():
    return [legset.CanonicalLeg(EVT, "spread", -1.5, "home"),
            legset.CanonicalLeg(EVT, "total", 8.5, "over"),
            legset.CanonicalLeg(EVT, "ml", None, "home")]


class FakeService:
    """Duck-typed SGPService: .books + .price_on_demand."""
    def __init__(self, books=("draftkings", "fanduel"), fair=0.14,
                 fail_books=(), latency=0.0):
        self.books = tuple(books)
        self.fair = fair
        self.fail_books = set(fail_books)
        self.latency = latency
        self.calls = []          # (book, hash) log
        self.lock = threading.Lock()

    def price_on_demand(self, book, game, legs):
        with self.lock:
            self.calls.append((book, legset.leg_set_hash(legs)))
        if self.latency:
            time.sleep(self.latency)
        if book in self.fail_books:
            return None
        return OnDemandBookResult(book=book, fair=self.fair, route="partition",
                                  n_cells_priced=8, latency_sec=0.1)


class FakeClock:
    def __init__(self):
        self.t = 1000.0
    def __call__(self):
        return self.t


def test_ensure_fetch_dedups_inflight_and_drain_lands_result():
    clock = FakeClock()
    svc = FakeService()
    eng = OnDemandEngine(svc, now_fn=clock, autostart=False)
    legs = _legs()
    h = legset.leg_set_hash(legs)
    eng.ensure_fetch(h, GAME, legs)
    eng.ensure_fetch(h, GAME, legs)          # dedup while queued
    assert eng._queue_len() == 1
    eng._drain_once()
    fairs = eng.lookup(h)
    assert fairs == {"draftkings": 0.14, "fanduel": 0.14}
    # both books called exactly once (shared fetch, not per-RFQ duplication)
    assert sorted(b for b, _ in svc.calls) == ["draftkings", "fanduel"]


def test_lookup_dies_after_fresh_window_and_refeeds():
    clock = FakeClock()
    svc = FakeService()
    eng = OnDemandEngine(svc, now_fn=clock, autostart=False)
    legs = _legs()
    h = legset.leg_set_hash(legs)
    eng.ensure_fetch(h, GAME, legs)
    eng._drain_once()
    assert eng.lookup(h) is not None
    clock.t += QUOTE_FRESH_SEC + 0.1
    assert eng.lookup(h) is None                      # dead for quoting
    eng.ensure_fetch(h, GAME, legs)                   # poll re-feeds
    assert eng._queue_len() == 1
    eng._drain_once()
    assert eng.lookup(h) is not None


def test_failed_books_dropped_and_all_failed_means_no_lookup():
    clock = FakeClock()
    svc = FakeService(fail_books=("fanduel",))
    eng = OnDemandEngine(svc, now_fn=clock, autostart=False)
    legs = _legs()
    h = legset.leg_set_hash(legs)
    eng.ensure_fetch(h, GAME, legs)
    eng._drain_once()
    assert eng.lookup(h) == {"draftkings": 0.14}
    svc2 = FakeService(fail_books=("draftkings", "fanduel"))
    eng2 = OnDemandEngine(svc2, now_fn=clock, autostart=False)
    eng2.ensure_fetch(h, GAME, legs)
    eng2._drain_once()
    assert eng2.lookup(h) is None                     # negative landing
    # a failed landing must not satisfy the next tick: poll re-feeds
    eng2.ensure_fetch(h, GAME, legs)
    assert eng2._queue_len() == 1


def test_lookup_results_exposes_routes_for_research():
    clock = FakeClock()
    eng = OnDemandEngine(FakeService(), now_fn=clock, autostart=False)
    legs = _legs()
    h = legset.leg_set_hash(legs)
    eng.ensure_fetch(h, GAME, legs)
    eng._drain_once()
    res = eng.lookup_results(h)
    assert res["draftkings"].route == "partition"


def test_worker_thread_processes_queue():
    svc = FakeService()
    eng = OnDemandEngine(svc)                          # real thread, real clock
    legs = _legs()
    h = legset.leg_set_hash(legs)
    eng.ensure_fetch(h, GAME, legs)
    deadline = time.monotonic() + 5.0
    while eng.lookup(h) is None and time.monotonic() < deadline:
        time.sleep(0.02)
    assert eng.lookup(h) is not None


def test_worker_survives_service_exceptions():
    class ExplodingService(FakeService):
        def price_on_demand(self, book, game, legs):
            raise RuntimeError("boom")
    eng = OnDemandEngine(ExplodingService())
    legs = _legs()
    h = legset.leg_set_hash(legs)
    eng.ensure_fetch(h, GAME, legs)
    deadline = time.monotonic() + 5.0
    while eng._inflight_len() and time.monotonic() < deadline:
        time.sleep(0.02)
    assert eng.lookup(h) is None                       # failed, not crashed
    # engine still serves new work after the explosion
    eng._service = FakeService()
    eng.ensure_fetch(h, GAME, legs)
    deadline = time.monotonic() + 5.0
    while eng.lookup(h) is None and time.monotonic() < deadline:
        time.sleep(0.02)
    assert eng.lookup(h) is not None


def test_refetch_now_lands_all_hashes_or_fails():
    svc = FakeService()
    eng = OnDemandEngine(svc)
    legs = _legs()
    h = legset.leg_set_hash(legs)
    ok = eng.refetch_now([(h, GAME, legs)], deadline_sec=5.0)
    assert ok is True
    assert eng.lookup(h) is not None
    # deadline overrun -> False (slow service)
    slow = FakeService(latency=0.3)
    eng2 = OnDemandEngine(slow)
    ok2 = eng2.refetch_now([(h, GAME, legs)], deadline_sec=0.05)
    assert ok2 is False


def test_refetch_now_awaits_inflight_rather_than_duplicating():
    svc = FakeService(latency=0.2)
    eng = OnDemandEngine(svc)
    legs = _legs()
    h = legset.leg_set_hash(legs)
    eng.ensure_fetch(h, GAME, legs)          # feed fetch in flight (0.2s/book)
    time.sleep(0.05)
    ok = eng.refetch_now([(h, GAME, legs)], deadline_sec=5.0)
    assert ok is True
    # 2 books × 1 fetch — the confirm lane did NOT fire a duplicate fetch
    assert len(svc.calls) == 2


def test_no_poll_no_fetch_invariant():
    """An RFQ absent from the poll generates zero further fetches: the engine
    only ever fetches via ensure_fetch/refetch_now — nothing self-schedules."""
    clock = FakeClock()
    svc = FakeService()
    eng = OnDemandEngine(svc, now_fn=clock, autostart=False)
    legs = _legs()
    h = legset.leg_set_hash(legs)
    eng.ensure_fetch(h, GAME, legs)
    eng._drain_once()
    n = len(svc.calls)
    clock.t += 500                            # way past freshness
    eng._drain_once()                         # nothing queued
    assert len(svc.calls) == n and eng._queue_len() == 0


def test_refetch_now_never_reuses_a_pre_entry_result():
    """Blocker fix (adversarial review #1): a feed result that is still
    FRESH when the accept arrives must NOT satisfy the confirm last look —
    the fill gate re-fetches unless it can await an in-flight landing."""
    svc = FakeService()
    eng = OnDemandEngine(svc)
    legs = _legs()
    h = legset.leg_set_hash(legs)
    eng.ensure_fetch(h, GAME, legs)
    deadline = time.monotonic() + 5.0
    while eng.lookup(h) is None and time.monotonic() < deadline:
        time.sleep(0.02)
    assert eng.lookup(h) is not None          # fresh feed result in store
    calls_before = len(svc.calls)
    ok = eng.refetch_now([(h, GAME, legs)], deadline_sec=5.0)
    assert ok is True
    assert len(svc.calls) == calls_before + len(svc.books)   # re-fetched live


def test_ensure_fetch_reaps_wedged_inflight_and_restarts_worker():
    clock = FakeClock()
    svc = FakeService()
    eng = OnDemandEngine(svc, now_fn=clock, autostart=False)
    legs = _legs()
    h = legset.leg_set_hash(legs)
    eng.ensure_fetch(h, GAME, legs)
    # simulate a worker death mid-job: job pulled off the queue, never finished
    with eng._lock:
        eng._queue.clear()
    assert eng._inflight_len() == 1
    # within the reap window the dedup holds
    assert eng.ensure_fetch(h, GAME, legs) is False
    # after the reap bound, the wedged flight is cleared and we re-enqueue
    clock.t += 301
    assert eng.ensure_fetch(h, GAME, legs) is True
    assert eng._queue_len() == 1
