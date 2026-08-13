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


# --------------------------------------------------------------------- #
# Issue #50 item 3: per-call budgets on the live path                     #
# --------------------------------------------------------------------- #

class UnevenService(FakeService):
    """Per-book latency map — models one hung book among fast ones."""

    def __init__(self, latency_by_book, **kwargs):
        super().__init__(**kwargs)
        self.latency_by_book = dict(latency_by_book)

    def price_on_demand(self, book, game, legs):
        with self.lock:
            self.calls.append((book, legset.leg_set_hash(legs)))
        time.sleep(self.latency_by_book.get(book, 0.0))
        if book in self.fail_books:
            return None
        return OnDemandBookResult(book=book, fair=self.fair, route="partition",
                                  n_cells_priced=8, latency_sec=0.1)


def _await_landing(eng, h, timeout=10.0):
    deadline = time.monotonic() + timeout
    while eng.landed_at(h) is None and time.monotonic() < deadline:
        time.sleep(0.02)
    return eng.landed_at(h) is not None


def test_hung_book_is_ignored_and_fast_books_land_within_budget():
    """One book over the on-demand budget must not stall the combo to the
    75s feed deadline — the fast books' results land, the straggler is
    dropped. This is the ticket's headline engine behavior."""
    svc = UnevenService({"draftkings": 0.0, "fanduel": 2.0})
    svc.on_demand_deadline_sec = 0.3
    eng = OnDemandEngine(svc)
    legs = _legs()
    h = legset.leg_set_hash(legs)
    t0 = time.monotonic()
    eng.ensure_fetch(h, GAME, legs)
    assert _await_landing(eng, h)
    wall = time.monotonic() - t0
    assert wall < 1.5, f"combo waited {wall:.2f}s for a hung book"
    assert eng.lookup(h) == {"draftkings": svc.fair}


def test_engine_prefers_on_demand_deadline_over_sweep_deadline():
    """The 75s sweep deadline (per_book_deadline_sec) is sized for the
    background path; the live path must read the service's dedicated
    on_demand_deadline_sec instead."""
    svc = UnevenService({"draftkings": 0.0, "fanduel": 1.5})
    svc.per_book_deadline_sec = 75.0
    svc.on_demand_deadline_sec = 0.2
    eng = OnDemandEngine(svc)
    legs = _legs()
    h = legset.leg_set_hash(legs)
    t0 = time.monotonic()
    eng.ensure_fetch(h, GAME, legs)
    assert _await_landing(eng, h)
    assert time.monotonic() - t0 < 1.0
    assert eng.lookup(h) == {"draftkings": svc.fair}


# --------------------------------------------------------------------- #
# Issue #50 item 3 (second finding): multi-game jobs drain concurrently   #
# --------------------------------------------------------------------- #

def test_multi_game_jobs_drain_concurrently_not_serially():
    """A cross-game RFQ enqueues one job per game. The old single-consumer
    worker drained them strictly serially — game 2 waited out ALL of game
    1, including its hung straggler. Now jobs overlap: game 1's hung book
    must not delay game 2's healthy books (per-book calls pipeline behind
    the #40 gate; other books run freely)."""
    svc = UnevenService({"draftkings": 0.0, "fanduel": 3.0})
    svc.on_demand_deadline_sec = 0.5
    eng = OnDemandEngine(svc)
    legs_a = _legs()
    legs_b = _legs()[:2]
    h_a = legset.leg_set_hash(legs_a)
    h_b = legset.leg_set_hash(legs_b)
    assert h_a != h_b
    t0 = time.monotonic()
    eng.ensure_fetch(h_a, GAME, legs_a)
    eng.ensure_fetch(h_b, GAME, legs_b)
    assert _await_landing(eng, h_a) and _await_landing(eng, h_b)
    wall = time.monotonic() - t0
    # serial drain: job 2 starts only after job 1 fully lands -> >= 2x the
    # hung book (or 2x deadline). Concurrent: both land ~one deadline.
    assert wall < 2.0, f"jobs drained serially: {wall:.2f}s"
    assert eng.lookup(h_a) == {"draftkings": svc.fair}
    assert eng.lookup(h_b) == {"draftkings": svc.fair}


def test_per_book_calls_never_overlap_across_jobs():
    """#40 pacing invariant survives concurrent jobs: at most ONE
    on-demand pricing call in flight per book, ever (Novig 403s under
    bursts — the gate is what makes concurrent jobs rate-limit-safe)."""
    class TrackingService(FakeService):
        def __init__(self, **kw):
            super().__init__(**kw)
            self.active = {}
            self.max_active = {}

        def price_on_demand(self, book, game, legs):
            with self.lock:
                self.active[book] = self.active.get(book, 0) + 1
                self.max_active[book] = max(self.max_active.get(book, 0),
                                            self.active[book])
            time.sleep(0.15)
            with self.lock:
                self.active[book] -= 1
            return OnDemandBookResult(book=book, fair=self.fair,
                                      route="partition", n_cells_priced=8,
                                      latency_sec=0.1)

    svc = TrackingService()
    svc.on_demand_deadline_sec = 5.0
    eng = OnDemandEngine(svc)
    hashes = []
    for i in range(3):
        legs = _legs()[:2] + [legset.CanonicalLeg(EVT, "total",
                                                  7.5 + i, "over")]
        h = legset.leg_set_hash(legs)
        hashes.append(h)
        eng.ensure_fetch(h, GAME, legs)
    for h in hashes:
        assert _await_landing(eng, h)
        assert eng.lookup(h) is not None
    assert svc.max_active, "no calls recorded"
    assert all(v == 1 for v in svc.max_active.values()), svc.max_active


def test_job_flood_past_the_concurrency_cap_fully_drains():
    """More jobs than _MAX_CONCURRENT_JOBS: the semaphore backpressures the
    dispatcher, every slot is released, and every job still lands."""
    svc = FakeService(latency=0.2)
    svc.on_demand_deadline_sec = 5.0
    eng = OnDemandEngine(svc)
    hashes = []
    for i in range(6):
        legs = _legs()[:2] + [legset.CanonicalLeg(EVT, "total",
                                                  7.5 + i, "over")]
        h = legset.leg_set_hash(legs)
        hashes.append(h)
        eng.ensure_fetch(h, GAME, legs)
    assert len(set(hashes)) == 6
    for h in hashes:
        assert _await_landing(eng, h), "a job never landed — leaked slot?"
        assert eng.lookup(h) is not None


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


# --------------------------------------------------------------------- #
# Issue #86: F5-winner tie conversion at the lookup choke point          #
# --------------------------------------------------------------------- #

class F5FakeService:
    """Duck-typed service whose results carry #86 semantics/tie fields."""
    def __init__(self, book_results):
        self.book_results = dict(book_results)   # book -> OnDemandBookResult
        self.books = tuple(self.book_results)

    def price_on_demand(self, book, game, legs):
        return self.book_results.get(book)


def _f5_result(book, fair, p_tie_book=None, semantics="push_2way"):
    return OnDemandBookResult(book=book, fair=fair, route="partition",
                              n_cells_priced=8, latency_sec=0.1,
                              f5_semantics=semantics, p_tie_book=p_tie_book)


def _f5_legs():
    return [legset.CanonicalLeg(EVT, "ml", None, "home", "F5"),
            legset.CanonicalLeg(EVT, "total", 4.5, "over", "F5")]


def _land(svc):
    eng = OnDemandEngine(svc, now_fn=FakeClock(), autostart=False)
    legs = _f5_legs()
    h = legset.leg_set_hash(legs)
    eng.ensure_fetch(h, GAME, legs)
    eng._drain_once()
    return eng, h


def test_lookup_converts_f5_winner_fairs_with_median_anchor():
    """Three conditional fairs, two anchors (FD 0.15, MGM 0.13): the flight
    median 0.14 converts EVERY book — the anchorless one included — and the
    converted (not raw) fairs are what consensus sees."""
    svc = F5FakeService({
        "fanduel": _f5_result("fanduel", 0.62, p_tie_book=0.15),
        "betmgm": _f5_result("betmgm", 0.64, p_tie_book=0.13),
        "novig": _f5_result("novig", 0.63),          # no 3-way at NV
    })
    eng, h = _land(svc)
    fairs = eng.lookup(h)
    assert fairs is not None
    assert fairs["fanduel"] == pytest.approx(0.62 * (1 - 0.14))
    assert fairs["betmgm"] == pytest.approx(0.64 * (1 - 0.14))
    assert fairs["novig"] == pytest.approx(0.63 * (1 - 0.14))


def test_lookup_returns_empty_when_no_tie_anchor_landed():
    """F5-winner flight with zero anchors: nothing is convertible, lookup
    serves {} (consensus declines too_few_books) — NEVER raw conditional
    fairs, and NEVER None (which would read as 'nothing landed' and mask
    the flight from research emission)."""
    svc = F5FakeService({
        "novig": _f5_result("novig", 0.63),
        "draftkings": _f5_result("draftkings", 0.62),
    })
    eng, h = _land(svc)
    assert eng.lookup(h) == {}
    # research reads see the RAW results (provenance stays intact)
    raw = eng.lookup_results(h)
    assert raw is not None and raw["novig"].fair == 0.63


def test_lookup_passes_non_f5_flights_through_unconverted():
    svc = F5FakeService({
        "fanduel": OnDemandBookResult(book="fanduel", fair=0.14,
                                      route="partition", n_cells_priced=8,
                                      latency_sec=0.1),
    })
    eng = OnDemandEngine(svc, now_fn=FakeClock(), autostart=False)
    legs = _legs()
    h = legset.leg_set_hash(legs)
    eng.ensure_fetch(h, GAME, legs)
    eng._drain_once()
    assert eng.lookup(h) == {"fanduel": 0.14}


def test_lookup_excludes_unconvertible_semantics_defensively():
    """A result labeled with an unknown semantics string converts to None
    and drops out; the labeled books still serve."""
    svc = F5FakeService({
        "fanduel": _f5_result("fanduel", 0.62, p_tie_book=0.14),
        "novig": _f5_result("novig", 0.63, semantics="mystery_mode"),
    })
    eng, h = _land(svc)
    fairs = eng.lookup(h)
    assert set(fairs) == {"fanduel"}
    assert fairs["fanduel"] == pytest.approx(0.62 * (1 - 0.14))


def test_mixed_semantics_agree_after_conversion_not_before():
    """Issue #86 spec case: a push-2way book (conditional 0.62) and a
    hypothetical 3-way book (unconditional 0.5332) DISAGREE by ~9pts raw —
    raw sigma_z would blow the gate on two views that are actually
    identical. Converted, both land at 0.5332 and consensus passes."""
    from kalshi_mlb_mm import router

    svc = F5FakeService({
        "fanduel": _f5_result("fanduel", 0.62, p_tie_book=0.14),
        "betmgm": _f5_result("betmgm", 0.62 * (1 - 0.14),
                             semantics="three_way_unconditional"),
    })
    eng, h = _land(svc)
    fairs = eng.lookup(h)
    assert fairs["fanduel"] == pytest.approx(fairs["betmgm"])
    cons, reason = router.consensus(fairs, min_books=2, sigma_z_max=0.07)
    assert reason == "ok"
    assert cons.fair == pytest.approx(0.62 * (1 - 0.14))
    # the RAW pair would have been declined for dispersion — the mix is
    # only quotable BECAUSE conversion ran before the gate
    raw = {b: r.fair for b, r in eng.lookup_results(h).items()}
    _, raw_reason = router.consensus(raw, min_books=2, sigma_z_max=0.07)
    assert raw_reason == "consensus_dispersion"


def test_tie_anchor_accessor_reports_median_and_none():
    svc = F5FakeService({
        "fanduel": _f5_result("fanduel", 0.62, p_tie_book=0.15),
        "betmgm": _f5_result("betmgm", 0.64, p_tie_book=0.13),
    })
    eng, h = _land(svc)
    assert eng.tie_anchor(h) == pytest.approx(0.14)
    svc2 = F5FakeService({"novig": _f5_result("novig", 0.63)})
    eng2, h2 = _land(svc2)
    assert eng2.tie_anchor(h2) is None
