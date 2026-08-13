"""Phase 2 on-demand pricing engine: queue / worker / result store.

The open-RFQ poll is the feed driver (spec §4.5): the discovery tick calls
`ensure_fetch` whenever an on-demand sub-combo has no fresh result, and
`lookup` only ever serves results younger than QUOTE_FRESH_SEC — so quotes
are always backed by a fetch triggered by the RFQ being priced, fetching
stops the moment an RFQ leaves the poll, and nothing here self-schedules.
Quorum quoting (issue #55): landings are INCREMENTAL — each book's result
is servable the moment it returns, so the poll quotes on the fast books
while stragglers are still in flight (see _Flight for the semantics).

Threading model: one daemon consumer thread drains a FIFO of combo jobs;
each job fans out to every book concurrently (thread per book, mirroring
SGPService.refresh). `refetch_now` (confirm-tick last look) runs on the
CALLER's thread so it can honour the confirm window deadline — it dedups
against an in-flight feed fetch by awaiting its landing instead of
duplicating it. Fail-safe: worker exceptions are caught; a dead worker is
restarted lazily on the next enqueue; all state is in-memory and ephemeral.
"""
from __future__ import annotations

import logging
import statistics
import threading
import time
from collections import deque
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor
from concurrent.futures import wait as futures_wait
from dataclasses import dataclass, field

log = logging.getLogger(__name__)

# A result may back a NEW quote only this long after landing (one 2s
# discovery tick + fetch-landing jitter). Correctness guard, not a knob.
QUOTE_FRESH_SEC = 15.0

# Store entries older than this are pruned (kept briefly past freshness
# only for research comparison / logging).
_STORE_RETENTION_SEC = 300.0

# In-flight entries older than this are reaped (a worker that died mid-job
# can no longer produce a usable result for them; without reaping, the
# dedup would suppress re-fetches of that combo forever).
_INFLIGHT_REAP_SEC = 300.0

# Fallback wall bound for FEED fetches when the service carries no
# on_demand_deadline_sec (duck-typed test services). The real budget comes
# from SGPService.on_demand_deadline_sec (issue #50): a book still running
# past it is dropped and the fast books' results land. The confirm lane
# passes its own confirm-window deadline.
_FEED_FETCH_DEADLINE_SEC = 75.0

@dataclass
class _Flight:
    """One leg-set's stored fetch state (issue #55: incremental landings).

    landed_at   monotonic time of the MOST RECENT landing (each book's
                arrival slides it forward; completion stamps it once more).
                Freshness is measured from here, so the worst-case serve
                window is unchanged from the pre-quorum atomic landing:
                job deadline + QUOTE_FRESH_SEC after fetch start.
    complete    False while the job is still fanning out. landed_empty
                (-> live_fetch_timeout) fires ONLY on complete+empty —
                a zero-book PARTIAL stays "pending", exactly as before.
    dropped_books  books whose call was still running at the job deadline
                (#50 budget). Failed books (returned None in time) are NOT
                dropped — the distinction feeds #38's health picture.
    """
    landed_at: float
    results: dict = field(default_factory=dict)
    complete: bool = False
    dropped_books: tuple = ()


# Issue #50: multi-game RFQs enqueue one job per game; the old single
# consumer drained them strictly serially, so game N's fetch waited out
# games 1..N-1 in full. Jobs run on up to ON_DEMAND_MAX_CONCURRENT_JOBS
# concurrent daemon threads (config knob, 2026-08-11: 4 slots capped combo
# throughput at ~50/min while option-B door survivors arrive at ~170/min —
# a job holds its slot for the full flight, including the slow books).
# Raising jobs does NOT raise per-book pressure: the per-book gates below
# still serialize to ONE pricing call in flight per book — more jobs just
# keep every book's serial lane full instead of idle.
_MAX_CONCURRENT_JOBS = 4   # fallback when config is absent (tests)


class OnDemandEngine:
    def __init__(self, service, now_fn=time.monotonic, autostart: bool = True,
                 max_concurrent_jobs: int | None = None):
        self._service = service
        self._now = now_fn
        self._autostart = autostart
        self._lock = threading.Lock()
        self._wake = threading.Event()
        self._queue: deque = deque()        # (hash, GameRef, legs)
        self._inflight: dict[str, float] = {}   # hash -> enqueue time
        self._landed: dict[str, threading.Event] = {}   # hash -> landing signal
        self._store: dict[str, _Flight] = {}            # hash -> _Flight
        self._worker: threading.Thread | None = None
        if max_concurrent_jobs is None:
            try:
                from kalshi_mlb_mm import config as _config
                max_concurrent_jobs = _config.ON_DEMAND_MAX_CONCURRENT_JOBS
            except (ImportError, AttributeError):
                max_concurrent_jobs = _MAX_CONCURRENT_JOBS
        self.max_concurrent_jobs = max_concurrent_jobs
        # Caps concurrent jobs; acquired by the dispatcher, released by the
        # job thread — the dispatcher blocks (backpressure) at the cap.
        self._job_slots = threading.Semaphore(max_concurrent_jobs)
        # Issue #40 pacing invariant: ONE on-demand pricing call in flight
        # per book, even with concurrent jobs (Novig 403s at ~26 rapid
        # calls; a 3-job burst of 8-cell partitions would blow past it).
        # Jobs therefore pipeline per book while their other books run
        # freely. Lazily created; guarded by self._lock.
        self._book_gates: dict[str, threading.Semaphore] = {}

    # ------------------------------------------------------------- #
    # Feed side (discovery tick)                                     #
    # ------------------------------------------------------------- #

    def ensure_fetch(self, hash_: str, game, legs) -> bool:
        """Queue a fetch for this sub-combo unless one is already in flight.
        Idempotent per flight — the poll calls this every tick while pending.
        Returns True iff this call newly enqueued (lets the caller emit the
        on_demand_requested research event once per flight, not per tick)."""
        if self._autostart:
            # Unconditionally, BEFORE the in-flight dedup: a worker that died
            # mid-job leaves its hash in-flight, and if all subsequent poll
            # traffic is for that hash the early-return would otherwise never
            # reach a restart.
            self._ensure_worker()
        now = self._now()
        with self._lock:
            # Reap wedged in-flight entries (worker died mid-job): a flight
            # older than the reap bound can't still be producing a usable
            # (<=15s-fresh) result, so clearing it is always safe.
            stale = [h for h, t in self._inflight.items()
                     if now - t > _INFLIGHT_REAP_SEC]
            for h in stale:
                del self._inflight[h]
                ev = self._landed.pop(h, None)
                if ev is not None:
                    ev.set()
                # Drop the dead flight's PARTIAL landing too: a later flight
                # for this hash must never merge with (or resurrect) book
                # fairs from a worker that died mid-job.
                ent = self._store.get(h)
                if ent is not None and not ent.complete:
                    del self._store[h]
            if hash_ in self._inflight:
                return False
            self._inflight[hash_] = now
            self._landed[hash_] = threading.Event()
            self._queue.append((hash_, game, legs))
        self._wake.set()
        return True

    def landed_at(self, hash_: str) -> float | None:
        """Monotonic time of the most recent landing (fresh or not), or
        None. Research-event dedup only — quoting freshness lives in lookup."""
        with self._lock:
            ent = self._store.get(hash_)
        return ent.landed_at if ent else None

    def landed_empty(self, hash_: str) -> bool:
        """True iff the last fetch COMPLETED within QUOTE_FRESH_SEC with ZERO
        book results (every book over budget or declining). The live quote
        path (#54) declines on this — never a cache fallback — and while the
        empty landing is fresh the poll does not re-feed, so a dead slate
        retries on the ~QUOTE_FRESH_SEC cadence instead of every tick.
        Quorum (#55): a zero-book PARTIAL is not a landing at all — the poll
        keeps the RFQ pending until the flight completes."""
        with self._lock:
            ent = self._store.get(hash_)
        if ent is None or not ent.complete:
            return False
        return (not ent.results
                and (self._now() - ent.landed_at) <= QUOTE_FRESH_SEC)

    def result_age_sec(self, hash_: str) -> float | None:
        """Seconds since the most recent landing (fresh or not), or None.
        Research payloads only (#54 live-trace proof): quoting freshness
        stays in lookup."""
        with self._lock:
            ent = self._store.get(hash_)
        if ent is None:
            return None
        return self._now() - ent.landed_at

    def completed_fetch_age_sec(self, hash_: str) -> float | None:
        """Seconds since the last COMPLETED landing that produced >=1 book
        result; None for never-fetched, still-in-flight, or zero-book
        completions. The #57 post-fill cooldown gate reads this: a complete,
        priced flight is proof the books were actually re-asked (and at
        least one answered) — an in-flight partial or an all-declined
        landing is not."""
        with self._lock:
            ent = self._store.get(hash_)
        if ent is None or not ent.complete or not ent.results:
            return None
        return self._now() - ent.landed_at

    def dropped_books(self, hash_: str) -> tuple:
        """Books dropped at the job deadline on the last COMPLETED flight
        (#55 item 4 observability; a chronically-dropped book is a de facto
        non-participant). Empty while in flight or when nothing is stored."""
        with self._lock:
            ent = self._store.get(hash_)
        return ent.dropped_books if ent else ()

    def lookup(self, hash_: str):
        """{book: fair} if landed within QUOTE_FRESH_SEC, else None.

        Issue #86: this is the single choke point every consumer reads
        (quote tick, confirm last-look, book-drift breaker), so F5-winner
        conversion lives HERE. A flight whose results carry an
        ``f5_semantics`` label priced a combo containing an F5-winner leg —
        its push-2way book fairs are CONDITIONAL on no-tie and convert via
        fair_value.f5_winner_fair_from_book with the flight's median tie
        anchor. Zero anchors landed -> {} (nothing convertible; consensus
        declines too_few_books) — NEVER the raw conditional fairs, and
        never None (which would read as "nothing landed" and hide the
        flight from research emission). ``lookup_results`` keeps serving
        the raw results for research provenance."""
        results = self.lookup_results(hash_)
        if not results:
            return None
        if not any(getattr(r, "f5_semantics", None) for r in results.values()):
            return {b: r.fair for b, r in results.items()}
        from kalshi_common import fair_value
        p_tie = self._median_tie_anchor(results)
        out = {}
        for b, r in results.items():
            fair = fair_value.f5_winner_fair_from_book(
                r.fair, p_tie, getattr(r, "f5_semantics", None))
            if fair is not None:
                out[b] = fair
        return out

    @staticmethod
    def _median_tie_anchor(results) -> float | None:
        """Median of the landed books' own 3-way tie anchors (FD/MGM), or
        None. One shared anchor per flight keeps the sigma_z dispersion
        gate measuring conditional-fair disagreement, not tie-estimate
        noise."""
        anchors = [r.p_tie_book for r in results.values()
                   if getattr(r, "p_tie_book", None) is not None]
        if not anchors:
            return None
        return statistics.median(anchors)

    def tie_anchor(self, hash_: str) -> float | None:
        """The flight's median P(tie) anchor while fresh (research trace:
        quote_priced carries it as per-game provenance), else None."""
        results = self.lookup_results(hash_)
        if not results:
            return None
        return self._median_tie_anchor(results)

    def lookup_results(self, hash_: str):
        """{book: OnDemandBookResult} while fresh, else None (research read).
        Quorum (#55): PARTIAL landings serve too — the moment one book has
        returned, the tick can see it; the #20 consensus gate (not the
        engine) decides whether the set is thick enough to quote."""
        with self._lock:
            ent = self._store.get(hash_)
        if ent is None:
            return None
        if (self._now() - ent.landed_at) > QUOTE_FRESH_SEC or not ent.results:
            return None
        return dict(ent.results)

    # ------------------------------------------------------------- #
    # Confirm-tick last look (priority lane, caller's thread)        #
    # ------------------------------------------------------------- #

    def refetch_now(self, jobs, deadline_sec: float) -> bool:
        """Synchronously land results for EVERY (hash, game, legs) job that
        arrive AFTER this call was entered, within deadline_sec. A result
        already in the store — however fresh — never satisfies the last
        look ("we never confirm on a previously-fetched number"); the only
        non-duplicating shortcut is awaiting an already-in-flight feed fetch,
        whose landing is also after entry. True iff every hash lands in
        time; False on any overrun/failure (caller voids — fail-safe)."""
        entry = self._now()
        deadline = entry + deadline_sec
        try:
            for hash_, game, legs in jobs:
                with self._lock:
                    inflight = hash_ in self._inflight
                    event = self._landed.get(hash_)
                if inflight and event is not None:
                    remaining = deadline - self._now()
                    if remaining <= 0 or not event.wait(timeout=remaining):
                        return False
                else:
                    if self._now() >= deadline:
                        return False
                    self._fetch_combo(hash_, game, legs,
                                      deadline=deadline)
                landed = self.landed_at(hash_)
                if landed is None or landed < entry or self.lookup(hash_) is None:
                    return False
            return True
        except Exception as e:                      # never raise into confirm
            log.warning("on_demand refetch_now error: %s", e)
            return False

    # ------------------------------------------------------------- #
    # Internals                                                      #
    # ------------------------------------------------------------- #

    def _ensure_worker(self) -> None:
        with self._lock:
            if self._worker is not None and self._worker.is_alive():
                return
            self._worker = threading.Thread(target=self._run, daemon=True,
                                            name="on-demand-worker")
            self._worker.start()

    def _run(self) -> None:
        while True:
            self._wake.wait(timeout=1.0)
            self._wake.clear()
            while self._dispatch_once():
                pass

    def _queue_len(self) -> int:
        with self._lock:
            return len(self._queue)

    def _inflight_len(self) -> int:
        with self._lock:
            return len(self._inflight)

    def _feed_deadline_sec(self) -> float:
        """The live path's per-job wall budget (issue #50): the service's
        on_demand_deadline_sec — NEVER per_book_deadline_sec, which is the
        75s background-sweep budget. A quote is worthless at sweep speed."""
        budget = getattr(self._service, "on_demand_deadline_sec", None)
        try:
            return float(budget)
        except (TypeError, ValueError):
            return _FEED_FETCH_DEADLINE_SEC

    def _dispatch_once(self) -> bool:
        """Move one queued job onto its own daemon thread. Returns False
        when the queue is empty. Blocks at _MAX_CONCURRENT_JOBS in flight
        (backpressure on the dispatcher, never a dropped job)."""
        with self._lock:
            if not self._queue:
                return False
            hash_, game, legs = self._queue.popleft()
        self._job_slots.acquire()
        try:
            threading.Thread(target=self._process_job,
                             args=(hash_, game, legs),
                             daemon=True, name="on-demand-job").start()
        except Exception as e:
            # Thread creation failing (fd/thread exhaustion) must not leak
            # the slot — four leaks would wedge the dispatcher forever.
            self._job_slots.release()
            log.warning("on_demand job spawn failed for %s: %s", hash_[:12], e)
            self._finish(hash_, {})
        return True

    def _process_job(self, hash_: str, game, legs) -> None:
        """One job, on its own thread. Never raises; always releases its
        slot and always lands SOMETHING for the hash (a failed job lands
        {} so the poll's dedup can re-feed instead of wedging)."""
        try:
            self._fetch_combo(hash_, game, legs,
                              deadline=self._now() + self._feed_deadline_sec())
        except Exception as e:
            log.warning("on_demand fetch error for %s: %s", hash_[:12], e)
            self._finish(hash_, {})
        finally:
            self._job_slots.release()

    def _drain_once(self) -> bool:
        """Synchronous single-job drain — the autostart=False test seam
        (production traffic flows through _dispatch_once). Never raises."""
        with self._lock:
            if not self._queue:
                return False
            hash_, game, legs = self._queue.popleft()
        try:
            self._fetch_combo(hash_, game, legs,
                              deadline=self._now() + self._feed_deadline_sec())
        except Exception as e:
            log.warning("on_demand fetch error for %s: %s", hash_[:12], e)
            self._finish(hash_, {})
        return True

    def _fetch_combo(self, hash_: str, game, legs, deadline=None) -> None:
        """Fan out to every book concurrently; land each result AS IT
        RETURNS (quorum quoting, #55) so the poll can quote on the fast
        books while stragglers are still in flight. Books still running at
        the deadline are dropped from this flight and recorded as such.

        Partial landings are feed-path only (_land_partial no-ops unless
        the hash is in-flight); correctness of the final landing rests on
        the local `results` dict, so the confirm lane's DIRECT call (hash
        not in-flight) still lands its full result set at _finish."""
        books = tuple(self._service.books)
        results = {}
        dropped: list = []
        pool = ThreadPoolExecutor(max_workers=max(1, len(books)),
                                  thread_name_prefix="on-demand-book")
        try:
            futs = {pool.submit(self._price_book_safe, b, game, legs,
                                deadline): b
                    for b in books}
            pending = set(futs)
            while pending:
                timeout = None
                if deadline is not None:
                    timeout = max(0.0, deadline - self._now())
                done, pending = futures_wait(pending, timeout=timeout,
                                             return_when=FIRST_COMPLETED)
                if not done:            # deadline hit — the rest are dropped
                    dropped = sorted(futs[f] for f in pending)
                    break
                for fut in done:
                    book = futs[fut]
                    try:
                        r = fut.result()
                    except Exception:
                        r = None
                    if r is not None:
                        results[book] = r
                        self._land_partial(hash_, book, r)
        finally:
            pool.shutdown(wait=False, cancel_futures=True)
        self._finish(hash_, results, dropped_books=dropped)

    def _book_gate(self, book) -> threading.Semaphore:
        with self._lock:
            gate = self._book_gates.get(book)
            if gate is None:
                gate = self._book_gates[book] = threading.Semaphore(1)
            return gate

    def _price_book_safe(self, book, game, legs, deadline=None):
        """One book's pricing call, gated to one-in-flight per book (#40).

        The gate wait is bounded by the job's deadline: a thread that
        cannot acquire in time returns None instead of firing a wire call
        whose result nobody will collect. The poll re-feeds the combo, so
        a skipped book costs one tick, not a rate-limit strike."""
        gate = self._book_gate(book)
        timeout = None
        if deadline is not None:
            timeout = max(0.0, deadline - self._now())
        if not gate.acquire(timeout=timeout):
            return None
        try:
            return self._service.price_on_demand(book, game, legs)
        except Exception as e:
            log.warning("on_demand %s error: %s", book, e)
            return None
        finally:
            gate.release()

    def _land_partial(self, hash_: str, book: str, result) -> None:
        """Land ONE book's result mid-flight (#55): lookup serves the
        partial set immediately. No-op once the flight is no longer
        in-flight — a book finishing past the deadline (or during a
        confirm-lane direct fetch) must never mutate a landed entry, which
        keeps refetch_now's landed-after-entry guarantee airtight."""
        now = self._now()
        with self._lock:
            if hash_ not in self._inflight:
                return
            prev = self._store.get(hash_)
            # A COMPLETE previous entry belongs to an older flight — start
            # fresh rather than mixing two flights' fairs.
            results = (dict(prev.results)
                       if prev is not None and not prev.complete else {})
            results[book] = result
            self._store[hash_] = _Flight(landed_at=now, results=results)

    def _finish(self, hash_: str, results: dict | None = None, *,
                dropped_books=()) -> None:
        """Finalize a flight: merge over any partials already landed (the
        failure paths call with no results and must never clobber them),
        mark it complete, release the in-flight dedup, wake refetch waiters."""
        now = self._now()
        with self._lock:
            prev = self._store.get(hash_)
            merged = (dict(prev.results)
                      if prev is not None and not prev.complete else {})
            if results:
                merged.update(results)
            self._store[hash_] = _Flight(landed_at=now, results=merged,
                                         complete=True,
                                         dropped_books=tuple(dropped_books))
            self._inflight.pop(hash_, None)
            event = self._landed.pop(hash_, None)
            # prune stale research leftovers
            for k in [k for k, f in self._store.items()
                      if now - f.landed_at > _STORE_RETENTION_SEC]:
                del self._store[k]
        if event is not None:
            event.set()
