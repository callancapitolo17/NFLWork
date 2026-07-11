"""Phase 2 on-demand pricing engine: queue / worker / result store.

The open-RFQ poll is the feed driver (spec §4.5): the discovery tick calls
`ensure_fetch` whenever an on-demand sub-combo has no fresh result, and
`lookup` only ever serves results younger than QUOTE_FRESH_SEC — so quotes
are always backed by a fetch triggered by the RFQ being priced, fetching
stops the moment an RFQ leaves the poll, and nothing here self-schedules.

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
import threading
import time
from collections import deque
from concurrent.futures import ThreadPoolExecutor

log = logging.getLogger(__name__)

# A result may back a NEW quote only this long after landing (one 2s
# discovery tick + fetch-landing jitter). Correctness guard, not a knob.
QUOTE_FRESH_SEC = 15.0

# Store entries older than this are pruned (kept briefly past freshness
# only for research comparison / logging).
_STORE_RETENTION_SEC = 300.0


class OnDemandEngine:
    def __init__(self, service, now_fn=time.monotonic, autostart: bool = True):
        self._service = service
        self._now = now_fn
        self._autostart = autostart
        self._lock = threading.Lock()
        self._wake = threading.Event()
        self._queue: deque = deque()        # (hash, GameRef, legs)
        self._inflight: set[str] = set()    # queued or currently fetching
        self._landed: dict[str, threading.Event] = {}   # hash -> landing signal
        self._store: dict = {}              # hash -> (landed_at, {book: OnDemandBookResult})
        self._worker: threading.Thread | None = None

    # ------------------------------------------------------------- #
    # Feed side (discovery tick)                                     #
    # ------------------------------------------------------------- #

    def ensure_fetch(self, hash_: str, game, legs) -> bool:
        """Queue a fetch for this sub-combo unless one is already in flight.
        Idempotent per flight — the poll calls this every tick while pending.
        Returns True iff this call newly enqueued (lets the caller emit the
        on_demand_requested research event once per flight, not per tick)."""
        with self._lock:
            if hash_ in self._inflight:
                return False
            self._inflight.add(hash_)
            self._landed[hash_] = threading.Event()
            self._queue.append((hash_, game, legs))
        self._wake.set()
        if self._autostart:
            self._ensure_worker()
        return True

    def landed_at(self, hash_: str) -> float | None:
        """Monotonic landing time of the stored result (fresh or not), or
        None. Research-event dedup only — quoting freshness lives in lookup."""
        with self._lock:
            ent = self._store.get(hash_)
        return ent[0] if ent else None

    def lookup(self, hash_: str):
        """{book: fair} if landed within QUOTE_FRESH_SEC, else None."""
        results = self.lookup_results(hash_)
        if not results:
            return None
        return {b: r.fair for b, r in results.items()}

    def lookup_results(self, hash_: str):
        """{book: OnDemandBookResult} while fresh, else None (research read)."""
        with self._lock:
            ent = self._store.get(hash_)
        if ent is None:
            return None
        landed_at, results = ent
        if (self._now() - landed_at) > QUOTE_FRESH_SEC or not results:
            return None
        return dict(results)

    # ------------------------------------------------------------- #
    # Confirm-tick last look (priority lane, caller's thread)        #
    # ------------------------------------------------------------- #

    def refetch_now(self, jobs, deadline_sec: float) -> bool:
        """Synchronously land fresh results for EVERY (hash, game, legs) job
        within deadline_sec. Awaits an already-in-flight feed fetch for a
        hash instead of duplicating it. True iff all hashes are fresh at
        return; False on any overrun/failure (caller voids — fail-safe)."""
        deadline = self._now() + deadline_sec
        try:
            for hash_, game, legs in jobs:
                with self._lock:
                    inflight = hash_ in self._inflight
                    event = self._landed.get(hash_)
                if inflight and event is not None:
                    remaining = deadline - self._now()
                    if remaining <= 0 or not event.wait(timeout=remaining):
                        return False
                elif self.lookup(hash_) is None:
                    if self._now() >= deadline:
                        return False
                    self._fetch_combo(hash_, game, legs,
                                      deadline=deadline)
                if self.lookup(hash_) is None:
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
            while self._drain_once():
                pass

    def _queue_len(self) -> int:
        with self._lock:
            return len(self._queue)

    def _inflight_len(self) -> int:
        with self._lock:
            return len(self._inflight)

    def _drain_once(self) -> bool:
        """Process one queued job. Returns False when the queue is empty.
        Never raises (a raise would kill the consumer thread)."""
        with self._lock:
            if not self._queue:
                return False
            hash_, game, legs = self._queue.popleft()
        try:
            self._fetch_combo(hash_, game, legs)
        except Exception as e:
            log.warning("on_demand fetch error for %s: %s", hash_[:12], e)
            self._finish(hash_, {})
        return True

    def _fetch_combo(self, hash_: str, game, legs, deadline=None) -> None:
        """Fan out to every book concurrently; land whatever survives."""
        books = tuple(self._service.books)
        results = {}
        pool = ThreadPoolExecutor(max_workers=max(1, len(books)),
                                  thread_name_prefix="on-demand-book")
        try:
            futs = {b: pool.submit(self._price_book_safe, b, game, legs)
                    for b in books}
            for b, fut in futs.items():
                timeout = None
                if deadline is not None:
                    timeout = max(0.0, deadline - self._now())
                try:
                    r = fut.result(timeout=timeout)
                except Exception:
                    r = None
                if r is not None:
                    results[b] = r
        finally:
            pool.shutdown(wait=False, cancel_futures=True)
        self._finish(hash_, results)

    def _price_book_safe(self, book, game, legs):
        try:
            return self._service.price_on_demand(book, game, legs)
        except Exception as e:
            log.warning("on_demand %s error: %s", book, e)
            return None

    def _finish(self, hash_: str, results: dict) -> None:
        now = self._now()
        with self._lock:
            self._store[hash_] = (now, results)
            self._inflight.discard(hash_)
            event = self._landed.pop(hash_, None)
            # prune stale research leftovers
            for k in [k for k, (t, _) in self._store.items()
                      if now - t > _STORE_RETENTION_SEC]:
                del self._store[k]
        if event is not None:
            event.set()
