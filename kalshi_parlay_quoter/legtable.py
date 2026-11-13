"""In-memory table of leg fairs, so pricing an RFQ is a dict lookup.

Why (2026-11-13): ~475 in-scope RFQs a second is ~1,900 leg prices a second,
far more than one runner call per RFQ could carry, and they re-ask the same
few thousand legs. The table holds every (ticker, side) seen in the last
`leg_ttl_sec`; Refresher asks the runner for all of them every couple of
seconds in chunks (pricer.fairs: board only, no Kalshi read). The work then
grows with distinct legs, not with RFQs.

A fair is served only while it is younger than `max_age_sec`, counted from
when its refresh call was SENT. Every served leg carries priced_at_ms (the
runner's stamp) and cache_age_ms, and fair_age_ms grows by the cache age, so
the log never shows a reused fair as fresh. Table legs are not Kalshi-checked
(kalshi_* fields are null): a would-quote off the table is confirmed through
pricer.price before it is logged as confirmed.
"""
import logging
import threading
import time
from typing import Callable

from kalshi_parlay_quoter.pricer import (FAIRS_CHUNK_LEGS, LegsNotReady, PriceResult,
                                         PricerError)

log = logging.getLogger("kalshi_parlay_quoter.legtable")

LegKey = tuple[str, str]          # (market_ticker, side)


class LegTable:
    def __init__(self, *, max_age_sec: float, leg_ttl_sec: float,
                 clock: Callable[[], float] = time.monotonic):
        self._max_age = max_age_sec
        self._leg_ttl = leg_ttl_sec
        self._clock = clock
        self._lock = threading.Lock()
        self._last_asked: dict[LegKey, float] = {}
        # key -> (leg as the runner answered, runner priced_at_ms, sent_at)
        self._fairs: dict[LegKey, tuple[dict, int | None, float]] = {}
        self.generation = 0        # bumps on every store(); waiters retry on change

    def __len__(self) -> int:
        return len(self._last_asked)

    def want(self, keys: list[LegKey]) -> None:
        """Keep these legs refreshed for leg_ttl_sec from now."""
        now = self._clock()
        with self._lock:
            for key in keys:
                self._last_asked[key] = now

    def wanted(self) -> list[LegKey]:
        """Legs asked for within leg_ttl_sec; older ones are forgotten."""
        cutoff = self._clock() - self._leg_ttl
        with self._lock:
            expired = [k for k, asked in self._last_asked.items() if asked < cutoff]
            for key in expired:
                del self._last_asked[key]
                self._fairs.pop(key, None)
            return list(self._last_asked)

    def store(self, keys: list[LegKey], legs: list[dict],
              priced_at_ms: int | None, sent_at: float) -> None:
        with self._lock:
            for key, leg in zip(keys, legs):
                if key in self._last_asked:
                    self._fairs[key] = (leg, priced_at_ms, sent_at)
            self.generation += 1

    def missing(self, keys: list[LegKey]) -> list[LegKey]:
        """The legs with no fair younger than max_age_sec."""
        cutoff = self._clock() - self._max_age
        with self._lock:
            return [k for k in keys
                    if (held := self._fairs.get(k)) is None or held[2] < cutoff]

    def price(self, rfq_id: str, keys: list[LegKey]) -> PriceResult:
        """decide.decide's price_fn. Raises LegsNotReady when any leg has no
        fresh fair. elapsed_ms is 0: nothing was fetched."""
        now = self._clock()
        cutoff = now - self._max_age
        out = []
        with self._lock:
            held = [self._fairs.get(k) for k in keys]
        stale = [k for k, h in zip(keys, held) if h is None or h[2] < cutoff]
        if stale:
            raise LegsNotReady(f"no fair under {self._max_age:.0f}s old for "
                               f"{len(stale)} of {len(keys)} legs, e.g. {stale[0][0]}")
        for leg, priced_at_ms, sent_at in held:
            age_ms = max(0.0, (now - sent_at) * 1000.0)
            aged = dict(leg)
            aged["priced_at_ms"] = priced_at_ms
            aged["cache_age_ms"] = age_ms
            if isinstance(leg.get("fair_age_ms"), (int, float)):
                aged["fair_age_ms"] = leg["fair_age_ms"] + age_ms
            out.append(aged)
        stamps = [leg["priced_at_ms"] for leg in out if leg["priced_at_ms"] is not None]
        return PriceResult(out, min(stamps) if stamps else None, 0.0)


class Refresher:
    """One thread: every interval_sec, re-price every wanted leg.

    fetch(keys) -> PriceResult for at most FAIRS_CHUNK_LEGS keys, may raise
    PricerError. A failed chunk ends the pass (the runner is down or busy;
    the next pass retries) and its legs simply age out of the table."""

    def __init__(self, table: LegTable, fetch: Callable[[list[LegKey]], PriceResult], *,
                 interval_sec: float, clock: Callable[[], float] = time.monotonic):
        self._table = table
        self._fetch = fetch
        self._interval = interval_sec
        self._clock = clock
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None
        self._lock = threading.Lock()
        self._failed_passes = 0
        self._last_pass_ms = 0.0

    def refresh_once(self) -> bool:
        """One pass over the wanted legs. True when every chunk landed."""
        keys = self._table.wanted()
        started = self._clock()
        ok = True
        for i in range(0, len(keys), FAIRS_CHUNK_LEGS):
            chunk = keys[i:i + FAIRS_CHUNK_LEGS]
            sent_at = self._clock()
            try:
                result = self._fetch(chunk)
            except PricerError as e:
                log.warning("[legtable] refresh failed at leg %d of %d: %s", i, len(keys), e)
                ok = False
                break
            self._table.store(chunk, result.legs, result.priced_at_ms, sent_at)
        with self._lock:
            self._last_pass_ms = (self._clock() - started) * 1000.0
            if not ok:
                self._failed_passes += 1
        return ok

    def take_stats(self) -> dict[str, int]:
        """Failed passes since the last call, and the last pass's duration."""
        with self._lock:
            failed, self._failed_passes = self._failed_passes, 0
            return {"table_refresh_failed": failed,
                    "gauge_table_refresh_ms": int(self._last_pass_ms)}

    def start(self) -> None:
        if self._thread is None:
            self._thread = threading.Thread(target=self._run, name="leg-table", daemon=True)
            self._thread.start()

    def stop(self) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=2.0)

    def _run(self) -> None:
        while not self._stop.is_set():
            started = self._clock()
            try:
                self.refresh_once()
            except Exception:
                log.exception("[legtable] refresh pass crashed")
            self._stop.wait(max(0.1, self._interval - (self._clock() - started)))
