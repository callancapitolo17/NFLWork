"""In-memory table of leg fairs, so pricing an RFQ is a dict lookup.

Why (2026-11-13): ~475 in-scope RFQs a second is ~1,900 leg prices a second,
far more than one runner call per RFQ could carry, and they re-ask the same
few thousand legs. The table holds every (ticker, side) seen in the last
`leg_ttl_sec`; Refresher asks the runner for all of them every couple of
seconds in chunks (pricer.fairs: board only, no Kalshi read). The work then
grows with distinct legs, not with RFQs.

A leg with no fair yet is priced first and between passes (Refresher), so
a new leg waits a fraction of a second, not a whole pass. A leg the runner
refuses outright is found by halving its chunk and declines
`runner_rejected` for REJECT_RETRY_SEC, so one bad ticker cannot stop
the other legs' refresh.

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
                                         PricerError, PricerRejected)

log = logging.getLogger("kalshi_parlay_quoter.legtable")

LegKey = tuple[str, str]          # (market_ticker, side)
NEW_LEG_POLL_SEC = 0.25


REJECT_RETRY_SEC = 300.0


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
        self._rejected: dict[LegKey, float] = {}    # key -> when rejected
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
        """Legs to refresh: asked for within leg_ttl_sec (older ones are
        forgotten) and not rejected, those with no fair yet first."""
        cutoff = self._clock() - self._leg_ttl
        with self._lock:
            expired = [k for k, asked in self._last_asked.items() if asked < cutoff]
            for key in expired:
                del self._last_asked[key]
                self._fairs.pop(key, None)
                self._rejected.pop(key, None)
            # A rejection is retried after REJECT_RETRY_SEC, so a refusal that
            # was the runner's fault (not the leg's) heals without a restart.
            retry_before = self._clock() - REJECT_RETRY_SEC
            for key in [k for k, at in self._rejected.items() if at < retry_before]:
                del self._rejected[key]
                self._fairs.pop(key, None)
            live = [k for k in self._last_asked if k not in self._rejected]
            return ([k for k in live if k not in self._fairs]
                    + [k for k in live if k in self._fairs])

    def unpriced(self) -> list[LegKey]:
        """Wanted legs with no fair at all yet."""
        with self._lock:
            return [k for k in self._last_asked
                    if k not in self._fairs and k not in self._rejected]

    def reject(self, key: LegKey, detail: str) -> None:
        """The runner refuses this leg: it declines runner_rejected, at any
        age, and is asked again after REJECT_RETRY_SEC."""
        with self._lock:
            if key not in self._last_asked:
                return
            self._rejected[key] = self._clock()
            self._fairs[key] = ({"market_ticker": key[0], "side": key[1], "ok": False,
                                 "reason": "runner_rejected", "detail": detail},
                                None, self._clock())
            self.generation += 1

    def store(self, keys: list[LegKey], legs: list[dict],
              priced_at_ms: int | None, sent_at: float) -> None:
        with self._lock:
            for key, leg in zip(keys, legs):
                if key in self._last_asked and key not in self._rejected:
                    self._fairs[key] = (leg, priced_at_ms, sent_at)
            self.generation += 1

    def missing(self, keys: list[LegKey]) -> list[LegKey]:
        """The legs with no fair younger than max_age_sec."""
        cutoff = self._clock() - self._max_age
        with self._lock:
            return [k for k in keys if self._is_stale(k, self._fairs.get(k), cutoff)]

    def _is_stale(self, key: LegKey, held, cutoff: float) -> bool:
        if held is None:
            return True
        return held[2] < cutoff and key not in self._rejected

    def price(self, rfq_id: str, keys: list[LegKey]) -> PriceResult:
        """decide.decide's price_fn. Raises LegsNotReady when any leg has no
        fresh fair. elapsed_ms is 0: nothing was fetched."""
        now = self._clock()
        cutoff = now - self._max_age
        out = []
        with self._lock:
            held = [self._fairs.get(k) for k in keys]
            stale = [k for k, h in zip(keys, held) if self._is_stale(k, h, cutoff)]
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
    """One thread: every interval_sec, re-price every wanted leg; between
    passes, every NEW_LEG_POLL_SEC, price the legs that have no fair yet.

    fetch(keys) -> PriceResult for at most FAIRS_CHUNK_LEGS keys, may raise
    PricerError. A chunk the runner could not answer (down, busy, timeout)
    ends the pass: the next pass retries and its legs simply age out of the
    table. A chunk the runner REJECTED (PricerRejected) is halved until the
    refused leg is alone, that leg is rejected in the table, and the pass
    goes on."""

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
        self._rejected_legs = 0
        self._last_pass_ms = 0.0

    def _price_chunk(self, chunk: list[LegKey]) -> None:
        """Fetch and store one chunk. Raises PricerError when the runner
        could not answer; a rejected chunk is halved down to the bad leg."""
        sent_at = self._clock()
        try:
            result = self._fetch(chunk)
        except PricerRejected as e:
            if len(chunk) == 1:
                log.warning("[legtable] the runner rejects leg %s:%s: %s", *chunk[0], e)
                self._table.reject(chunk[0], str(e))
                with self._lock:
                    self._rejected_legs += 1
                return
            half = len(chunk) // 2
            self._price_chunk(chunk[:half])
            self._price_chunk(chunk[half:])
            return
        self._table.store(chunk, result.legs, result.priced_at_ms, sent_at)

    def _price(self, keys: list[LegKey]) -> bool:
        for i in range(0, len(keys), FAIRS_CHUNK_LEGS):
            try:
                self._price_chunk(keys[i:i + FAIRS_CHUNK_LEGS])
            except PricerError as e:
                log.warning("[legtable] refresh failed at leg %d of %d: %s", i, len(keys), e)
                return False
        return True

    def refresh_once(self) -> bool:
        """One pass over the wanted legs, unpriced ones first. True when
        every chunk landed."""
        started = self._clock()
        ok = self._price(self._table.wanted())
        with self._lock:
            self._last_pass_ms = (self._clock() - started) * 1000.0
            if not ok:
                self._failed_passes += 1
        return ok

    def price_new_legs(self) -> bool:
        """Price the legs with no fair yet. True when there were none or all landed."""
        return self._price(self._table.unpriced())

    def take_stats(self) -> dict[str, int]:
        """Failed passes and newly rejected legs since the last call, and
        the last pass's duration."""
        with self._lock:
            failed, self._failed_passes = self._failed_passes, 0
            rejected, self._rejected_legs = self._rejected_legs, 0
            return {"table_refresh_failed": failed, "table_legs_rejected": rejected,
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
                answered = self.refresh_once()
                while not self._stop.wait(NEW_LEG_POLL_SEC) and self._clock() - started < self._interval:
                    # A runner that just failed is not asked again before the next pass.
                    if answered:
                        answered = self.price_new_legs()
            except Exception:
                log.exception("[legtable] refresh pass crashed")
                self._stop.wait(self._interval)
