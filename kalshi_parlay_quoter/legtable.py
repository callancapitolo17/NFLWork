"""In-memory table of leg fairs, so pricing an RFQ is a dict lookup.

Why (2026-11-13): ~475 in-scope RFQs a second is ~1,900 leg prices a second,
far more than one runner call per RFQ could carry, and they re-ask the same
few thousand legs. The table holds every (ticker, side) seen in the last
`leg_ttl_sec` and prices an RFQ from memory, so the work grows with distinct
legs, not with RFQs.

A leg is priced when it changes, not on a clock. A leg's fair moves only
when Unabated's snapshot for its league lands on the runner (about once a
minute). Refresher reads each league's snapshot time off the runner's
/health every couple of seconds, which the runner answers without touching
the thread its taker runs on, and asks for fairs (pricer.fairs: board only,
no Kalshi read) only for:
  - legs with no fair yet (first, and between polls, so a new leg waits a
    fraction of a second),
  - legs whose league's snapshot changed since they were priced,
  - legs priced over `max_reuse_sec` ago (a backstop for anything a snapshot
    time does not show).
The first version re-priced every leg every 2 s; ~97% of those calls
returned the same number and they loaded the taker's thread (2026-11-13: the
taker froze 40 minutes after the quoter went live).

A fair is served only while (a) the runner's snapshot times were read under
`max_age_sec` ago, (b) the leg was priced on its league's current snapshot,
or that snapshot changed under `max_age_sec` ago and the re-price is on its
way, and (c) the leg was priced under `max_reuse_sec + max_age_sec` ago.
Every served leg carries priced_at_ms (the runner's stamp) and cache_age_ms,
and fair_age_ms grows by the cache age, so the log never shows a reused fair
as fresh. Table legs are not Kalshi-checked (kalshi_* fields are null).

A leg the runner refuses outright is found by halving its chunk and declines
`runner_rejected` for REJECT_RETRY_SEC, so one bad ticker cannot stop the
other legs' refresh.
"""
import logging
import threading
import time
from typing import Callable

from kalshi_parlay_quoter.pricer import (FAIRS_CHUNK_LEGS, FairsResult, LegsNotReady,
                                         PriceResult, PricerError, PricerRejected)

log = logging.getLogger("kalshi_parlay_quoter.legtable")

LegKey = tuple[str, str]          # (market_ticker, side)
NEW_LEG_POLL_SEC = 0.25
REJECT_RETRY_SEC = 300.0
# Between two chunks of one pass: lets the runner's thread take other work
# when a whole league (a few thousand legs) is re-priced at once.
CHUNK_PAUSE_SEC = 0.05


class LegTable:
    def __init__(self, *, max_age_sec: float, leg_ttl_sec: float, max_reuse_sec: float,
                 league_of: Callable[[LegKey], str | None],
                 clock: Callable[[], float] = time.monotonic):
        self._max_age = max_age_sec
        self._leg_ttl = leg_ttl_sec
        self._max_reuse = max_reuse_sec
        self._league_of = league_of
        self._clock = clock
        self._lock = threading.Lock()
        self._last_asked: dict[LegKey, float] = {}
        # key -> (leg as the runner answered, runner priced_at_ms, sent_at,
        #         its league's snapshot time when it was priced)
        self._fairs: dict[LegKey, tuple[dict, int | None, float, int | None]] = {}
        self._rejected: dict[LegKey, float] = {}    # key -> when rejected
        # league label -> the runner's snapshot time (ms), and when we first saw it.
        self._versions: dict[str, int] = {}
        self._changed_at: dict[str, float] = {}
        self._versions_read_at: float | None = None
        self.generation = 0        # bumps on every store(); waiters retry on change

    def __len__(self) -> int:
        return len(self._last_asked)

    def want(self, keys: list[LegKey]) -> None:
        """Keep these legs priced for leg_ttl_sec from now."""
        now = self._clock()
        with self._lock:
            for key in keys:
                self._last_asked[key] = now

    def note_versions(self, versions: dict[str, int], *, read: bool = True) -> None:
        """The runner's snapshot time per league. A time only moves forward:
        /health can be a few seconds behind a fairs answer.
          read  True for a /health read (it also proves the runner is up)."""
        now = self._clock()
        with self._lock:
            for league, version in versions.items():
                if version > self._versions.get(league, -1):
                    self._versions[league] = version
                    self._changed_at[league] = now
            if read:
                self._versions_read_at = now

    def _forget_expired(self, now: float) -> None:
        cutoff = now - self._leg_ttl
        for key in [k for k, asked in self._last_asked.items() if asked < cutoff]:
            del self._last_asked[key]
            self._fairs.pop(key, None)
            self._rejected.pop(key, None)
        # A rejection is retried after REJECT_RETRY_SEC, so a refusal that
        # was the runner's fault (not the leg's) heals without a restart.
        retry_before = now - REJECT_RETRY_SEC
        for key in [k for k, at in self._rejected.items() if at < retry_before]:
            del self._rejected[key]
            self._fairs.pop(key, None)

    def _on_current_snapshot(self, key: LegKey, held) -> bool:
        return held[3] == self._versions.get(self._league_of(key))

    def due(self) -> list[LegKey]:
        """Legs to price now, those with no fair first: asked for within
        leg_ttl_sec (older ones are forgotten), not rejected, and unpriced,
        priced on an older snapshot of their league, or priced over
        max_reuse_sec ago."""
        now = self._clock()
        with self._lock:
            self._forget_expired(now)
            unpriced, changed = [], []
            for key in self._last_asked:
                if key in self._rejected:
                    continue
                held = self._fairs.get(key)
                if held is None:
                    unpriced.append(key)
                elif not self._on_current_snapshot(key, held) or now - held[2] > self._max_reuse:
                    changed.append(key)
            return unpriced + changed

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
                                None, self._clock(), None)
            self.generation += 1

    def store(self, keys: list[LegKey], legs: list[dict], priced_at_ms: int | None,
              sent_at: float, versions: dict[str, int] | None = None) -> None:
        """Fairs the runner just answered, and the snapshot times they were
        priced on (none given: the times already known)."""
        if versions:
            self.note_versions(versions, read=False)
        with self._lock:
            for key, leg in zip(keys, legs):
                if key in self._last_asked and key not in self._rejected:
                    self._fairs[key] = (leg, priced_at_ms, sent_at,
                                        self._versions.get(self._league_of(key)))
            self.generation += 1

    def _not_servable(self, key: LegKey, held, now: float) -> bool:
        if held is None:
            return True
        if key in self._rejected:
            return False
        if self._versions_read_at is None or now - self._versions_read_at > self._max_age:
            return True     # the runner has not been heard from: nothing is known fresh
        if now - held[2] > self._max_reuse + self._max_age:
            return True
        if self._on_current_snapshot(key, held):
            return False
        changed_at = self._changed_at.get(self._league_of(key))
        return changed_at is None or now - changed_at > self._max_age

    def missing(self, keys: list[LegKey]) -> list[LegKey]:
        """The legs with no fair that may be served now."""
        now = self._clock()
        with self._lock:
            return [k for k in keys if self._not_servable(k, self._fairs.get(k), now)]

    def price(self, rfq_id: str, keys: list[LegKey]) -> PriceResult:
        """decide.decide's price_fn. Raises LegsNotReady when any leg has no
        servable fair. elapsed_ms is 0: nothing was fetched."""
        now = self._clock()
        out = []
        with self._lock:
            held = [self._fairs.get(k) for k in keys]
            stale = [k for k, h in zip(keys, held) if self._not_servable(k, h, now)]
        if stale:
            raise LegsNotReady(f"no current fair for {len(stale)} of {len(keys)} legs, "
                               f"e.g. {stale[0][0]}")
        for leg, priced_at_ms, sent_at, _version in held:
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
    """One thread: every interval_sec, read the runner's snapshot times and
    price the legs that are due (LegTable.due); between those passes, every
    NEW_LEG_POLL_SEC, price the legs that have no fair yet.

    read_versions() -> {league: snapshot ms}, may raise PricerError.
    fetch(keys) -> FairsResult for at most FAIRS_CHUNK_LEGS keys, may raise
    PricerError. A read or a chunk the runner could not answer (down, busy,
    timeout) ends the pass: the next pass retries, and fairs stop being
    served once the snapshot times are max_age_sec old. A chunk the runner
    REJECTED (PricerRejected) is halved until the refused leg is alone, that
    leg is rejected in the table, and the pass goes on."""

    def __init__(self, table: LegTable, fetch: Callable[[list[LegKey]], FairsResult],
                 read_versions: Callable[[], dict[str, int]], *,
                 interval_sec: float, clock: Callable[[], float] = time.monotonic,
                 chunk_pause_sec: float = CHUNK_PAUSE_SEC):
        self._table = table
        self._fetch = fetch
        self._read_versions = read_versions
        self._interval = interval_sec
        self._clock = clock
        self._chunk_pause = chunk_pause_sec
        self._stop = threading.Event()
        self._paused = threading.Event()
        self._thread: threading.Thread | None = None
        self._lock = threading.Lock()
        self._failed_passes = 0
        self._rejected_legs = 0
        self._legs_priced = 0
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
        self._table.store(chunk, result.legs, result.priced_at_ms, sent_at, result.league_loaded_at)
        with self._lock:
            self._legs_priced += len(chunk)

    def _price(self, keys: list[LegKey]) -> bool:
        for i in range(0, len(keys), FAIRS_CHUNK_LEGS):
            if self._paused.is_set() or self._stop.is_set():
                return True     # a pause ends the pass at the next chunk
            if i and self._chunk_pause:
                self._stop.wait(self._chunk_pause)
            try:
                self._price_chunk(keys[i:i + FAIRS_CHUNK_LEGS])
            except PricerError as e:
                log.warning("[legtable] refresh failed at leg %d of %d: %s", i, len(keys), e)
                return False
        return True

    def refresh_once(self) -> bool:
        """One pass: the runner's snapshot times, then every due leg, unpriced
        ones first. True when the read and every chunk landed."""
        started = self._clock()
        try:
            self._table.note_versions(self._read_versions())
            ok = self._price(self._table.due())
        except PricerError as e:
            log.warning("[legtable] the runner's snapshot times were not read: %s", e)
            ok = False
        with self._lock:
            self._last_pass_ms = (self._clock() - started) * 1000.0
            if not ok:
                self._failed_passes += 1
        return ok

    def price_new_legs(self) -> bool:
        """Price the legs with no fair yet. True when there were none or all landed."""
        return self._price(self._table.unpriced())

    def take_stats(self) -> dict[str, int]:
        """Failed passes, newly rejected legs and legs priced since the last
        call, and the last pass's duration."""
        with self._lock:
            failed, self._failed_passes = self._failed_passes, 0
            rejected, self._rejected_legs = self._rejected_legs, 0
            priced, self._legs_priced = self._legs_priced, 0
            return {"table_refresh_failed": failed, "table_legs_rejected": rejected,
                    "table_legs_priced": priced,
                    "gauge_table_refresh_ms": int(self._last_pass_ms)}

    def start(self) -> None:
        if self._thread is None:
            self._thread = threading.Thread(target=self._run, name="leg-table", daemon=True)
            self._thread.start()

    def stop(self) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=2.0)

    def set_paused(self, paused: bool) -> None:
        """While paused, nothing is asked of the runner (the stop file must
        take the quoter's whole load off it, not only the feed's)."""
        if paused:
            self._paused.set()
        else:
            self._paused.clear()

    def _run(self) -> None:
        while not self._stop.is_set():
            if self._paused.is_set():
                self._stop.wait(NEW_LEG_POLL_SEC)
                continue
            started = self._clock()
            try:
                answered = self.refresh_once()
                while (not self._stop.wait(NEW_LEG_POLL_SEC) and not self._paused.is_set()
                       and self._clock() - started < self._interval):
                    # A runner that just failed is not asked again before the next pass.
                    if answered:
                        answered = self.price_new_legs()
            except Exception:
                log.exception("[legtable] refresh pass crashed")
                self._stop.wait(self._interval)
