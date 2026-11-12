"""Dry-run parlay quoter (#283): RFQ feed -> scope rule -> runner leg pricer
-> decision log. NO ORDERS: nothing in this package can send a quote.

    python -m kalshi_parlay_quoter.main

Feed: kalshi_mlb_mm's WebSocketRFQSource with its door widened to the
leagues in scope.LEAGUES and no leg cap. Every RFQ id gets exactly one
decision row, except those the door drops (other leagues, below the size
floor), which are counted per minute in door_counts.
"""
import logging
import signal
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone

from kalshi_common import auth_client
from kalshi_parlay_quoter import config, decide, pricer, scope
from kalshi_parlay_quoter.db import DecisionLog

log = logging.getLogger("kalshi_parlay_quoter")


class _RateLimit:
    """At most `per_sec` acquisitions per rolling second; never blocks."""

    def __init__(self, per_sec: float):
        self._per_sec = per_sec
        self._stamps: list[float] = []
        self._lock = threading.Lock()

    def try_acquire(self) -> bool:
        now = time.monotonic()
        with self._lock:
            self._stamps = [t for t in self._stamps if now - t < 1.0]
            if len(self._stamps) >= self._per_sec:
                return False
            self._stamps.append(now)
            return True


class Quoter:
    def __init__(self, source, dlog: DecisionLog, *, price_fn=None,
                 workers: int = config.PRICER_WORKERS,
                 seen_ttl_sec: float = config.SEEN_TTL_SEC,
                 market_fetch_rate: float = config.MAX_MARKET_FETCHES_PER_SEC,
                 max_backlog: int = config.MAX_BACKLOG,
                 min_cost_usd: float = config.MIN_RFQ_TARGET_COST_USD):
        self.source = source
        self.dlog = dlog
        self.price_fn = price_fn or pricer.CachingPricer(
            lambda rid, legs: pricer.price(config.RUNNER_URL, rid, legs,
                                           timeout=config.PRICER_TIMEOUT_SEC),
            ttl_sec=config.PRICE_CACHE_SEC, backoff_sec=config.PRICER_BACKOFF_SEC)
        self._max_backlog = max_backlog
        self._min_cost = min_cost_usd
        self._below_floor = 0
        self._pending = 0
        self._pool = ThreadPoolExecutor(max_workers=workers,
                                        thread_name_prefix="pricer")
        self._seen: dict[str, float] = {}
        self._seen_ttl = seen_ttl_sec
        self._last_prune = time.monotonic()
        self._fetch_limit = _RateLimit(market_fetch_rate)
        self._lock = threading.Lock()
        # Local counters for the per-minute door row.
        self._logged = 0
        self._other_league = 0
        self._door_base = dict(getattr(source, "ingestion_drops", {}) or {})
        self._malformed_base = getattr(source, "malformed_count", 0)
        self.decisions: dict[str, int] = {}

    # -- per RFQ -------------------------------------------------------------

    def _legs_for(self, rfq: dict):
        legs = rfq.get("mve_selected_legs")
        if isinstance(legs, list) and legs:
            return legs, None
        ticker = rfq.get("market_ticker")
        if not ticker:
            return None, "no_combo_ticker"
        if not self._fetch_limit.try_acquire():
            return None, "market_fetch_budget"
        market = self.source.get_market(ticker)
        if not market:
            return None, "market_fetch_failed"
        legs = market.get("mve_selected_legs")
        return (legs, None) if isinstance(legs, list) and legs else (None, "no_mve_legs")

    @staticmethod
    def _backlogged(rid, legs):
        raise pricer.PricerBacklog("")

    def _handle_pooled(self, rfq: dict, seen_at: datetime) -> None:
        try:
            self.handle(rfq, seen_at)
        finally:
            with self._lock:
                self._pending -= 1

    def handle(self, rfq: dict, seen_at: datetime | None = None,
               price_fn=None) -> None:
        """seen_at is when tick() first saw the RFQ, so a pricer backlog
        shows up as queue_ms instead of hiding inside seen_at."""
        try:
            seen_at = seen_at or datetime.now(timezone.utc)
            legs, why = self._legs_for(rfq)
            if legs is not None and not scope.in_league(legs):
                # REST fallback is unfiltered; count like the WS door does.
                with self._lock:
                    self._other_league += 1
                return
            row, leg_rows = decide.decide(rfq, legs, price_fn or self.price_fn,
                                          seen_at=seen_at)
            if row["n_legs"] and row["reason"] in scope.SCOPE_REASONS:
                # Out-of-scope legs say nothing a margin needs, and at
                # ~100 RFQs/s they were most of the DB; reason_detail keeps
                # the series / game that failed.
                leg_rows = []
            row["queue_ms"] = (datetime.now(timezone.utc) - seen_at).total_seconds() * 1000.0
            if legs is None and why:
                row["reason_detail"] = why
            self.dlog.add(row, leg_rows)
            key = row["reason"] or row["decision"]
            with self._lock:
                self._logged += 1
                self.decisions[key] = self.decisions.get(key, 0) + 1
        except Exception:
            log.exception("[quoter] failed handling rfq %s", rfq.get("id"))

    # -- loop ----------------------------------------------------------------

    def tick(self) -> int:
        """Submit every not-yet-decided open RFQ. Returns how many."""
        now = time.monotonic()
        submitted = 0
        for rfq in self.source.poll():
            rid = rfq.get("id")
            if not rid or rid in self._seen:
                continue
            self._seen[rid] = now
            if self._below_contract_floor(rfq):
                with self._lock:
                    self._below_floor += 1
                continue
            seen_at = datetime.now(timezone.utc)
            with self._lock:
                full = self._pending >= self._max_backlog
                if not full:
                    self._pending += 1
            if full:
                # Decided here, without the pricer: scope still runs, so an
                # out-of-scope RFQ keeps its real reason.
                self.handle(rfq, seen_at, price_fn=self._backlogged)
            else:
                self._pool.submit(self._handle_pooled, rfq, seen_at)
            submitted += 1
        if now - self._last_prune > 60 or len(self._seen) > 200_000:
            cutoff = now - self._seen_ttl
            self._seen = {k: t for k, t in self._seen.items() if t >= cutoff}
            self._last_prune = now
        return submitted

    def _below_contract_floor(self, rfq: dict) -> bool:
        """The WS door's $ floor skips contract-sized RFQs (no dollar target),
        which were 63% of in-league RFQs on 2026-11-07. A contract pays at
        most $1, so contracts < floor means the RFQ's cost is under the floor
        whatever the price."""
        if self._min_cost <= 0:
            return False
        try:
            if float(rfq.get("target_cost_dollars") or 0) > 0:
                return False          # dollar-sized: the WS door judged it
            contracts = float(rfq.get("contracts_fp") or 0)
        except (TypeError, ValueError):
            return False
        return 0 < contracts < self._min_cost

    def door_row(self, minute: datetime) -> dict:
        drops = dict(getattr(self.source, "ingestion_drops", {}) or {})
        delta = {k: drops.get(k, 0) - self._door_base.get(k, 0) for k in drops}
        self._door_base = drops
        malformed = getattr(self.source, "malformed_count", 0)
        mal_delta, self._malformed_base = malformed - self._malformed_base, malformed
        with self._lock:
            logged, other, floor = self._logged, self._other_league, self._below_floor
            self._logged = self._other_league = self._below_floor = 0
        return {"minute": minute, "logged": logged,
                "backlog": self._pending,
                "other_league": other + delta.get("non_mlb", 0),
                "below_floor": floor + delta.get("below_floor", 0),
                "too_many_legs": delta.get("too_many_legs", 0),
                "malformed": mal_delta,
                "feed_mode": getattr(self.source, "_mode", None)}

    def shutdown(self, cancel_pending: bool = False) -> None:
        """cancel_pending drops queued (not running) RFQs, so SIGTERM ends
        within Docker's 10 s grace instead of draining the queue; the first
        VM deploy was SIGKILLed (exit 137) doing that."""
        self._pool.shutdown(wait=True, cancel_futures=cancel_pending)


def _build_source():
    from kalshi_mlb_mm import ws_rfq_source
    from kalshi_mlb_mm.rfq_source import RestRFQSource
    return ws_rfq_source.WebSocketRFQSource(
        rest_source=RestRFQSource(),
        door_leg_prefixes=scope.LEAGUE_PREFIXES,
        door_max_legs=10**6,                       # no leg cap (Cal, 2026-11-07)
        door_min_cost_usd=config.MIN_RFQ_TARGET_COST_USD,
        mirror_ttl_sec=config.MIRROR_TTL_SEC)


def run() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s %(message)s")
    if not (config.KALSHI_API_KEY_ID and config.KALSHI_PRIVATE_KEY_PATH):
        raise SystemExit("KALSHI_API_KEY_ID and KALSHI_PRIVATE_KEY_PATH are required")
    auth_client.configure(config.KALSHI_API_KEY_ID, config.KALSHI_PRIVATE_KEY_PATH,
                          config.KALSHI_BASE_URL, config.PROJECT_ROOT)
    config.DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    dlog = DecisionLog(config.DB_PATH)
    source = _build_source()
    quoter = Quoter(source, dlog)
    stop = threading.Event()
    for sig in (signal.SIGTERM, signal.SIGINT):
        signal.signal(sig, lambda *_: stop.set())

    source.start()
    log.info("[quoter] DRY RUN: no orders. db=%s runner=%s floor=$%.0f",
             config.DB_PATH, config.RUNNER_URL, config.MIN_RFQ_TARGET_COST_USD)
    last_flush = time.monotonic()
    minute = datetime.now(timezone.utc).replace(second=0, microsecond=0)
    while not stop.is_set():
        if config.STOP_PATH.exists():
            log.info("[quoter] stop file %s present; idling", config.STOP_PATH)
            stop.wait(5)
            continue
        try:
            quoter.tick()
        except Exception:
            log.exception("[quoter] tick failed")
        now_min = datetime.now(timezone.utc).replace(second=0, microsecond=0)
        if now_min > minute:
            row = quoter.door_row(minute)
            dlog.add_door(row)
            log.info("[quoter] minute %s: %s | decisions so far %s",
                     minute.strftime("%H:%M"), row, quoter.decisions)
            minute = now_min
        if time.monotonic() - last_flush >= config.FLUSH_SEC:
            try:
                dlog.flush()
            except Exception:
                log.exception("[quoter] decision-log flush failed; batch dropped")
            last_flush = time.monotonic()
        stop.wait(config.TICK_SEC)

    log.info("[quoter] stopping")
    source.stop()
    quoter.shutdown(cancel_pending=True)
    dlog.close()


if __name__ == "__main__":
    run()
