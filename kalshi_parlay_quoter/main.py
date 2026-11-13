"""Dry-run parlay quoter (#283): RFQ feed -> leg table -> decision log.
NO ORDERS: nothing in this package can send a quote.

    python -m kalshi_parlay_quoter.main

Threads
  rfq-feed    feed.RFQFeed reads Kalshi's RFQ websocket; the door (scope rule
              + size floor) leaves only RFQs we could price, on a queue.
  leg-table   legtable.Refresher re-prices every leg seen in the last 10 min
              off the runner's /rfq_fairs.json, every 2 s.
  main        Quoter.step: takes RFQs off the queue and decides each from the
              leg table (memory only), then flushes the decision log.
  confirm-*   a capped number of would-quote parlays a second are re-priced
              through the runner's Kalshi-checked /rfq_price.json.

Side effects: appends to DB_PATH (db.DecisionLog: rfq_decisions, rfq_legs,
minute_counts). Every frame and decision is counted per minute; confirmed
parlays are stored in full and STORED_SAMPLE_PCT of the rest.
"""
import logging
import queue
import signal
import threading
import time
import zlib
from collections import Counter, deque
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from typing import Callable

from kalshi_common import auth_client
from kalshi_parlay_quoter import config, decide, feed, legtable, pricer
from kalshi_parlay_quoter.db import DecisionLog

log = logging.getLogger("kalshi_parlay_quoter")

# RFQs decided per step() before the loop flushes, so a full queue cannot starve the log.
MAX_RFQS_PER_STEP = 5000
SAMPLE_BUCKETS = 10_000


def stored_in_sample(rfq_id: str, sample_pct: float) -> bool:
    """Whether an unconfirmed RFQ's row is stored: a fixed hash of its id, so
    the sample is the same on every replay and unrelated to arrival order."""
    return zlib.crc32(rfq_id.encode()) % SAMPLE_BUCKETS < sample_pct * SAMPLE_BUCKETS / 100.0


def leg_keys_of(rfq: dict) -> list[legtable.LegKey]:
    return [(str(leg.get("market_ticker") or ""), str(leg.get("side") or "").lower())
            for leg in rfq.get("mve_selected_legs") or []]


class Quoter:
    """Decides the RFQs the feed kept. step() is called from one thread; the
    confirm workers only call _record (locked) and DecisionLog.add (locked)."""

    def __init__(self, rfqs: queue.Queue, table: legtable.LegTable, dlog: DecisionLog, *,
                 confirm_fn: Callable[[str, list[legtable.LegKey]], pricer.PriceResult],
                 confirms_per_sec: float = config.CONFIRMS_PER_SEC,
                 confirm_workers: int = config.CONFIRM_WORKERS,
                 pending_max_wait_sec: float = config.PENDING_MAX_WAIT_SEC,
                 max_backlog: int = config.MAX_BACKLOG,
                 stored_sample_pct: float = config.STORED_SAMPLE_PCT,
                 clock: Callable[[], float] = time.monotonic):
        self._rfqs = rfqs
        self._table = table
        self._dlog = dlog
        self._confirm_fn = confirm_fn
        self._confirm_gap_sec = 1.0 / confirms_per_sec if confirms_per_sec > 0 else None
        self._confirm_workers = confirm_workers
        self._next_confirm_at = 0.0
        self._confirms_in_flight = 0
        self._pool = ThreadPoolExecutor(max_workers=confirm_workers, thread_name_prefix="confirm")
        self._pending_max_wait = pending_max_wait_sec
        self._max_backlog = max_backlog
        self._stored_sample_pct = stored_sample_pct
        self._clock = clock
        # RFQs waiting for a leg's first fair: (rfq, seen_at, waiting since).
        self._pending: deque[tuple[dict, datetime, float]] = deque()
        self._table_generation = table.generation
        self._lock = threading.Lock()
        self._counts: Counter = Counter()

    # -- one RFQ ---------------------------------------------------------

    def _record(self, row: dict, leg_rows: list[dict], seen_at: datetime) -> None:
        """Count the decision; store it when confirmed or in the sample."""
        row["queue_ms"] = (datetime.now(timezone.utc) - seen_at).total_seconds() * 1000.0
        key = row["reason"] or row["decision"]
        if row["decision"] == "would_quote" and not row["confirmed"]:
            key = "would_quote_unconfirmed"
        with self._lock:
            self._counts[key] += 1
        if row["confirmed"] or stored_in_sample(str(row["rfq_id"]), self._stored_sample_pct):
            self._dlog.add(row, leg_rows)
            with self._lock:
                self._counts["stored"] += 1

    def _take_confirm_slot(self) -> bool:
        if self._confirm_gap_sec is None:
            return False
        now = self._clock()
        with self._lock:
            if now < self._next_confirm_at or self._confirms_in_flight >= self._confirm_workers:
                return False
            self._next_confirm_at = now + self._confirm_gap_sec
            self._confirms_in_flight += 1
            return True

    def _confirm(self, rfq: dict, seen_at: datetime) -> None:
        try:
            row, leg_rows = decide.decide(rfq, rfq.get("mve_selected_legs"),
                                          self._confirm_fn, seen_at=seen_at)
            row["confirmed"] = True
            self._record(row, leg_rows, seen_at)
        except Exception:
            log.exception("[quoter] confirm failed for rfq %s", rfq.get("id"))
            with self._lock:
                self._counts["confirm_crashed"] += 1
        finally:
            with self._lock:
                self._confirms_in_flight -= 1

    def _decide(self, rfq: dict, seen_at: datetime, price_fn) -> None:
        row, leg_rows = decide.decide(rfq, rfq.get("mve_selected_legs"), price_fn,
                                      seen_at=seen_at)
        if row["decision"] == "would_quote" and self._take_confirm_slot():
            self._pool.submit(self._confirm, rfq, seen_at)
            return
        self._record(row, leg_rows, seen_at)

    @staticmethod
    def _backlogged(rfq_id, legs):
        raise pricer.PricerBacklog("")

    def handle(self, rfq: dict, seen_at: datetime, waiting_since: float | None = None) -> None:
        """Decide one RFQ, or park it until the next table refresh when a leg
        has no fair yet."""
        keys = leg_keys_of(rfq)
        self._table.want(keys)
        if not self._table.missing(keys):
            self._decide(rfq, seen_at, self._table.price)
            return
        now = self._clock()
        if waiting_since is None:
            if len(self._pending) >= self._max_backlog:
                self._decide(rfq, seen_at, self._backlogged)
                return
            self._pending.append((rfq, seen_at, now))
        elif now - waiting_since >= self._pending_max_wait:
            # table.price raises LegsNotReady -> decline pricer_error, with the leg named.
            self._decide(rfq, seen_at, self._table.price)
        else:
            self._pending.append((rfq, seen_at, waiting_since))

    # -- loop ------------------------------------------------------------

    def _retry_pending(self) -> None:
        for _ in range(len(self._pending)):
            rfq, seen_at, waiting_since = self._pending.popleft()
            self._handle_safely(rfq, seen_at, waiting_since)

    def _handle_safely(self, rfq: dict, seen_at: datetime, waiting_since: float | None = None) -> None:
        try:
            self.handle(rfq, seen_at, waiting_since)
        except Exception:
            log.exception("[quoter] failed handling rfq %s", rfq.get("id"))
            with self._lock:
                self._counts["handle_crashed"] += 1

    def step(self, wait_sec: float = 0.0) -> int:
        """Decide what is on the queue (waiting up to wait_sec for the first
        RFQ) and retry parked RFQs after a table refresh or when the oldest
        has waited out. Returns how many RFQs came off the queue."""
        taken = 0
        try:
            rfq = self._rfqs.get(timeout=wait_sec) if wait_sec > 0 else self._rfqs.get_nowait()
            while True:
                taken += 1
                self._handle_safely(rfq, datetime.now(timezone.utc))
                if taken >= MAX_RFQS_PER_STEP:
                    break
                rfq = self._rfqs.get_nowait()
        except queue.Empty:
            pass
        oldest_expired = bool(self._pending) and (
            self._clock() - self._pending[0][2] >= self._pending_max_wait)
        if self._table.generation != self._table_generation or oldest_expired:
            self._table_generation = self._table.generation
            self._retry_pending()
        return taken

    def take_counts(self) -> dict[str, int]:
        """Decision counts since the last call, plus gauges (current sizes)."""
        with self._lock:
            counts, self._counts = self._counts, Counter()
        out = dict(counts)
        out["gauge_pending"] = len(self._pending)
        out["gauge_queue"] = self._rfqs.qsize()
        out["gauge_table_legs"] = len(self._table)
        return out

    def shutdown(self) -> None:
        """Queued confirms are dropped so SIGTERM ends inside Docker's 10 s
        grace; one already at the runner finishes (its timeout is 3 s)."""
        self._pool.shutdown(wait=True, cancel_futures=True)


def run() -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s %(message)s")
    if not (config.KALSHI_API_KEY_ID and config.KALSHI_PRIVATE_KEY_PATH):
        raise SystemExit("KALSHI_API_KEY_ID and KALSHI_PRIVATE_KEY_PATH are required")
    auth_client.configure(config.KALSHI_API_KEY_ID, config.KALSHI_PRIVATE_KEY_PATH,
                          config.KALSHI_BASE_URL, config.PROJECT_ROOT)
    config.DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    dlog = DecisionLog(config.DB_PATH)
    rfq_feed = feed.RFQFeed(
        feed.kalshi_connect(config.KALSHI_WS_URL, auth_client.ws_auth_headers),
        min_cost_usd=config.MIN_RFQ_TARGET_COST_USD, queue_max=config.FEED_QUEUE_MAX,
        idle_timeout_sec=config.FEED_IDLE_TIMEOUT_SEC)
    table = legtable.LegTable(max_age_sec=config.TABLE_MAX_AGE_SEC,
                              leg_ttl_sec=config.TABLE_LEG_TTL_SEC)
    refresher = legtable.Refresher(
        table,
        lambda keys: pricer.fairs(config.RUNNER_URL, keys, timeout=config.TABLE_FETCH_TIMEOUT_SEC),
        interval_sec=config.TABLE_REFRESH_SEC)
    quoter = Quoter(
        rfq_feed.rfqs, table, dlog,
        confirm_fn=lambda rfq_id, keys: pricer.price(config.RUNNER_URL, rfq_id, keys,
                                                     timeout=config.CONFIRM_TIMEOUT_SEC))
    stop = threading.Event()
    for sig in (signal.SIGTERM, signal.SIGINT):
        signal.signal(sig, lambda *_: stop.set())

    rfq_feed.start()
    refresher.start()
    log.info("[quoter] DRY RUN: no orders. db=%s runner=%s floor=$%.0f",
             config.DB_PATH, config.RUNNER_URL, config.MIN_RFQ_TARGET_COST_USD)
    last_flush = time.monotonic()
    minute = datetime.now(timezone.utc).replace(second=0, microsecond=0)
    paused = False
    while not stop.is_set():
        stop_file = config.STOP_PATH.exists()
        if stop_file != paused:
            paused = stop_file
            rfq_feed.set_paused(paused)
            log.info("[quoter] stop file %s %s", config.STOP_PATH,
                     "present: feed closed, idling" if paused else "gone: resuming")
        if paused:
            stop.wait(5)
            continue
        quoter.step(wait_sec=0.25)
        now_min = datetime.now(timezone.utc).replace(second=0, microsecond=0)
        if now_min > minute:
            counts = {**rfq_feed.take_counts(), **refresher.take_stats(), **quoter.take_counts()}
            dlog.add_counts(minute, counts)
            log.info("[quoter] minute %s: %s", minute.strftime("%H:%M"),
                     dict(sorted(counts.items())))
            minute = now_min
            try:
                dlog.flush()
                dlog.checkpoint()
            except Exception:
                log.exception("[quoter] decision-log checkpoint failed")
            last_flush = time.monotonic()
        if time.monotonic() - last_flush >= config.FLUSH_SEC:
            try:
                dlog.flush()
            except Exception:
                log.exception("[quoter] decision-log flush failed; batch dropped")
            last_flush = time.monotonic()

    log.info("[quoter] stopping")
    rfq_feed.stop()
    refresher.stop()
    quoter.shutdown()
    dlog.close()


if __name__ == "__main__":
    run()
