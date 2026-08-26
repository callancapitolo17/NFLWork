"""The leg-surface ingest loop (issue #96).

Inputs: Kalshi (slate) + the books (prices). Outputs: a populated
``LegSurface`` in memory, mirrored to ``kalshi_mlb_mm_surface.duckdb``.
Side effects: live HTTP to Kalshi and the books; DELETE+INSERT of that DB's
three tables. NEVER touches the maker's state or market DBs, and never calls
into the RFQ path or the on-demand engine.

Shape: one thread per (book, route), each on its OWN cadence and its own
clock, plus a slate thread and a maintenance thread. A slow book must not
drag a fast one — DK's whole-slate scrape is 28.4s while BetMGM's structure
pass is 2.58s — and each row carries its own ``built_at`` so the staleness
gate (#99) judges books independently.

The SGPService here is a SECOND, PRIVATE instance with
``single_leg_structure_fair=True`` and a structure TTL equal to the cadence.
The maker's on-demand engine keeps its own service and its 420s structure
cache untouched; sharing one would make the surface's freshness requirement
silently change the same-game path's caching.
"""
from __future__ import annotations

import logging
import threading
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone

from kalshi_mlb_mm import config
from kalshi_mlb_mm.leg_surface import db, singles, slate, structure
from kalshi_mlb_mm.leg_surface.devig import ExclusionCounts
from kalshi_mlb_mm.leg_surface.store import LegSurface

log = logging.getLogger(__name__)


@dataclass
class PassResult:
    """What one book's pass did — the surface_refresh_log row."""
    book: str
    route: str
    started_at: datetime
    duration_sec: float
    games_attempted: int
    games_priced: int
    rungs_priced: int
    legs_written: int
    counts: ExclusionCounts = field(default_factory=ExclusionCounts)
    error_class: str | None = None


class SurfaceIngest:
    """Owns the slate, the per-book workers and the DuckDB mirror.

    Lifecycle: ``start()`` then ``stop()``. Construction does no I/O, so the
    maker (#98) can build it before deciding to run it.
    """

    def __init__(self, *, surface: LegSurface | None = None,
                 service=None, books_structure=None, books_singles=None):
        self.surface = surface or LegSurface()
        self._service = service
        self._owns_service = service is None
        self._books_structure = tuple(
            books_structure if books_structure is not None
            else config.SURFACE_BOOKS_STRUCTURE)
        self._books_singles = tuple(
            books_singles if books_singles is not None
            else config.SURFACE_BOOKS_SINGLES)
        self._slate: list = []
        self._slate_lock = threading.Lock()
        self._slate_ready = threading.Event()
        self._running = threading.Event()
        self._threads: list[threading.Thread] = []

    # ---------------------------------------------------------------- #
    # slate
    # ---------------------------------------------------------------- #
    def current_slate(self) -> list:
        with self._slate_lock:
            return list(self._slate)

    def refresh_slate(self) -> int:
        games = slate.discover_slate(
            min_minutes=config.SURFACE_GAME_MIN_MINUTES,
            max_hours=config.SURFACE_GAME_MAX_HOURS)
        if not games:
            # Keep the previous slate: an empty return is far more often a
            # Kalshi blip than a real off-day, and dropping the slate would
            # blank every book's rows on the next pass. A genuine off-day
            # ages out through built_at instead.
            log.warning("surface slate: 0 games discovered — keeping previous")
            return 0
        with self._slate_lock:
            self._slate = games
        self._slate_ready.set()
        try:
            db.record_slate(games)
        except Exception as e:
            log.warning("surface slate: mirror write failed: %s", e)
        return len(games)

    # ---------------------------------------------------------------- #
    # per-book passes
    # ---------------------------------------------------------------- #
    def _ensure_service(self):
        if self._service is not None:
            return self._service
        from kalshi_common.sgp_service import SGPService
        # structure_ttl_sec == the cadence: the structure IS the odds on this
        # route, so a cache older than one pass would serve a stale price as
        # a fresh one (the mistake kalshi_rfi avoids by passing 0.0).
        # health_db_path=None: no writes to the market DB the pricing path
        # reads — the refresh log is this loop's observability.
        self._service = SGPService(
            books=self._books_structure, health_db_path=None,
            structure_ttl_sec=config.SURFACE_CADENCE_DEFAULT_SEC,
            single_leg_structure_fair=True)
        return self._service

    def structure_pass(self, book: str) -> PassResult:
        started_at = datetime.now(timezone.utc)
        t0 = time.monotonic()
        games = self.current_slate()
        error_class = None
        rows, counts, games_priced = [], ExclusionCounts(), 0
        try:
            rows, counts, games_priced = structure.run_pass(
                book, self._ensure_service(), games,
                band_min=config.SURFACE_OVERROUND_MIN,
                band_max=config.SURFACE_OVERROUND_MAX,
                max_req_per_sec=config.SURFACE_MAX_REQ_PER_SEC_PER_BOOK)
        except Exception as e:
            # A worker must survive anything one book does. The pass lands in
            # the log with its error_class and the book keeps its previous
            # rows, which then age out under the staleness gate.
            error_class = type(e).__name__
            log.error("surface %s structure pass failed: %s", book, e)
        return self._finish(book, structure.ROUTE, started_at, t0, games,
                            rows, counts, games_priced, error_class)

    def singles_pass(self, book: str) -> PassResult:
        started_at = datetime.now(timezone.utc)
        t0 = time.monotonic()
        games = self.current_slate()
        error_class = None
        rows, counts, games_priced = [], ExclusionCounts(), 0
        try:
            rows, counts, games_priced = singles.run_pass(
                book, games,
                owned_markets=set(config.SURFACE_SINGLES_MARKETS.get(book, ())),
                band_min=config.SURFACE_OVERROUND_MIN,
                band_max=config.SURFACE_OVERROUND_MAX,
                tolerance_min=config.SURFACE_START_TOLERANCE_MIN)
        except Exception as e:
            error_class = type(e).__name__
            log.error("surface %s singles pass failed: %s", book, e)
        return self._finish(book, singles.ROUTE, started_at, t0, games,
                            rows, counts, games_priced, error_class)

    def _finish(self, book, route, started_at, t0, games, rows, counts,
                games_priced, error_class) -> PassResult:
        """Publish a pass and log it. A FAILED pass publishes nothing — the
        book keeps its previous rows and ages out — while a pass that simply
        found nothing publishes an empty slice, so a book that stopped
        offering a market does not linger at its last price forever."""
        if error_class is None:
            self.surface.publish(book, route, rows)
            try:
                db.replace_slice_rows(book, route, rows)
            except Exception as e:
                log.warning("surface %s/%s: mirror write failed: %s",
                            book, route, e)
        result = PassResult(
            book=book, route=route, started_at=started_at,
            duration_sec=time.monotonic() - t0, games_attempted=len(games),
            games_priced=games_priced, rungs_priced=len(rows) // 2,
            legs_written=len(rows), counts=counts, error_class=error_class)
        try:
            db.record_pass(result)
        except Exception as e:
            log.warning("surface %s/%s: refresh-log write failed: %s",
                        book, route, e)
        log.info("surface %s/%s: %.1fs %d/%d games, %d legs "
                 "(crossed=%d overround=%d one_sided=%d unresolved=%d "
                 "unmatched=%d ambiguous=%d)%s",
                 book, route, result.duration_sec, games_priced,
                 len(games), len(rows), counts.crossed, counts.overround,
                 counts.one_sided, counts.unresolved, counts.game_unmatched,
                 counts.game_ambiguous,
                 f" ERROR={error_class}" if error_class else "")
        return result

    # ---------------------------------------------------------------- #
    # threads
    # ---------------------------------------------------------------- #
    def _slate_loop(self):
        while self._running.is_set():
            try:
                self.refresh_slate()
            except Exception as e:
                log.error("surface slate refresh failed: %s", e)
            self._running.wait(config.SURFACE_SLATE_REFRESH_SEC)

    def _book_loop(self, book: str, route: str, cadence_sec: float):
        # Every worker waits for the first slate: a pass over an empty slate
        # would publish an empty slice and log a misleading zero.
        self._slate_ready.wait()
        while self._running.is_set():
            t0 = time.monotonic()
            if route == structure.ROUTE:
                self.structure_pass(book)
            else:
                self.singles_pass(book)
            # Cadence is a FLOOR, not a deadline: a pass that overran its
            # cadence (DK's 28.4s slate scrape against a 60s cadence, or any
            # book having a slow minute) starts the next one immediately
            # rather than piling up.
            self._running.wait(max(0.0, cadence_sec - (time.monotonic() - t0)))

    def _maintenance_loop(self):
        while self._running.is_set():
            self._running.wait(config.SURFACE_DB_FLUSH_SEC)
            if not self._running.is_set():
                return
            try:
                dropped = db.prune_refresh_log(
                    config.SURFACE_LOG_RETENTION_HOURS)
                if dropped:
                    log.debug("surface: pruned %d refresh-log rows", dropped)
            except Exception as e:
                log.warning("surface: refresh-log prune failed: %s", e)

    def worker_specs(self) -> list[tuple[str, str, float]]:
        """(book, route, cadence_sec) for every worker this ingest runs."""
        specs = [(b, structure.ROUTE, config.SURFACE_CADENCE_DEFAULT_SEC)
                 for b in self._books_structure]
        specs += [(b, singles.ROUTE,
                   config.SURFACE_CADENCE_SINGLES_SEC.get(
                       b, config.SURFACE_CADENCE_DEFAULT_SEC))
                  for b in self._books_singles]
        return specs

    def start(self) -> None:
        db.init_database()
        self._running.set()
        self._spawn("surface-slate", self._slate_loop)
        for book, route, cadence in self.worker_specs():
            self._spawn(f"surface-{book}-{route}", self._book_loop,
                        book, route, cadence)
        self._spawn("surface-maint", self._maintenance_loop)
        log.info("surface ingest started: %d workers", len(self._threads) - 2)

    def _spawn(self, name, target, *args):
        t = threading.Thread(target=target, args=args, name=name, daemon=True)
        t.start()
        self._threads.append(t)

    def stop(self, timeout_sec: float = 10.0) -> None:
        self._running.clear()
        self._slate_ready.set()      # release any worker still waiting
        for t in self._threads:
            t.join(timeout=timeout_sec)
        self._threads.clear()
        if self._owns_service and self._service is not None:
            self._service.close()
            self._service = None
        log.info("surface ingest stopped")
