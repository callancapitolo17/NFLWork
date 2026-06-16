"""In-process SGP pricing service shared by the Kalshi MLB bots.

Replaces the subprocess-per-cycle scraper model for the taker
(kalshi_mlb_rfq) and maker (kalshi_mlb_mm): one persistent HTTP client
per book held across cycles (no per-cycle TLS handshake), the four book
orchestrators run concurrently under a per-book deadline, and slow
structure fetches (event lists, DK's 2MB selection-id payload) are
TTL-cached. Prices are NEVER cached — every refresh() re-prices.

The dashboard keeps the CLI-shim subprocess path and never touches this.
"""
from __future__ import annotations
import logging
import sys
import time
from concurrent.futures import ThreadPoolExecutor, TimeoutError as FutureTimeout
from pathlib import Path

from mlb_sgp._shared import PricedRow, TargetLine, TTLCache

log = logging.getLogger(__name__)

_REPO_ROOT = Path(__file__).resolve().parents[1]
_MLB_SGP_DIR = _REPO_ROOT / "mlb_sgp"

DEFAULT_BOOKS = ("draftkings", "fanduel", "prophetx", "novig")

# TTLs for structure fetches (market *structure*, never prices).
EVENTS_TTL_SEC = 300       # game list churns ~daily
STRUCTURE_TTL_SEC = 180    # sel-ids / runners churn when a book re-mains a line


class _BookState:
    __slots__ = ("client", "failures", "last_success", "caches")

    def __init__(self):
        self.client = None          # persistent per-book HTTP client
        self.failures = 0           # consecutive refresh failures
        self.last_success = None    # now_fn timestamp of last success
        self.caches = {}            # name -> TTLCache


class SGPService:
    """Holds persistent book clients; prices targets at all due books
    concurrently. See refresh() for the result contract."""

    MAX_FAILURES_BEFORE_REINIT = 3

    def __init__(
        self,
        books: tuple[str, ...] = DEFAULT_BOOKS,
        per_book_deadline_sec: float = 75.0,
        min_refresh_sec: dict[str, float] | None = None,
        runners: dict | None = None,   # test seam: {book: callable(targets)->rows}
        now_fn=time.monotonic,
    ):
        self.books = tuple(books)
        self.per_book_deadline_sec = per_book_deadline_sec
        self.min_refresh_sec = dict(min_refresh_sec or {})
        self._runners = runners
        self._now = now_fn
        self._state = {b: _BookState() for b in self.books}
        # The orchestrators lazily import the legacy scraper modules by
        # top-level name (`from scraper_draftkings_sgp import ...`),
        # which only resolves with mlb_sgp/ itself on sys.path — true
        # for CLI runs (cwd=mlb_sgp/) but not for the bots (cwd=repo
        # root). Make it resolvable here, once.
        if str(_MLB_SGP_DIR) not in sys.path:
            sys.path.insert(0, str(_MLB_SGP_DIR))

    # ------------------------------------------------------------------ #
    # Public API                                                          #
    # ------------------------------------------------------------------ #

    def refresh(self, targets: list[TargetLine]) -> dict[str, list[PricedRow] | None]:
        """Price `targets` at every due book concurrently.

        Returns {book: result} for books ATTEMPTED this call:
          list[PricedRow] (possibly empty) -> success
          None                             -> failure or deadline timeout
        Books skipped by min_refresh_sec are ABSENT from the dict —
        callers must not clear their previously-written rows.
        """
        due = [b for b in self.books if self._due(b)]
        results: dict[str, list[PricedRow] | None] = {}
        if not due:
            return results
        # Fresh pool per refresh: a hung book thread from a previous
        # cycle must not occupy a worker slot forever. shutdown(wait=
        # False) lets a still-running (timed-out) thread finish in the
        # background; its client gets torn down via the failure path.
        pool = ThreadPoolExecutor(max_workers=len(due),
                                  thread_name_prefix="sgp-book")
        try:
            futs = {b: pool.submit(self._run_book_safe, b, targets) for b in due}
            # All futures are already running concurrently (one worker per
            # due book), so per-book `remaining` is wall-bounded by the
            # deadline, not summed across books.
            wall_deadline = time.monotonic() + self.per_book_deadline_sec
            for b, fut in futs.items():
                remaining = max(0.0, wall_deadline - time.monotonic())
                timed_out = False
                try:
                    rows = fut.result(timeout=remaining)
                except FutureTimeout:
                    # The worker thread is still running and still holds
                    # this book's curl_cffi session (not thread-safe).
                    # Drop the client NOW so the next cycle builds a fresh
                    # session the lingering thread no longer shares —
                    # don't wait for the 3-strike teardown.
                    rows = None
                    timed_out = True
                except Exception:
                    rows = None
                self._book_done(b, rows, timed_out=timed_out)
                results[b] = rows
        finally:
            pool.shutdown(wait=False, cancel_futures=True)
        return results

    def close(self):
        """Drop all persistent clients (sessions close on GC)."""
        for st in self._state.values():
            st.client = None
            st.caches = {}

    # ------------------------------------------------------------------ #
    # Internals                                                           #
    # ------------------------------------------------------------------ #

    def _due(self, book: str) -> bool:
        min_s = self.min_refresh_sec.get(book, 0)
        st = self._state[book]
        if not min_s or st.last_success is None:
            return True
        return (self._now() - st.last_success) >= min_s

    def _book_done(self, book: str, rows, timed_out: bool = False) -> None:
        st = self._state[book]
        if rows is None:
            st.failures += 1
            log.warning("sgp_service: %s failed (consecutive=%d, timeout=%s)",
                        book, st.failures, timed_out)
            if timed_out or st.failures >= self.MAX_FAILURES_BEFORE_REINIT:
                log.warning("sgp_service: %s client torn down for reinit "
                            "(timeout=%s)", book, timed_out)
                st.client = None
                st.caches = {}
                st.failures = 0
        else:
            st.failures = 0
            st.last_success = self._now()

    def _run_book_safe(self, book: str, targets):
        """Runs in a worker thread. Returns rows or None — never raises
        (a raise would surface as a generic future error and lose the
        book attribution in logs)."""
        try:
            if self._runners is not None:
                return self._runners[book](targets)
            return self._run_book(book, targets)
        except Exception as e:
            log.warning("sgp_service: %s runner error: %s", book, e)
            return None

    def _run_book(self, book: str, targets):
        st = self._state[book]
        if book == "draftkings":
            from mlb_sgp import draftkings as mod
            if st.client is None:
                from mlb_sgp.dk_client import DraftKingsClient
                st.client = DraftKingsClient()
                st.caches = {"events": TTLCache(EVENTS_TTL_SEC),
                             "structure": TTLCache(STRUCTURE_TTL_SEC)}
            return mod.price_sgps(targets, periods=("FG",), client=st.client,
                                  fetchers=self._dk_fetchers(st))
        if book == "fanduel":
            from mlb_sgp import fanduel as mod
            if st.client is None:
                from mlb_sgp.fd_client import FanDuelClient
                st.client = FanDuelClient()
                st.caches = {"events": TTLCache(EVENTS_TTL_SEC),
                             "structure": TTLCache(STRUCTURE_TTL_SEC)}
            return mod.price_sgps(targets, periods=("FG",), client=st.client,
                                  fetchers=self._fd_fetchers(st))
        if book == "prophetx":
            from mlb_sgp import prophetx as mod
            if st.client is None:
                from mlb_sgp.prophetx_client import ProphetXClient
                st.client = ProphetXClient()
            return mod.price_sgps(targets, periods=("FG",), client=st.client)
        if book == "novig":
            from mlb_sgp import novig as mod
            if st.client is None:
                from mlb_sgp.novig_client import NovigClient
                st.client = NovigClient()
            return mod.price_sgps(targets, periods=("FG",), client=st.client)
        raise ValueError(f"unknown book {book!r}")

    @staticmethod
    def _dk_fetchers(st: _BookState) -> dict:
        import scraper_draftkings_sgp as legacy
        ev, struct = st.caches["events"], st.caches["structure"]
        return {
            "fetch_dk_events": lambda session: ev.get_or_fetch(
                "dk_events", lambda: legacy.fetch_dk_events(session)),
            "fetch_main_market_nums": lambda session, eid: struct.get_or_fetch(
                ("nums", eid), lambda: legacy.fetch_main_market_nums(session, eid)),
            "fetch_selection_ids": lambda session, eid, nums, verbose=False:
                struct.get_or_fetch(
                    ("selids", eid),
                    lambda: legacy.fetch_selection_ids(session, eid, nums, verbose)),
        }

    @staticmethod
    def _fd_fetchers(st: _BookState) -> dict:
        import scraper_fanduel_sgp as legacy
        ev, struct = st.caches["events"], st.caches["structure"]
        return {
            "fetch_fd_events": lambda session: ev.get_or_fetch(
                "fd_events", lambda: legacy.fetch_fd_events(session)),
            "fetch_event_runners": lambda session, eid, h, a: struct.get_or_fetch(
                ("runners", eid),
                lambda: legacy.fetch_event_runners(session, eid, h, a)),
        }
