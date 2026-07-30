"""Run-time alerting when an SGP book stops answering (issue #37).

WHY THIS EXISTS
    #38 gave the bots a MEMORY of every fetch (``sgp_fetch_health``). This
    module is the alarm bell on top of it. The failure it prevents already
    happened: the taker's feed died on 2026-06-28 and went unnoticed for
    four weeks, because a dead book produced no sound at all.

    Deliberately NOT a daily coverage-audit cron (user decision, do not
    relitigate): a cron is too slow for a book that dies mid-slate, and too
    noisy while the bots are intentionally down — which they are, for weeks
    at a time. Alerting is run-time, inside the bot loop.

WHAT IT KEYS ON — STREAKS, NOT DATA AGE
    An age-based rule ("no rows for 10 minutes") would false-fire by design
    once #57 slows the sweep to background structure-warming. Consecutive
    failed fetches are invariant to cadence.

THE OUTCOME VOCABULARY IS FIVE VALUES, NOT FOUR
    ok / empty / transport_error / timeout / error. ``empty`` means the book
    ANSWERED and had no markets: an off-day, a thin slate, or (on-demand)
    "this book won't price that combo" — which is ``_Verdict``'s DEFAULT for
    a plain ``None``. Treating it as a failure pages on every quiet morning
    and on every RFQ for a game a book doesn't list. The predicate is
    ``outcome NOT IN ('ok','empty')``, identical to the shipped
    ``current_failure_streak`` query in ``fetch_health_queries.sql``.

ABSENCE IS NOT FAILURE
    A book skipped by ``min_refresh_sec`` writes no health row and makes no
    observation here. A gap means "we didn't ask", and state simply persists.

WHY IN-MEMORY AND NOT A QUERY OVER sgp_fetch_health
    Two reasons. (1) The table is a LOSSY view of the signal: ``flush()``
    swallows failures and retains its buffer, so rows land up to
    ``HEALTH_FLUSH_MIN_INTERVAL_SEC`` late and, under a stuck lock, later
    still — an alerter reading the table would under-count exactly when the
    bot is under stress. (2) Reproduced on duckdb 1.4.4: while ANY
    read-write connection to a file is open in this process, an in-process
    READ-ONLY connect to the same file raises ``ConnectionException``
    immediately — not an IOException, not retryable. The maker holds rw
    handles on the market DB to flush health and upsert odds, and
    ``kalshi_mlb_mm/main.py``'s narrow ``except (IOException,
    CatalogException)`` guards would not catch it: it would fall to their
    outer fail-safe, return None, and COST A QUOTE. An alerting reader on a
    timer thread is precisely what would open that window. The monitor
    dashboard reads the table instead — it is a separate PROCESS, so it hits
    the ordinary cross-process file lock, which ``kalshi_mlb_monitor`` already
    retries.

SIDE EFFECTS
    ``observe()``: none — an O(1) counter update under a lock. Safe on the
    quote path and on sweep worker threads.
    ``dispatch()``: logs, and spawns a short-lived daemon thread per alert to
    run the notifier (``osascript``, up to 10s). Call it from the bot's tick
    loop. It touches no database, so the hand-off is safe.
"""
from __future__ import annotations

import logging
import subprocess
import threading
from dataclasses import dataclass

log = logging.getLogger(__name__)

# The one predicate. Kept in lockstep with `outcome IN ('ok','empty')` in
# fetch_health_queries.sql::current_failure_streak by a test.
HEALTHY_OUTCOMES = frozenset({"ok", "empty"})

DEFAULT_PATHS = ("sweep", "on_demand")
DEFAULT_STREAK_THRESHOLD = 3
DEFAULT_MIN_HEALTHY_BOOKS = 2

NOTIFY_TIMEOUT_SEC = 10.0


def is_failure(outcome: str | None) -> bool:
    """True when a fetch outcome means the book did not answer."""
    return outcome not in HEALTHY_OUTCOMES


def notify_macos(title: str, message: str) -> None:
    """Fire a macOS notification. Never raises.

    Same mechanism as ``coverage_audit/audit.py::notify``, copied rather
    than imported: a trading bot importing the coverage-audit package to
    send a toast is a worse coupling than eight duplicated lines, and #37
    explicitly keeps SGP books out of that audit's registry.
    """
    try:
        subprocess.run(
            ["osascript", "-e",
             f'display notification "{_osa(message)}" with title "{_osa(title)}"'],
            check=False, capture_output=True, timeout=NOTIFY_TIMEOUT_SEC)
    except Exception:
        pass


def _osa(text: str) -> str:
    """Escape for embedding in an AppleScript string literal.

    ``error_class`` and book names are ours, but a stray quote would produce
    a silently-failing osascript — i.e. exactly the silent failure this
    ticket exists to remove.
    """
    return str(text).replace("\\", "\\\\").replace('"', '\\"')


@dataclass(frozen=True)
class BookAlert:
    """One emitted alert. Returned by ``dispatch()`` for tests and callers."""
    kind: str               # book_degraded | book_recovered
                            # consensus_dark | consensus_restored
    message: str
    notify: bool
    book: str | None = None
    path: str | None = None
    streak: int = 0
    error_class: str | None = None


class BookHealthAlerter:
    """Consecutive-failure tracker with per-incident alert suppression.

    Parameters
    ----------
    label
        Which bot is speaking ("MM" / "RFQ") — it lands in every message, so
        a toast is unambiguous when both bots run.
    streak_threshold
        Rule A: consecutive failures on one (book, path) before alerting.
        One 403 is noise; three in a row is a dead book.
    min_healthy_books
        Rule B: below this many healthy books the bot cannot form consensus
        and is quoting blind. Pass the bot's REAL gate (maker
        ``MIN_AGREEING_BOOKS``, taker ``MIN_BOOK_COUNT_FOR_BLEND``) rather
        than a second knob that can silently disagree with it.
    paths
        Which fetch paths count. Post-#53 every quote is priced by an
        on-demand fetch and the sweep is background warming, so #57 flips
        this to ``("on_demand",)`` in CONFIG, not code.
    notifier
        ``callable(title, message)``. Injected for tests.
    notify_async
        Run the notifier on a daemon thread (default). The maker ticks 4x/s
        and ``osascript`` can block for seconds.
    """

    def __init__(self, *, label: str,
                 streak_threshold: int = DEFAULT_STREAK_THRESHOLD,
                 min_healthy_books: int = DEFAULT_MIN_HEALTHY_BOOKS,
                 paths: tuple[str, ...] = DEFAULT_PATHS,
                 notifier=notify_macos,
                 notify_async: bool = True,
                 enabled: bool = True):
        self.label = label
        self.streak_threshold = max(1, int(streak_threshold))
        self.min_healthy_books = int(min_healthy_books)
        self.paths = tuple(paths)
        self.notifier = notifier
        self.notify_async = notify_async
        self.enabled = enabled
        self._lock = threading.Lock()
        self._streaks: dict[tuple[str, str], int] = {}
        self._degraded: set[tuple[str, str]] = set()
        self._pending: list[BookAlert] = []
        self._dark = False

    # ------------------------------------------------------------------ #
    # Observation — hot path, called per fetch                            #
    # ------------------------------------------------------------------ #

    def observe(self, *, book: str, path: str, outcome: str,
                error_class: str | None = None) -> None:
        """Record one LOGICAL fetch's outcome. O(1), no I/O, never raises.

        Called from the sweep's result loop (main thread) and from
        ``price_on_demand`` on arbitrary RFQ threads — hence the lock. One
        call per logical fetch, never per wire attempt: #34's
        ``request_with_retry`` has already absorbed its retries by the time
        an outcome exists.
        """
        if not self.enabled or path not in self.paths:
            return
        try:
            key = (book, path)
            with self._lock:
                if is_failure(outcome):
                    streak = self._streaks.get(key, 0) + 1
                    self._streaks[key] = streak
                    if (streak >= self.streak_threshold
                            and key not in self._degraded):
                        self._degraded.add(key)
                        self._pending.append(self._degraded_alert(
                            book, path, streak, error_class))
                    return
                self._streaks[key] = 0
                if key in self._degraded:
                    self._degraded.discard(key)
                    self._pending.append(BookAlert(
                        kind="book_recovered", book=book, path=path,
                        notify=False,
                        message=(f"[{self.label}] SGP book recovered: "
                                 f"{book}/{path} answered again")))
        except Exception:   # pragma: no cover — must never break a fetch
            pass

    # ------------------------------------------------------------------ #
    # Dispatch — tick thread only                                         #
    # ------------------------------------------------------------------ #

    def dispatch(self) -> list[BookAlert]:
        """Emit everything queued since the last call. Never raises.

        Rule B is evaluated HERE rather than in ``observe`` on purpose:
        ``refresh()`` records its books one at a time, so a mid-loop reading
        of "how many books are healthy" can catch a transient. Once per tick,
        on a settled state, cannot.
        """
        if not self.enabled:
            return []
        try:
            with self._lock:
                alerts = self._pending
                self._pending = []
                capacity = self._capacity_alert_locked()
            if capacity is not None:
                alerts.append(capacity)
            for alert in alerts:
                self._emit(alert)
            return alerts
        except Exception:   # pragma: no cover — alerting must never halt a bot
            log.exception("book-health dispatch failed")
            return []

    def snapshot(self) -> dict:
        """Current view, for tests and debugging.

        ``streaks`` is nested book -> path -> int; ``healthy_books`` lists
        observed books with no degraded path; ``dark`` is Rule B's state.
        """
        with self._lock:
            streaks: dict[str, dict[str, int]] = {}
            for (book, path), n in self._streaks.items():
                streaks.setdefault(book, {})[path] = n
            return {
                "streaks": streaks,
                "degraded": sorted(f"{b}/{p}" for b, p in self._degraded),
                "healthy_books": self._healthy_books_locked(),
                "known_books": sorted({b for b, _ in self._streaks}),
                "dark": self._dark,
            }

    # ------------------------------------------------------------------ #
    # Internals                                                           #
    # ------------------------------------------------------------------ #

    def _degraded_alert(self, book, path, streak, error_class) -> BookAlert:
        detail = f" (last: {error_class})" if error_class else ""
        return BookAlert(
            kind="book_degraded", book=book, path=path, streak=streak,
            error_class=error_class, notify=True,
            message=(f"[{self.label}] SGP book degraded: {book}/{path} — "
                     f"{streak} consecutive failed fetches{detail}"))

    def _healthy_books_locked(self) -> list[str]:
        """Books observed at least once with no degraded considered-path.

        A book we have never fetched is neither healthy nor unhealthy — it is
        unknown, and must not be counted either way (that is what keeps a
        never-due book from reading as a capacity loss).
        """
        known = {b for b, _ in self._streaks}
        unhealthy = {b for b, _ in self._degraded}
        return sorted(known - unhealthy)

    def _capacity_alert_locked(self) -> BookAlert | None:
        """Rule B: are we below the consensus floor? Called under the lock.

        The ``len(known) >= min_healthy_books`` precondition is what stops a
        false alarm on the process's FIRST observation, where exactly one
        book is known and "1 healthy < 2" would otherwise be true.
        """
        known = {b for b, _ in self._streaks}
        healthy = self._healthy_books_locked()
        if len(known) < self.min_healthy_books:
            return None
        if len(healthy) < self.min_healthy_books:
            if self._dark:
                return None
            self._dark = True
            degraded = ", ".join(sorted(f"{b}/{p}" for b, p in self._degraded))
            return BookAlert(
                kind="consensus_dark", notify=True,
                message=(f"[{self.label}] SGP consensus DARK: "
                         f"{len(healthy)} healthy book(s), need "
                         f"{self.min_healthy_books} — quoting effectively "
                         f"blind. Degraded: {degraded}"))
        if self._dark:
            self._dark = False
            return BookAlert(
                kind="consensus_restored", notify=True,
                message=(f"[{self.label}] SGP consensus restored: "
                         f"{len(healthy)} healthy books "
                         f"({', '.join(healthy)})"))
        return None

    def _emit(self, alert: BookAlert) -> None:
        if alert.notify:
            log.error("%s", alert.message)
        else:
            # Recovery is log-only: a "DK is back" toast is noise you cannot
            # act on, and Rule B's restore already covers the case that
            # changes whether the bot can quote.
            log.warning("%s", alert.message)
        if not alert.notify or self.notifier is None:
            return
        title = f"Kalshi {self.label}: SGP book health"
        if self.notify_async:
            threading.Thread(target=self._notify_safe,
                             args=(title, alert.message),
                             name="book-health-notify", daemon=True).start()
        else:
            self._notify_safe(title, alert.message)

    def _notify_safe(self, title: str, message: str) -> None:
        try:
            self.notifier(title, message)
        except Exception as exc:
            log.warning("book-health notifier failed: %s", exc)


def build_alerter(*, label: str, enabled: bool, streak_threshold: int,
                  min_healthy_books: int,
                  paths: tuple[str, ...]) -> BookHealthAlerter | None:
    """Config -> alerter, or None when disabled.

    Both bots build theirs this way so the "disabled means no object at all"
    convention (and therefore a fully inert ``SGPService``) lives in one place.
    """
    if not enabled:
        return None
    return BookHealthAlerter(label=label, streak_threshold=streak_threshold,
                             min_healthy_books=min_healthy_books, paths=paths)
