"""Per-book SGP fetch-health history (issue #38).

WHY THIS EXISTS
    Nothing in this repo remembered whether a book answered. #33 made a dead
    book raise instead of returning ``[]``; #35 made a silently-broken parser
    countable. Both are LOG-only — and logs are truncated, rotated, and
    unqueryable. The taker's feed died on 2026-06-28 and went unnoticed for
    four weeks because no one could ask "when did DraftKings last return a
    price?". This module answers that question in SQL.

WHAT IT WRITES
    One row per book per LOGICAL fetch into ``sgp_fetch_health``, in the same
    market DB (and therefore under the same write lock) as ``mlb_sgp_odds``:

        fetched_at TIMESTAMPTZ, book, path, outcome, rows_or_prices,
        duration_sec, error_class, counters_json

    "Logical fetch" means the whole operation the caller asked for, including
    every retry #34 made inside it. A fetch that retried twice is ONE row.

SAFETY INVARIANT
    A health write must never reach pricing. ``record()`` only appends to an
    in-memory buffer; ``flush()`` swallows every exception and RETAINS the
    buffer so the next flush retries. A locked DB costs us history, never a
    quote. This mirrors ``kalshi_mlb_rfq/research.py`` — but the state lives
    on the instance, not in module globals, so two services (or two tests)
    can never contaminate each other.

SIDE EFFECTS
    ``flush()`` connects read-write to the configured DB path, CREATEs the
    table if absent, INSERTs the buffered rows, and (at most once per
    ``prune_interval_sec``) DELETEs rows older than ``retention_days``.
    ``record()`` touches no disk. A recorder built with ``db_path=None`` is
    inert end to end and never creates a file.
"""
from __future__ import annotations

import json
import logging
import threading
import time
from datetime import datetime, timezone

import duckdb

from mlb_sgp._shared import COUNTER_NAMES, FetchCountersSnapshot

log = logging.getLogger(__name__)

# 'error' is NOT in the ticket's four-value vocabulary. It is here because a
# TypeError in our own parse code is not a dead endpoint: folding it into
# 'transport_error' would point a fixer at the network, which is exactly the
# misdiagnosis #33 exists to remove. Anything outside {ok, empty} is a
# failure for streak purposes (#37).
OUTCOMES = frozenset({"ok", "empty", "transport_error", "timeout", "error"})
PATHS = frozenset({"sweep", "on_demand", "warming"})
# Issue #50: on-demand rows tag whether the call paid a wire structure fetch
# (cold) or was fully cache-served (warm). Sweep/warming rows carry NULL —
# the split only means something on the quote path.
CACHE_STATES = frozenset({"cold", "warm"})

DEFAULT_BUFFER_MAX = 5_000       # ~1 hour of maker traffic; drops oldest
DEFAULT_RETENTION_DAYS = 30
DEFAULT_PRUNE_INTERVAL_SEC = 3600.0

# Deliberately no index. Measured through this writer: 1M rows (~30 days at
# post-#53 quote volume) is ~21 MB and scans in milliseconds — an ART index
# would only add write cost on the path we keep off the quote budget.
SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS sgp_fetch_health (
    fetched_at     TIMESTAMPTZ,
    book           VARCHAR,
    path           VARCHAR,
    outcome        VARCHAR,
    rows_or_prices INTEGER,
    duration_sec   DOUBLE,
    error_class    VARCHAR,
    counters_json  VARCHAR,
    cache_state    VARCHAR
);
"""

# Pre-#50 tables lack cache_state; ALTER is a no-op once applied. Runs on
# every flush right after SCHEMA_SQL — same idempotent-DDL pattern, so a
# live DB migrates on the first post-upgrade flush with no manual step.
_MIGRATE_SQL = """
ALTER TABLE sgp_fetch_health ADD COLUMN IF NOT EXISTS cache_state VARCHAR;
"""

_INSERT_SQL = """
INSERT INTO sgp_fetch_health
    (fetched_at, book, path, outcome, rows_or_prices, duration_sec,
     error_class, counters_json, cache_state)
VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
"""


def counters_to_json(snapshot: FetchCountersSnapshot | None) -> str | None:
    """Serialize a #35 counter snapshot for the ``counters_json`` column.

    Every counter is stored even when zero: a reader that has to COALESCE a
    missing key will eventually read "absent" as 0 and call a dead parser
    healthy. Storage is not the constraint — see the measurement above.
    """
    if snapshot is None:
        return None
    payload = {name: getattr(snapshot, name) for name in COUNTER_NAMES}
    if snapshot.by_stage:
        payload["by_stage"] = dict(snapshot.by_stage)
    return json.dumps(payload)


def transport_error_class(exc) -> str:
    """A stable, groupable label for a ``BookTransportError``.

    ``BookTransportError:events:403`` says "the event listing 403'd" four
    weeks later without cross-referencing bot.log; the stage alone would not
    distinguish DK's Akamai 403 from Novig's DNS failure.
    """
    parts = [type(exc).__name__, getattr(exc, "stage", "unknown")]
    status = getattr(exc, "status_code", None)
    if status is not None:
        parts.append(str(status))
    return ":".join(parts)


class OnDemandCoverageCounter:
    """Per-book outcome tally over the on_demand path (#81).

    The maker's sweep — and with it `scrape_done`'s periodic book counts —
    is gone; this counter is the research-side heir: SGPService owns one
    UNCONDITIONALLY (never the alerter, which `build_alerter` nulls out
    when alerting is disabled) and feeds it from the same call site as
    every health row, and the bot's tick drains it into an
    `on_demand_coverage` research event. Only on_demand outcomes count —
    that is the path that prices quotes; warming (and the taker's sweep)
    say nothing about coverage.

    `observe()` is O(1), lock-protected and never raises — it runs on
    arbitrary RFQ threads. `snapshot_and_reset()` returns a deep copy.
    """

    def __init__(self):
        self._lock = threading.Lock()
        self._counts: dict[str, dict[str, int]] = {}

    def observe(self, *, book: str, path: str, outcome: str) -> None:
        if path != "on_demand":
            return
        try:
            with self._lock:
                by_outcome = self._counts.setdefault(book, {})
                by_outcome[outcome] = by_outcome.get(outcome, 0) + 1
        except Exception:   # pragma: no cover — must never break a fetch
            pass

    def snapshot_and_reset(self) -> dict[str, dict[str, int]]:
        with self._lock:
            snap = {book: dict(by_outcome)
                    for book, by_outcome in self._counts.items()}
            self._counts.clear()
        return snap


class FetchHealthRecorder:
    """Buffered, failure-proof writer for ``sgp_fetch_health``.

    Parameters
    ----------
    db_path
        The market DB the fetch itself writes to. ``None`` disables the
        recorder entirely (dashboard CLI path, and every test that doesn't
        ask for health) — no buffer, no file.
    """

    def __init__(self, db_path: str | None, *,
                 buffer_max: int = DEFAULT_BUFFER_MAX,
                 retention_days: int = DEFAULT_RETENTION_DAYS,
                 prune_interval_sec: float = DEFAULT_PRUNE_INTERVAL_SEC,
                 now_fn=time.monotonic):
        self.db_path = str(db_path) if db_path else None
        self.buffer_max = buffer_max
        self.retention_days = retention_days
        self.prune_interval_sec = prune_interval_sec
        self._now = now_fn
        self._lock = threading.Lock()
        self._buffer: list[tuple] = []
        self._last_prune: float | None = None
        self._last_warn = 0.0

    @property
    def enabled(self) -> bool:
        return self.db_path is not None

    def pending(self) -> int:
        with self._lock:
            return len(self._buffer)

    # -------------------------------------------------------------- #
    # Write path                                                      #
    # -------------------------------------------------------------- #

    def record(self, *, book: str, path: str, outcome: str,
               rows_or_prices: int | None = None,
               duration_sec: float | None = None,
               error_class: str | None = None,
               counters: FetchCountersSnapshot | None = None,
               cache_state: str | None = None) -> None:
        """Buffer one fetch's health. O(1), touches no disk, NEVER raises.

        Called from sweep worker threads, the on-demand caller thread, and
        the refresh thread — hence the lock.

        ``cache_state`` (#50): 'cold' | 'warm' on on-demand rows, None
        elsewhere — cold/warm latency percentiles are one GROUP BY away.
        """
        if not self.enabled:
            return
        # Warn, don't raise: a vocabulary typo silently poisons #37's streak
        # queries, but it must not be able to break a fetch.
        if outcome not in OUTCOMES or path not in PATHS:
            log.warning("sgp_fetch_health: unknown outcome=%r / path=%r for "
                        "book=%s", outcome, path, book)
        if cache_state is not None and cache_state not in CACHE_STATES:
            log.warning("sgp_fetch_health: unknown cache_state=%r for "
                        "book=%s", cache_state, book)
        try:
            row = (datetime.now(timezone.utc), book, path, outcome,
                   None if rows_or_prices is None else int(rows_or_prices),
                   None if duration_sec is None else float(duration_sec),
                   error_class, counters_to_json(counters), cache_state)
            with self._lock:
                self._buffer.append(row)
                if len(self._buffer) > self.buffer_max:
                    del self._buffer[:len(self._buffer) - self.buffer_max]
        except Exception:   # pragma: no cover — record must never break pricing
            pass

    def flush(self) -> int:
        """Write the buffer in one batch. Returns rows written (0 on failure).

        Never raises. On failure the buffer is RETAINED (capped by record) so
        the next flush retries — a locked market DB during a sweep upsert is
        expected, not exceptional.
        """
        if not self.enabled:
            return 0
        with self._lock:
            batch = list(self._buffer)
        if not batch:
            return 0
        try:
            # retries=1: a locked DB must not stall a quote tick. The maker
            # flushes ~4x/sec; the next one picks these rows straight up.
            con = _connect(self.db_path, retries=1)
            try:
                con.execute(SCHEMA_SQL)
                con.execute(_MIGRATE_SQL)
                con.executemany(_INSERT_SQL, batch)
                if self._prune_due():
                    self._prune(con)
                    self._last_prune = self._now()
            finally:
                con.close()
        except Exception as exc:
            self._warn_rate_limited(exc)
            return 0
        with self._lock:
            del self._buffer[:len(batch)]   # only drop what actually landed
        return len(batch)

    # -------------------------------------------------------------- #
    # Retention                                                       #
    # -------------------------------------------------------------- #

    def _prune_due(self) -> bool:
        if self._last_prune is None:
            return True
        return (self._now() - self._last_prune) >= self.prune_interval_sec

    def _prune(self, con) -> None:
        con.execute(
            "DELETE FROM sgp_fetch_health "
            "WHERE fetched_at < now() - (INTERVAL 1 DAY * ?)",
            [self.retention_days])

    def _warn_rate_limited(self, exc: BaseException) -> None:
        now = time.time()
        if now - self._last_warn >= 60.0:
            log.warning("sgp_fetch_health flush failed (%d rows retained): %s",
                        self.pending(), exc)
            self._last_warn = now


def _connect(db_path: str, retries: int = 1):
    """Open the market DB read-write, retrying only on a lock conflict."""
    last_err = None
    for attempt in range(max(1, retries)):
        try:
            return duckdb.connect(db_path)
        except duckdb.IOException as e:
            last_err = e
            if attempt < retries - 1:
                time.sleep(0.05 * (2 ** attempt))
    raise last_err
