# DuckDB Write-Race Fix Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate the chronic DuckDB single-writer lock race between the dashboard's in-process trades poll, its scrape subprocesses, and any external writer.

**Architecture:** Keep the dashboard self-scheduling (user preference). Three layered mitigations: (1) batch `fetch_trades()` writes into one connection per cycle, cutting the in-process lock-acquisition rate ~100×. (2) Add exponential-backoff retry to `write_connection()` so any remaining race self-heals. (3) Reap scrape subprocesses via `subprocess.run` in daemon threads so zombie children stop blocking DuckDB's live-PID check. Plus an audit step to flag any scraper that holds the lock across network I/O.

**Tech Stack:** Python 3.14, DuckDB, Dash, pytest.

**Spec:** `docs/superpowers/specs/2026-04-23-duckdb-write-race-design.md`

---

## File Structure

Files created/modified and their responsibilities:

- **Modify** `nfl_draft/lib/db.py` — relocate `_is_lock_error` helper here from `queries.py`; add module-level backoff constants; add retry loop inside `write_connection()`.
- **Modify** `nfl_draft/lib/queries.py` — import `_is_lock_error` from `db.py` instead of defining it locally. No behavioral change.
- **Modify** `nfl_draft/scrapers/kalshi.py` — refactor `fetch_trades()` so all `INSERT OR IGNORE INTO kalshi_trades` and `kalshi_poll_state` upserts execute on a single `write_connection()` opened at the top of the function.
- **Modify** `kalshi_draft/app.py` — replace fire-and-forget `subprocess.Popen` calls with `subprocess.run` inside daemon threads, in both `_run_scrape_once` and the `__main__` startup block.
- **Modify** `nfl_draft/README.md` — add a "Concurrency model" subsection describing the single-writer constraint, retry layer, and trades-batch optimization.
- **Create** `nfl_draft/tests/unit/test_db_retry.py` — retry-layer unit tests.
- **Create** `nfl_draft/tests/unit/test_kalshi_trades_batching.py` — single-connection invariant for `fetch_trades()`.

No other files. Other scraper call sites are audited in Task 5 and only modified if the audit finds lock-across-IO issues.

---

### Task 1: Create worktree and feature branch

**Files:**
- None (infra only)

- [ ] **Step 1: Verify you are on `main` and the tree is clean**

Run: `git -C /Users/callancapitolo/NFLWork status --short && git -C /Users/callancapitolo/NFLWork branch --show-current`

Expected: only the pre-existing untracked entries (`.playwright-mcp/`, `docs/superpowers/plans/2026-04-22-bookmaker-no-popup.md`, `docs/superpowers/plans/2026-04-23-nfl-draft-staleness-filter.md`), and `main` for the branch.

- [ ] **Step 2: Create the worktree**

Run:
```bash
git -C /Users/callancapitolo/NFLWork worktree add \
  /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race \
  -b fix/duckdb-write-race
```

Expected: `Preparing worktree ...` and `HEAD is now at ...` messages.

- [ ] **Step 3: Verify the worktree**

Run: `git -C /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race branch --show-current`

Expected: `fix/duckdb-write-race`.

**From Task 2 onward, run every command from inside `/Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race`.**

---

### Task 2: Relocate `_is_lock_error` from queries.py to db.py

**Files:**
- Modify: `nfl_draft/lib/db.py` — append helper.
- Modify: `nfl_draft/lib/queries.py` — delete local definition; import from `db.py`.

Rationale: one helper, two call sites. The retry layer in Task 3 needs the same signature detector that `_safe_read` uses. Relocating first keeps Task 3's diff focused on behavior change.

- [ ] **Step 1: Append helper + constants to `nfl_draft/lib/db.py` above the `write_connection` function**

In `nfl_draft/lib/db.py`, after the existing imports (line 7) and `DB_PATH` definition (line 9), insert:

```python
# -- Lock-contention retry layer ----------------------------------------------
# DuckDB allows exactly one writer across processes. The dashboard runs an
# in-process trades poll plus scrape subprocesses; any concurrent writer
# (manual CLI, cron) can race their write windows. `write_connection()`
# retries with exponential backoff when it hits a lock error, so transient
# collisions self-heal instead of crashing the caller.

_INITIAL_BACKOFF_SEC = 0.1
_MAX_BACKOFF_SEC = 5.0
_DEFAULT_RETRY_TIMEOUT_SEC = 30.0


def _is_lock_error(err: BaseException) -> bool:
    """True if this exception matches the DuckDB lock-contention signature.

    DuckDB surfaces lock errors as IOException with messages like
    'Could not set lock on file' / 'Conflicting lock is held'. Also
    includes 'file handle conflict' for the in-process collision mode
    that happens when two readers race to open the same file.
    """
    msg = str(err).lower()
    return any(marker in msg for marker in (
        "lock on file",
        "conflicting lock",
        "could not set lock",
        "i/o error",
        "io error",
        "file handle conflict",
    ))
```

- [ ] **Step 2: Replace the private definition in `nfl_draft/lib/queries.py` with an import**

In `nfl_draft/lib/queries.py`, after the existing imports (lines 27-33), add this import near the other `nfl_draft.lib.db` import:

Find (around line 33):
```python
from nfl_draft.lib.db import read_connection
```

Replace with:
```python
from nfl_draft.lib.db import read_connection, _is_lock_error
```

Then delete the local `def _is_lock_error(...)` block and its preceding docstring-style comment (roughly lines 56-72, the block that starts with `def _is_lock_error` and ends with the closing `))` of the return). The surrounding `QueryLocked`, `_LOCKED`, and `_safe_read` stay unchanged — they call `_is_lock_error` by name, which now resolves to the imported symbol.

- [ ] **Step 3: Run the full nfl_draft test suite to confirm no behavioral regression**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest nfl_draft/tests/ -q
```

Expected: same pass/fail count as `main` (154 passed, 1 flaky `test_writer_and_reader_coexist`, 1 skipped). No new failures.

- [ ] **Step 4: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
git add nfl_draft/lib/db.py nfl_draft/lib/queries.py && \
git commit -m "$(cat <<'EOF'
refactor(nfl_draft): share _is_lock_error between db.py and queries.py

Moves the lock-signature detector from queries.py (private helper used by
_safe_read) into db.py (public helper) so the upcoming write_connection
retry layer can reuse it. No behavior change — queries.py now imports the
symbol from db.py.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

Expected: one new commit on `fix/duckdb-write-race`.

---

### Task 3: Add retry-with-backoff to `write_connection()`

**Files:**
- Create: `nfl_draft/tests/unit/test_db_retry.py`
- Modify: `nfl_draft/lib/db.py` (the `write_connection` contextmanager)

TDD order: write the three failing tests first, confirm they fail, then implement.

- [ ] **Step 1: Create `nfl_draft/tests/unit/test_db_retry.py`**

Create the file with:

```python
"""Unit tests for write_connection's retry-on-lock-conflict layer.

These tests monkeypatch duckdb.connect so we can inject synthetic
lock-signature errors without needing a second writer process.
"""
import time

import duckdb
import pytest

from nfl_draft.lib import db as db_module


class _FakeConn:
    """Stand-in returned by a successful duckdb.connect() mock."""
    def close(self):
        pass


def _lock_error():
    return duckdb.IOException(
        'IO Error: Could not set lock on file "fake.duckdb": '
        'Conflicting lock is held in Python (PID 9999)'
    )


def test_write_connection_retries_on_lock_error(monkeypatch):
    """If duckdb.connect raises a lock error twice then succeeds, the retry
    loop absorbs the first two attempts and yields the third connection.
    The caller never sees the transient error."""
    call_count = {"n": 0}
    conn = _FakeConn()

    def _fake_connect(path):
        call_count["n"] += 1
        if call_count["n"] < 3:
            raise _lock_error()
        return conn

    monkeypatch.setattr(duckdb, "connect", _fake_connect)
    # Speed up the test: backoff sleep can be very small in tests.
    monkeypatch.setattr(db_module, "_INITIAL_BACKOFF_SEC", 0.001)
    monkeypatch.setattr(db_module, "_MAX_BACKOFF_SEC", 0.001)

    with db_module.write_connection(timeout_sec=5.0) as yielded:
        assert yielded is conn

    assert call_count["n"] == 3


def test_write_connection_gives_up_at_timeout(monkeypatch):
    """If the lock never clears, write_connection re-raises after timeout_sec."""
    def _always_lock(path):
        raise _lock_error()

    monkeypatch.setattr(duckdb, "connect", _always_lock)
    monkeypatch.setattr(db_module, "_INITIAL_BACKOFF_SEC", 0.01)
    monkeypatch.setattr(db_module, "_MAX_BACKOFF_SEC", 0.01)

    started = time.monotonic()
    with pytest.raises(duckdb.IOException):
        with db_module.write_connection(timeout_sec=0.2):
            pytest.fail("Should never enter the with-block body")
    elapsed = time.monotonic() - started

    # Ran for at least the timeout, but not absurdly longer.
    assert 0.15 <= elapsed <= 2.0


def test_write_connection_propagates_non_lock_errors_immediately(monkeypatch):
    """Non-lock errors bypass the retry loop and propagate on first call."""
    call_count = {"n": 0}

    def _other_error(path):
        call_count["n"] += 1
        raise duckdb.Error("schema mismatch: expected column X")

    monkeypatch.setattr(duckdb, "connect", _other_error)
    monkeypatch.setattr(db_module, "_INITIAL_BACKOFF_SEC", 0.001)

    started = time.monotonic()
    with pytest.raises(duckdb.Error, match="schema mismatch"):
        with db_module.write_connection(timeout_sec=5.0):
            pytest.fail("Should never enter the with-block body")
    elapsed = time.monotonic() - started

    # Exactly one attempt, essentially instant.
    assert call_count["n"] == 1
    assert elapsed < 0.1
```

- [ ] **Step 2: Run the new tests — confirm they fail against the current `write_connection`**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest \
  nfl_draft/tests/unit/test_db_retry.py -v
```

Expected: all three tests FAIL. Specifically:
- `test_write_connection_retries_on_lock_error` fails because the current `write_connection` has no retry — the first `_lock_error()` propagates immediately.
- `test_write_connection_gives_up_at_timeout` will raise on the first attempt (fine assertion-wise) but the `timeout_sec` keyword doesn't exist — TypeError.
- `test_write_connection_propagates_non_lock_errors_immediately` fails on the `timeout_sec` keyword too.

All three failures are expected; this is the red phase.

- [ ] **Step 3: Replace `write_connection` in `nfl_draft/lib/db.py`**

In `nfl_draft/lib/db.py`, replace the current `write_connection` block (lines 12-19):

```python
@contextmanager
def write_connection() -> Iterator[duckdb.DuckDBPyConnection]:
    """Short-lived write connection. Lock held for milliseconds."""
    con = duckdb.connect(str(DB_PATH))
    try:
        yield con
    finally:
        con.close()
```

Replace with:

```python
@contextmanager
def write_connection(
    timeout_sec: float = _DEFAULT_RETRY_TIMEOUT_SEC,
) -> Iterator[duckdb.DuckDBPyConnection]:
    """Short-lived write connection with retry-on-lock-contention.

    DuckDB allows only one writer across processes. When the dashboard is
    holding the write lock (trades-poll, subprocess scrape, bet-log), a
    concurrent writer from another process (manual CLI, cron) races the
    open. Instead of failing on the first collision, retry with
    exponential backoff for up to ``timeout_sec`` seconds. Non-lock errors
    (schema mismatch, typos, permissions) propagate immediately.
    """
    import time
    deadline = time.monotonic() + timeout_sec
    backoff = _INITIAL_BACKOFF_SEC
    while True:
        try:
            con = duckdb.connect(str(DB_PATH))
            break
        except (duckdb.IOException, duckdb.Error) as e:
            if not _is_lock_error(e) or time.monotonic() >= deadline:
                raise
            time.sleep(backoff)
            backoff = min(backoff * 2, _MAX_BACKOFF_SEC)
    try:
        yield con
    finally:
        con.close()
```

Leave `read_connection` (lines 22-36) unchanged — callbacks keep the fast-fail sentinel semantics.

- [ ] **Step 4: Run the unit tests — confirm all three pass**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest \
  nfl_draft/tests/unit/test_db_retry.py -v
```

Expected: 3 passed, 0 failed.

- [ ] **Step 5: Run the full suite to confirm no regression**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest nfl_draft/tests/ -q
```

Expected: 157 passed (154 from main + 3 new retry tests), 1 flaky `test_writer_and_reader_coexist`, 1 skipped. No new failures besides the known flaky one.

- [ ] **Step 6: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
git add nfl_draft/lib/db.py nfl_draft/tests/unit/test_db_retry.py && \
git commit -m "$(cat <<'EOF'
fix(nfl_draft): retry write_connection() on DuckDB lock contention

DuckDB allows one writer across processes. The dashboard's in-process
trades poll and its subprocess scrape regularly collide with external
writers (manual CLI, cron). Retrying with exponential backoff (100ms ->
5s cap, 30s total budget) lets callers ride through transient contention
instead of crashing. Non-lock errors propagate immediately.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 4: Batch `fetch_trades()` writes into a single connection

**Files:**
- Create: `nfl_draft/tests/unit/test_kalshi_trades_batching.py`
- Modify: `nfl_draft/scrapers/kalshi.py` (the `fetch_trades` function body)

TDD order: write the single-connection invariant test first, confirm it fails, then refactor.

- [ ] **Step 1: Create `nfl_draft/tests/unit/test_kalshi_trades_batching.py`**

Create the file with:

```python
"""Unit test: fetch_trades must open exactly ONE write_connection per call,
regardless of how many tickers or batches it processes.

Before this fix, fetch_trades opened a write_connection per ticker-batch
AND per series poll_state update — dozens to hundreds per 15-second cycle.
That rate made the dashboard a near-permanent DuckDB writer, starving
external writers (CLI, cron).
"""
from datetime import datetime
from unittest.mock import MagicMock

import pytest


class _ConnCounter:
    """Counts how many times write_connection() is entered as a context manager."""
    def __init__(self, fake_conn):
        self.fake_conn = fake_conn
        self.enter_count = 0

    def __call__(self):
        return self

    def __enter__(self):
        self.enter_count += 1
        return self.fake_conn

    def __exit__(self, *args):
        return False


def test_fetch_trades_opens_single_write_connection(monkeypatch):
    """With two series and several tickers, one batch returned per ticker,
    fetch_trades must open exactly one write connection for the whole call."""
    # Import inside the test so monkeypatching resolves against the module
    # scope fetch_trades actually imports from.
    from nfl_draft.scrapers import kalshi as kalshi_mod

    # Stub the legacy discover call and per-series ticker enumeration.
    monkeypatch.setattr(
        kalshi_mod.legacy_fetcher, "discover_draft_series",
        lambda: [{"series_ticker": "S1"}, {"series_ticker": "S2"}],
    )
    monkeypatch.setattr(
        kalshi_mod, "_tickers_for_series",
        lambda rcon, st: ["T1", "T2"] if st == "S1" else ["T3"],
    )

    # read_connection stays a normal no-op context manager over a MagicMock
    # so the poll-state cursor query runs harmlessly.
    class _ReadCtx:
        def __enter__(self):
            m = MagicMock()
            m.execute.return_value.fetchall.return_value = []
            return m
        def __exit__(self, *a): return False
    monkeypatch.setattr(kalshi_mod, "read_connection", lambda: _ReadCtx())

    # Each /markets/trades call returns one fake batch then stops.
    def _fake_public_request(path):
        if "cursor=" in path:
            return None
        return {"trades": [{"fake": True}], "cursor": None}
    monkeypatch.setattr(kalshi_mod, "public_request", _fake_public_request)

    now = datetime.now()
    def _fake_parse(raw):
        if not raw:
            return []
        TradeRow = kalshi_mod.TradeRow
        return [TradeRow(
            trade_id=f"id-{id(raw)}", ticker="T?",
            side="yes", price_cents=50, count=1,
            traded_at=now, fetched_at=now,
        )]
    monkeypatch.setattr(kalshi_mod, "parse_trades_response", _fake_parse)

    # The key assertion: swap write_connection with the counter.
    fake_conn = MagicMock()
    counter = _ConnCounter(fake_conn)
    monkeypatch.setattr(kalshi_mod, "write_connection", counter)

    kalshi_mod.fetch_trades()

    assert counter.enter_count == 1, (
        f"fetch_trades opened write_connection {counter.enter_count} times; "
        "must be exactly 1 regardless of ticker or batch count."
    )


def test_fetch_trades_inserts_all_batches_into_shared_connection(monkeypatch):
    """The shared write connection must receive INSERT OR IGNORE calls for
    every batch plus the poll_state upsert for every series."""
    from nfl_draft.scrapers import kalshi as kalshi_mod

    monkeypatch.setattr(
        kalshi_mod.legacy_fetcher, "discover_draft_series",
        lambda: [{"series_ticker": "S1"}],
    )
    monkeypatch.setattr(
        kalshi_mod, "_tickers_for_series",
        lambda rcon, st: ["T1"],
    )

    class _ReadCtx:
        def __enter__(self):
            m = MagicMock()
            m.execute.return_value.fetchall.return_value = []
            return m
        def __exit__(self, *a): return False
    monkeypatch.setattr(kalshi_mod, "read_connection", lambda: _ReadCtx())

    def _fake_public_request(path):
        if "cursor=" in path:
            return None
        return {"trades": [{"fake": True}], "cursor": None}
    monkeypatch.setattr(kalshi_mod, "public_request", _fake_public_request)

    now = datetime.now()
    TradeRow = __import__(
        "nfl_draft.scrapers.kalshi", fromlist=["TradeRow"]
    ).TradeRow
    monkeypatch.setattr(
        kalshi_mod, "parse_trades_response",
        lambda raw: [TradeRow(
            trade_id="id-1", ticker="T1",
            side="yes", price_cents=50, count=1,
            traded_at=now, fetched_at=now,
        )] if raw else [],
    )

    fake_conn = MagicMock()
    class _CtxWrap:
        def __enter__(self): return fake_conn
        def __exit__(self, *a): return False
    monkeypatch.setattr(kalshi_mod, "write_connection", lambda: _CtxWrap())

    kalshi_mod.fetch_trades()

    exec_calls = [call.args[0] for call in fake_conn.execute.call_args_list]
    assert any("kalshi_trades" in sql for sql in exec_calls), (
        "Expected INSERT OR IGNORE INTO kalshi_trades; got SQL: " + repr(exec_calls)
    )
    assert any("kalshi_poll_state" in sql for sql in exec_calls), (
        "Expected kalshi_poll_state upsert; got SQL: " + repr(exec_calls)
    )
```

- [ ] **Step 2: Run the new tests — confirm they fail against the current `fetch_trades`**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest \
  nfl_draft/tests/unit/test_kalshi_trades_batching.py -v
```

Expected:
- `test_fetch_trades_opens_single_write_connection` FAILS with `counter.enter_count == 4` (or similar >1): current code opens one per batch per ticker (3 tickers × 1 batch = 3), plus one per series poll_state (2 series = 2), total ~5 depending on exact fixture behavior. The failure message will show the actual count.
- `test_fetch_trades_inserts_all_batches_into_shared_connection` FAILS because the current `fetch_trades` creates a new connection per batch — the test's single `fake_conn` never receives the poll_state upsert (that goes to a different, unsubstituted connection instance in real code, but under monkeypatch the SAME fake_conn is reused because our patched `write_connection` always yields it; the actual failure mode may be that SQL strings are present anyway — which is FINE, this assertion passes in the current code). **Accept either: the first test failing is the key red signal.**

If only the first test fails and the second passes (because of the shared `fake_conn`), that's acceptable — the first test is the strict invariant we need.

- [ ] **Step 3: Refactor `fetch_trades()` in `nfl_draft/scrapers/kalshi.py`**

In `nfl_draft/scrapers/kalshi.py`, locate `fetch_trades` (starts at line 528) and replace its body from line 560 (`now_local = datetime.now()`) through line 639 (the `except Exception as e: print(f"  [trades] poll_state update failed...")`) with a single-connection version.

**Before** (current lines 560-640):
```python
    now_local = datetime.now()

    for s in series_list:
        series_ticker = s["series_ticker"]
        last_local = cursor_by_series.get(series_ticker)
        min_ts = _local_to_epoch_seconds(last_local) if last_local else None
        max_traded_at = last_local

        for ticker in tickers_by_series.get(series_ticker, []):
            try:
                cursor = None
                while True:
                    path = f"/markets/trades?ticker={ticker}&limit=100"
                    if min_ts is not None:
                        path += f"&min_ts={min_ts}"
                    if cursor:
                        path += f"&cursor={cursor}"

                    raw = public_request(path)
                    if not raw:
                        break

                    batch = parse_trades_response(raw)
                    if batch:
                        with write_connection() as wcon:
                            for t in batch:
                                wcon.execute(
                                    "INSERT OR IGNORE INTO kalshi_trades "
                                    "VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                                    [
                                        t.trade_id, t.ticker, t.side,
                                        t.price_cents, t.count,
                                        t.count * t.price_cents * 0.01,
                                        t.traded_at, t.fetched_at,
                                    ],
                                )
                        ingested.extend(batch)
                        batch_max = max(t.traded_at for t in batch)
                        if max_traded_at is None or batch_max > max_traded_at:
                            max_traded_at = batch_max

                    cursor = raw.get("cursor")
                    if not cursor or not raw.get("trades"):
                        break
            except Exception as e:
                print(f"  [trades] ticker={ticker} failed: {e}")
                traceback.print_exc()
                continue

        try:
            with write_connection() as wcon:
                if max_traded_at is not None and max_traded_at != last_local:
                    wcon.execute(
                        "INSERT INTO kalshi_poll_state "
                        "(series_ticker, last_traded_at, last_polled_at) "
                        "VALUES (?, ?, ?) "
                        "ON CONFLICT (series_ticker) DO UPDATE SET "
                        "last_traded_at=excluded.last_traded_at, "
                        "last_polled_at=excluded.last_polled_at",
                        [series_ticker, max_traded_at, now_local],
                    )
                else:
                    wcon.execute(
                        "INSERT INTO kalshi_poll_state "
                        "(series_ticker, last_traded_at, last_polled_at) "
                        "VALUES (?, ?, ?) "
                        "ON CONFLICT (series_ticker) DO UPDATE SET "
                        "last_polled_at=excluded.last_polled_at",
                        [series_ticker, last_local, now_local],
                    )
        except Exception as e:
            print(f"  [trades] poll_state update failed for {series_ticker}: {e}")

    return ingested
```

**After** — open one connection at the top; share it across all loops; close once at exit:

```python
    now_local = datetime.now()

    # One write connection for the entire fetch. The retry layer in
    # write_connection() absorbs any lock-contention at open time; all
    # inserts/upserts below share this single lock acquisition. This drops
    # the dashboard's trades-poll lock-acquisition rate from
    # O(tickers * batches) per 15-sec cycle to exactly 1.
    with write_connection() as wcon:
        for s in series_list:
            series_ticker = s["series_ticker"]
            last_local = cursor_by_series.get(series_ticker)
            min_ts = _local_to_epoch_seconds(last_local) if last_local else None
            max_traded_at = last_local

            for ticker in tickers_by_series.get(series_ticker, []):
                try:
                    cursor = None
                    while True:
                        path = f"/markets/trades?ticker={ticker}&limit=100"
                        if min_ts is not None:
                            path += f"&min_ts={min_ts}"
                        if cursor:
                            path += f"&cursor={cursor}"

                        raw = public_request(path)
                        if not raw:
                            break

                        batch = parse_trades_response(raw)
                        if batch:
                            for t in batch:
                                wcon.execute(
                                    "INSERT OR IGNORE INTO kalshi_trades "
                                    "VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                                    [
                                        t.trade_id, t.ticker, t.side,
                                        t.price_cents, t.count,
                                        t.count * t.price_cents * 0.01,
                                        t.traded_at, t.fetched_at,
                                    ],
                                )
                            ingested.extend(batch)
                            batch_max = max(t.traded_at for t in batch)
                            if max_traded_at is None or batch_max > max_traded_at:
                                max_traded_at = batch_max

                        cursor = raw.get("cursor")
                        if not cursor or not raw.get("trades"):
                            break
                except Exception as e:
                    print(f"  [trades] ticker={ticker} failed: {e}")
                    traceback.print_exc()
                    continue

            # poll_state upsert for this series, on the shared connection.
            try:
                if max_traded_at is not None and max_traded_at != last_local:
                    wcon.execute(
                        "INSERT INTO kalshi_poll_state "
                        "(series_ticker, last_traded_at, last_polled_at) "
                        "VALUES (?, ?, ?) "
                        "ON CONFLICT (series_ticker) DO UPDATE SET "
                        "last_traded_at=excluded.last_traded_at, "
                        "last_polled_at=excluded.last_polled_at",
                        [series_ticker, max_traded_at, now_local],
                    )
                else:
                    wcon.execute(
                        "INSERT INTO kalshi_poll_state "
                        "(series_ticker, last_traded_at, last_polled_at) "
                        "VALUES (?, ?, ?) "
                        "ON CONFLICT (series_ticker) DO UPDATE SET "
                        "last_polled_at=excluded.last_polled_at",
                        [series_ticker, last_local, now_local],
                    )
            except Exception as e:
                print(f"  [trades] poll_state update failed for {series_ticker}: {e}")

    return ingested
```

Two substantive differences from the original:
1. The outer `with write_connection() as wcon:` wraps everything.
2. The two inner `with write_connection() as wcon:` blocks are gone — their bodies now run directly against the outer `wcon`.

Key invariant preserved: the per-ticker `try/except` still isolates ticker failures, and the per-series `try/except` still isolates poll-state failures.

- [ ] **Step 4: Run the batching tests — confirm both pass**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest \
  nfl_draft/tests/unit/test_kalshi_trades_batching.py -v
```

Expected: both PASS. Specifically, `enter_count == 1`.

- [ ] **Step 5: Full suite sanity**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest nfl_draft/tests/ -q
```

Expected: 159 passed (+2 from Task 4), 1 flaky, 1 skipped.

- [ ] **Step 6: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
git add nfl_draft/scrapers/kalshi.py \
         nfl_draft/tests/unit/test_kalshi_trades_batching.py && \
git commit -m "$(cat <<'EOF'
perf(nfl_draft): batch fetch_trades() writes into one connection

Before: fetch_trades opened a fresh write_connection per ticker-batch plus
another per series for the poll_state upsert. On a 15-second cadence with
~350 ticker-batches per cycle, that was hundreds of DuckDB lock
acquisitions per minute — effectively making the dashboard a permanent
writer from another process's perspective.

After: one write_connection for the whole fetch. All INSERT OR IGNOREs
and poll_state upserts execute on the shared connection. Lock
acquisitions drop from O(tickers * batches) to exactly 1 per cycle.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 5: Audit every `write_connection()` call site for lock-across-IO

**Files:**
- Read-only audit; writes only if findings require fixes.

Goal: confirm no caller holds the DuckDB write lock across slow I/O (network request, sleep, user-visible operation). The retry layer from Task 3 bounds at 30 sec; a long-held writer could starve everyone else.

Authoritative call-site list (from `grep write_connection` at plan-writing time):
- `nfl_draft/lib/seed.py:16` — `seed.run()` at start of every scrape.
- `nfl_draft/lib/db.py:180` — inside `init_schema()`.
- `nfl_draft/lib/quarantine.py:86` — `write_or_quarantine()` called after every scraper.
- `nfl_draft/scrapers/kalshi.py:382` — `_write_legacy_kalshi_odds`.
- `nfl_draft/scrapers/kalshi.py:587, 616` — inside `fetch_trades()` (these are gone after Task 4; verify).
- `kalshi_draft/app.py:1381` — bet-log submission callback.

Plus the raw `duckdb.connect` in `kalshi_draft/db.py:45` — a legacy shim. Not via `write_connection`, but still takes the write lock. Flag it too.

- [ ] **Step 1: Read each call site and classify**

For each file/line above, read ~30 lines around the `with write_connection() as con:` block. For each site, answer ONE question:

> Does any network call, user input, `time.sleep`, or non-DuckDB I/O happen INSIDE the `with` block?

Record findings in a table. Example shape (fill in as you audit):

| File:line                                         | Classification | Notes |
|---------------------------------------------------|----------------|-------|
| `nfl_draft/lib/seed.py:16`                        | CLEAN          | only CREATE TABLE IF NOT EXISTS + INSERT OR REPLACE |
| `nfl_draft/lib/db.py:180` (init_schema)           | CLEAN          | only schema DDL |
| `nfl_draft/lib/quarantine.py:86`                  | ?              | check whether DNS resolution / API call happens inside |
| `nfl_draft/scrapers/kalshi.py:382` (legacy write) | ?              | |
| `kalshi_draft/app.py:1381` (bet log)              | ?              | |
| `kalshi_draft/db.py:45` (raw connect)             | ?              | check what callers do with this |

Replace `?` with `CLEAN` or `DIRTY` after reading. DIRTY = slow I/O happens inside the `with`.

- [ ] **Step 2: If all sites are CLEAN, commit the audit findings as a docstring comment in the spec**

Add a paragraph to `docs/superpowers/specs/2026-04-23-duckdb-write-race-design.md` under a new `## Audit findings (implementation)` section at the bottom, summarizing the per-site classification. Commit as:

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
git add docs/superpowers/specs/2026-04-23-duckdb-write-race-design.md && \
git commit -m "$(cat <<'EOF'
docs(spec): record audit findings for DuckDB write-race fix

All write_connection() call sites hold the lock only across in-memory or
local-DB operations — no network I/O, sleeps, or user input inside any
`with write_connection()` block. The retry layer's 30-second budget is
sufficient for all currently-observed writer windows.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 3: If any site is DIRTY, STOP and escalate**

If you found lock-across-IO at any site, do NOT silently "fix" it by pulling I/O out of the block — that's a non-trivial refactor that deserves its own plan. Report back with:
- The dirty site(s).
- A description of what I/O happens inside.
- A recommendation on whether to fix here or defer.

Wait for controller guidance before proceeding.

---

### Task 6: Reap scrape subprocesses via `subprocess.run` in daemon threads

**Files:**
- Modify: `kalshi_draft/app.py` (two blocks: `_run_scrape_once` around lines 38-50; the `__main__` startup block at lines 1513-1528)

- [ ] **Step 1: Replace `_run_scrape_once` body**

In `kalshi_draft/app.py`, replace the current `_run_scrape_once` (approx lines 38-50):

```python
def _run_scrape_once() -> None:
    """Fire a detached scrape subprocess. Never raises."""
    try:
        repo_root = Path(__file__).resolve().parent.parent
        subprocess.Popen(
            [sys.executable, "-m", "nfl_draft.run", "--mode", "scrape", "--book", "all"],
            cwd=str(repo_root),
            stdout=open("/tmp/nfl_draft_periodic_scrape.log", "a"),
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
    except Exception as e:
        print(f"[periodic] scrape launch failed: {e}")
```

With:

```python
def _run_scrape_once() -> None:
    """Fire a scrape subprocess in a background thread; auto-reaps on exit.

    The thread blocks on subprocess.run() — when the scrape finishes (or
    fails), run() reaps the child. That keeps `ps` clean and stops DuckDB
    from reporting zombie PIDs as stale lock-holders. Never raises.
    """
    repo_root = Path(__file__).resolve().parent.parent
    def _runner():
        try:
            subprocess.run(
                [sys.executable, "-m", "nfl_draft.run", "--mode", "scrape", "--book", "all"],
                cwd=str(repo_root),
                stdout=open("/tmp/nfl_draft_periodic_scrape.log", "a"),
                stderr=subprocess.STDOUT,
                check=False,
            )
        except Exception as e:
            print(f"[periodic] scrape run failed: {e}")
    threading.Thread(
        target=_runner, daemon=True,
        name="nfl_draft-scrape-runner",
    ).start()
```

The `start_new_session=True` flag is gone — no longer needed, since the runner thread itself is the owner and will wait/reap.

- [ ] **Step 2: Replace the startup-scrape block in `__main__`**

In `kalshi_draft/app.py`, in the `__main__` block (approx lines 1513-1528), replace:

```python
if __name__ == "__main__":
    import os
    # Kick off a background scrape so the dashboard has fresh data
    # within ~2 min of launch. Detached so it doesn't block startup.
    repo_root = Path(__file__).resolve().parent.parent
    try:
        subprocess.Popen(
            [sys.executable, "-m", "nfl_draft.run", "--mode", "scrape", "--book", "all"],
            cwd=str(repo_root),
            stdout=open("/tmp/nfl_draft_startup_scrape.log", "a"),
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        print("[startup] background scrape kicked off — see /tmp/nfl_draft_startup_scrape.log")
    except Exception as e:
        print(f"[startup] could not launch scrape: {e}")
```

With:

```python
if __name__ == "__main__":
    import os
    # Kick off a background scrape so the dashboard has fresh data within
    # ~2 min of launch. subprocess.run inside a daemon thread so the scrape
    # child is reaped when it exits (no zombies) without blocking startup.
    repo_root = Path(__file__).resolve().parent.parent
    def _startup_scrape_runner():
        try:
            subprocess.run(
                [sys.executable, "-m", "nfl_draft.run", "--mode", "scrape", "--book", "all"],
                cwd=str(repo_root),
                stdout=open("/tmp/nfl_draft_startup_scrape.log", "a"),
                stderr=subprocess.STDOUT,
                check=False,
            )
        except Exception as e:
            print(f"[startup] scrape run failed: {e}")
    try:
        threading.Thread(
            target=_startup_scrape_runner, daemon=True,
            name="nfl_draft-startup-scrape",
        ).start()
        print("[startup] background scrape kicked off — see /tmp/nfl_draft_startup_scrape.log")
    except Exception as e:
        print(f"[startup] could not launch scrape: {e}")
```

Note: `threading` is already imported at the top of the file (used for `_periodic_scrape_loop` etc.). If your local check shows otherwise, import it at the top of the file alongside the other stdlib imports.

- [ ] **Step 3: Import smoke test**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -c "import kalshi_draft.app"
```

Expected: no output, exit code 0.

- [ ] **Step 4: Full-suite regression**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest nfl_draft/tests/ -q
```

Expected: same pass count as after Task 4 (159 passed + the known flaky + 1 skipped). No new failures.

- [ ] **Step 5: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
git add kalshi_draft/app.py && \
git commit -m "$(cat <<'EOF'
fix(app): reap scrape subprocess with subprocess.run in a daemon thread

subprocess.Popen (fire-and-forget) left zombie children until the
dashboard exited. DuckDB reports those PIDs as stale lock-holders in its
conflicting-lock error, masking the real writer. Replacing with
subprocess.run inside a daemon thread auto-waits and reaps the child.

Applies to both _run_scrape_once (periodic 15-min tick) and the __main__
startup-scrape block. start_new_session=True is dropped; no longer needed
since the thread owns and waits the child.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 7: Document the concurrency model in the README

**Files:**
- Modify: `nfl_draft/README.md`

- [ ] **Step 1: Locate insertion point**

In `nfl_draft/README.md`, the file already has a Setup section starting at line 52. Insert the new "Concurrency model" subsection immediately before `## Setup` (i.e., right after the existing spec-link paragraph at line 48-50). A blank line between the new subsection and `## Setup` keeps rendering clean.

- [ ] **Step 2: Insert the new subsection**

Insert this block before `## Setup`:

```markdown
## Concurrency model

DuckDB allows exactly **one writer across processes** at a time. The
dashboard (`kalshi_draft/app.py`) is the primary writer:

- **Trades poll** — in-process daemon thread, every 15 s, calls
  `nfl_draft.scrapers.kalshi.fetch_trades()`. This function opens a
  single `write_connection()` for the entire poll (all
  `INSERT OR IGNORE INTO kalshi_trades` and `kalshi_poll_state` upserts
  share one lock acquisition).
- **Periodic scrape** — every 15 min, the dashboard spawns
  `python -m nfl_draft.run --mode scrape --book all` as a subprocess
  via `subprocess.run` inside a daemon thread (auto-reaped on exit).
- **Startup scrape** — a one-shot subprocess fired the same way at
  dashboard launch.
- **Bet logging** — short, user-triggered write from the dashboard UI.

External writers (manual `python -m nfl_draft.run ...`, cron-fired
scrapes) race the dashboard for the write lock. To absorb transient
collisions, `write_connection()` retries on lock-signature errors with
exponential backoff (100 ms → 5 s cap, 30 s total budget). Non-lock
errors propagate immediately. Readers (`read_connection()`) keep the
fast-fail sentinel behavior (`QueryLocked`) — callbacks re-render on
the next interval tick rather than blocking on a lock.

See `docs/superpowers/specs/2026-04-23-duckdb-write-race-design.md`
for the full design rationale.
```

- [ ] **Step 3: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
git add nfl_draft/README.md && \
git commit -m "$(cat <<'EOF'
docs(nfl_draft): document concurrency model

Explains the single-writer constraint, the retry layer in
write_connection(), and where the dashboard's writer paths live.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 8: Pre-merge review + manual smoke test + merge + cleanup

**Files:**
- None beyond the diff. Manual verification only.

- [ ] **Step 1: Full suite sanity run**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest nfl_draft/tests/ -v
```

Expected: all tests pass, same pass count as Task 6 plus any audit tests. No new failures beyond the known flaky `test_writer_and_reader_coexist`.

- [ ] **Step 2: Produce the full diff and exec checklist**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
git log main..HEAD --oneline && \
git diff main..HEAD --stat && \
git diff main..HEAD
```

Review against CLAUDE.md pre-merge checklist:
- **Data integrity**: no schema changes; all writes preserve existing semantics.
- **Resource safety**: retry layer has a bounded deadline; daemon threads are daemon (process exit kills them); `subprocess.run` reaps children.
- **Edge cases**: non-lock errors propagate immediately (no silent swallow); backoff capped so we don't sleep longer than `_MAX_BACKOFF_SEC`.
- **Dead code**: `_is_lock_error` duplicate removed from `queries.py`; `start_new_session=True` dropped from both subprocess sites.
- **Security**: no logs/secrets touched.
- **Docs**: README + spec updated.

- [ ] **Step 3: Manual smoke test — confirm the fix against a live dashboard**

This is the only way to prove the race is actually resolved. Do it before asking for merge approval.

1. From `main` (before merging), kill the current dashboard:
   ```bash
   DASH_PID=$(pgrep -f "kalshi_draft/app.py" | head -1)
   if [ -n "$DASH_PID" ]; then kill "$DASH_PID"; sleep 3; fi
   ```
2. Start the dashboard from the worktree:
   ```bash
   cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
   nohup /Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python kalshi_draft/app.py \
     >> /tmp/nfl_draft_dashboard.log 2>&1 &
   disown
   sleep 5
   ```
3. Wait ~30 seconds for the trades-poll to start firing, then trigger a manual scrape from another shell context:
   ```bash
   /Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m nfl_draft.run \
     --mode scrape --book kalshi 2>&1 | grep -E "lock|mapped|ERROR" | head -10
   ```
4. Expected: `[scrape] kalshi: mapped=<N>` with `N > 0` and no lock-error traceback. If a transient lock error appears but the retry layer resolves it, that's acceptable — the final line must show successful `mapped=<N>`.
5. Confirm zombie check:
   ```bash
   DASH_PID=$(pgrep -f "kalshi_draft/app.py" | head -1)
   ps -o pid,ppid,state -A | awk -v p=$DASH_PID '$2==p'
   ```
   Expected: no rows with `ZN` (zombie) state after a periodic-scrape tick completes (may need to wait up to 15 min for the first tick, or invoke `_run_scrape_once` via an interactive Python session targeting the dashboard).
6. Record the outcome before Step 4.

- [ ] **Step 4: Ask the user to approve the merge**

Report the commit list, diff summary, test results, and smoke-test outcome to the user. Do NOT merge without an explicit "yes".

- [ ] **Step 5: After approval, merge to main via fast-forward**

If the feature branch is behind `main` (main advanced during implementation), first rebase the feature branch onto `main`:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
git fetch origin && git rebase main
```
Resolve any conflicts following the same rule used in the earlier asymmetric-flag merge: prefer the feature-branch changes for files only this branch touched, keep main's changes for files this branch did not touch. After rebase, re-run the full suite.

Then fast-forward merge:
```bash
git -C /Users/callancapitolo/NFLWork checkout main && \
git -C /Users/callancapitolo/NFLWork merge --ff-only fix/duckdb-write-race && \
git -C /Users/callancapitolo/NFLWork log --oneline -10
```

Expected: fast-forward succeeds; the worktree commits now show on main.

- [ ] **Step 6: Push to remote (only after explicit user approval in Step 4)**

Run:
```bash
git -C /Users/callancapitolo/NFLWork push origin main
```

Expected: `main -> main` update reported cleanly.

- [ ] **Step 7: Remove the worktree and branch**

Run:
```bash
git -C /Users/callancapitolo/NFLWork worktree remove /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
git -C /Users/callancapitolo/NFLWork branch -d fix/duckdb-write-race && \
git -C /Users/callancapitolo/NFLWork worktree list
```

Expected: worktree gone; branch deleted.

- [ ] **Step 8: Final sanity run from main and restart the live dashboard**

```bash
cd /Users/callancapitolo/NFLWork && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest nfl_draft/tests/ -q
```

Expected: all tests pass.

Then restart the dashboard one more time from main so the live process is running the merged code:
```bash
DASH_PID=$(pgrep -f "kalshi_draft/app.py" | head -1)
if [ -n "$DASH_PID" ]; then kill "$DASH_PID"; sleep 3; fi
nohup /Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python \
  /Users/callancapitolo/NFLWork/kalshi_draft/app.py \
  >> /tmp/nfl_draft_dashboard.log 2>&1 &
disown
sleep 5
pgrep -f "kalshi_draft/app.py"
tail -n 5 /tmp/nfl_draft_dashboard.log
```

Expected: new PID listed; log shows `Dash is running on http://127.0.0.1:8090/`.
