# DuckDB Write-Race Fix (Dashboard Self-Scheduling)

**Date:** 2026-04-23
**Status:** Approved for implementation
**Scope:** `nfl_draft/lib/db.py`, `nfl_draft/lib/queries.py` (helper relocation),
`nfl_draft/scrapers/kalshi.py`, `kalshi_draft/app.py`, associated tests,
`nfl_draft/README.md`.

## Motivation

DuckDB allows only one writer across processes at a time. The NFL-draft
dashboard (`kalshi_draft/app.py`) runs **two independent writers** on
tight cadences, plus spawns scrape subprocesses that never get reaped:

1. `_periodic_trades_loop` calls `kalshi.fetch_trades()` in-process
   every 15 seconds. `fetch_trades()` opens a fresh `write_connection()`
   per ticker-batch and again per series poll-state update — dozens to
   hundreds of short-lived writer-lock acquisitions per cycle.
2. `_periodic_scrape_loop` fires `subprocess.Popen(...)` every 15 min
   (plus a one-shot startup scrape at boot). These subprocesses are
   never `wait()`-ed, so they linger as zombies after exit. DuckDB
   still reports those PIDs as lock-holders.
3. External writers (manual CLI `nfl_draft.run --mode scrape`, future
   cron-fired scrapes) race against #1 constantly. Any micro-window
   overlap produces a hard `Conflicting lock is held by PID <dashboard>`
   error — the writer just dies and that cycle's data is lost.

Observed impact today: `draft_odds` rows for `book='kalshi'` were 86
minutes stale; the recent `MAX_AGE_MINUTES=20` staleness filter hid
them from the Cross-Book Grid. The race is not a new bug — the new
filter just made an old race visible.

User preference (option A in brainstorming): keep the dashboard
self-scheduling. Fix the race without restructuring to a cron-only
writer.

## Design

### Change 1 — `write_connection()` retry-with-backoff

In `nfl_draft/lib/db.py`, replace the current bare-`duckdb.connect` call
with an exponential-backoff retry when the error matches the lock
signature.

```python
_INITIAL_BACKOFF_SEC = 0.1
_MAX_BACKOFF_SEC = 5.0
_DEFAULT_TIMEOUT_SEC = 30.0


def _is_lock_error(err: BaseException) -> bool:
    """True if this exception matches the DuckDB lock-contention signature.

    DuckDB surfaces lock errors as IOException with messages like
    'Could not set lock on file' / 'Conflicting lock is held'. Shared
    with queries.py read-path (single definition).
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


@contextmanager
def write_connection(timeout_sec: float = _DEFAULT_TIMEOUT_SEC) -> Iterator[duckdb.DuckDBPyConnection]:
    """Short-lived write connection with retry-on-lock-contention.

    DuckDB allows only one writer across processes. When the dashboard is
    holding the write lock (trades-poll or subprocess scrape), a concurrent
    writer from another process (manual CLI, cron) races. Instead of dying
    on the first collision, we retry with exponential backoff for up to
    ``timeout_sec`` seconds. Non-lock errors propagate immediately.
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

`queries.py` keeps its `_safe_read` fast-fail sentinel — Dash callbacks
should not block on the interval tick; the next tick will re-render.
Only write paths retry. The `_is_lock_error` helper moves from
`queries.py` into `db.py` (single source of truth); `queries.py` imports
it from there.

### Change 2 — batch `fetch_trades()` writes into a single connection

In `nfl_draft/scrapers/kalshi.py::fetch_trades`, open one
`write_connection()` at the top of the function body and thread that
`wcon` through all inner loops. All `INSERT OR IGNORE INTO kalshi_trades`
statements and the `kalshi_poll_state` upsert execute on that single
connection. Close it once at function exit.

Current (simplified):
```python
for series_ticker in series_list:
    for ticker in tickers_by_series.get(series_ticker, []):
        while True:
            batch = parse_trades_response(raw)
            if batch:
                with write_connection() as wcon:
                    for t in batch:
                        wcon.execute("INSERT OR IGNORE INTO kalshi_trades ...")
            ...
    with write_connection() as wcon:
        wcon.execute("INSERT INTO kalshi_poll_state ... ON CONFLICT DO UPDATE ...")
```

New (simplified):
```python
with write_connection() as wcon:
    for series_ticker in series_list:
        for ticker in tickers_by_series.get(series_ticker, []):
            while True:
                batch = parse_trades_response(raw)
                if batch:
                    for t in batch:
                        wcon.execute("INSERT OR IGNORE INTO kalshi_trades ...")
                ...
        wcon.execute("INSERT INTO kalshi_poll_state ... ON CONFLICT DO UPDATE ...")
```

The `try/except Exception as e` around each ticker stays so one bad
ticker cannot take down the whole poll.

Lock acquisitions drop from O(tickers × batches) per cycle to exactly
one. Retry-on-lock (Change 1) absorbs contention when another writer
happens to be first.

### Change 3 — reap scrape subprocesses

In `kalshi_draft/app.py`:

- `_run_scrape_once` (currently fires `subprocess.Popen(...)` and
  returns immediately) becomes a daemon-thread wrapper around
  `subprocess.run(...)`. `run()` auto-waits and reaps on exit.
- The one-shot startup-scrape block (~line 1519) gets the same
  treatment.

New skeleton:
```python
def _run_scrape_once() -> None:
    """Fire a detached scrape subprocess. Auto-reaps. Never raises."""
    def _runner():
        try:
            subprocess.run(
                [sys.executable, "-m", "nfl_draft.run", "--mode", "scrape", "--book", "all"],
                cwd=str(Path(__file__).resolve().parent.parent),
                stdout=open("/tmp/nfl_draft_periodic_scrape.log", "a"),
                stderr=subprocess.STDOUT,
                check=False,
            )
        except Exception as e:
            print(f"[periodic] scrape run failed: {e}")
    threading.Thread(target=_runner, daemon=True, name="nfl_draft-scrape-runner").start()
```

Startup block mirrors this pattern. The `start_new_session=True` flag
from the old Popen is no longer needed — a daemon thread running
`subprocess.run` exits cleanly when the dashboard exits, and the
subprocess it spawns is a normal child that the runner thread waits on.

### Interaction of changes

After these three changes, the lock-contention picture becomes:

- The **trades poll** opens one writer connection per 15-second cycle
  (from Change 2). Total time holding the lock per cycle: milliseconds.
- The **subprocess scrape** still spawns an external Python process,
  but only every 15 minutes. When that external process writes, it
  races the trades poll exactly once per cycle (if at all). Change 1
  makes the external writer retry until the trades poll releases;
  DuckDB's per-process internal mutex handles the case where both are
  the dashboard.
- Zombies cleared by Change 3 — `ps` and DuckDB's live-PID check agree
  again.

## Tests

### Unit — `tests/unit/test_db.py` (new file OR extend existing)

- `test_write_connection_retries_on_lock_error`: monkeypatch
  `duckdb.connect` to raise a lock-signature `IOException` twice then
  succeed; assert `write_connection()` yields a connection and the mock
  was called three times.
- `test_write_connection_gives_up_at_timeout`: monkeypatch
  `duckdb.connect` to always raise a lock-signature `IOException`;
  call `write_connection(timeout_sec=0.5)`; assert the exception
  propagates after approximately the timeout.
- `test_write_connection_propagates_non_lock_errors_immediately`:
  monkeypatch `duckdb.connect` to raise a non-lock error; assert it
  propagates without any sleep/retry.

### Unit — `tests/unit/test_kalshi_trades_batching.py` (new file)

- `test_fetch_trades_opens_single_write_connection`: mock
  `fetch_trades` dependencies (discover_draft_series, public_request,
  parse_trades_response, etc.) so the function executes with a few
  fake tickers and batches. Monkeypatch `write_connection` with a
  call-counter. Assert the counter reached exactly 1 regardless of
  ticker count.
- `test_fetch_trades_writes_same_rows_as_before`: after the
  single-connection refactor, the number and content of inserted rows
  must match the pre-refactor behavior for the same fake input.

### Integration — manual smoke (documented in spec, not automated)

- Restart dashboard with new code. From a second shell run
  `python -m nfl_draft.run --mode scrape --book kalshi` while the
  trades-poll is live. Expect success (no lock error).
- `ps --ppid <dashboard-pid>` after a periodic-scrape tick shows no
  zombie child.

### Existing-test sanity

- Run full `nfl_draft/tests/` suite; all currently-green tests must
  stay green. The pre-existing flaky
  `test_writer_and_reader_coexist` stays out of scope (known
  concurrent-DuckDB intermittent; not caused or fixed by this work).

## Documentation

- `nfl_draft/README.md`: add a short "Concurrency model" subsection under
  an existing architecture paragraph. Points:
  - DuckDB is single-writer across processes.
  - Dashboard is the primary writer (in-process trades poll + spawned
    scrape subprocess).
  - `write_connection()` retries with backoff on lock contention so
    external/CLI/cron writers self-heal.
  - `fetch_trades()` batches all inserts into one connection.
- `nfl_draft/CLAUDE.md`: does not currently exist; no update needed.
- Do **not** add a design doc pointer into the README body — the
  `docs/superpowers/specs/` file is the canonical design.

## Version control

- **Branch:** `fix/duckdb-write-race`
- **Worktree:** `.worktrees/duckdb-write-race` (created off `main`).
- **Commits (planned, rough order):**
  1. `refactor(nfl_draft): share _is_lock_error between db.py and queries.py`
  2. `fix(nfl_draft): retry write_connection() on DuckDB lock contention`
  3. `perf(nfl_draft): batch fetch_trades() writes into one connection`
  4. `fix(app): reap scrape subprocess with subprocess.run in a daemon thread`
  5. `docs(nfl_draft): document concurrency model`
  A plan will settle the exact commit structure; 3–5 commits total.
- **Merge:** pre-merge executive review per CLAUDE.md, then explicit
  user approval before merging to `main`. Worktree + branch cleanup
  after merge.

## Non-goals

- Moving scraping to cron-only (option B — rejected per user choice).
- Bringing `fetch_draft_odds()` in-process (keeps the subprocess design).
- Replacing DuckDB with a multi-writer database.
- Any UI changes.
- Any change to the read path (`read_connection`, `_safe_read`).
  Callbacks keep the fail-fast sentinel behavior.
- Any change to the flag rules from the earlier asymmetric-outlier-flag
  work.
