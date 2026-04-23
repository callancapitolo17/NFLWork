# DuckDB Write-Race Fix — In-Process Writer Consolidation

**Date:** 2026-04-23
**Status:** Approved for implementation (revised from earlier retry-and-batch draft)
**Scope:** `kalshi_draft/app.py`, `nfl_draft/scrapers/kalshi.py`, `nfl_draft/README.md`,
associated tests.

## Motivation

DuckDB allows only one writer across processes at a time. The NFL-draft
dashboard (`kalshi_draft/app.py`) currently spawns **scrape subprocesses**
from both its periodic 15-min loop and its one-shot startup block. Those
subprocesses run `python -m nfl_draft.run --mode scrape --book all`, each
opening its own DuckDB write connection. Meanwhile an **in-process
trades-poll daemon thread** opens dozens of short-lived write connections
per 15-second cycle (once per ticker-batch, plus one per series for the
poll-state upsert). External writers (manual CLI, cron) compete with both.

Observed impact today: `draft_odds` rows for `book='kalshi'` stalled for
86 minutes because the Kalshi scraper kept losing the race to the
trades-poll's lock acquisitions. Once the new `MAX_AGE_MINUTES=20`
staleness filter landed, Kalshi vanished from the Cross-Book Grid.

## Design

### Core idea

Eliminate **all subprocess-spawned writers** so every write comes from the
same Python process (the dashboard). DuckDB allows multiple connections
from within a single process and serializes them on its internal mutex —
no cross-process race, no retry layer needed.

Three concrete changes:

### Change 1 — Replace scrape subprocesses with in-process calls

In `kalshi_draft/app.py`, `_run_scrape_once` currently does:

```python
subprocess.Popen(
    [sys.executable, "-m", "nfl_draft.run", "--mode", "scrape", "--book", "all"],
    ...
)
```

Replace with a direct function call in a daemon thread:

```python
def _run_scrape_once() -> None:
    """Fire the scrape in-process in a background thread. Never raises."""
    def _runner():
        try:
            from nfl_draft.run import run_scrape
            run_scrape("all")
        except Exception as e:
            print(f"[periodic] scrape failed: {e}")
            traceback.print_exc()
    threading.Thread(
        target=_runner, daemon=True,
        name="nfl_draft-scrape-runner",
    ).start()
```

Apply the same pattern to the one-shot startup-scrape block in `__main__`
(~line 1519). The existing `_periodic_scrape_loop` that calls
`_run_scrape_once()` every 15 min stays unchanged — it keeps calling
the new helper, which now spawns an in-process thread instead of a
subprocess. Zero subprocesses survive.

**Crash isolation.** `run_scrape` already wraps each per-book scrape in
`try/except`, so one bad scraper doesn't abort the rest. The outer
`_runner` adds a final try/except so an uncaught error in `run_scrape`
itself (e.g. `seed.run()` fails) logs and returns rather than killing
the thread pool.

**Log destination.** Old code redirected subprocess stdout/stderr to
`/tmp/nfl_draft_periodic_scrape.log` and `/tmp/nfl_draft_startup_scrape.log`.
In-process output goes to whatever the dashboard's own stdout is —
typically `/tmp/nfl_draft_dashboard.log` per current `nohup` invocation.
Acceptable collapse: scrape output interleaves with dashboard output in
one log. If users want separate log files later, that's a follow-up.

### Change 2 — Batch `fetch_trades()` writes into one connection

Same as the earlier draft. `fetch_trades` in `nfl_draft/scrapers/kalshi.py`
currently opens `write_connection()` per ticker-batch inside nested loops
(~line 587) and again per series for poll-state (~line 616). Refactor to
open one `write_connection()` at the top of the function; thread the
`wcon` through the inner loops; close once at exit.

Rationale still holds even after Change 1: within a single process,
DuckDB serializes connections via its internal mutex. Fewer connection
opens = less contention with the periodic scrape's own writes (which
also run in this process now). Also less CPU spent on DuckDB
connect/close churn on the 15-second cadence.

### Change 3 — No retry layer needed

The earlier draft added exponential-backoff retry to `write_connection()`.
Under the in-process architecture, cross-process contention is gone, so
retry is unnecessary for the dashboard's normal operation.

**CLI trade-off.** External writers (`python -m nfl_draft.run --mode scrape`
from a shell) will fail with a DuckDB lock error if the dashboard is
running. This is intentional: running the CLI writer concurrently with
the dashboard was always the ambiguous case. Document the constraint in
the README: "Stop the dashboard before running the CLI scraper, or wait
for the dashboard's next in-process cycle to fetch the data."

No retry layer is added in this change. If a future need arises (e.g.
CLI-with-dashboard becomes common), revisit as a separate spec.

### What's NOT removed

- `_periodic_scrape_loop` — still fires every 15 min, still in-process.
- `_periodic_trades_loop` — still fires every 15 s, still in-process.
- `_trades_lock` — still guards the trades poll against overlapping
  ticks. Kept.
- Dashboard's own `write_connection()` calls for bet logging
  (`kalshi_draft/app.py:1381`) — unchanged.

## Tests

### Unit — `tests/unit/test_kalshi_trades_batching.py` (new)

Same two tests as the earlier draft:
- `test_fetch_trades_opens_single_write_connection`: counter-based
  invariant; exactly 1 open regardless of ticker count.
- `test_fetch_trades_inserts_all_batches_into_shared_connection`:
  verifies the shared connection receives the expected SQL statements.

### Unit — `tests/unit/test_app_scrape_helpers.py` (new)

- `test_run_scrape_once_does_not_spawn_subprocess`: monkeypatch
  `subprocess.Popen` and `subprocess.run` to raise; confirm
  `_run_scrape_once` executes without hitting either.
- `test_run_scrape_once_invokes_run_scrape_in_thread`: monkeypatch
  `nfl_draft.run.run_scrape` with a counter; call
  `_run_scrape_once`; assert the counter increments to 1 (possibly
  after a short `thread.join`-style wait).

### Integration — manual smoke (documented, not automated)

- Restart dashboard with the new code.
- From a second shell: try `python -m nfl_draft.run --mode scrape --book kalshi`.
  Expected: lock error (documented as the tradeoff). This is NOT a failure
  — it's the intentional constraint.
- Stop the dashboard. Re-run the CLI. Expected: success (single writer).
- Restart the dashboard. Watch `draft_odds` freshness for Kalshi after
  the first periodic tick (15 min) or wait for the startup-scrape tick.
  Expected: all venues within staleness filter; no zombies in `ps`.

### Existing-test sanity

- Full `nfl_draft/tests/` suite stays green.

## Documentation

`nfl_draft/README.md` gets a new "Concurrency model" subsection before
`## Setup`. Points:

- DuckDB is single-writer across processes.
- The dashboard is the sole writer: its in-process daemon threads
  handle trades polling (15 s) and venue scraping (15 min), plus one
  startup scrape. No subprocesses.
- External writers (CLI `python -m nfl_draft.run`, cron) **must not run
  concurrently with the dashboard** — the DuckDB lock prevents it.
  Stop the dashboard first, or rely on its own in-process schedule.
- `fetch_trades()` opens a single `write_connection()` per cycle; all
  inserts and poll-state upserts share it.

## Version control

- **Branch:** `fix/duckdb-write-race`
- **Worktree:** `.worktrees/duckdb-write-race`
- **Planned commits (rough):**
  1. `refactor(app): run periodic scrape in-process instead of subprocess`
  2. `refactor(app): run startup scrape in-process instead of subprocess`
  3. `perf(nfl_draft): batch fetch_trades() writes into one connection`
  4. `docs(nfl_draft): document single-writer concurrency model`
- **Merge:** pre-merge executive review + live smoke test + explicit user
  approval before merging to `main`. Worktree + branch cleanup after
  merge.

## Non-goals

- No retry layer in `write_connection()`. (Removed from earlier draft.)
- No subprocess zombie-reap fix. (Obsolete — no subprocesses exist.)
- No IPC layer for CLI-write-through-dashboard. (YAGNI; separate spec if
  it ever becomes needed.)
- No read-path changes.
- No changes to flag rules, staleness filter, or any UI.

## Why this supersedes the earlier retry-and-batch draft

The earlier draft added retry-with-backoff on `write_connection()` and
kept the subprocess scrape design. That absorbed the race but did not
eliminate it — under heavy load the retry budget could still time out,
and subprocess zombies still needed reaping. The in-process
consolidation fixes the root cause (cross-process writers) instead,
making the retry layer and the zombie-reap fix both obsolete.

The earlier draft's commits (in `git log` history) are not cherry-picked
forward. This spec replaces it entirely.
