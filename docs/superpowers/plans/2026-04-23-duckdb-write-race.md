# DuckDB Write-Race Fix Implementation Plan (In-Process Consolidation)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate the DuckDB single-writer race by consolidating all writers into the dashboard process — no subprocess scrapes anywhere.

**Architecture:** Replace both scrape-subprocess call sites in `kalshi_draft/app.py` (periodic 15-min tick + one-shot startup) with direct in-process calls to `nfl_draft.run.run_scrape("all")` inside a daemon thread. Batch `fetch_trades()` writes into one shared `write_connection()` to reduce within-process DuckDB connect/close churn. Document the resulting constraint (no concurrent CLI scraper) in the README.

**Tech Stack:** Python 3.14, DuckDB, Dash, pytest.

**Spec:** `docs/superpowers/specs/2026-04-23-duckdb-write-race-design.md`

---

## File Structure

Files created/modified and their responsibilities:

- **Modify** `kalshi_draft/app.py` — replace the `subprocess.Popen` call in `_run_scrape_once` (lines ~38-50) and the `__main__` startup-scrape block (lines ~1513-1528) with `threading.Thread` wrappers around a direct `run_scrape("all")` call. Remove `subprocess`/`sys` imports that become unused.
- **Modify** `nfl_draft/scrapers/kalshi.py` — refactor `fetch_trades()` (lines ~528-641) so one `write_connection()` opens at the top of the function body, the inner loops share it, and it closes once at exit.
- **Modify** `nfl_draft/README.md` — add a "Concurrency model" subsection before `## Setup` explaining the single-writer constraint.
- **Create** `nfl_draft/tests/unit/test_kalshi_trades_batching.py` — single-connection invariant for `fetch_trades()`.
- **Create** `nfl_draft/tests/unit/test_app_scrape_helpers.py` — verify `_run_scrape_once` uses in-process threads, not subprocesses.

No other files. `nfl_draft/lib/db.py` and `nfl_draft/lib/queries.py` are untouched; there's no retry layer in this design.

---

### Task 1: Create worktree and feature branch

**Files:**
- None (infra only)

- [ ] **Step 1: Verify you are on `main` and the tree is clean enough**

Run:
```bash
git -C /Users/callancapitolo/NFLWork status --short && \
git -C /Users/callancapitolo/NFLWork branch --show-current
```

Expected: only the pre-existing untracked entries (e.g. `.playwright-mcp/`, stray plan drafts), and `main` for the current branch.

- [ ] **Step 2: Create the worktree**

Run:
```bash
git -C /Users/callancapitolo/NFLWork worktree add \
  /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race \
  -b fix/duckdb-write-race
```

Expected: `Preparing worktree ...` and `HEAD is now at ...` messages.

- [ ] **Step 3: Verify the worktree + branch**

Run:
```bash
git -C /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race branch --show-current
```

Expected: `fix/duckdb-write-race`.

**From Task 2 onward, run every command from inside `/Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race`.**

---

### Task 2: Replace periodic scrape subprocess with in-process thread

**Files:**
- Create: `nfl_draft/tests/unit/test_app_scrape_helpers.py`
- Modify: `kalshi_draft/app.py` (approx lines 38-50 — the `_run_scrape_once` function)

TDD order: write the two failing tests first, confirm they fail, then flip the implementation.

- [ ] **Step 1: Create `nfl_draft/tests/unit/test_app_scrape_helpers.py`**

Create the file with:

```python
"""Unit tests: _run_scrape_once must run in-process in a daemon thread,
not via a subprocess. Under the in-process writer consolidation, any
subprocess scrape would re-introduce cross-process DuckDB lock races."""
import subprocess
import threading
import time

import pytest


def _import_app():
    """Import kalshi_draft.app lazily so test collection doesn't start
    the dashboard's daemon threads before monkeypatching is in place."""
    import kalshi_draft.app as app_mod
    return app_mod


def test_run_scrape_once_does_not_spawn_subprocess(monkeypatch):
    """If _run_scrape_once ever falls back to subprocess, the patched
    Popen/run will raise and this test will fail — acting as a tripwire."""
    app_mod = _import_app()

    def _boom(*a, **kw):
        raise AssertionError("subprocess usage is forbidden under in-process consolidation")
    monkeypatch.setattr(subprocess, "Popen", _boom)
    monkeypatch.setattr(subprocess, "run", _boom)

    # Also stub run_scrape so it doesn't actually hit the network.
    import nfl_draft.run as nfl_run
    monkeypatch.setattr(nfl_run, "run_scrape", lambda book: None)

    app_mod._run_scrape_once()
    # Give the daemon thread a moment to actually call run_scrape.
    time.sleep(0.2)


def test_run_scrape_once_invokes_run_scrape_in_thread(monkeypatch):
    """The in-process refactor must delegate to nfl_draft.run.run_scrape('all')."""
    app_mod = _import_app()

    called = {"count": 0, "arg": None, "thread_name": None}

    def _fake_run_scrape(book):
        called["count"] += 1
        called["arg"] = book
        called["thread_name"] = threading.current_thread().name

    import nfl_draft.run as nfl_run
    monkeypatch.setattr(nfl_run, "run_scrape", _fake_run_scrape)

    app_mod._run_scrape_once()

    # Wait up to ~1s for the thread to run.
    for _ in range(20):
        if called["count"] >= 1:
            break
        time.sleep(0.05)

    assert called["count"] == 1
    assert called["arg"] == "all"
    # Must run in a NON-main thread so it doesn't block the caller.
    assert called["thread_name"] != "MainThread"


def test_run_scrape_once_survives_run_scrape_exception(monkeypatch):
    """If run_scrape raises, _run_scrape_once must not propagate the
    exception into the caller's thread (the _periodic_scrape_loop)."""
    app_mod = _import_app()

    def _boom(book):
        raise RuntimeError("simulated scrape failure")

    import nfl_draft.run as nfl_run
    monkeypatch.setattr(nfl_run, "run_scrape", _boom)

    # Must not raise:
    app_mod._run_scrape_once()
    time.sleep(0.2)
```

- [ ] **Step 2: Run the new tests — confirm they fail against the current subprocess-based code**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest \
  nfl_draft/tests/unit/test_app_scrape_helpers.py -v
```

Expected:
- `test_run_scrape_once_does_not_spawn_subprocess` FAILS with `AssertionError: subprocess usage is forbidden ...` because current code calls `subprocess.Popen`.
- `test_run_scrape_once_invokes_run_scrape_in_thread` FAILS because current code spawns a subprocess instead of calling `run_scrape` in-process — `called["count"]` remains 0.
- `test_run_scrape_once_survives_run_scrape_exception` likely PASSES trivially because the current subprocess path doesn't raise on subprocess.Popen failure (it catches Exception). That's fine; it's a forward-looking test that will remain green after the refactor.

The first two failures are the red-phase target.

- [ ] **Step 3: Replace `_run_scrape_once` in `kalshi_draft/app.py`**

In `kalshi_draft/app.py`, find the current `_run_scrape_once`:

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

Replace with:

```python
def _run_scrape_once() -> None:
    """Fire the scrape in-process in a background daemon thread. Never raises.

    Under in-process writer consolidation (see
    docs/superpowers/specs/2026-04-23-duckdb-write-race-design.md), all
    DuckDB writes happen in the dashboard process so DuckDB's internal
    mutex serializes them. A subprocess scrape would re-introduce the
    cross-process lock race we're fixing.
    """
    def _runner():
        try:
            from nfl_draft.run import run_scrape
            run_scrape("all")
        except Exception as e:
            import traceback
            print(f"[periodic] scrape failed: {e}")
            traceback.print_exc()
    threading.Thread(
        target=_runner, daemon=True,
        name="nfl_draft-scrape-runner",
    ).start()
```

Note: `threading` is already imported at the top of the file (line ~31). `traceback` is imported lazily inside the except so the import cost is only paid when something fails.

- [ ] **Step 4: Run the scrape-helper tests — confirm all three pass**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest \
  nfl_draft/tests/unit/test_app_scrape_helpers.py -v
```

Expected: 3 passed, 0 failed.

- [ ] **Step 5: Full-suite sanity**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest nfl_draft/tests/ -q
```

Expected: same pass count as `main` plus the three new tests (so ~157 passed), 1 known flaky `test_writer_and_reader_coexist` may or may not fail depending on timing (known pre-existing issue), 1 skipped.

- [ ] **Step 6: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
git add kalshi_draft/app.py nfl_draft/tests/unit/test_app_scrape_helpers.py && \
git commit -m "$(cat <<'EOF'
refactor(app): run periodic scrape in-process instead of subprocess

_run_scrape_once used to fire subprocess.Popen(nfl_draft.run --mode scrape
--book all). Under DuckDB's one-writer-per-process constraint, that
subprocess raced the dashboard's own in-process trades-poll and any other
writer, producing chronic "Conflicting lock is held" errors and zombie
children that DuckDB still reported as live lock-holders.

Replace with a daemon thread that calls nfl_draft.run.run_scrape('all')
directly. All writes now happen in the dashboard process, serialized by
DuckDB's internal mutex.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: Replace startup-scrape subprocess with in-process thread

**Files:**
- Modify: `kalshi_draft/app.py` (approx lines 1513-1528, the `__main__` block)

The periodic tick is fixed; the one-shot startup scrape still uses `subprocess.Popen`. Same pattern, separate concern.

- [ ] **Step 1: Replace the startup-scrape block**

In `kalshi_draft/app.py`, inside `if __name__ == "__main__":` locate:

```python
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

Replace with:

```python
    # Kick off a background scrape so the dashboard has fresh data within
    # ~2 min of launch. Runs in-process in a daemon thread so DuckDB's
    # internal mutex serializes it with the trades-poll and any future
    # periodic-scrape tick. No subprocess = no cross-process lock race.
    def _startup_scrape_runner():
        try:
            from nfl_draft.run import run_scrape
            run_scrape("all")
        except Exception as e:
            import traceback
            print(f"[startup] scrape failed: {e}")
            traceback.print_exc()
    try:
        threading.Thread(
            target=_startup_scrape_runner, daemon=True,
            name="nfl_draft-startup-scrape",
        ).start()
        print("[startup] background scrape kicked off (in-process)")
    except Exception as e:
        print(f"[startup] could not launch scrape: {e}")
```

The `repo_root = Path(...)` line above this block can also be removed if no other code in the `__main__` block uses it; leave it alone if it does.

- [ ] **Step 2: Check for now-unused imports**

If both subprocess usages are gone, `subprocess` and `sys` may be unused at the top of `kalshi_draft/app.py`. Check:

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
grep -n "^import subprocess\|^import sys\|subprocess\.\|sys\." kalshi_draft/app.py | head -30
```

If `sys` is still used anywhere (e.g. `sys.path.insert`, `sys.executable`), leave the `import sys` alone. If `subprocess` is no longer used anywhere in the file, remove `import subprocess` from the top.

- [ ] **Step 3: Import smoke test**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -c "import kalshi_draft.app"
```

Expected: no output, exit code 0. If there's an ImportError (e.g. removed an import that was still needed), restore it and retry.

- [ ] **Step 4: Full-suite sanity**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest nfl_draft/tests/ -q
```

Expected: same count as Task 2's green state.

- [ ] **Step 5: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
git add kalshi_draft/app.py && \
git commit -m "$(cat <<'EOF'
refactor(app): run startup scrape in-process instead of subprocess

The one-shot startup scrape in __main__ used the same subprocess.Popen
pattern as the periodic tick. Converting it to an in-process daemon
thread keeps the dashboard as the sole writer from the very first
second after launch — no race window during bootstrap.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 4: Batch `fetch_trades()` writes into a single connection

**Files:**
- Create: `nfl_draft/tests/unit/test_kalshi_trades_batching.py`
- Modify: `nfl_draft/scrapers/kalshi.py` (the `fetch_trades` function body, approx lines 528-641)

Still valuable post-Change-1: within a single process, DuckDB serializes connections via its internal mutex. Fewer connection opens = less CPU churn and less within-process contention with the periodic scrape.

- [ ] **Step 1: Create `nfl_draft/tests/unit/test_kalshi_trades_batching.py`**

Create the file with:

```python
"""Unit test: fetch_trades must open exactly ONE write_connection per call,
regardless of how many tickers or batches it processes.

Before this fix, fetch_trades opened a write_connection per ticker-batch
AND per series poll_state update — dozens to hundreds per 15-second cycle.
Batching into one connection reduces DuckDB connect/close churn and
within-process contention with the periodic scrape."""
from datetime import datetime
from unittest.mock import MagicMock


class _ConnCounter:
    """Stand-in for write_connection that counts context-manager entries."""
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
    """With two series and three tickers total (one batch each),
    fetch_trades must open exactly one write connection for the whole call."""
    from nfl_draft.scrapers import kalshi as kalshi_mod

    monkeypatch.setattr(
        kalshi_mod.legacy_fetcher, "discover_draft_series",
        lambda: [{"series_ticker": "S1"}, {"series_ticker": "S2"}],
    )
    monkeypatch.setattr(
        kalshi_mod, "_tickers_for_series",
        lambda rcon, st: ["T1", "T2"] if st == "S1" else ["T3"],
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
    TradeRow = kalshi_mod.TradeRow

    def _fake_parse(raw):
        if not raw:
            return []
        return [TradeRow(
            trade_id=f"id-{id(raw)}", ticker="T?",
            side="yes", price_cents=50, count=1,
            traded_at=now, fetched_at=now,
        )]
    monkeypatch.setattr(kalshi_mod, "parse_trades_response", _fake_parse)

    fake_conn = MagicMock()
    counter = _ConnCounter(fake_conn)
    monkeypatch.setattr(kalshi_mod, "write_connection", counter)

    kalshi_mod.fetch_trades()

    assert counter.enter_count == 1, (
        f"fetch_trades opened write_connection {counter.enter_count} times; "
        "must be exactly 1 regardless of ticker or batch count."
    )


def test_fetch_trades_writes_trades_and_poll_state_on_shared_conn(monkeypatch):
    """The single write connection must receive both INSERT OR IGNORE INTO
    kalshi_trades and the kalshi_poll_state upsert."""
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
    TradeRow = kalshi_mod.TradeRow
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

    exec_sqls = [call.args[0] for call in fake_conn.execute.call_args_list]
    assert any("kalshi_trades" in sql for sql in exec_sqls), (
        "Expected INSERT OR IGNORE INTO kalshi_trades; got SQL: " + repr(exec_sqls)
    )
    assert any("kalshi_poll_state" in sql for sql in exec_sqls), (
        "Expected kalshi_poll_state upsert; got SQL: " + repr(exec_sqls)
    )
```

- [ ] **Step 2: Run the new tests — confirm the first fails and the second may pass trivially**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest \
  nfl_draft/tests/unit/test_kalshi_trades_batching.py -v
```

Expected:
- `test_fetch_trades_opens_single_write_connection` FAILS with `counter.enter_count == N` where N > 1 (current code opens one per batch plus one per series). This is the red-phase target.
- `test_fetch_trades_writes_trades_and_poll_state_on_shared_conn` may PASS coincidentally because the monkeypatched `write_connection` keeps handing back the same `fake_conn`; the test is a forward-looking invariant that must still hold after the refactor.

- [ ] **Step 3: Refactor `fetch_trades()` in `nfl_draft/scrapers/kalshi.py`**

In `nfl_draft/scrapers/kalshi.py`, locate the `fetch_trades` function (starts at line 528). Replace its body from `now_local = datetime.now()` through the end of the `for s in series_list:` loop (currently ends around line 639) with a single-connection version. Show the full replacement below (from line 560 through end of the outer loop):

**Delete** the current lines 560-639 and **replace** with:

```python
    now_local = datetime.now()

    # One write connection for the entire fetch. DuckDB's in-process
    # mutex serializes this with the dashboard's other writers (trades
    # poll, periodic scrape, bet log). Dropping from one-open-per-batch
    # to one-open-per-fetch cuts CPU + mutex churn dramatically on the
    # 15-second cadence.
    with write_connection() as wcon:
        for s in series_list:
            series_ticker = s["series_ticker"]
            last_local = cursor_by_series.get(series_ticker)
            # Kalshi min_ts is EPOCH SECONDS (not ISO) - verified live.
            min_ts = _local_to_epoch_seconds(last_local) if last_local else None

            max_traded_at = last_local  # track advance

            for ticker in tickers_by_series.get(series_ticker, []):
                # Per-ticker try/except so one bad ticker can't break the poll.
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
                            # Track max traded_at for cursor advance.
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
            # Only write when we actually saw new trades (idempotent:
            # don't regress the traded_at cursor on empty windows).
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
                    # Still record that we polled (observability), but
                    # don't regress the traded_at cursor.
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

Two substantive differences from the pre-refactor code:
1. Outer `with write_connection() as wcon:` wraps the whole series loop.
2. The two inner `with write_connection() as wcon:` blocks are gone — their bodies now run directly against the outer `wcon`.

Key invariants preserved: per-ticker `try/except` (one bad ticker doesn't take down the poll), per-series `try/except` (one bad poll-state upsert doesn't take down other series), cursor advancement semantics unchanged.

- [ ] **Step 4: Run the batching tests — both must pass**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest \
  nfl_draft/tests/unit/test_kalshi_trades_batching.py -v
```

Expected: 2 passed. The counter assertion passes at `enter_count == 1`.

- [ ] **Step 5: Full-suite sanity**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest nfl_draft/tests/ -q
```

Expected: pass count from Task 3 + 2 new batching tests. No new failures.

- [ ] **Step 6: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
git add nfl_draft/scrapers/kalshi.py \
         nfl_draft/tests/unit/test_kalshi_trades_batching.py && \
git commit -m "$(cat <<'EOF'
perf(nfl_draft): batch fetch_trades() writes into one connection

Before: fetch_trades opened a fresh write_connection per ticker-batch
AND another per series for the poll_state upsert — on a 15-second
cadence that churns DuckDB's connect/close path constantly. With the
in-process consolidation landed in the prior commits, the dashboard is
already the sole writer; the trades-poll's own connection churn is the
last remaining source of intra-process contention with the periodic
scrape's writes. One connection per fetch drops it to the minimum.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 5: Document the concurrency model in the README

**Files:**
- Modify: `nfl_draft/README.md`

- [ ] **Step 1: Locate insertion point**

In `nfl_draft/README.md`, find the line containing `## Setup`. The new "Concurrency model" subsection goes immediately before it. There should already be a paragraph ending with a spec link (around line 48-50 in the current main branch); insert the new subsection after that paragraph and before the blank line that precedes `## Setup`.

- [ ] **Step 2: Insert the new subsection**

Insert this block immediately before `## Setup`:

```markdown
## Concurrency model

DuckDB allows exactly **one writer across processes** at a time. The
dashboard (`kalshi_draft/app.py`) is the **sole writer**. All writes —
venue scraping, Kalshi trade polling, bet logging — run inside the
dashboard process, serialized by DuckDB's internal mutex. No
subprocesses are spawned for scraping anywhere in normal operation.

Writers in the dashboard process:

- **Trades poll** — in-process daemon thread, every 15 s, calls
  `nfl_draft.scrapers.kalshi.fetch_trades()`. All `INSERT OR IGNORE INTO
  kalshi_trades` statements and `kalshi_poll_state` upserts for the
  cycle share a single `write_connection()`.
- **Periodic scrape** — in-process daemon thread, every 15 min, calls
  `nfl_draft.run.run_scrape("all")` directly (not a subprocess).
- **Startup scrape** — a one-shot daemon-thread scrape fired once on
  dashboard launch so data is fresh within the first few minutes.
- **Bet logging** — short, user-triggered write from the dashboard UI.

**Constraint**: running `python -m nfl_draft.run --mode scrape ...` from
the shell **while the dashboard is up** will fail with a DuckDB lock
error. That's intentional — DuckDB's single-writer guarantee enforces
it, and by design the dashboard is the canonical writer. If you need
to run the CLI scraper ad-hoc, stop the dashboard first:

```bash
pkill -f "kalshi_draft/app.py"
python -m nfl_draft.run --mode scrape --book all
# then restart the dashboard
```

See `docs/superpowers/specs/2026-04-23-duckdb-write-race-design.md`
for the full design rationale.
```

- [ ] **Step 3: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
git add nfl_draft/README.md && \
git commit -m "$(cat <<'EOF'
docs(nfl_draft): document single-writer concurrency model

Explains that the dashboard process is the sole DuckDB writer and
enumerates the in-process writer paths. Documents the intentional
constraint that CLI scrapes cannot run concurrently with the dashboard,
and how to handle that case.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 6: Pre-merge review + live smoke test + merge + cleanup

**Files:**
- None beyond the diff. Manual verification only.

- [ ] **Step 1: Full-suite sanity run**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pytest nfl_draft/tests/ -v
```

Expected: all tests pass, 5 new (3 scrape-helper + 2 batching), pass count = main's 154 + 5 = 159. Known flaky `test_writer_and_reader_coexist` may or may not fail (pre-existing DuckDB-concurrency intermittent; not caused or fixed by this work).

- [ ] **Step 2: Produce the full diff and exec checklist**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
git log main..HEAD --oneline && \
git diff main..HEAD --stat && \
git diff main..HEAD
```

Review against CLAUDE.md pre-merge checklist:
- **Data integrity**: no schema changes; trades-poll semantics preserved (per-ticker/per-series error isolation, cursor advance, `INSERT OR IGNORE` dedupe).
- **Resource safety**: daemon threads die with the process; `write_connection` context manager still closes the connection on exception.
- **Edge cases**: `run_scrape` failure in the in-process thread is caught by `_runner`'s try/except — the periodic-scrape-loop thread keeps running.
- **Dead code**: if `import subprocess` was removed from `kalshi_draft/app.py`, confirm no remaining call sites.
- **Security**: no logs/secrets touched. Note: scrape stdout now interleaves with dashboard stdout in `/tmp/nfl_draft_dashboard.log` (no longer goes to `/tmp/nfl_draft_periodic_scrape.log` / `/tmp/nfl_draft_startup_scrape.log`). Acceptable per spec.
- **Docs**: README + spec updated.

- [ ] **Step 3: Live smoke test — confirm the fix against a running dashboard**

This is the only way to prove the race is gone. Do this before requesting merge approval.

1. Kill the currently-running dashboard (if any):
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
3. Wait 30-60 seconds for the startup scrape thread to make meaningful progress, then confirm no subprocesses were spawned:
   ```bash
   DASH_PID=$(pgrep -f "kalshi_draft/app.py" | head -1)
   pgrep -f "nfl_draft.run" | head -5
   ps -eo pid,ppid,state,command | awk -v p=$DASH_PID '$2==p'
   ```
   Expected: `pgrep -f "nfl_draft.run"` returns nothing (no subprocess exists). The dashboard's child-process listing (second command) shows no entries at all — no zombies, no workers.
4. Confirm the CLI constraint is enforced:
   ```bash
   /Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m nfl_draft.run \
     --mode scrape --book kalshi 2>&1 | grep -E "lock|mapped" | head -3
   ```
   Expected: `Conflicting lock is held ...` error (INTENTIONAL — documented tradeoff). If this succeeds, the dashboard is NOT holding its lock consistently and something is wrong with the in-process consolidation.
5. Wait for the startup scrape thread to finish (check dashboard log for `[scrape] TOTAL: mapped=...`), then verify Kalshi `draft_odds` is fresh:
   ```bash
   tail -n 40 /tmp/nfl_draft_dashboard.log | grep -E "\[scrape\]|\[startup\]"
   /Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -c "
   import duckdb
   con = duckdb.connect('/Users/callancapitolo/NFLWork/nfl_draft/nfl_draft.duckdb', read_only=True)
   r = con.execute('''SELECT book, MAX(fetched_at), DATEDIFF('minute', MAX(fetched_at), NOW()) FROM draft_odds GROUP BY book ORDER BY 3''').fetchall()
   for row in r: print(f'  {row[0]:12s}  latest={row[1]}  ({row[2]} min ago)')
   con.close()
   "
   ```
   Expected: `[scrape] TOTAL: mapped=...` appears in the log; all 7 venues are within a few minutes of `NOW()`. If the read_only connection itself errors with a lock message, the dashboard's writer is holding the file — run the read again after a few seconds; DuckDB briefly blocks read-only opens during active writes.

- [ ] **Step 4: Report smoke-test results and ask the user to approve merge**

Report to the user:
- Commit list (from `git log main..HEAD --oneline`).
- Diff stat.
- Test suite pass count.
- Live smoke test outcome (no subprocesses, CLI correctly rejected, all venues fresh after startup scrape).

Do NOT merge without an explicit "yes".

- [ ] **Step 5: After approval, rebase if needed and fast-forward merge**

If `main` advanced during implementation:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/duckdb-write-race && \
git fetch origin && git rebase main
```
Resolve conflicts following the same rule used for the asymmetric-flag merge: prefer this branch's changes for files only this branch touched; keep main's for files this branch did not touch. Re-run the full suite after rebase.

Then fast-forward:
```bash
git -C /Users/callancapitolo/NFLWork checkout main && \
git -C /Users/callancapitolo/NFLWork merge --ff-only fix/duckdb-write-race && \
git -C /Users/callancapitolo/NFLWork log --oneline -10
```

Expected: fast-forward succeeds; the worktree commits are now on main.

- [ ] **Step 6: Push to remote (only after explicit user approval in Step 4)**

```bash
git -C /Users/callancapitolo/NFLWork push origin main
```

Expected: `main -> main` update reported cleanly.

- [ ] **Step 7: Remove the worktree and branch**

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

Restart the live dashboard from main so the running process is on the merged code (if a dashboard was running on the branch during the smoke test, it's already fine — but re-confirm):
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

Expected: new PID listed; log shows `Dash is running on http://127.0.0.1:8090/`; startup-scrape thread kicks off.
