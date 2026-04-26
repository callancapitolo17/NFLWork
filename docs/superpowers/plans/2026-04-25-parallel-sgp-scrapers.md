# Parallel SGP Scrapers Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Cut wall-clock time of the MLB SGP refresh block in `mlb_correlated_parlay.R` by running the 4 book scrapers (DK / FD / PX / NV) in parallel instead of sequentially, surface per-book timings via log files, and bump Novig's pricing concurrency to match peers.

**Architecture:** R orchestrator currently calls `system2(..., wait=TRUE)` four times back-to-back, so wall-clock = sum of all four scrapes. Switch to `parallel::mclapply` to fork 4 R workers that each invoke one scraper concurrently — wall-clock becomes the *max* of the four. Because all 4 scrapers write to the same `mlb.duckdb` file (and DuckDB allows only one writer at a time), `mlb_sgp/db.py` gets a small retry-with-backoff helper that handles transient lock contention at `clear_source` and `upsert_sgp_odds`. Per-scraper stdout/stderr is redirected to `mlb_sgp/logs/*.log` so timings are visible. Novig's `PARALLEL_PRICING` is bumped from 3 → 4 to match the other scrapers.

**Tech Stack:** R (`parallel` package, `system2`), Python 3 (DuckDB, `curl_cffi`), DuckDB write locking semantics.

---

## Worktree & Version Control

- **Branch:** `feat/parallel-sgp-scrapers`
- **Worktree:** `/Users/callancapitolo/NFLWork/.worktrees/parallel-sgp-scrapers` (create via `git worktree add`)
- **DuckDB note:** Tasks 1 and 2 use a throwaway `/tmp/test_sgp.duckdb` for unit tests — they do **not** need a copy of `mlb.duckdb`. Task 3 (R orchestrator change) is reviewed in worktree but **end-to-end verification (Task 5) runs from `main` after merge** to avoid copying or symlinking the live `mlb.duckdb`. (Per CLAUDE.md: never symlink DuckDB files.)
- **Commits:** One commit per task. Conventional-commit style messages.
- **Merge:** After all tasks pass and the user has approved, merge to `main` and clean up the worktree + branch.

## Files Touched

- **Modify:** `mlb_sgp/db.py` — add `_connect_with_retry()`, route write `connect()` calls through it
- **Modify:** `mlb_sgp/scraper_novig_sgp.py` — bump `PARALLEL_PRICING` from 3 to 4 (line 92)
- **Modify:** `Answer Keys/mlb_correlated_parlay.R` — replace lines 441-455 with parallel launch + log file routing
- **Modify:** `mlb_sgp/README.md` — document parallel execution and log file location
- **Create:** `mlb_sgp/test_db_retry.py` — unit tests for the retry helper
- **Create:** `mlb_sgp/logs/.gitkeep` — placeholder so the log directory exists in git
- **Modify:** `.gitignore` — ignore `mlb_sgp/logs/*.log`

## Documentation

- `mlb_sgp/README.md` — add a "Concurrency & Logs" section describing the parallel scrape, the log location, and the DuckDB retry behavior. Updated in the same merge to `main`.

---

### Task 1: Add retry-with-backoff helper in `mlb_sgp/db.py`

**Files:**
- Create: `mlb_sgp/test_db_retry.py`
- Modify: `mlb_sgp/db.py`

**Why:** When the 4 scrapers run in parallel they will occasionally collide at `clear_source()` and `upsert_sgp_odds()` because DuckDB only allows one writer per file. A small exponential-backoff retry around `duckdb.connect()` makes those collisions invisible to the scrapers.

- [ ] **Step 1: Write failing tests for the retry helper**

Create `mlb_sgp/test_db_retry.py`:

```python
"""Unit tests for db._connect_with_retry."""
import unittest
from unittest.mock import patch, MagicMock

import duckdb

import db


class TestConnectWithRetry(unittest.TestCase):
    def test_succeeds_on_first_try_when_no_lock(self):
        fake_con = MagicMock()
        with patch.object(db.duckdb, "connect", return_value=fake_con) as m:
            con = db._connect_with_retry("/tmp/test_sgp.duckdb")
        self.assertIs(con, fake_con)
        self.assertEqual(m.call_count, 1)

    def test_retries_then_succeeds_on_lock_error(self):
        fake_con = MagicMock()
        attempts = {"n": 0}

        def fake_connect(*args, **kwargs):
            attempts["n"] += 1
            if attempts["n"] < 3:
                raise duckdb.IOException("Could not set lock on file")
            return fake_con

        with patch.object(db.duckdb, "connect", side_effect=fake_connect):
            con = db._connect_with_retry(
                "/tmp/test_sgp.duckdb", base_delay=0.001
            )
        self.assertIs(con, fake_con)
        self.assertEqual(attempts["n"], 3)

    def test_raises_on_non_lock_error(self):
        with patch.object(db.duckdb, "connect",
                          side_effect=ValueError("bad path")):
            with self.assertRaises(ValueError):
                db._connect_with_retry("/bad/path", base_delay=0.001)

    def test_gives_up_after_max_attempts(self):
        with patch.object(db.duckdb, "connect",
                          side_effect=duckdb.IOException("lock")):
            with self.assertRaises(duckdb.IOException):
                db._connect_with_retry(
                    "/tmp/test_sgp.duckdb",
                    max_attempts=3,
                    base_delay=0.001,
                )


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run the test, confirm it fails**

Run:
```bash
cd /Users/callancapitolo/NFLWork/mlb_sgp && ./venv/bin/python -m unittest test_db_retry -v
```
Expected: `AttributeError: module 'db' has no attribute '_connect_with_retry'` (or similar — the helper doesn't exist yet).

- [ ] **Step 3: Add the retry helper to `db.py`**

Edit `mlb_sgp/db.py`. Update the imports near the top (lines 10-12):

```python
import duckdb
import random
import time
from pathlib import Path
from datetime import datetime
```

Then insert this helper immediately after the `MLB_DB = ...` line (currently line 20):

```python
# DuckDB allows only one writer per file. When the four SGP scrapers run in
# parallel they will occasionally collide at connect() — this helper retries
# with exponential backoff + jitter so the contention is invisible to callers.
def _connect_with_retry(db_path, *, read_only=False,
                        max_attempts=10, base_delay=0.1, max_delay=1.5):
    last_err = None
    for attempt in range(max_attempts):
        try:
            return duckdb.connect(db_path, read_only=read_only)
        except duckdb.IOException as e:
            msg = str(e).lower()
            if "lock" not in msg and "in use" not in msg:
                raise
            last_err = e
            if attempt == max_attempts - 1:
                break
            delay = min(base_delay * (2 ** attempt), max_delay)
            time.sleep(delay + random.uniform(0, 0.05))
    raise last_err
```

- [ ] **Step 4: Route the write `connect()` calls through the helper**

In `mlb_sgp/db.py`, replace the three `duckdb.connect(db_path)` write calls (currently in `ensure_table`, `clear_source`, `upsert_sgp_odds`) with `_connect_with_retry(db_path)`. Leave the read-only call in `get_sgp_odds` alone — readers don't contend with each other.

`ensure_table` (line ~39):
```python
con = _connect_with_retry(db_path)
```

`clear_source` (line ~54):
```python
con = _connect_with_retry(db_path)
```

`upsert_sgp_odds` (line ~76):
```python
con = _connect_with_retry(db_path)
```

- [ ] **Step 5: Run tests, confirm they pass**

Run:
```bash
cd /Users/callancapitolo/NFLWork/mlb_sgp && ./venv/bin/python -m unittest test_db_retry -v
```
Expected: 4 tests pass.

- [ ] **Step 6: Commit**

```bash
git add mlb_sgp/db.py mlb_sgp/test_db_retry.py
git commit -m "$(cat <<'EOF'
feat(mlb_sgp): add retry-with-backoff for DuckDB write contention

Wrap duckdb.connect() write calls in _connect_with_retry() with exponential
backoff + jitter (10 attempts, 0.1s → 1.5s cap). Read connections are
unaffected. This is a prerequisite for parallel scraper execution: with 4
scrapers writing to mlb.duckdb concurrently, transient lock collisions at
clear_source/upsert_sgp_odds would otherwise crash a scraper.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: Bump Novig `PARALLEL_PRICING` from 3 → 4

**Files:**
- Modify: `mlb_sgp/scraper_novig_sgp.py:92`

**Why:** All other scrapers run pricing at 4–6 workers. Novig was the conservative outlier. If Novig is on the critical path of the parallel scrape, bumping to 4 trims its pricing phase by ~25%. We can revert if Novig starts returning 401s.

- [ ] **Step 1: Edit the constant**

In `mlb_sgp/scraper_novig_sgp.py`, change line 92 from:
```python
PARALLEL_PRICING = 3
```
to:
```python
PARALLEL_PRICING = 4
```

- [ ] **Step 2: Smoke-test the Novig scraper standalone**

Run:
```bash
cd /Users/callancapitolo/NFLWork/mlb_sgp && ./venv/bin/python scraper_novig_sgp.py 2>&1 | tail -30
```
Expected: scraper completes without 401/403 errors, prints `Wrote N Novig SGP odds in T.Ts total`. Compare T to a recent run if possible — should be ≤ previous.

If you see 401/403 errors or a flood of "Retry error", revert to 3 and stop.

- [ ] **Step 3: Commit**

```bash
git add mlb_sgp/scraper_novig_sgp.py
git commit -m "$(cat <<'EOF'
perf(mlb_sgp): bump Novig PARALLEL_PRICING from 3 to 4

Brings Novig in line with DK (6), FD (4), and ProphetX (4). Smoke-tested
standalone — no auth errors observed.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: Parallelize the 4 scraper launches in R + route stdout to log files

**Files:**
- Modify: `Answer Keys/mlb_correlated_parlay.R:441-455`
- Create: `mlb_sgp/logs/.gitkeep`
- Modify: `.gitignore`

**Why:** Today the 4 `system2()` calls run sequentially with `wait = TRUE`, so the SGP refresh wall-clock is the *sum* of the four scrapes. `parallel::mclapply` forks 4 R workers; each invokes one scraper with `wait = TRUE`; wall-clock collapses to the *max* of the four. Per-scraper stdout/stderr goes to `mlb_sgp/logs/<scraper>.log` so we can finally see timing prints. The orchestrator prints a one-line elapsed-time summary per scraper after they all finish.

- [ ] **Step 1: Create the log directory placeholder**

```bash
mkdir -p /Users/callancapitolo/NFLWork/mlb_sgp/logs
touch /Users/callancapitolo/NFLWork/mlb_sgp/logs/.gitkeep
```

- [ ] **Step 2: Add `mlb_sgp/logs/*.log` to `.gitignore`**

In `/Users/callancapitolo/NFLWork/.gitignore`, append (preserving any existing trailing newline):

```
# SGP scraper per-book log files (overwritten each run)
mlb_sgp/logs/*.log
```

- [ ] **Step 3: Replace the SGP refresh block in the R orchestrator**

In `Answer Keys/mlb_correlated_parlay.R`, replace lines 441–455 (the `cat("Refreshing SGP odds...` block through the matching `else { ... }` clause) with:

```r
cat("Refreshing SGP odds (DK + FD + PX + NV) in parallel...\n")
sgp_scraper_dir <- file.path(path.expand("~"), "NFLWork", "mlb_sgp")
sgp_venv_python <- file.path(sgp_scraper_dir, "venv", "bin", "python")
sgp_log_dir     <- file.path(sgp_scraper_dir, "logs")
dir.create(sgp_log_dir, showWarnings = FALSE, recursive = TRUE)

if (file.exists(sgp_venv_python)) {
  sgp_scrapers <- c(
    "scraper_draftkings_sgp.py",
    "scraper_fanduel_sgp.py",
    "scraper_prophetx_sgp.py",
    "scraper_novig_sgp.py"
  )

  sgp_t0 <- Sys.time()
  # mclapply forks 4 R workers; each one runs system2(wait=TRUE) on one
  # scraper, so wall-clock = max(per-scraper time) instead of the sum.
  # mc.preschedule = FALSE: each scraper is its own job (no batching).
  sgp_results <- parallel::mclapply(sgp_scrapers, function(scr) {
    log_file <- file.path(sgp_log_dir,
                          sub("\\.py$", ".log", scr))
    t0 <- Sys.time()
    rc <- system2(
      sgp_venv_python,
      args   = file.path(sgp_scraper_dir, scr),
      wait   = TRUE,
      stdout = log_file,
      stderr = log_file
    )
    list(scraper = scr, exit_code = rc,
         elapsed = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  }, mc.cores = 4, mc.preschedule = FALSE)
  sgp_wall <- as.numeric(difftime(Sys.time(), sgp_t0, units = "secs"))

  # Per-scraper summary: elapsed seconds + non-zero exit codes are loud.
  for (res in sgp_results) {
    if (inherits(res, "try-error")) {
      cat(sprintf("  [SGP] FORK ERROR: %s\n", as.character(res)))
      next
    }
    status <- if (res$exit_code == 0) "ok" else
              sprintf("EXIT %d", res$exit_code)
    cat(sprintf("  [SGP] %-28s %6.1fs  %s  (log: %s)\n",
                res$scraper, res$elapsed, status,
                file.path("mlb_sgp/logs", sub("\\.py$", ".log", res$scraper))))
  }
  cat(sprintf("  [SGP] Wall clock: %.1fs\n", sgp_wall))
} else {
  cat("  SGP scraper venv not found — skipping. Run: cd mlb_sgp && python -m venv venv && pip install curl_cffi duckdb\n")
}
```

- [ ] **Step 4: Sanity-check the R syntax**

Run:
```bash
Rscript -e 'parse("/Users/callancapitolo/NFLWork/Answer Keys/mlb_correlated_parlay.R")' 2>&1 | tail -5
```
Expected: no parse errors (the command is silent on success).

- [ ] **Step 5: Commit**

```bash
git add "Answer Keys/mlb_correlated_parlay.R" mlb_sgp/logs/.gitkeep .gitignore
git commit -m "$(cat <<'EOF'
perf(mlb): run 4 SGP scrapers in parallel via mclapply

Wall clock for the SGP refresh block was sum(DK + FD + PX + NV); now it's
max() of the four. Each scraper's stdout/stderr is captured to
mlb_sgp/logs/<scraper>.log so per-book timing prints are finally visible,
and the orchestrator prints a one-line summary (elapsed + exit code) per
scraper plus an overall wall-clock line. Log files overwrite each run
(gitignored).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 4: Update `mlb_sgp/README.md`

**Files:**
- Modify: `mlb_sgp/README.md`

**Why:** Per CLAUDE.md doc discipline, architectural changes ship with their docs. New behavior to document: parallel scrape, log file location, DuckDB retry.

- [ ] **Step 1: Read the existing README to find a good insertion point**

```bash
grep -n "^#" /Users/callancapitolo/NFLWork/mlb_sgp/README.md
```

- [ ] **Step 2: Add a new "Concurrency & Logs" section**

Append the following section to `mlb_sgp/README.md` (insert before any "Troubleshooting" section if one exists, otherwise at the end):

```markdown
## Concurrency & Logs

The four SGP scrapers (DK, FD, ProphetX, Novig) are launched in parallel by
`Answer Keys/mlb_correlated_parlay.R` via `parallel::mclapply` (4 forked R
workers, one per scraper). Wall-clock time of the SGP refresh block is the
*max* of the four scrape times rather than the sum.

**Per-scraper logs:** stdout + stderr from each scraper are captured to
`mlb_sgp/logs/<scraper>.log` (overwritten each run, gitignored). The
orchestrator also prints a one-line summary per scraper (elapsed seconds,
exit code, log path) plus an overall wall-clock line — check the R console
output to see which book is the slowest on a given run.

**DuckDB write contention:** All four scrapers write to `Answer Keys/mlb.duckdb`,
which DuckDB only allows one writer to open at a time. `db.py` wraps every
write `connect()` call in `_connect_with_retry()` (exponential backoff +
jitter, up to 10 attempts) so transient lock collisions between scrapers are
invisible to callers. Read connections are not retried (DuckDB allows
unlimited concurrent readers).
```

- [ ] **Step 3: Commit**

```bash
git add mlb_sgp/README.md
git commit -m "$(cat <<'EOF'
docs(mlb_sgp): document parallel scrape, log files, and DB retry

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 5: Pre-merge review and end-to-end verification

**Files:** none modified — review-only.

**Why:** Per CLAUDE.md, every feature branch gets a pre-merge executive review of the full diff, plus end-to-end verification before merging to `main`.

- [ ] **Step 1: Review the full diff**

```bash
git diff main..HEAD --stat
git diff main..HEAD
```

Walk the executive-review checklist from CLAUDE.md:
- Data integrity: `mlb_sgp_odds` writes still tagged by `source`, no double-writes
- Resource safety: `mlb_sgp/db.py` still uses `try/finally` to close connections
- Edge cases: zero-game days (mclapply with one empty result), Novig 401 floods
- Dead code: no leftover sequential `system2` calls in the R orchestrator
- Log/disk hygiene: log files overwrite (no unbounded growth), gitignored
- Security: no API keys logged

Document any findings as ISSUES TO FIX vs ACCEPTABLE RISKS in the conversation.

- [ ] **Step 2: Merge to `main` (only after explicit user approval)**

Per CLAUDE.md: never merge without the user's go-ahead. After approval:

```bash
cd /Users/callancapitolo/NFLWork
git checkout main
git merge --no-ff feat/parallel-sgp-scrapers -m "Merge feat/parallel-sgp-scrapers: parallel SGP refresh + DB retry + per-book logs"
```

- [ ] **Step 3: End-to-end pipeline run from `main`**

The real verification — run the full MLB correlated-parlay pipeline once and confirm:
1. R console prints the new `[SGP]` summary lines (one per scraper + wall-clock).
2. `mlb_sgp/logs/scraper_*.log` files exist and contain the `Wrote N XXX SGP odds in T.Ts total` lines.
3. Wall-clock for the SGP block is roughly the max of the 4 individual times (not the sum). On a typical run that's a 50–70% reduction.
4. `mlb_sgp_odds` table in `mlb.duckdb` has rows for all 4 sources (`draftkings_direct`, `fanduel_direct`, `prophetx_direct`, `novig_direct`) — confirming no scraper silently failed under contention.

Quick check after running:
```bash
duckdb "/Users/callancapitolo/NFLWork/Answer Keys/mlb.duckdb" \
  "SELECT source, COUNT(*) FROM mlb_sgp_odds GROUP BY source ORDER BY 1"
```

- [ ] **Step 4: Clean up worktree and branch**

```bash
git worktree remove /Users/callancapitolo/NFLWork/.worktrees/parallel-sgp-scrapers
git branch -d feat/parallel-sgp-scrapers
```
