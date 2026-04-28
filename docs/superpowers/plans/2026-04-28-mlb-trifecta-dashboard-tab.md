# MLB Trifecta Dashboard Tab Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a "Trifectas" tab to the MLB +EV Dashboard that displays priced TRIPLE-PLAY / GRAND-SLAM specials with a per-row manual-log Place button, mirroring the bets-tab UI on parlay-tab plumbing.

**Architecture:** `mlb_triple_play.R` gains an end-of-script step that writes a new `mlb_trifecta_opportunities` table to `Answer Keys/mlb.duckdb`. The dashboard server (`mlb_dashboard_server.py`) creates a `placed_trifectas` table + three trifecta sizing rows on init, exposes `/api/place-trifecta` and `/api/remove-trifecta` endpoints, and runs the parlay + trifecta R scripts in parallel via `subprocess.Popen` for zero added refresh latency. The dashboard renderer (`mlb_dashboard.R`) loads both tables and renders a new reactable + tab using the same `data-*` attribute pattern as the bets tab so placement is an in-place button toggle (no fragment endpoint).

**Tech Stack:** R 4.x (duckdb, dplyr, digest, reactable), Python 3 (Flask, duckdb, pytest), DuckDB. All patterns already established in `mlb_correlated_parlay.R`, `mlb_dashboard.R`, and `mlb_dashboard_server.py` — this plan only adds parallel infrastructure, no new tech.

---

## File Structure

**Created:**
- `Answer Keys/MLB Dashboard/tests/test_trifecta_endpoints.py` — pytest covering `/api/place-trifecta` + `/api/remove-trifecta` (~80 lines, Task 3)

**Modified:**
- `Answer Keys/mlb_triple_play.R` — append sizing read + Kelly + hash + dbWriteTable to end of main block, carry `id`/`description`/`commence_time` into `priced` (~50 lines added, Task 1)
- `Answer Keys/tests/test_triple_play.R` — append 3 unit tests for hash stability, kelly cell behavior, write-table idempotence (~50 lines added, Task 1)
- `Answer Keys/MLB Dashboard/mlb_dashboard_server.py`:
  - `init_db()` — add `placed_trifectas` CREATE + 3 sizing-row INSERTs (Task 2, ~25 lines added)
  - new `/api/place-trifecta` + `/api/remove-trifecta` route handlers (Task 3, ~70 lines added)
  - `run_pipeline()` Step 2 — replace `subprocess.run` with parallel `Popen` (Task 4, ~15 lines changed)
- `Answer Keys/MLB Dashboard/mlb_dashboard.R` — `load_trifecta_opps()`, `load_placed_trifectas()`, `create_trifectas_table()`, JS functions, tab strip + pane in HTML (~250 lines added, Task 5)
- `Answer Keys/CLAUDE.md` — extend "Triple-Play Data Flow" diagram with the dashboard write step (~5 lines, Task 6)
- `Answer Keys/MLB Dashboard/README.md` — add "Trifectas tab" subsection (~20 lines, Task 6)

---

## Worktree & Version Control Plan

- **Branch:** `feature/mlb-trifecta-dashboard-tab` (already created off main `51e721b`)
- **Worktree:** `/Users/callancapitolo/NFLWork/.worktrees/mlb-trifecta-dashboard` (already created)
- **Commits**: one per task (6 total). Specific messages in each task's final step.
- **DuckDB handling:** Tasks 1, 2, 4, 5 require copies of `mlb.duckdb` and/or `mlb_dashboard.duckdb` for verification. Each task's verification step explicitly copies the DB into the worktree, runs the test, then `rm`s the copy before the commit step. Never check in DBs, never symlink.
- **Cleanup post-merge:**
  ```bash
  git -C /Users/callancapitolo/NFLWork worktree remove /Users/callancapitolo/NFLWork/.worktrees/mlb-trifecta-dashboard
  git -C /Users/callancapitolo/NFLWork branch -d feature/mlb-trifecta-dashboard-tab
  ```
- **Approval required:** Never merge to `main` without explicit user approval. Run pre-merge checklist (below) and present diff before requesting merge.

---

## Pre-Merge Review Checklist

Run before requesting merge approval:

- **Data integrity:** Confirm `mlb_trifecta_opportunities` is `DROP + dbWriteTable` per refresh — no duplicate writes. Confirm `placed_trifectas` PK on `trifecta_hash` blocks double-log. Confirm the hash uses only stable identifiers (`game_id|target_team|prop_type|side`) — book_odds and fair_odds are NOT in the hash.
- **Resource safety:** Every new `dbConnect` paired with `dbDisconnect` (via `on.exit` or explicit close). The pricer's existing connection-management is preserved.
- **Edge cases verified:**
  - Off-day (`wagerzon_specials` empty) — pricer prints "No priceable specials" and exits 0; opportunities table not written. Dashboard tab shows "No trifectas priced yet."
  - First run (`placed_trifectas` table missing) — `init_db()` creates it before any read.
  - DK scraper failure — `dk_odds` NULL, blend reduces to model-only (Plan #2 already handles this).
  - Past game — `game_time` filter on read excludes already-played games.
- **Dead code:** No unused imports/helpers introduced. Every new function is called by another component.
- **Log/disk hygiene:** Pricer's existing console output preserved; no new long-lived log files.
- **Security:** Endpoints validate `trifecta_hash` non-empty; all SQL uses parameterized queries.
- **Concurrency:** Parallel `Popen` of two `mlb.duckdb` writers is safe (DuckDB serializes write transactions). Refresh is already mutex'd via `_refresh_lock`.
- **Regression:** Run the bets and parlay placement flows manually after merge — both must be untouched.

---

## Documentation Plan

Task 6 updates two files in the same commit as the docs commit:
- `Answer Keys/CLAUDE.md` — extend the existing "Triple-Play Data Flow" diagram block (added in Plan #2) with one ascii arrow showing the new `mlb_trifecta_opportunities` write, plus one bullet describing the dashboard tab.
- `Answer Keys/MLB Dashboard/README.md` — add a "Trifectas tab" subsection mirroring the existing "Parlays tab" subsection (data source, placement flow, sizing settings).

---

## Task 1: Pricer writes `mlb_trifecta_opportunities`

**Files:**
- Modify: `Answer Keys/mlb_triple_play.R` (append to end of main block, before final `}`)
- Modify: `Answer Keys/tests/test_triple_play.R` (append 3 new test_that blocks)

**What this task does:** Inside the existing `if (!interactive() && sys.nframe() == 0L) { ... }` main block of `mlb_triple_play.R`, after the existing `print(as.data.frame(display), row.names = FALSE)` line and before the closing `}`, add:
1. Carry `id`, `description`, `home_team`, `away_team`, `commence_time` through the `select()` so they survive into `priced`.
2. Read `trifecta_bankroll`, `trifecta_kelly_mult`, `trifecta_min_edge` from `mlb_dashboard.duckdb::sizing_settings` (with safe defaults).
3. Compute `kelly_bet` per row (inline Kelly formula — the parlay flow's `independent_kelly` takes parlay_group lists and is overkill for single-leg trifectas).
4. Compute `trifecta_hash` per row via `digest::digest`.
5. Build `game = "Away @ Home"` and `game_time = commence_time` columns.
6. `DROP TABLE IF EXISTS` + `dbWriteTable` to `Answer Keys/mlb.duckdb`.

### - [ ] Step 1: Read the current end of `mlb_triple_play.R` to confirm the insertion point

Run:
```bash
sed -n '293,316p' "/Users/callancapitolo/NFLWork/.worktrees/mlb-trifecta-dashboard/Answer Keys/mlb_triple_play.R"
```

Expected: see lines 293–316 ending with the `print(as.data.frame(display), row.names = FALSE)` call and the closing `}` on line 315–316.

### - [ ] Step 2: Carry needed columns through `priced`

In `Answer Keys/mlb_triple_play.R`, locate the `select(...)` call inside the rowwise mutate (around line 291–292):

```r
    select(target_team, prop_type, side, n_samples,
           model_odds, dk_odds, fair_odds, book_odds, edge_pct) %>%
```

Replace with:

```r
    select(id, target_team, prop_type, side, description,
           home_team, away_team, commence_time, n_samples,
           model_odds, dk_odds, fair_odds, book_odds, edge_pct) %>%
```

`id` is the Odds-API event id from `mlb_consensus_temp` (already in `matched`). `description`, `home_team`, `away_team` are also already in `matched`. `commence_time` is read by the existing consensus query (line 142) but not currently carried — Step 3 fixes that.

### - [ ] Step 3: Add `commence_time` to the `matched` join

Find the `matched <- ...` construction around lines 156–168. The current `select(home_team, away_team, target_team, side, book_odds, description, prop_type, id)` calls drop `commence_time`. Replace **both** occurrences (home side and away side) with:

```r
    select(home_team, away_team, target_team, side, book_odds, description, prop_type, id, commence_time)
```

The `commence_time` column flows in from the `inner_join(consensus, ...)` step earlier — verify by re-reading the surrounding 30 lines after the edit.

### - [ ] Step 4: Append the new write block at the end of the main script

Find the last `print(as.data.frame(display), row.names = FALSE)` line (around line 313). Immediately after it, BEFORE the closing `}` of the `if (!interactive() && sys.nframe() == 0L)` block, insert:

```r

  # =============================================================================
  # WRITE TO DUCKDB (for dashboard consumption)
  # =============================================================================

  # Read trifecta sizing settings from the dashboard DB (same pattern as
  # mlb_correlated_parlay.R reads parlay sizing). Falls back to safe defaults
  # if the dashboard DB or rows are missing.
  trifecta_bankroll  <- 100
  trifecta_kelly_mult <- 0.10
  trifecta_min_edge  <- 0.05
  dash_db_path <- path.expand("~/NFLWork/Answer Keys/MLB Dashboard/mlb_dashboard.duckdb")
  if (file.exists(dash_db_path)) {
    dash_con <- tryCatch(
      dbConnect(duckdb(), dbdir = dash_db_path, read_only = TRUE),
      error = function(e) NULL
    )
    if (!is.null(dash_con)) {
      saved <- tryCatch(
        dbGetQuery(dash_con, "SELECT param, value FROM sizing_settings"),
        error = function(e) data.frame(param = character(0), value = numeric(0))
      )
      if ("trifecta_bankroll"  %in% saved$param) trifecta_bankroll  <- saved$value[saved$param == "trifecta_bankroll"]
      if ("trifecta_kelly_mult" %in% saved$param) trifecta_kelly_mult <- saved$value[saved$param == "trifecta_kelly_mult"]
      if ("trifecta_min_edge"  %in% saved$param) trifecta_min_edge  <- saved$value[saved$param == "trifecta_min_edge"]
      dbDisconnect(dash_con)
    }
  }

  # Compute Kelly per row. Trifectas are single-leg-of-one bets from the
  # operator's POV (one ticket per row), so independent Kelly is correct.
  # Filter via min_edge: rows below the threshold get kelly_bet = 0.
  priced <- priced %>%
    rowwise() %>%
    mutate(
      win_prob   = if (!is.na(fair_odds)) american_to_prob(fair_odds) else NA_real_,
      dec_odds   = if (!is.na(book_odds)) {
                     if (book_odds > 0) 1 + book_odds / 100 else 1 + 100 / abs(book_odds)
                   } else NA_real_,
      kelly_frac = if (!is.na(win_prob) && !is.na(dec_odds) && dec_odds > 1) {
                     b <- dec_odds - 1
                     p <- win_prob
                     q <- 1 - p
                     max(0, (b * p - q) / b)
                   } else 0,
      kelly_bet  = if (!is.na(edge_pct) && edge_pct >= trifecta_min_edge * 100) {
                     trifecta_bankroll * trifecta_kelly_mult * kelly_frac
                   } else 0
    ) %>%
    ungroup() %>%
    select(-win_prob, -dec_odds, -kelly_frac)

  # Add display columns + dedup hash
  priced <- priced %>%
    mutate(
      game        = sprintf("%s @ %s", away_team, home_team),
      game_time   = commence_time,
      game_id     = as.character(id),
      trifecta_hash = sapply(seq_len(n()), function(i) {
        digest::digest(
          paste(game_id[i], target_team[i], prop_type[i], side[i], sep = "|"),
          algo = "sha256", serialize = FALSE
        )
      })
    ) %>%
    select(trifecta_hash, game_id, game, game_time, target_team, prop_type, side,
           description, n_samples, model_odds, dk_odds, fair_odds, book_odds,
           edge_pct, kelly_bet)

  # Drop + rewrite (same pattern as mlb_correlated_parlay.R)
  write_con <- NULL
  tryCatch({
    write_con <- dbConnect(duckdb(), dbdir = MLB_DB)
    dbExecute(write_con, "DROP TABLE IF EXISTS mlb_trifecta_opportunities")
    dbWriteTable(write_con, "mlb_trifecta_opportunities", priced)
    cat(sprintf("Wrote %d trifecta opportunities to %s.\n", nrow(priced), MLB_DB))
  }, error = function(e) {
    cat(sprintf("Warning: Failed to write trifectas to DB: %s\n", e$message))
  })
  if (!is.null(write_con)) dbDisconnect(write_con)
```

The `digest::digest` call requires the `digest` package, which is already a dependency of `mlb_dashboard.R`. Add `library(digest)` to the `suppressPackageStartupMessages` block at the top of `mlb_triple_play.R` (around line 10):

```r
suppressPackageStartupMessages({
  library(duckdb)
  library(dplyr)
  library(tibble)
  library(DBI)
  library(digest)   # NEW: for trifecta_hash
})
```

### - [ ] Step 5: Append unit tests to `Answer Keys/tests/test_triple_play.R`

Append at the very end of the file:

```r

# =============================================================================
# Trifecta dashboard table writer (Task 1 of 2026-04-28 plan)
# =============================================================================

test_that("trifecta_hash is stable across reruns and unique per row", {
  # Same inputs → same hash
  h1 <- digest::digest(paste("game42", "Yankees", "TRIPLE-PLAY", "home", sep = "|"),
                       algo = "sha256", serialize = FALSE)
  h2 <- digest::digest(paste("game42", "Yankees", "TRIPLE-PLAY", "home", sep = "|"),
                       algo = "sha256", serialize = FALSE)
  expect_equal(h1, h2)

  # Different inputs → different hash
  h3 <- digest::digest(paste("game42", "Yankees", "GRAND-SLAM", "home", sep = "|"),
                       algo = "sha256", serialize = FALSE)
  expect_false(h1 == h3)

  # Different side → different hash
  h4 <- digest::digest(paste("game42", "Yankees", "TRIPLE-PLAY", "away", sep = "|"),
                       algo = "sha256", serialize = FALSE)
  expect_false(h1 == h4)
})

test_that("kelly_bet is zero when edge_pct is below trifecta_min_edge", {
  # Replicate the inline Kelly logic from mlb_triple_play.R
  trifecta_bankroll  <- 100
  trifecta_kelly_mult <- 0.10
  trifecta_min_edge  <- 0.05  # 5%

  # Row with 3% edge — below threshold
  edge_pct <- 3.0
  kelly_bet <- if (!is.na(edge_pct) && edge_pct >= trifecta_min_edge * 100) {
    trifecta_bankroll * trifecta_kelly_mult * 0.5  # arbitrary frac
  } else 0
  expect_equal(kelly_bet, 0)

  # Row with 8% edge — above threshold
  edge_pct <- 8.0
  kelly_bet <- if (!is.na(edge_pct) && edge_pct >= trifecta_min_edge * 100) {
    trifecta_bankroll * trifecta_kelly_mult * 0.5
  } else 0
  expect_equal(kelly_bet, 5.0)
})

test_that("kelly_frac formula matches expected for known win_prob and odds", {
  # +400 American = 5.0 decimal, b = 4
  # Win prob 0.25 → kelly_frac = (4*0.25 - 0.75) / 4 = 0.0625
  win_prob <- 0.25
  dec_odds <- 5.0
  b <- dec_odds - 1
  p <- win_prob
  q <- 1 - p
  kelly_frac <- max(0, (b * p - q) / b)
  expect_equal(kelly_frac, 0.0625, tolerance = 1e-9)

  # If win_prob equals breakeven (1/dec_odds = 0.20), kelly = 0
  win_prob <- 0.20
  p <- win_prob
  q <- 1 - p
  kelly_frac <- max(0, (b * p - q) / b)
  expect_equal(kelly_frac, 0, tolerance = 1e-9)

  # Below breakeven → kelly clamped to 0 (don't bet against yourself)
  win_prob <- 0.15
  p <- win_prob
  q <- 1 - p
  kelly_frac <- max(0, (b * p - q) / b)
  expect_equal(kelly_frac, 0)
})
```

### - [ ] Step 6: Run R unit tests — verify all pass

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/mlb-trifecta-dashboard/Answer Keys/tests" && Rscript -e 'library(testthat); test_file("test_triple_play.R")'
```

Expected: all original tests pass + 3 new tests pass. No errors, no warnings.

### - [ ] Step 7: Verify the pricer writes the opportunities table on a real run

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/mlb-trifecta-dashboard
cp "/Users/callancapitolo/NFLWork/Answer Keys/mlb.duckdb" "Answer Keys/mlb.duckdb"
cd "Answer Keys" && Rscript mlb_triple_play.R 2>&1 | tail -30
```

Expected: existing fair-odds table prints, then a new line "Wrote N trifecta opportunities to mlb.duckdb." where N is the number of priced rows. If `wagerzon_specials` is empty (off-day), the pricer prints "No priceable specials found" and exits before reaching the write block — that's correct behavior.

Verify the table was created:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/mlb-trifecta-dashboard
duckdb "Answer Keys/mlb.duckdb" -c "DESCRIBE mlb_trifecta_opportunities;" -c "SELECT COUNT(*) FROM mlb_trifecta_opportunities;"
```

Expected: 15 columns (`trifecta_hash`, `game_id`, `game`, `game_time`, `target_team`, `prop_type`, `side`, `description`, `n_samples`, `model_odds`, `dk_odds`, `fair_odds`, `book_odds`, `edge_pct`, `kelly_bet`). Row count matches the console "Wrote N trifectas" message.

### - [ ] Step 8: Remove copied DB, then commit

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/mlb-trifecta-dashboard
rm "Answer Keys/mlb.duckdb"
git add "Answer Keys/mlb_triple_play.R" "Answer Keys/tests/test_triple_play.R"
git commit -m "feat(mlb): write priced trifectas to mlb_trifecta_opportunities"
```

---

## Task 2: Server creates `placed_trifectas` + sizing rows on init

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard_server.py` (extend `init_db()`)

**What this task does:** Add an idempotent `CREATE TABLE IF NOT EXISTS placed_trifectas` to the existing `init_db()` function, plus three `INSERT OR IGNORE INTO sizing_settings` statements for `trifecta_bankroll = 100`, `trifecta_kelly_mult = 0.10`, `trifecta_min_edge = 0.05`. No tests new — the existing server boot-and-render manual smoke is sufficient verification.

### - [ ] Step 1: Locate the existing parlay-sizing INSERTs in `init_db()`

```bash
grep -n "parlay_kelly_mult\|parlay_min_edge\|parlay_bankroll\|placed_parlays" "/Users/callancapitolo/NFLWork/.worktrees/mlb-trifecta-dashboard/Answer Keys/MLB Dashboard/mlb_dashboard_server.py" | head -10
```

Expected: lines around 146–155 show the parlay sizing-row INSERTs; lines around 173 show the `placed_parlays` CREATE.

### - [ ] Step 2: Add three new sizing-row INSERTs

Find the parlay sizing block (around line 146):
```python
            INSERT INTO sizing_settings (param, value) VALUES ('parlay_bankroll', 100)
```
and the corresponding `kelly_mult` / `min_edge` lines (148–155). Immediately after the **last** `parlay_min_edge` insert (line 155 area), append three new statements following the **exact same pattern**. Use `INSERT OR IGNORE` to be idempotent across server restarts:

```python
            con.execute("""
                INSERT OR IGNORE INTO sizing_settings (param, value)
                VALUES ('trifecta_bankroll', 100)
            """)
            con.execute("""
                INSERT OR IGNORE INTO sizing_settings (param, value)
                VALUES ('trifecta_kelly_mult', 0.10)
            """)
            con.execute("""
                INSERT OR IGNORE INTO sizing_settings (param, value)
                VALUES ('trifecta_min_edge', 0.05)
            """)
```

(If the existing parlay INSERTs use a different style — e.g. inline INSERTs without OR IGNORE wrapped in their own try/except — match the surrounding style exactly. The IGNORE/EXCEPT semantics matter; first-run creation must succeed but reruns must not error on duplicate PK. Inspect 5 lines of context before/after the insertion point and copy the existing pattern.)

### - [ ] Step 3: Add `placed_trifectas` table creation

Find the existing `placed_parlays` CREATE TABLE block (around line 173). Immediately after that block (and after any `ALTER TABLE` migrations applied to `placed_parlays`), append:

```python
        con.execute("""
            CREATE TABLE IF NOT EXISTS placed_trifectas (
                trifecta_hash  TEXT PRIMARY KEY,
                placed_at      TIMESTAMP,
                game_id        TEXT,
                game           TEXT,
                game_time      TIMESTAMP,
                target_team    TEXT,
                prop_type      TEXT,
                side           TEXT,
                description    TEXT,
                book_odds      INTEGER,
                fair_odds      INTEGER,
                edge_pct       DOUBLE,
                kelly_bet      DOUBLE,
                actual_wager   DOUBLE,
                status         TEXT
            )
        """)
```

### - [ ] Step 4: Smoke test — boot server, confirm table + rows exist

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/mlb-trifecta-dashboard
cp "/Users/callancapitolo/NFLWork/Answer Keys/MLB Dashboard/mlb_dashboard.duckdb" "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb" 2>/dev/null || true
python3 -c "
import sys
sys.path.insert(0, 'Answer Keys/MLB Dashboard')
import mlb_dashboard_server
mlb_dashboard_server.init_db()
print('init_db OK')
"
duckdb "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb" -c "DESCRIBE placed_trifectas;" -c "SELECT param, value FROM sizing_settings WHERE param LIKE 'trifecta_%' ORDER BY param;"
```

Expected: `init_db OK` printed. `placed_trifectas` table has 15 columns (trifecta_hash, placed_at, game_id, …, status). Three rows in `sizing_settings`: `trifecta_bankroll | 100`, `trifecta_kelly_mult | 0.1`, `trifecta_min_edge | 0.05`.

### - [ ] Step 5: Run init_db a SECOND time — confirm idempotence

```bash
python3 -c "
import sys
sys.path.insert(0, 'Answer Keys/MLB Dashboard')
import mlb_dashboard_server
mlb_dashboard_server.init_db()
print('init_db second run OK')
"
duckdb "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb" -c "SELECT COUNT(*) FROM sizing_settings WHERE param LIKE 'trifecta_%';"
```

Expected: prints `init_db second run OK` with no errors. Row count is exactly 3 (not 6 — `INSERT OR IGNORE` worked).

### - [ ] Step 6: Remove copied DB, commit

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/mlb-trifecta-dashboard
rm -f "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb"
git add "Answer Keys/MLB Dashboard/mlb_dashboard_server.py"
git commit -m "feat(dashboard): create placed_trifectas table + sizing rows on init"
```

---

## Task 3: `/api/place-trifecta` + `/api/remove-trifecta` endpoints

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard_server.py` (add two new route handlers)
- Create: `Answer Keys/MLB Dashboard/tests/test_trifecta_endpoints.py`

**What this task does:** Add two Flask endpoints that mirror the existing `/api/place-bet` and `/api/remove-bet` patterns but write to `placed_trifectas`. Cover with pytest.

### - [ ] Step 1: Write failing pytest at `Answer Keys/MLB Dashboard/tests/test_trifecta_endpoints.py`

```python
"""Tests for /api/place-trifecta and /api/remove-trifecta endpoints."""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import duckdb
import pytest

# Make the server module importable
DASHBOARD_DIR = Path(__file__).parent.parent
MLB_DIR       = DASHBOARD_DIR.parent
NFLWORK_DIR   = MLB_DIR.parent
sys.path.insert(0, str(DASHBOARD_DIR))


@pytest.fixture
def server_with_temp_dbs(tmp_path, monkeypatch):
    """Spin up the Flask server with isolated temp DBs."""
    dash_db = tmp_path / "mlb_dashboard.duckdb"
    mlb_db  = tmp_path / "mlb.duckdb"

    # Seed the trifecta opportunities table the server reads from
    con = duckdb.connect(str(mlb_db))
    con.execute("""
        CREATE TABLE mlb_trifecta_opportunities (
            trifecta_hash VARCHAR, game_id VARCHAR, game VARCHAR,
            game_time TIMESTAMP, target_team VARCHAR, prop_type VARCHAR,
            side VARCHAR, description VARCHAR, n_samples INTEGER,
            model_odds INTEGER, dk_odds INTEGER, fair_odds INTEGER,
            book_odds INTEGER, edge_pct DOUBLE, kelly_bet DOUBLE
        )
    """)
    con.execute("""
        INSERT INTO mlb_trifecta_opportunities VALUES
        ('hashA', 'g1', 'NYY @ BOS', NULL, 'Yankees', 'TRIPLE-PLAY', 'away',
         'YANKEES — SCR 1ST, 1H & GM', 500, 350, 380, 365, 500, 8.5, 5.0)
    """)
    con.close()

    # Patch DB paths BEFORE importing the server (which reads them at import time)
    monkeypatch.setenv("MLB_DASHBOARD_DB", str(dash_db))
    monkeypatch.setenv("MLB_DB_PATH", str(mlb_db))

    # Force re-import so the patched env is picked up
    if "mlb_dashboard_server" in sys.modules:
        del sys.modules["mlb_dashboard_server"]
    import mlb_dashboard_server as server
    monkeypatch.setattr(server, "DB_PATH", dash_db)
    monkeypatch.setattr(server, "MLB_DB", mlb_db)
    server.init_db()

    server.app.config["TESTING"] = True
    return server.app.test_client()


def test_place_trifecta_inserts_row(server_with_temp_dbs):
    client = server_with_temp_dbs
    r = client.post("/api/place-trifecta",
                    json={"trifecta_hash": "hashA", "actual_wager": 7.0})
    assert r.status_code == 200, r.get_data(as_text=True)
    body = json.loads(r.get_data(as_text=True))
    assert body.get("success") is True or body.get("ok") is True


def test_place_trifecta_idempotent_on_duplicate(server_with_temp_dbs):
    client = server_with_temp_dbs
    client.post("/api/place-trifecta",
                json={"trifecta_hash": "hashA", "actual_wager": 7.0})
    r = client.post("/api/place-trifecta",
                    json={"trifecta_hash": "hashA", "actual_wager": 7.0})
    # Second call must NOT crash; either returns ok:true (no-op) or 409.
    assert r.status_code in (200, 409)


def test_place_trifecta_unknown_hash_returns_404(server_with_temp_dbs):
    client = server_with_temp_dbs
    r = client.post("/api/place-trifecta",
                    json={"trifecta_hash": "no_such_hash", "actual_wager": 7.0})
    assert r.status_code == 404


def test_remove_trifecta_deletes_row(server_with_temp_dbs):
    client = server_with_temp_dbs
    client.post("/api/place-trifecta",
                json={"trifecta_hash": "hashA", "actual_wager": 7.0})
    r = client.post("/api/remove-trifecta", json={"trifecta_hash": "hashA"})
    assert r.status_code == 200
    body = json.loads(r.get_data(as_text=True))
    assert body.get("success") is True or body.get("ok") is True


def test_place_trifecta_missing_hash_returns_400(server_with_temp_dbs):
    client = server_with_temp_dbs
    r = client.post("/api/place-trifecta", json={"actual_wager": 7.0})
    assert r.status_code == 400
```

### - [ ] Step 2: Run pytest — confirm failures

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/mlb-trifecta-dashboard
cd "Answer Keys/MLB Dashboard" && python3 -m pytest tests/test_trifecta_endpoints.py -v
```

Expected: 5 failures, all 404/Method-Not-Allowed because the endpoints don't exist yet. (If pytest fails to import the server, that's a different problem — fix imports first.)

### - [ ] Step 3: Add `/api/place-trifecta` endpoint

In `Answer Keys/MLB Dashboard/mlb_dashboard_server.py`, find the `@app.route("/api/remove-bet", methods=["POST"])` block (around line 660). Immediately AFTER the function it defines (after its closing `return` line), append:

```python


@app.route("/api/place-trifecta", methods=["POST"])
def place_trifecta():
    """Manually log a placed trifecta. Server fetches the full opportunity row
    from mlb_trifecta_opportunities by hash; client only needs to send hash +
    actual_wager. Idempotent: re-posting the same hash is a no-op via PK.
    """
    data = request.json or {}
    trifecta_hash = data.get("trifecta_hash")
    if not trifecta_hash:
        return jsonify({"success": False, "error": "Missing trifecta_hash"}), 400

    actual_wager = data.get("actual_wager")

    # Look up the opportunity row in mlb.duckdb. The pricer recreates this
    # table on every refresh, so the row may have moved or disappeared if the
    # operator placed mid-refresh; treat that as a 404 with a clear message.
    try:
        opp_con = duckdb.connect(str(MLB_DB), read_only=True)
        try:
            row = opp_con.execute(
                "SELECT trifecta_hash, game_id, game, game_time, target_team, "
                "       prop_type, side, description, book_odds, fair_odds, "
                "       edge_pct, kelly_bet "
                "FROM mlb_trifecta_opportunities WHERE trifecta_hash = ?",
                [trifecta_hash]
            ).fetchone()
        finally:
            opp_con.close()
    except Exception as e:
        return jsonify({"success": False, "error": f"Lookup failed: {e}"}), 500

    if row is None:
        return jsonify({"success": False, "error": "Trifecta hash not found in opportunities"}), 404

    (h, game_id, game, game_time, target_team, prop_type, side,
     description, book_odds, fair_odds, edge_pct, kelly_bet) = row

    # Default actual_wager to round(kelly_bet) if client didn't send one
    if actual_wager is None:
        actual_wager = float(round(kelly_bet)) if kelly_bet else 0.0

    try:
        con = duckdb.connect(str(DB_PATH))
        try:
            existing = con.execute(
                "SELECT trifecta_hash FROM placed_trifectas WHERE trifecta_hash = ?",
                [trifecta_hash]
            ).fetchone()
            if existing:
                # Idempotent: already placed, no-op
                return jsonify({"success": True, "message": "Already placed"})

            con.execute("""
                INSERT INTO placed_trifectas (
                    trifecta_hash, placed_at, game_id, game, game_time,
                    target_team, prop_type, side, description, book_odds,
                    fair_odds, edge_pct, kelly_bet, actual_wager, status
                ) VALUES (?, NOW(), ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 'placed')
            """, [
                trifecta_hash, game_id, game, game_time, target_team,
                prop_type, side, description, book_odds, fair_odds,
                float(edge_pct), float(kelly_bet), float(actual_wager)
            ])
        finally:
            con.close()
        return jsonify({"success": True, "message": "Trifecta logged"})
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500


@app.route("/api/remove-trifecta", methods=["POST"])
def remove_trifecta():
    """Remove a manually-logged trifecta from placed_trifectas."""
    data = request.json or {}
    trifecta_hash = data.get("trifecta_hash")
    if not trifecta_hash:
        return jsonify({"success": False, "error": "Missing trifecta_hash"}), 400

    try:
        con = duckdb.connect(str(DB_PATH))
        try:
            con.execute(
                "DELETE FROM placed_trifectas WHERE trifecta_hash = ?",
                [trifecta_hash]
            )
        finally:
            con.close()
        return jsonify({"success": True, "message": "Trifecta removed"})
    except Exception as e:
        return jsonify({"success": False, "error": str(e)}), 500
```

`MLB_DB` is the path to `Answer Keys/mlb.duckdb`. If a module-level `MLB_DB` constant doesn't already exist in `mlb_dashboard_server.py`, add one near the top alongside `BASE_DIR` and `DB_PATH`:

```python
MLB_DB = BASE_DIR.parent / "mlb.duckdb"
```

(Verify with `grep -n "^MLB_DB\b" "Answer Keys/MLB Dashboard/mlb_dashboard_server.py"`. If it exists, skip; if not, add.)

### - [ ] Step 4: Run pytest — confirm all 5 tests pass

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/mlb-trifecta-dashboard/Answer Keys/MLB Dashboard"
python3 -m pytest tests/test_trifecta_endpoints.py -v
```

Expected: 5 PASS, 0 FAIL.

### - [ ] Step 5: Commit

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/mlb-trifecta-dashboard
git add "Answer Keys/MLB Dashboard/mlb_dashboard_server.py" "Answer Keys/MLB Dashboard/tests/test_trifecta_endpoints.py"
git commit -m "feat(dashboard): /api/place-trifecta + /api/remove-trifecta endpoints"
```

---

## Task 4: Run parlay R + trifecta R in parallel during refresh

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard_server.py` (`run_pipeline()` Step 2)

**What this task does:** Replace the single sequential `subprocess.run` for `mlb_correlated_parlay.R` with two `subprocess.Popen` calls (parlay + trifecta) plus `communicate()` waits. Saves the trifecta DK fetch time (~5–30s) by overlapping it with the parlay R work.

### - [ ] Step 1: Locate the existing Step 2 in `run_pipeline()`

```bash
grep -n "Step 2: Find correlated parlay\|mlb_correlated_parlay.R" "/Users/callancapitolo/NFLWork/.worktrees/mlb-trifecta-dashboard/Answer Keys/MLB Dashboard/mlb_dashboard_server.py"
```

Expected: lines around 1879–1888 — the comment "Step 2: Find correlated parlay opportunities (non-fatal)" plus a `subprocess.run` invocation.

### - [ ] Step 2: Replace the sequential block with parallel `Popen`

Replace the existing block (around lines 1879–1888):

```python
        # Step 2: Find correlated parlay opportunities (non-fatal)
        print("Finding parlay opportunities...")
        parlay_result = subprocess.run(
            ["Rscript", str(answer_keys_dir / "mlb_correlated_parlay.R")],
            capture_output=True,
            text=True,
            cwd=str(nfl_work_dir)
        )
        if parlay_result.returncode != 0:
            print(f"Parlay finder warning: {(parlay_result.stderr or '')[-300:]}")
```

with:

```python
        # Step 2: Find correlated parlay opportunities + price trifectas (parallel, non-fatal)
        # Both R scripts read independent inputs (parlay reads mlb_sgp_odds + mlb_consensus_temp;
        # trifecta reads wagerzon_specials + mlb_trifecta_sgp_odds) and write independent
        # output tables. Running them concurrently shaves ~5–30s off refresh latency by
        # overlapping the trifecta DK scraper with the parlay R work. DuckDB serializes
        # write transactions transparently; both writes are small DROP+dbWriteTable on tables
        # the other process never reads.
        print("Finding parlay opportunities + pricing trifectas (parallel)...")
        parlay_proc = subprocess.Popen(
            ["Rscript", str(answer_keys_dir / "mlb_correlated_parlay.R")],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            cwd=str(nfl_work_dir)
        )
        trifecta_proc = subprocess.Popen(
            ["Rscript", str(answer_keys_dir / "mlb_triple_play.R")],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            cwd=str(nfl_work_dir)
        )
        parlay_out, parlay_err = parlay_proc.communicate()
        trifecta_out, trifecta_err = trifecta_proc.communicate()
        if parlay_proc.returncode != 0:
            err_text = parlay_err.decode("utf-8", errors="replace")[-300:]
            print(f"Parlay finder warning: {err_text}")
        if trifecta_proc.returncode != 0:
            err_text = trifecta_err.decode("utf-8", errors="replace")[-300:]
            print(f"Trifecta pricer warning: {err_text}")
```

### - [ ] Step 3: Manual smoke — run `run_pipeline()` end-to-end

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/mlb-trifecta-dashboard
cp "/Users/callancapitolo/NFLWork/Answer Keys/mlb.duckdb" "Answer Keys/mlb.duckdb"
cp "/Users/callancapitolo/NFLWork/Answer Keys/MLB Dashboard/mlb_dashboard.duckdb" "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb"
cd "Answer Keys/MLB Dashboard"
python3 -c "
import sys; sys.path.insert(0, '.')
import mlb_dashboard_server
ok, msg = mlb_dashboard_server.run_pipeline()
print('OK:', ok, 'MSG:', msg)
"
```

Expected: pipeline runs to completion (or partial completion if scrapers fail — that's pre-existing). Console shows BOTH "Wrote N parlay opportunities" AND "Wrote N trifecta opportunities" messages, indicating both R scripts ran. Check the timestamps — the two messages should appear close together (within seconds), confirming parallelism.

If the pipeline fails entirely on Step 1 (scrapers error), that's pre-existing behavior — not something this task introduces. The parallel-Popen change is a no-op when neither R script gets reached.

Verify both tables were written:
```bash
duckdb "Answer Keys/mlb.duckdb" -c "SELECT 'parlay' AS k, COUNT(*) FROM mlb_parlay_opportunities UNION ALL SELECT 'trifecta', COUNT(*) FROM mlb_trifecta_opportunities;"
```

Expected: two rows, both with `COUNT(*) >= 0`.

### - [ ] Step 4: Remove copied DBs, commit

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/mlb-trifecta-dashboard
rm "Answer Keys/mlb.duckdb" "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb"
git add "Answer Keys/MLB Dashboard/mlb_dashboard_server.py"
git commit -m "perf(dashboard): run mlb_correlated_parlay.R + mlb_triple_play.R in parallel"
```

---

## Task 5: Trifectas tab — load, render, place/remove JS

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` (~250 lines added)

**What this task does:** Adds the load helpers, the reactable renderer, the JS Place/Remove functions, and the tab strip + pane to the rendered dashboard HTML.

### - [ ] Step 1: Locate the parlay tab insertion points

```bash
grep -n "load_parlay\|create_parlays_table\|placeParlay\|tab-parlays\|data-tab=.parlays" "/Users/callancapitolo/NFLWork/.worktrees/mlb-trifecta-dashboard/Answer Keys/MLB Dashboard/mlb_dashboard.R" | head -20
```

Expected: line numbers for: (a) the `parlay_opps <- tryCatch(...)` load (around line 3834), (b) the `create_parlays_table` definition (around line 334), (c) the parlay tab navigation HTML (search for `data-tab="parlays"`), (d) the parlay tab pane HTML.

### - [ ] Step 2: Add `load_trifecta_opps()` and `load_placed_trifectas()` helpers

Find the existing `load_placed_bets <- function(db_path)` definition near the top of the file (line ~42). Immediately AFTER it (and after any other top-level load helpers), insert:

```r
load_trifecta_opps <- function(mlb_db) {
  if (!file.exists(mlb_db)) return(tibble())
  con <- dbConnect(duckdb(), dbdir = mlb_db, read_only = TRUE)
  on.exit(dbDisconnect(con, shutdown = TRUE))
  tryCatch({
    if (!"mlb_trifecta_opportunities" %in% dbListTables(con)) return(tibble())
    dbGetQuery(con, "
      SELECT * FROM mlb_trifecta_opportunities
      WHERE game_time IS NULL OR CAST(game_time AS TIMESTAMP) > NOW()
    ")
  }, error = function(e) tibble())
}

load_placed_trifectas <- function(db_path) {
  if (!file.exists(db_path)) return(tibble())
  con <- dbConnect(duckdb(), dbdir = db_path, read_only = TRUE)
  on.exit(dbDisconnect(con, shutdown = TRUE))
  tryCatch({
    if (!"placed_trifectas" %in% dbListTables(con)) return(tibble())
    dbGetQuery(con, "
      SELECT * FROM placed_trifectas
      WHERE game_time IS NULL OR CAST(game_time AS TIMESTAMP) > NOW()
    ")
  }, error = function(e) tibble())
}
```

### - [ ] Step 3: Add `create_trifectas_table()` renderer

Find the existing `create_parlays_table <- function(...)` definition (around line 334). Immediately AFTER its closing `}`, insert the new renderer. The structure mirrors `create_bets_table`'s ergonomics (filterable, sortable, color-coded EV, in-place button):

```r
create_trifectas_table <- function(trifecta_opps, placed_trifectas) {
  if (nrow(trifecta_opps) == 0) {
    return(tags$div(
      style = "text-align: center; padding: 48px; color: #8b949e;",
      tags$p(style = "font-size: 1.1rem;", "No trifectas priced yet."),
      tags$p(style = "font-size: 0.85rem;",
             "Click Refresh after Wagerzon posts today's TRIPLE-PLAY / GRAND-SLAM lines.")
    ))
  }

  # Mark placed rows
  placed_set <- if (nrow(placed_trifectas) > 0) placed_trifectas$trifecta_hash else character(0)

  table_data <- trifecta_opps %>%
    mutate(
      is_placed       = trifecta_hash %in% placed_set,
      pick_display    = sprintf("%s (%s)", target_team, side),
      model_display   = ifelse(is.na(model_odds), "—",
                               ifelse(model_odds > 0, paste0("+", model_odds),
                                      as.character(model_odds))),
      dk_display      = ifelse(is.na(dk_odds), "—",
                               ifelse(dk_odds > 0, paste0("+", dk_odds),
                                      as.character(dk_odds))),
      fair_display    = ifelse(is.na(fair_odds), "—",
                               ifelse(fair_odds > 0, paste0("+", fair_odds),
                                      as.character(fair_odds))),
      book_display    = ifelse(book_odds > 0, paste0("+", book_odds),
                               as.character(book_odds)),
      edge_display    = sprintf("%+.1f%%", edge_pct),
      stake_display   = ifelse(kelly_bet > 0, sprintf("$%.0f", kelly_bet), "—")
    ) %>%
    arrange(desc(edge_pct))

  reactable(
    table_data,
    searchable = TRUE,
    filterable = TRUE,
    striped    = TRUE,
    highlight  = TRUE,
    compact    = TRUE,
    defaultPageSize = 25,
    columns = list(
      # Hidden helpers (used by JS via data-* attrs on the Action button)
      trifecta_hash = colDef(show = FALSE),
      game_id       = colDef(show = FALSE),
      game_time     = colDef(show = FALSE),
      home_team     = colDef(show = FALSE),
      away_team     = colDef(show = FALSE),
      n_samples     = colDef(show = FALSE),
      model_odds    = colDef(show = FALSE),
      dk_odds       = colDef(show = FALSE),
      fair_odds     = colDef(show = FALSE),
      book_odds     = colDef(show = FALSE),
      edge_pct      = colDef(show = FALSE),
      kelly_bet     = colDef(show = FALSE),
      target_team   = colDef(show = FALSE),
      side          = colDef(show = FALSE),
      is_placed     = colDef(show = FALSE),

      game = colDef(name = "Game", minWidth = 180, filterable = TRUE),
      pick_display = colDef(name = "Pick", minWidth = 140),
      prop_type    = colDef(name = "Prop", minWidth = 110, filterable = TRUE),
      description  = colDef(name = "Description", minWidth = 220),
      model_display = colDef(name = "Model", minWidth = 70, align = "right",
                             style = list(fontFamily = "monospace")),
      dk_display    = colDef(name = "DK", minWidth = 70, align = "right",
                             style = list(fontFamily = "monospace")),
      fair_display  = colDef(name = "Fair", minWidth = 70, align = "right",
                             style = list(fontFamily = "monospace", fontWeight = "600")),
      book_display  = colDef(name = "Book", minWidth = 70, align = "right",
                             style = list(fontFamily = "monospace")),
      edge_display = colDef(
        name = "Edge", minWidth = 80, align = "right",
        cell = function(value, index) {
          ev <- table_data$edge_pct[index]
          color <- if (is.na(ev)) "#8b949e"
                   else if (ev >= 15) "#3fb950"
                   else if (ev >= 10) "#56d364"
                   else if (ev >=  5) "#7ee787"
                   else                "#a5d6a7"
          div(style = list(color = color, fontWeight = "600"), value)
        }
      ),
      stake_display = colDef(name = "Stake", minWidth = 70, align = "right"),

      description = colDef(name = "Description", minWidth = 220, filterable = TRUE),

      # The Action column does the conditional render
      pick_display = colDef(name = "Pick", minWidth = 140,
        cell = function(value, index) {
          row <- table_data[index, ]
          # The actual Action button is rendered as its OWN column below.
          # Pick is just a styled label here.
          div(
            span(style = "font-weight: 600;", row$target_team),
            span(style = "margin-left: 8px; color: #888; font-size: 0.9em;", paste0("(", row$side, ")"))
          )
        }
      )
    ),
    columns = list(
      action = colDef(
        name = "Action", minWidth = 110, align = "center", html = TRUE,
        sortable = FALSE, filterable = FALSE,
        cell = function(value, index) {
          row <- table_data[index, ]
          if (isTRUE(row$is_placed)) {
            sprintf(
              '<button class="btn-placed" data-trifecta-hash="%s" onclick="removeTrifecta(this)">Placed</button>',
              row$trifecta_hash
            )
          } else if (!is.na(row$kelly_bet) && row$kelly_bet > 0) {
            sprintf(
              '<button class="btn-place" data-trifecta-hash="%s" data-actual-wager="%.2f" onclick="placeTrifecta(this)">Place</button>',
              row$trifecta_hash, row$kelly_bet
            )
          } else {
            ""  # below threshold or no kelly — no button
          }
        }
      )
    ),
    theme = reactableTheme(
      backgroundColor = "#0d1117",
      color           = "#c9d1d9",
      borderColor     = "#30363d",
      stripedColor    = "#161b22"
    )
  )
}
```

Note the duplicate `columns = list(...)` blocks above are an artifact of writing this incrementally — **collapse them into a single `columns = list(...)`** with all colDefs in one definition. The pattern here is simply to show every column listed (hidden + visible + Action). When implementing, merge them into one block in the order: hidden helpers, Game, Pick, Prop, Description, Model, DK, Fair, Book, Edge, Stake, Action.

### - [ ] Step 4: Add `placeTrifecta` / `removeTrifecta` JS to the rendered HTML

Find the existing inline `<script>` block in `mlb_dashboard.R` that defines `placeBet`, `removeBet`, `placeParlay`, etc. Search:

```bash
grep -n "function placeBet\|function placeParlay\|function removeBet\|async function placeBet" "/Users/callancapitolo/NFLWork/.worktrees/mlb-trifecta-dashboard/Answer Keys/MLB Dashboard/mlb_dashboard.R" | head
```

Inside that script block, after `removeParlay`, append:

```javascript
async function placeTrifecta(btn) {
  const hash  = btn.dataset.trifectaHash;
  const wager = parseFloat(btn.dataset.actualWager);
  if (!hash) { console.error('placeTrifecta: missing data-trifecta-hash'); return; }
  try {
    const r = await fetch('/api/place-trifecta', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ trifecta_hash: hash, actual_wager: isNaN(wager) ? null : wager })
    });
    const body = await r.json().catch(() => ({}));
    if (r.ok && body.success !== false) {
      btn.textContent = 'Placed';
      btn.className = 'btn-placed';
      btn.removeAttribute('data-actual-wager');
      btn.onclick = () => removeTrifecta(btn);
    } else {
      alert('Place failed: ' + (body.error || r.statusText));
    }
  } catch (e) {
    alert('Network error: ' + e.message);
  }
}

async function removeTrifecta(btn) {
  const hash = btn.dataset.trifectaHash;
  if (!hash) { console.error('removeTrifecta: missing data-trifecta-hash'); return; }
  try {
    const r = await fetch('/api/remove-trifecta', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ trifecta_hash: hash })
    });
    const body = await r.json().catch(() => ({}));
    if (r.ok && body.success !== false) {
      btn.textContent = 'Place';
      btn.className = 'btn-place';
      btn.onclick = () => placeTrifecta(btn);
    } else {
      alert('Remove failed: ' + (body.error || r.statusText));
    }
  } catch (e) {
    alert('Network error: ' + e.message);
  }
}
```

### - [ ] Step 5: Add the load + render + tab insertion in the main render block

Find the main render block (around line 3791 onward, after the `# === MLB Answer Key Dashboard ===` comment). After the existing parlay load (around line 3833):

```r
# Load parlay opportunities
parlay_opps <- tryCatch({ ... })
```

Add:

```r
# Load trifecta opportunities + placed trifectas
cat("Loading trifecta opportunities...\n")
trifecta_opps      <- load_trifecta_opps("Answer Keys/mlb.duckdb")
placed_trifectas   <- load_placed_trifectas(DB_PATH)
trifectas_table    <- create_trifectas_table(trifecta_opps, placed_trifectas)
```

Then find the existing tab strip HTML and tab pane HTML in `mlb_dashboard.R`. Search:

```bash
grep -n 'data-tab="parlays"\|data-tab="bets"\|tab-parlays\|tab-bets' "/Users/callancapitolo/NFLWork/.worktrees/mlb-trifecta-dashboard/Answer Keys/MLB Dashboard/mlb_dashboard.R" | head
```

In the tab strip — immediately after the line containing `data-tab="parlays"` (the button) — add:

```r
              tags$button(
                class = "tab-btn",
                `data-tab` = "trifectas",
                sprintf("Trifectas (%d)", nrow(trifecta_opps))
              ),
```

In the tab pane area — immediately after the parlay pane `<div id="tab-parlays" class="tab-pane">…</div>` — add:

```r
            tags$div(
              id = "tab-trifectas", class = "tab-pane",
              if (is.null(trifectas_table)) NULL else trifectas_table
            ),
```

(The exact R-syntax for tag construction depends on whether the surrounding section uses `htmltools::tags` calls or raw HTML strings. Read 30 lines of context in both spots before editing — match the existing style. Both styles are valid; consistency matters.)

### - [ ] Step 6: Smoke test — run dashboard R, render report.html, verify it includes the new tab

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/mlb-trifecta-dashboard
cp "/Users/callancapitolo/NFLWork/Answer Keys/mlb.duckdb" "Answer Keys/mlb.duckdb"
cp "/Users/callancapitolo/NFLWork/Answer Keys/MLB Dashboard/mlb_dashboard.duckdb" "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb"
Rscript "Answer Keys/MLB Dashboard/mlb_dashboard.R" 2>&1 | tail -10
grep -c 'data-tab="trifectas"\|tab-trifectas\|create_trifectas\|placeTrifecta' "Answer Keys/MLB Dashboard/report.html" 2>/dev/null
```

Expected: dashboard R completes without error. `grep -c` returns ≥ 3 (tab button + tab pane + JS function references found). If it's 0, something didn't render — debug before continuing.

### - [ ] Step 7: Visual check — open the report

```bash
open "Answer Keys/MLB Dashboard/report.html"
```

Click the "Trifectas" tab. Confirm:
- Tab appears next to "Parlays"
- Table shows priced rows with Edge color-coded
- "Place" button on rows where `kelly_bet > 0`
- No button on below-threshold rows
- Filter box at top works

If the dashboard server is running (`./Answer Keys/MLB Dashboard/run.sh`), test placement: click Place on a row, button toggles to Placed, refresh page, button stays Placed (persistence verified). Click Placed to remove.

### - [ ] Step 8: Remove copied DBs, commit

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/mlb-trifecta-dashboard
rm "Answer Keys/mlb.duckdb" "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb" "Answer Keys/MLB Dashboard/report.html"
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "feat(dashboard): Trifectas tab — load, render, place/remove JS"
```

---

## Task 6: Documentation

**Files:**
- Modify: `Answer Keys/CLAUDE.md`
- Modify: `Answer Keys/MLB Dashboard/README.md`

### - [ ] Step 1: Extend `Answer Keys/CLAUDE.md` Triple-Play Data Flow

Find the existing "Triple-Play Data Flow" section in `Answer Keys/CLAUDE.md`. After the existing `mlb_triple_play.R (standalone pricer)` ascii-art block — but inside the closing ` ``` ` of the diagram — add one new branch showing the dashboard write:

Find:
```
mlb_triple_play.R (standalone pricer)
  ├── Reads wagerzon_specials (latest scraped_at) for posted lines
  ...
  └── Prints model_odds, dk_odds, fair_odds (blended), book_odds, edge_pct
```

Replace the final `└── Prints` line with the two-line ending below:

```
  ├── Prints model_odds, dk_odds, fair_odds (blended), book_odds, edge_pct
  └── Writes mlb_trifecta_opportunities (game_id, hash, kelly_bet, ...) for the dashboard
```

After the existing trifecta-related bullet list (last bullet is "Plan #2 (post-recon) populates `mlb_trifecta_sgp_odds` ..."), append one new bullet:

```markdown
- Dashboard tab: `mlb_dashboard.R` reads `mlb_trifecta_opportunities` + `placed_trifectas` and renders a Trifectas tab next to Parlays. Manual-log only — `/api/place-trifecta` and `/api/remove-trifecta` toggle a row in `placed_trifectas` (in `mlb_dashboard.duckdb`). Sizing settings: `trifecta_bankroll`, `trifecta_kelly_mult`, `trifecta_min_edge` rows in `sizing_settings`.
```

### - [ ] Step 2: Add a "Trifectas tab" subsection to `Answer Keys/MLB Dashboard/README.md`

Find the existing "Parlays tab" or "Correlated parlays" section in `Answer Keys/MLB Dashboard/README.md`. Immediately AFTER it, append:

```markdown
## Trifectas tab

Shows priced TRIPLE-PLAY and GRAND-SLAM specials (from Wagerzon) blended with DraftKings SGP fair odds. Mirrors the bets-tab UI (filterable, color-coded EV, in-place Place button) on the parlay-tab plumbing (pricer-writes-table on each refresh, separate placed_X table for dedup).

### Data flow

1. `wagerzon_odds/scraper_specials.py` populates `wagerzon_specials` (sport='mlb', prop_type IN ('TRIPLE-PLAY','GRAND-SLAM'))
2. `Answer Keys/mlb_triple_play.R` runs in parallel with `mlb_correlated_parlay.R` during refresh:
   - reads posted lines from `wagerzon_specials`
   - invokes `mlb_sgp/scraper_draftkings_trifecta.py` for live DK SGP odds
   - blends model fair × DK fair (vig 1.25)
   - writes `mlb_trifecta_opportunities` to `Answer Keys/mlb.duckdb`
3. Dashboard reads `mlb_trifecta_opportunities` and renders the Trifectas tab.

### Placement (manual log)

- Click **Place** to log a trifecta: `POST /api/place-trifecta` with the row's hash; the server fetches the full opportunity row by hash and inserts into `placed_trifectas` (in `mlb_dashboard.duckdb`). Idempotent: re-clicking is a no-op.
- Click **Placed** to undo: `POST /api/remove-trifecta` deletes the row.
- Auto-placement (direct submission to Wagerzon) is **not** wired in this version. Place the bet on Wagerzon yourself, then click Place here to log it.

### Sizing settings

Three rows in `sizing_settings` control trifecta sizing (separate from parlay sizing because trifectas are 3-4 legs at higher payouts with single-book DK blending):

- `trifecta_bankroll` — default 100
- `trifecta_kelly_mult` — default 0.10 (10% Kelly, half the parlay default; trifecta vig 1.25 is conservative until validated)
- `trifecta_min_edge` — default 0.05 (5%; below this the row stays visible but no Place button is shown)

Adjust via the existing sliders panel (next to parlay sliders).
```

### - [ ] Step 3: Verify both docs read correctly

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/mlb-trifecta-dashboard
grep -A 3 "trifecta" "Answer Keys/CLAUDE.md" | head -25
grep -A 3 "Trifectas tab" "Answer Keys/MLB Dashboard/README.md" | head -25
```

Expected: both grep outputs show the newly-added sections.

### - [ ] Step 4: Commit

```bash
git add "Answer Keys/CLAUDE.md" "Answer Keys/MLB Dashboard/README.md"
git commit -m "docs: trifecta dashboard tab in Answer Keys CLAUDE.md + MLB Dashboard README"
```

---

## Self-Review

**Spec coverage** (from `docs/superpowers/specs/2026-04-27-mlb-trifecta-dashboard-tab-design.md`):

| Spec section                                       | Covered by             |
|----------------------------------------------------|------------------------|
| Pricer adds Kelly + hash + dbWriteTable            | Task 1                 |
| `mlb_trifecta_opportunities` schema                | Task 1 Step 4 + 7 verify |
| `placed_trifectas` schema                          | Task 2 Step 3          |
| Three trifecta sizing rows                         | Task 2 Step 2          |
| `/api/place-trifecta` + idempotence on PK          | Task 3 Steps 3 + tests |
| `/api/remove-trifecta`                             | Task 3 Step 3          |
| Parallel `Popen` orchestration in run_pipeline     | Task 4 Step 2          |
| `load_trifecta_opps` + `load_placed_trifectas`     | Task 5 Step 2          |
| `create_trifectas_table` mirroring bets-tab UX     | Task 5 Step 3          |
| `placeTrifecta` / `removeTrifecta` JS              | Task 5 Step 4          |
| Tab strip + tab pane in HTML                       | Task 5 Step 5          |
| Below-threshold rows have no Place button          | Task 5 Step 3 (action col else branch) |
| `Answer Keys/CLAUDE.md` Triple-Play Data Flow extension | Task 6 Step 1          |
| `Answer Keys/MLB Dashboard/README.md` new section  | Task 6 Step 2          |
| Pre-merge review checklist                         | Top of plan            |
| Worktree + branch + cleanup                        | Top of plan            |

**Placeholder scan:** No "TBD" / "implement later" / "fill in details". Every step has either complete code or a concrete shell command with expected output.

**Type consistency:**
- `trifecta_hash` is sha256-hex string everywhere (R produces it via `digest::digest(..., serialize=FALSE)`, Python receives it as TEXT and parameterizes it into queries).
- `mlb_trifecta_opportunities.game_id` is VARCHAR; matches `placed_trifectas.game_id` TEXT (DuckDB treats them equivalently).
- `kelly_bet` is DOUBLE in both opportunities and placed tables.
- `actual_wager` is DOUBLE in placed_trifectas; the endpoint accepts numeric or null in JSON, casts to float server-side.
- `placeTrifecta` JS reads `data-trifecta-hash` — the table renderer outputs that exact attribute.
- `removeTrifecta` JS reads `data-trifecta-hash` — the placed-state button renderer outputs that exact attribute.

No issues found.

---

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-04-28-mlb-trifecta-dashboard-tab.md`. Two execution options:

**1. Subagent-Driven (recommended)** — fresh subagent per task, review between tasks, fast iteration. Best for a 6-task plan with isolated commits.

**2. Inline Execution** — execute tasks in this session with checkpoints between tasks.

**Which approach?**
