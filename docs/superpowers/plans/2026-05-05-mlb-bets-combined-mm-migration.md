# mlb_bets_combined → mlb_mm.duckdb Migration Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move `mlb_bets_combined` from `Answer Keys/mlb.duckdb` to `Answer Keys/mlb_mm.duckdb` so the MLB dashboard's bet loader is no longer blocked by the pipeline's long write lock on `mlb.duckdb`. Finishes the partial migration documented in `2026-04-30-mlb-mm-separation.md` (its "layer 3 / out of scope" note).

**Architecture:** `MLB.R` already opens a brief `con_mm` for `mlb_game_samples` / `mlb_samples_meta` and disconnects immediately to release the lock. We mirror that pattern at the bottom of `MLB.R` for `mlb_bets_combined`. The dashboard reader switches its `dbdir` from `mlb.duckdb` to `mlb_mm.duckdb` and gains a `tryCatch` wrapper for defense-in-depth (matching the parlay/trifecta loaders directly below it). `mlb_bets_combined` is regenerated on every pipeline run, so there is no historical data to migrate — the next pipeline run after merge populates the new location.

**Tech Stack:** R 4.5 (duckdb, DBI), DuckDB 1.x, Python 3.14 (Flask dashboard server — read-only impact: it never touches `mlb_bets_combined`).

---

## File Structure

| File | Change | Reason |
|---|---|---|
| `Answer Keys/MLB Answer Key/MLB.R` | Modify lines ~844-846 | Move write of `mlb_bets_combined` from `con_mlb` to a fresh `con_mm` |
| `Answer Keys/MLB Dashboard/mlb_dashboard.R` | Modify lines ~4467-4481 | Read from `mlb_mm.duckdb` instead of `mlb.duckdb`; wrap in `tryCatch` |
| `Answer Keys/CLAUDE.md` | Modify line 50-51 | Move `mlb_bets_combined` listing from `mlb.duckdb` row to `mlb_mm.duckdb` row |
| `Answer Keys/MLB Dashboard/README.md` | Modify lines 8, 67, 178 | Update data flow diagram, source description, troubleshooting hint |
| `Answer Keys/MLB Answer Key/README.md` | Modify line 110 | Update schema table to show new DB location |
| (one-time, no file change) | `DROP TABLE` in `mlb.duckdb` | Remove the now-orphaned table from the old DB |

No new files. No new dependencies. No schema changes — same columns, just a different file.

---

## Risks & Verification Strategy

**Risk 1: First dashboard render after merge but before first pipeline run.**
- The writer change won't actually create `mlb_bets_combined` in `mlb_mm.duckdb` until the next pipeline run.
- Mitigation: the existing `if (!"mlb_bets_combined" %in% tables)` guard already handles missing-table → empty tibble. We add `tryCatch` around the `dbConnect` itself so a missing-DB or locked-DB also degrades to empty bets list rather than crashing.

**Risk 2: `mlb_mm.duckdb` could itself be locked when a long write happens.**
- `con_mm` writes in `MLB.R` are brief (DROP + dbWriteTable + disconnect). Unlike `con_mlb`, no long-held write lock. Lock-conflict probability on `mlb_mm.duckdb` is dramatically lower.
- Mitigation: the new `tryCatch` covers this case too.

**Risk 3: A stale orphan table in `mlb.duckdb` could confuse future debugging.**
- Mitigation: Task 4 explicitly drops it.

**TDD note (transparency):** This change is structurally a config move, not new logic. There is no existing R unit test infrastructure for `MLB.R` or the dashboard's main render path — both run as scripts. Rather than introduce a fragile integration-test scaffold for a one-line `dbdir` change, verification is done by:
1. Running the actual pipeline once and checking the table moved (Task 6).
2. Running the dashboard once and confirming bets render (Task 6).
3. Reproducing the original race (start `MLB.R` in background, render dashboard during the write window) and confirming the dashboard now succeeds (Task 7).

---

## Version Control & Worktree

- **Branch:** `feature/mlb-bets-combined-mm-migration` (already created)
- **Worktree path:** `.worktrees/mlb-bets-mm/` (already created)
- **Commit structure:** One commit per task (writer change, reader change, doc updates, orphan-drop). Final integration-test smoke run does not commit unless it surfaces fixes.
- **Merge plan:** After all tasks pass and explicit user approval, merge feature branch to `main` with `--no-ff`. After merge, remove worktree (`git worktree remove .worktrees/mlb-bets-mm`) and delete branch (`git branch -d feature/mlb-bets-combined-mm-migration`).
- **No `git stash` between branches.** All work stays on this worktree.

---

## Documentation Updates Required

In the same merge to `main`:
- `Answer Keys/CLAUDE.md` — table-location list (Task 5a)
- `Answer Keys/MLB Dashboard/README.md` — data flow + troubleshooting (Task 5b)
- `Answer Keys/MLB Answer Key/README.md` — schema table (Task 5c)

Memory file to update after merge:
- `~/.claude-personal/projects/-Users-callancapitolo-NFLWork/memory/mlb_mm_separation.md` — note that the "/refresh path still vulnerable" caveat from 2026-04-30 is now resolved for `mlb_bets_combined`.

---

## Tasks

### Task 1: Move the writer in `MLB.R`

**Files:**
- Modify: `Answer Keys/MLB Answer Key/MLB.R:840-851`

**Why:** `con_mlb` is held open across the whole script and write-locks `mlb.duckdb` for minutes. We open a fresh, short-lived `con_mm` only for this final write, mirroring the pattern already used at lines 342-353 for `mlb_game_samples`.

- [ ] **Step 1: Read the existing block to confirm exact current code**

```bash
sed -n '840,851p' "Answer Keys/MLB Answer Key/MLB.R"
```

Expected output (current state):
```r
# =============================================================================
# PHASE 8: SAVE TO DUCKDB
# =============================================================================

dbExecute(con_mlb, "DROP TABLE IF EXISTS mlb_bets_combined")
dbWriteTable(con_mlb, "mlb_bets_combined", all_bets_combined)
dbDisconnect(con_mlb)
on.exit(NULL)

timer$mark("save_bets")
cat(sprintf("Saved %d bets to mlb_bets_combined table.\n", nrow(all_bets_combined)))
```

- [ ] **Step 2: Replace the block with the new writer**

Replace lines 844-846 (the three lines that write `mlb_bets_combined` and disconnect `con_mlb`) with:

```r
# Disconnect the long-held mlb.duckdb writer FIRST so it cannot block the
# brief mlb_mm.duckdb write below (no other process should ever wait on
# mlb.duckdb during this window).
dbDisconnect(con_mlb)

# Write mlb_bets_combined to mlb_mm.duckdb (moved from mlb.duckdb on
# 2026-05-05 to free the dashboard's bet loader from the pipeline's long
# write lock on mlb.duckdb — finishes the migration started 2026-04-30).
con_bets <- duckdb_connect_retry("mlb_mm.duckdb")
dbExecute(con_bets, "DROP TABLE IF EXISTS mlb_bets_combined")
dbWriteTable(con_bets, "mlb_bets_combined", all_bets_combined)
dbDisconnect(con_bets)
```

Use Edit with the unique three-line `old_string`:
```r
dbExecute(con_mlb, "DROP TABLE IF EXISTS mlb_bets_combined")
dbWriteTable(con_mlb, "mlb_bets_combined", all_bets_combined)
dbDisconnect(con_mlb)
```

- [ ] **Step 3: Verify the change reads correctly**

```bash
sed -n '840,860p' "Answer Keys/MLB Answer Key/MLB.R"
```

Expected: the new `con_bets` block appears, `con_mlb` is disconnected before the `con_bets` block opens, the trailing `on.exit(NULL)` and `timer$mark("save_bets")` lines remain.

- [ ] **Step 4: Commit**

```bash
git add "Answer Keys/MLB Answer Key/MLB.R"
git commit -m "$(cat <<'EOF'
refactor(mlb): write mlb_bets_combined to mlb_mm.duckdb

Moves the final pipeline write from the long-locked mlb.duckdb to the
short-locked mlb_mm.duckdb, so the dashboard's bet loader can read it
even while a future pipeline run is mid-write on mlb.duckdb. Disconnects
con_mlb first so the brief con_bets write never overlaps the long lock.

Finishes the partial mlb_mm separation from 2026-04-30.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 2: Repoint the dashboard reader

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R:4467-4481`

**Why:** The dashboard's only reader of `mlb_bets_combined` must follow the data. Also wrap the connect+read in `tryCatch` so a transient lock or missing DB no longer crashes the whole render — matches the pattern already used by the parlay/trifecta loaders immediately below.

- [ ] **Step 1: Read existing block**

```bash
sed -n '4465,4485p' "Answer Keys/MLB Dashboard/mlb_dashboard.R"
```

Expected (current):
```r
cat("=== MLB Answer Key Dashboard ===\n\n")

# Load bets from duckdb (saved by MLBCombine.R or equivalent)
cat("Loading bets from database...\n")
con <- dbConnect(duckdb(), dbdir = "Answer Keys/mlb.duckdb", read_only = TRUE)

# Check if table exists
tables <- dbListTables(con)
if (!"mlb_bets_combined" %in% tables) {
  all_bets <- tibble()
} else {
  all_bets <- dbGetQuery(con, "SELECT * FROM mlb_bets_combined") %>%
    filter(is.na(pt_start_time) | pt_start_time > Sys.time())
}
dbDisconnect(con, shutdown = TRUE)

cat(sprintf("Loaded %d bets\n", nrow(all_bets)))
```

- [ ] **Step 2: Replace with mlb_mm.duckdb + tryCatch wrapper**

Use Edit. `old_string`:
```r
# Load bets from duckdb (saved by MLBCombine.R or equivalent)
cat("Loading bets from database...\n")
con <- dbConnect(duckdb(), dbdir = "Answer Keys/mlb.duckdb", read_only = TRUE)

# Check if table exists
tables <- dbListTables(con)
if (!"mlb_bets_combined" %in% tables) {
  all_bets <- tibble()
} else {
  all_bets <- dbGetQuery(con, "SELECT * FROM mlb_bets_combined") %>%
    filter(is.na(pt_start_time) | pt_start_time > Sys.time())
}
dbDisconnect(con, shutdown = TRUE)
```

`new_string`:
```r
# Load bets from duckdb. Source moved from mlb.duckdb → mlb_mm.duckdb
# on 2026-05-05 so the dashboard is no longer blocked by the pipeline's
# long write lock on mlb.duckdb. tryCatch mirrors the parlay/trifecta
# loaders below: a transient lock or missing DB degrades to empty bets
# rather than crashing the whole render.
cat("Loading bets from database...\n")
all_bets <- tryCatch({
  con <- dbConnect(duckdb(), dbdir = "Answer Keys/mlb_mm.duckdb", read_only = TRUE)
  tables <- dbListTables(con)
  result <- if (!"mlb_bets_combined" %in% tables) {
    tibble()
  } else {
    dbGetQuery(con, "SELECT * FROM mlb_bets_combined") %>%
      filter(is.na(pt_start_time) | pt_start_time > Sys.time())
  }
  dbDisconnect(con, shutdown = TRUE)
  result
}, error = function(e) {
  cat(sprintf("  Warning: could not load mlb_bets_combined (%s) — rendering with empty bets\n", e$message))
  tibble()
})
```

- [ ] **Step 3: Verify**

```bash
sed -n '4465,4495p' "Answer Keys/MLB Dashboard/mlb_dashboard.R"
```

Expected: the new block appears, the `cat(sprintf("Loaded %d bets\n", ...))` call below it still works (it reads `nrow(all_bets)`).

- [ ] **Step 4: Commit**

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "$(cat <<'EOF'
fix(mlb-dashboard): read mlb_bets_combined from mlb_mm.duckdb + degrade gracefully

Repoints the dashboard's bet loader at mlb_mm.duckdb (where the writer
moved in the previous commit) and wraps it in tryCatch so a transient
DuckDB file lock — the errno 35 crash root cause — no longer halts the
whole render. Matches the resilience pattern already used by the parlay
and trifecta loaders directly below.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 3: Drop the orphan table from `mlb.duckdb`

**Files:** none (one-time DB cleanup via `Rscript -e`)

**Why:** Otherwise a stale `mlb_bets_combined` lingers in `mlb.duckdb` forever, confusing future debugging ("which one is current?"). This is a one-shot cleanup, intentionally not a script in the repo.

- [ ] **Step 1: Wait until any in-flight pipeline finishes**

```bash
ps aux | grep -E "MLB\.R|run\.py mlb" | grep -v grep
```

Expected: no matching process. If a run is in progress, wait for it to finish (running this DROP while the pipeline holds the write lock will fail with the same errno 35 we are fixing).

- [ ] **Step 2: Drop the table**

Run from the worktree root:
```bash
Rscript -e 'library(duckdb); con <- dbConnect(duckdb(), dbdir = "Answer Keys/mlb.duckdb"); dbExecute(con, "DROP TABLE IF EXISTS mlb_bets_combined"); print(dbListTables(con)); dbDisconnect(con, shutdown = TRUE)'
```

Expected: `mlb_bets_combined` is no longer in the printed table list. `mlb_team_dict`, `mlb_consensus_temp`, `mlb_odds_temp`, etc. should still be present.

- [ ] **Step 3: No commit** — this is a DB-state change, not a source change.

---

### Task 4: Update `Answer Keys/CLAUDE.md`

**Files:**
- Modify: `Answer Keys/CLAUDE.md` lines around 50-51

**Why:** The architecture doc must reflect where the table actually lives, otherwise the next time someone (you or a future agent) reads it, they'll be misled and reintroduce the bug.

- [ ] **Step 1: Show current state**

```bash
sed -n '48,53p' "Answer Keys/CLAUDE.md"
```

Expected (current):
```
- `cbb_mm.duckdb` — CBB MM export (predictions, game samples). Separate to avoid lock contention.
- `mlb.duckdb` — MLB pipeline output (`mlb_bets_combined`, `mlb_team_dict`, `mlb_consensus_temp`, `mlb_odds_temp`, `mlb_trifecta_sgp_odds`, scraper raw odds, historical PBP joins). Consumer-facing tables moved to `mlb_mm.duckdb`.
- `mlb_mm.duckdb` — MLB MM consumer tables (`mlb_game_samples`, `mlb_samples_meta`, `mlb_sgp_odds`, `mlb_parlay_lines`, `mlb_parlay_opportunities`, `mlb_trifecta_opportunities`). Separate from `mlb.duckdb` so the pipeline's write lock on `mlb.duckdb` never blocks the RFQ bot or the dashboard's bet-log loaders. Mirrors the CBB pattern.
```

- [ ] **Step 2: Edit both lines**

Use Edit twice:

`old_string` (line 50, mlb.duckdb description):
```
- `mlb.duckdb` — MLB pipeline output (`mlb_bets_combined`, `mlb_team_dict`, `mlb_consensus_temp`, `mlb_odds_temp`, `mlb_trifecta_sgp_odds`, scraper raw odds, historical PBP joins). Consumer-facing tables moved to `mlb_mm.duckdb`.
```
`new_string`:
```
- `mlb.duckdb` — MLB pipeline working tables (`mlb_team_dict`, `mlb_consensus_temp`, `mlb_odds_temp`, `mlb_trifecta_sgp_odds`, scraper raw odds, historical PBP joins). Held under a long write lock during a pipeline run; nothing the dashboard or RFQ bot needs lives here anymore.
```

`old_string` (line 51, mlb_mm.duckdb description):
```
- `mlb_mm.duckdb` — MLB MM consumer tables (`mlb_game_samples`, `mlb_samples_meta`, `mlb_sgp_odds`, `mlb_parlay_lines`, `mlb_parlay_opportunities`, `mlb_trifecta_opportunities`). Separate from `mlb.duckdb` so the pipeline's write lock on `mlb.duckdb` never blocks the RFQ bot or the dashboard's bet-log loaders. Mirrors the CBB pattern.
```
`new_string`:
```
- `mlb_mm.duckdb` — MLB MM consumer tables (`mlb_bets_combined`, `mlb_game_samples`, `mlb_samples_meta`, `mlb_sgp_odds`, `mlb_parlay_lines`, `mlb_parlay_opportunities`, `mlb_trifecta_opportunities`). Separate from `mlb.duckdb` so the pipeline's write lock on `mlb.duckdb` never blocks the RFQ bot or the dashboard. `mlb_bets_combined` moved here on 2026-05-05; the others moved 2026-04-30. Mirrors the CBB pattern.
```

- [ ] **Step 3: Verify**

```bash
grep -n "mlb_bets_combined" "Answer Keys/CLAUDE.md"
```

Expected: exactly one hit, on the `mlb_mm.duckdb` line.

- [ ] **Step 4: Commit (deferred — bundled with Task 5)**

---

### Task 5: Update READMEs

**Files:**
- Modify: `Answer Keys/MLB Dashboard/README.md` (lines 8, 67, 178)
- Modify: `Answer Keys/MLB Answer Key/README.md` (line 110)

- [ ] **Step 1: Show current state**

```bash
grep -n "mlb_bets_combined\|mlb\.duckdb" "Answer Keys/MLB Dashboard/README.md" "Answer Keys/MLB Answer Key/README.md"
```

- [ ] **Step 2: Update `Answer Keys/MLB Dashboard/README.md` line 8 (data flow diagram)**

`old_string`:
```
run.py mlb (pipeline)  →  mlb.duckdb/mlb_bets_combined             (dashboard refresh)
```
`new_string`:
```
run.py mlb (pipeline)  →  mlb_mm.duckdb/mlb_bets_combined          (dashboard refresh)
```

- [ ] **Step 3: Update `Answer Keys/MLB Dashboard/README.md` line 67 (source description)**

`old_string`:
```
- Pipeline bets read from `Answer Keys/mlb.duckdb` (`mlb_bets_combined`) and `Answer Keys/mlb_mm.duckdb` (`mlb_parlay_opportunities`, `mlb_trifecta_opportunities`)
```
`new_string`:
```
- Pipeline bets read from `Answer Keys/mlb_mm.duckdb` (`mlb_bets_combined`, `mlb_parlay_opportunities`, `mlb_trifecta_opportunities`) — all consumer tables now in one DB to avoid contention with the pipeline's long write lock on `mlb.duckdb`
```

- [ ] **Step 4: Update `Answer Keys/MLB Dashboard/README.md` line 178 (troubleshooting)**

`old_string`:
```
- **No bets shown** — check `mlb_bets_combined` has rows; if empty, check MLB.R output in pipeline logs
```
`new_string`:
```
- **No bets shown** — check `mlb_bets_combined` has rows in `mlb_mm.duckdb`; if empty, check MLB.R output in pipeline logs
```

- [ ] **Step 5: Update `Answer Keys/MLB Answer Key/README.md` line 110 (schema table)**

`old_string`:
```
| `mlb.duckdb` | `mlb_bets_combined` | Pipeline output (daily +EV bets) |
```
`new_string`:
```
| `mlb_mm.duckdb` | `mlb_bets_combined` | Pipeline output (daily +EV bets) — moved from `mlb.duckdb` on 2026-05-05 to avoid lock contention |
```

- [ ] **Step 6: Verify all four README updates**

```bash
grep -n "mlb_bets_combined" "Answer Keys/MLB Dashboard/README.md" "Answer Keys/MLB Answer Key/README.md"
```

Expected: all hits should now reference `mlb_mm.duckdb`, none should reference `mlb.duckdb` next to `mlb_bets_combined`.

- [ ] **Step 7: Commit Tasks 4 + 5 together**

```bash
git add "Answer Keys/CLAUDE.md" "Answer Keys/MLB Dashboard/README.md" "Answer Keys/MLB Answer Key/README.md"
git commit -m "$(cat <<'EOF'
docs(mlb): record mlb_bets_combined location move to mlb_mm.duckdb

Updates the architecture doc and the dashboard / answer-key READMEs so
they reflect where the table actually lives after the writer + reader
changes.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

### Task 6: Smoke test — pipeline run + dashboard render

**Files:** none — this is verification, not a code change.

- [ ] **Step 1: Confirm no in-flight pipeline before starting**

```bash
ps aux | grep -E "MLB\.R|run\.py mlb" | grep -v grep
```

Expected: empty. If something is running, wait.

- [ ] **Step 2: Run the pipeline once from the worktree**

⚠️ **DuckDB symlink rule:** Per `CLAUDE.md`, never symlink `.duckdb` files into a worktree. Run the pipeline from the **main repo path**, not the worktree, so writes land in the canonical DBs. We test the worktree's *code* by pointing the main-repo run at the worktree's MLB.R via the working directory:

```bash
cd /Users/callancapitolo/NFLWork
python "Answer Keys/run.py" mlb 2>&1 | tail -40
```

Wait, that runs main's MLB.R. The cleanest verification path: merge to a temp local branch first... no, simpler: just check out the feature branch in the main repo briefly for the test, then switch back. Easier still: copy the modified files into main, test, revert. **Simplest of all:** run the pipeline normally from `main` only AFTER merge (per CLAUDE.md preferred testing path: "or better yet, test from `main` after merging"). For pre-merge verification, do a **read-only** check:

Pre-merge static checks (these are sufficient given the change's scope — single-line `dbdir` swap):
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/mlb-bets-mm
# Sanity: Rscript can parse the modified files
Rscript -e 'parse("Answer Keys/MLB Answer Key/MLB.R")'
Rscript -e 'parse("Answer Keys/MLB Dashboard/mlb_dashboard.R")'
```

Expected: both commands exit 0 with no parse errors.

- [ ] **Step 3: Confirm post-merge integration test plan with user**

Tell the user: "Pre-merge checks pass. The real integration test (run pipeline → confirm `mlb_bets_combined` in `mlb_mm.duckdb` → render dashboard) needs to run from `main` after merge per the CLAUDE.md no-symlink rule. Want me to merge now and run the integration test on main?"

If user says yes:
- Merge (Task 8 below).
- From `main`, run:
  ```bash
  python "Answer Keys/run.py" mlb 2>&1 | tail -40
  ```
- Verify table moved:
  ```bash
  Rscript -e 'library(duckdb); con <- dbConnect(duckdb(), dbdir = "Answer Keys/mlb_mm.duckdb", read_only = TRUE); print("mlb_bets_combined" %in% dbListTables(con)); print(dbGetQuery(con, "SELECT COUNT(*) AS n FROM mlb_bets_combined")); dbDisconnect(con, shutdown = TRUE)'
  Rscript -e 'library(duckdb); con <- dbConnect(duckdb(), dbdir = "Answer Keys/mlb.duckdb", read_only = TRUE); print("mlb_bets_combined" %in% dbListTables(con)); dbDisconnect(con, shutdown = TRUE)'
  ```
  Expected: `TRUE` for `mlb_mm.duckdb`, row count > 0; `FALSE` for `mlb.duckdb`.
- Render dashboard via `Rscript "Answer Keys/MLB Dashboard/mlb_dashboard.R" --parlay-fragment` (quick mode) and confirm exit status 0.

If user prefers pre-merge testing, we'd need to copy files manually into main without committing — fragile and error-prone. Not recommended.

---

### Task 7: Concurrency repro — confirm the bug is actually fixed

**Files:** none — this is the canonical "did we fix the right thing" check.

- [ ] **Step 1: Start a pipeline run (which holds the write lock on `mlb.duckdb`)**

In one terminal, on `main` (post-merge):
```bash
cd /Users/callancapitolo/NFLWork
python "Answer Keys/run.py" mlb 2>&1 | tail -50
```

- [ ] **Step 2: While it's running, render the dashboard**

In a second terminal:
```bash
cd /Users/callancapitolo/NFLWork
Rscript "Answer Keys/MLB Dashboard/mlb_dashboard.R" 2>&1 | tail -30
```

Expected: exit status 0. The "Loading bets from database..." line is followed by either a successful row count (if the previous pipeline run already populated `mlb_mm.duckdb`) or the new graceful warning ("Warning: could not load mlb_bets_combined (...) — rendering with empty bets") if the current run hasn't written yet. **Crucially: no errno 35 crash, no `Execution halted`.**

- [ ] **Step 3: Report results to user**

Paste the relevant output lines into chat.

---

### Task 8: Pre-merge review + merge

- [ ] **Step 1: Run the executive engineer review checklist on the full diff**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/mlb-bets-mm
git diff main..HEAD --stat
git diff main..HEAD
```

Confirm:
- **Data integrity:** writer is DROP+CREATE on every run (idempotent); reader handles missing-table; tryCatch handles missing-DB.
- **Resource safety:** new `con_bets` opens, writes, disconnects in a 3-line block; reader's tryCatch path also disconnects. No leaks.
- **Edge cases:** first dashboard render after deploy but before first pipeline run → empty `all_bets` (already handled by existing `if (!"mlb_bets_combined" %in% tables)` guard); off-season (no rows from MLB.R) → still works.
- **Dead code:** none introduced.
- **Log/disk hygiene:** no new files, no unbounded growth.
- **Security:** no secrets touched.

- [ ] **Step 2: Ask user for explicit merge approval**

Per CLAUDE.md: "Never merge to `main` without explicit user approval, even if tests pass."

- [ ] **Step 3: On approval, merge**

```bash
cd /Users/callancapitolo/NFLWork
git merge --no-ff feature/mlb-bets-combined-mm-migration
```

- [ ] **Step 4: Run Task 6+7 integration tests on `main` and report results to user**

- [ ] **Step 5: On all-green, clean up worktree + branch**

```bash
git worktree remove .worktrees/mlb-bets-mm
git branch -d feature/mlb-bets-combined-mm-migration
```

- [ ] **Step 6: Update memory**

Update `~/.claude-personal/projects/-Users-callancapitolo-NFLWork/memory/mlb_mm_separation.md` to remove the "/refresh path still vulnerable" caveat for `mlb_bets_combined` (the others — `mlb_consensus_temp`, etc. — are still on `mlb.duckdb` but the dashboard doesn't read them).

---

## Summary of Commits

1. `refactor(mlb): write mlb_bets_combined to mlb_mm.duckdb` — Task 1
2. `fix(mlb-dashboard): read mlb_bets_combined from mlb_mm.duckdb + degrade gracefully` — Task 2
3. `docs(mlb): record mlb_bets_combined location move to mlb_mm.duckdb` — Tasks 4 + 5

Plus one merge commit on `main`.
