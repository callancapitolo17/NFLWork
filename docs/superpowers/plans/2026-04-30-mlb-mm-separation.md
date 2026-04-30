# MLB MM Separation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Split MLB pipeline output across two DuckDB files so the RFQ bot's reads and the dashboard's bet-log reads stop colliding with the pipeline's write lock on `mlb.duckdb`. Mirrors the proven CBB pattern (`cbb.duckdb` + `cbb_mm.duckdb`).

**Architecture:** Introduce a new file `Answer Keys/mlb_mm.duckdb` that holds only the six tables consumers actually need. Each writer (MLB.R, mlb_correlated_parlay.R, mlb_triple_play.R, the SGP scrapers) writes those tables to `mlb_mm.duckdb` instead of `mlb.duckdb`. The RFQ bot and the dashboard's parlay/trifecta loaders point at `mlb_mm.duckdb`. Everything else (historical, scraper raw data, dashboard refresh) keeps using `mlb.duckdb`. No new logic, no atomic copy, no snapshot — just a second DuckDB file living alongside the first.

**Tech Stack:** R (`duckdb`), Python (`duckdb` package), DuckDB ≥ 0.10.

---

## Table-to-Consumer Mapping

| Table | Writer | Consumer | Layer |
|---|---|---|---|
| `mlb_game_samples` | `MLB.R:339-340` | RFQ bot (`main.py:101`) | 1 |
| `mlb_samples_meta` | `MLB.R:342-345` | RFQ bot (`main.py:132-135`) | 1 |
| `mlb_parlay_lines` | `mlb_correlated_parlay.R:426-427` | RFQ bot (`main.py:120-129`) | 1 |
| `mlb_sgp_odds` | `mlb_sgp/db.py:125` | RFQ bot + `mlb_correlated_parlay.R` | 1 |
| `mlb_parlay_opportunities` | `mlb_correlated_parlay.R:921-923` | Dashboard (`mlb_dashboard_server.py:1014` + `mlb_dashboard.R:4081,4162`) | 2 |
| `mlb_trifecta_opportunities` | `mlb_triple_play.R:417-418` | Dashboard (`mlb_dashboard_server.py` + `mlb_dashboard.R:65-67`) | 2 |

Tables that **stay in `mlb.duckdb`**: `mlb_team_dict`, `mlb_consensus_temp`, `mlb_odds_temp`, `mlb_bets_combined`, plus all scraper raw-odds tables and historical PBP joins.

## Version Control

- **Branch:** `feature/mlb-mm-separation` (already created)
- **Worktree:** `/Users/callancapitolo/NFLWork/.worktrees/mlb-mm-separation` (already created)
- **Commit style:** one commit per task. Conventional-commit prefixes: `feat(mlb-mm)`, `refactor(mlb-mm)`, `docs(mlb-mm)`, `chore(mlb-mm)`.
- **Merge to main:** only after pre-merge review (Task 9) AND user approval. Per project rules: never merge without explicit OK.
- **Worktree cleanup:** after merge, `git worktree remove .worktrees/mlb-mm-separation && git branch -d feature/mlb-mm-separation`.

## Worktree Lifecycle

- **Setup:** worktree already created at `.worktrees/mlb-mm-separation` on branch `feature/mlb-mm-separation`.
- **DB testing inside worktree:** intentionally minimal. The pipeline reads `pbp.duckdb` (2.5 GB) and depends on live scraper output, so full pipeline runs happen on `main` after merge. In-worktree verification is limited to: (1) R/Python syntax checks; (2) reading the modified files for visual review.
- **Pre-merge:** run all syntax checks, push branch (no PR needed for personal repo), do executive review of `git diff main..HEAD`.
- **After merge:** run pipeline ON MAIN once, verify `mlb_mm.duckdb` has all six expected tables, restart bot with new config, click Place Bet to confirm fix.
- **Cleanup:** remove worktree + branch immediately after merge.

## Documentation

The following docs are touched in **Task 8** as part of the same merge:

- `Answer Keys/CLAUDE.md` — add `mlb_mm.duckdb` to the DuckDB Databases section, mirroring the `cbb_mm.duckdb` line that already exists.
- `kalshi_mlb_rfq/README.md` — update any reference that says the bot reads `mlb.duckdb` to say `mlb_mm.duckdb`.
- `Answer Keys/MLB Dashboard/README.md` — update the data-flow diagram that shows pipeline outputs to include the new file.

---

## Task 1: MLB.R writes bot-facing tables to `mlb_mm.duckdb`

**Why:** `mlb_game_samples` and `mlb_samples_meta` are read only by the bot. Moving them to a separate file means MLB.R's long write window on `mlb.duckdb` no longer matters to the bot.

**Files:**
- Modify: `Answer Keys/MLB Answer Key/MLB.R:339-345`

- [ ] **Step 1.1: Read the current write block**

```r
# Existing code at MLB.R:339-345 (uses con_mlb which targets mlb.duckdb)
dbExecute(con_mlb, "DROP TABLE IF EXISTS mlb_game_samples")
dbWriteTable(con_mlb, "mlb_game_samples", sample_rows)

dbExecute(con_mlb, "DROP TABLE IF EXISTS mlb_samples_meta")
dbExecute(con_mlb, "CREATE TABLE mlb_samples_meta (generated_at TIMESTAMP)")
dbExecute(con_mlb, sprintf(
  "INSERT INTO mlb_samples_meta VALUES ('%s')",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
```

- [ ] **Step 1.2: Replace with a separate-file write**

Open a new connection to `mlb_mm.duckdb`, write the two tables there, close it. Leave `con_mlb` (the `mlb.duckdb` connection) alone for the other writes that stay.

```r
# MLB.R:339-345 — replacement
con_mm <- duckdb_connect_retry("mlb_mm.duckdb")
dbExecute(con_mm, "DROP TABLE IF EXISTS mlb_game_samples")
dbWriteTable(con_mm, "mlb_game_samples", sample_rows)

dbExecute(con_mm, "DROP TABLE IF EXISTS mlb_samples_meta")
dbExecute(con_mm, "CREATE TABLE mlb_samples_meta (generated_at TIMESTAMP)")
dbExecute(con_mm, sprintf(
  "INSERT INTO mlb_samples_meta VALUES ('%s')",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
dbDisconnect(con_mm)
```

- [ ] **Step 1.3: Syntax check**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/mlb-mm-separation/Answer Keys/MLB Answer Key"
Rscript -e 'parse(file = "MLB.R"); cat("OK\n")'
```

Expected: `OK`. (No execution — just a parse check, since full execution needs `pbp.duckdb` and live scraper data.)

- [ ] **Step 1.4: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/mlb-mm-separation"
git add "Answer Keys/MLB Answer Key/MLB.R"
git commit -m "refactor(mlb-mm): MLB.R writes mlb_game_samples + mlb_samples_meta to mlb_mm.duckdb"
```

---

## Task 2: SGP scrapers write `mlb_sgp_odds` to `mlb_mm.duckdb`

**Why:** `mlb_sgp_odds` is the bot's third hot-path table. The four SGP scrapers all use `mlb_sgp/db.py` as their write helper — single point of change.

**Files:**
- Modify: `mlb_sgp/db.py:22`

- [ ] **Step 2.1: Read current MLB_DB definition**

```python
# mlb_sgp/db.py:22
MLB_DB = _REPO_ROOT / "Answer Keys" / "mlb.duckdb"
```

- [ ] **Step 2.2: Replace with mlb_mm.duckdb**

```python
# mlb_sgp/db.py:22 — replacement
MLB_DB = _REPO_ROOT / "Answer Keys" / "mlb_mm.duckdb"
```

That single-line change reroutes all four SGP scrapers (DK, FD, ProphetX, Novig) since they all import `MLB_DB` from this module. The `CREATE TABLE IF NOT EXISTS mlb_sgp_odds` block at line 45 will run on first connect to the new file and create the schema there.

- [ ] **Step 2.3: Syntax check**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/mlb-mm-separation"
python3 -c "import ast; ast.parse(open('mlb_sgp/db.py').read()); print('OK')"
```

Expected: `OK`.

- [ ] **Step 2.4: Commit**

```bash
git add mlb_sgp/db.py
git commit -m "refactor(mlb-mm): mlb_sgp_odds writes go to mlb_mm.duckdb"
```

---

## Task 3: `mlb_correlated_parlay.R` writes `mlb_parlay_lines` and reads `mlb_sgp_odds` from `mlb_mm.duckdb`

**Why:** `mlb_parlay_lines` is a bot-facing table; `mlb_sgp_odds` was just moved in Task 2 — the parlay finder reads it as input. Both need to point at the new file. The script's *output* table (`mlb_parlay_opportunities`) moves in Task 4 (Layer 2).

**Files:**
- Modify: `Answer Keys/mlb_correlated_parlay.R:424-432` (write of `mlb_parlay_lines`)
- Modify: any read of `mlb_sgp_odds` in the same file (need to grep)

- [ ] **Step 3.1: Find all references to mlb.duckdb in the script**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/mlb-mm-separation"
grep -nE 'MLB_DB|mlb\.duckdb|mlb_sgp_odds|mlb_parlay_lines' "Answer Keys/mlb_correlated_parlay.R"
```

Expected: shows the `MLB_DB` constant definition, the `staging_con` block (lines 424-432), the `write_con` block for `mlb_parlay_opportunities` (lines 919-928), and any `dbReadTable`/`dbGetQuery` for `mlb_sgp_odds`.

- [ ] **Step 3.2: Add a parallel `MLB_MM_DB` constant**

Find the existing `MLB_DB` definition (use grep output from Step 3.1 to locate it — it's near the top of the file). Immediately after it, add:

```r
# Bot- and dashboard-facing tables live in a separate DB to avoid lock contention
# with the long write window on mlb.duckdb. Mirrors the CBB pattern.
MLB_MM_DB <- gsub("mlb\\.duckdb$", "mlb_mm.duckdb", MLB_DB)
```

- [ ] **Step 3.3: Update mlb_parlay_lines write to use MLB_MM_DB**

Replace lines 424-432:

```r
# OLD:
staging_con <- dbConnect(duckdb(), dbdir = MLB_DB)
# ...
dbExecute(staging_con, "DROP TABLE IF EXISTS mlb_parlay_lines")
dbWriteTable(staging_con, "mlb_parlay_lines", staging)
# ...
dbDisconnect(staging_con)

# NEW:
staging_con <- dbConnect(duckdb(), dbdir = MLB_MM_DB)
# ...
dbExecute(staging_con, "DROP TABLE IF EXISTS mlb_parlay_lines")
dbWriteTable(staging_con, "mlb_parlay_lines", staging)
# ...
dbDisconnect(staging_con)
```

(Only the `dbConnect(... dbdir = ...)` line changes; the rest is the same.)

- [ ] **Step 3.4: Update any read of `mlb_sgp_odds` in this script to use MLB_MM_DB**

Use the grep output from Step 3.1 to locate each `mlb_sgp_odds` reference. For each `dbConnect` opened to read `mlb_sgp_odds`, change `dbdir = MLB_DB` → `dbdir = MLB_MM_DB`. If the same connection reads tables that are still in `mlb.duckdb`, split it: open one connection per file. (DuckDB connections are cheap.)

- [ ] **Step 3.5: Syntax check**

```bash
Rscript -e 'parse(file = "Answer Keys/mlb_correlated_parlay.R"); cat("OK\n")'
```

Expected: `OK`.

- [ ] **Step 3.6: Commit**

```bash
git add "Answer Keys/mlb_correlated_parlay.R"
git commit -m "refactor(mlb-mm): correlated_parlay reads sgp_odds + writes parlay_lines via mlb_mm.duckdb"
```

---

## Task 4: `mlb_correlated_parlay.R` writes `mlb_parlay_opportunities` to `mlb_mm.duckdb` (Layer 2)

**Why:** This is the dashboard's bet-log hot-path table. Moving it to `mlb_mm.duckdb` is what fixes the 500 you saw at 21:18.

**Files:**
- Modify: `Answer Keys/mlb_correlated_parlay.R:919-928`

- [ ] **Step 4.1: Update the write_con destination**

```r
# OLD (lines 921-922):
write_con <- dbConnect(duckdb(), dbdir = MLB_DB)
dbExecute(write_con, "DROP TABLE IF EXISTS mlb_parlay_opportunities")

# NEW:
write_con <- dbConnect(duckdb(), dbdir = MLB_MM_DB)
dbExecute(write_con, "DROP TABLE IF EXISTS mlb_parlay_opportunities")
```

- [ ] **Step 4.2: Syntax check**

```bash
Rscript -e 'parse(file = "Answer Keys/mlb_correlated_parlay.R"); cat("OK\n")'
```

- [ ] **Step 4.3: Commit**

```bash
git add "Answer Keys/mlb_correlated_parlay.R"
git commit -m "refactor(mlb-mm): mlb_parlay_opportunities writes to mlb_mm.duckdb"
```

---

## Task 5: `mlb_triple_play.R` writes `mlb_trifecta_opportunities` to `mlb_mm.duckdb` (Layer 2)

**Why:** Same as Task 4, for the trifecta tab.

**Files:**
- Modify: `Answer Keys/mlb_triple_play.R:87-100` (cleanup block) + `:417-418` (main write)

- [ ] **Step 5.1: Find existing MLB_DB references**

```bash
grep -nE 'MLB_DB|mlb\.duckdb|mlb_trifecta_opportunities' "Answer Keys/mlb_triple_play.R"
```

- [ ] **Step 5.2: Add MLB_MM_DB constant**

Same pattern as Task 3, immediately after the existing `MLB_DB` definition:

```r
MLB_MM_DB <- gsub("mlb\\.duckdb$", "mlb_mm.duckdb", MLB_DB)
```

- [ ] **Step 5.3: Update the cleanup block (line 87-100)**

The cleanup block opens a connection to drop the table when no triple-play candidates are found. Change its `dbConnect(..., dbdir = MLB_DB)` to `dbdir = MLB_MM_DB`.

- [ ] **Step 5.4: Update the main write block (line 417-418)**

Same change: `dbConnect(..., dbdir = MLB_DB)` → `dbdir = MLB_MM_DB` for the connection that wraps these two lines:

```r
dbExecute(write_con, "DROP TABLE IF EXISTS mlb_trifecta_opportunities")
dbWriteTable(write_con, "mlb_trifecta_opportunities", priced)
```

- [ ] **Step 5.5: Syntax check**

```bash
Rscript -e 'parse(file = "Answer Keys/mlb_triple_play.R"); cat("OK\n")'
```

- [ ] **Step 5.6: Commit**

```bash
git add "Answer Keys/mlb_triple_play.R"
git commit -m "refactor(mlb-mm): mlb_trifecta_opportunities writes to mlb_mm.duckdb"
```

---

## Task 6: RFQ bot reads from `mlb_mm.duckdb`

**Why:** Now that the writers all target the new file, point the consumer at it. One-line config change.

**Files:**
- Modify: `kalshi_mlb_rfq/config.py:82`

- [ ] **Step 6.1: Update ANSWER_KEY_DB**

```python
# kalshi_mlb_rfq/config.py:82 — OLD
ANSWER_KEY_DB = PROJECT_ROOT / "Answer Keys" / "mlb.duckdb"

# NEW
ANSWER_KEY_DB = PROJECT_ROOT / "Answer Keys" / "mlb_mm.duckdb"
```

- [ ] **Step 6.2: Syntax check**

```bash
python3 -c "import ast; ast.parse(open('kalshi_mlb_rfq/config.py').read()); print('OK')"
```

- [ ] **Step 6.3: Commit**

```bash
git add kalshi_mlb_rfq/config.py
git commit -m "refactor(mlb-mm): RFQ bot reads from mlb_mm.duckdb"
```

---

## Task 7: Dashboard reads parlay/trifecta opportunities from `mlb_mm.duckdb`

**Why:** Layer 2 reader side — completes the bet-log fix.

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard_server.py` (add `MLB_MM_DB` constant, update `_load_parlay_row` + the trifecta loader)
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R:65-67, 4081, 4161-4163`

- [ ] **Step 7.1: Add `MLB_MM_DB` constant in mlb_dashboard_server.py**

Find the existing `MLB_DB` definition at line 51:

```python
MLB_DB = REPO_ROOT / "Answer Keys" / "mlb.duckdb"
```

Add immediately after:

```python
# Bot- and dashboard-bet-log-facing tables live in mlb_mm.duckdb so they're
# never blocked by the pipeline's long write lock on mlb.duckdb.
MLB_MM_DB = REPO_ROOT / "Answer Keys" / "mlb_mm.duckdb"
```

- [ ] **Step 7.2: Update `_load_parlay_row` (line 1014)**

```python
# OLD:
con = duckdb.connect(str(MLB_DB), read_only=True)

# NEW:
con = duckdb.connect(str(MLB_MM_DB), read_only=True)
```

- [ ] **Step 7.3: Locate and update the trifecta row loader**

Find the trifecta equivalent of `_load_parlay_row`:

```bash
grep -n "mlb_trifecta_opportunities\|_load_trifecta\|trifecta_hash" "Answer Keys/MLB Dashboard/mlb_dashboard_server.py"
```

For every site that opens `MLB_DB` to read `mlb_trifecta_opportunities`, change to `MLB_MM_DB`. Leave reads of any other table (e.g. `mlb_bets_combined`) on `MLB_DB`.

- [ ] **Step 7.4: Update mlb_dashboard.R reads of these two tables**

Find every read of `mlb_parlay_opportunities` and `mlb_trifecta_opportunities`:

```bash
grep -nE 'mlb_parlay_opportunities|mlb_trifecta_opportunities' "Answer Keys/MLB Dashboard/mlb_dashboard.R"
```

For each, ensure the surrounding `dbConnect` opens `mlb_mm.duckdb`. Where the connection is shared with reads of *other* tables (e.g. `mlb_bets_combined` reads stay on `mlb.duckdb`), split the connection so each table is read from the file that holds it.

Specifically:
- Line 65-67: trifectas tab — open a connection to `mlb_mm.duckdb` here
- Line 4081: switch the `mlb_parlay_opportunities` query to a connection on `mlb_mm.duckdb`
- Line 4161-4163: same as 4081

- [ ] **Step 7.5: Syntax check both files**

```bash
python3 -c "import ast; ast.parse(open('Answer Keys/MLB Dashboard/mlb_dashboard_server.py').read()); print('py OK')"
Rscript -e 'parse(file = "Answer Keys/MLB Dashboard/mlb_dashboard.R"); cat("R OK\n")'
```

- [ ] **Step 7.6: Commit**

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard_server.py" "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "refactor(mlb-mm): dashboard reads parlay/trifecta opps from mlb_mm.duckdb"
```

---

## Task 8: Update documentation

**Why:** Project rule: every architectural change updates the docs in the same merge.

**Files:**
- Modify: `Answer Keys/CLAUDE.md`
- Modify: `kalshi_mlb_rfq/README.md`
- Modify: `Answer Keys/MLB Dashboard/README.md`

- [ ] **Step 8.1: Update Answer Keys/CLAUDE.md DuckDB Databases section**

Find the existing list (search for `### DuckDB Databases`). Add a new line, mirroring the `cbb_mm.duckdb` line that's already there:

```markdown
- `mlb_mm.duckdb` — MLB MM consumer tables (samples, sgp_odds, parlay_lines, samples_meta, parlay_opportunities, trifecta_opportunities). Separate from `mlb.duckdb` so the pipeline's write lock on `mlb.duckdb` never blocks the RFQ bot or the dashboard's bet-log loaders.
```

- [ ] **Step 8.2: Update kalshi_mlb_rfq/README.md**

Find any reference to `mlb.duckdb` and update to `mlb_mm.duckdb`. Add a sentence explaining the separation:

> The bot reads `Answer Keys/mlb_mm.duckdb` (read-only). This file holds the six consumer-facing tables and is written by `MLB.R`, `mlb_correlated_parlay.R`, `mlb_triple_play.R`, and the SGP scrapers — separate from `mlb.duckdb` (the pipeline's main write target) so lock contention can't block the bot's cache refreshes.

- [ ] **Step 8.3: Update Answer Keys/MLB Dashboard/README.md**

Locate the data-flow diagram showing `run.py mlb (pipeline) → mlb.duckdb/mlb_bets_combined`. Add a parallel arrow:

```
run.py mlb (pipeline)  →  mlb.duckdb/mlb_bets_combined           (dashboard refresh)
                       →  mlb_mm.duckdb/mlb_parlay_opportunities  (place-bet loader)
                       →  mlb_mm.duckdb/mlb_trifecta_opportunities
```

- [ ] **Step 8.4: Commit**

```bash
git add "Answer Keys/CLAUDE.md" kalshi_mlb_rfq/README.md "Answer Keys/MLB Dashboard/README.md"
git commit -m "docs(mlb-mm): document mlb_mm.duckdb consumer-isolation pattern"
```

---

## Task 9: Pre-merge executive review

**Why:** Project rule — every feature branch gets reviewed against the executive-engineer checklist before merging to main.

- [ ] **Step 9.1: Diff review**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/mlb-mm-separation"
git diff main..HEAD --stat
git diff main..HEAD
```

- [ ] **Step 9.2: Run the project checklist**

For each item, confirm or note as ACCEPTABLE RISK:

- **Data integrity:** No table is silently lost — every write of the moved tables now lands in `mlb_mm.duckdb`. The OLD copies still in `mlb.duckdb` will go stale. (See Task 10 cleanup.)
- **Resource safety:** Every new `dbConnect` is paired with `dbDisconnect` on the same code path. No leaked connections.
- **Edge cases:**
  - First-run / `mlb_mm.duckdb` doesn't exist yet → `duckdb_connect_retry` creates it. ✓
  - Bot starts before pipeline runs once → connect succeeds, table reads return empty/error → existing retry logic handles. ✓
  - Pipeline runs concurrently from /refresh and bot's `_run_pipeline` → both write to `mlb_mm.duckdb` → DuckDB lock + retry → same behavior as before. ✓
- **Dead code:** No new functions or branches that aren't used.
- **Log/disk hygiene:** Adds one new ~1MB DuckDB file. Negligible.
- **Security:** No secrets touched.

- [ ] **Step 9.3: Document findings**

In the conversation (not a file), produce: ISSUES TO FIX vs ACCEPTABLE RISKS list. Get explicit user approval before proceeding to Task 10.

---

## Task 10: Merge to main + live cutover

**Why:** All worktree tests are syntax-only. Real verification requires running the pipeline against live data on main. This task is the operational handover.

**Approval gate:** the user must explicitly OK the merge before any of these steps runs.

- [ ] **Step 10.1: Stop the RFQ bot**

```bash
ps aux | grep "kalshi_mlb_rfq.main" | grep -v grep   # find PID
kill <PID>
```

The bot does not auto-restart (no launchd entry for it). Note: per memory, `Bash(kill:*)` is allowlisted but flag the bot kill in chat first.

- [ ] **Step 10.2: Merge worktree branch into main**

```bash
cd "/Users/callancapitolo/NFLWork"
git checkout main
git merge --no-ff feature/mlb-mm-separation -m "feat(mlb-mm): split MLB pipeline output into mlb.duckdb + mlb_mm.duckdb"
```

- [ ] **Step 10.3: Run pipeline once to populate mlb_mm.duckdb**

```bash
cd "/Users/callancapitolo/NFLWork/Answer Keys"
python3 run.py mlb 2>&1 | tail -40
Rscript mlb_correlated_parlay.R 2>&1 | tail -20
Rscript mlb_triple_play.R 2>&1 | tail -20
```

Expected: each exits 0; `mlb_mm.duckdb` now exists.

- [ ] **Step 10.4: Verify mlb_mm.duckdb has all six tables**

```bash
python3 -c "
import duckdb
con = duckdb.connect('Answer Keys/mlb_mm.duckdb', read_only=True)
print(con.execute('SHOW TABLES').fetchall())
for t in ['mlb_game_samples','mlb_samples_meta','mlb_parlay_lines','mlb_sgp_odds','mlb_parlay_opportunities','mlb_trifecta_opportunities']:
    n = con.execute(f'SELECT COUNT(*) FROM {t}').fetchone()[0]
    print(f'{t}: {n} rows')
"
```

Expected: all six tables present, non-zero row counts (except sgp_odds and trifecta_opportunities, which can be 0 if no qualifying games).

- [ ] **Step 10.5: Restart RFQ bot**

```bash
cd "/Users/callancapitolo/NFLWork"
PYTHONPATH=. ./kalshi_mlb_rfq/venv/bin/python -u -m kalshi_mlb_rfq.main >> kalshi_mlb_rfq/bot.log 2>&1 &
```

- [ ] **Step 10.6: Verify bot cache_refresh succeeds**

```bash
tail -20 "/Users/callancapitolo/NFLWork/kalshi_mlb_rfq/bot.log"
```

Look for a line like:

```
cache_refresh: 11 games, 128 sgp_odds rows, 11 parlay_lines, samples gen_at=...
```

- [ ] **Step 10.7: Restart dashboard server (so it picks up MLB_MM_DB)**

The dashboard is a Flask app; restart however it's normally started.

- [ ] **Step 10.8: Manual smoke test — Place a parlay**

Open the dashboard, click Place Bet on a parlay, confirm: no 500, ticket recorded normally.

- [ ] **Step 10.9: Worktree + branch cleanup**

```bash
cd "/Users/callancapitolo/NFLWork"
git worktree remove .worktrees/mlb-mm-separation
git branch -d feature/mlb-mm-separation
```

---

## Task 11 (deferred — separate cleanup PR): Drop stale tables from `mlb.duckdb`

**Why:** After cutover, the moved tables still exist as stale copies in `mlb.duckdb`. They're not read by anything but they're confusing. A one-time `DROP TABLE` cleans them up.

**Defer this** until layers 1+2 have been running cleanly for a day. No code changes — just a one-shot SQL cleanup.

```sql
-- Run once after layers 1+2 are stable:
DROP TABLE IF EXISTS mlb_game_samples;
DROP TABLE IF EXISTS mlb_samples_meta;
DROP TABLE IF EXISTS mlb_parlay_lines;
DROP TABLE IF EXISTS mlb_sgp_odds;
DROP TABLE IF EXISTS mlb_parlay_opportunities;
DROP TABLE IF EXISTS mlb_trifecta_opportunities;
```

Not part of this merge.

---

## Pre-mortem: vulnerabilities considered

1. **Bootstrap race.** First MLB.R run after merge creates `mlb_mm.duckdb`. ✓ (`duckdb_connect_retry` creates the file on connect.)
2. **Bot reads stale tables.** During the migration window, if the bot is reading from `mlb.duckdb` before the config is flipped, it reads the LAST-good data the writer left there. → Solved by the explicit `stop bot → merge → populate → flip config (already in commit) → restart` ordering in Task 10.
3. **`mlb_correlated_parlay.R` reads `mlb_sgp_odds`.** The script needs to read `mlb_sgp_odds` from the new file. → Handled in Task 3 step 3.4.
4. **Dashboard `/refresh` chain has multiple R scripts that touch the moved tables.** All such writes/reads are explicitly enumerated in Tasks 3-7. ✓
5. **`/refresh` still spawns `MLB.R`** which still write-locks `mlb.duckdb` for ~30-90s. The dashboard's `mlb_dashboard.R` still reads other tables from `mlb.duckdb` (`mlb_bets_combined`, etc.) and would still 500 if a bot-spawned MLB.R is mid-write at refresh time. → OUT OF SCOPE for layers 1+2. Tracked separately as "/refresh retry-on-lock wrapper" or "stop dashboard from running its own MLB.R."
6. **Two parallel writers.** `mlb_correlated_parlay.R` and `mlb_triple_play.R` are launched in parallel by `/refresh`. After Task 4+5 they both write to `mlb_mm.duckdb`. They write *different tables*. DuckDB serializes via the file lock → tryCatch retry on each → same behavior as before with `mlb.duckdb`. ✓

## Self-review checklist

- ✅ Spec coverage: each table in the mapping has an explicit task that moves it.
- ✅ No placeholders: every R/Python edit shows the exact code.
- ✅ Type consistency: `MLB_MM_DB` constant pattern used identically across `mlb_correlated_parlay.R`, `mlb_triple_play.R`, and `mlb_dashboard_server.py`.
- ✅ Worktree lifecycle, version control, documentation, pre-merge review all explicit.
- ✅ Operational cutover (Task 10) covers stop-bot → merge → populate → restart.
