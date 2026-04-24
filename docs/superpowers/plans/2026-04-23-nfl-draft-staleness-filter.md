# NFL Draft Dashboard — 20-Minute Staleness Filter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the 2-hour `MAX_AGE_HOURS` staleness filter in the NFL Draft dashboard with a 20-minute `MAX_AGE_MINUTES` flat threshold, so dead venues are hidden within one pre-draft scrape cycle instead of lingering for hours.

**Architecture:** One module-level constant in `nfl_draft/lib/queries.py` gates two SQL queries (`cross_book_grid`, `kalshi_tooltip_data`). The change is a rename + unit-swap (hours → minutes) with a new value (20). Two integration tests import the constant by name, so they must be updated in lockstep. README documents the threshold for future readers.

**Tech Stack:** Python 3, DuckDB, pytest, Dash (unchanged — only SQL query functions are touched).

**Spec:** `docs/superpowers/specs/2026-04-23-nfl-draft-staleness-filter-design.md`

---

## Version Control

- **Feature branch:** `feature/staleness-20min`, branched from `main`.
- **Worktree:** `/Users/callancapitolo/NFLWork/.worktrees/staleness-20min` (follows existing convention under `.worktrees/`).
- **DuckDB safety:** This change only touches `queries.py` and tests (tests use `tmp_path` fixtures, no live DB access). No symlinks, no shared-DB concerns.
- **Commits:** One commit with source + tests + README together. The spec allows same-commit or immediately-after; single-commit keeps the merge atomic.
- **Merge gate:** Pre-merge review (`git diff main..HEAD`) + explicit user approval before fast-forward merge to `main`.
- **Cleanup:** Remove worktree and delete feature branch immediately after merge.

---

## File Structure

Files touched in this plan:

- **Modify:** `nfl_draft/lib/queries.py` — rename constant, rewrite comment block, update two SQL interval clauses.
- **Modify:** `nfl_draft/tests/integration/test_dashboard_queries.py` — import rename, timedelta unit-swap, two docstring updates.
- **Modify:** `nfl_draft/README.md` — state the current 20-min threshold in the existing staleness-filter sentence.

No new files.

---

## Documentation Updates

- `nfl_draft/README.md` — add the threshold value to the existing FanDuel/staleness sentence (Task 5).
- Spec doc (`docs/superpowers/specs/2026-04-23-nfl-draft-staleness-filter-design.md`) — already committed on `main` at `a15f32c`, no further changes needed.
- No `CLAUDE.md` updates needed — the constant and its reasoning live in `queries.py` comments.

---

## Task 0: Create Worktree and Feature Branch

**Files:** None (setup only).

- [ ] **Step 1: Confirm you are on `main` at the repo root**

```bash
cd /Users/callancapitolo/NFLWork
git branch --show-current
```

Expected: `main`

- [ ] **Step 2: Create the worktree on a new feature branch**

```bash
git worktree add /Users/callancapitolo/NFLWork/.worktrees/staleness-20min -b feature/staleness-20min main
```

Expected: `Preparing worktree (new branch 'feature/staleness-20min')` then `HEAD is now at <sha>`

- [ ] **Step 3: Switch into the worktree for all remaining tasks**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/staleness-20min
git branch --show-current
```

Expected: `feature/staleness-20min`. All remaining file paths below are **relative to this worktree root** unless otherwise noted.

---

## Task 1: Update the Stale-Snapshot Test to Use `MAX_AGE_MINUTES`

This is TDD step one: rewrite the test that imports the constant, run it, and watch it fail because the constant name no longer matches. That failure is the signal that drives the code change in Task 2.

**Files:**
- Modify: `nfl_draft/tests/integration/test_dashboard_queries.py:757-791` (the `test_kalshi_tooltip_data_excludes_stale_snapshots` function, specifically the import on line 766 and the `timedelta` line on 767).

- [ ] **Step 1: Open the test file and locate the import + timedelta lines**

Open `nfl_draft/tests/integration/test_dashboard_queries.py`. Lines 766–767 currently read:

```python
    from nfl_draft.lib.queries import MAX_AGE_HOURS
    stale = datetime.now() - timedelta(hours=MAX_AGE_HOURS + 1)
```

- [ ] **Step 2: Replace those two lines with the minute-scoped version**

Change them to:

```python
    from nfl_draft.lib.queries import MAX_AGE_MINUTES
    stale = datetime.now() - timedelta(minutes=MAX_AGE_MINUTES + 5)
```

The `+ 5` cushion keeps the row comfortably past the cutoff (25 min when threshold is 20), mirroring the intent of the old `+ 1` hour.

- [ ] **Step 3: Run the test and confirm it fails with an ImportError**

```bash
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest \
  nfl_draft/tests/integration/test_dashboard_queries.py::test_kalshi_tooltip_data_excludes_stale_snapshots -v
```

Expected: FAIL. Specifically an `ImportError: cannot import name 'MAX_AGE_MINUTES' from 'nfl_draft.lib.queries'` (because the code change hasn't happened yet).

If the test passes or fails for any other reason, stop and diagnose before continuing.

- [ ] **Step 4: Do NOT commit yet**

We want the source + test changes in a single commit. Hold.

---

## Task 2: Rename the Constant and Update SQL in `queries.py`

TDD step two: make the failing test pass by introducing `MAX_AGE_MINUTES = 20` and swapping the SQL intervals from hours to minutes.

**Files:**
- Modify: `nfl_draft/lib/queries.py:36-42` (constant + comment block), `nfl_draft/lib/queries.py:121` (cross_book_grid SQL), `nfl_draft/lib/queries.py:371` (kalshi_tooltip_data SQL).

- [ ] **Step 1: Replace the comment block and constant**

Open `nfl_draft/lib/queries.py`. Lines 36–42 currently read:

```python
# Any row older than MAX_AGE_HOURS is excluded from Cross-Book Grid /
# ev_candidates. This prevents venues that have silently stopped scraping
# (e.g. FD regression 2026-04-19) from polluting the grid with 24h-old
# prices that get flagged as "edges" when in reality the feed just died.
# Tune via this constant — 2h is a conservative default for pre-draft
# cadence (scrapers run every few minutes).
MAX_AGE_HOURS = 2
```

Replace with:

```python
# Any row older than MAX_AGE_MINUTES is excluded from Cross-Book Grid,
# +EV Candidates, and the Kalshi tooltip. This prevents venues that have
# silently stopped scraping (e.g. FD regression 2026-04-19) from polluting
# the grid with hours-old prices that get flagged as "edges" when the feed
# has actually died. 20 minutes = pre-draft scrape cadence (15 min from
# crontab.pre) + 5-minute cushion for a late cron run; it is also ~10x the
# draft-day cadence (2 min from crontab.draft), so dead venues drop out
# within ~20 min on draft day instead of the prior ~2h.
MAX_AGE_MINUTES = 20
```

- [ ] **Step 2: Update the `cross_book_grid` SQL interval**

At `nfl_draft/lib/queries.py:121`, change:

```python
                  WHERE fetched_at > NOW() - INTERVAL '{MAX_AGE_HOURS} hours'
```

to:

```python
                  WHERE fetched_at > NOW() - INTERVAL '{MAX_AGE_MINUTES} minutes'
```

Then update the enclosing f-string interpolation so the new constant name is actually passed in. Search for the `.format(` or f-string preamble of this query block (it lives inside `cross_book_grid`, lines 89–160ish). Confirm the SQL is built with the `MAX_AGE_HOURS` name and swap that reference to `MAX_AGE_MINUTES`.

**How to verify:** after the edit, `grep -n "MAX_AGE_HOURS" nfl_draft/lib/queries.py` should return **no matches**.

- [ ] **Step 3: Update the `kalshi_tooltip_data` SQL interval**

At `nfl_draft/lib/queries.py:371`, change:

```python
                  WHERE fetch_time > NOW() - INTERVAL '{MAX_AGE_HOURS} hours'
```

to:

```python
                  WHERE fetch_time > NOW() - INTERVAL '{MAX_AGE_MINUTES} minutes'
```

Again confirm the format / f-string feeder references `MAX_AGE_MINUTES`.

- [ ] **Step 4: Re-run the stale-snapshot test and verify it now PASSES**

```bash
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest \
  nfl_draft/tests/integration/test_dashboard_queries.py::test_kalshi_tooltip_data_excludes_stale_snapshots -v
```

Expected: PASS.

- [ ] **Step 5: Confirm `MAX_AGE_HOURS` is fully gone from the `nfl_draft/` tree**

```bash
grep -rn "MAX_AGE_HOURS" nfl_draft/
```

Expected: no matches. If any remain (including in tests), fix them — they will break the next test run.

- [ ] **Step 6: Do NOT commit yet**

Docstring touch-ups + README still ahead.

---

## Task 3: Update the Two `MAX_AGE_HOURS` Docstring References

Both docstrings mention the old constant name. Update them so a future reader searching for `MAX_AGE_MINUTES` finds them.

**Files:**
- Modify: `nfl_draft/tests/integration/test_dashboard_queries.py:389`, `nfl_draft/tests/integration/test_dashboard_queries.py:758`.

- [ ] **Step 1: Update docstring at line 389**

Currently:

```python
    """Rows older than MAX_AGE_HOURS must be filtered out of cross_book_grid.
```

Change to:

```python
    """Rows older than MAX_AGE_MINUTES must be filtered out of cross_book_grid.
```

- [ ] **Step 2: Update docstring at line 758**

Currently:

```python
    """Snapshots older than MAX_AGE_HOURS must be filtered out so the tooltip
    never shows a price next to a cell whose underlying scrape is stale."""
```

Change to:

```python
    """Snapshots older than MAX_AGE_MINUTES must be filtered out so the tooltip
    never shows a price next to a cell whose underlying scrape is stale."""
```

- [ ] **Step 3: Grep once more to confirm zero references to the old name**

```bash
grep -rn "MAX_AGE_HOURS" nfl_draft/
```

Expected: no matches.

---

## Task 4: Run the Full Dashboard Query Test File

Acceptance criterion #3: `pytest nfl_draft/tests/integration/test_dashboard_queries.py` passes. Run it and confirm no regressions — the 3-day-old row in `test_cross_book_grid_excludes_stale_rows` is still stale under 20 minutes, so that test should keep passing without modification.

**Files:** None (verification step).

- [ ] **Step 1: Run the whole file**

```bash
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest \
  nfl_draft/tests/integration/test_dashboard_queries.py -v
```

Expected: ALL PASS. Pay particular attention to:
- `test_cross_book_grid_excludes_stale_rows` — PASS (3-day-old row is still stale under 20-min threshold)
- `test_kalshi_tooltip_data_excludes_stale_snapshots` — PASS (now that the constant exists)

- [ ] **Step 2: If anything fails, stop and diagnose**

Do not proceed to the README / commit steps until the full file is green. Common failure modes to look for:
- An `f"..."` vs `.format(...)` mismatch where the SQL still tries to interpolate `MAX_AGE_HOURS`.
- A missed reference somewhere in the file (Task 2 step 5's grep should have caught this).

---

## Task 5: Update the README's Staleness-Filter Sentence

**Files:**
- Modify: `nfl_draft/README.md` — the existing sentence about "the dashboard's staleness filter hides FD rows while they're gone" (appears in the "Current venue status" block near the top).

- [ ] **Step 1: Read the current line**

Open `nfl_draft/README.md` and find the sentence. It currently reads approximately:

> the dashboard's staleness filter hides FD rows while they're gone

- [ ] **Step 2: Append the threshold value**

Replace that clause with:

> the dashboard's 20-minute staleness filter hides FD rows while they're gone (see `MAX_AGE_MINUTES` in `nfl_draft/lib/queries.py`)

Keep the surrounding sentence structure intact — the goal is only to document the current threshold and point the reader to the source of truth. Do not rewrite other parts of the README.

- [ ] **Step 3: Verify the change reads cleanly**

Skim the surrounding paragraph to make sure the added clause still flows in context. No other README sections need edits.

---

## Task 6: Smoke-Test the Dashboard Against the Live DuckDB

Acceptance criterion #5: dashboard boots cleanly and the affected tabs render. The worktree shares no DuckDB with `main` (we're only reading, not writing), but we'll read against the live `nfl_draft/nfl_draft.duckdb` by pointing the dashboard at the absolute path.

**Files:** None (smoke test).

- [ ] **Step 1: Start the dashboard from the worktree**

From inside the worktree (`/Users/callancapitolo/NFLWork/.worktrees/staleness-20min`):

```bash
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m kalshi_draft.app
```

(or whatever the existing `run.sh` invocation is in the repo — check `kalshi_draft/run.sh` if unsure).

Expected: Dash prints `Dash is running on http://127.0.0.1:8090/` (or similar port). No stack traces.

- [ ] **Step 2: Load the four portal tabs in a browser**

Visit `http://127.0.0.1:8090/` and click through:
- Cross-Book Grid — grid renders, no errors in the terminal.
- +EV Candidates — candidates list renders.
- Trade Tape — renders.
- Bet Log — renders.

- [ ] **Step 3: Hover a Kalshi cell in Cross-Book Grid**

Tooltip should appear with buy/sell/last. If the tooltip is empty but the cell has a value, something is off with the `kalshi_tooltip_data` SQL — stop and diagnose.

- [ ] **Step 4: Stop the dashboard**

Ctrl-C the process in the terminal.

---

## Task 7: Commit Source + Tests + README in a Single Commit

**Files:** All three modified files.

- [ ] **Step 1: Sanity-check the diff**

```bash
git diff --stat main..HEAD
git diff main..HEAD
```

Expected: exactly three files changed — `nfl_draft/lib/queries.py`, `nfl_draft/tests/integration/test_dashboard_queries.py`, `nfl_draft/README.md`. No stray edits.

- [ ] **Step 2: Stage the three files explicitly**

```bash
git add nfl_draft/lib/queries.py \
        nfl_draft/tests/integration/test_dashboard_queries.py \
        nfl_draft/README.md
```

Do not use `git add -A` — that could pick up unrelated untracked files.

- [ ] **Step 3: Create the commit**

```bash
git commit -m "$(cat <<'EOF'
feat(nfl_draft): tighten dashboard staleness filter to 20 minutes

Replaces the 2-hour MAX_AGE_HOURS filter in queries.py with a 20-minute
MAX_AGE_MINUTES flat threshold. Dead venues (e.g. FD on 2026-04-19) now
drop off Cross-Book Grid / +EV Candidates / Kalshi tooltip within one
pre-draft scrape cycle instead of lingering for hours.

Covers two integration tests that import the constant by name, plus a
README note pointing readers at the new constant.

Spec: docs/superpowers/specs/2026-04-23-nfl-draft-staleness-filter-design.md

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

- [ ] **Step 4: Confirm commit landed**

```bash
git log --oneline -1
```

Expected: one new commit on `feature/staleness-20min` with the message above.

---

## Task 8: Pre-Merge Executive Review

**Files:** None (review checklist).

- [ ] **Step 1: Print the full diff against `main` and read it end-to-end**

```bash
git diff main..HEAD
```

- [ ] **Step 2: Walk the checklist**

For each item, confirm PASS or flag it:

1. **Data integrity:** no writes are involved; the change is read-only. PASS by construction.
2. **Resource safety:** no new DB connections opened; existing `read_connection` context manager handles cleanup. PASS.
3. **Edge cases:**
   - Empty DB (first run, no rows) → SQL still returns nothing, grid renders empty. PASS.
   - Off-season (scrapers not running) → all rows older than 20 min, grid renders empty. This is the desired behavior, documented in the spec's Known Tradeoff section.
   - Timezone: `NOW()` in DuckDB and `datetime.now()` in tests both use the same local-time reference. PASS.
4. **Dead code:** no new flags, functions, or imports introduced — only rename. PASS.
5. **Log/disk hygiene:** no logging changes. PASS.
6. **Security:** no secrets touched. PASS.

- [ ] **Step 3: Present findings to the user**

Summarize the review in chat: "Pre-merge review complete, X items checked, Y issues found (or none). Requesting approval to merge."

- [ ] **Step 4: Wait for explicit user approval before merging**

Do NOT run `git merge` or any push command until the user says "go."

---

## Task 9: Fast-Forward Merge to `main`

**Files:** None (merge step). **Only proceed after user approval in Task 8 Step 4.**

- [ ] **Step 1: Return to the main working tree**

```bash
cd /Users/callancapitolo/NFLWork
git branch --show-current
```

Expected: `main`.

- [ ] **Step 2: Pull latest (defensive, in case anything was pushed in the meantime)**

```bash
git fetch origin
git status
```

Expected: "Your branch is up to date with 'origin/main'" (or a warning that you are ahead of origin, which is fine).

- [ ] **Step 3: Fast-forward merge**

```bash
git merge --ff-only feature/staleness-20min
```

Expected: `Fast-forward` followed by the file summary. If this fails with "not possible to fast-forward," stop and investigate — it means `main` has diverged.

- [ ] **Step 4: Confirm the merge landed on `main`**

```bash
git log --oneline -3
grep -n "MAX_AGE_MINUTES" nfl_draft/lib/queries.py
```

Expected: feat commit visible; grep returns at least three matches (constant + two SQL references).

---

## Task 10: Clean Up the Worktree and Feature Branch

**Files:** None (cleanup step).

- [ ] **Step 1: Remove the worktree**

```bash
cd /Users/callancapitolo/NFLWork
git worktree remove /Users/callancapitolo/NFLWork/.worktrees/staleness-20min
```

Expected: no output; `ls .worktrees/` shows the directory is gone.

- [ ] **Step 2: Delete the feature branch**

```bash
git branch -d feature/staleness-20min
```

Expected: `Deleted branch feature/staleness-20min (was <sha>).`

- [ ] **Step 3: Verify cleanup**

```bash
git worktree list
git branch --list feature/staleness-20min
```

Expected: no worktree entry for `staleness-20min`; branch list empty for that name.

---

## Acceptance-Criteria Verification (post-merge)

Run these commands from `/Users/callancapitolo/NFLWork` on `main` after the merge lands. Each maps to one spec acceptance criterion:

1. **Constant renamed, old name gone:**
   ```bash
   grep -rn "MAX_AGE_HOURS" nfl_draft/
   grep -n "MAX_AGE_MINUTES" nfl_draft/lib/queries.py
   ```
   Expected: first command returns nothing; second returns `MAX_AGE_MINUTES = 20`.

2. **SQL uses minute intervals:**
   ```bash
   grep -n "INTERVAL" nfl_draft/lib/queries.py
   ```
   Expected: two `INTERVAL '{MAX_AGE_MINUTES} minutes'` clauses.

3. **Tests green:**
   ```bash
   /Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest \
     nfl_draft/tests/integration/test_dashboard_queries.py -v
   ```
   Expected: all PASS.

4. **README mentions threshold:**
   ```bash
   grep -n "20-minute\|MAX_AGE_MINUTES" nfl_draft/README.md
   ```
   Expected: at least one match.

5. **Dashboard renders:** already confirmed in Task 6. No need to re-run unless you suspect a merge anomaly.

---

## Self-Review Notes (for plan author)

- **Spec coverage:** every numbered acceptance criterion maps to a task step.
  - AC1 → Task 2 Step 1 + Task 2 Step 5 grep.
  - AC2 → Task 2 Steps 2–3.
  - AC3 → Task 4.
  - AC4 → Task 5.
  - AC5 → Task 6.
- **Placeholder scan:** no TBDs, no "implement later", no "similar to above" shortcuts. Every code block shows the exact edit.
- **Type consistency:** constant name `MAX_AGE_MINUTES` used consistently across all tasks (Task 1 import, Task 2 definition, Task 3 docstrings, Task 5 README, verification step).
