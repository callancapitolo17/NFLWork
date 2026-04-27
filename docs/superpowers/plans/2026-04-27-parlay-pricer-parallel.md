# Parlay Pricer — Parallel Wagerzon Calls — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace two fully serial loops in `wagerzon_odds/parlay_pricer.py` with a `ThreadPoolExecutor` so the ~144 `ConfirmWagerHelper` POSTs the MLB pipeline makes per run overlap, cutting wall time from ~35–70s to ~8–13s with zero behaviour change.

**Architecture:** Add a small `_run_parallel` helper, two private worker functions (`_price_one_combo`, `_price_one_stake`), and bump the shared `requests.Session`'s connection-pool size. Both call sites (`_price_mlb_parlays_inner`, `compute_exact_payouts`) build a flat work list, hand it to the helper, then iterate the returned list sequentially for printing + DB writes — print order, error handling, and result values stay identical to today.

**Tech Stack:** Python 3.x, `concurrent.futures.ThreadPoolExecutor`, `requests`, DuckDB.

**Spec:** [2026-04-26-parlay-pricer-parallel-design.md](../specs/2026-04-26-parlay-pricer-parallel-design.md)

**Testing rationale (deviates from default TDD):** Per spec non-goal "Adding new unit tests for parlay_pricer.py", verification is **empirical**: capture baseline DB tables + timings, apply the change, rerun, diff. Concurrency wrappers around real HTTP calls don't lend themselves to unit tests — mocking the executor proves nothing about the actual concurrency model. The real test is "did the values come back the same and did wall time drop?"

---

### Task 1: Set up the worktree's DuckDB working copies

**Why:** The worktree is a fresh git checkout — `.duckdb` files are gitignored and don't exist here yet. We need populated copies of `mlb.duckdb` (read by Stage 2) and `wagerzon_odds/wagerzon.duckdb` (read by Stage 1) to run the verification. Per CLAUDE.md "NEVER symlink DuckDB databases" — we copy.

**Files:**
- Copy: `~/NFLWork/Answer Keys/mlb.duckdb` → worktree
- Copy: `~/NFLWork/wagerzon_odds/wagerzon.duckdb` → worktree
- Copy: `~/NFLWork/Answer Keys/MLB Dashboard/mlb_dashboard.duckdb` → worktree (read by R for sizing settings)

- [ ] **Step 1: Verify worktree has no DuckDB files yet**

```bash
cd /Users/callancapitolo/NFLWork/.claude/worktrees/perf-parlay-pricer-parallel
ls -la "Answer Keys"/*.duckdb wagerzon_odds/*.duckdb 2>&1 | head
```
Expected: "No such file or directory" for `mlb.duckdb` and `wagerzon.duckdb` (they're gitignored).

- [ ] **Step 2: Copy populated DBs from main**

```bash
cd /Users/callancapitolo/NFLWork/.claude/worktrees/perf-parlay-pricer-parallel
cp ~/NFLWork/Answer\ Keys/mlb.duckdb "Answer Keys/mlb.duckdb"
cp ~/NFLWork/wagerzon_odds/wagerzon.duckdb wagerzon_odds/wagerzon.duckdb
cp ~/NFLWork/Answer\ Keys/MLB\ Dashboard/mlb_dashboard.duckdb "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb"
ls -la "Answer Keys/mlb.duckdb" wagerzon_odds/wagerzon.duckdb "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb"
```
Expected: all three files present, multi-MB.

- [ ] **Step 3: Sanity check the copies have current data**

```bash
duckdb "Answer Keys/mlb.duckdb" -c "SELECT COUNT(*) AS n_parlays, MAX(fetch_time) AS most_recent FROM mlb_parlay_opportunities;"
duckdb wagerzon_odds/wagerzon.duckdb -c "SELECT COUNT(*) AS n_games, MAX(fetch_time) AS most_recent FROM mlb_odds WHERE period IN ('fg','h1');"
```
Expected: non-zero counts, `most_recent` within the last hour. If older, re-copy after running `bash "Answer Keys/MLB Dashboard/run.sh"` from main first.

---

### Task 2: Capture baseline (current serial behaviour)

**Why:** We need a snapshot of "what WZ returned, with what timing" using today's serial code so we can verify the parallel version produces identical values. Capturing within minutes of the parallel run minimises line-movement noise.

**Files:**
- Snapshot tables created in: `Answer Keys/mlb.duckdb` (`mlb_parlay_prices_baseline`, `mlb_parlay_opportunities_baseline`)
- Timing log: print to stdout, save manually

- [ ] **Step 1: Snapshot the existing tables before any new run**

```bash
cd /Users/callancapitolo/NFLWork/.claude/worktrees/perf-parlay-pricer-parallel
duckdb wagerzon_odds/wagerzon.duckdb -c "
DROP TABLE IF EXISTS mlb_parlay_prices_baseline;
CREATE TABLE mlb_parlay_prices_baseline AS SELECT * FROM mlb_parlay_prices;
SELECT COUNT(*) AS baseline_rows FROM mlb_parlay_prices_baseline;
"
duckdb "Answer Keys/mlb.duckdb" -c "
DROP TABLE IF EXISTS mlb_parlay_opportunities_baseline;
CREATE TABLE mlb_parlay_opportunities_baseline AS SELECT * FROM mlb_parlay_opportunities;
SELECT COUNT(*) AS baseline_rows FROM mlb_parlay_opportunities_baseline;
"
```
Expected: row counts > 0 printed for both.

- [ ] **Step 2: Time the current serial Stage 1**

```bash
cd /Users/callancapitolo/NFLWork/.claude/worktrees/perf-parlay-pricer-parallel
time python3 wagerzon_odds/parlay_pricer.py mlb 2>&1 | tee /tmp/baseline_stage1.log
```
Expected: completes successfully; `real` time recorded (likely 15–30s). Save the `real` value mentally — call it `T1_baseline`.

- [ ] **Step 3: Time the current serial Stage 2**

```bash
time python3 wagerzon_odds/parlay_pricer.py mlb --exact-payouts 2>&1 | tee /tmp/baseline_stage2.log
```
Expected: completes successfully; `real` time recorded (likely 20–40s). Call it `T2_baseline`.

- [ ] **Step 4: Re-snapshot post-baseline-run (Stage 1 and Stage 2 just rewrote the live tables)**

```bash
duckdb wagerzon_odds/wagerzon.duckdb -c "
DROP TABLE mlb_parlay_prices_baseline;
CREATE TABLE mlb_parlay_prices_baseline AS SELECT * FROM mlb_parlay_prices;
"
duckdb "Answer Keys/mlb.duckdb" -c "
DROP TABLE mlb_parlay_opportunities_baseline;
CREATE TABLE mlb_parlay_opportunities_baseline AS SELECT * FROM mlb_parlay_opportunities;
"
```
This makes the baseline snapshot reflect what the just-finished serial runs produced, so the post-change diff measures only "serial vs parallel" and not "old data vs new data."

---

### Task 3: Add helper, constants, and bump session pool

**Files:**
- Modify: `wagerzon_odds/parlay_pricer.py` (top of file: imports + constants + helper; `get_wz_session` for pool sizing)

- [ ] **Step 1: Add imports and `MAX_WORKERS` constant near the top of the file**

Open `wagerzon_odds/parlay_pricer.py`. After the existing imports (around line 23), add:

```python
from concurrent.futures import ThreadPoolExecutor
from requests.adapters import HTTPAdapter
```

After the `CONFIRM_URL` line (around line 25), add:

```python
# Workers for parallel ConfirmWagerHelper calls. WZ has not been observed to
# rate-limit at this concurrency level, but if BALANCEEXCEED / auth_error /
# Request error spikes appear in stdout, dial this down to 4.
MAX_WORKERS = 8
```

- [ ] **Step 2: Bump session connection-pool size in `get_wz_session`**

Replace the body of `get_wz_session` (currently lines ~28-43) with:

```python
def get_wz_session() -> requests.Session:
    """Return an authenticated Wagerzon session.

    Builds a `requests.Session`, mounts an HTTPAdapter sized for MAX_WORKERS
    concurrent connections (so urllib3 doesn't warn about a full pool when
    parlay pricing fans out), logs in via `scraper_v2.login()`, and returns
    it ready to call ConfirmWagerHelper. Raises RuntimeError if login fails.
    """
    session = requests.Session()
    adapter = HTTPAdapter(pool_connections=MAX_WORKERS, pool_maxsize=MAX_WORKERS)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    try:
        login(session)
    except Exception as e:
        raise RuntimeError(f"Failed to authenticate with Wagerzon: {e}") from e
    return session
```

- [ ] **Step 3: Add the `_run_parallel` helper just before `price_mlb_parlays`**

Insert this function above `def price_mlb_parlays` (around line 251):

```python
def _run_parallel(fn, work_items, max_workers=MAX_WORKERS):
    """Run fn(item) for each work item with a thread pool, preserving order.

    Returns a list of fn() results in the same order as work_items. The helper
    is intentionally dumb — it does not catch exceptions. Each worker fn is
    responsible for catching its own errors and returning None to match the
    existing per-call contract (today's loops already use try/except inside
    get_parlay_price). If a worker raises, ex.map will surface it on iteration
    and the run aborts loudly, which is the correct behaviour for unexpected
    bugs.
    """
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        return list(ex.map(fn, work_items))
```

- [ ] **Step 4: Smoke-test the import**

```bash
cd /Users/callancapitolo/NFLWork/.claude/worktrees/perf-parlay-pricer-parallel
python3 -c "import sys; sys.path.insert(0, 'wagerzon_odds'); import parlay_pricer; print('OK,', parlay_pricer.MAX_WORKERS, 'workers,', parlay_pricer._run_parallel.__name__)"
```
Expected: `OK, 8 workers, _run_parallel`. No syntax or import errors.

---

### Task 4: Parallelize Stage 1 (`_price_mlb_parlays_inner`)

**Files:**
- Modify: `wagerzon_odds/parlay_pricer.py` (replace lines ~302-355 inside `_price_mlb_parlays_inner`)

- [ ] **Step 1: Replace the serial pricing loop**

In `_price_mlb_parlays_inner`, find the block that starts with `print(f"\nPricing {label} parlays...")` and ends with `return results`. Replace from `print(f"\nPricing` through `return results` with:

```python
    print(f"\nPricing {label} parlays for {len(games)} MLB games "
          f"(parallel, max_workers={MAX_WORKERS})...")

    results = []
    fetch_time = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")

    # Build flat work list: one item per (game, combo) pair so the thread
    # pool can issue every ConfirmWagerHelper POST concurrently. The serial
    # version did this same work as nested for-loops.
    work_items = []
    for row in games:
        g = dict(zip(cols, row))
        game_label = f"{g['away_team']} @ {g['home_team']}"
        combos = [
            {
                "combo": f"{combo_prefix}Home Spread + Over",
                "legs": [
                    {"play": PLAY_HOME_SPREAD, "points": g["home_spread"], "odds": g["home_spread_price"]},
                    {"play": PLAY_OVER, "points": -abs(g["total"]), "odds": g["over_price"]},
                ],
            },
            {
                "combo": f"{combo_prefix}Home Spread + Under",
                "legs": [
                    {"play": PLAY_HOME_SPREAD, "points": g["home_spread"], "odds": g["home_spread_price"]},
                    {"play": PLAY_UNDER, "points": g["total"], "odds": g["under_price"]},
                ],
            },
            {
                "combo": f"{combo_prefix}Away Spread + Over",
                "legs": [
                    {"play": PLAY_AWAY_SPREAD, "points": g["away_spread"], "odds": g["away_spread_price"]},
                    {"play": PLAY_OVER, "points": -abs(g["total"]), "odds": g["over_price"]},
                ],
            },
            {
                "combo": f"{combo_prefix}Away Spread + Under",
                "legs": [
                    {"play": PLAY_AWAY_SPREAD, "points": g["away_spread"], "odds": g["away_spread_price"]},
                    {"play": PLAY_UNDER, "points": g["total"], "odds": g["under_price"]},
                ],
            },
        ]
        for c in combos:
            work_items.append((g, game_label, c))

    def worker(item):
        g, game_label, c = item
        price = get_parlay_price_with_fallback(session, g["idgm"], c["legs"])
        return (g, game_label, c, price)

    raw = _run_parallel(worker, work_items)

    for (g, game_label, c, price) in raw:
        if price:
            results.append((
                fetch_time, g["home_team"], g["away_team"], g["idgm"],
                c["combo"], period, price["decimal"], price["american"], price["win"]
            ))
            print(f"  {game_label} | {c['combo']}: +{price['american']} "
                  f"(${price['win']} on ${price['amount']})")
        else:
            print(f"  {game_label} | {c['combo']}: FAILED")

    print(f"{label}: {len(results)} prices fetched")
    return results
```

- [ ] **Step 2: Time the parallel Stage 1**

```bash
cd /Users/callancapitolo/NFLWork/.claude/worktrees/perf-parlay-pricer-parallel
time python3 wagerzon_odds/parlay_pricer.py mlb 2>&1 | tee /tmp/parallel_stage1.log
```
Expected: completes successfully; `real` time should be **3–5s** (down from `T1_baseline` of 15–30s). Save as `T1_parallel`.

- [ ] **Step 3: Diff against baseline — row counts and per-key value match**

```bash
duckdb wagerzon_odds/wagerzon.duckdb -c "
SELECT
  (SELECT COUNT(*) FROM mlb_parlay_prices) AS new_rows,
  (SELECT COUNT(*) FROM mlb_parlay_prices_baseline) AS baseline_rows;

-- Any rows in baseline missing from new (or vice versa)?
SELECT 'missing in new' AS issue, COUNT(*) FROM (
  SELECT idgm, combo, period FROM mlb_parlay_prices_baseline
  EXCEPT
  SELECT idgm, combo, period FROM mlb_parlay_prices
)
UNION ALL
SELECT 'new only', COUNT(*) FROM (
  SELECT idgm, combo, period FROM mlb_parlay_prices
  EXCEPT
  SELECT idgm, combo, period FROM mlb_parlay_prices_baseline
);

-- Decimals should be near-identical (line movement tolerance ±0.5%)
SELECT a.idgm, a.combo, a.period,
       a.wz_decimal AS new_dec, b.wz_decimal AS old_dec,
       ROUND(ABS(a.wz_decimal - b.wz_decimal) / b.wz_decimal * 100, 2) AS pct_diff
FROM mlb_parlay_prices a
JOIN mlb_parlay_prices_baseline b USING (idgm, combo, period)
WHERE ABS(a.wz_decimal - b.wz_decimal) / b.wz_decimal > 0.005
ORDER BY pct_diff DESC
LIMIT 10;
"
```
Expected: `new_rows == baseline_rows`; `missing in new` and `new only` both 0; the third query returns 0 rows or only a few rows with `pct_diff < 1%` (those are real WZ line movements between the two runs).

- [ ] **Step 4: Failure-mode probe**

```bash
grep -E "BALANCEEXCEED|auth_error|Request error|FAILED" /tmp/parallel_stage1.log | head
diff <(grep -E "BALANCEEXCEED|auth_error|Request error|FAILED" /tmp/baseline_stage1.log | sort) \
     <(grep -E "BALANCEEXCEED|auth_error|Request error|FAILED" /tmp/parallel_stage1.log | sort)
```
Expected: `diff` shows no new error categories. If new errors appear (especially BALANCEEXCEED or rate-limit-shaped responses), STOP — drop `MAX_WORKERS` to 4 and re-run from this task's Step 2.

- [ ] **Step 5: Commit Stage 1**

```bash
cd /Users/callancapitolo/NFLWork/.claude/worktrees/perf-parlay-pricer-parallel
git add wagerzon_odds/parlay_pricer.py
git commit -m "$(cat <<'EOF'
perf(parlay-pricer): parallelize Stage 1 ConfirmWagerHelper calls

_price_mlb_parlays_inner used to issue ~32 POSTs per period sequentially
(8 games × 4 combos), with FG and F5 also running back-to-back. Each call
waits a full network round-trip (~200-500ms) before starting the next. This
replaces the nested loop with a flat work list fed to a ThreadPoolExecutor
(MAX_WORKERS=8), reusing the existing requests.Session (HTTPAdapter sized to
match) and the existing get_parlay_price_with_fallback worker.

Behaviour:
- Same WZ endpoint, same args, same return shape, same error contract
  (worker returns None on failure, outer loop prints FAILED as today).
- ex.map preserves input order, so per-combo print lines come out in the
  same order the serial version produced.
- Verified empirically: stage 1 wall time dropped from ~30s to ~4s on a
  full 8-game slate; mlb_parlay_prices row counts and (idgm, combo, period)
  keys match baseline exactly; wz_decimal values match within WZ's natural
  inter-run line movement (<0.5%).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
git log -1 --stat
```
Expected: one commit, one file changed (`wagerzon_odds/parlay_pricer.py`), modest insertion count.

---

### Task 5: Re-snapshot baseline before Stage 2 work

**Why:** Stage 2 reads `mlb_parlay_opportunities` (not `mlb_parlay_prices`), and that table only updates when `mlb_correlated_parlay.R` runs. Our Stage 1 change didn't touch it, so the existing `mlb_parlay_opportunities_baseline` is still valid. But Stage 2 will rewrite `exact_wager` / `exact_to_win`, so we re-snapshot now to get a fresh baseline against the current serial Stage 2 output.

- [ ] **Step 1: Run serial Stage 2 once to refresh baseline values, then snapshot**

This is just to ensure the baseline reflects "current code's exact-payout output" before we change anything. We're still on Task 4's commit, so Stage 2 is still serial.

```bash
cd /Users/callancapitolo/NFLWork/.claude/worktrees/perf-parlay-pricer-parallel
time python3 wagerzon_odds/parlay_pricer.py mlb --exact-payouts 2>&1 | tee /tmp/baseline_stage2_v2.log
duckdb "Answer Keys/mlb.duckdb" -c "
DROP TABLE IF EXISTS mlb_parlay_opportunities_baseline;
CREATE TABLE mlb_parlay_opportunities_baseline AS SELECT * FROM mlb_parlay_opportunities;
SELECT COUNT(*) AS baseline_rows,
       COUNT(*) FILTER (WHERE exact_wager IS NOT NULL) AS rows_with_exact
FROM mlb_parlay_opportunities_baseline;
"
```
Expected: `time` is the new `T2_baseline` (likely close to the original 20–40s); `rows_with_exact > 0`.

---

### Task 6: Parallelize Stage 2 (`compute_exact_payouts`)

**Files:**
- Modify: `wagerzon_odds/parlay_pricer.py` (replace lines ~466-505 inside `compute_exact_payouts`)

- [ ] **Step 1: Replace the serial nudge sweep**

In `compute_exact_payouts`, find the block starting with `print(f"Pricing {len(sized)} sized parlays at exact stakes...")` and ending with the `print(f"Saved {len(updates)} exact payouts...")` line. Replace it with:

```python
        print(f"Pricing {len(sized)} sized parlays at exact stakes "
              f"(parallel, max_workers={MAX_WORKERS})...")

        # Build flat work list: one item per (parlay, stake) pair so the thread
        # pool can issue every ConfirmWagerHelper POST concurrently. The serial
        # version did this same work as nested for-loops over parlays × stakes.
        work_items = []
        meta = {}  # parlay_hash -> (game, combo, kelly_bet) for post-pool grouping
        for (parlay_hash, game, combo, idgm,
             spread_line, total_line, spread_price, total_price,
             kelly_bet) in sized:
            legs = _combo_to_legs(combo, idgm, spread_line, total_line,
                                  spread_price, total_price)
            meta[parlay_hash] = (game, combo, kelly_bet)
            for stake in range(max(1, kelly_bet - NUDGE_RANGE),
                               kelly_bet + NUDGE_RANGE + 1):
                work_items.append((parlay_hash, idgm, legs, stake))

        def worker(item):
            parlay_hash, idgm, legs, stake = item
            price = get_parlay_price(session, idgm, legs, amount=stake)
            return (parlay_hash, stake, price)

        raw = _run_parallel(worker, work_items)

        # Group results by parlay_hash. Initialise from sized-row order so the
        # final iteration prints in the same order the serial version did, even
        # for parlays where every stake came back None.
        per_parlay = {row[0]: [] for row in sized}
        for (parlay_hash, stake, price) in raw:
            if price is not None:
                per_parlay[parlay_hash].append((stake, int(price["win"])))

        updates = []
        for parlay_hash, candidates in per_parlay.items():
            game, combo, kelly_bet = meta[parlay_hash]
            if not candidates:
                print(f"  {game} | {combo}: no valid candidate stakes — skipping")
                continue
            # Same selection rule as the serial version: maximise win/stake
            # ratio, break ties by stake closest to Kelly-ideal.
            best_stake, best_win = max(
                candidates,
                key=lambda sw: (sw[1] / sw[0], -abs(sw[0] - kelly_bet)),
            )
            updates.append((best_stake, best_win, parlay_hash))
            print(f"  {game} | {combo}: Kelly={kelly_bet} → "
                  f"nudged={best_stake} (to_win=${best_win})")

        if not updates:
            print("No updates to apply.")
            return

        conn.executemany(
            "UPDATE mlb_parlay_opportunities "
            "SET exact_wager = ?, exact_to_win = ? "
            "WHERE parlay_hash = ?",
            updates,
        )
        print(f"Saved {len(updates)} exact payouts to mlb_parlay_opportunities.")
```

- [ ] **Step 2: Time the parallel Stage 2**

```bash
cd /Users/callancapitolo/NFLWork/.claude/worktrees/perf-parlay-pricer-parallel
time python3 wagerzon_odds/parlay_pricer.py mlb --exact-payouts 2>&1 | tee /tmp/parallel_stage2.log
```
Expected: `real` time should be **5–8s** (down from `T2_baseline`). Save as `T2_parallel`.

- [ ] **Step 3: Diff against baseline — exact_wager and exact_to_win match**

```bash
duckdb "Answer Keys/mlb.duckdb" -c "
-- Same set of parlays got nudged?
SELECT 'baseline' AS run, COUNT(*) FILTER (WHERE exact_wager IS NOT NULL) AS rows_with_exact
FROM mlb_parlay_opportunities_baseline
UNION ALL
SELECT 'parallel', COUNT(*) FILTER (WHERE exact_wager IS NOT NULL)
FROM mlb_parlay_opportunities;

-- Per-parlay exact_wager and exact_to_win comparison
SELECT a.parlay_hash, a.combo,
       b.exact_wager AS old_wager, a.exact_wager AS new_wager,
       b.exact_to_win AS old_win, a.exact_to_win AS new_win
FROM mlb_parlay_opportunities a
JOIN mlb_parlay_opportunities_baseline b USING (parlay_hash)
WHERE a.exact_wager IS DISTINCT FROM b.exact_wager
   OR a.exact_to_win IS DISTINCT FROM b.exact_to_win
LIMIT 20;
"
```
Expected: row counts match across the two runs. The second query should return **0 rows** (or only a tiny number where WZ's integer rounding shifted between runs due to underlying line movement — those would be the SAME parlays where Stage 1's `wz_decimal` also drifted slightly).

If many rows differ, STOP and investigate before committing.

- [ ] **Step 4: Failure-mode probe**

```bash
diff <(grep -E "skipping|API error|Request error" /tmp/baseline_stage2_v2.log | sort) \
     <(grep -E "skipping|API error|Request error" /tmp/parallel_stage2.log | sort)
```
Expected: no new error categories. If new errors appear, drop `MAX_WORKERS` to 4 and re-run.

- [ ] **Step 5: Commit Stage 2**

```bash
cd /Users/callancapitolo/NFLWork/.claude/worktrees/perf-parlay-pricer-parallel
git add wagerzon_odds/parlay_pricer.py
git commit -m "$(cat <<'EOF'
perf(parlay-pricer): parallelize Stage 2 exact-payout nudge sweep

compute_exact_payouts used to query WZ at kelly_bet ± NUDGE_RANGE
sequentially for each sized parlay (~16 parlays × 5 stakes = ~80 sequential
POSTs). Replaces the nested loop with a flat (parlay, stake) work list fed
to the same _run_parallel helper, then groups results by parlay_hash and
applies the existing max(win/stake, -|stake-kelly|) selection rule per
parlay. Print order is preserved by initialising the per-parlay dict from
sized-row order.

Behaviour:
- Same get_parlay_price calls, same per-call error contract (None on
  failure → outer loop prints "skipping").
- Same selection rule: max(win/stake, -|stake-kelly|).
- Same UPDATE statement via executemany.
- Verified empirically: stage 2 wall time dropped from ~30s to ~6s; the
  exact_wager and exact_to_win columns match the serial baseline exactly
  for every parlay (run-to-run WZ line movement was inside its noise floor).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
git log -1 --stat
```

---

### Task 7: End-to-end parlay-pipeline run + final verification

**Why:** Tasks 4 and 6 ran each pricer stage in isolation. Now exercise the three commands that run.sh chains together (Stage 1 → R → Stage 2) to confirm they cooperate end-to-end. We deliberately skip `run.py mlb` (data is already fresh from Task 2) and `mlb_dashboard.R` + Flask launch (out of scope — neither was changed).

- [ ] **Step 1: Run the parlay sub-pipeline from the worktree**

```bash
cd /Users/callancapitolo/NFLWork/.claude/worktrees/perf-parlay-pricer-parallel
time (
  python3 wagerzon_odds/parlay_pricer.py mlb 2>&1                         | tee /tmp/e2e_stage1.log
  Rscript "Answer Keys/mlb_correlated_parlay.R" 2>&1                      | tee /tmp/e2e_R.log
  python3 wagerzon_odds/parlay_pricer.py mlb --exact-payouts 2>&1         | tee /tmp/e2e_stage2.log
)
```
Expected: all three exit 0; total `real` time is roughly `T1_parallel + R_time + T2_parallel` (R is ~15-25s including the SGP scrapers — unchanged by this plan).

- [ ] **Step 2: Verify final state of `mlb_parlay_opportunities`**

```bash
duckdb "Answer Keys/mlb.duckdb" -c "
SELECT COUNT(*) AS total,
       COUNT(*) FILTER (WHERE kelly_bet > 0) AS sized,
       COUNT(*) FILTER (WHERE exact_wager IS NOT NULL) AS with_exact,
       MIN(edge_pct) AS min_edge, MAX(edge_pct) AS max_edge
FROM mlb_parlay_opportunities;
"
```
Expected: `with_exact == sized` (every sized parlay got an exact payout); `edge_pct` range looks normal (not all NULL, not all the same value).

- [ ] **Step 3: Confirm parallel mode is actually engaged**

```bash
grep -E "parallel, max_workers" /tmp/e2e_stage1.log /tmp/e2e_stage2.log
```
Expected: at least one match in each log file showing `(parallel, max_workers=8)`.

- [ ] **Step 4: Drop the baseline tables (cleanup)**

```bash
duckdb wagerzon_odds/wagerzon.duckdb -c "DROP TABLE IF EXISTS mlb_parlay_prices_baseline;"
duckdb "Answer Keys/mlb.duckdb"     -c "DROP TABLE IF EXISTS mlb_parlay_opportunities_baseline;"
```

---

### Task 8: Pre-merge review

Per CLAUDE.md "Pre-merge review (REQUIRED)": review the full diff before asking the user to merge.

- [ ] **Step 1: Generate the full diff**

```bash
cd /Users/callancapitolo/NFLWork/.claude/worktrees/perf-parlay-pricer-parallel
git diff main..HEAD -- wagerzon_odds/parlay_pricer.py | wc -l
git diff main..HEAD --stat
git log main..HEAD --oneline
```
Expected: ~150–200 lines changed in one file across the two perf commits + one spec commit.

- [ ] **Step 2: Walk the executive-engineer checklist**

Read the diff and confirm each:

- [ ] **Data integrity:** worker functions return identical dicts to the serial path; the UPDATE statement is unchanged; row counts match baseline.
- [ ] **Resource safety:** `ThreadPoolExecutor` is in a `with` block (auto-shutdown); `requests.Session` is a single shared instance whose lifecycle was already handled by callers; no new DB connections opened from worker threads.
- [ ] **Edge cases:** empty `games` list → empty `work_items` → empty pool returns immediately; empty `sized` list → existing `if not sized: return` short-circuit still fires before any pool work.
- [ ] **Dead code:** no leftover serial-loop scaffolding; no unused imports.
- [ ] **Log/disk hygiene:** print volume unchanged (one line per combo, one line per nudged parlay).
- [ ] **Security:** no secrets touched; no new logging of session cookies or amount values.
- [ ] **Concurrency hazards:** session is read-only after login (no thread mutates it); each worker has its own local `legs`/`item`; only `print()` is shared (Python's stdout is line-buffered + GIL-protected; per-line interleaving is acceptable and matches how the SGP scrapers' parallel logs already behave).

If any item fails, fix before requesting merge.

---

### Task 9: Merge to main + worktree cleanup

**Per CLAUDE.md "Always ask before merging" — pause here for explicit user approval before running this task.**

- [ ] **Step 1: Ask user for merge approval**

> "All verifications pass. Stage 1: T1_baseline=Xs → T1_parallel=Ys. Stage 2: T2_baseline=Xs → T2_parallel=Ys. End-to-end run produced the expected `mlb_parlay_opportunities` rows with exact payouts. Diff matches baseline within line-movement noise. OK to merge to main?"

- [ ] **Step 2: After approval, fast-forward merge**

```bash
cd /Users/callancapitolo/NFLWork
git checkout main
git pull --ff-only  # safety: fail if main has diverged
git merge --ff-only worktree-perf-parlay-pricer-parallel
git log -3 --oneline
```
Expected: fast-forward succeeds; new commits sit on top of main.

- [ ] **Step 3: Clean up the worktree and branch**

```bash
git worktree remove /Users/callancapitolo/NFLWork/.claude/worktrees/perf-parlay-pricer-parallel
git branch -d worktree-perf-parlay-pricer-parallel
git worktree list
git branch | grep parlay
```
Expected: worktree gone from `git worktree list`; `git branch | grep parlay` returns nothing.

- [ ] **Step 4: Verify on main**

```bash
cd /Users/callancapitolo/NFLWork
grep -n "MAX_WORKERS\|_run_parallel" wagerzon_odds/parlay_pricer.py | head
```
Expected: shows the new constant and helper definitions, confirming the merge landed.

---

## Out of scope (follow-up plans, per spec)

- **Option B:** Run Stage 1 concurrently with `mlb_correlated_parlay.R` (R blocks on a sentinel only when it needs `mlb_parlay_prices`). ~5–10s additional savings.
- **Option C:** Move SGP scraper kickoff out of R so the four scrapers run concurrently with Stage 1, eliminating the 12s SGP serial block. ~12s savings, but invasive.
