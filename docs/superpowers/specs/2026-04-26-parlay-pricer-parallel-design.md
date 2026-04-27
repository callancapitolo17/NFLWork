# Parlay Pricer — Parallel Wagerzon Calls

**Date:** 2026-04-26
**Worktree / branch:** `worktree-perf-parlay-pricer-parallel`
**Files touched:** `wagerzon_odds/parlay_pricer.py` (only)

## Problem

`bash "Answer Keys/MLB Dashboard/run.sh"` has gotten noticeably slower since the
combined-parlay + exact-payout features landed. Tracing the pipeline shows two
fully serial loops in `wagerzon_odds/parlay_pricer.py` issuing
`ConfirmWagerHelper` POSTs one at a time:

| Stage | Function | Calls per run | Approx wall time |
|---|---|---|---|
| 1 | `_price_mlb_parlays_inner` (FG + F5, called twice from `price_mlb_parlays`) | ~64 (8 games × 4 combos × 2 periods) | 15–30s |
| 2 | `compute_exact_payouts` (NUDGE_RANGE sweep) | ~80 (≈16 sized parlays × 5 stakes) | 20–40s |

Together these account for roughly 35–70s of the run. Every other step in the
pipeline either is already parallel (the four SGP scrapers via `mclapply`) or
is fast (sample load, R compute, dashboard HTML).

## Goal

Parallelize both loops with a thread pool so the WZ network round-trips
overlap. Aim for an order-of-magnitude reduction in pricer wall time without
changing any computed values.

## Non-goals (deferred to follow-up plans)

- Running Stage 1 concurrently with `mlb_correlated_parlay.R` (option B in the
  trace report).
- Moving the four SGP scraper kickoffs out of the R script so they overlap
  Stage 1 (option C).
- Changing how Stage 2 chooses the best stake (`max(... key=...)` logic stays
  identical).
- Adding new unit tests for `parlay_pricer.py`. The change is structural and is
  validated by before/after timing + row-level diff (see Verification).

## Design

### Concurrency model

- Python `concurrent.futures.ThreadPoolExecutor` with `max_workers=8`.
- One shared `requests.Session` reused across all worker threads. `requests`
  documents its `Session` object as safe to share across threads for the
  connection-pool / cookie-jar use-case we have here, and the WZ ASP.NET
  session cookie (`ASP.NET_SessionId`) is a single-user identity that all
  threads should reuse.
- The session's HTTPS adapter gets `pool_maxsize=8` and
  `pool_connections=8` so we don't see urllib3's "Connection pool is full"
  warning under concurrent load.
- Worker function: a thin wrapper around the existing `get_parlay_price` /
  `get_parlay_price_with_fallback`. Same retry-on-amount-fallback ladder, same
  `try/except` returning `None` on individual failures.

### Shared helper

```python
def _run_parallel(fn, work_items, max_workers=8):
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

`ex.map` preserves input order, which keeps the existing print order
deterministic without extra sorting.

`max_workers` is a module-level constant so the dial-down to 4 (per the
verification fallback below) is a single-line change rather than a
multi-call-site edit.

### Stage 1 — `_price_mlb_parlays_inner`

Today: nested `for game in games: for combo in combos: get_parlay_price_with_fallback(...)`.

Change:

1. Build a flat list of `WorkItem(period, game_dict, combo_spec)` tuples for
   the period being priced (FG or F5).
2. Hand the list to `_run_parallel(_price_one_combo, work_items)`.
3. Iterate the returned list in order to print the per-combo line and append
   to `results`.

`price_mlb_parlays` (which calls `_price_mlb_parlays_inner` for FG then F5)
keeps its current shape: FG and F5 still execute sequentially as two pool
runs. Merging the two periods into one pool would save another small slice of
time but complicates the `period` filter; out of scope for this change.

### Stage 2 — `compute_exact_payouts`

Today: `for parlay in sized: for stake in range(...): get_parlay_price(...)`.

Change:

1. Build a flat list of `WorkItem(parlay_hash, kelly_bet, stake, idgm, legs)`
   tuples for every (parlay × stake) pair (≈80 items).
2. `_run_parallel(_price_one_stake, work_items)` returns one result per item.
3. Group results by `parlay_hash` (preserves input order), apply the existing
   `max(candidates, key=lambda sw: (sw[1] / sw[0], -abs(sw[0] - kelly_bet)))`
   selection per parlay, print the per-parlay summary line, and append to the
   `updates` list. The `executemany` UPDATE step is unchanged.

### Error handling

Each worker wrapper preserves today's "log + return None" contract:

- `get_parlay_price` failure (HTTP error, timeout, malformed JSON, WZ error
  key) → returns `None`. Outer code skips with the same "FAILED" / "skipping"
  print as today.
- A single bad combo or stake never aborts the pool — each task is
  independent.
- If WZ rate-limits us under 8-way concurrency (we've never tested this
  empirically), individual calls will surface as `None` and be loud in the
  output. We will dial `max_workers` down to 4 if that happens during
  verification.

### Print-order considerations

`ex.map` returns results in input order, so the existing per-combo / per-parlay
print lines come out in the same order they do today. No interleaved output.

## Verification (during build)

Run from the worktree:

1. **Baseline (main):** record wall time of current `parlay_pricer.py` and
   `parlay_pricer.py --exact-payouts` against a live WZ session.
2. **After change:** rerun both with the same DB state and confirm:
   - Stage 1 wall time drops from ~15–30s to ~3–5s.
   - Stage 2 wall time drops from ~20–40s to ~5–8s.
   - `mlb_parlay_prices` row count is identical.
   - `mlb_parlay_opportunities.exact_wager` and `exact_to_win` values match
     pre-change values (to within the WZ stake-window tie-breaking that
     already exists; identical inputs should yield identical outputs).
3. **Failure-mode probe:** scan stdout for any new `BALANCEEXCEED`,
   `auth_error`, `Request error:` lines compared to baseline. Investigate
   before merging if any appear.

## Version control

- **Worktree:** `.claude/worktrees/perf-parlay-pricer-parallel/`
- **Branch:** `worktree-perf-parlay-pricer-parallel`
- **Commits (planned):**
  1. Add `_run_parallel` helper + parallelize Stage 1
     (`_price_mlb_parlays_inner`).
  2. Parallelize Stage 2 (`compute_exact_payouts`).
- **Merge gate:** verification step 2 above passes; no new error lines from
  step 3.
- **Cleanup after merge to `main`:** `git worktree remove
  .claude/worktrees/perf-parlay-pricer-parallel` and
  `git branch -d worktree-perf-parlay-pricer-parallel`.

## Documentation

`wagerzon_odds/README.md` mentions `parlay_pricer.py` only by name; no usage
or perf claims to update. `Answer Keys/MLB Dashboard/README.md` describes the
pipeline but does not quote timings. **No README updates required.** A short
sentence on the pricer's concurrency goes in the `compute_exact_payouts` /
`price_mlb_parlays` docstrings instead.

## Out of scope / follow-ups

- **Option B (R + Stage 1 concurrent in `run.sh`):** moves the WZ Stage 1 call
  to the background, has R block on a sentinel only when it needs
  `mlb_parlay_prices`. ~5–10s additional savings. Separate plan.
- **Option C (move SGP kickoff out of R):** lets the four SGP scrapers run
  concurrently with Stage 1, eliminating the 12s SGP serial block. ~12s
  savings but invasive — splits responsibilities across `run.sh`, the pricer,
  and the R script. Separate plan.
- Persistent connection-pool warmup (TLS handshake amortization) — likely
  noise once we're already parallel.
