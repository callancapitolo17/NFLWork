# Issue #99 — leg surface: staleness gate + pre-quote constituent freshness veto

**Epic:** #94. **Depends on:** #96 (ingest loop), #98 (router reads the surface).
**Base:** `main` @ e25c882. **Branch:** `feature/issue-99-surface-staleness-gate`
in worktree `.claude/worktrees/upbeat-solomon-69a657`.

## What is broken today

#98 shipped surface pricing with **no age bound**. `LegSurface.book_fairs`
returns whatever is published, and a failed ingest pass deliberately publishes
NOTHING and keeps the book's previous rows (an empty slice would blank a live
book on a blip). So a book that goes dark **holds its last prices forever**, and
the maker will happily quote off them. The age gate is the mechanism that turns
"dark book" into "countable decline".

## Guard 1 — age gate

A quote may only use surface rows younger than `SURFACE_MAX_AGE_SEC = 30`
(user decision 2026-08-25, reaffirmed in #99's comment).

**Where.** Above the store, never inside it. `LegSurface.book_fairs` keeps
returning rows of any age on purpose — filtering there would make the decline
counts unreadable (same reason #98 didn't do it). The gate lives in a small
callable in `main.py` that wraps the store and **records what it dropped**:

```
class _SurfaceAgeGate:          # the router's `surface_fairs` argument
    __call__(game_legs) -> {book: fair}      # fresh rows only
    used:        {leg_set_hash: {book: age_sec}}
    excluded:    {leg_set_hash: {book: age_sec}}   # dropped by age
```

It is callable, so every existing `surface_fairs=...` call site is unchanged.
Constructed per pricing attempt so its `now` and its records are attributable
to one RFQ.

**Decline naming — three different things, three different reasons:**

| situation | reason |
|---|---|
| no book published this leg at all | `surface_too_few_books` (#98, unchanged) |
| rows existed but age-excluded ⇒ consensus thin | **`surface_stale`** (new) |
| enough fresh rows, but they disagree | `surface_dispersion` (#98, unchanged) |

The classification is exact, not a re-read: the gate recorded the exclusions
during the very call that produced the gate reason. If the router returns
`surface_too_few_books` **and** the gate excluded ≥1 row by age → `surface_stale`.

**DK-excluded-by-age is counted distinctly from "DK had no row".** Three
instruments:
1. `quote_priced.surface_games` gains `excluded_by_age: {book: age_sec}` and
   `oldest_used_age_sec` per game; `books` continues to list only what was USED.
   So query 18 keeps measuring used rows, and the exclusions are a sibling key.
2. A periodic `surface_age_summary` research event + INFO log every
   `COVERAGE_SUMMARY_SEC` (300s) — per book, `used` vs `excluded_by_age` counts.
   Exact precedent: #81's `on_demand_coverage`.
3. A **startup WARNING** naming every configured book whose cadence alone is
   ≥ `SURFACE_MAX_AGE_SEC`, i.e. structurally excluded. DraftKings (60s singles,
   21–28s scrape) trips this on every boot. Loud, deliberate, not silent.

**Where the gate applies:**

| call site | gated? | why |
|---|---|---|
| discovery tick pricing | **yes** | the ticket |
| confirm last look (`_confirm_tick`) | **yes** | this is the FILL moment — the strictest place. Void reason `voided_surface_stale`, split out of `voided_no_fresh_books` |
| risk sweep drift (`_current_consensus_fair`) | **yes** | a stale row would produce a fake drift number; gated it goes quiet, and `constituent_jump` + tipoff still cover the resting quote |
| post-fill cooldown (`_surface_refreshed_since`) | **no** | it already requires `built_at` > `filled_at`, which is strictly stronger in the direction that matters. If those rows later go stale, the quote path's gate declines anyway. Documented, not changed |
| same-game / live on-demand path | **no** | out of scope. `QUOTE_FRESH_SEC` is that path's rule and is untouched |

**Rollback:** `SURFACE_MAX_AGE_SEC=0` disables the gate (documented explicitly as
"no age bound"), which is byte-for-byte #98 behaviour. Config only.

## Guard 2 — pre-quote constituent freshness veto

Before quoting, ask Kalshi's own constituent single-leg markets whether the
market moved **since the surface row was built**. If it did, the book's cached
number is stale by construction — refuse.

**The baseline problem.** We have Kalshi *now* (the #17/#23 snapshot). We need
Kalshi *at `built_at`*. That requires history, and history must cost zero API
calls. So: a new pure in-memory `ConstituentTape`
(`kalshi_mlb_mm/constituent_tape.py`) that **remembers reads we already make**:

* discovery tick — `_leg_market_prices(legs)` (the #17 quote-time snapshot)
* confirm tick — `_leg_market_prices(legs)` (the last look re-read)
* risk sweep — `singles.fetch_market_prices(...)` (the 10s poll)

`record(prices, observed_at)` stores devigged P(YES) per ticker;
`price_at_or_before(ticker, when)` answers the veto. Bounded by
`CONSTITUENT_TAPE_RETENTION_SEC` (180s) and a per-ticker point cap, pruned on
write. **Zero extra Kalshi API calls** — the tape is only a memory of existing ones.

Coverage is good in practice because the discovery tick re-enters every open RFQ
every ~2s and snapshots its legs each time it reaches the quote step, so an
actively-quoted leg is sampled far faster than the ≤30s row age we compare against.

**The check** (discovery tick, immediately after the leg snapshot, beside the
existing `corr_sanity` gate that consumes the same snapshot). For each
surface-routed game group of the combo:

```
ticker    = the group's leg -> market_ticker  (from legset.parse_leg per raw leg)
p_now     = singles.devigged_yes(snapshot[ticker])
oldest    = min built_at over the books that BACKED this group   (from the gate's records)
baseline  = tape.price_at_or_before(ticker, oldest)
if baseline is None            -> no_baseline, fail-open, counted
if |p_now - baseline| > SURFACE_CONSTITUENT_MOVE_THRESHOLD -> decline
```

Decline reason: **`surface_constituent_moved`**.

**Design decisions to call out:**

* **Combo-level, not per-row.** The issue says "refuse the row". By the time we
  hold current Kalshi, the fair has already run the margin / size-gate /
  exposure-cap / hysteresis chain; refusing individual rows would mean
  re-running all of it. Declining the whole combo is fail-closed, far simpler,
  and forfeits little — the age gate already bounds the oldest backing row to
  30s. Taking the snapshot *before* pricing instead is rejected: it would spend
  2–3 Kalshi GETs on every RFQ that then declines on a price gate, which at this
  RFQ volume is a large, permanent increase in Kalshi load.
* **`min(built_at)` over backing books** — the veto asks whether *any* input we
  used predates the market's last move; the oldest row is the weak link.
* **No baseline ⇒ fail-open, counted.** Same contract as `jumped_tickers`
  ("absence of signal is not a jump"), and unavoidable without extra calls.
  Counted in `surface_constituent_check` so the miss rate is measurable.
* **Not repeated at confirm.** #17's `singles_moved` already voids on ANY
  one-tick move of ANY leg since the quote snapshot — strictly stricter than
  this veto. Adding it there would be dead code.
* **Threshold** `SURFACE_CONSTITUENT_MOVE_THRESHOLD` defaults to
  `CONSTITUENT_JUMP_THRESHOLD`'s value (0.03), same units (|Δ devigged P(YES)|).
  Kept equal until the emitted delta distribution says otherwise — a decline is
  cheaper than #23's cancel, so it may want tightening later.
* **Rollback:** `SURFACE_CONSTITUENT_VETO_ENABLED=false` (no tape recording, no
  veto). Config only.

**Research:** `surface_constituent_check` (emitted for every surface-routed
combo reaching the veto, mirroring `corr_sanity_check`): per leg `ticker`,
`p_now`, `p_baseline`, `delta`, `baseline_age_sec`, `oldest_row_age_sec`,
`verdict` ∈ {ok, moved, no_baseline, unreadable}. That is the tuning dataset.

## Cadence — the decision #96 deferred

**Keep `SURFACE_CADENCE_DEFAULT_SEC = 20`.**

Measured (#98's live run, to be re-verified in this ticket's live run): structure
books p50 7–8s, p95 19s, max 20.6s. Against a 30s gate that is ~10s of headroom.
One skipped pass lands a book at ~40s → that book drops for one cycle; with three
structure books, `MIN_AGREEING_BOOKS=2` still clears. Cutting to 15s costs +33%
requests (~350–690k/day, per the README cost table) to buy 5s of headroom we have
no evidence we need, and ProphetX already 403s us. The new
`surface_age_summary` counters are the instrument that will justify a change.

**DraftKings cannot clear a 30s gate at any cadence** — its slate scrape alone is
21–28s, so its rows are 30–90s old. It is deliberately excluded, warned at
startup, and counted per-book. #96's README predicted exactly this.

## Files

| file | change |
|---|---|
| `kalshi_mlb_mm/constituent_tape.py` | **new** — pure bounded tape, no I/O |
| `kalshi_mlb_mm/config.py` | `SURFACE_MAX_AGE_SEC`, `SURFACE_CONSTITUENT_VETO_ENABLED`, `SURFACE_CONSTITUENT_MOVE_THRESHOLD`, `CONSTITUENT_TAPE_RETENTION_SEC`; rewrite the "Caveat until #99" comment |
| `kalshi_mlb_mm/main.py` | `_SurfaceAgeGate`; `_surface_enabled()` split from the lookup; veto + tape feeding at 3 existing read sites; new decline/void reasons; `surface_games` trace gains exclusions; periodic `surface_age_summary`; startup warning |
| `kalshi_mlb_mm/leg_surface/store.py` | comment only — `book_fairs` stays age-blind (say the gate now exists and lives above) |
| `kalshi_mlb_mm/research_queries.sql` | query 19 gains the new reasons; new queries for age exclusions + veto deltas |
| `kalshi_mlb_mm/report.py` | staleness section reports age-excluded counts |
| `kalshi_mlb_mm/README.md` | replace "Caveat until #99" with the shipped gates; cadence decision; new env rows |
| `CLAUDE.md` | maker blurb: #99 shipped |
| `tests/` | unit tests for the tape + the gate |

## Version control

* Worktree `.claude/worktrees/upbeat-solomon-69a657`, branch
  `feature/issue-99-surface-staleness-gate` (created before any file was written).
* Commits: (1) config + tape + unit tests, (2) main.py gates, (3) observability
  (research/report/queries), (4) docs — README + CLAUDE.md in the same merge.
* Live acceptance run BEFORE the pre-merge review. Pre-merge review against the
  CLAUDE.md checklist, findings presented, then explicit user approval to merge.
* Cleanup after merge: `git worktree remove` + `git branch -d`. `kalshi_mlb_mm/.env`
  is copied in for the live run and **deleted** before cleanup.

## Live acceptance (real books, bot NOT started as a daemon)

Harness precedent: #98's run drove the real `main._discovery_tick` in dry-run
against `RestRFQSource` with a counting engine wrapper.

1. Run the real `SurfaceIngest` against live books until the surface is warm.
2. Drive `_discovery_tick` in dry-run over real open RFQs. Assert: zero
   single-leg on-demand flights (the #98 invariant still holds), and record the
   per-book used/excluded-by-age split — expect DK ~always excluded, structure
   books ~always used.
3. **Age gate fires:** freeze a book's ingest (stop its worker) and re-tick;
   assert its rows get excluded and, once too few remain, a `surface_stale`
   decline appears — not `surface_too_few_books`.
4. **Veto fires:** with a real warm surface, inject a moved Kalshi price for one
   constituent (the tape and snapshot are both in-process) and assert
   `surface_constituent_moved`, per the issue's acceptance criterion.
5. Re-measure query 18's age percentiles to confirm the 30s/20s pairing.
6. Confirm `sgp_fetch_health` shows no new book traffic from any of this.
