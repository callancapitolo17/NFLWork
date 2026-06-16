# MLB SGP Scraper Speedup — Design

**Date:** 2026-06-08
**Branch:** `worktree-sgp-scraper-speedup`
**Status:** Draft for review

---

## Review Pack

**What we're building** — The SGP price-fetching layer is the speed
bottleneck for both Kalshi bots: each refresh cycle takes ~160s against a
60s target, so book prices are ~3 minutes stale by the time the bots quote
or accept. We're making the scrapers faster three ways: parallelizing the
per-target pricing calls (DK/FD currently price one target at a time),
replacing the subprocess-per-cycle model with a long-lived in-process
`SGPService` that both bots share, and caching slow-but-stable fetches
(event lists, DK's 2MB selection-ID dictionary) so only the actual prices
are refetched every cycle.

**Key decisions**

1. **In-process `SGPService` for the bots; CLI shims stay for the
   dashboard.** Alternative rejected: persistent worker *processes* with
   an IPC channel. The orchestrators are already pure library functions
   (`price_sgps: targets → rows`), so calling them in-process deletes
   complexity (no subprocess spawn, no DB round-trip for results) instead
   of adding it. The dashboard keeps its existing subprocess path and is
   untouched.
2. **Architecture changes cover all 4 books and both bots; probe effort
   focuses on DK + ProphetX.** Alternative rejected: DK+PX-only changes.
   Both bots flow through one shared module (`kalshi_common/sgp_runner.py`),
   so the service rebuild serves maker and taker at once. The concurrency
   probe focuses where the wall-clock is: DK (sequential today) and PX
   (2-wide today). FD/NV get sensible defaults, not probe campaigns.
3. **Asymmetric aggression: DK aggressive, PX gentle.** DK's
   `calculateBets` is a pure price calculation — worst case is a
   recoverable rate-limit/IP block. PX's endpoint submits *real RFQs* to
   market makers; hundreds of rapid RFQs are a detectable footprint
   (adversarial-review vector E3) with account-level risk. The probe ramps
   DK hard and PX cautiously, stopping at the first degradation signal.
4. **Cache the selection-ID dictionaries (TTL ~3–5 min), never the
   prices.** The sel-id dict is "which bets exist" — stable for minutes.
   Prices come fresh from `calculateBets`/`implyBets`/RFQ every cycle
   regardless, so the cache has zero freshness cost on the numbers the
   bots trade on.
5. **Phased delivery with a probe first.** Phase 0 measures each book's
   concurrency ceiling before we hard-code worker counts. Alternative
   rejected: guessing ceilings and shipping — a wrong guess either leaves
   speed on the table or triggers a DK block during the season.

**Risks / push back here** *(both resolved 2026-06-08 — decisions below)*

- **PX RFQ footprint — RESOLVED: probe gently to 4–6, ship at ceiling;
  per-book cadence is the real dial.** Parallelism does not change RFQ
  *volume* per cycle (~1,100 RFQs either way — it only compresses the
  window); the volume increase comes from faster cycles (~2.5× more
  RFQs/day once cycles hit 60s), which happens regardless of PX's worker
  count. Mitigation: `SGPService` supports per-book refresh intervals
  from day one (cheap staleness check per book), so PX can be refreshed
  every 2–3 cycles later without rearchitecting if footprint becomes a
  concern. Not enabled by default in v1 — measure first.
- **DK probe Akamai risk — RESOLVED: run ~7–8am PT (no games in
  progress, hours of recovery buffer before first pitch), hard call
  budget ~150 calls for the full DK ramp (noise vs ~6,000 calls in one
  normal cycle), and a post-probe health check (one vanilla
  `calculateBets`) so a tripped defense is detected immediately, not at
  game time. If blocked at level N, ship the last clean level.
- **Bot restart required.** Both bots must restart to pick up the
  in-process service (Phase 2). Restart gotchas are documented (SIGTERM
  starvation, must run from main repo cwd) but it's still an operational
  touch on live-money processes.
- **Scope check:** this spec does NOT build the top-down bivariate-Poisson
  pricer (separate future project) and does NOT touch Wagerzon scraping,
  `mlb_parlay_lines`, the R blender, or any dashboard rendering.

**Worth understanding** (opt-in)

1. **Thread pools for I/O-bound work** — like R's `parallel::mclapply`,
   which this repo already uses to run the 4 scrapers side by side. Python
   threads don't speed up *computation* (the GIL serializes Python
   bytecode), but they're perfect when each task spends its time *waiting
   on the network*: while one request waits, another sends. That's why
   `ThreadPoolExecutor(max_workers=N)` over HTTP calls gives near-linear
   speedup up to the server's tolerance.
2. **Connection reuse / TLS handshake cost** — every fresh HTTPS
   connection pays a multi-round-trip setup (TCP + TLS + Akamai's TLS
   fingerprint check). A persistent session pays it once and reuses the
   socket. Subprocess-per-cycle forces a fresh handshake every cycle ×
   every book — pure overhead.
3. **TTL caching** — "time-to-live": keep a fetched value and reuse it
   until it's older than X, then refetch. The art is choosing what's safe
   to cache (market *structure*, stable for minutes) vs. what never is
   (market *prices*).

---

## 1. Problem statement

Observed from `kalshi_mlb_rfq/bot.log` (2026-06-03): every `sgp_cycle`
takes 153–178s wall-clock against `SGP_REFRESH_SEC = 60`. The bots
therefore trade on book fairs that are ~2.5–3 minutes old. Both the taker
(`kalshi_mlb_rfq`) and the maker (`kalshi_mlb_mm`) consume SGP odds through
the same `kalshi_common/sgp_runner.py::sgp_cycle`, so both inherit the
staleness.

### Where the time goes (per cycle, last observed: 495 targets, 15 games)

| Cost | Mechanism |
|---|---|
| **DK sequential target loop** (~146s) | `mlb_sgp/draftkings.py::price_sgps` iterates ~365 surviving targets one at a time; only the 4 combos *within* a target run in parallel. ~365 × ~400ms RTT. |
| **PX 2-wide target loop** | `PX_TARGET_PARALLELISM = 2`; 1,118 rows priced through a 2-lane funnel. |
| **Subprocess-per-cycle** | `sgp_runner.run_scrapers` spawns 4 fresh Python processes every tick: 4× interpreter boot, 4× TLS handshake, 4× event-list fetch, DK's 2MB sel-id fetch per game — all redone every cycle. |
| **DB round-trip for results** | Targets written to `mlb_target_lines`, rows read back from `mlb_sgp_odds` with lock-retry logic, purely to pass data between processes on the same machine. |

FD prices few rows (96 — strict line matching filters most targets) and NV
is already 4-wide; neither is the long pole.

### Current state of target-level parallelism

| Book | Target-level parallelism | Combo-level | Long pole? |
|---|---|---|---|
| DraftKings | **none (sequential)** | 4-wide | **yes** |
| ProphetX | 2-wide | 4-wide | **yes** |
| FanDuel | none (sequential) | 4-wide | no (few targets survive) |
| Novig | 4-wide | 4-wide | no |

## 2. Goals / non-goals

**Goals**
- SGP refresh cycle ≤ 60s for the current line surface (~500 targets,
  15 games), with headroom for growth.
- One shared fetch architecture serving taker, maker, and dashboard.
- Empirically measured (not guessed) concurrency ceilings for DK and PX.
- No behavior change in priced output: same rows, same `sgp_decimal`
  values, same source labels.

**Non-goals**
- Top-down internal SGP pricer (bivariate Poisson) — separate project.
- Any change to Wagerzon scraping, `mlb_parlay_lines`, the R blender, or
  dashboard rendering.
- WebSocket/streaming price feeds.
- Changing what the bots *do* with prices (gates, Kelly, RFQ logic).

## 3. Design

### 3.1 Component overview

```
                 kalshi_common/sgp_runner.py
                 ┌──────────────────────────────────────┐
 taker main ────▶│  SGPService  (NEW)                   │
 maker main ────▶│    .refresh(targets) → list[PricedRow]│
                 │                                      │
                 │  Holds across cycles:                │
                 │   dk_client · fd_client ·            │
                 │   px_client · nv_client  (persistent)│
                 │                                      │
                 │  Per refresh():                      │
                 │   ├─ 4 books concurrently (1 thread  │
                 │   │   per book, bounded deadline)    │
                 │   ├─ TTL caches: event lists,        │
                 │   │   sel-id / runner dicts          │
                 │   └─ auto-reinit a failed client     │
                 └──────────────────────────────────────┘
                            │ in-process calls
                 mlb_sgp/{draftkings,fanduel,prophetx,novig}.py
                    price_sgps(targets, client=...)
                    └─ DK & FD gain target-level ThreadPoolExecutor
                       (worker counts from Phase 0 probe)

 DASHBOARD (unchanged path): mlb_correlated_parlay.R
   └─ spawns scraper_*_sgp.py CLI shims → price_sgps() → mlb_mm.duckdb
      (inherits the Phase 1 parallelism win for free; no other change)
```

### 3.2 Phase 0 — concurrency probe harness

Standalone script `mlb_sgp/probe_concurrency.py` (kept in the repo —
rerun when books change defenses; it is a tool, not a temp file).

- For each concurrency level N in a ramp schedule, fire N simultaneous
  *representative* pricing calls (DK: `calculateBets` on real sel-id pairs;
  PX: parlay RFQs on real legs), then record per-level: p50/p95 latency,
  HTTP status mix (200 / 422 / 429 / 403), and error bodies.
- **DK ramp:** 2 → 4 → 8 → 12 → 16 → 24, stop at first degradation
  (latency knee, 429s, or Akamai challenge), back off one level.
- **PX ramp:** 2 → 3 → 4 → 6, stop at *any* anomaly (slower offers,
  empty offer ladders, 4xx). PX calls are real RFQs — the probe uses the
  minimum stake and a small call budget per level (~10 calls), and we
  accept a coarser ceiling estimate in exchange for a small footprint.
- Output: a printed table + recommended `DK_TARGET_PARALLELISM` /
  `PX_TARGET_PARALLELISM` values. Run manually; never from the bots.
- Cooldown pauses between levels so a triggered defense at level N
  doesn't contaminate the read at N+1.
- **Run window:** ~7–8am PT (no MLB games in progress; hours of recovery
  buffer before first pitch). **Hard budgets:** ≤150 total DK calls,
  ≤40 total PX RFQs across the whole ramp. **Post-probe health check:**
  one vanilla `calculateBets` (DK) and one minimum-stake RFQ (PX) after
  the ramp completes, so a tripped defense is detected immediately.

### 3.3 Phase 1 — target-level parallelism in the orchestrators

`mlb_sgp/draftkings.py` and `mlb_sgp/fanduel.py`: wrap the per-target
pricing block (currently `for t in filtered_targets:` body) in a
`ThreadPoolExecutor`, mirroring the existing pattern in `novig.py` /
`prophetx.py` (`_price_one_target` + `as_completed`). Worker counts come
from module-level constants set by the probe:

- `DK_TARGET_PARALLELISM` (new; probe-derived, expected 8–16)
- `FD_TARGET_PARALLELISM` (new; conservative default 4 — FD is not a
  long pole, no probe campaign)
- `PX_TARGET_PARALLELISM` (existing; raised per probe, expected 3–4)
- NV unchanged at 4.

Within each target, the existing 4-wide combo parallelism stays. Total
in-flight requests per book = `TARGET_PARALLELISM × 4`; the probe measures
at this granularity (simultaneous *requests*, not targets) so the product
is what's validated.

Notes:
- The shared `curl_cffi` session is already used concurrently by the
  combo-level pool (proven in production); we are widening, not
  introducing, multi-threaded session use.
- Integer-fallback targets (`try_integer_fallback_*`, 8 internal calls)
  run inside the same worker — no nested pools.
- Order-independence: rows are keyed by (game, period, combo, lines), so
  `as_completed` ordering does not affect output content. The regression
  tests compare row sets, not order.

### 3.4 Phase 2 — `SGPService` (in-process, both bots)

New class.

> **Implementation note (as shipped):** `SGPService` lives in its own
> module `kalshi_common/sgp_service.py` (re-exported via
> `sgp_runner.py`), not inside `sgp_runner.py` as sketched below.

```python
class SGPService:
    def __init__(self, books=("draftkings", "fanduel", "prophetx", "novig"),
                 per_book_deadline_sec=75): ...
    def refresh(self, targets: list[TargetLine]) -> list[PricedRow]:
        # 1 thread per book; each calls mlb_sgp.<book>.price_sgps(
        #    targets, periods=("FG",), client=self._client(book))
        # future.result(timeout=per_book_deadline) per book;
        # a late/hung book contributes [] this cycle and its client is
        # torn down + lazily rebuilt next cycle.
    def close(self): ...
```

- **Persistent clients.** Created lazily on first use, held across
  cycles. A book whose calls start failing (3 consecutive, matching the
  existing auto-reinit heuristic) gets its client rebuilt.
- **Per-book refresh intervals.** `SGPService` carries an optional
  `min_refresh_sec` per book (default 0 = every cycle). On `refresh()`,
  a book whose last successful fetch is younger than its interval is
  skipped and its previous rows are not re-emitted (callers already
  treat `mlb_sgp_odds` freshness via `fetch_time`). This is the PX
  footprint dial — present from day one, not enabled in v1.
- **Bot integration.** `sgp_cycle()` keeps its signature but gains a
  service-backed implementation: enumerate targets → `service.refresh()`
  → write rows to the bot's market DB (same `mlb_sgp_odds` schema — the
  research DB, dashboards-on-bot-DB, and `read_priced_rows` consumers
  keep working). `mlb_target_lines` is still written (tipoff gating and
  debugging read it). What disappears: `run_scrapers` subprocess spawning
  and the results round-trip through the DB.
- Both `kalshi_mlb_rfq/main.py` and `kalshi_mlb_mm/main.py` construct one
  `SGPService` at startup and pass it to `sgp_cycle`. The taker's
  `kalshi_mlb_rfq/sgp_runner.py` shim re-exports unchanged.
- **Threading containment:** `refresh()` is called from the bots'
  existing SGP cadence position in the main loop (synchronous, like
  `sgp_cycle` today) — the win is eliminating spawn + handshake + DB
  round-trip and adding within-call parallelism, not changing the loop's
  concurrency model. The per-book deadline replaces
  `SGP_SCRAPER_TIMEOUT_SEC` as the runaway guard.
- The CLI shims (`scraper_*_sgp.py`) remain for the dashboard and are
  not modified beyond what Phase 1 does to the orchestrators they call.

### 3.5 Phase 3 — TTL caching of structure fetches

Inside `SGPService` (not inside the orchestrators — the dashboard path
must stay stateless):

> **Implementation note (as shipped):** structure caching is implemented
> for **DK and FD only** in v1. PX and NV get the persistent-client win
> but not structure caching — their per-cycle structure fetches are only
> 1–2 cheap calls, so the PX/NV row in the table below is deferred.

| Cached value | Fetch today | TTL | Rationale |
|---|---|---|---|
| DK event list (`fetch_dk_events`) | every cycle | 5 min | Game list changes ~daily; live filtering already excludes started games downstream. |
| DK per-game `main_market_nums` + `sel_ids` (2MB/game) | every cycle | 3 min | Selection IDs are market *structure*; only churn when DK re-mains a line. Stale entry worst case: a sel-id pair 404s/422s → combo skipped that cycle → refetched on expiry. |
| FD event list + per-game runners | every cycle | 3–5 min | Same structure-vs-price argument. |
| PX/NV event + market objects | every cycle | 3–5 min | Same. |

Implementation: the orchestrators' `price_sgps` gain optional injected
fetcher hooks (default = current direct fetch, used by dashboard shims;
`SGPService` injects cache-wrapped versions). Prices are **never**
cached.

Line-move staleness guard: a cached sel-id dict can briefly lack a *new*
line Kalshi just listed; the bot's existing per-accept staleness and
line-move gates already cover the trading risk, and the entry refreshes
within TTL.

### 3.6 Error handling

- **Per-book isolation:** one book hanging or erroring cannot delay the
  others (independent threads + per-book deadline) or crash the cycle
  (exceptions caught per book, logged, book contributes []).
- **Session death:** auto-reinit after 3 consecutive failures (port of
  existing scraper behavior into the service).
- **Probe safety:** probe is manual-only, has per-level call budgets, and
  hard-stops on 403/Akamai challenge for DK and on any anomaly for PX.
- **Rate-limit regression in production:** if a book starts 429ing at the
  shipped parallelism, the orchestrator already treats failed calls as
  "no price" (row simply absent). Mitigation is config: lower the
  parallelism constant. No new automatic backoff machinery in v1.

### 3.7 Testing & verification

- **Existing regression gate:** `mlb_sgp/tests/test_sgp_regression.py`
  golden-baseline comparison must pass after Phase 1 (parallelism must
  not change priced values).
- **New unit tests:**
  - DK/FD orchestrators with a mock client: parallel path returns the
    same row set as the sequential path on a fixture of targets
    (including integer-fallback targets and per-target failures).
  - `SGPService.refresh`: per-book deadline enforcement (slow mock book
    → [] for that book, others unaffected); client reinit after 3
    failures; TTL cache hit/miss/expiry behavior.
- **Live verification (pre-merge):**
  1. Run probe; record ceilings in the spec/README.
  2. One supervised `SGPService.refresh()` against live books from the
     worktree with a copied market DB; compare row counts + spot-check
     ~20 `sgp_decimal` values against a concurrent CLI-shim run.
  3. Timed before/after: full cycle wall-clock on a real slate.
  4. Dashboard parity: run the dashboard's SGP refresh path once
     (CLI shims) and confirm `mlb_sgp_odds` in `mlb_mm.duckdb` populates
     normally.
- **Success criterion:** cycle ≤ 60s on a ~15-game slate; no increase in
  per-book error rates at shipped parallelism over a 1-hour soak.

## 4. Version control

- **Branch / worktree:** `worktree-sgp-scraper-speedup` at
  `.claude/worktrees/sgp-scraper-speedup` (already created).
- **Commit structure:** one commit per phase —
  (0) probe harness, (1) DK/FD parallelism + regression-test pass,
  (2) `SGPService` + bot integration, (3) TTL caching. Spec committed
  first.
- **DuckDB:** never symlink; live verification uses copied DB files
  inside the worktree or runs from `main` post-merge.
- **Merge:** pre-merge executive review of `git diff main..HEAD`, then
  explicit user approval before merging to `main`. Worktree + branch
  deleted immediately after merge.
- **Bot restart** (from main repo cwd, not the worktree — restart-gotchas
  memory) after merge so both bots pick up the service.

## 5. Documentation updates (same merge)

- `mlb_sgp/README.md` — new "Concurrency" section: per-book parallelism
  constants, probe harness usage, measured ceilings.
- `kalshi_mlb_rfq/README.md` + `kalshi_mlb_mm/README.md` — SGP cadence
  section: subprocess model → in-process `SGPService`.
- Root `CLAUDE.md` — one-line update to the `kalshi_common/` bullet
  (`sgp_runner` now hosts `SGPService`).

## 6. Out of scope / future work

- Top-down bivariate-Poisson SGP pricer (separate spec; this work makes
  its calibration inputs fresher).
- Persistent-worker *processes* for the dashboard path.
- Automatic adaptive backoff / dynamic parallelism tuning.
- WebSocket or push-based book feeds.
