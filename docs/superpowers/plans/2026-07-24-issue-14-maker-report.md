# Daily State-of-the-Maker Report (issue #14) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A standalone read-only script `kalshi_mlb_mm/report.py` that aggregates the maker bot's three DuckDBs into a printed markdown "state of the maker" report: RFQ funnel + on-demand coverage, quotable universe, staleness, demand curve, settlement P&L, health.

**Architecture:** One module with a pure-compute layer (functions that query the DBs and return plain dicts/lists — unit-testable) and a thin render layer (dict → markdown string). All DB access goes through one lock-safe `_read()` helper copied from the monitor's pattern (`read_only=True`, 3 retries, `LOCKED` sentinel; missing file → `MISSING` sentinel; missing table → empty result). No section can crash the report — empty/missing/locked each degrade to a note line.

**Tech Stack:** Python stdlib + duckdb only. No bot-loop imports (`main.py` is never imported); only `kalshi_mlb_mm.config` (for default DB paths + `MIN_AGREEING_BOOKS`/`TARGET_ROI` context).

## Global Constraints

- **Read-only, always**: every `duckdb.connect` in report.py uses `read_only=True`. The report NEVER writes any DB.
- **Fresh-DB safe**: missing DB files, missing tables, and zero rows must all produce a readable report (acceptance criterion).
- **JSON access**: `json_extract_string(payload, '$.key')` — NEVER `->>` (documented gotcha at top of `kalshi_mlb_mm/research_queries.sql`).
- **Timestamps**: state + research DB columns are TIMESTAMPTZ → bind aware `datetime.now(timezone.utc)` cutoffs. Market DB (`mlb_sgp_odds.fetch_time`) is naive UTC → bind `now.replace(tzinfo=None)`.
- **Markout content dropped**: #13 was descoped 2026-07-24. Section 5 is settlement P&L only (from #12's `settlements` + `fills.realized_pnl`).
- **No bot restart, no monitor changes** — standalone script only.
- Tests: `python -m pytest kalshi_mlb_mm/tests/test_report.py -v`; full suite `python -m pytest kalshi_mlb_mm/tests/ kalshi_common/tests/` before review.

## Version control / worktree / docs

- Worktree: `.claude/worktrees/issue-14-maker-report`, branch `worktree-issue-14-maker-report`, based on main `28ac097` (includes #12). Already created.
- Files: create `kalshi_mlb_mm/report.py`, `kalshi_mlb_mm/tests/test_report.py`; modify `kalshi_mlb_mm/README.md` (new "Reports" section, same merge).
- Commit per task (`feat(mm): ...` / `test(mm): ...`); pre-merge executive review; **no merge without explicit user approval**; delete worktree + branch after merge.

## Data sources (established by code research)

| Question | Source |
|---|---|
| RFQs seen / in-scope | ~~state `seen_rfqs`~~ **CORRECTED during live smoke test:** `seen_rfqs` only records out-of-scope + live-quoted RFQs (main.py:1124 writes inside the quote transaction) — the funnel is derived rfq-level from `quote_decisions` instead |
| Decision funnel + skip reasons | state `quote_decisions` (`decision`, `reason`, `observed_at`, `rfq_id`) — vocabulary at `main.py:682-1541` |
| Grid vs on-demand vs skipped path | rfq-level: `last_decision LIKE 'out_of_scope%'` = skipped; any decision `reason='on_demand_pending'` = on-demand path; other in-scope = grid (approximation, footnoted: an on-demand RFQ whose fetch lands within its first 2s tick never logs `on_demand_pending`) |
| On-demand fetch latency/success | research `events`: pair `on_demand_requested` ↔ `on_demand_result` on `json_extract_string(payload,'$.leg_set_hash')`; per-book `latency_sec` inside `on_demand_result` payload `books` dict (keys via `json_keys`); hash with request but no result = failed flight |
| Quotable universe | market `mlb_sgp_odds`: per-day distinct `(game_id, combo, spread_line, total_line)` having `count(DISTINCT bookmaker) >= MIN_AGREEING_BOOKS`; per-book participation + last-row age |
| Staleness at quote | research: ASOF join `quote_priced` events to latest prior `scrape_done` event, age = ts diff |
| Staleness at accept | research `confirm_singles_check` payload `quote_age_sec` |
| Demand curve | state `quote_decisions` (quotes: margin = `blended_fair - yes_bid` / `(1-blended_fair) - no_bid`) vs `fills` (accepts: margin = side-adjusted `blended_fair_at_quote - price`) |
| P&L | state `fills.realized_pnl` + `settlements` (reuse research_queries.sql §7-8 math) |
| Health | voids/confirms + `sweep_cancel` reasons + `circuit_breaker` + `halted_high_void_rate` (state); `halt_event`, `reconcile_done` outcomes (research); phantom = fills `contracts=0 AND reconciled` |

---

### Task 1: Plumbing — lock-safe reader, sentinels, CLI skeleton, empty-report smoke test

**Files:**
- Create: `kalshi_mlb_mm/report.py`
- Test: `kalshi_mlb_mm/tests/test_report.py`

**Interfaces (produced, used by all later tasks):**
- `LOCKED`, `MISSING` — module-level sentinels (falsy classes, like monitor `queries.py:26-33`)
- `_read(db_path, sql, params=None, retries=3) -> list[tuple] | LOCKED | MISSING` — read-only + retry; `CatalogException` (missing table) → `[]`
- `Windows = namedtuple("Windows", "now h24 d7 d7_naive h24_naive")` — precomputed cutoffs
- `build_report(state_db, research_db, market_db, now=None) -> str` — assembles all sections
- `main(argv=None)` — argparse: `--state-db/--research-db/--market-db` (default `config.DB_PATH`/`config.RESEARCH_DB_PATH`/`config.MARKET_DB`), prints `build_report(...)`
- Test helper `_setup_dbs(tmp_path, monkeypatch) -> (state, research, market)` — state: monkeypatch `cfg.DB_PATH` + `importlib.reload(db)` + `db.init_database()` (convention from `tests/test_settlement.py:19-25`); research: monkeypatch `cfg.RESEARCH_DB_PATH` + `research.init_research_db()`; market: `mlb_sgp.db.ensure_table(db_path=...)`

**Steps:**
- [ ] Write failing tests: `test_empty_dbs_report_runs` (fresh schemas → `build_report` returns str containing "State of the Maker" and "no fills yet", no raise), `test_missing_db_files` (paths that don't exist → report still returns, contains "not found"), `test_read_locked_returns_sentinel` (hold a write conn on the state DB, `_read` → `LOCKED`).
- [ ] Run: `python -m pytest kalshi_mlb_mm/tests/test_report.py -v` → FAIL (module missing).
- [ ] Implement report.py skeleton: docstring (inputs/outputs/side-effects: READ-ONLY, writes nothing), sentinels, `_read`, `Windows`, `build_report` with placeholder sections, `main`, `if __name__ == "__main__"`.
- [ ] Tests pass; commit `feat(mm): report.py skeleton — lock-safe read-only plumbing (#14)`.

### Task 2: Section 1 — RFQ funnel (24h + 7d) + path split (grid / on-demand / skipped)

**Interfaces:**
- `funnel_counts(state_db, since) -> dict` — keys: `seen, in_scope, quoted, accepted, confirmed, filled, decisions` (list of (decision, reason, n)), `paths` = `{skipped_out_of_scope, on_demand, grid}` each `{rfqs, quoted}`
- Constants: `QUOTED_DECISIONS = ("quoted", "dry_run_quote")`, `ACCEPT_DECISIONS = ("confirmed", "voided_no_legs", "voided_singles_moved", "voided_no_fresh_books", "voided_last_look")`
- `render_funnel(d24, d7) -> str`

**Steps:**
- [ ] Failing test `test_funnel_math`: seed 5 `seen_rfqs` (3 in-scope), decisions: 2×`skipped/out_of_scope`, rfq A `skipped/on_demand_pending` then `quoted`, rfq B `quoted` then `confirmed`, rfq C `skipped/no_fair`. Assert seen=5, in_scope=3, quoted=2, accepted=1, confirmed=1; paths: skipped=2 rfqs, on_demand=1 rfq (1 quoted), grid=2 rfqs (1 quoted).
- [ ] Implement (rfq-level path CTE over `quote_decisions` + `seen_rfqs.last_decision`); render 24h/7d side-by-side markdown table + full decision/reason breakdown table + approximation footnote.
- [ ] Tests pass; commit.

### Task 3: Section 1b — on-demand coverage: fetch success rate + latency

**Interfaces:**
- `on_demand_stats(research_db, since) -> dict` — `{flights_requested, hashes_requested, hashes_landed, success_rate, wall_latency_p50/p95, per_book: [(book, landings, latency_p50, latency_p95)]}`

**Steps:**
- [ ] Failing test: seed research `events` — `on_demand_requested` for hashes h1,h2,h3; `on_demand_result` for h1 (books dk latency 1.2 / fd 2.0, +5s wall) and h2 (dk 0.8, +3s). Assert hashes_requested=3, landed=2, success 2/3; dk latency p50 = 1.0; wall p50 = 4.0.
- [ ] Implement: hash pairing via min(result.ts) − max(request.ts ≤ result.ts); per-book via `json_keys(payload, '$.books')` unnest + `json_extract_string(payload, '$.books.' || book || '.latency_sec')` cast DOUBLE; note in output that confirm-lane re-fetch latency is not persisted (voids show as `voided_no_fresh_books` in section 6).
- [ ] Tests pass; commit.

### Task 4: Section 2 — quotable universe

**Interfaces:**
- `universe_stats(market_db, since_naive) -> dict` — `{per_day: [(date, combos_passing, combos_total)], per_book: [(book, days_present, distinct_combos, last_age_min)]}` using `config.MIN_AGREEING_BOOKS`

**Steps:**
- [ ] Failing test: seed `mlb_sgp_odds` across 2 days — day1: combo X with 3 books, combo Y with 1 book; day2: combo X with 2 books. Assert per_day = [(d1, 1, 2), (d2, 1, 1)]; per-book rows for all books.
- [ ] Implement (`GROUP BY game_id, combo, spread_line, total_line` — DuckDB GROUP BY treats NULLs as equal, so ML×total's NULL spread_line is safe); render.
- [ ] Tests pass; commit.

### Task 5: Section 3 — staleness

**Interfaces:**
- `staleness_stats(research_db, since) -> dict` — `{quote_age: {n, p50, p90, p95, max}, accept_age: {n, p50, p90, p95, max}}` (seconds)

**Steps:**
- [ ] Failing test: seed `scrape_done` at T, T+120; `quote_priced` at T+30 (age 30) and T+150 (age 30); `confirm_singles_check` with `quote_age_sec` 5 and 45. Assert quote_age p50=30, accept n=2, p50=25.
- [ ] Implement: ASOF JOIN (`FROM q ASOF LEFT JOIN s ON q.ts >= s.ts`) for quote-time age; `CAST(json_extract_string(payload,'$.quote_age_sec') AS DOUBLE)` for accept-time; `quantile_cont` for percentiles; quotes with no prior scrape_done reported as `unknown_age` count.
- [ ] Tests pass; commit.

### Task 6: Section 4 — demand curve (quotes vs accepts by margin × fair band)

**Interfaces:**
- `demand_stats(state_db, since) -> dict` — `{cells: [(margin_bucket, fair_band, quotes, accepts)], quote_sides, fill_count}`
- Buckets: margin `(-inf,1c),[1,2),[2,3),[3,4),[4,5),[5c,inf)`; fair band `<.10,[.10,.25),[.25,.50),[.50,.75),[.75,.90),>=.90`
- Quote side-margins from `quote_decisions` rows in `QUOTED_DECISIONS` (two rows per decision: yes margin `blended_fair - yes_bid` at fair `blended_fair`; no margin `(1-blended_fair) - no_bid` at fair `1-blended_fair`; NULL bids skipped). Accepts from `fills` with `contracts > 0`: side-adjusted fair = `blended_fair_at_quote` if `side_held='yes'` else `1-blended_fair_at_quote`; margin = side-adjusted fair − `price`.

**Steps:**
- [ ] Failing test: seed one quoted decision (fair .60, yes_bid .57, no_bid .37) → yes side lands (3c,[.50,.75)), no side ([2,3)c,[.25,.50)); one fill (side yes, price .57, fair_at_quote .60) → accept in (3c,[.50,.75)). Assert those cells.
- [ ] Implement + render pivot table (margin rows × fair-band cols, `q/a` cells) with note "feeds #26".
- [ ] Tests pass; commit.

### Task 7: Section 5 — settlement P&L (markouts descoped)

**Interfaces:**
- `pnl_stats(state_db) -> dict` — all-time: `{fills, contracts, settled_fills, unsettled_fills, realized_pnl_total, avg_quoted_margin_per_ct, realized_per_ct, by_result: [(result, fills, pnl)]}`

**Steps:**
- [ ] Failing tests: (a) fills exist, none settled → dict flags `settled_fills=0`, render says "fills exist but none settled yet"; (b) seeded settled fill (contracts 10, price .55, blended_fair_at_quote .58, fee .01, realized_pnl from win) + `settlements` row → totals + per-result row + quoted-vs-realized per-contract match hand-computed values (query-8 math: quoted margin/ct = `fair_at_quote − price`; realized/ct = `sum(realized_pnl)/sum(contracts)`); (c) empty fills → "no fills yet".
- [ ] Implement (phantom fills `contracts=0` excluded); render.
- [ ] Tests pass; commit.

### Task 8: Section 6 — health

**Interfaces:**
- `health_stats(state_db, research_db, since) -> dict` — `{voids_by_decision, confirmed, void_rate, sweep_cancels_by_reason, circuit_breakers, void_rate_halts, phantom_fills, reconcile_outcomes}`

**Steps:**
- [ ] Failing test: seed decisions (2 confirmed, 1 `voided_singles_moved`, 1 `sweep_cancel/tipoff`, 1 `halted_high_void_rate`), research `halt_event` (fire), `reconcile_done` (`outcome='phantom'`), phantom fill row (`contracts=0, reconciled=TRUE`). Assert void_rate = 1/3, each count.
- [ ] Implement (void_rate = voids / (voids + confirmed), guard div-zero); render.
- [ ] Tests pass; commit.

### Task 9: README + full suite + review

**Files:** Modify `kalshi_mlb_mm/README.md`

**Steps:**
- [ ] Add "Reports" section: what the report answers, `python -m kalshi_mlb_mm.report` usage (run from main repo cwd for live DBs; `--state-db` etc. overrides), read-only guarantee, section list, the on-demand-path-classification approximation, "markouts descoped (#13) — settlement P&L is the adverse-selection signal".
- [ ] Run `python -m kalshi_mlb_mm.report` against the live DBs (read-only — safe) as a smoke test; eyeball output.
- [ ] Full suite: `python -m pytest kalshi_mlb_mm/tests/ kalshi_common/tests/` → all green.
- [ ] Commit; pre-merge executive review (data integrity / resource safety / edge cases / dead code / hygiene / security); present findings; **ask user for merge approval**.

## Self-review notes

- Spec coverage: all 6 issue sections mapped (funnel T2+T3, universe T4, staleness T5, demand T6, P&L T7, health T8); acceptance criteria: fresh-DB safety (T1), funnel-math + out-of-scope-share fixture tests (T2), read-only + lock-safe (T1), README (T9).
- Type consistency: all sections consume `_read` sentinels; renders check `is LOCKED / is MISSING` before shaping.
- Known approximations, stated in report output: grid-vs-on-demand path split heuristic; quote-time staleness proxied via scrape_done cadence (per-book age not persisted at quote time); confirm-lane re-fetch latency not persisted.
