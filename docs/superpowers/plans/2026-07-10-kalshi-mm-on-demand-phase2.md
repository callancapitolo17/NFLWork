# Kalshi MLB MM Phase 2 — On-Demand Pricing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Price and quote any same-game leg combo the six sportsbooks will price, via live on-demand SGP queries at RFQ time (spec: `docs/superpowers/specs/2026-07-08-kalshi-mlb-mm-on-demand-pricing-design.md`, rev 5).

**Architecture:** Pure math (partition devig + correlation transfer with Fréchet bounds) in `kalshi_common/fair_value.py`; per-book leg resolution extracted from the existing orchestrator closures in `mlb_sgp/*`; a single-combo entry point on `SGPService`; a queue/worker/result-store engine in `kalshi_mlb_mm/on_demand.py`; a one-argument seam in `router.subcombo_fair`; discovery/confirm integration + instrumentation in `main.py`. The open-RFQ poll drives the fetch feed.

**Tech Stack:** Python 3, pytest, pandas, curl_cffi sessions (existing per-book clients), DuckDB (state), no new dependencies.

## Global Constraints

- Grid pricing must stay byte-identical: `kalshi_mlb_mm/tests/test_router_integration.py` (2-leg == median(_book_fairs) to 1e-9) must pass UNMODIFIED after every task.
- No new config keys. Module constants only: `QUOTE_FRESH_SEC = 15` (on_demand.py), `PARTITION_OVERROUND_PER_LEG = 0.12` (fair_value.py).
- Fail-safe everywhere: any doubt → return None → no quote. Never raise into the trading loop.
- Existing orchestrator paths (`price_sgps`) byte-identical — new functions only, no edits to existing function bodies unless a task explicitly says so.
- Lone singles stay excluded: the RFQ-level `len(canon) >= 2` rule in main.py is retained.
- Run the FULL suite (`python3 -m pytest kalshi_common/tests kalshi_mlb_mm/tests mlb_sgp/tests -q`) before every commit; all green.
- Commit after every task with the message given in the task.
- Cell ordering convention (used by legset, fair_value, sgp_service — must agree): cell index `i` in `range(2**n)`; bit `j` of `i` (LSB = leg 0) says leg `j` is FLIPPED to its opposite side; cell 0 = the target (all chosen sides).

---

### Task 1: `legset.enumerate_partition` (pure)

**Files:**
- Modify: `kalshi_common/legset.py` (append)
- Test: `kalshi_common/tests/test_legset.py` (append)

**Interfaces:**
- Produces: `enumerate_partition(game_legs: list[CanonicalLeg]) -> list[list[CanonicalLeg]]`, cell ordering per Global Constraints; `flip_leg(leg: CanonicalLeg) -> CanonicalLeg`.

- [ ] Step 1: failing tests — flip semantics (spread/ml home↔away, total over↔under, line unchanged), 2^N count and ordering (cell 0 == input, cell index bit j flips leg j), determinism, same-market-both-sides input still enumerates (impossible cells are the caller's problem).
- [ ] Step 2: run, verify fail (`AttributeError`/`ImportError`).
- [ ] Step 3: implement:

```python
_FLIP = {"home": "away", "away": "home", "over": "under", "under": "over"}

def flip_leg(leg: CanonicalLeg) -> CanonicalLeg:
    return CanonicalLeg(leg.game_id, leg.market_type, leg.line, _FLIP[leg.side])

def enumerate_partition(game_legs):
    n = len(game_legs)
    out = []
    for i in range(2 ** n):
        out.append([flip_leg(l) if (i >> j) & 1 else l
                    for j, l in enumerate(game_legs)])
    return out
```

- [ ] Step 4: tests pass; full suite green.
- [ ] Step 5: commit `feat(mm-p2): legset.enumerate_partition — 2^N side-combination cells`.

### Task 2: `fair_value.devig_partition` (pure, Route A)

**Files:**
- Modify: `kalshi_common/fair_value.py` (append; add `PARTITION_OVERROUND_PER_LEG = 0.12`)
- Test: `kalshi_common/tests/test_fair_value_on_demand.py` (create)

**Interfaces:**
- Consumes: `_probit_devig_n(raw_probs)` (existing).
- Produces: `devig_partition(cell_decimals: list[float], n_legs: int) -> float | None` — index 0 = target cell.

- [ ] Step 1: failing tests — (a) REGRESSION BRIDGE: for a 4-cell grid, result equals `devig_book` on the same rows to 1e-12 (build a tiny DataFrame with 4 combos, call `devig_book(rows, combo=<cell0>)`, compare); (b) overround > 1 + 0.12·n → None; (c) overround < 1.0 → None; (d) any cell decimal ≤ 1.0 or None/NaN → None; (e) wrong cell count (≠ 2^n) → None.
- [ ] Step 2: run, verify fail.
- [ ] Step 3: implement:

```python
PARTITION_OVERROUND_PER_LEG = 0.12

def devig_partition(cell_decimals, n_legs):
    if not cell_decimals or len(cell_decimals) != 2 ** n_legs:
        return None
    try:
        raw = [1.0 / float(d) for d in cell_decimals]
    except (TypeError, ValueError, ZeroDivisionError):
        return None
    if any(not (0.0 < p < 1.0) for p in raw):
        return None
    s = sum(raw)
    if not (1.0 <= s <= 1.0 + PARTITION_OVERROUND_PER_LEG * n_legs):
        return None
    return float(_probit_devig_n(raw)[0])
```

- [ ] Step 4: tests pass; full suite green.
- [ ] Step 5: commit `feat(mm-p2): devig_partition — N-aware overround gate + probit over 2^N cells`.

### Task 3: `fair_value.fair_by_correlation_transfer` + route selector (pure, Route B)

**Files:**
- Modify: `kalshi_common/fair_value.py` (append)
- Test: `kalshi_common/tests/test_fair_value_on_demand.py` (append)

**Interfaces:**
- Produces:
  - `devig_two_way(dec_a: float, dec_b: float) -> tuple[float, float] | None` — exact 2-way probit devig of a single market's two sides (reuse `_probit_devig_n` on 2 cells).
  - `fair_by_correlation_transfer(sgp_decimal: float, singles: list[tuple[float, float]]) -> float | None` — `singles[i] = (vigged_implied_chosen_side, devigged_fair_chosen_side)`; returns `∏fair × (1/sgp_decimal)/∏vigged`, then Fréchet gate `max(0.0, sum(fairs) - (n-1)) <= result <= min(fairs)`, and `0 < result < 1`, else None.
  - `book_on_demand_fair(cells: list[float] | None, sgp_decimal: float | None, singles: list | None, n_legs: int) -> tuple[float, str] | None` — try `devig_partition(cells, n_legs)` when `cells` is complete → `(fair, "partition")`; else `fair_by_correlation_transfer` → `(fair, "transfer")`; else None.
- [ ] Step 1: failing tests — vig cancellation (synthetic book: known fair joint × uniform margin on SGP and singles → recovers fair within 1e-6); Fréchet upper rejection (nested totals: singles fairs .55/.30, transfer result > .30 → None); Fréchet lower rejection; high-correlation multiplier > 2 ACCEPTED when within Fréchet (RL+ML stack case); route selector A-preferred, B-fallback, None when neither.
- [ ] Step 2: run, verify fail.  Step 3: implement per interface.  Step 4: green + full suite.
- [ ] Step 5: commit `feat(mm-p2): Route B correlation transfer with Fréchet bounds + route selector`.

### Task 4: shared on-demand types

**Files:**
- Modify: `mlb_sgp/_shared.py` (append)
- Test: `mlb_sgp/tests/test_on_demand_types.py` (create)

**Interfaces (produced — every later task consumes these EXACT names):**

```python
@dataclass(frozen=True)
class ResolvedLeg:
    ref: object                       # book-specific descriptor, opaque
    opposite_ref: object | None       # None => one-sided at this book
    single_decimal: float | None      # chosen-side single odds from structure
    opposite_decimal: float | None

@dataclass(frozen=True)
class OnDemandBookResult:
    book: str
    fair: float
    route: str        # "partition" | "transfer"
    n_cells_priced: int
    latency_sec: float
```

- [ ] Steps: failing test (construct, frozen, defaults) → implement → green → commit `feat(mm-p2): shared ResolvedLeg/OnDemandBookResult types`.

### Tasks 5–10: per-book `resolve_legs` + `price_selection_set` (one task/commit per book)

Common shape for each book task (books in order: DK, FD, PX, NV, MGM, CZR):

**Files (per book B):**
- Modify: `mlb_sgp/<B>.py` (append two functions; do NOT edit `price_sgps`)
- Modify (FD only): `mlb_sgp/scraper_fanduel_sgp.py::fetch_event_runners` — additionally capture each runner's odds into the returned tuples (extend `(marketId, selectionId)` → `(marketId, selectionId, decimal_or_None)`; update its two existing consumers in `fanduel.py` to unpack 3-tuples — behavior otherwise unchanged; verify runner odds field name against a recorded fixture: expected `winRunnerOdds.trueOdds.decimalOdds.decimalOdds` per the shape at `scraper_fanduel_sgp.py:482`)
- Test: `mlb_sgp/tests/test_<B>_on_demand.py` (create; fixtures under `mlb_sgp/tests/fixtures/` following the existing fixture pattern in that directory)

**Interfaces (produced, identical across books):**
- `resolve_legs(structure, legs: list[CanonicalLeg], home_team: str, away_team: str) -> list[ResolvedLeg] | None` — `structure` is that book's existing per-event structure object (DK: `fetch_selection_ids(...)["fg"]` dict; FD: `fetch_event_runners(...)["fg"]`; PX: market list from `fetch_event_markets`; NV: `EventLegs` fg bucket; MGM: `parse_markets(...)["FG"]`; CZR: `parse_markets(event)["FG"]`). None if any leg's CHOSEN side is unresolvable; a missing OPPOSITE side yields `opposite_ref=None`.
- `price_selection_set(client, refs: list) -> float | None` — one wire call; a 1-element list prices a single.

**Key-mapping facts (extract, don't reinvent — mirror these exact sites):**
- DK: spreads keyed `(sign, abs_line, participant)` with home participant `"1"`, away `"3"` — mirror `scraper_draftkings_sgp.py:194-195` for sign construction; totals keyed `("O"|"U", line)`; ML bucket `{"1": [...], "3": [...]}` (`draftkings.py:455-461`). Wire call: generalize `post_calculate_bets(session, sel_ids: list)` from `scraper_draftkings_trifecta.py:104` (move/import, don't duplicate). Singles: structure has NO odds → `single_decimal=None`; Route B singles priced via 1-leg `price_selection_set` (2 calls/leg) inside Task 11, DK only. 422/`combinabilityRestrictions` → None (mirror `calculate_sgp` handling at `scraper_draftkings_sgp.py:613-636`).
- FD: runners keyed `("home"|"away", abs_line)` / `("O"|"U", line)` / `moneyline` dict (`scraper_fanduel_sgp.py:303` docstring); wire = implyBets with N `betLegs` (generalize the 2-leg body at `scraper_fanduel_sgp.py:452`).
- PX: legs from market tree via `_selection_to_leg` (`prophetx.py:695`) / `_pick_ml_selection` (`prophetx.py:478`); singles odds via `_safe_leg_decimal` pattern (`sel["odds"]`, american → decimal); wire = `client.submit_parlay_rfq(legs)` → `(offer, used_fallback)`; decimal from the offer node as in the existing orchestrator; used_fallback=True → treat as priced (teaser filter already applied).
- NV: fg bucket keys `home_spread/away_spread/over/under/home_ml/away_ml` per line (`scraper_novig_sgp.py:440-459`); legs carry `available` implied prob → `single_decimal = 1/available`; wire = `client.submit_parlay([uuid, ...])["decimal"]`.
- MGM: `spreads: {home_signed_line: {"home": (mid,oid,dec), "away": ...}}`, `totals: {line: {"over": ..., "under": ...}}`, `moneyline` (`betmgm.py:95` docstring); wire = `client.price_picks(fixture_id, [(mid,oid), ...])["decimal"]`.
- CZR: `spreads: {home_line: {"home": leg, "away": leg}}`, `totals`, `moneyline` — leg dicts ready for `/bets/details`; AWAY spread line must be NEGATED (`caesars.py:42-50`); singles from `_price_d(leg)` (`caesars.py:105`); wire = `client.price_combo(legs)["decimal"]`.

**Steps per book:** failing fixture tests (chosen-side hit; chosen-side miss → None; opposite-side miss → `opposite_ref=None`; singles odds captured where the book carries them; canonical→book line sign conventions incl. away-negation cases) → implement → green + full suite → commit `feat(mm-p2): <book> on-demand leg resolution + N-leg pricing`.

### Task 11: `SGPService.price_on_demand`

**Files:**
- Modify: `kalshi_common/sgp_service.py` (append method + per-book dispatch)
- Test: `kalshi_common/tests/test_sgp_service_on_demand.py` (create)

**Interfaces:**
- Consumes: Tasks 1–10 (`enumerate_partition` ordering = bitmask; `book_on_demand_fair`; per-book `resolve_legs`/`price_selection_set`; existing `_BookState`, `_book_done`, structure fetchers/TTL caches).
- Produces: `price_on_demand(book: str, game: GameRef, legs: list[CanonicalLeg]) -> OnDemandBookResult | None` with `GameRef = namedtuple("GameRef", "game_id home_team away_team commence_time")` (define in `_shared.py` here).
- Logic: resolve event (reuse each book's existing event-match helpers, cached) → structure (TTL cache) → `resolve_legs` → if `len(legs) <= 3` and all `opposite_ref` present: price cells by bitmask product over (ref, opposite_ref) — cell 0 first, stop on first missing → else/fallback: price target refs once; singles = structure odds where present, else (DK) 1-leg calls both sides; `devig_two_way` per leg (vig-fallback haircut when one-sided single, mirror `devig_book`'s `<4 sides` fallback with the book's existing `*_VIG_FALLBACK`) → `book_on_demand_fair`. Failure accounting: success/failure feeds `_book_done` exactly like `refresh`. NO cross-lock with the sweep; runs on the caller's thread.
- [ ] Steps: failing tests with fake resolver/pricer injected via a `_on_demand_hooks` test seam dict (mirroring the `runners` seam) — Route A happy path (8 cells), missing cell → Route B, one-sided → Route B, DK singles-via-1-leg path, book failure increments failures → green + full suite → commit `feat(mm-p2): SGPService.price_on_demand single-combo entry point`.

### Task 12: `OnDemandEngine`

**Files:**
- Create: `kalshi_mlb_mm/on_demand.py`
- Test: `kalshi_mlb_mm/tests/test_on_demand_engine.py`

**Interfaces:**
- Consumes: `SGPService.price_on_demand`, `legset.leg_set_hash`.
- Produces (main.py consumes exactly these):
  - `OnDemandEngine(service, now_fn=time.monotonic)`
  - `ensure_fetch(hash_: str, game: GameRef, legs: list[CanonicalLeg]) -> None` (idempotent while in flight)
  - `lookup(hash_: str) -> dict[str, float] | None` (fresh ≤ `QUOTE_FRESH_SEC=15` only; also `lookup_results(hash_)` returning `{book: OnDemandBookResult}` for research)
  - `refetch_now(jobs: list[tuple[str, GameRef, list[CanonicalLeg]]], deadline_sec: float) -> bool` (priority lane; awaits in-flight rather than duplicating; True iff ALL hashes land fresh in time)
- Internals: unbounded FIFO + in-flight set; ONE consumer thread (daemon), lazily (re)started on enqueue; per combo, thread-per-book fan-out with `per_book_deadline_sec`; at most one combo in flight per book (books busy with the current combo — single consumer gives this for free); store `{hash: (landed_at, {book: OnDemandBookResult})}`; single `threading.Lock`; every worker exception caught+logged.
- [ ] Steps: failing tests with a fake service (controllable latency/results, fake clock) — enqueue/dedup, landing → lookup fresh, 15s death, feed invariant (ensure_fetch after death re-enqueues; no poll → no fetches), refetch_now all-land / partial / deadline overrun / awaits in-flight, worker-death lazy restart, no stale served ever → green + full suite → commit `feat(mm-p2): OnDemandEngine — queue/worker/result store, poll-driven feed`.

### Task 13: router seam + `consensus_detail`

**Files:**
- Modify: `kalshi_mlb_mm/router.py` (`subcombo_fair`/`combo_fair` gain `on_demand_fairs=None` kw-only arg; new `consensus_detail`)
- Test: `kalshi_mlb_mm/tests/test_router_on_demand.py` (create)

- [ ] Steps: failing tests — on_demand route with injected lookup → consensus over returned fairs; lookup None/missing → None; grids ignore the arg; `consensus_detail(book_fairs, min_books, band) -> tuple[float, list[str]] | None` agrees with `consensus` on the fair to 1e-12 → implement (spec §4.6 code) → **run `test_router_integration.py` UNMODIFIED — must pass** → full suite → commit `feat(mm-p2): router on_demand seam (grid path byte-identical) + consensus_detail`.

### Task 14: main.py integration + instrumentation

**Files:**
- Modify: `kalshi_mlb_mm/main.py`
- Test: `kalshi_mlb_mm/tests/test_main_on_demand.py` (create; follow `test_main_smoke.py` monkeypatch style)

Changes (each with a test):
1. `_priceable_in_phase1` → `_priceable`: allow `on_demand`; keep `len(canon) >= 2`.
2. First: READ `_resolve_game_for_legs` and the target-lines/game cache to build `GameRef` (home/away/commence from `mlb_target_lines` via existing cache; if absent for a game → treat as unresolvable, skip `no_game`).
3. Discovery: per on-demand game sub-combo — `engine.lookup(hash)`; fresh → price via `router.combo_fair(..., on_demand_fairs=engine.lookup)` (hysteresis/replace path untouched); else `engine.ensure_fetch(...)` + `_log_decision("skipped", reason="on_demand_pending")` + research `on_demand_requested` (first enqueue per hash per flight only).
4. Confirm: combos with on-demand games → `deadline = 30s window remainder − 5s buffer` → `engine.refetch_now(jobs, deadline)`; False → void via existing `voided_no_fresh_books` block; True → `cur_fair = router.combo_fair(..., on_demand_fairs=engine.lookup)`.
5. Risk sweep: pass `on_demand_fairs=engine.lookup` in its `combo_fair` call (works when fresh, no-ops otherwise — zero structural change).
6. Fill research payload: add `consensus_books` (via `consensus_detail`) + per-book `route` (from `lookup_results`).
7. Out-of-scope instrumentation: store `legs_json` + granular `last_decision` (`out_of_scope_non_mlb`, `out_of_scope_lone_single`, `out_of_scope_unparseable`) in the `seen_rfqs` insert at `main.py:469-476`.
8. Research event `on_demand_result` emitted by the discovery tick when a landing is first consumed (per-book route/fair/latency/cells, `route_gap` where both routes free).
- [ ] Steps: failing tests (fake engine + fake gateway; pending→landed→quoted; confirm void on refetch False; instrumentation rows) → implement → green + full suite → commit `feat(mm-p2): wire on-demand engine into discovery/confirm/risk + instrumentation`.

### Task 15: end-to-end integration test

**Files:**
- Test: `kalshi_mlb_mm/tests/test_on_demand_e2e.py` (create)

- [ ] Simulated loop with fake service (fake clock): RFQ appears → pending → landing → quoted; fair moves → replace beyond hysteresis / hold within; result ages → re-fetch fires only while RFQ in poll; RFQ leaves poll → zero further fetches; accept → refetch_now → confirm/void both paths; multi-game (grid + on-demand mix) product pricing. Commit `test(mm-p2): on-demand end-to-end loop invariants`.

### Task 16: docs + memory + finish

**Files:**
- Modify: `kalshi_mlb_mm/README.md` (architecture, two routes, feed model, failure modes, rollout/retreat playbook)
- Modify: `CLAUDE.md` (maker bullet: Phase 2 summary)
- Memory: `kalshi_mm_generalized_combos.md` status update

- [ ] Full suite; spec cross-check (every §4 item has code + tests); commit `docs(mm-p2): Phase 2 on-demand pricing docs`. Then: live REPL verification per book (one exotic combo per route; DK 1-leg calculateBets check; FD odds-field check; PX 1-leg RFQ check) — findings recorded in the PR/summary, NOT merged on failure. Adversarial review. STOP — user approval required before merge/push.
