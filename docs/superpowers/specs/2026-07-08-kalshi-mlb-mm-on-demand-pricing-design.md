# Kalshi MLB MM — Phase 2: On-Demand Pricing of Arbitrary Leg Combos

**Date:** 2026-07-08
**Status:** Draft — awaiting user review
**Depends on:** Phase 1 generalized combos (merged to main 2026-07-08, `2e9a0e8`);
spec at `docs/superpowers/specs/2026-06-28-kalshi-mlb-mm-generalized-combos-design.md`

---

## Review Pack

**What we're building** — The maker currently skips any same-game leg combo
that isn't one of the two pre-scraped 2-leg grids (spread×total, ml×total).
Phase 2 makes the quoting engine fully general: when an RFQ contains a novel
same-game shape (3+ legs, spread+ml, total+total, …), the bot live-queries the
six sportsbooks' SGP endpoints for that exact leg set, devigs, builds
consensus, and quotes it — asynchronously, so the 2-second RFQ discovery loop
never blocks on a book. Goal per user direction: "quote any combination of
legs given we have access to it through our sportsbook scraping."

**Key decisions**

1. **Async fetch, quote next tick — never block the 2s discovery loop.**
   A novel shape enqueues a background fetch (~10–20s) and the RFQ is skipped
   this tick; a later tick prices it from cache. Rejected: blocking the
   discovery tick on live HTTP. Evidence async works: the bot already first
   acts on RFQs at median age 15s (p90 29s) and those quotes succeed — RFQs
   rest minutes.
2. **Full 2^N-partition devig, no fallback.** Fetch every side-combination of
   the leg set from each book (3 legs → 8 cells), probit-devig across the
   partition, read the target cell — the exact generalization of Phase 1's
   4-cell grid rule ("full grid or nothing"). Rejected: single-price + vig
   haircut (guessy, inconsistent with Phase 1). Leg cap: 3 (config).
3. **All 6 books participate** (user decision). Caesars (WAF fragility) and
   ProphetX (an on-demand query is a real RFQ on their exchange) were offered
   as exclusions and declined; both get failure-isolation + measurement.
4. **Plain 2-book consensus, same as grids** (user decision). A DK+Novig pair
   counts even though Novig mirrors DK. The stricter independence rule was
   offered and declined; fills get a `dk_nv_only` research flag so the risk is
   measured, not silent.
5. **Router stays pure; grid pricing byte-identical.** `subcombo_fair` gains
   one optional `on_demand_fairs=None` lookup argument. Default None
   reproduces Phase 1 exactly; the 1e-9 regression lock is untouched. All
   stateful machinery (queue, worker, cache) lives in a new
   `kalshi_mlb_mm/on_demand.py`.

**Risks / push back here**

- **Target flow is currently zero.** Live measurement (June 9 → July 8): 0
  novel same-game shapes in ~1,100 classified RFQ markets (exhaustive last-7-
  days + 1,500-market random sample; 95% upper bound ≈ 0.3% of out-of-scope
  flow). You chose to build the general engine anyway; instrumentation ships
  in the same merge so usage is visible from day one. If Kalshi's MVE builder
  structurally can't produce these shapes, this code may never fire.
- **Caesars blast radius.** On-demand shares Caesars' WAF-token session with
  the 60s grid sweep. If on-demand traffic trips the WAF, the *grid* pipeline
  loses Caesars too. Mitigation: per-book failure isolation + a per-book
  on-demand kill switch (`ON_DEMAND_BOOKS` config), but the residual risk is
  yours by choice.
- **DK+Novig-only consensus on correlation-heavy shapes** is the known
  "fill our fair-error tail" vector from the adversarial review — exotic
  combos are where one book's correlation model being wrong hurts most, and a
  DK+NV pair is ~one independent source. Accepted; measured via research flag.
- **In-scope quoting bottleneck is elsewhere.** All-time skip reasons:
  `no_fair` 5,431 vs 518 quoted. If the goal is more fills, book
  freshness/consensus tuning on already-supported shapes is the larger lever.

**Worth understanding** (opt-in)

1. *Joint-distribution partitions.* A 2-leg combo's 4-cell grid is a 2×2
   contingency table whose cells sum to 1. N legs → `expand.grid` of each
   leg's two sides = 2^N mutually exclusive cells. Devigging over the full
   partition is what removes the book's margin without assuming where it hides.
2. *TTL caches.* Like R's `memoise::memoise(..., cache = cache_mem(max_age))`:
   a fetched price is reused until its age exceeds the freshness bound, then
   it's as if it were never fetched. Two different TTLs here: "may I quote off
   this?" (180s) vs "should I re-fetch it?" (60s, only for combos with open
   quotes).
3. *Producer/consumer worker.* The discovery tick *produces* fetch jobs onto a
   queue and moves on; a background thread *consumes* them and fills the
   cache. The two sides only communicate through the queue and cache — no
   shared mutable pricing state, which is what keeps the 2s tick safe.

---

## 1. Context and evidence

### 1.1 Where Phase 1 stops

`legset.classify_subcombo` routes each game's legs to `single`,
`grid_spread_total`, `grid_ml_total`, `on_demand`, or `unpriceable`.
`router.subcombo_fair` returns `None` for `on_demand` (`router.py:130`), and
`main._priceable_in_phase1` (`main.py:250`) marks the whole RFQ out-of-scope —
it never even reaches the research firehose (`rfq_received` fires post-filter,
`main.py:478`), and `seen_rfqs.legs_json` is NULL for out-of-scope rows. Phase
2 removes the pricing gap *and* the measurement blind spot.

### 1.2 Live flow measurement (2026-07-08, read-only)

| Sample | Markets classified | Novel same-game shapes |
|---|---|---|
| Exhaustive: out-of-scope RFQ markets, last 7 days | 177 (189 RFQs) | 0 |
| Random 1,500 of all-time out-of-scope (381,315 markets since 2026-06-09) | 911 (589 delisted → unfetchable) | 0 |
| All 8,140 in-scope `rfq_received` RFQs | 8,140 | 0 (by scope construction) |

Out-of-scope composition: ~85% contain non-MLB legs (cross-category MVE
collections carry every sport), ~15% are shapes Phase 1 now prices (seen
before its deploy). Zero lone singles in the sample. RFQ flow overall
collapsed from ~200–280k/week (June, multi-sport) to ~376 last week.

RFQ resting behavior (research DB, n=8,140): age at our first action — median
15.2s, p75 23s, p90 28.5s; 2.7% first seen at ≥60s. 979 `quote_priced` events
at those ages. Conclusion: a 10–20s async fetch delay forfeits only the tail.

### 1.3 Feasibility (per-book wire audit, 2026-07-08)

The wire layer mostly exists; the resolution layer does not.

| Book | Arbitrary-N wire call | Leg→selection lookup exists? | On-demand caveats |
|---|---|---|---|
| DraftKings | `post_calculate_bets(session, sel_ids)` (`scraper_draftkings_trifecta.py:104`) | yes — sel-id dicts for spread/total/ML per event (`scraper_draftkings_sgp.py:420,576`) | server rejects some pairings (`combinabilityRestrictions`, 422); 2MB structure payload (TTL-cached) |
| FanDuel | `price_combo` implyBets — body is a list; needs 2→N loop (`scraper_fanduel_sgp.py:452`) | yes — runners cache spread/total/ML (`fanduel.py:409`) | PerimeterX token is long-lived but manual-refresh |
| ProphetX | `submit_parlay_rfq(legs: list)` (`prophetx_client.py:138`) | yes — market tree + ML selections (`prophetx.py:478,695`) | **each query = a real RFQ on their exchange** (footprint); teaser-tier filter |
| Novig | `submit_parlay(outcome_ids: list)` (`novig_client.py:157`) | yes — GraphQL market tree (`scraper_novig_sgp.py:440`) | mirrors DK prices (weak independence) |
| BetMGM | `price_picks(fixture, legs: list)` (`betmgm_client.py:172`) | yes — `parse_markets` incl. ML (`betmgm.py:95`) | declines un-comboable picks silently (no pricing group) |
| Caesars | `price_combo(legs: list)` (`caesars_client.py:228`) | yes — exact-name classify (`caesars.py:114`) | **AWS-WAF token (~1.5s mint, ~4min TTL), 3-wide cap** — shared with grid sweep |

What's missing (all books): the lookups live in per-book closures inside each
orchestrator's `price_sgps`, keyed differently per book; there is no
`resolve(CanonicalLeg) → selection` function, no single-combo entry point on
`SGPService` (its only public method sweeps all targets), and no
representation for exotic-shape results (`PricedRow.combo` is a fixed name
string; `mlb_sgp_odds` only has `spread_line`/`total_line` columns).

## 2. Goal and non-goals

**Goal.** Any RFQ where every game's sub-combo is either a Phase 1 route or an
on-demand-priceable shape (≤ `MAX_ON_DEMAND_LEGS` legs that ≥2 books will
price as a full partition) gets a two-sided quote. Cross-game multiplication,
gates, sizing, and lifecycle are unchanged.

**Non-goals.**
- No new market types: legs are still MLB FG spread / total / moneyline
  (`KXMLBSPREAD-`/`KXMLBTOTAL-`/`KXMLBGAME-`). Player props, F5, etc. stay
  unpriceable.
- No persistence of on-demand prices into `mlb_sgp_odds` (in-memory cache +
  research events only; schema untouched).
- No change to taker (`kalshi_mlb_rfq`) — maker only.
- Lone singles stay excluded (unchanged Phase 1 scope decision).

## 3. Architecture

```
                       main.py discovery tick (2s, never blocks)
                       ────────────────────────────────────────
 RFQ → parse → partition by game → per game:
   grid shape ──────────────► price from _SGP_ODDS (Phase 1, unchanged)
   on_demand ─┬─ cache fresh ► price from OnDemandCache ─► consensus ─► quote
              └─ cache miss ─► enqueue(leg_set_hash) ─► skip: "on_demand_pending"

 kalshi_mlb_mm/on_demand.py (new)                 kalshi_common / mlb_sgp
 ───────────────────────────────                  ──────────────────────
 FetchQueue (dedup by hash, cap)                  SGPService.price_partition(
   │  background worker thread                        book, game, legs)
   ▼                                                  │ per-book thread, reuses
 for each book in ON_DEMAND_BOOKS ──────────────────► │ persistent client +
   fetch 2^N cells → sanity gate → probit devig       │ structure TTL caches
   ▼                                                  ▼
 OnDemandCache[hash] = {book: fair, fetched_at}   resolve_legs(event, legs)
   ▲                                                  → book selection IDs
   │ re-fetch every SGP_REFRESH_SEC for combos
   │ with an OPEN quote (keeps last-look alive)
 router.subcombo_fair(..., on_demand_fairs)  ← pure lookup, injected by main
```

Data flow for one novel RFQ (worked example — 3-leg: Home -1.5, Over 8.5,
Home ML):

1. Tick T: classify → `on_demand`, no cache → enqueue, skip
   (`on_demand_pending`).
2. Worker (runs ~T+0 to T+15s): for each book, resolve the 3 legs to
   selection IDs from the (cached) event structure, then price all 8 cells:

   | cell | legs queried |
   |---|---|
   | target | Home -1.5, Over 8.5, Home ML |
   | 2 | Home -1.5, Over 8.5, Away ML |
   | 3 | Home -1.5, Under 8.5, Home ML |
   | … | … all 8 sign combinations |

   Probit-devig the 8 implied probabilities → target-cell fair per book.
3. Tick T+16s (RFQ still open — median rest ≫ this): cache hit →
   `consensus()` over per-book fairs → gates → quote as normal.
4. Open-quote maintenance: the combo's hash joins the refresh set; the worker
   re-prices it every `SGP_REFRESH_SEC` so the risk sweep, circuit breaker,
   and confirm-tick last look see fresh fairs. Quote closed → hash drops out.

## 4. Components

### 4.1 `kalshi_common/legset.py` — partition enumeration (pure)

`enumerate_partition(game_legs) -> list[list[CanonicalLeg]]`: the 2^N lists of
legs covering every side combination (each leg flipped through its two sides;
spread/ml flip home↔away, total flips over↔under). First element = the target
(input) combination. Pure, fully unit-testable. Line values never change —
only sides — so leg→selection resolution is shared across cells.

### 4.2 `kalshi_common/fair_value.py` — partition devig (pure)

`devig_partition(cell_decimals: list[float]) -> float | None`: probit devig
(`_probit_devig_n`, already n-ary) over the 2^N raw implied probabilities;
returns the first (target) cell's devigged probability. Sanity gate before
devig: `1.0 <= sum(1/d) <= PARTITION_OVERROUND_MAX` (default 1.35) — catches
books returning nonsense for degenerate cells (e.g. total+total shapes where
`Over 8.5 & Under 7.5` is impossible); any missing/insane cell → return None
(book dropped). This is the same "full partition or nothing" rule as Phase 1's
4-cell requirement, generalized.

### 4.3 `mlb_sgp/*` — per-book leg resolution + single-combo pricing

Per book, two additions to the existing module (no new scrapers):

- `resolve_legs(structure, legs: list[CanonicalLeg]) -> list[SelectionRef] | None`
  — maps canonical legs to that book's selection descriptors using the SAME
  structure objects the orchestrators already fetch and cache (DK sel-id
  dicts, FD runners, PX market tree, NV UUIDs, MGM `(mid, oid)`, CZR leg
  dicts). None if any leg is unresolvable at the book (line not offered).
  This is extraction/refactor of logic that exists in closures — TDD'd against
  recorded structure fixtures, then live-verified per book.
- `price_selection_set(client, refs) -> float | None` — one wire call, decimal
  odds or None. Mostly exists; FD needs the 2→N loop, DK uses the trifecta
  N-leg call.

Existing orchestrator paths (`price_sgps`) are untouched — the grid sweep and
its outputs stay byte-identical.

### 4.4 `kalshi_common/sgp_service.py` — single-combo entry point

`SGPService.price_partition(book, game_key, legs) -> dict[cell_idx, float] | None`
— runs on the calling (worker) thread, reuses the book's persistent client and
structure TTL caches (a cold game costs one structure fetch, then 2^N price
calls). Respects the same failure accounting as `refresh()` (a failing book
increments `failures`, 3 strikes → client reinit) so on-demand and grid-sweep
health share one recovery path. Concurrency guard: on-demand calls to a book
serialize with that book's grid-sweep runner via a per-book lock —
`curl_cffi` sessions are not thread-safe (documented in `sgp_service.py`), and
today's code never uses one session from two threads; this preserves that.

### 4.5 `kalshi_mlb_mm/on_demand.py` — the stateful engine (new)

- `OnDemandEngine(service, books, config)` owning:
  - **Queue**: deduped by `leg_set_hash`; bounded (`ON_DEMAND_QUEUE_MAX`,
    default 8; overflow → drop + research event, fail-safe).
  - **Worker**: one background thread consuming the queue; per combo, fans out
    to the 6 books concurrently (thread per book, same pattern and deadline
    discipline as `SGPService.refresh`); writes
    `cache[hash] = {book: fair, fetched_at, legs}`.
  - **Cache**: in-memory dict; entry quote-eligible while
    `age <= MAX_BOOK_STALENESS_SEC` (180s, same constant as grids). Negative
    results cache too (`no_books_priced` until +`ON_DEMAND_NEGATIVE_TTL_SEC`,
    default 300s) so a shape no book offers doesn't get re-fetched every tick.
  - **Refresh set**: hashes of combos with an open `live_quotes` row; worker
    re-prices them every `SGP_REFRESH_SEC` (60s). Everything else ages out.
  - **Rate budget**: global `ON_DEMAND_MAX_COMBOS_PER_MIN` (default 6) — with
    the 3-leg cap that bounds worst-case on-demand traffic at 48
    requests/book/min, comparable to one grid sweep.
  - `lookup(hash) -> {book: fair} | None` — the pure read handed to the
    router. Thread-safety: single lock around cache mutations; lookups copy.
- Crash-safety: worker exceptions are caught and logged; a dead worker is
  restarted lazily on next enqueue; cache loss on bot restart is acceptable
  (re-fetch on demand).

### 4.6 `kalshi_mlb_mm/router.py` — the seam (Phase 1 lock intact)

```python
def subcombo_fair(game_id, game_legs, sgp_df, min_books, band,
                  on_demand_fairs=None):
    ...
    elif route == "on_demand" and on_demand_fairs is not None:
        book_fairs = on_demand_fairs(legset.leg_set_hash(game_legs)) or {}
    else:   # "on_demand" without lookup, or "unpriceable"
        return None
    return consensus(book_fairs, min_books, band)
```

`combo_fair` threads the same optional argument through. Grid and single
routes untouched; with `on_demand_fairs=None` behavior is byte-identical, so
the `router 2-leg == median(_book_fairs) to 1e-9` regression test passes
unmodified. Consensus: same `consensus()`, same `MIN_AGREEING_BOOKS=2`, same
band (user decision — no independence adjustment).

### 4.7 `kalshi_mlb_mm/main.py` — integration + instrumentation

- `_priceable_in_phase1` → `_priceable(canon)`: `on_demand` routes become
  in-scope when `ON_DEMAND_ENABLED` and the game's leg count ≤
  `MAX_ON_DEMAND_LEGS`.
- Discovery tick: on-demand games with a fresh cache price through the normal
  path (all existing gates apply unchanged); cache miss → `engine.enqueue()` +
  skip reason `on_demand_pending` (a new reason value; the monitor reads
  reason vocabularies from data, so it appears automatically).
- Confirm tick / risk sweep: no structural change — they call
  `router.combo_fair(..., on_demand_fairs=engine.lookup)`; a stale cache
  yields None → existing `voided_no_fresh_books` / quote-pull behavior. The
  refresh set keeps open-quote combos fresh so voids stay rare.
- Fill research event gains `consensus_books` (which books agreed) so
  DK+NV-only fills are countable (accepted-risk measurement).
- **Instrumentation (measurement blind spot fix):** out-of-scope MLB RFQs now
  store `legs_json` and a shape class in `seen_rfqs`, with distinct
  `last_decision` values (`out_of_scope_non_mlb`, `out_of_scope_lone_single`,
  `out_of_scope_shape_cap`, …). New research events: `on_demand_requested`,
  `on_demand_result` (per-book cell coverage, fairs, duration, consensus
  outcome).

### 4.8 Config additions (`kalshi_mlb_mm/config.py`)

```
ON_DEMAND_ENABLED=true            # master kill switch
ON_DEMAND_BOOKS=<all 6>           # per-book opt-out (user chose all 6)
MAX_ON_DEMAND_LEGS=3              # 2^N cap: 3 legs = 8 cells/book
ON_DEMAND_QUEUE_MAX=8
ON_DEMAND_MAX_COMBOS_PER_MIN=6
ON_DEMAND_NEGATIVE_TTL_SEC=300
PARTITION_OVERROUND_MAX=1.35
```

Freshness deliberately reuses `MAX_BOOK_STALENESS_SEC` and `SGP_REFRESH_SEC` —
one staleness regime for the whole bot.

## 5. Failure modes (all fail toward "don't quote")

| Failure | Behavior |
|---|---|
| Book can't resolve a leg (line not offered) | book returns None → dropped from consensus |
| Book prices some cells but not all 2^N | partition incomplete → book dropped |
| Partition overround insane (degenerate shape) | sanity gate → book dropped |
| < MIN_AGREEING_BOOKS survive | consensus None → `no_fair` skip |
| Fetch in flight when RFQ re-seen | still `on_demand_pending`; no duplicate enqueue |
| No book prices the shape | negative cache 300s → no re-fetch churn |
| Cache stale at confirm-time last look | `voided_no_fresh_books` (existing rule) |
| Worker thread dies | lazy restart on next enqueue; meanwhile all on_demand → skip |
| Queue full / rate budget exhausted | drop with research event; RFQ simply not quoted |
| Caesars WAF trips | book 3-strike reinit (shared with grid sweep — accepted risk) |
| Bot restart | cache empty → shapes re-fetch on demand |

## 6. Testing (TDD throughout)

1. **Pure math first**: `enumerate_partition` (side flipping, ordering,
   hash stability), `devig_partition` (4-cell case must equal existing
   `devig_book` output on the same grid — a second regression bridge to
   Phase 1; sanity-gate edges; degenerate-cell None).
2. **Resolution layer**: per-book `resolve_legs` against recorded structure
   fixtures (the repos' existing fixture pattern); unresolvable-leg → None.
3. **Engine**: fake `SGPService` (constructor-injected, mirrors the existing
   `runners` test seam) — queue dedup, TTL expiry, negative cache, refresh
   set, rate budget, worker-death restart. No network in tests.
4. **Router**: injected `on_demand_fairs` — consensus paths, None paths, and
   the untouched-grid regression suite (must pass unmodified, 1e-9 lock).
5. **Integration**: discovery→pending→cached→quoted flow with everything
   faked; confirm-tick staleness void.
6. **Live verification (manual, before merge)**: one real exotic combo priced
   per book from a REPL, logged; then `--dry-run` soak with synthetic
   self-created RFQs (we can create RFQs on Kalshi as the taker does) to
   exercise the full path end-to-end.

## 7. Rollout & measurement

1. Merge → restart maker **from main repo cwd** ([[kalshi_mlb_rfq_restart_gotchas]]
   applies to the maker too) in `--dry-run`; confirm `on_demand_*` research
   events flow on synthetic RFQs.
2. Live-enable. Watch (monitor + research DB): on-demand RFQ count/shapes
   (finally measurable), per-book partition-completion rate, fetch latency,
   `dk_nv_only` consensus share, Caesars WAF failure rate vs grid-sweep
   health, void rate on on-demand quotes.
3. Pre-registered retreat triggers: Caesars grid-sweep failures attributable
   to on-demand → drop CZR from `ON_DEMAND_BOOKS`; on-demand void rate
   materially above grid void rate → raise freshness strictness or disable.

## 8. Version control & worktree

- Worktree `worktree-kalshi-mm-on-demand-phase2` (this spec is its first
  commit); branch of the same name off local main.
- Commit sequence (one reviewable unit each): spec → pure math
  (legset/fair_value) → per-book resolvers (one commit per book, live-verified)
  → SGPService entry point → engine → router seam + main integration +
  instrumentation → docs.
- Full test suite (1,042 + new) green before merge; adversarial review before
  merge; **no merge or push without explicit user approval**; worktree + branch
  cleaned up immediately after merge.

## 9. Documentation (same merge as code)

- `kalshi_mlb_mm/README.md`: on-demand architecture section, config table,
  failure modes, rollout/retreat playbook.
- `NFLWork/CLAUDE.md`: maker bullet — Phase 2 on-demand pricing summary.
- `kalshi_mlb_monitor/README.md`: only if a monitor tweak is needed (reason
  vocabulary is data-driven, so likely none).
- Memory file `kalshi_mm_generalized_combos.md`: Phase 2 status update.
