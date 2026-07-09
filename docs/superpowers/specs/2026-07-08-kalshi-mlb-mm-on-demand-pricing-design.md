# Kalshi MLB MM — Phase 2: On-Demand Pricing of Arbitrary Leg Combos

**Date:** 2026-07-08 (rev 2, 2026-07-09 — live-per-RFQ pricing, one-sided-leg
route, limitation configs removed per user review)
**Status:** Draft — awaiting user review
**Depends on:** Phase 1 generalized combos (merged to main 2026-07-08, `2e9a0e8`);
spec at `docs/superpowers/specs/2026-06-28-kalshi-mlb-mm-generalized-combos-design.md`

---

## Review Pack

**What we're building** — The maker currently skips any same-game leg combo
that isn't one of the two pre-scraped 2-leg grids (spread×total, ml×total).
Phase 2 makes the quoting engine fully general: when an RFQ contains a novel
same-game shape (any leg count, any market mix), the bot live-queries the six
sportsbooks' SGP endpoints for that exact leg set, devigs, builds consensus,
and quotes — asynchronously, so the 2-second RFQ discovery loop never blocks
on a book, and **every quote is backed by a fetch triggered by that specific
RFQ** (no price reuse, no staleness window). Goal per user direction: "quote
any combination of legs given we have access to it through our sportsbook
scraping," with no artificial leg-count or rate limitations.

**Key decisions**

1. **Async fetch, quote on landing — never block the 2s discovery loop.**
   A novel shape enqueues a background fetch (~10–20s); the RFQ is quoted
   within one tick of the price landing. Rejected: blocking the discovery tick
   on live HTTP. Evidence async works: the bot already first acts on RFQs at
   median age 15s (p90 29s) and those quotes succeed — RFQs rest minutes.
2. **Live-per-RFQ pricing, no reuse** (user requirement, rev 2). A fetched
   price may back a NEW quote only within ~15s of landing (one tick + jitter);
   a later RFQ on the same combo triggers a fresh fetch. Resting quotes are
   re-fetched every 60s and pulled immediately if a refresh fails; the
   confirm-time last look is a synchronous live re-fetch (can't complete →
   void). Rejected: the rev-1 TTL cache with 180s quote-eligibility — user:
   "we can't quote stale prices."
3. **Two devig routes, no leg cap** (rev 2). Route A (gold standard, ≤3
   legs, all cells available): fetch the full 2^N side-combination partition,
   probit-devig, read the target cell — the exact generalization of Phase 1's
   4-cell rule. Route B (any N, and any leg one-sided at the book):
   book-implied correlation transfer — one SGP call/book; fair =
   ∏(devigged singles) × [SGP implied / ∏(vigged singles)], singles devigged
   from the structure fetch we already cache. Rejected: partition-only with a
   hard leg cap (rev 1) — it made one-sided lines and 4+ legs refusal paths.
4. **All 6 books participate** (user decision, rev 1). Caesars (WAF fragility)
   and ProphetX (an on-demand query is a real RFQ on their exchange) were
   offered as exclusions and declined; both get failure isolation +
   measurement, and `ON_DEMAND_BOOKS` remains the retreat lever.
5. **Plain 2-book consensus, same as grids** (user decision, rev 1). A
   DK+Novig pair counts even though Novig mirrors DK; fills carry a
   `consensus_books` research field so the risk is measured, not silent.
6. **Router stays pure; grid pricing byte-identical.** `subcombo_fair` gains
   one optional `on_demand_fairs=None` lookup argument. Default None
   reproduces Phase 1 exactly; the 1e-9 regression lock is untouched. All
   stateful machinery lives in a new `kalshi_mlb_mm/on_demand.py`.

**Risks / push back here**

- **Target flow is currently zero.** Live measurement (June 9 → July 8): 0
  novel same-game shapes in ~1,100 classified RFQ markets (exhaustive last-7-
  days + 1,500-market random sample; 95% upper bound ≈ 0.3% of out-of-scope
  flow). You chose to build the general engine anyway; instrumentation ships
  in the same merge so usage is visible from day one.
- **Unbounded quoting needs bounded fetching.** With rate-limit configs
  removed, the only pacing is per-book serialization (one request stream per
  book — forced anyway by `curl_cffi` thread-unsafety). In a June-style RFQ
  flood (~35k/day) the fetch queue would grow and prices would land after
  RFQs expire — quotes forgone, never stale, and the books see at most one
  serialized stream per book. That's the graceful-degradation story; if you'd
  rather shed load explicitly under flood, say so and a backstop returns.
- **Caesars blast radius.** On-demand shares Caesars' WAF-token session with
  the 60s grid sweep; tripping the WAF costs the grid pipeline Caesars too.
  Mitigation: failure isolation + `ON_DEMAND_BOOKS` retreat lever. Residual
  risk accepted by user (rev 1).
- **Route B rests on a vig-cancellation assumption** (the book's SGP margin ≈
  compounded leg vig, so it divides out of the correlation ratio). Extra SGP
  margin beyond that inflates Route B fairs slightly. Defense: consensus band
  across books + a partition-vs-ratio cross-check research metric on every
  ≤3-leg combo (both routes computable there), so the assumption is measured
  from day one.
- **DK+Novig-only consensus on correlation-heavy shapes** is the known
  "fill our fair-error tail" vector from the adversarial review. Accepted;
  measured via `consensus_books`.
- **In-scope quoting bottleneck is elsewhere.** All-time skip reasons:
  `no_fair` 5,431 vs 518 quoted. If the goal is more fills, book
  freshness/consensus tuning on already-supported shapes is the larger lever.

**Worth understanding** (opt-in)

1. *Joint-distribution partitions.* A 2-leg combo's 4-cell grid is a 2×2
   contingency table whose cells sum to 1; N legs → `expand.grid` of sides =
   2^N cells. Devigging over the full partition removes the book's margin
   without assuming where it hides — that's why it's the preferred route.
2. *Ratio estimators.* Route B is a classic ratio trick: a bias present in
   both numerator (vigged SGP) and denominator (product of vigged singles)
   approximately cancels, leaving the correlation multiplier. Like comparing
   two measurements from the same miscalibrated scale.
3. *Producer/consumer worker.* The discovery tick *produces* fetch jobs onto a
   queue and moves on; a background thread *consumes* them and posts results.
   The two sides communicate only through the queue and the results store —
   no shared mutable pricing state, which is what keeps the 2s tick safe.

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
before its deploy). RFQ flow overall collapsed from ~200–280k/week (June,
multi-sport) to ~376 last week.

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
on-demand shape that ≥2 books will price gets a two-sided quote backed by
prices fetched for that RFQ. No leg-count cap, no per-shape refusal rules —
the only refusals are fail-safe ones (books won't price it / consensus fails /
price can't be kept live). Cross-game multiplication, gates, sizing, and
lifecycle are unchanged.

**Non-goals.**
- No new market types: legs are still MLB FG spread / total / moneyline
  (`KXMLBSPREAD-`/`KXMLBTOTAL-`/`KXMLBGAME-`). Player props, F5, etc. stay
  unpriceable (the books' structure caches don't carry them; separate project).
- No persistence of on-demand prices into `mlb_sgp_odds` (results live in the
  engine + research events; schema untouched).
- No change to taker (`kalshi_mlb_rfq`) — maker only.
- Lone singles stay excluded (unchanged Phase 1 scope decision).

## 3. Architecture

```
                      main.py discovery tick (2s, never blocks)
                      ────────────────────────────────────────
 RFQ → parse → partition by game → per game:
   grid shape ─────────────► price from _SGP_ODDS (Phase 1, unchanged)
   on_demand ─┬ live result (≤15s old, THIS combo) ► consensus ► quote
              └ none/aged ► ensure fetch queued ► skip: "on_demand_pending"
                              (re-fetch per RFQ — results are never reused)

 kalshi_mlb_mm/on_demand.py (new)                kalshi_common / mlb_sgp
 ───────────────────────────────                 ──────────────────────
 FetchQueue (dedup in-flight, FIFO,              SGPService.price_on_demand(
   unbounded — pacing = per-book                     book, game, legs)
   serialization, never refusal)                     │ reuses persistent client
   │ background worker                               │ + structure TTL caches
   ▼                                                 ▼
 per book (concurrent, serialized           resolve_legs(structure, legs)
 against that book's grid sweep):                → book selection IDs
   Route A (≤3 legs, all cells offered):
     fetch 2^N cells → sanity gate → probit devig over partition
   Route B (any N, or any one-sided leg):
     fetch 1 SGP price → × correlation transfer vs structure singles
   ▼
 ResultStore[hash] = {book: fair, route, landed_at}   (fresh ≤15s for NEW
   ▲                                                   quotes; then dead)
   │ open-quote maintenance: re-fetch every 60s; any failed
   │ refresh → PULL the quote immediately
 router.subcombo_fair(..., on_demand_fairs)  ← pure lookup, injected by main
```

**Worked example** — 3-leg RFQ: Home -1.5, Over 8.5, Home ML:

1. Tick T: classify → `on_demand`, no live result → enqueue, skip
   (`on_demand_pending`).
2. Worker (~T+0 to T+15s), per book: resolve the 3 legs to selection IDs from
   the (cached) event structure. All 6 sides offered → Route A — price all 8
   cells:

   | cell | legs queried |
   |---|---|
   | target | Home -1.5, Over 8.5, Home ML |
   | 2 | Home -1.5, Over 8.5, Away ML |
   | 3 | Home -1.5, Under 8.5, Home ML |
   | … | … all 8 sign combinations |

   → probit devig over the 8 → target-cell fair.
   **One-sided variant:** the book has pulled Away +1.5 → cells needing it are
   unpriceable → Route B instead: one SGP call for the target combo, singles
   devigged from the structure (e.g. SGP implied 16.7%; vigged singles 52% ×
   51% × 60%; multiplier 16.7%/15.9% = 1.05; fair = 49%×48%×57% × 1.05 =
   14.1%). If the *single* is one-sided too, its devig uses the book's
   existing vig-fallback haircut (`DK_VIG_FALLBACK` etc.); still
   consensus-gated downstream.
3. Tick ≈ landing+≤2s (RFQ still open — median rest ≫ this): live result →
   `consensus()` over per-book fairs → all existing gates → quote.
4. While the quote rests: the combo joins the maintenance set — full re-fetch
   every `SGP_REFRESH_SEC` (60s); reprice through the normal hysteresis /
   circuit-breaker path; **a failed refresh pulls the quote immediately** (no
   stale resting exposure beyond one refresh interval).
5. On acceptance: confirm-tick last look triggers a **synchronous live
   re-fetch** of this combo (worker priority lane); if it can't complete and
   pass the drift check inside the confirm window → void. We never confirm on
   a previously-fetched number.

## 4. Components

### 4.1 `kalshi_common/legset.py` — partition enumeration (pure)

`enumerate_partition(game_legs) -> list[list[CanonicalLeg]]`: the 2^N lists of
legs covering every side combination (spread/ml flip home↔away, total flips
over↔under); first element = the target combination. Used by Route A (N ≤ 3 —
beyond that Route B is used, so 2^N never exceeds 8 requests). Pure, fully
unit-testable. Lines never change — only sides — so leg→selection resolution
is shared across cells.

### 4.2 `kalshi_common/fair_value.py` — two devig routes (pure)

- `devig_partition(cell_decimals: list[float]) -> float | None` — Route A:
  probit devig (`_probit_devig_n`, already n-ary) over the 2^N raw implied
  probabilities; returns the target cell's devigged probability. Sanity gate:
  `1.0 <= sum(1/d) <= PARTITION_OVERROUND_MAX` (constant 1.35) — catches books
  returning nonsense for degenerate cells (e.g. total+total shapes where
  `Over 8.5 & Under 7.5` is impossible). Any missing/insane cell → None; the
  caller then tries Route B.
- `fair_by_correlation_transfer(sgp_decimal, singles) -> float | None` — Route
  B: `singles` is per-leg `(vigged_implied_of_chosen_side, devigged_fair)`
  from the book's own structure payload (all six structure fetches carry
  two-sided single odds for spread/total/ML). Returns
  `∏ fair_i × (1/sgp_decimal) / ∏ vigged_i`. Assumption (measured, see §7):
  the book's SGP margin ≈ compounded leg vig, so it cancels in the ratio.
  Sanity gate: resulting multiplier within [0.5, 2.0], else None.
- Route selection per (book, combo): A if N ≤ 3 and all 2^N cells priced,
  else B, else book drops out. On every combo where BOTH routes are
  computable, both are computed and their gap emitted to research
  (`route_gap`) — a free live experiment on the Route B assumption.

### 4.3 `mlb_sgp/*` — per-book leg resolution + single-combo pricing

Per book, two additions to the existing module (no new scrapers):

- `resolve_legs(structure, legs: list[CanonicalLeg]) -> list[SelectionRef] | None`
  — maps canonical legs to that book's selection descriptors using the SAME
  structure objects the orchestrators already fetch and cache (DK sel-id
  dicts, FD runners, PX market tree, NV UUIDs, MGM `(mid, oid)`, CZR leg
  dicts), and also returns each leg's two-sided single odds for Route B. None
  for any unresolvable leg (that side not offered → drives the one-sided
  routing above). Extraction/refactor of logic that exists in closures —
  TDD'd against recorded structure fixtures, then live-verified per book.
- `price_selection_set(client, refs) -> float | None` — one wire call, decimal
  odds or None. Mostly exists; FD needs the 2→N loop, DK uses the trifecta
  N-leg call.

Existing orchestrator paths (`price_sgps`) are untouched — the grid sweep and
its outputs stay byte-identical.

### 4.4 `kalshi_common/sgp_service.py` — single-combo entry point

`SGPService.price_on_demand(book, game_key, legs) -> OnDemandBookResult | None`
(cells or single SGP price + singles, per the route) — runs on the calling
(worker) thread, reuses the book's persistent client and structure TTL caches
(a cold game costs one structure fetch, then the price calls). Shares
`refresh()`'s failure accounting (3 strikes → client reinit) so on-demand and
grid-sweep health share one recovery path. Concurrency: a per-book lock
serializes on-demand calls with that book's grid-sweep runner — `curl_cffi`
sessions are not thread-safe (documented in `sgp_service.py`); this per-book
serialization is also the system's natural fetch pacing.

### 4.5 `kalshi_mlb_mm/on_demand.py` — the stateful engine (new)

- `OnDemandEngine(service, books)` owning:
  - **Queue**: FIFO, unbounded, deduped on in-flight `leg_set_hash` (two RFQs
    on the same combo in the same window share one fetch — that's
    concurrency-dedup, not price reuse). No drop path: under load, fetches
    land later and RFQs that expired meanwhile simply go unquoted — degraded
    means *fewer quotes*, never *staler quotes*.
  - **Worker**: one background thread consuming the queue; per combo, fans out
    to `ON_DEMAND_BOOKS` concurrently (thread per book, same deadline
    discipline as `SGPService.refresh`). A **priority lane** serves confirm-
    time last-look re-fetches ahead of discovery fetches.
  - **ResultStore**: `hash → {book: fair, route, landed_at}`. A result may
    back a NEW quote only while `age <= QUOTE_FRESH_SEC` (constant 15s — one
    discovery tick + landing jitter; the price is seconds old when quoted).
    After that it is dead for quoting; it is kept briefly only for research
    comparison. There is deliberately NO quote-eligibility TTL to tune.
  - **Maintenance set**: hashes of combos with an open `live_quotes` row —
    full re-fetch every `SGP_REFRESH_SEC` (60s). A refresh that fails (books
    < consensus minimum) → the engine flags the hash and the risk sweep pulls
    the quote on its next pass (≤10s later). Quote closed → hash drops out.
  - `lookup(hash) -> {book: fair} | None` — returns per-book fairs only while
    fresh (for new quotes) or maintained (for open-quote re-pricing); the pure
    read handed to the router. Single lock around store mutations; lookups
    copy.
- Crash-safety: worker exceptions caught and logged; dead worker restarts
  lazily on next enqueue; restart loses nothing durable (results are
  per-RFQ ephemeral by design).

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

- `_priceable_in_phase1` → `_priceable(canon)`: `on_demand` routes are
  in-scope whenever `ON_DEMAND_ENABLED` (no leg-count condition).
- Discovery tick: on-demand games with a live (≤15s) result price through the
  normal path (all existing gates unchanged); otherwise
  `engine.ensure_fetch()` + skip reason `on_demand_pending` (new reason value;
  the monitor reads reason vocabularies from data, so it appears
  automatically).
- Confirm tick: for combos containing an on-demand game, the last look calls
  `engine.refetch_now(hash, deadline)` (priority lane, synchronous); on
  timeout or consensus failure → existing `voided_no_fresh_books` path. Grid
  combos keep today's behavior.
- Risk sweep: unchanged code path via `combo_fair(...,
  on_demand_fairs=engine.lookup)`; the maintenance set guarantees those
  lookups are ≤60s old or the quote is being pulled.
- Fill research event gains `consensus_books` and `route` per book.
- **Instrumentation (measurement blind-spot fix):** out-of-scope MLB RFQs now
  store `legs_json` and a shape class in `seen_rfqs`, with distinct
  `last_decision` values (`out_of_scope_non_mlb`, `out_of_scope_lone_single`,
  …). New research events: `on_demand_requested`, `on_demand_result`
  (per-book route, cell coverage, fairs, latency, consensus outcome,
  `route_gap` when both routes computed).

### 4.8 Config additions (`kalshi_mlb_mm/config.py`)

```
ON_DEMAND_ENABLED=true      # master kill switch (operational, not a limit)
ON_DEMAND_BOOKS=<all 6>     # per-book retreat lever (user chose all 6)
```

Deliberately nothing else (user decision, rev 2): no leg cap, no queue cap,
no rate budget, no tunable TTLs. `QUOTE_FRESH_SEC=15` and
`PARTITION_OVERROUND_MAX=1.35` are module constants (correctness guards, not
knobs); freshness cadence reuses `SGP_REFRESH_SEC`; pacing is per-book
serialization.

## 5. Failure modes (all fail toward "don't quote" / "fewer quotes", never "staler quotes")

| Failure | Behavior |
|---|---|
| Leg's side not offered at a book (one-sided line) | Route A impossible at that book → Route B (single SGP call + correlation transfer) |
| Single itself one-sided at the book | Route B single devig uses the book's vig-fallback haircut; consensus still gates |
| Book prices some cells but not all 2^N | Route A → Route B fallback; if the SGP call also fails, book drops out |
| Partition overround / multiplier insane | sanity gate → book drops out |
| Book rejects the combo (DK combinability, MGM no-group) | book drops out |
| < MIN_AGREEING_BOOKS survive | consensus None → `no_fair` skip |
| Fetch in flight when RFQ re-seen | shared fetch (dedup); still `on_demand_pending` |
| Result older than 15s when RFQ (re-)seen | dead for quoting → fresh fetch enqueued |
| RFQ flood (June-style) | queue grows; prices land late; expired RFQs go unquoted; books see one serialized stream each |
| Open-quote refresh fails | quote pulled within one risk-sweep pass (≤10s) |
| Confirm-time re-fetch can't complete in window | `voided_no_fresh_books` — never confirm on an old number |
| Worker thread dies | lazy restart on next enqueue; meanwhile all on_demand → skip |
| Caesars WAF trips | book 3-strike reinit (shared with grid sweep — accepted risk; retreat lever = `ON_DEMAND_BOOKS`) |
| Bot restart | engine state empty → everything re-fetches per RFQ (by design) |

## 6. Testing (TDD throughout)

1. **Pure math first**: `enumerate_partition` (side flipping, ordering, hash
   stability); `devig_partition` (4-cell case must equal existing `devig_book`
   on the same grid — a second regression bridge to Phase 1; sanity-gate
   edges); `fair_by_correlation_transfer` (vig-cancellation on synthetic
   books, multiplier sanity bounds, one-sided-single haircut path); route
   selection (A↔B fallbacks).
2. **Resolution layer**: per-book `resolve_legs` against recorded structure
   fixtures (repo's existing fixture pattern); unresolvable-side → drives
   Route B; singles extraction both sides.
3. **Engine**: fake `SGPService` (constructor-injected, mirrors the existing
   `runners` test seam) — in-flight dedup, 15s freshness death, maintenance
   refresh + pull-on-fail flag, priority lane, worker-death restart, flood
   behavior (queue grows, no drops, no stale results served). No network.
4. **Router**: injected `on_demand_fairs` — consensus paths, None paths, and
   the untouched-grid regression suite (must pass unmodified, 1e-9 lock).
5. **Integration**: discovery→pending→landed→quoted flow with everything
   faked; confirm-tick synchronous re-fetch success/timeout/void; risk-sweep
   pull on failed maintenance.
6. **Live verification (manual, before merge)**: one real exotic combo priced
   per book per route from a REPL, logged; then `--dry-run` soak with
   synthetic self-created RFQs (as the taker does) to exercise the full path
   end-to-end, including a deliberately one-sided leg.

## 7. Rollout & measurement

1. Merge → restart maker **from main repo cwd** (restart-gotchas rule applies
   to the maker) in `--dry-run`; confirm `on_demand_*` research events flow on
   synthetic RFQs.
2. Live-enable. Watch (monitor + research DB): on-demand RFQ count/shapes
   (finally measurable), per-book route mix (A vs B) and partition-completion
   rate, `route_gap` distribution (validates the Route B vig-cancellation
   assumption), fetch latency, `consensus_books` composition (DK+NV-only
   share), Caesars WAF failure rate vs grid-sweep health, void rate and
   quote-pull rate on on-demand quotes.
3. Pre-registered retreat triggers: Caesars grid-sweep failures attributable
   to on-demand → drop CZR from `ON_DEMAND_BOOKS`; `route_gap` shows Route B
   systematically rich → haircut or disable Route B (that's a code decision,
   surfaced with data, not a config).

## 8. Version control & worktree

- Worktree `worktree-kalshi-mm-on-demand-phase2` (this spec is its first
  commit); branch of the same name off local main.
- Commit sequence (one reviewable unit each): spec → pure math
  (legset/fair_value, both routes) → per-book resolvers (one commit per book,
  live-verified) → SGPService entry point → engine → router seam + main
  integration + instrumentation → docs.
- Full test suite (1,042 + new) green before merge; adversarial review before
  merge; **no merge or push without explicit user approval**; worktree +
  branch cleaned up immediately after merge.

## 9. Documentation (same merge as code)

- `kalshi_mlb_mm/README.md`: on-demand architecture, the two devig routes,
  liveness rules, failure modes, rollout/retreat playbook.
- `NFLWork/CLAUDE.md`: maker bullet — Phase 2 on-demand pricing summary.
- `kalshi_mlb_monitor/README.md`: only if a monitor tweak is needed (reason
  vocabulary is data-driven, so likely none).
- Memory file: Phase 2 status update.
