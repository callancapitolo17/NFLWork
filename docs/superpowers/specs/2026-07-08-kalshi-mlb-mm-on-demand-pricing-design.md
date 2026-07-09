# Kalshi MLB MM — Phase 2: On-Demand Pricing of Arbitrary Leg Combos

**Date:** 2026-07-08 (rev 4, 2026-07-09 — RFQ-poll-driven live feed replaces
both rev 2's maintenance loop and rev 3's rest-untouched model; rev 3 removed
on-demand configs; rev 2 added live-per-RFQ pricing + one-sided-leg route)
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
on a book. Liveness model: **while an RFQ is open it is on a live feed** —
re-fetched every ~15s and re-quoted through the existing hysteresis/replace
logic, exactly how grid quotes already track the 60s sweep at 2s cadence —
and a fill is confirmed only after a synchronous live re-fetch. Goal per user
direction: "quote any combination of legs given we have access to it through
our sportsbook scraping," with no artificial limitations or switches.

**Key decisions**

1. **Async fetch, quote on landing — never block the 2s discovery loop.**
   A novel shape enqueues a background fetch (~10–20s); the RFQ is quoted
   within one tick of the price landing. Rejected: blocking the discovery tick
   on live HTTP. Evidence async works: the bot already first acts on RFQs at
   median age 15s (p90 29s) and those quotes succeed.
2. **The RFQ poll drives a continuous live feed** (rev 4). Kalshi RFQs are
   competitive auctions: quotes exist only while their RFQ is open (a quote
   dies with its RFQ), and the discovery tick already re-prices every open
   RFQ every 2s, replacing our quote when fair moves beyond hysteresis
   (`main.py:632-649`). On-demand combos plug into that loop: while the RFQ
   is in the poll, a re-fetch fires whenever the last result is ≥15s old;
   each landing re-quotes via the existing replace logic; when the RFQ leaves
   the poll, fetching stops naturally. No separate maintenance subsystem
   (rejected, rev 2 — redundant with the poll) and no rest-untouched model
   (rejected, rev 3 — incoherent with the discovery loop and left drift
   uncovered for the RFQ's lifetime). Residual drift exposure = one feed
   interval (~15–25s), backstopped by the confirm-time live re-fetch.
3. **Two devig routes, no leg cap** (rev 2). Route A (gold standard, ≤3 legs,
   all cells offered): fetch the full 2^N side-combination partition,
   probit-devig, read the target cell — the exact generalization of Phase 1's
   4-cell rule. Route B (any N, and any leg one-sided at the book):
   book-implied correlation transfer — one SGP call/book; fair =
   ∏(devigged singles) × [SGP implied / ∏(vigged singles)]. A lone SGP price
   cannot be devigged by itself (fair and margin are entangled in one
   number); the partition and the book's own singles are the only two margin
   references available — hence exactly these two routes. The singles cost
   nothing: they ride along in the structure fetch already required to
   resolve legs → selection IDs (TTL-cached, usually warm from the grid
   sweep). Rejected: partition-only with a leg cap (rev 1) — it made
   one-sided lines and 4+ legs refusal paths.
4. **All 6 books, always on, zero new config** (user decisions, revs 1+3).
   No switches: the on-demand path is a permanent part of the pricing engine,
   books = the SGP service's book set. The existing bot-wide kill file covers
   emergencies; dropping a misbehaving book is a one-line code change.
   `QUOTE_FRESH_SEC=15` and `PARTITION_OVERROUND_MAX=1.35` are module
   constants (correctness guards, not knobs).
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
- **Feed cadence is ~15–25s, not 2s.** Between fetch landings an on-demand
  quote can drift off a fast-moving market (grids have the same gap vs their
  60s sweep, defended identically: margin + hysteresis + circuit-breaker +
  confirm last look, which here is a live re-fetch). If books reprice faster
  than the feed on news (lineup scratch), the last look is the backstop and a
  void is the cost.
- **No rate limiting anywhere.** Pacing is per-book serialization only (one
  request stream per book — forced anyway by `curl_cffi` thread-unsafety).
  Steady-state cost: each open on-demand RFQ ≈ 2–4 fetch rounds/min for its
  open life. In a June-style RFQ flood (~35k/day) the fetch queue grows and
  prices land late — quotes forgone or lagged, never artificially refused. If
  you'd rather shed load explicitly under flood, say so and a backstop
  returns.
- **Caesars blast radius.** On-demand shares Caesars' WAF-token session with
  the 60s grid sweep; tripping the WAF costs the grid pipeline Caesars too.
  With no config lever, retreat is a one-line code change shipped on
  evidence. Residual risk accepted by user.
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

1. *One price can't be devigged alone.* A quoted price is fair probability ×
   margin, entangled in one number. Separating them needs a reference that
   reveals the margin's size: either a full partition (cells must sum to 1,
   the excess IS the margin — Route A) or the book's own two-sided singles
   (per-leg margin computed exactly, assumed to compound into the SGP's —
   Route B). Like correcting an instrument's bias by measuring a known
   standard on the same instrument.
2. *RFQs are auctions, not order books.* A quote is an answer to one
   creator's open request, competing with other makers' answers; it exists
   only while the request does. The maker's job while the auction runs is to
   keep its answer current (re-quote on movement); its protection at the
   moment of truth is the last look. Both already exist in this bot — Phase 2
   just puts on-demand combos on the same treadmill.
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
on-demand shape that ≥2 books will price gets a two-sided quote that stays
current for as long as the RFQ is open. No leg-count cap, no per-shape refusal
rules, no switches — the only refusals are fail-safe ones (books won't price
it / consensus fails / price isn't fresh at action time). Cross-game
multiplication, gates, sizing, and lifecycle are unchanged.

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
 every open RFQ, every tick → parse → partition by game → per game:
   grid shape ─────────────► price from _SGP_ODDS (Phase 1, unchanged)
   on_demand ─┬ live result (≤15s) ► consensus ► quote / hysteresis-replace
              └ aged or absent ───► ensure fetch queued ► skip this tick
                                     (the open-RFQ poll IS the feed driver:
                                      re-fetch fires each time the result
                                      ages past 15s, for as long as the RFQ
                                      stays in the poll; RFQ closes → feed
                                      stops naturally)

 kalshi_mlb_mm/on_demand.py (new)                kalshi_common / mlb_sgp
 ───────────────────────────────                 ──────────────────────
 FetchQueue (dedup in-flight, FIFO,              SGPService.price_on_demand(
   unbounded — pacing = per-book                     book, game, legs)
   serialization, never refusal)                     │ reuses persistent client
   │ background worker                               │ + structure TTL caches
   ▼                                                 ▼
 per book (concurrent, serialized           resolve_legs(structure, legs)
 against that book's grid sweep):                → book selection IDs
   Route A (≤3 legs, all cells offered):           + two-sided single odds
     fetch 2^N cells → sanity gate → probit devig    (ride along, no extra
   Route B (any N, or any one-sided leg):             request — Route B input)
     fetch 1 SGP price → × correlation transfer
   ▼
 ResultStore[hash] = {book: fair, route, landed_at}   (quotable ≤15s,
                                                       then re-fetch)
 accept ─► confirm tick: synchronous live re-fetch (priority lane)
           + drift check → confirm, else VOID
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
   51% × 60% = 15.9%; multiplier 16.7/15.9 = 1.05; fair = 49%×48%×57% × 1.05
   = 14.1%). If the *single* is one-sided too, its devig uses the book's
   existing vig-fallback haircut (`DK_VIG_FALLBACK` etc.); still
   consensus-gated downstream.
3. Tick ≈ landing+≤2s (RFQ still open — median rest ≫ this): live result →
   `consensus()` over per-book fairs → all existing gates → quote.
4. While the RFQ stays open, the discovery tick keeps re-pricing it (that is
   existing behavior for all quotes): the result ages past 15s → re-fetch →
   landing re-prices → the existing hysteresis/replace logic
   (`main.py:632-649`) leaves the quote if fair is unchanged or replaces it
   if it moved. Our quote tracks the books at feed cadence (~15–25s) the same
   way grid quotes track the 60s sweep at 2s cadence. The RFQ closes (creator
   accepted someone, expired, cancelled) → it leaves the poll → the feed for
   it stops; our quote's terminal status lands via the confirm tick.
5. If OUR quote is accepted: confirm-tick last look triggers a **synchronous
   live re-fetch** of this combo (worker priority lane); drift check against
   the fresh consensus; can't complete in the confirm window or consensus
   fails → void (`voided_no_fresh_books`). We never confirm on a
   previously-fetched number.

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
  read out of the book's own structure payload — no extra requests or latency
  (all six structure fetches carry two-sided single odds for
  spread/total/ML, and the structure is already needed for leg resolution).
  Returns `∏ fair_i × (1/sgp_decimal) / ∏ vigged_i`. Assumption (measured,
  see §7): the book's SGP margin ≈ compounded leg vig, so it cancels in the
  ratio. Sanity gate: resulting multiplier within [0.5, 2.0], else None.
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

- `OnDemandEngine(service)` owning:
  - **Queue**: FIFO, unbounded, deduped on in-flight `leg_set_hash` (two RFQs
    on the same combo in the same window share one fetch — that's
    concurrency-dedup, not price reuse). No drop path: under load, fetches
    land later and RFQs that expired meanwhile simply go unquoted — degraded
    means *fewer or laggier quotes*, never *staler quotes*.
  - **Worker**: one background thread consuming the queue; per combo, fans out
    to the service's books concurrently (thread per book, same deadline
    discipline as `SGPService.refresh`). A **priority lane** serves confirm-
    time last-look re-fetches ahead of feed fetches.
  - **ResultStore**: `hash → {book: fair, route, landed_at}`. A result is
    quotable while `age <= QUOTE_FRESH_SEC` (constant 15s); past that, the
    next discovery tick that still sees the RFQ enqueues a re-fetch — the
    open-RFQ poll is the feed driver, so fetching stops the moment the RFQ
    leaves the poll. No TTL knobs, no standing per-quote refresh loop.
  - `lookup(hash) -> {book: fair} | None` — fresh results only; the pure read
    handed to the router. Single lock around store mutations; lookups copy.
  - `refetch_now(hash, deadline) -> {book: fair} | None` — synchronous
    priority re-fetch for the confirm-tick last look.
- Crash-safety: worker exceptions caught and logged; dead worker restarts
  lazily on next enqueue; restart loses nothing durable (results are
  ephemeral by design).

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

- `_priceable_in_phase1` → `_priceable(canon)`: `on_demand` routes are always
  in-scope (no switch, no leg-count condition).
- Discovery tick: on-demand games with a live (≤15s) result price through the
  normal path — including the existing hysteresis/replace of an already-open
  quote (`main.py:632-649`), which is what keeps our quote current at feed
  cadence; result aged/absent → `engine.ensure_fetch()` + skip reason
  `on_demand_pending` (new reason value; the monitor reads reason
  vocabularies from data, so it appears automatically). Because the tick
  re-visits every open RFQ, this single rule yields the continuous feed —
  there is no separate scheduler.
- Confirm tick: for combos containing an on-demand game, the last look calls
  `engine.refetch_now(hash, deadline)` (priority lane, synchronous); on
  timeout or consensus failure → existing `voided_no_fresh_books` path. Grid
  combos keep today's behavior.
- Risk sweep: **zero changes.** Tipoff and global books-stale pulls apply to
  all quotes; the drift-since-quote check works whenever a fresh on-demand
  result happens to be in the store and no-ops otherwise (`main.py:941`
  requires `cur_med is not None`) — drift defense for on-demand lives in the
  feed's hysteresis/replace plus the confirm last look.
- Fill research event gains `consensus_books` and `route` per book.
- **Instrumentation (measurement blind-spot fix):** out-of-scope MLB RFQs now
  store `legs_json` and a shape class in `seen_rfqs`, with distinct
  `last_decision` values (`out_of_scope_non_mlb`, `out_of_scope_lone_single`,
  …). New research events: `on_demand_requested`, `on_demand_result`
  (per-book route, cell coverage, fairs, latency, consensus outcome,
  `route_gap` when both routes computed).

### 4.8 Config

**None.** No new config keys (user decision, rev 3). `QUOTE_FRESH_SEC=15` and
`PARTITION_OVERROUND_MAX=1.35` are module constants (correctness guards, not
knobs); pacing is per-book serialization; the existing bot-wide kill file
(`config.KILL_FILE`) remains the emergency stop; removing a misbehaving book
from on-demand is a one-line code change shipped on evidence.

## 5. Failure modes (all fail toward "fewer/laggier quotes", never "staler quotes")

| Failure | Behavior |
|---|---|
| Leg's side not offered at a book (one-sided line) | Route A impossible at that book → Route B (single SGP call + correlation transfer) |
| Single itself one-sided at the book | Route B single devig uses the book's vig-fallback haircut; consensus still gates |
| Book prices some cells but not all 2^N | Route A → Route B fallback; if the SGP call also fails, book drops out |
| Partition overround / multiplier insane | sanity gate → book drops out |
| Book rejects the combo (DK combinability, MGM no-group) | book drops out |
| < MIN_AGREEING_BOOKS survive | consensus None → `no_fair` skip; open quote stays until feed or last look resolves it |
| Fetch in flight when RFQ re-seen | shared fetch (dedup); still `on_demand_pending` |
| Result ages past 15s while RFQ open | next tick enqueues re-fetch (the feed); landing re-quotes via hysteresis/replace |
| Fair moves between feed landings (~15–25s window) | next landing replaces the quote; if accepted inside the window, confirm-time live re-fetch + drift check → void if bad |
| RFQ closes (someone else taken / expired) | leaves the poll → feed stops; quote's terminal status via confirm tick |
| RFQ flood (June-style) | queue grows; prices land late; expired RFQs go unquoted; books see one serialized stream each |
| Confirm-time re-fetch can't complete in window | `voided_no_fresh_books` — never confirm on an old number |
| Worker thread dies | lazy restart on next enqueue; meanwhile all on_demand → skip |
| Caesars WAF trips | book 3-strike reinit (shared with grid sweep — accepted risk; retreat = code change) |
| Bot restart | engine state empty → open RFQs re-fetch on next poll (by design) |

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
   `runners` test seam) — in-flight dedup, 15s freshness death, priority-lane
   `refetch_now` (success/deadline), worker-death restart, flood behavior
   (queue grows, no drops, no stale results served), **feed invariants**:
   an RFQ present in the poll re-fetches each time its result ages past 15s;
   an RFQ absent from the poll generates zero further fetches.
4. **Router**: injected `on_demand_fairs` — consensus paths, None paths, and
   the untouched-grid regression suite (must pass unmodified, 1e-9 lock).
5. **Integration**: discovery→pending→landed→quoted flow with everything
   faked; feed re-quote via hysteresis/replace on a moved fair (and no
   replace inside hysteresis); confirm-tick synchronous re-fetch
   success/timeout/void; feed stop when the RFQ leaves the poll.
6. **Live verification (manual, before merge)**: one real exotic combo priced
   per book per route from a REPL, logged; then `--dry-run` soak with
   synthetic self-created RFQs (as the taker does) to exercise the full path
   end-to-end, including a deliberately one-sided leg and a fair moved
   mid-RFQ.

## 7. Rollout & measurement

1. Merge → restart maker **from main repo cwd** (restart-gotchas rule applies
   to the maker) in `--dry-run`; confirm `on_demand_*` research events flow on
   synthetic RFQs.
2. Live-enable. Watch (monitor + research DB): on-demand RFQ count/shapes
   (finally measurable), per-book route mix (A vs B) and partition-completion
   rate, `route_gap` distribution (validates the Route B vig-cancellation
   assumption), feed latency (enqueue→landing) and effective re-quote
   cadence, `consensus_books` composition (DK+NV-only share), Caesars WAF
   failure rate vs grid-sweep health, void rate on on-demand quotes vs grid
   quotes.
3. Pre-registered retreat triggers (all code changes shipped on evidence, no
   configs): Caesars grid-sweep failures attributable to on-demand → drop CZR
   from the on-demand book set; `route_gap` shows Route B systematically rich
   → haircut or disable Route B; feed cadence proves too slow on fast markets
   (void cluster) → shorten `QUOTE_FRESH_SEC` / widen on-demand margin.

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
  the RFQ-poll-driven feed model, failure modes, rollout/retreat playbook.
- `NFLWork/CLAUDE.md`: maker bullet — Phase 2 on-demand pricing summary.
- `kalshi_mlb_monitor/README.md`: only if a monitor tweak is needed (reason
  vocabulary is data-driven, so likely none).
- Memory file: Phase 2 status update.
