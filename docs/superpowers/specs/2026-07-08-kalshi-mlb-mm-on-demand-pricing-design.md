# Kalshi MLB MM — Phase 2: On-Demand Pricing of Arbitrary Leg Combos

**Date:** 2026-07-08 (rev 5, 2026-07-10 — adversarial-review fixes: per-book
singles sourcing corrected, Fréchet sanity bounds, N-aware overround gate, no
sweep cross-lock, confirm deadline budget, multi-game refetch; rev 4
RFQ-poll-driven live feed; rev 3 removed configs; rev 2 live-per-RFQ +
one-sided route)
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
   (rejected, rev 3 — incoherent with the discovery loop). Residual drift
   exposure = one feed interval, backstopped by the confirm-time live
   re-fetch.
3. **Two devig routes, no leg cap** (rev 2; sourcing corrected rev 5).
   Route A (gold standard, ≤3 legs, all cells offered): fetch the full 2^N
   side-combination partition, probit-devig, read the target cell — the exact
   generalization of Phase 1's 4-cell rule. Route B (any N, and any leg
   one-sided at the book): book-implied correlation transfer — one SGP
   call/book; fair = ∏(devigged singles) × [SGP implied / ∏(vigged singles)].
   A lone SGP price cannot be devigged by itself (fair and margin are
   entangled in one number); the partition and the book's own singles are the
   only two margin references available — hence exactly these two routes.
   Singles sourcing is per-book (rev 5): PX/NV/MGM/CZR market trees carry
   two-sided single odds (free); FD's structure payload carries them but the
   current parser discards them — capture is extended (free); DK's payload
   has selection IDs only — singles priced via 1-leg `calculateBets` calls
   (2 per leg), made only when Route B is actually taken. Rejected:
   partition-only with a leg cap (rev 1).
4. **All 6 books, always on, zero new config** (user decisions, revs 1+3).
   No switches: the on-demand path is a permanent part of the pricing engine,
   books = the SGP service's book set. The existing bot-wide kill file covers
   emergencies; dropping a misbehaving book is a one-line code change.
   Module constants (correctness guards, not knobs): `QUOTE_FRESH_SEC=15`,
   `PARTITION_OVERROUND_PER_LEG=0.12`.
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
- **Feed cadence is ~15–25s between sweeps — and worse across them.** The
  60s grid sweep runs INLINE in the main loop and blocks it for ~30–40s
  (existing architecture); during that window no RFQs are discovered,
  re-quoted, or re-fetch-enqueued — for grids and on-demand alike. Worst-case
  on-demand re-quote gap ≈ sweep length + one fetch (~50–60s). Defenses are
  the same as grids': margin, hysteresis, circuit breaker, and the confirm
  last look (here a live re-fetch). If books reprice faster on news (lineup
  scratch), voids are the cost.
- **Accept bursts can void.** The confirm-time live re-fetch is synchronous
  in the confirm tick and budgeted against the ~30s confirm window; several
  simultaneous accepts on on-demand combos are re-fetched sequentially, and
  later ones that run out of window are voided (never confirmed on an old
  number). This is the same failure the old in-confirm positions-retry hit
  (removed as N5) — reintroduced here deliberately, but deadline-budgeted and
  fail-safe. At observed accept rates (few per week) this is theoretical.
- **No rate limiting anywhere.** Pacing is one on-demand combo in flight per
  book. Steady-state cost: each open on-demand RFQ ≈ 2–4 fetch rounds/min for
  its open life — and at ProphetX every fetch is a REAL RFQ on their
  exchange, so an open combo shows as 2–4 PX RFQs/min until it closes. In a
  June-style flood (~35k RFQs/day) the fetch queue grows and prices land late
  — quotes forgone or lagged, never artificially refused. If you'd rather
  shed load explicitly under flood, say so and a backstop returns.
- **Caesars blast radius.** On-demand shares Caesars' WAF-token session with
  the 60s grid sweep; tripping the WAF costs the grid pipeline Caesars too.
  Note CZR will usually price via Route B anyway (8 cells × 30s timeouts
  rarely fit the per-book deadline). Retreat is a one-line code change.
- **Route B rests on a vig-cancellation assumption** (the book's SGP margin ≈
  compounded leg vig, so it divides out of the correlation ratio). Extra SGP
  margin beyond that inflates Route B fairs slightly. Defenses: Fréchet
  bounds (rev 5) reject logically impossible fairs, the consensus band
  rejects outliers, and a partition-vs-ratio `route_gap` research metric is
  emitted wherever both routes come free — so the assumption is measured.
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
2. *Fréchet bounds.* Whatever the correlation, a joint probability can never
   exceed its smallest component (P(A∩B) ≤ min) nor fall below
   max(0, ΣP − (n−1)). These are exact, assumption-free limits — the right
   sanity gate for an estimated joint fair, unlike an arbitrary multiplier
   cap that would reject legitimately high-correlation combos (e.g. run line
   + moneyline stacks).
3. *RFQs are auctions, not order books.* A quote is an answer to one
   creator's open request, competing with other makers' answers; it exists
   only while the request does. The maker's job while the auction runs is to
   keep its answer current (re-quote on movement); its protection at the
   moment of truth is the last look. Both already exist in this bot — Phase 2
   just puts on-demand combos on the same treadmill.

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

### 1.3 Feasibility (per-book wire audit, 2026-07-08; singles column rev 5)

The wire layer mostly exists; the resolution layer does not.

| Book | Arbitrary-N wire call | Leg→selection lookup | Singles odds for Route B | On-demand caveats |
|---|---|---|---|---|
| DraftKings | `post_calculate_bets(session, sel_ids)` (`scraper_draftkings_trifecta.py:104`) | sel-id dicts: spreads `(sign, abs_line, participant)`, totals `("O"/"U", line)`, ML `{"1","3"}` (`scraper_draftkings_sgp.py:420,576`) | **not in structure** (regex-extracted IDs only) → 1-leg `calculateBets` ×2 per leg, Route-B-only | rejects some pairings (`combinabilityRestrictions`, 422); 2MB structure payload (TTL-cached) |
| FanDuel | `price_combo` implyBets — body is a list; needs 2→N loop (`scraper_fanduel_sgp.py:452`) | runners `(side, line) → (marketId, selectionId)` (`scraper_fanduel_sgp.py:303`) | in the same event-page payload but currently DISCARDED → extend `fetch_event_runners` capture (verify runner odds field live) | PerimeterX token long-lived but manual-refresh |
| ProphetX | `submit_parlay_rfq(legs: list)` (`prophetx_client.py:138`) | market tree + ML selections (`prophetx.py:478,695`) | in market tree (`sel["odds"]`, cf. `_safe_leg_decimal`) | **each query = a real RFQ on their exchange**; teaser-tier filter |
| Novig | `submit_parlay(outcome_ids: list)` (`novig_client.py:157`) | GraphQL market tree → home/away/over/under/ml legs (`scraper_novig_sgp.py:440`) | in market tree (legs carry `available` implied prob) | mirrors DK prices (weak independence) |
| BetMGM | `price_picks(fixture, legs: list)` (`betmgm_client.py:172`) | `parse_markets` → `(market_id, option_id, decimal_odds)` incl. ML (`betmgm.py:95`) | in structure (leg tuple carries decimal) | declines un-comboable picks silently (no pricing group) |
| Caesars | `price_combo(legs: list)` (`caesars_client.py:228`) | `parse_markets` → leg dicts with price (`caesars.py:114`; away-spread line negation `caesars.py:42`) | in structure (`sel.price.d`) | **AWS-WAF token (~1.5s mint, ~4min TTL), 3-wide cap** — shared with grid sweep |

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
multiplication, gates, sizing, and lifecycle are unchanged — including the
RFQ-level ≥2-total-legs rule (lone singles stay excluded; the cap that Phase 2
removes is the per-game SHAPE restriction, not that rule).

**Non-goals.**
- No new market types: legs are still MLB FG spread / total / moneyline
  (`KXMLBSPREAD-`/`KXMLBTOTAL-`/`KXMLBGAME-`). Player props, F5, etc. stay
  unpriceable (the books' structure caches don't carry them; separate project).
- No persistence of on-demand prices into `mlb_sgp_odds` (results live in the
  engine + research events; schema untouched).
- No change to taker (`kalshi_mlb_rfq`) — maker only.

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
   unbounded — pacing = one combo                    book, game, legs)
   in flight per book, never refusal)                │ reuses persistent client
   │ background worker                               │ + structure TTL caches
   ▼                                                 │ (concurrent with the
 per book (concurrent across books):                 │  grid sweep — matches
   Route A (≤3 legs, all cells offered):             │  existing in-book
     fetch 2^N cells → sanity gate → probit devig    │  shared-session use)
   Route B (any N, or any one-sided leg):            ▼
     fetch 1 SGP price → × correlation          resolve_legs(structure, legs)
     transfer vs singles                             → book selection IDs
   ▼                                                 → singles odds (or DK
 ResultStore[hash] = {book: fair, route, landed_at}    1-leg-call fallback)
                                        (quotable ≤15s, then re-fetch)
 accept ─► confirm tick: synchronous live re-fetch of ALL the combo's
           on-demand games (priority lane, budgeted against the confirm
           window) + drift check → confirm, else VOID
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
   devigged per §4.2's sourcing (e.g. SGP implied 16.7%; vigged singles 52% ×
   51% × 60% = 15.9%; multiplier 16.7/15.9 = 1.05; fair = 49%×48%×57% × 1.05
   = 14.1%). If the *single* is one-sided too, its devig uses the book's
   existing vig-fallback haircut (`DK_VIG_FALLBACK` etc.); still gated by
   Fréchet bounds + consensus downstream.
3. Tick ≈ landing+≤2s (RFQ still open — median rest ≫ this): live result →
   `consensus()` over per-book fairs → all existing gates → quote.
4. While the RFQ stays open, the discovery tick keeps re-pricing it (existing
   behavior): the result ages past 15s → re-fetch → landing re-prices → the
   existing hysteresis/replace logic (`main.py:632-649`) leaves the quote if
   fair is unchanged or replaces it if it moved. Feed cadence ~15–25s between
   grid sweeps; during the inline ~30–40s sweep each minute the whole loop
   (grids included) pauses. A combo spanning several on-demand games needs
   all of them fresh in the same tick — their re-fetches co-trigger on the
   same ticks, so landings converge after a round or two.
5. If OUR quote is accepted: confirm-tick last look re-fetches ALL the
   combo's on-demand games synchronously (worker priority lane), budgeted
   against the remaining confirm window; drift check against the fresh
   consensus; overrun or consensus failure → void (`voided_no_fresh_books`).
   We never confirm on a previously-fetched number.

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
  probabilities; returns the target cell's devigged probability. Sanity gate
  (rev 5, N-aware — SGP vig compounds per leg, so a fixed cap would reject
  healthy 8-cell partitions): `1.0 <= sum(1/d) <= 1.0 + 0.12·N`
  (`PARTITION_OVERROUND_PER_LEG = 0.12`; N=2 → 1.24, N=3 → 1.36). Any
  missing/insane cell → None; the caller then tries Route B.
- `fair_by_correlation_transfer(sgp_decimal, singles) -> float | None` — Route
  B: `singles` is per-leg `(vigged_implied_of_chosen_side, devigged_fair)`.
  Returns `∏ fair_i × (1/sgp_decimal) / ∏ vigged_i`. Assumption (measured,
  see §7): the book's SGP margin ≈ compounded leg vig, so it cancels in the
  ratio. Sanity gate (rev 5): **Fréchet bounds** — the result must satisfy
  `max(0, Σ fair_i − (n−1)) ≤ fair ≤ min(fair_i)`, else None. (Exact,
  assumption-free joint-probability limits; replaces rev 4's arbitrary
  multiplier cap, which wrongly rejected legitimate high-correlation stacks
  like run line + moneyline where the multiplier exceeds 2.)
- Singles sourcing (rev 5, per book): PX (`sel["odds"]`), NV (`available`),
  MGM (leg decimal), CZR (`price.d`) come free from the structure already
  fetched for leg resolution. FD: same event-page payload carries runner
  odds; `fetch_event_runners` is extended to capture them (field name
  verified live at implementation). DK: structure has IDs only → each leg's
  two sides priced via 1-leg `calculateBets` (2 calls/leg), executed ONLY
  when Route B is actually taken at DK.
- Route selection per (book, combo): A if N ≤ 3 and all 2^N cells priced,
  else B, else book drops out. `route_gap` (A-vs-B comparison) is emitted to
  research ONLY where both routes come free (never extra DK calls just for
  research).

### 4.3 `mlb_sgp/*` — per-book leg resolution + single-combo pricing

Per book, two additions to the existing module (no new scrapers):

- `resolve_legs(structure, legs: list[CanonicalLeg]) -> list[SelectionRef] | None`
  — maps canonical legs to that book's selection descriptors using the SAME
  structure objects the orchestrators already fetch and cache, returning per
  leg: chosen-side ref, opposite-side ref (partition + singles devig), and
  singles odds where the structure carries them (§4.2). None for any
  unresolvable chosen side; a missing OPPOSITE side resolves with
  `opposite=None` (that's what routes the book to Route B). Extraction/
  refactor of logic that exists in closures — TDD'd against recorded
  structure fixtures, then live-verified per book.
- `price_selection_set(client, refs) -> float | None` — one wire call, decimal
  odds or None. Mostly exists; FD needs the 2→N loop, DK uses the trifecta
  N-leg call. A single ref prices a single (used for DK Route B singles).

Existing orchestrator paths (`price_sgps`) are untouched — the grid sweep and
its outputs stay byte-identical.

### 4.4 `kalshi_common/sgp_service.py` — single-combo entry point

`SGPService.price_on_demand(book, game_key, legs) -> OnDemandBookResult | None`
(cells or single SGP price + singles, per the route) — runs on the calling
(worker) thread, reuses the book's persistent client and structure TTL caches
(a cold game costs one structure fetch, then the price calls). Shares
`refresh()`'s failure accounting (3 strikes → client reinit) so on-demand and
grid-sweep health share one recovery path. Concurrency (rev 5): on-demand
runs CONCURRENTLY with the grid sweep on the same session — this matches
existing practice (each book's sweep already fans one session across a
thread pool, e.g. `draftkings.py:410` at 8–32 wide), so no cross-lock; a
cross-lock would stall on-demand fetches ~30–75s during every sweep and
break the feed cadence. Pacing: at most ONE on-demand combo in flight per
book (the engine's per-book width), so a book sees the sweep plus one
serialized on-demand stream.

### 4.5 `kalshi_mlb_mm/on_demand.py` — the stateful engine (new)

- `OnDemandEngine(service)` owning:
  - **Queue**: FIFO, unbounded, deduped on in-flight `leg_set_hash` (two RFQs
    on the same combo in the same window share one fetch — that's
    concurrency-dedup, not price reuse). No drop path: under load, fetches
    land later and RFQs that expired meanwhile simply go unquoted — degraded
    means *fewer or laggier quotes*, never *staler quotes*.
  - **Worker**: one background thread consuming the queue; per combo, fans out
    to the service's books concurrently (thread per book, same deadline
    discipline as `SGPService.refresh`); at most one combo in flight per book.
    Known limit (accepted): a slow combo delays the next in queue — at
    observed volumes (≈0 on-demand RFQs/week) a single consumer is right;
    widening to a small pool is a follow-up if usage appears.
  - **ResultStore**: `hash → {book: fair, route, landed_at}`. A result is
    quotable while `age <= QUOTE_FRESH_SEC` (constant 15s); past that, the
    next discovery tick that still sees the RFQ enqueues a re-fetch — the
    open-RFQ poll is the feed driver, so fetching stops the moment the RFQ
    leaves the poll. No TTL knobs, no standing per-quote refresh loop.
  - `lookup(hash) -> {book: fair} | None` — fresh results only; the pure read
    handed to the router. Single lock around store mutations; lookups copy.
  - `refetch_now(hashes, deadline) -> bool` (rev 5: plural) — synchronous
    priority re-fetch of ALL of a combo's on-demand game sub-combos for the
    confirm-tick last look; returns success only if every hash lands fresh
    within `deadline`. If a feed fetch for a hash is already in flight, its
    landing is awaited (within the deadline) rather than duplicated.
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
band (user decision — no independence adjustment). For the `consensus_books`
research field, a pure sibling `consensus_detail()` returns
`(fair, agreeing_books)`; `consensus()` itself is untouched.

### 4.7 `kalshi_mlb_mm/main.py` — integration + instrumentation

- `_priceable_in_phase1` → `_priceable(canon)`: `on_demand` routes are always
  in-scope. The RFQ-level `len(canon) >= 2` rule stays (lone singles remain
  excluded — see §2); what's removed is only the per-game shape restriction.
- Discovery tick: on-demand games with a live (≤15s) result price through the
  normal path — including the existing hysteresis/replace of an already-open
  quote (`main.py:632-649`), which is what keeps our quote current at feed
  cadence; result aged/absent → `engine.ensure_fetch()` (per on-demand GAME
  sub-combo) + skip reason `on_demand_pending` (new reason value; the monitor
  reads reason vocabularies from data, so it appears automatically). Because
  the tick re-visits every open RFQ, this single rule yields the continuous
  feed — there is no separate scheduler. Note `on_demand_pending` decisions
  repeat per tick while pending — same volume behavior as today's repeated
  `no_fair` skips (precedented; research DB already absorbs ~84
  decisions/RFQ).
- Confirm tick: for combos containing on-demand games, the last look calls
  `engine.refetch_now(hashes, deadline)` with `deadline` = remaining confirm
  window minus a safety buffer for the confirm API call itself; on overrun or
  consensus failure → existing `voided_no_fresh_books` path. Accept bursts
  re-fetch sequentially; later accepts that run out of window void
  (fail-safe, disclosed in Review Pack). Grid combos keep today's behavior.
- Risk sweep: **zero changes.** Tipoff and global books-stale pulls apply to
  all quotes; the drift-since-quote check works whenever a fresh on-demand
  result happens to be in the store and no-ops otherwise (`main.py:941`
  requires `cur_med is not None`) — drift defense for on-demand lives in the
  feed's hysteresis/replace plus the confirm last look.
- Fill research event gains `consensus_books` and per-book `route`.
- **Instrumentation (measurement blind-spot fix):** out-of-scope MLB RFQs now
  store `legs_json` and a shape class in `seen_rfqs`, with distinct
  `last_decision` values (`out_of_scope_non_mlb`, `out_of_scope_lone_single`,
  …). New research events: `on_demand_requested`, `on_demand_result`
  (per-book route, cell coverage, fairs, latency, consensus outcome,
  `route_gap` where free).

### 4.8 Config

**None.** No new config keys (user decision, rev 3). Module constants
(correctness guards, not knobs): `QUOTE_FRESH_SEC = 15`,
`PARTITION_OVERROUND_PER_LEG = 0.12`; Route B is bounded by Fréchet limits
(no constant needed). Pacing is one on-demand combo in flight per book; the
existing bot-wide kill file (`config.KILL_FILE`) remains the emergency stop;
removing a misbehaving book from on-demand is a one-line code change shipped
on evidence.

## 5. Failure modes (all fail toward "fewer/laggier quotes", never "staler quotes")

| Failure | Behavior |
|---|---|
| Leg's side not offered at a book (one-sided line) | Route A impossible at that book → Route B (single SGP call + correlation transfer) |
| Single itself one-sided at the book | Route B single devig uses the book's vig-fallback haircut; Fréchet + consensus still gate |
| Book prices some cells but not all 2^N | Route A → Route B fallback; if the SGP call also fails, book drops out |
| Partition overround outside N-aware bound / fair outside Fréchet bounds | sanity gate → book drops out |
| Book rejects the combo (DK combinability, MGM no-group) | book drops out |
| < MIN_AGREEING_BOOKS survive | consensus None → `no_fair` skip; open quote stays until feed or last look resolves it |
| Fetch in flight when RFQ re-seen | shared fetch (dedup); still `on_demand_pending` |
| Result ages past 15s while RFQ open | next tick enqueues re-fetch (the feed); landing re-quotes via hysteresis/replace |
| Fair moves between feed landings | next landing replaces the quote; if accepted inside the window, confirm-time live re-fetch + drift check → void if bad |
| Grid sweep blocking the main loop (~30–40s/min) | no discovery/re-quotes for anything during the window (existing architecture); on-demand worst-case re-quote gap ≈ sweep + fetch |
| Multi-on-demand-game combo, landings staggered | joint freshness required at one tick; co-triggered re-fetches converge within a round or two |
| Accept burst on on-demand combos | sequential confirm re-fetches; later accepts that exhaust the confirm window void (never confirm stale) |
| RFQ closes (someone else taken / expired) | leaves the poll → feed stops; quote's terminal status via confirm tick |
| RFQ flood (June-style) | queue grows; prices land late; expired RFQs go unquoted; books see sweep + one serialized on-demand stream each |
| Worker thread dies | lazy restart on next enqueue; meanwhile all on_demand → skip |
| Caesars WAF trips | book 3-strike reinit (shared with grid sweep — accepted risk; retreat = code change) |
| Bot restart | engine state empty → open RFQs re-fetch on next poll (by design) |

## 6. Testing (TDD throughout)

1. **Pure math first**: `enumerate_partition` (side flipping, ordering, hash
   stability, same-market-both-sides edge → impossible target cell);
   `devig_partition` (4-cell case must equal existing `devig_book` on the
   same grid — a second regression bridge to Phase 1; N-aware overround gate
   edges); `fair_by_correlation_transfer` (vig-cancellation on synthetic
   books, Fréchet rejection cases incl. nested totals, one-sided-single
   haircut path); route selection (A↔B fallbacks).
2. **Resolution layer**: per-book `resolve_legs` against recorded structure
   fixtures (repo's existing fixture pattern); unresolvable chosen side →
   None; missing opposite side → `opposite=None` (Route B path); singles
   extraction both sides (and FD odds capture).
3. **Engine**: fake `SGPService` (constructor-injected, mirrors the existing
   `runners` test seam) — in-flight dedup, 15s freshness death, per-book
   width 1, priority-lane `refetch_now` (multi-hash success / partial /
   deadline overrun / await-in-flight), worker-death restart, flood behavior
   (queue grows, no drops, no stale results served), **feed invariants**: an
   RFQ present in the poll re-fetches each time its result ages past 15s; an
   RFQ absent from the poll generates zero further fetches.
4. **Router**: injected `on_demand_fairs` — consensus paths, None paths,
   `consensus_detail` agreement bookkeeping, and the untouched-grid
   regression suite (must pass unmodified, 1e-9 lock).
5. **Integration**: discovery→pending→landed→quoted flow with everything
   faked; feed re-quote via hysteresis/replace on a moved fair (and no
   replace inside hysteresis); confirm-tick multi-game re-fetch
   success/timeout/void; feed stop when the RFQ leaves the poll.
6. **Live verification (manual, before merge)**: one real exotic combo priced
   per book per route from a REPL, logged; then a `--dry-run` soak with
   synthetic self-created RFQs (as the taker does) to exercise the full path
   end-to-end, including a deliberately one-sided leg and a fair moved
   mid-RFQ. **Dry-run only** — never live-quote-and-confirm our own RFQs
   (self-trade/wash concerns); live enablement waits for organic flow.

## 7. Rollout & measurement

1. Merge → restart maker **from main repo cwd** (restart-gotchas rule applies
   to the maker) in `--dry-run`; confirm `on_demand_*` research events flow on
   synthetic RFQs.
2. Live-enable. Watch (monitor + research DB): on-demand RFQ count/shapes
   (finally measurable), per-book route mix (A vs B) and partition-completion
   rate, `route_gap` distribution (validates the Route B vig-cancellation
   assumption), feed latency (enqueue→landing) and effective re-quote
   cadence, `consensus_books` composition (DK+NV-only share), Caesars WAF
   failure rate vs grid-sweep health, PX on-demand RFQ volume, void rate on
   on-demand quotes vs grid quotes.
3. Pre-registered retreat triggers (all code changes shipped on evidence, no
   configs): Caesars grid-sweep failures attributable to on-demand → drop CZR
   from the on-demand book set; `route_gap` shows Route B systematically rich
   → haircut or disable Route B; feed cadence proves too slow on fast markets
   (void cluster) → shorten `QUOTE_FRESH_SEC` / widen on-demand margin;
   sustained on-demand volume → widen the engine's worker pool.

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
