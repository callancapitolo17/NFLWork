# Kalshi MLB MM — Generalized N-Leg Combo Pricing

**Status:** Design — awaiting user review
**Date:** 2026-06-28
**Branch:** `worktree-kalshi-mm-generalized-combos`
**Supersedes scope of:** `docs/superpowers/specs/2026-05-26-kalshi-mlb-mm-design.md` (the fixed 2-leg grid maker)

---

## Review Pack

**What we're building**

Today the Kalshi MLB maker only prices *exactly* 2-leg full-game combos of one fixed shape (spread×total or ml×total). We're turning it into a true RFQ maker that prices **any combination of any number of legs** drawn from the three core game-line markets — spread, moneyline, total. Cross-game legs are priced as independent (multiply); same-game correlated legs are priced by forwarding the exact leg set to each sportsbook's own SGP builder on-demand, devigging across the leg-set's outcome corners, and taking a multi-book consensus.

**Key decisions**

1. **Books-as-oracle, on-demand, for novel same-game combos.** When an in-scope RFQ arrives, forward its exact legs live to each book's SGP endpoint, devig, consensus. *Rejected:* building the internal bivariate-Poisson model now — it's the long-term moat but more work, and the books already price correlation on demand. The internal model is deferred behind the same `price_legset` interface (Phase 3).
2. **Sided cross-product devig (no reduction logic).** To devig an arbitrary combo we enumerate its `2^k` outcome corners by flipping each leg's side, price the ones the book will price, and probit-devig the set. *Rejected:* our own per-axis redundant-leg/contradiction reduction — the book's pricer already handles redundancy and impossible corners natively, so reduction code would be redundant with the book and pure risk surface.
3. **Hash-keyed canonical leg-set as the identity.** A combo is identified by `sha1` of its sorted, normalized legs. On-demand prices live in a new `mlb_sgp_legset_odds` table keyed by `(leg_set_hash, bookmaker, source)`. *Rejected:* extending the fixed `spread_line/total_line` columns — they don't generalize to N legs and a NULL-sentinel scheme gets ugly fast.
4. **Keep the existing grid table as a pre-scraped cache.** `mlb_sgp_odds` (and its 60s scrape) stays; exact 2-leg spread×total / ml×total RFQs still hit it for free, only novel shapes go on-demand. *Rejected:* a hard migration of the grid into the hash-keyed store — the MLB dashboard + R pipeline read `mlb_sgp_odds` with its current columns, so migrating is large blast radius for zero correctness gain. Folding the grid in is noted as future cleanup.
5. **Per-game consensus, then multiply across games.** Each game's sub-combo gets its own median-of-agreeing-books fair; independent games' fairs multiply. *Rejected:* requiring one book to cover all games then multiply-per-book — that needlessly drops combos when books cover different slates.

**Risks / push back here**

- **DK ≡ Novig is ~one source, and we're keeping `MIN_AGREEING_BOOKS=2`** (your call). A flexible bot leans harder on consensus to validate correlated prices, so at 2 a DK+Novig "agreement" is false confidence. Documented as an accepted risk with an independence-aware-consensus ready-lever — but this is the top pricing-integrity exposure of the project. Flag if you'd rather raise the bar now.
- **On-demand latency vs. an unknown RFQ window.** Kalshi exposes no RFQ TTL; we don't actually know how long we have to respond. The design fails safe (skip if too few books price in time) and *measures* the real window, but if it turns out to be very short, on-demand same-game coverage may be thin until we pre-warm popular leg-sets. Acceptable for a measurement phase?
- **Footprint.** More combo shapes = more distinct quotes = a larger leaked fair surface and more book API calls per RFQ. Mitigated by a short-TTL leg-set cache and price-once-per-leg-set, but it is a real increase over the 2-leg bot.
- **Corner cap is a coverage choice.** Capping at ≤4 same-game legs (≤16 corners) is a cost/coverage tradeoff, not a correctness one. If you want deeper same-game combos priced, the cap moves — at more book calls per RFQ.

**Worth understanding** (opt-in)

- **MECE partition + devig.** "Mutually exclusive and exhaustive" — every possible game outcome falls in exactly one corner of the cross-product, so the corners' implied probabilities sum to `1 + vig`. Dividing each by that sum removes the vig. This is the same idea as today's 4-cell grid, just generalized to `2^k` corners. (In R terms: it's like `prop.table()` over a complete contingency table — you can only normalize to 1 if the cells partition the whole space.)
- **Content-addressing via a hash.** Using `sha1(canonical_legs)` as the row key means identical combos always collide to the same key regardless of leg order — the hash *is* the identity. (Like using a `digest::digest()` of a sorted key vector instead of a hand-built composite string key.)

---

## 1. Background & problem

The current maker (`kalshi_mlb_mm/`) is a reactive RFQ daemon: it polls open Kalshi RFQs, and for each one it looks up a **pre-scraped 4-cell grid** in an in-memory cache (`_SGP_ODDS`, refreshed every 60s), devigs the cell matching the RFQ's legs, takes a book-consensus median, and quotes at `p / (1 + TARGET_ROI)`.

This works *only* because the two supported shapes — spread×total and ml×total — each form a clean Home/Away × Over/Under 2×2 partition that devigs cleanly via `devig_book`. `combo_descriptor` (`kalshi_common/leg_types.py`) recognizes exactly those two shapes and returns `None` for everything else. That `None` is the deliberate fail-safe that has kept the bot correct: anything it can't represent, it doesn't quote.

We want to price **any** combo of spread/ml/total legs across **any** number of legs and games — `spread+total+ML` same game, `ML+ML` cross-game, `spread+spread` across two games, `total+total`, etc. — while preserving the fail-safe.

### What the recon established (2026-06-28)

- **Kalshi RFQ payloads carry no TTL/expiry field.** RFQ lifetime is a server-side constant, undocumented in the repo. The maker currently never prices on-demand — it reads the 60s cache. On-demand per-RFQ pricing is new behavior.
- **All six book SGP endpoints accept an arbitrary leg list** (`selectionsForYourBet` / `betLegs` / `outcome_ids` / `legs`). Novig, ProphetX, BetMGM, Caesars client wrappers already take a list; DK and FD wrappers hardcode two args but POST a list internally (signature-widening only).
- **The fixed-grid assumption lives only in the orchestrators** (`kalshi_common/sgp_runner.py::enumerate_kalshi_targets` + each `mlb_sgp/<book>.py::price_sgps` / `_price_one_target`), not in the clients. A single combo POST is ~1–10s.

## 2. Goals / non-goals

**Goals**
- Price arbitrary N-leg combos of {spread, moneyline, total} across one or more MLB games.
- Cross-game = independent multiply; same-game = book-correlated, devigged, consensus.
- Preserve the fail-safe: anything we cannot price soundly → no quote.
- Keep all existing risk gates (tipoff, freshness, circuit-breaker, last-look, cooldown, void-rate halt, reconciliation) unchanged.
- Lay the foundation so an internal correlation model can be slotted in behind one interface later.

**Non-goals (v1)**
- Internal bivariate-Poisson model (Phase 3, deferred behind `price_legset`).
- Periods other than FG (matches current maker scope; the architecture doesn't preclude F5/F7 later).
- Independence-aware consensus / DK≡Novig de-duplication (documented accepted risk + ready-lever).
- Player props or any market outside spread/ml/total (normalizer rejects → skip).
- Taker (`kalshi_mlb_rfq`) wiring — maker is the goal.

## 3. Architecture: the leg-set router

The maker becomes a router. The discovery loop, risk gates, quote submission, confirm, and reconciliation are all unchanged — only the *fair-value computation* between "RFQ is in scope" and "we have a fair probability" is replaced.

```
Incoming RFQ (mve_selected_legs = arbitrary list)
        │
        ▼
[1] NORMALIZE      parse each leg → typed CanonicalLeg
                   reject any untypeable leg → None (skip)        ← fail-safe
        │
        ▼
[2] PARTITION BY GAME   {gameA: [...], gameB: [...]}
        │
        ▼  per game sub-combo:
[3] CLASSIFY ROUTE
      • exact 2-leg grid shape  → cached grid (existing fast path)
      • else (novel same-game / single) → on-demand book SGP
      • untypeable / over cap   → None (skip whole combo)
        │
        ▼
[4] PRICE EACH GAME → per-game devigged fair per book
        │
        ▼
[5] PER-GAME CONSENSUS   median of agreeing books (band, MIN_AGREEING_BOOKS)
        │
        ▼
[6] COMBINE ACROSS GAMES   combo_fair = Π gameFair_i   (independence)
        │
        ▼
quote at p/(1+TARGET_ROI)  → existing gates + submission unchanged
```

**Fail-safe rule (unchanged contract):** if *any* game sub-combo cannot be priced to consensus, the whole combo returns `None` and we do not quote.

## 4. Canonical leg-set + hash schema

### 4.1 CanonicalLeg

A leg is normalized to a tuple:

```
CanonicalLeg = (game_id, market_type, period, line, side)
   market_type ∈ {"spread", "total", "ml"}
   period      = "FG"   (v1)
   line        = signed home-perspective number; None for ml
   side        = home-perspective outcome ("home"/"away" for spread/ml,
                 "over"/"under" for total) — Kalshi YES/NO is resolved to this
                 home-perspective orientation at parse time so the same economic
                 leg always canonicalizes identically regardless of which side
                 the ticker expressed
```

Parsing reuses the existing Kalshi ticker decoders in `kalshi_common/leg_types.py`
(`_leg_dict_to_typed`, `_moneyline_side`, `_event_codes_from_legs`,
`_MLB_CODE_TO_TEAM`). Any leg whose market_type can't be resolved → the whole
combo is unpriceable (skip).

### 4.2 Canonicalization + hash

- Sort legs by `(game_id, market_type, line, side)`.
- Serialize to compact JSON.
- `leg_set_hash = sha1(json)`.

Identical combos collide to the same hash regardless of RFQ leg order. The hash is the content-addressed identity used for caching, dedup, and exposure tracking.

### 4.3 Storage

New table in the maker's market DB (`kalshi_mlb_mm_market.duckdb`):

```sql
CREATE TABLE IF NOT EXISTS mlb_sgp_legset_odds (
    leg_set_hash   VARCHAR,      -- sha1 of canonical legs
    game_id        VARCHAR,      -- the single game these legs belong to
    legs_json      VARCHAR,      -- canonical leg-set JSON (audit / re-pricing)
    bookmaker      VARCHAR,
    source         VARCHAR,
    fair_prob      DOUBLE,       -- devigged correlated fair for this leg-set, this book
    raw_decimal    DOUBLE,       -- the book's vigged combo decimal (audit)
    n_corners      INTEGER,      -- corners priced in the devig (audit)
    fetch_time     TIMESTAMP,
    PRIMARY KEY (leg_set_hash, bookmaker, source)
);
```

`mlb_sgp_odds` (the grid) is **unchanged** and continues to serve the cached
fast-path and the dashboard. Each `mlb_sgp_legset_odds` row is per-game (a
multi-game combo never lands here whole — it's decomposed first), so the hash
covers a single game's leg subset.

## 5. The pricing engine

### 5.1 Cross-game = independent multiply

Legs in different games share no outcome, so their joint probability is the
product of per-game probabilities. After partitioning by game we price each
game's sub-combo independently, consensus each, and multiply:

```
combo_fair = Π_game  consensus_fair(game_sub_combo)
```

If any game's sub-combo fails consensus, the product is undefined → skip. This
is the math the archived `archive/odds-api-maker-expansion` prototype already
used for cross-game; we keep it and route same-game pieces through the oracle.

### 5.2 Same-game = sided cross-product devig (the correctness core)

A book's SGP endpoint returns the **vigged** correlated decimal for the exact
leg set. To get a fair probability we must remove the vig. We do this with the
**sided cross-product**, the rigorous generalization of today's 4-cell grid:

1. **Enumerate corners.** For a `k`-leg same-game sub-combo, flip each leg's
   side to produce the `2^k` outcome corners. These corners are MECE by
   construction — every game outcome lands in exactly one.
2. **Price each corner** through the book's on-demand SGP endpoint. Impossible
   corners (e.g. `home -1.5 = yes` AND `home ML = no`) the book declines to
   price; their true probability is 0, so they fall out with no special-case
   logic. *This is how redundancy and contradictions are handled — by the book,
   not by us.*
3. **Probit-devig** the corners that priced. Their raw implied probs sum to
   `1 + vig`; `devig_book`'s n-way probit (`kalshi_common/fair_value.py`)
   normalizes them to fair probabilities summing to 1.
4. **Read the target corner** — the one matching the RFQ's stated sides — as
   the book's devigged correlated fair for this leg-set.

Worked example — `home -1.5 + over 8.5` (`k=2`, today's grid, for continuity):

```
corners (flip each side):
  (-1.5, over)  (-1.5, under)  (+1.5, over)  (+1.5, under)
price all 4 → Σ implied = 1 + vig → probit-devig → take (-1.5, over)
```

Worked example — `home -1.5 + over 8.5 + home ML` (`k=3`):

```
8 corners. The impossible ones — any corner with "-1.5 yes" but "ML no",
or "-1.5 no" but "+... " contradictions — the book won't price. The ~6
priceable corners devig to fairs summing to 1; we read the target corner.
The redundant ML leg required zero handling from us.
```

**Why not cheaper devig methods** (recorded so we don't relitigate):
- *Single price + per-leg-vig decomposition* (1 POST/book): assumes the book
  applies vig per-leg multiplicatively. If it doesn't, we mis-devig and get
  picked off. Correctness > latency here.
- *Consensus-only, no devig*: vig is one-sided overround; it does **not** cancel
  across books. We'd systematically overprice our fair and bleed.

### 5.3 Corner cap + fail-safe

- Cap same-game legs at `MAX_SAME_GAME_LEGS` (default 4 → ≤16 corners). Beyond
  the cap → skip (fail-safe). This is a cost knob, not a correctness limit.
- If fewer than `MIN_AGREEING_BOOKS` books fully price all corners of a game's
  partition within the deadline → that game sub-combo fails → skip the combo.

### 5.4 Per-game consensus

Identical to today: median of the supplied book fairs, keep books within
`±BOOK_CONSENSUS_BAND` of the median, require `≥ MIN_AGREEING_BOOKS` survivors,
fair = median of survivors. Applied **per game sub-combo** (not per whole combo)
so different-slate book coverage doesn't sink a multi-game combo.

### 5.5 Sanity guard (generalized SANITY_MULT_RATIO)

For each book's corner-set, compare the book's correlated target-corner implied
prob to the **naive independent multiply** of that book's own devigged single
legs. If `book_correlated / naive_independent` is outside `[1/1.5, 1.5]`,
quarantine *that book's* price for *that combo* (don't let it into consensus).
This is the generalization of the ProphetX F5-Over defense — book SGP builders
have bugs, and a flexible bot ingests far more of them.

## 6. On-demand mechanics

- **`price_legset(legs)` interface** (`kalshi_common/sgp_oracle.py`): takes a
  canonical same-game leg-set, returns `{bookmaker: fair_prob}` (post-sanity,
  pre-consensus). Internally: enumerate corners → fan out corner×book POSTs
  concurrently under an overall deadline → devig per book → sanity-filter.
  The internal model (Phase 3) implements the same signature as another
  "book" and slots into consensus.
- **Concurrency + deadline.** Corner×book calls run concurrently behind one
  `LEGSET_PRICE_DEADLINE_SEC`. Gather what returns; missing books just reduce
  the consensus count (and may trip the `MIN_AGREEING_BOOKS` fail-safe).
- **Short-TTL leg-set cache.** Cache `mlb_sgp_legset_odds` rows by
  `leg_set_hash` for `LEGSET_CACHE_TTL_SEC`; reuse across RFQs so a burst on the
  same combo doesn't refire the corner fan-out or grow our footprint.
- **Client wrapper widening.** DK/FD wrappers gain a generic `legs: list`
  signature (bodies already POST a list). The other four already accept lists.
- **RFQ-window measurement.** Log `created_ts → first-quote latency` and watch
  for `quote_expired` / `rfq_closed` errors to empirically learn the window the
  docs don't give us.

## 7. Risk / exposure / footprint

- **Exposure key generalizes** from `combo_market_ticker` to `leg_set_hash`
  (per-combo cap, cooldown, and in-flight accounting all key on the hash). The
  Kalshi `market_ticker` is still stored for the actual quote/position calls.
- **All existing gates unchanged** — they act on the quote/fill, which is
  identical regardless of how the fair was computed.
- **Accepted risks (measured, not prevented in v1):**
  - DK≡Novig as ~one source at `MIN_AGREEING_BOOKS=2` — heightened for
    correlated combos; ready-lever = independence-aware consensus.
  - Larger leaked fair surface + more correlated exposure across more shapes —
    same posture as maker v1, now over more combos.
  - On-demand latency vs. unknown RFQ window — fail-safe + measure.

## 8. Module structure

| Unit | Job | Notes |
|---|---|---|
| `kalshi_common/legset.py` | normalizer: parse → CanonicalLeg, canonicalize+hash, partition by game, classify route, enumerate corners | pure functions, network-free, unit-testable |
| `kalshi_common/sgp_oracle.py` | `price_legset(legs)` → per-book devigged fair via generic client endpoints + sided-cross-product devig + sanity guard | reusable by taker later; Phase-3 model slots in here |
| `mlb_sgp/dk_client.py`, `fd_client.py` | widen price wrappers to `legs: list` | bodies already POST a list |
| `kalshi_mlb_mm/db.py` | add `mlb_sgp_legset_odds` table + helpers | market DB |
| `kalshi_mlb_mm/fairs.py` / `main.py::_book_fairs` | route via normalizer; multiply across games; consensus per game | replaces direct grid lookup |
| `kalshi_common/leg_types.py::combo_descriptor` | thin shim over the normalizer for the cached 2-leg fast-path | preserves the existing grid hit |

## 9. Phasing

One spec, staged commits:

1. **Phase 1 — Foundation, no new network behavior.** Normalizer + canonical/hash
   + cross-game multiply + existing grid fast-path. Novel same-game shapes
   fail-safe skip. End-to-end testable with the current 60s cache only.
2. **Phase 2 — On-demand.** `price_legset` + sided-cross-product devig + sanity
   guard + leg-set cache + `mlb_sgp_legset_odds` + DK/FD wrapper widening. This
   is where live per-RFQ book calls turn on.
3. **Phase 3 — Deferred.** Internal bivariate-Poisson model implementing
   `price_legset` as an additional "book." Out of scope for this spec; the
   interface is built to receive it.

## 10. Testing

- **Unit (network-free):** `legset.py` — parsing, canonicalization, hash
  stability across leg order, game partitioning, corner enumeration (corner
  count = `2^k`), route classification, cap/skip behavior.
- **Devig math:** sided-cross-product devig on synthetic corner sets — sum-to-one
  invariant, target-corner extraction, probit parity with the existing
  `devig_book` on the 4-cell case (must reproduce today's grid number exactly).
- **Sanity guard:** book price outside `[1/1.5, 1.5]` of naive multiply → quarantined.
- **Cross-game:** two-game combo fair = product of per-game consensus fairs;
  any game failing consensus → skip.
- **Fail-safe:** untypeable leg, over-cap legs, too-few books → `None`/skip.
- **Regression:** the existing 2-leg spread×total and ml×total RFQs still price
  through the cached fast path with byte-identical fairs (no behavior change for
  today's flow).
- **Dry-run + live-small:** dry-run validates routing/pricing with no exchange
  writes; live-small measures fill-vs-fair across the new shapes.

## 11. Version control / worktree / docs

- **Branch / worktree:** `worktree-kalshi-mm-generalized-combos` (already created).
  Implement + test in the worktree; merge to `main` only on explicit approval;
  clean up worktree + branch after merge.
- **DuckDB safety:** never symlink DBs into the worktree; copy if needed, or test
  from `main` after merge (per repo rules). The new table is created idempotently
  via `CREATE TABLE IF NOT EXISTS` on startup.
- **Docs updated in the merge commit:**
  - `kalshi_mlb_mm/README.md` — pricing section rewritten for the router + on-demand path.
  - `mlb_sgp/README.md` — note the generic-leg-list pricing path + DK/FD wrapper change.
  - Memory `kalshi_mlb_mm_coverage_expansion` — append the generalization.
- **Pre-merge review:** executive-engineer review of the full diff per repo
  checklist (data integrity, resource safety, edge cases, dead code, log hygiene,
  security) before requesting merge approval.

## 12. Open questions for review

1. `MIN_AGREEING_BOOKS=2` confirmed despite DK≡Novig? (accepted-risk path assumed)
2. `MAX_SAME_GAME_LEGS=4` an acceptable v1 coverage cap?
3. `LEGSET_PRICE_DEADLINE_SEC` / `LEGSET_CACHE_TTL_SEC` starting values — propose
   `3s` deadline, `15s` cache TTL; refine once the RFQ window is measured.
