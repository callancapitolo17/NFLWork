# Kalshi MLB MM (Maker) — Full RFQ System Audit

**Date:** 2026-07-19
**Scope:** `kalshi_mlb_mm/` + shared `kalshi_common/` pricing path (the RFQ maker), benchmarked
against combo-maker practice on Polymarket / prediction-exchange MM literature.
**Branch:** `claude/mlb-rfq-audit-8t6fq4`
**Deliverable:** this audit + GitHub epic with tickets (see §6).

---

## 1. Executive summary

The maker is a well-instrumented, fail-safe-first v1 whose stated purpose is to *measure*
whether a fixed quoted margin survives adverse selection. The audit finds:

1. **The measurement loop is broken.** `fills.realized_pnl` is written as `NULL` and no
   settlement sweep exists — the spec's `settlements` table (§7 of the design spec) was never
   built. The single question v1 exists to answer ("does 3% survive adverse selection?")
   cannot currently be answered from the data the bot writes. There is also no post-fill
   markout tracking, so adverse selection can't even be estimated pre-settlement.
2. **Four correctness bugs** in the live path, two of them exchange-facing:
   a cross-game quote can rest past the *earlier* game's first pitch (risk sweep checks only
   the primary game); the quote-replace flow can orphan a live exchange quote when a
   re-submit fails; last-look re-prices from the same stale in-memory cache the quote was
   priced from (README claims a fresh pull — the code does not do one); accept-time fill
   size parsing defaults to 1 contract for dollar-denominated RFQs.
3. **Structural pricing leaks**: constant-ROI margin gives near-zero absolute cushion on
   longshot combos (≈0.3¢ at p=0.10) exactly where combo fair error is largest; the
   consensus band is absolute (±2¢) regardless of price level; DK≡Novig can form a false
   2-book consensus; cross-game products compound per-game fair error under a constant
   margin; and the bot never uses Kalshi's own real-time single-leg markets as an anchor or
   steam detector, despite them being faster than the 150–165s book scrape cycle.
4. **Speed is the binding constraint**, consistent with the prior Polymarket MM research
   (issues #1–#7): quotes rest on data up to ~3 minutes old with a 60s cooldown that expires
   *before* the data refreshes, so the same stale price can be re-hit.

None of this says the strategy is dead — 0 fills to date means nothing has been lost — but
until the measurement loop is fixed and the correctness bugs are closed, scaling or further
margin-tightening would be flying blind.

---

## 2. System as-built (reference)

- **Flow:** poll open RFQs (2s) → scope filter (Phase-1 shapes: 2-leg same-game
  spread×total / ml×total grids, cross-game products of those + marginalized singles) →
  book-consensus fair (median of agreeing books, ±2¢ band, ≥2 books) → two-sided quote at
  `p/(1+3%) − maker_fee`, floored to the $0.001 grid → confirm loop with last-look →
  reconcile sweep vs `/portfolio/positions`.
- **Data:** own SGP scrape of 6 books every 60s (a full cycle is ~150–165s end-to-end;
  staleness window widened to 180s to compensate), into
  `kalshi_mlb_mm_market.duckdb::mlb_sgp_odds`.
- **Defenses:** margin → consensus gate → freshness auto-pull → book-move circuit breaker
  (3¢, per-tick + per-quote drift) → tipoff blackout (5 min) → last-look (2¢ drift) →
  reconciliation; plus void-rate halt, per-creator halt (inert — Kalshi anonymizes
  creator ids), per-combo cap/cooldown, per-game/daily/per-fill dollar caps.
- **Status:** live-small measurement phase; 13 quotes in ~8 days at 5% → 0 accepts →
  tightened to 3% (2026-06-18) and `MIN_AGREEING_BOOKS` 3→2 to raise quotable surface
  (~5 → ~40-47 tuples).

---

## 3. Findings — correctness bugs (fix before anything else)

### B1. Cross-game quotes survive the earlier game's first pitch — `main.py:924-927`
`_risk_sweep_tick` re-checks tipoff using only `live_quotes.game_id` (the *primary* game,
which is just the first key of a dict partition). Discovery gates on the earliest commence
time across all games at quote time, but a resting quote is only ever swept against the
primary game's clock. A spread×total(game A) × ml(game B) quote where B starts first rests
into B's live window — in-play pickoff on a stale pre-game fair is the single worst fill
available to a counterparty. **Fix:** persist all game_ids per quote (reuse the
`fill_games` pattern, e.g. a `quote_games` map or re-derive from `seen_rfqs.legs_json`) and
sweep against the earliest commence time.

### B2. Quote-replace can orphan a live exchange quote — `main.py:639-649`
On a beyond-hysteresis re-price the old row is marked `status='replaced'` *before*
`submit_quote` is called. If the submit fails (`qid None` — 4xx/5xx/timeout), Kalshi never
auto-cancels the prior quote (that only happens on successful resubmission), but the DB now
says `replaced`, so the confirm loop stops polling it and the risk sweep ignores it. Result:
an untracked live quote resting on the exchange at a price we can no longer manage. An
accept on it would silently void (30s window) — "safe" but it burns void-rate budget and is
invisible. **Fix:** submit first, mark `replaced` only on success; on failure, explicitly
`cancel()` the old quote and mark it accordingly.

### B3. Last-look re-prices from the same stale cache — `main.py:707-708`
`_confirm_tick` computes `cur_fair` from the in-memory `_SGP_ODDS` — the exact snapshot the
quote was priced from (books refresh every ~150–165s; the confirm loop runs every 2s). For
any accept that arrives within the same scrape cycle, `cur_fair ≈ prev_fair` by
construction and the drift check is a no-op precisely in the fast-pickoff window it exists
to guard. The README's claim ("recomputes fair from a fresh book pull") does not match the
code. **Fix:** on accept, trigger a targeted on-demand re-price (single-combo scrape via
`SGPService`, and/or a Kalshi single-leg midpoint sanity check — see P-5) before
confirming; keep the hard-fail semantics.

### B4. Accept-time fill size defaults to 1 contract — `main.py:720`
`contracts = int(float((q or {}).get("contracts", 1) or 1))`. RFQs are ~85%
dollar-denominated and Kalshi's live payloads use fixed-point strings (`contracts_fp`) —
the quote-status response very likely does too. A 100-contract fill recorded as 1 contract
under-counts every cap for up to 30s until the reconcile sweep corrects it (N8's
conservative treatment covers the *caps* partially, but the `fills` row itself and the
notify are wrong). **Fix:** parse `contracts_fp` first (mirror `_rfq_requested_contracts`),
and treat an unparseable size as reconcile-blocking, not 1.

### B5. Transient market-fetch failure permanently blacklists a combo — `main.py:459-465`
If `get_market()` returns None (any HTTP blip), `legs=None` → `in_scope=False` is cached in
`_SCOPE_CACHE` forever (no TTL, no retry, unbounded growth). Every future RFQ on that
ticker is silently skipped until restart. **Fix:** only cache negative verdicts when the
market was actually fetched and decoded; add a TTL / size bound to the cache.

---

## 4. Findings — structural leaks (pricing & adverse selection)

### P-1. Constant-ROI margin ≈ no cushion on longshots
Cushion in probability points is `p·r/(1+r)` ≈ 0.029·p. At fair p=0.50 that's ~1.5¢/side;
at p=0.10 it's ~0.3¢; at `MIN_FAIR_PROB`=0.05 it's ~0.15¢. Combo fair error (book
dispersion + devig error + correlation error) is *not* proportional to p — it's roughly
constant-to-growing in absolute terms at the longshot end. So the bot's thinnest absolute
protection sits exactly where its fair is least reliable, and cross-category MVE flow skews
longshot. A patient counterparty only needs to hit the low-p quotes. **Fix:** floor the
margin in probability points (e.g. `margin = max(ROI-margin, k·σ_books + c)`), or price the
spread in z-space (probit) so it scales with distributional uncertainty; at minimum raise
`MIN_FAIR_PROB` until the floor exists.

### P-2. Consensus band is absolute (±2¢) at every price level
±2¢ at p=0.50 is a tight 4% relative test; ±2¢ at p=0.08 tolerates a 25% relative
disagreement — and then the *median* of books that disagree by 25% becomes the fair we
quote 0.2¢ around. Same class of bug as P-1: the gate weakens exactly where quoting is most
dangerous. **Fix:** make the band relative (e.g. ±k% of median) or defined in z/logit
space; scale `BOOK_CONSENSUS_BAND` and the quoted margin off the *measured cross-book σ*
per combo instead of constants.

### P-3. DK≡Novig false consensus at MIN_AGREEING_BOOKS=2
Documented as the "top pricing-integrity exposure" in the generalized-combos spec and still
open: Novig mirrors DK, so a DK+Novig pair passes the 2-book gate as ~one independent
source. **Fix:** independence-aware consensus — book "source groups" (DK+Novig = 1), gate
on distinct groups; log group composition per quote so fills can be attributed to weak- vs
strong-consensus quotes.

### P-4. Cross-game products compound error under a constant margin
`combo_fair = Π game_fair` with independent per-game errors: relative error grows ~√N (or
linearly if biases correlate — same books, same devig method, same direction). The 3% ROI
margin does not grow with N. A 3-game combo priced from 2-book consensuses on each game
carries far more fair uncertainty than a single-game grid but gets the same cushion.
**Fix:** per-game-count margin scaling (e.g. `(1+r)^N` pricing or an additive per-game
z-buffer), and/or require stronger consensus (more books, tighter band) for multi-game
combos.

### P-5. Kalshi's own single-leg markets are unused
The v1.1 correlation-premium gate (design spec §13) remains deferred, but the bigger prize
is latency: KXMLBSPREAD/KXMLBTOTAL/KXMLBGAME orderbooks move in real time, while our books
refresh every ~2.5 min. Kalshi singles give us, for free: (a) a real-time marginal anchor
to sanity-check every combo fair (`combo_fair` must respect Fréchet bounds and stay near
`corr_premium × Π marginals`); (b) a sub-second steam detector — if a constituent single's
mid jumps while our books sleep, pull the combo quotes *now* (this is exactly the
"mid-jump circuit breaker" pattern from issue #1, applied to constituents); (c) the
denominator for the correlation-premium gate. **Fix:** subscribe (WS or fast REST) to the
constituent singles of every combo we're resting quotes on; wire a constituent-move
circuit breaker + Fréchet/premium sanity gate.

### P-6. Symmetric two-sided quoting into (likely) asymmetric flow
Both sides always get the same margin, but MVE RFQ flow is plausibly bimodal: retail parlay
buyers (buy YES on longshots — flow we *want* to fill) and sharps picking off stale quotes
(flow we don't). The two sides of our quote face completely different toxicity, and we
can't currently even see which side fills. **Fix:** first measure (fill side × p-level ×
markout, from the research firehose once fills exist); then skew per-side margin by
side-toxicity (the reservation-price / inventory-skew pattern from issue #3 applies once
positions accumulate).

### P-7. Re-pickoff window: cooldown (60s) < data refresh (~150–165s)
`COMBO_COOLDOWN_SEC=60` expires before the books that priced the previous fill have
refreshed even once. A counterparty who just picked us off can post a new RFQ on the same
combo and hit the *same stale price* again. **Fix:** cooldown until *fresh post-fill data*
exists for the combo (cooled_until ≥ next successful scrape covering it), not a fixed 60s;
or block re-quoting a combo at a price within ε of the just-filled price until fair moves.

### P-8. No margin/competitiveness experimentation framework
5%→3% was a manual global change judged on 13 quotes / 0 fills. Competitor quotes are
invisible (403), so the *only* estimator of the demand curve is our own fill rate vs
margin. A static global margin learns nothing. **Fix:** per-combo-family margin ladder with
deliberate randomization (e.g. margin drawn from {2.5%, 3%, 4%} per quote, logged),
so fill-rate-vs-margin and markout-vs-margin curves accumulate even at low volume.

---

## 5. Findings — measurement, risk controls, ops

### M-1. Settlement sweep missing → `realized_pnl` never populated (spec §7 gap)
The primary v1 deliverable is uncomputable. **Fix:** settlement sweep — poll
`/markets/{ticker}` for `settled`/result on tickers with open fills, write `realized_pnl`
(and a `settlements` audit table), plus a P&L attribution split: margin captured vs fair
drift (quote→confirm) vs settlement vs fee.

### M-2. No post-fill markout tracking
Binary settlement P&L converges slowly (huge variance per fill); consensus-fair markout at
+60s/+5m/+30m after each fill estimates adverse selection with far less data — this is the
same markout-EWMA pattern poly-maker uses (issue #2) and what the `fills` table's
`fair_at_confirm` column starts but doesn't finish. **Fix:** post-fill markout job
(re-price the combo at fixed horizons, store in a `fill_markouts` table), feeding both the
validation report and (later) a toxicity brake.

### M-3. Funnel telemetry exists but has no aggregated readout
`research_queries.sql` is a good start (void breakdown, book participation, latency), but
there's no periodic "state of the maker" report: quotable-universe size, skip reasons,
quote rate, margin distribution, staleness-age at quote, and (once fills happen) fill rate
and markout. The `kalshi_mlb_monitor` dashboard shows operational state, not strategy
learning. **Fix:** scheduled report (extend monitor or a daily script) so the
measurement phase actually produces a readable answer.

### R-1. Open-quote exposure not counted in per-game / daily caps
N7 counts in-flight quotes only in the *per-combo* cap. 25 open quotes across different
combos of the same game can, in a burst (one sweep by one counterparty), fill to
25 × $50 = $1,250 against a $50 per-game / $375 daily cap. Worst-case open exposure should
gate discovery for *all* caps, not just per-combo.

### R-2. No exit / hedge path
Every fill rides to settlement. When post-fill info arrives (constituent steam, scratch),
the bot can neither trade out of the MVE position nor hedge delta with Kalshi singles.
Acceptable at $50/fill; unacceptable at scale. **Fix (later phase):** hedging module —
compute per-leg deltas of held combos, offset the dominant leg in the singles orderbook
when markout breaches a threshold.

### O-1. Hot-loop DuckDB connect-per-call
`_commence_time`, `_resolve_game_for_legs`, and several per-RFQ lookups open a fresh DuckDB
connection each call, per RFQ, per 2s tick — lock churn against the scraper writes and the
monitor's readers. Cache game resolution (it's static per day) and commence times in
memory with TTL.

### O-2. Timestamp discipline
`fills.filled_at` etc. are naive `TIMESTAMP` compared against naive-local `now()` in the
reconcile sweep, with tz-normalization comments papering over it — this repo's own rule is
TIMESTAMPTZ UTC everywhere. Migrate the maker schema before the data matters.

### O-3. Still-unverified exchange semantics
Maker fee = 25% of taker (assumed); quote-status response shape / `accepted_side`
semantics (side-held inference); RFQ TTL. All flagged in the README; all still open —
the accept-observed firehose event will answer them on the first real accept, but the
fee assumption also silently shapes every quoted price today.

### O-4. Phase 2 (on-demand novel shapes) unbuilt — coverage unknown
Scope is still the two grid families (+ cross-game products). The share of RFQ flow
skipped as `out_of_scope`/`on_demand` is loggable today but not reported (see M-3); if it
dominates, Phase 2 (spec'd 2026-06-28: sided cross-product devig, `price_legset`,
leg-set cache) is the volume unlock.

---

## 6. Industry research — how the best combo makers do it

*(Web research pass 2026-07-19 + the 2026-07-18 Polymarket MM research already encoded in
issues #1–#7. Sources are search-extracted where docs are proxy-blocked; the biggest
non-public gaps — Kalshi DMM obligations, Polymarket combo maker penalty terms, any firm's
internal combo pricing — are flagged as such.)*

### 6.1 Polymarket Combos (the closest comparable)
Polymarket launched sports parlays ("Combos", CFTC-self-certified CAOCs, ~May 2026) as an
**RFQ auction, not a CLOB**: user request → connected makers get `RFQ_REQUEST` over
authenticated WebSocket → **makers have 400ms to quote** → best quote wins → user has 10s
to accept → makers with **Last Look enabled get a 1s confirm/decline window** for final
risk checks. Nobody rests combo liquidity anywhere; competitive combo making is a
sub-second quote SLA plus a last-look risk check. Their CLOB LP program (not combos) pays
daily rewards via a quadratic proximity-to-mid score with two-sided boost; makers pay zero
fees and receive 25% (15% sports) of taker fees as rebates.

### 6.2 Kalshi RFQ mechanics that matter to us
- Quotes are implicitly **full-RFQ-size**; either side may be 0 (**one-sided quoting is
  allowed** — decline a side without skipping the RFQ); rejected if yes_bid+no_bid > $1.
- Two-step execution: requester accepts → **maker confirms within 30s — but only 2s in
  High Volatility Markets**. The confirm step is effectively free last-look; no published
  penalties for declining (DMM obligations are private bilateral agreements). Our confirm
  path must complete in well under 2s to be HVM-safe.
- Maker fee = 25% of taker fee (`0.07·C·(1−C)` base); combos charged one fee on the
  combined contract, not per leg. Multiplier for parlays needs re-verification against the
  7.7.26 fee schedule.
- SIG (first institutional MM on Kalshi, owns Nellie Analytics) and Jump Trading both make
  prediction markets with options-style statistical MM; nothing about their combo pricing
  is public.

### 6.3 Literature → our situation
- **Glosten–Milgrom:** quote around the fair *conditional on being hit*, not the
  unconditional fair; every accept is information — shift fair toward the hit side and
  widen that side.
- **Avellaneda–Stoikov:** reservation price `r = fair − q·γ·σ²(T−t)` (inventory skew) and
  spread `∝ γσ²(T−t)` — margin should scale with belief volatility and time-to-event, not
  sit at a flat 3%. Probability-price volatility is state-dependent (explodes near event
  start / news, compresses in quiet hours), so flat margin under-prices the dangerous
  windows and over-prices the safe ones.
- **Last-look literature (FX):** momentum-conditioned rejection (did the price move
  between quote and execution?) is the single best discriminator of toxic flow — which is
  exactly what our B3-bugged last-look fails to implement.
- **SGP correlation pricing:** industry standard is empirical joints where dense + Gaussian
  copula over devigged marginals elsewhere; same-game correlations commonly 0.3–0.5+; SGP
  hold is 3–5× straight hold, so book SGP prices embed large correlation-inclusive vig —
  consensus-devig is sound but devig-method error concentrates at longshot prices.
- **Sportsbook/exchange practice:** books *suspend* markets on pitcher scratches and
  lineup news rather than trying to out-price them; in-play exchange makers get farmed by
  faster event data. A maker is stale whenever its slowest input is.

### 6.4 Ranked actionable techniques (top of the research brief's 18)
1. Real accept-time last-look re-price (fresh inputs + Kalshi singles), momentum-keyed.
2. Staleness-aware quote fading (widen/decline as input age grows; per-book age tracking).
3. Event-window halts: lineup posting, probable-pitcher changes, pre-tipoff.
4. Volatility/time-scaled margin replacing flat ROI.
5. Inventory-skewed reservation price (fill_games ledger already provides q).
6. Post-fill fair updating + hit-side widening (Glosten–Milgrom ratchet).
7. Markout-based toxicity scoring by shape/size/time bucket, fed back into quoting.
8. Longshot devig validation + tail shading below ~10¢.
9. Copula upgrade path for correlation (grid > copula > Fréchet).
10. Cross-book combo-price dispersion as a per-quote uncertainty input (also softens the
    DK≡Novig problem more gracefully than a count gate).
11. One-sided declines instead of full skips.
12. Size-conditional margin (large RFQs are more informed).
13. Quiet-period tightening to win auctions when inputs are fresh.
14. Kalshi single-leg jump detection ("we are the stale ones" detector).

---

## 7. Epic + tickets

Tracked on GitHub: epic issue **kalshi_mlb_mm: RFQ maker improvement epic** with one issue
per ticket, linked as sub-issues, sequenced in four milestones:

**M0 — Trustworthy measurement** *(you cannot iterate on what you cannot measure)*
- T1. Settlement sweep: populate `realized_pnl` + `settlements` table + P&L attribution (M-1)
- T2. Post-fill markout tracking (+60s/+5m/+30m) + toxicity readout by shape/size/time (M-2)
- T3. Daily "state of the maker" report: funnel, skip-reason & out-of-scope share,
  staleness-age at quote, margin vs fill-rate (M-3, O-4 measurement)

**M1 — Correctness** *(close exchange-facing bugs before more volume)*
- T4. Cross-game tipoff sweep bug (B1)
- T5. Quote-replace orphan on failed resubmit (B2)
- T6. Real last-look: fresh targeted re-price + Kalshi-singles check inside the confirm
  window; HVM-safe (<2s) path (B3)
- T7. Correctness batch: accept-time `contracts_fp` parsing + scope-cache transient
  blacklist + hot-loop DB connection caching + TIMESTAMPTZ migration (B4, B5, O-1, O-2)

**M2 — Pricing integrity** *(make the quoted margin mean what it says)*
- T8. Uncertainty-scaled margin: probability-point floor + dispersion & time-to-game
  scaling + per-game-count scaling for cross-game combos (P-1, P-4, techniques 4/10/12)
- T9. Consensus overhaul: relative/z-space band + DK≡Novig independence groups (P-2, P-3)
- T10. Cooldown ≥ data refresh + same-price re-quote block + re-request pattern detection
  (P-7)
- T11. Open-quote worst-case exposure counted in per-game and daily caps (R-1)

**M3 — Speed & competitiveness** *(turn a safe bot into a competitive one)*
- T12. Kalshi singles anchor: constituent-jump circuit breaker + correlation-premium /
  Fréchet sanity gate (P-5, spec §13 v1.1)
- T13. MLB event-window halts: lineup posting + probable-pitcher-change detection from
  schedule diffs (research technique 3)
- T14. One-sided quoting + per-side toxicity asymmetry + post-fill fair nudge (P-6)
- T15. Margin experimentation ladder: randomized per-quote margin, demand-curve estimation
  (P-8)
- T16. Phase 2 on-demand pricing for novel shapes — gated on T3 showing material
  out-of-scope share (O-4)
- T17. Inventory skew + delta hedging via Kalshi singles — later phase, needs fills (R-2)

Dependencies: T2/T3 unblock T8/T14/T15 calibration; T3 gates T16; T12 feeds T6's singles
check; T17 waits for realized fill data from M0.
