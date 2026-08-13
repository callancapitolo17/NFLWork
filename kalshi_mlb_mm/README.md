# Kalshi MLB MM (Maker) Bot

Independent maker daemon that listens for others' RFQs on the Kalshi cross-category MVE collection, prices MLB combos against a book-consensus fair value, and provides two-sided quotes at a fixed 5% ROI margin. Coexists with the taker (`kalshi_mlb_rfq/`) as a separate OS process with no runtime dependency on it.

**Combo scope (Phase 1).** Pricing flows through `kalshi_common.legset` + `kalshi_mlb_mm.router`, which generalize the original fixed 2-leg path:

- **2-leg same-game grids** — spread×total and moneyline×total (both FG), priced from each book's 4-cell devig grid (unchanged math).
- **Cross-game combos** — each game's sub-combo is priced independently and the per-game fairs are **multiplied** (independence assumption); a single leg within a game is marginalized out of that game's grid.
- **On-demand same-game shapes (Phase 2)** — 3-leg, spread+ml, total+total and any other novel shape route to `on_demand` and are priced by live book queries at RFQ time (see **On-demand pricing** below). Lone single-leg RFQs remain out of scope (skipped fail-safe).
- **F5 (first 5 innings) legs — issues #84/#85, epic #82.** `CanonicalLeg` carries a `period` field (`"FG"` default / `"F5"`); `KXMLBF5SPREAD` / `KXMLBF5TOTAL` legs parse with the same suffix grammar and sign semantics as FG (live-verified 2026-08-11: total suffix `-7` ⇔ floor_strike 6.5; spread `-LAD3` ⇔ "LAD −2.5 first 5"). Any F5-containing multi-leg set routes `on_demand` — the 2-leg grids stay **FG-only** (the grid tables are full-game surfaces), and `period` joins the leg-set hash + duplicate-market guard key (FG −1.5 with F5 −1.5 is not a contradiction). Since #85 the six on-demand book hooks are **period-aware**: each `build_structure` returns a normalized `{"FG": bucket, "F5": bucket}` structure (DK/FD re-keyed from their `fg`/`f5` fetch output, MGM/CZR pass `parse_markets` through, Novig builds both periods from one raw tree via `build_line_structure(period=...)` — its `SPREAD_1H`/`TOTAL_1H`/`MONEY_1H` types ARE first-5, confirmed 2026-08-12 by line scale vs FD's FG total), each `resolve_legs` picks the bucket per leg from `leg.period` (missing bucket → clean whole-book decline), and ProphetX dispatches on a `(market_type, period)` name-alias map (`_OD_MARKET_NAMES`; PX listed **zero** F5 markets on the 2026-08-12 recon board, so PX cleanly declines F5 until it re-lists them). Per-book F5 status: DK ✅ (subcategory 15628), FD ✅, MGM ✅, CZR ✅, Novig ✅ (`_1H` types), PX ⚠️ name-mapped but board-absent. `SGPService.price_on_demand` still fails closed on any **unknown period** (`ON_DEMAND_PERIODS`). The firehose `on_demand_result` payload and `quote_priced.live_games` trace carry a `periods` list since #85.
- **F5-winner (`KXMLBF5`) pricing via ±0.5 run-line re-encoding — issue #86.** Kalshi's F5 winner family is 3-way ({TEAM_A, TEAM_B, TIE} tickers per game) and the team markets are **unconditional — ties lose** ("all team strikes will resolve to No", live rules 2026-08-12). Books' F5 *moneyline* is the wrong row to read: it's conditional push-2-way (recon implied sums 1.01–1.07 across DK/FD/MGM/NV; a tie refunds), so quoting it unconverted would overprice Kalshi YES by ~P(tie) ≈ 12–15% — a craftable pick-off. The books' F5 *run line at ±0.5* is the **right** row: "team wins F5" ⟺ "team −0.5 covers" (win by ≥1; tie loses) and a team market's NO side ⟺ "other team +0.5" (win **or tie** — the tie mass a NO holder wins, which a naive ml-complement encoding would silently drop). Books price that market two-sided with **no push**, so `legset.parse_leg` simply **re-encodes `KXMLBF5-{TEAM}` legs (both sides) as F5 spread legs at ±0.5** and the entire #85 pipeline — resolvers, partition devig, consensus, confirm — prices them with zero tie modeling, zero conversion math, zero new book surface. **TIE-side legs stay unpriceable** (the one leg ±0.5 can't express; no book prices it — declined `out_of_scope_f5_tie_leg`). A **contradiction guard** in `classify_subcombo` backstops the encoding: opposed spread/total legs whose margin/total bounds are jointly empty (e.g. home −0.5 + away −0.5 — both F5 winners — or Over 8.5 + Under 4.5, per period) are `unpriceable`; pre-#86 the home+away winner pair was caught by the duplicate-market guard's `("ml","F5",None)` collision, and the new guard also closes the pre-existing FG hole (home −1.5 + away −1.5 slipped past the same-key dup check). A book missing the ±0.5 line for a game simply declines that flight (fail-closed; the F5 run-line main is typically ±0.5 — the #85 live smoke priced F5 spread −0.5 at 4 books).
- **RFI (`KXMLBRFI`) pricing via 1st-inning-total re-encoding — issue #87.** `KXMLBRFI` is ONE binary market per event (market ticker = event ticker, no team/line suffix) whose rules read "if either team scores a run in the first inning" with `floor_strike 1 / greater_or_equal` — YES is **exactly** "1st-inning total runs ≥ 1" (live-verified 2026-08-13). So `legset.parse_leg` re-encodes RFI legs as **total legs at 0.5 with `period="I1"`** (yes→over, no→under) — the #86 pattern: books price 1st-inning 0.5 totals two-sided with no push, and the whole #85 period-aware pipeline prices them with zero new pricing machinery. `ON_DEMAND_PERIODS` admits `"I1"`; the dup guard kills a YES+NO pair per game (`("total","I1",0.5)` collision); I1 bounds never interact with FG/F5 totals (different variables). Per-book I1 status (live recon 2026-08-13, all combinability tests = a real 2-leg SGP price of {I1 O0.5 + FG ML}): **DK ✅** "Runs - 1st Inning" (name-seeded into an `i1` bucket from the parlays payload — the ≤2.0-line range heuristic deliberately skips inning totals; priced 3.16), **FD ✅** "1st Inning 0.5 Runs" (runners are bare Over/Under with `handicap 0`, so the line is FIXED at 0.5 from the market name via `_FIXED_LINE_MARKETS`; priced 3.35), **Novig ✅** `FIRST_INNING_TOTAL` type (strike 0.5; priced 3.33), **MGM ✅** "Will there be a run in the 1st inning?" Yes/No mapped Yes→over / No→under at 0.5 (its "1st Inning Runs" sibling is 3-way 0/1/2+ and is never read; priced 4.5), **PX ⚠️** "1st Inning Total Runs" name-mapped but the board carried outcome definitions with no priced lines — declines cleanly until PX carries I1 liquidity, **CZR ⏸** unmapped (AWS-WAF auth outage blocked name recon 2026-08-13; no `I1` bucket emitted → clean whole-book decline; follow-up when auth recovers). Live smoke (`RFI_LIVE_SMOKE=1 … test_rfi_live_smoke.py`): 3 books priced {RFI yes + home ML} with σ_z 0.043 (< 0.07 gate).

Legs are grouped into games by `legset.game_id_of` — the **game code** inside
the event ticker (`KXMLBTOTAL-26JUN102105MILATH` → `26JUN102105MILATH`), not
the ticker itself. Kalshi gives every market family its own event, so keying on
the whole ticker split a game's spread and total legs into two partitions and
made the same-game grids unreachable (**issue #71**, fixed 2026-08-04). The
symptom was silent: same-game combos were priced as two independent games and
multiplied, and the routing census read as "no same-game demand". Re-running
that census over the 738 in-scope RFQs recorded to date, post-fix:

| route | sub-combos |
|---|---|
| `single` | 1611 |
| `grid_spread_total` | 50 |
| `grid_ml_total` | 36 |
| `on_demand` | 1 |

Spread lines are **signed home-perspective throughout** (**issue #70**, fixed
2026-08-04): a home-team margin ticker (`…-COL2` with COL home) canonicalizes
to `spread_line = -(N-0.5)`, an away-team margin ticker (`…-MIA2` with MIA
away) to `+(N-0.5)` — matching `mlb_sgp_odds.spread_line`, where the
`Home Spread` cell is always the home team at the signed line and `Away
Spread` its negation. Before the fix both teams' tickers collapsed onto
`-(N-0.5)`, so an away-favourite leg (away −1.5, a hard event) was priced
with the away **+1.5** cell (an easy event) — fair systematically too high on
58% of in-scope spread legs. The maker's scrape cycle now passes
`both_teams=True` (as the taker always did), scraping **both** teams' margin
grids per game — roughly double the spread target lines per cycle; without
the positive-line grids an away-favourite leg would decline `too_few_books`.
Research-firehose rows written before 2026-08-04 encode away favourites with
a negative `spread_line`; rows after are signed (convention epoch — mind the
boundary in longitudinal queries).

The read side is **book-agnostic**: on-demand flights price whichever of the 6 books' scraper modules answer (the same per-book resolve/price hooks the scrapers use), so a new book's moneyline rows are priced with no maker-side change.

**Spec:** `docs/superpowers/specs/2026-05-26-kalshi-mlb-mm-design.md`
**Plan:** `docs/superpowers/plans/2026-05-27-kalshi-mlb-mm-maker-bot.md`

## Quick start

```bash
# Run from the main repo root (NOT from a worktree — same requirement as the taker)
cd ~/NFLWork

# 1. Create a venv and install deps
python3 -m venv kalshi_mlb_mm/venv
./kalshi_mlb_mm/venv/bin/pip install -r kalshi_mlb_mm/requirements.txt

# 2. Copy env template and fill credentials
cp kalshi_mlb_mm/.env.example kalshi_mlb_mm/.env
# Edit kalshi_mlb_mm/.env: set KALSHI_API_KEY_ID, KALSHI_PRIVATE_KEY_PATH, KALSHI_USER_ID

# 3. Dry-run first (computes + logs quotes, makes NO exchange writes)
./kalshi_mlb_mm/venv/bin/python -m kalshi_mlb_mm.main --dry-run
# Inspect quote_decisions in kalshi_mlb_mm/kalshi_mlb_mm.duckdb to verify pricing looks sane

# 4. Live mode (backgrounded)
./kalshi_mlb_mm/venv/bin/python -u -m kalshi_mlb_mm.main >> kalshi_mlb_mm/bot.log 2>&1 &
tail -f kalshi_mlb_mm/bot.log
```

The bot must run from the main repo root (`~/NFLWork`), not from inside the package directory or a worktree. The same restriction applies to the taker — see the taker README and `kalshi_mlb_rfq/README.md` for the reason (`.env` loading and the `kalshi_draft/auth.py` path resolver both depend on the working directory).

The dry-run is a wiring smoke-test: it decodes open RFQs, detects in-scope combos, computes prices, and logs everything to `quote_decisions` — but makes no `POST /communications/quotes` calls. It is not the validation mechanism; live-small is.

## Stopping

```bash
# Graceful SIGTERM — cancels all open quotes before exiting
kill $(pgrep -f "kalshi_mlb_mm.main")

# Kill switch — stops quoting but keeps the process alive
touch kalshi_mlb_mm/.kill
# Resume quoting:
rm kalshi_mlb_mm/.kill
```

## Architecture

REST-polling daemon, single process. Eight timed sub-loops:

| Loop | Cadence | Job |
|---|---|---|
| Discovery + quote | 2s | Poll open RFQs → scope filter → price → submit or refresh quote |
| Confirm | 2s | Poll open quotes → on `accepted`, last-look gate → confirm or void; an explicit 404 on the quote GET (RFQ expired → Kalshi deleted the quote) closes the row as `expired` — transient errors leave it open |
| Risk sweep | 10s | Kill-switch, books-unhealthy auto-pull (#38 health-dark), tipoff cancel, constituent-jump + drift-since-quote cancel |
| Reconcile sweep | 30s | Verify recorded fill side/size against Kalshi `/portfolio/positions` (live only) |
| Settlement sweep | 600s | Poll `GET /markets/{ticker}` for combos with unsettled fills → write `fills.realized_pnl` + `settlements` audit row (live only; issue #12) |
| Structure warming | 120s | #50: structure-only pass (no pricing calls) keeping every book's events/structure TTL caches + the CZR WAF token hot, so live fetches never pay a cold start |
| Target-line refresh | 300s | #81: Kalshi MVE enumeration + Odds API schedule → `mlb_target_lines` (game resolution, tipoff gating and warming read it). **Zero book requests** |
| Coverage summary | 300s | #81: drain the service's per-book on-demand outcome tally into an `on_demand_coverage` research event — the record of which books actually answer live fetches |

There is **no background SGP sweep** (#81 deleted #57's demoted remnant): the
maker's only book traffic is on-demand flights triggered by live RFQs plus
the structure-warming pass. `sgp_fetch_health` proves it — post-#81 rows
carry only the `on_demand` and `warming` paths (`book_requests_per_day` in
`kalshi_common/fetch_health_queries.sql` is the before/after readout).

**Transport seam.** All exchange I/O goes through two interfaces: `RFQSource.poll()` / `RFQSource.get_market()` and `QuoteGateway.submit_quote()` / `.confirm()` / `.cancel()` (v1: `RestQuoteGateway`). RFQ discovery is **WS-only** (issue #56, user decision 2026-08-07 — no mode switch, mirroring #54's live-only decision; rollback is a git revert): `WebSocketRFQSource` (`ws_rfq_source.py`) runs a daemon reader thread subscribed to Kalshi's `communications` WS channel (`rfq_created` / `rfq_deleted`) and mirrors it into a stateful open-RFQ set, so `poll()` keeps the **level-triggered** contract the discovery tick depends on (all open RFQs every tick — that re-entry drives on-demand re-fetch, quorum refinement, and pulls). Payloads are stored verbatim (`contracts_fp` stays a string), so downstream is byte-compatible with the old REST poll. `RestRFQSource` (the 2s poll of `GET /communications/rfqs?status=open`) is not a mode: it survives only as the WS source's automatic fallback and gap-fill path, and `get_market()` stays pure REST. **Ingestion filter + mirror TTL (2026-08-10):** the `communications` channel is the entire exchange's RFQ firehose (~5–10k `rfq_created`/min at Monday-evening peak) and Kalshi sends no `rfq_deleted` on expiry, so an unfiltered mirror leaks without bound (48k entries in 5 min observed). The reader drops any RFQ whose legs aren't all `KXMLB*` at the door (no mirror entry, no research event, no downstream decision row; legless frames pass fail-open to the tick's budgeted `get_market` fallback), logs a count every 25k drops, and `poll()` evicts entries older than `MAX_RFQ_AGE_SEC` on an arrival clock (`mirror_ttl_sec` ctor override for tests). Gap-fill snapshots pass through the same filter.

**WS fallback semantics (never silently deaf).** Liveness = any WS frame (server pings every ~10s count) within `WS_HEARTBEAT_TIMEOUT_SEC`. The watchdog runs inside `poll()` on the tick thread — a wedged reader thread still trips it. Trip → ERROR log + one notify + `rfq_source_transition` research event, and `poll()` transparently serves the REST source until the feed recovers, which flips back the same loud way. The reader thread applies the same staleness test itself: a subscribed socket with no frame for a full heartbeat timeout is torn down and reconnected — a TCP connection killed without a close/error frame (network blip, NAT drop; live incident 2026-08-12 left the reader deaf 13h) raises only recv timeouts, which the poll()-side watchdog cannot recover from on its own. Reconnects use backoff+jitter and re-subscribe; readiness is gated on the server's `subscribed` ack (an unacknowledged or rejected subscribe forces a reconnect — a live socket with a dead subscription passes the ping watchdog, the one deafness heartbeats can't see); each (re)connect performs ONE REST gap-fill that reconciles the mirror wholesale — additions AND removals, so a missed `rfq_deleted` can never leave a zombie RFQ being re-fetched forever (events buffered during the gap-fill re-apply on top, so nothing lands lost). Repeated connect failures (`WS_CONNECT_ALERT_STREAK`) alert once per outage. Every WS `rfq_created` emits `rfq_seen_latency` (payload `created_ts` vs receipt) to the firehose, so the discovery-latency win over the 2s poll is measured per RFQ.

**Shared math.** Fair value, EV calc (including `maker_fee_per_contract`), authenticated HTTP, SGP orchestration, and leg-typing helpers all live in `kalshi_common/`. Both bots import from there; the taker's original files are now one-line re-export shims (behavior unchanged).

**SGP pricing (2026-06, rescoped by #81).** The bot prices SGPs **in-process** via `kalshi_common/sgp_service.py::SGPService` — persistent per-book HTTP clients, TTL-cached structure fetches, prices never cached. Since #81 the maker uses the service ONLY for on-demand flights (`price_on_demand`, budgeted by `ON_DEMAND_DEADLINE_SEC`) and structure warming — `service.refresh()`/`sgp_cycle` (the full-slate scrape and its `SGP_SCRAPER_TIMEOUT_SEC` per-book budget) remain in `kalshi_common` for the taker's sweep only.

**State DB.** The bot writes `kalshi_mlb_mm/kalshi_mlb_mm.duckdb` (quotes, fills, positions, decisions, `fill_games`, and `quote_games`). The `fill_games` table maps each fill to **every** game its combo touches (one row per game); the per-game exposure cap reads a fanned-out `fills⋈fill_games` join so a cross-game combo's full stake counts against each of its games, while the daily cap and P&L keep reading `fills` (one row per combo) and are never double-counted. The per-game and daily caps additionally count **open quotes' worst-case exposure** (issue #22): each `live_quotes` row freezes `worst_exposure_usd` at quote time, `quote_games` fans it out per game (same attribution rule as `fill_games`, written in the same transaction as the quote insert), and both cap gates sum `fills + open quotes` so a burst sweep of resting quotes can no longer fill to multiples of the caps before the first fill registers. A `NULL` `worst_exposure_usd` (pre-migration row) is counted at the per-fill cap, mirroring the N8 unreconciled-fill convention; the RFQ currently being re-quoted is excluded from the sums (its resting quote is superseded, not stacked). The sibling `kalshi_mlb_mm/kalshi_mlb_mm_market.duckdb` holds SGP-line and SGP-odds data (same pattern as the taker's `kalshi_mlb_rfq_market.duckdb`). The v1-hardening pass removed the model component of the blend, so there is no longer a read-only dependency on `Answer Keys/mlb_mm.duckdb`. All timestamp columns are `TIMESTAMPTZ` UTC (repo rule); existing naive-local DBs are migrated in place on startup by the idempotent `MIGRATE_SQL` (instants preserved), and the monitor normalizes reads so both frames render identically.

**Quote replace ordering.** When a re-priced quote moves beyond `QUOTE_HYSTERESIS`, the bot submits the NEW quote first and only then marks the old row `replaced` — Kalshi auto-cancels the prior quote only on a successful resubmit, so a failed submit leaves the old row `open` (still confirm-polled and risk-swept) and the next tick retries.

## Pricing

### Live pricing (issue #54) — the only quote path

EVERY in-scope sub-combo — 2-leg grids, each game of a cross-game combo (including its single legs), novel shapes — rides the on-demand engine: the RFQ landing triggers a fetch of exactly that leg set at all six books, and the quote is priced from those fresh per-book fairs through the same #20 consensus gate and #19 margin machinery (the math is input-agnostic; `σ_books` comes from the live set). **There is no cache in the quote path — and since #81 no sweep at all** (deliberately no cached-pricing mode or fallback — user decision 2026-08-05; rollback is a git revert). A fetch that lands with **zero** books declines with `live_fetch_timeout` — never a silent cache fallback — and is not re-fetched until the empty landing ages out (~`QUOTE_FRESH_SEC`), so a dead slate retries on a ~15s cadence instead of hammering per tick. A live consensus with too few books declines `live_too_few_books` (named so the monitor separates a thin live slate from pre-#54 cache-era data). Confirm last look re-fetches **every** sub-combo live (see below). Every `quote_priced` research event carries a per-game `live_games` trace (leg-set hash, result age, per-book fair/route/latency) — combined with `on_demand_requested`/`on_demand_result` timestamps, the research DB proves each quoted fair traces to a fetch initiated after its RFQ landed. `test_live_pricing.py` pins all of this, including "cache present but never consulted".

Why live is right for this flow: ~40 in-scope RFQs/day × 6 books is a few hundred book hits/day vs the sweep's thousands (~10× less fingerprint surface), and quote-time staleness → ~0 — the adverse-selection window the staleness gates only bounded. DraftKings is never load-bearing: its price call is bimodal (0.8–1.6s usual, 8–22s stalls) and the `ON_DEMAND_DEADLINE_SEC=10` budget drops it; consensus needs `MIN_AGREEING_BOOKS=2` of the remaining books (warm p95s: FD 0.77s, PX 1.55s, MGM 1.66s, CZR 2.37s, NV 4.00s — `mlb_sgp/README.md`).

### Removed sweep (issue #81; demoted first in #57)

#57 demoted the full-slate sweep to a 300s background record; #81 deleted it outright (user decision 2026-08-09 — the measurement window the ticket's gate asked for was never produced because the bots were off, and the code audit showed the sweep record had NO remaining consumer: the monitor reads `sgp_fetch_health`, never `mlb_sgp_odds`; `research_queries.sql` reads fills/settlements; the firehose's `on_demand_result` events already record every flight's per-book fairs). What replaced its two residual jobs: `mlb_target_lines` upkeep is its own 300s Kalshi-only loop arm (`target_line_cycle`, zero book requests — it also inherits the game-metadata cache invalidation, since new games appear there now), and the per-pass book counts became the periodic `on_demand_coverage` research event fed by an unconditional per-book outcome tally on the service (deliberately NOT on the #37 alerter, which is None when alerting is disabled). The maker's `mlb_sgp_odds` table simply stops being written; the market DB's live tables are `mlb_target_lines` + `sgp_fetch_health`. Everything that used to gate on sweep freshness was already deleted in #57 (no discovery top gate; risk-sweep mass-cancel keys on #38 failure streaks via `BookHealthAlerter.consensus_dark()`; drift breaker prices from live engine fairs; post-fill cooldown clears on a **targeted** on-demand fetch — below). Rollback is a git revert.

### Consensus fair (shared math)

Fair value is the median of the devigged book fairs, gated on cross-book agreement. `router.combo_fair_detail` parses the RFQ into canonical legs (`legset`), partitions them by game, prices each game's sub-combo, and **multiplies** the per-game fairs (cross-game independence). For each **2-leg same-game grid** (family + cell resolved from the legs — spread×total keyed by game × **signed** spread_line × total_line (#70: negative = home margin market, positive = away margin market), moneyline×total keyed by game × total_line with `spread_line` NULL) we:

1. Take every book's devigged fair for the exact leg set from the RFQ's live on-demand fetch (#54/#81: the fetch prices the full 2^N partition per book — the full 4-cell grid for a 2-leg combo, no fallback — and devigs it in `kalshi_common.fair_value`; the router's grid-frame path survives only as the test oracle).
2. **Consensus gate (issue #20 — z-space dispersion threshold).** Transform each book fair via `norm.ppf` (probit); decline the combo (`too_few_books`) if fewer than `MIN_AGREEING_BOOKS` books priced it, or (`consensus_dispersion`) if the sample stddev of the z-values exceeds `SIGMA_Z_MAX`. There is NO outlier removal — a dissenting book is as likely the informed one (news mid-propagation) as a broken scrape, so genuine disagreement declines the quote instead of outvoting the dissenter. Constant width in z-space = the same amount of disagreement at every price level: the tolerated absolute gap tightens automatically at the tails (~2¢ at p=0.50 → ~0.6¢ at p=0.08), where the old absolute ±2¢ band tolerated 25% relative disagreement.
3. Fair = median of **all** books that priced the combo. The gate also reports `sigma_pts` — the stddev of the book fairs in probability points — which feeds the margin floor (below).

A **single leg** (within a cross-game combo) is priced by marginalizing that game's grid — summing the two devigged cells along the free axis. **Cross-game** fairs are the product of the per-game consensus fairs (combo `sigma_pts` propagates as relative variances adding: `σ_combo = fair·sqrt(Σ(σ_g/f_g)²)`); if any game fails the gate or is `unpriceable`, the combo is not quoted (fail-safe), with the first failing game's reason logged in `quote_decisions`. Novel same-game shapes (`on_demand` — 3+ legs, spread+ml, total+total, …) price via live book queries — see **On-demand pricing** below. The v1.1 explicit correlation-premium gate is deferred — see spec section 13.

### Margin (issue #19 — uncertainty-scaled)

Per side (p = side win-prob, N = number of games in the combo):

```
margin_pts = max( p·(1 − (1+TARGET_ROI)^−N),  MIN_MARGIN_PTS + K_SIGMA·σ_books )
bid        = (p − margin_pts) − maker_fee,  floored to the $0.001 grid
```

- **ROI part, compounded per game:** one `(1+r)` divisor per game — the same structure sportsbooks use to build multi-leg parlay hold (per-leg vig multiplies), because cross-game products multiply per-game fairs whose errors compound. At N=1 this is exactly the old constant-ROI price `p/(1+r)`.
- **Probability-point floor:** the constant-ROI cushion ≈ `0.029·p` collapses at longshot fairs (~0.3¢ at p=0.10) while absolute fair error does not shrink with p — so the floor keeps an absolute minimum cushion, widened by `K_SIGMA·σ_books` when the books visibly disagree. `MIN_MARGIN_PTS` covers error σ cannot see (shared devig/correlation bias, mirrored books). This mirrors how books load vig onto longshots (favorite–longshot distribution) and how options MMs widen deep-OTM quotes.
- **Quorum add-on (issue #55):** a quote whose thinnest per-game consensus is **exactly 2 books** carries one extra floor term — `floor = MIN_MARGIN_PTS + K_SIGMA·σ_books + QUORUM_MARGIN_ADDON` — because a 2-book σ is a 2-sample stddev, too noisy to price the visible disagreement alone. The add-on drops out automatically when a straggler book lands and the quote refines (see quorum quoting below); it is one term inside this framework, not a parallel margin system.
- The stepped-fee/grid-floor guard is unchanged: the bid steps down until `bid + fee(bid) ≤ p − margin_pts`, so the **realized** margin is ≥ target by construction; the `yes_bid + no_bid < 1` and positivity guards still apply (an unquotable side declines the RFQ).
- Margin components (`sigma_pts`, `floor_pts`, `roi_pts_*`, `margin_pts_*`, `n_games`) are logged in the research firehose `quote_priced` payload — `quote_decisions` columns are unchanged.

Interplay of the two tickets: **moderate** dispersion widens the margin (#19); dispersion **past `SIGMA_Z_MAX`** kills the quote (#20). One uncertainty number drives both.

## On-demand pricing (Phase 2)

Spec: `docs/superpowers/specs/2026-07-08-kalshi-mlb-mm-on-demand-pricing-design.md` (rev 5).

Since #54 (above) this engine prices **all** in-scope shapes — grids and
cross-game singles included, not just novel shapes: every quote is a live
query of the six books' SGP endpoints for that exact leg set, at RFQ time.
Always on — no config switch; the bot-wide kill file remains the emergency
stop, and dropping a misbehaving book is a one-line code change.

**Feed model — the open-RFQ poll drives everything.** The 2s discovery tick
never blocks on a book. An on-demand shape with no fresh result enqueues a
background fetch (`OnDemandEngine`, `kalshi_mlb_mm/on_demand.py`) and skips
with reason `on_demand_pending`; the fetch (all books concurrent, each book
capped at `ON_DEMAND_DEADLINE_SEC` — a straggler is dropped, never waited
on) lands per-book fairs in an in-memory store, and the next tick prices
them through the normal consensus + gates + hysteresis/replace path.

**Quorum quoting (issue #55) — landings are incremental.** Each book's
result becomes servable the moment it returns (`_Flight` in
`on_demand.py`), so the tick quotes as soon as **2 fresh books** pass the
#20 dispersion gate — RFQ→first-quote ≈ the 2nd-fastest book's response
(~1.6s at the #50 latency table) + one ≤2s tick, never bounded by DK's
tail. Freshness semantics for a partially-landed flight: `landed_at`
slides forward with each arrival (worst-case serve window unchanged:
job deadline + `QUOTE_FRESH_SEC`); a zero-book **partial** stays
`on_demand_pending`; `live_fetch_timeout` fires only when a **completed**
flight landed empty; a book finishing after its flight completed is
discarded (post-completion immutability keeps the confirm lane's
landed-after-entry guarantee). Stragglers then refine through the
existing machinery as they land: within `QUOTE_HYSTERESIS` → hold;
beyond → cancel/replace (#16 orphan-safety: submit first, mark
`replaced` only on success); a straggler that busts the dispersion gate
**pulls** the resting quote (decision `pulled`, reason
`quorum_dispersion_bust`) — the resting price came from less information
than we now have, and what we have says the books disagree. Exactly-2-book
quotes carry `QUORUM_MARGIN_ADDON` (see Margin above), so expect one
refine-replace per quote when the 3rd book lands and the add-on drops out.

**Pass latency IS quote latency (issue #89).** The first live measurement
of #55 showed 18.5s median RFQ→quote on RFQs whose median lifetime is 10s
— not because quorum landings weren't served (they were; the engine tests
pin it) but because the discovery pass itself ran 8–12s: the in-scope
gate path opened 4–5 DB connections **per RFQ** (~17ms each on the live
state DB), so the RFQ re-visit cadence was the pass duration, not
`DISCOVERY_SEC` (2s). The fix loads one `_PassSnapshot` per pass (open
quotes + `quote_games`, cooldowns + last fills, per-combo and per-creator
fill aggregates — ONE read connection) and answers every per-RFQ gate
from memory; quotes submitted mid-pass fold back into the snapshot so the
cap gates keep their read-your-writes semantics. A 200-RFQ pass drops
from ~8.4s of connection churn to ~0.25s. The pass-summary log line now
carries `pass_sec` and fires unconditionally whenever a pass exceeds
`DISCOVERY_SEC` — a slow pass is a quote-latency regression even when
nothing was deferred. Regression tests:
`test_discovery_pass_constant_connections_for_in_scope_rfqs` (O(1)
connections) and `test_staggered_landings_quote_on_next_tick_not_deadline`
(quote fires on the 2-book partial, never at the flight deadline).
Multi-game RFQs enqueue one job per game; jobs run on up to 4 concurrent
daemon threads (issue #50 — previously strictly serial) while the #40
pacing invariant still holds: a per-book gate keeps at most ONE pricing
call in flight per book, so concurrent jobs pipeline per book instead of
bursting into a rate limit. A background structure-warming pass
(`STRUCTURE_WARM_SEC`) keeps every book's events/structure caches and
Caesars' WAF token hot so the live fetch never pays cold discovery. A result may
back a NEW quote only within `QUOTE_FRESH_SEC` (15s, module constant) of
landing — after that the next tick that still sees the RFQ re-fetches. When
the RFQ leaves the poll, fetching stops; nothing self-schedules, and there is
no price reuse across RFQs.

**Two devig routes** (`kalshi_common/fair_value.py`), selected per (book,
combo) by `SGPService.price_on_demand`:

- **Route A — full partition** (≤3 legs, all sides offered): price all 2^N
  side-combination cells (`legset.enumerate_partition`; cell 0 = target,
  bit j of cell i flips leg j), probit-devig across the partition
  (`devig_partition`), read the target cell. Overround gate scales with leg
  count (`1 ≤ Σ(1/dec) ≤ 1 + 0.25·N`, live-calibrated — real FD 2-leg partitions sum to ~1.28). No fallback within the route — any
  missing/insane cell abandons to Route B.
- **Route B — correlation transfer** (any N, or any one-sided leg): one SGP
  call; `fair = ∏(devigged singles) × [SGP implied / ∏(vigged singles)]`
  (the book's SGP margin ≈ compounded leg vig cancels in the ratio). Singles
  come free from the structure fetch at FD/PX/NV/MGM/CZR; DK's structure has
  IDs only, so DK singles are priced as 1-leg `calculateBets` calls (only
  when Route B is actually taken). Result gated by exact **Fréchet bounds**
  (`max(0, Σp−(n−1)) ≤ fair ≤ min(p)`) — assumption-free joint-probability
  limits.

Consensus is the same `MIN_AGREEING_BOOKS=2` + band rule as grids. Per-book
`resolve_legs` / `price_selection_set` live in each `mlb_sgp/<book>.py`,
reusing the orchestrators' own lookup structures (`_priceable_in_phase1`
and `test_router_integration.py` are the retained Phase 1 regression
oracles).

**Accepted quotes get a live last look**, two layers (issue #17):

1. **Singles-move veto (all combo shapes).** At quote time the bot snapshots
   the raw Kalshi odds of every leg (each leg is its own Kalshi singles
   market; `live_quotes.leg_prices_json`) — no snapshot, no quote. On accept
   it re-reads those markets live (one GET per leg, <~1s — fast enough for
   the 2s High-Volatility confirm window) and voids (`voided_singles_moved`)
   if ANY leg's bid or ask moved even one tick. Rationale: over the seconds
   between quote and accept the legs' correlation is structural and static —
   the marginals carry the pickoff signal. Missing baseline / failed read /
   unparseable snapshot all void ("can't verify ⇒ don't confirm").
2. **Live re-price (EVERY sub-combo, #54).** Accepted quotes get the
   synchronous priority re-fetch (budgeted `CONFIRM_REFETCH_BUDGET_SEC=20s`)
   for every game they touch — the quote's fair came from the engine, whose
   result has aged out by accept time, and the cache must not stand in. The
   confirm gate honours `refetch_now`'s landed-after-entry guarantee (a
   failed re-fetch voids even if a pre-entry result still sits in the
   15s-fresh store) — we never confirm on a previously-fetched number. The
   void/veto rate is the live-pricing health signal: re-checking a
   seconds-old number should collapse it vs the cache era (research query
   9); if it does not, live pricing is not actually working.

**Failure direction is always "fewer/laggier quotes", never "staler
quotes"**: unresolvable leg → book drops; incomplete partition → Route B →
drop; sanity/Fréchet/consensus failure → no quote; worker death → lazy
restart, meanwhile skip; RFQ flood → queue grows, late landings, expired
RFQs go unquoted (pacing = one on-demand combo in flight per book).

**Observability**: research events `on_demand_requested` (once per fetch
flight), `on_demand_result` (once per landing — with incremental landings
that is once per book arrival, so per-book arrival order and latency are
recorded; payload also carries `dropped_books`, the books over their #50
budget — a chronically-dropped book is a de facto non-participant for
#38), `quote_priced` with the per-game `live_games` trace plus (#55)
`quorum_size`, `quorum_addon_pts`, `refine_action`
(initial/hold/replace) and `resting_quote_id`, `quote_pulled` on a
dispersion-bust pull, and fill payloads gain per-book `route` +
`consensus_books` (DK+Novig-only fills are a named risk metric). New
`quote_decisions` reasons: `live_fetch_timeout`, `live_too_few_books`,
and decision `pulled` / reason `quorum_dispersion_bust` (the monitor
reads reason vocabularies from data — no dashboard change).
Out-of-scope RFQs now store `legs_json` and granular reasons
(`out_of_scope_non_mlb` / `out_of_scope_lone_single` /
`out_of_scope_unparseable` since #84; `out_of_scope_f5_winner` was
emitted #84→#86 and became `out_of_scope_f5_tie_leg` — TIE legs only —
when #86 made team legs priceable via ±0.5 re-encoding) in
`seen_rfqs` — on-demand demand is finally
measurable (it was 0 in ~1,100 classified markets over the 30 days before
this shipped).

**Rollout / retreat playbook**: restart from main repo cwd in `--dry-run`
first; watch route mix, `route_gap`, per-book completion, void rate
on-demand-vs-grid, Caesars WAF health, PX on-demand volume (each PX fetch is
a real RFQ on their exchange). Retreat triggers (all code changes, no
config): drop CZR from on-demand if its WAF failures spill into the warming
pass; haircut/disable Route B if `route_gap` shows it systematically rich;
reinstate a pull loop only if on-demand voids threaten the global void-rate
halt.

The model (a fraction-of-sample-paths estimate driven by `mlb_game_samples` from the R answer-key pipeline) was removed in the v1 hardening pass: it was being medianed out of the blend, carried documented bias on certain combo families ([[mlb_parlay_edge_overestimation]]), and added a soft dependency on the R pipeline.

For each side of a two-sided quote:

```
p        = side win-probability (fair for YES, 1 - fair for NO)
raw_bid  = p / (1 + TARGET_ROI)            # = p / 1.05
bid      = raw_bid - maker_fee_per_contract(raw_bid)   # one refinement step
bid      = round_down_to_grid(bid)         # $0.001 (deci-cent) grid, round DOWN

assert yes_bid + no_bid < 1               # sum ≈ 0.94 at fair=0.50, always valid
```

The 5% is a *quoted/expected* ROI — what we actually realize is what the validation phase measures. Maker fee = 25% of the taker fee on the same quadratic base (`maker_fee_per_contract` in `kalshi_common/ev_calc.py`). **Verify the exact charge on the first real fill** — the assumption is strongly implied by the Kalshi fee schedule but not yet confirmed against a real fill.

## Reports

`report.py` (issue #14) aggregates the three bot DBs into a printed markdown
"state of the maker" — the measurement-phase readout that the raw
`quote_decisions` / research-firehose rows can't answer at a glance:

```bash
# From the main repo root — reads the live DBs READ-ONLY (lock-safe, never writes)
./kalshi_mlb_mm/venv/bin/python -m kalshi_mlb_mm.report

# Or point at explicit DB copies
./kalshi_mlb_mm/venv/bin/python -m kalshi_mlb_mm.report \
    --state-db path/to/kalshi_mlb_mm.duckdb \
    --research-db path/to/kalshi_mlb_mm_research.duckdb
```

Sections: **1** RFQ funnel (24h + 7d: seen → in-scope → quoted → accepted →
confirmed → filled, full decision/reason breakdown, and the grid /
on-demand / out-of-scope path split), **1b** on-demand fetch success rate +
latency (paired `on_demand_requested`/`on_demand_result` events), **2**
quotable universe (combos passing `MIN_AGREEING_BOOKS` per day + per-book
participation), **3** staleness (book-data age at quote, quote age at
accept), **4** demand curve (quotes vs accepts by margin × fair-prob band —
feeds #26), **5** settlement P&L (quoted vs realized margin per contract
from the #12 sweep; **markouts descoped with #13** — settlement P&L is the
adverse-selection signal), **6** health (void rate, sweep cancels, halts,
phantom fills, reconcile outcomes).

Caveats printed in the report itself: the funnel is derived from
`quote_decisions` (`seen_rfqs` only records out-of-scope and live-quoted
RFQs); the grid-vs-on-demand split is a heuristic (an on-demand RFQ whose
fetch lands within its first discovery tick counts as grid); quote-time
staleness is exact since #81 — the max sub-fetch age from each
`quote_priced` event's `live_games` trace. Fresh/empty/locked DBs degrade
to notes — the report never crashes and never writes.

## Per-book fetch health

Every SGP fetch this bot makes — every on-demand live fetch and every
warming pass (#81: those are the only paths left) — writes one row to
`sgp_fetch_health` in
`kalshi_mlb_mm_market.duckdb` (book, path, outcome, duration, error class,
#35 counters). It is the memory that makes "DraftKings has been 403ing since
Tuesday" answerable. Buffered and flushed once per tick; a locked DB costs
history, never a quote. 30-day retention.

```bash
duckdb kalshi_mlb_mm/kalshi_mlb_mm_market.duckdb -readonly \
    < kalshi_common/fetch_health_queries.sql
```

Schema, outcome vocabulary and the units caveat (**always group by `path`**)
are documented in `mlb_sgp/README.md` → "Fetch-health history".

**Run-time alerting (#37).** `check_book_health()` runs every tick off
in-memory streaks fed by the same two call sites that write those rows.
`BOOK_ALERT_STREAK` (default 3) consecutive non-`ok`/non-`empty` outcomes on
one (book, path) → one ERROR log + one macOS notification, once per incident,
re-armed only when the book answers again. If healthy books drop below
`MIN_AGREEING_BOOKS` — the same gate that decides whether this bot can quote
at all — it alerts "quoting effectively blind". `empty` is a healthy answer
and a book skipped by `min_refresh_sec` is not a failure; see
`mlb_sgp/README.md` → "Run-time book-health alerting" for the full semantics,
the known silent-parser gap, and why this reads memory instead of the DB.
Knobs: `BOOK_ALERT_ENABLED`, `BOOK_ALERT_STREAK`, `BOOK_ALERT_PATHS`.

## Expired-quote outcome labels

Almost every quote expires unfilled — but *why* matters for margin tuning.
`expiry_outcome.py` (own sweep arm, `EXPIRY_OUTCOME_SWEEP_SEC` = 600s,
runs in dry-run too) answers it per quote from the combo market's **public
trade tape** (`GET /markets/trades?ticker=...` — live-verified 2026-08-13;
MVE-combo quirk: `yes_price`/`no_price`/`count` come back `None`, so only
`created_time` + `taker_side` are used):

- `competitor_traded_in_window` — the market printed a trade inside our
  quote's resting window (`submitted_at`..`closed_at`, both converted to
  UTC): a competing maker won the RFQ — **we're being outpriced**.
- `traded_shortly_after` — first print landed within
  `EXPIRY_OUTCOME_GRACE_SEC` (300s) after our quote closed: a near-miss.
- `no_trade` — nothing printed: the creator never executed with anyone —
  **cutting margin would donate edge**.

Each quote is labeled exactly once, after its grace window has fully
elapsed, into `quote_expiry_outcomes` in the state DB (PK `quote_id` is the
dedupe marker; also stores `n_trades_in_window`,
`first_trade_after_close_sec`, and the UTC window bounds) plus one buffered
`quote_expiry_outcome` research event. One tape fetch per distinct ticker
per sweep, capped at `EXPIRY_OUTCOME_MAX_TICKERS_PER_SWEEP`; failures are
logged and retried next sweep — nothing here touches the quote path.
Quotes WE cancelled (risk pulls, tipoff, breakers) are excluded unless
`EXPIRY_OUTCOME_INCLUDE_CANCELLED=true`. The daily label ratio is query 11
in `research_queries.sql`.

## Knobs

All knobs are overridable via `kalshi_mlb_mm/.env` or environment variables. Defaults come from `kalshi_mlb_mm/config.py`.

| Knob | Default | Purpose |
|---|---|---|
| `BANKROLL` | `500.0` | Master risk dial — raise this one number to scale everything |
| `DAILY_EXPOSURE_CAP_PCT` | `0.75` | Daily hard stop as a fraction of BANKROLL ($375 at default). Counts today's fills PLUS open quotes' worst-case exposure |
| `MAX_GAME_EXPOSURE_PCT` | `0.10` | Per-game exposure cap as fraction of BANKROLL ($50 at default). Counts fills (`fill_games`) PLUS open quotes (`quote_games`) touching the game |
| `MAX_FILL_EXPOSURE_PCT` | `0.10` | Per-fill dollar cap as fraction of BANKROLL ($50 at default). Quote-or-skip only — the RFQ creator fixes fill size; this is the only lever. |
| `MAX_OPEN_QUOTES` | `25` | Cap on simultaneously resting quotes (well under Kalshi's 100 limit) |
| `TARGET_ROI` | `0.03` | Quoted margin — the ROI part of the per-side margin, compounded per game: `p·(1 − (1+TARGET_ROI)^−n_games)` |
| `MIN_MARGIN_PTS` | `0.01` | Probability-point margin floor (issue #19) — minimum absolute cushion per side, covering fair error the books' visible disagreement can't measure |
| `K_SIGMA` | `1.0` | Margin widening per unit of cross-book disagreement: floor = `MIN_MARGIN_PTS + K_SIGMA·σ_books` |
| `QUORUM_MARGIN_ADDON` | `0.01` | Extra probability-point floor term on quotes whose thinnest per-game consensus is exactly 2 books (issue #55) — a 2-sample σ underprices visible disagreement; drops out when a straggler book lands |
| `FAIR_DRIFT_TOLERANCE` | `0.02` | Last-look: void confirm if fair drifted >2¢ against filled side since quote time |
| `BOOK_MOVE_CB_THRESHOLD` | `0.03` | Circuit breaker: cancel a combo's quotes if the live book fair jumps >3¢ between consecutive pricings (per-tick) or if drift since quote exceeds this (per-quote risk sweep) |
| `TIPOFF_CANCEL_MIN` | `5` | Pull quotes this many minutes before first pitch |
| `QUOTE_HYSTERESIS` | `0.005` | Don't replace a resting quote unless fair moved more than ½¢. Also the ε of the post-fill `same_price_block` |
| `COMBO_COOLDOWN_SEC` | `60` | Hard FLOOR of the post-fill per-combo cooldown. The combo additionally stays cooled until every game it touches has a completed post-fill TARGETED fetch (#57), then until consensus fair moves off the filled fair (defense item 5) |
| `MAX_COMBO_EXPOSURE_USD` | `50.0` | Per-combo concentration cap (H8/N7): fills + in-flight open quotes on one combo may not exceed this |
| `SIGMA_Z_MAX` | `0.07` | Consensus gate (issue #20): max sample stddev of the books' fairs in z-space (`norm.ppf`); above it the combo is declined (`consensus_dispersion`) — no outlier removal. 0.07 ≈ continuity with the old ±2¢ band at p=0.50 |
| `MIN_AGREEING_BOOKS` | `2` | Consensus gate: minimum number of books with a fair for the combo (no longer "band survivors" — there is no band) |
| `CORR_SANITY_FRECHET_ENABLED` | `true` | Reject a quote whose fair violates the Fréchet bounds implied by the live Kalshi marginals (issue #23) — parameter-free, so it gates by default |
| `CORR_SANITY_PREMIUM_ENABLED` | `false` | Also reject on a correlation premium outside the band. Ships **log-only**: the band is a guess until the `corr_sanity_check` firehose shows the real distribution |
| `CORR_PREMIUM_MIN` / `CORR_PREMIUM_MAX` | `0.5` / `2.0` | Band for `combo_fair / Π(live marginals)` — deliberately wide, since same-game stacks legitimately run ~2× |
| `CONSTITUENT_JUMP_THRESHOLD` | `0.03` | Devigged move on a constituent single, since quote placement, that pulls every resting quote on that game (issue #23) |
| `CONSTITUENT_POLL_BUDGET_SEC` | `1.0` | Wall-clock ceiling on the constituent poll. All ticks share one thread, so this protects the 2s HVM confirm window; the poll order rotates so no ticker is starved |
| `KALSHI_WS_URL` | `wss://api.elections.kalshi.com/trade-api/ws/v2` | WS endpoint (signed handshake via `kalshi_common.auth_client.ws_auth_headers()`) |
| `WS_HEARTBEAT_TIMEOUT_SEC` | `30` | Watchdog: no WS frame (pings included) for this long → REST fallback, loudly; the reader thread also tears down + reconnects a subscribed socket this stale (silent TCP death). Kalshi pings ~10s, so 30 = three missed heartbeats |
| `WS_RECONNECT_BASE_SEC` / `WS_RECONNECT_MAX_SEC` | `1` / `30` | Reconnect backoff bounds (exponential + jitter) |
| `WS_CONNECT_ALERT_STREAK` | `3` | Consecutive connect failures before the one-per-outage notification |
| `DISCOVERY_SEC` | `2` | Discovery + quote loop cadence (seconds) |
| `SCOPE_FETCH_BUDGET_PER_TICK` | `40` | Max REST market fetches (scope checks) per discovery tick. Scope normally resolves free from the RFQ's own `mve_selected_legs`; this bounds the `get_market` fallback for payloads missing that field, so an exchange-wide open-RFQ backlog (4k+ tickers, 2026-08-10 incident) can't starve the loop's other arms or wedge SIGTERM; excess tickers wait for a later tick |
| `MAX_RFQ_AGE_SEC` | `30` | Skip RFQs older than this before any work. An RFQ resting 10+ min un-quoted is near expiry or already picked over — quoting late is adverse selection, and classifying it burns lap budget. Missing/unparseable `created_ts` passes (fail-open); skips are counted in the periodic pass-summary log, not per-RFQ decision rows |
| `MIN_RFQ_TARGET_COST_USD` | `250` | Door-filter dollar floor (option B, 2026-08-11: NO fetch ceiling by user decision — the door gates are the only limit on book traffic). Applies only to dollar-denominated RFQs; contracts-denominated ones meet the tick's size gate instead |
| `MAX_RFQ_LEG_COUNT` | `3` | Door-filter leg cap: full-devig rigor stops at 3 legs and the 4-8-leg flood is algo basket spam |
| `DISCOVERY_PASS_BUDGET_SEC` | `20` | Wall-clock cap on one discovery pass. The WS mirror can serve the whole exchange's open-RFQ set (9k+ at Monday-evening peak, 2026-08-10); the lap processes what it can, defers the rest, and resumes at a rotating cursor next tick so tail RFQs never starve. Always makes progress (≥1 RFQ) |
| `CONFIRM_SEC` | `2` | Confirm loop cadence (seconds) |
| `RISK_SWEEP_SEC` | `10` | Risk sweep cadence (seconds) |
| `RECONCILE_SWEEP_SEC` | `30` | Fill side/size reconciliation cadence against Kalshi positions (seconds) |
| `SETTLEMENT_SWEEP_SEC` | `600` | Settlement sweep cadence (seconds) — populates `realized_pnl` once markets settle; only matters hours post-game |
| `EXPIRY_OUTCOME_SWEEP_SEC` | `600` | Expired-quote outcome labeler cadence (seconds) — labels each expired unfilled quote from the combo's public trade tape: competitor traded in our window vs no trade at all (the margin-tuning ratio) |
| `EXPIRY_OUTCOME_GRACE_SEC` | `300` | Trades landing within this many seconds AFTER our quote closed label `traded_shortly_after`; labeling waits until this window has elapsed so each quote is labeled exactly once |
| `EXPIRY_OUTCOME_MAX_TICKERS_PER_SWEEP` | `25` | Tape fetches per labeler sweep — only bounds the catch-up burst after downtime |
| `EXPIRY_OUTCOME_INCLUDE_CANCELLED` | `false` | Also label quotes WE cancelled (risk pulls, tipoff, breakers). Off by default — pulls are our own decisions and muddy the headline ratio |
| `TARGET_LINE_REFRESH_SEC` | `300` | #81: `mlb_target_lines` refresh cadence (Kalshi MVE enumeration + Odds API schedule; zero book requests). Sets how fast a NEW game becomes quotable |
| `COVERAGE_SUMMARY_SEC` | `300` | #81: cadence of the `on_demand_coverage` research event (per-book live-fetch outcome tally since the last summary; idle windows emit nothing) |
| `STRUCTURE_WARM_BUDGET_SEC` | `360.0` | #81: wall budget for one warming pass (pre-#81 warming rode the sweep's per-book deadline, which the live env had raised to 360 — this keeps that proven value). A book still running at the budget is dropped with a warming-path timeout health row |
| `STRUCTURE_WARM_SEC` | `120` | Structure-only warming cadence (issue #50): every book's events/structure TTL caches + Caesars' WAF token are re-warmed with ZERO pricing calls, so an RFQ never pays cold-structure discovery. Keep under `STRUCTURE_TTL_SEC` (180) and the CZR token TTL (240) |
| `FLIGHT_HORIZON_HOURS` | `6` | Fly flights only when every game in the combo starts within this many hours — books price SGP combos near game time, so far-out flights waste fetches on guaranteed too_few_books. 0 disables |
| `ON_DEMAND_MAX_CONCURRENT_JOBS` | `16` | Concurrent pricing jobs in the on-demand engine (was hard-coded 4, which capped throughput at ~50 combos/min vs ~170/min option-B demand). Per-book pressure unchanged at any value — per-book gates still serialize to one call in flight per book |
| `ON_DEMAND_DEADLINE_SEC` | `10.0` | Per-book wall budget for LIVE (on-demand) pricing fetches (issue #50). A book still running at the cap is dropped; the fast books' results land. Sized so warm Novig (p95 ~9s) barely fits |

## Defense hierarchy (stale-quote / adverse-selection risk)

Resting quotes are priced off books that lag reality (books refresh every 60s). An informed counterparty can cross a stale quote before our data reacts. Defenses, in decreasing priority:

1. **Margin (primary, continuous coverage).** The per-side cushion `max(ROI part, MIN_MARGIN_PTS + K_SIGMA·σ_books)` absorbs *continuous* fair drift — a cent or two of movement between scrapes — and, since issue #19, keeps an absolute floor at longshot fairs (where the old flat-ROI cushion collapsed to ~0.3¢) and widens with cross-book disagreement and game count. This is the primary defense and fires on every fill. It does NOT cover discrete events (scratch / postponement / steam move).

2. **Book-consensus gate (issue #20 — z-space dispersion threshold).** Before quoting, we require `>= MIN_AGREEING_BOOKS` books to have priced the combo AND their fairs to agree within `SIGMA_Z_MAX` in probit space. No outlier removal: one loudly-dissenting book declines the whole combo — it may be the informed one, and quoting the stale majority's median is adverse selection by construction (a false decline is ~free; a false quote is not). This gate measures whether the books agree with EACH OTHER; the explicit correlation-premium gate that measures their consensus against an OUTSIDE anchor (spec section 13) is now implemented as defense item 3 (issue #23).

3. **Correlation sanity vs Kalshi's own singles (issue #23).** Every leg of a combo is its own real-time 2-way Kalshi market. Before submitting, the devigged marginals of all legs are read — the same snapshot item 9's veto takes, so **no extra API calls** — and the combo fair is checked against them: it must respect the **Fréchet bounds** `max(0, Σp − (n−1)) ≤ fair ≤ min(p)` and its implied **correlation premium** `fair / Πp` must sit inside `[CORR_PREMIUM_MIN, CORR_PREMIUM_MAX]`. This is the only defense that can catch the books being *jointly* wrong: item 2 only sees whether they agree with **each other**, and tightly-agreeing books get a *thin* margin under item 1 — so a shared error is exactly the case with the least cushion and no other detector. Fréchet gates by default (`corr_sanity_frechet`); the premium band ships log-only (`CORR_SANITY_PREMIUM_ENABLED=false`) because same-game stacks legitimately price well above 1× (run line + moneyline ~2×) and a guessed band would decline real business. Every quote emits a `corr_sanity_check` research event carrying marginals, baseline, premium and both bounds — the evidence for where the band belongs, and the calibration dataset for a Phase-3 correlation model. Degenerate books (`yes_ask=100`, empty 0–100, crossed) yield no marginals: the miss is logged and the quote proceeds on book consensus alone rather than declining on missing information.

4. **Constituent-jump circuit breaker (issue #23).** Our books refresh every ~150–165s; the constituent singles trade in real time, which makes them the fastest evidence that the market has left a resting quote behind. Every risk sweep (`RISK_SWEEP_SEC`, 10s) the bot polls the **distinct** constituent tickers of all resting quotes and compares each against the placement baseline in `live_quotes.leg_prices_json`. A devigged move past `CONSTITUENT_JUMP_THRESHOLD` cancels **every** resting quote touching that game (`constituent_jump`), not just the quote whose own leg moved. Evaluated *before* the item 6 book-drift check, so the faster and more specific signal wins the logged reason. The jump is unconditional (since #54): every quote is priced from a live fetch, so the fair was fresh at placement and ANY post-placement jump means the market moved after us. The pre-#54 `book_quiet` guard — which additionally required our book consensus to have stayed put, excusing cache-priced quotes catching up to their own stale data — was deleted along with cache pricing (there is no case left for it to distinguish). An unreadable constituent is deliberately **not** treated as a jump: a transient API failure must not flush the resting book, and item 9's veto still fails closed on any resulting accept. **Rate load: one `GET /markets/{ticker}` per DISTINCT constituent ticker per sweep** — bounded by (open quotes × legs), far smaller in practice because quotes on the same game share legs, and exactly zero when flat or when the books-unhealthy pull is cancelling everything anyway. The poll is additionally capped by `CONSTITUENT_POLL_BUDGET_SEC` of wall clock: every tick runs on ONE thread, so an unbounded poll inside the sweep would delay the confirm tick, and Kalshi allows only **2s to confirm in High Volatility Markets** — blowing a confirm window is a far worse failure than polling fewer tickers this pass. When the budget bites, the iteration order rotates each sweep so no ticker is starved; a large resting book is therefore covered over several passes rather than all at once. Wide *and* instant coverage is what the WebSocket feed is for — the transport sits behind a small interface precisely so that swap is one adapter.

5. **Books-unhealthy auto-pull (#38 health-dark, re-keyed by #57).** The old rule read sweep-row age, which under the slow background sweep is meaningless. Now: the shared `BookHealthAlerter` tracks consecutive live-fetch failures per book (the same streaks #37 alerts on, `on_demand` path only), and when fewer than `MIN_AGREEING_BOOKS` books are healthy (`consensus_dark()`) the risk sweep cancels every open quote (`books_unhealthy`) — resting quotes can no longer be re-verified, so stop resting risk. Cadence-invariant, and fail-safe on startup: books never observed are unknown, not unhealthy, so a fresh process can't false-pull. Blind → no live quotes, same principle, live-keyed signal. Caveat: `BOOK_ALERT_ENABLED=false` disables the alerter entirely and therefore this pull too — the #23 constituent-jump/tipoff/drift breakers still run, but leave alerting on in live sessions.

6. **Book-move circuit breaker (discrete events).** Two layers: (a) per-tick — if a re-pricing shows a live book-fair jump greater than `BOOK_MOVE_CB_THRESHOLD` for a combo vs its prior pricing, the bot immediately cancels that combo's resting quotes (does not wait for the next discovery tick). (b) per-quote in the risk sweep — if current book consensus has drifted more than `BOOK_MOVE_CB_THRESHOLD` from the `book_fair_at_quote` stored when the quote was placed, the quote is cancelled. The per-quote sweep catches gradual drift the per-tick threshold misses (e.g., five 1¢ moves across ticks).

7. **Post-fill cooldown + same-price block (issue #21, re-pointed by #57).** After a fill, the combo is cooled for `COMBO_COOLDOWN_SEC` (hard floor) — and stays cooled past the floor until **every game the combo touches has a completed post-fill TARGETED on-demand fetch** (`in_cooldown_awaiting_refresh`). The confirm tick queues that fetch the instant the fill is recorded (research event `on_demand_requested` with `trigger="post_fill"`), so it lands within seconds and the gate normally clears right at the 60s floor; if the landing is lost (process restart, engine store prune) the cooldown gate lazily re-queues it — self-healing, never wedged, and never waiting on a sweep pass. "Landed" means the engine flight COMPLETED with ≥1 book fair after `filled_at` (`OnDemandEngine.completed_fetch_age_sec`). Once the books have been re-asked, the combo is still skipped (`same_price_block`) while the **live** consensus fair remains within `QUOTE_HYSTERESIS` of the fair the fill transacted against (`fair_at_confirm`, falling back to `blended_fair_at_quote`) — quote prices are fair ± a fixed ROI margin, so unchanged fair ⇔ the identical just-picked-off price. Both skip reasons are logged distinctly in `quote_decisions` for the funnel report. The refresh check fails CLOSED (no engine, missing landings, comparison errors keep the combo cooled).

8. **Tipoff blackout.** `TIPOFF_CANCEL_MIN` pulls all quotes for a game before first pitch. A cross-game combo is swept against the **earliest** first pitch across ALL its games (re-derived from `seen_rfqs.legs_json` each sweep — same legset path as discovery); any game the sweep can't resolve cancels the quote (fail-safe, never fail open).

9. **Last-look backstop (discrete events).** On accept, the bot first runs the **singles-move veto** (issue #17): the raw Kalshi odds of every leg, snapshotted at quote time, are re-read live (<~1s) and ANY one-tick move on any leg's bid or ask voids the fill (`voided_singles_moved`). Previously the last look re-priced grid combos from the same ~150s scrape cache the quote came from — `cur_fair == prev_fair` by construction, a no-op exactly in the fast-pickoff window. Combos that pass the veto still run the existing gates: (a) cannot re-price (no fresh books / failed on-demand re-fetch → `voided_no_fresh_books`), (b) no longer +EV (`price + fee >= current_fair`), (c) fair drifted past `FAIR_DRIFT_TOLERANCE` (→ `voided_last_look`). "Can't verify ⇒ don't confirm" is intentional throughout. Every accept emits a `confirm_singles_check` research event (quote age at accept, per-leg snapshot-vs-fresh odds, moved flag) — both the gate's save-rate and what a non-zero tolerance would have passed are measurable from day one. Non-confirms are abusive behavior Kalshi can throttle; the zero-tolerance veto's void rate is a watch item — if crowd noise voids too many clean accepts, the firehose deltas are the tuning dataset.

10. **Position reconciliation.** After every confirm we call `/portfolio/positions` for the combo ticker and trust Kalshi as the source of truth; if the confirm response's side or size disagrees, the `fills` row is written with the reconciled values and a `[position_mismatch]` warning is printed.

11. **Measure.** The `fills` table records `book_fair_at_quote`, `blended_fair_at_quote`, `fair_at_confirm`, and — via the **settlement sweep** (`settlement.py`, every `SETTLEMENT_SWEEP_SEC`; issue #12) — `realized_pnl` per fill: for each combo with reconciled unsettled fills the sweep polls `GET /markets/{ticker}` and, once Kalshi reports a yes/no result, writes `realized_pnl = contracts × ((won ? 1−price : −price) − fee)` plus a `settlements` audit row holding the raw market payload (the settlement response shape is unverified until real data lands). Only `reconciled=TRUE` fills settle, so a later side/size correction can never invalidate a written P&L; the sweep is idempotent and API failures just retry next pass. On the first settled ticker ever it also best-effort checks the maker-fee assumption against the `/portfolio` endpoints (`[FEE_MISMATCH]` logs loudly). Each settled fill emits a `settlement_recorded` research event, and `research_queries.sql` queries 7–8 decompose per-contract P&L into quoted margin vs fair drift (quote→confirm) vs settlement-vs-fair vs fee — the measurement-phase headline. If pickoffs swamp the margin, the honest conclusion is that making is not viable at this data latency → improve data speed (WebSocket feeds, faster scrape cadence) before scaling.

**Crash safety.** An unconfirmed quote that the process never confirms voids automatically — no surprise open position on crash. Graceful SIGTERM cancels all open quotes before exit. The main loop uses a 250ms sleep between sub-loop checks so SIGTERM is never starved (directly addresses the taker's known SIGTERM-starvation issue with its 640s RFQ refresh block).

## Open items to verify on first real fill

These are noted in spec §10 and remain unconfirmed until live fills happen:

1. **Maker fee** — confirm the actual charge matches 25% of taker fee (`maker_fee_per_contract` in `kalshi_common/ev_calc.py`). The settlement sweep automates a first pass: on the first settled ticker it fetches `/portfolio/fills` + `/portfolio/settlements`, compares any fee field found against our recorded fees (`[FEE_MISMATCH]` warning on disagreement), and captures the raw payloads in a `fee_verification` research event for manual inspection either way. Adjust the rounding in `ev_calc.py` if it differs.
2. **Quote-status polling shape** — exact fields on `GET /communications/quotes/{id}` (`status`, `accepted_side`, `contracts`). The confirm path in `main.py::_confirm_tick` infers `side_held` from `accepted_side` (`"no" if accepted_side == "yes" else "yes"`). Verify this field name and semantics on the first accepted quote.
3. **`get_competitors`** — competing quotes on an RFQ are stubbed out in v1 (future competitive pricing). Confirm the API shape if needed later.

## Accepted risks (v1)

Eight adversarial vectors are accepted for v1, measured not prevented. See spec §12 for the full table. The key ones: our quotes leak our fair surface (vector #1); fills are disproportionately from combos we underpriced (vector #2, the quiet killer — if per-combo fair error exceeds the margin, a patient sharp grinds us down regardless of other gates); the fav −1.5 + over family has documented model bias (vector #3). The validation dataset (`fills` + `realized_pnl` at settlement) is built specifically to measure fill-vs-fair error and answer "is 5% enough?"
