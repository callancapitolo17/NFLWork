# Kalshi MLB MM (Maker) Bot

Independent maker daemon that listens for others' RFQs on the Kalshi cross-category MVE collection, prices MLB combos against a book-consensus fair value, and provides two-sided quotes at a fixed 5% ROI margin. Coexists with the taker (`kalshi_mlb_rfq/`) as a separate OS process with no runtime dependency on it.

**Combo scope (Phase 1).** Pricing flows through `kalshi_common.legset` + `kalshi_mlb_mm.router`, which generalize the original fixed 2-leg path:

- **2-leg same-game grids** — spread×total and moneyline×total (both FG), priced from each book's 4-cell devig grid (unchanged math).
- **Cross-game combos** — each game's sub-combo is priced independently and the per-game fairs are **multiplied** (independence assumption); a single leg within a game is marginalized out of that game's grid.
- **On-demand same-game shapes (Phase 2)** — 3-leg, spread+ml, total+total and any other novel shape route to `on_demand` and are priced by live book queries at RFQ time (see **On-demand pricing** below). Lone single-leg RFQs remain out of scope (skipped fail-safe).

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

The read side is **book-agnostic**: it consumes whatever the 6 SGP scrapers write to `mlb_sgp_odds`, so all books' moneyline rows are priced with no maker-side change.

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

REST-polling daemon, single process. Six timed sub-loops:

| Loop | Cadence | Job |
|---|---|---|
| Discovery + quote | 2s | Poll open RFQs → scope filter → price → submit or refresh quote |
| Confirm | 2s | Poll open quotes → on `accepted`, last-look gate → confirm or void |
| Risk sweep | 10s | Kill-switch, books-unhealthy auto-pull (#38 health-dark), tipoff cancel, constituent-jump + drift-since-quote cancel |
| Reconcile sweep | 30s | Verify recorded fill side/size against Kalshi `/portfolio/positions` (live only) |
| Settlement sweep | 600s | Poll `GET /markets/{ticker}` for combos with unsettled fills → write `fills.realized_pnl` + `settlements` audit row (live only; issue #12) |
| Structure warming | 120s | #50: structure-only pass (no pricing calls) keeping every book's events/structure TTL caches + the CZR WAF token hot, so live fetches never pay a cold start |
| SGP scrape | 300s | **Background only (#57):** full-slate pricing pass → `kalshi_mlb_mm_market.duckdb` for research/monitor continuity + `mlb_target_lines` upkeep. Quoting NEVER reads it |

**Transport seam.** All exchange I/O goes through two interfaces: `RFQSource.poll()` / `RFQSource.get_market()` and `QuoteGateway.submit_quote()` / `.confirm()` / `.cancel()` (v1: `RestQuoteGateway`). RFQ discovery is **WS-only** (issue #56, user decision 2026-08-07 — no mode switch, mirroring #54's live-only decision; rollback is a git revert): `WebSocketRFQSource` (`ws_rfq_source.py`) runs a daemon reader thread subscribed to Kalshi's `communications` WS channel (`rfq_created` / `rfq_deleted`) and mirrors it into a stateful open-RFQ set, so `poll()` keeps the **level-triggered** contract the discovery tick depends on (all open RFQs every tick — that re-entry drives on-demand re-fetch, quorum refinement, and pulls). Payloads are stored verbatim (`contracts_fp` stays a string), so downstream is byte-compatible with the old REST poll. `RestRFQSource` (the 2s poll of `GET /communications/rfqs?status=open`) is not a mode: it survives only as the WS source's automatic fallback and gap-fill path, and `get_market()` stays pure REST.

**WS fallback semantics (never silently deaf).** Liveness = any WS frame (server pings every ~10s count) within `WS_HEARTBEAT_TIMEOUT_SEC`. The watchdog runs inside `poll()` on the tick thread — a wedged reader thread still trips it. Trip → ERROR log + one notify + `rfq_source_transition` research event, and `poll()` transparently serves the REST source until the feed recovers, which flips back the same loud way. Reconnects use backoff+jitter and re-subscribe; readiness is gated on the server's `subscribed` ack (an unacknowledged or rejected subscribe forces a reconnect — a live socket with a dead subscription passes the ping watchdog, the one deafness heartbeats can't see); each (re)connect performs ONE REST gap-fill that reconciles the mirror wholesale — additions AND removals, so a missed `rfq_deleted` can never leave a zombie RFQ being re-fetched forever (events buffered during the gap-fill re-apply on top, so nothing lands lost). Repeated connect failures (`WS_CONNECT_ALERT_STREAK`) alert once per outage. Every WS `rfq_created` emits `rfq_seen_latency` (payload `created_ts` vs receipt) to the firehose, so the discovery-latency win over the 2s poll is measured per RFQ.

**Shared math.** Fair value, EV calc (including `maker_fee_per_contract`), authenticated HTTP, SGP orchestration, and leg-typing helpers all live in `kalshi_common/`. Both bots import from there; the taker's original files are now one-line re-export shims (behavior unchanged).

**SGP pricing (2026-06).** The bot prices SGPs **in-process** via `kalshi_common/sgp_service.py::SGPService` — no subprocess per cycle. The service holds persistent per-book HTTP clients reused across cycles (no per-cycle TLS handshake) and prices the four books concurrently under a per-book deadline (`SGP_SCRAPER_TIMEOUT_SEC`). DK/FD structure fetches (event lists, selection-id dicts) are TTL-cached; prices are never cached, and a failed or timed-out book keeps its prior rows. The old subprocess-per-cycle model is retained as a rollback hatch — calling `sgp_cycle` without `service=`.

**State DB.** The bot writes `kalshi_mlb_mm/kalshi_mlb_mm.duckdb` (quotes, fills, positions, decisions, `fill_games`, and `quote_games`). The `fill_games` table maps each fill to **every** game its combo touches (one row per game); the per-game exposure cap reads a fanned-out `fills⋈fill_games` join so a cross-game combo's full stake counts against each of its games, while the daily cap and P&L keep reading `fills` (one row per combo) and are never double-counted. The per-game and daily caps additionally count **open quotes' worst-case exposure** (issue #22): each `live_quotes` row freezes `worst_exposure_usd` at quote time, `quote_games` fans it out per game (same attribution rule as `fill_games`, written in the same transaction as the quote insert), and both cap gates sum `fills + open quotes` so a burst sweep of resting quotes can no longer fill to multiples of the caps before the first fill registers. A `NULL` `worst_exposure_usd` (pre-migration row) is counted at the per-fill cap, mirroring the N8 unreconciled-fill convention; the RFQ currently being re-quoted is excluded from the sums (its resting quote is superseded, not stacked). The sibling `kalshi_mlb_mm/kalshi_mlb_mm_market.duckdb` holds SGP-line and SGP-odds data (same pattern as the taker's `kalshi_mlb_rfq_market.duckdb`). The v1-hardening pass removed the model component of the blend, so there is no longer a read-only dependency on `Answer Keys/mlb_mm.duckdb`. All timestamp columns are `TIMESTAMPTZ` UTC (repo rule); existing naive-local DBs are migrated in place on startup by the idempotent `MIGRATE_SQL` (instants preserved), and the monitor normalizes reads so both frames render identically.

**Quote replace ordering.** When a re-priced quote moves beyond `QUOTE_HYSTERESIS`, the bot submits the NEW quote first and only then marks the old row `replaced` — Kalshi auto-cancels the prior quote only on a successful resubmit, so a failed submit leaves the old row `open` (still confirm-polled and risk-swept) and the next tick retries.

## Pricing

### Live pricing (issue #54) — the only quote path

EVERY in-scope sub-combo — 2-leg grids, each game of a cross-game combo (including its single legs), novel shapes — rides the on-demand engine: the RFQ landing triggers a fetch of exactly that leg set at all six books, and the quote is priced from those fresh per-book fairs through the same #20 consensus gate and #19 margin machinery (the math is input-agnostic; `σ_books` comes from the live set). **The sweep cache is never in the quote path** (there is deliberately no cached-pricing mode or fallback — user decision 2026-08-05; rollback is a git revert). A fetch that lands with **zero** books declines with `live_fetch_timeout` — never a silent cache fallback — and is not re-fetched until the empty landing ages out (~`QUOTE_FRESH_SEC`), so a dead slate retries on a ~15s cadence instead of hammering per tick. A live consensus with too few books declines `live_too_few_books` (named so the monitor separates a thin live slate from pre-#54 cache-era data). Confirm last look re-fetches **every** sub-combo live (see below). Every `quote_priced` research event carries a per-game `live_games` trace (leg-set hash, result age, per-book fair/route/latency) — combined with `on_demand_requested`/`on_demand_result` timestamps, the research DB proves each quoted fair traces to a fetch initiated after its RFQ landed. `test_live_pricing.py` pins all of this, including "cache present but never consulted".

Why live is right for this flow: ~40 in-scope RFQs/day × 6 books is a few hundred book hits/day vs the sweep's thousands (~10× less fingerprint surface), and quote-time staleness → ~0 — the adverse-selection window the staleness gates only bounded. DraftKings is never load-bearing: its price call is bimodal (0.8–1.6s usual, 8–22s stalls) and the `ON_DEMAND_DEADLINE_SEC=10` budget drops it; consensus needs `MIN_AGREEING_BOOKS=2` of the remaining books (warm p95s: FD 0.77s, PX 1.55s, MGM 1.66s, CZR 2.37s, NV 4.00s — `mlb_sgp/README.md`).

### Demoted sweep (issue #57)

With every quote priced live, the full-slate sweep's residual jobs are (a) the market-wide `mlb_sgp_odds` record the monitor dashboard and research queries read, and (b) refreshing `mlb_target_lines` (leg→game resolution + the warming pass's game list). It runs at `SGP_REFRESH_SEC=300` for exactly those two, while #50's structure warming stays on its own hot 120s knob — warming, not sweeping, is what keeps live fetches fast. Everything that used to *gate* on sweep freshness was deleted, not slowed: the discovery tick has no sweep-freshness top gate (an empty `_SGP_ODDS` is the normal state between passes), the risk sweep's mass-cancel keys on #38 live-fetch failure streaks instead of data age (`books_unhealthy` via `BookHealthAlerter.consensus_dark()`), the drift breaker prices from live engine fairs, and the post-fill cooldown clears on a **targeted** on-demand fetch (below). `BOOK_ALERT_PATHS` counts only `on_demand`, so a slow — or entirely dead — sweep can never trip a #37 alert. The `book_requests_per_day` queries in `kalshi_common/fetch_health_queries.sql` measure the resulting request-volume drop from #38's health rows; full sweep removal is #81, decision-gated on that live data.

### Consensus fair (shared math)

Fair value is the median of the devigged book fairs, gated on cross-book agreement. `router.combo_fair_detail` parses the RFQ into canonical legs (`legset`), partitions them by game, prices each game's sub-combo, and **multiplies** the per-game fairs (cross-game independence). For each **2-leg same-game grid** (family + cell resolved from the legs — spread×total keyed by game × **signed** spread_line × total_line (#70: negative = home margin market, positive = away margin market), moneyline×total keyed by game × total_line with `spread_line` NULL) we:

1. Pull every book's 4-cell grid from `mlb_sgp_odds` (filtered by combo family, require all 4 cells — no fallback).
2. Devig each book's 4-way grid to a single combo fair (`devig_book` in `kalshi_common.fair_value`).
3. **Consensus gate (issue #20 — z-space dispersion threshold).** Transform each book fair via `norm.ppf` (probit); decline the combo (`too_few_books`) if fewer than `MIN_AGREEING_BOOKS` books priced it, or (`consensus_dispersion`) if the sample stddev of the z-values exceeds `SIGMA_Z_MAX`. There is NO outlier removal — a dissenting book is as likely the informed one (news mid-propagation) as a broken scrape, so genuine disagreement declines the quote instead of outvoting the dissenter. Constant width in z-space = the same amount of disagreement at every price level: the tolerated absolute gap tightens automatically at the tails (~2¢ at p=0.50 → ~0.6¢ at p=0.08), where the old absolute ±2¢ band tolerated 25% relative disagreement.
4. Fair = median of **all** books that priced the combo. The gate also reports `sigma_pts` — the stddev of the book fairs in probability points — which feeds the margin floor (below).

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
The live RFQ→first-quote latency vs the #50 table is a **post-restart
measurement deliverable** — the bots were off when this shipped.
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
reusing the orchestrators' own lookup structures — the grid-sweep path is
untouched (`_priceable_in_phase1` and `test_router_integration.py` are the
retained Phase 1 regression oracles).

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
`out_of_scope_unparseable`) in `seen_rfqs` — on-demand demand is finally
measurable (it was 0 in ~1,100 classified markets over the 30 days before
this shipped).

**Rollout / retreat playbook**: restart from main repo cwd in `--dry-run`
first; watch route mix, `route_gap`, per-book completion, void rate
on-demand-vs-grid, Caesars WAF health, PX on-demand volume (each PX fetch is
a real RFQ on their exchange). Retreat triggers (all code changes, no
config): drop CZR from on-demand if its WAF failures spill into the grid
sweep; haircut/disable Route B if `route_gap` shows it systematically rich;
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
    --research-db path/to/kalshi_mlb_mm_research.duckdb \
    --market-db path/to/kalshi_mlb_mm_market.duckdb
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
staleness is the gap to the latest prior `scrape_done` event (per-book ages
aren't persisted at quote time). Fresh/empty/locked DBs degrade to notes —
the report never crashes and never writes.

## Per-book fetch health

Every SGP fetch this bot makes — the periodic sweep and every on-demand live
fetch — writes one row to `sgp_fetch_health` in
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
| `BOOK_MOVE_CB_THRESHOLD` | `0.03` | Circuit breaker: cancel a combo's quotes if book fair jumps >3¢ between scrapes (per-tick) or if drift since quote exceeds this (per-quote risk sweep) |
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
| `WS_HEARTBEAT_TIMEOUT_SEC` | `30` | Watchdog: no WS frame (pings included) for this long → REST fallback, loudly. Kalshi pings ~10s, so 30 = three missed heartbeats |
| `WS_RECONNECT_BASE_SEC` / `WS_RECONNECT_MAX_SEC` | `1` / `30` | Reconnect backoff bounds (exponential + jitter) |
| `WS_CONNECT_ALERT_STREAK` | `3` | Consecutive connect failures before the one-per-outage notification |
| `DISCOVERY_SEC` | `2` | Discovery + quote loop cadence (seconds) |
| `CONFIRM_SEC` | `2` | Confirm loop cadence (seconds) |
| `RISK_SWEEP_SEC` | `10` | Risk sweep cadence (seconds) |
| `RECONCILE_SWEEP_SEC` | `30` | Fill side/size reconciliation cadence against Kalshi positions (seconds) |
| `SETTLEMENT_SWEEP_SEC` | `600` | Settlement sweep cadence (seconds) — populates `realized_pnl` once markets settle; only matters hours post-game |
| `SGP_REFRESH_SEC` | `300` | Background full-slate sweep cadence (#57) — research/monitor continuity + `mlb_target_lines` upkeep only; quoting never reads it. The taker keeps its own 60s value under the same knob name |
| `SGP_SCRAPER_TIMEOUT_SEC` | `90` | Per-book deadline passed to `SGPService` (seconds) — a book exceeding it contributes nothing that cycle and its client is rebuilt |
| `STRUCTURE_WARM_SEC` | `120` | Structure-only warming cadence (issue #50): every book's events/structure TTL caches + Caesars' WAF token are re-warmed with ZERO pricing calls, so an RFQ never pays cold-structure discovery. Keep under `STRUCTURE_TTL_SEC` (180) and the CZR token TTL (240) |
| `ON_DEMAND_DEADLINE_SEC` | `10.0` | Per-book wall budget for LIVE (on-demand) pricing fetches (issue #50). A book still running at the cap is dropped; the fast books' results land. Sized so warm Novig (p95 ~9s) barely fits — the 75s sweep budget never applies to the quote path |

## Defense hierarchy (stale-quote / adverse-selection risk)

Resting quotes are priced off books that lag reality (books refresh every 60s). An informed counterparty can cross a stale quote before our data reacts. Defenses, in decreasing priority:

1. **Margin (primary, continuous coverage).** The per-side cushion `max(ROI part, MIN_MARGIN_PTS + K_SIGMA·σ_books)` absorbs *continuous* fair drift — a cent or two of movement between scrapes — and, since issue #19, keeps an absolute floor at longshot fairs (where the old flat-ROI cushion collapsed to ~0.3¢) and widens with cross-book disagreement and game count. This is the primary defense and fires on every fill. It does NOT cover discrete events (scratch / postponement / steam move).

2. **Book-consensus gate (issue #20 — z-space dispersion threshold).** Before quoting, we require `>= MIN_AGREEING_BOOKS` books to have priced the combo AND their fairs to agree within `SIGMA_Z_MAX` in probit space. No outlier removal: one loudly-dissenting book declines the whole combo — it may be the informed one, and quoting the stale majority's median is adverse selection by construction (a false decline is ~free; a false quote is not). This gate measures whether the books agree with EACH OTHER; the explicit correlation-premium gate that measures their consensus against an OUTSIDE anchor (spec section 13) is now implemented as defense item 3 (issue #23).

3. **Correlation sanity vs Kalshi's own singles (issue #23).** Every leg of a combo is its own real-time 2-way Kalshi market. Before submitting, the devigged marginals of all legs are read — the same snapshot item 9's veto takes, so **no extra API calls** — and the combo fair is checked against them: it must respect the **Fréchet bounds** `max(0, Σp − (n−1)) ≤ fair ≤ min(p)` and its implied **correlation premium** `fair / Πp` must sit inside `[CORR_PREMIUM_MIN, CORR_PREMIUM_MAX]`. This is the only defense that can catch the books being *jointly* wrong: item 2 only sees whether they agree with **each other**, and tightly-agreeing books get a *thin* margin under item 1 — so a shared error is exactly the case with the least cushion and no other detector. Fréchet gates by default (`corr_sanity_frechet`); the premium band ships log-only (`CORR_SANITY_PREMIUM_ENABLED=false`) because same-game stacks legitimately price well above 1× (run line + moneyline ~2×) and a guessed band would decline real business. Every quote emits a `corr_sanity_check` research event carrying marginals, baseline, premium and both bounds — the evidence for where the band belongs, and the calibration dataset for a Phase-3 correlation model. Degenerate books (`yes_ask=100`, empty 0–100, crossed) yield no marginals: the miss is logged and the quote proceeds on book consensus alone rather than declining on missing information.

4. **Constituent-jump circuit breaker (issue #23).** Our books refresh every ~150–165s; the constituent singles trade in real time, which makes them the fastest evidence that the market has left a resting quote behind. Every risk sweep (`RISK_SWEEP_SEC`, 10s) the bot polls the **distinct** constituent tickers of all resting quotes and compares each against the placement baseline in `live_quotes.leg_prices_json`. A devigged move past `CONSTITUENT_JUMP_THRESHOLD` cancels **every** resting quote touching that game (`constituent_jump`), not just the quote whose own leg moved. Evaluated *before* the item 6 book-drift check, so the faster and more specific signal wins the logged reason. The jump is unconditional (since #54): every quote is priced from a live fetch, so the fair was fresh at placement and ANY post-placement jump means the market moved after us. The pre-#54 `book_quiet` guard — which additionally required our book consensus to have stayed put, excusing cache-priced quotes catching up to their own stale data — was deleted along with cache pricing (there is no case left for it to distinguish). An unreadable constituent is deliberately **not** treated as a jump: a transient API failure must not flush the resting book, and item 9's veto still fails closed on any resulting accept. **Rate load: one `GET /markets/{ticker}` per DISTINCT constituent ticker per sweep** — bounded by (open quotes × legs), far smaller in practice because quotes on the same game share legs, and exactly zero when flat or when the books-unhealthy pull is cancelling everything anyway. The poll is additionally capped by `CONSTITUENT_POLL_BUDGET_SEC` of wall clock: every tick runs on ONE thread, so an unbounded poll inside the sweep would delay the confirm tick, and Kalshi allows only **2s to confirm in High Volatility Markets** — blowing a confirm window is a far worse failure than polling fewer tickers this pass. When the budget bites, the iteration order rotates each sweep so no ticker is starved; a large resting book is therefore covered over several passes rather than all at once. Wide *and* instant coverage is what the WebSocket feed is for — the transport sits behind a small interface precisely so that swap is one adapter.

5. **Books-unhealthy auto-pull (#38 health-dark, re-keyed by #57).** The old rule read sweep-row age, which under the slow background sweep is meaningless. Now: the shared `BookHealthAlerter` tracks consecutive live-fetch failures per book (the same streaks #37 alerts on, `on_demand` path only), and when fewer than `MIN_AGREEING_BOOKS` books are healthy (`consensus_dark()`) the risk sweep cancels every open quote (`books_unhealthy`) — resting quotes can no longer be re-verified, so stop resting risk. Cadence-invariant, and fail-safe on startup: books never observed are unknown, not unhealthy, so a fresh process can't false-pull. Blind → no live quotes, same principle, live-keyed signal. Caveat: `BOOK_ALERT_ENABLED=false` disables the alerter entirely and therefore this pull too — the #23 constituent-jump/tipoff/drift breakers still run, but leave alerting on in live sessions.

6. **Book-move circuit breaker (discrete events).** Two layers: (a) per-tick — if a scrape shows a book-fair jump greater than `BOOK_MOVE_CB_THRESHOLD` for a combo vs the prior scrape, the bot immediately cancels that combo's resting quotes (does not wait for the next discovery tick). (b) per-quote in the risk sweep — if current book consensus has drifted more than `BOOK_MOVE_CB_THRESHOLD` from the `book_fair_at_quote` stored when the quote was placed, the quote is cancelled. The per-quote sweep catches gradual drift the per-tick threshold misses (e.g., five 1¢ moves across ticks).

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
