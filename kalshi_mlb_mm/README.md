# Kalshi MLB MM (Maker) Bot

Independent maker daemon that listens for others' RFQs on the Kalshi cross-category MVE collection, prices MLB combos against a book-consensus fair value, and provides two-sided quotes at a fixed 5% ROI margin. Coexists with the taker (`kalshi_mlb_rfq/`) as a separate OS process with no runtime dependency on it.

**Combo scope (Phase 1).** Pricing flows through `kalshi_common.legset` + `kalshi_mlb_mm.router`, which generalize the original fixed 2-leg path:

- **2-leg same-game grids** — spread×total and moneyline×total (both FG), priced from each book's 4-cell devig grid (unchanged math).
- **Cross-game combos** — each game's sub-combo is priced independently and the per-game fairs are **multiplied** (independence assumption); a single leg within a game is marginalized out of that game's grid.
- **On-demand same-game shapes (Phase 2)** — 3-leg, spread+ml, total+total and any other novel shape route to `on_demand` and are priced by live book queries at RFQ time (see **On-demand pricing** below). Lone single-leg RFQs remain out of scope (skipped fail-safe).

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
| Risk sweep | 10s | Kill-switch, book-staleness auto-pull, tipoff cancel, drift-since-quote cancel |
| Reconcile sweep | 30s | Verify recorded fill side/size against Kalshi `/portfolio/positions` (live only) |
| Settlement sweep | 600s | Poll `GET /markets/{ticker}` for combos with unsettled fills → write `fills.realized_pnl` + `settlements` audit row (live only; issue #12) |
| SGP scrape | 60s | Own scraper cadence → `kalshi_mlb_mm_market.duckdb` (sibling market DB) |

**Transport seam.** All exchange I/O goes through two interfaces: `RFQSource.poll()` / `RFQSource.get_market()` (v1: `RestRFQSource`, REST poll of `GET /communications/rfqs?status=open`) and `QuoteGateway.submit_quote()` / `.confirm()` / `.cancel()` (v1: `RestQuoteGateway`). A WebSocket adapter is a drop-in replacement behind these interfaces — no pricing or risk code changes.

**Shared math.** Fair value, EV calc (including `maker_fee_per_contract`), authenticated HTTP, SGP orchestration, and leg-typing helpers all live in `kalshi_common/`. Both bots import from there; the taker's original files are now one-line re-export shims (behavior unchanged).

**SGP pricing (2026-06).** The bot prices SGPs **in-process** via `kalshi_common/sgp_service.py::SGPService` — no subprocess per cycle. The service holds persistent per-book HTTP clients reused across cycles (no per-cycle TLS handshake) and prices the four books concurrently under a per-book deadline (`SGP_SCRAPER_TIMEOUT_SEC`). DK/FD structure fetches (event lists, selection-id dicts) are TTL-cached; prices are never cached, and a failed or timed-out book keeps its prior rows. The old subprocess-per-cycle model is retained as a rollback hatch — calling `sgp_cycle` without `service=`.

**State DB.** The bot writes `kalshi_mlb_mm/kalshi_mlb_mm.duckdb` (quotes, fills, positions, decisions, `fill_games`, and `quote_games`). The `fill_games` table maps each fill to **every** game its combo touches (one row per game); the per-game exposure cap reads a fanned-out `fills⋈fill_games` join so a cross-game combo's full stake counts against each of its games, while the daily cap and P&L keep reading `fills` (one row per combo) and are never double-counted. The per-game and daily caps additionally count **open quotes' worst-case exposure** (issue #22): each `live_quotes` row freezes `worst_exposure_usd` at quote time, `quote_games` fans it out per game (same attribution rule as `fill_games`, written in the same transaction as the quote insert), and both cap gates sum `fills + open quotes` so a burst sweep of resting quotes can no longer fill to multiples of the caps before the first fill registers. A `NULL` `worst_exposure_usd` (pre-migration row) is counted at the per-fill cap, mirroring the N8 unreconciled-fill convention; the RFQ currently being re-quoted is excluded from the sums (its resting quote is superseded, not stacked). The sibling `kalshi_mlb_mm/kalshi_mlb_mm_market.duckdb` holds SGP-line and SGP-odds data (same pattern as the taker's `kalshi_mlb_rfq_market.duckdb`). The v1-hardening pass removed the model component of the blend, so there is no longer a read-only dependency on `Answer Keys/mlb_mm.duckdb`. All timestamp columns are `TIMESTAMPTZ` UTC (repo rule); existing naive-local DBs are migrated in place on startup by the idempotent `MIGRATE_SQL` (instants preserved), and the monitor normalizes reads so both frames render identically.

**Quote replace ordering.** When a re-priced quote moves beyond `QUOTE_HYSTERESIS`, the bot submits the NEW quote first and only then marks the old row `replaced` — Kalshi auto-cancels the prior quote only on a successful resubmit, so a failed submit leaves the old row `open` (still confirm-polled and risk-swept) and the next tick retries.

## Pricing

Fair value is the median of the devigged book fairs, gated on cross-book agreement. `router.combo_fair_detail` parses the RFQ into canonical legs (`legset`), partitions them by game, prices each game's sub-combo, and **multiplies** the per-game fairs (cross-game independence). For each **2-leg same-game grid** (family + cell resolved from the legs — spread×total keyed by game × spread_line × total_line, moneyline×total keyed by game × total_line with `spread_line` NULL) we:

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
- The stepped-fee/grid-floor guard is unchanged: the bid steps down until `bid + fee(bid) ≤ p − margin_pts`, so the **realized** margin is ≥ target by construction; the `yes_bid + no_bid < 1` and positivity guards still apply (an unquotable side declines the RFQ).
- Margin components (`sigma_pts`, `floor_pts`, `roi_pts_*`, `margin_pts_*`, `n_games`) are logged in the research firehose `quote_priced` payload — `quote_decisions` columns are unchanged.

Interplay of the two tickets: **moderate** dispersion widens the margin (#19); dispersion **past `SIGMA_Z_MAX`** kills the quote (#20). One uncertainty number drives both.

## On-demand pricing (Phase 2)

Spec: `docs/superpowers/specs/2026-07-08-kalshi-mlb-mm-on-demand-pricing-design.md` (rev 5).

Any same-game shape the two pre-scraped grids can't price is queried **live**
at the six books' SGP endpoints for that exact leg set, at RFQ time. Always
on — no config switch; the bot-wide kill file remains the emergency stop, and
dropping a misbehaving book is a one-line code change.

**Feed model — the open-RFQ poll drives everything.** The 2s discovery tick
never blocks on a book. An on-demand shape with no fresh result enqueues a
background fetch (`OnDemandEngine`, `kalshi_mlb_mm/on_demand.py`) and skips
with reason `on_demand_pending`; the fetch (~10–20s, all books concurrent)
lands per-book fairs in an in-memory store, and the next tick prices them
through the normal consensus + gates + hysteresis/replace path. A result may
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
2. **On-demand re-price (novel shapes only).** Combos with on-demand games
   additionally get the synchronous priority re-fetch (budgeted
   `CONFIRM_REFETCH_BUDGET_SEC=20s`) — those shapes have no cache fairs, so
   without a live fetch they cannot be re-priced for the EV/drift gate at
   all. The confirm gate honours `refetch_now`'s landed-after-entry guarantee
   (a failed re-fetch voids even if a pre-entry result still sits in the
   15s-fresh store) — we never confirm on a previously-fetched number.

**Failure direction is always "fewer/laggier quotes", never "staler
quotes"**: unresolvable leg → book drops; incomplete partition → Route B →
drop; sanity/Fréchet/consensus failure → no quote; worker death → lazy
restart, meanwhile skip; RFQ flood → queue grows, late landings, expired
RFQs go unquoted (pacing = one on-demand combo in flight per book).

**Observability**: research events `on_demand_requested` (once per fetch
flight), `on_demand_result` (once per landing; per-book route/fair/latency,
`route_gap` where both routes came free), fill payloads gain per-book
`route` + `consensus_books` (DK+Novig-only fills are a named risk metric).
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
| `FAIR_DRIFT_TOLERANCE` | `0.02` | Last-look: void confirm if fair drifted >2¢ against filled side since quote time |
| `MAX_BOOK_STALENESS_SEC` | `60` | Withhold and pull quotes if book odds older than this |
| `BOOK_MOVE_CB_THRESHOLD` | `0.03` | Circuit breaker: cancel a combo's quotes if book fair jumps >3¢ between scrapes (per-tick) or if drift since quote exceeds this (per-quote risk sweep) |
| `TIPOFF_CANCEL_MIN` | `5` | Pull quotes this many minutes before first pitch |
| `QUOTE_HYSTERESIS` | `0.005` | Don't replace a resting quote unless fair moved more than ½¢. Also the ε of the post-fill `same_price_block` |
| `COMBO_COOLDOWN_SEC` | `60` | Hard FLOOR of the post-fill per-combo cooldown. The combo additionally stays cooled until every game it touches has a post-fill SGP scrape, then until consensus fair moves off the filled fair (defense item 5) |
| `MAX_COMBO_EXPOSURE_USD` | `50.0` | Per-combo concentration cap (H8/N7): fills + in-flight open quotes on one combo may not exceed this |
| `SIGMA_Z_MAX` | `0.07` | Consensus gate (issue #20): max sample stddev of the books' fairs in z-space (`norm.ppf`); above it the combo is declined (`consensus_dispersion`) — no outlier removal. 0.07 ≈ continuity with the old ±2¢ band at p=0.50 |
| `MIN_AGREEING_BOOKS` | `2` | Consensus gate: minimum number of books with a fair for the combo (no longer "band survivors" — there is no band) |
| `DISCOVERY_SEC` | `2` | Discovery + quote loop cadence (seconds) |
| `CONFIRM_SEC` | `2` | Confirm loop cadence (seconds) |
| `RISK_SWEEP_SEC` | `10` | Risk sweep cadence (seconds) |
| `RECONCILE_SWEEP_SEC` | `30` | Fill side/size reconciliation cadence against Kalshi positions (seconds) |
| `SETTLEMENT_SWEEP_SEC` | `600` | Settlement sweep cadence (seconds) — populates `realized_pnl` once markets settle; only matters hours post-game |
| `SGP_REFRESH_SEC` | `60` | SGP scrape cadence (seconds) |
| `SGP_SCRAPER_TIMEOUT_SEC` | `90` | Per-book deadline passed to `SGPService` (seconds) — a book exceeding it contributes nothing that cycle and its client is rebuilt |

## Defense hierarchy (stale-quote / adverse-selection risk)

Resting quotes are priced off books that lag reality (books refresh every 60s). An informed counterparty can cross a stale quote before our data reacts. Defenses, in decreasing priority:

1. **Margin (primary, continuous coverage).** The per-side cushion `max(ROI part, MIN_MARGIN_PTS + K_SIGMA·σ_books)` absorbs *continuous* fair drift — a cent or two of movement between scrapes — and, since issue #19, keeps an absolute floor at longshot fairs (where the old flat-ROI cushion collapsed to ~0.3¢) and widens with cross-book disagreement and game count. This is the primary defense and fires on every fill. It does NOT cover discrete events (scratch / postponement / steam move).

2. **Book-consensus gate (issue #20 — z-space dispersion threshold).** Before quoting, we require `>= MIN_AGREEING_BOOKS` books to have priced the combo AND their fairs to agree within `SIGMA_Z_MAX` in probit space. No outlier removal: one loudly-dissenting book declines the whole combo — it may be the informed one, and quoting the stale majority's median is adverse selection by construction (a false decline is ~free; a false quote is not). This is v1's only correlation defense — the v1.1 explicit correlation-premium gate (where Kalshi singles serve as the marginal anchor) is documented in spec section 13 and deferred.

3. **Freshness gate + auto-pull (discrete events).** Before submitting any quote and inside the risk sweep, the bot checks that fresh book odds exist (`_SGP_ODDS` non-empty within `MAX_BOOK_STALENESS_SEC`). The instant books go stale or a scrape fails, all open quotes are cancelled. Blind → no live quotes.

4. **Book-move circuit breaker (discrete events).** Two layers: (a) per-tick — if a scrape shows a book-fair jump greater than `BOOK_MOVE_CB_THRESHOLD` for a combo vs the prior scrape, the bot immediately cancels that combo's resting quotes (does not wait for the next discovery tick). (b) per-quote in the risk sweep — if current book consensus has drifted more than `BOOK_MOVE_CB_THRESHOLD` from the `book_fair_at_quote` stored when the quote was placed, the quote is cancelled. The per-quote sweep catches gradual drift the per-tick threshold misses (e.g., five 1¢ moves across ticks).

5. **Post-fill cooldown + same-price block (issue #21).** After a fill, the combo is cooled for `COMBO_COOLDOWN_SEC` (hard floor) — and stays cooled past the floor until **every game the combo touches has a post-fill SGP scrape** (`in_cooldown_awaiting_refresh`): the 60s floor expires before the ~150–165s scrape cycle, so until fresh data lands a re-quote would re-post the exact books that priced the picked-off fill. Once books have refreshed, the combo is still skipped (`same_price_block`) while consensus fair remains within `QUOTE_HYSTERESIS` of the fair the fill transacted against (`fair_at_confirm`, falling back to `blended_fair_at_quote`) — quote prices are fair ± a fixed ROI margin, so unchanged fair ⇔ the identical just-picked-off price. Both skip reasons are logged distinctly in `quote_decisions` for the funnel report. The refresh check fails CLOSED (missing games / NaT / comparison errors keep the combo cooled).

6. **Tipoff blackout.** `TIPOFF_CANCEL_MIN` pulls all quotes for a game before first pitch. A cross-game combo is swept against the **earliest** first pitch across ALL its games (re-derived from `seen_rfqs.legs_json` each sweep — same legset path as discovery); any game the sweep can't resolve cancels the quote (fail-safe, never fail open).

7. **Last-look backstop (discrete events).** On accept, the bot first runs the **singles-move veto** (issue #17): the raw Kalshi odds of every leg, snapshotted at quote time, are re-read live (<~1s) and ANY one-tick move on any leg's bid or ask voids the fill (`voided_singles_moved`). Previously the last look re-priced grid combos from the same ~150s scrape cache the quote came from — `cur_fair == prev_fair` by construction, a no-op exactly in the fast-pickoff window. Combos that pass the veto still run the existing gates: (a) cannot re-price (no fresh books / failed on-demand re-fetch → `voided_no_fresh_books`), (b) no longer +EV (`price + fee >= current_fair`), (c) fair drifted past `FAIR_DRIFT_TOLERANCE` (→ `voided_last_look`). "Can't verify ⇒ don't confirm" is intentional throughout. Every accept emits a `confirm_singles_check` research event (quote age at accept, per-leg snapshot-vs-fresh odds, moved flag) — both the gate's save-rate and what a non-zero tolerance would have passed are measurable from day one. Non-confirms are abusive behavior Kalshi can throttle; the zero-tolerance veto's void rate is a watch item — if crowd noise voids too many clean accepts, the firehose deltas are the tuning dataset.

8. **Position reconciliation.** After every confirm we call `/portfolio/positions` for the combo ticker and trust Kalshi as the source of truth; if the confirm response's side or size disagrees, the `fills` row is written with the reconciled values and a `[position_mismatch]` warning is printed.

9. **Measure.** The `fills` table records `book_fair_at_quote`, `blended_fair_at_quote`, `fair_at_confirm`, and — via the **settlement sweep** (`settlement.py`, every `SETTLEMENT_SWEEP_SEC`; issue #12) — `realized_pnl` per fill: for each combo with reconciled unsettled fills the sweep polls `GET /markets/{ticker}` and, once Kalshi reports a yes/no result, writes `realized_pnl = contracts × ((won ? 1−price : −price) − fee)` plus a `settlements` audit row holding the raw market payload (the settlement response shape is unverified until real data lands). Only `reconciled=TRUE` fills settle, so a later side/size correction can never invalidate a written P&L; the sweep is idempotent and API failures just retry next pass. On the first settled ticker ever it also best-effort checks the maker-fee assumption against the `/portfolio` endpoints (`[FEE_MISMATCH]` logs loudly). Each settled fill emits a `settlement_recorded` research event, and `research_queries.sql` queries 7–8 decompose per-contract P&L into quoted margin vs fair drift (quote→confirm) vs settlement-vs-fair vs fee — the measurement-phase headline. If pickoffs swamp the margin, the honest conclusion is that making is not viable at this data latency → improve data speed (WebSocket feeds, faster scrape cadence) before scaling.

**Crash safety.** An unconfirmed quote that the process never confirms voids automatically — no surprise open position on crash. Graceful SIGTERM cancels all open quotes before exit. The main loop uses a 250ms sleep between sub-loop checks so SIGTERM is never starved (directly addresses the taker's known SIGTERM-starvation issue with its 640s RFQ refresh block).

## Open items to verify on first real fill

These are noted in spec §10 and remain unconfirmed until live fills happen:

1. **Maker fee** — confirm the actual charge matches 25% of taker fee (`maker_fee_per_contract` in `kalshi_common/ev_calc.py`). The settlement sweep automates a first pass: on the first settled ticker it fetches `/portfolio/fills` + `/portfolio/settlements`, compares any fee field found against our recorded fees (`[FEE_MISMATCH]` warning on disagreement), and captures the raw payloads in a `fee_verification` research event for manual inspection either way. Adjust the rounding in `ev_calc.py` if it differs.
2. **Quote-status polling shape** — exact fields on `GET /communications/quotes/{id}` (`status`, `accepted_side`, `contracts`). The confirm path in `main.py::_confirm_tick` infers `side_held` from `accepted_side` (`"no" if accepted_side == "yes" else "yes"`). Verify this field name and semantics on the first accepted quote.
3. **`get_competitors`** — competing quotes on an RFQ are stubbed out in v1 (future competitive pricing). Confirm the API shape if needed later.

## Accepted risks (v1)

Eight adversarial vectors are accepted for v1, measured not prevented. See spec §12 for the full table. The key ones: our quotes leak our fair surface (vector #1); fills are disproportionately from combos we underpriced (vector #2, the quiet killer — if per-combo fair error exceeds the margin, a patient sharp grinds us down regardless of other gates); the fav −1.5 + over family has documented model bias (vector #3). The validation dataset (`fills` + `realized_pnl` at settlement) is built specifically to measure fill-vs-fair error and answer "is 5% enough?"
