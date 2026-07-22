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

REST-polling daemon, single process. Four timed sub-loops:

| Loop | Cadence | Job |
|---|---|---|
| Discovery + quote | 2s | Poll open RFQs → scope filter → price → submit or refresh quote |
| Confirm | 2s | Poll open quotes → on `accepted`, last-look gate → confirm or void |
| Risk sweep | 10s | Kill-switch, book-staleness auto-pull, tipoff cancel, drift-since-quote cancel |
| SGP scrape | 60s | Own scraper cadence → `kalshi_mlb_mm_market.duckdb` (sibling market DB) |

**Transport seam.** All exchange I/O goes through two interfaces: `RFQSource.poll()` / `RFQSource.get_market()` (v1: `RestRFQSource`, REST poll of `GET /communications/rfqs?status=open`) and `QuoteGateway.submit_quote()` / `.confirm()` / `.cancel()` (v1: `RestQuoteGateway`). A WebSocket adapter is a drop-in replacement behind these interfaces — no pricing or risk code changes.

**Shared math.** Fair value, EV calc (including `maker_fee_per_contract`), authenticated HTTP, SGP orchestration, and leg-typing helpers all live in `kalshi_common/`. Both bots import from there; the taker's original files are now one-line re-export shims (behavior unchanged).

**SGP pricing (2026-06).** The bot prices SGPs **in-process** via `kalshi_common/sgp_service.py::SGPService` — no subprocess per cycle. The service holds persistent per-book HTTP clients reused across cycles (no per-cycle TLS handshake) and prices the four books concurrently under a per-book deadline (`SGP_SCRAPER_TIMEOUT_SEC`). DK/FD structure fetches (event lists, selection-id dicts) are TTL-cached; prices are never cached, and a failed or timed-out book keeps its prior rows. The old subprocess-per-cycle model is retained as a rollback hatch — calling `sgp_cycle` without `service=`.

**State DB.** The bot writes `kalshi_mlb_mm/kalshi_mlb_mm.duckdb` (quotes, fills, positions, decisions, `fill_games`, and `quote_games`). The `fill_games` table maps each fill to **every** game its combo touches (one row per game); the per-game exposure cap reads a fanned-out `fills⋈fill_games` join so a cross-game combo's full stake counts against each of its games, while the daily cap and P&L keep reading `fills` (one row per combo) and are never double-counted. The per-game and daily caps additionally count **open quotes' worst-case exposure** (issue #22): each `live_quotes` row freezes `worst_exposure_usd` at quote time, `quote_games` fans it out per game (same attribution rule as `fill_games`, written in the same transaction as the quote insert), and both cap gates sum `fills + open quotes` so a burst sweep of resting quotes can no longer fill to multiples of the caps before the first fill registers. A `NULL` `worst_exposure_usd` (pre-migration row) is counted at the per-fill cap, mirroring the N8 unreconciled-fill convention; the RFQ currently being re-quoted is excluded from the sums (its resting quote is superseded, not stacked). The sibling `kalshi_mlb_mm/kalshi_mlb_mm_market.duckdb` holds SGP-line and SGP-odds data (same pattern as the taker's `kalshi_mlb_rfq_market.duckdb`). The v1-hardening pass removed the model component of the blend, so there is no longer a read-only dependency on `Answer Keys/mlb_mm.duckdb`. All timestamp columns are `TIMESTAMPTZ` UTC (repo rule); existing naive-local DBs are migrated in place on startup by the idempotent `MIGRATE_SQL` (instants preserved), and the monitor normalizes reads so both frames render identically.

**Quote replace ordering.** When a re-priced quote moves beyond `QUOTE_HYSTERESIS`, the bot submits the NEW quote first and only then marks the old row `replaced` — Kalshi auto-cancels the prior quote only on a successful resubmit, so a failed submit leaves the old row `open` (still confirm-polled and risk-swept) and the next tick retries.

## Pricing

Fair value is the median of *book-consensus-agreeing* devigged book fairs. `router.combo_fair` parses the RFQ into canonical legs (`legset`), partitions them by game, prices each game's sub-combo, and **multiplies** the per-game fairs (cross-game independence). For each **2-leg same-game grid** (family + cell resolved from the legs — spread×total keyed by game × spread_line × total_line, moneyline×total keyed by game × total_line with `spread_line` NULL) we:

1. Pull every book's 4-cell grid from `mlb_sgp_odds` (filtered by combo family, require all 4 cells — no fallback).
2. Devig each book's 4-way grid to a single combo fair (`devig_book` in `kalshi_common.fair_value`).
3. Compute the median across books, then keep only books within `±BOOK_CONSENSUS_BAND` of that median.
4. If `>= MIN_AGREEING_BOOKS` survive, the fair is the median of the survivors. Otherwise we do not quote.

A **single leg** (within a cross-game combo) is priced by marginalizing that game's grid — summing the two devigged cells along the free axis. **Cross-game** fairs are the product of the per-game consensus fairs; if any game fails consensus or is `unpriceable`, `combo_fair` returns `None` and we do not quote (fail-safe). Novel same-game shapes (`on_demand` — 3+ legs, spread+ml, total+total, …) price via live book queries — see **On-demand pricing** below. This is the v1 correlation defense (mirrors the MLB answer-key dashboard's consensus-band pattern). The v1.1 explicit correlation-premium gate is deferred — see spec section 13.

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

**Accepted quotes get a live last look**: the confirm tick synchronously
re-fetches every on-demand game **and every 2-leg grid game** in the combo
(issue #17 — a grid sub-combo is just a 2²-cell partition to
`price_on_demand`, so grids re-price on a live fetch instead of the ~150s
scrape cache; priority lane, all jobs share the
`CONFIRM_REFETCH_BUDGET_SEC=20s` budget against the ~30s confirm window).
A failed or late re-fetch leaves the lookup stale → `voided_no_fresh_books`.
We never confirm on a previously-fetched number. Single-leg sub-combos of
cross-game combos still last-look from the cache marginal — issue #23 adds
a Kalshi-singles fast check at that seam.

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

## Knobs

All knobs are overridable via `kalshi_mlb_mm/.env` or environment variables. Defaults come from `kalshi_mlb_mm/config.py`.

| Knob | Default | Purpose |
|---|---|---|
| `BANKROLL` | `500.0` | Master risk dial — raise this one number to scale everything |
| `DAILY_EXPOSURE_CAP_PCT` | `0.75` | Daily hard stop as a fraction of BANKROLL ($375 at default). Counts today's fills PLUS open quotes' worst-case exposure |
| `MAX_GAME_EXPOSURE_PCT` | `0.10` | Per-game exposure cap as fraction of BANKROLL ($50 at default). Counts fills (`fill_games`) PLUS open quotes (`quote_games`) touching the game |
| `MAX_FILL_EXPOSURE_PCT` | `0.10` | Per-fill dollar cap as fraction of BANKROLL ($50 at default). Quote-or-skip only — the RFQ creator fixes fill size; this is the only lever. |
| `MAX_OPEN_QUOTES` | `25` | Cap on simultaneously resting quotes (well under Kalshi's 100 limit) |
| `TARGET_ROI` | `0.05` | Quoted margin — the `p / (1 + TARGET_ROI)` divisor in pricing |
| `FAIR_DRIFT_TOLERANCE` | `0.02` | Last-look: void confirm if fair drifted >2¢ against filled side since quote time |
| `MAX_BOOK_STALENESS_SEC` | `60` | Withhold and pull quotes if book odds older than this |
| `BOOK_MOVE_CB_THRESHOLD` | `0.03` | Circuit breaker: cancel a combo's quotes if book fair jumps >3¢ between scrapes (per-tick) or if drift since quote exceeds this (per-quote risk sweep) |
| `TIPOFF_CANCEL_MIN` | `5` | Pull quotes this many minutes before first pitch |
| `QUOTE_HYSTERESIS` | `0.005` | Don't replace a resting quote unless fair moved more than ½¢. Also the ε of the post-fill `same_price_block` |
| `COMBO_COOLDOWN_SEC` | `60` | Hard FLOOR of the post-fill per-combo cooldown. The combo additionally stays cooled until every game it touches has a post-fill SGP scrape, then until consensus fair moves off the filled fair (defense item 5) |
| `MAX_COMBO_EXPOSURE_USD` | `50.0` | Per-combo concentration cap (H8/N7): fills + in-flight open quotes on one combo may not exceed this |
| `BOOK_CONSENSUS_BAND` | `0.02` | v1 correlation defense: max distance from per-combo book median (fair-prob units) for a book to count as "agreeing"; outliers are discarded |
| `MIN_AGREEING_BOOKS` | `3` | v1 correlation defense: minimum number of books that must agree before we quote |
| `DISCOVERY_SEC` | `2` | Discovery + quote loop cadence (seconds) |
| `CONFIRM_SEC` | `2` | Confirm loop cadence (seconds) |
| `RISK_SWEEP_SEC` | `10` | Risk sweep cadence (seconds) |
| `SGP_REFRESH_SEC` | `60` | SGP scrape cadence (seconds) |
| `SGP_SCRAPER_TIMEOUT_SEC` | `90` | Per-book deadline passed to `SGPService` (seconds) — a book exceeding it contributes nothing that cycle and its client is rebuilt |

## Defense hierarchy (stale-quote / adverse-selection risk)

Resting quotes are priced off books that lag reality (books refresh every 60s). An informed counterparty can cross a stale quote before our data reacts. Defenses, in decreasing priority:

1. **Margin (primary, continuous coverage).** The ~2.5–3¢ per-side cushion at 5% ROI absorbs *continuous* fair drift — a cent or two of movement between scrapes. This is the primary defense and fires on every fill. It does NOT cover discrete events (scratch / postponement / steam move).

2. **Book-consensus gate (correlation defense).** Before quoting, we require `>= MIN_AGREEING_BOOKS` books within `±BOOK_CONSENSUS_BAND` of the per-combo book median; outliers are discarded and the fair is the median of survivors. A single rogue book cannot anchor our quote. This is v1's only correlation defense — the v1.1 explicit correlation-premium gate (where Kalshi singles serve as the marginal anchor) is documented in spec section 13 and deferred.

3. **Freshness gate + auto-pull (discrete events).** Before submitting any quote and inside the risk sweep, the bot checks that fresh book odds exist (`_SGP_ODDS` non-empty within `MAX_BOOK_STALENESS_SEC`). The instant books go stale or a scrape fails, all open quotes are cancelled. Blind → no live quotes.

4. **Book-move circuit breaker (discrete events).** Two layers: (a) per-tick — if a scrape shows a book-fair jump greater than `BOOK_MOVE_CB_THRESHOLD` for a combo vs the prior scrape, the bot immediately cancels that combo's resting quotes (does not wait for the next discovery tick). (b) per-quote in the risk sweep — if current book consensus has drifted more than `BOOK_MOVE_CB_THRESHOLD` from the `book_fair_at_quote` stored when the quote was placed, the quote is cancelled. The per-quote sweep catches gradual drift the per-tick threshold misses (e.g., five 1¢ moves across ticks).

5. **Post-fill cooldown + same-price block (issue #21).** After a fill, the combo is cooled for `COMBO_COOLDOWN_SEC` (hard floor) — and stays cooled past the floor until **every game the combo touches has a post-fill SGP scrape** (`in_cooldown_awaiting_refresh`): the 60s floor expires before the ~150–165s scrape cycle, so until fresh data lands a re-quote would re-post the exact books that priced the picked-off fill. Once books have refreshed, the combo is still skipped (`same_price_block`) while consensus fair remains within `QUOTE_HYSTERESIS` of the fair the fill transacted against (`fair_at_confirm`, falling back to `blended_fair_at_quote`) — quote prices are fair ± a fixed ROI margin, so unchanged fair ⇔ the identical just-picked-off price. Both skip reasons are logged distinctly in `quote_decisions` for the funnel report. The refresh check fails CLOSED (missing games / NaT / comparison errors keep the combo cooled).

6. **Tipoff blackout.** `TIPOFF_CANCEL_MIN` pulls all quotes for a game before first pitch. A cross-game combo is swept against the **earliest** first pitch across ALL its games (re-derived from `seen_rfqs.legs_json` each sweep — same legset path as discovery); any game the sweep can't resolve cancels the quote (fail-safe, never fail open).

7. **Last-look backstop (discrete events).** On accept, the bot re-prices the combo on a **live book fetch** — on-demand games AND 2-leg grid games (issue #17; grids previously re-read the ~150s scrape cache, making the drift check a no-op exactly in the fast-pickoff window) — and voids the confirm if (a) we cannot re-price live (fetch failed/late, too few books → `voided_no_fresh_books`), (b) the filled side is no longer +EV (`price + fee >= current_fair`), or (c) fair drifted past `FAIR_DRIFT_TOLERANCE` (→ `voided_last_look`). The "can't re-price ⇒ don't confirm" rule is intentional: silently falling back to the stored fair or the scrape cache would neuter the drift check. `fills.fair_at_confirm` stores the fresh fair, and every accept emits a `confirm_fresh_check` research event (quote age at accept, per-book fetch latency, stale-cache vs fresh fair delta) so the gate's save-rate is measurable. Non-confirms are abusive behavior Kalshi can throttle — do not lean on this gate.

8. **Position reconciliation.** After every confirm we call `/portfolio/positions` for the combo ticker and trust Kalshi as the source of truth; if the confirm response's side or size disagrees, the `fills` row is written with the reconciled values and a `[position_mismatch]` warning is printed.

9. **Measure.** The `fills` table records `book_fair_at_quote`, `blended_fair_at_quote`, `fair_at_confirm`, and (at settlement) `realized_pnl` per fill. The primary deliverable of v1 is computing whether the 5% margin survives the adverse-selection tail. If pickoffs swamp the margin, the honest conclusion is that making is not viable at this data latency → improve data speed (WebSocket feeds, faster scrape cadence) before scaling.

**Crash safety.** An unconfirmed quote that the process never confirms voids automatically — no surprise open position on crash. Graceful SIGTERM cancels all open quotes before exit. The main loop uses a 250ms sleep between sub-loop checks so SIGTERM is never starved (directly addresses the taker's known SIGTERM-starvation issue with its 640s RFQ refresh block).

## Open items to verify on first real fill

These are noted in spec §10 and remain unconfirmed until live fills happen:

1. **Maker fee** — confirm the actual charge matches 25% of taker fee (`maker_fee_per_contract` in `kalshi_common/ev_calc.py`). Adjust rounding there if it differs from the observed amount.
2. **Quote-status polling shape** — exact fields on `GET /communications/quotes/{id}` (`status`, `accepted_side`, `contracts`). The confirm path in `main.py::_confirm_tick` infers `side_held` from `accepted_side` (`"no" if accepted_side == "yes" else "yes"`). Verify this field name and semantics on the first accepted quote.
3. **`get_competitors`** — competing quotes on an RFQ are stubbed out in v1 (future competitive pricing). Confirm the API shape if needed later.

## Accepted risks (v1)

Eight adversarial vectors are accepted for v1, measured not prevented. See spec §12 for the full table. The key ones: our quotes leak our fair surface (vector #1); fills are disproportionately from combos we underpriced (vector #2, the quiet killer — if per-combo fair error exceeds the margin, a patient sharp grinds us down regardless of other gates); the fav −1.5 + over family has documented model bias (vector #3). The validation dataset (`fills` + `realized_pnl` at settlement) is built specifically to measure fill-vs-fair error and answer "is 5% enough?"
