# Kalshi MLB MM (Maker) Bot

Independent maker daemon that listens for others' RFQs on the Kalshi cross-category MVE collection, prices **arbitrary MLB combos** (moneyline / spread / total, single- or cross-game, N legs) against a book-consensus fair value, and provides two-sided quotes at a fixed 5% ROI margin. Coexists with the taker (`kalshi_mlb_rfq/`) as a separate OS process with no runtime dependency on it.

**Coverage (2026-06 expansion).** The maker originally quoted only single-game 2-leg spread×total combos — which, when the live RFQ flow was measured, turned out to be ~6 of every ~75 all-MLB combos (most flow is cross-game moneyline parlays). Coverage now includes **moneyline (`KXMLBGAME`)**, single legs, and **cross-game / N-leg parlays** of any mix of spread/total/moneyline. Player props (`KXMLBHR` / `KXMLBHIT` / `KXMLBTB` / `KXMLBKS` / `KXMLBHRR`) are the remaining un-covered MLB leg type — they are correctly classified **out of scope** (a prop leg drops the whole combo, so we never mis-price one) pending a player-name matching layer. See "Pricing" below.

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
| SGP scrape | 60s | Own scraper cadence → `kalshi_mlb_mm_market.duckdb`: SGP spread×total grid (`mlb_sgp_odds`) **plus** moneyline/spread/total singles (`mlb_singles_odds`, one Odds API request, all books) |

**Transport seam.** All exchange I/O goes through two interfaces: `RFQSource.poll()` / `RFQSource.get_market()` (v1: `RestRFQSource`, REST poll of `GET /communications/rfqs?status=open`) and `QuoteGateway.submit_quote()` / `.confirm()` / `.cancel()` (v1: `RestQuoteGateway`). A WebSocket adapter is a drop-in replacement behind these interfaces — no pricing or risk code changes.

**Shared math.** Fair value, EV calc (including `maker_fee_per_contract`), authenticated HTTP, SGP orchestration, and leg-typing helpers all live in `kalshi_common/`. Both bots import from there; the taker's original files are now one-line re-export shims (behavior unchanged).

**SGP pricing (2026-06).** The bot prices SGPs **in-process** via `kalshi_common/sgp_service.py::SGPService` — no subprocess per cycle. The service holds persistent per-book HTTP clients reused across cycles (no per-cycle TLS handshake) and prices the four books concurrently under a per-book deadline (`SGP_SCRAPER_TIMEOUT_SEC`). DK/FD structure fetches (event lists, selection-id dicts) are TTL-cached; prices are never cached, and a failed or timed-out book keeps its prior rows. The old subprocess-per-cycle model is retained as a rollback hatch — calling `sgp_cycle` without `service=`.

**State DB.** The bot writes `kalshi_mlb_mm/kalshi_mlb_mm.duckdb` (quotes, fills, positions, decisions). The sibling `kalshi_mlb_mm/kalshi_mlb_mm_market.duckdb` holds SGP-line and SGP-odds data (same pattern as the taker's `kalshi_mlb_rfq_market.duckdb`). The v1-hardening pass removed the model component of the blend, so there is no longer a read-only dependency on `Answer Keys/mlb_mm.duckdb`.

## Pricing

`combo_pricer.combo_fair` prices an arbitrary MLB combo by **grouping legs by game**, pricing each game's sub-combo to ONE consensus fair, then **multiplying the independent games together** (games are statistically independent, so the parlay fair is the product of the per-game fairs):

1. **Group legs by game** (legs of the same game share an event-ticker suffix).
2. **Price each game group** to a single consensus fair:
   - **Single leg** (one moneyline / spread / total) → 2-way singles devig (`devig_two_way`) across the Odds API book universe in `mlb_singles_odds`.
   - **Same-game spread+total pair** → the 4-cell SGP-grid devig (`devig_book` over `mlb_sgp_odds`), which captures within-game correlation. The grid label is now derived from the *actual* leg sides (e.g. `Away Spread + Under`) — the previous code hardcoded `Home Spread + Over`, mis-pricing every other side combination.
   - **Any other same-game shape** (3-leg, spread+moneyline, …) → unpriceable → the whole combo is dropped. We never assume independence between correlated same-game legs.
3. **Per-group consensus gate** (v1 correlation defense): within each game's book universe, take the median, keep books within `±BOOK_CONSENSUS_BAND`, and require `>= MIN_AGREEING_BOOKS` survivors — else the group (and the combo) is not quoted.
4. **Combo fair = product of the per-group consensus fairs.** A combo is quotable only if *every* game group reached consensus. The combo's agreeing-book count is the minimum across groups (weakest link).

**Book data sources.** Same-game spread+total correlation comes from the maker's own SGP scrapers (`mlb_sgp_odds`, 4 books). Moneyline and single-leg spread/total fairs come from the Odds API h2h/spreads/totals markets (`mlb_singles_odds`, ~7–11 books per game), fetched once per SGP cycle. A cross-game moneyline parlay therefore prices entirely off the singles universe; a same-game spread+total off the SGP grid; mixes use whichever source each group needs.

This is the v1 correlation defense (mirrors the MLB answer-key dashboard's consensus-band pattern). The v1.1 explicit correlation-premium gate is deferred — see spec section 13.

**Multi-game caps (v1 simplification).** Per-game exposure and tipoff gates are checked across *all* games a combo spans (a cross-game combo is pulled if *any* of its games is within the tipoff blackout). Fill exposure attributes to the combo's first game — a documented measurement-phase simplification.

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
| `DAILY_EXPOSURE_CAP_PCT` | `0.75` | Daily hard stop as a fraction of BANKROLL ($375 at default) |
| `MAX_GAME_EXPOSURE_PCT` | `0.10` | Per-game exposure cap as fraction of BANKROLL ($50 at default) |
| `MAX_FILL_EXPOSURE_PCT` | `0.10` | Per-fill dollar cap as fraction of BANKROLL ($50 at default). Quote-or-skip only — the RFQ creator fixes fill size; this is the only lever. |
| `MAX_OPEN_QUOTES` | `25` | Cap on simultaneously resting quotes (well under Kalshi's 100 limit) |
| `TARGET_ROI` | `0.05` | Quoted margin — the `p / (1 + TARGET_ROI)` divisor in pricing |
| `FAIR_DRIFT_TOLERANCE` | `0.02` | Last-look: void confirm if fair drifted >2¢ against filled side since quote time |
| `MAX_BOOK_STALENESS_SEC` | `60` | Withhold and pull quotes if book odds older than this |
| `BOOK_MOVE_CB_THRESHOLD` | `0.03` | Circuit breaker: cancel a combo's quotes if book fair jumps >3¢ between scrapes (per-tick) or if drift since quote exceeds this (per-quote risk sweep) |
| `TIPOFF_CANCEL_MIN` | `5` | Pull quotes this many minutes before first pitch |
| `QUOTE_HYSTERESIS` | `0.005` | Don't replace a resting quote unless fair moved more than ½¢ |
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

5. **Tipoff blackout.** `TIPOFF_CANCEL_MIN` pulls all quotes for a game before first pitch.

6. **Last-look backstop (discrete, but limited scope).** On accept, the bot recomputes fair from a fresh book pull and voids the confirm if (a) we cannot re-price (no fresh books, blend fails), (b) the filled side is no longer +EV (`price + fee >= current_fair`), or (c) fair drifted past `FAIR_DRIFT_TOLERANCE`. The "can't re-price ⇒ don't confirm" rule is intentional: silently falling back to the stored fair would neuter the drift check. Non-confirms are abusive behavior Kalshi can throttle — do not lean on this gate.

7. **Position reconciliation.** After every confirm we call `/portfolio/positions` for the combo ticker and trust Kalshi as the source of truth; if the confirm response's side or size disagrees, the `fills` row is written with the reconciled values and a `[position_mismatch]` warning is printed.

8. **Measure.** The `fills` table records `book_fair_at_quote`, `blended_fair_at_quote`, `fair_at_confirm`, and (at settlement) `realized_pnl` per fill. The primary deliverable of v1 is computing whether the 5% margin survives the adverse-selection tail. If pickoffs swamp the margin, the honest conclusion is that making is not viable at this data latency → improve data speed (WebSocket feeds, faster scrape cadence) before scaling.

**Crash safety.** An unconfirmed quote that the process never confirms voids automatically — no surprise open position on crash. Graceful SIGTERM cancels all open quotes before exit. The main loop uses a 250ms sleep between sub-loop checks so SIGTERM is never starved (directly addresses the taker's known SIGTERM-starvation issue with its 640s RFQ refresh block).

## Open items to verify on first real fill

These are noted in spec §10 and remain unconfirmed until live fills happen:

1. **Maker fee** — confirm the actual charge matches 25% of taker fee (`maker_fee_per_contract` in `kalshi_common/ev_calc.py`). Adjust rounding there if it differs from the observed amount.
2. **Quote-status polling shape** — exact fields on `GET /communications/quotes/{id}` (`status`, `accepted_side`, `contracts`). The confirm path in `main.py::_confirm_tick` infers `side_held` from `accepted_side` (`"no" if accepted_side == "yes" else "yes"`). Verify this field name and semantics on the first accepted quote.
3. **`get_competitors`** — competing quotes on an RFQ are stubbed out in v1 (future competitive pricing). Confirm the API shape if needed later.

## Accepted risks (v1)

Eight adversarial vectors are accepted for v1, measured not prevented. See spec §12 for the full table. The key ones: our quotes leak our fair surface (vector #1); fills are disproportionately from combos we underpriced (vector #2, the quiet killer — if per-combo fair error exceeds the margin, a patient sharp grinds us down regardless of other gates); the fav −1.5 + over family has documented model bias (vector #3). The validation dataset (`fills` + `realized_pnl` at settlement) is built specifically to measure fill-vs-fair error and answer "is 5% enough?"
