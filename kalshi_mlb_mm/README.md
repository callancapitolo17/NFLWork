# Kalshi MLB MM (Maker) Bot

Independent maker daemon that listens for others' RFQs on the Kalshi cross-category MVE collection, prices 2-leg spread×total MLB combos using the same model+book blended fair value as the taker, and provides two-sided quotes at a fixed 5% ROI margin. Coexists with the taker (`kalshi_mlb_rfq/`) as a separate OS process with no runtime dependency on it.

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

REST-polling daemon, single process. Five timed sub-loops:

| Loop | Cadence | Job |
|---|---|---|
| Discovery + quote | 2s | Poll open RFQs → scope filter → price → submit or refresh quote |
| Confirm | 2s | Poll open quotes → on `accepted`, last-look gate → confirm or void |
| Risk sweep | 10s | Kill-switch, samples-staleness auto-pull, tipoff cancel |
| SGP scrape | 60s | Own scraper cadence → `kalshi_mlb_mm_market.duckdb` (sibling market DB) |
| Samples refresh | 600s | Reload `Answer Keys/mlb_mm.duckdb::mlb_game_samples` read-only |

**Transport seam.** All exchange I/O goes through two interfaces: `RFQSource.poll()` / `RFQSource.get_market()` (v1: `RestRFQSource`, REST poll of `GET /communications/rfqs?status=open`) and `QuoteGateway.submit_quote()` / `.confirm()` / `.cancel()` (v1: `RestQuoteGateway`). A WebSocket adapter is a drop-in replacement behind these interfaces — no pricing or risk code changes.

**Shared math.** Fair value, EV calc (including `maker_fee_per_contract`), authenticated HTTP, SGP orchestration, and leg-typing helpers all live in `kalshi_common/`. Both bots import from there; the taker's original files are now one-line re-export shims (behavior unchanged).

**State DB.** The bot writes `kalshi_mlb_mm/kalshi_mlb_mm.duckdb` (quotes, fills, positions, decisions). The sibling `kalshi_mlb_mm/kalshi_mlb_mm_market.duckdb` holds SGP-line and SGP-odds data (same pattern as the taker's `kalshi_mlb_rfq_market.duckdb`). The answer-key DB (`Answer Keys/mlb_mm.duckdb`) is read-only.

## Pricing

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
| `MAX_RFQ_CONTRACTS` | `5` | Decline any RFQ requesting more contracts than this (~$5/fill max) |
| `MAX_OPEN_QUOTES` | `25` | Cap on simultaneously resting quotes (well under Kalshi's 100 limit) |
| `TARGET_ROI` | `0.05` | Quoted margin — the `p / (1 + TARGET_ROI)` divisor in pricing |
| `FAIR_DRIFT_TOLERANCE` | `0.02` | Last-look: void confirm if fair drifted >2¢ against filled side since quote time |
| `MAX_PREDICTION_STALENESS_SEC` | `600` | Withhold and pull all quotes if model samples older than this |
| `MAX_BOOK_STALENESS_SEC` | `60` | Withhold and pull quotes if book odds older than this |
| `BOOK_MOVE_CB_THRESHOLD` | `0.03` | Circuit breaker: cancel a game's quotes if book fair jumps >3¢ between scrapes |
| `TIPOFF_CANCEL_MIN` | `5` | Pull quotes this many minutes before first pitch |
| `QUOTE_HYSTERESIS` | `0.005` | Don't replace a resting quote unless fair moved more than ½¢ |
| `DISCOVERY_SEC` | `2` | Discovery + quote loop cadence (seconds) |
| `CONFIRM_SEC` | `2` | Confirm loop cadence (seconds) |
| `RISK_SWEEP_SEC` | `10` | Risk sweep cadence (seconds) |
| `SGP_REFRESH_SEC` | `60` | SGP scrape cadence (seconds) |
| `SGP_SCRAPER_TIMEOUT_SEC` | `90` | Per-scraper kill deadline (seconds) |
| `SAMPLES_REFRESH_SEC` | `600` | Model samples reload cadence (seconds) |
| `MIN_BOOK_COUNT_FOR_BLEND` | `2` | Minimum number of books required to produce a blended fair |

## Defense hierarchy (stale-quote / adverse-selection risk)

Resting quotes are priced off inputs that lag reality — books refresh every 60s, model samples every 600s. An informed counterparty can cross a stale quote before our data reacts. Defenses, in decreasing priority:

1. **Margin (primary, continuous coverage).** The ~2.5–3¢ per-side cushion at 5% ROI absorbs *continuous* fair drift — a cent or two of movement between scrapes. This is the primary defense and fires on every fill. It does NOT cover discrete events (scratch / postponement / steam move).

2. **Freshness gate + auto-pull (discrete events).** Before submitting any quote and inside the risk sweep, the bot checks that both model samples and book odds are within their staleness thresholds. The instant either goes stale or a scrape fails, all open quotes are cancelled. Blind → no live quotes.

3. **Book-move circuit breaker (discrete events).** If a scrape shows a book-fair jump greater than `BOOK_MOVE_CB_THRESHOLD` for a game vs the prior scrape, the bot immediately cancels that game's quotes — it does not wait for the next 2s discovery tick. A sudden move is the signal that unmodeled news (scratch, postponement, sharp steam) just landed.

4. **Tipoff blackout.** `TIPOFF_CANCEL_MIN` pulls all quotes for a game before first pitch. Lineup- and weather-driven moves are caught indirectly by the circuit breaker once they hit the books.

5. **Last-look backstop (discrete, but limited scope).** On accept, the bot recomputes fair and voids the confirm if the filled side is no longer +EV (`price + fee >= current_fair`) or fair drifted past `FAIR_DRIFT_TOLERANCE`. This is an explicit backstop, not the first line of defense: it only catches information the pipeline already absorbed but hadn't yet caused a re-quote. It is blind to faster-than-our-data pickoffs. Non-confirms are abusive behavior Kalshi can throttle — do not lean on this gate.

6. **Measure.** The `fills` + `settlements` tables record `model_fair_at_quote`, `book_fair_at_quote`, `blended_fair_at_quote`, `fair_at_confirm`, and `realized_pnl` per fill. The primary deliverable of v1 is computing whether the 5% margin survives the adverse-selection tail. If pickoffs swamp the margin, the honest conclusion is that making is not viable at this data latency → improve data speed (WebSocket feeds, faster scrape cadence) before scaling.

**Crash safety.** An unconfirmed quote that the process never confirms voids automatically — no surprise open position on crash. Graceful SIGTERM cancels all open quotes before exit. The main loop uses a 250ms sleep between sub-loop checks so SIGTERM is never starved (directly addresses the taker's known SIGTERM-starvation issue with its 640s RFQ refresh block).

## Open items to verify on first real fill

These are noted in spec §10 and remain unconfirmed until live fills happen:

1. **Maker fee** — confirm the actual charge matches 25% of taker fee (`maker_fee_per_contract` in `kalshi_common/ev_calc.py`). Adjust rounding there if it differs from the observed amount.
2. **Quote-status polling shape** — exact fields on `GET /communications/quotes/{id}` (`status`, `accepted_side`, `contracts`). The confirm path in `main.py::_confirm_tick` infers `side_held` from `accepted_side` (`"no" if accepted_side == "yes" else "yes"`). Verify this field name and semantics on the first accepted quote.
3. **`get_competitors`** — competing quotes on an RFQ are stubbed out in v1 (future competitive pricing). Confirm the API shape if needed later.

## Accepted risks (v1)

Eight adversarial vectors are accepted for v1, measured not prevented. See spec §12 for the full table. The key ones: our quotes leak our fair surface (vector #1); fills are disproportionately from combos we underpriced (vector #2, the quiet killer — if per-combo fair error exceeds the margin, a patient sharp grinds us down regardless of other gates); the fav −1.5 + over family has documented model bias (vector #3). The validation dataset (`fills` + `realized_pnl` at settlement) is built specifically to measure fill-vs-fair error and answer "is 5% enough?"
