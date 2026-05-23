# Kalshi MLB RFQ Bot

Autonomous taker daemon that auto-RFQs MLB game-line same-game-parlays on Kalshi (cross-category MVE collection), evaluates maker quotes against a model + sportsbook blended fair value, and auto-accepts +EV quotes within a complete safety scaffold.

**Spec:** `docs/superpowers/specs/2026-04-27-kalshi-mlb-rfq-bot-design.md`
**Plan:** `docs/superpowers/plans/2026-04-29-kalshi-mlb-rfq-bot.md`

## Quick start

```bash
cd ~/NFLWork/kalshi_mlb_rfq

# 1. Copy env template, fill credentials
cp .env.example .env
# Edit .env: set KALSHI_API_KEY_ID, KALSHI_PRIVATE_KEY_PATH, KALSHI_USER_ID

# 2. Install deps
python3 -m venv venv
./venv/bin/pip install -r requirements.txt

# 3. Make sure the MLB pipeline has run at least once
cd "../Answer Keys" && python run.py --sport mlb && cd ../kalshi_mlb_rfq

# 4. Dry-run validation (foreground, no real accepts)
./venv/bin/python -m kalshi_mlb_rfq.main --dry-run

# 5. Live mode (backgrounded)
./venv/bin/python -u -m kalshi_mlb_rfq.main >> bot.log 2>&1 &
tail -f bot.log
```

## Stopping the bot

```bash
# Graceful (cancels all live RFQs first via SIGTERM)
kill $(pgrep -f "kalshi_mlb_rfq.main")

# Emergency kill switch (file-based; bot stays alive but stops accepting)
touch .kill
# Resume:
rm .kill
```

## Status

```bash
./venv/bin/python -m kalshi_mlb_rfq.dashboard
```

## Architecture

See spec §3. TL;DR: standalone daemon, reads `Answer Keys/mlb_mm.duckdb` (samples + sgp_odds) read-only, writes `kalshi_mlb_rfq.duckdb`. Continuous priority-queue pipeline of up to 80 in-flight RFQs.

## Devigging

`fair_value.devig_book` uses probit (additive z-shift) devigging on the 4-cell SGP joint distribution (Home/Away Spread × Over/Under). When fewer than 4 sides are visible (interpolated combos with partial coverage), falls back to `(1/decimal) / (1 + vig_fallback)` heuristic. See `docs/superpowers/specs/2026-05-11-probit-devig-design.md`.

Requires `scipy>=1.10` (added to `requirements.txt`). Run `pip install -r requirements.txt` after pulling this change.

## Data Sources

The bot reads `Answer Keys/mlb_mm.duckdb` (read-only). This file holds the six consumer-facing tables and is written by `MLB.R`, `mlb_correlated_parlay.R`, `mlb_triple_play.R`, and the SGP scrapers — separate from `mlb.duckdb` (the pipeline's main write target) so lock contention can't block the bot's cache refreshes.

Four internal loops:
- **RFQ refresh @ 30s** — enumerate, score, prioritize, submit/cancel RFQs
- **Quote poll @ 2s** — fetch quotes, evaluate gates, accept under lock
- **Risk sweep @ 10s** — tipoff cancels, kill-switch check
- **Pipeline refresh @ 600s** — re-run answer key

## Sizing (two-RFQ per-side Kelly, 2026-05-23)

The bot now sizes each combo with **two RFQs**, one per side (YES/NO), each
sized via `target_cost_dollars` for that side's Kelly count at its
worst-acceptable price. Replaces the prior hardcoded `target_cost_dollars=$1`.

**Why two RFQs.** Kalshi's RFQ has one size knob per request, and REST
accept is all-or-nothing. With a single RFQ you must either under-size
the dominant side (min budget) or over-bet the opposite side on rare
mispricings (max budget). Two RFQs let each side carry its own correctly
Kelly-sized budget. The maker quotes each RFQ independently; the bot
accepts the one whose +EV side matches the RFQ's `intended_side` and
declines the mismatch.

**How sizing is computed.** Per candidate, the bot computes:
1. `worst_yes_ask`, `worst_no_ask` — binary search via `ev_calc` for the
   highest price on each side that still meets `KELLY_CREATE_EV_FLOOR_PCT`
   (defaults to `MIN_EV_PCT=0.05`) after fees.
2. `kelly_yes_n`, `kelly_no_n` — Kelly's contract count at each worst
   price, using the side-appropriate fair (`fair` for YES, `1-fair` for
   NO) and existing-position correlation on the same game.
3. `target_cost_dollars = kelly_n × worst_ask` per side. Skip the side if
   either is zero (no acceptable price after fees, or Kelly says don't
   bet).

**Per-RFQ accept gate.** Each RFQ is stored in `live_rfqs` with an
`intended_side` column. At quote-evaluation time, the bot computes
`chosen_side` (the +EV side per math invariant) and declines if it
doesn't match `intended_side` (`decision='declined_side_mismatch'`). NULL
`intended_side` on legacy rows is treated as "no constraint" for
backward compatibility.

**Audit columns in `live_rfqs`:**
- `intended_side` — "yes" or "no" (or NULL for legacy rows)
- `kelly_yes_n_at_submit`, `kelly_no_n_at_submit`
- `worst_yes_ask_at_submit`, `worst_no_ask_at_submit`
- `target_cost_dollars_at_submit` — the actual budget sent to Kalshi

## Knobs

See `.env.example`. Most relevant:
- `BANKROLL`, `KELLY_FRACTION` — sizing inputs
- `MIN_EV_PCT` — accept threshold (also default for sizing's EV floor)
- `KELLY_CREATE_EV_FLOOR_PCT` — EV floor used to find worst-acceptable
  prices at create time. Defaults to `MIN_EV_PCT`. Raise to size more
  conservatively, lower to be more aggressive
- `MAX_GAME_EXPOSURE_PCT`, `DAILY_EXPOSURE_CAP_USD` — exposure caps
- `MAX_PREDICTION_STALENESS_SEC` — accept-gate staleness threshold (10 min default)
- `MIN_FILL_RATIO`, `FILL_RATIO_WINDOW` — adverse-selection halt

## Safety scaffold (per-accept gates)

Every quote acceptance must pass ALL of:
- Min EV after fee · Fair-value bounds · Prediction staleness
- Tipoff window (5 min before first pitch) · Line-move check · Per-game exposure cap (% of bankroll)
- Daily exposure cap · Kill-switch off · Inverse-combo not held · 2-source gate (model + ≥1 book)
- Per-combo cooldown (30s post-fill) · Positions API health · Fill-ratio halt (rolling 50 attempts)

Acceptance is serialized via `ACCEPT_LOCK` so concurrent quotes Kelly-size against fresh DB state.

## Troubleshooting

- **No combos surfacing:** check that `mlb_game_samples` and `mlb_sgp_odds` are populated. Run the MLB pipeline first.
- **All quotes declining `declined_stale_predictions`:** rerun the pipeline; samples are over 10 minutes old.
- **Bot halted on `fill_ratio_collapse`:** investigate — makers are walking on accepts at a rate that suggests adverse selection. Check `quote_log` for the maker `creator_id`s causing it.
- **`mint_combo_ticker` failing with 400:** the MVE collection ticker may have changed or one of the leg market_tickers doesn't exist. Re-run `mlb_sgp/recon_kalshi_mlb_rfq.py` for a fresh probe.

## Accept semantics (symmetric NO-side support, 2026-05-22)

Kalshi's REST `PUT /communications/quotes/{id}/accept` is **all-or-nothing** — body is `{"accepted_side": "yes"|"no"}` and Kalshi returns 204 on success. There is no per-accept contract count parameter.

Sizing happens **upstream at RFQ creation** via `target_cost_dollars` (currently flat $1 diagnostic). With `rest_remainder: False`, the LP's quote is fill-or-kill at the requested size — accepting fills exactly what was requested. Real-size sizing per-contract awaits the FIX migration (separate plan).

Historical bug fixed 2026-05-21: from first commit (2026-05-02) through 2026-05-21, this code sent body `{"contracts": N}`, which Kalshi rejected as `invalid_parameters` on every accept. The bug was masked because the failure looked identical to "lost the race." The walk-diagnostic columns added 2026-05-20 finally exposed it — every walked row showed `accept_response_body = "invalid_parameters"` AND `rfq_terminal_status = "open"` (no competitor had filled the RFQ; we just couldn't accept it).

### Side semantics — `accepted_side` is INVERTED from intuition

The `accepted_side` field names the side of the LP's two-sided quote we're accepting — and accepting their bid on a side means *they buy that side from us*, leaving us LONG the OPPOSITE side. Verified empirically 2026-05-21 from the first real fill:

| Field value | What Kalshi does | What we end up holding | Bot sends this when |
|---|---|---|---|
| `accepted_side="yes"` | LP buys YES from us at `yes_bid` | **LONG NO** at `1 − yes_bid` (≈ `no_ask`) | `chosen_side="no"` (we want to buy NO) |
| `accepted_side="no"` | LP buys NO from us at `no_bid` | **LONG YES** at `1 − no_bid` (≈ `yes_ask`) | `chosen_side="yes"` (we want to buy YES) |

### Side selection — symmetric EV evaluation

Per quote, the bot computes EV on both sides after fees:

```
ev_yes = post_fee_ev_buy_yes(fair, no_bid)   # pay yes_ask = 1 − no_bid, payout = fair
ev_no  = post_fee_ev_buy_no (fair, yes_bid)  # pay no_ask  = 1 − yes_bid, payout = 1 − fair
```

It then routes the quote based on which side(s) clear `MIN_EV_PCT`:

| Case | Action |
|---|---|
| Only YES eligible | `chosen_side="yes"` → accept with `accepted_side="no"` |
| Only NO  eligible | `chosen_side="no"`  → accept with `accepted_side="yes"` |
| Neither eligible  | `decline_ev` |
| **Both eligible** | `decline_math_invariant` + alert (impossible after fees IRL; see below) |

**Why both-sides-+EV is a guard, not a feature.** For any LP making money on a two-sided quote, `yes_bid + no_bid < 1`, which forces `yes_ask + no_ask > 1`. Adding fees: `yes_ask + no_ask + fees > 1`. The both-+EV condition requires `yes_ask + no_ask + fees < 1` — contradiction. If the guard ever fires it means the fee model or fair-value pipeline has drifted; the bot logs `MATH_INVARIANT_BROKEN` and refuses to trade against a broken model. See `tests/kalshi_mlb_rfq/test_ev_calc.py::test_math_invariant_both_sides_cannot_be_positive` for the parametric proof across realistic LP spreads.

**Cross-LP arb falls out for free.** Different LPs are independent — nothing forces `yes_ask_A + no_ask_B > 1`. When LP_A's YES and LP_B's NO both clear our gate against the same fair, the math says `yes_ask_A + no_ask_B + fees < 1` — combined cost < $1 guaranteed payout = profitable hedge. Since `_evaluate_quote` runs per-quote, the bot picks both up naturally without any special-case code.

### Hedge-formation diagnostic (non-blocking)

When an accept creates a position on the side opposite a row we already hold on the same combo, the `quote_log` row is tagged:

| Column | Meaning |
|---|---|
| `hedge_added` | TRUE when a same-combo opposite-side position exists at accept time |
| `hedge_original_side` | `"yes"` or `"no"` — which side we already held |
| `hedge_original_price` | The held row's `weighted_price` |
| `hedge_new_price` | The ask we just paid on the new fill |
| `hedge_current_fair` | Fair value at time of new fill (useful for fair-drift analysis) |
| `hedge_projected_net` | `1 − held_price − new_price` (before fees), the net P&L on the matched contracts at settlement |

The diagnostic does **not** block the accept. The forward-looking math says each individual +EV accept is correct even when it creates a hedge — the sunk YES is sunk, and adding NO when NO is +EV against current fair improves expected P&L from here. The risk we measure (not prevent) is the case where fair has drifted by less than `2 × fees` between the two fills, which makes the *cumulative* combined position a small loss. If after 1-2 weeks of data the hedge-pattern P&L is meaningfully negative, we'll revisit with a "block when projected_net < threshold" rule.

Useful query:
```sql
-- Aggregate hedge P&L pattern (need >= 1-2 weeks of data)
SELECT
  COUNT(*) AS hedge_count,
  AVG(hedge_projected_net) AS avg_projected_net,
  SUM(hedge_projected_net) AS total_projected_net
FROM quote_log
WHERE decision = 'accepted' AND hedge_added = TRUE;
```

### Per-side cooldown

`combo_cooldown` PK is `(leg_set_hash, side)` — after the bot accepts YES on combo X it cools down YES on that combo, but NO on the same combo remains eligible for its own RFQ. This is intentional: they're genuinely different positions in our book.

### Schema v2 (2026-05-22)

Three tables changed shape to support NO-side accepts:

| Table | Change |
|---|---|
| `positions` | PK: `(combo_market_ticker, side)`. Existing rows backfill to `side='yes'`. Opposite-side fills no longer corrupt the upsert path. |
| `combo_cooldown` | PK: `(leg_set_hash, side)`. Same migration. |
| `quote_log` | New columns: `chosen_side`, `ev_yes_pct`, `ev_no_pct`, `hedge_added`, `hedge_original_side`, `hedge_original_price`, `hedge_new_price`, `hedge_current_fair`, `hedge_projected_net`. |

The v2 migration in `db._migrate_v2_side_columns` is idempotent, transactional, and backs up the DB to `<DB_PATH>.bak.v2_side_columns.<timestamp>` before the destructive PK swap. Runs automatically on `db.init_database()` at bot startup.

Partial accepts exist at the protocol level via the FIX interface (`OrderQty` on `35=UA`), but the REST endpoint doesn't expose it. If partial accepts become necessary, a FIX-protocol migration would be required (separate plan).

## Walk diagnostics

Every row in `quote_log` carries walk-diagnostic context so we can tell apart **"we were too slow"** vs **"we were too cheap"** when a quote we tried to accept walked.

**Columns** (in addition to the core decision fields):

| Column | Meaning |
|---|---|
| `competitor_count` | Number of OTHER open quotes on the same RFQ at poll time |
| `best_competitor_no_bid_dollars` | Max `no_bid_dollars` among competitors — the better quote we could have tried instead |
| `accept_response_body` | Kalshi's error body on walk (e.g. `quote_expired`, `rfq_closed`) — null on success |
| `rfq_terminal_status` | Status of the RFQ after walk, from a follow-up `get_rfq` call |
| `quote_first_seen_at` | When our `poll_quotes` returned this quote |
| `accept_attempted_at` | Right before we sent the PUT /accept |
| `accept_response_at` | Right after Kalshi responded to our /accept |

**Latency derivation:** `accept_response_at - quote_first_seen_at` = end-to-end round-trip. `accept_response_at - accept_attempted_at` = pure Kalshi PUT latency. The gap between them is internal bot processing.

**Useful queries:**

```sql
-- How many walks were races vs price issues, last 24h?
SELECT
  CASE WHEN best_competitor_no_bid_dollars > no_bid_dollars THEN 'better_quote_existed'
       WHEN competitor_count = 0 THEN 'sole_quote_walked'
       ELSE 'we_were_best' END AS bucket,
  COUNT(*) AS n,
  AVG((accept_response_at - quote_first_seen_at) * 1000) AS avg_ms
FROM quote_log
WHERE decision = 'failed_quote_walked' AND observed_at > now() - INTERVAL 24 HOUR
GROUP BY bucket;

-- Latency distribution on walks vs accepts
SELECT decision,
       percentile_cont(0.5) WITHIN GROUP (ORDER BY epoch(accept_response_at - quote_first_seen_at) * 1000) AS p50_ms,
       percentile_cont(0.95) WITHIN GROUP (ORDER BY epoch(accept_response_at - quote_first_seen_at) * 1000) AS p95_ms
FROM quote_log
WHERE accept_response_at IS NOT NULL
GROUP BY decision;

-- Walks by Kalshi error reason
SELECT accept_response_body, COUNT(*) FROM quote_log
WHERE decision = 'failed_quote_walked' GROUP BY 1 ORDER BY 2 DESC;
```

## SGP cadence loop (line-source pivot, 2026-05-13)

The bot drives its own SGP scrape cadence independent of the MLB dashboard.

**Data flow on each SGP tick (every `SGP_REFRESH_SEC`, default 60s):**

1. Bot enumerates open Kalshi MVE markets per MLB game (every `(spread, total)` tuple Kalshi lists).
2. Bot rewrites `mlb_target_lines` in `kalshi_mlb_rfq_market.duckdb` (sibling to state DB).
3. Bot spawns the four scrapers (`mlb_sgp/scraper_{draftkings,fanduel,prophetx,novig}_sgp.py`) with `MLB_SGP_DB_PATH=<market DB>` and `MLB_SGP_PERIODS=FG` env overrides.
4. Scrapers read `mlb_target_lines`, price every tuple at their respective book, write back to `mlb_sgp_odds` in the bot's market DB with new `spread_line`/`total_line` columns.
5. Bot reloads `_SGP_ODDS_CACHE` from the bot market DB.

**Edge surface:** any Kalshi MVE combo with ≥2 books priced at the matching (spread, total). Off-line combos (only 1 book) are dropped — the bot does not bet model-only or single-book candidates.

**Schedule source:** game IDs and team metadata come from `Answer Keys/mlb.duckdb::mlb_odds_temp` (read-only). The bot has no dependency on Wagerzon-derived `mlb_parlay_lines` anymore.

**Cold start:** the first SGP cycle runs synchronously before `main_loop` enters its tick. Bot blocks ~60-90s on startup.

**Config:**
- `SGP_REFRESH_SEC` (default 60) — SGP cadence interval
- `SGP_SCRAPER_TIMEOUT_SEC` (default 90) — per-scraper kill deadline
- `BOT_MARKET_DB` (default `kalshi_mlb_rfq_market.duckdb` in this package) — sibling market DB
- `MIN_BOOK_COUNT_FOR_BLEND` (default 2) — drop-candidate threshold
