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

## Knobs

See `.env.example`. Most relevant:
- `BANKROLL`, `KELLY_FRACTION` — sizing
- `MIN_EV_PCT` — accept threshold
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

## Accept semantics (corrected 2026-05-21)

Kalshi's REST `PUT /communications/quotes/{id}/accept` is **all-or-nothing** — body is `{"accepted_side": "yes"|"no"}` and Kalshi returns 204 on success. There is no per-accept contract count parameter.

Sizing happens **upstream at RFQ creation** via `target_cost_dollars`. With `rest_remainder: False`, the LP's quote is fill-or-kill at the requested size — accepting fills exactly what was requested.

Historical bug: from first commit (2026-05-02) through 2026-05-21, this code sent body `{"contracts": N}`, which Kalshi rejected as `invalid_parameters` on every accept. The bug was masked because the failure looked identical to "lost the race." The walk-diagnostic columns added 2026-05-20 finally exposed it — every walked row showed `accept_response_body = "invalid_parameters"` AND `rfq_terminal_status = "open"` (no competitor had filled the RFQ; we just couldn't accept it).

### Side semantics — `accepted_side` is INVERTED from intuition

The `accepted_side` field names the side of the LP's two-sided quote we're accepting — and accepting their bid on a side means *they buy that side from us*, leaving us LONG the OPPOSITE side. Verified empirically 2026-05-21 from the first real fill:

| Field value | What Kalshi does | What we end up holding |
|---|---|---|
| `accepted_side="yes"` | LP buys YES from us at `yes_bid` | **LONG NO** at `1 − yes_bid` (≈ `no_ask`) |
| `accepted_side="no"` | LP buys NO from us at `no_bid` | **LONG YES** at `1 − no_bid` (≈ `yes_ask`) |

The bot's EV gate is `ev_calc.post_fee_ev_buy_yes(fair, no_bid)` — it's evaluating whether to BUY YES at `1 − no_bid`. To actually open that position, the bot must send `accepted_side="no"`. The first ever fill landed at `no_price=$0.969` (1 NO contract) when the bot sent `"yes"`, costing a small −EV position before the bot was stopped. Single-contract scope thanks to the $1 RFQ default.

Partial accepts exist at the protocol level via the FIX interface (`OrderQty` on `35=UA`), but the REST endpoint doesn't expose it. If partial accepts become necessary, a FIX-protocol migration would be required.

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
