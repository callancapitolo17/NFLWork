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
- `MIN_EV_PCT`, `MAX_QUOTE_DEVIATION` — accept thresholds
- `MAX_GAME_EXPOSURE_PCT`, `DAILY_EXPOSURE_CAP_USD` — exposure caps
- `MAX_PREDICTION_STALENESS_SEC` — accept-gate staleness threshold (10 min default)
- `MIN_FILL_RATIO`, `FILL_RATIO_WINDOW` — adverse-selection halt

## Safety scaffold (per-accept gates)

Every quote acceptance must pass ALL of:
- Min EV after fee · Fair-value bounds · Sanity bound (quote vs fair) · Prediction staleness
- Tipoff window (5 min before first pitch) · Line-move check · Per-game exposure cap (% of bankroll)
- Daily exposure cap · Kill-switch off · Inverse-combo not held · 2-source gate (model + ≥1 book)
- Per-combo cooldown (30s post-fill) · Positions API health · Fill-ratio halt (rolling 50 attempts)

Acceptance is serialized via `ACCEPT_LOCK` so concurrent quotes Kelly-size against fresh DB state.

## Troubleshooting

- **No combos surfacing:** check that `mlb_game_samples` and `mlb_sgp_odds` are populated. Run the MLB pipeline first.
- **All quotes declining `declined_stale_predictions`:** rerun the pipeline; samples are over 10 minutes old.
- **Bot halted on `fill_ratio_collapse`:** investigate — makers are walking on accepts at a rate that suggests adverse selection. Check `quote_log` for the maker `creator_id`s causing it.
- **`mint_combo_ticker` failing with 400:** the MVE collection ticker may have changed or one of the leg market_tickers doesn't exist. Re-run `mlb_sgp/recon_kalshi_mlb_rfq.py` for a fresh probe.
