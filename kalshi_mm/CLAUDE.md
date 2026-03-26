# Kalshi Market Maker — AI Context

## What This Is
Automated market maker for CBB 1st-half markets on Kalshi. Posts bid/ask quotes, monitors fills, manages risk. Supports spreads, totals, moneylines, and race-to-10 props.

## Architecture
- `main.py` — Bot loop: fetch prices → compute quotes → place orders → monitor
- `quoter.py` — Pricing engine: uses oracle (Bookmaker/Bet105 scrapers) to set fair value, then applies spread
- `orders.py` — Order management: cancel stale, place new, track fills
- `risk.py` — Position limits, exposure tracking, P&L computation
- `taker.py` — Taker mode: crosses the spread when oracle shows strong edge
- `scan_book.py` — Scans all available markets for opportunities
- `config.py` — Configuration (spread width, position limits, update interval)
- `db.py` — DuckDB persistence (orders, fills, positions)

## Critical Details

### Oracle Dependency
Fair value comes from offshore scraper databases (Bookmaker, Bet105). If those scrapers haven't run recently, the oracle is stale and quotes will be wrong. Predictions and game samples are read from `cbb_mm.duckdb` (separate from `cbb.duckdb` to avoid lock contention during pipeline runs).

### DuckDB File
`kalshi_mm.duckdb` — NEVER symlink. WAL files must be co-located with the database.

### Tipoff Safety
Bot cancels all resting orders and stops quoting/taking 1 minute before game tipoff (`TIPOFF_PULLBACK_MIN`). `sweep_tipoff_cancel()` runs every quote cycle (10s), independent of market matching — it sweeps ALL entries in `resting_by_ticker` using stored `_commence_time`. Cancels use `batch_cancel()` and only remove local tracking after confirmed API success. If cancel fails, the sweep retries next cycle. If `commence_time` is missing or unparseable, the bot refuses to quote (fail-safe).

### Line Move Cancel
When the line monitor detects a sharp move, all resting orders for affected games are immediately batch-cancelled using stored `_home_team`/`_away_team` metadata in `resting_by_ticker`. This avoids stale quotes while the prediction pipeline refreshes.

### Batch Order Placement
New orders are batched via `POST /portfolio/orders/batched` (max 20 per call) instead of placed individually. This reduces first-cycle startup from ~30 min to ~2 min for ~710 orders. Amends remain individual API calls. Cancels within the quote cycle are also batched via `batch_cancel()`. The `batch_place()` function in `orders.py` handles chunking and per-order failure logging.

### Budget-Based Maker Sizing
Maker sizes use a budget approach instead of Kelly matrix optimization. Markets are grouped by game + type (spreads, totals, moneyline). Each group gets a budget = best standalone Kelly fraction × kelly_mult × bankroll, capped at MAX_GAME_TYPE_EXPOSURE_PCT. The budget is distributed within the group by Kelly weight (proportional to each ticker's standalone Kelly fraction). Bids and asks share one budget per group to prevent cross-side directional doubling. Existing Kalshi positions are subtracted from the budget before distributing. When over budget, the bot posts orders to unwind. Kelly matrix sizing (`batch_kelly_sizes_for_game`) is retained for the taker only.

**Quoter-price sizing:** The quote cycle runs the quoter first (prices only) on all markets, then passes those prices to `compute_maker_sizes(quoter_prices=...)`. Dollar-to-contract conversion uses the quoter's actual posting price — not the Kalshi book price — so that `contracts × cost_per_contract` never exceeds the dollar budget. The quoter determines the price, Kelly determines the size at that price.

### Taker Cooldowns
The taker uses a simple 10s per-ticker tactical cooldown to prevent hammering the same contract. Cross-market correlation (same game, same event) is handled entirely by Kelly conditional sizing — `kelly.clear_positions_cache()` is called after each fill so subsequent takes see updated positions. No event-level or game-level cooldowns.

### Sample Cache Pre-warm
At startup and after every prediction refresh, `kelly.prewarm_sample_cache()` bulk-loads all game simulation samples in one DuckDB query (~15K rows, <1MB). This eliminates the ~15 min cold-cache Kelly cycle on first boot. `clear_sample_cache()` is now an alias that clears and reloads.

### Sample Column Layout
`kelly.py` defines `SAMPLE_COLUMNS = ["home_margin_h1", "total_h1", "first_to_10_h1"]` with named column indices. The cache dynamically discovers which columns exist in `cbb_game_samples` and NaN-fills missing ones. This allows the R pipeline to add new columns without breaking the running bot.

### Race-to-10 Markets
- Series: `KXNCAAMBFIRST10` — 2 contracts per event (one per team)
- Title format: "Will [Team] be the first to reach 10 points?"
- Fair value: `mean(first_to_10_h1)` from resampled historical games in `cbb_game_samples`
- `first_to_10_h1` = 1 (home first), 0 (away first), NA (neither/push)
- Prediction rows: `race_to_10_h1` in `cbb_raw_predictions` table
- Backfill: `Rscript "Acquire CBB Data.R" --backfill-first10` to populate existing `cbb_pbp_v2` rows

### Prediction Auto-Reload
The main loop checks `cbb_prediction_meta.updated_at` every 30s. If predictions were updated externally (manual pipeline run, cron), the bot reloads predictions, re-matches markets, and rewarms the sample cache — no restart needed.

### Shared DB Write Connection
`db.py` uses a single shared write connection (`get_write_conn()`) opened lazily on first use and closed on shutdown (`close_write_conn()`). All write functions use it instead of opening/closing their own connections. This eliminates ~1,200 connection open/close cycles per quote cycle. Init functions (`init_database`, `init_taker_tables`) keep their own connections since they run before the shared connection is available. Batch versions exist for hot-path writes: `batch_save_resting_orders()`, `batch_remove_resting_orders()`, `flush_quote_log()`.

### Positions: Kalshi API is Source of Truth
Kelly sizing and the exposure cap read positions from the Kalshi API (`GET /portfolio/positions`), NOT the local DuckDB. This prevents position drift when fills are missed by `poll_for_fills` (e.g., fully-filled orders that disappear between poll cycles). The local DB is still updated on fills for audit logging, but is never used for sizing decisions. Positions are cached per cycle with a 5s TTL to avoid hammering the API from the taker's 1s poll loop.

### Manual Order Preservation
Bot-placed orders are tagged with `client_order_id` prefixed `mm_` (configured via `config.BOT_ORDER_PREFIX`). All cancel logic — phantom detection, `_nuke_and_adopt()`, shutdown kill switch — checks this prefix and skips orders without it. This lets users place manual sell orders on Kalshi without the bot cancelling them. The `is_bot_order(order)` helper in `orders.py` is the single source of truth for this check, using `(order.get("client_order_id") or "")` to handle null values from the API.

### Kalshi API Auth
Uses `KALSHI_API_KEY` and `KALSHI_PRIVATE_KEY` from `.env`. Keys are RSA — the private key file path goes in `.env`.

## When Making Changes
- Test with minimum position sizes first
- The bot runs continuously — changes require restart
- Check risk.py limits before increasing position sizes
