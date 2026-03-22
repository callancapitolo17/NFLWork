# Kalshi Market Maker — AI Context

## What This Is
Automated market maker for CBB 1st-half markets on Kalshi. Posts bid/ask quotes, monitors fills, manages risk.

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
Fair value comes from offshore scraper databases (Bookmaker, Bet105). If those scrapers haven't run recently, the oracle is stale and quotes will be wrong.

### DuckDB File
`kalshi_mm.duckdb` — NEVER symlink. WAL files must be co-located with the database.

### Tipoff Safety
Bot cancels all resting orders and stops quoting/taking 1 minute before game tipoff (`TIPOFF_PULLBACK_MIN`). `sweep_tipoff_cancel()` runs every quote cycle (10s), independent of market matching — it sweeps ALL entries in `resting_by_ticker` using stored `_commence_time`. Cancels use `batch_cancel()` and only remove local tracking after confirmed API success. If cancel fails, the sweep retries next cycle. If `commence_time` is missing or unparseable, the bot refuses to quote (fail-safe).

### Line Move Cancel
When the line monitor detects a sharp move, all resting orders for affected games are immediately batch-cancelled using stored `_home_team`/`_away_team` metadata in `resting_by_ticker`. This avoids stale quotes while the prediction pipeline refreshes.

### Batch Order Placement
New orders are batched via `POST /portfolio/orders/batched` (max 20 per call) instead of placed individually. This reduces first-cycle startup from ~30 min to ~2 min for ~710 orders. Amends remain individual API calls. Cancels within the quote cycle are also batched via `batch_cancel()`. The `batch_place()` function in `orders.py` handles chunking and per-order failure logging.

### Batch Kelly Sizing
Kelly sizes are pre-computed for all markets grouped by game before the main quote loop. Instead of calling `conditional_kelly_sizes()` twice per market (~610 calls), `batch_kelly_sizes_for_game()` computes all bid+ask sizes for a game in one matrix solve (~57 calls). This reduces the Kelly loop from ~15 min to ~1.5 min.

### Taker Cooldowns
The taker uses a simple 10s per-ticker tactical cooldown to prevent hammering the same contract. Cross-market correlation (same game, same event) is handled entirely by Kelly conditional sizing — `kelly.clear_positions_cache()` is called after each fill so subsequent takes see updated positions. No event-level or game-level cooldowns.

### Sample Cache Pre-warm
At startup and after every prediction refresh, `kelly.prewarm_sample_cache()` bulk-loads all game simulation samples in one DuckDB query (~15K rows, <1MB). This eliminates the ~15 min cold-cache Kelly cycle on first boot. `clear_sample_cache()` is now an alias that clears and reloads.

### Prediction Auto-Reload
The main loop checks `cbb_prediction_meta.updated_at` every 30s. If predictions were updated externally (manual pipeline run, cron), the bot reloads predictions, re-matches markets, and rewarms the sample cache — no restart needed.

### Shared DB Write Connection
`db.py` uses a single shared write connection (`get_write_conn()`) opened lazily on first use and closed on shutdown (`close_write_conn()`). All write functions use it instead of opening/closing their own connections. This eliminates ~1,200 connection open/close cycles per quote cycle. Init functions (`init_database`, `init_taker_tables`) keep their own connections since they run before the shared connection is available. Batch versions exist for hot-path writes: `batch_save_resting_orders()`, `batch_remove_resting_orders()`, `flush_quote_log()`.

### Positions: Kalshi API is Source of Truth
Kelly sizing and the exposure cap read positions from the Kalshi API (`GET /portfolio/positions`), NOT the local DuckDB. This prevents position drift when fills are missed by `poll_for_fills` (e.g., fully-filled orders that disappear between poll cycles). The local DB is still updated on fills for audit logging, but is never used for sizing decisions. Positions are cached per cycle with a 5s TTL to avoid hammering the API from the taker's 1s poll loop.

### Kalshi API Auth
Uses `KALSHI_API_KEY` and `KALSHI_PRIVATE_KEY` from `.env`. Keys are RSA — the private key file path goes in `.env`.

## When Making Changes
- Test with minimum position sizes first
- The bot runs continuously — changes require restart
- Check risk.py limits before increasing position sizes
