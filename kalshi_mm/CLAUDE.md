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
Bot cancels all resting orders and stops quoting/taking 1 minute before game tipoff (`TIPOFF_PULLBACK_MIN`). `sweep_tipoff_cancel()` runs every quote cycle (10s), independent of market matching — it sweeps ALL entries in `resting_by_ticker` using stored `_commence_time`. If `commence_time` is missing or unparseable, the bot refuses to quote (fail-safe).

### Kalshi API Auth
Uses `KALSHI_API_KEY` and `KALSHI_PRIVATE_KEY` from `.env`. Keys are RSA — the private key file path goes in `.env`.

## When Making Changes
- Test with minimum position sizes first
- The bot runs continuously — changes require restart
- Check risk.py limits before increasing position sizes
