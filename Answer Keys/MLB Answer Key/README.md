# MLB Answer Key

Consensus betting odds model for MLB. Aggregates moneyline and totals from multiple sportsbooks, devigs, and weights by historical accuracy.

## Pipeline

1. **Data Acquisition** (`Acquire New MLB Data.R`): Fetches live odds via Odds API, loads historical outcomes from DuckDB
2. **Consensus Building** (`MLB Answer Key 2.0.R`): Devigs each book's odds, computes weighted average (sharper books weighted higher)
3. **Output**: `mlb_odds` table with consensus devigged home/away odds and totals

## Key Files

| File | Purpose |
|------|---------|
| `MLB Answer Key 2.0.R` | Main consensus model (current) |
| `MLB Answer Key.R` | Legacy version with mean-matching algorithm |
| `Acquire New MLB Data.R` | Data fetching and preprocessing |
| `MLB Odds API Scraping.R` | OddsAPI scraper |
| `MLB Updated API Pitch Data.R` | Pitch data from API |
| `Back Testing Odds.R` | Backtesting framework |
| `Consensus Betting History.R` | Historical consensus tracking |
| `run_fetch.sh` | Friday-only cron scheduler |

## Scheduling

`run_fetch.sh` runs on Fridays only (date-gated). Logs to `~/Library/Logs/fetch-mlb-odds/`.

## Data Storage

Historical odds and outcomes stored in `pbp.duckdb` (table: `mlb_betting_pbp`).
