# MLB Answer Key

Pricing engine for MLB First-5-Innings (F5) markets. Generates fair odds by matching current game lines to historically similar games, then compares predictions to offshore sportsbook odds to identify +EV opportunities.

## Architecture

```
run.py mlb (orchestrator)
  ├── [first]    Sharp scrapers (bookmaker, bet105) → their DuckDB files
  ├── [parallel] Other scrapers (wagerzon, hoop88, bfa) + MLB.R
  ├── Sentinel: .scrapers_done_mlb signals scrapers complete
  └── MLB.R reads scraper DBs, generates fair prices, writes pipeline output
```

## Markets Priced

| Market | Period | Odds API Key |
|--------|--------|-------------|
| Moneyline | F5 | `h2h_1st_5_innings` |
| Totals | F5 | `totals_1st_5_innings` |
| Spreads | F5 | `spreads_1st_5_innings` |

## Key Files

| File | Purpose |
|------|---------|
| `MLB.R` | Main merged pipeline (8 phases) |
| `Acquire New MLB Data.R` | Historical odds + PBP fetching |
| `clv_compute.py` | Post-game CLV computation |
| `run_mlb_daily.sh` | Daily acquisition scheduler |

## Setup

### Prerequisites
- R with packages: `data.table`, `oddsapiR`, `duckdb`, `dplyr`, `tidyr`, `lubridate`, `httr`, `jsonlite`
- Python 3.10+ with `flask`, `duckdb`, `requests`, `scipy`
- `ODDS_API_KEY` in `~/.Renviron`
- Historical data in `Answer Keys/pbp.duckdb` (table: `mlb_betting_pbp`)

### Running

```bash
# Full pipeline (scrapers + predictions + dashboard)
cd "Answer Keys/MLB Dashboard"
bash run.sh

# Pipeline only (no dashboard)
cd "Answer Keys"
python3 run.py mlb

# R script standalone (no scrapers)
cd "Answer Keys"
Rscript "MLB Answer Key/MLB.R"
```

### Daily Data Acquisition

```bash
# Manual run
cd "Answer Keys/MLB Answer Key"
bash run_mlb_daily.sh

# Logs to ~/Library/Logs/mlb-daily-acquire/ with 30-day rotation
```

### CLV Computation

```bash
cd "Answer Keys/MLB Answer Key"
python3 clv_compute.py
```

## DuckDB Tables

| Database | Table | Purpose |
|----------|-------|---------|
| `pbp.duckdb` | `mlb_betting_pbp` | Historical games with inning-by-inning outcomes (12,719 games) |
| `mlb.duckdb` | `mlb_bets_combined` | Pipeline output (daily +EV bets) |
| `mlb.duckdb` | `mlb_team_dict` | Team name dictionary |
| `mlb_dashboard.duckdb` | `placed_bets` | Bet placement tracking |
| `mlb_dashboard.duckdb` | `bet_clv` | Post-game CLV results |

## Key Design Decisions

1. **Moneyline-based matching** (`use_spread_line = FALSE`): MLB uses ML probability (not spread) as the primary consensus metric for sample generation.
2. **Sharp-only consensus**: Uses `SHARP_BOOKS` from Tools.R (Pinnacle 1.1, Bookmaker 1.1, LowVig/Circa/Bet105 1.0). Rec books excluded.
3. **F5 only (initial)**: Full-game markets deferred as follow-up work.
4. **Inning columns are cumulative**: `game_home_margin_inning_inning_5` = home margin through first 5 innings.
