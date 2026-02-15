# NFL +EV Betting Dashboard

Interactive dashboard for viewing and managing positive expected value bets from the NFL Answer Key pipeline.

## Features

- **+EV Bet Display**: Shows all bets sorted by expected value
- **Multi-Book Support**: Displays bets from API books (DraftKings, FanDuel, etc.) and offshore scrapers (Wagerzon, Hoop88, BFA)
- **Bet Tracking**: Mark bets as placed, persisted in local database
- **Correlation Warnings**: Alerts when a bet correlates with one you've already placed
  - Red = high correlation (90%+, same market group)
  - Yellow = medium correlation (60-80%, cross-group like ML ↔ spreads)
  - Hover to see the specific correlated bet
- **Filtering**: Filter by game, market type, book, or search
- **One-Click Refresh**: Runs full pipeline (scrapers + R predictions)

## Quick Start

```bash
# First time: Generate predictions
cd ~/NFLWork
python "Answer Keys/run.py" nfl

# Start the dashboard
cd "Answer Keys/NFL Dashboard"
./run.sh
```

Dashboard opens at http://127.0.0.1:8081

**Refresh button** runs the full pipeline automatically:
1. Scrapers (wagerzon, hoop88, bfa) run in parallel
2. R sample generation runs in parallel
3. NFLCombine.R merges everything and finds +EV bets
4. Dashboard reloads with fresh data

## File Structure

```
NFL Dashboard/
├── nfl_dashboard.R           # Generates HTML from duckdb data
├── nfl_dashboard_server.py   # Flask server (serves HTML + API)
├── nfl_dashboard.duckdb      # Placed bets tracking
├── market_relationships.json # Correlation rules
├── report.html               # Generated dashboard
├── run.sh                    # Startup script
├── lib/                      # Static JS/CSS assets
└── README.md
```

## How It Works

1. **Data Source**: Reads `nfl_bets_combined` table from `pbp.duckdb` (saved by NFLCombine.R)
2. **HTML Generation**: R script creates interactive reactable with all bets
3. **Flask Server**: Serves HTML and handles bet placement API
4. **Refresh**: Calls `python run.py nfl` then regenerates HTML

## API Endpoints

| Endpoint | Method | Description |
|----------|--------|-------------|
| `/` | GET | Serve the dashboard |
| `/refresh` | POST | Run full pipeline and regenerate |
| `/api/place-bet` | POST | Mark a bet as placed |
| `/api/remove-bet` | POST | Remove from placed status |
| `/api/placed-bets` | GET | Get all placed bets |
| `/api/exposure` | GET | Get exposure summary by game |

## Correlation Groups

The dashboard warns about correlated bets to prevent over-exposure:

| Group | Markets |
|-------|---------|
| h2h | ML Q1-Q4, ML H1-H2 |
| h2h_3way | ML-3 Q1-Q4, ML-3 H1-H2 |
| spreads | Spread Q1-Q4, Spread H1-H2, alternates |
| totals | Total Q1-Q4, Total H1-H2, alternates |
| team_totals | Team totals by period |

**Cross-group correlations:**
- ML ↔ Spreads: 75%
- ML ↔ ML-3: 90%
- Totals ↔ Team Totals: 60%
- Spreads ↔ Team Totals: 50%

## Manual Setup

### Requirements
- Python 3.x with flask, duckdb
- R with: tidyverse, duckdb, reactable, htmltools, htmlwidgets, jsonlite, digest

### Install
```bash
cd "Answer Keys/NFL Dashboard"
python3 -m venv venv
source venv/bin/activate
pip install flask duckdb pandas numpy
```

### Run
```bash
source venv/bin/activate
python3 nfl_dashboard_server.py
```
