# NFL Draft Prediction Market Dashboard

Interactive Dash dashboard for tracking NFL Draft prediction markets on Kalshi. Auto-discovers all draft series, tracks price history, detects cross-market edges, and compares to mock draft consensus.

> **Status (April 2026):** The dashboard has been extended with the NFL Draft EV
> portal (see [`../nfl_draft/README.md`](../nfl_draft/README.md)). New tabs —
> Cross-Book Grid, +EV Candidates, Trade Tape, Bet Log — render data for Kalshi
> plus five sportsbooks (DK, FD, Bookmaker, Wagerzon, Hoop88). The database pointer
> has been repointed at `nfl_draft/nfl_draft.duckdb`; the legacy
> `kalshi_draft.duckdb` has been migrated in and retired. Dashboard port moved
> from 8083 to 8090 (override via `NFL_DRAFT_DASHBOARD_PORT`).

## Features

- **Auto-Discovery**: Finds all NFL Draft series on Kalshi dynamically (no hardcoded tickers)
- **Market Overview**: Heatmap + sortable/filterable table of all draft markets
- **Price History**: Time-series charts showing how odds move over time
- **Edge Detection**: Sum-to-1 violations, cross-market inconsistencies, spread analysis
- **Consensus Comparison**: Tankathon big board vs Kalshi implied probabilities
- **Portfolio Tracking**: Positions, resting orders, P&L (requires API credentials)

## Setup

```bash
cd kalshi_draft

# Create virtual environment
python3 -m venv venv
./venv/bin/pip install -r requirements.txt

# (Optional) Set up Kalshi API credentials for portfolio tracking
cp .env.example .env
# Edit .env with your credentials
```

## Usage

```bash
# Full pipeline: fetch data + start dashboard
./run.sh

# Fetch data only (no server)
./run.sh --fetch-only

# Run individual components
./venv/bin/python fetcher.py        # Fetch odds from Kalshi
./venv/bin/python edge_detector.py  # Run edge detection
./venv/bin/python consensus.py      # Scrape mock draft consensus
./venv/bin/python app.py            # Start dashboard server
```

Dashboard runs at `http://127.0.0.1:8090` (override via `NFL_DRAFT_DASHBOARD_PORT`).

## Architecture

```
fetcher.py ──→ DuckDB ──→ app.py (Dash)
                 ↑
edge_detector.py
consensus.py
```

- **fetcher.py**: Queries Kalshi `/series` and `/markets` endpoints, stores snapshots
- **edge_detector.py**: Analyzes latest odds for mispricing
- **consensus.py**: Scrapes Tankathon big board, fuzzy-matches to Kalshi players
- **db.py**: DuckDB schema + query helpers
- **auth.py**: RSA-PSS signing for authenticated Kalshi endpoints
- **app.py**: Dash application with 5 tabs

## Kalshi Market Series

Auto-discovered series (as of Feb 2026):

| Series | Description | Markets |
|--------|-------------|---------|
| KXNFLDRAFT1 | #1 overall pick (by player) | ~48 |
| KXNFLDRAFT1ST | Which team picks #1 | ~32 |
| KXNFLDRAFTPICK | Draft pick position | ~152 |
| KXNFLDRAFTTOP | Top N draft range | ~114 |

## Edge Types Detected

1. **Sum Violations**: Markets within a mutually-exclusive event should sum to ~100%
2. **Cross-Market**: Player probabilities should be consistent across series (e.g., #1 pick prob <= top 5 prob)
3. **Spread Opportunities**: Wide bid/ask spreads with stale last-trade prices

## Data Storage

All data stored in `nfl_draft/nfl_draft.duckdb` (migrated from `kalshi_draft.duckdb`):
- `kalshi_odds`: Historical odds snapshots (append-only; formerly `draft_odds`)
- `draft_series`: Discovered series metadata
- `market_info`: Cached market details
- `positions` / `resting_orders`: Portfolio snapshots
- `consensus_board`: Mock draft rankings
- `detected_edges`: Computed edges
