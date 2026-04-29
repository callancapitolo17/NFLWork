# NFLWork

Quantitative sports betting edge detection and execution platform.

## What This Does

Scrapes live odds from 6+ sportsbooks, generates statistically-backed fair prices using historical game matching, identifies +EV opportunities, and tracks closing line value (CLV) to validate edge.

## Directory Structure

```
NFLWork/
├── Answer Keys/              # Core pricing models (CBB, NFL, MLB, CFB, EPL, March Madness)
│   ├── CBB Answer Key/       #   College basketball model + backtests
│   ├── CBB Dashboard/        #   Interactive dashboard (Flask + R-generated HTML)
│   ├── NFL Answer Key/       #   NFL model + backtests
│   ├── NFL Dashboard/        #   NFL interactive dashboard
│   ├── MLB Answer Key/       #   MLB consensus model
│   ├── College Football Answer Key/
│   ├── EPL Answer Key/       #   Soccer (English Premier League)
│   ├── Tools.R               #   Shared R utilities (devig, odds conversion, sampling)
│   ├── canonical_match.py    #   Team name resolution (offshore → canonical)
│   ├── run.py                #   Pipeline orchestrator (parallel scrapers + R)
│   ├── parlay.R              #   Fair parlay pricing from samples
│   └── props.R               #   Fair prop pricing from play-by-play
│
├── bet_logger/               # Bet history scrapers → Google Sheets
├── bet_placer/               # Automated bet placement via Playwright
│
├── bet105_odds/              # Bet105 odds scraper (WebSocket)
├── bfa_odds/                 # BFA Gaming odds scraper (REST, no auth)
├── bookmaker_odds/           # Bookmaker.eu odds scraper (curl_cffi)
├── hoop88_odds/              # Hoop88 odds scraper (REST + JWT)
├── kalshi_odds/              # Kalshi prediction market scraper (REST)
├── wagerzon_odds/            # Wagerzon odds scraper (ASP.NET)
│
├── hoop88_correlation/       # Correlated parlay edge finder
├── sharp_analyst/            # AI-powered sharp line movement detector
│
├── kalshi_mm/                # Kalshi market maker bot (CBB 1H)
├── kalshi_mlb_rfq/           # Kalshi MLB SGP RFQ taker bot (cross-category MVE)
├── kalshi_draft/             # Kalshi NFL draft prediction markets
├── nfl_draft/                # Cross-venue NFL Draft EV portal (Kalshi + 5 books)
├── kalshi coaching/          # Kalshi coaching hire markets
│
├── March Madness/            # Tournament simulator (Monte Carlo)
└── docs/                     # Agent & plugin specifications
```

## Architecture

### Data Flow

```
Odds API + 6 Offshore Scrapers  →  DuckDB (per-book)
        ↓
   run.py (orchestrator)        →  Parallel scraper execution
        ↓
   Answer Key (R)               →  Historical game matching → fair prices
        ↓
   Dashboard (Flask + HTML)     →  +EV bets, Kelly sizing, CLV tracking
        ↓
   bet_placer (Playwright)      →  Automated bet placement
        ↓
   bet_logger                   →  Google Sheets P&L tracking
```

### Key Concepts

- **Answer Key Algorithm**: Matches current game lines to historically similar games, builds balanced samples, prices arbitrary markets from actual outcomes
- **Devigging**: Removes bookmaker vig to get true implied probabilities before comparing across books
- **CLV (Closing Line Value)**: Tracks whether bets beat the closing line — the ultimate edge metric
- **Kelly Sizing**: Fractional Kelly criterion (25-50%) for position sizing based on edge magnitude

## Tech Stack

| Layer | Technology |
|-------|-----------|
| Odds scraping | Python (Playwright, requests, curl_cffi, websockets) |
| Statistical models | R (cbbdata, hoopR, Rcpp for performance) |
| Data storage | DuckDB (one database per scraper/tool) |
| Dashboards | R (reactable → static HTML) + Flask (API + serving) |
| Bet placement | Python (Playwright browser automation) |
| Bet tracking | Python → Google Sheets API |

## Shared Conventions

- **Team names**: All scrapers resolve to canonical names via `canonical_match.py` (Python) or `resolve_offshore_teams()` (R in Tools.R)
- **DuckDB schema**: 18-column standard across all odds scrapers (fetch_time, sport_key, game_id, game_date, game_time, away_team, home_team, market, period, away_spread, away_spread_price, home_spread, home_spread_price, total, over_price, under_price, away_ml, home_ml)
- **Environment**: All credentials in `bet_logger/.env` (shared across scrapers)
- **No temp files**: DuckDB tables for all shared state; no CSV/RDS intermediates

## Quick Start

```bash
# Run the CBB pipeline (scrapers + model + dashboard)
cd "Answer Keys/CBB Dashboard"
bash run.sh

# Run a single odds scraper
cd bet105_odds && python scraper.py cbb

# Log bet history to Google Sheets
cd bet_logger && bash run_all_scrapers.sh
```

## NFL Draft Portal

See [`nfl_draft/README.md`](nfl_draft/README.md) for the cross-venue EV portal
(Kalshi + DK/FD/Bookmaker/Wagerzon/Hoop88). Dashboard runs at http://127.0.0.1:8090/.

