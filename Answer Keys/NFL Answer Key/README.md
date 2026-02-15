# NFL Answer Key

Generates derivative market predictions (quarters, halves) for NFL games by building consensus lines and sampling from historical game data. Compares predictions against multiple bookmakers (API + offshore scrapers).

## Quick Start

```bash
# Run the full pipeline (scrapers + predictions)
cd ~/NFLWork
python "Answer Keys/run.py" nfl

# Or use the dashboard
cd "Answer Keys/NFL Dashboard"
./run.sh
```

## Architecture

The pipeline runs in two phases:

### Phase 1: Parallel Execution
Runs simultaneously:
- **Scrapers** - wagerzon, hoop88, bfa (fetch offshore odds)
- **NFLPrepare.R** - Load data, build consensus, generate samples

### Phase 2: Combine
- **NFLCombine.R** - Load samples + all scraped odds, find +EV bets

```
run.py nfl
    │
    ├── [parallel] wagerzon scraper → wagerzon.duckdb
    ├── [parallel] hoop88 scraper   → hoop88.duckdb
    ├── [parallel] bfa scraper      → bfa.duckdb
    ├── [parallel] NFLPrepare.R     → pbp.duckdb (samples)
    │
    └── [sequential] NFLCombine.R
            │
            ├── Load samples from pbp.duckdb
            ├── Load odds from all scrapers
            ├── Generate predictions for all markets
            ├── Find +EV bets across all books
            └── Save to nfl_bets_combined table
```

## Files

| File | Purpose |
|------|---------|
| `NFLPrepare.R` | Phase 1: Load data, build consensus, generate samples, save to duckdb |
| `NFLCombine.R` | Phase 2: Load samples + scraped odds, generate predictions, find edge |
| `NFLAnswerKey2.0.R` | Legacy standalone script (deprecated) |

## Prerequisites

**Database:** `pbp.duckdb` in the Answer Keys directory containing:
- `nfl_betting_pbp` - Play-by-play betting data (2020+)
- `nfl_pre_20_betting_history` - Historical betting data (1999-2019)
- `nfl_weights` - Book weights for consensus building

**API:** The Odds API key set as environment variable `ODDS_API_KEY`

**Scrapers:** Python venvs set up in:
- `~/NFLWork/wagerzon_odds/`
- `~/NFLWork/hoop88_odds/`
- `~/NFLWork/bfa_odds/`

## How It Works

### 1. Build Consensus Lines
Fetches current odds from multiple books, weights by book sharpness, and calculates consensus spread/total for each game.

### 2. Generate Historical Samples
For each upcoming game, samples ~5% of historical games with similar spread/total profiles. Uses distance-weighted sampling based on how close historical games match current consensus.

### 3. Scrape Offshore Odds
Parallel scrapers fetch current odds from wagerzon, hoop88, and bfa aggregators.

### 4. Predict Derivative Markets
Uses the historical sample to calculate probabilities for:
- **Moneylines** - Q1, Q2, Q3, Q4, H1, H2 (2-way and 3-way)
- **Spreads** - Q1, Q2, Q3, Q4, H1, H2 + alternates
- **Totals** - Q1, Q2, Q3, Q4, H1, H2 + alternates
- **Team Totals** - Q1, Q2, Q3, Q4, H1, H2 + alternates

### 5. Find Edge
Compares model probabilities to odds from all books (API + offshore), calculates EV, and sizes bets using fractional Kelly criterion.

## Output

`nfl_bets_combined` table in duckdb with columns:
| Column | Description |
|--------|-------------|
| `market` | Market name (e.g., `h2h_q1`, `spreads_h1`) |
| `bet_on` | Team name or Over/Under/Tie |
| `line` | Spread or total line (if applicable) |
| `prob` | Model probability |
| `ev` | Expected value |
| `bet_size` | Kelly-sized bet amount |
| `odds` | American odds |
| `bookmaker_key` | Book for this bet (draftkings, wagerzon, hoop88, etc.) |

## Configuration

In `NFLPrepare.R`:
```r
bankroll   <- 100      # Bankroll for Kelly sizing
kelly_mult <- 0.25     # Fractional Kelly (25%)
N <- round(nrow(DT) * 0.05, 0)  # Sample size (~5% of historical data)
```

## Dashboard

See `NFL Dashboard/README.md` for the interactive betting dashboard that:
- Displays all +EV bets from the pipeline
- Tracks placed bets
- Warns about correlated market exposure
- Refreshes with one click (runs full pipeline)

## Supporting Scripts

| Script | Purpose |
|--------|---------|
| `Acquire New NFL Data.R` | Update play-by-play data from nflverse |
| `Acquire Historical Derivative Odds.R` | Fetch historical Q/H odds from The Odds API |
| `Historical Betting History 1999-2019.R` | Process pre-2020 betting data |
| `Consensus Betting History.R` | Build consensus lines for historical games |
| `All_Quarters_Backtest.R` | Backtest derivative market predictions |
