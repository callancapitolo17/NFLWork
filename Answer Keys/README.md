# Answer Keys

Core pricing engine for sports betting edge detection. Generates statistically-backed fair odds by matching current game lines to historical outcomes.

## How It Works

The "Answer Key" algorithm:
1. For each upcoming game, compute the Pythagorean distance between its pregame line and every historical game
2. Select the N closest matches (default ~500 games)
3. Iteratively balance the sample so the mean spread/total converges to the target line (±1)
4. Price any derivative market (1H, quarters, team totals, alt lines) directly from sample outcomes
5. Compare fair price to posted odds across books → flag +EV opportunities

## Directory Structure

```
Answer Keys/
├── Tools.R                  # Shared utilities (devig, odds conversion, Rcpp sampling)
├── canonical_match.py       # Team name resolution (offshore → Odds API canonical)
├── run.py                   # Pipeline orchestrator (async parallel scrapers + R)
├── parlay.R                 # Fair parlay pricing from sample data
├── props.R                  # Fair prop pricing from play-by-play data
│
├── CBB Answer Key/          # College basketball (primary active model)
├── CBB Dashboard/           # CBB interactive dashboard (Flask + HTML)
├── NFL Answer Key/          # NFL model
├── NFL Dashboard/           # NFL interactive dashboard
├── MLB Answer Key/          # MLB consensus model
├── College Football Answer Key/  # CFB model
├── EPL Answer Key/          # Soccer model
│
├── diagnostics/             # Deep diagnostic tools
└── tests/                   # Test suite (answer key validation, Rcpp)
```

## Key Files

### Tools.R — Shared Utilities
- `odds_to_prob()` / `prob_to_odds()`: American ↔ implied probability
- `devig_american()` / `devig_american_3way()`: Remove vigorish
- `pick_consensus_line()`: Weighted consensus across books
- `resolve_offshore_teams()`: R-side team name standardization
- `run_answer_key_sample()`: Core sampling with optional Rcpp acceleration
- `pipeline_timer()`: Performance instrumentation

### canonical_match.py — Team Name Resolution
Two-layer matching:
1. **Dictionary lookup**: ESPN team data with stripped punctuation
2. **Game-level fallback**: Substring + mascot-free matching

Used by all offshore scrapers to map sportsbook team names → Odds API canonical names.

### run.py — Pipeline Orchestrator
- Fetches canonical game list from Odds API (cached if <2 hrs old)
- Launches all scrapers in parallel (async)
- Signals R via sentinel file (`.scrapers_done_<sport>`) when scrapers complete
- Reports per-scraper timing

### parlay.R — Parlay Pricing
```r
# CLI usage
Rscript parlay.R "1H Away spread" "FG Over total"
```
Computes joint probability across legs from historical sample data.

### props.R — Prop Pricing
```r
# Interactive R usage
prop_odds("pass_touchdown", "> 2.5", period = "Half1")
```
Joins samples to play-by-play data, evaluates conditions, returns fair odds.

## Data Storage

All pipeline data lives in DuckDB databases (not CSV files):
- `cbb.duckdb` — CBB pipeline data (odds, samples, results)
- `pbp.duckdb` — Play-by-play data across sports
- `cbb_dashboard.duckdb` — Dashboard state (placed bets, settings, CLV)
- `cbb_20XX.duckdb` — Historical CBB data by season

## Running the Pipeline

```bash
# Full CBB pipeline (scrapers + model + dashboard)
python run.py cbb

# Just the R answer key (after scrapers have run)
Rscript "CBB Answer Key/CBB.R"

# Dashboard only
cd "CBB Dashboard" && bash run.sh
```
