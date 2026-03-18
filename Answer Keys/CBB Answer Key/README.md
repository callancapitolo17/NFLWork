# CBB Answer Key

Generates derivative market predictions (halves, team totals) for college basketball games by building consensus lines and sampling from historical data. Compares predictions against multiple bookmakers (API + offshore scrapers). Uses conditional Kelly criterion to size bets while accounting for already-placed positions.

## Quick Start

```bash
# Run the full pipeline (scrapers + predictions)
cd ~/NFLWork
python "Answer Keys/run.py" cbb

# Or use the dashboard
cd "Answer Keys/CBB Dashboard"
./run.sh
```

## Architecture

The pipeline runs in three phases orchestrated by `run.py`:

### Phase 1: Sharp Scrapers First
- **bookmaker, bet105** — Run first (in parallel with each other). Must complete before R starts Phase 2 so sharp data is available for consensus.

### Phase 2: Parallel Execution
Runs simultaneously after sharp scrapers finish:
- **Other scrapers** — wagerzon, hoop88, bfa, kalshi (fetch rec/offshore odds)
- **CBB.R Phases 1-5** — Load data, build sharp consensus, generate samples, build API predictions

### Phase 3: Combine + Size
- **CBB.R Phases 6-7.5** — Wait for scrapers, load offshore odds, find +EV bets, apply correlation-adjusted Kelly sizing

```
run.py cbb
    │
    ├── [first]    bookmaker scraper → bookmaker.duckdb  (sharp)
    ├── [first]    bet105 scraper    → bet105.duckdb     (sharp)
    │
    ├── [parallel] wagerzon scraper  → wagerzon.duckdb
    ├── [parallel] hoop88 scraper    → hoop88.duckdb
    ├── [parallel] bfa scraper       → bfa.duckdb
    ├── [parallel] kalshi scraper    → kalshi.duckdb
    ├── [parallel] CBB.R (Phases 1-5)
    │       ├── Phase 1: Load historical PBP data
    │       ├── Phase 2: Fetch odds, load sharp scrapers, build sharp-only consensus
    │       ├── Phase 3: Build team name dictionary
    │       ├── Phase 4: Generate samples + fetch derivative odds
    │       └── Phase 5: Build predictions (API bets)
    │
    └── [sequential] CBB.R (Phases 6-7.5)
            ├── Phase 6: Wait for scrapers, compare to offshore
            ├── Phase 7: Combine all bets, deduplicate
            └── Phase 7.5: Correlation-adjusted Kelly sizing
```

### Consensus Strategy

**Historical (Phase 1):** All-book median, 0.5 prob hardcode. Used for sample building from PBP data.

**Live (Phase 2):** Sharp books only. Defined in `SHARP_BOOKS` (`Tools.R`):

| Book | Weight | Source |
|------|--------|--------|
| Pinnacle | 1.1 | Odds API |
| Bookmaker.eu | 1.1 | Scraper |
| LowVig | 1.0 | Odds API |
| Circa | 1.0 | Odds API |
| Bet105 | 1.0 | Scraper |

Rec books (DraftKings, FanDuel, etc.) get weight=0 — present in data but excluded from consensus. Games with no sharp book coverage are dropped. Fallback: if zero sharp books exist in the entire dataset, all books weighted equally.

## Correlation-Adjusted Kelly Sizing (Phase 7.5)

### Problem
Standard Kelly sizes each bet independently, but correlated bets (e.g., Under +75 and Team Tot Under +37 in the same game) compound exposure. Worse, when bets are already placed, standard Kelly recommends from-scratch fractions that ignore existing exposure.

### Solution: Conditional Kelly

Phase 7.5 uses `adjust_kelly_for_correlation()` (defined in `Tools.R`) to handle two cases:

**Case 1: No placed bets in correlation group**
Standard multivariate Kelly: `f* = Sigma^{-1} * mu`

**Case 2: Placed bets exist in correlation group (Conditional Kelly)**
```
f_new* = Sigma_nn^{-1} * (mu_new - Sigma_np * f_placed)
                                    ^^^^^^^^^^^^^^^^
                                    penalty term
```

The penalty term reduces new recommendations proportionally to:
- How much money is already wagered (`f_placed`)
- How correlated the new bet is with placed positions (`Sigma_np`)

**Examples:**
- Placed Under +75.5 ($132), pipeline recommends Under +75 → ρ ≈ 0.99 → penalty ≈ full → **$0**
- Placed Under +75.5 ($132), pipeline recommends Alt Spread +1 → ρ ≈ 0.3 → small penalty → **$128**
- No placed bets → penalty = 0 → standard Kelly (unchanged)

### Regularization
Ridge regularization (`Sigma + 0.01 * I`) prevents matrix singularity from near-duplicate lines (e.g., +75.5 vs +75 have ρ ≈ 0.99).

### Fallback
If the multivariate solve fails, a per-bet fallback scales each bet by average pairwise correlation. Bets with ρ > 0.90 to any placed bet are set to $0.

## Files

| File | Purpose |
|------|---------|
| `CBB.R` | Main pipeline: consensus → samples → predictions → offshore compare → Kelly sizing |
| `Acquire CBB Data.R` | Update play-by-play data from hoopR |
| `Acquire CBB Derivative Odds v2.R` | Fetch historical derivative odds from The Odds API |
| `build_betting_pbp.R` | Process raw PBP into betting-ready format |
| `CBB_Backtest.R` | Backtest derivative market predictions |
| `CBB_Backtest_Correlation.R` | Backtest correlation adjustment impact |
| `CBB_Parameter_Sweep.R` | Sweep Kelly multiplier / bankroll parameters |
| `run_cbb_daily.sh` | Cron wrapper for daily pipeline runs |

## Prerequisites

**Database:** `cbb.duckdb` in the Answer Keys directory containing:
- `cbb_betting_pbp` - Play-by-play betting data
- `cbb_closing_odds` - Historical closing odds (used for historical consensus)

**Dashboard DB:** `cbb_dashboard.duckdb` containing:
- `placed_bets` - Already-placed bets (read during Phase 7.5 for conditional Kelly)
- `book_settings` - Which books to show
- `sizing_settings` - Bankroll and Kelly multiplier

**API:** The Odds API key set as environment variable `ODDS_API_KEY`

**Scrapers:** Python venvs set up in:
- `~/NFLWork/bookmaker_odds/` (sharp — runs before R)
- `~/NFLWork/bet105_odds/` (sharp — runs before R)
- `~/NFLWork/wagerzon_odds/`
- `~/NFLWork/hoop88_odds/`
- `~/NFLWork/bfa_odds/`
- `~/NFLWork/kalshi_odds/`

## How It Works

### 1. Build Consensus Lines
Fetches current odds from Odds API + sharp scraper DuckDBs. Only sharp books (Pinnacle, Bookmaker, LowVig, Circa, Bet105) participate in consensus — rec books are zeroed out. `scraper_to_odds_api_format()` converts scraper data to Odds API format, then `pick_consensus_line()` picks the line with highest total sharp weight.

### 2. Generate Historical Samples
For each upcoming game, samples historical games with similar spread/total profiles. Uses distance-weighted sampling based on how close historical games match current consensus.

### 3. Predict Derivative Markets
Uses the historical sample to calculate probabilities for:
- **Spreads** - H1 + alternates
- **Totals** - H1 + alternates
- **Team Totals** - H1 home/away + alternates
- **Moneylines** - H1

### 4. Find Edge
Compares model probabilities to odds from all books (API + offshore), calculates EV, and sizes bets using fractional Kelly criterion.

### 5. Correlation-Adjusted Sizing
Groups bets by game, builds covariance matrices from historical samples, and applies conditional Kelly when placed bets exist. See [Correlation-Adjusted Kelly Sizing](#correlation-adjusted-kelly-sizing-phase-75) above.

## Output

`cbb_bets_combined` table in `cbb.duckdb` with columns:

| Column | Description |
|--------|-------------|
| `market` | Market name (e.g., `spreads_h1`, `totals_h1`) |
| `bet_on` | Team name or Over/Under |
| `line` | Spread or total line |
| `prob` | Model probability |
| `ev` | Expected value |
| `bet_size` | Kelly-sized bet amount (post-correlation adjustment) |
| `correlation_adj` | Ratio of adjusted to original bet size |
| `odds` | American odds |
| `bookmaker_key` | Book for this bet |

## Dashboard

See `CBB Dashboard/` for the interactive betting dashboard that:
- Displays all +EV bets from the pipeline
- Tracks placed bets with fill status (placed/partial/open)
- Applies conditional Kelly to account for existing exposure
- Shows correlation adjustment tooltips
- Refreshes with one click (runs full pipeline)
