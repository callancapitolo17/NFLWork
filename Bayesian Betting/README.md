# Bayesian Betting

Bayesian MLB outcome modeling using RStan and historical play-by-play data from Retrosheet.

## Files

| File | Description |
|------|-------------|
| `Bayesian MLB Betting.R` | Loads Retrosheet MLB data (all2023.csv), maps team abbreviations for all 30 MLB teams, merges with player rosters, and fits Bayesian models using RStan for game outcome prediction. |

## Dependencies

- `rstan`, `baseballr`, `tidyverse`
- Retrosheet play-by-play data files (e.g., `all2023.csv`)
