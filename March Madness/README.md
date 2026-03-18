# March Madness

NCAA tournament simulator using Bayesian power ratings and Monte Carlo simulation.

## What It Does

- Fetches bracket structure and results from ESPN API
- Pulls power ratings from 4 sources (BPI, KenPom, Torvik, EvanMiya)
- Simulates 10,000 tournament outcomes via residual sampling
- Outputs championship probabilities and Survivor Value rankings

## Key Files

| File | Purpose |
|------|---------|
| `shared.R` | Central utilities: team name resolution, rating fetchers, `simulate_game()` |
| `Basketball Model.R` | Pre-tournament full bracket simulation (68 teams, all 63 games) |
| `Dynamic Simulator.R` | Mid-tournament simulation with actual results (remaining games only) |
| `End of Round Simulator.R` | End-of-round specific simulation |
| `espn_bracket.R` | ESPN API scraper for bracket structure and results |

## Game Simulation Method

```
predicted_margin = rating_A - rating_B
actual_margin = predicted_margin + residual (sampled from historical distribution)
winner = team with positive margin
```

Residual sampling uses actual historical variance rather than assumed distributions.

## Survivor Value

For survivor pool strategy:
```
Survivor Value = P(win current round) × (1 - P(bust in future rounds))
```

Ranks teams by combined safety + longevity, not just win probability.

## Team Name Resolution

Centralized `TEAM_NAME_FIXES` mapping in `shared.R` covers all 68 tournament teams across all rating sources. Uses hoopR's `get_standard_team()` as fallback.

## Usage

```r
# Pre-tournament: full bracket simulation
source("Basketball Model.R")

# Mid-tournament: update with actual results
source("Dynamic Simulator.R")
```
