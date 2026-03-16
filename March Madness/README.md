# March Madness

NCAA tournament bracket simulators using Monte Carlo methods. Scrapes bracket data from ncaa.com and runs thousands of simulations to project tournament outcomes.

## Files

| File | Description |
|------|-------------|
| `Basketball Model.R` | Main tournament simulator. Scrapes bracket from ncaa.com, loads power ratings, handles First Four play-in games, runs Monte Carlo simulations for full bracket projections. |
| `Dynamic Simulator.R` | Alternative simulation approach with dynamic updating |
| `End of Round Simulator.R` | Simulator that picks up from a specific round |
| `team_name_mapping.csv` | Standardized team name mapping for joining data sources |

## Dependencies

- `rvest`, `httr`, `jsonlite`, `googlesheets4`, `tidyverse`
