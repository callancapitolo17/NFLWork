# JS Projects

Formula 1 racing analysis. Despite the directory name, this contains R-based F1 analytics (originally for a job trial project).

## Files

| File | Description |
|------|-------------|
| `ALT Sports Trial.R` | F1 driver analysis - loads race_info, race_sessions, session_competitors data. Identifies regular drivers (75%+ participation), builds gt() tables with 538 theme. |
| `Callan Capitolo ALT Sports Analysis.Rmd` | R Markdown version of the F1 analysis for presentation |

## Data Files

| File | Description |
|------|-------------|
| `f1_2024.csv` | 2024 F1 season data |
| `f1_2025.csv` | 2025 F1 season data |
| `f1_odds_24.csv` | 2024 F1 betting odds |
| `race_info_2024.csv` | Race metadata |
| `race_sessions_2024.csv` | Session-level data |
| `session_competitors_2024.csv` | Competitor session results (12.6MB) |
| `start_list.csv` | Race start list |

## Dependencies

- `tidyverse`, `gt`, `gtExtras`
