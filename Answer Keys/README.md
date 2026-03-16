# Answer Keys

Betting consensus line systems for multiple sports. These scripts build "answer keys" - historical records of consensus betting lines derived by de-vigging odds from multiple sportsbooks.

## Core Utilities

| File | Description |
|------|-------------|
| `Tools.R` | Shared utility functions: `odds_to_prob`, `devig_american`, `american_prob`, `logit`/`invlogit`, `logloss`, `pick_consensus_line`. Foundation for all answer key work. |

## NFL Answer Key

| File | Description |
|------|-------------|
| `NFLAnswerKey.R` | Loads 1999-2024 NFL schedules with historical odds. De-vigs spread and over/under to true probabilities. Implements distance index for consensus line building. |
| `NFLAnswerKey2.0.R` | Enhanced version with improved consensus logic |
| `Acquire New NFL Data.R` | Pulls live NFL odds via oddsapiR, stores in duckdb, uses fuzzyjoin for team matching |
| `Consensus Betting History.R` | Constructs historical betting records from consensus lines |

## MLB Answer Key

| File | Description |
|------|-------------|
| `MLB Answer Key.R` | Loads 2015-2024 MLB historical odds via oddsapiR. De-vigs moneylines from multiple books. |
| `MLB Odds API Scraping.R` | Scrapes and processes MLB betting odds via OddsAPI |
| `MLB API Data Analysis.R` | Pitch-level MLB data analysis |
| `MLB Updated API Pitch Data.R` | Pulls pitch-level data via baseballr `get_pbp_mlb()` (2015-2024), saves as RDS |
| `MLB Flat Odds.csv` | Large historical MLB odds database (32.8MB) |

## College Football Answer Key

| File | Description |
|------|-------------|
| `CFB Answer Key.R` | Uses cfbfastR to query CFBD API for college football data |

## EPL Answer Key

| File | Description |
|------|-------------|
| `EPL Answer Key.R` | Consolidates per-season EPL CSV files into `PremBettingHistory.csv`. Parses 1X2 markets (home/draw/away) from bookmaker JSON. |
| `PremBettingHistory.csv` | Consolidated EPL betting history (1.2MB) |
| `epl_XXXX-XXXX.csv` | Per-season EPL match + odds data (9 seasons, 2015-16 through 2024-25) |

## Dependencies

- `oddsapiR` (requires API key), `duckdb`, `fuzzyjoin`, `baseballr`, `cfbfastR`, `tidyverse`
