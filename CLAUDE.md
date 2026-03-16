# CLAUDE.md - NFLWork Repository Guide

## README Maintenance Policy (CRITICAL)

**Every Claude Code session that modifies files MUST update the relevant README.md files before completing.**

- If you add, rename, or significantly modify a script, update the corresponding README.md.
- If you create a new subdirectory, create a README.md inside it describing its purpose and contents.
- If you change how a script works (new data sources, new outputs, new dependencies), reflect that in the README.
- The root README.md must always reflect the current state of the repository.
- Subdirectory READMEs must stay in sync with their contents.
- **Do not consider a task complete until all affected READMEs are up to date.**

## Project Overview

This is a sports analytics repository focused primarily on NFL quantitative analysis, sports betting modeling, and multi-sport simulation. The codebase is predominantly R with one Python scraper.

**Owner:** Callan Capitolo

## Repository Structure

```
NFLWork/
├── *.R / *.Rmd / *.py          # Root-level analysis scripts (see categories below)
├── Answer Keys/                 # Betting answer key systems (NFL, MLB, CFB, EPL)
├── Bayesian Betting/            # Bayesian MLB modeling with RStan
├── EPAbyWeek/                   # Shiny app: cumulative EPA tracking by week
├── EfficiencyLandscape/         # Shiny app: NFL offensive efficiency scatter plots
├── JS Projects/                 # Formula 1 racing analysis
├── March Madness/               # NCAA tournament bracket simulator
├── NFL Draft/                   # Mock draft consensus tracker and analysis
├── OddsHarvester/               # (Currently empty) Odds scraping project
├── TeamSimilarity/              # Shiny app: historical team similarity via clustering
├── transactions.csv             # Personal betting transaction history
├── Clean Win Totals - Sheet1.csv
└── AnalystDataUseCase.xlsx      # SDFC ticket sales data
```

## Script Categories

### NFL Core Analytics
- `NFL Play by Play.R` - Offensive efficiency landscape (EPA vs success rate)
- `EPA by Week.R` - Cumulative EPA tracking with plotly
- `NFLEloModel.R` - Bayesian Elo ratings via RStan
- `Play Prediction.R` - Play-level prediction features (1999-2024)
- `NFL Regression.R` - Fourth-down decision logistic regression
- `49ers Analysis.R` / `49ers Breakdown.R` - Team-specific deep dives
- `Home vs Road.R` - Home/away offensive efficiency splits
- `Redzon Analysis.R` - Red zone efficiency analysis
- `NFL Viz Tutorial.R` - Visualization tutorial

### Scouting & Team Profiling
- `Scouting Report.R` - Comprehensive offensive scouting (2006-2023, 50+ metrics)
- `New Scouting Report.R` - Enhanced scouting with FTN charting data (2024)
- `Scouting 2.0.R` - Latest scouting iteration
- `Similarity Index.R` - Cross-era team similarity profiles
- `RadarPlot.R` / `RadarPlot2.0.R` - Spider/radar chart team comparisons

### Sports Betting & Modeling
- `TeaserSGP.R` - Teaser and SGP cover rate analysis
- `Updated Teaser.R` - Advanced teaser modeling with logistic regression (2010-2024)
- `SuperBowl Props.R` - Super Bowl prop modeling with odds utilities
- `Bets Analysis.R` - Personal betting ROI/EV tracking
- `Bovada Scraper.py` - Selenium-based bet history scraper

### Weekly/Seasonal Analysis
- `Weekly Deep Dives.R` - Weekly NFL analysis with FTN + participation data
- `Weekly Breakouts.Rmd` - R Markdown weekly reports
- `Offseason Deep Dives.R` - Offseason exploratory analysis
- `One Off Analysis.R` - Ad-hoc analytical explorations

### Multi-Sport
- `MLB Simulator.R` / `NBA Simulator.R` / `WNBA Simulator.R` - Season/series simulators
- `WildCard Simulator.R` / `Updated Series Simulator.R` - Playoff simulators
- `Soccer Plots.R` / `Soccer Pizza Plots.R` - EPL player/team viz (understatr)
- `SDFC.R` - San Diego FC ticket sales analysis
- `NFL Draft Analysis.R` - Historical draft position trends
- `Longest Reception Model.R` - WR longest reception ML prediction (RF, XGBoost)
- `NFLBigDataBowl.R` - NFL Big Data Bowl 2025 data loading

## Key Technology Stack

| Category | Packages |
|----------|----------|
| NFL Data | nflfastR, nflreadr, nflplotR |
| MLB Data | baseballr, oddsapiR |
| CFB Data | cfbfastR |
| Visualization | ggplot2, plotly, gt, gtExtras, patchwork |
| Web Apps | shiny (3 deployed apps on shinyapps.io) |
| Modeling | rstan, caret, randomForest, xgboost |
| Data Wrangling | tidyverse, data.table, duckdb, fuzzyjoin |
| Web Scraping | rvest, httr, selenium (Python) |

## Shiny Apps (Deployed)

1. **EPAbyWeek** - Interactive cumulative EPA explorer (shinyapps.io/callan-capitolo)
2. **EfficiencyLandscape** - NFL efficiency scatter plots
3. **TeamSimilarity** - Historical team comparison tool with DBSCAN clustering

## Build & Run Notes

- R scripts are standalone; run in RStudio with `NFLWork.Rproj`
- Most scripts auto-install missing packages or assume tidyverse + nflfastR ecosystem
- Shiny apps can be run locally via `shiny::runApp()` in their directories
- `Bovada Scraper.py` requires Python + Selenium + ChromeDriver
- Answer Key scripts require an OddsAPI key for live data pulls
- Some scripts reference FTN charting data (requires FTN subscription)

## Data Sources

- **nflfastR** - NFL play-by-play data (1999-present)
- **baseballr** - MLB Statcast and game data
- **oddsapiR** - Historical and live betting odds
- **FTN Data** - Charting data (personnel, play design, coverage)
- **Retrosheet** - Historical MLB play-by-play
- **understatr** - EPL expected goals and shot data
- **NCAA** - Tournament bracket data (web scraped)

## Sensitive Files

- `transactions.csv` - Personal betting history (do not share publicly)
- API keys may appear in some scripts (OddsAPI, CFBD) - sanitize before sharing
- `.DS_Store` - macOS metadata (should be in .gitignore)
