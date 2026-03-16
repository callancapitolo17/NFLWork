# NFLWork

A sports analytics repository focused on NFL quantitative analysis, sports betting modeling, and multi-sport simulation. Built primarily in R with interactive Shiny dashboards deployed on shinyapps.io.

## What's Inside

### NFL Analytics
| Script | Description |
|--------|-------------|
| `NFL Play by Play.R` | Offensive efficiency landscape - EPA vs success rate scatter plots with team logos |
| `EPA by Week.R` | Cumulative EPA tracking across weeks with interactive plotly visualization |
| `NFLEloModel.R` | Bayesian Elo rating system using RStan with home field advantage modeling |
| `Play Prediction.R` | Play-level feature engineering for prediction models (1999-2024 pbp) |
| `NFL Regression.R` | Logistic regression on fourth-down decision-making (2018-2022) |
| `49ers Analysis.R` | 49ers offensive EPA/play, success rates, CPOE by game segment |
| `49ers Breakdown.R` | Play-action, RPO, and screen efficiency comparison with FTN charting |
| `Home vs Road.R` | Home vs away offensive efficiency splits for 2023 season |
| `Redzon Analysis.R` | Red zone EPA/play vs non-red zone efficiency comparison |
| `Weekly Deep Dives.R` | Comprehensive weekly analysis with FTN charting + participation data |
| `Weekly Breakouts.Rmd` | R Markdown weekly reports with formatted tables |
| `Offseason Deep Dives.R` | Multi-team offseason exploratory analysis |
| `One Off Analysis.R` | Ad-hoc analytical explorations |

### Scouting & Team Profiling
| Script | Description |
|--------|-------------|
| `Scouting Report.R` | Offensive scouting profiles with 50+ metrics (2006-2023) |
| `New Scouting Report.R` | Enhanced scouting with FTN charting data for 2024 |
| `Scouting 2.0.R` | Latest scouting report iteration |
| `Similarity Index.R` | Cross-era team similarity using comprehensive metric profiles |
| `RadarPlot.R` / `RadarPlot2.0.R` | Spider chart team comparison visualizations |

### Sports Betting
| Script | Description |
|--------|-------------|
| `TeaserSGP.R` | Teaser and same-game parlay cover rate analysis |
| `Updated Teaser.R` | Advanced teaser modeling with logistic regression (2010-2024) |
| `SuperBowl Props.R` | Super Bowl prop modeling with odds conversion utilities |
| `Bets Analysis.R` | Personal betting ROI and expected value tracking |
| `Bovada Scraper.py` | Selenium-based Bovada bet history scraper |

### Multi-Sport
| Script | Description |
|--------|-------------|
| `MLB Simulator.R` | MLB season and series outcome simulator |
| `NBA Simulator.R` | NBA season and series outcome simulator |
| `WNBA Simulator.R` | WNBA tournament simulator |
| `WildCard Simulator.R` | NFL wildcard round simulator |
| `Updated Series Simulator.R` | Playoff series outcome simulator |
| `Soccer Plots.R` / `Soccer Pizza Plots.R` | EPL player analysis with understatr |
| `SDFC.R` | San Diego FC ticket sales geographic segmentation |
| `NFL Draft Analysis.R` | Historical draft position and age trends (1999-present) |
| `Longest Reception Model.R` | ML model (RF, XGBoost) predicting WR longest reception |
| `NFLBigDataBowl.R` | Data loading functions for NFL Big Data Bowl 2025 |

## Subdirectories

| Directory | Description |
|-----------|-------------|
| [Answer Keys/](Answer%20Keys/) | Betting consensus line systems for NFL, MLB, CFB, and EPL with odds de-vigging and historical validation |
| [Bayesian Betting/](Bayesian%20Betting/) | Bayesian MLB modeling using RStan and Retrosheet data |
| [EPAbyWeek/](EPAbyWeek/) | Shiny app for interactive cumulative EPA exploration by week |
| [EfficiencyLandscape/](EfficiencyLandscape/) | Shiny app for NFL offensive efficiency scatter plots |
| [JS Projects/](JS%20Projects/) | Formula 1 racing analysis with driver participation and odds data |
| [March Madness/](March%20Madness/) | NCAA tournament bracket simulator with Monte Carlo methods |
| [NFL Draft/](NFL%20Draft/) | 2025 mock draft big board consensus tracker |
| [OddsHarvester/](OddsHarvester/) | Odds scraping project (in progress) |
| [TeamSimilarity/](TeamSimilarity/) | Shiny app for finding historically similar NFL teams via DBSCAN clustering |

## Tech Stack

- **R** (primary) - tidyverse, nflfastR, ggplot2, plotly, shiny, gt, rstan, caret, xgboost
- **Python** - selenium for web scraping
- **Data sources** - nflfastR, baseballr, oddsapiR, FTN charting, Retrosheet, understatr

## Setup

1. Open `NFLWork.Rproj` in RStudio
2. Most scripts auto-install required packages
3. Shiny apps: `shiny::runApp("EPAbyWeek")`, `shiny::runApp("EfficiencyLandscape")`, `shiny::runApp("TeamSimilarity")`
4. Python scraper: requires `selenium` and ChromeDriver
