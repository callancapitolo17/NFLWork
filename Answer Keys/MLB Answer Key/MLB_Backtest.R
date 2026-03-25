# MLB Answer Key Backtest
# Leave-one-out backtest of F5 derivative markets against historical data
# Uses run_derivative_backtest() from Tools.R

setwd("~/NFLWork/Answer Keys")
suppressPackageStartupMessages({
  library(data.table)
  library(duckdb)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(DBI)
})
source("Tools.R")

cat("=== MLB BACKTEST ===\n")

# Load historical data (same as MLB.R Phase 1)
con <- dbConnect(duckdb(), dbdir = "pbp.duckdb", read_only = TRUE)

DT <- dbGetQuery(con, "
  SELECT *
  FROM mlb_betting_pbp
  WHERE game_home_margin_inning_inning_5 IS NOT NULL
") %>%
  rename(
    game_home_margin_period_F5 = game_home_margin_inning_inning_5,
    game_total_period_F5       = game_total_inning_inning_5,
    game_home_margin_period_FG = home_margin,
    game_total_period_FG       = total_final_score,
    home_ml_odds  = consensus_devig_home_odds,
    away_ml_odds  = consensus_devig_away_odds,
    over_odds     = consensus_devig_over_odds,
    under_odds    = consensus_devig_under_odds,
    actual_cover  = home_winner
  ) %>%
  mutate(
    actual_over = ifelse(game_total_period_FG > total_line, 1, 0),
    # Backtest needs a game_id column
    game_id = game_pk
  )

dbDisconnect(con)

cat(sprintf("Loaded %d games for backtest.\n", nrow(DT)))

# Get MLB-specific config
config <- get_sport_backtest_config("mlb")

# Run backtest (sample 500 games for speed, or all for full test)
args <- commandArgs(trailingOnly = TRUE)
n_test <- if (length(args) > 0) as.integer(args[1]) else 500

set.seed(42)
test_ids <- sample(unique(DT$game_id), min(n_test, length(unique(DT$game_id))))

cat(sprintf("Testing %d games (of %d total).\n", length(test_ids), length(unique(DT$game_id))))

results <- run_derivative_backtest(
  pbp_data     = DT,
  sport_config = config,
  sample_pct   = 0.02,
  test_game_ids = test_ids,
  verbose      = TRUE
)

cat("\n=== BACKTEST COMPLETE ===\n")
