# CBB Derivative Backtest
# Tests H1 and H2 spreads, totals, moneylines, and team totals

setwd("~/NFLWork/Answer Keys")
library(data.table)
library(duckdb)
library(dplyr)
library(tidyverse)
library(DBI)
source("Tools.R")

# =============================================================================
# LOAD AND PREPARE DATA
# =============================================================================

cat("Loading CBB data...\n")
con <- dbConnect(duckdb(), dbdir = "cbb.duckdb")

# Get PBP outcomes
pbp <- dbGetQuery(con, "SELECT * FROM cbb_betting_pbp")

# Get closing odds (take consensus across books)
closing_odds <- dbGetQuery(con, "
  SELECT
    id as odds_id,
    home_spread,
    total_line,
    spread_home_odds,
    spread_away_odds,
    tot_over_odds,
    tot_under_odds
  FROM cbb_closing_odds
  WHERE spread_home_odds IS NOT NULL
    AND spread_away_odds IS NOT NULL
    AND tot_over_odds IS NOT NULL
    AND tot_under_odds IS NOT NULL
")

dbDisconnect(con)

cat("PBP rows:", nrow(pbp), "\n")
cat("Closing odds rows:", nrow(closing_odds), "\n")

# =============================================================================
# CALCULATE CONSENSUS LINES AND PROBABILITIES
# =============================================================================

# Devig function
american_to_prob <- function(odds) {
  ifelse(odds > 0, 100 / (odds + 100), abs(odds) / (abs(odds) + 100))
}

# Calculate devigged probabilities for each book
closing_odds <- closing_odds %>%
  mutate(
    # Implied probs
    implied_home_cover = american_to_prob(spread_home_odds),
    implied_away_cover = american_to_prob(spread_away_odds),
    implied_over = american_to_prob(tot_over_odds),
    implied_under = american_to_prob(tot_under_odds),
    # Devigged probs
    spread_juice = implied_home_cover + implied_away_cover,
    total_juice = implied_over + implied_under,
    devig_home_cover = implied_home_cover / spread_juice,
    devig_away_cover = implied_away_cover / spread_juice,
    devig_over = implied_over / total_juice,
    devig_under = implied_under / total_juice
  )

# Take consensus (average) across all books for each game
consensus <- closing_odds %>%
  group_by(odds_id) %>%
  summarize(
    home_spread = median(home_spread, na.rm = TRUE),
    total_line = median(total_line, na.rm = TRUE),
    consensus_home_cover_prob = mean(devig_home_cover, na.rm = TRUE),
    consensus_over_prob = mean(devig_over, na.rm = TRUE),
    n_books = n(),
    .groups = "drop"
  )

cat("Consensus lines computed for", nrow(consensus), "games\n")

# Join PBP with consensus
DT <- pbp %>%
  inner_join(consensus, by = "odds_id") %>%
  filter(
    !is.na(home_spread),
    !is.na(total_line),
    !is.na(consensus_home_cover_prob),
    !is.na(consensus_over_prob),
    !is.na(game_home_margin_h1),
    !is.na(game_home_margin_h2)
  ) %>%
  # Add actual_cover and actual_over columns required by balance_sample
  mutate(
    # Home covers if margin > -spread
    actual_cover = ifelse(game_home_margin_fg > -home_spread, 1, 0),
    # Over hits if total > line
    actual_over = ifelse(game_total_fg > total_line, 1, 0)
  ) %>%
  as.data.table()

cat("Joined dataset rows:", nrow(DT), "\n\n")

# =============================================================================
# RUN BACKTEST
# =============================================================================

# Configuration - adjust these for speed vs accuracy tradeoff
N_TEST_GAMES <- 500  # Number of games to test (use NULL for all games)
SAMPLE_PCT <- 0.05   # Use 5% sample instead of 10% (faster)

# Randomly sample test game IDs if N_TEST_GAMES is set
if (!is.null(N_TEST_GAMES) && N_TEST_GAMES < nrow(DT)) {
  set.seed(42)  # For reproducibility
  all_game_ids <- unique(DT$game_id)
  test_ids <- sample(all_game_ids, min(N_TEST_GAMES, length(all_game_ids)))
  cat("Testing on random sample of", length(test_ids), "games (out of", nrow(DT), ")\n\n")
} else {
  test_ids <- NULL  # Test all games
  cat("Testing on all", nrow(DT), "games\n\n")
}

# Create CBB-specific config
cbb_config <- list(
  sport = "cbb",
  margin_col_prefix = "game_home_margin_",
  total_col_prefix = "game_total_",
  home_score_prefix = "home_",
  away_score_prefix = "away_",
  parent_spread_col = "home_spread",
  parent_total_col = "total_line",
  consensus_home_odds_col = "consensus_home_cover_prob",
  consensus_over_odds_col = "consensus_over_prob",
  periods = list(
    h1 = list(suffix = "h1", display = "H1"),
    h2 = list(suffix = "h2", display = "H2")
  )
)

# Run the backtest - uses FULL dataset for training, but only tests subset
results <- run_derivative_backtest(
  pbp_data = DT,
  sport_config = cbb_config,
  kelly_mult = 0.25,
  bankroll = 1000,
  sample_pct = SAMPLE_PCT,
  test_game_ids = test_ids,
  verbose = TRUE
)

# =============================================================================
# ADDITIONAL ANALYSIS
# =============================================================================

if (!is.null(results)) {
  cat("\n===========================================\n")
  cat("MARKET TYPE COMPARISON\n")
  cat("===========================================\n\n")

  type_summary <- results$by_market %>%
    group_by(market_type) %>%
    summarize(
      n_markets = n(),
      total_predictions = sum(n_predictions),
      avg_logloss = mean(algo_logloss),
      avg_cal_error = mean(calibration_error),
      avg_abs_cal_error = mean(abs_cal_error),
      .groups = "drop"
    )

  print(type_summary)

  cat("\n===========================================\n")
  cat("PERIOD COMPARISON\n")
  cat("===========================================\n\n")

  period_summary <- results$by_market %>%
    group_by(period) %>%
    summarize(
      n_markets = n(),
      total_predictions = sum(n_predictions),
      avg_logloss = mean(algo_logloss),
      avg_cal_error = mean(calibration_error),
      .groups = "drop"
    )

  print(period_summary)

  # =============================================================================
  # CALIBRATION PLOT
  # =============================================================================

  cat("\n===========================================\n")
  cat("CALIBRATION ANALYSIS\n")
  cat("===========================================\n\n")

  # Calculate calibration by probability bins
  cal_by_bin <- results$predictions %>%
    mutate(prob_bin = cut(algo_prob, breaks = seq(0, 1, 0.1), include.lowest = TRUE)) %>%
    group_by(prob_bin) %>%
    summarize(
      n = n(),
      avg_predicted = mean(algo_prob),
      avg_actual = mean(actual_outcome),
      std_error = sd(actual_outcome) / sqrt(n()),
      .groups = "drop"
    ) %>%
    filter(!is.na(prob_bin))

  cat("Calibration by probability bin:\n")
  print(cal_by_bin)

  # Perfect calibration = predicted equals actual
  cat("\nCalibration quality (slope should be ~1, intercept ~0):\n")
  cal_model <- lm(avg_actual ~ avg_predicted, data = cal_by_bin)
  cat("  Slope:", round(coef(cal_model)[2], 3), "\n")
  cat("  Intercept:", round(coef(cal_model)[1], 3), "\n")
  cat("  R-squared:", round(summary(cal_model)$r.squared, 3), "\n")

  # =============================================================================
  # SUMMARY BY MARKET TYPE AND SIDE
  # =============================================================================

  cat("\n===========================================\n")
  cat("DETAILED RESULTS BY BET SIDE\n")
  cat("===========================================\n\n")

  by_side <- results$predictions %>%
    group_by(market_type, bet_side) %>%
    summarize(
      n = n(),
      avg_predicted = mean(algo_prob),
      avg_actual = mean(actual_outcome),
      cal_error = avg_predicted - avg_actual,
      win_rate = mean(actual_outcome) * 100,
      .groups = "drop"
    )

  print(by_side)
}

cat("\n===========================================\n")
cat("BACKTEST COMPLETE\n")
cat("===========================================\n")
