# All Quarters Backtest
# Tests Q1, Q2, Q3, Q4, H1 markets for both 2-way and 3-way moneylines

setwd("~/NFLWork/Answer Keys")
library(data.table)
library(duckdb)
library(dplyr)
library(tidyverse)
library(DBI)
source("Tools.R")

KELLY_MULT <- 0.25
BANKROLL <- 1000

# =============================================================================
# 3-WAY HELPER FUNCTIONS
# =============================================================================

# devig_american_3way() is sourced from Tools.R (probit-based, returns
# data.frame with p_home/p_away/p_tie columns).

# Compute EV for 3-way (same formula, just accepts American odds)
calc_ev_3way <- function(pred_prob, book_odds) {
  decimal_odds <- ifelse(book_odds > 0, 1 + book_odds / 100, 1 + 100 / abs(book_odds))
  pred_prob * decimal_odds - 1
}

# Load data
cat("Loading data...\n")
con <- dbConnect(duckdb(), dbdir = "pbp.duckdb")

DT_2020_plus <- dbGetQuery(con, "SELECT * FROM nfl_betting_pbp") %>%
  rename(
    home_spread_odds = "consensus_devig_home_odds",
    away_spread_odds = "consensus_devig_away_odds",
    over_odds = "consensus_devig_over_odds",
    under_odds = "consensus_devig_under_odds"
  ) %>%
  mutate(
    actual_over = ifelse(total_final_score > total_line, 1, 0),
    actual_cover = ifelse(home_margin > -home_spread, 1, 0)
  ) %>%
  as.data.table()

pre_20 <- dbGetQuery(con, "SELECT * FROM nfl_pre_20_betting_history") %>%
  rename(home_spread = "spread") %>%
  mutate(
    actual_over = ifelse(total_final_score > total_line, 1, 0),
    actual_cover = ifelse(home_margin > -home_spread, 1, 0)
  ) %>%
  as.data.table()

DT_all <- bind_rows(DT_2020_plus, pre_20 %>% select(any_of(names(DT_2020_plus)))) %>%
  as.data.table()

derivative_odds <- dbGetQuery(con, "SELECT * FROM nfl_derivative_closing_odds")
dbDisconnect(con)

cat("Markets available:", paste(unique(derivative_odds$market), collapse = ", "), "\n")
cat("Rows per market:\n")
print(derivative_odds %>% count(market))

disp <- compute_dispersion(DT_2020_plus, moneyline = FALSE)
ss <- disp$ss
st <- disp$st
N <- round(nrow(DT_all) * 0.10, 0)

# Helper functions
# devig_american() (2-way pair) is sourced from Tools.R (probit-based,
# returns data.frame with p1/p2 columns).

calc_ev <- function(pred_prob, book_odds) {
  decimal_odds <- ifelse(book_odds > 0, 1 + book_odds / 100, 1 + 100 / abs(book_odds))
  pred_prob * decimal_odds - 1
}

kelly_stake <- function(ev, book_odds, bankroll, kelly_mult = 0.25) {
  if (ev <= 0) return(0)
  decimal_odds <- ifelse(book_odds > 0, 1 + book_odds / 100, 1 + 100 / abs(book_odds))
  b <- decimal_odds - 1
  p <- (ev + 1) / decimal_odds
  q <- 1 - p
  kelly_full <- max(0, (b * p - q) / b)
  bankroll * kelly_mult * kelly_full
}

# Market to period column mapping
# 2-way markets use the standard column names
# 3-way markets use the same period columns but need separate probability columns
market_to_period <- list(
  # 2-way markets
  "h2h_q1" = list(period = "1", col = "game_home_margin_period_1", type = "2way"),
  "h2h_q2" = list(period = "2", col = "game_home_margin_period_2", type = "2way"),
  "h2h_q3" = list(period = "3", col = "game_home_margin_period_3", type = "2way"),
  "h2h_q4" = list(period = "4", col = "game_home_margin_period_4", type = "2way"),
  "h2h_h1" = list(period = "Half1", col = "game_home_margin_period_Half1", type = "2way"),
  # 3-way markets (same periods, different prediction columns)
  "h2h_3_way_q1" = list(period = "1", col = "game_home_margin_period_1", type = "3way"),
  "h2h_3_way_q2" = list(period = "2", col = "game_home_margin_period_2", type = "3way"),
  "h2h_3_way_q3" = list(period = "3", col = "game_home_margin_period_3", type = "3way"),
  "h2h_3_way_q4" = list(period = "4", col = "game_home_margin_period_4", type = "3way"),
  "h2h_3_way_h1" = list(period = "Half1", col = "game_home_margin_period_Half1", type = "3way")
)

# Get all test games
all_test_games <- unique(derivative_odds$game_id)
all_test_games <- all_test_games[all_test_games %in% DT_2020_plus$game_id]

cat("\n=========================================\n")
cat("ALL QUARTERS BACKTEST\n")
cat("=========================================\n\n")

cat("Test games:", length(all_test_games), "\n")
cat("Sample size N:", N, "\n\n")

cat("Running predictions for all games...\n")

all_predictions <- list()

for (i in seq_along(all_test_games)) {
  if (i %% 50 == 0) cat("  ", i, "/", length(all_test_games), "\n")

  test_game_id <- all_test_games[i]
  test_game <- DT_2020_plus[DT_2020_plus$game_id == test_game_id, ][1, ]

  game_odds <- derivative_odds %>% filter(game_id == test_game_id)
  if (nrow(game_odds) == 0) next

  DT_train <- DT_all[DT_all$game_id != test_game_id, ]

  parent_spread <- test_game$home_spread
  parent_total <- test_game$total_line
  target_cover <- test_game$home_spread_odds
  target_over <- test_game$over_odds

  if (is.na(parent_spread) || is.na(parent_total) || is.na(target_cover) || is.na(target_over)) next

  # Step 1: Run Answer Key to get balanced sample (once per game)
  sample_result <- tryCatch({
    run_answer_key_sample(
      id = test_game_id, parent_spread = parent_spread, parent_total = parent_total,
      target_cover = target_cover, target_over = target_over,
      DT = DT_train, ss = ss, st = st, N = N,
      use_spread_line = TRUE
    )
  }, error = function(e) NULL)

  if (is.null(sample_result)) next

  # Step 2: Generate moneyline predictions from the sample
  predictions <- predict_moneyline_from_sample(
    sample_result$sample,
    margin_col = "game_home_margin_period"
  )

  # Process each market
  for (market in names(market_to_period)) {
    market_info <- market_to_period[[market]]
    period_suffix <- market_info$period
    pred_col <- market_info$col
    market_type <- market_info$type

    # Get market odds
    market_odds <- game_odds %>% filter(market == !!market)
    if (nrow(market_odds) == 0) next

    # Get actual margin for this period
    actual_margin <- test_game[[pred_col]]
    if (is.na(actual_margin)) next

    if (market_type == "2way") {
      # 2-way: exclude ties
      if (actual_margin == 0) next

      algo_col <- paste0("game_home_margin_period_", period_suffix)
      if (!(algo_col %in% names(predictions))) next

      pred_prob <- predictions[[algo_col]]
      if (is.na(pred_prob)) next

      actual_outcome <- ifelse(actual_margin > 0, 1, 0)

      best_home_odds <- max(market_odds$odds_home, na.rm = TRUE)
      best_away_odds <- max(market_odds$odds_away, na.rm = TRUE)
      if (is.infinite(best_home_odds) || is.infinite(best_away_odds)) next

      devigged <- devig_american(best_home_odds, best_away_odds)
      ev_home <- calc_ev(pred_prob, best_home_odds)
      ev_away <- calc_ev(1 - pred_prob, best_away_odds)

      # Home bet
      all_predictions[[length(all_predictions) + 1]] <- list(
        game_id = test_game_id, market = market, market_type = "2way",
        bet_side = "home", spread = parent_spread,
        algo_prob = pred_prob, book_prob = devigged$p1,
        book_odds = best_home_odds, ev = ev_home,
        actual_outcome = actual_outcome
      )

      # Away bet
      all_predictions[[length(all_predictions) + 1]] <- list(
        game_id = test_game_id, market = market, market_type = "2way",
        bet_side = "away", spread = parent_spread,
        algo_prob = 1 - pred_prob, book_prob = devigged$p2,
        book_odds = best_away_odds, ev = ev_away,
        actual_outcome = 1 - actual_outcome
      )

    } else if (market_type == "3way") {
      # 3-way: include ties as a distinct outcome
      algo_col_home <- paste0("game_home_margin_period_", period_suffix, "_3way_home")
      algo_col_away <- paste0("game_home_margin_period_", period_suffix, "_3way_away")
      algo_col_tie <- paste0("game_home_margin_period_", period_suffix, "_3way_tie")

      if (!(algo_col_home %in% names(predictions))) next

      pred_prob_home <- predictions[[algo_col_home]]
      pred_prob_away <- predictions[[algo_col_away]]
      pred_prob_tie <- predictions[[algo_col_tie]]

      if (is.na(pred_prob_home) || is.na(pred_prob_away) || is.na(pred_prob_tie)) next

      # Actual outcome: 1 = home win, 0 = away win, 2 = tie
      actual_outcome_home <- ifelse(actual_margin > 0, 1, 0)
      actual_outcome_away <- ifelse(actual_margin < 0, 1, 0)
      actual_outcome_tie <- ifelse(actual_margin == 0, 1, 0)

      # Get odds (3-way markets should have odds_tie column)
      best_home_odds <- max(market_odds$odds_home, na.rm = TRUE)
      best_away_odds <- max(market_odds$odds_away, na.rm = TRUE)
      best_tie_odds <- if ("odds_tie" %in% names(market_odds)) max(market_odds$odds_tie, na.rm = TRUE) else NA

      if (is.infinite(best_home_odds) || is.infinite(best_away_odds) || is.na(best_tie_odds) || is.infinite(best_tie_odds)) next

      devigged <- devig_american_3way(best_home_odds, best_away_odds, best_tie_odds)
      ev_home <- calc_ev_3way(pred_prob_home, best_home_odds)
      ev_away <- calc_ev_3way(pred_prob_away, best_away_odds)
      ev_tie <- calc_ev_3way(pred_prob_tie, best_tie_odds)

      # Home bet
      all_predictions[[length(all_predictions) + 1]] <- list(
        game_id = test_game_id, market = market, market_type = "3way",
        bet_side = "home", spread = parent_spread,
        algo_prob = pred_prob_home, book_prob = devigged$p_home,
        book_odds = best_home_odds, ev = ev_home,
        actual_outcome = actual_outcome_home
      )

      # Away bet
      all_predictions[[length(all_predictions) + 1]] <- list(
        game_id = test_game_id, market = market, market_type = "3way",
        bet_side = "away", spread = parent_spread,
        algo_prob = pred_prob_away, book_prob = devigged$p_away,
        book_odds = best_away_odds, ev = ev_away,
        actual_outcome = actual_outcome_away
      )

      # Tie bet
      all_predictions[[length(all_predictions) + 1]] <- list(
        game_id = test_game_id, market = market, market_type = "3way",
        bet_side = "tie", spread = parent_spread,
        algo_prob = pred_prob_tie, book_prob = devigged$p_tie,
        book_odds = best_tie_odds, ev = ev_tie,
        actual_outcome = actual_outcome_tie
      )
    }
  }
}

results <- map_dfr(all_predictions, as_tibble)
cat("\nTotal predictions:", nrow(results), "\n\n")

# =============================================================================
# RESULTS BY MARKET
# =============================================================================

cat("=========================================\n")
cat("RESULTS BY MARKET\n")
cat("=========================================\n\n")

# Calculate log-loss for each prediction
results <- results %>%
  mutate(
    algo_logloss = -(actual_outcome * log(pmax(algo_prob, 1e-10)) + (1 - actual_outcome) * log(pmax(1 - algo_prob, 1e-10))),
    book_logloss = -(actual_outcome * log(pmax(book_prob, 1e-10)) + (1 - actual_outcome) * log(pmax(1 - book_prob, 1e-10)))
  )

# Summary by market and market_type
by_market <- results %>%
  group_by(market, market_type) %>%
  summarize(
    n = n(),
    n_games = n_distinct(game_id),
    algo_logloss = mean(algo_logloss, na.rm = TRUE),
    book_logloss = mean(book_logloss, na.rm = TRUE),
    algo_better = algo_logloss < book_logloss,
    logloss_diff = book_logloss - algo_logloss,
    avg_algo_prob = mean(algo_prob),
    avg_book_prob = mean(book_prob),
    avg_actual = mean(actual_outcome),
    algo_cal_error = mean(algo_prob) - mean(actual_outcome),
    book_cal_error = mean(book_prob) - mean(actual_outcome),
    .groups = "drop"
  ) %>%
  arrange(market_type, market)

cat("Log-loss comparison (lower is better):\n\n")
print(by_market %>% select(market, market_type, n, n_games, algo_logloss, book_logloss, algo_better, logloss_diff))

cat("\n\nCalibration (closer to 0 is better):\n\n")
print(by_market %>% select(market, market_type, avg_algo_prob, avg_book_prob, avg_actual, algo_cal_error, book_cal_error))

# =============================================================================
# ROI BY MARKET (EV > 0)
# =============================================================================

cat("\n=========================================\n")
cat("ROI BY MARKET (EV > 0 filter)\n")
cat("=========================================\n\n")

roi_by_market <- results %>%
  filter(ev > 0) %>%
  mutate(
    stake = map2_dbl(ev, book_odds, ~kelly_stake(.x, .y, BANKROLL, KELLY_MULT)),
    decimal_odds = ifelse(book_odds > 0, 1 + book_odds / 100, 1 + 100 / abs(book_odds)),
    pnl = ifelse(actual_outcome == 1, stake * (decimal_odds - 1), -stake)
  ) %>%
  group_by(market, market_type) %>%
  summarize(
    n_bets = n(),
    n_games = n_distinct(game_id),
    avg_ev = mean(ev) * 100,
    total_staked = sum(stake),
    total_pnl = sum(pnl),
    roi = total_pnl / total_staked * 100,
    win_rate = mean(actual_outcome) * 100,
    pred_win_rate = mean(algo_prob) * 100,
    cal_gap = win_rate - pred_win_rate,
    .groups = "drop"
  ) %>%
  arrange(market_type, market)

print(roi_by_market)

# =============================================================================
# OVERALL SUMMARY
# =============================================================================

cat("\n=========================================\n")
cat("OVERALL SUMMARY\n")
cat("=========================================\n\n")

# Total across all markets, by market type
overall_by_type <- results %>%
  group_by(market_type) %>%
  summarize(
    total_predictions = n(),
    total_games = n_distinct(game_id),
    algo_logloss = mean(algo_logloss, na.rm = TRUE),
    book_logloss = mean(book_logloss, na.rm = TRUE),
    algo_wins = algo_logloss < book_logloss,
    .groups = "drop"
  )

cat("Summary by market type:\n")
print(overall_by_type)

# Total across all markets
overall <- results %>%
  summarize(
    total_predictions = n(),
    total_games = n_distinct(game_id),
    algo_logloss = mean(algo_logloss, na.rm = TRUE),
    book_logloss = mean(book_logloss, na.rm = TRUE),
    algo_wins = algo_logloss < book_logloss
  )

cat("\nTotal predictions:", overall$total_predictions, "\n")
cat("Total games:", overall$total_games, "\n")
cat("Overall algo log-loss:", round(overall$algo_logloss, 4), "\n")
cat("Overall book log-loss:", round(overall$book_logloss, 4), "\n")
cat("Winner:", ifelse(overall$algo_wins, "ALGO", "BOOK"), "\n")

# ROI summary
ev_bets <- results %>% filter(ev > 0)
if (nrow(ev_bets) > 0) {
  ev_bets <- ev_bets %>%
    mutate(
      stake = map2_dbl(ev, book_odds, ~kelly_stake(.x, .y, BANKROLL, KELLY_MULT)),
      decimal_odds = ifelse(book_odds > 0, 1 + book_odds / 100, 1 + 100 / abs(book_odds)),
      pnl = ifelse(actual_outcome == 1, stake * (decimal_odds - 1), -stake)
    )

  cat("\nOverall ROI (EV > 0 bets):\n")
  cat("  Total bets:", nrow(ev_bets), "\n")
  cat("  Total staked: $", round(sum(ev_bets$stake), 0), "\n")
  cat("  Total P&L: $", round(sum(ev_bets$pnl), 0), "\n")
  cat("  ROI:", round(sum(ev_bets$pnl) / sum(ev_bets$stake) * 100, 1), "%\n")
}

# =============================================================================
# WHICH QUARTER IS MOST PREDICTABLE?
# =============================================================================

cat("\n=========================================\n")
cat("MARKET RANKINGS\n")
cat("=========================================\n\n")

cat("By log-loss improvement over book (algo - book, negative = algo better):\n")
print(by_market %>% arrange(logloss_diff) %>% select(market, market_type, logloss_diff, algo_better))

cat("\nBy ROI:\n")
print(roi_by_market %>% arrange(desc(roi)) %>% select(market, market_type, n_bets, roi, total_pnl))

# =============================================================================
# 2-WAY vs 3-WAY COMPARISON
# =============================================================================

cat("\n=========================================\n")
cat("2-WAY vs 3-WAY COMPARISON\n")
cat("=========================================\n\n")

# Compare performance by market type
type_comparison <- roi_by_market %>%
  group_by(market_type) %>%
  summarize(
    n_markets = n(),
    total_bets = sum(n_bets),
    total_staked = sum(total_staked),
    total_pnl = sum(total_pnl),
    weighted_roi = total_pnl / total_staked * 100,
    avg_roi = mean(roi),
    .groups = "drop"
  )

cat("Overall performance by market type:\n")
print(type_comparison)

cat("\n=========================================\n")
cat("BACKTEST COMPLETE\n")
cat("=========================================\n")
