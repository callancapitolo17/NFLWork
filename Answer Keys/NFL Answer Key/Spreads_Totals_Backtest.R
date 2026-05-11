# Spreads and Totals Backtest
# Tests Q1, Q2, Q3, Q4, H1, H2 spreads and totals markets

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
# HELPER FUNCTIONS
# =============================================================================

# devig_american() is sourced from Tools.R (probit-based, returns
# data.frame with p1/p2 columns).

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

# Market to period mapping for spreads and totals
market_to_period <- list(
  # Spreads
  "spreads_q1" = list(period = "1", margin_col = "game_home_margin_period_1"),
  "spreads_q2" = list(period = "2", margin_col = "game_home_margin_period_2"),
  "spreads_q3" = list(period = "3", margin_col = "game_home_margin_period_3"),
  "spreads_q4" = list(period = "4", margin_col = "game_home_margin_period_4"),
  "spreads_h1" = list(period = "Half1", margin_col = "game_home_margin_period_Half1"),
  "spreads_h2" = list(period = "Half2", margin_col = "game_home_margin_period_Half2"),
  # Alternate spreads
  "alternate_spreads_q1" = list(period = "1", margin_col = "game_home_margin_period_1"),
  "alternate_spreads_q2" = list(period = "2", margin_col = "game_home_margin_period_2"),
  "alternate_spreads_q3" = list(period = "3", margin_col = "game_home_margin_period_3"),
  "alternate_spreads_q4" = list(period = "4", margin_col = "game_home_margin_period_4"),
  "alternate_spreads_h1" = list(period = "Half1", margin_col = "game_home_margin_period_Half1"),
  "alternate_spreads_h2" = list(period = "Half2", margin_col = "game_home_margin_period_Half2"),
  # Totals
  "totals_q1" = list(period = "1", total_col = "game_total_period_1"),
  "totals_q2" = list(period = "2", total_col = "game_total_period_2"),
  "totals_q3" = list(period = "3", total_col = "game_total_period_3"),
  "totals_q4" = list(period = "4", total_col = "game_total_period_4"),
  "totals_h1" = list(period = "Half1", total_col = "game_total_period_Half1"),
  "totals_h2" = list(period = "Half2", total_col = "game_total_period_Half2"),
  # Alternate totals
  "alternate_totals_q1" = list(period = "1", total_col = "game_total_period_1"),
  "alternate_totals_q2" = list(period = "2", total_col = "game_total_period_2"),
  "alternate_totals_q3" = list(period = "3", total_col = "game_total_period_3"),
  "alternate_totals_q4" = list(period = "4", total_col = "game_total_period_4"),
  "alternate_totals_h1" = list(period = "Half1", total_col = "game_total_period_Half1"),
  "alternate_totals_h2" = list(period = "Half2", total_col = "game_total_period_Half2")
)

# =============================================================================
# LOAD DATA
# =============================================================================

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

# Filter to only spreads and totals
derivative_odds <- derivative_odds %>%
  filter(grepl("spread|total", market))

cat("\nSpreads/totals rows:\n")
print(derivative_odds %>% count(market))

disp <- compute_dispersion(DT_2020_plus, moneyline = FALSE)
ss <- disp$ss
st <- disp$st
N <- round(nrow(DT_all) * 0.10, 0)

# Get all test games
all_test_games <- unique(derivative_odds$game_id)
all_test_games <- all_test_games[all_test_games %in% DT_2020_plus$game_id]

cat("\n=========================================\n")
cat("SPREADS AND TOTALS BACKTEST\n")
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

  sample <- sample_result$sample

  # Process each market
  for (market_name in names(market_to_period)) {
    market_info <- market_to_period[[market_name]]
    period <- market_info$period

    # Get market odds
    market_odds <- game_odds %>% filter(market == !!market_name)
    if (nrow(market_odds) == 0) next

    # Determine if this is a spread or total market
    is_spread <- grepl("spread", market_name)
    is_alt <- grepl("alternate", market_name)

    if (is_spread) {
      # SPREADS: predict P(home margin > -spread)
      margin_col <- market_info$margin_col
      if (!(margin_col %in% names(sample))) next

      margins <- sample[[margin_col]]
      actual_margin <- test_game[[margin_col]]
      if (is.na(actual_margin)) next

      # Get all unique lines in this market
      unique_lines <- unique(market_odds$line)

      for (line in unique_lines) {
        line_odds <- market_odds %>% filter(line == !!line)
        if (nrow(line_odds) == 0) next

        best_home_odds <- max(line_odds$odds_home, na.rm = TRUE)
        best_away_odds <- max(line_odds$odds_away, na.rm = TRUE)

        if (is.infinite(best_home_odds) || is.infinite(best_away_odds)) next
        if (is.na(best_home_odds) || is.na(best_away_odds)) next

        # Calculate prediction: P(home margin > -line | not push)
        # For spread, line is from home's perspective (negative = favorite)
        # Home covers if margin > -line
        home_covers <- margins > -line
        pushes <- margins == -line
        non_pushes <- !pushes

        if (sum(non_pushes) == 0) next

        home_cover_prob <- mean(home_covers[non_pushes])
        away_cover_prob <- 1 - home_cover_prob

        # Actual outcome (exclude pushes in backtest)
        if (actual_margin == -line) next  # Push - skip
        actual_home_cover <- ifelse(actual_margin > -line, 1, 0)

        devigged <- devig_american(best_home_odds, best_away_odds)
        ev_home <- calc_ev(home_cover_prob, best_home_odds)
        ev_away <- calc_ev(away_cover_prob, best_away_odds)

        # Home cover bet
        all_predictions[[length(all_predictions) + 1]] <- list(
          game_id = test_game_id, market = market_name, market_type = ifelse(is_alt, "spreads_alt", "spreads"),
          bet_side = "home", line = line, period = period,
          algo_prob = home_cover_prob, book_prob = devigged$p1,
          book_odds = best_home_odds, ev = ev_home,
          actual_outcome = actual_home_cover
        )

        # Away cover bet
        all_predictions[[length(all_predictions) + 1]] <- list(
          game_id = test_game_id, market = market_name, market_type = ifelse(is_alt, "spreads_alt", "spreads"),
          bet_side = "away", line = line, period = period,
          algo_prob = away_cover_prob, book_prob = devigged$p2,
          book_odds = best_away_odds, ev = ev_away,
          actual_outcome = 1 - actual_home_cover
        )
      }

    } else {
      # TOTALS: predict P(total > line)
      total_col <- market_info$total_col
      if (!(total_col %in% names(sample))) next

      totals <- sample[[total_col]]
      actual_total <- test_game[[total_col]]
      if (is.na(actual_total)) next

      # Get all unique lines in this market
      unique_lines <- unique(market_odds$line)

      for (line in unique_lines) {
        line_odds <- market_odds %>% filter(line == !!line)
        if (nrow(line_odds) == 0) next

        best_over_odds <- max(line_odds$odds_over, na.rm = TRUE)
        best_under_odds <- max(line_odds$odds_under, na.rm = TRUE)

        if (is.infinite(best_over_odds) || is.infinite(best_under_odds)) next
        if (is.na(best_over_odds) || is.na(best_under_odds)) next

        # Calculate prediction: P(total > line | not push)
        overs <- totals > line
        pushes <- totals == line
        non_pushes <- !pushes

        if (sum(non_pushes) == 0) next

        over_prob <- mean(overs[non_pushes])
        under_prob <- 1 - over_prob

        # Actual outcome (exclude pushes in backtest)
        if (actual_total == line) next  # Push - skip
        actual_over <- ifelse(actual_total > line, 1, 0)

        devigged <- devig_american(best_over_odds, best_under_odds)
        ev_over <- calc_ev(over_prob, best_over_odds)
        ev_under <- calc_ev(under_prob, best_under_odds)

        # Over bet
        all_predictions[[length(all_predictions) + 1]] <- list(
          game_id = test_game_id, market = market_name, market_type = ifelse(is_alt, "totals_alt", "totals"),
          bet_side = "over", line = line, period = period,
          algo_prob = over_prob, book_prob = devigged$p1,
          book_odds = best_over_odds, ev = ev_over,
          actual_outcome = actual_over
        )

        # Under bet
        all_predictions[[length(all_predictions) + 1]] <- list(
          game_id = test_game_id, market = market_name, market_type = ifelse(is_alt, "totals_alt", "totals"),
          bet_side = "under", line = line, period = period,
          algo_prob = under_prob, book_prob = devigged$p2,
          book_odds = best_under_odds, ev = ev_under,
          actual_outcome = 1 - actual_over
        )
      }
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

# Summary by market
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
# SPREADS vs TOTALS COMPARISON
# =============================================================================

cat("\n=========================================\n")
cat("SPREADS vs TOTALS COMPARISON\n")
cat("=========================================\n\n")

# Compare performance by market type
type_comparison <- roi_by_market %>%
  mutate(broad_type = ifelse(grepl("spread", market_type), "spreads", "totals")) %>%
  group_by(broad_type) %>%
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

# Main vs Alternate comparison
alt_comparison <- roi_by_market %>%
  mutate(is_alt = grepl("alt", market_type)) %>%
  group_by(is_alt) %>%
  summarize(
    n_markets = n(),
    total_bets = sum(n_bets),
    total_staked = sum(total_staked),
    total_pnl = sum(total_pnl),
    weighted_roi = total_pnl / total_staked * 100,
    avg_roi = mean(roi),
    .groups = "drop"
  ) %>%
  mutate(line_type = ifelse(is_alt, "Alternate", "Main"))

cat("\nMain vs Alternate lines:\n")
print(alt_comparison)

cat("\n=========================================\n")
cat("BACKTEST COMPLETE\n")
cat("=========================================\n")
