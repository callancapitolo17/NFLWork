# CBB Derivative Backtest v2
# Uses actual book derivative odds for proper algo vs book comparison

setwd("~/NFLWork/Answer Keys")
library(data.table)
library(duckdb)
library(dplyr)
library(tidyverse)
library(DBI)
source("Tools.R")

# =============================================================================
# CONFIGURATION
# =============================================================================

N_TEST_GAMES <- 1000  # Number of games to test (NULL for all)
SAMPLE_PCT <- 0.05    # 5% sample for speed
KELLY_MULT <- 0.25
BANKROLL <- 1000

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

american_to_prob <- function(odds) {
  ifelse(odds > 0, 100 / (odds + 100), abs(odds) / (abs(odds) + 100))
}

devig_pair <- function(odds1, odds2) {
  p1 <- american_to_prob(odds1)
  p2 <- american_to_prob(odds2)
  total <- p1 + p2
  list(prob1 = p1 / total, prob2 = p2 / total)
}

calc_logloss <- function(prob, actual) {
  -(actual * log(pmax(prob, 1e-10)) + (1 - actual) * log(pmax(1 - prob, 1e-10)))
}

# =============================================================================
# LOAD DATA
# =============================================================================

cat("Loading data...\n")
con <- dbConnect(duckdb(), dbdir = "cbb.duckdb")

# Load PBP with consensus lines
pbp <- dbGetQuery(con, "SELECT * FROM cbb_betting_pbp")

closing_odds <- dbGetQuery(con, "
  SELECT id as odds_id, home_spread, total_line,
         spread_home_odds, spread_away_odds, tot_over_odds, tot_under_odds
  FROM cbb_closing_odds
  WHERE spread_home_odds IS NOT NULL
")

# Load derivative odds
deriv_odds <- dbGetQuery(con, "SELECT * FROM cbb_derivative_closing_odds")

dbDisconnect(con)

cat("PBP rows:", nrow(pbp), "\n")
cat("Derivative odds rows:", nrow(deriv_odds), "\n")

# =============================================================================
# PREPARE CONSENSUS DATA
# =============================================================================

closing_odds <- closing_odds %>%
  mutate(
    implied_home = american_to_prob(spread_home_odds),
    implied_away = american_to_prob(spread_away_odds),
    implied_over = american_to_prob(tot_over_odds),
    implied_under = american_to_prob(tot_under_odds),
    spread_juice = implied_home + implied_away,
    total_juice = implied_over + implied_under,
    devig_home = implied_home / spread_juice,
    devig_over = implied_over / total_juice
  )

consensus <- closing_odds %>%
  group_by(odds_id) %>%
  summarize(
    home_spread = median(home_spread, na.rm = TRUE),
    total_line = median(total_line, na.rm = TRUE),
    consensus_home_cover_prob = mean(devig_home, na.rm = TRUE),
    consensus_over_prob = mean(devig_over, na.rm = TRUE),
    .groups = "drop"
  )

# Join PBP with consensus
DT <- pbp %>%
  inner_join(consensus, by = "odds_id") %>%
  filter(!is.na(home_spread), !is.na(total_line),
         !is.na(game_home_margin_h1), !is.na(game_home_margin_h2)) %>%
  mutate(
    actual_cover = ifelse(game_home_margin_fg > -home_spread, 1, 0),
    actual_over = ifelse(game_total_fg > total_line, 1, 0)
  ) %>%
  as.data.table()

cat("Games with full data:", nrow(DT), "\n")

# =============================================================================
# PREPARE DERIVATIVE ODDS (PIVOT TO WIDE FORMAT)
# =============================================================================

cat("\nPreparing derivative odds...\n")

# Get best odds for spreads_h1
spreads_h1 <- deriv_odds %>%
  filter(market == "spreads_h1") %>%
  mutate(
    side = ifelse(outcome_name == home_team, "home", "away"),
    line = ifelse(side == "home", outcome_point, -outcome_point)
  ) %>%
  group_by(event_id, line) %>%
  summarize(
    home_team = first(home_team),
    away_team = first(away_team),
    best_home_odds = max(outcome_price[side == "home"], na.rm = TRUE),
    best_away_odds = max(outcome_price[side == "away"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(best_home_odds), is.finite(best_away_odds))

# Get best odds for totals_h1
totals_h1 <- deriv_odds %>%
  filter(market == "totals_h1") %>%
  mutate(side = tolower(outcome_name)) %>%
  group_by(event_id, outcome_point) %>%
  summarize(
    home_team = first(home_team),
    best_over_odds = max(outcome_price[side == "over"], na.rm = TRUE),
    best_under_odds = max(outcome_price[side == "under"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(line = outcome_point) %>%
  filter(is.finite(best_over_odds), is.finite(best_under_odds))

# Get best odds for h2h_h1 (moneylines)
ml_h1 <- deriv_odds %>%
  filter(market == "h2h_h1") %>%
  mutate(side = ifelse(outcome_name == home_team, "home", "away")) %>%
  group_by(event_id) %>%
  summarize(
    home_team = first(home_team),
    best_home_odds = max(outcome_price[side == "home"], na.rm = TRUE),
    best_away_odds = max(outcome_price[side == "away"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(best_home_odds), is.finite(best_away_odds))

cat("H1 spreads lines:", nrow(spreads_h1), "\n")
cat("H1 totals lines:", nrow(totals_h1), "\n")
cat("H1 moneylines:", nrow(ml_h1), "\n")

# =============================================================================
# RUN BACKTEST
# =============================================================================

cat("\n===========================================\n")
cat("RUNNING BACKTEST\n")
cat("===========================================\n\n")

# Get games that have derivative odds
games_with_derivs <- unique(c(spreads_h1$event_id, totals_h1$event_id, ml_h1$event_id))

# Match to our PBP data via odds_id
DT_test <- DT %>% filter(odds_id %in% games_with_derivs)
cat("Games with derivative odds:", nrow(DT_test), "\n")

# Sample if needed
if (!is.null(N_TEST_GAMES) && N_TEST_GAMES < nrow(DT_test)) {
  set.seed(42)
  test_ids <- sample(unique(DT_test$odds_id), N_TEST_GAMES)
  cat("Testing on random sample of", length(test_ids), "games\n\n")
} else {
  test_ids <- unique(DT_test$odds_id)
  cat("Testing on all", length(test_ids), "games\n\n")
}

# Compute dispersion for sampling
disp <- compute_dispersion(DT, moneyline = FALSE, spread_col = "home_spread", total_col = "total_line")
ss <- disp$ss
st <- disp$st
N <- round(nrow(DT) * SAMPLE_PCT, 0)

cat("Sample size N:", N, "\n\n")

# Run predictions
all_results <- list()

for (i in seq_along(test_ids)) {
  if (i %% 50 == 0) cat("Progress:", i, "/", length(test_ids), "\n")

  test_id <- test_ids[i]
  test_game <- DT_test[DT_test$odds_id == test_id, ][1, ]

  # Get parent lines
  parent_spread <- test_game$home_spread
  parent_total <- test_game$total_line
  target_cover <- test_game$consensus_home_cover_prob
  target_over <- test_game$consensus_over_prob

  if (is.na(parent_spread) || is.na(parent_total)) next
  if (is.na(target_cover) || is.na(target_over)) next

  # Leave-one-out training data
  DT_train <- DT[DT$odds_id != test_id, ]

  # Run answer key sampling
  sample_result <- tryCatch({
    run_answer_key_sample(
      id = test_id,
      parent_spread = parent_spread,
      parent_total = parent_total,
      target_cover = target_cover,
      target_over = target_over,
      DT = DT_train,
      ss = ss, st = st, N = N,
      use_spread_line = TRUE
    )
  }, error = function(e) NULL)

  if (is.null(sample_result)) next

  sample <- sample_result$sample

  # =========== H1 SPREADS ===========
  game_spreads <- spreads_h1 %>% filter(event_id == test_id)
  if (nrow(game_spreads) > 0) {
    margins <- sample$game_home_margin_h1
    actual_margin <- test_game$game_home_margin_h1

    for (j in 1:nrow(game_spreads)) {
      line <- game_spreads$line[j]
      home_odds <- game_spreads$best_home_odds[j]
      away_odds <- game_spreads$best_away_odds[j]

      # Algo prediction: P(home covers) = P(margin > -line)
      non_push <- margins != -line
      if (sum(non_push) < 10) next

      algo_home <- mean(margins[non_push] > -line)
      algo_away <- 1 - algo_home

      # Book implied probs
      book <- devig_pair(home_odds, away_odds)

      # Actual outcome
      if (actual_margin == -line) next  # Push
      actual_home <- ifelse(actual_margin > -line, 1, 0)

      all_results[[length(all_results) + 1]] <- tibble(
        game_id = test_id, market = "spreads_h1", line = line,
        side = "home", algo_prob = algo_home, book_prob = book$prob1,
        book_odds = home_odds, actual = actual_home
      )
      all_results[[length(all_results) + 1]] <- tibble(
        game_id = test_id, market = "spreads_h1", line = line,
        side = "away", algo_prob = algo_away, book_prob = book$prob2,
        book_odds = away_odds, actual = 1 - actual_home
      )
    }
  }

  # =========== H1 TOTALS ===========
  game_totals <- totals_h1 %>% filter(event_id == test_id)
  if (nrow(game_totals) > 0) {
    totals <- sample$game_total_h1
    actual_total <- test_game$game_total_h1

    for (j in 1:nrow(game_totals)) {
      line <- game_totals$line[j]
      over_odds <- game_totals$best_over_odds[j]
      under_odds <- game_totals$best_under_odds[j]

      # Algo prediction
      non_push <- totals != line
      if (sum(non_push) < 10) next

      algo_over <- mean(totals[non_push] > line)
      algo_under <- 1 - algo_over

      # Book implied probs
      book <- devig_pair(over_odds, under_odds)

      # Actual outcome
      if (actual_total == line) next
      actual_over <- ifelse(actual_total > line, 1, 0)

      all_results[[length(all_results) + 1]] <- tibble(
        game_id = test_id, market = "totals_h1", line = line,
        side = "over", algo_prob = algo_over, book_prob = book$prob1,
        book_odds = over_odds, actual = actual_over
      )
      all_results[[length(all_results) + 1]] <- tibble(
        game_id = test_id, market = "totals_h1", line = line,
        side = "under", algo_prob = algo_under, book_prob = book$prob2,
        book_odds = under_odds, actual = 1 - actual_over
      )
    }
  }

  # =========== H1 MONEYLINE ===========
  game_ml <- ml_h1 %>% filter(event_id == test_id)
  if (nrow(game_ml) > 0) {
    margins <- sample$game_home_margin_h1
    actual_margin <- test_game$game_home_margin_h1

    home_odds <- game_ml$best_home_odds[1]
    away_odds <- game_ml$best_away_odds[1]

    # Algo prediction (exclude ties)
    non_tie <- margins != 0
    if (sum(non_tie) >= 10) {
      algo_home <- mean(margins[non_tie] > 0)
      algo_away <- 1 - algo_home

      # Book implied probs
      book <- devig_pair(home_odds, away_odds)

      # Actual outcome
      if (actual_margin != 0) {
        actual_home <- ifelse(actual_margin > 0, 1, 0)

        all_results[[length(all_results) + 1]] <- tibble(
          game_id = test_id, market = "h2h_h1", line = NA_real_,
          side = "home", algo_prob = algo_home, book_prob = book$prob1,
          book_odds = home_odds, actual = actual_home
        )
        all_results[[length(all_results) + 1]] <- tibble(
          game_id = test_id, market = "h2h_h1", line = NA_real_,
          side = "away", algo_prob = algo_away, book_prob = book$prob2,
          book_odds = away_odds, actual = 1 - actual_home
        )
      }
    }
  }
}

# Combine results
results <- bind_rows(all_results)
cat("\nTotal predictions:", nrow(results), "\n")

# =============================================================================
# CALCULATE METRICS
# =============================================================================

results <- results %>%
  mutate(
    algo_logloss = calc_logloss(algo_prob, actual),
    book_logloss = calc_logloss(book_prob, actual),
    edge = algo_prob - book_prob
  )

# =============================================================================
# RESULTS
# =============================================================================

cat("\n===========================================\n")
cat("RESULTS BY MARKET\n")
cat("===========================================\n\n")

by_market <- results %>%
  group_by(market) %>%
  summarize(
    n = n(),
    n_games = n_distinct(game_id),
    algo_logloss = mean(algo_logloss),
    book_logloss = mean(book_logloss),
    logloss_diff = book_logloss - algo_logloss,
    algo_wins = algo_logloss < book_logloss,
    avg_edge = mean(edge),
    .groups = "drop"
  )

print(by_market)

cat("\n===========================================\n")
cat("OVERALL COMPARISON\n")
cat("===========================================\n\n")

overall <- results %>%
  summarize(
    total_predictions = n(),
    total_games = n_distinct(game_id),
    algo_logloss = mean(algo_logloss),
    book_logloss = mean(book_logloss),
    logloss_diff = book_logloss - algo_logloss,
    algo_wins = algo_logloss < book_logloss
  )

cat("Total predictions:", overall$total_predictions, "\n")
cat("Total games:", overall$total_games, "\n")
cat("Algo log-loss:", round(overall$algo_logloss, 4), "\n")
cat("Book log-loss:", round(overall$book_logloss, 4), "\n")
cat("Difference (book - algo):", round(overall$logloss_diff, 4), "\n")
cat("Winner:", ifelse(overall$algo_wins, "ALGO", "BOOK"), "\n")

cat("\n===========================================\n")
cat("CALIBRATION BY PROBABILITY BIN\n")
cat("===========================================\n\n")

calibration <- results %>%
  mutate(prob_bin = cut(algo_prob, breaks = seq(0, 1, 0.1), include.lowest = TRUE)) %>%
  group_by(prob_bin) %>%
  summarize(
    n = n(),
    avg_algo = mean(algo_prob),
    avg_book = mean(book_prob),
    avg_actual = mean(actual),
    algo_cal_err = avg_algo - avg_actual,
    book_cal_err = avg_book - avg_actual,
    .groups = "drop"
  )

print(calibration)

cat("\n===========================================\n")
cat("BACKTEST COMPLETE\n")
cat("===========================================\n")
