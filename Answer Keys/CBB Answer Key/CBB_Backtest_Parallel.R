# CBB Derivative Backtest - Parallel Version with ROI
# Uses parallel processing for speed + calculates ROI on +EV bets

setwd("~/NFLWork/Answer Keys")
library(data.table)
library(duckdb)
library(dplyr)
library(tidyverse)
library(DBI)
library(parallel)
source("Tools.R")

# =============================================================================
# CONFIGURATION
# =============================================================================

N_TEST_GAMES <- 1000    # Number of games to test (NULL for all)
SAMPLE_PCT <- 0.05      # 5% sample for answer key
KELLY_MULT <- 0.25      # Fractional Kelly
BANKROLL <- 1000        # Starting bankroll for ROI simulation
N_CORES <- detectCores() - 1  # Leave 1 core free

cat("Using", N_CORES, "cores for parallel processing\n\n")

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

calc_ev <- function(algo_prob, book_odds) {
  decimal_odds <- ifelse(book_odds > 0, 1 + book_odds / 100, 1 + 100 / abs(book_odds))
  algo_prob * decimal_odds - 1
}

kelly_stake <- function(ev, book_odds, bankroll, kelly_mult) {
  if (ev <= 0) return(0)
  decimal_odds <- ifelse(book_odds > 0, 1 + book_odds / 100, 1 + 100 / abs(book_odds))
  b <- decimal_odds - 1
  p <- (ev + 1) / decimal_odds
  q <- 1 - p
  kelly_full <- max(0, (b * p - q) / b)
  bankroll * kelly_mult * kelly_full
}

# =============================================================================
# LOAD DATA (done once in main process)
# =============================================================================

cat("Loading data...\n")
con <- dbConnect(duckdb(), dbdir = "cbb.duckdb")

pbp <- dbGetQuery(con, "SELECT * FROM cbb_betting_pbp")

closing_odds <- dbGetQuery(con, "
  SELECT id as odds_id, home_spread, total_line,
         spread_home_odds, spread_away_odds, tot_over_odds, tot_under_odds
  FROM cbb_closing_odds
  WHERE spread_home_odds IS NOT NULL
")

deriv_odds <- dbGetQuery(con, "SELECT * FROM cbb_derivative_closing_odds")

dbDisconnect(con)

cat("PBP rows:", nrow(pbp), "\n")
cat("Derivative odds rows:", nrow(deriv_odds), "\n")

# =============================================================================
# PREPARE DATA
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

# Prepare derivative odds (pivot to wide format)
cat("\nPreparing derivative odds...\n")

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

# Get test games
games_with_derivs <- unique(c(spreads_h1$event_id, totals_h1$event_id, ml_h1$event_id))
DT_test <- DT %>% filter(odds_id %in% games_with_derivs)
cat("\nGames with derivative odds:", nrow(DT_test), "\n")

# Sample if needed
if (!is.null(N_TEST_GAMES) && N_TEST_GAMES < nrow(DT_test)) {
  set.seed(42)
  test_ids <- sample(unique(DT_test$odds_id), N_TEST_GAMES)
  cat("Testing on random sample of", length(test_ids), "games\n")
} else {
  test_ids <- unique(DT_test$odds_id)
  cat("Testing on all", length(test_ids), "games\n")
}

# Compute dispersion
disp <- compute_dispersion(DT, moneyline = FALSE, spread_col = "home_spread", total_col = "total_line")
ss <- disp$ss
st <- disp$st
N <- round(nrow(DT) * SAMPLE_PCT, 0)

cat("Sample size N:", N, "\n\n")

# =============================================================================
# DEFINE WORKER FUNCTION (processes a single game)
# =============================================================================

process_single_game <- function(test_id) {
  results <- list()

  test_game <- DT_test[DT_test$odds_id == test_id, ][1, ]

  parent_spread <- test_game$home_spread
  parent_total <- test_game$total_line
  target_cover <- test_game$consensus_home_cover_prob
  target_over <- test_game$consensus_over_prob

  if (is.na(parent_spread) || is.na(parent_total)) return(NULL)
  if (is.na(target_cover) || is.na(target_over)) return(NULL)

  DT_train <- DT[DT$odds_id != test_id, ]

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

  if (is.null(sample_result)) return(NULL)

  sample <- sample_result$sample

  # H1 SPREADS
  game_spreads <- spreads_h1[spreads_h1$event_id == test_id, ]
  if (nrow(game_spreads) > 0) {
    margins <- sample$game_home_margin_h1
    actual_margin <- test_game$game_home_margin_h1

    for (j in 1:nrow(game_spreads)) {
      line <- game_spreads$line[j]
      home_odds <- game_spreads$best_home_odds[j]
      away_odds <- game_spreads$best_away_odds[j]

      non_push <- margins != -line
      if (sum(non_push) < 10) next

      algo_home <- mean(margins[non_push] > -line)
      algo_away <- 1 - algo_home

      book <- devig_pair(home_odds, away_odds)

      if (actual_margin == -line) next
      actual_home <- ifelse(actual_margin > -line, 1, 0)

      results[[length(results) + 1]] <- data.frame(
        game_id = test_id, market = "spreads_h1", line = line,
        side = "home", algo_prob = algo_home, book_prob = book$prob1,
        book_odds = home_odds, actual = actual_home, stringsAsFactors = FALSE
      )
      results[[length(results) + 1]] <- data.frame(
        game_id = test_id, market = "spreads_h1", line = line,
        side = "away", algo_prob = algo_away, book_prob = book$prob2,
        book_odds = away_odds, actual = 1 - actual_home, stringsAsFactors = FALSE
      )
    }
  }

  # H1 TOTALS
  game_totals <- totals_h1[totals_h1$event_id == test_id, ]
  if (nrow(game_totals) > 0) {
    totals <- sample$game_total_h1
    actual_total <- test_game$game_total_h1

    for (j in 1:nrow(game_totals)) {
      line <- game_totals$line[j]
      over_odds <- game_totals$best_over_odds[j]
      under_odds <- game_totals$best_under_odds[j]

      non_push <- totals != line
      if (sum(non_push) < 10) next

      algo_over <- mean(totals[non_push] > line)
      algo_under <- 1 - algo_over

      book <- devig_pair(over_odds, under_odds)

      if (actual_total == line) next
      actual_over <- ifelse(actual_total > line, 1, 0)

      results[[length(results) + 1]] <- data.frame(
        game_id = test_id, market = "totals_h1", line = line,
        side = "over", algo_prob = algo_over, book_prob = book$prob1,
        book_odds = over_odds, actual = actual_over, stringsAsFactors = FALSE
      )
      results[[length(results) + 1]] <- data.frame(
        game_id = test_id, market = "totals_h1", line = line,
        side = "under", algo_prob = algo_under, book_prob = book$prob2,
        book_odds = under_odds, actual = 1 - actual_over, stringsAsFactors = FALSE
      )
    }
  }

  # H1 MONEYLINE
  game_ml <- ml_h1[ml_h1$event_id == test_id, ]
  if (nrow(game_ml) > 0) {
    margins <- sample$game_home_margin_h1
    actual_margin <- test_game$game_home_margin_h1

    home_odds <- game_ml$best_home_odds[1]
    away_odds <- game_ml$best_away_odds[1]

    non_tie <- margins != 0
    if (sum(non_tie) >= 10 && actual_margin != 0) {
      algo_home <- mean(margins[non_tie] > 0)
      algo_away <- 1 - algo_home

      book <- devig_pair(home_odds, away_odds)
      actual_home <- ifelse(actual_margin > 0, 1, 0)

      results[[length(results) + 1]] <- data.frame(
        game_id = test_id, market = "h2h_h1", line = NA_real_,
        side = "home", algo_prob = algo_home, book_prob = book$prob1,
        book_odds = home_odds, actual = actual_home, stringsAsFactors = FALSE
      )
      results[[length(results) + 1]] <- data.frame(
        game_id = test_id, market = "h2h_h1", line = NA_real_,
        side = "away", algo_prob = algo_away, book_prob = book$prob2,
        book_odds = away_odds, actual = 1 - actual_home, stringsAsFactors = FALSE
      )
    }
  }

  if (length(results) == 0) return(NULL)
  do.call(rbind, results)
}

# =============================================================================
# RUN PARALLEL BACKTEST
# =============================================================================

cat("===========================================\n")
cat("RUNNING PARALLEL BACKTEST\n")
cat("===========================================\n\n")

start_time <- Sys.time()

# Create cluster
cl <- makeCluster(N_CORES)

# Export data and functions to workers
clusterExport(cl, c("DT", "DT_test", "spreads_h1", "totals_h1", "ml_h1",
                    "ss", "st", "N", "devig_pair", "process_single_game",
                    "american_to_prob"))

# Load required packages and source Tools.R on workers
clusterEvalQ(cl, {
  library(data.table)
  library(dplyr)
  setwd("~/NFLWork/Answer Keys")
  source("Tools.R")
})

cat("Starting parallel processing on", N_CORES, "cores...\n")
cat("Processing", length(test_ids), "games\n\n")

# Run in parallel with error capture
results_list <- parLapplyLB(cl, test_ids, function(id) {
  tryCatch(
    process_single_game(id),
    error = function(e) {
      list(error = TRUE, id = id, msg = conditionMessage(e))
    }
  )
})

stopCluster(cl)

end_time <- Sys.time()
elapsed <- as.numeric(difftime(end_time, start_time, units = "mins"))

cat("Completed in", round(elapsed, 1), "minutes\n\n")

# Check for errors
errors <- sapply(results_list, function(x) {
  if (is.list(x) && !is.null(x$error)) TRUE else FALSE
})
n_errors <- sum(errors)
cat("Games with errors:", n_errors, "\n")
if (n_errors > 0) {
  error_msgs <- results_list[errors]
  cat("Sample errors:\n")
  for (i in 1:min(5, length(error_msgs))) {
    cat("  ID:", error_msgs[[i]]$id, "- Error:", error_msgs[[i]]$msg, "\n")
  }
}

# Filter to successful results only
valid_results <- results_list[!errors & !sapply(results_list, is.null)]
cat("Successful games:", length(valid_results), "\n")

# Combine results
results <- bind_rows(valid_results)
cat("Total predictions:", nrow(results), "\n")

if (nrow(results) == 0) {
  cat("\nNo valid results to analyze. Exiting.\n")
  quit(save = "no", status = 1)
}

# =============================================================================
# CALCULATE METRICS
# =============================================================================

results <- results %>%
  mutate(
    algo_logloss = calc_logloss(algo_prob, actual),
    book_logloss = calc_logloss(book_prob, actual),
    edge = algo_prob - book_prob,
    ev = calc_ev(algo_prob, book_odds)
  )

# =============================================================================
# LOG-LOSS RESULTS
# =============================================================================

cat("\n===========================================\n")
cat("LOG-LOSS RESULTS BY MARKET\n")
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
    .groups = "drop"
  )

print(by_market)

cat("\n===========================================\n")
cat("OVERALL LOG-LOSS COMPARISON\n")
cat("===========================================\n\n")

overall_ll <- results %>%
  summarize(
    total_predictions = n(),
    total_games = n_distinct(game_id),
    algo_logloss = mean(algo_logloss),
    book_logloss = mean(book_logloss),
    logloss_diff = book_logloss - algo_logloss,
    algo_wins = algo_logloss < book_logloss
  )

cat("Total predictions:", overall_ll$total_predictions, "\n")
cat("Total games:", overall_ll$total_games, "\n")
cat("Algo log-loss:", round(overall_ll$algo_logloss, 4), "\n")
cat("Book log-loss:", round(overall_ll$book_logloss, 4), "\n")
cat("Difference (book - algo):", round(overall_ll$logloss_diff, 4), "\n")
cat("Log-loss Winner:", ifelse(overall_ll$algo_wins, "ALGO", "BOOK"), "\n")

# =============================================================================
# ROI ANALYSIS (on +EV bets)
# =============================================================================

cat("\n===========================================\n")
cat("ROI ANALYSIS (+EV BETS)\n")
cat("===========================================\n\n")

# Filter to +EV bets only
ev_bets <- results %>%
  filter(ev > 0) %>%
  mutate(
    stake = mapply(kelly_stake, ev, book_odds,
                   MoreArgs = list(bankroll = BANKROLL, kelly_mult = KELLY_MULT)),
    decimal_odds = ifelse(book_odds > 0, 1 + book_odds / 100, 1 + 100 / abs(book_odds)),
    pnl = ifelse(actual == 1, stake * (decimal_odds - 1), -stake)
  )

cat("Total +EV bets:", nrow(ev_bets), "\n")
cat("Percentage of all bets:", round(nrow(ev_bets) / nrow(results) * 100, 1), "%\n\n")

if (nrow(ev_bets) > 0) {
  # ROI by market
  roi_by_market <- ev_bets %>%
    group_by(market) %>%
    summarize(
      n_bets = n(),
      n_games = n_distinct(game_id),
      avg_ev = mean(ev) * 100,
      total_staked = sum(stake),
      total_pnl = sum(pnl),
      roi = total_pnl / total_staked * 100,
      win_rate = mean(actual) * 100,
      avg_algo_prob = mean(algo_prob) * 100,
      .groups = "drop"
    )

  cat("ROI by Market:\n")
  print(roi_by_market)

  # Overall ROI
  cat("\n--- Overall ROI ---\n")
  cat("Total bets:", nrow(ev_bets), "\n")
  cat("Total staked: $", round(sum(ev_bets$stake), 0), "\n")
  cat("Total P&L: $", round(sum(ev_bets$pnl), 0), "\n")
  cat("ROI:", round(sum(ev_bets$pnl) / sum(ev_bets$stake) * 100, 2), "%\n")
  cat("Win rate:", round(mean(ev_bets$actual) * 100, 1), "%\n")
  cat("Avg predicted win rate:", round(mean(ev_bets$algo_prob) * 100, 1), "%\n")

  # ROI by EV bucket
  cat("\n--- ROI by EV Bucket ---\n")
  ev_buckets <- ev_bets %>%
    mutate(ev_bucket = cut(ev * 100, breaks = c(0, 2, 5, 10, 20, Inf),
                           labels = c("0-2%", "2-5%", "5-10%", "10-20%", "20%+"))) %>%
    group_by(ev_bucket) %>%
    summarize(
      n_bets = n(),
      avg_ev = mean(ev) * 100,
      total_staked = sum(stake),
      total_pnl = sum(pnl),
      roi = total_pnl / total_staked * 100,
      win_rate = mean(actual) * 100,
      .groups = "drop"
    )
  print(ev_buckets)
}

# =============================================================================
# CALIBRATION
# =============================================================================

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
