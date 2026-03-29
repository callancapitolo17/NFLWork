# MLB F5 ROI Backtest
# Leave-one-out backtest comparing algo F5 predictions to actual book F5 closing odds
# Parallel version based on CBB_Backtest_Parallel.R
#
# Requires: pbp.duckdb/mlb_f5_odds_history (from fetch_mlb_f5_odds.py)
#
# Usage:
#   Rscript "MLB Answer Key/MLB_ROI_Backtest.R"          # default 500 games
#   Rscript "MLB Answer Key/MLB_ROI_Backtest.R" 2000     # test 2000 games

setwd("~/NFLWork/Answer Keys")
suppressPackageStartupMessages({
  library(data.table)
  library(duckdb)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(DBI)
  library(parallel)
})
source("Tools.R")

# =============================================================================
# CONFIGURATION
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
N_TEST_GAMES <- if (length(args) > 0) as.integer(args[1]) else 500
SAMPLE_PCT   <- 0.02       # 2% sample (matches production MLB.R)
KELLY_MULT   <- 0.25       # Fractional Kelly
BANKROLL     <- 1000       # Starting bankroll for ROI sim
EV_THRESHOLD <- 0.05       # 5% min EV
N_CORES      <- detectCores() - 1

cat("=== MLB F5 ROI BACKTEST ===\n")
cat("Cores:", N_CORES, "| Test games:", N_TEST_GAMES, "\n\n")

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
# LOAD DATA
# =============================================================================

cat("Loading data...\n")
con <- dbConnect(duckdb(), dbdir = "pbp.duckdb", read_only = TRUE)
on.exit(dbDisconnect(con))

# Historical PBP (same rename as MLB.R Phase 1)
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
    game_id = game_pk
  ) %>%
  as.data.table()

cat("  PBP games:", nrow(DT), "\n")

# F5 historical book odds (from fetch_mlb_f5_odds.py)
f5_odds_count <- dbGetQuery(con, "
  SELECT COUNT(DISTINCT game_pk) as n
  FROM mlb_f5_odds_history
")$n

if (f5_odds_count == 0) {
  cat("\nERROR: No F5 historical odds found in mlb_f5_odds_history.\n")
  cat("Run fetch_mlb_f5_odds.py first to acquire F5 closing odds.\n")
  quit(save = "no", status = 1)
}

f5_raw <- dbGetQuery(con, "SELECT * FROM mlb_f5_odds_history")
dbDisconnect(con)
on.exit()  # clear the on.exit since we just disconnected

cat("  F5 odds records:", nrow(f5_raw), "covering", f5_odds_count, "games\n")

# =============================================================================
# PREPARE F5 CLOSING ODDS (pivot to wide best-book format)
# =============================================================================

cat("\nPreparing F5 odds...\n")

# --- F5 Spreads: best home/away spread odds per game/line ---
spreads_f5 <- f5_raw %>%
  filter(market == "spreads_1st_5_innings") %>%
  # Join to get home_team for this game_pk
  inner_join(
    DT %>% select(game_pk, home_team, away_team) %>% distinct(),
    by = "game_pk"
  ) %>%
  mutate(
    side = ifelse(key == home_team, "home",
           ifelse(key == away_team, "away", NA_character_)),
    # Normalize line to home perspective (home point stays, away point flips)
    home_line = ifelse(side == "home", point, -point)
  ) %>%
  filter(!is.na(side)) %>%
  group_by(game_pk, home_line) %>%
  summarize(
    best_home_odds = max(price[side == "home"], na.rm = TRUE),
    best_away_odds = max(price[side == "away"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(line = home_line) %>%
  filter(is.finite(best_home_odds), is.finite(best_away_odds))

# --- F5 Totals: best over/under odds per game/line ---
totals_f5 <- f5_raw %>%
  filter(market == "totals_1st_5_innings") %>%
  mutate(side = tolower(key)) %>%
  group_by(game_pk, point) %>%
  summarize(
    best_over_odds = max(price[side == "over"], na.rm = TRUE),
    best_under_odds = max(price[side == "under"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(line = point) %>%
  filter(is.finite(best_over_odds), is.finite(best_under_odds))

# --- F5 Moneylines: best home/away ML per game ---
ml_f5 <- f5_raw %>%
  filter(market == "h2h_1st_5_innings") %>%
  inner_join(
    DT %>% select(game_pk, home_team, away_team) %>% distinct(),
    by = "game_pk"
  ) %>%
  mutate(
    side = ifelse(key == home_team, "home",
           ifelse(key == away_team, "away", NA_character_))
  ) %>%
  filter(!is.na(side)) %>%
  group_by(game_pk) %>%
  summarize(
    best_home_odds = max(price[side == "home"], na.rm = TRUE),
    best_away_odds = max(price[side == "away"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(best_home_odds), is.finite(best_away_odds))

cat("  F5 spread lines:", nrow(spreads_f5), "\n")
cat("  F5 total lines:", nrow(totals_f5), "\n")
cat("  F5 moneylines:", nrow(ml_f5), "\n")

# =============================================================================
# SELECT TEST GAMES
# =============================================================================

# Games that have F5 derivative odds
games_with_derivs <- unique(c(spreads_f5$game_pk, totals_f5$game_pk, ml_f5$game_pk))
DT_test <- DT[DT$game_pk %in% games_with_derivs, ]
cat("\nGames with F5 odds:", nrow(DT_test), "\n")

if (nrow(DT_test) == 0) {
  cat("ERROR: No games matched between PBP and F5 odds.\n")
  quit(save = "no", status = 1)
}

# Sample test games
if (!is.null(N_TEST_GAMES) && N_TEST_GAMES < nrow(DT_test)) {
  set.seed(42)
  test_ids <- sample(unique(DT_test$game_pk), N_TEST_GAMES)
  cat("Testing random sample of", length(test_ids), "games\n")
} else {
  test_ids <- unique(DT_test$game_pk)
  cat("Testing all", length(test_ids), "games\n")
}

# Compute dispersion (MLB uses moneyline-based matching)
disp <- compute_dispersion(DT, moneyline = TRUE)
ss <- disp$ss
st <- disp$st
N <- round(nrow(DT) * SAMPLE_PCT, 0)

cat("Sample size N:", N, "\n\n")

# =============================================================================
# WORKER FUNCTION (processes a single game)
# =============================================================================

process_single_game <- function(test_id) {
  results <- list()

  test_game <- DT_test[DT_test$game_pk == test_id, ][1, ]

  # Parent lines — MLB uses moneyline probability, not spread
  parent_spread <- test_game$home_ml_odds   # this IS the probability (use_spread_line = FALSE)
  parent_total  <- test_game$total_line
  target_cover  <- test_game$home_ml_odds   # same as parent_spread for MLB
  target_over   <- test_game$over_odds

  if (is.na(parent_spread) || is.na(parent_total)) return(NULL)
  if (is.na(target_cover) || is.na(target_over)) return(NULL)

  # Leave-one-out
  DT_train <- DT[DT$game_pk != test_id, ]

  sample_result <- tryCatch({
    run_answer_key_sample(
      id = test_id,
      parent_spread = parent_spread,
      parent_total = parent_total,
      target_cover = target_cover,
      target_over = target_over,
      DT = DT_train,
      ss = ss, st = st, N = N,
      use_spread_line = FALSE  # CRITICAL for MLB: match on ML probability
    )
  }, error = function(e) NULL)

  if (is.null(sample_result)) return(NULL)

  sample <- sample_result$sample

  # =========== F5 SPREADS ===========
  game_spreads <- spreads_f5[spreads_f5$game_pk == test_id, ]
  if (nrow(game_spreads) > 0) {
    margins <- sample$game_home_margin_period_F5
    actual_margin <- test_game$game_home_margin_period_F5

    for (j in 1:nrow(game_spreads)) {
      line <- game_spreads$line[j]
      home_odds <- game_spreads$best_home_odds[j]
      away_odds <- game_spreads$best_away_odds[j]

      # P(home covers) = P(margin > -line), excluding pushes
      non_push <- margins != -line
      if (sum(non_push) < 10) next

      algo_home <- mean(margins[non_push] > -line)
      algo_away <- 1 - algo_home

      book <- devig_pair(home_odds, away_odds)

      if (actual_margin == -line) next  # Push
      actual_home <- ifelse(actual_margin > -line, 1, 0)

      results[[length(results) + 1]] <- data.frame(
        game_id = test_id, market = "spreads_f5", line = line,
        side = "home", algo_prob = algo_home, book_prob = book$prob1,
        book_odds = home_odds, actual = actual_home, stringsAsFactors = FALSE
      )
      results[[length(results) + 1]] <- data.frame(
        game_id = test_id, market = "spreads_f5", line = line,
        side = "away", algo_prob = algo_away, book_prob = book$prob2,
        book_odds = away_odds, actual = 1 - actual_home, stringsAsFactors = FALSE
      )
    }
  }

  # =========== F5 TOTALS ===========
  game_totals <- totals_f5[totals_f5$game_pk == test_id, ]
  if (nrow(game_totals) > 0) {
    totals_vec <- sample$game_total_period_F5
    actual_total <- test_game$game_total_period_F5

    for (j in 1:nrow(game_totals)) {
      line <- game_totals$line[j]
      over_odds <- game_totals$best_over_odds[j]
      under_odds <- game_totals$best_under_odds[j]

      non_push <- totals_vec != line
      if (sum(non_push) < 10) next

      algo_over <- mean(totals_vec[non_push] > line)
      algo_under <- 1 - algo_over

      book <- devig_pair(over_odds, under_odds)

      if (actual_total == line) next  # Push
      actual_over <- ifelse(actual_total > line, 1, 0)

      results[[length(results) + 1]] <- data.frame(
        game_id = test_id, market = "totals_f5", line = line,
        side = "over", algo_prob = algo_over, book_prob = book$prob1,
        book_odds = over_odds, actual = actual_over, stringsAsFactors = FALSE
      )
      results[[length(results) + 1]] <- data.frame(
        game_id = test_id, market = "totals_f5", line = line,
        side = "under", algo_prob = algo_under, book_prob = book$prob2,
        book_odds = under_odds, actual = 1 - actual_over, stringsAsFactors = FALSE
      )
    }
  }

  # =========== F5 MONEYLINE ===========
  game_ml <- ml_f5[ml_f5$game_pk == test_id, ]
  if (nrow(game_ml) > 0) {
    margins <- sample$game_home_margin_period_F5
    actual_margin <- test_game$game_home_margin_period_F5

    home_odds <- game_ml$best_home_odds[1]
    away_odds <- game_ml$best_away_odds[1]

    # F5 CAN tie (unlike FG with extra innings), so exclude ties
    non_tie <- margins != 0
    if (sum(non_tie) >= 10 && actual_margin != 0) {
      algo_home <- mean(margins[non_tie] > 0)
      algo_away <- 1 - algo_home

      book <- devig_pair(home_odds, away_odds)
      actual_home <- ifelse(actual_margin > 0, 1, 0)

      results[[length(results) + 1]] <- data.frame(
        game_id = test_id, market = "h2h_f5", line = NA_real_,
        side = "home", algo_prob = algo_home, book_prob = book$prob1,
        book_odds = home_odds, actual = actual_home, stringsAsFactors = FALSE
      )
      results[[length(results) + 1]] <- data.frame(
        game_id = test_id, market = "h2h_f5", line = NA_real_,
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

cl <- makeCluster(N_CORES)

clusterExport(cl, c("DT", "DT_test", "spreads_f5", "totals_f5", "ml_f5",
                    "ss", "st", "N", "devig_pair", "process_single_game",
                    "american_to_prob"))

clusterEvalQ(cl, {
  library(data.table)
  library(dplyr)
  setwd("~/NFLWork/Answer Keys")
  source("Tools.R")
})

cat("Processing", length(test_ids), "games on", N_CORES, "cores...\n\n")

results_list <- parLapplyLB(cl, test_ids, function(id) {
  tryCatch(
    process_single_game(id),
    error = function(e) {
      list(error = TRUE, id = id, msg = conditionMessage(e))
    }
  )
})

stopCluster(cl)

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat("Completed in", round(elapsed, 1), "minutes\n\n")

# Check errors
errors <- sapply(results_list, function(x) {
  is.list(x) && !is.null(x$error)
})
n_errors <- sum(errors)
cat("Games with errors:", n_errors, "\n")
if (n_errors > 0) {
  error_msgs <- results_list[errors]
  for (i in 1:min(5, length(error_msgs))) {
    cat("  ID:", error_msgs[[i]]$id, "- Error:", error_msgs[[i]]$msg, "\n")
  }
}

valid_results <- results_list[!errors & !sapply(results_list, is.null)]
cat("Successful games:", length(valid_results), "\n")

results <- bind_rows(valid_results)
cat("Total predictions:", nrow(results), "\n")

if (nrow(results) == 0) {
  cat("\nNo valid results. Check that F5 odds and PBP data overlap.\n")
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
cat("LOG-LOSS: ALGO vs BOOK (by market)\n")
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

cat("\n--- Overall ---\n")
overall <- results %>%
  summarize(
    n = n(), n_games = n_distinct(game_id),
    algo_ll = mean(algo_logloss), book_ll = mean(book_logloss)
  )
cat("Predictions:", overall$n, "| Games:", overall$n_games, "\n")
cat("Algo log-loss:", round(overall$algo_ll, 4), "\n")
cat("Book log-loss:", round(overall$book_ll, 4), "\n")
cat("Winner:", ifelse(overall$algo_ll < overall$book_ll, "ALGO", "BOOK"),
    "(diff:", round(overall$book_ll - overall$algo_ll, 4), ")\n")

# =============================================================================
# ROI ANALYSIS (+EV bets)
# =============================================================================

cat("\n===========================================\n")
cat("ROI ANALYSIS (+EV BETS, threshold =", EV_THRESHOLD * 100, "%)\n")
cat("===========================================\n\n")

ev_bets <- results %>%
  filter(ev > EV_THRESHOLD) %>%
  mutate(
    stake = mapply(kelly_stake, ev, book_odds,
                   MoreArgs = list(bankroll = BANKROLL, kelly_mult = KELLY_MULT)),
    decimal_odds = ifelse(book_odds > 0, 1 + book_odds / 100, 1 + 100 / abs(book_odds)),
    pnl = ifelse(actual == 1, stake * (decimal_odds - 1), -stake)
  )

cat("+EV bets:", nrow(ev_bets), "of", nrow(results), "total",
    "(", round(nrow(ev_bets) / nrow(results) * 100, 1), "%)\n\n")

if (nrow(ev_bets) > 0) {
  roi_by_market <- ev_bets %>%
    group_by(market) %>%
    summarize(
      n_bets = n(),
      n_games = n_distinct(game_id),
      avg_ev = round(mean(ev) * 100, 1),
      total_staked = round(sum(stake), 0),
      total_pnl = round(sum(pnl), 0),
      roi = round(total_pnl / total_staked * 100, 2),
      win_rate = round(mean(actual) * 100, 1),
      .groups = "drop"
    )

  cat("--- By Market ---\n")
  print(roi_by_market)

  cat("\n--- Overall ---\n")
  cat("Total bets:", nrow(ev_bets), "\n")
  cat("Total staked: $", round(sum(ev_bets$stake), 0), "\n")
  cat("Total P&L: $", round(sum(ev_bets$pnl), 0), "\n")
  cat("ROI:", round(sum(ev_bets$pnl) / sum(ev_bets$stake) * 100, 2), "%\n")
  cat("Win rate:", round(mean(ev_bets$actual) * 100, 1), "%\n")
  cat("Avg predicted prob:", round(mean(ev_bets$algo_prob) * 100, 1), "%\n")

  # ROI by EV bucket
  cat("\n--- By EV Bucket ---\n")
  ev_buckets <- ev_bets %>%
    mutate(ev_bucket = cut(ev * 100, breaks = c(0, 2, 5, 10, 20, Inf),
                           labels = c("0-2%", "2-5%", "5-10%", "10-20%", "20%+"))) %>%
    group_by(ev_bucket) %>%
    summarize(
      n_bets = n(),
      avg_ev = round(mean(ev) * 100, 1),
      total_staked = round(sum(stake), 0),
      total_pnl = round(sum(pnl), 0),
      roi = round(total_pnl / total_staked * 100, 2),
      win_rate = round(mean(actual) * 100, 1),
      .groups = "drop"
    )
  print(ev_buckets)
}

# =============================================================================
# CALIBRATION
# =============================================================================

cat("\n===========================================\n")
cat("CALIBRATION (10% bins)\n")
cat("===========================================\n\n")

calibration <- results %>%
  mutate(prob_bin = cut(algo_prob, breaks = seq(0, 1, 0.1), include.lowest = TRUE)) %>%
  group_by(prob_bin) %>%
  summarize(
    n = n(),
    avg_algo = round(mean(algo_prob), 3),
    avg_book = round(mean(book_prob), 3),
    avg_actual = round(mean(actual), 3),
    algo_err = round(mean(algo_prob) - mean(actual), 3),
    book_err = round(mean(book_prob) - mean(actual), 3),
    .groups = "drop"
  )

print(calibration)

cat("\n===========================================\n")
cat("MLB F5 ROI BACKTEST COMPLETE\n")
cat("===========================================\n")
