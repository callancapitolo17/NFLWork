# CBB Derivative Backtest - BetOnline Parameter Sweep
# Tests different sample sizes and EV thresholds using ONLY BetOnline odds

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

BOOKMAKER <- "betonlineag"
N_TEST_GAMES <- NULL  # All games
SAMPLE_PCTS <- c(0.02, 0.05, 0.075, 0.10)  # 2%, 5%, 7.5%, 10%
EV_THRESHOLDS <- c(0, 0.02, 0.05, 0.10)    # 0%, 2%, 5%, 10%
KELLY_MULT <- 0.25
BANKROLL <- 1000
N_CORES <- detectCores() - 1
RANDOM_SEED <- 42

cat("===========================================\n")
cat("CBB PARAMETER SWEEP - BETONLINE ONLY\n")
cat("===========================================\n\n")
cat("Bookmaker:", BOOKMAKER, "\n")
cat("Test games:", ifelse(is.null(N_TEST_GAMES), "All", N_TEST_GAMES), "\n")
cat("Sample sizes:", paste0(SAMPLE_PCTS * 100, "%", collapse = ", "), "\n")
cat("EV thresholds:", paste0(EV_THRESHOLDS * 100, "%", collapse = ", "), "\n")
cat("Cores:", N_CORES, "\n\n")

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
# LOAD DATA (once for all runs)
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

# FILTER TO BETONLINE ONLY
deriv_odds <- dbGetQuery(con, sprintf("
  SELECT * FROM cbb_derivative_closing_odds
  WHERE bookmaker_key = '%s'
", BOOKMAKER))

dbDisconnect(con)

cat("PBP rows:", nrow(pbp), "\n")
cat("BetOnline derivative odds rows:", nrow(deriv_odds), "\n")

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

# Prepare BetOnline derivative odds
cat("\nPreparing BetOnline derivative odds...\n")

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

# Get test games (same for all runs)
games_with_derivs <- unique(c(spreads_h1$event_id, totals_h1$event_id, ml_h1$event_id))
DT_test <- DT %>% filter(odds_id %in% games_with_derivs)
cat("\nGames with BetOnline derivative odds:", nrow(DT_test), "\n")

set.seed(RANDOM_SEED)
if (!is.null(N_TEST_GAMES) && N_TEST_GAMES < nrow(DT_test)) {
  test_ids <- sample(unique(DT_test$odds_id), N_TEST_GAMES)
} else {
  test_ids <- unique(DT_test$odds_id)
}
cat("Testing on", length(test_ids), "games (same across all sample sizes)\n")

# Compute dispersion (same for all)
disp <- compute_dispersion(DT, moneyline = FALSE, spread_col = "home_spread", total_col = "total_line")
ss <- disp$ss
st <- disp$st

# =============================================================================
# WORKER FUNCTION
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
      ss = ss, st = st, N = N_SAMPLE,
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
# RUN BACKTEST FOR EACH SAMPLE SIZE
# =============================================================================

all_results <- list()

for (sample_pct in SAMPLE_PCTS) {
  cat("\n===========================================\n")
  cat("SAMPLE SIZE:", sample_pct * 100, "%\n")
  cat("===========================================\n")

  N_SAMPLE <- round(nrow(DT) * sample_pct, 0)
  cat("N_SAMPLE:", N_SAMPLE, "\n")

  start_time <- Sys.time()

  # Create cluster
  cl <- makeCluster(N_CORES)

  # Export data and functions to workers
  clusterExport(cl, c("DT", "DT_test", "spreads_h1", "totals_h1", "ml_h1",
                      "ss", "st", "N_SAMPLE", "devig_pair", "process_single_game",
                      "american_to_prob"))

  # Load packages and source Tools.R on workers
  clusterEvalQ(cl, {
    library(data.table)
    library(dplyr)
    setwd("~/NFLWork/Answer Keys")
    source("Tools.R")
  })

  cat("Starting parallel processing...\n")

  results_list <- parLapplyLB(cl, test_ids, function(id) {
    tryCatch(
      process_single_game(id),
      error = function(e) list(error = TRUE, id = id, msg = conditionMessage(e))
    )
  })

  stopCluster(cl)

  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "mins"))
  cat("Completed in", round(elapsed, 1), "minutes\n")

  # Check for errors
  errors <- sapply(results_list, function(x) {
    if (is.list(x) && !is.null(x$error)) TRUE else FALSE
  })
  n_errors <- sum(errors)
  cat("Errors:", n_errors, "\n")

  # Combine valid results
  valid_results <- results_list[!errors & !sapply(results_list, is.null)]
  results <- bind_rows(valid_results)

  if (nrow(results) == 0) {
    cat("No valid results for this sample size. Skipping.\n")
    next
  }

  # Add EV column
  results <- results %>%
    mutate(
      ev = calc_ev(algo_prob, book_odds),
      sample_pct = sample_pct
    )

  cat("Predictions:", nrow(results), "\n")

  all_results[[as.character(sample_pct)]] <- results
}

# Combine all results
full_results <- bind_rows(all_results)
cat("\n\nTotal predictions across all sample sizes:", nrow(full_results), "\n")

# =============================================================================
# APPLY EV THRESHOLDS AND GENERATE SUMMARY
# =============================================================================

cat("\n===========================================\n")
cat("BETONLINE PARAMETER SWEEP RESULTS\n")
cat("===========================================\n\n")

summary_rows <- list()

for (sample_pct in SAMPLE_PCTS) {
  results <- full_results %>% filter(sample_pct == !!sample_pct)

  for (ev_thresh in EV_THRESHOLDS) {
    # Filter to bets above threshold
    ev_bets <- results %>%
      filter(ev > ev_thresh) %>%
      mutate(
        stake = mapply(kelly_stake, ev, book_odds,
                       MoreArgs = list(bankroll = BANKROLL, kelly_mult = KELLY_MULT)),
        decimal_odds = ifelse(book_odds > 0, 1 + book_odds / 100, 1 + 100 / abs(book_odds)),
        pnl = ifelse(actual == 1, stake * (decimal_odds - 1), -stake)
      )

    if (nrow(ev_bets) == 0) {
      summary_rows[[length(summary_rows) + 1]] <- tibble(
        sample_pct = sample_pct * 100,
        ev_thresh = ev_thresh * 100,
        n_bets = 0, win_rate = NA, roi = NA,
        roi_spreads = NA, roi_totals = NA, roi_ml = NA
      )
      next
    }

    # Overall stats
    n_bets <- nrow(ev_bets)
    win_rate <- mean(ev_bets$actual) * 100
    total_staked <- sum(ev_bets$stake)
    total_pnl <- sum(ev_bets$pnl)
    roi <- if (total_staked > 0) total_pnl / total_staked * 100 else NA

    # By market
    calc_market_roi <- function(df, mkt) {
      mkt_bets <- df %>% filter(market == mkt)
      if (nrow(mkt_bets) == 0 || sum(mkt_bets$stake) == 0) return(NA)
      sum(mkt_bets$pnl) / sum(mkt_bets$stake) * 100
    }

    roi_spreads <- calc_market_roi(ev_bets, "spreads_h1")
    roi_totals <- calc_market_roi(ev_bets, "totals_h1")
    roi_ml <- calc_market_roi(ev_bets, "h2h_h1")

    summary_rows[[length(summary_rows) + 1]] <- tibble(
      sample_pct = sample_pct * 100,
      ev_thresh = ev_thresh * 100,
      n_bets = n_bets,
      win_rate = round(win_rate, 1),
      roi = round(roi, 2),
      roi_spreads = round(roi_spreads, 2),
      roi_totals = round(roi_totals, 2),
      roi_ml = round(roi_ml, 2)
    )
  }
}

summary_table <- bind_rows(summary_rows)

cat("SUMMARY TABLE\n")
cat("-------------\n")
print(summary_table, n = 20)

# =============================================================================
# BEST PARAMETERS
# =============================================================================

cat("\n===========================================\n")
cat("BEST PARAMETERS FOR BETONLINE\n")
cat("===========================================\n\n")

best_overall <- summary_table %>%
  filter(!is.na(roi)) %>%
  arrange(desc(roi)) %>%
  head(3)

cat("Top 3 by Overall ROI:\n")
print(best_overall)

best_totals <- summary_table %>%
  filter(!is.na(roi_totals)) %>%
  arrange(desc(roi_totals)) %>%
  head(3)

cat("\nTop 3 by Totals ROI:\n")
print(best_totals)

# Compare to all-books results
cat("\n===========================================\n")
cat("COMPARISON: BETONLINE vs ALL BOOKS\n")
cat("===========================================\n")
cat("(Load parameter_sweep_full_summary.csv to compare)\n")

# =============================================================================
# SAVE RESULTS
# =============================================================================

saveRDS(full_results, "CBB Answer Key/betonline_sweep_predictions.rds")
write.csv(summary_table, "CBB Answer Key/betonline_sweep_summary.csv", row.names = FALSE)

cat("\n===========================================\n")
cat("BETONLINE SWEEP COMPLETE\n")
cat("===========================================\n")
cat("Full predictions saved to: betonline_sweep_predictions.rds\n")
cat("Summary table saved to: betonline_sweep_summary.csv\n")
