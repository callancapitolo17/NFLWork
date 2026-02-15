# CBB Backtest - BetOnline Only
# Tests algorithm performance using only BetOnline odds
# Uses optimal parameters: 2% sample, 5% EV threshold

setwd("~/NFLWork/Answer Keys")
library(data.table)
library(duckdb)
library(dplyr)
library(tidyverse)
library(DBI)
library(parallel)
source("Tools.R")

# =============================================================================
# CONFIGURATION (Optimal parameters from sweep)
# =============================================================================

SAMPLE_PCT <- 0.02      # 2% sample
EV_THRESHOLD <- 0.05    # 5% EV threshold
KELLY_MULT <- 0.25
BANKROLL <- 1000
N_CORES <- detectCores() - 1
BOOKMAKER <- "betonlineag"

cat("===========================================\n")
cat("CBB BACKTEST - BETONLINE ONLY\n")
cat("===========================================\n\n")
cat("Sample size:", SAMPLE_PCT * 100, "%\n")
cat("EV threshold:", EV_THRESHOLD * 100, "%\n")
cat("Bookmaker:", BOOKMAKER, "\n")
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
# LOAD DATA
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

# Load ONLY BetOnline derivative odds
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

# Prepare BetOnline derivative odds (no need to find "best" - only one book)
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
    home_odds = max(outcome_price[side == "home"], na.rm = TRUE),
    away_odds = max(outcome_price[side == "away"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(home_odds), is.finite(away_odds))

totals_h1 <- deriv_odds %>%
  filter(market == "totals_h1") %>%
  mutate(side = tolower(outcome_name)) %>%
  group_by(event_id, outcome_point) %>%
  summarize(
    home_team = first(home_team),
    over_odds = max(outcome_price[side == "over"], na.rm = TRUE),
    under_odds = max(outcome_price[side == "under"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(line = outcome_point) %>%
  filter(is.finite(over_odds), is.finite(under_odds))

ml_h1 <- deriv_odds %>%
  filter(market == "h2h_h1") %>%
  mutate(side = ifelse(outcome_name == home_team, "home", "away")) %>%
  group_by(event_id) %>%
  summarize(
    home_team = first(home_team),
    home_odds = max(outcome_price[side == "home"], na.rm = TRUE),
    away_odds = max(outcome_price[side == "away"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(home_odds), is.finite(away_odds))

cat("H1 spreads lines:", nrow(spreads_h1), "\n")
cat("H1 totals lines:", nrow(totals_h1), "\n")
cat("H1 moneylines:", nrow(ml_h1), "\n")

# Get test games
games_with_derivs <- unique(c(spreads_h1$event_id, totals_h1$event_id, ml_h1$event_id))
DT_test <- DT %>% filter(odds_id %in% games_with_derivs)
cat("\nGames with BetOnline derivative odds:", nrow(DT_test), "\n")

test_ids <- unique(DT_test$odds_id)
cat("Testing on all", length(test_ids), "games\n")

# Compute dispersion
disp <- compute_dispersion(DT, moneyline = FALSE, spread_col = "home_spread", total_col = "total_line")
ss <- disp$ss
st <- disp$st
N_SAMPLE <- round(nrow(DT) * SAMPLE_PCT, 0)
cat("Sample size N:", N_SAMPLE, "\n\n")

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
      home_odds <- game_spreads$home_odds[j]
      away_odds <- game_spreads$away_odds[j]

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
      over_odds <- game_totals$over_odds[j]
      under_odds <- game_totals$under_odds[j]

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

    home_odds <- game_ml$home_odds[1]
    away_odds <- game_ml$away_odds[1]

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
cat("RUNNING BETONLINE BACKTEST\n")
cat("===========================================\n\n")

start_time <- Sys.time()

cl <- makeCluster(N_CORES)

clusterExport(cl, c("DT", "DT_test", "spreads_h1", "totals_h1", "ml_h1",
                    "ss", "st", "N_SAMPLE", "devig_pair", "process_single_game",
                    "american_to_prob"))

clusterEvalQ(cl, {
  library(data.table)
  library(dplyr)
  setwd("~/NFLWork/Answer Keys")
  source("Tools.R")
})

cat("Starting parallel processing on", N_CORES, "cores...\n")
cat("Processing", length(test_ids), "games\n\n")

results_list <- parLapplyLB(cl, test_ids, function(id) {
  tryCatch(
    process_single_game(id),
    error = function(e) list(error = TRUE, id = id, msg = conditionMessage(e))
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
cat("Errors:", n_errors, "\n")

# Combine valid results
valid_results <- results_list[!errors & !sapply(results_list, is.null)]
results <- bind_rows(valid_results)
cat("Total predictions:", nrow(results), "\n")

if (nrow(results) == 0) {
  cat("\nNo valid results. Exiting.\n")
  quit(save = "no", status = 1)
}

# Add EV
results <- results %>%
  mutate(ev = calc_ev(algo_prob, book_odds))

# =============================================================================
# FILTER TO EV THRESHOLD AND CALCULATE ROI
# =============================================================================

cat("\n===========================================\n")
cat("RESULTS (EV >=", EV_THRESHOLD * 100, "%)\n")
cat("===========================================\n\n")

ev_bets <- results %>%
  filter(ev >= EV_THRESHOLD) %>%
  mutate(
    stake = mapply(kelly_stake, ev, book_odds,
                   MoreArgs = list(bankroll = BANKROLL, kelly_mult = KELLY_MULT)),
    decimal_odds = ifelse(book_odds > 0, 1 + book_odds / 100, 1 + 100 / abs(book_odds)),
    pnl = ifelse(actual == 1, stake * (decimal_odds - 1), -stake)
  )

cat("Total +EV bets (>=5%):", nrow(ev_bets), "\n")
cat("Percentage of all predictions:", round(nrow(ev_bets) / nrow(results) * 100, 1), "%\n\n")

if (nrow(ev_bets) > 0) {
  # By market
  by_market <- ev_bets %>%
    group_by(market) %>%
    summarize(
      n_bets = n(),
      n_games = n_distinct(game_id),
      win_rate = mean(actual) * 100,
      avg_ev = mean(ev) * 100,
      total_staked = sum(stake),
      total_pnl = sum(pnl),
      roi = total_pnl / total_staked * 100,
      .groups = "drop"
    )

  cat("BY MARKET:\n")
  print(by_market)

  # Overall
  cat("\n--- OVERALL ---\n")
  cat("Total bets:", nrow(ev_bets), "\n")
  cat("Total staked: $", round(sum(ev_bets$stake), 0), "\n")
  cat("Total P&L: $", round(sum(ev_bets$pnl), 0), "\n")
  cat("ROI:", round(sum(ev_bets$pnl) / sum(ev_bets$stake) * 100, 2), "%\n")
  cat("Win rate:", round(mean(ev_bets$actual) * 100, 1), "%\n")

  # Calibration
  cat("\n--- CALIBRATION ---\n")
  cat("Predicted win rate:", round(mean(ev_bets$algo_prob) * 100, 1), "%\n")
  cat("Actual win rate:", round(mean(ev_bets$actual) * 100, 1), "%\n")
  cat("Calibration error:", round((mean(ev_bets$algo_prob) - mean(ev_bets$actual)) * 100, 1), "pp\n")
}

cat("\n===========================================\n")
cat("BETONLINE BACKTEST COMPLETE\n")
cat("===========================================\n")
