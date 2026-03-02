# CBB Backtest: Correlation-Adjusted Kelly vs Independent Kelly
# Compares two sizing strategies on the same predictions to validate
# whether correlation adjustment improves risk-adjusted returns.

setwd("~/NFLWork/Answer Keys")
library(data.table)
library(duckdb)
library(dplyr)
library(tidyverse)
library(DBI)
library(parallel)
# Source Tools.R from worktree (has bet_to_leg, multivariate_kelly, etc.)
WORKTREE_TOOLS <- "~/NFLWork/.claude/worktrees/correlation-kelly/Answer Keys/Tools.R"
source(WORKTREE_TOOLS)

# =============================================================================
# CONFIGURATION
# =============================================================================

N_TEST_GAMES <- NULL       # NULL = all games with derivative odds
SAMPLE_PCT   <- 0.02       # 2% sample (matches production)
KELLY_MULT   <- 0.25       # Fractional Kelly
BANKROLL     <- 1000       # Static bankroll for sizing
EV_THRESHOLD <- 0.05       # 5% min EV
N_CORES      <- detectCores() - 1

cat("Correlation Kelly Backtest\n")
cat("Using", N_CORES, "cores\n\n")

# =============================================================================
# HELPER FUNCTIONS (local copies matching existing backtest pattern)
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

calc_kelly <- function(ev, book_odds, bankroll, kelly_mult) {
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

# Add column aliases expected by evaluate_leg() (period naming convention)
DT[, game_home_margin_period_Half1 := game_home_margin_h1]
DT[, game_home_margin_period_Half2 := game_home_margin_h2]
DT[, game_total_period_Half1 := game_total_h1]
DT[, game_total_period_Half2 := game_total_h2]
DT[, home_margin := game_home_margin_fg]
DT[, total_final_score := game_total_fg]

cat("Games with full data:", nrow(DT), "\n")

# =============================================================================
# PREPARE DERIVATIVE ODDS (H1 + H2)
# =============================================================================

cat("\nPreparing derivative odds (H1 + H2)...\n")

prep_spread_odds <- function(odds_df, mkt) {
  odds_df %>%
    filter(market == mkt) %>%
    mutate(
      side = ifelse(outcome_name == home_team, "home", "away"),
      line = ifelse(side == "home", outcome_point, -outcome_point)
    ) %>%
    group_by(event_id, line) %>%
    summarize(
      home_team = first(home_team), away_team = first(away_team),
      best_home_odds = max(outcome_price[side == "home"], na.rm = TRUE),
      best_away_odds = max(outcome_price[side == "away"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(is.finite(best_home_odds), is.finite(best_away_odds)) %>%
    mutate(market = mkt)
}

prep_total_odds <- function(odds_df, mkt) {
  odds_df %>%
    filter(market == mkt) %>%
    mutate(side = tolower(outcome_name)) %>%
    group_by(event_id, outcome_point) %>%
    summarize(
      home_team = first(home_team),
      best_over_odds = max(outcome_price[side == "over"], na.rm = TRUE),
      best_under_odds = max(outcome_price[side == "under"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(line = outcome_point) %>%
    filter(is.finite(best_over_odds), is.finite(best_under_odds)) %>%
    mutate(market = mkt)
}

prep_ml_odds <- function(odds_df, mkt) {
  odds_df %>%
    filter(market == mkt) %>%
    mutate(side = ifelse(outcome_name == home_team, "home", "away")) %>%
    group_by(event_id) %>%
    summarize(
      home_team = first(home_team),
      best_home_odds = max(outcome_price[side == "home"], na.rm = TRUE),
      best_away_odds = max(outcome_price[side == "away"], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(is.finite(best_home_odds), is.finite(best_away_odds)) %>%
    mutate(market = mkt)
}

spreads_h1 <- prep_spread_odds(deriv_odds, "spreads_h1")
spreads_h2 <- prep_spread_odds(deriv_odds, "spreads_h2")
totals_h1  <- prep_total_odds(deriv_odds, "totals_h1")
totals_h2  <- prep_total_odds(deriv_odds, "totals_h2")
ml_h1      <- prep_ml_odds(deriv_odds, "h2h_h1")
ml_h2      <- prep_ml_odds(deriv_odds, "h2h_h2")

cat("H1 spreads:", nrow(spreads_h1), " H2 spreads:", nrow(spreads_h2), "\n")
cat("H1 totals:", nrow(totals_h1), " H2 totals:", nrow(totals_h2), "\n")
cat("H1 ML:", nrow(ml_h1), " H2 ML:", nrow(ml_h2), "\n")

# Get test games (must have at least one derivative market)
games_with_derivs <- unique(c(
  spreads_h1$event_id, spreads_h2$event_id,
  totals_h1$event_id, totals_h2$event_id,
  ml_h1$event_id, ml_h2$event_id
))
DT_test <- DT %>% filter(odds_id %in% games_with_derivs)
cat("Games with derivative odds:", nrow(DT_test), "\n")

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
# WORKER FUNCTION (processes a single game, returns predictions + sample)
# =============================================================================

process_single_game <- function(test_id) {
  preds <- list()

  test_game <- DT_test[DT_test$odds_id == test_id, ][1, ]
  parent_spread <- test_game$home_spread
  parent_total <- test_game$total_line
  target_cover <- test_game$consensus_home_cover_prob
  target_over <- test_game$consensus_over_prob
  game_date <- test_game$game_date

  if (is.na(parent_spread) || is.na(parent_total)) return(NULL)
  if (is.na(target_cover) || is.na(target_over)) return(NULL)

  DT_train <- DT[DT$odds_id != test_id, ]

  sample_result <- tryCatch({
    run_answer_key_sample(
      id = test_id, parent_spread = parent_spread, parent_total = parent_total,
      target_cover = target_cover, target_over = target_over,
      DT = DT_train, ss = ss, st = st, N = N, use_spread_line = TRUE
    )
  }, error = function(e) NULL)

  if (is.null(sample_result)) return(NULL)
  samp <- sample_result$sample

  # --- Helper to add spread predictions ---
  add_spreads <- function(game_odds, margin_col, actual_margin_col) {
    if (nrow(game_odds) == 0) return()
    margins <- samp[[margin_col]]
    actual_margin <- test_game[[actual_margin_col]]
    if (is.na(actual_margin)) return()

    for (j in 1:nrow(game_odds)) {
      line <- game_odds$line[j]
      non_push <- margins != -line
      if (sum(non_push) < 10) next
      algo_home <- mean(margins[non_push] > -line)
      book <- devig_pair(game_odds$best_home_odds[j], game_odds$best_away_odds[j])
      if (actual_margin == -line) next
      actual_home <- ifelse(actual_margin > -line, 1, 0)
      mkt <- game_odds$market[j]

      preds[[length(preds) + 1]] <<- data.frame(
        game_id = test_id, game_date = game_date, market = mkt, line = line,
        side = "home", algo_prob = algo_home, book_prob = book$prob1,
        book_odds = game_odds$best_home_odds[j], actual = actual_home,
        home_team = test_game$home_team, away_team = test_game$away_team,
        stringsAsFactors = FALSE
      )
      preds[[length(preds) + 1]] <<- data.frame(
        game_id = test_id, game_date = game_date, market = mkt, line = line,
        side = "away", algo_prob = 1 - algo_home, book_prob = book$prob2,
        book_odds = game_odds$best_away_odds[j], actual = 1 - actual_home,
        home_team = test_game$home_team, away_team = test_game$away_team,
        stringsAsFactors = FALSE
      )
    }
  }

  # --- Helper to add total predictions ---
  add_totals <- function(game_odds, total_col, actual_total_col) {
    if (nrow(game_odds) == 0) return()
    totals <- samp[[total_col]]
    actual_total <- test_game[[actual_total_col]]
    if (is.na(actual_total)) return()

    for (j in 1:nrow(game_odds)) {
      line <- game_odds$line[j]
      non_push <- totals != line
      if (sum(non_push) < 10) next
      algo_over <- mean(totals[non_push] > line)
      book <- devig_pair(game_odds$best_over_odds[j], game_odds$best_under_odds[j])
      if (actual_total == line) next
      actual_over <- ifelse(actual_total > line, 1, 0)
      mkt <- game_odds$market[j]

      preds[[length(preds) + 1]] <<- data.frame(
        game_id = test_id, game_date = game_date, market = mkt, line = line,
        side = "over", algo_prob = algo_over, book_prob = book$prob1,
        book_odds = game_odds$best_over_odds[j], actual = actual_over,
        home_team = test_game$home_team, away_team = test_game$away_team,
        stringsAsFactors = FALSE
      )
      preds[[length(preds) + 1]] <<- data.frame(
        game_id = test_id, game_date = game_date, market = mkt, line = line,
        side = "under", algo_prob = 1 - algo_over, book_prob = book$prob2,
        book_odds = game_odds$best_under_odds[j], actual = 1 - actual_over,
        home_team = test_game$home_team, away_team = test_game$away_team,
        stringsAsFactors = FALSE
      )
    }
  }

  # --- Helper to add ML predictions ---
  add_ml <- function(game_odds, margin_col, actual_margin_col) {
    if (nrow(game_odds) == 0) return()
    margins <- samp[[margin_col]]
    actual_margin <- test_game[[actual_margin_col]]
    if (is.na(actual_margin) || actual_margin == 0) return()

    non_tie <- margins != 0
    if (sum(non_tie) < 10) return()
    algo_home <- mean(margins[non_tie] > 0)
    book <- devig_pair(game_odds$best_home_odds[1], game_odds$best_away_odds[1])
    actual_home <- ifelse(actual_margin > 0, 1, 0)
    mkt <- game_odds$market[1]

    preds[[length(preds) + 1]] <<- data.frame(
      game_id = test_id, game_date = game_date, market = mkt, line = NA_real_,
      side = "home", algo_prob = algo_home, book_prob = book$prob1,
      book_odds = game_odds$best_home_odds[1], actual = actual_home,
      home_team = test_game$home_team, away_team = test_game$away_team,
      stringsAsFactors = FALSE
    )
    preds[[length(preds) + 1]] <<- data.frame(
      game_id = test_id, game_date = game_date, market = mkt, line = NA_real_,
      side = "away", algo_prob = 1 - algo_home, book_prob = book$prob2,
      book_odds = game_odds$best_away_odds[1], actual = 1 - actual_home,
      home_team = test_game$home_team, away_team = test_game$away_team,
      stringsAsFactors = FALSE
    )
  }

  # H1 markets
  add_spreads(spreads_h1[spreads_h1$event_id == test_id, ], "game_home_margin_h1", "game_home_margin_h1")
  add_totals(totals_h1[totals_h1$event_id == test_id, ], "game_total_h1", "game_total_h1")
  add_ml(ml_h1[ml_h1$event_id == test_id, ], "game_home_margin_h1", "game_home_margin_h1")

  # H2 markets
  add_spreads(spreads_h2[spreads_h2$event_id == test_id, ], "game_home_margin_h2", "game_home_margin_h2")
  add_totals(totals_h2[totals_h2$event_id == test_id, ], "game_total_h2", "game_total_h2")
  add_ml(ml_h2[ml_h2$event_id == test_id, ], "game_home_margin_h2", "game_home_margin_h2")

  if (length(preds) == 0) return(NULL)

  list(
    predictions = do.call(rbind, preds),
    sample = samp
  )
}

# =============================================================================
# RUN PARALLEL BACKTEST
# =============================================================================

cat("===========================================\n")
cat("RUNNING PARALLEL BACKTEST (H1 + H2)\n")
cat("===========================================\n\n")

start_time <- Sys.time()

cl <- makeCluster(N_CORES)
clusterExport(cl, c("DT", "DT_test", "spreads_h1", "spreads_h2",
                     "totals_h1", "totals_h2", "ml_h1", "ml_h2",
                     "ss", "st", "N", "devig_pair", "process_single_game",
                     "american_to_prob", "WORKTREE_TOOLS"))
clusterEvalQ(cl, {
  library(data.table)
  library(dplyr)
  setwd("~/NFLWork/Answer Keys")
  source(WORKTREE_TOOLS)
})

cat("Processing", length(test_ids), "games...\n")

results_list <- parLapplyLB(cl, test_ids, function(id) {
  tryCatch(
    process_single_game(id),
    error = function(e) list(error = TRUE, id = id, msg = conditionMessage(e))
  )
})

stopCluster(cl)

elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat("Completed in", round(elapsed, 1), "minutes\n\n")

# Separate predictions and samples
errors <- sapply(results_list, function(x) is.list(x) && !is.null(x$error))
n_errors <- sum(errors)
cat("Games with errors:", n_errors, "\n")

valid <- results_list[!errors & !sapply(results_list, is.null)]
cat("Successful games:", length(valid), "\n")

all_preds <- bind_rows(lapply(valid, `[[`, "predictions"))
cat("Total predictions:", nrow(all_preds), "\n")

# Build named list of samples keyed by game_id
samples_list <- setNames(
  lapply(valid, function(x) list(sample = x$sample)),
  sapply(valid, function(x) x$predictions$game_id[1])
)

# =============================================================================
# PHASE 2: FILTER TO +EV BETS, COMPUTE INDEPENDENT KELLY
# =============================================================================

cat("\n===========================================\n")
cat("COMPUTING INDEPENDENT KELLY STAKES\n")
cat("===========================================\n\n")

ev_bets <- all_preds %>%
  mutate(
    ev = calc_ev(algo_prob, book_odds),
    decimal_odds = ifelse(book_odds > 0, 1 + book_odds / 100, 1 + 100 / abs(book_odds))
  ) %>%
  filter(ev > EV_THRESHOLD) %>%
  mutate(
    stake_independent = mapply(calc_kelly, ev, book_odds,
                               MoreArgs = list(bankroll = BANKROLL, kelly_mult = KELLY_MULT))
  )

cat("Total +EV bets:", nrow(ev_bets), "\n")
cat("Games with +EV bets:", n_distinct(ev_bets$game_id), "\n")

# Count games with multiple bets
multi_bet_games <- ev_bets %>%
  group_by(game_id) %>%
  filter(n() >= 2) %>%
  ungroup()
cat("Bets on multi-bet games:", nrow(multi_bet_games), "\n")
cat("Multi-bet games:", n_distinct(multi_bet_games$game_id), "\n\n")

# =============================================================================
# PHASE 3: APPLY CORRELATION ADJUSTMENT
# =============================================================================

cat("===========================================\n")
cat("COMPUTING CORRELATION-ADJUSTED STAKES\n")
cat("===========================================\n\n")

# Start with independent stakes as default
ev_bets$stake_corr <- ev_bets$stake_independent
ev_bets$correlation_adj <- 1.0

game_ids <- unique(ev_bets$game_id)
n_mv <- 0L; n_fb <- 0L; n_skip <- 0L

# Verify correlation functions are available
stopifnot("bet_to_leg must be loaded" = exists("bet_to_leg"))
stopifnot("multivariate_kelly must be loaded" = exists("multivariate_kelly"))

for (gid in game_ids) {
  idx <- which(ev_bets$game_id == gid)
  if (length(idx) < 2) next

  if (is.null(samples_list[[gid]])) { n_skip <- n_skip + 1L; next }
  game_sample <- samples_list[[gid]]$sample

  bets_group <- ev_bets[idx, ]

  # Build bet rows compatible with bet_to_leg()
  # bet_to_leg expects: market, bet_on, home_team, away_team, line
  bet_rows <- bets_group %>%
    mutate(
      bet_on = case_when(
        side == "home" ~ home_team,
        side == "away" ~ away_team,
        side == "over" ~ "Over",
        side == "under" ~ "Under",
        TRUE ~ side
      ),
      prob = algo_prob,
      odds = book_odds
    )

  # Try multivariate Kelly
  mv_result <- tryCatch(
    multivariate_kelly(bet_rows, game_sample),
    error = function(e) NULL
  )

  if (!is.null(mv_result)) {
    mv_sizes <- mv_result$f_star * KELLY_MULT * BANKROLL
    mv_sizes <- round(mv_sizes, 2)
    original_sizes <- bets_group$stake_independent
    mv_sizes <- pmin(mv_sizes, original_sizes * 1.5)
    mv_sizes <- pmax(mv_sizes, 0)

    ev_bets$stake_corr[idx] <- mv_sizes
    ev_bets$correlation_adj[idx] <- ifelse(original_sizes > 0, mv_sizes / original_sizes, 1.0)
    n_mv <- n_mv + 1L
  } else {
    # Fallback: per-bet average rho
    n_bets <- nrow(bets_group)
    outcome_matrix <- matrix(NA_real_, nrow = nrow(game_sample), ncol = n_bets)
    leg_failed <- FALSE
    for (i in seq_len(n_bets)) {
      leg <- tryCatch(bet_to_leg(bet_rows[i, ]), error = function(e) NULL)
      if (is.null(leg)) { leg_failed <- TRUE; break }
      outcomes <- tryCatch(evaluate_leg(game_sample, leg), error = function(e) NULL)
      if (is.null(outcomes) || length(outcomes) == 0) { leg_failed <- TRUE; break }
      outcome_matrix[, i] <- as.numeric(outcomes)
    }
    if (leg_failed) { n_skip <- n_skip + 1L; next }

    complete_rows <- complete.cases(outcome_matrix)
    if (sum(complete_rows) < 30) { n_skip <- n_skip + 1L; next }
    R <- cor(outcome_matrix[complete_rows, , drop = FALSE])
    if (any(is.na(R))) { n_skip <- n_skip + 1L; next }

    for (i in seq_len(n_bets)) {
      avg_rho_i <- mean(R[i, -i])
      denom <- 1 + (n_bets - 1) * avg_rho_i
      scale_i <- if (denom <= 0) 1.5 else min(1.5, 1 / sqrt(denom))
      scale_i <- max(scale_i, 0)
      ev_bets$stake_corr[idx[i]] <- round(ev_bets$stake_independent[idx[i]] * scale_i, 2)
      ev_bets$correlation_adj[idx[i]] <- scale_i
    }
    n_fb <- n_fb + 1L
  }
}

cat(sprintf("Correlation adjustment: multivariate=%d, fallback=%d, skipped=%d\n", n_mv, n_fb, n_skip))

adjusted <- ev_bets$correlation_adj[ev_bets$correlation_adj != 1.0]
if (length(adjusted) > 0) {
  cat(sprintf("Adjusted bets: %d, avg scale: %.3f, range: [%.3f, %.3f]\n\n",
              length(adjusted), mean(adjusted), min(adjusted), max(adjusted)))
}

# =============================================================================
# PHASE 4: COMPUTE PNL FOR BOTH STRATEGIES
# =============================================================================

ev_bets <- ev_bets %>%
  mutate(
    pnl_independent = ifelse(actual == 1, stake_independent * (decimal_odds - 1), -stake_independent),
    pnl_corr = ifelse(actual == 1, stake_corr * (decimal_odds - 1), -stake_corr),
    is_multi_bet = game_id %in% multi_bet_games$game_id
  )

# =============================================================================
# PHASE 5: COMPARISON METRICS
# =============================================================================

cat("===========================================\n")
cat("STRATEGY COMPARISON: INDEPENDENT vs CORRELATION-ADJUSTED\n")
cat("===========================================\n\n")

# Overall comparison
summary_independent <- ev_bets %>%
  summarize(
    strategy = "Independent Kelly",
    n_bets = n(),
    total_staked = round(sum(stake_independent)),
    total_pnl = round(sum(pnl_independent)),
    roi = round(sum(pnl_independent) / sum(stake_independent) * 100, 2),
    avg_stake = round(mean(stake_independent), 2),
    win_rate = round(mean(actual) * 100, 1)
  )

summary_corr <- ev_bets %>%
  summarize(
    strategy = "Corr-Adjusted Kelly",
    n_bets = n(),
    total_staked = round(sum(stake_corr)),
    total_pnl = round(sum(pnl_corr)),
    roi = round(sum(pnl_corr) / sum(stake_corr) * 100, 2),
    avg_stake = round(mean(stake_corr), 2),
    win_rate = round(mean(actual) * 100, 1)
  )

cat("--- OVERALL ---\n")
print(bind_rows(summary_independent, summary_corr))

# By market
cat("\n--- BY MARKET ---\n")
by_market <- ev_bets %>%
  group_by(market) %>%
  summarize(
    n_bets = n(),
    roi_independent = round(sum(pnl_independent) / sum(stake_independent) * 100, 2),
    roi_corr = round(sum(pnl_corr) / sum(stake_corr) * 100, 2),
    roi_diff = round(roi_corr - roi_independent, 2),
    avg_corr_adj = round(mean(correlation_adj), 3),
    .groups = "drop"
  ) %>%
  arrange(desc(n_bets))
print(by_market)

# Multi-bet games vs single-bet games
cat("\n--- MULTI-BET GAMES vs SINGLE-BET GAMES ---\n")
by_multi <- ev_bets %>%
  group_by(is_multi_bet) %>%
  summarize(
    label = ifelse(first(is_multi_bet), "Multi-bet games", "Single-bet games"),
    n_bets = n(),
    n_games = n_distinct(game_id),
    staked_ind = round(sum(stake_independent)),
    pnl_ind = round(sum(pnl_independent)),
    roi_ind = round(sum(pnl_independent) / sum(stake_independent) * 100, 2),
    staked_corr = round(sum(stake_corr)),
    pnl_corr = round(sum(pnl_corr)),
    roi_corr = round(sum(pnl_corr) / sum(stake_corr) * 100, 2),
    avg_adj = round(mean(correlation_adj), 3),
    .groups = "drop"
  )
print(by_multi %>% select(-is_multi_bet))

# Correlation adjustment distribution
cat("\n--- CORRELATION ADJUSTMENT DISTRIBUTION ---\n")
adj_dist <- ev_bets %>%
  filter(correlation_adj != 1.0) %>%
  mutate(adj_bucket = cut(correlation_adj,
    breaks = c(0, 0.5, 0.7, 0.9, 1.0, 1.5),
    labels = c("<0.50", "0.50-0.70", "0.70-0.90", "0.90-1.00", ">1.00"),
    include.lowest = TRUE
  )) %>%
  group_by(adj_bucket) %>%
  summarize(
    n_bets = n(),
    avg_ev = round(mean(ev) * 100, 1),
    roi_ind = round(sum(pnl_independent) / sum(stake_independent) * 100, 2),
    roi_corr = round(sum(pnl_corr) / sum(stake_corr) * 100, 2),
    .groups = "drop"
  )
print(adj_dist)

# ROI by number of bets per game
cat("\n--- ROI BY BETS PER GAME ---\n")
bets_per_game <- ev_bets %>%
  group_by(game_id) %>%
  mutate(n_bets_in_game = n()) %>%
  ungroup() %>%
  mutate(bets_bucket = ifelse(n_bets_in_game >= 4, "4+", as.character(n_bets_in_game))) %>%
  group_by(bets_bucket) %>%
  summarize(
    n_bets = n(),
    n_games = n_distinct(game_id),
    roi_ind = round(sum(pnl_independent) / sum(stake_independent) * 100, 2),
    roi_corr = round(sum(pnl_corr) / sum(stake_corr) * 100, 2),
    roi_diff = round(roi_corr - roi_ind, 2),
    avg_adj = round(mean(correlation_adj), 3),
    .groups = "drop"
  )
print(bets_per_game)

# =============================================================================
# PHASE 6: ROLLING BANKROLL SIMULATION
# =============================================================================

cat("\n===========================================\n")
cat("ROLLING BANKROLL SIMULATION\n")
cat("===========================================\n\n")

STARTING_BANKROLL <- 1000

# Order bets chronologically (game_date, then game_id for tie-breaking)
sim <- ev_bets %>% arrange(game_date, game_id)

# Flag best-EV bet per game (for the "best EV only" strategy)
sim <- sim %>%
  group_by(game_id) %>%
  mutate(is_best_ev = ev == max(ev)) %>%
  # Break ties: keep first occurrence only
  mutate(is_best_ev = is_best_ev & cumsum(is_best_ev) == 1) %>%
  ungroup()

cat("Best-EV-per-game bets:", sum(sim$is_best_ev), "\n\n")

# --- Simulate day-by-day: size all bets on a given date using that day's bankroll ---
# (Within a day, bets are simultaneous so bankroll doesn't update intra-day)
sim_dates <- sort(unique(sim$game_date))

bankroll_ind <- STARTING_BANKROLL
bankroll_corr <- STARTING_BANKROLL
bankroll_best <- STARTING_BANKROLL

# Track daily bankroll snapshots
daily_log <- data.frame(
  date = as.Date(character()),
  bankroll_ind = numeric(), bankroll_corr = numeric(), bankroll_best = numeric(),
  daily_pnl_ind = numeric(), daily_pnl_corr = numeric(), daily_pnl_best = numeric(),
  daily_staked_ind = numeric(), daily_staked_corr = numeric(), daily_staked_best = numeric(),
  n_bets = integer(), n_bets_best = integer()
)

sim$rolling_stake_ind <- 0
sim$rolling_stake_corr <- 0
sim$rolling_stake_best <- 0
sim$rolling_pnl_ind <- 0
sim$rolling_pnl_corr <- 0
sim$rolling_pnl_best <- 0

for (d in seq_along(sim_dates)) {
  today <- sim_dates[d]
  day_idx <- which(sim$game_date == today)

  # Size bets as fraction of current bankroll (using the same Kelly fractions)
  for (i in day_idx) {
    ev_i <- sim$ev[i]
    dec_i <- sim$decimal_odds[i]
    b <- dec_i - 1
    edge_frac <- ev_i / b

    # Independent: all bets, independent Kelly
    kelly_ind <- max(edge_frac * KELLY_MULT * bankroll_ind, 0)

    # Corr-adjusted: all bets, scaled by correlation adjustment
    kelly_corr <- max(kelly_ind * (bankroll_corr / bankroll_ind) * sim$correlation_adj[i], 0)
    # Fix: size against corr bankroll, not ind bankroll
    kelly_corr <- max(edge_frac * KELLY_MULT * bankroll_corr * sim$correlation_adj[i], 0)

    # Best-EV-only: only bet if this is the best EV for the game
    kelly_best <- 0
    if (sim$is_best_ev[i]) {
      kelly_best <- max(edge_frac * KELLY_MULT * bankroll_best, 0)
    }

    sim$rolling_stake_ind[i] <- round(kelly_ind, 2)
    sim$rolling_stake_corr[i] <- round(kelly_corr, 2)
    sim$rolling_stake_best[i] <- round(kelly_best, 2)

    sim$rolling_pnl_ind[i] <- ifelse(sim$actual[i] == 1,
      kelly_ind * (dec_i - 1), -kelly_ind)
    sim$rolling_pnl_corr[i] <- ifelse(sim$actual[i] == 1,
      kelly_corr * (dec_i - 1), -kelly_corr)
    sim$rolling_pnl_best[i] <- ifelse(sim$actual[i] == 1,
      kelly_best * (dec_i - 1), -kelly_best)
  }

  day_pnl_ind <- sum(sim$rolling_pnl_ind[day_idx])
  day_pnl_corr <- sum(sim$rolling_pnl_corr[day_idx])
  day_pnl_best <- sum(sim$rolling_pnl_best[day_idx])

  bankroll_ind <- bankroll_ind + day_pnl_ind
  bankroll_corr <- bankroll_corr + day_pnl_corr
  bankroll_best <- bankroll_best + day_pnl_best

  daily_log <- rbind(daily_log, data.frame(
    date = today,
    bankroll_ind = bankroll_ind, bankroll_corr = bankroll_corr, bankroll_best = bankroll_best,
    daily_pnl_ind = day_pnl_ind, daily_pnl_corr = day_pnl_corr, daily_pnl_best = day_pnl_best,
    daily_staked_ind = sum(sim$rolling_stake_ind[day_idx]),
    daily_staked_corr = sum(sim$rolling_stake_corr[day_idx]),
    daily_staked_best = sum(sim$rolling_stake_best[day_idx]),
    n_bets = length(day_idx),
    n_bets_best = sum(sim$is_best_ev[day_idx])
  ))
}

# --- Compute metrics ---
cat("--- FINAL BANKROLL ---\n")
cat(sprintf("Independent:    $%.2f  (%.1f%% return)\n",
    bankroll_ind, (bankroll_ind / STARTING_BANKROLL - 1) * 100))
cat(sprintf("Corr-Adjusted:  $%.2f  (%.1f%% return)\n",
    bankroll_corr, (bankroll_corr / STARTING_BANKROLL - 1) * 100))
cat(sprintf("Best-EV-Only:   $%.2f  (%.1f%% return)\n",
    bankroll_best, (bankroll_best / STARTING_BANKROLL - 1) * 100))

# Max drawdown
calc_max_drawdown <- function(bankroll_series) {
  peak <- cummax(bankroll_series)
  drawdown <- (peak - bankroll_series) / peak
  max(drawdown) * 100
}

dd_ind <- calc_max_drawdown(daily_log$bankroll_ind)
dd_corr <- calc_max_drawdown(daily_log$bankroll_corr)
dd_best <- calc_max_drawdown(daily_log$bankroll_best)

cat(sprintf("\n--- MAX DRAWDOWN ---\n"))
cat(sprintf("Independent:    %.2f%%\n", dd_ind))
cat(sprintf("Corr-Adjusted:  %.2f%%\n", dd_corr))
cat(sprintf("Best-EV-Only:   %.2f%%\n", dd_best))

# Sharpe ratio (daily PnL / daily staked, annualized)
daily_roi_ind <- daily_log$daily_pnl_ind / pmax(daily_log$daily_staked_ind, 1)
daily_roi_corr <- daily_log$daily_pnl_corr / pmax(daily_log$daily_staked_corr, 1)
daily_roi_best <- daily_log$daily_pnl_best / pmax(daily_log$daily_staked_best, 1)

sharpe_ind <- mean(daily_roi_ind) / sd(daily_roi_ind) * sqrt(252)
sharpe_corr <- mean(daily_roi_corr) / sd(daily_roi_corr) * sqrt(252)
sharpe_best <- mean(daily_roi_best) / sd(daily_roi_best) * sqrt(252)

cat(sprintf("\n--- SHARPE RATIO (annualized on daily ROI) ---\n"))
cat(sprintf("Independent:    %.3f\n", sharpe_ind))
cat(sprintf("Corr-Adjusted:  %.3f\n", sharpe_corr))
cat(sprintf("Best-EV-Only:   %.3f\n", sharpe_best))

# Total staked and ROI
total_staked_ind <- sum(sim$rolling_stake_ind)
total_staked_corr <- sum(sim$rolling_stake_corr)
total_staked_best <- sum(sim$rolling_stake_best)
total_pnl_ind <- sum(sim$rolling_pnl_ind)
total_pnl_corr <- sum(sim$rolling_pnl_corr)
total_pnl_best <- sum(sim$rolling_pnl_best)

cat(sprintf("\n--- ROLLING TOTALS ---\n"))
cat(sprintf("Independent:    Staked $%.0f, PnL $%.0f, ROI %.2f%%\n",
    total_staked_ind, total_pnl_ind, total_pnl_ind / total_staked_ind * 100))
cat(sprintf("Corr-Adjusted:  Staked $%.0f, PnL $%.0f, ROI %.2f%%\n",
    total_staked_corr, total_pnl_corr, total_pnl_corr / total_staked_corr * 100))
cat(sprintf("Best-EV-Only:   Staked $%.0f, PnL $%.0f, ROI %.2f%%\n",
    total_staked_best, total_pnl_best, total_pnl_best / total_staked_best * 100))

# Worst single-day loss
worst_ind <- min(daily_log$daily_pnl_ind)
worst_corr <- min(daily_log$daily_pnl_corr)
worst_best <- min(daily_log$daily_pnl_best)
worst_ind_date <- daily_log$date[which.min(daily_log$daily_pnl_ind)]
worst_corr_date <- daily_log$date[which.min(daily_log$daily_pnl_corr)]
worst_best_date <- daily_log$date[which.min(daily_log$daily_pnl_best)]

cat(sprintf("\n--- WORST SINGLE-DAY LOSS ---\n"))
cat(sprintf("Independent:    $%.2f on %s\n", worst_ind, worst_ind_date))
cat(sprintf("Corr-Adjusted:  $%.2f on %s\n", worst_corr, worst_corr_date))
cat(sprintf("Best-EV-Only:   $%.2f on %s\n", worst_best, worst_best_date))

# Longest losing streak (consecutive negative-PnL days)
longest_streak <- function(pnl) {
  runs <- rle(pnl < 0)
  max(runs$lengths[runs$values], 0)
}
streak_ind <- longest_streak(daily_log$daily_pnl_ind)
streak_corr <- longest_streak(daily_log$daily_pnl_corr)
streak_best <- longest_streak(daily_log$daily_pnl_best)

cat(sprintf("\n--- LONGEST LOSING STREAK (days) ---\n"))
cat(sprintf("Independent:    %d days\n", streak_ind))
cat(sprintf("Corr-Adjusted:  %d days\n", streak_corr))
cat(sprintf("Best-EV-Only:   %d days\n", streak_best))

# Bankroll below starting at any point?
below_start_ind <- sum(daily_log$bankroll_ind < STARTING_BANKROLL)
below_start_corr <- sum(daily_log$bankroll_corr < STARTING_BANKROLL)
below_start_best <- sum(daily_log$bankroll_best < STARTING_BANKROLL)
cat(sprintf("\n--- DAYS BELOW STARTING BANKROLL ---\n"))
cat(sprintf("Independent:    %d / %d days (%.1f%%)\n",
    below_start_ind, nrow(daily_log), below_start_ind / nrow(daily_log) * 100))
cat(sprintf("Corr-Adjusted:  %d / %d days (%.1f%%)\n",
    below_start_corr, nrow(daily_log), below_start_corr / nrow(daily_log) * 100))
cat(sprintf("Best-EV-Only:   %d / %d days (%.1f%%)\n",
    below_start_best, nrow(daily_log), below_start_best / nrow(daily_log) * 100))

cat("\n===========================================\n")
cat("BACKTEST COMPLETE\n")
cat("===========================================\n")
