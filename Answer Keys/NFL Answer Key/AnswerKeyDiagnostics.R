# Answer Key Algorithm Diagnostics
# Purpose: Verify that mean_match and balance_sample produce correct outputs
# Run this after loading your data and sourcing Tools.R

setwd("~/NFLWork/Answer Keys")
library(data.table)
library(duckdb)
library(dplyr)
library(tidyverse)
library(DBI)
source("Tools.R")

# =============================================================================
# STEP 1: Load Data (same as NFLAnswerKey2.0.R)
# =============================================================================

con <- dbConnect(duckdb(), dbdir = "pbp.duckdb") 

DT <- dbGetQuery(con, "SELECT * FROM nfl_betting_pbp") %>% 
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

dbDisconnect(con)

# Compute dispersion for scaling
disp <- compute_dispersion(DT, moneyline = FALSE)
ss <- disp$ss
st <- disp$st

cat("=== DATA OVERVIEW ===\n")
cat("Total games in DT:", nrow(DT), "\n")
cat("Spread range (5th-95th percentile):", ss, "\n")
cat("Total range (5th-95th percentile):", st, "\n")
cat("Spread distribution:\n")
print(summary(DT$home_spread))
cat("Total distribution:\n")
print(summary(DT$total_line))
cat("\n")

# =============================================================================
# STEP 2: Diagnostic Functions
# =============================================================================

#' Diagnose mean_match output
diagnose_mean_match <- function(mm_result, parent_spread, parent_total, 
                                 use_spread_line = FALSE) {
  dt <- mm_result$dt
  inc <- dt[included == TRUE]
  spread_col <- if (use_spread_line) "home_spread" else "home_ml_odds"
  
  tibble(
    # Sample size
    n_included = nrow(inc),
    
    # Spread statistics
    target_spread = parent_spread,
    sample_mean_spread = mean(inc[[spread_col]], na.rm = TRUE),
    sample_median_spread = median(inc[[spread_col]], na.rm = TRUE),
    sample_sd_spread = sd(inc[[spread_col]], na.rm = TRUE),
    spread_mean_error = mean(inc[[spread_col]], na.rm = TRUE) - parent_spread,
    spread_median_error = median(inc[[spread_col]], na.rm = TRUE) - parent_spread,
    
    # Total statistics
    target_total = parent_total,
    sample_mean_total = mean(inc$total_line, na.rm = TRUE),
    sample_median_total = median(inc$total_line, na.rm = TRUE),
    sample_sd_total = sd(inc$total_line, na.rm = TRUE),
    total_mean_error = mean(inc$total_line, na.rm = TRUE) - parent_total,
    total_median_error = median(inc$total_line, na.rm = TRUE) - parent_total,
    
    # Index statistics (lower = better match)
    min_index = min(inc$index, na.rm = TRUE),
    max_index = max(inc$index, na.rm = TRUE),
    median_index = median(inc$index, na.rm = TRUE),
    
    # Outcome rates (pre-balancing)
    cover_rate = mean(inc$actual_cover, na.rm = TRUE),
    over_rate = mean(inc$actual_over, na.rm = TRUE)
  )
}

#' Diagnose balance_sample output
diagnose_balance_sample <- function(bal_result, target_cover, target_over) {
  dt <- bal_result$dt
  inc <- dt[included == TRUE]
  n <- nrow(inc)
  
  tibble(
    # Sample size
    final_N = bal_result$final_N,
    actual_included = n,
    
    # Cover statistics
    target_cover = target_cover,
    expected_covers = round(target_cover * n),
    actual_covers = sum(inc$actual_cover, na.rm = TRUE),
    actual_cover_rate = mean(inc$actual_cover, na.rm = TRUE),
    cover_error = bal_result$cover_error,
    cover_rate_error = mean(inc$actual_cover, na.rm = TRUE) - target_cover,
    
    # Over statistics
    target_over = target_over,
    expected_overs = round(target_over * n),
    actual_overs = sum(inc$actual_over, na.rm = TRUE),
    actual_over_rate = mean(inc$actual_over, na.rm = TRUE),
    over_error = bal_result$over_error,
    over_rate_error = mean(inc$actual_over, na.rm = TRUE) - target_over,
    
    # Index statistics (should still be reasonably low)
    min_index = min(inc$index, na.rm = TRUE),
    max_index = max(inc$index, na.rm = TRUE),
    median_index = median(inc$index, na.rm = TRUE),
    
    # Check for pathological samples
    converged = abs(bal_result$cover_error) <= 1 & abs(bal_result$over_error) <= 1
  )
}

#' Run full diagnostic on a single parent line
run_full_diagnostic <- function(DT, parent_spread, parent_total, 
                                 target_cover, target_over,
                                 ss, st, N,
                                 max_iter_mean = 500, tol_mean = 0.005, tol_error = 1,
                                 use_spread_line = TRUE,
                                 verbose = TRUE) {
  
  if (verbose) {
    cat("\n=== DIAGNOSTIC FOR PARENT LINE ===\n")
    cat("Parent Spread:", parent_spread, "\n")
    cat("Parent Total:", parent_total, "\n")
    cat("Target Cover Rate:", target_cover, "\n")
    cat("Target Over Rate:", target_over, "\n")
    cat("Sample Size N:", N, "\n")
    cat("Scaling: ss =", ss, ", st =", st, "\n")
    cat("\n")
  }
  
  # Step 1: Mean Match
  if (verbose) cat("--- Running mean_match ---\n")
  mm_result <- mean_match(
    dt = DT,
    N = N,
    parent_spread = parent_spread,
    parent_total = parent_total,
    ss = ss,
    st = st,
    max_iter_mean = max_iter_mean,
    tol_mean = tol_mean,
    use_spread_line = use_spread_line
  )
  
  mm_diag <- diagnose_mean_match(mm_result, parent_spread, parent_total, use_spread_line)
  
  if (verbose) {
    cat("Mean Match Results:\n")
    cat("  Sample size:", mm_diag$n_included, "\n")
    cat("  Spread - Target:", mm_diag$target_spread, 
        "Mean:", round(mm_diag$sample_mean_spread, 3),
        "Median:", round(mm_diag$sample_median_spread, 3),
        "Error:", round(mm_diag$spread_mean_error, 4), "\n")
    cat("  Total - Target:", mm_diag$target_total,
        "Mean:", round(mm_diag$sample_mean_total, 3),
        "Median:", round(mm_diag$sample_median_total, 3),
        "Error:", round(mm_diag$total_mean_error, 4), "\n")
    cat("  Pre-balance Cover Rate:", round(mm_diag$cover_rate, 4), "\n")
    cat("  Pre-balance Over Rate:", round(mm_diag$over_rate, 4), "\n")
    cat("\n")
  }
  
  # Step 2: Balance Sample
  if (verbose) cat("--- Running balance_sample ---\n")
  bal_result <- balance_sample(
    dt = mm_result$dt,
    N = N,
    target_cover = target_cover,
    target_over = target_over,
    tol_error = tol_error
  )
  
  bal_diag <- diagnose_balance_sample(bal_result, target_cover, target_over)
  
  if (verbose) {
    cat("Balance Sample Results:\n")
    cat("  Final N:", bal_diag$final_N, "(started with", N, ")\n")
    cat("  Cover - Target:", round(bal_diag$target_cover, 4),
        "Actual:", round(bal_diag$actual_cover_rate, 4),
        "Error:", bal_diag$cover_error, "games\n")
    cat("  Over - Target:", round(bal_diag$target_over, 4),
        "Actual:", round(bal_diag$actual_over_rate, 4),
        "Error:", bal_diag$over_error, "games\n")
    cat("  Converged:", bal_diag$converged, "\n")
    cat("\n")
  }
  
  # Step 3: Check spread/total drift after balancing
  inc_after <- bal_result$dt[included == TRUE]
  spread_col <- if (use_spread_line) "home_spread" else "home_ml_odds"
  
  post_balance_diag <- tibble(
    post_mean_spread = mean(inc_after[[spread_col]], na.rm = TRUE),
    post_median_spread = median(inc_after[[spread_col]], na.rm = TRUE),
    post_mean_total = mean(inc_after$total_line, na.rm = TRUE),
    post_median_total = median(inc_after$total_line, na.rm = TRUE),
    spread_drift = mean(inc_after[[spread_col]], na.rm = TRUE) - parent_spread,
    total_drift = mean(inc_after$total_line, na.rm = TRUE) - parent_total
  )
  
  if (verbose) {
    cat("--- Post-Balance Spread/Total Check ---\n")
    cat("  Spread Drift (mean):", round(post_balance_diag$spread_drift, 4), "\n")
    cat("  Total Drift (mean):", round(post_balance_diag$total_drift, 4), "\n")
    cat("\n")
  }
  
  # Return all diagnostics
  list(
    mm_result = mm_result,
    mm_diag = mm_diag,
    bal_result = bal_result,
    bal_diag = bal_diag,
    post_balance_diag = post_balance_diag
  )
}

#' Debug a specific prediction with full traceability
#' 
#' @param DT Historical data
#' @param parent_spread The spread for the game
#' @param parent_total The total for the game
#' @param target_cover De-vigged probability home covers spread
#' @param target_over De-vigged probability game goes over
#' @param margin_col Prefix for margin columns (e.g., "game_home_margin_period")
#' @param periods Which periods to analyze (e.g., c(1,2,3,4) for quarters)
#' @param home_team Optional team name for display
#' @param away_team Optional team name for display
#' 
#' @return List with full diagnostic info
debug_prediction <- function(DT, parent_spread, parent_total, 
                              target_cover, target_over,
                              ss, st, N,
                              margin_col = "game_home_margin_period",
                              periods = c(1, 2, 3, 4),
                              home_team = "Home",
                              away_team = "Away",
                              use_spread_line = TRUE) {
  
  cat("=============================================================\n")
  cat("DEBUG PREDICTION:", home_team, "vs", away_team, "\n")
  cat("=============================================================\n\n")
  
  # --- Input Parameters ---
  cat("--- INPUT PARAMETERS ---\n")
  cat("Parent Spread:", parent_spread, 
      ifelse(parent_spread < 0, "(home favored)", "(home underdog)"), "\n")
  cat("Parent Total:", parent_total, "\n")
  cat("Target Cover Prob:", round(target_cover, 4), 
      sprintf("(%.1f%% home covers)", target_cover * 100), "\n")
  cat("Target Over Prob:", round(target_over, 4), 
      sprintf("(%.1f%% goes over)", target_over * 100), "\n")
  cat("Sample Size N:", N, "\n")
  cat("Scaling - ss:", round(ss, 3), "st:", round(st, 3), "\n\n")
  
  # --- Step 1: Mean Match ---
  cat("--- STEP 1: MEAN MATCHING ---\n")
  mm <- mean_match(DT, N, parent_spread, parent_total, ss, st, 
                   max_iter_mean = 500, tol_mean = 0.005, 
                   use_spread_line = use_spread_line)
  
  inc_mm <- mm$dt[included == TRUE]
  spread_col <- if (use_spread_line) "home_spread" else "home_ml_odds"
  
  mm_stats <- tibble(
    sample_n = nrow(inc_mm),
    mean_spread = mean(inc_mm[[spread_col]], na.rm = TRUE),
    median_spread = median(inc_mm[[spread_col]], na.rm = TRUE),
    mean_total = mean(inc_mm$total_line, na.rm = TRUE),
    median_total = median(inc_mm$total_line, na.rm = TRUE),
    pre_balance_cover_rate = mean(inc_mm$actual_cover, na.rm = TRUE),
    pre_balance_over_rate = mean(inc_mm$actual_over, na.rm = TRUE)
  )
  
  cat("Sample size:", mm_stats$sample_n, "\n")
  cat("Spread - Mean:", round(mm_stats$mean_spread, 3), 
      "Median:", round(mm_stats$median_spread, 3),
      "Target:", parent_spread, "\n")
  cat("Total - Mean:", round(mm_stats$mean_total, 3),
      "Median:", round(mm_stats$median_total, 3),
      "Target:", parent_total, "\n")
  cat("Pre-balance cover rate:", round(mm_stats$pre_balance_cover_rate, 4), 
      sprintf("(%.1f%%)", mm_stats$pre_balance_cover_rate * 100), "\n")
  cat("Pre-balance over rate:", round(mm_stats$pre_balance_over_rate, 4),
      sprintf("(%.1f%%)", mm_stats$pre_balance_over_rate * 100), "\n\n")
  
  # --- Step 2: Balance Sample ---
  cat("--- STEP 2: BALANCE SAMPLE ---\n")
  bal <- balance_sample(mm$dt, N, target_cover, target_over, tol_error = 1)
  
  inc_bal <- bal$dt[included == TRUE]
  
  bal_stats <- tibble(
    final_n = nrow(inc_bal),
    final_cover_rate = mean(inc_bal$actual_cover, na.rm = TRUE),
    final_over_rate = mean(inc_bal$actual_over, na.rm = TRUE),
    cover_error = bal$cover_error,
    over_error = bal$over_error,
    post_mean_spread = mean(inc_bal[[spread_col]], na.rm = TRUE),
    post_median_spread = median(inc_bal[[spread_col]], na.rm = TRUE),
    post_mean_total = mean(inc_bal$total_line, na.rm = TRUE),
    post_median_total = median(inc_bal$total_line, na.rm = TRUE)
  )
  
  cat("Final sample size:", bal_stats$final_n, "\n")
  cat("Cover rate - Target:", round(target_cover, 4), 
      "Actual:", round(bal_stats$final_cover_rate, 4),
      "Error:", bal_stats$cover_error, "games\n")
  cat("Over rate - Target:", round(target_over, 4),
      "Actual:", round(bal_stats$final_over_rate, 4),
      "Error:", bal_stats$over_error, "games\n")
  cat("Post-balance spread - Mean:", round(bal_stats$post_mean_spread, 3),
      "Median:", round(bal_stats$post_median_spread, 3), "\n")
  cat("Post-balance total - Mean:", round(bal_stats$post_mean_total, 3),
      "Median:", round(bal_stats$post_median_total, 3), "\n\n")
  
  # --- Step 3: Period-by-Period Analysis ---
  cat("--- STEP 3: PERIOD-BY-PERIOD PREDICTIONS ---\n")
  
  period_results <- list()
  
  for (p in periods) {
    col_name <- paste0(margin_col, "_", p)
    
    if (!col_name %in% names(inc_bal)) {
      cat("WARNING: Column", col_name, "not found!\n")
      next
    }
    
    margins <- inc_bal[[col_name]]
    
    # Calculate stats
    home_wins <- sum(margins > 0, na.rm = TRUE)
    away_wins <- sum(margins < 0, na.rm = TRUE)
    ties <- sum(margins == 0, na.rm = TRUE)
    valid_games <- sum(margins != 0, na.rm = TRUE)
    
    home_win_prob <- home_wins / valid_games
    
    period_results[[paste0("period_", p)]] <- tibble(
      period = p,
      home_wins = home_wins,
      away_wins = away_wins,
      ties = ties,
      valid_games = valid_games,
      home_win_prob = home_win_prob,
      away_win_prob = 1 - home_win_prob
    )
    
    cat(sprintf("Period %d: Home wins %d/%d (%.1f%%), Away wins %d/%d (%.1f%%), Ties: %d\n",
                p, home_wins, valid_games, home_win_prob * 100,
                away_wins, valid_games, (1 - home_win_prob) * 100, ties))
  }
  
  period_df <- bind_rows(period_results)
  
  # --- Step 4: Sanity Checks ---
  cat("\n--- SANITY CHECKS ---\n")
  
  # Check 1: Full game margin consistency
  if ("home_margin" %in% names(inc_bal)) {
    fg_home_wins <- sum(inc_bal$home_margin > 0, na.rm = TRUE)
    fg_valid <- sum(inc_bal$home_margin != 0, na.rm = TRUE)
    fg_home_win_pct <- fg_home_wins / fg_valid
    
    cat("Full game home win rate:", round(fg_home_win_pct, 4), 
        sprintf("(%.1f%%)\n", fg_home_win_pct * 100))
    
    # For an underdog (+spread), full game win rate should be < 50%
    if (parent_spread > 0 && fg_home_win_pct > 0.50) {
      cat("⚠️  WARNING: Home is underdog (+", parent_spread, 
          ") but wins full game ", round(fg_home_win_pct * 100, 1), 
          "% - seems inconsistent!\n", sep = "")
    }
    # For a favorite (-spread), full game win rate should be > 50%
    if (parent_spread < 0 && fg_home_win_pct < 0.50) {
      cat("⚠️  WARNING: Home is favorite (", parent_spread,
          ") but wins full game only ", round(fg_home_win_pct * 100, 1),
          "% - seems inconsistent!\n", sep = "")
    }
  }
  
  # Check 2: Period win rates vs spread
  cat("\nExpected pattern for home ", 
      ifelse(parent_spread > 0, "UNDERDOG", "FAVORITE"), ":\n", sep = "")
  if (parent_spread > 0) {
    cat("  - Quarter win probabilities should generally be < 50%\n")
  } else {
    cat("  - Quarter win probabilities should generally be > 50%\n")
  }
  
  inconsistent_periods <- c()
  for (p in periods) {
    prob <- period_df$home_win_prob[period_df$period == p]
    if (length(prob) > 0) {
      if (parent_spread > 0 && prob > 0.52) {
        inconsistent_periods <- c(inconsistent_periods, p)
      }
      if (parent_spread < -3 && prob < 0.48) {
        inconsistent_periods <- c(inconsistent_periods, p)
      }
    }
  }
  
  if (length(inconsistent_periods) > 0) {
    cat("⚠️  Potentially inconsistent periods:", paste(inconsistent_periods, collapse = ", "), "\n")
  } else {
    cat("✓ All period predictions appear directionally consistent\n")
  }
  
  cat("\n=============================================================\n\n")
  
  # Return everything for further inspection
  list(
    inputs = tibble(
      home_team = home_team,
      away_team = away_team,
      parent_spread = parent_spread,
      parent_total = parent_total,
      target_cover = target_cover,
      target_over = target_over,
      N = N
    ),
    mean_match_stats = mm_stats,
    balance_stats = bal_stats,
    period_predictions = period_df,
    sample = inc_bal  # The actual sample for manual inspection
  )
}

# =============================================================================
# STEP 3: Run Diagnostics on Test Cases
# =============================================================================

N <- round(nrow(DT) * 0.1, 0)
cat("Using N =", N, "games per sample\n\n")

# --- Test Case 1: Pickem with average total ---
cat("========================================\n")
cat("TEST CASE 1: Pickem (spread = 0), Total = 44\n")
cat("========================================\n")

test1 <- run_full_diagnostic(
  DT = DT,
  parent_spread = 0,
  parent_total = 44,
  target_cover = 0.50,  # Pickem should be 50/50
  target_over = 0.50,
  ss = ss, st = st, N = N,
  use_spread_line = TRUE
)

# --- Test Case 2: Favorite with higher total ---
cat("========================================\n")
cat("TEST CASE 2: Home Favorite -7, Total = 47\n")
cat("========================================\n")

test2 <- run_full_diagnostic(
  DT = DT,
  parent_spread = 3,
  parent_total = 48,
  target_cover = 0.51,  # Slight favorite edge
  target_over = 0.50,
  ss = ss, st = st, N = N,
  use_spread_line = TRUE
)

# --- Test Case 3: Underdog with lower total ---
cat("========================================\n")
cat("TEST CASE 3: Home Underdog +3, Total = 48\n")
cat("========================================\n")

test3 <- run_full_diagnostic(
  DT = DT,
  parent_spread = 3,
  parent_total = 48,
  target_cover = 0.5,
  target_over = 0.5,
  ss = ss, st = st, N = N,
  use_spread_line = TRUE
)

# --- Test Case 4: Extreme case ---
cat("========================================\n")
cat("TEST CASE 4: Heavy Favorite -14, Total = 52\n")
cat("========================================\n")

test4 <- run_full_diagnostic(
  DT = DT,
  parent_spread = -14,
  parent_total = 52,
  target_cover = 0.52,
  target_over = 0.52,
  ss = ss, st = st, N = N,
  use_spread_line = TRUE
)

# =============================================================================
# STEP 4: Summary Table
# =============================================================================

cat("\n========================================\n")
cat("SUMMARY OF ALL TEST CASES\n")
cat("========================================\n\n")

summary_table <- bind_rows(
  test1$bal_diag %>% mutate(test_case = "Pickem/44", parent_spread = 0, parent_total = 44),
  test2$bal_diag %>% mutate(test_case = "Fav-7/47", parent_spread = -7, parent_total = 47),
  test3$bal_diag %>% mutate(test_case = "Dog+3/41", parent_spread = 3, parent_total = 41),
  test4$bal_diag %>% mutate(test_case = "Fav-14/52", parent_spread = -14, parent_total = 52)
) %>%
  select(test_case, parent_spread, parent_total, final_N, 
         target_cover, actual_cover_rate, cover_error,
         target_over, actual_over_rate, over_error, converged)

print(summary_table)

# Add post-balance spread/total drift
drift_table <- bind_rows(
  test1$post_balance_diag %>% mutate(test_case = "Pickem/44"),
  test2$post_balance_diag %>% mutate(test_case = "Fav-7/47"),
  test3$post_balance_diag %>% mutate(test_case = "Dog+3/41"),
  test4$post_balance_diag %>% mutate(test_case = "Fav-14/52")
) %>%
  select(test_case, spread_drift, total_drift)

cat("\nPost-Balance Spread/Total Drift:\n")
print(drift_table)

# =============================================================================
# STEP 5: Visual Distribution Check (Optional)
# =============================================================================

# Check the distribution of included games vs. parent line
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  
  # Use test1 as example
  inc_games <- test1$bal_result$dt[included == TRUE]
  
  p1 <- ggplot(inc_games, aes(x = home_spread)) +
    geom_histogram(binwidth = 0.5, fill = "steelblue", alpha = 0.7) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1) +
    labs(title = "Test 1: Sample Spread Distribution",
         subtitle = "Red line = target spread (0)",
         x = "Home Spread", y = "Count") +
    theme_minimal()
  
  p2 <- ggplot(inc_games, aes(x = total_line)) +
    geom_histogram(binwidth = 0.5, fill = "darkgreen", alpha = 0.7) +
    geom_vline(xintercept = 44, color = "red", linetype = "dashed", size = 1) +
    labs(title = "Test 1: Sample Total Distribution",
         subtitle = "Red line = target total (44)",
         x = "Total Line", y = "Count") +
    theme_minimal()
  
  print(p1)
  print(p2)
}

# =============================================================================
# STEP 6: Feustel Median Check (Critical Validation)
# =============================================================================

cat("\n========================================\n")
cat("FEUSTEL MEDIAN CHECK\n")
cat("(Feustel says to match MEDIAN, not mean)\n")
cat("========================================\n\n")

median_check <- bind_rows(
  test1$mm_diag %>% mutate(test_case = "Pickem/44"),
  test2$mm_diag %>% mutate(test_case = "Fav-7/47"),
  test3$mm_diag %>% mutate(test_case = "Dog+3/41"),
  test4$mm_diag %>% mutate(test_case = "Fav-14/52")
) %>%
  select(test_case, target_spread, sample_mean_spread, sample_median_spread,
         spread_mean_error, spread_median_error,
         target_total, sample_mean_total, sample_median_total,
         total_mean_error, total_median_error)

print(median_check)

cat("\n*** NOTE: If spread_median_error or total_median_error is large,")
cat("\n*** consider switching mean_match to use median instead of mean. ***\n")

# =============================================================================
# STEP 7: Debug Specific Predictions
# =============================================================================
# Use this section to debug specific games that look suspicious

#' Debug the prediction pipeline for a specific game
#' This traces through run_means_for_id to see exactly what's happening
#' 
debug_prediction_pipeline <- function(DT, parent_spread, parent_total,
                                       target_cover, target_over,
                                       ss, st, N,
                                       margin_col = "game_home_margin_period",
                                       use_spread_line = TRUE,
                                       home_team = "Home",
                                       away_team = "Away") {
  
  cat("=============================================================\n")
  cat("PIPELINE DEBUG:", home_team, "vs", away_team, "\n")
  cat("=============================================================\n\n")
  
  # --- Step 0: Verify column names exist ---
  cat("--- STEP 0: VERIFY DATA COLUMNS ---\n")
  margin_cols_in_dt <- names(DT)[startsWith(names(DT), margin_col)]
  cat("Margin columns in DT:\n")
  print(margin_cols_in_dt)
  
  if (length(margin_cols_in_dt) == 0) {
    cat("ERROR: No columns starting with '", margin_col, "' found!\n", sep = "")
    cat("Available columns:\n")
    print(names(DT))
    return(NULL)
  }
  
  # Check for NAs in margin columns
  for (col in margin_cols_in_dt) {
    na_count <- sum(is.na(DT[[col]]))
    cat(sprintf("  %s: %d NAs (%.1f%%)\n", col, na_count, 100 * na_count / nrow(DT)))
  }
  cat("\n")
  
  # --- Step 1: Mean Match ---
  cat("--- STEP 1: MEAN MATCHING ---\n")
  cat("Parent Spread:", parent_spread, ifelse(parent_spread > 0, "(underdog)", "(favorite)"), "\n")
  cat("Parent Total:", parent_total, "\n")
  cat("Target Cover:", target_cover, "\n")
  cat("Target Over:", target_over, "\n")
  cat("N:", N, "\n\n")
  
  mm <- mean_match(DT, N, parent_spread, parent_total, ss, st,
                   max_iter_mean = 500, tol_mean = 0.005,
                   use_spread_line = use_spread_line)
  
  inc_mm <- mm$dt[included == TRUE]
  spread_col <- if (use_spread_line) "home_spread" else "home_ml_odds"
  
  cat("Mean match sample size:", nrow(inc_mm), "\n")
  cat("Sample mean spread:", round(mean(inc_mm[[spread_col]]), 3), "\n")
  cat("Sample mean total:", round(mean(inc_mm$total_line), 3), "\n\n")
  
  # --- Step 2: Balance Sample ---
  cat("--- STEP 2: BALANCE SAMPLE ---\n")
  bal <- balance_sample(mm$dt, N, target_cover, target_over, tol_error = 1)
  inc_bal <- bal$dt[included == TRUE]
  
  cat("Balanced sample size:", nrow(inc_bal), "\n")
  cat("Final cover rate:", round(mean(inc_bal$actual_cover), 4), "\n")
  cat("Final over rate:", round(mean(inc_bal$actual_over), 4), "\n\n")
  
  # --- Step 3: What run_means_for_id produces ---
  cat("--- STEP 3: RUN_MEANS_FOR_ID OUTPUT ---\n")
  cat("This is what run_means_for_id returns:\n\n")
  
  # Replicate exactly what run_means_for_id does
  bets_summary <- inc_bal %>%
    summarise(across(starts_with(margin_col), ~ sum(.x > 0, na.rm = TRUE) / sum(.x != 0, na.rm = TRUE))) %>%
    ungroup()
  
  print(bets_summary)
  cat("\n")
  
  # --- Step 4: Manual verification ---
  cat("--- STEP 4: MANUAL VERIFICATION ---\n")
  cat("Manually counting from balanced sample:\n\n")
  
  for (col in margin_cols_in_dt) {
    if (col %in% names(inc_bal)) {
      vals <- inc_bal[[col]]
      home_wins <- sum(vals > 0, na.rm = TRUE)
      away_wins <- sum(vals < 0, na.rm = TRUE)
      ties <- sum(vals == 0, na.rm = TRUE)
      valid <- sum(vals != 0, na.rm = TRUE)
      home_pct <- home_wins / valid
      
      cat(sprintf("%s:\n", col))
      cat(sprintf("  Home wins: %d, Away wins: %d, Ties: %d, Valid: %d\n", 
                  home_wins, away_wins, ties, valid))
      cat(sprintf("  Home win prob: %.4f (%.1f%%)\n", home_pct, home_pct * 100))
      
      # Sanity check
      if (parent_spread > 0 && home_pct > 0.52) {
        cat("  ⚠️  WARNING: Underdog showing >52% - suspicious!\n")
      }
      if (parent_spread < -3 && home_pct < 0.48) {
        cat("  ⚠️  WARNING: Favorite showing <48% - suspicious!\n")
      }
      cat("\n")
    }
  }
  
  # --- Step 5: Check full game outcome ---
  cat("--- STEP 5: FULL GAME CONSISTENCY CHECK ---\n")
  
  if ("home_margin" %in% names(inc_bal)) {
    fg_wins <- sum(inc_bal$home_margin > 0, na.rm = TRUE)
    fg_losses <- sum(inc_bal$home_margin < 0, na.rm = TRUE)
    fg_valid <- sum(inc_bal$home_margin != 0, na.rm = TRUE)
    fg_pct <- fg_wins / fg_valid
    
    cat("Full game home wins:", fg_wins, "/", fg_valid, 
        sprintf("(%.1f%%)\n", fg_pct * 100))
    
    # Compare to quarter averages
    avg_q_pct <- mean(sapply(margin_cols_in_dt, function(col) {
      if (col %in% names(inc_bal)) {
        vals <- inc_bal[[col]]
        sum(vals > 0, na.rm = TRUE) / sum(vals != 0, na.rm = TRUE)
      } else NA
    }), na.rm = TRUE)
    
    cat("Avg quarter home win %:", sprintf("%.1f%%\n", avg_q_pct * 100))
    cat("Full game home win %:", sprintf("%.1f%%\n", fg_pct * 100))
    
    if (abs(avg_q_pct - fg_pct) > 0.05) {
      cat("⚠️  Large gap between quarter avg and full game!\n")
    }
  }
  
  # --- Step 6: Sample some actual games ---
  cat("\n--- STEP 6: SAMPLE GAMES FROM BALANCED SET ---\n")
  cat("First 10 games in balanced sample:\n\n")
  
  sample_cols <- c("home_spread", "total_line", "home_margin", "actual_cover", margin_cols_in_dt[1:min(2, length(margin_cols_in_dt))])
  sample_cols <- sample_cols[sample_cols %in% names(inc_bal)]
  
  print(inc_bal %>% select(all_of(sample_cols)) %>% head(10))
  
  cat("\n=============================================================\n\n")
  
  # Return for further inspection
  list(
    inputs = tibble(
      home_team = home_team,
      away_team = away_team,
      parent_spread = parent_spread,
      parent_total = parent_total,
      target_cover = target_cover,
      target_over = target_over
    ),
    bets_summary = bets_summary,
    sample = inc_bal
  )
}

# Example: Run this to debug a specific prediction
# 
# result <- debug_prediction_pipeline(
#   DT = DT,
#   parent_spread = 3,      # +3 underdog
#   parent_total = 48,
#   target_cover = 0.50,
#   target_over = 0.50,
#   ss = ss, st = st, N = N,
#   margin_col = "game_home_margin_period",
#   use_spread_line = TRUE,
#   home_team = "Patriots",
#   away_team = "Bills"
# )
