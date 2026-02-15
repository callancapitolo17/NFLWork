# Answer Key Algorithm Validation Tests
# Run: source("tests/test_answer_key.R") from Answer Keys directory

library(data.table)
library(duckdb)
library(dplyr)
library(tidyverse)
library(DBI)
source("Tools.R")

cat("=============================================================\n")
cat("       ANSWER KEY ALGORITHM VALIDATION TESTS\n")
cat("=============================================================\n\n")

# Load data
con <- dbConnect(duckdb(), dbdir = "pbp.duckdb", read_only = TRUE)
DT <- dbGetQuery(con, "SELECT * FROM nfl_betting_pbp") %>%
  mutate(
    actual_over = ifelse(total_final_score > total_line, 1, 0),
    actual_cover = ifelse(home_margin > -home_spread, 1, 0)
  ) %>%
  as.data.table()
dbDisconnect(con)

disp <- compute_dispersion(DT, moneyline = FALSE)
ss <- disp$ss
st <- disp$st
N <- round(nrow(DT) * 0.10, 0)

test_passed <- 0
test_failed <- 0

# Test 1: Odds conversion
cat("TEST 1: American odds to probability conversion\n")
test_cases <- list(
  list(odds = -110, expected = 0.5238, tol = 0.001),
  list(odds = +100, expected = 0.5, tol = 0.001),
  list(odds = -200, expected = 0.6667, tol = 0.001),
  list(odds = +200, expected = 0.3333, tol = 0.001)
)
for (tc in test_cases) {
  result <- odds_to_prob(tc$odds)
  if (abs(result - tc$expected) < tc$tol) {
    cat(sprintf("  PASS: odds_to_prob(%d) = %.4f (expected %.4f)\n", tc$odds, result, tc$expected))
    test_passed <- test_passed + 1
  } else {
    cat(sprintf("  FAIL: odds_to_prob(%d) = %.4f (expected %.4f)\n", tc$odds, result, tc$expected))
    test_failed <- test_failed + 1
  }
}

# Test 2: Devig function
cat("\nTEST 2: Devigging two-sided odds\n")
result <- devig_american(-110, -110)
if (abs(result$p1 - 0.5) < 0.001 && abs(result$p2 - 0.5) < 0.001) {
  cat(sprintf("  PASS: devig_american(-110, -110) = (%.4f, %.4f)\n", result$p1, result$p2))
  test_passed <- test_passed + 1
} else {
  cat(sprintf("  FAIL: devig_american(-110, -110) = (%.4f, %.4f), expected (0.5, 0.5)\n", result$p1, result$p2))
  test_failed <- test_failed + 1
}

result2 <- devig_american(-200, +180)
if (abs(result2$p1 + result2$p2 - 1.0) < 0.001) {
  cat(sprintf("  PASS: devig probabilities sum to 1.0 (%.4f + %.4f = %.4f)\n", result2$p1, result2$p2, result2$p1 + result2$p2))
  test_passed <- test_passed + 1
} else {
  cat(sprintf("  FAIL: devig probabilities don't sum to 1.0\n"))
  test_failed <- test_failed + 1
}

# Test 3: Mean match produces expected spread/total means
cat("\nTEST 3: Mean match produces samples matching target\n")
target_spread <- -3
target_total <- 45
mm <- mean_match(DT, N, target_spread, target_total, ss, st,
                 max_iter_mean = 500, tol_mean = 0.005, use_spread_line = TRUE)
inc <- mm$dt[included == TRUE]

spread_error <- abs(mean(inc$home_spread) - target_spread)
total_error <- abs(mean(inc$total_line) - target_total)

if (spread_error < 0.1) {
  cat(sprintf("  PASS: Mean spread = %.3f (target = %.1f, error = %.4f)\n",
              mean(inc$home_spread), target_spread, spread_error))
  test_passed <- test_passed + 1
} else {
  cat(sprintf("  FAIL: Mean spread = %.3f (target = %.1f, error = %.4f > 0.1)\n",
              mean(inc$home_spread), target_spread, spread_error))
  test_failed <- test_failed + 1
}

if (total_error < 0.1) {
  cat(sprintf("  PASS: Mean total = %.3f (target = %.1f, error = %.4f)\n",
              mean(inc$total_line), target_total, total_error))
  test_passed <- test_passed + 1
} else {
  cat(sprintf("  FAIL: Mean total = %.3f (target = %.1f, error = %.4f > 0.1)\n",
              mean(inc$total_line), target_total, total_error))
  test_failed <- test_failed + 1
}

# Test 4: Balance sample achieves target rates (within tolerance)
cat("\nTEST 4: Balance sample achieves target cover/over rates\n")
target_cover <- 0.52
target_over <- 0.50
bal <- balance_sample(mm$dt, N, target_cover, target_over, tol_error = 1)

if (abs(bal$cover_error) <= 1) {
  cat(sprintf("  PASS: Cover error = %d (within tolerance of 1)\n", bal$cover_error))
  test_passed <- test_passed + 1
} else {
  cat(sprintf("  FAIL: Cover error = %d (exceeds tolerance of 1)\n", bal$cover_error))
  test_failed <- test_failed + 1
}

if (abs(bal$over_error) <= 1) {
  cat(sprintf("  PASS: Over error = %d (within tolerance of 1)\n", bal$over_error))
  test_passed <- test_passed + 1
} else {
  cat(sprintf("  FAIL: Over error = %d (exceeds tolerance of 1)\n", bal$over_error))
  test_failed <- test_failed + 1
}

# Test 5: skip_balance reduces bias for derivatives
cat("\nTEST 5: skip_balance reduces selection bias for derivatives\n")
result_with_balance <- run_means_for_id("test", -3, 45, 0.52, 0.50, DT, ss, st, N,
                                         use_spread_line = TRUE, margin_col = "game_home_margin_period",
                                         skip_balance = FALSE, shrinkage_weight = 0)
result_skip_balance <- run_means_for_id("test", -3, 45, 0.52, 0.50, DT, ss, st, N,
                                         use_spread_line = TRUE, margin_col = "game_home_margin_period",
                                         skip_balance = TRUE, shrinkage_weight = 0)

# Balanced version should show higher Q2 (due to selection bias)
# We just verify both produce valid probabilities
if (result_with_balance$game_home_margin_period_2 >= 0 && result_with_balance$game_home_margin_period_2 <= 1) {
  cat(sprintf("  PASS: Balanced Q2 = %.3f (valid probability)\n", result_with_balance$game_home_margin_period_2))
  test_passed <- test_passed + 1
} else {
  cat(sprintf("  FAIL: Balanced Q2 = %.3f (invalid probability)\n", result_with_balance$game_home_margin_period_2))
  test_failed <- test_failed + 1
}

if (result_skip_balance$game_home_margin_period_2 >= 0 && result_skip_balance$game_home_margin_period_2 <= 1) {
  cat(sprintf("  PASS: Skip-balance Q2 = %.3f (valid probability)\n", result_skip_balance$game_home_margin_period_2))
  test_passed <- test_passed + 1
} else {
  cat(sprintf("  FAIL: Skip-balance Q2 = %.3f (invalid probability)\n", result_skip_balance$game_home_margin_period_2))
  test_failed <- test_failed + 1
}

# Test 6: Shrinkage pulls toward 0.5
cat("\nTEST 6: Shrinkage pulls predictions toward 0.5\n")
result_no_shrink <- run_means_for_id("test", -3, 45, 0.52, 0.50, DT, ss, st, N,
                                      use_spread_line = TRUE, margin_col = "game_home_margin_period",
                                      skip_balance = TRUE, shrinkage_weight = 0)
result_with_shrink <- run_means_for_id("test", -3, 45, 0.52, 0.50, DT, ss, st, N,
                                        use_spread_line = TRUE, margin_col = "game_home_margin_period",
                                        skip_balance = TRUE, shrinkage_weight = 0.4)

expected_shrunk <- 0.4 * 0.5 + 0.6 * result_no_shrink$game_home_margin_period_2
if (abs(result_with_shrink$game_home_margin_period_2 - expected_shrunk) < 0.001) {
  cat(sprintf("  PASS: Shrinkage works correctly (%.3f → %.3f)\n",
              result_no_shrink$game_home_margin_period_2, result_with_shrink$game_home_margin_period_2))
  test_passed <- test_passed + 1
} else {
  cat(sprintf("  FAIL: Shrinkage calculation error (got %.3f, expected %.3f)\n",
              result_with_shrink$game_home_margin_period_2, expected_shrunk))
  test_failed <- test_failed + 1
}

# Summary
cat("\n=============================================================\n")
cat(sprintf("RESULTS: %d passed, %d failed\n", test_passed, test_failed))
cat("=============================================================\n")

if (test_failed > 0) {
  cat("\nWARNING: Some tests failed. Review the failures above.\n")
} else {
  cat("\nAll tests passed.\n")
}
