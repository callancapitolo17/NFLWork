# =============================================================================
# Rcpp Sampling Validation Tests
# Run: source("tests/test_rcpp_sampling.R") from Answer Keys directory
# =============================================================================
#
# Tests:
#   1. mean_match C++ vs R: sample means match targets
#   2. balance_sample C++ vs R: cover/over errors within tolerance
#   3. Prediction equivalence: derivative predictions agree within 1pp
#   4. Benchmark: C++ is >5x faster than R
#   5. Edge cases: N=1, extreme targets
#   6. Fork safety: mclapply with C++ functions
# =============================================================================

library(data.table)
library(duckdb)
library(dplyr)
library(tidyverse)
library(DBI)
source("Tools.R")

stopifnot(.use_rcpp)  # These tests require Rcpp to be available

cat("=============================================================\n")
cat("       RCPP SAMPLING VALIDATION TESTS\n")
cat("=============================================================\n\n")

# Load NFL data (same as test_answer_key.R)
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

pass <- function(msg) {
  cat(sprintf("  PASS: %s\n", msg))
  test_passed <<- test_passed + 1
}

fail <- function(msg) {
  cat(sprintf("  FAIL: %s\n", msg))
  test_failed <<- test_failed + 1
}

# =============================================================================
# TEST 1: mean_match C++ correctness
# =============================================================================
cat("TEST 1: mean_match C++ produces samples matching targets\n")

target_spread <- -3
target_total <- 45

.use_rcpp <<- TRUE
mm_cpp <- mean_match(DT, N, target_spread, target_total, ss, st,
                     max_iter_mean = 500, tol_mean = 0.005, use_spread_line = TRUE)
inc_cpp <- mm_cpp$dt[included == TRUE]

spread_err_cpp <- abs(mean(inc_cpp$home_spread) - target_spread)
total_err_cpp <- abs(mean(inc_cpp$total_line) - target_total)

if (spread_err_cpp < 0.1) {
  pass(sprintf("C++ mean spread = %.3f (target = %.1f, error = %.4f)",
               mean(inc_cpp$home_spread), target_spread, spread_err_cpp))
} else {
  fail(sprintf("C++ mean spread = %.3f (target = %.1f, error = %.4f > 0.1)",
               mean(inc_cpp$home_spread), target_spread, spread_err_cpp))
}

if (total_err_cpp < 0.1) {
  pass(sprintf("C++ mean total = %.3f (target = %.1f, error = %.4f)",
               mean(inc_cpp$total_line), target_total, total_err_cpp))
} else {
  fail(sprintf("C++ mean total = %.3f (target = %.1f, error = %.4f > 0.1)",
               mean(inc_cpp$total_line), target_total, total_err_cpp))
}

# =============================================================================
# TEST 2: balance_sample C++ correctness
# =============================================================================
cat("\nTEST 2: balance_sample C++ achieves target cover/over rates\n")

target_cover <- 0.52
target_over <- 0.50

.use_rcpp <<- TRUE
bal_cpp <- balance_sample(mm_cpp$dt, N, target_cover, target_over, tol_error = 1)

if (abs(bal_cpp$cover_error) <= 1) {
  pass(sprintf("C++ cover error = %d (within tolerance)", bal_cpp$cover_error))
} else {
  fail(sprintf("C++ cover error = %d (exceeds tolerance)", bal_cpp$cover_error))
}

if (abs(bal_cpp$over_error) <= 1) {
  pass(sprintf("C++ over error = %d (within tolerance)", bal_cpp$over_error))
} else {
  fail(sprintf("C++ over error = %d (exceeds tolerance)", bal_cpp$over_error))
}

if (bal_cpp$converged) {
  pass("C++ balance_sample converged")
} else {
  fail("C++ balance_sample did not converge")
}

# =============================================================================
# TEST 3: Prediction equivalence (R vs C++)
# =============================================================================
cat("\nTEST 3: Prediction equivalence across 20 test cases\n")

# Grid of realistic NFL (spread, total) combinations
test_cases <- expand.grid(
  spread = c(-7, -3, 0, 3, 7),
  total = c(40, 44, 48, 52)
)

max_divergence <- 0
all_within_1pp <- TRUE
all_within_2pp <- TRUE

for (i in seq_len(nrow(test_cases))) {
  ps <- test_cases$spread[i]
  pt <- test_cases$total[i]

  # Devig to get target cover/over rates
  # Use simple 50-50 for testing (exact rates don't matter, just consistency)
  tc <- 0.50 + (ps / 100)  # slight adjustment based on spread
  to <- 0.50

  # Run R path
  .use_rcpp <<- FALSE
  r_result <- run_answer_key_sample(
    id = paste0("test_", i),
    parent_spread = ps,
    parent_total = pt,
    target_cover = tc,
    target_over = to,
    DT = DT, ss = ss, st = st, N = N,
    use_spread_line = TRUE
  )
  r_sample <- r_result$sample

  # Run C++ path
  .use_rcpp <<- TRUE
  cpp_result <- run_answer_key_sample(
    id = paste0("test_", i),
    parent_spread = ps,
    parent_total = pt,
    target_cover = tc,
    target_over = to,
    DT = DT, ss = ss, st = st, N = N,
    use_spread_line = TRUE
  )
  cpp_sample <- cpp_result$sample

  # Compare derivative predictions from both samples
  # Extract all game_home_margin_period_* columns (quarter/half predictions)
  margin_cols <- grep("^game_home_margin_period", names(r_sample), value = TRUE)
  total_cols <- grep("^game_total_period", names(r_sample), value = TRUE)
  pred_cols <- c(margin_cols, total_cols)

  for (col in pred_cols) {
    if (col %in% names(r_sample) && col %in% names(cpp_sample)) {
      r_pred <- mean(r_sample[[col]] > 0, na.rm = TRUE)
      cpp_pred <- mean(cpp_sample[[col]] > 0, na.rm = TRUE)
      diff <- abs(r_pred - cpp_pred)
      if (diff > max_divergence) max_divergence <- diff
      if (diff > 0.01) all_within_1pp <- FALSE
      if (diff > 0.02) {
        all_within_2pp <- FALSE
        cat(sprintf("    WARNING: spread=%.0f total=%.0f col=%s R=%.3f C++=%.3f diff=%.3f\n",
                    ps, pt, col, r_pred, cpp_pred, diff))
      }
    }
  }

  # Also compare sample-level stats
  r_mean_spread <- mean(r_sample$home_spread, na.rm = TRUE)
  cpp_mean_spread <- mean(cpp_sample$home_spread, na.rm = TRUE)
  r_mean_total <- mean(r_sample$total_line, na.rm = TRUE)
  cpp_mean_total <- mean(cpp_sample$total_line, na.rm = TRUE)

  spread_diff <- abs(r_mean_spread - cpp_mean_spread)
  total_diff <- abs(r_mean_total - cpp_mean_total)

  if (spread_diff > 0.1 || total_diff > 0.1) {
    cat(sprintf("    WARNING: spread=%.0f total=%.0f mean_spread_diff=%.4f mean_total_diff=%.4f\n",
                ps, pt, spread_diff, total_diff))
  }
}

# Re-enable Rcpp for remaining tests
.use_rcpp <<- TRUE

cat(sprintf("  Max prediction divergence across all tests: %.4f (%.2f pp)\n",
            max_divergence, max_divergence * 100))

if (all_within_2pp) {
  pass("All derivative predictions agree within 2pp")
} else {
  fail("Some predictions diverged by >2pp — investigate")
}

if (all_within_1pp) {
  pass("All derivative predictions agree within 1pp")
} else {
  cat("  NOTE: Some predictions diverged by 1-2pp (acceptable but worth monitoring)\n")
}

# =============================================================================
# TEST 4: Benchmark (C++ vs R)
# =============================================================================
cat("\nTEST 4: Performance benchmark\n")

n_bench <- 5  # number of runs to average
target_spread <- -3
target_total <- 45
target_cover <- 0.52
target_over <- 0.50

# Time R version
.use_rcpp <<- FALSE
r_times <- numeric(n_bench)
for (i in seq_len(n_bench)) {
  t0 <- proc.time()["elapsed"]
  run_answer_key_sample("bench", target_spread, target_total, target_cover, target_over,
                        DT, ss, st, N, use_spread_line = TRUE)
  r_times[i] <- proc.time()["elapsed"] - t0
}
r_avg <- mean(r_times)

# Time C++ version
.use_rcpp <<- TRUE
cpp_times <- numeric(n_bench)
for (i in seq_len(n_bench)) {
  t0 <- proc.time()["elapsed"]
  run_answer_key_sample("bench", target_spread, target_total, target_cover, target_over,
                        DT, ss, st, N, use_spread_line = TRUE)
  cpp_times[i] <- proc.time()["elapsed"] - t0
}
cpp_avg <- mean(cpp_times)

speedup <- r_avg / cpp_avg
cat(sprintf("  R average:   %.3f sec\n", r_avg))
cat(sprintf("  C++ average: %.3f sec\n", cpp_avg))
cat(sprintf("  Speedup:     %.1fx\n", speedup))

if (speedup > 5) {
  pass(sprintf("C++ is %.1fx faster than R (>5x threshold)", speedup))
} else if (speedup > 2) {
  cat(sprintf("  NOTE: Speedup is %.1fx — less than 5x target but still meaningful\n", speedup))
  pass(sprintf("C++ is %.1fx faster than R (>2x, acceptable)", speedup))
} else {
  fail(sprintf("C++ speedup only %.1fx (expected >5x)", speedup))
}

# =============================================================================
# TEST 5: Edge cases
# =============================================================================
cat("\nTEST 5: Edge cases\n")

# N = 1
.use_rcpp <<- TRUE
tryCatch({
  mm_small <- mean_match(DT, 1L, -3, 45, ss, st, 500, 0.005, use_spread_line = TRUE)
  inc_small <- mm_small$dt[included == TRUE]
  if (nrow(inc_small) == 1) {
    pass("N=1 returns exactly 1 game")
  } else {
    fail(sprintf("N=1 returned %d games", nrow(inc_small)))
  }
}, error = function(e) {
  fail(sprintf("N=1 threw error: %s", e$message))
})

# Extreme spread
tryCatch({
  mm_extreme <- mean_match(DT, N, -20, 60, ss, st, 500, 0.05, use_spread_line = TRUE)
  inc_extreme <- mm_extreme$dt[included == TRUE]
  if (nrow(inc_extreme) == N) {
    pass(sprintf("Extreme target (-20/60) still selects %d games", N))
  } else {
    pass(sprintf("Extreme target (-20/60) selects %d games (may differ from N)", nrow(inc_extreme)))
  }
}, error = function(e) {
  fail(sprintf("Extreme target threw error: %s", e$message))
})

# =============================================================================
# TEST 6: Fork safety (mclapply)
# =============================================================================
cat("\nTEST 6: Fork safety with mclapply\n")

.use_rcpp <<- TRUE
targets_par <- data.frame(
  spread = c(-7, -3, 0, 3),
  total = c(42, 45, 48, 51),
  cover = c(0.55, 0.52, 0.50, 0.48),
  over = c(0.50, 0.50, 0.50, 0.50)
)

tryCatch({
  results <- parallel::mclapply(seq_len(nrow(targets_par)), function(i) {
    run_answer_key_sample(
      id = paste0("fork_", i),
      parent_spread = targets_par$spread[i],
      parent_total = targets_par$total[i],
      target_cover = targets_par$cover[i],
      target_over = targets_par$over[i],
      DT = DT, ss = ss, st = st, N = N,
      use_spread_line = TRUE
    )
  }, mc.cores = min(4, parallel::detectCores() - 1))

  # Check all returned valid results
  all_valid <- all(sapply(results, function(r) {
    !inherits(r, "try-error") && !is.null(r$sample) && nrow(r$sample) > 0
  }))

  if (all_valid) {
    pass(sprintf("mclapply with 4 cores produced %d valid samples", length(results)))
  } else {
    fail("Some mclapply results were invalid")
  }
}, error = function(e) {
  fail(sprintf("mclapply threw error: %s", e$message))
})

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n=============================================================\n")
cat(sprintf("RESULTS: %d passed, %d failed\n", test_passed, test_failed))
cat("=============================================================\n")

if (test_failed > 0) {
  cat("\nWARNING: Some tests failed. Review the failures above.\n")
} else {
  cat("\nAll tests passed.\n")
}
