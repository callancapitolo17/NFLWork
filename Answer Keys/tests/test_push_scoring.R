# =============================================================================
# Push scoring in balance_sample
# Run: Rscript tests/test_push_scoring.R  (from the Answer Keys directory)
#
# Pushes on whole-number lines are scored 0.5 (half cover / half over). Targets
# are no-vig probabilities that exclude pushes, so scoring a push as 0 biased
# balanced samples toward home/over. Uses synthetic data only (no DB needed).
# =============================================================================

suppressMessages({ library(Rcpp); library(data.table) })
sourceCpp("src/sampling.cpp")

score <- function(x, line) ifelse(x > line, 1, ifelse(x == line, 0.5, 0))
failed <- 0
check <- function(ok, msg) {
  cat(sprintf("  %s: %s\n", if (ok) "PASS" else "FAIL", msg))
  if (!ok) failed <<- failed + 1
}

# 1. Scoring helper semantics
check(identical(score(c(8, 7, 6), 7), c(1, 0.5, 0)), "win / push / loss -> 1 / 0.5 / 0")
check(is.na(score(NA_real_, 7)), "NA outcome stays NA")

# 2. Pool where the true no-push over rate is 50% but 15% of games push.
#    Balancing to a 50% target must leave non-push overs/unders even.
set.seed(1)
M <- 4000
pool <- sample(c(1, 0.5, 0), M, replace = TRUE, prob = c(0.425, 0.15, 0.425))
cov  <- sample(c(1, 0), M, replace = TRUE)
N <- 1000
res <- balance_sample_cpp(seq_len(M), cov, pool, N, 0.5, 0.5, 2L)
inc <- pool[res$included_indices]
np <- inc[inc != 0.5]
check(res$converged, "balance converges with pushes present")
check(abs(mean(np) - 0.5) < 0.01,
      sprintf("non-push over share %.3f ~= 0.5 target", mean(np)))

# Old scoring (push = 0) on the same pool skews toward over
old <- ifelse(pool == 1, 1, 0)
res_old <- balance_sample_cpp(seq_len(M), cov, old, N, 0.5, 0.5, 2L)
inc_old <- pool[res_old$included_indices]
np_old <- inc_old[inc_old != 0.5]
check(mean(np_old) > 0.52,
      sprintf("old push=0 scoring biased over (%.3f), confirming the bug", mean(np_old)))

# 3. Without pushes, results are identical to the old integer behaviour
pool01 <- sample(c(1, 0), M, replace = TRUE, prob = c(0.45, 0.55))
r1 <- balance_sample_cpp(seq_len(M), cov, pool01, N, 0.52, 0.5, 2L)
check(abs(r1$over_error) <= 2 && abs(r1$cover_error) <= 2, "0/1 data still balances within tol")

if (failed > 0) stop(failed, " test(s) failed")
cat("All push scoring tests passed\n")
