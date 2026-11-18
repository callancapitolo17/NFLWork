# =============================================================================
# Push scoring in balance_sample
# Run: Rscript tests/test_push_scoring.R  (from the Answer Keys directory)
#
# A push on a whole-number line is stored as 0.5 in actual_cover / actual_over.
# Balance targets are no-vig probabilities that exclude pushes, so
# balance_sample weights a push by the target: the balanced sample's
# non-push cover/over share must land on the target. Synthetic data only.
# =============================================================================

suppressMessages(library(data.table))
source("Tools.R")
stopifnot(.use_rcpp)

failed <- 0
check <- function(ok, msg) {
  cat(sprintf("  %s: %s\n", if (ok) "PASS" else "FAIL", msg))
  if (!ok) failed <<- failed + 1
}

M <- 5000
N <- 1000
TOL <- 2L
PUSH_SHARE <- 0.15

make_pool <- function(seed) {
  set.seed(seed)
  data.table(
    actual_cover = sample(c(1, 0.5, 0), M, replace = TRUE,
                          prob = c(0.5 - PUSH_SHARE / 2, PUSH_SHARE, 0.5 - PUSH_SHARE / 2)),
    actual_over  = sample(c(1, 0.5, 0), M, replace = TRUE,
                          prob = c(0.5 - PUSH_SHARE / 2, PUSH_SHARE, 0.5 - PUSH_SHARE / 2)),
    included     = c(rep(TRUE, N), rep(FALSE, M - N))
  )
}

non_push_share <- function(x) mean(x[x != 0.5])

run_both_paths <- function(pool, target_cover, target_over) {
  .use_rcpp <<- TRUE
  cpp <- balance_sample(pool, N, target_cover, target_over, TOL)
  .use_rcpp <<- FALSE
  r <- balance_sample(pool, N, target_cover, target_over, TOL)
  .use_rcpp <<- TRUE
  list(cpp = cpp, r = r)
}

# 1. Non-push share hits the target off 0.5, on both paths, and both agree
for (targets in list(c(0.5, 0.5), c(0.60, 0.40), c(0.45, 0.58))) {
  pool <- make_pool(1)
  res <- run_both_paths(pool, targets[1], targets[2])
  label <- sprintf("targets %.2f / %.2f", targets[1], targets[2])
  for (path in c("cpp", "r")) {
    inc <- res[[path]]$dt[included == TRUE]
    check(res[[path]]$converged, sprintf("%s %s converges", path, label))
    cover_gap <- abs(non_push_share(inc$actual_cover) - targets[1])
    over_gap  <- abs(non_push_share(inc$actual_over)  - targets[2])
    check(cover_gap < 0.004 && over_gap < 0.004,
          sprintf("%s %s non-push shares %.4f / %.4f", path, label,
                  non_push_share(inc$actual_cover), non_push_share(inc$actual_over)))
  }
  check(identical(res$cpp$dt$included, res$r$dt$included),
        sprintf("C++ and R pick the same sample (%s)", label))
  check(identical(names(res$r$dt), names(pool)) &&
          identical(res$r$dt$actual_cover, pool$actual_cover),
        sprintf("R fallback returns the stored 0 / 0.5 / 1 columns unchanged (%s)", label))
}

# 1b. C++ and R agree across many small random pools
mismatches <- 0
for (seed in 1:200) {
  set.seed(seed)
  small <- data.table(
    actual_cover = sample(c(1, 0.5, 0), 400, replace = TRUE, prob = c(0.4, 0.2, 0.4)),
    actual_over  = sample(c(1, 0.5, 0), 400, replace = TRUE, prob = c(0.4, 0.2, 0.4)),
    included     = c(rep(TRUE, 100), rep(FALSE, 300))
  )
  .use_rcpp <<- TRUE
  a <- balance_sample(small, 100L, 0.55, 0.45, 1L)
  .use_rcpp <<- FALSE
  b <- balance_sample(small, 100L, 0.55, 0.45, 1L)
  .use_rcpp <<- TRUE
  if (!identical(a$dt$included, b$dt$included)) mismatches <- mismatches + 1
}
check(mismatches == 0, sprintf("C++ and R pick the same sample on 200 random pools (%d differ)", mismatches))

# 2. The old push = 0 scoring biases the same pool toward cover/over
pool <- make_pool(1)
old <- copy(pool)[, `:=`(actual_cover = as.numeric(actual_cover == 1),
                         actual_over  = as.numeric(actual_over == 1))]
inc_old <- pool[balance_sample(old, N, 0.5, 0.5, TOL)$dt$included]
check(non_push_share(inc_old$actual_over) > 0.52,
      sprintf("old push = 0 scoring biased over (%.3f), confirming the bug",
              non_push_share(inc_old$actual_over)))

# 3. Data with no pushes balances as before (pushes are the only change)
pool01 <- make_pool(2)[, `:=`(actual_cover = round(actual_cover), actual_over = round(actual_over))]
res01 <- balance_sample(pool01, N, 0.52, 0.5, TOL)
check(res01$converged && abs(res01$cover_error) <= TOL && abs(res01$over_error) <= TOL,
      "0/1 data still balances within tol")

if (failed > 0) stop(failed, " test(s) failed")
cat("All push scoring tests passed\n")
