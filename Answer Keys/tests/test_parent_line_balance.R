# Parent-line balancing tests (no database needed)
# Run from the Answer Keys directory:
#   Rscript tests/test_parent_line_balance.R            (C++ sampler)
#   Rscript tests/test_parent_line_balance.R pure-r     (pure-R sampler)

suppressMessages(library(data.table))
source("Tools.R")
if ("pure-r" %in% commandArgs(trailingOnly = TRUE)) .use_rcpp <- FALSE
cat("sampler:", if (.use_rcpp) "C++" else "pure R", "\n")

ok <- 0; bad <- 0
check <- function(cond, msg) {
  if (isTRUE(cond)) { cat("  PASS:", msg, "\n"); ok <<- ok + 1 }
  else              { cat("  FAIL:", msg, "\n"); bad <<- bad + 1 }
}

# --- 1. grade_vs_parent_line ---
dt <- data.table(home_spread = c(-1, -7, -3.5, -2), total_line = c(150, 130, 145, 140),
                 home_margin = c(3, 5, -2, NA), total_final_score = c(147, 140, 145, 150),
                 actual_cover = c(1L, 0L, 0L, 1L), actual_over = c(0L, 1L, 0L, 1L))
g <- grade_vs_parent_line(dt, -3.5, 145, use_spread_line = TRUE)
check(identical(g$actual_cover, c(0, 1, 0)), "spread cover graded vs parent -3.5")
check(identical(g$actual_over,  c(1, 0, 0.5)), "over graded vs parent 145, push = 0.5")
check(nrow(g) == 3 && identical(dt$actual_cover, c(1L, 0L, 0L, 1L)),
      "drops NA-score games and leaves the input untouched")
g <- grade_vs_parent_line(dt, 0.6, 145, use_spread_line = FALSE)
check(identical(g$actual_cover, c(1, 1, 0)), "moneyline mode: cover = home win")

g <- grade_vs_parent_line(data.table(game_home_margin_fg = c(4, -4), game_total_fg = c(150, 140)),
                          -3.5, 145, TRUE)
check(identical(g$actual_cover, c(1, 0)) && identical(g$actual_over, c(1, 0)),
      "finds game_home_margin_fg / game_total_fg")
g <- grade_vs_parent_line(data.table(game_home_margin_period_FG = c(1, -1), game_total_period_FG = c(9, 7)),
                          0.55, 8, FALSE)
check(identical(g$actual_cover, c(1, 0)) && identical(g$actual_over, c(1, 0)),
      "finds MLB game_home_margin_period_FG / game_total_period_FG")
err <- tryCatch({ grade_vs_parent_line(data.table(hm = 1, tot = 2), -3.5, 145, TRUE); NULL },
                error = function(e) e)
check(!is.null(err), "stops when the score columns are missing")

# --- 2. End to end: the sample hits the target on the PARENT line ---
# Integer final scores, no ties, half-point parent lines: no pushes. The
# stored actual_cover/actual_over are coin flips unrelated to the scores, so a
# regression to balancing on them lands near the natural rate, not the target.
set.seed(1)
M <- 8000
hist <- data.table(game_date = as.Date("2020-01-01") + sample(0:1500, M, TRUE),
                   home_spread = round(rnorm(M, -3.5, 4) * 2) / 2,
                   total_line  = round(rnorm(M, 145, 8) * 2) / 2)
hist[, home_margin := round(-home_spread + rnorm(M, 0, 11))]
hist[home_margin == 0, home_margin := 1]
hist[, total_final_score := round(total_line + rnorm(M, 0, 15))]
hist[, home_ml_odds := plogis(-home_spread / 6)]
hist[, actual_cover := rbinom(M, 1, 0.5)]
hist[, actual_over  := rbinom(M, 1, 0.5)]
ss <- sd(hist$home_spread); st <- sd(hist$total_line)

r <- run_answer_key_sample("t", -6.5, 150.5, 0.40, 0.60, hist, ss, st, N = 500, use_spread_line = TRUE)
s <- r$sample
cov <- mean(s$home_margin > 6.5); ov <- mean(s$total_final_score > 150.5)
check(abs(cov - 0.40) <= 0.004, sprintf("spread mode: covers -6.5 at %.3f (target 0.400)", cov))
check(abs(ov - 0.60) <= 0.004,  sprintf("spread mode: over 150.5 at %.3f (target 0.600)", ov))

# Control: balancing the stored columns misses both parent targets.
mm <- mean_match(hist, 500L, -6.5, 150.5, ss, st, 500, 0.005, TRUE)
own <- balance_sample(mm$dt, 500L, 0.40, 0.60, 1L)$dt[included == TRUE]
oc <- mean(own$home_margin > 6.5); oo <- mean(own$total_final_score > 150.5)
check(abs(oc - 0.40) > 0.04 && abs(oo - 0.60) > 0.04,
      sprintf("control: stored-column balancing gives cover %.3f / over %.3f", oc, oo))

r <- run_answer_key_sample("t", 0.70, 150.5, 0.70, 0.50, hist, ss, st, N = 500, use_spread_line = FALSE)
s <- r$sample
win <- mean(s$home_margin > 0)
check(abs(win - 0.70) <= 0.004, sprintf("moneyline mode: home wins %.3f (target 0.700)", win))

cat(sprintf("\n%d passed, %d failed\n", ok, bad))
if (bad > 0) quit(status = 1)
