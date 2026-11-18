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
check(identical(g$actual_cover, c(0L, 1L, 0L)), "spread cover graded vs parent -3.5")
check(identical(g$actual_over,  c(1L, 0L, 0L)), "over graded vs parent 145")
check(nrow(g) == 3 && identical(dt$actual_cover, c(1L, 0L, 0L, 1L)),
      "drops NA-score games and leaves the input untouched")
g <- grade_vs_parent_line(dt, 0.6, 145, use_spread_line = FALSE)
check(identical(g$actual_cover, c(1L, 1L, 0L)), "moneyline mode: cover = home win")

g <- grade_vs_parent_line(data.table(game_home_margin_fg = c(4, -4), game_total_fg = c(150, 140)),
                          -3.5, 145, TRUE)
check(identical(g$actual_cover, c(1L, 0L)) && identical(g$actual_over, c(1L, 0L)),
      "finds game_home_margin_fg / game_total_fg")
g <- grade_vs_parent_line(data.table(game_home_margin_period_FG = c(1, -1), game_total_period_FG = c(9, 7)),
                          0.55, 8, FALSE)
check(identical(g$actual_cover, c(1L, 0L)) && identical(g$actual_over, c(1L, 0L)),
      "finds MLB game_home_margin_period_FG / game_total_period_FG")
err <- tryCatch({ grade_vs_parent_line(data.table(hm = 1, tot = 2), -3.5, 145, TRUE); NULL },
                error = function(e) e)
check(!is.null(err), "stops when the score columns are missing")

# --- 2. End to end: the sample hits the target on the PARENT line ---
# Synthetic history where own-line and parent-line grading disagree: games on
# bigger own spreads cover them less often, so balancing own-line cover leaves
# the parent cover rate off target.
set.seed(1)
M <- 8000
hist <- data.table(game_date = as.Date("2020-01-01") + sample(0:1500, M, TRUE),
                   home_spread = round(rnorm(M, -3.5, 4) * 2) / 2,
                   total_line  = round(rnorm(M, 145, 8) * 2) / 2)
hist[, home_margin := round(-home_spread * 0.6 + rnorm(M, 0, 11)) + 0.5]
hist[, total_final_score := round(total_line * 0.9 + 14.5 + rnorm(M, 0, 15)) + 0.5]
hist[, home_ml_odds := plogis(-home_spread / 6)]
hist[, actual_cover := as.integer(home_margin > -home_spread)]
hist[, actual_over  := as.integer(total_final_score > total_line)]
ss <- sd(hist$home_spread); st <- sd(hist$total_line)

r <- run_answer_key_sample("t", -6.5, 150, 0.50, 0.50, hist, ss, st, N = 500, use_spread_line = TRUE)
s <- r$sample
cov <- mean(s$home_margin > 6.5); ov <- mean(s$total_final_score > 150)
check(abs(cov - 0.50) <= 0.004, sprintf("spread mode: covers -6.5 at %.3f (target 0.500)", cov))
check(abs(ov - 0.50) <= 0.004,  sprintf("spread mode: over 150 at %.3f (target 0.500)", ov))

# Control: balancing the stored own-line columns misses the parent cover, so
# the spread check above would catch a regression (the over grading is pinned
# by the unit checks in section 1).
mm <- mean_match(hist, 500L, -6.5, 150, ss, st, 500, 0.005, TRUE)
b <- balance_sample(mm$dt, 500L, 0.50, 0.50, 1L)
own <- b$dt[included == TRUE]
check(abs(mean(own$home_margin > 6.5) - 0.50) > 0.015,
      sprintf("control: own-line balancing covers -6.5 at %.3f (off target)", mean(own$home_margin > 6.5)))

r <- run_answer_key_sample("t", 0.70, 150, 0.70, 0.50, hist, ss, st, N = 500, use_spread_line = FALSE)
s <- r$sample
win <- mean(s$home_margin > 0)
check(abs(win - 0.70) <= 0.004, sprintf("moneyline mode: home wins %.3f (target 0.700)", win))

cat(sprintf("\n%d passed, %d failed\n", ok, bad))
if (bad > 0) quit(status = 1)
