# Parent-line balancing tests (no database needed)
# Run: Rscript tests/test_parent_line_balance.R  from the Answer Keys directory
#
# balance_sample must grade each historical game against the PARENT line, not
# the game's own closing line (Feustel: the sample covers the parent at the
# no-vig rate, so its median margin lands on the parent).

suppressMessages(library(data.table))
source("Tools.R")

ok <- 0; bad <- 0
check <- function(cond, msg) {
  if (isTRUE(cond)) { cat("  PASS:", msg, "\n"); ok <<- ok + 1 }
  else              { cat("  FAIL:", msg, "\n"); bad <<- bad + 1 }
}

# --- 1. Grading against the parent line ---
dt <- data.table(home_spread = c(-1, -7, -3.5), total_line = c(150, 130, 145),
                 home_margin = c(3, 5, -2), total_final_score = c(147, 140, 145),
                 actual_cover = c(1L, 0L, 0L), actual_over = c(0L, 1L, 0L))
set_parent_line_outcomes(dt, parent_spread = -3.5, parent_total = 145, use_spread_line = TRUE)
check(identical(dt$actual_cover, c(0L, 1L, 0L)), "spread cover graded vs parent -3.5")
check(identical(dt$actual_over,  c(1L, 0L, 0L)), "over graded vs parent 145")

set_parent_line_outcomes(dt, parent_spread = 0.6, parent_total = 145, use_spread_line = FALSE)
check(identical(dt$actual_cover, c(1L, 1L, 0L)), "moneyline mode: cover = home win")

# Alternate column names (CBB backtests, MLB)
dt2 <- data.table(game_home_margin_fg = c(4, -4), game_total_fg = c(150, 140))
set_parent_line_outcomes(dt2, -3.5, 145, TRUE)
check(identical(dt2$actual_cover, c(1L, 0L)) && identical(dt2$actual_over, c(1L, 0L)),
      "auto-detects game_home_margin_fg / game_total_fg")

# Missing score columns: keep stored outcomes, warn
dt3 <- data.table(actual_cover = c(1L, 0L), actual_over = c(0L, 1L))
w <- tryCatch({ set_parent_line_outcomes(dt3, -3.5, 145, TRUE); NULL }, warning = function(w) w)
check(!is.null(w), "warns when no score columns exist")

# --- 2. End to end: sample covers the parent at the target rate ---
# Synthetic history where each game's own line differs from the parent, so
# own-line balancing and parent-line balancing disagree.
set.seed(1)
M <- 6000
hist <- data.table(game_date = as.Date("2020-01-01") + sample(0:1500, M, TRUE),
                   home_spread = round(rnorm(M, -3.5, 4) * 2) / 2,
                   total_line  = round(rnorm(M, 145, 8) * 2) / 2)
hist[, home_margin := round(-home_spread + rnorm(M, 0, 11))]
hist[, total_final_score := round(total_line + rnorm(M, 0, 15))]
hist[, home_ml_odds := 0.5]
hist[, actual_cover := as.integer(home_margin > -home_spread)]
hist[, actual_over  := as.integer(total_final_score > total_line)]

r <- run_answer_key_sample("t", -3.5, 145, 0.52, 0.5, hist,
                           ss = sd(hist$home_spread), st = sd(hist$total_line),
                           N = 500, use_spread_line = TRUE)
s <- r$sample
cov <- mean(s$home_margin > 3.5); ov <- mean(s$total_final_score > 145)
check(abs(cov - 0.52) <= 0.004, sprintf("sample covers -3.5 at %.3f (target 0.520)", cov))
check(abs(ov - 0.50) <= 0.004,  sprintf("sample goes over 145 at %.3f (target 0.500)", ov))

cat(sprintf("\n%d passed, %d failed\n", ok, bad))
if (bad > 0) quit(status = 1)
