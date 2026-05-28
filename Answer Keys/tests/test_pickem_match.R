# Tests for derive_pickem_american() — pick'em DNB math.
# Run from worktree root: Rscript "Answer Keys/tests/test_pickem_match.R"
suppressWarnings(suppressMessages({
  source("Answer Keys/Tools.R")          # devig_american, devig_american_3way, prob_to_american
  source("Answer Keys/MLB Answer Key/odds_screen.R")
}))

# 2-way: a winner market already excludes ties.
# DK 1st-3-innings COL +140 / LA -180. Probit 2-way devig -> matched probs.
r <- derive_pickem_american(home_raw = -180, away_raw = 140, tie_raw = NA)
stopifnot(is.list(r))
stopifnot(all(c("home_raw_dnb","away_raw_dnb",
                "home_fair_dnb","away_fair_dnb") %in% names(r)))
stopifnot(r$home_raw_dnb == -180)   # raw unchanged for 2-way
stopifnot(r$away_raw_dnb == 140)
# Fair: devigged probs sum to 1, so the two Americans round-trip to probs summing to 1.
p_home <- if (r$home_fair_dnb < 0) -r$home_fair_dnb / (-r$home_fair_dnb + 100) else 100 / (r$home_fair_dnb + 100)
p_away <- if (r$away_fair_dnb < 0) -r$away_fair_dnb / (-r$away_fair_dnb + 100) else 100 / (r$away_fair_dnb + 100)
stopifnot(abs((p_home + p_away) - 1) < 1e-6)
cat("OK 2-way path\n")

# 3-way: synthetic Home -150 / Tie +400 / Away +180.
r3 <- derive_pickem_american(home_raw = -150, away_raw = 180, tie_raw = 400)
stopifnot(!is.na(r3$home_fair_dnb), !is.na(r3$away_fair_dnb))
p_h_fair <- if (r3$home_fair_dnb < 0) -r3$home_fair_dnb / (-r3$home_fair_dnb + 100) else 100 / (r3$home_fair_dnb + 100)
p_a_fair <- if (r3$away_fair_dnb < 0) -r3$away_fair_dnb / (-r3$away_fair_dnb + 100) else 100 / (r3$away_fair_dnb + 100)
stopifnot(abs((p_h_fair + p_a_fair) - 1) < 1e-6)
# Raw DNB: sum to 1 and home favored (home was -150).
p_h_raw <- if (r3$home_raw_dnb < 0) -r3$home_raw_dnb / (-r3$home_raw_dnb + 100) else 100 / (r3$home_raw_dnb + 100)
p_a_raw <- if (r3$away_raw_dnb < 0) -r3$away_raw_dnb / (-r3$away_raw_dnb + 100) else 100 / (r3$away_raw_dnb + 100)
stopifnot(abs((p_h_raw + p_a_raw) - 1) < 1e-6)
stopifnot(p_h_raw > p_a_raw)
cat("OK 3-way path\n")

# Degenerate inputs.
stopifnot(is.na(derive_pickem_american(NA, 140, NA)$home_fair_dnb))
cat("OK degenerate inputs\n")
