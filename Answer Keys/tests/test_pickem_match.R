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

# 3-way canonicalization: a wide row with odds_tie (line == NA) -> three
# h2h_3way rows (home, away, Tie) with the right period.
library(tibble)
raw3 <- tibble(
  market = "h2h_3way_1st_5_innings",
  home_team = "LA Dodgers", away_team = "COL Rockies",
  line = NA_real_,
  odds_home = -150, odds_away = 180, odds_tie = 400,
  odds_over = NA_real_, odds_under = NA_real_,
  fetch_time = Sys.time()
)
lookup3 <- tibble(id = "g1", home_team = "LA Dodgers", away_team = "COL Rockies")
canon <- scraper_to_canonical(raw3, lookup3, book_name = "test")
stopifnot(nrow(canon) == 3)
stopifnot(all(canon$market == "h2h_3way"))
stopifnot(all(canon$period == "F5"))
stopifnot(setequal(canon$side, c("LA Dodgers", "COL Rockies", "Tie")))
cat("OK 3-way canonicalization\n")

# === Task 6: matcher line-0 pick'em branch ===
# Bet: PHI 0 spread, period F3. Book DK: 2-way h2h_1st_3_innings winner.
bets <- tibble(
  bet_row_id = "b1", id = "g1", home_team = "SD Padres", away_team = "PHI Phillies",
  market = "spreads_1st_3_innings", market_type = "spreads", period = "F3",
  line = 0, bet_on = "PHI Phillies", opposite_side = "SD Padres",
  pt_start_time = as.POSIXct(NA, tz = "UTC")
)
dk_canonical <- normalize_book_odds_frame(tibble(
  game_id = "g1", market_name = "h2h_1st_3_innings",
  bet_on = c("SD Padres", "PHI Phillies"),
  line = c(NA_real_, NA_real_), american_odds = c(-180L, 140L),
  fetch_time = Sys.time()
))
out <- expand_bets_to_book_prices(bets, list(dk = dk_canonical))
stopifnot(nrow(out) == 2)
stopifnot(all(out$line == 0))
stopifnot(all(out$is_exact_line))
stopifnot(all(!is.na(out$derived_fair_odds)))
pick <- out[out$side == "pick", ]; opp <- out[out$side == "opposite", ]
stopifnot(pick$american_odds == 140)   # PHI (away) raw
stopifnot(opp$american_odds == -180)   # SD (home) raw
cat("OK matcher 2-way pick'em\n")

# Same bet, book only has a 3-way h2h_3way winner.
wz_canonical <- normalize_book_odds_frame(tibble(
  game_id = "g1", market_name = "h2h_3way_1st_3_innings",
  bet_on = c("SD Padres", "PHI Phillies", "Tie"),
  line = c(NA_real_, NA_real_, NA_real_), american_odds = c(-150L, 180L, 400L),
  fetch_time = Sys.time()
))
out2 <- expand_bets_to_book_prices(bets, list(wz = wz_canonical))
stopifnot(nrow(out2) == 2)
stopifnot(all(out2$is_exact_line))
stopifnot(all(!is.na(out2$derived_fair_odds)))
cat("OK matcher 3-way pick'em\n")

# A NON-zero spread bet must NOT enter the pick'em branch (regression).
bets_alt <- tibble(
  bet_row_id = "b2", id = "g1", home_team = "SD Padres", away_team = "PHI Phillies",
  market = "alternate_spreads_fg", market_type = "alternate_spreads", period = "FG",
  line = -2.5, bet_on = "PHI Phillies", opposite_side = "SD Padres",
  pt_start_time = as.POSIXct(NA, tz = "UTC")
)
fake_book <- normalize_book_odds_frame(tibble(
  game_id = "g1", market_name = "alternate_spreads_fg",
  bet_on = "PHI Phillies", line = -2.5, american_odds = 240L, fetch_time = Sys.time()
))
out3 <- expand_bets_to_book_prices(bets_alt, list(test = fake_book))
stopifnot(nrow(out3) >= 1)
stopifnot(all(is.na(out3$derived_fair_odds)))
cat("OK non-zero spread unaffected\n")
