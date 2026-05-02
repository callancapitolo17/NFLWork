# Answer Keys/tests/test_compare_alts_to_samples.R
# Fixture-driven tests for compare_alts_to_samples derivative-h2h + odd/even branches.
# Run from "Answer Keys/" directory: Rscript -e 'testthat::test_file("tests/test_compare_alts_to_samples.R")'

library(testthat)
library(dplyr)
library(tibble)
library(data.table)
source("../Tools.R")

# Helper: build a synthetic samples list with one game.
# margins_f3 / margins_f7 / totals_fg are vectors of 1000 sim outcomes.
make_synthetic_samples <- function(
    game_id = "test-game-1",
    margins_f3 = NULL,
    margins_f7 = NULL,
    totals_fg = NULL
) {
  n <- 1000L
  if (is.null(margins_f3)) margins_f3 <- rep(0L, n)
  if (is.null(margins_f7)) margins_f7 <- rep(0L, n)
  if (is.null(totals_fg)) totals_fg <- rep(8L, n)
  stopifnot(length(margins_f3) == n, length(margins_f7) == n, length(totals_fg) == n)
  samp <- data.frame(
    game_home_margin_period_F3 = margins_f3,
    game_total_period_F3       = abs(margins_f3) + 6L,  # plausible total
    game_home_margin_period_F5 = margins_f3,            # placeholder
    game_total_period_F5       = abs(margins_f3) + 7L,
    game_home_margin_period_F7 = margins_f7,
    game_total_period_F7       = abs(margins_f7) + 8L,
    game_home_margin_period_FG = margins_f7,
    game_total_period_FG       = totals_fg
  )
  list(setNames(list(list(sample = samp)), game_id))[[1]]
}

make_consensus <- function(game_id = "test-game-1") {
  tibble(
    id = game_id,
    home_team = "Test Home",
    away_team = "Test Away",
    commence_time = as.POSIXct("2026-05-02 19:00:00", tz = "UTC")
  )
}

test_that("compare_alts_to_samples emits h2h_1st_3_innings bets when home wins F3 60% of the time", {
  # 600 home wins (margin > 0), 300 away wins (margin < 0), 100 ties
  margins <- c(rep(2L, 600), rep(-2L, 300), rep(0L, 100))
  samples <- make_synthetic_samples(margins_f3 = margins)
  consensus <- make_consensus()

  # Wagerzon-shaped row: away_ml -110 / home_ml -110 (so devigged fair = 50/50)
  # Model says home wins 60% of non-push → big edge on home.
  # game_date/game_time required by nearest_game_match() (called inside
  # compare_alts_to_samples) — these come from the standard 18-column scraper schema.
  offshore <- tibble(
    bookmaker_key = "wagerzon",
    home_team = "Test Home",
    away_team = "Test Away",
    market = "h2h_1st_3_innings",
    line = NA_real_,
    odds_away = -110L,
    odds_home = -110L,
    odds_over = NA_integer_,
    odds_under = NA_integer_,
    home_spread = NA_real_,
    away_spread = NA_real_,
    game_date = "2026-05-02",
    game_time = "19:00",
    commence_time = as.POSIXct("2026-05-02 19:00:00", tz = "UTC")
  )

  bets <- compare_alts_to_samples(
    samples = samples,
    offshore_odds = offshore,
    consensus_odds = consensus,
    bankroll = 100,
    kelly_mult = 0.25,
    ev_threshold = 0.02
  )

  expect_true(nrow(bets) >= 1)
  home_bet <- bets[bets$bet_on == "Test Home", ]
  expect_equal(nrow(home_bet), 1)
  # 600/(600+300) = 0.6667 — exclude pushes from denominator
  expect_equal(home_bet$prob, 600 / 900, tolerance = 1e-9)
  expect_equal(home_bet$market, "h2h_1st_3_innings")
})
