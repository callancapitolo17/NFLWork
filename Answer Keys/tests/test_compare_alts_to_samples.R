# Answer Keys/tests/test_compare_alts_to_samples.R
# Fixture-driven tests for compare_alts_to_samples derivative-h2h + odd/even branches.
# Run from "Answer Keys/" directory: Rscript -e 'testthat::test_file("tests/test_compare_alts_to_samples.R")'

library(testthat)
library(dplyr)
library(tibble)
library(data.table)
source("../Tools.R")

# Helper: build a synthetic samples list with one game.
# Every axis (F3/F7/FG margin and F3/F7/FG total) is independently overridable
# so downstream tests (Task 2 = derivative h2h, Task 3 = odd/even) can dial in
# whatever distribution they need without surprise coupling.
make_synthetic_samples <- function(
    game_id    = "test-game-1",
    margins_f3 = NULL,
    margins_f7 = NULL,
    margins_fg = NULL,
    totals_f3  = NULL,
    totals_f7  = NULL,
    totals_fg  = NULL
) {
  n <- 1000L
  if (is.null(margins_f3)) margins_f3 <- rep(0L, n)
  if (is.null(margins_f7)) margins_f7 <- rep(0L, n)
  if (is.null(margins_fg)) margins_fg <- margins_f7   # default: FG margin mirrors F7
  if (is.null(totals_f3))  totals_f3  <- abs(margins_f3) + 6L
  if (is.null(totals_f7))  totals_f7  <- abs(margins_f7) + 8L
  if (is.null(totals_fg))  totals_fg  <- rep(8L, n)
  stopifnot(
    length(margins_f3) == n, length(margins_f7) == n, length(margins_fg) == n,
    length(totals_f3) == n,  length(totals_f7) == n,  length(totals_fg) == n
  )
  samp <- data.frame(
    game_home_margin_period_F3 = margins_f3,
    game_total_period_F3       = totals_f3,
    game_home_margin_period_F5 = margins_f3,            # placeholder — F5 not exercised in this file
    game_total_period_F5       = abs(margins_f3) + 7L,  # placeholder
    game_home_margin_period_F7 = margins_f7,
    game_total_period_F7       = totals_f7,
    game_home_margin_period_FG = margins_fg,
    game_total_period_FG       = totals_fg
  )
  setNames(list(list(sample = samp)), game_id)
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
  # game_date/game_time are required by nearest_game_match() (called inside
  # compare_alts_to_samples). NOTE: real get_wagerzon_odds() output does NOT
  # carry a commence_time column — the consensus_odds side supplies it through
  # the inner join inside nearest_game_match. Adding commence_time to the
  # scraper side would create suffixed `commence_time.x` / `commence_time.y`
  # after the join and break the post-join `row$commence_time` lookup.
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
    game_time = "19:00"
  )

  bets <- compare_alts_to_samples(
    samples = samples,
    offshore_odds = offshore,
    consensus_odds = consensus,
    bankroll = 100,
    kelly_mult = 0.25,
    ev_threshold = 0.02
  )

  expect_equal(nrow(bets), 1L)
  home_bet <- bets[bets$bet_on == "Test Home", ]
  expect_equal(nrow(home_bet), 1)
  # 600/(600+300) = 0.6667 — exclude pushes from denominator
  expect_equal(home_bet$prob, 600 / 900, tolerance = 1e-9)
  expect_equal(home_bet$market, "h2h_1st_3_innings")
})

test_that("compare_alts_to_samples emits h2h_1st_7_innings bets using F7 margins", {
  # 700 home wins, 250 away wins, 50 ties — at F7
  margins <- c(rep(3L, 700), rep(-3L, 250), rep(0L, 50))
  samples <- make_synthetic_samples(margins_f7 = margins)
  consensus <- make_consensus()

  offshore <- tibble(
    bookmaker_key = "wagerzon",
    home_team = "Test Home",
    away_team = "Test Away",
    market = "h2h_1st_7_innings",
    line = NA_real_,
    odds_away = +200L,   # implied 33.3% — way too long if home is 73.7%
    odds_home = -150L,   # implied 60.0% — undervaluing home
    odds_over = NA_integer_,
    odds_under = NA_integer_,
    home_spread = NA_real_,
    away_spread = NA_real_,
    game_date = "2026-05-02",
    game_time = "19:00"
  )

  bets <- compare_alts_to_samples(
    samples = samples,
    offshore_odds = offshore,
    consensus_odds = consensus,
    bankroll = 100,
    kelly_mult = 0.25,
    ev_threshold = 0.02
  )

  home_bet <- bets[bets$bet_on == "Test Home" & bets$market == "h2h_1st_7_innings", ]
  expect_equal(nrow(home_bet), 1)
  expect_equal(home_bet$prob, 700 / 950, tolerance = 1e-9)
})

test_that("compare_alts_to_samples returns no bets when h2h row has NA odds", {
  margins <- c(rep(2L, 600), rep(-2L, 400))
  samples <- make_synthetic_samples(margins_f3 = margins)
  consensus <- make_consensus()

  offshore <- tibble(
    bookmaker_key = "wagerzon",
    home_team = "Test Home",
    away_team = "Test Away",
    market = "h2h_1st_3_innings",
    line = NA_real_,
    odds_away = NA_integer_,   # missing odds — must skip cleanly
    odds_home = NA_integer_,
    odds_over = NA_integer_,
    odds_under = NA_integer_,
    home_spread = NA_real_,
    away_spread = NA_real_,
    game_date = "2026-05-02",
    game_time = "19:00"
  )

  bets <- compare_alts_to_samples(
    samples = samples, offshore_odds = offshore, consensus_odds = consensus,
    bankroll = 100, kelly_mult = 0.25, ev_threshold = 0.02
  )

  expect_equal(nrow(bets), 0)
})

test_that("compare_alts_to_samples emits spreads_1st_3_innings bets — verifies suffix fix unblocks the spread branch for F3", {
  # Home covers -1.5 in 700/1000 sims → 70% cover. Wagerzon -110/-110 = 50/50 fair → big home edge.
  # Encode this with margins: 700 sims at margin=2 (home wins by 2, covers -1.5),
  # 300 sims at margin=-3 (away wins, home doesn't cover).
  margins <- c(rep(2L, 700), rep(-3L, 300))
  samples <- make_synthetic_samples(margins_f3 = margins)
  consensus <- make_consensus()

  offshore <- tibble(
    bookmaker_key = "wagerzon",
    home_team = "Test Home",
    away_team = "Test Away",
    game_date = "2026-05-02",
    game_time = "19:00",
    market = "spreads_1st_3_innings",
    line = -1.5,
    odds_away = -110L,
    odds_home = -110L,
    odds_over = NA_integer_,
    odds_under = NA_integer_,
    home_spread = -1.5,
    away_spread = 1.5
  )

  bets <- compare_alts_to_samples(
    samples = samples,
    offshore_odds = offshore,
    consensus_odds = consensus,
    bankroll = 100,
    kelly_mult = 0.25,
    ev_threshold = 0.02
  )

  # Home -1.5 covers in 700/1000 sims (margin > -(-1.5) = 1.5)
  home_bet <- bets[bets$bet_on == "Test Home" & bets$market == "spreads_1st_3_innings", ]
  expect_equal(nrow(home_bet), 1)
  expect_equal(home_bet$prob, 700 / 1000, tolerance = 1e-9)
})

test_that("compare_alts_to_samples emits odd_even_runs bet when totals are odd 60% of the time", {
  totals <- c(rep(7L, 600), rep(8L, 400))   # 60% odd, 40% even
  samples <- make_synthetic_samples(totals_fg = totals)
  consensus <- make_consensus()

  # Wagerzon convention: away_ml = ODD price, home_ml = EVEN price
  # (scraper_v2.py:446-453 — vtm = "TOTAL RUNS ODD", htm = "TOTAL RUNS EVEN")
  offshore <- tibble(
    bookmaker_key = "wagerzon",
    home_team = "Test Home",
    away_team = "Test Away",
    game_date = "2026-05-02",
    game_time = "19:00",
    market = "odd_even_runs",
    line = NA_real_,
    odds_away = -110L,   # ODD — fair if 52.4%
    odds_home = -110L,   # EVEN — fair if 52.4%
    odds_over = NA_integer_,
    odds_under = NA_integer_,
    home_spread = NA_real_,
    away_spread = NA_real_
  )

  bets <- compare_alts_to_samples(
    samples = samples, offshore_odds = offshore, consensus_odds = consensus,
    bankroll = 100, kelly_mult = 0.25, ev_threshold = 0.02
  )

  expect_equal(nrow(bets), 1L)             # exactly one bet emitted (Odd only; Even is -EV)
  odd_bet <- bets[bets$bet_on == "Odd" & bets$market == "odd_even_runs", ]
  expect_equal(nrow(odd_bet), 1)
  expect_equal(odd_bet$prob, 0.60, tolerance = 1e-9)
})

test_that("american_prob_3way devigs three American odds via sum-and-normalize", {
  # CLE @ ATH 3-way from the recon: away=+105, home=+130, draw=+475
  # Implied: 100/205, 100/230, 100/575 = 0.4878, 0.4348, 0.1739; sum=1.0965
  # Devigged: 0.4448, 0.3966, 0.1586
  result <- american_prob_3way(105, 130, 475)
  expect_equal(result$p_away, 0.4448, tolerance = 1e-3)
  expect_equal(result$p_home, 0.3966, tolerance = 1e-3)
  expect_equal(result$p_draw, 0.1586, tolerance = 1e-3)
  expect_equal(result$p_away + result$p_home + result$p_draw, 1.0, tolerance = 1e-9)
})

test_that("american_prob_3way handles negative odds", {
  # Heavy home favorite: away=+200, home=-150, draw=+500
  # Implied: 0.3333, 0.6, 0.1667; sum=1.1
  # Devigged: 0.3030, 0.5455, 0.1515
  result <- american_prob_3way(200, -150, 500)
  expect_equal(result$p_away, 0.3030, tolerance = 1e-3)
  expect_equal(result$p_home, 0.5455, tolerance = 1e-3)
  expect_equal(result$p_draw, 0.1515, tolerance = 1e-3)
})

test_that("american_prob_3way returns NA for any NA input", {
  result <- american_prob_3way(NA_integer_, 130, 475)
  expect_true(is.na(result$p_away))
  expect_true(is.na(result$p_home))
  expect_true(is.na(result$p_draw))
})
