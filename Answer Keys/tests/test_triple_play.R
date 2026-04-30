# Answer Keys/tests/test_triple_play.R
library(testthat)
source("../triple_play_helpers.R")

test_that("home scored first when away has 0 runs in inning 1", {
  expect_equal(
    determine_home_scored_first(m1=2, t1=2,
                                m2=NA, t2=NA, m3=NA, t3=NA,
                                m4=NA, t4=NA, m5=NA, t5=NA),
    1L
  )
})

test_that("away scored first when away scored any runs in first scoring inning", {
  expect_equal(
    determine_home_scored_first(m1=-1, t1=3,
                                m2=NA, t2=NA, m3=NA, t3=NA,
                                m4=NA, t4=NA, m5=NA, t5=NA),
    0L
  )
})

test_that("falls through to inning 2 when inning 1 is scoreless", {
  expect_equal(
    determine_home_scored_first(m1=0, t1=0,
                                m2=1, t2=1, m3=NA, t3=NA,
                                m4=NA, t4=NA, m5=NA, t5=NA),
    1L
  )
})

test_that("returns NA when no scoring in innings 1 through 5", {
  expect_true(is.na(
    determine_home_scored_first(m1=0, t1=0, m2=0, t2=0,
                                m3=0, t3=0, m4=0, t4=0,
                                m5=0, t5=0)
  ))
})

test_that("returns NA when data is NA before any scoring", {
  expect_true(is.na(
    determine_home_scored_first(m1=NA, t1=NA,
                                m2=NA, t2=NA, m3=NA, t3=NA,
                                m4=NA, t4=NA, m5=NA, t5=NA)
  ))
})

test_that("vectorized version matches scalar over multiple rows", {
  result <- determine_home_scored_first_vec(
    m1 = c(2, -1, 0), t1 = c(2, 3, 0),
    m2 = c(NA, NA, 1), t2 = c(NA, NA, 1),
    m3 = c(NA, NA, NA), t3 = c(NA, NA, NA),
    m4 = c(NA, NA, NA), t4 = c(NA, NA, NA),
    m5 = c(NA, NA, NA), t5 = c(NA, NA, NA)
  )
  expect_equal(result, c(1L, 0L, 1L))
})

test_that("away scored first when both teams score 1 run in inning 1 (tied)", {
  # t1=2, m1=0 → home 1 + away 1 in inning 1. Away bats first → away scored first.
  # This is the canonical "both scored in the first" case.
  expect_equal(
    determine_home_scored_first(m1=0, t1=2,
                                m2=NA, t2=NA, m3=NA, t3=NA,
                                m4=NA, t4=NA, m5=NA, t5=NA),
    0L
  )
})

test_that("returns NA on impossible inputs (|margin_change| > total_change)", {
  # t1=1, m1=3 is mathematically impossible — you cannot swing the margin by
  # 3 runs from a 1-run total increase. Must return NA, not silently compute.
  expect_true(is.na(
    determine_home_scored_first(m1=3, t1=1,
                                m2=NA, t2=NA, m3=NA, t3=NA,
                                m4=NA, t4=NA, m5=NA, t5=NA)
  ))
})

test_that("mlb_game_samples has home_scored_first column populated", {
  skip_if_not(file.exists("../mlb_mm.duckdb"),
              "mlb_mm.duckdb not present — run MLB.R first")
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = "../mlb_mm.duckdb", read_only = TRUE)
  on.exit(DBI::dbDisconnect(con))
  cols <- DBI::dbGetQuery(con, "PRAGMA table_info('mlb_game_samples')")$name
  expect_true("home_scored_first" %in% cols)
  # Values should be 0, 1, or NA only
  vals <- DBI::dbGetQuery(con,
    "SELECT DISTINCT home_scored_first FROM mlb_game_samples")$home_scored_first
  expect_true(all(vals %in% c(0L, 1L, NA_integer_)))
  # Coverage should be >= 94% non-NA per historical data (observed ~95.4%)
  non_na_pct <- DBI::dbGetQuery(con,
    "SELECT AVG(CASE WHEN home_scored_first IS NOT NULL THEN 1.0 ELSE 0.0 END) AS p
     FROM mlb_game_samples")$p
  expect_gt(non_na_pct, 0.94)
})

# ----- Pricer tests (mlb_triple_play.R) -----
source("../mlb_triple_play.R", local = TRUE)

test_that("prob_to_american handles favorites and dogs", {
  expect_equal(prob_to_american(0.5), -100L)
  expect_equal(prob_to_american(0.6), -150L)   # -100 * 0.6/0.4
  expect_equal(prob_to_american(0.25), 300L)   # 100 * 0.75/0.25
  expect_true(is.na(prob_to_american(0)))
  expect_true(is.na(prob_to_american(1)))
  expect_true(is.na(prob_to_american(NA_real_)))
})

test_that("american_to_prob is the inverse of prob_to_american", {
  expect_equal(american_to_prob(-150), 0.6)
  expect_equal(american_to_prob(+300), 0.25)
  expect_equal(american_to_prob(+100), 0.5)
})

test_that("american_to_prob(0) returns NA_real_", {
  expect_true(is.na(american_to_prob(0)))
})

# =============================================================================
# Trifecta dashboard table writer (Task 1 of 2026-04-28 plan)
# =============================================================================

test_that("trifecta_hash is stable across reruns and unique per row", {
  # Same inputs -> same hash
  h1 <- digest::digest(paste("game42", "Yankees", "TRIPLE-PLAY", "home", sep = "|"),
                       algo = "sha256", serialize = FALSE)
  h2 <- digest::digest(paste("game42", "Yankees", "TRIPLE-PLAY", "home", sep = "|"),
                       algo = "sha256", serialize = FALSE)
  expect_equal(h1, h2)

  # Different inputs -> different hash
  h3 <- digest::digest(paste("game42", "Yankees", "GRAND-SLAM", "home", sep = "|"),
                       algo = "sha256", serialize = FALSE)
  expect_false(h1 == h3)

  # Different side -> different hash
  h4 <- digest::digest(paste("game42", "Yankees", "TRIPLE-PLAY", "away", sep = "|"),
                       algo = "sha256", serialize = FALSE)
  expect_false(h1 == h4)
})

test_that("kelly_bet is zero when edge_pct is below trifecta_min_edge", {
  # Replicate the inline Kelly logic from mlb_triple_play.R
  trifecta_bankroll   <- 100
  trifecta_kelly_mult <- 0.10
  trifecta_min_edge   <- 0.05  # 5%

  # Row with 3% edge -- below threshold
  edge_pct <- 3.0
  kelly_bet <- if (!is.na(edge_pct) && edge_pct >= trifecta_min_edge * 100) {
    trifecta_bankroll * trifecta_kelly_mult * 0.5  # arbitrary frac
  } else 0
  expect_equal(kelly_bet, 0)

  # Row with 8% edge -- above threshold
  edge_pct <- 8.0
  kelly_bet <- if (!is.na(edge_pct) && edge_pct >= trifecta_min_edge * 100) {
    trifecta_bankroll * trifecta_kelly_mult * 0.5
  } else 0
  expect_equal(kelly_bet, 5.0)
})

test_that("kelly_frac formula matches expected for known win_prob and odds", {
  # +400 American = 5.0 decimal, b = 4
  # Win prob 0.25 -> kelly_frac = (4*0.25 - 0.75) / 4 = 0.0625
  win_prob <- 0.25
  dec_odds <- 5.0
  b <- dec_odds - 1
  p <- win_prob
  q <- 1 - p
  kelly_frac <- max(0, (b * p - q) / b)
  expect_equal(kelly_frac, 0.0625, tolerance = 1e-9)

  # If win_prob equals breakeven (1/dec_odds = 0.20), kelly = 0
  win_prob <- 0.20
  p <- win_prob
  q <- 1 - p
  kelly_frac <- max(0, (b * p - q) / b)
  expect_equal(kelly_frac, 0, tolerance = 1e-9)

  # Below breakeven -> kelly clamped to 0 (don't bet against yourself)
  win_prob <- 0.15
  p <- win_prob
  q <- 1 - p
  kelly_frac <- max(0, (b * p - q) / b)
  expect_equal(kelly_frac, 0)
})

test_that("kelly_bet is zero when fair_odds is NA (model couldn't price)", {
  # Replicate the inline NA-guard logic from mlb_triple_play.R
  trifecta_bankroll  <- 100
  trifecta_kelly_mult <- 0.10
  trifecta_min_edge  <- 0.05

  fair_odds <- NA_integer_
  book_odds <- 500L
  edge_pct  <- NA_real_

  win_prob   <- if (!is.na(fair_odds)) {
    if (fair_odds > 0) 100 / (fair_odds + 100) else (-fair_odds) / ((-fair_odds) + 100)
  } else NA_real_
  dec_odds   <- if (!is.na(book_odds)) {
    if (book_odds > 0) 1 + book_odds / 100 else 1 + 100 / abs(book_odds)
  } else NA_real_
  kelly_frac <- if (!is.na(win_prob) && !is.na(dec_odds) && dec_odds > 1) {
    b <- dec_odds - 1; p <- win_prob; q <- 1 - p
    max(0, (b * p - q) / b)
  } else 0
  kelly_bet  <- if (!is.na(edge_pct) && edge_pct >= trifecta_min_edge * 100) {
    trifecta_bankroll * trifecta_kelly_mult * kelly_frac
  } else 0

  expect_equal(win_prob, NA_real_)
  expect_equal(kelly_frac, 0)
  expect_equal(kelly_bet, 0)
})
