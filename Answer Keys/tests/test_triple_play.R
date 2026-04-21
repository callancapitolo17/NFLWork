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
  skip_if_not(file.exists("../mlb.duckdb"),
              "mlb.duckdb not present — run MLB.R first")
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = "../mlb.duckdb", read_only = TRUE)
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

test_that("compute_triple_play_fair (home) returns joint P of 3 legs", {
  # 10 synthetic rows. home_triple = scored_first==1 AND margin_f5>0 AND home_margin>0.
  # Hits: rows 1 (1,1,2), 2 (1,2,3), 4 (1,3,4), 6 (1,2,5), 9 (1,4,6) → 5/10 = 0.5
  samples <- data.frame(
    home_margin       = c( 2,  3, -1,  4,  1,  5, -2,  2,  6, -3),
    home_margin_f5    = c( 1,  2,  0,  3,  0,  2, -1,  1,  4, -2),
    home_scored_first = c( 1L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 1L, 0L)
  )
  expect_equal(compute_triple_play_fair(samples, side = "home"), 0.5)
})

test_that("compute_triple_play_fair (away) returns joint P for away side", {
  # away_triple = scored_first==0 AND margin_f5<0 AND home_margin<0
  # Hits: rows 1, 2, 4, 6, 9 → 5/10 = 0.5
  samples <- data.frame(
    home_margin       = c(-2, -3,  1, -4,  1, -5,  2, -2, -6,  3),
    home_margin_f5    = c(-1, -2,  0, -3,  0, -2,  1, -1, -4,  2),
    home_scored_first = c( 0L, 0L, 1L, 0L, 0L, 0L, 1L, 1L, 0L, 1L)
  )
  expect_equal(compute_triple_play_fair(samples, side = "away"), 0.5)
})

test_that("compute_triple_play_fair excludes NA home_scored_first rows", {
  # 6 rows — 3 valid with 1 hit, 3 NA. Fair should be 1/3, not 1/6.
  samples <- data.frame(
    home_margin       = c( 2,  3,  1,  4,  5,  6),
    home_margin_f5    = c( 1,  2,  0,  3,  4,  5),
    home_scored_first = c( 1L, 0L, 1L, NA, NA, NA)
  )
  expect_equal(compute_triple_play_fair(samples, side = "home"), 1/3)
})

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
