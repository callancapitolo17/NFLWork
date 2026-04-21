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
  # Coverage should be >= 90% non-NA per historical data
  non_na_pct <- DBI::dbGetQuery(con,
    "SELECT AVG(CASE WHEN home_scored_first IS NOT NULL THEN 1.0 ELSE 0.0 END) AS p
     FROM mlb_game_samples")$p
  expect_gt(non_na_pct, 0.90)
})
