library(testthat)
source("../triple_play_helpers.R")

test_that("home and away both score in 1st inning", {
  # Top of 1st: away scores 1 (margin -1, total 1).
  # Bottom of 1st: home scores 2 (margin +1 cumulative, total 3).
  result <- determine_inning_1_scoring(m1 = 1, t1 = 3)
  expect_equal(result$home, 1L)
  expect_equal(result$away, 1L)
})

test_that("only away scores in 1st inning", {
  # Away scores 2, home scores 0. margin = -2, total = 2.
  result <- determine_inning_1_scoring(m1 = -2, t1 = 2)
  expect_equal(result$home, 0L)
  expect_equal(result$away, 1L)
})

test_that("only home scores in 1st inning", {
  # Away scores 0, home scores 1. margin = +1, total = 1.
  result <- determine_inning_1_scoring(m1 = 1, t1 = 1)
  expect_equal(result$home, 1L)
  expect_equal(result$away, 0L)
})

test_that("neither team scores in 1st inning", {
  result <- determine_inning_1_scoring(m1 = 0, t1 = 0)
  expect_equal(result$home, 0L)
  expect_equal(result$away, 0L)
})

test_that("NA inputs propagate to NA outputs", {
  result <- determine_inning_1_scoring(m1 = NA, t1 = NA)
  expect_true(is.na(result$home))
  expect_true(is.na(result$away))
})

test_that("impossible inputs return NA", {
  # |margin| > total is impossible — corrupt data.
  result <- determine_inning_1_scoring(m1 = 5, t1 = 1)
  expect_true(is.na(result$home))
  expect_true(is.na(result$away))
})

test_that("vec wrapper returns parallel vectors", {
  out <- determine_inning_1_scoring_vec(
    m1 = c(1, -2, 0, NA),
    t1 = c(3,  2, 0, NA)
  )
  expect_equal(out$home, c(1L, 0L, 0L, NA_integer_))
  expect_equal(out$away, c(1L, 1L, 0L, NA_integer_))
})
