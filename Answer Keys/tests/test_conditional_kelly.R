library(testthat)
source("../conditional_kelly.R")

test_that("conditional Kelly returns non-negative residuals", {
  res <- conditional_kelly_residuals(
    p_a = 0.50, d_a = 2.10,
    p_b = 0.50, d_b = 2.10,
    s_combo = 50, d_combo = 4.30,
    bankroll = 1000, kelly_mult = 1.0
  )
  expect_true(res$s_a >= 0)
  expect_true(res$s_b >= 0)
})

test_that("conditional Kelly residual <= unconstrained Kelly", {
  full_kelly_a <- ((0.55 * 2.10) - 1) / 1.10 * 1000  # ~$45
  res <- conditional_kelly_residuals(
    p_a = 0.55, d_a = 2.10,
    p_b = 0.55, d_b = 2.10,
    s_combo = 30, d_combo = 4.30,
    bankroll = 1000, kelly_mult = 1.0
  )
  expect_true(res$s_a <= full_kelly_a + 0.01)
  expect_true(res$s_b <= full_kelly_a + 0.01)
})

test_that("conditional Kelly with s_combo = 0 matches independent Kelly", {
  full_kelly_a <- ((0.55 * 2.10) - 1) / 1.10 * 1000
  res <- conditional_kelly_residuals(
    p_a = 0.55, d_a = 2.10,
    p_b = 0.55, d_b = 2.10,
    s_combo = 0, d_combo = 4.30,
    bankroll = 1000, kelly_mult = 1.0
  )
  expect_equal(res$s_a, full_kelly_a, tolerance = 1.0)
  expect_equal(res$s_b, full_kelly_a, tolerance = 1.0)
})

test_that("conditional Kelly with kelly_mult halves the residual", {
  full <- conditional_kelly_residuals(
    p_a = 0.55, d_a = 2.10, p_b = 0.55, d_b = 2.10,
    s_combo = 30, d_combo = 4.30,
    bankroll = 1000, kelly_mult = 1.0
  )
  half <- conditional_kelly_residuals(
    p_a = 0.55, d_a = 2.10, p_b = 0.55, d_b = 2.10,
    s_combo = 30, d_combo = 4.30,
    bankroll = 1000, kelly_mult = 0.5
  )
  expect_equal(half$s_a, full$s_a * 0.5, tolerance = 0.5)
})

test_that("high-edge legs do not crash the optimizer", {
  res <- conditional_kelly_residuals(
    p_a = 0.85, d_a = 2.10,
    p_b = 0.85, d_b = 2.10,
    s_combo = 50, d_combo = 4.30,
    bankroll = 1000
  )
  expect_true(is.finite(res$s_a) && res$s_a >= 0)
  expect_true(is.finite(res$s_b) && res$s_b >= 0)
})

test_that("negative-edge legs return zero residuals", {
  res <- conditional_kelly_residuals(
    p_a = 0.40, d_a = 2.10,  # implied 47.6%, true 40% -> negative EV
    p_b = 0.40, d_b = 2.10,
    s_combo = 0, d_combo = 4.30,
    bankroll = 1000
  )
  expect_equal(res$s_a, 0)
  expect_equal(res$s_b, 0)
})

test_that("kelly_mult = 0 yields zero residuals", {
  res <- conditional_kelly_residuals(
    p_a = 0.55, d_a = 2.10, p_b = 0.55, d_b = 2.10,
    s_combo = 30, d_combo = 4.30,
    bankroll = 1000, kelly_mult = 0
  )
  expect_equal(res$s_a, 0)
  expect_equal(res$s_b, 0)
})

test_that("zero bankroll returns zero residuals", {
  res <- conditional_kelly_residuals(
    p_a = 0.55, d_a = 2.10, p_b = 0.55, d_b = 2.10,
    s_combo = 0, d_combo = 4.30,
    bankroll = 0
  )
  expect_equal(res$s_a, 0)
  expect_equal(res$s_b, 0)
})
