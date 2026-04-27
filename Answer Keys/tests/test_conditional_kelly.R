library(testthat)
source("../conditional_kelly.R")
source("../Tools.R")

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

test_that("compute_combined_parlay_pricing multiplies fair decimal odds", {
  res <- compute_combined_parlay_pricing(
    fair_dec_a = 4.32, fair_dec_b = 4.85,
    wz_dec = 22.10
  )
  expect_equal(res$joint_fair_dec, 4.32 * 4.85, tolerance = 0.001)
  expect_equal(res$joint_fair_prob, 1 / (4.32 * 4.85), tolerance = 0.0001)
})

test_that("compute_combined_parlay_pricing computes joint edge correctly", {
  res <- compute_combined_parlay_pricing(
    fair_dec_a = 4.32, fair_dec_b = 4.85,
    wz_dec = 22.10
  )
  expected_edge <- (1 / (4.32 * 4.85)) * 22.10 - 1
  expect_equal(res$joint_edge, expected_edge, tolerance = 0.0001)
})

test_that("compute_combined_parlay_pricing returns kelly stake using parlay_kelly_mult", {
  res <- compute_combined_parlay_pricing(
    fair_dec_a = 4.32, fair_dec_b = 4.85,
    wz_dec = 22.10,
    bankroll = 1000, kelly_mult = 0.5
  )
  expect_true(res$kelly_stake >= 0)
  expect_true(res$kelly_stake <= 500)  # never more than half bankroll at Kelly_mult = 0.5
})

test_that("compute_combined_parlay_pricing floors kelly_stake at 0 on negative edge", {
  # Two coinflips combined, book pays only 3.0 (true fair = 4.0) -> -EV bet
  res <- compute_combined_parlay_pricing(
    fair_dec_a = 2.0, fair_dec_b = 2.0,
    wz_dec = 3.0,
    bankroll = 1000, kelly_mult = 0.5
  )
  expect_true(res$joint_edge < 0)
  expect_equal(res$kelly_stake, 0)
})
