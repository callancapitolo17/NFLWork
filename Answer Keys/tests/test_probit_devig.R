# Answer Keys/tests/test_probit_devig.R
library(testthat)
source("../Tools.R")

# ---------------------------------------------------------------------------
# 2-way devig_american
# ---------------------------------------------------------------------------

test_that("devig_american(-110, -110) gives 50/50 (symmetric case)", {
  result <- devig_american(-110, -110)
  expect_equal(result$p1, 0.5, tolerance = 1e-9)
  expect_equal(result$p2, 0.5, tolerance = 1e-9)
})

test_that("devig_american(-200, +180) sums to 1.0", {
  result <- devig_american(-200, 180)
  expect_equal(result$p1 + result$p2, 1.0, tolerance = 1e-9)
})

test_that("devig_american(-200, +180) matches frozen probit values", {
  # Computed by reference probit math: p_raw = (0.6667, 0.3571);
  # z = (qnorm(0.6667), qnorm(0.3571)); c = -(z1+z2)/2; output = pnorm(z+c)
  result <- devig_american(-200, 180)
  expect_equal(result$p1, 0.6548, tolerance = 0.001)
  expect_equal(result$p2, 0.3452, tolerance = 0.001)
})

test_that("devig_american(-1000, +800) tail case matches frozen probit values", {
  # p_raw = (0.9091, 0.1111); probit gives favorite ~0.8995, dog ~0.1005
  result <- devig_american(-1000, 800)
  expect_equal(result$p1, 0.8995, tolerance = 0.001)
  expect_equal(result$p2, 0.1005, tolerance = 0.001)
})

test_that("devig_american(0, +100) returns NA on invalid input", {
  result <- devig_american(0, 100)
  expect_true(is.na(result$p1))
  expect_true(is.na(result$p2))
})

test_that("devig_american(NA, -110) returns NA on NA input", {
  result <- devig_american(NA, -110)
  expect_true(is.na(result$p1))
  expect_true(is.na(result$p2))
})

# ---------------------------------------------------------------------------
# Negative test: probit != multiplicative at tails
# ---------------------------------------------------------------------------

test_that("probit output diverges from multiplicative on tail case", {
  result <- devig_american(-1000, 800)

  # Compute multiplicative output for comparison
  p1_raw <- 1000 / 1100
  p2_raw <- 100 / 900
  mult_p1 <- p1_raw / (p1_raw + p2_raw)

  # Must differ by at least 0.001 — catches accidental method regression
  expect_gt(abs(result$p1 - mult_p1), 0.001)
})

# ---------------------------------------------------------------------------
# 3-way devig_american_3way
# ---------------------------------------------------------------------------

test_that("devig_american_3way sums to 1.0", {
  result <- devig_american_3way(150, 200, 400)
  expect_equal(result$p_home + result$p_away + result$p_tie, 1.0, tolerance = 1e-9)
})

test_that("devig_american_3way falls back to 2-way when tie is NA", {
  result <- devig_american_3way(-110, -110, NA)
  expect_equal(result$p_home, 0.5, tolerance = 1e-9)
  expect_equal(result$p_away, 0.5, tolerance = 1e-9)
  expect_true(is.na(result$p_tie))
})

test_that("devig_american_3way invalid input returns NA", {
  result <- devig_american_3way(0, 150, 200)
  expect_true(is.na(result$p_home))
  expect_true(is.na(result$p_away))
  expect_true(is.na(result$p_tie))
})
