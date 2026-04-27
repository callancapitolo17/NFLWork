# Answer Keys/tests/test_dk_blend.R
library(testthat)
source("../parse_legs.R")

test_that("blend returns model_prob when DK decimal is NA", {
  expect_equal(blend_dk_with_model(0.30, NA_real_, 1.10), 0.30)
})

test_that("blend averages model and DK devigged when both present", {
  # DK SGP +200 → decimal 3.0 → implied 0.3333... Devigged at vig=1.10 → 0.30303
  # Model 0.30. Expected blend = mean(0.30, 0.30303) ≈ 0.30152
  expect_equal(blend_dk_with_model(0.30, 3.0, 1.10),
               (0.30 + (1/3.0) / 1.10) / 2,
               tolerance = 1e-9)
})

test_that("blend treats DK decimal of 0 or negative as missing", {
  expect_equal(blend_dk_with_model(0.30, 0,    1.10), 0.30)
  expect_equal(blend_dk_with_model(0.30, -1.0, 1.10), 0.30)
})

test_that("blend returns DK alone when model is NA", {
  # Edge case: if compute_prop_fair returned NA (e.g. zero valid samples),
  # DK is the only available signal.
  expect_equal(blend_dk_with_model(NA_real_, 3.0, 1.10),
               (1/3.0) / 1.10,
               tolerance = 1e-9)
})

test_that("blend is NA when both inputs are NA / missing", {
  expect_true(is.na(blend_dk_with_model(NA_real_, NA_real_, 1.10)))
  expect_true(is.na(blend_dk_with_model(NA_real_, 0,        1.10)))
})
