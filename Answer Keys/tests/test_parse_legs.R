# Answer Keys/tests/test_parse_legs.R
library(testthat)
source("../parse_legs.R")

test_that("parse_fraction handles integer input", {
  expect_equal(parse_fraction("3"), 3)
})

test_that("parse_fraction handles decimal input", {
  expect_equal(parse_fraction("2.5"), 2.5)
})

test_that("parse_fraction handles ½ unicode half", {
  expect_equal(parse_fraction("2½"), 2.5)
})

test_that("parse_fraction handles ¼ ¾ unicode quarters", {
  expect_equal(parse_fraction("1¼"), 1.25)
  expect_equal(parse_fraction("3¾"), 3.75)
})

test_that("parse_fraction returns NA on unparseable input", {
  expect_true(is.na(parse_fraction("abc")))
  expect_true(is.na(parse_fraction("")))
})
