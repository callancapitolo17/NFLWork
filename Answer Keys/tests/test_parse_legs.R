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

test_that("parse_legs handles triple-play (3 legs)", {
  expect_equal(
    parse_legs("GIANTS TRIPLE-PLAY (SCR 1ST, 1H & GM)"),
    list(
      list(type = "scores_first"),
      list(type = "wins_period", period = "F5"),
      list(type = "wins_period", period = "FG")
    )
  )
})

test_that("parse_legs handles grand slam with unicode fraction", {
  expect_equal(
    parse_legs("GIANTS GRAND-SLAM (SCR 1ST, 1H, GM & SCR U2½)"),
    list(
      list(type = "scores_first"),
      list(type = "wins_period", period = "F5"),
      list(type = "wins_period", period = "FG"),
      list(type = "team_total_under", line = 2.5)
    )
  )
})

test_that("parse_legs handles grand slam with decimal fraction", {
  expect_equal(
    parse_legs("GIANTS GRAND-SLAM (SCR 1ST, 1H, GM & SCR U2.5)"),
    list(
      list(type = "scores_first"),
      list(type = "wins_period", period = "F5"),
      list(type = "wins_period", period = "FG"),
      list(type = "team_total_under", line = 2.5)
    )
  )
})

test_that("parse_legs handles F3 and F7 period tokens", {
  expect_equal(
    parse_legs("GIANTS MEGA (SCR 1ST, F3, 1H, F7 & GM)"),
    list(
      list(type = "scores_first"),
      list(type = "wins_period", period = "F3"),
      list(type = "wins_period", period = "F5"),
      list(type = "wins_period", period = "F7"),
      list(type = "wins_period", period = "FG")
    )
  )
})

test_that("parse_legs handles SCR O<N> team-total over", {
  expect_equal(
    parse_legs("GIANTS BIG (SCR 1ST, GM & SCR O4.5)"),
    list(
      list(type = "scores_first"),
      list(type = "wins_period", period = "FG"),
      list(type = "team_total_over", line = 4.5)
    )
  )
})

test_that("parse_legs returns NULL and warns on unknown token", {
  expect_warning(
    result <- parse_legs("GIANTS MYSTERY (SCR 1ST, XYZ & GM)"),
    "Unknown leg token"
  )
  expect_null(result)
})

test_that("parse_legs returns NULL on description with no parenthetical", {
  expect_null(parse_legs("GIANTS TRIPLE-PLAY"))
})
