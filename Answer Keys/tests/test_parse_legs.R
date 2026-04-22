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

test_that("eval_leg scores_first home matches home_scored_first == 1L", {
  samples <- data.frame(
    home_scored_first = c(1L, 0L, 1L, NA)
  )
  team_runs <- c(NA, NA, NA, NA)
  opp_runs  <- c(NA, NA, NA, NA)
  result <- eval_leg(list(type = "scores_first"),
                     samples, side = "home",
                     team_runs = team_runs, opp_runs = opp_runs)
  expect_equal(result, c(TRUE, FALSE, TRUE, NA))
})

test_that("eval_leg scores_first away matches home_scored_first == 0L", {
  samples <- data.frame(home_scored_first = c(1L, 0L, 1L))
  result <- eval_leg(list(type = "scores_first"),
                     samples, side = "away",
                     team_runs = c(NA, NA, NA), opp_runs = c(NA, NA, NA))
  expect_equal(result, c(FALSE, TRUE, FALSE))
})

test_that("eval_leg wins_period F5 home checks home_margin_f5 > 0", {
  samples <- data.frame(home_margin_f5 = c(2, 0, -1, 3))
  result <- eval_leg(list(type = "wins_period", period = "F5"),
                     samples, side = "home",
                     team_runs = c(NA, NA, NA, NA), opp_runs = c(NA, NA, NA, NA))
  expect_equal(result, c(TRUE, FALSE, FALSE, TRUE))
})

test_that("eval_leg wins_period F5 away checks home_margin_f5 < 0 (strict)", {
  samples <- data.frame(home_margin_f5 = c(2, 0, -1, -3))
  result <- eval_leg(list(type = "wins_period", period = "F5"),
                     samples, side = "away",
                     team_runs = c(NA, NA, NA, NA), opp_runs = c(NA, NA, NA, NA))
  expect_equal(result, c(FALSE, FALSE, TRUE, TRUE))
})

test_that("eval_leg wins_period FG uses home_margin column", {
  samples <- data.frame(home_margin = c(1, -1, 0))
  result <- eval_leg(list(type = "wins_period", period = "FG"),
                     samples, side = "home",
                     team_runs = c(NA, NA, NA), opp_runs = c(NA, NA, NA))
  expect_equal(result, c(TRUE, FALSE, FALSE))
})

test_that("eval_leg wins_period F3 uses home_margin_f3 column", {
  samples <- data.frame(home_margin_f3 = c(2, -1, 0))
  result <- eval_leg(list(type = "wins_period", period = "F3"),
                     samples, side = "home",
                     team_runs = c(NA, NA, NA), opp_runs = c(NA, NA, NA))
  expect_equal(result, c(TRUE, FALSE, FALSE))
})

test_that("eval_leg wins_period F7 uses home_margin_f7 column", {
  samples <- data.frame(home_margin_f7 = c(2, -1, 0))
  result <- eval_leg(list(type = "wins_period", period = "F7"),
                     samples, side = "home",
                     team_runs = c(NA, NA, NA), opp_runs = c(NA, NA, NA))
  expect_equal(result, c(TRUE, FALSE, FALSE))
})

test_that("eval_leg team_total_under uses team_runs (side-dependent)", {
  samples <- data.frame(x = c(1, 2, 3))  # unused
  result <- eval_leg(list(type = "team_total_under", line = 2.5),
                     samples, side = "home",
                     team_runs = c(1, 2, 3), opp_runs = c(5, 5, 5))
  expect_equal(result, c(TRUE, TRUE, FALSE))
})

test_that("eval_leg team_total_over uses team_runs", {
  samples <- data.frame(x = c(1, 2, 3))
  result <- eval_leg(list(type = "team_total_over", line = 2.5),
                     samples, side = "home",
                     team_runs = c(1, 2, 3), opp_runs = c(5, 5, 5))
  expect_equal(result, c(FALSE, FALSE, TRUE))
})

test_that("eval_leg unknown type raises error", {
  samples <- data.frame(x = 1)
  expect_error(
    eval_leg(list(type = "nonexistent"),
             samples, side = "home",
             team_runs = c(1), opp_runs = c(1)),
    "Unknown leg type"
  )
})
