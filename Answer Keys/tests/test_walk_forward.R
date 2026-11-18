# Answer Keys/tests/test_walk_forward.R
# Regression guard for Tools.R::walk_forward_training_set — the rule that a
# backtest prices each game only off games played BEFORE it, as live does.
#
# Background (2026-11-17): every backtest trained on all games except the test
# game (leave-one-out), so later seasons leaked into the sample.
#
# Run from "Answer Keys/" directory:
#   Rscript -e 'testthat::test_file("tests/test_walk_forward.R")'

library(testthat)
library(data.table)
source("../Tools.R")

.history <- function(game_dates) {
  data.table(game_id = seq_along(game_dates), game_date = as.Date(game_dates))
}

test_that("only games dated strictly before the test game are kept", {
  DT <- .history(c("2024-01-01", "2024-01-02", "2024-01-03", "2024-01-04"))
  train <- walk_forward_training_set(DT, as.Date("2024-01-04"), min_train_share = 0)
  expect_equal(train$game_id, 1:3)
})

test_that("same-day games are excluded, including the test game itself", {
  DT <- .history(c("2024-01-01", "2024-01-02", "2024-01-02", "2024-01-03"))
  train <- walk_forward_training_set(DT, as.Date("2024-01-02"), min_train_share = 0)
  expect_equal(train$game_id, 1L)
})

test_that("later games never enter the training set", {
  DT <- .history(c("2023-12-01", "2024-02-01", "2025-02-01"))
  train <- walk_forward_training_set(DT, as.Date("2024-01-15"), min_train_share = 0)
  expect_true(all(train$game_date < as.Date("2024-01-15")))
  expect_equal(nrow(train), 1)
})

test_that("returns NULL when less than min_train_share of history precedes the test game", {
  DT <- .history(sprintf("2024-01-%02d", 1:10))
  expect_null(walk_forward_training_set(DT, as.Date("2024-01-05"), min_train_share = 0.5))  # 4 of 10
  expect_equal(nrow(walk_forward_training_set(DT, as.Date("2024-01-06"), min_train_share = 0.5)), 5)
})

test_that("default share is WALK_FORWARD_MIN_TRAIN_SHARE", {
  DT <- .history(sprintf("2024-01-%02d", 1:10))
  expect_null(walk_forward_training_set(DT, as.Date("2024-01-02")))
  expect_equal(nrow(walk_forward_training_set(DT, as.Date("2024-01-10"))), 9)
})

test_that("POSIXct and character test dates work like Date", {
  DT <- .history(c("2024-01-01", "2024-01-02", "2024-01-03"))
  as_posix <- walk_forward_training_set(DT, as.POSIXct("2024-01-03 19:00:00", tz = "UTC"), min_train_share = 0)
  as_char  <- walk_forward_training_set(DT, "2024-01-03", min_train_share = 0)
  expect_equal(as_posix$game_id, 1:2)
  expect_equal(as_char$game_id, 1:2)
})

test_that("rows with an NA game_date are never used for training", {
  DT <- data.table(game_id = 1:3, game_date = as.Date(c("2024-01-01", NA, "2024-01-02")))
  train <- walk_forward_training_set(DT, as.Date("2024-01-05"), min_train_share = 0)
  expect_equal(train$game_id, c(1L, 3L))
})

test_that("fails loudly without a game_date column or with a bad test date", {
  expect_error(walk_forward_training_set(data.table(game_id = 1:3), as.Date("2024-01-01")),
               "expected a game_date column")
  DT <- .history(c("2024-01-01", "2024-01-02"))
  expect_error(walk_forward_training_set(DT, NA), "expected one non-NA test_game_date")
  expect_error(walk_forward_training_set(DT, as.Date(c("2024-01-01", "2024-01-02"))),
               "expected one non-NA test_game_date")
})
