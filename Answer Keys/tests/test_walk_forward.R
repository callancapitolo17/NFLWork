# Answer Keys/tests/test_walk_forward.R
# Regression guard for Tools.R::walk_forward_training_set / has_walk_forward_history:
# a backtest prices each game only off games dated before it, as live does.
#
# Run from "Answer Keys/" directory:
#   Rscript -e 'testthat::test_file("tests/test_walk_forward.R")'

library(testthat)
library(data.table)
source("../Tools.R")

.history <- function(game_dates) {
  data.table(game_id = seq_along(game_dates), game_date = as.Date(game_dates),
             total_line = 140 + seq_along(game_dates))
}

test_that("training set is DT's rows dated strictly before the test game, unchanged", {
  DT <- .history(c("2024-01-01", "2024-01-02", "2024-01-03", "2024-01-04"))
  train <- walk_forward_training_set(DT, as.Date("2024-01-04"), min_train_share = 0)
  expect_equal(train, DT[1:3])
})

test_that("same-day games are excluded, including the test game itself", {
  DT <- .history(c("2024-01-01", "2024-01-02", "2024-01-02", "2024-01-03"))
  train <- walk_forward_training_set(DT, as.Date("2024-01-02"), min_train_share = 0)
  expect_equal(train, DT[1])
})

test_that("NULL when less than min_train_share of history precedes the test game", {
  DT <- .history(sprintf("2024-01-%02d", 1:10))
  expect_null(walk_forward_training_set(DT, as.Date("2024-01-05"), min_train_share = 0.5))  # 4 of 10
  expect_equal(walk_forward_training_set(DT, as.Date("2024-01-06"), min_train_share = 0.5), DT[1:5])
  expect_null(walk_forward_training_set(DT, as.Date("2024-01-02")))  # default share is 0.5
})

test_that("NULL when no game precedes the test game, even with share 0 or an empty DT", {
  DT <- .history(c("2024-01-05", "2024-01-06"))
  expect_null(walk_forward_training_set(DT, as.Date("2024-01-05"), min_train_share = 0))
  expect_null(walk_forward_training_set(DT[0], as.Date("2024-01-05")))
})

test_that("rows with an NA game_date are never used for training", {
  DT <- data.table(game_id = 1:3, game_date = as.Date(c("2024-01-01", NA, "2024-01-02")))
  train <- walk_forward_training_set(DT, as.Date("2024-01-05"), min_train_share = 0)
  expect_equal(train, DT[c(1, 3)])
})

test_that("works on a plain data.frame", {
  DF <- as.data.frame(.history(c("2024-01-01", "2024-01-02", "2024-01-03")))
  train <- walk_forward_training_set(DF, as.Date("2024-01-03"), min_train_share = 0)
  expect_equal(train, DF[1:2, ])
})

test_that("IDate and POSIXct game_date columns work like Date", {
  DT <- .history(c("2024-01-01", "2024-01-02", "2024-01-03"))
  as_idate <- copy(DT)[, game_date := as.IDate(game_date)]
  as_posix <- copy(DT)[, game_date := as.POSIXct(paste(game_date, "12:00:00"), tz = "UTC")]
  expect_equal(walk_forward_training_set(as_idate, as.IDate("2024-01-03"), min_train_share = 0), as_idate[1:2])
  expect_equal(walk_forward_training_set(as_posix, as_posix$game_date[3], min_train_share = 0), as_posix[1:2])
})

test_that("has_walk_forward_history agrees with walk_forward_training_set per test date", {
  DT <- .history(sprintf("2024-01-%02d", 1:10))
  test_dates <- as.Date(c("2024-01-01", "2024-01-05", "2024-01-06", "2024-02-01", NA))
  expect_equal(has_walk_forward_history(DT, test_dates), c(FALSE, FALSE, TRUE, TRUE, FALSE))
  for (k in which(!is.na(test_dates))) {
    expect_equal(is.null(walk_forward_training_set(DT, test_dates[k])),
                 !has_walk_forward_history(DT, test_dates[k]))
  }
})

test_that("fails loudly on a missing game_date column, character dates, or a bad test date", {
  expect_error(walk_forward_training_set(data.table(game_id = 1:3), as.Date("2024-01-01")),
               "expected a game_date column")
  char_dates <- data.table(game_id = 1:2, game_date = c("01/02/2024", "01/03/2024"))
  expect_error(walk_forward_training_set(char_dates, as.Date("2024-01-04")), "as Date or POSIXct")
  DT <- .history(c("2024-01-01", "2024-01-02"))
  expect_error(walk_forward_training_set(DT, "2024-01-03"), "as Date or POSIXct")
  expect_error(walk_forward_training_set(DT, NA), "expected one non-NA test_game_date")
  expect_error(walk_forward_training_set(DT, as.Date(c("2024-01-01", "2024-01-02"))),
               "expected one non-NA test_game_date")
})
