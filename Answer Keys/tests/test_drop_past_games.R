# Answer Keys/tests/test_drop_past_games.R
# Regression guard for Tools.R::.drop_past_games — the gate that keeps a
# scraper's stale snapshot (yesterday's slate, still in the DuckDB) out of
# today's odds screen. Every per-book table carries game_start_time as
# TIMESTAMPTZ UTC; the filter must key on that column directly.
#
# Background (2026-09-01): after the 2026-05-22 timezone standardization the
# filter still looked for the old game_date/game_time columns, warned
# "skipping past-game filter", and returned ALL rows — ~100 of 110 Wagerzon
# rows were a two-month-old slate.
#
# Run from "Answer Keys/" directory:
#   Rscript -e 'testthat::test_file("tests/test_drop_past_games.R")'

library(testthat)
source("../Tools.R")

.frame_with_start <- function(game_start_time) {
  data.frame(
    home_team       = sprintf("H%d", seq_along(game_start_time)),  # sprintf keeps length 0, paste0 would not
    away_team       = sprintf("A%d", seq_along(game_start_time)),
    fetch_time      = rep(Sys.time(), length(game_start_time)),  # fresh on every row
    game_start_time = game_start_time,
    stringsAsFactors = FALSE
  )
}

test_that("a POSIXct game_start_time in the past is dropped, a future one kept", {
  now <- Sys.time()
  raw <- .frame_with_start(as.POSIXct(c(now + 3600, now - 60 * 24 * 3600), tz = "UTC"))
  out <- .drop_past_games(raw, source_label = "test")
  expect_equal(nrow(out), 1)
  expect_equal(out$home_team, "H1")
  expect_equal(attr(out$game_start_time, "tzone"), "UTC")
})

test_that("5-minute grace keeps a game that started 4 minutes ago, drops 6 minutes ago", {
  now <- Sys.time()
  raw <- .frame_with_start(as.POSIXct(c(now - 4 * 60, now - 6 * 60), tz = "UTC"))
  out <- .drop_past_games(raw, source_label = "test")
  expect_equal(out$home_team, "H1")
})

test_that("a non-UTC POSIXct compares on the instant, not the wall clock", {
  # A game 1h in the future expressed in New York time must still be kept.
  future_ny <- as.POSIXct(Sys.time() + 3600, tz = "America/New_York")
  out <- .drop_past_games(.frame_with_start(future_ny), source_label = "test")
  expect_equal(nrow(out), 1)
  expect_equal(attr(out$game_start_time, "tzone"), "UTC")
})

test_that("a character ISO 8601 UTC game_start_time is parsed and filtered", {
  future_iso <- format(Sys.time() + 86400, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  raw <- .frame_with_start(c(future_iso, "2020-01-01T20:10:00Z"))
  out <- .drop_past_games(raw, source_label = "test")
  expect_equal(nrow(out), 1)
  expect_s3_class(out$game_start_time, "POSIXct")
})

test_that("NA game_start_time rows are dropped with a warning", {
  raw <- .frame_with_start(as.POSIXct(c(Sys.time() + 3600, NA), tz = "UTC"))
  expect_warning(out <- .drop_past_games(raw, source_label = "test"),
                 "1/2 rows have NA game_start_time")
  expect_equal(nrow(out), 1)
})

test_that("all-past frame warns that the book will be invisible and returns 0 rows", {
  raw <- .frame_with_start(as.POSIXct(Sys.time() - c(3600, 7200), tz = "UTC"))
  expect_warning(out <- .drop_past_games(raw, source_label = "test"),
                 "all 2 rows are past games")
  expect_equal(nrow(out), 0)
})

test_that("missing game_start_time column (pre-migration table) fails closed: warns, returns 0 rows", {
  raw <- data.frame(home_team = "H1", game_date = "06/28", game_time = "19:05",
                    stringsAsFactors = FALSE)
  expect_warning(out <- .drop_past_games(raw, source_label = "test"),
                 "no game_start_time column")
  expect_equal(nrow(out), 0)
  expect_equal(names(out), names(raw))
})

test_that("empty and NULL inputs pass through", {
  expect_null(.drop_past_games(NULL, source_label = "test"))
  empty <- .frame_with_start(as.POSIXct(character(), tz = "UTC"))
  expect_equal(nrow(.drop_past_games(empty, source_label = "test")), 0)
})
