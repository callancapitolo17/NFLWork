# Answer Keys/tests/test_market_edge.R
# Unit tests for find_market_edges() — leave-one-out market-consensus edges.
# Run from "Answer Keys/":
#   Rscript -e 'testthat::test_file("tests/test_market_edge.R")'

library(testthat)
library(dplyr)
library(tibble)
source("../Tools.R")
source("../MLB Answer Key/odds_screen.R")
source("../MLB Answer Key/market_edge.R")

FT <- as.POSIXct("2026-06-05 18:00:00", tz = "UTC")
NOW <- as.POSIXct("2026-06-05 18:05:00", tz = "UTC")   # 5 min later

# canonical book row (matches scraper_to_canonical output shape)
crow <- function(game_id, market, period, side, line, american_odds, fetch_time = FT) {
  tibble(game_id = game_id, market = market, period = period, side = side,
         line = line, american_odds = as.integer(american_odds), fetch_time = fetch_time)
}

# A totals market at one book = Over + Under at the same line.
tot <- function(game_id, book_side_odds, line = 8.5, ft = FT) {
  bind_rows(
    crow(game_id, "totals", "F5", "Over",  line, book_side_odds[1], ft),
    crow(game_id, "totals", "F5", "Under", line, book_side_odds[2], ft)
  )
}

# A spread market at one book = the two teams at +L / -L (signed per team).
spr <- function(game_id, home_team, away_team, home_odds, away_odds,
                home_line = -1.5, ft = FT) {
  bind_rows(
    crow(game_id, "spreads", "F5", home_team,  home_line, home_odds, ft),
    crow(game_id, "spreads", "F5", away_team, -home_line, away_odds, ft)
  )
}

# A moneyline (h2h) market at one book = the two teams, line NA.
ml <- function(game_id, home_team, away_team, home_odds, away_odds, ft = FT) {
  bind_rows(
    crow(game_id, "h2h", "F5", home_team, NA_real_, home_odds, ft),
    crow(game_id, "h2h", "F5", away_team, NA_real_, away_odds, ft)
  )
}

test_that("a soft book beating the others flags as a market edge", {
  books <- list(
    pinnacle = tot("g1", c(-110, -110)),
    bfa      = tot("g1", c(-110, -110)),
    wagerzon = tot("g1", c(+110, -130))
  )
  out <- find_market_edges(books, now = NOW)
  over <- out %>% filter(bet_on == "Over")
  expect_equal(nrow(over), 1)
  expect_equal(over$bookmaker_key, "wagerzon")
  expect_equal(over$edge_source, "market")
  expect_true(over$ev >= 0.02)
  expect_equal(over$market_type, "totals")
  expect_equal(over$period, "F5")
  expect_equal(over$line, 8.5)
  expect_equal(over$n_books, 3L)
})

test_that("prices in line with the market produce no flag", {
  books <- list(
    pinnacle = tot("g1", c(-110, -110)),
    bfa      = tot("g1", c(-108, -112)),
    wagerzon = tot("g1", c(-109, -111))
  )
  out <- find_market_edges(books, now = NOW)
  expect_equal(nrow(out), 0)
})

test_that("two books (edge + 1 other) is enough to flag", {
  books <- list(
    pinnacle = tot("g1", c(-110, -110)),
    wagerzon = tot("g1", c(+110, -130))
  )
  out <- find_market_edges(books, now = NOW)
  expect_true(any(out$bet_on == "Over" & out$bookmaker_key == "wagerzon"))
})

test_that("a single book (no other to compare) is skipped", {
  books <- list(wagerzon = tot("g1", c(+150, -180)))
  out <- find_market_edges(books, now = NOW)
  expect_equal(nrow(out), 0)
})

test_that("a one-sided book is excluded from the consensus", {
  books <- list(
    pinnacle = tot("g1", c(-110, -110)),
    bfa      = tot("g1", c(-108, -112)),
    wagerzon = crow("g1", "totals", "F5", "Over", 8.5, +104)
  )
  out <- find_market_edges(books, now = NOW)
  expect_false(any(out$bookmaker_key == "wagerzon"))
})

test_that("stale quotes are dropped before devigging", {
  stale_ft <- as.POSIXct("2026-06-05 17:00:00", tz = "UTC")  # 65 min before NOW
  books <- list(
    pinnacle = tot("g1", c(-110, -110)),
    wagerzon = tot("g1", c(+104, -124), ft = stale_ft)
  )
  out <- find_market_edges(books, now = NOW, staleness_min = 30)
  expect_equal(nrow(out), 0)
})

test_that("bet_row_id matches the shared helper for the same wager", {
  books <- list(
    pinnacle = tot("g1", c(-110, -110)),
    wagerzon = tot("g1", c(+110, -130))
  )
  out <- find_market_edges(books, now = NOW) %>% filter(bet_on == "Over")
  expected <- compute_bet_row_id("g1", "totals", "F5", 8.5, "Over")
  expect_equal(out$bet_row_id, expected)
})

test_that("a soft book flags on a SPREAD (signed-line pairing)", {
  books <- list(
    pinnacle = spr("g1", "NYY", "BOS", -110, -110),
    wagerzon = spr("g1", "NYY", "BOS", +110, -130)
  )
  out <- find_market_edges(books, now = NOW)
  hit <- out %>% filter(market_type == "spreads", bet_on == "NYY")
  expect_equal(nrow(hit), 1)
  expect_equal(hit$bookmaker_key, "wagerzon")
  expect_equal(hit$line, -1.5)
  expect_true(hit$ev >= 0.02)
})

test_that("a soft book flags on a MONEYLINE (NA-line grouping)", {
  books <- list(
    pinnacle = ml("g1", "NYY", "BOS", -110, -110),
    wagerzon = ml("g1", "NYY", "BOS", +110, -130)
  )
  out <- find_market_edges(books, now = NOW)
  hit <- out %>% filter(market_type == "h2h", bet_on == "NYY")
  expect_equal(nrow(hit), 1)
  expect_equal(hit$bookmaker_key, "wagerzon")
  expect_true(is.na(hit$line))
  expect_true(hit$ev >= 0.02)
})
