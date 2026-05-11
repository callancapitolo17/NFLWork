# Answer Keys/tests/test_odds_screen.R
# Tests for expand_bets_to_book_prices — long-format expansion of bets
# into per-book pill rows with +/-1 unit nearest-line matching.
#
# Run from "Answer Keys/" directory:
#   Rscript -e 'testthat::test_file("tests/test_odds_screen.R")'

library(testthat)
library(dplyr)
library(tibble)
source("../MLB Answer Key/odds_screen.R")

# A canonical 1-bet input frame matching the columns mlb_bets_combined
# carries: id, market, line, bet_on, market_type, plus the helper's
# computed bet_row_id (the caller hashes it in).
make_bet_row <- function(game_id = "g1",
                         market = "totals_1st_5_innings",
                         line = 5.5,
                         bet_on = "Over",
                         market_type = "totals") {
  tibble(
    bet_row_id  = "fakehash",
    game_id     = game_id,
    market      = market,
    market_type = market_type,
    period      = "F5",
    line        = line,
    bet_on      = bet_on,
    pick_side   = "pick"
  )
}

# Per-book odds frame matching the shape scrapers + game_odds expose.
make_book_odds <- function(rows) {
  bind_rows(rows)
}

book_row <- function(game_id, market, period, side, line, american_odds,
                     fetch_time = as.POSIXct("2026-05-11 12:00:00", tz = "UTC")) {
  tibble(
    game_id = game_id, market = market, period = period, side = side,
    line = line, american_odds = american_odds, fetch_time = fetch_time
  )
}

test_that("exact-line quote at one book emits two rows (pick + opposite)", {
  bets <- make_bet_row()
  book_odds <- list(wagerzon = bind_rows(
    book_row("g1", "totals", "F5", "Over",  5.5, +120),
    book_row("g1", "totals", "F5", "Under", 5.5, -140)
  ))
  out <- expand_bets_to_book_prices(bets, book_odds)
  expect_equal(nrow(out), 2)
  expect_setequal(out$side, c("pick", "opposite"))
  expect_true(all(out$is_exact_line))
  expect_equal(out$line_quoted[out$side == "pick"], 5.5)
  expect_equal(out$american_odds[out$side == "pick"], 120)
})

test_that("missing book emits zero rows for that book", {
  bets <- make_bet_row()
  book_odds <- list(wagerzon = bind_rows(
    book_row("g1", "totals", "F5", "Over", 5.5, +120)
  ))
  out <- expand_bets_to_book_prices(bets, book_odds)
  expect_false("hoop88" %in% out$bookmaker)
})

test_that("nearest-line within +/-1 unit is taken with is_exact_line=FALSE", {
  bets <- make_bet_row(line = 5.5)
  book_odds <- list(bfa = bind_rows(
    book_row("g1", "totals", "F5", "Over",  5.0, -115),
    book_row("g1", "totals", "F5", "Under", 5.0, -105)
  ))
  out <- expand_bets_to_book_prices(bets, book_odds)
  expect_equal(nrow(out), 2)
  expect_false(any(out$is_exact_line))
  expect_equal(unique(out$line_quoted), 5.0)
})

test_that("lines outside +/-1 unit emit no row", {
  bets <- make_bet_row(line = 5.5)
  book_odds <- list(bfa = bind_rows(
    book_row("g1", "totals", "F5", "Over", 3.5, -115)  # 2 units off
  ))
  out <- expand_bets_to_book_prices(bets, book_odds)
  expect_equal(nrow(out), 0)
})

test_that("when both 0.5-unit and 1.0-unit lines exist, the closer one wins", {
  bets <- make_bet_row(line = 5.5)
  book_odds <- list(bfa = bind_rows(
    book_row("g1", "totals", "F5", "Over", 5.0, -115),
    book_row("g1", "totals", "F5", "Over", 4.5, -150),
    book_row("g1", "totals", "F5", "Under", 5.0, -105),
    book_row("g1", "totals", "F5", "Under", 4.5, +130)
  ))
  out <- expand_bets_to_book_prices(bets, book_odds)
  expect_equal(unique(out$line_quoted), 5.0)
})

test_that("equidistant tiebreak prefers the line worse for the bettor (over)", {
  bets <- make_bet_row(line = 5.5, bet_on = "Over")
  book_odds <- list(bfa = bind_rows(
    book_row("g1", "totals", "F5", "Over", 5.0, +110),
    book_row("g1", "totals", "F5", "Over", 6.0, -130),
    book_row("g1", "totals", "F5", "Under", 5.0, -130),
    book_row("g1", "totals", "F5", "Under", 6.0, +110)
  ))
  out <- expand_bets_to_book_prices(bets, book_odds) %>%
    filter(side == "pick")
  expect_equal(out$line_quoted, 6.0)
})

test_that("bet_row_id is preserved on every emitted row", {
  bets <- make_bet_row()
  book_odds <- list(wagerzon = bind_rows(
    book_row("g1", "totals", "F5", "Over", 5.5, +120),
    book_row("g1", "totals", "F5", "Under", 5.5, -140)
  ))
  out <- expand_bets_to_book_prices(bets, book_odds)
  expect_true(all(out$bet_row_id == "fakehash"))
})

test_that("multiple books on the same bet expand independently", {
  bets <- make_bet_row()
  book_odds <- list(
    wagerzon = bind_rows(
      book_row("g1", "totals", "F5", "Over",  5.5, +120),
      book_row("g1", "totals", "F5", "Under", 5.5, -140)
    ),
    bfa = bind_rows(
      book_row("g1", "totals", "F5", "Over",  5.0, -115),
      book_row("g1", "totals", "F5", "Under", 5.0, -105)
    )
  )
  out <- expand_bets_to_book_prices(bets, book_odds)
  expect_setequal(out$bookmaker, c("wagerzon", "bfa"))
  expect_equal(nrow(out), 4)
})

test_that("normalize_book_odds_frame splits market into (market, period)", {
  raw <- tibble(
    game_id = "g1",
    market_name = c("totals_1st_5_innings", "spreads_1st_3_innings",
                    "h2h", "alternate_totals"),
    bet_on = c("Over", "Home", "Home", "Over"),
    line = c(5.5, -1.5, NA, 8.5),
    american_odds = c(-110L, +120L, -150L, +105L),
    fetch_time = rep(as.POSIXct("2026-05-11 12:00", tz="UTC"), 4)
  )
  out <- normalize_book_odds_frame(raw)
  expect_setequal(out$market, c("totals", "spreads", "h2h", "alternate_totals"))
  expect_setequal(out$period, c("F5", "F3", "FG", "FG"))
  expect_equal(out$side, c("Over", "Home", "Home", "Over"))
})

test_that("normalize_book_odds_frame handles F7 and alts", {
  raw <- tibble(
    game_id = "g1",
    market_name = c("totals_1st_7_innings", "alternate_spreads_1st_5_innings"),
    bet_on = c("Under", "Away"),
    line = c(6.5, 1.5),
    american_odds = c(-115L, -120L),
    fetch_time = rep(as.POSIXct("2026-05-11 12:00", tz="UTC"), 2)
  )
  out <- normalize_book_odds_frame(raw)
  expect_equal(out$market, c("totals", "alternate_spreads"))
  expect_equal(out$period, c("F7", "F5"))
})
