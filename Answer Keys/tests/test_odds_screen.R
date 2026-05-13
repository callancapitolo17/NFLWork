# Answer Keys/tests/test_odds_screen.R
# Tests for expand_bets_to_book_prices — long-format expansion of bets
# into per-book pill rows with +/-1 unit nearest-line matching.
#
# Run from "Answer Keys/" directory:
#   Rscript -e 'testthat::test_file("tests/test_odds_screen.R")'

library(testthat)
library(dplyr)
library(tibble)
library(jsonlite)
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

test_that("lines outside the display tolerance emit no row", {
  # LINE_MATCH_TOLERANCE = 3.0 (raised from 1.0 on 2026-05-13 to expose more
  # book quotes with amber line tags). Anything > 3.0 units off should still
  # be dropped — a bet at 5.5 with a book line at 1.5 (4 units off) qualifies.
  bets <- make_bet_row(line = 5.5)
  book_odds <- list(bfa = bind_rows(
    book_row("g1", "totals", "F5", "Over", 1.5, -115)  # 4 units off
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

test_that("empty book_odds_by_book list returns empty schema tibble", {
  bets <- make_bet_row()
  out <- expand_bets_to_book_prices(bets, list())
  expect_equal(nrow(out), 0)
  # Verify columns exist even when empty
  expect_true(all(c("bet_row_id", "game_id", "market", "period", "side",
                    "bookmaker", "line", "line_quoted", "is_exact_line",
                    "american_odds", "fetch_time") %in% names(out)))
})

test_that("equidistant tiebreak prefers the line worse for the Under bettor", {
  bets <- make_bet_row(line = 5.5, bet_on = "Under")
  book_odds <- list(bfa = bind_rows(
    book_row("g1", "totals", "F5", "Over", 5.0, +110),
    book_row("g1", "totals", "F5", "Over", 6.0, -130),
    book_row("g1", "totals", "F5", "Under", 5.0, -130),
    book_row("g1", "totals", "F5", "Under", 6.0, +110)
  ))
  out <- expand_bets_to_book_prices(bets, book_odds) %>%
    filter(side == "pick")
  expect_equal(out$line_quoted, 5.0)  # lower line is harder for Under bettor
})

test_that("normalize_book_odds_frame returns empty-schema tibble on empty input", {
  raw <- tibble(game_id = character(), market_name = character(),
                bet_on = character(), line = numeric(),
                american_odds = integer(),
                fetch_time = as.POSIXct(character()))
  out <- normalize_book_odds_frame(raw)
  expect_equal(nrow(out), 0)
  expect_true(all(c("game_id", "market", "period", "side", "line",
                    "american_odds", "fetch_time") %in% names(out)))
})

test_that("expand_bets_to_book_prices accepts 'id' as an alias for 'game_id'", {
  bets <- make_bet_row() %>%
    rename(game_id_alt_renamed_via_id = game_id) %>%
    rename(id = game_id_alt_renamed_via_id)
  # Drop the period column to also exercise the period-derivation fallback
  bets <- bets %>% select(-period)
  book_odds <- list(wagerzon = bind_rows(
    book_row("g1", "totals", "F5", "Over",  5.5, +120),
    book_row("g1", "totals", "F5", "Under", 5.5, -140)
  ))
  out <- expand_bets_to_book_prices(bets, book_odds)
  expect_equal(nrow(out), 2)
  expect_setequal(out$side, c("pick", "opposite"))
})

# =============================================================================
# parse_prefetched_to_long tests
# =============================================================================

test_that("parse_prefetched_to_long emits one row per outcome", {
  fake_json <- toJSON(list(
    id = "evt1",
    bookmakers = list(
      list(
        key = "draftkings",
        last_update = "2026-05-11T12:00:00Z",
        markets = list(
          list(
            key = "totals_1st_5_innings",
            last_update = "2026-05-11T12:00:00Z",
            outcomes = list(
              list(name = "Over", price = -110, point = 5.5),
              list(name = "Under", price = -110, point = 5.5)
            )
          )
        )
      )
    )
  ), auto_unbox = TRUE, pretty = FALSE)
  prefetched <- tibble(event_id = "evt1", json_response = fake_json)
  out <- parse_prefetched_to_long(prefetched, bookmaker_keys = c("draftkings"))
  expect_equal(nrow(out), 2)
  expect_setequal(out$outcomes_name, c("Over", "Under"))
  expect_equal(unique(out$bookmaker_key), "draftkings")
  expect_equal(unique(out$market_key), "totals_1st_5_innings")
  expect_false(any(is.na(out$fetch_time)))
})

test_that("parse_prefetched_to_long filters out non-listed bookmakers", {
  fake_json <- toJSON(list(
    id = "evt1",
    bookmakers = list(
      list(key = "draftkings", last_update = "2026-05-11T12:00:00Z",
           markets = list(list(key = "h2h", last_update = "2026-05-11T12:00:00Z",
             outcomes = list(list(name = "Home", price = -110, point = NULL),
                             list(name = "Away", price = -110, point = NULL))))),
      list(key = "betmgm", last_update = "2026-05-11T12:00:00Z",
           markets = list(list(key = "h2h", last_update = "2026-05-11T12:00:00Z",
             outcomes = list(list(name = "Home", price = -105, point = NULL),
                             list(name = "Away", price = -115, point = NULL)))))
    )
  ), auto_unbox = TRUE, null = "null")
  prefetched <- tibble(event_id = "evt1", json_response = fake_json)
  out <- parse_prefetched_to_long(prefetched, bookmaker_keys = c("draftkings", "fanduel"))
  expect_equal(nrow(out), 2)
  expect_true(all(out$bookmaker_key == "draftkings"))
})

test_that("parse_prefetched_to_long emits warning on malformed JSON", {
  prefetched <- tibble(event_id = "evt1", json_response = "{ malformed json")
  expect_warning(
    out <- parse_prefetched_to_long(prefetched, bookmaker_keys = c("draftkings")),
    "JSON parse failed"
  )
  expect_equal(nrow(out), 0)
})

test_that("parse_prefetched_to_long handles NA point on moneyline", {
  fake_json <- toJSON(list(
    id = "evt1",
    bookmakers = list(list(
      key = "draftkings", last_update = "2026-05-11T12:00:00Z",
      markets = list(list(
        key = "h2h", last_update = "2026-05-11T12:00:00Z",
        outcomes = list(list(name = "Home", price = -110, point = NULL),
                        list(name = "Away", price = -110, point = NULL))
      ))
    ))
  ), auto_unbox = TRUE, null = "null")
  prefetched <- tibble(event_id = "evt1", json_response = fake_json)
  out <- parse_prefetched_to_long(prefetched, bookmaker_keys = c("draftkings"))
  expect_equal(nrow(out), 2)
  expect_true(all(is.na(out$outcomes_point)))
})

test_that("parse_prefetched_to_long returns empty schema on empty input", {
  out <- parse_prefetched_to_long(
    tibble(event_id = character(), json_response = character()),
    bookmaker_keys = c("draftkings"))
  expect_equal(nrow(out), 0)
  expect_true(all(c("id", "bookmaker_key", "market_key", "outcomes_name",
                    "outcomes_price", "outcomes_point", "fetch_time")
                  %in% names(out)))
})

test_that("expand_bets_to_book_prices emits opposite row for spread bets", {
  bets <- tibble(
    bet_row_id   = "row1",
    id           = "g1",
    home_team    = "New York Yankees",
    away_team    = "Boston Red Sox",
    market       = "spreads",
    market_type  = "spreads",
    period       = "FG",
    line         = -1.5,
    bet_on       = "New York Yankees"
  )

  # WZ has -1.5 on home and +1.5 on away.
  wz <- tibble(
    game_id       = "g1",
    market        = "spreads",
    period        = "FG",
    side          = c("New York Yankees", "Boston Red Sox"),
    line          = c(-1.5, 1.5),
    american_odds = c(-120L, 100L),
    fetch_time    = as.POSIXct("2026-05-13 12:00", tz = "UTC")
  )

  out <- expand_bets_to_book_prices(bets, list(wz = wz))

  expect_equal(nrow(out), 2)  # one pick row + one opposite row
  expect_setequal(out$side, c("pick", "opposite"))

  pick <- out[out$side == "pick", ]
  opp  <- out[out$side == "opposite", ]
  expect_equal(pick$american_odds, -120L)
  expect_equal(pick$line_quoted, -1.5)
  expect_equal(opp$american_odds,  100L)
  expect_equal(opp$line_quoted,    1.5)
})

test_that("expand_bets_to_book_prices emits opposite row for moneyline bets", {
  bets <- tibble(
    bet_row_id   = "row2",
    id           = "g2",
    home_team    = "New York Yankees",
    away_team    = "Boston Red Sox",
    market       = "h2h",
    market_type  = "h2h",
    period       = "FG",
    line         = NA_real_,
    bet_on       = "Boston Red Sox"
  )

  wz <- tibble(
    game_id       = "g2",
    market        = "h2h",
    period        = "FG",
    side          = c("Boston Red Sox", "New York Yankees"),
    line          = c(NA_real_, NA_real_),
    american_odds = c(150L, -170L),
    fetch_time    = as.POSIXct("2026-05-13 12:00", tz = "UTC")
  )

  out <- expand_bets_to_book_prices(bets, list(wz = wz))
  expect_equal(nrow(out), 2)
  expect_setequal(out$side, c("pick", "opposite"))
  expect_equal(out[out$side == "pick", ]$american_odds,  150L)
  expect_equal(out[out$side == "opposite", ]$american_odds, -170L)
})
