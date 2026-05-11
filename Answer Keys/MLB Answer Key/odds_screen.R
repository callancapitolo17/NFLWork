# Answer Keys/MLB Answer Key/odds_screen.R
# Pure helpers for the bets-tab odds-screen pipeline path:
#   - normalize_book_odds_frame()    -- raw scraper rows -> canonical shape
#   - expand_bets_to_book_prices()   -- +/-1 unit nearest-line matching
#
# Sourced by MLB.R between the dedup step and the table write.
# Tested by tests/test_odds_screen.R.

library(dplyr)
library(tibble)
library(stringr)

LINE_MATCH_TOLERANCE <- 1.0  # max abs(line_quoted - model_line) we'll emit

#' Normalize a raw scraper / Odds API frame to the canonical shape consumed
#' by expand_bets_to_book_prices.
normalize_book_odds_frame <- function(raw) {
  if (nrow(raw) == 0) {
    return(tibble(game_id = character(), market = character(),
                  period = character(), side = character(),
                  line = numeric(), american_odds = integer(),
                  fetch_time = as.POSIXct(character())))
  }
  raw %>%
    mutate(
      period = case_when(
        str_detect(market_name, "_1st_3_innings$") ~ "F3",
        str_detect(market_name, "_1st_5_innings$") ~ "F5",
        str_detect(market_name, "_1st_7_innings$") ~ "F7",
        TRUE                                       ~ "FG"
      ),
      market = str_replace(market_name, "_1st_[357]_innings$", "")
    ) %>%
    transmute(
      game_id, market, period,
      side = bet_on,
      line, american_odds = as.integer(american_odds),
      fetch_time
    )
}

#' Pick the closest matching line from a candidate set of book quotes.
.pick_closest_line <- function(candidates, model_line, bet_on) {
  if (nrow(candidates) == 0) {
    return(tibble(line_quoted = numeric(), american_odds = integer(),
                  is_exact_line = logical()))
  }
  c2 <- candidates %>%
    mutate(.dist = abs(line - model_line)) %>%
    filter(.dist <= LINE_MATCH_TOLERANCE)
  if (nrow(c2) == 0) {
    return(tibble(line_quoted = numeric(), american_odds = integer(),
                  is_exact_line = logical()))
  }
  min_dist <- min(c2$.dist)
  c3 <- c2 %>% filter(.dist == min_dist)
  if (nrow(c3) > 1) {
    # Equidistant tiebreaker: prefer the line worse for the bettor.
    # Over side: pick higher line; Under side: pick lower line.
    # Spread favorite (negative line): pick lower (more negative); dog: pick higher.
    pick_high <- grepl("^(Over|over)", bet_on) || (model_line < 0)
    if (pick_high) c3 <- c3 %>% slice_max(line, n = 1)
    else c3 <- c3 %>% slice_min(line, n = 1)
  }
  row <- c3[1, , drop = FALSE]
  tibble(
    line_quoted   = row$line,
    american_odds = row$american_odds,
    is_exact_line = abs(row$line - model_line) < 1e-9
  )
}

#' Map a bet's bet_on string to the "pick" and "opposite" side labels used
#' in the per-book odds frames.
.sides_for_bet <- function(bet_on, market_type) {
  if (market_type == "totals") {
    list(pick = bet_on,
         opposite = if (grepl("^Over", bet_on)) "Under" else "Over")
  } else if (market_type %in% c("spreads", "moneyline")) {
    list(pick = bet_on, opposite = NA_character_)  # caller wires up opp
  } else {
    list(pick = bet_on, opposite = NA_character_)
  }
}

#' Main entry point.
expand_bets_to_book_prices <- function(bets, book_odds_by_book) {
  if (nrow(bets) == 0) {
    return(tibble(bet_row_id = character(), game_id = character(),
                  market = character(), period = character(),
                  side = character(), bookmaker = character(),
                  line = numeric(), line_quoted = numeric(),
                  is_exact_line = logical(), american_odds = integer(),
                  fetch_time = as.POSIXct(character())))
  }

  out_rows <- vector("list", length(book_odds_by_book) * nrow(bets) * 2)
  k <- 0L

  for (i in seq_len(nrow(bets))) {
    bet <- bets[i, , drop = FALSE]
    sides <- .sides_for_bet(bet$bet_on, bet$market_type)
    side_labels <- c(pick = sides$pick, opposite = sides$opposite)

    for (book_name in names(book_odds_by_book)) {
      book_frame <- book_odds_by_book[[book_name]]
      if (is.null(book_frame) || nrow(book_frame) == 0) next

      for (slot in names(side_labels)) {
        side_value <- side_labels[[slot]]
        if (is.na(side_value)) next

        candidates <- book_frame %>%
          filter(game_id == bet$game_id,
                 market  == bet$market_type,
                 period  == bet$period,
                 side    == side_value)

        chosen <- .pick_closest_line(candidates, bet$line, bet$bet_on)
        if (nrow(chosen) == 0) next

        ft <- candidates %>%
          filter(abs(line - chosen$line_quoted) < 1e-9) %>%
          slice_head(n = 1) %>%
          pull(fetch_time)

        k <- k + 1L
        out_rows[[k]] <- tibble(
          bet_row_id    = bet$bet_row_id,
          game_id       = bet$game_id,
          market        = bet$market_type,
          period        = bet$period,
          side          = slot,
          bookmaker     = book_name,
          line          = bet$line,
          line_quoted   = chosen$line_quoted,
          is_exact_line = chosen$is_exact_line,
          american_odds = as.integer(chosen$american_odds),
          fetch_time    = ft
        )
      }
    }
  }

  if (k == 0L) {
    return(tibble(bet_row_id = character(), game_id = character(),
                  market = character(), period = character(),
                  side = character(), bookmaker = character(),
                  line = numeric(), line_quoted = numeric(),
                  is_exact_line = logical(), american_odds = integer(),
                  fetch_time = as.POSIXct(character())))
  }
  bind_rows(out_rows[1:k])
}
