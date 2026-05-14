# Answer Keys/MLB Answer Key/odds_screen.R
# Pure helpers for the bets-tab odds-screen pipeline path:
#   - normalize_book_odds_frame()        -- raw scraper rows -> canonical shape
#   - expand_bets_to_book_prices()       -- +/-1 unit nearest-line matching
#   - scraper_to_canonical()             -- wide scraper frame -> canonical
#   - odds_api_to_canonical()            -- Odds API long frame -> canonical
#   - parse_prefetched_to_long()         -- raw JSON cache -> long frame
#
# Sourced by MLB.R between the dedup step and the table write.
# Tested by tests/test_odds_screen.R.

library(dplyr)
library(tibble)
library(stringr)
library(lubridate)
library(jsonlite)

LINE_MATCH_TOLERANCE <- 3.0  # max abs(line_quoted - model_line) we'll emit
# Raised from 1.0 -> 3.0 on 2026-05-13. With 1.0, F7 totals where the
# model bets Over 5.5 but DK posts 7.5 (a 2.0 spread normal for F7)
# showed no DK pill at all. book_pill.R already renders the pill amber
# with a line tag like "O7.5" when is_exact_line is FALSE, so wider
# tolerance exposes more book quotes with their actual line clearly
# marked. The pick book is still required to be on the exact line per
# the data-bug guard in MLB.R, so this affects display only.

# Internal helpers for period/market-type derivation.
# Extracted here so normalize_book_odds_frame and expand_bets_to_book_prices
# share the same logic without duplication.

.derive_period <- function(market_name) {
  case_when(
    str_detect(market_name, "_1st_3_innings$") ~ "F3",
    str_detect(market_name, "_1st_5_innings$") ~ "F5",
    str_detect(market_name, "_1st_7_innings$") ~ "F7",
    # Alt-market convention used by compare_alts_to_samples: <market>_<period>
    # (e.g. alternate_totals_fg, alternate_spreads_f5). Recognize all four.
    str_detect(market_name, "_f3$")            ~ "F3",
    str_detect(market_name, "_f5$")            ~ "F5",
    str_detect(market_name, "_f7$")            ~ "F7",
    str_detect(market_name, "_fg$")            ~ "FG",
    TRUE                                       ~ "FG"
  )
}

.derive_market_type <- function(market_name) {
  market_name %>%
    str_replace("_1st_[357]_innings$", "") %>%
    str_replace("_(fg|f3|f5|f7)$", "")
}

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
      period = .derive_period(market_name),
      market = .derive_market_type(market_name)
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
  # Moneyline / h2h have no line on either side — model_line is NA and the
  # book quotes have NA line. abs(NA - NA) is NA, which dplyr::filter would
  # drop, so handle this case explicitly: any NA-line candidate is an
  # "exact" match.
  if (is.na(model_line)) {
    c_ml <- candidates %>% filter(is.na(line))
    if (nrow(c_ml) == 0) {
      return(tibble(line_quoted = numeric(), american_odds = integer(),
                    is_exact_line = logical()))
    }
    row <- c_ml[1, , drop = FALSE]
    return(tibble(
      line_quoted   = NA_real_,
      american_odds = row$american_odds,
      is_exact_line = TRUE
    ))
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
    pick_high <- grepl("^Over", bet_on, ignore.case = TRUE) || (model_line < 0)
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
#'
#' For totals, opposite is derived automatically (Over <-> Under).
#' For spreads, moneyline, and h2h, the opposite side is the OTHER team —
#' resolved from home_team / away_team if present, else from the explicit
#' opposite_side column, else NA (no opposite-side row emitted).
.sides_for_bet <- function(bet_on, market_type,
                            opposite_side = NA_character_,
                            home_team = NA_character_,
                            away_team = NA_character_) {
  if (market_type == "totals" || grepl("^alternate_totals", market_type) ||
      grepl("^totals_1st_[357]_innings$", market_type)) {
    list(pick = bet_on,
         opposite = if (grepl("^Over", bet_on, ignore.case = TRUE)) "Under" else "Over")
  } else {
    # Spreads / alt_spreads / h2h / moneyline: opposite is the other team.
    derived_opp <- if (!is.na(home_team) && !is.na(away_team)) {
      if (identical(bet_on, home_team)) away_team
      else if (identical(bet_on, away_team)) home_team
      else NA_character_
    } else NA_character_
    final_opp <- if (!is.na(derived_opp)) derived_opp else opposite_side
    list(pick = bet_on, opposite = final_opp)
  }
}

#' Main entry point.
#'
#' @param bets         A data frame of bets. Must contain: bet_row_id, game_id
#'                     (or `id` as an alias), market_type, period (or `market`
#'                     as a fallback to derive it from), line, bet_on.
#'                     Optional: opposite_side (required for spreads/moneyline
#'                     to emit the opposite-side row).
#' @param book_odds_by_book Named list of normalized book frames (output of
#'                     normalize_book_odds_frame). Each frame may also use `id`
#'                     instead of `game_id`; both are accepted.
expand_bets_to_book_prices <- function(bets, book_odds_by_book) {
  if (nrow(bets) == 0) {
    return(tibble(bet_row_id = character(), game_id = character(),
                  market = character(), period = character(),
                  side = character(), bookmaker = character(),
                  line = numeric(), line_quoted = numeric(),
                  is_exact_line = logical(), american_odds = integer(),
                  fetch_time = as.POSIXct(character())))
  }

  # --- Defensive column renames for bets frame ---
  # all_bets_combined uses `id` (Odds API game id), not `game_id`.
  if (!"game_id" %in% names(bets) && "id" %in% names(bets)) {
    bets <- bets %>% rename(game_id = id)
  }

  # all_bets_combined has `market` (e.g. "totals_1st_5_innings") but no
  # separate `period` or `market_type` columns. Derive them if absent.
  if (!"period" %in% names(bets)) {
    bets <- bets %>% mutate(period = .derive_period(market))
  }
  if (!"market_type" %in% names(bets)) {
    bets <- bets %>% mutate(market_type = .derive_market_type(market))
  }

  # --- Defensive column renames for each book frame ---
  book_odds_by_book <- lapply(book_odds_by_book, function(b) {
    if (!is.null(b) && !"game_id" %in% names(b) && "id" %in% names(b)) {
      b %>% rename(game_id = id)
    } else b
  })

  out_rows <- vector("list", length(book_odds_by_book) * nrow(bets) * 2)
  k <- 0L

  for (i in seq_len(nrow(bets))) {
    bet <- bets[i, , drop = FALSE]
    opposite_col_value <- if ("opposite_side" %in% names(bet)) bet$opposite_side else NA_character_
    home_team_value    <- if ("home_team"     %in% names(bet)) bet$home_team     else NA_character_
    away_team_value    <- if ("away_team"     %in% names(bet)) bet$away_team     else NA_character_
    sides <- .sides_for_bet(bet$bet_on, bet$market_type, opposite_col_value,
                            home_team = home_team_value, away_team = away_team_value)
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

        # Spread/alt_spread: away-side candidates carry the negated line
        # (see scraper_to_canonical()). The opposite slot must compare
        # against -bet$line, not bet$line, or every same-line book gets
        # flagged as a line mismatch. Totals share one line across both
        # sides; moneyline has line=NA.
        is_spread_market <- bet$market_type %in% c("spreads", "alternate_spreads")
        effective_line <- if (slot == "opposite" && is_spread_market && !is.na(bet$line)) {
          -bet$line
        } else {
          bet$line
        }
        chosen <- .pick_closest_line(candidates, effective_line, bet$bet_on)
        if (nrow(chosen) == 0) next

        # Find the fetch_time for the chosen quote. Moneyline / h2h have
        # NA lines on both sides; match on is.na() in that case.
        ft <- if (is.na(chosen$line_quoted)) {
          candidates %>% filter(is.na(line)) %>%
            slice_head(n = 1) %>% pull(fetch_time)
        } else {
          candidates %>%
            filter(abs(line - chosen$line_quoted) < 1e-9) %>%
            slice_head(n = 1) %>% pull(fetch_time)
        }

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

# =============================================================================
# PREFETCH PARSER & CANONICAL CONVERTERS
# (moved from MLB.R so they can be sourced and tested independently)
# =============================================================================

# Small NULL-coalesce helper (matches rlang's %||% behaviour without the dep).
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Parse the prefetched_odds frame (event_id + raw JSON from the Odds API
#' /events/{id}/odds endpoint) into the long format expected by
#' odds_api_to_canonical().
#'
#' @param prefetched_odds tibble with columns event_id, json_response (char)
#' @param bookmaker_keys  character vector of book keys to retain
#' @return tibble: id, bookmaker_key, market_key, outcomes_name,
#'         outcomes_price, outcomes_point, fetch_time
parse_prefetched_to_long <- function(prefetched_odds, bookmaker_keys) {
  if (is.null(prefetched_odds) || nrow(prefetched_odds) == 0) {
    return(tibble(
      id = character(), bookmaker_key = character(), market_key = character(),
      outcomes_name = character(), outcomes_price = integer(),
      outcomes_point = numeric(), fetch_time = as.POSIXct(character(), tz = "UTC")
    ))
  }

  out <- list()
  for (i in seq_len(nrow(prefetched_odds))) {
    js <- prefetched_odds$json_response[i]
    if (is.na(js) || !nzchar(js)) next
    parsed <- tryCatch(
      jsonlite::fromJSON(js, simplifyVector = FALSE),
      error = function(e) {
        warning(sprintf(
          "[parse_prefetched] JSON parse failed for event %s: %s",
          prefetched_odds$event_id[i], conditionMessage(e)))
        NULL
      }
    )
    if (is.null(parsed) || is.null(parsed$bookmakers)) next
    eid <- parsed$id %||% prefetched_odds$event_id[i]

    for (bm in parsed$bookmakers) {
      bk <- bm$key
      if (is.null(bk) || !(bk %in% bookmaker_keys)) next
      bm_lu <- bm$last_update  # ISO8601 string, or NULL

      for (mk in bm$markets %||% list()) {
        mkey <- mk$key
        if (is.null(mkey)) next
        # Prefer market-level last_update (more granular); fall back to bookmaker-level.
        ft_str <- mk$last_update %||% bm_lu
        ft_val <- if (!is.null(ft_str)) {
          tryCatch(ymd_hms(ft_str, tz = "UTC", quiet = TRUE),
                   error = function(e) Sys.time())
        } else Sys.time()
        if (is.na(ft_val)) ft_val <- Sys.time()

        for (oc in mk$outcomes %||% list()) {
          out[[length(out) + 1]] <- tibble(
            id             = as.character(eid),
            bookmaker_key  = bk,
            market_key     = mkey,
            outcomes_name  = oc$name %||% NA_character_,
            outcomes_price = if (!is.null(oc$price)) as.integer(oc$price) else NA_integer_,
            outcomes_point = if (!is.null(oc$point)) as.numeric(oc$point) else NA_real_,
            fetch_time     = ft_val
          )
        }
      }
    }
  }

  if (length(out) == 0) {
    return(tibble(
      id = character(), bookmaker_key = character(), market_key = character(),
      outcomes_name = character(), outcomes_price = integer(),
      outcomes_point = numeric(), fetch_time = as.POSIXct(character(), tz = "UTC")
    ))
  }
  bind_rows(out)
}

#' Convert a wide-format scraper frame to canonical shape for
#' normalize_book_odds_frame(). Scraper frames have odds_home/odds_away/
#' odds_over/odds_under as separate columns; we pivot to long.
#'
#' @param raw    Wide scraper tibble (wagerzon/hoop88/bfa/bookmaker/bet105)
#' @param lookup tibble with columns id, home_team, away_team — used to join
#'               scraper team names to Odds API game IDs. Typically derived
#'               from mlb_odds in MLB.R.
#' @return normalized tibble (output of normalize_book_odds_frame), or NULL
scraper_to_canonical <- function(raw, lookup) {
  if (is.null(raw) || nrow(raw) == 0) return(NULL)
  if (!"fetch_time" %in% names(raw)) raw$fetch_time <- as.POSIXct(NA, tz = "UTC")

  # Replace NA fetch_times with now for consistency with odds_api_to_canonical.
  if (!is.null(raw$fetch_time)) {
    raw$fetch_time[is.na(raw$fetch_time)] <- Sys.time()
  }

  # Join to Odds API game IDs. Scraper frames have home_team + away_team
  # already resolved to canonical names by resolve_offshore_teams().
  joined <- raw %>%
    inner_join(lookup, by = c("home_team", "away_team")) %>%
    rename(game_id = id)

  if (nrow(joined) == 0) return(NULL)

  rows <- list()

  for (i in seq_len(nrow(joined))) {
    row <- joined[i, , drop = FALSE]
    mkt <- row$market
    ft  <- row$fetch_time
    gid <- row$game_id

    # Spread side (home / away)
    if (!is.na(row$odds_home) && !is.na(row$line)) {
      rows[[length(rows) + 1]] <- tibble(
        game_id = gid, market_name = mkt,
        bet_on = row$home_team, line = as.numeric(row$line),
        american_odds = as.integer(row$odds_home), fetch_time = ft
      )
    }
    if (!is.na(row$odds_away) && !is.na(row$line)) {
      rows[[length(rows) + 1]] <- tibble(
        game_id = gid, market_name = mkt,
        bet_on = row$away_team,
        line = -as.numeric(row$line),   # away spread is negated home spread
        american_odds = as.integer(row$odds_away), fetch_time = ft
      )
    }
    # Totals (Over / Under)
    if (!is.na(row$odds_over) && !is.na(row$line)) {
      rows[[length(rows) + 1]] <- tibble(
        game_id = gid, market_name = mkt,
        bet_on = "Over", line = as.numeric(row$line),
        american_odds = as.integer(row$odds_over), fetch_time = ft
      )
    }
    if (!is.na(row$odds_under) && !is.na(row$line)) {
      rows[[length(rows) + 1]] <- tibble(
        game_id = gid, market_name = mkt,
        bet_on = "Under", line = as.numeric(row$line),
        american_odds = as.integer(row$odds_under), fetch_time = ft
      )
    }
    # Moneyline (h2h) — odds_home/odds_away with no line.
    if (!is.na(row$odds_home) && is.na(row$line)) {
      rows[[length(rows) + 1]] <- tibble(
        game_id = gid, market_name = mkt,
        bet_on = row$home_team, line = NA_real_,
        american_odds = as.integer(row$odds_home), fetch_time = ft
      )
    }
    if (!is.na(row$odds_away) && is.na(row$line)) {
      rows[[length(rows) + 1]] <- tibble(
        game_id = gid, market_name = mkt,
        bet_on = row$away_team, line = NA_real_,
        american_odds = as.integer(row$odds_away), fetch_time = ft
      )
    }
  }

  if (length(rows) == 0) return(NULL)
  result <- bind_rows(rows)
  normalize_book_odds_frame(result)
}

#' Convert an Odds API long-format frame (one row per outcome) to canonical
#' shape. Expected columns: id, bookmaker_key, market_key, outcomes_name,
#' outcomes_price, outcomes_point. Optionally honors a fetch_time column;
#' falls back to Sys.time() when absent so staleness chips don't show stale.
#'
#' @param raw Odds API long tibble
#' @return normalized tibble (output of normalize_book_odds_frame), or NULL
odds_api_to_canonical <- function(raw) {
  if (is.null(raw) || nrow(raw) == 0) return(NULL)
  ft <- if ("fetch_time" %in% names(raw)) raw$fetch_time else rep(Sys.time(), nrow(raw))
  # Replace any NA fetch_times with now (NA breaks staleness chip rendering).
  ft <- as.POSIXct(ft, tz = "UTC")
  ft[is.na(ft)] <- Sys.time()

  result <- raw %>%
    transmute(
      game_id      = id,
      market_name  = market_key,
      bet_on       = outcomes_name,
      line         = as.numeric(outcomes_point),
      american_odds = as.integer(outcomes_price),
      fetch_time   = ft
    )
  normalize_book_odds_frame(result)
}
