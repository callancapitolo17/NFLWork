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
    # H2 = MLB 2nd halves (innings 6-9) — emitted by the BKM scraper from
    # league 503. Not currently modeled, but recognized here so H2 rows
    # are tagged explicitly instead of silently defaulting to FG (which
    # would mis-classify them as full-game data).
    str_detect(market_name, "_h2$")            ~ "H2",
    TRUE                                       ~ "FG"
  )
}

.derive_market_type <- function(market_name) {
  market_name %>%
    str_replace("_1st_[357]_innings$", "") %>%
    str_replace("_(fg|f3|f5|f7|h2)$", "")
}

#' Canonical, book-agnostic bet identity hash shared by the model path
#' (MLB.R) and the market-edge path (market_edge.R). Hashes
#' (id, base_market_type, period, normalized_line, bet_on) so that the same
#' wager always maps to the same id regardless of the verbose market string
#' or the alternate_/main labeling. `market_type` may be the verbose model
#' string OR a canonical type; callers should pass the canonical type, but the
#' alternate_ prefix is stripped here defensively.
compute_bet_row_id <- function(id, market_type, period, line, bet_on) {
  base_mt   <- gsub("^alternate_", "", market_type)
  line_norm <- ifelse(is.na(line), "", as.character(round(as.numeric(line), 1)))
  key <- paste(id, base_mt, period, line_norm, bet_on, sep = "|")
  vapply(key, function(s) digest::digest(s, algo = "md5"), character(1),
         USE.NAMES = FALSE)
}

# When matching a bet to per-book candidates, treat alt and main as the
# same bucket. Different scrapers use different conventions:
#   - WZ/BKM/Bet105 label alt rows "alternate_spreads"/"alternate_totals"
#   - DK/FD (via get_dk_odds/get_fd_odds in Tools.R) collapse alt rows into
#     "spreads"/"totals" with the alt line carried in the `line` column
# The closest-line picker decides which row wins, so we just union the
# possible labels here. Moneyline (h2h) has no alt variant.
.related_market_types <- function(mt) {
  if (mt == "spreads"           || mt == "alternate_spreads")
    return(c("spreads", "alternate_spreads"))
  if (mt == "totals"            || mt == "alternate_totals")
    return(c("totals", "alternate_totals"))
  if (mt == "h2h" || mt == "h2h_3way")
    return(c("h2h", "h2h_3way"))
  mt
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
    # Over total: pick_high = TRUE -> slice_max -> higher total (harder to go over).
    # Under total: pick_high = FALSE -> slice_min -> lower total (harder to go under).
    # Spread dog (model_line > 0): pick_high = TRUE -> slice_max -> larger +N
    #   (less cushion = harder to cover).
    # Spread favorite (model_line < 0): pick_high = FALSE -> slice_min -> more
    #   negative line (must cover by more = harder).
    if (grepl("^Over", bet_on, ignore.case = TRUE)) {
      pick_high <- TRUE
    } else if (grepl("^Under", bet_on, ignore.case = TRUE)) {
      pick_high <- FALSE
    } else {
      # It's a spread; use line sign to determine dog vs favorite
      pick_high <- model_line > 0
    }
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
                  derived_fair_odds = numeric(),
                  fetch_time = as.POSIXct(character()),
                  game_start_time = as.POSIXct(character(), tz = "UTC")))
  }

  # --- Defensive column renames for bets frame ---
  # all_bets_combined uses `id` (Odds API game id), not `game_id`.
  if (!"game_id" %in% names(bets) && "id" %in% names(bets)) {
    bets <- bets %>% rename(game_id = id)
  }

  # Always re-derive market_type from the canonical `market` column.
  # Production callers (MLB.R::wz_alt_bets etc.) pre-set market_type to a
  # stale bare value like "totals" for alt bets that have
  # market="alternate_totals_fg" — the join then fails to match book frames
  # where market="alternate_totals". Re-deriving unconditionally keeps the
  # join honest: .derive_market_type("alternate_totals_fg") == "alternate_totals"
  # matches the book frame, and for non-alt markets the function is idempotent
  # (e.g. .derive_market_type("spreads") == "spreads") so nothing changes.
  bets <- bets %>%
    mutate(market_type = .derive_market_type(market))

  # Period — derive when the column is absent OR when individual rows are
  # NA. Both happen in production: compare_alts_to_samples emits frames
  # without a period column, then bind_rows in MLB.R merges them with
  # non-alt frames that do carry period, filling the alt rows with NA.
  # The previous "absent column" guard never fired because the column
  # technically existed. Row-level if_else preserves explicitly-set
  # period values (e.g. "F5" set by callers / format_bets_table) and
  # fills NAs from the canonical `market` column.
  if (!"period" %in% names(bets)) {
    bets <- bets %>% mutate(period = .derive_period(market))
  } else {
    bets <- bets %>% mutate(period = if_else(is.na(period),
                                             .derive_period(market),
                                             period))
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
    bet_game_start <- if ("pt_start_time" %in% names(bet)) bet$pt_start_time else as.POSIXct(NA, tz = "UTC")
    opposite_col_value <- if ("opposite_side" %in% names(bet)) bet$opposite_side else NA_character_
    home_team_value    <- if ("home_team"     %in% names(bet)) bet$home_team     else NA_character_
    away_team_value    <- if ("away_team"     %in% names(bet)) bet$away_team     else NA_character_
    sides <- .sides_for_bet(bet$bet_on, bet$market_type, opposite_col_value,
                            home_team = home_team_value, away_team = away_team_value)
    side_labels <- c(pick = sides$pick, opposite = sides$opposite)

    for (book_name in names(book_odds_by_book)) {
      book_frame <- book_odds_by_book[[book_name]]
      if (is.null(book_frame) || nrow(book_frame) == 0) next

      # Pick'em (line-0 spread): derive draw-no-bet from the book's period
      # winner instead of falling back to its nearest run line. Priority:
      # 2-way h2h > 3-way h2h_3way > true 0-handicap spread.
      is_pickem <- bet$market_type %in% c("spreads", "alternate_spreads") &&
                   !is.na(bet$line) && bet$line == 0
      if (is_pickem) {
        winners <- book_frame %>%
          filter(game_id == bet$game_id, period == bet$period,
                 market %in% c("h2h", "h2h_3way", "spreads", "alternate_spreads"))
        winners <- winners %>% filter(
          market %in% c("h2h", "h2h_3way") |
          (market %in% c("spreads", "alternate_spreads") & !is.na(line) & abs(line) < 1e-9)
        )
        if (nrow(winners) == 0) next   # cell -> "—"

        h2h_rows <- winners %>% filter(market == "h2h")
        h3w_rows <- winners %>% filter(market == "h2h_3way")
        sp0_rows <- winners %>% filter(market %in% c("spreads", "alternate_spreads"))
        use_3way <- nrow(h2h_rows) < 2 && nrow(h3w_rows) >= 3
        use_sp0  <- nrow(h2h_rows) < 2 && nrow(h3w_rows) < 3 && nrow(sp0_rows) >= 2

        if (use_3way) {
          home_row <- h3w_rows %>% filter(side == home_team_value) %>% slice_head(n = 1)
          away_row <- h3w_rows %>% filter(side == away_team_value) %>% slice_head(n = 1)
          tie_row  <- h3w_rows %>% filter(side == "Tie")           %>% slice_head(n = 1)
          if (nrow(home_row) == 0 || nrow(away_row) == 0 || nrow(tie_row) == 0) next
          dnb <- derive_pickem_american(home_row$american_odds, away_row$american_odds, tie_row$american_odds)
          ft  <- home_row$fetch_time
        } else if (use_sp0) {
          home_row <- sp0_rows %>% filter(side == home_team_value) %>% slice_head(n = 1)
          away_row <- sp0_rows %>% filter(side == away_team_value) %>% slice_head(n = 1)
          if (nrow(home_row) == 0 || nrow(away_row) == 0) next
          dnb <- derive_pickem_american(home_row$american_odds, away_row$american_odds, NA)
          ft  <- home_row$fetch_time
        } else {
          home_row <- h2h_rows %>% filter(side == home_team_value) %>% slice_head(n = 1)
          away_row <- h2h_rows %>% filter(side == away_team_value) %>% slice_head(n = 1)
          if (nrow(home_row) == 0 || nrow(away_row) == 0) next
          dnb <- derive_pickem_american(home_row$american_odds, away_row$american_odds, NA)
          ft  <- home_row$fetch_time
        }
        if (is.na(dnb$home_raw_dnb) || is.na(dnb$away_raw_dnb)) next

        bet_is_home <- identical(bet$bet_on, home_team_value)
        pick_raw  <- if (bet_is_home) dnb$home_raw_dnb  else dnb$away_raw_dnb
        opp_raw   <- if (bet_is_home) dnb$away_raw_dnb  else dnb$home_raw_dnb
        pick_fair <- if (bet_is_home) dnb$home_fair_dnb else dnb$away_fair_dnb
        opp_fair  <- if (bet_is_home) dnb$away_fair_dnb else dnb$home_fair_dnb

        for (slot in c("pick", "opposite")) {
          american <- if (slot == "pick") pick_raw  else opp_raw
          fair_odd <- if (slot == "pick") pick_fair else opp_fair
          if (is.na(american)) next
          k <- k + 1L
          out_rows[[k]] <- tibble(
            bet_row_id      = bet$bet_row_id,
            game_id         = bet$game_id,
            market          = bet$market_type,
            period          = bet$period,
            side            = slot,
            bookmaker       = book_name,
            line            = bet$line,
            line_quoted     = 0,
            is_exact_line   = TRUE,
            american_odds   = as.integer(round(american)),
            derived_fair_odds = as.numeric(fair_odd),
            fetch_time      = ft,
            game_start_time = bet_game_start
          )
        }
        next   # skip the run-line slot loop for pick'em bets
      }

      for (slot in names(side_labels)) {
        side_value <- side_labels[[slot]]
        if (is.na(side_value)) next

        candidates <- book_frame %>%
          filter(game_id == bet$game_id,
                 market  %in% .related_market_types(bet$market_type),
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
          bet_row_id        = bet$bet_row_id,
          game_id           = bet$game_id,
          market            = bet$market_type,
          period            = bet$period,
          side              = slot,
          bookmaker         = book_name,
          line              = bet$line,
          line_quoted       = chosen$line_quoted,
          is_exact_line     = chosen$is_exact_line,
          american_odds     = as.integer(chosen$american_odds),
          derived_fair_odds = NA_real_,
          fetch_time        = ft,
          game_start_time   = bet_game_start
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
                  derived_fair_odds = numeric(),
                  fetch_time = as.POSIXct(character()),
                  game_start_time = as.POSIXct(character(), tz = "UTC")))
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
#' @param book_name Optional book label used in coverage-drop log messages
#'               (e.g. "dk", "fd", "wagerzon"). Defaults to "<book>".
#' @return normalized tibble (output of normalize_book_odds_frame), or NULL
scraper_to_canonical <- function(raw, lookup, book_name = NULL) {
  if (is.null(raw) || nrow(raw) == 0) return(NULL)
  if (!"fetch_time" %in% names(raw)) raw$fetch_time <- as.POSIXct(NA, tz = "UTC")

  # Replace NA fetch_times with now for consistency with odds_api_to_canonical.
  if (!is.null(raw$fetch_time)) {
    raw$fetch_time[is.na(raw$fetch_time)] <- Sys.time()
  }

  # Coverage diagnostic — surface rows that would be silently dropped by inner_join.
  # The #1 invisible failure mode of the bets tab is unmapped team-name pairs.
  dropped <- anti_join(raw, lookup, by = c("home_team", "away_team"))
  if (nrow(dropped) > 0) {
    dropped_pairs <- dropped %>%
      distinct(home_team, away_team) %>%
      mutate(pair = paste(home_team, "vs", away_team)) %>%
      pull(pair)
    message(sprintf(
      "[scraper_to_canonical] %s: %d row(s) across %d game(s) dropped (no game_id match): %s",
      book_name %||% "<book>",
      nrow(dropped),
      length(dropped_pairs),
      paste(dropped_pairs, collapse = "; ")
    ))
  }

  # Join to Odds API game IDs. Scraper frames have home_team + away_team
  # already resolved to canonical names by resolve_offshore_teams().
  joined <- raw %>%
    inner_join(lookup, by = c("home_team", "away_team")) %>%
    rename(game_id = id)

  # Surface silent drops: if any raw rows did not match the lookup, log the
  # unmatched team pairs so a roster rename or scraper drift surfaces in the
  # run output instead of silently producing fewer pills on the dashboard.
  n_dropped <- nrow(raw) - nrow(joined)
  if (n_dropped > 0) {
    unmatched <- raw %>%
      anti_join(lookup, by = c("home_team", "away_team")) %>%
      distinct(home_team, away_team)
    warning(sprintf(
      "[scraper_to_canonical] dropped %d row(s) for %d unmatched team pair(s): %s",
      n_dropped, nrow(unmatched),
      paste(paste(unmatched$away_team, "@", unmatched$home_team), collapse = "; ")
    ))
  }

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
    # 3-way moneyline (h2h_3way_*): emit a "Tie" row alongside home/away.
    if (!is.na(row$odds_home) && is.na(row$line) &&
        "odds_tie" %in% names(row) && !is.na(row$odds_tie)) {
      rows[[length(rows) + 1]] <- tibble(
        game_id = gid, market_name = mkt,
        bet_on = "Tie", line = NA_real_,
        american_odds = as.integer(row$odds_tie), fetch_time = ft
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
  ft <- as.POSIXct(ft, tz = "UTC")

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

#' Derive the draw-no-bet (pick'em) American odds for a book's period winner.
#'
#' Two source shapes:
#'  - 2-way winner (no tie outcome): the market already excludes ties, so the
#'    devigged probabilities ARE the DNB probabilities. Raw = input unchanged.
#'  - 3-way winner: probit-devig home/away/tie, drop the tie, renormalize.
#'
#' Depends on Tools.R (devig_american, devig_american_3way, prob_to_american).
#'
#' @param home_raw American odds, home side (numeric, e.g. -180)
#' @param away_raw American odds, away side
#' @param tie_raw  American odds for the tie outcome, or NA for 2-way
#' @return list with home_raw_dnb, away_raw_dnb, home_fair_dnb, away_fair_dnb
derive_pickem_american <- function(home_raw, away_raw, tie_raw = NA) {
  na_result <- list(home_raw_dnb = NA_real_, away_raw_dnb = NA_real_,
                    home_fair_dnb = NA_real_, away_fair_dnb = NA_real_)
  if (is.na(home_raw) || is.na(away_raw)) return(na_result)

  # Guarded prob -> American: reuse Tools.R::prob_to_american for valid probs.
  to_amer <- function(p) if (is.na(p) || p <= 0 || p >= 1) NA_real_ else prob_to_american(p)

  if (is.na(tie_raw)) {
    # 2-way: raw_dnb = input; fair_dnb = probit 2-way devig.
    # devig_american(odd1, odd2) -> df with $p1 (corresponds to odd1), $p2 (odd2).
    devig <- devig_american(away_raw, home_raw)   # p1 = away, p2 = home
    if (is.null(devig) || any(is.na(c(devig$p1, devig$p2)))) {
      return(list(home_raw_dnb = home_raw, away_raw_dnb = away_raw,
                  home_fair_dnb = NA_real_, away_fair_dnb = NA_real_))
    }
    return(list(
      home_raw_dnb  = home_raw,
      away_raw_dnb  = away_raw,
      home_fair_dnb = to_amer(devig$p2),
      away_fair_dnb = to_amer(devig$p1)
    ))
  }
  # 3-way: devig home/away/tie -> drop the tie -> renormalize home/away.
  devig3 <- devig_american_3way(home_raw, away_raw, tie_raw)
  if (is.null(devig3) || any(is.na(c(devig3$p_home, devig3$p_away)))) return(na_result)
  denom_f <- devig3$p_home + devig3$p_away
  if (!is.finite(denom_f) || denom_f <= 0) return(na_result)
  # Raw DNB: same drop-tie-renormalize on RAW implied probs (no devig).
  implied <- function(american) {
    if (is.na(american)) return(NA_real_)
    if (american < 0) -american / (-american + 100) else 100 / (american + 100)
  }
  q_h <- implied(home_raw); q_a <- implied(away_raw); denom_r <- q_h + q_a
  list(
    home_raw_dnb  = to_amer(q_h / denom_r),
    away_raw_dnb  = to_amer(q_a / denom_r),
    home_fair_dnb = to_amer(devig3$p_home / denom_f),
    away_fair_dnb = to_amer(devig3$p_away / denom_f)
  )
}
