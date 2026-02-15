# Acquire Historical Derivative Odds from The Odds API
# Purpose: Fetch Q1, Q2, Q3, Q4, H1 moneyline closing odds for NFL games
#
# API Documentation: https://the-odds-api.com/historical-odds-data/
#
# IMPORTANT: Historical derivative markets are only available from May 2023 onwards.
#
# Cost: 10 quota units per region per market per event
# With 3 regions (us, us2, eu) and 5 markets:
# Cost per game = 10 × 3 × 5 = 150 quota units

setwd("~/NFLWork/Answer Keys")
library(httr)
library(jsonlite)
library(dplyr)
library(purrr)
library(lubridate)
library(tidyr)
library(duckdb)
library(DBI)

# =============================================================================
# CONFIGURATION
# =============================================================================

api_key <- Sys.getenv("ODDS_API_KEY")
if (api_key == "") {
  stop("ODDS_API_KEY environment variable not set. Set it with Sys.setenv(ODDS_API_KEY = 'your_key')")
}

# All available derivative markets

# 2-way moneyline markets: home/away (ties excluded)
MARKETS_ML_2WAY <- c("h2h_q1", "h2h_q2", "h2h_q3", "h2h_q4", "h2h_h1", "h2h_h2")

# 3-way moneyline markets: home/away/tie (ties included as distinct outcome)
MARKETS_ML_3WAY <- c("h2h_3_way_q1", "h2h_3_way_q2", "h2h_3_way_q3", "h2h_3_way_q4", "h2h_3_way_h1", "h2h_3_way_h2")

# Spreads markets: home/away with point spread
MARKETS_SPREADS <- c("spreads_q1", "spreads_q2", "spreads_q3", "spreads_q4", "spreads_h1", "spreads_h2")
MARKETS_ALT_SPREADS <- c("alternate_spreads_q1", "alternate_spreads_q2", "alternate_spreads_q3", "alternate_spreads_q4", "alternate_spreads_h1", "alternate_spreads_h2")

# Totals markets: over/under with point total
MARKETS_TOTALS <- c("totals_q1", "totals_q2", "totals_q3", "totals_q4", "totals_h1", "totals_h2")
MARKETS_ALT_TOTALS <- c("alternate_totals_q1", "alternate_totals_q2", "alternate_totals_q3", "alternate_totals_q4", "alternate_totals_h1", "alternate_totals_h2")

# Market groups
ALL_ML_MARKETS <- c(MARKETS_ML_2WAY, MARKETS_ML_3WAY)
ALL_SPREADS_MARKETS <- c(MARKETS_SPREADS, MARKETS_ALT_SPREADS)
ALL_TOTALS_MARKETS <- c(MARKETS_TOTALS, MARKETS_ALT_TOTALS)

# Default: all markets
ALL_MARKETS <- c(ALL_ML_MARKETS, ALL_SPREADS_MARKETS, ALL_TOTALS_MARKETS)

REGIONS <- "us,us2,eu"
MIN_DATE <- as.Date("2023-05-01")
API_DELAY <- 0.5

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Check API quota
check_quota <- function() {
  res <- GET(
    url = "https://api.the-odds-api.com/v4/sports",
    query = list(apiKey = api_key)
  )

  remaining <- as.numeric(headers(res)$`x-requests-remaining`)
  used <- headers(res)$`x-requests-used`

  cat("API Quota Status:\n")
  cat("  Remaining:", remaining, "\n")
  cat("  Used:", used, "\n")

  remaining
}

#' Get historical events for a specific date
get_historical_events <- function(date, sport_key = "americanfootball_nfl") {
  snapshot <- paste0(format(date, "%Y-%m-%d"), "T14:00:00Z")

  res <- GET(
    url = paste0("https://api.the-odds-api.com/v4/historical/sports/", sport_key, "/events"),
    query = list(
      apiKey = api_key,
      date = snapshot,
      dateFormat = "iso"
    )
  )

  if (http_error(res)) {
    warning(paste("Failed to get events for date:", date))
    return(tibble())
  }

  parsed <- fromJSON(content(res, "text"), flatten = TRUE)

  if (is.null(parsed$data) || length(parsed$data) == 0) {
    return(tibble())
  }

  as_tibble(parsed$data) %>%
    mutate(query_date = date)
}

#' Get historical odds for a specific event
get_historical_event_odds <- function(
    event_id,
    commence_time,
    markets,
    regions = REGIONS,
    sport_key = "americanfootball_nfl"
) {
  if (is.character(commence_time)) {
    commence_time <- ymd_hms(commence_time, tz = "UTC")
  }
  snapshot <- format(commence_time - minutes(15), "%Y-%m-%dT%H:%M:%SZ")

  res <- GET(
    url = paste0("https://api.the-odds-api.com/v4/historical/sports/", sport_key,
                 "/events/", event_id, "/odds"),
    query = list(
      apiKey = api_key,
      date = snapshot,
      regions = regions,
      markets = paste(markets, collapse = ","),
      oddsFormat = "american",
      dateFormat = "iso"
    )
  )

  if (http_error(res)) {
    warning(paste("Failed to get odds for event:", event_id))
    return(tibble())
  }

  parsed <- fromJSON(content(res, "text"), flatten = TRUE)

  if (is.null(parsed$data) || length(parsed$data) == 0) {
    return(tibble())
  }

  odds_data <- parsed$data

  if (is.null(odds_data$bookmakers) || length(odds_data$bookmakers) == 0) {
    return(tibble())
  }

  tibble(
    event_id = event_id,
    snapshot_time = parsed$timestamp,
    data = list(odds_data)
  )
}

#' Flatten odds data from API response
#' Handles moneylines (2-way, 3-way), spreads, and totals
flatten_derivative_odds <- function(raw_odds, home_team, away_team) {
  if (nrow(raw_odds) == 0) return(tibble())

  tryCatch({
    odds_data <- raw_odds$data[[1]]

    if (is.null(odds_data$bookmakers) || nrow(odds_data$bookmakers) == 0) {
      return(tibble())
    }

    bookmakers <- odds_data$bookmakers

    expanded <- bookmakers %>%
      rename(bookmaker_key = key, bookmaker_title = title) %>%
      unnest(markets, names_sep = "_") %>%
      rename(market = markets_key) %>%
      unnest(markets_outcomes, names_sep = "_")

    # Determine market type based on market name
    all_results <- list()

    # Process moneyline markets (h2h_*, h2h_3_way_*)
    ml_markets <- expanded %>% filter(grepl("^h2h", market))
    if (nrow(ml_markets) > 0) {
      ml_result <- flatten_moneyline_odds(ml_markets, home_team, away_team, raw_odds, odds_data)
      if (nrow(ml_result) > 0) all_results[[length(all_results) + 1]] <- ml_result
    }

    # Process spreads markets (spreads_*, alternate_spreads_*)
    spreads_markets <- expanded %>% filter(grepl("spread", market))
    if (nrow(spreads_markets) > 0) {
      spreads_result <- flatten_spreads_odds(spreads_markets, home_team, away_team, raw_odds, odds_data)
      if (nrow(spreads_result) > 0) all_results[[length(all_results) + 1]] <- spreads_result
    }

    # Process totals markets (totals_*, alternate_totals_*)
    totals_markets <- expanded %>% filter(grepl("total", market))
    if (nrow(totals_markets) > 0) {
      totals_result <- flatten_totals_odds(totals_markets, home_team, away_team, raw_odds, odds_data)
      if (nrow(totals_result) > 0) all_results[[length(all_results) + 1]] <- totals_result
    }

    if (length(all_results) == 0) return(tibble())

    bind_rows(all_results)
  }, error = function(e) {
    warning(paste("Error flattening odds:", e$message))
    tibble()
  })
}

#' Flatten moneyline odds (2-way and 3-way)
flatten_moneyline_odds <- function(expanded, home_team, away_team, raw_odds, odds_data) {
  is_3way <- any(grepl("3_way", expanded$market))

  final <- expanded %>%
    mutate(
      side = case_when(
        markets_outcomes_name == home_team ~ "home",
        markets_outcomes_name == away_team ~ "away",
        markets_outcomes_name == "Draw" | markets_outcomes_name == "Tie" ~ "tie",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(side)) %>%
    select(bookmaker_key, market, side, price = markets_outcomes_price)

  if (is_3way) {
    final <- final %>%
      pivot_wider(names_from = side, values_from = price, names_prefix = "odds_",
                  values_fn = \(x) x[1]) %>%
      mutate(
        event_id = raw_odds$event_id[1],
        snapshot_time = raw_odds$snapshot_time[1],
        home_team = home_team,
        away_team = away_team,
        commence_time = odds_data$commence_time,
        market_type = "moneyline_3way",
        line = NA_real_
      )
  } else {
    final <- final %>%
      filter(side != "tie") %>%
      pivot_wider(names_from = side, values_from = price, names_prefix = "odds_",
                  values_fn = \(x) x[1]) %>%
      mutate(
        event_id = raw_odds$event_id[1],
        snapshot_time = raw_odds$snapshot_time[1],
        home_team = home_team,
        away_team = away_team,
        commence_time = odds_data$commence_time,
        market_type = "moneyline_2way",
        line = NA_real_
      )
  }

  final
}

#' Flatten spreads odds
#' Structure: outcomes have name (team), price, and point (spread value)
#'
#' IMPORTANT: The API returns spread points from each team's perspective:
#'   - Home team outcome has point = home_spread (e.g., -1.5 if home favored)
#'   - Away team outcome has point = away_spread (e.g., +1.5 if away getting points)
#'
#' We normalize to HOME perspective: home_line = home_point = -away_point
flatten_spreads_odds <- function(expanded, home_team, away_team, raw_odds, odds_data) {
  # Check if point column exists
  if (!"markets_outcomes_point" %in% names(expanded)) {
    return(tibble())
  }

  final <- expanded %>%
    mutate(
      side = case_when(
        markets_outcomes_name == home_team ~ "home",
        markets_outcomes_name == away_team ~ "away",
        TRUE ~ NA_character_
      ),
      raw_point = as.numeric(markets_outcomes_point),
      # Normalize point to HOME perspective
      # Home team point is already from home perspective
      # Away team point needs to be negated to get home perspective
      # e.g., away +1.5 → home -1.5
      home_line = case_when(
        side == "home" ~ raw_point,
        side == "away" ~ -raw_point,
        TRUE ~ NA_real_
      )
    ) %>%
    filter(!is.na(side)) %>%
    select(bookmaker_key, market, side, price = markets_outcomes_price, home_line)

  # For spreads, each row is a unique combination of book + market + home_line
  # We pivot so each line gets odds_home and odds_away
  # Use values_fn to handle duplicates (take first value)
  final <- final %>%
    pivot_wider(
      id_cols = c(bookmaker_key, market, home_line),
      names_from = side,
      values_from = price,
      names_prefix = "odds_",
      values_fn = \(x) x[1]
    ) %>%
    rename(line = home_line) %>%
    mutate(
      event_id = raw_odds$event_id[1],
      snapshot_time = raw_odds$snapshot_time[1],
      home_team = home_team,
      away_team = away_team,
      commence_time = odds_data$commence_time,
      market_type = ifelse(grepl("alternate", market), "spreads_alt", "spreads")
    )

  final
}

#' Flatten totals odds
#' Structure: outcomes have name (Over/Under), price, and point (total line)
flatten_totals_odds <- function(expanded, home_team, away_team, raw_odds, odds_data) {
  # Check if point column exists
  if (!"markets_outcomes_point" %in% names(expanded)) {
    return(tibble())
  }

  final <- expanded %>%
    mutate(
      side = case_when(
        markets_outcomes_name == "Over" ~ "over",
        markets_outcomes_name == "Under" ~ "under",
        TRUE ~ NA_character_
      ),
      point = as.numeric(markets_outcomes_point)
    ) %>%
    filter(!is.na(side)) %>%
    select(bookmaker_key, market, side, price = markets_outcomes_price, point)

  # Pivot so each line gets odds_over and odds_under
  # Use values_fn to handle duplicates (take first value)
  final <- final %>%
    pivot_wider(
      id_cols = c(bookmaker_key, market, point),
      names_from = side,
      values_from = price,
      names_prefix = "odds_",
      values_fn = \(x) x[1]
    ) %>%
    rename(line = point) %>%
    mutate(
      event_id = raw_odds$event_id[1],
      snapshot_time = raw_odds$snapshot_time[1],
      home_team = home_team,
      away_team = away_team,
      commence_time = odds_data$commence_time,
      market_type = ifelse(grepl("alternate", market), "totals_alt", "totals")
    )

  final
}

# =============================================================================
# MAIN ACQUISITION FUNCTION
# =============================================================================

#' Acquire derivative odds for NFL games
#'
#' @param markets Vector of markets to fetch (default: all 5 markets)
#' @param start_date Start date (default: 2023-05-01)
#' @param end_date End date (default: today)
#' @param max_games Maximum games to fetch (NULL for all)
#' @param append If TRUE, append to existing table; if FALSE, overwrite
#' @return Combined tibble of all derivative odds
acquire_derivative_odds <- function(
    markets = ALL_MARKETS,
    start_date = MIN_DATE,
    end_date = Sys.Date(),
    max_games = NULL,
    append = FALSE
) {
  con <- dbConnect(duckdb(), dbdir = "pbp.duckdb")
  on.exit(dbDisconnect(con))

  # Get games from database
  games_needed <- dbGetQuery(con, sprintf("
    SELECT DISTINCT
      game_id,
      home_team,
      away_team,
      game_date,
      game_start_time
    FROM nfl_betting_pbp
    WHERE game_date >= '%s'
      AND game_date <= '%s'
    ORDER BY game_date
  ", start_date, end_date))

  cat("Found", nrow(games_needed), "games in date range\n")

  if (!is.null(max_games) && max_games < nrow(games_needed)) {
    games_needed <- games_needed[1:max_games, ]
    cat("Limited to first", max_games, "games\n")
  }

  # Get event IDs
  unique_dates <- unique(as.Date(games_needed$game_date))
  cat("Fetching event IDs for", length(unique_dates), "dates...\n")

  all_events <- map_dfr(unique_dates, function(d) {
    events <- get_historical_events(d)
    Sys.sleep(API_DELAY)
    events
  }, .progress = TRUE)

  cat("Found", nrow(all_events), "events from API\n")

  # Match events to games
  games_with_events <- games_needed %>%
    mutate(game_date = as.Date(game_date)) %>%
    inner_join(
      all_events %>%
        mutate(
          game_date = as.Date(query_date),
          commence_time = ymd_hms(commence_time, tz = "UTC")
        ) %>%
        select(event_id = id, home_team, away_team, commence_time, game_date),
      by = c("home_team", "away_team", "game_date")
    )

  cat("Matched", nrow(games_with_events), "games to event IDs\n")

  if (nrow(games_with_events) == 0) {
    warning("No games matched to events.")
    return(tibble())
  }

  # Fetch odds
  cat("\nFetching derivative odds for", nrow(games_with_events), "games...\n")
  cat("Markets:", paste(markets, collapse = ", "), "\n")
  cat("Estimated quota:", nrow(games_with_events) * length(markets) * 3 * 10, "units\n\n")

  all_odds <- list()

  for (i in 1:nrow(games_with_events)) {
    game <- games_with_events[i, ]

    if (i %% 25 == 0) {
      cat("Progress:", i, "/", nrow(games_with_events),
          "(", game$home_team, "vs", game$away_team, ")\n")
    }

    odds <- get_historical_event_odds(
      event_id = game$event_id,
      commence_time = game$commence_time,
      markets = markets
    )

    if (nrow(odds) > 0) {
      flat_odds <- flatten_derivative_odds(odds, game$home_team, game$away_team)
      if (nrow(flat_odds) > 0) {
        flat_odds$game_id <- game$game_id
        all_odds[[length(all_odds) + 1]] <- flat_odds
      }
    }

    Sys.sleep(API_DELAY)

    # Checkpoint every 50 games
    if (i %% 50 == 0 && length(all_odds) > 0) {
      combined <- bind_rows(all_odds)
      dbWriteTable(con, "nfl_derivative_odds_temp", combined, overwrite = TRUE)
      cat("Checkpoint saved:", nrow(combined), "rows\n")
    }
  }

  if (length(all_odds) == 0) {
    warning("No derivative odds fetched")
    return(tibble())
  }

  new_odds <- bind_rows(all_odds)

  # Save to database
  if (append && dbExistsTable(con, "nfl_derivative_closing_odds")) {
    existing_odds <- dbGetQuery(con, "SELECT * FROM nfl_derivative_closing_odds")
    final_odds <- bind_rows(existing_odds, new_odds) %>%
      distinct()  # Remove duplicates
  } else {
    final_odds <- new_odds
  }

  dbWriteTable(con, "nfl_derivative_closing_odds", final_odds, overwrite = TRUE)

  cat("\n=== ACQUISITION COMPLETE ===\n")
  cat("New rows:", nrow(new_odds), "\n")
  cat("Total rows:", nrow(final_odds), "\n")
  cat("Unique games:", n_distinct(final_odds$game_id), "\n")
  cat("\nRows by market:\n")
  print(final_odds %>% count(market))

  final_odds
}

# =============================================================================
# USAGE
# =============================================================================

cat("=== NFL DERIVATIVE ODDS ACQUISITION ===\n\n")

remaining <- check_quota()

cat("\nAvailable market groups:\n")
cat("  ALL_ML_MARKETS      - Moneylines (2-way and 3-way)\n")
cat("  ALL_SPREADS_MARKETS - Spreads (main and alternate)\n")
cat("  ALL_TOTALS_MARKETS  - Totals (main and alternate)\n")
cat("  ALL_MARKETS         - All markets combined\n\n")

cat("Usage examples:\n")
cat("  # Fetch all markets for all games:\n")
cat("  all_odds <- acquire_derivative_odds()\n\n")
cat("  # Fetch only spreads markets:\n")
cat("  odds <- acquire_derivative_odds(markets = ALL_SPREADS_MARKETS)\n\n")
cat("  # Fetch only totals markets:\n")
cat("  odds <- acquire_derivative_odds(markets = ALL_TOTALS_MARKETS)\n\n")
cat("  # Fetch spreads and totals together:\n")
cat("  odds <- acquire_derivative_odds(markets = c(ALL_SPREADS_MARKETS, ALL_TOTALS_MARKETS))\n\n")
cat("  # Test with 5 games:\n")
cat("  test <- acquire_derivative_odds(max_games = 5)\n\n")
cat("  # Append new markets to existing data:\n")
cat("  odds <- acquire_derivative_odds(markets = ALL_SPREADS_MARKETS, append = TRUE)\n")
