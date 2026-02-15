# Acquire Historical Derivative Odds from The Odds API for CBB
# Purpose: Fetch H1, H2 spreads, totals, moneylines, and team totals for CBB games
#
# API Documentation: https://the-odds-api.com/historical-odds-data/
#
# IMPORTANT: Historical derivative markets are only available from May 2023 onwards.
#
# Cost: 10 quota units per region per market per event
# With 3 regions (us, us2, eu) and multiple markets:
# Cost per game depends on number of markets requested

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

SPORT_KEY <- "basketball_ncaab"
DB_PATH <- "cbb.duckdb"

# CBB only has halves, no quarters
# 2-way moneyline markets: home/away (ties excluded - rare in basketball)
MARKETS_ML_2WAY <- c("h2h_h1", "h2h_h2")

# 3-way moneyline markets: home/away/tie (ties included as distinct outcome)
MARKETS_ML_3WAY <- c("h2h_3_way_h1", "h2h_3_way_h2")

# Spreads markets: home/away with point spread
MARKETS_SPREADS <- c("spreads_h1", "spreads_h2")
MARKETS_ALT_SPREADS <- c("alternate_spreads_h1", "alternate_spreads_h2")

# Totals markets: over/under with point total
MARKETS_TOTALS <- c("totals_h1", "totals_h2")
MARKETS_ALT_TOTALS <- c("alternate_totals_h1", "alternate_totals_h2")

# Team totals
MARKETS_TEAM_TOTALS <- c("team_totals_h1", "team_totals_h2")
MARKETS_ALT_TEAM_TOTALS <- c("alternate_team_totals_h1", "alternate_team_totals_h2")

# Market groups
ALL_ML_MARKETS <- c(MARKETS_ML_2WAY, MARKETS_ML_3WAY)
ALL_SPREADS_MARKETS <- c(MARKETS_SPREADS, MARKETS_ALT_SPREADS)
ALL_TOTALS_MARKETS <- c(MARKETS_TOTALS, MARKETS_ALT_TOTALS)
ALL_TEAM_TOTALS_MARKETS <- c(MARKETS_TEAM_TOTALS, MARKETS_ALT_TEAM_TOTALS)

# Recommended: Start with main markets (not alternates) to save quota
MAIN_MARKETS <- c(MARKETS_ML_2WAY, MARKETS_SPREADS, MARKETS_TOTALS, MARKETS_TEAM_TOTALS)

# All markets (expensive)
ALL_MARKETS <- c(ALL_ML_MARKETS, ALL_SPREADS_MARKETS, ALL_TOTALS_MARKETS, ALL_TEAM_TOTALS_MARKETS)

REGIONS <- "us,us2,eu"
MIN_DATE <- as.Date("2023-05-01")  # Historical derivatives available from May 2023
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
get_historical_events <- function(date, sport_key = SPORT_KEY) {
  # Use noon Eastern for snapshot to catch games that day
  snapshot <- paste0(format(date, "%Y-%m-%d"), "T17:00:00Z")

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
    sport_key = SPORT_KEY
) {
  if (is.character(commence_time)) {
    commence_time <- ymd_hms(commence_time, tz = "UTC")
  }
  # Snapshot 15 min before tip-off for closing odds
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
#' Handles moneylines (2-way, 3-way), spreads, totals, and team totals
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
    totals_markets <- expanded %>% filter(grepl("^totals|^alternate_totals", market))
    if (nrow(totals_markets) > 0) {
      totals_result <- flatten_totals_odds(totals_markets, home_team, away_team, raw_odds, odds_data)
      if (nrow(totals_result) > 0) all_results[[length(all_results) + 1]] <- totals_result
    }

    # Process team totals markets
    team_totals_markets <- expanded %>% filter(grepl("team_total", market))
    if (nrow(team_totals_markets) > 0) {
      tt_result <- flatten_team_totals_odds(team_totals_markets, home_team, away_team, raw_odds, odds_data)
      if (nrow(tt_result) > 0) all_results[[length(all_results) + 1]] <- tt_result
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
flatten_spreads_odds <- function(expanded, home_team, away_team, raw_odds, odds_data) {
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
      home_line = case_when(
        side == "home" ~ raw_point,
        side == "away" ~ -raw_point,
        TRUE ~ NA_real_
      )
    ) %>%
    filter(!is.na(side)) %>%
    select(bookmaker_key, market, side, price = markets_outcomes_price, home_line)

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
flatten_totals_odds <- function(expanded, home_team, away_team, raw_odds, odds_data) {
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

#' Flatten team totals odds
flatten_team_totals_odds <- function(expanded, home_team, away_team, raw_odds, odds_data) {
  if (!"markets_outcomes_point" %in% names(expanded)) {
    return(tibble())
  }

  # Team totals have team name embedded in description
  final <- expanded %>%
    mutate(
      is_home_team = grepl(home_team, markets_outcomes_description, fixed = TRUE),
      is_away_team = grepl(away_team, markets_outcomes_description, fixed = TRUE),
      team_side = case_when(
        is_home_team ~ "home",
        is_away_team ~ "away",
        TRUE ~ NA_character_
      ),
      ou_side = case_when(
        markets_outcomes_name == "Over" ~ "over",
        markets_outcomes_name == "Under" ~ "under",
        TRUE ~ NA_character_
      ),
      point = as.numeric(markets_outcomes_point)
    ) %>%
    filter(!is.na(team_side), !is.na(ou_side)) %>%
    select(bookmaker_key, market, team_side, ou_side, price = markets_outcomes_price, point)

  if (nrow(final) == 0) return(tibble())

  # Pivot to get home_over, home_under, away_over, away_under columns
  final <- final %>%
    pivot_wider(
      id_cols = c(bookmaker_key, market, point),
      names_from = c(team_side, ou_side),
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
      market_type = ifelse(grepl("alternate", market), "team_totals_alt", "team_totals")
    )

  final
}

# =============================================================================
# MAIN ACQUISITION FUNCTION
# =============================================================================

#' Acquire derivative odds for CBB games
#'
#' @param markets Vector of markets to fetch (default: main markets only)
#' @param start_date Start date (default: 2023-05-01)
#' @param end_date End date (default: today)
#' @param max_games Maximum games to fetch (NULL for all)
#' @param append If TRUE, append to existing table; if FALSE, overwrite
#' @return Combined tibble of all derivative odds
acquire_derivative_odds <- function(
    markets = MAIN_MARKETS,
    start_date = MIN_DATE,
    end_date = Sys.Date(),
    max_games = NULL,
    append = FALSE
) {
  con <- dbConnect(duckdb(), dbdir = DB_PATH)
  on.exit(dbDisconnect(con))

  # Get games from database that have odds
  games_needed <- dbGetQuery(con, sprintf("
    SELECT DISTINCT
      o.id as odds_id,
      o.home_team,
      o.away_team,
      o.commence_time,
      DATE(o.commence_time) as game_date
    FROM cbb_closing_odds o
    WHERE DATE(o.commence_time) >= '%s'
      AND DATE(o.commence_time) <= '%s'
    ORDER BY o.commence_time
  ", start_date, end_date))

  cat("Found", nrow(games_needed), "games in date range\n")

  if (!is.null(max_games) && max_games < nrow(games_needed)) {
    games_needed <- games_needed[1:max_games, ]
    cat("Limited to first", max_games, "games\n")
  }

  # Get unique dates
  unique_dates <- unique(as.Date(games_needed$game_date))
  cat("Fetching event IDs for", length(unique_dates), "dates...\n")

  # Fetch events for each date
  all_events <- map_dfr(unique_dates, function(d) {
    events <- get_historical_events(d)
    Sys.sleep(API_DELAY)
    events
  }, .progress = TRUE)

  cat("Found", nrow(all_events), "events from API\n")

  # Match events to games by team names and date
  games_with_events <- games_needed %>%
    mutate(game_date = as.Date(game_date)) %>%
    inner_join(
      all_events %>%
        mutate(
          game_date = as.Date(query_date),
          commence_time_api = ymd_hms(commence_time, tz = "UTC")
        ) %>%
        select(event_id = id, home_team, away_team, commence_time_api, game_date),
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
      commence_time = game$commence_time_api,
      markets = markets
    )

    if (nrow(odds) > 0) {
      flat_odds <- flatten_derivative_odds(odds, game$home_team, game$away_team)
      if (nrow(flat_odds) > 0) {
        flat_odds$odds_id <- game$odds_id
        all_odds[[length(all_odds) + 1]] <- flat_odds
      }
    }

    Sys.sleep(API_DELAY)

    # Checkpoint every 50 games
    if (i %% 50 == 0 && length(all_odds) > 0) {
      combined <- bind_rows(all_odds)
      dbWriteTable(con, "cbb_derivative_odds_temp", combined, overwrite = TRUE)
      cat("Checkpoint saved:", nrow(combined), "rows\n")
    }
  }

  if (length(all_odds) == 0) {
    warning("No derivative odds fetched")
    return(tibble())
  }

  new_odds <- bind_rows(all_odds)

  # Save to database
  if (append && dbExistsTable(con, "cbb_derivative_closing_odds")) {
    existing_odds <- dbGetQuery(con, "SELECT * FROM cbb_derivative_closing_odds")
    final_odds <- bind_rows(existing_odds, new_odds) %>%
      distinct()
  } else {
    final_odds <- new_odds
  }

  dbWriteTable(con, "cbb_derivative_closing_odds", final_odds, overwrite = TRUE)

  cat("\n=== ACQUISITION COMPLETE ===\n")
  cat("New rows:", nrow(new_odds), "\n")
  cat("Total rows:", nrow(final_odds), "\n")
  cat("Unique games (odds_id):", n_distinct(final_odds$odds_id), "\n")
  cat("\nRows by market:\n")
  print(final_odds %>% count(market))

  final_odds
}

# =============================================================================
# USAGE
# =============================================================================

cat("=== CBB DERIVATIVE ODDS ACQUISITION ===\n\n")

remaining <- check_quota()

cat("\nAvailable market groups:\n")
cat("  MAIN_MARKETS        - Main H1/H2 markets (recommended start)\n")
cat("  ALL_ML_MARKETS      - Moneylines (2-way and 3-way)\n")
cat("  ALL_SPREADS_MARKETS - Spreads (main and alternate)\n")
cat("  ALL_TOTALS_MARKETS  - Totals (main and alternate)\n")
cat("  ALL_TEAM_TOTALS_MARKETS - Team totals (main and alternate)\n")
cat("  ALL_MARKETS         - All markets combined (expensive!)\n\n")

cat("Usage examples:\n")
cat("  # Fetch main markets (spreads, totals, MLs, team totals) for all games:\n")
cat("  all_odds <- acquire_derivative_odds()\n\n")
cat("  # Fetch only H1/H2 spreads:\n")
cat("  odds <- acquire_derivative_odds(markets = MARKETS_SPREADS)\n\n")
cat("  # Test with 10 games:\n")
cat("  test <- acquire_derivative_odds(max_games = 10)\n\n")
cat("  # Append new markets to existing data:\n")
cat("  odds <- acquire_derivative_odds(markets = ALL_TOTALS_MARKETS, append = TRUE)\n\n")
cat("  # Only games from 2024 season:\n")
cat("  odds <- acquire_derivative_odds(start_date = as.Date('2023-11-01'), end_date = as.Date('2024-04-15'))\n")
