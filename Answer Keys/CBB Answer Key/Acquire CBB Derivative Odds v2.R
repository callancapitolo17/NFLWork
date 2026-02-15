# Acquire CBB Derivative Odds - v2 (Uses existing event IDs)
# Much faster: Uses event IDs we already have from cbb_closing_odds

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
  stop("ODDS_API_KEY environment variable not set")
}

SPORT_KEY <- "basketball_ncaab"
DB_PATH <- "cbb.duckdb"

# Main derivative markets (no alternates to save quota)
MARKETS <- c("h2h_h1", "h2h_h2", "spreads_h1", "spreads_h2",
             "totals_h1", "totals_h2", "team_totals_h1", "team_totals_h2")

REGIONS <- "us,us2,eu"
API_DELAY <- 0.5

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

get_historical_event_odds <- function(event_id, commence_time, markets, regions = REGIONS) {
  if (is.character(commence_time)) {
    commence_time <- ymd_hms(commence_time, tz = "UTC")
  }
  snapshot <- format(commence_time - minutes(15), "%Y-%m-%dT%H:%M:%SZ")

  res <- GET(
    url = paste0("https://api.the-odds-api.com/v4/historical/sports/", SPORT_KEY,
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
    return(NULL)
  }

  parsed <- fromJSON(content(res, "text"), flatten = TRUE)

  if (is.null(parsed$data) || length(parsed$data) == 0) {
    return(NULL)
  }

  parsed$data
}

flatten_derivative_odds <- function(odds_data, event_id, home_team, away_team) {
  if (is.null(odds_data$bookmakers) || length(odds_data$bookmakers) == 0) {
    return(tibble())
  }

  tryCatch({
    bookmakers <- odds_data$bookmakers

    expanded <- bookmakers %>%
      rename(bookmaker_key = key, bookmaker_title = title) %>%
      unnest(markets, names_sep = "_") %>%
      rename(market = markets_key) %>%
      unnest(markets_outcomes, names_sep = "_")

    # Add point column if missing
    if (!"markets_outcomes_point" %in% names(expanded)) {
      expanded$markets_outcomes_point <- NA_real_
    }

    result <- expanded %>%
      mutate(
        event_id = event_id,
        home_team = home_team,
        away_team = away_team,
        commence_time = odds_data$commence_time,
        outcome_name = markets_outcomes_name,
        outcome_price = markets_outcomes_price,
        outcome_point = as.numeric(markets_outcomes_point)
      ) %>%
      select(event_id, home_team, away_team, commence_time,
             bookmaker_key, market, outcome_name, outcome_price, outcome_point)

    result
  }, error = function(e) {
    tibble()
  })
}

# =============================================================================
# MAIN ACQUISITION
# =============================================================================

cat("===========================================\n")
cat("CBB DERIVATIVE ODDS ACQUISITION (v2)\n")
cat("===========================================\n\n")

# Check quota
res <- GET("https://api.the-odds-api.com/v4/sports", query = list(apiKey = api_key))
cat("API Quota Remaining:", headers(res)$`x-requests-remaining`, "\n\n")

# Get unique events from existing data (post May 2023 for derivatives)
con <- dbConnect(duckdb(), dbdir = DB_PATH)

events <- dbGetQuery(con, "
  SELECT DISTINCT
    id as event_id,
    home_team,
    away_team,
    commence_time
  FROM cbb_closing_odds
  WHERE commence_time >= '2023-11-01'
  ORDER BY commence_time
")

cat("Found", nrow(events), "unique events with existing IDs\n")
cat("Markets:", paste(MARKETS, collapse = ", "), "\n")
cat("Estimated quota:", nrow(events) * length(MARKETS) * 3 * 10, "units\n\n")

# Check for already acquired data
if (dbExistsTable(con, "cbb_derivative_closing_odds")) {
  existing <- dbGetQuery(con, "SELECT DISTINCT event_id FROM cbb_derivative_closing_odds")
  events <- events %>% filter(!event_id %in% existing$event_id)
  cat("Skipping", nrow(existing), "already acquired events\n")
  cat("Remaining:", nrow(events), "events to fetch\n\n")
}

if (nrow(events) == 0) {
  cat("All events already acquired!\n")
  dbDisconnect(con)
  stop("Nothing to do")
}

# Fetch odds
all_odds <- list()

for (i in 1:nrow(events)) {
  event <- events[i, ]

  if (i %% 25 == 0) {
    cat("Progress:", i, "/", nrow(events),
        "(", event$home_team, "vs", event$away_team, ")\n")
  }

  odds_data <- get_historical_event_odds(
    event_id = event$event_id,
    commence_time = event$commence_time,
    markets = MARKETS
  )

  if (!is.null(odds_data)) {
    flat <- flatten_derivative_odds(odds_data, event$event_id,
                                    event$home_team, event$away_team)
    if (nrow(flat) > 0) {
      all_odds[[length(all_odds) + 1]] <- flat
    }
  }

  Sys.sleep(API_DELAY)

  # Checkpoint every 100 games
  if (i %% 100 == 0 && length(all_odds) > 0) {
    combined <- bind_rows(all_odds)

    if (dbExistsTable(con, "cbb_derivative_closing_odds")) {
      dbAppendTable(con, "cbb_derivative_closing_odds", combined)
    } else {
      dbWriteTable(con, "cbb_derivative_closing_odds", combined)
    }

    cat("Checkpoint:", nrow(combined), "new rows saved\n")
    all_odds <- list()  # Clear buffer
  }
}

# Save final batch
if (length(all_odds) > 0) {
  combined <- bind_rows(all_odds)
  if (dbExistsTable(con, "cbb_derivative_closing_odds")) {
    dbAppendTable(con, "cbb_derivative_closing_odds", combined)
  } else {
    dbWriteTable(con, "cbb_derivative_closing_odds", combined)
  }
  cat("Final batch:", nrow(combined), "rows saved\n")
}

# Summary
total <- dbGetQuery(con, "SELECT COUNT(*) as n FROM cbb_derivative_closing_odds")
unique_events <- dbGetQuery(con, "SELECT COUNT(DISTINCT event_id) as n FROM cbb_derivative_closing_odds")
cat("\n===========================================\n")
cat("ACQUISITION COMPLETE\n")
cat("Total rows:", total$n, "\n")
cat("Unique events:", unique_events$n, "\n")
cat("===========================================\n")

dbDisconnect(con)
