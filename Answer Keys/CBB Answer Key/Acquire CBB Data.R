setwd("~/NFLWork/Answer Keys")
options(warn = 1)

# Package loading with auto-install
packages <- c("data.table", "duckdb", "dplyr", "purrr", "lubridate", "DBI", "httr", "jsonlite", "tidyverse")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE))
}

source("Tools.R")

# Configuration
SPORT_KEY <- "basketball_ncaab"
MARKETS <- "totals,spreads"
REGIONS <- "us,eu,us2"
DB_PATH <- "cbb.duckdb"
TABLE_NAME <- "cbb_closing_odds"
CHECKPOINT_INTERVAL <- 50  # Save to DB every N games
REQUEST_TIMEOUT <- 30  # Timeout in seconds for API requests
MAX_RETRIES <- 3  # Max retries per request

# Wrapper for GET with timeout and retry
safe_GET <- function(url, query, retries = MAX_RETRIES) {
  for (attempt in 1:retries) {
    res <- tryCatch({
      GET(url, query = query, timeout(REQUEST_TIMEOUT))
    }, error = function(e) {
      if (attempt < retries) {
        message(sprintf("    Request failed (attempt %d/%d): %s. Retrying...", attempt, retries, e$message))
        Sys.sleep(2)  # Wait before retry
      }
      return(NULL)
    })
    if (!is.null(res)) return(res)
  }
  return(NULL)
}

# Helper: Get historical events for a date
get_cbb_events_for_date <- function(date) {
  snapshot <- paste0(format(date, "%Y-%m-%d"), "T11:59:59Z")
  res <- safe_GET(
    url = paste0("https://api.the-odds-api.com/v4/historical/sports/", SPORT_KEY, "/events"),
    query = list(
      apiKey = Sys.getenv("ODDS_API_KEY"),
      date = snapshot,
      dateFormat = "iso"
    )
  )

  if (is.null(res) || http_status(res)$category != "Success") {
    warning(paste("Failed to get events for", date))
    return(tibble())
  }

  parsed <- fromJSON(content(res, "text"), flatten = TRUE)
  if (length(parsed) == 0 || is.null(parsed$data) || length(parsed$data) == 0) {
    return(tibble())
  }

  as_tibble(parsed$data) %>%
    select(id, home_team, away_team, commence_time) %>%
    distinct(id, .keep_all = TRUE) %>%
    mutate(
      commence_time = ymd_hms(commence_time, tz = "UTC"),
      snapshot_date = date
    )
}

# Helper: Fetch odds for a single event (spreads + totals only)
get_cbb_event_odds <- function(event_id, commence_time) {
  snapshot <- format(commence_time - minutes(15), "%Y-%m-%dT%H:%M:%SZ")

  res <- safe_GET(
    url = paste0("https://api.the-odds-api.com/v4/historical/sports/", SPORT_KEY, "/odds"),
    query = list(
      apiKey = Sys.getenv("ODDS_API_KEY"),
      date = snapshot,
      eventIds = event_id,
      regions = REGIONS,
      markets = MARKETS,
      oddsFormat = "american",
      dateFormat = "iso"
    )
  )

  if (is.null(res) || http_status(res)$category != "Success") {
    warning(paste("Request failed for", event_id, "at snapshot:", snapshot))
    return(tibble())
  }

  parsed <- fromJSON(content(res, as = "text"), flatten = TRUE)
  if (length(parsed$data) == 0) {
    return(tibble())
  }

  as_tibble(parsed$data)
}

# Flatten odds to one row per game per bookmaker
flatten_cbb_odds <- function(history_df) {
  clean_df <- history_df %>%
    as_tibble() %>%
    filter(map_lgl(bookmakers, ~ length(.x) > 0))

  if (nrow(clean_df) == 0) return(tibble())

  flat_odds <- clean_df %>%
    unnest_longer(bookmakers) %>%
    unnest_wider(bookmakers, names_sep = "_") %>%
    unnest_longer(bookmakers_markets) %>%
    unnest_wider(bookmakers_markets, names_sep = "_") %>%
    unnest_longer(bookmakers_markets_outcomes) %>%
    unnest_wider(bookmakers_markets_outcomes, names_sep = "_") %>%
    mutate(
      commence_time = ymd_hms(commence_time, tz = "UTC"),
      bookmaker_update = ymd_hms(bookmakers_last_update, tz = "UTC")
    ) %>%
    rename(
      bookmaker_key = bookmakers_key,
      bookmaker_title = bookmakers_title,
      market_key = bookmakers_markets_key,
      outcome_name = bookmakers_markets_outcomes_name,
      closing_odds = bookmakers_markets_outcomes_price
    ) %>%
    unnest_wider(outcome_name, names_sep = "_") %>%
    unnest_wider(closing_odds, names_sep = "_") %>%
    unnest_wider(bookmakers_markets_outcomes_point, names_sep = "_") %>%
    mutate(
      market_type = case_when(
        outcome_name_1 == "Over" ~ "totals",
        !is.na(bookmakers_markets_outcomes_point_1) ~ "spread",
        TRUE ~ "other"
      ),
      home_odds = case_when(
        market_type == "totals" ~ NA_real_,
        market_type == "spread" & home_team == outcome_name_1 ~ closing_odds_1,
        market_type == "spread" & away_team == outcome_name_1 ~ closing_odds_2,
        TRUE ~ NA_real_
      ),
      away_odds = case_when(
        market_type == "totals" ~ NA_real_,
        market_type == "spread" & away_team == outcome_name_1 ~ closing_odds_1,
        market_type == "spread" & home_team == outcome_name_1 ~ closing_odds_2,
        TRUE ~ NA_real_
      ),
      home_spread = case_when(
        market_type != "spread" ~ NA_real_,
        home_team == outcome_name_1 ~ bookmakers_markets_outcomes_point_1,
        TRUE ~ bookmakers_markets_outcomes_point_2
      ),
      away_spread = case_when(
        market_type != "spread" ~ NA_real_,
        away_team == outcome_name_1 ~ bookmakers_markets_outcomes_point_1,
        TRUE ~ bookmakers_markets_outcomes_point_2
      )
    ) %>%
    group_by(id, commence_time, home_team, away_team, bookmaker_key, bookmaker_title, bookmaker_update) %>%
    summarise(
      home_spread = first(home_spread[market_type == "spread"]),
      away_spread = first(away_spread[market_type == "spread"]),
      spread_home_odds = first(home_odds[market_type == "spread"]),
      spread_away_odds = first(away_odds[market_type == "spread"]),
      total_line = first(bookmakers_markets_outcomes_point_1[market_type == "totals"]),
      tot_over_odds = first(closing_odds_1[market_type == "totals"]),
      tot_under_odds = first(closing_odds_2[market_type == "totals"]),
      .groups = "drop"
    )

  flat_odds
}

# Main acquisition function
acquire_cbb_odds <- function(start_date, end_date, resume = TRUE) {
  message(sprintf("Acquiring CBB odds from %s to %s", start_date, end_date))

  # Connect to DB
  con <- dbConnect(duckdb(), dbdir = DB_PATH)
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

  # Get existing records for deduplication
  existing_ids <- tryCatch({
    dbGetQuery(con, sprintf("SELECT DISTINCT id || bookmaker_key as unique_id FROM %s", TABLE_NAME)) %>%
      pull(unique_id)
  }, error = function(e) character(0))

  message(sprintf("Found %d existing records", length(existing_ids)))

  # Generate date sequence
  dates <- seq(as.Date(start_date), as.Date(end_date), by = "day")

  all_odds <- tibble()
  games_processed <- 0

  for (i in seq_along(dates)) {
    date <- dates[i]
    message(sprintf("[%d/%d] Processing %s...", i, length(dates), date))

    # Get events for this date
    events <- get_cbb_events_for_date(date)

    if (nrow(events) == 0) {
      message(sprintf("  No events found for %s", date))
      next
    }

    message(sprintf("  Found %d events", nrow(events)))

    # Fetch odds for each event
    for (j in seq_len(nrow(events))) {
      event <- events[j, ]

      # Skip if already in DB (check if any existing record starts with this event ID)
      if (any(startsWith(existing_ids, event$id))) {
        next
      }

      odds_raw <- get_cbb_event_odds(event$id, event$commence_time)

      if (nrow(odds_raw) == 0) {
        next
      }

      odds_flat <- flatten_cbb_odds(odds_raw)

      if (nrow(odds_flat) > 0) {
        all_odds <- bind_rows(all_odds, odds_flat)
        games_processed <- games_processed + 1
      }

      # Checkpoint: save to DB periodically
      if (games_processed > 0 && games_processed %% CHECKPOINT_INTERVAL == 0) {
        message(sprintf("  Checkpoint: saving %d records to DB", nrow(all_odds)))

        # Filter out duplicates before writing
        new_odds <- all_odds %>%
          mutate(unique_id = paste0(id, bookmaker_key)) %>%
          filter(!(unique_id %in% existing_ids)) %>%
          select(-unique_id)

        if (nrow(new_odds) > 0) {
          dbWriteTable(con, TABLE_NAME, new_odds, append = TRUE)
          existing_ids <- c(existing_ids, paste0(new_odds$id, new_odds$bookmaker_key))
        }

        all_odds <- tibble()  # Clear buffer
      }

      # Small delay to avoid rate limiting
      Sys.sleep(0.1)
    }
  }

  # Final write of remaining records
  if (nrow(all_odds) > 0) {
    new_odds <- all_odds %>%
      mutate(unique_id = paste0(id, bookmaker_key)) %>%
      filter(!(unique_id %in% existing_ids)) %>%
      select(-unique_id)

    if (nrow(new_odds) > 0) {
      message(sprintf("Final write: %d records", nrow(new_odds)))
      dbWriteTable(con, TABLE_NAME, new_odds, append = TRUE)
    }
  }

  # Report
  total_records <- dbGetQuery(con, sprintf("SELECT COUNT(*) as n FROM %s", TABLE_NAME))$n
  message(sprintf("Complete. Total records in %s: %d", TABLE_NAME, total_records))
}

# Run for a test date range (uncomment to execute)
# acquire_cbb_odds("2025-01-15", "2025-01-15")

# Full backfill by season (uncomment to execute)
# acquire_cbb_odds("2020-11-01", "2021-04-15")  # 2020-21 season
# acquire_cbb_odds("2021-11-01", "2022-04-15")  # 2021-22 season
# acquire_cbb_odds("2022-11-01", "2023-04-15")  # 2022-23 season
# acquire_cbb_odds("2023-11-01", "2024-04-15")  # 2023-24 season
# acquire_cbb_odds("2024-11-01", "2025-04-15")  # 2024-25 season

# ============================================================================
# PLAY-BY-PLAY DATA ACQUISITION
# ============================================================================

# Load hoopR for CBB play-by-play data
if (!requireNamespace("hoopR", quietly = TRUE)) install.packages("hoopR")
library(hoopR)

# Acquire PBP data and compute game margins by half
acquire_cbb_pbp <- function(seasons = 2021:2026) {
  message(sprintf("Acquiring CBB PBP data for seasons: %s", paste(seasons, collapse = ", ")))

  # Load PBP data from hoopR
  message("Loading PBP data from hoopR (this may take a few minutes)...")
  pbp <- load_mbb_pbp(seasons)
  message(sprintf("Loaded %d play-by-play records", nrow(pbp)))

  # Compute scores by half for each game
  # Half 1 and 2 are regulation; 3+ are overtime periods
  message("Computing game margins by half...")

  game_by_half <- pbp %>%
    filter(half <= 2) %>%  # Only regulation halves
    group_by(game_id, half) %>%
    summarize(
      home_team = first(home_team_name),
      away_team = first(away_team_name),
      game_date = as.Date(first(game_date)),
      # Get max score at end of each half
      home_score_end = max(home_score, na.rm = TRUE),
      away_score_end = max(away_score, na.rm = TRUE),
      # Get min score at start of each half (for half-specific scoring)
      home_score_start = min(home_score, na.rm = TRUE),
      away_score_start = min(away_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      # Points scored in this half specifically
      home_half_score = home_score_end - home_score_start,
      away_half_score = away_score_end - away_score_start,
      # Margin for this half
      half_margin = home_half_score - away_half_score,
      half_total = home_half_score + away_half_score
    )

  # Pivot to wide format - one row per game
  game_margins <- game_by_half %>%
    select(game_id, game_date, home_team, away_team, half, half_margin, half_total,
           home_score_end, away_score_end) %>%
    pivot_wider(
      id_cols = c(game_id, game_date, home_team, away_team),
      names_from = half,
      values_from = c(half_margin, half_total, home_score_end, away_score_end),
      names_glue = "{.value}_h{half}"
    ) %>%
    # Compute full game stats (from end of half 2)
    mutate(
      game_home_margin_h1 = half_margin_h1,
      game_home_margin_h2 = half_margin_h2,
      game_home_margin_fg = home_score_end_h2 - away_score_end_h2,
      game_total_h1 = half_total_h1,
      game_total_h2 = half_total_h2,
      game_total_fg = home_score_end_h2 + away_score_end_h2,
      home_final_score = home_score_end_h2,
      away_final_score = away_score_end_h2,
      home_winner = case_when(
        home_final_score > away_final_score ~ 1L,
        home_final_score < away_final_score ~ 0L,
        TRUE ~ NA_integer_
      ),
      # For Answer Key: actual cover/over indicators will be computed when joining with odds
      actual_cover = NA_integer_,
      actual_over = NA_integer_
    ) %>%
    select(
      game_id, game_date, home_team, away_team,
      game_home_margin_h1, game_home_margin_h2, game_home_margin_fg,
      game_total_h1, game_total_h2, game_total_fg,
      home_final_score, away_final_score, home_winner
    )

  message(sprintf("Computed margins for %d games", nrow(game_margins)))

  # Connect to DB and check existing data
  con <- dbConnect(duckdb(), dbdir = DB_PATH)
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

  # Get existing game IDs for deduplication
  existing_game_ids <- tryCatch({
    dbGetQuery(con, "SELECT DISTINCT game_id FROM cbb_pbp") %>% pull(game_id)
  }, error = function(e) character(0))

  # Filter to new games only
  new_games <- game_margins %>%
    filter(!(game_id %in% existing_game_ids))

  message(sprintf("Found %d existing games, %d new games to add",
                  length(existing_game_ids), nrow(new_games)))

  if (nrow(new_games) > 0) {
    dbWriteTable(con, "cbb_pbp", new_games, append = TRUE)
    message(sprintf("Added %d new games to cbb_pbp", nrow(new_games)))
  }

  # Report final count
  total_games <- dbGetQuery(con, "SELECT COUNT(*) as n FROM cbb_pbp")$n
  message(sprintf("Complete. Total games in cbb_pbp: %d", total_games))
}

# Run PBP acquisition (uncomment to execute)
# acquire_cbb_pbp(2021:2026)
