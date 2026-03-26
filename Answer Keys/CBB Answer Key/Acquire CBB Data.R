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
# CLI: --daily mode (dynamic gap fill)
# ============================================================================
cli_args <- commandArgs(trailingOnly = TRUE)
if ("--daily" %in% cli_args) {
  message("Running in --daily mode: filling gap from last acquired date to yesterday")

  con <- dbConnect(duckdb(), dbdir = DB_PATH)
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

  last_date <- tryCatch({
    dbGetQuery(con, sprintf("SELECT MAX(DATE(commence_time)) as d FROM %s", TABLE_NAME))$d
  }, error = function(e) NA)

  dbDisconnect(con, shutdown = TRUE)

  if (is.na(last_date)) {
    # No data yet — start from current season
    today <- Sys.Date()
    start_date <- if (as.integer(format(today, "%m")) >= 11) {
      as.Date(paste0(format(today, "%Y"), "-11-01"))
    } else {
      as.Date(paste0(as.integer(format(today, "%Y")) - 1, "-11-01"))
    }
  } else {
    start_date <- as.Date(last_date)
  }

  end_date <- Sys.Date() - 1

  if (start_date > end_date) {
    message("Already up to date — nothing to fetch.")
  } else {
    acquire_cbb_odds(as.character(start_date), as.character(end_date))
  }

  # Exit after daily odds run — don't load hoopR section below
  quit(save = "no", status = 0)
}

# ============================================================================
# CLI: --daily-pbp mode (hoopR PBP to cbb_pbp_v2, matching cbbpy schema)
# ============================================================================
# Helper: extract first-to-10 from raw hoopR PBP data
# Returns tibble with game_id + first_to_10_h1 (1 = home, 0 = away, NA = neither reached 10)
extract_first_to_10_h1 <- function(pbp) {
  # Filter to H1 plays with valid running scores
  h1 <- pbp %>%
    filter(half == 1, !is.na(home_score), !is.na(away_score)) %>%
    arrange(game_id, id) %>%
    group_by(game_id) %>%
    summarize(
      # Find first play where home reaches 10
      home_first_10_idx = {
        idx <- which(home_score >= 10)
        if (length(idx) > 0) min(idx) else NA_integer_
      },
      # Find first play where away reaches 10
      away_first_10_idx = {
        idx <- which(away_score >= 10)
        if (length(idx) > 0) min(idx) else NA_integer_
      },
      .groups = "drop"
    ) %>%
    mutate(
      first_to_10_h1 = case_when(
        is.na(home_first_10_idx) & is.na(away_first_10_idx) ~ NA_integer_,
        is.na(away_first_10_idx) ~ 1L,   # only home reached 10
        is.na(home_first_10_idx) ~ 0L,   # only away reached 10
        home_first_10_idx < away_first_10_idx ~ 1L,   # home reached first
        away_first_10_idx < home_first_10_idx ~ 0L,   # away reached first
        TRUE ~ NA_integer_  # tied (same play) — treat as NA
      ),
      game_id = as.character(game_id)
    ) %>%
    select(game_id, first_to_10_h1)

  h1
}

# Helper: extract race-to-X for full game (any threshold)
# Returns tibble with game_id + first_to_{threshold}_fg (1 = home, 0 = away, NA = neither/tie)
extract_race_to_fg <- function(pbp, threshold = 10) {
  col_name <- paste0("first_to_", threshold, "_fg")
  fg <- pbp %>%
    filter(!is.na(home_score), !is.na(away_score)) %>%  # NO half filter — full game
    arrange(game_id, id) %>%
    group_by(game_id) %>%
    summarize(
      home_first_idx = { idx <- which(home_score >= threshold); if (length(idx) > 0) min(idx) else NA_integer_ },
      away_first_idx = { idx <- which(away_score >= threshold); if (length(idx) > 0) min(idx) else NA_integer_ },
      .groups = "drop"
    ) %>%
    mutate(
      !!col_name := case_when(
        is.na(home_first_idx) & is.na(away_first_idx) ~ NA_integer_,
        is.na(away_first_idx) ~ 1L,
        is.na(home_first_idx) ~ 0L,
        home_first_idx < away_first_idx ~ 1L,
        away_first_idx < home_first_idx ~ 0L,
        TRUE ~ NA_integer_  # simultaneous (tie)
      ),
      game_id = as.character(game_id)
    ) %>%
    select(game_id, all_of(col_name))
  fg
}

if ("--daily-pbp" %in% cli_args) {
  if (!requireNamespace("hoopR", quietly = TRUE)) install.packages("hoopR")
  library(hoopR)

  message("Running in --daily-pbp mode: filling PBP gap using hoopR")

  # Determine current season (Nov = new season starts)
  today <- Sys.Date()
  current_season <- if (as.integer(format(today, "%m")) >= 11) {
    as.integer(format(today, "%Y")) + 1L
  } else {
    as.integer(format(today, "%Y"))
  }

  con <- dbConnect(duckdb(), dbdir = DB_PATH)
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

  # Check existing data in cbb_pbp_v2
  existing_ids <- tryCatch({
    dbGetQuery(con, "SELECT DISTINCT game_id FROM cbb_pbp_v2") %>% pull(game_id)
  }, error = function(e) character(0))
  message(sprintf("Existing games in cbb_pbp_v2: %d", length(existing_ids)))

  # Load PBP from hoopR for current season
  message(sprintf("Loading hoopR PBP for season %d...", current_season))
  pbp <- load_mbb_pbp(current_season)
  message(sprintf("Loaded %d PBP records", nrow(pbp)))

  # Compute per-half scores (regulation halves 1 & 2)
  game_by_half <- pbp %>%
    filter(half <= 2) %>%
    group_by(game_id, half) %>%
    summarize(
      home_team = paste(first(home_team_name), first(home_team_mascot)),
      away_team = paste(first(away_team_name), first(away_team_mascot)),
      game_date = as.character(as.Date(first(game_date))),
      home_score_end = max(home_score, na.rm = TRUE),
      away_score_end = max(away_score, na.rm = TRUE),
      home_score_start = min(home_score, na.rm = TRUE),
      away_score_start = min(away_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      home_half_score = home_score_end - home_score_start,
      away_half_score = away_score_end - away_score_start
    )

  # Compute OT info
  ot_games <- pbp %>%
    filter(half > 2) %>%
    group_by(game_id) %>%
    summarize(
      went_to_ot = 1L,
      ot_home_end = max(home_score, na.rm = TRUE),
      ot_away_end = max(away_score, na.rm = TRUE),
      .groups = "drop"
    )

  # End of regulation scores (end of half 2)
  reg_end <- game_by_half %>%
    filter(half == 2) %>%
    select(game_id, reg_home_end = home_score_end, reg_away_end = away_score_end)

  ot_info <- ot_games %>%
    left_join(reg_end, by = "game_id") %>%
    mutate(
      home_ot_score = as.integer(ot_home_end - reg_home_end),
      away_ot_score = as.integer(ot_away_end - reg_away_end),
      game_home_margin_ot = as.integer(home_ot_score - away_ot_score),
      game_total_ot = as.integer(home_ot_score + away_ot_score)
    ) %>%
    select(game_id, went_to_ot, home_ot_score, away_ot_score, game_home_margin_ot, game_total_ot)

  # Extract first-to-10 from raw PBP before pivoting
  first10 <- extract_first_to_10_h1(pbp)
  message(sprintf("First-to-10 H1: computed for %d games", nrow(first10)))

  # Extract race-to-X full game for thresholds 10, 20, 40
  race10_fg <- extract_race_to_fg(pbp, 10)
  race20_fg <- extract_race_to_fg(pbp, 20)
  race40_fg <- extract_race_to_fg(pbp, 40)
  message(sprintf("Race-to-X FG: computed for %d games", nrow(race10_fg)))

  # Pivot to wide format - one row per game
  game_margins <- game_by_half %>%
    select(game_id, game_date, home_team, away_team, half,
           home_half_score, away_half_score, home_score_end, away_score_end) %>%
    pivot_wider(
      id_cols = c(game_id, game_date, home_team, away_team),
      names_from = half,
      values_from = c(home_half_score, away_half_score, home_score_end, away_score_end),
      names_glue = "{.value}_h{half}"
    ) %>%
    # Final scores = end of regulation if no OT, or max overall
    left_join(ot_info, by = "game_id") %>%
    left_join(first10, by = "game_id") %>%
    left_join(race10_fg, by = "game_id") %>%
    left_join(race20_fg, by = "game_id") %>%
    left_join(race40_fg, by = "game_id") %>%
    mutate(
      home_h1_score = as.integer(home_half_score_h1),
      away_h1_score = as.integer(away_half_score_h1),
      game_home_margin_h1 = as.integer(home_h1_score - away_h1_score),
      game_total_h1 = as.integer(home_h1_score + away_h1_score),
      home_h2_score = as.integer(home_half_score_h2),
      away_h2_score = as.integer(away_half_score_h2),
      game_home_margin_h2 = as.integer(home_h2_score - away_h2_score),
      game_total_h2 = as.integer(home_h2_score + away_h2_score),
      # OT defaults
      went_to_ot = coalesce(went_to_ot, 0L),
      home_ot_score = coalesce(home_ot_score, 0L),
      away_ot_score = coalesce(away_ot_score, 0L),
      game_home_margin_ot = coalesce(game_home_margin_ot, 0L),
      game_total_ot = coalesce(game_total_ot, 0L),
      first_to_10_h1 = coalesce(first_to_10_h1, NA_integer_),
      first_to_10_fg = coalesce(first_to_10_fg, NA_integer_),
      first_to_20_fg = coalesce(first_to_20_fg, NA_integer_),
      first_to_40_fg = coalesce(first_to_40_fg, NA_integer_),
      # Final scores (regulation end + OT)
      home_final_score = as.integer(home_score_end_h2 + home_ot_score),
      away_final_score = as.integer(away_score_end_h2 + away_ot_score),
      game_home_margin_fg = as.integer(home_final_score - away_final_score),
      game_total_fg = as.integer(home_final_score + away_final_score),
      home_winner = case_when(
        home_final_score > away_final_score ~ 1L,
        home_final_score < away_final_score ~ 0L,
        TRUE ~ NA_integer_
      ),
      game_id = as.character(game_id)
    ) %>%
    select(
      game_id, game_date, home_team, away_team,
      home_h1_score, away_h1_score, game_home_margin_h1, game_total_h1,
      home_h2_score, away_h2_score, game_home_margin_h2, game_total_h2,
      home_ot_score, away_ot_score, game_home_margin_ot, game_total_ot,
      home_final_score, away_final_score, game_home_margin_fg, game_total_fg,
      home_winner, went_to_ot, first_to_10_h1,
      first_to_10_fg, first_to_20_fg, first_to_40_fg
    )

  # Filter to completed games only (must have both halves with valid scores)
  # Without this, in-progress games get inserted with NAs and are never corrected
  # (PRIMARY KEY dedup would skip them on future runs)
  game_margins <- game_margins %>%
    filter(!is.na(home_h2_score) & !is.na(home_final_score))

  # Filter to new games only (dedup against existing cbb_pbp_v2)
  new_games <- game_margins %>%
    filter(!(game_id %in% existing_ids))

  message(sprintf("Computed %d games, %d are new", nrow(game_margins), nrow(new_games)))

  if (nrow(new_games) > 0) {
    # Ensure table exists
    tryCatch(
      dbGetQuery(con, "SELECT 1 FROM cbb_pbp_v2 LIMIT 0"),
      error = function(e) {
        dbExecute(con, "CREATE TABLE cbb_pbp_v2 (
          game_id VARCHAR PRIMARY KEY, game_date VARCHAR,
          home_team VARCHAR, away_team VARCHAR,
          home_h1_score INTEGER, away_h1_score INTEGER,
          game_home_margin_h1 INTEGER, game_total_h1 INTEGER,
          home_h2_score INTEGER, away_h2_score INTEGER,
          game_home_margin_h2 INTEGER, game_total_h2 INTEGER,
          home_ot_score INTEGER, away_ot_score INTEGER,
          game_home_margin_ot INTEGER, game_total_ot INTEGER,
          home_final_score INTEGER, away_final_score INTEGER,
          game_home_margin_fg INTEGER, game_total_fg INTEGER,
          home_winner INTEGER, went_to_ot INTEGER,
          first_to_10_h1 INTEGER,
          first_to_10_fg INTEGER, first_to_20_fg INTEGER, first_to_40_fg INTEGER
        )")
      }
    )
    dbWriteTable(con, "cbb_pbp_v2", new_games, append = TRUE)
    message(sprintf("Added %d new games to cbb_pbp_v2", nrow(new_games)))
  }

  total <- dbGetQuery(con, "SELECT COUNT(*) as n FROM cbb_pbp_v2")$n
  message(sprintf("Complete. Total games in cbb_pbp_v2: %d", total))

  quit(save = "no", status = 0)
}

# ============================================================================
# CLI: --backfill-first10 mode (add first_to_10_h1 to existing cbb_pbp_v2 rows)
# ============================================================================
if ("--backfill-first10" %in% cli_args) {
  if (!requireNamespace("hoopR", quietly = TRUE)) install.packages("hoopR")
  library(hoopR)

  message("Running in --backfill-first10 mode: adding first_to_10_h1 to existing games")

  # Use --db <path> if provided, otherwise default to parent Answer Keys/cbb.duckdb
  # (the main PBP database, not the local CBB Answer Key copy)
  db_idx <- which(cli_args == "--db")
  backfill_db <- if (length(db_idx) > 0 && db_idx < length(cli_args)) {
    cli_args[db_idx + 1]
  } else {
    DB_PATH  # Answer Keys/cbb.duckdb (same as DB_PATH after setwd)
  }
  message(sprintf("Database: %s", backfill_db))

  con <- dbConnect(duckdb(), dbdir = backfill_db)
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

  # Add column if it doesn't exist
  tryCatch(
    dbExecute(con, "ALTER TABLE cbb_pbp_v2 ADD COLUMN first_to_10_h1 INTEGER"),
    error = function(e) message("Column first_to_10_h1 already exists")
  )

  # Find games missing first_to_10_h1
  missing <- dbGetQuery(con, "
    SELECT game_id FROM cbb_pbp_v2 WHERE first_to_10_h1 IS NULL
  ") %>% pull(game_id)
  message(sprintf("Games missing first_to_10_h1: %d", length(missing)))

  if (length(missing) == 0) {
    message("Nothing to backfill.")
    quit(save = "no", status = 0)
  }

  # Load PBP for all relevant seasons
  seasons <- 2021:2026
  message(sprintf("Loading hoopR PBP for seasons %s...", paste(seasons, collapse = ", ")))
  pbp <- load_mbb_pbp(seasons)
  message(sprintf("Loaded %d PBP records", nrow(pbp)))

  # Extract first-to-10
  first10 <- extract_first_to_10_h1(pbp)
  message(sprintf("Computed first-to-10 for %d games", nrow(first10)))

  # Update matching rows
  to_update <- first10 %>% filter(game_id %in% missing, !is.na(first_to_10_h1))
  message(sprintf("Updating %d games with first_to_10_h1 values", nrow(to_update)))

  if (nrow(to_update) > 0) {
    # Write to temp table, then UPDATE
    dbWriteTable(con, "tmp_first10", to_update, overwrite = TRUE)
    updated <- dbExecute(con, "
      UPDATE cbb_pbp_v2 SET first_to_10_h1 = t.first_to_10_h1
      FROM tmp_first10 t WHERE cbb_pbp_v2.game_id = t.game_id
    ")
    dbExecute(con, "DROP TABLE IF EXISTS tmp_first10")
    message(sprintf("Updated %d rows", updated))
  }

  # Report
  filled <- dbGetQuery(con, "
    SELECT COUNT(*) as n FROM cbb_pbp_v2 WHERE first_to_10_h1 IS NOT NULL
  ")$n
  total <- dbGetQuery(con, "SELECT COUNT(*) as n FROM cbb_pbp_v2")$n
  message(sprintf("Complete. %d/%d games have first_to_10_h1", filled, total))

  quit(save = "no", status = 0)
}

# ============================================================================
# CLI: --backfill-race-fg <threshold> mode (add first_to_{N}_fg to existing cbb_pbp_v2 rows)
# ============================================================================
if ("--backfill-race-fg" %in% cli_args) {
  if (!requireNamespace("hoopR", quietly = TRUE)) install.packages("hoopR")
  library(hoopR)

  threshold <- as.integer(cli_args[which(cli_args == "--backfill-race-fg") + 1])
  if (is.na(threshold)) stop("Usage: --backfill-race-fg <threshold>  (e.g. 10, 20, 40)")
  col_name <- paste0("first_to_", threshold, "_fg")

  message(sprintf("Running in --backfill-race-fg mode: adding %s to existing games", col_name))

  # Use --db <path> if provided, otherwise default
  db_idx <- which(cli_args == "--db")
  backfill_db <- if (length(db_idx) > 0 && db_idx < length(cli_args)) {
    cli_args[db_idx + 1]
  } else {
    DB_PATH
  }
  message(sprintf("Database: %s", backfill_db))

  con <- dbConnect(duckdb(), dbdir = backfill_db)
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE)

  # Add column if missing
  tryCatch(
    dbExecute(con, sprintf("ALTER TABLE cbb_pbp_v2 ADD COLUMN %s INTEGER", col_name)),
    error = function(e) message(sprintf("Column %s already exists", col_name))
  )

  # Find games with NULL
  missing <- dbGetQuery(con, sprintf(
    "SELECT DISTINCT game_id FROM cbb_pbp_v2 WHERE %s IS NULL", col_name
  ))
  message(sprintf("Found %d games missing %s", nrow(missing), col_name))

  if (nrow(missing) == 0) {
    message("Nothing to backfill")
    quit(save = "no")
  }

  # Load PBP for each season and extract
  all_results <- tibble()
  for (yr in 2021:2026) {
    message(sprintf("Loading hoopR PBP for %d...", yr))
    tryCatch({
      pbp <- load_mbb_pbp(yr)
      if (nrow(pbp) > 0) {
        result <- extract_race_to_fg(pbp, threshold)
        # Only keep games that are in our missing list
        result <- result %>% filter(game_id %in% missing$game_id)
        all_results <- bind_rows(all_results, result)
        message(sprintf("  %d: %d games extracted", yr, nrow(result)))
      }
    }, error = function(e) message(sprintf("  %d: skipped (%s)", yr, e$message)))
  }

  # Update via temp table
  if (nrow(all_results) > 0) {
    dbWriteTable(con, "tmp_race_fg", all_results, overwrite = TRUE)
    updated <- dbExecute(con, sprintf("
      UPDATE cbb_pbp_v2 SET %s = t.%s
      FROM tmp_race_fg t WHERE cbb_pbp_v2.game_id = t.game_id
    ", col_name, col_name))
    dbExecute(con, "DROP TABLE IF EXISTS tmp_race_fg")
    message(sprintf("Updated %d rows for %s", updated, col_name))
  }

  # Summary
  filled <- dbGetQuery(con, sprintf(
    "SELECT COUNT(*) n FROM cbb_pbp_v2 WHERE %s IS NOT NULL", col_name
  ))
  total <- dbGetQuery(con, "SELECT COUNT(*) n FROM cbb_pbp_v2")
  message(sprintf("Coverage: %d/%d games have %s (%.1f%%)",
                  filled$n, total$n, col_name, 100 * filled$n / total$n))
  quit(save = "no")
}

# ============================================================================
# PLAY-BY-PLAY DATA ACQUISITION (hoopR - legacy function for manual use)
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
      home_team = paste(first(home_team_name), first(home_team_mascot)),
      away_team = paste(first(away_team_name), first(away_team_mascot)),
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
