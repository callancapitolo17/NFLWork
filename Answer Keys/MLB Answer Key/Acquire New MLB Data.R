setwd("~/NFLWork/Answer Keys")
options(warn = 1)
    if (!requireNamespace("data.table", quietly = TRUE)) {
      install.packages("data.table")
    }
    suppressPackageStartupMessages(library(
      data.table,
      quietly = TRUE,
      warn.conflicts = FALSE
    ))
    if (!requireNamespace("oddsapiR", quietly = TRUE)) {
      install.packages("oddsapiR")
    }
    suppressPackageStartupMessages(library(
      oddsapiR,
      quietly = TRUE,
      warn.conflicts = FALSE
    ))
    if (!requireNamespace("fuzzyjoin", quietly = TRUE)) {
      install.packages("fuzzyjoin")
    }
    suppressPackageStartupMessages(library(
      fuzzyjoin,
      quietly = TRUE,
      warn.conflicts = FALSE
    ))
    if (!requireNamespace("duckdb", quietly = TRUE)) {
      install.packages("duckdb")
    }
    suppressPackageStartupMessages(library(
      duckdb,
      quietly = TRUE,
      warn.conflicts = FALSE
    ))
    if (!requireNamespace("dplyr", quietly = TRUE)) {
      install.packages("dplyr")
    }
    suppressPackageStartupMessages(library(
      dplyr,
      quietly = TRUE,
      warn.conflicts = FALSE
    ))
    if (!requireNamespace("purrr", quietly = TRUE)) {
      install.packages("purrr")
    }
    suppressPackageStartupMessages(library(
      purrr,
      quietly = TRUE,
      warn.conflicts = FALSE
    ))
    if (!requireNamespace("baseballr", quietly = TRUE)) {
      install.packages("baseballr")
    }
    suppressPackageStartupMessages(library(
      baseballr,
      quietly = TRUE,
      warn.conflicts = FALSE
    ))
    if (!requireNamespace("lubridate", quietly = TRUE)) {
      install.packages("lubridate")
    }
    suppressPackageStartupMessages(library(
      lubridate,
      quietly = TRUE,
      warn.conflicts = FALSE
    ))
    if (!requireNamespace("DBI", quietly = TRUE)) {
      install.packages("DBI")
    }
    suppressPackageStartupMessages(library(
      DBI,
      quietly = TRUE,
      warn.conflicts = FALSE
    ))
    if (!requireNamespace("httr", quietly = TRUE)) {
      install.packages("httr")
    }
    suppressPackageStartupMessages(library(
      httr,
      quietly = TRUE,
      warn.conflicts = FALSE
    ))
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      install.packages("jsonlite")
    }
    suppressPackageStartupMessages(library(
      jsonlite,
      quietly = TRUE,
      warn.conflicts = FALSE
    ))
    source("Tools.R")

# ============================================================================
# CLI: --daily runs odds only, --daily-pbp runs PBP only, no flags runs both
# ============================================================================
cli_args <- commandArgs(trailingOnly = TRUE)
run_odds <- length(cli_args) == 0 || "--daily" %in% cli_args
run_pbp  <- length(cli_args) == 0 || "--daily-pbp" %in% cli_args

# ============================================================================
# Section 1: Acquire Closing Odds from Odds API
# ============================================================================
if (run_odds) {
  message("=== Acquiring closing odds ===")
  tryCatch(
    {
    con <- dbConnect(duckdb(), dbdir = "pbp.duckdb")
    on.exit(tryCatch(dbDisconnect(con, shutdown = TRUE), error = function(e) NULL), add = TRUE)

    mlb_db_dates <- dbGetQuery(
      con,
      "
  SELECT commence_time as date
  FROM mlb_betting_history
"
    ) %>%
      mutate(date = as.Date(date)) %>%
      unique() %>%
      pull()

    betting_db_ids <- dbGetQuery(
      con,
      "
  SELECT id,commence_time, bookmaker_key
  FROM mlb_betting_history
"
    ) %>%
      mutate(unique_db_book_id = paste0(id, bookmaker_key)) %>%
      pull(unique_db_book_id)

    dbDisconnect(con, shutdown = TRUE)
    on.exit(NULL)

    # T-2 safety buffer: only fetch odds for games that finished at least
    # 2 days ago. Guarantees no live/in-progress games get captured even if
    # the schedule API has stale Final flags.
    cutoff_date <- Sys.Date() - 2

    sched <- map_dfr(2020:year(Sys.Date()), mlb_schedule)
    mlb_dates <- sched %>%
      filter(status_abstract_game_state == "Final") %>%
      mutate(
        game_datetime_utc = ymd_hms(game_date, tz = "UTC"),
        date = as.Date(game_datetime_utc)
      ) %>%
      filter(
        !series_description %in% c("Exhibition", "Spring Training") &
          status_coded_game_state == "F"
      ) %>%
      filter(date >= as.Date("2020-06-06")) %>%
      filter(date <= cutoff_date) %>%
      filter(!date %in% mlb_db_dates) %>%
      pull(date) %>%
      unique() %>%
      sort()

    if (length(mlb_dates) == 0L) {
      message("No new completed MLB dates to process.")
    } else {
      event_list <- map_dfr(mlb_dates, function(day) {
        snapshot <- paste0(format(day, "%Y-%m-%d"), "T14:59:59Z")
        res <- GET(
          url = "https://api.the-odds-api.com/v4/historical/sports/baseball_mlb/events",
          query = list(
            apiKey = Sys.getenv("ODDS_API_KEY"),
            date = snapshot,
            dateFormat = "iso"
          )
        )
        stop_for_status(res)
        parsed <- fromJSON(content(res, "text"), flatten = TRUE)

        if (length(parsed) == 0) {
          return(tibble())
        }

        as_tibble(parsed) %>%
          mutate(snapshot_date = day)
      })
      event_list_clean <- event_list$data %>%
        select(id, home_team, away_team, commence_time) %>%
        distinct(id, .keep_all = TRUE) %>%
        mutate(commence_time = ymd_hms(commence_time, tz = "UTC"))

      history_df <- map2_dfr(
        event_list_clean$id,
        event_list_clean$commence_time,
        get_event_odds_by_id
      )

      library(tidyr)

      clean_history_df <- history_df %>%
        as_tibble() %>%
        filter(map_lgl(bookmakers, ~ length(.x) > 0))

      flat_odds <- clean_history_df %>%
        unnest_longer(bookmakers) %>%
        unnest_wider(bookmakers, names_sep = "_") %>%
        unnest_longer(bookmakers_markets) %>%
        unnest_wider(bookmakers_markets, names_sep = "_") %>%
        unnest_longer(bookmakers_markets_outcomes) %>%
        unnest_wider(bookmakers_markets_outcomes, names_sep = "_") %>%
        mutate(
          commence_time = ymd_hms(commence_time, tz = "UTC"),
          bookmaker_update = ymd_hms(bookmakers_last_update, tz = "UTC"),
          market_update = ymd_hms(bookmakers_markets_last_update, tz = "UTC")
        ) %>%
        rename(
          bookmaker_key = bookmakers_key,
          bookmaker_title = bookmakers_title,
          market_key = bookmakers_markets_key,
          outcome_name = bookmakers_markets_outcomes_name,
          closing_odds = bookmakers_markets_outcomes_price
        ) %>%
        unnest_wider(outcome_name, names_sep = c("_")) %>%
        unnest_wider(closing_odds, names_sep = c("_")) %>%
        unnest_wider(bookmakers_markets_outcomes_point, names_sep = c("_")) %>%
        mutate(
          market_type = if_else(outcome_name_1 == "Over", "totals", "moneyline"),
          home_odds = ifelse(
            market_type == "totals",
            NA,
            ifelse(
              market_type == "moneyline" & home_team == outcome_name_1,
              closing_odds_1,
              closing_odds_2
            )
          ),
          away_odds = ifelse(
            market_type == "totals",
            NA,
            ifelse(
              market_type == "moneyline" & away_team == outcome_name_1,
              closing_odds_1,
              closing_odds_2
            )
          ),
        ) %>%
        group_by(
          id,
          commence_time,
          home_team,
          away_team,
          bookmaker_key,
          bookmaker_title,
          bookmaker_update
        ) %>%
        summarise(
          ml_home_odds = first(home_odds[market_type == "moneyline"]),
          ml_away_odds = first(away_odds[market_type == "moneyline"]),
          # Books may post ML without totals (e.g. BetAnything on some games).
          # Empty subscript needs explicit handling — if_else(empty > 0, ...)
          # returns size 0 and crashes summarise's size-1 requirement.
          total_line = {
            pts <- bookmakers_markets_outcomes_point_1[market_type == "totals"]
            if (length(pts) > 0 && !is.na(pts[1]) && pts[1] > 0) pts[1] else NA_real_
          },
          tot_over_odds = first(closing_odds_1[market_type == "totals"]),
          tot_under_odds = first(closing_odds_2[market_type == "totals"]),
          .groups = "drop"
        )

      new_betting_history <- flat_odds %>%
        filter(commence_time < with_tz(Sys.time(), "UTC")) %>%
        mutate(unique_book_id = paste0(id, bookmaker_key)) %>%
        mutate(
          new_id = ifelse(unique_book_id %in% betting_db_ids, FALSE, TRUE)
        ) %>%
        filter(new_id == TRUE)

      con <- dbConnect(duckdb(), dbdir = "pbp.duckdb")
      on.exit(tryCatch(dbDisconnect(con, shutdown = TRUE), error = function(e) NULL), add = TRUE)

      duckdb_register(con, "new_rows", new_betting_history)

      dbExecute(
        con,
        "
INSERT INTO mlb_betting_history AS t (
  id, commence_time, home_team, away_team,
  bookmaker_key, bookmaker_title, bookmaker_update,
  ml_home_odds, ml_away_odds, total_line,
  tot_over_odds, tot_under_odds
)
SELECT
  n.id, n.commence_time, n.home_team, n.away_team,
  n.bookmaker_key, n.bookmaker_title, n.bookmaker_update,
  n.ml_home_odds, n.ml_away_odds, n.total_line,
  n.tot_over_odds, n.tot_under_odds
FROM new_rows n
WHERE NOT EXISTS (
  SELECT 1 FROM mlb_betting_history t2
  WHERE t2.id = n.id AND t2.bookmaker_key = n.bookmaker_key
);
"
      )

      message(sprintf("Inserted %d new odds rows.", nrow(new_betting_history)))
      dbDisconnect(con, shutdown = TRUE)
      on.exit(NULL)
    }
  },
  error = function(e) {
    message("ERROR in odds acquisition: ", conditionMessage(e))
    if (!run_pbp) quit(save = "no", status = 1)
  }
  )
}

# ============================================================================
# Section 2: Acquire PBP Data
# ============================================================================
if (run_pbp) {
  message("=== Acquiring PBP data ===")
  tryCatch(
    {
    con <- dbConnect(duckdb(), dbdir = "pbp.duckdb")
    on.exit(tryCatch(dbDisconnect(con, shutdown = TRUE), error = function(e) NULL), add = TRUE)

    # Cast game_date to Date server-side so comparisons are type-safe regardless
    # of the column's underlying storage format (currently VARCHAR 'YYYY-MM-DD').
    db_game_dates <- dbGetQuery(
      con,
      "
  SELECT DISTINCT CAST(game_date AS DATE) AS game_date
  FROM mlb_pbp_all
  WHERE year(CAST(game_date AS DATE)) >= 2020
"
    ) %>%
      pull(game_date) %>%
      as.Date() %>%
      unique()

    # Get column schema while connection is still open
    cols_name <- colnames(dbGetQuery(con, "SELECT * FROM mlb_pbp_all LIMIT 1"))
    dbDisconnect(con, shutdown = TRUE)
    on.exit(NULL)

    # T-2 safety buffer: only fetch PBP for games that finished at least
    # 2 days ago (mirrors the odds section).
    cutoff_date <- Sys.Date() - 2

    game_dates <- map_dfr(year(Sys.Date()), mlb_schedule) %>%
      filter(
        (!series_description %in% c("Exhibition", "Spring Training")),
        status_coded_game_state == "F"
      ) %>%
      mutate(date = as.Date(ymd_hms(game_date, tz = "UTC"))) %>%
      filter(date <= cutoff_date) %>%
      select(date, game_pk)

    # Both `date` and `db_game_dates` are Date objects — type-safe comparison.
    new_game_ids <- game_dates %>%
      filter(!date %in% db_game_dates) %>%
      pull(game_pk)

    if (length(new_game_ids) == 0L) {
      message("No new MLB games to fetch PBP for.")
    } else {
      message(sprintf("Fetching PBP for %d new games...", length(new_game_ids)))
      season_pbp <- map_dfr(new_game_ids, get_pbp_mlb)

      new_pbp <- season_pbp %>%
        select(all_of(cols_name))

      con <- dbConnect(duckdb(), dbdir = "pbp.duckdb")
      on.exit(tryCatch(dbDisconnect(con, shutdown = TRUE), error = function(e) NULL), add = TRUE)

      dbWithTransaction(con, {
        dbAppendTable(con, "mlb_pbp_all", new_pbp)
      })

      message(sprintf("Inserted %d PBP rows for %d games.", nrow(new_pbp), length(new_game_ids)))
      dbDisconnect(con, shutdown = TRUE)
      on.exit(NULL)
    }
  },
  error = function(e) {
    message("ERROR in PBP acquisition: ", conditionMessage(e))
    quit(save = "no", status = 1)
  }
  )
}

message("=== MLB data acquisition complete ===")
