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

tryCatch(
      {
    #Acquire New Data----
    con <- dbConnect(duckdb(), dbdir = "pbp.duckdb") #ensure in NFL work directory

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

    #this needs to go after odds are pulled
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
    sched <- map_dfr(2020:year(Sys.Date()), mlb_schedule) #grab MLB schedule
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
      mutate(new_date = ifelse(date %in% mlb_db_dates, FALSE, TRUE)) %>% #new addition to focus just on new dates
      filter(new_date == TRUE) %>% #new addition to focus just on new dates
      pull(date) %>%
      unique() %>%
      sort()

    if (length(mlb_dates) == 0L) {
      message("No new completed MLB dates to process. Exiting.")
      quit(save = "no", status = 0)
    }

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


    # Apply this function across all events
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
      # 1. one row per bookie
      unnest_longer(bookmakers) %>%
      unnest_wider(bookmakers, names_sep = "_") %>%
      # 2. one row per market (h2h, totals, …)
      unnest_longer(bookmakers_markets) %>%
      unnest_wider(bookmakers_markets, names_sep = "_") %>%
      # 3. one row per single outcome (e.g. “Detroit Tigers” @ –113)
      unnest_longer(bookmakers_markets_outcomes) %>%
      unnest_wider(bookmakers_markets_outcomes, names_sep = "_") %>%
      # 4. parse your datetimes
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
        total_line = if_else(
          bookmakers_markets_outcomes_point_1[market_type == "totals"] > 0,
          first(bookmakers_markets_outcomes_point_1[market_type == "totals"]),
          NA
        ),
        tot_over_odds = first(closing_odds_1[market_type == "totals"]),
        tot_under_odds = first(closing_odds_2[market_type == "totals"]),
        .groups = "drop"
      )

    # Define the weights for each bookmaker to create consensus total line
    odds_to_prob <- function(odds) {
      ifelse(odds > 0, 100 / (odds + 100), -odds / (-odds + 100))
    }

    #need to figure out how to insert into DB Might Need to Join?
    new_betting_history <- flat_odds %>%
      filter(commence_time < with_tz(Sys.time(), "UTC")) %>%
      mutate(unique_book_id = paste0(id, bookmaker_key)) %>%
      mutate(
        new_id = ifelse(unique_book_id %in% betting_db_ids, FALSE, TRUE)
      ) %>%
      filter(new_id == TRUE)
    con <- dbConnect(duckdb(), dbdir = "pbp.duckdb")

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

    dbDisconnect(con, shutdown = TRUE)
  },
  error = function(e) {
    message("ERROR: ", conditionMessage(e))
    quit(save = "no", status = 1)
  }
)

#Acquire PBP Data----
con <- dbConnect(duckdb(), dbdir = "pbp.duckdb") #ensure in correct working directory
db_game_dates <- dbGetQuery(
  con,
  "
  SELECT DISTINCT game_pk,game_date
  FROM mlb_pbp
  WHERE SEASON >= 2020
"
) %>%
  mutate(db_unique_id = paste0(game_pk, year(as.Date(game_date)))) %>%
  pull(game_date) %>%
  unique()
dbDisconnect(con, shutdown = TRUE)

game_dates <- map_dfr(year(Sys.Date()), mlb_schedule) %>%
  filter(
    (!series_description %in% c("Exhibition", "Spring Training")),
    status_coded_game_state == "F"
  ) %>%
  mutate(date = as.Date(ymd_hms(game_date, tz = "UTC"))) %>%
  select(date, game_pk) #grab MLB schedule

new_game_ids <- game_dates %>%
  mutate(
    new_id = ifelse(as.character(date) %in% db_game_dates, FALSE, TRUE)
  ) %>%
  filter(new_id == TRUE) %>%
  pull(game_pk)

season_pbp <- map_dfr(new_game_ids, get_pbp_mlb)

cols_name <- colnames(dbGetQuery(con, "SELECT * FROM mlb_pbp LIMIT 1"))

new_pbp <- season_pbp %>%
  select(cols_name)

con <- dbConnect(duckdb(), dbdir = "pbp.duckdb")

dbWithTransaction(con, {
  dbAppendTable(con, "mlb_pbp", new_pbp) # adds rows
})
dbDisconnect(con, shutdown = TRUE)


