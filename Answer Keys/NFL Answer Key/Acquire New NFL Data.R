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
library(nflfastR)
library(nflreadr)
library(tidyverse)
source("Tools.R")

tryCatch(
  {
    #Acquire New Data----
    con <- dbConnect(duckdb(), dbdir = "pbp.duckdb") #ensure in NFL work directory
    
    nfl_db_dates <- dbGetQuery(
      con,
      "
  SELECT commence_time as date
  FROM nfl_closing_odds
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
  FROM nfl_closing_odds
"
    ) %>%
      mutate(unique_db_book_id = paste0(id, bookmaker_key)) %>%
      pull(unique_db_book_id)
    
    dbDisconnect(con, shutdown = TRUE)
    nfl_sched <- fast_scraper_schedules(2020:2025)
    nfl_dates  <- nfl_sched %>% 
      mutate(
        kickoff_et  = ymd_hm(paste(gameday, gametime), tz = "America/New_York"),
        kickoff_utc = with_tz(kickoff_et, "UTC"),
        kickoff_date = as.Date(kickoff_utc)) %>% 
      pull(kickoff_date) %>% 
      unique() %>% 
      sort()
    if (length(nfl_dates) == 0L) {
      message("No new completed NFL dates to process. Exiting.")
      quit(save = "no", status = 0)
    }
    
    event_list <- map_dfr(nfl_dates, function(day) {
      snapshot <- paste0(format(day, "%Y-%m-%d"), "T11:59:59Z")
      res <- GET(
        url = "https://api.the-odds-api.com/v4/historical/sports/americanfootball_nfl/events",
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
      
      as_tibble(compact(parsed)) %>%
        mutate(snapshot_date = day)
    })
    event_list_clean <- event_list$data %>%
      select(id, home_team, away_team, commence_time) %>%
      distinct(id, .keep_all = TRUE) %>%
      mutate(commence_time = ymd_hms(commence_time, tz = "UTC"))
    
    
    # Apply this function across all events
    history_df <- map2_dfr(
      event_list_clean$id[1],
      event_list_clean$commence_time[1],
      get_event_odds_by_id,
      url = "https://api.the-odds-api.com/v4/historical/sports/americanfootball_nfl/odds"
    )
    
    library(tidyr)
    
    clean_history_df <- history_df %>%
      as_tibble() %>%
      filter(map_lgl(bookmakers, ~ length(.x) > 0))
    
    flat_odds <- clean_history_df %>%
      unnest_longer(bookmakers) %>% 
      unnest_wider (bookmakers,   names_sep = "_") %>% 
      # 2. one row per market (h2h, totals, …)
      unnest_longer(bookmakers_markets) %>% 
      unnest_wider (bookmakers_markets, names_sep = "_") %>% 
      # 3. one row per single outcome (e.g. “Detroit Tigers” @ –113)
      unnest_longer(bookmakers_markets_outcomes) %>% 
      unnest_wider (bookmakers_markets_outcomes, names_sep = "_") %>% 
      # 4. parse your datetimes
      mutate(
        commence_time       = ymd_hms(commence_time,               tz="UTC"),
        bookmaker_update    = ymd_hms(bookmakers_last_update,      tz="UTC"),
        market_update       = ymd_hms(bookmakers_markets_last_update, tz="UTC")
      ) %>%
      rename(
        bookmaker_key   = bookmakers_key,
        bookmaker_title = bookmakers_title,
        market_key      = bookmakers_markets_key,
        outcome_name    = bookmakers_markets_outcomes_name,
        closing_odds    = bookmakers_markets_outcomes_price
      ) %>% 
      unnest_wider(outcome_name, names_sep = c("_")) %>% 
      unnest_wider(closing_odds, names_sep = c("_")) %>% 
      unnest_wider(bookmakers_markets_outcomes_point, names_sep = c("_")) %>% 
      mutate(
        market_type = if_else(outcome_name_1 == "Over", "totals", ifelse(is.na(bookmakers_markets_outcomes_point_1) & outcome_name_1 != "Over","moneyline","spread" )),
        home_odds = ifelse(market_type == "totals",NA, ifelse(market_type %in% c("moneyline","spread") & home_team == outcome_name_1  , closing_odds_1,closing_odds_2)),
        away_odds = ifelse(market_type == "totals",NA, ifelse(market_type %in% c("moneyline","spread") & away_team == outcome_name_1, closing_odds_1,closing_odds_2)),
        home_spread = ifelse(market_type != "spread",NA,ifelse(market_type == "spread" & home_team == outcome_name_1, bookmakers_markets_outcomes_point_1, bookmakers_markets_outcomes_point_2)),
        away_spread = ifelse(market_type != "spread",NA,ifelse(market_type == "spread" & away_team == outcome_name_1, bookmakers_markets_outcomes_point_1, bookmakers_markets_outcomes_point_2))
      ) %>% 
      group_by(
        id, commence_time, home_team, away_team,
        bookmaker_key, bookmaker_title, bookmaker_update
      ) %>% 
      summarise(
        ml_home_odds = first(home_odds[market_type == "moneyline"]),
        ml_away_odds = first(away_odds[market_type == "moneyline"]),
        home_spread = first(home_spread[market_type == "spread"]),
        away_spread = first(away_spread[market_type == "spread"]),
        spread_home_odds = first(home_odds[market_type == "spread"]),
        spread_away_odds = first(away_odds[market_type == "spread"]),
        total_line = if_else(bookmakers_markets_outcomes_point_1[market_type=="totals"] >0,
                             first(bookmakers_markets_outcomes_point_1[market_type=="totals"]),
                             NA),
        tot_over_odds  = first(closing_odds_1[market_type == "totals"]),
        tot_under_odds = first(closing_odds_2[market_type == "totals"]),
        .groups = "drop"
      )
    
    #need to figure out how to insert into DB Might Need to Join?
    new_betting_history <- flat_odds %>%
      filter(commence_time < with_tz(Sys.time(), "UTC")) %>%
      mutate(unique_book_id = paste0(id, bookmaker_key)) %>%
      mutate(
        new_id = ifelse(unique_book_id %in% betting_db_ids, FALSE, TRUE)
      ) %>%
      filter(new_id == TRUE) %>% 
      select(-unique_book_id, -new_id)
    con <- dbConnect(duckdb(), dbdir = "pbp.duckdb")
    dbWriteTable(con, "nfl_closing_odds", new_betting_history, append = TRUE)

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
  SELECT DISTINCT game_id,game_date
  FROM nfl_pbp
  WHERE SEASON >= 2020
"
) %>%
  mutate(db_unique_id = paste0(game_id, year(as.Date(game_date)))) %>%
  pull(game_date) %>%
  unique()
dbDisconnect(con, shutdown = TRUE)

season_pbp <- load_pbp(max(nfl_sched$season))

new_pbp <- season_pbp %>% 
  mutate(date = as.Date(time_of_day),
         new_id = ifelse(as.character(date) %in% db_game_dates, FALSE, TRUE)) %>% 
  filter(new_id == TRUE)
  
con <- dbConnect(duckdb(), dbdir = "pbp.duckdb")
dbWriteTable(con, "nfl_pbp", new_pbp %>% select(-date,-new_id), append = TRUE)
dbDisconnect(con, shutdown = TRUE)


