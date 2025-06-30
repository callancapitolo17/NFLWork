library(httr)
library(jsonlite)
library(dplyr)
library(purrr)
library(lubridate)
library(baseballr)
sched <- map_dfr(2020:2025, mlb_schedule) 
mlb_dates  <- sched %>% 
  mutate(date = as.Date(date)) %>% 
  filter(!series_description %in% c("Exhibition","Spring Training") & status_coded_game_state == "F") %>% 
  filter(date >= as.Date("2020-06-06")) %>% 
  pull(date) %>% 
  unique() %>% 
  sort()


history_df <- map_dfr(mlb_dates, function(day) {
  snapshot <- paste0(format(day, "%Y-%m-%d"), "T23:59:59Z")
  res <- GET(
    "https://api.the-odds-api.com/v4/historical/sports/baseball_mlb/odds",
    query = list(
      apiKey     = Sys.getenv("ODDS_API_KEY"),
      date       = snapshot,
      regions    = "us,eu",
      markets    = "h2h,totals",
      oddsFormat = "american",
      dateFormat = "iso"
    )
  )
  stop_for_status(res) #make sure web request works or else throw an error
  parsed <- fromJSON(content(res, "text"), flatten = TRUE) #put the json into an oject in R that works and I can work with
  
  # if there were no games that day, just return an empty tibble
  if (length(parsed$data) == 0) {
    return(tibble())
  }
  
  # otherwise pull out the 'data' element and coerce that to a tibble
  as_tibble(parsed$data) %>%
    mutate(snapshot_date = day)  # tag each row with its snapshot date
})

clean_history_df <- history_df %>% 
  as_tibble()
flat_odds <- clean_history_df %>%
  # 1. one row per bookie
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
    market_type = if_else(outcome_name_1 == "Over", "totals", "moneyline")
  ) %>% 
  group_by(
    id, commence_time, home_team, away_team,
    bookmaker_key, bookmaker_title, bookmaker_update
  ) %>%
  summarise(
    ml_home_odds   = closing_odds_1[market_type == "moneyline"],
    ml_away_odds   = closing_odds_2[market_type == "moneyline"],
    ttotal_line = if_else(any(market_type=="totals"),
                          bookmakers_markets_outcomes_point_1[market_type=="totals"],
                          NA_real_),
    tot_over_odds  = closing_odds_1[market_type == "totals"],
    tot_under_odds = closing_odds_2[market_type == "totals"],
    .groups = "drop"
  )
