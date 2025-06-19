library(httr)
res <- GET(
  url = "https://api.the-odds-api.com/v4/historical/sports/basketball_nba/events/da359da99aa27e97d38f2df709343998/odds",
  query = list(
    apiKey  = Sys.getenv("ODDS_API_KEY"),
    date    = "2023-11-29T22:45:00Z",
    regions = "us",
    markets = "h2h"
  )
)
stop_for_status(res)
content(res, "parsed")

library(httr)
library(jsonlite)
library(purrr)
library(dplyr)
library(tidyr)
library(lubridate)

api_key  <- Sys.getenv("ODDS_API_KEY")
sport    <- "baseball_mlb"
regions  <- "us"
markets  <- "h2h"                # moneyline only
snapshot <- "2023-08-29T23:59:59Z"

# 1) Fetch the list of event IDs for that snapshot
evt_url <- sprintf("https://api.the-odds-api.com/v4/historical/sports/%s/events", sport)
evt_res <- GET(evt_url, query = list(apiKey = api_key, date = snapshot))
stop_for_status(evt_res)

# Parse into a tibble so we can pull the IDs easily
events_df <- fromJSON(content(evt_res, "text"), flatten=TRUE) %>% as_tibble()

# Extract the vector of IDs
event_ids <- pull(events_df, id)

# 2) Define the odds‐fetching function
odds_url <- sprintf("https://api.the-odds-api.com/v4/historical/sports/%s/events/%%s/odds", sport)
fetch_odds <- function(eid) {
  res <- GET(sprintf(odds_url, eid), query = list(
    apiKey     = api_key,
    date       = snapshot,
    regions    = regions,
    markets    = markets,
    oddsFormat = "american",
    dateFormat = "iso"
  ))
  stop_for_status(res)
  fromJSON(content(res, "text"), flatten=TRUE) %>% 
    as_tibble() 
}

# 3) Map over the ID vector (no set_names needed)
mlb_829 <- map_dfr(event_ids, fetch_odds, .id = "game_id")

# 4) Flatten to one row per outcome
mlb_829_flat <- mlb_829 %>%
  unnest(markets) %>%
  unnest(outcomes) %>%
  transmute(
    game_id      = game_id,
    home         = home_team,
    away         = away_team,
    commence     = ymd_hms(commence_time, tz="UTC"),
    market       = markets.key,       # always "h2h" here
    outcome      = outcomes.name,     # "Home" or "Away"
    closing_odds = outcomes.price
  )

# Inspect
print(mlb_829_flat)
