library(httr)
library(jsonlite)
library(dplyr)
library(purrr)
library(lubridate)
library(nflfastR)
nfl_sched <- fast_scraper_schedules(2020:2025)
nfl_dates  <- nfl_sched %>% 
  mutate(
    kickoff_et  = ymd_hm(paste(gameday, gametime), tz = "America/New_York"),
    kickoff_utc = with_tz(kickoff_et, "UTC"),
    kickoff_date = as.Date(kickoff_utc)) %>% 
  pull(kickoff_date) %>% 
  unique() %>% 
  sort()

# Loop over each day to collect event info
event_list <- map_dfr(as.Date(nfl_dates), function(day) {
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
  
  if (length(parsed) == 0) return(tibble())
  
  as_tibble(parsed) %>%
    mutate(snapshot_date = day)
})
event_list_clean <- event_list$data %>%
  select(id, home_team, away_team, commence_time) %>%
  distinct(id, .keep_all = TRUE) %>%
  mutate(commence_time = ymd_hms(commence_time, tz = "UTC"))


# Function to get odds for a single event via eventIds param
get_event_odds_by_id <- function(event_id, commence_time) {
  snapshot <- format(commence_time - minutes(9360), "%Y-%m-%dT%H:%M:%SZ")
  
  res <- GET(
    url = "https://api.the-odds-api.com/v4/historical/sports/americanfootball_nfl/odds",
    query = list(
      apiKey     = Sys.getenv("ODDS_API_KEY"),
      date       = snapshot,
      eventIds   = event_id,
      regions    = "us,us2,eu,us_ex",
      markets    = "h2h,totals,spreads",
      oddsFormat = "american",
      dateFormat = "iso"
    )
  )
  
  if (http_status(res)$category != "Success") {
    warning(paste("Request failed for", event_id, "at snapshot:", snapshot))
    return(tibble())
  }
  
  parsed <- fromJSON(content(res, as = "text"), flatten = TRUE)
  
  if (length(parsed$data) == 0) {
    return(tibble())
  }
  as_tibble(parsed$data)
}


# Apply this function across all events
nfl_history_df <- map2_dfr(
  event_list_clean$id,
  event_list_clean$commence_time,
  get_event_odds_by_id
)

library(tidyr)

nfl_clean_history_df <- nfl_history_df %>% 
  as_tibble() %>% 
  filter(map_lgl(bookmakers, ~ length(.x) > 0))

nfl_flat_odds <- nfl_clean_history_df %>%
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


