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

# Loop over each day to collect event info
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
  snapshot <- format(commence_time - minutes(15), "%Y-%m-%dT%H:%M:%SZ")
  
  res <- GET(
    url = "https://api.the-odds-api.com/v4/historical/sports/baseball_mlb/odds",
    query = list(
      apiKey     = Sys.getenv("ODDS_API_KEY"),
      date       = snapshot,
      eventIds   = event_id,
      regions    = "us,eu",
      markets    = "h2h,totals",
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
history_df <- map2_dfr(
  event_list_clean$id,
  event_list_clean$commence_time,
  get_event_odds_by_id
)

library(tidyr)

clean_history_df <- history_df %>% 
  as_tibble() %>% 
  filter(map_lgl(bookmakers, ~ length(.x) > 0))

# clean_history_df <- history_df %>%
#   # mutate(bookmakers = map(bookmakers, ~ as.list(as.data.frame(.x)))) %>% 
#   slice_head(n = 1069) %>% 
#   as_tibble()

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
    ml_home_odds   = first(closing_odds_1[market_type == "moneyline"]),
    ml_away_odds   = first(closing_odds_2[market_type == "moneyline"]),
    total_line = if_else(bookmakers_markets_outcomes_point_1[market_type=="totals"] >0,
                          first(bookmakers_markets_outcomes_point_1[market_type=="totals"]),
                          NA),
    tot_over_odds  = first(closing_odds_1[market_type == "totals"]),
    tot_under_odds = first(closing_odds_2[market_type == "totals"]),
    .groups = "drop"
  )

write.csv(flat_odds,"MLB Flat Odds.csv")
flat_odds <- read.csv("MLB Flat Odds.csv")
# Define the weights for each bookmaker to create consensus total line
book_weights <- tibble::tibble(
  bookmaker_key = c(
    "pinnacle", "circasports", "betonlineag", "lowvig", "matchbook", "bookmaker.eu", "coolbet",
    "everygame", "intertops", "unibet", "unibet_us", "unibet_it", "betmgm", "williamhill_us",
    "caesars", "fanduel", "draftkings", "pointsbetus", "betrivers", "bovada", "sport888", "superbook",
    "wynnbet", "sugarhouse", "betus", "mybookieag", "foxbet", "barstool", "fanatics", "gtbets",
    "tipico_de", "twinspires", "livescorebet_eu", "nordicbet", "betsson", "onexbet"
  ),
  weight = c(
    1.00, 0.95, 0.92, 0.90, 0.88, 0.88, 0.85, 0.80, 0.80, 0.80, 0.75, 0.75, 0.72, 0.72,
    0.70, 0.68, 0.68, 0.65, 0.60, 0.55, 0.52, 0.52, 0.50, 0.50, 0.48, 0.45, 0.45, 0.45, 0.40, 0.40,
    0.38, 0.38, 0.35, 0.35, 0.35, 0.38
  )
)
#Calculate probability from american odds
odds_to_prob <- function(odds) {
  ifelse(odds > 0, 100 / (odds + 100), -odds / (-odds + 100))
}
clean_flat_odds <- flat_odds %>%
  left_join(book_weights, by = "bookmaker_key") %>% 
  filter(if_all(c(ml_home_odds, ml_away_odds, tot_over_odds, tot_under_odds), ~ .x > -400)) %>% 
  filter(tot_over_odds > -200, tot_under_odds > -200) %>%  #having a total that high doesn't make sense could make to maybe lower
  mutate(
    prob_home = odds_to_prob(ml_home_odds),
    prob_away = odds_to_prob(ml_away_odds),
    prob_over = odds_to_prob(tot_over_odds),
    prob_under = odds_to_prob(tot_under_odds)
  ) %>% 
  group_by(id) %>% #as of now not de-vigging
  filter(
    between(
      prob_home, quantile(prob_home, 0.10, na.rm = TRUE), quantile(prob_home, 0.90, na.rm = TRUE) # removing moneyline outliers
    ),
      between(
        prob_over, quantile(prob_over, 0.10, na.rm = TRUE), quantile(prob_over, 0.90, na.rm = TRUE) # removing total outliers
      ))
consensus_ml <- clean_flat_odds %>% 
  group_by(id,commence_time) %>% 
  summarize(
    home_team = first(home_team),
    away_team = first(away_team),
    consensus_prob_home = median(prob_home, na.rm = TRUE),
    consensus_prob_away = median(prob_away, na.rm = TRUE)
  ) %>%
  ungroup()

 consensus_over <- clean_flat_odds %>% 
  group_by(id,total_line) %>% 
  mutate(total_weight = sum(weight, na.rm = T)) %>% 
  ungroup() %>% 
  group_by(id) %>% 
  filter(total_weight == max(total_weight, na.rm = TRUE)) %>% 
  group_by(id,total_line) %>% 
  summarize(consensus_over = median(prob_over, na.rm = T),
            consensus_under = median(prob_under,na.rm = T))
mlb_betting_history <- consensus_ml %>% inner_join(consensus_over, by = "id")  %>% 
  # mutate(start_time = with_tz(ymd_hms(commence_time, tz = "UTC"),tzone = "America/Los_Angeles"))
  mutate(commence_time = ymd_hms(commence_time, tz = "UTC")) %>% 
  mutate(game_date = as.Date(commence_time))

library(fuzzyjoin)

# test <- mlb_betting_history %>% difference_inner_join(game_by_inning %>% mutate(game_date = as.Date(game_date)), by = c("game_date", "home_team","away_team","commence_time" = "game_start_time"),
#                                                       max_dist = c(0, 0,0,as.difftime(1, units = "hours")),
#                                                       distance_col = "time_diff")  
library(data.table)
odds_dt <- as.data.table(mlb_betting_history)
inning_dt <- as.data.table(game_by_inning)

# Make sure commence_time and game_start_time are POSIXct!
odds_dt[, commence_time := ymd_hms(commence_time, tz = "UTC")]
inning_dt[, game_start_time := ymd_hms(game_start_time, tz = "UTC")]

# Standardize team names if needed (trimws, toupper, etc)
# odds_dt[, home_team := toupper(trimws(home_team))]
# inning_dt[, home_team := toupper(trimws(home_team))]

# Set keys for both tables
setkey(odds_dt, home_team, away_team, commence_time)
setkey(inning_dt, home_team, away_team, game_start_time)

# Rolling join: match to the nearest game_start_time within 1 hour (3600 seconds)
joined_dt <- inning_dt[odds_dt, 
                       on = .(home_team, away_team, game_start_time = commence_time), 
                       roll = "nearest",  # closest time (could use "forward" or "backward" for direction)
                       nomatch = 0L
]
