# 1) Install + load only what you need
for (pkg in c("httr","jsonlite","dplyr","tidyr","lubridate")) {
  if (!requireNamespace(pkg, quietly=TRUE)) install.packages(pkg)
}
library(httr)
library(jsonlite)
library(dplyr)
library(tidyr)
library(lubridate)
library(oddsapiR)

# 2) Build the “closing snapshot” timestamp for yesterday at 23:59:59 UTC
snapshot <- paste0(format(Sys.Date() - 1, "%Y-%m-%d"), "T23:59:59Z")

# 3) Hit the Odds API historical endpoint
res <- GET(
  "https://api.the-odds-api.com/v4/historical/sports/baseball_mlb/odds",
  query = list(
    apiKey      = Sys.getenv("ODDS_API_KEY"),  # make sure you’ve exported this
    date        = snapshot,
    regions     = "us",
    markets     = "h2h,totals",                # moneyline + totals
    oddsFormat  = "american",
    dateFormat  = "iso"
  )
)
stop_for_status(res)

# 4) Parse & flatten the nested JSON
history_df <- fromJSON(content(res, "text"), flatten = TRUE) %>%
  as_tibble()

library(dplyr)
library(tidyr)
library(lubridate)

flat_odds <- history_df$data %>%
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
    total_line     = bookmakers_markets_outcomes_point_1[market_type == "totals"],
    tot_over_odds  = closing_odds_1[market_type == "totals"],
    tot_under_odds = closing_odds_2[market_type == "totals"],
    .groups = "drop"
  )
  

