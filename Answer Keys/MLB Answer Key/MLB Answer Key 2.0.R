setwd("~/NFLWork/Answer Keys")
library(data.table)
library(oddsapiR)
library(fuzzyjoin)
library(duckdb)
library(dplyr)
library(purrr)
library(baseballr)
library(lubridate)
library(DBI)
library(httr)
library(jsonlite)
library(tidyverse)
source("Tools.R")
library(data.table)

# Answer Key----
con <- dbConnect(duckdb(), dbdir = "pbp.duckdb") 

DT <- dbGetQuery(con,"select * from mlb_betting_pbp") %>% 
  rename(
    home_ml_odds = "consensus_devig_home_odds",
    away_ml_odds = "consensus_devig_away_odds",
    over_odds = "consensus_devig_over_odds",
    under_odds = "consensus_devig_under_odds",
    actual_cover = "home_winner"
  ) %>%
  mutate(actual_over = ifelse(total_final_score > total_line, 1, 0)) %>% 
  as.data.table()

dbDisconnect(con)


disp <- compute_dispersion(DT)
ss <- disp$ss
st <- disp$st


# get game odds----
game_odds <- toa_sports_odds(
  sport_key = "baseball_mlb",
  regions = "us,eu",
  markets = "h2h,totals",
  odds_format = "american",
  date_format = "iso"
)
con <- dbConnect(duckdb(), dbdir = "pbp.duckdb") 
book_weights <- dbGetQuery(con, "
  SELECT *
  FROM mlb_weights")
dbDisconnect(con)

ml_consensus <- game_odds %>%
  filter(market_key == "h2h") %>%
  mutate(date = as.Date(commence_time)) %>% 
  mutate(odds_type = ifelse(home_team == outcomes_name, "home", ifelse(away_team == outcomes_name, "away", outcomes_name))) %>%
  pivot_wider(id_cols = c(id, commence_time, home_team, away_team, bookmaker_key,outcomes_point,date), names_from = odds_type, values_from = outcomes_price, names_prefix = "odds_") %>%
  mutate(as_tibble(american_prob(odds_home, odds_away))) %>% # consider devigging which can be done by switching devig_american()
  rename(prob_home = p1, prob_away = p2) %>%
  left_join(book_weights, by = "bookmaker_key") 

ml_consensus <- prepare_two_way_odds(
  game_odds    = game_odds,
  mkt_key      = "h2h",
  book_weights = book_weights,
  prob_fun     = devig_american,  # or devig_american if you want
  prob_names   = c("prob_home", "prob_away"),
  odds_names = c("odds_home", "odds_away")
) %>% 
  moneyline_consensus(
    game_id_col = "id",
    weight_col = "ml_weight",
    date_col = "date",
    time_col = "commence_time",
    market1 = "prob_home",
    market2 = "prob_away"
  )

consensus_total <- prepare_two_way_odds(
    game_odds    = game_odds,
    mkt_key      = "totals",
    book_weights = book_weights,
    prob_fun     = devig_american,
    prob_names   = c("prob_over", "prob_under"),
    odds_names = c("odds_Over", "odds_Under")
  ) %>%
  rename(total_line = outcomes_point) %>% 
  pick_consensus_line(
    game_id_col = "id",
    line_col = "total_line",
    weight_col = "tot_weight",
    date_col = "date",
    time_col = "commence_time",
    market1 = "prob_over",
    market2 = "prob_under"
  )
  

mlb_odds <- ml_consensus %>% inner_join(consensus_total, by = "id")

# Generate Side predictions----
targets <- mlb_odds %>%
  transmute(
    id,
    parent_spread = consensus_prob_home,
    parent_total  = total_line,
    target_cover  = consensus_prob_home,
    target_over   = consensus_over
  )

final_bets <- targets %>%
  mutate(res = pmap(
    list(id, parent_spread, parent_total, target_cover, target_over),
    ~ run_means_for_id(..1, ..2, ..3, ..4, ..5,
      DT = DT, ss = ss, st = st, N = N
    )
  )) %>%
  unnest(res) %>%
  inner_join(game_odds %>% group_by(id) %>%
    summarize(home_team = first(home_team), away_team = first(away_team), commence_time = first(commence_time)), by = "id") %>%
  select(home_team, away_team, commence_time, everything())

predictions <- final_bets %>%
  mutate(across(starts_with("game_home_margin_in"),
    ~ prob_to_american(.x),
    .names = "{.col}_american"
  )) %>%
  relocate(ends_with("_american"), .before = starts_with("game_home_margin_in"))


# 1) Get upcoming MLB events (you can pass regions here, but markets aren’t needed yet)
period     <- 5
bankroll   <- 200
kelly_mult <- 0.25

events <- get_events("baseball_mlb", regions = "us")

ml_results <- build_moneyline_market(
  DT            = DT,
  consensus_odds = mlb_odds,
  ss            = ss,
  st            = st,
  N             = N,
  period        = period,
  events        = events,
  market        = "h2h_1st_5_innings",
  bankroll      = bankroll,
  kelly_mult    = kelly_mult
)
my_bets <- ml_results$bets

total_results <- build_totals_market(
  DT             = DT,
  consensus_odds = mlb_odds,
  ss             = ss,
  st             = st,
  N              = N,
  period         = period,
  events         = events,
  market         = "alternate_totals_1st_5_innings",
  bankroll       = bankroll,
  kelly_mult     = kelly_mult
)
total_my_bets <- total_results$bets

spread_results <- build_spread_market(
  DT             = DT,
  consensus_odds = mlb_odds,
  ss             = ss,
  st             = st,
  N              = N,
  period         = period,
  events         = events,
  market         = "spreads_1st_5_innings",
  bankroll       = bankroll,
  kelly_mult     = kelly_mult
)
my_spread_bets <- spread_results$bets
