#TODO generate single sample answer key for each game and then use functions generates bets based on that
setwd("~/NFLWork/Answer Keys")
library(data.table)
library(oddsapiR)
library(fuzzyjoin)
library(duckdb)
library(dplyr)
library(purrr)
library(lubridate)
library(DBI)
library(httr)
library(jsonlite)
library(tidyverse)
source("Tools.R")
library(data.table)

# Answer Key----
con <- dbConnect(duckdb(), dbdir = "pbp.duckdb") 

DT <- dbGetQuery(con,"select * from nfl_betting_pbp") %>% 
  rename(
    home_spread_odds = "consensus_devig_home_odds",
    away_spread_odds = "consensus_devig_away_odds",
    over_odds = "consensus_devig_over_odds",
    under_odds = "consensus_devig_under_odds"
  ) %>%
  mutate(actual_over = ifelse(total_final_score > total_line, 1, 0),
         actual_cover = ifelse(home_margin > -home_spread, 1, 0)) %>% 
  as.data.table()


disp <- compute_dispersion(DT,moneyline =FALSE)
ss <- disp$ss
st <- disp$st


# get game odds----
game_odds <- toa_sports_odds(
  sport_key = "americanfootball_nfl",
  regions = "us,us2,eu,us_ex",
  markets = "spreads,totals",
  odds_format = "american",
  date_format = "iso"
)

book_weights <- dbGetQuery(con, "
  SELECT *
  FROM nfl_weights")
dbDisconnect(con)

consensus_spread <- prepare_two_way_odds(
  game_odds    = game_odds,
  mkt_key      = "spreads",
  book_weights = book_weights,
  prob_fun     = devig_american,
  prob_names   = c("prob_home", "prob_away"),
  odds_names = c("outcomes_price_home", "outcomes_price_away")
) %>% 
  select(-outcomes_point_away) %>% 
  rename(spread = outcomes_point_home) %>% 
  pick_consensus_line(
    game_id_col = "id",
    line_col = "spread",
    weight_col = "spread_weight",
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
  odds_names = c("outcomes_price_Over", "outcomes_price_Under")
) %>%
  select(-outcomes_point_Under) %>% 
  rename(total_line = outcomes_point_Over) %>% 
  pick_consensus_line(
    game_id_col = "id",
    line_col = "total_line",
    weight_col = "totals_weight",
    date_col = "date",
    time_col = "commence_time",
    market1 = "prob_over",
    market2 = "prob_under"
  )


nfl_odds <- consensus_spread %>% inner_join(consensus_total %>% ungroup() %>%  select(-home_team,-away_team,-date,-commence_time), by = "id") %>% 
  filter(if_all(everything(), ~ !is.na(.)))

# Generate Side predictions----
targets <- nfl_odds %>%
  transmute(
    id,
    parent_spread = spread,
    parent_total  = total_line,
    target_cover  = consensus_prob_home,
    target_over   = consensus_prob_over
  )


# 1) Get upcoming MLB events (you can pass regions here, but markets aren’t needed yet)
period     <- 1
bankroll   <- 50
kelly_mult <- 0.25
N <-  round(nrow(DT)*0.05,0)
tol_error <- 1

events <- get_events("americanfootball_nfl", regions = "us")

# Fetch ALL moneyline derivative markets at once
ml_results <- build_multi_moneyline_markets(
  DT              = DT,
  consensus_odds  = nfl_odds,
  ss              = ss,
  st              = st,
  N               = N,
  periods         = c(1, 2, 3, 4),  # Quarters 1-4
  events          = events,
  markets         = c("h2h_q1", "h2h_q2", "h2h_q3", "h2h_q4"),  # Corresponding markets
  sport_key       = "americanfootball_nfl",
  bankroll        = bankroll,
  kelly_mult      = kelly_mult,
  use_spread_line = TRUE,
  targets         = targets,
  margin_col      = "game_home_margin_period"
)

# View all bets across all markets
my_bets <- ml_results$bets

# View market summary to see which quarters have most edge
ml_results$markets_summary

# Filter to specific markets if needed
q1_bets <- my_bets %>% filter(market == "h2h_q1")
q2_bets <- my_bets %>% filter(market == "h2h_q2")

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
