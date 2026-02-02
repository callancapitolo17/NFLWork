# NFL Answer Key 2.0
# Uses the new architecture: run_answer_key_sample() once per game,
# then generate predictions for multiple market types from the same sample

setwd("~/NFLWork/Answer Keys")
library(data.table)
library(oddsapiR)
library(duckdb)
library(dplyr)
library(purrr)
library(lubridate)
library(DBI)
library(httr)
library(jsonlite)
library(tidyverse)
source("Tools.R")

# =============================================================================
# LOAD HISTORICAL DATA
# =============================================================================

con <- dbConnect(duckdb(), dbdir = "pbp.duckdb")

betting_20_plus <- dbGetQuery(con, "SELECT * FROM nfl_betting_pbp")
pre_20 <- dbGetQuery(con, "SELECT * FROM nfl_pre_20_betting_history") %>%
  rename(home_spread = "spread")

DT <- bind_rows(betting_20_plus, pre_20) %>%
  rename(
    home_spread_odds = "consensus_devig_home_odds",
    away_spread_odds = "consensus_devig_away_odds",
    over_odds = "consensus_devig_over_odds",
    under_odds = "consensus_devig_under_odds"
  ) %>%
  mutate(
    actual_over = ifelse(total_final_score > total_line, 1, 0),
    actual_cover = ifelse(home_margin > -home_spread, 1, 0)
  ) %>%
  as.data.table()

# Compute dispersion for distance weighting
disp <- compute_dispersion(DT, moneyline = FALSE)
ss <- disp$ss
st <- disp$st

# =============================================================================
# GET CURRENT GAME ODDS & BUILD CONSENSUS
# =============================================================================

game_odds <- toa_sports_odds(
  sport_key = "americanfootball_nfl",
  regions = "us,us2,eu,us_ex",
  markets = "spreads,totals",
  odds_format = "american",
  date_format = "iso"
)

book_weights <- dbGetQuery(con, "SELECT * FROM nfl_weights")
dbDisconnect(con)

# Build consensus spread
consensus_spread <- prepare_two_way_odds(
  game_odds    = game_odds,
  mkt_key      = "spreads",
  book_weights = book_weights,
  prob_fun     = devig_american,
  prob_names   = c("prob_home", "prob_away"),
  odds_names   = c("outcomes_price_home", "outcomes_price_away")
) %>%
  select(-outcomes_point_away) %>%
  rename(spread = outcomes_point_home) %>%
  pick_consensus_line(
    game_id_col = "id",
    line_col    = "spread",
    weight_col  = "spread_weight",
    date_col    = "date",
    time_col    = "commence_time",
    market1     = "prob_home",
    market2     = "prob_away"
  )

# Build consensus total
consensus_total <- prepare_two_way_odds(
  game_odds    = game_odds,
  mkt_key      = "totals",
  book_weights = book_weights,
  prob_fun     = devig_american,
  prob_names   = c("prob_over", "prob_under"),
  odds_names   = c("outcomes_price_Over", "outcomes_price_Under")
) %>%
  select(-outcomes_point_Under) %>%
  rename(total_line = outcomes_point_Over) %>%
  pick_consensus_line(
    game_id_col = "id",
    line_col    = "total_line",
    weight_col  = "totals_weight",
    date_col    = "date",
    time_col    = "commence_time",
    market1     = "prob_over",
    market2     = "prob_under"
  )

# Combine into single odds table
nfl_odds <- consensus_spread %>%
  inner_join(
    consensus_total %>% ungroup() %>% select(-home_team, -away_team, -date, -commence_time),
    by = "id"
  ) %>%
  filter(if_all(everything(), ~ !is.na(.)))

# =============================================================================
# SETUP PREDICTION PARAMETERS
# =============================================================================

# Targets for Answer Key algorithm
targets <- nfl_odds %>%
  transmute(
    id,
    parent_spread = spread,
    parent_total  = total_line,
    target_cover  = consensus_prob_home,
    target_over   = consensus_prob_over
  )

# Betting parameters
bankroll   <- 100
kelly_mult <- 0.25
N <- round(nrow(DT) * 0.05, 0)

# Get upcoming NFL events
events <- get_events("americanfootball_nfl", regions = "us")

# =============================================================================
# GENERATE SAMPLES (ONCE PER GAME)
# =============================================================================

# Generate samples for all games upfront - this is the expensive operation
# and should only run ONCE, then be reused for all market types
cat("Generating samples for all games...\n")
samples <- generate_all_samples(
  targets         = targets,
  DT              = DT,
  ss              = ss,
  st              = st,
  N               = N,
  use_spread_line = TRUE
)

# =============================================================================
# GENERATE PREDICTIONS FOR ALL DERIVATIVE MARKETS
# =============================================================================

# Base periods (Q1-Q4 + H1 + H2)
periods_base <- c("1", "2", "3", "4", "Half1", "Half2")

# --- MONEYLINES (Q1-Q4 + H1 + H2) ---
# No alternate lines for moneylines
ml_results <- build_moneylines_from_samples(
  samples         = samples,
  consensus_odds  = nfl_odds,
  events          = events,
  periods         = periods_base,
  markets         = c("h2h_q1", "h2h_q2", "h2h_q3", "h2h_q4", "h2h_h1", "h2h_h2"),
  sport_key       = "americanfootball_nfl",
  bankroll        = bankroll,
  kelly_mult      = kelly_mult,
  margin_col      = "game_home_margin_period"
)

ml_bets <- ml_results$bets
ml_results$markets_summary

# --- 3-WAY MONEYLINES (Q1-Q4 + H1 + H2) ---
# Includes tie as a distinct outcome
ml_3way_results <- build_3way_from_samples(
  samples         = samples,
  consensus_odds  = nfl_odds,
  events          = events,
  periods         = periods_base,
  markets         = c("h2h_3_way_q1", "h2h_3_way_q2", "h2h_3_way_q3", "h2h_3_way_q4", "h2h_3_way_h1", "h2h_3_way_h2"),
  sport_key       = "americanfootball_nfl",
  bankroll        = bankroll,
  kelly_mult      = kelly_mult,
  margin_col      = "game_home_margin_period"
)

ml_3way_bets <- ml_3way_results$bets
ml_3way_results$markets_summary

# --- TOTALS (Q1-Q4 + H1 + H2) + ALTERNATE TOTALS ---
# Include both main lines and alternate lines
totals_markets <- c(
  "totals_q1", "totals_q2", "totals_q3", "totals_q4", "totals_h1", "totals_h2",
  "alternate_totals_q1", "alternate_totals_q2", "alternate_totals_q3", "alternate_totals_q4", "alternate_totals_h1", "alternate_totals_h2"
)
totals_periods <- c(
  "1", "2", "3", "4", "Half1", "Half2",
  "1", "2", "3", "4", "Half1", "Half2"
)

total_results <- build_totals_from_samples(
  samples         = samples,
  consensus_odds  = nfl_odds,
  events          = events,
  periods         = totals_periods,
  markets         = totals_markets,
  sport_key       = "americanfootball_nfl",
  bankroll        = bankroll,
  kelly_mult      = kelly_mult
)

total_bets <- total_results$bets
total_results$markets_summary

# --- SPREADS (Q1-Q4 + H1 + H2) + ALTERNATE SPREADS ---
# Include both main lines and alternate lines
spreads_markets <- c(
  "spreads_q1", "spreads_q2", "spreads_q3", "spreads_q4", "spreads_h1", "spreads_h2",
  "alternate_spreads_q1", "alternate_spreads_q2", "alternate_spreads_q3", "alternate_spreads_q4", "alternate_spreads_h1", "alternate_spreads_h2"
)
spreads_periods <- c(
  "1", "2", "3", "4", "Half1", "Half2",
  "1", "2", "3", "4", "Half1", "Half2"
)

spread_results <- build_spreads_from_samples(
  samples         = samples,
  consensus_odds  = nfl_odds,
  events          = events,
  periods         = spreads_periods,
  markets         = spreads_markets,
  sport_key       = "americanfootball_nfl",
  bankroll        = bankroll,
  kelly_mult      = kelly_mult
)

spread_bets <- spread_results$bets
spread_results$markets_summary

# --- TEAM TOTALS (Q1-Q4 + H1 + H2) + ALTERNATE TEAM TOTALS ---
# Team-specific over/under lines
team_totals_markets <- c(
  "team_totals_q1", "team_totals_q2", "team_totals_q3", "team_totals_q4",
  "team_totals_h1", "team_totals_h2",
  "alternate_team_totals_q1", "alternate_team_totals_q2", "alternate_team_totals_q3", "alternate_team_totals_q4",
  "alternate_team_totals_h1", "alternate_team_totals_h2"
)
team_totals_periods <- c(
  "1", "2", "3", "4", "Half1", "Half2",
  "1", "2", "3", "4", "Half1", "Half2"
)

team_totals_results <- build_team_totals_from_samples(
  samples         = samples,
  consensus_odds  = nfl_odds,
  events          = events,
  periods         = team_totals_periods,
  markets         = team_totals_markets,
  sport_key       = "americanfootball_nfl",
  bankroll        = bankroll,
  kelly_mult      = kelly_mult
)

team_totals_bets <- team_totals_results$bets
team_totals_results$markets_summary

# =============================================================================
# COMBINE ALL BETS
# =============================================================================

all_bets <- bind_rows(
  ml_bets %>% mutate(market_type = "moneyline"),
  ml_3way_bets %>% mutate(market_type = "moneyline_3way"),
  total_bets %>% mutate(market_type = "totals"),
  spread_bets %>% mutate(market_type = "spreads"),
  team_totals_bets %>% mutate(market_type = "team_totals")
) %>%
  arrange(desc(ev))

# Summary by market type
all_bets %>%
  group_by(market_type) %>%
  summarise(
    n_bets = n(),
    total_stake = sum(bet_size),
    avg_ev = mean(ev),
    .groups = "drop"
  )

# Summary by main vs alternate lines
all_bets %>%
  group_by(market_type) %>%
  summarise(
    n_bets = n(),
    total_stake = sum(bet_size),
    avg_ev = mean(ev),
    .groups = "drop"
  )

# Top bets across all markets
all_bets %>%
  head(20)

# =============================================================================
# WAGERZON OFFSHORE ODDS COMPARISON
# =============================================================================

# Run Wagerzon comparison - scrapes fresh odds and compares to model predictions
wagerzon_results <- run_wagerzon_comparison(
  spread_results = spread_results,
  totals_results = total_results,
  ml_results     = ml_results,
  sport          = "nfl",
  bankroll       = bankroll,
  kelly_mult     = kelly_mult,
  ev_threshold   = 0.05,
  scrape_first   = TRUE  # Set to FALSE if odds were recently scraped
)

# Wagerzon-specific bets
wagerzon_bets <- wagerzon_results$bets

if (nrow(wagerzon_bets) > 0) {
  cat("\n=== TOP WAGERZON BETS ===\n")
  print(wagerzon_bets %>% head(10))
}

# =============================================================================
# COMBINE ALL BETS (API + WAGERZON)
# =============================================================================

all_bets_combined <- bind_rows(
  all_bets %>% mutate(source = "api"),
  wagerzon_bets %>% mutate(source = "wagerzon")
) %>%
  arrange(desc(ev))

# Final summary
cat("\n=== FINAL BETTING SUMMARY ===\n")
all_bets_combined %>%
  group_by(source, market_type) %>%
  summarise(
    n_bets = n(),
    total_stake = sum(bet_size),
    avg_ev = mean(ev),
    .groups = "drop"
  ) %>%
  print()
