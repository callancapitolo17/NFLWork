# NFL Prepare - Part 1 of parallel pipeline
# Runs in parallel with scrapers
# Generates samples and saves to DuckDB for NFLCombine.R

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

cat("=== NFL PREPARE: Starting sample generation ===\n")

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

book_weights <- dbGetQuery(con, "SELECT * FROM nfl_weights")

cat("Historical data loaded.\n")

# =============================================================================
# GET CURRENT GAME ODDS & BUILD CONSENSUS
# =============================================================================

cat("Fetching Odds API...\n")

game_odds <- toa_sports_odds(
  sport_key = "americanfootball_nfl",
  regions = "us,us2,eu,us_ex",
  markets = "spreads,totals",
  odds_format = "american",
  date_format = "iso"
)

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

cat(sprintf("Built consensus for %d games.\n", nrow(nfl_odds)))

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

cat("Generating samples for all games...\n")
samples <- generate_all_samples(
  targets         = targets,
  DT              = DT,
  ss              = ss,
  st              = st,
  N               = N,
  use_spread_line = TRUE
)

cat(sprintf("Generated %d samples.\n", length(samples)))

# =============================================================================
# SAVE TO DUCKDB (shared state, no temp files)
# =============================================================================

cat("Saving to DuckDB...\n")

# Flatten samples list to dataframe for storage
# Each element of samples is a list with $sample (dataframe), $final_N, etc.
# Extract the $sample dataframes and add game_id
samples_df <- samples %>%
  imap_dfr(~ .x$sample %>% mutate(game_id = .y))

# Save to DuckDB tables
dbWriteTable(con, "nfl_samples_temp", samples_df, overwrite = TRUE)
dbWriteTable(con, "nfl_odds_temp", nfl_odds, overwrite = TRUE)
dbWriteTable(con, "nfl_events_temp", events, overwrite = TRUE)

# Also save betting parameters
params_df <- tibble(
  param = c("bankroll", "kelly_mult"),
  value = c(bankroll, kelly_mult)
)
dbWriteTable(con, "nfl_params_temp", params_df, overwrite = TRUE)

# Save generation timestamp for freshness checking
dbExecute(con, "CREATE OR REPLACE TABLE nfl_samples_meta AS SELECT CURRENT_TIMESTAMP as generated_at")

dbDisconnect(con)

cat("=== NFL PREPARE: Complete. Samples saved to DuckDB. ===\n")
