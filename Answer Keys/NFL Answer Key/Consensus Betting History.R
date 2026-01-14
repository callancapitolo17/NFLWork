setwd("~/NFLWork/Answer Keys")
library(data.table)
library(duckdb)
library(dplyr)
library(lubridate)
library(DBI)
library(tidyverse)
source("Tools.R")

# 1) Load Data ----
con <- dbConnect(duckdb(), dbdir = "pbp.duckdb", read_only = TRUE)
pbp_all <- dbGetQuery(con, "SELECT * FROM nfl_pbp where season >=2020")
flat_odds <- dbGetQuery(con, "SELECT * FROM nfl_closing_odds")
dbDisconnect(con)

data("teams_colors_logos", package = "nflfastR")
team_lookup <- setNames(teams_colors_logos$team_name, teams_colors_logos$team_abbr)

# 2) Prep Box Scores ----
game_box_score <- pbp_all %>%
  group_by(game_id) %>%
  summarize(
    home_team = first(home_team),
    away_team = first(away_team),
    home_score = max(home_score, na.rm = TRUE),
    away_score = max(away_score, na.rm = TRUE),
    # Try time_of_day first, fall back to game_date if missing
    start_time = if_else(
      all(is.na(time_of_day)),
      NA_POSIXct_,
      ymd_hms(min(time_of_day, na.rm = TRUE), tz = "UTC")
    ),
    # Extract date from game_id as fallback (format: YYYY_WK_AWAY_HOME)
    game_date_from_id = as.Date(paste0(substr(game_id[1], 1, 4), "-01-01")) + 
      (as.numeric(substr(game_id[1], 6, 7)) - 1) * 7,
    .groups = "drop"
  ) %>%
  mutate(
    # Use game_date if available, otherwise approximate from game_id
    game_date = coalesce(as.Date(start_time), game_date_from_id),
    # For start_time, default to noon UTC if missing
    start_time = coalesce(start_time, as.POSIXct(paste(game_date, "12:00:00"), tz = "UTC")),
    total = home_score + away_score,
    home_winner = ifelse(home_score == away_score, NA, ifelse(home_score > away_score, 1, 0)),
    home_margin = home_score - away_score,
    home_team = unname(team_lookup[home_team]),
    away_team = unname(team_lookup[away_team])
  ) %>%
  select(-game_date_from_id)

# 3) Prep Odds ----
clean_odds <- flat_odds %>%
  mutate(
    devig_home_odds = devig_american(spread_home_odds, spread_away_odds)[[1]],
    devig_away_odds = devig_american(spread_home_odds, spread_away_odds)[[2]],
    devig_over_odds = devig_american(tot_over_odds, tot_under_odds)[[1]],
    devig_under_odds = devig_american(tot_over_odds, tot_under_odds)[[2]],
    game_date = as.Date(commence_time),
    commence_time = ymd_hms(commence_time, tz = "UTC")
  )

# 4) Join ----
eval_dt <- join_pbp_odds(
  pbp_dt = game_box_score %>% rename(game_start_time = start_time),
  odds_dt = clean_odds %>% rename(game_start_time = commence_time),
  join_cols = c("home_team", "away_team", "game_date", "game_start_time")
)

# 5) Spreads Consensus ----
spread_results <- build_market_consensus(
  eval_dt = eval_dt %>% mutate(spread_hit = ifelse(home_margin > -home_spread, 1, 0)),
  market_type = "spread",
  game_id_col = "game_id",
  outcome_col = "spread_hit",
  prob_col = "devig_home_odds",
  min_count_1yr = 0 #todo needs to be fixed
)

# 6) Totals Consensus ----
total_results <- build_market_consensus(
  eval_dt = eval_dt %>% mutate(over_hit = ifelse(total > total_line, 1, 0)),
  market_type = "totals",
  game_id_col = "game_id",
  outcome_col = "over_hit",
  prob_col = "devig_over_odds",
  min_count_1yr = 0
  
)

# 7) Build Final Consensus Lines ----
spread_consensus_df <- pick_consensus_line(
  df = eval_dt %>% inner_join(spread_results$weights, by = "bookmaker_key"),
  game_id_col = "game_id",
  line_col = "home_spread",
  weight_col = "spread_weight",
  date_col = "game_date",
  time_col = "game_start_time",
  market1 = "devig_home_odds",
  market2 = "devig_away_odds"
)

totals_consensus_df <- pick_consensus_line(
  df = eval_dt %>% inner_join(total_results$weights, by = "bookmaker_key"),
  game_id_col = "game_id",
  line_col = "total_line",
  weight_col = "totals_weight",
  date_col = "game_date",
  time_col = "game_start_time",
  market1 = "devig_over_odds",
  market2 = "devig_under_odds"
)

nfl_betting_history <- spread_consensus_df %>%
  inner_join(totals_consensus_df, by = c("game_id", "game_date", "game_start_time", "home_team", "away_team"))

# 8) Add PBP Periods & Save ----
game_by_period <- pbp_all %>%
  group_by(game_id, qtr) %>%
  summarize(
    home_team = first(home_team),
    away_team = first(away_team),
    # Try time_of_day first, fall back if missing
    period_start_time = if_else(
      all(is.na(time_of_day)),
      NA_POSIXct_,
      ymd_hms(min(time_of_day, na.rm = TRUE), tz = "UTC")
    ),
    # Extract date from game_id as fallback
    game_date_from_id = as.Date(paste0(substr(game_id[1], 1, 4), "-01-01")) + 
      (as.numeric(substr(game_id[1], 6, 7)) - 1) * 7,
    home_score_max = max(total_home_score, na.rm = TRUE),
    away_score_max = max(total_away_score, na.rm = TRUE),
    home_score_min = min(total_home_score, na.rm = TRUE),
    away_score_min = min(total_away_score, na.rm = TRUE),
    home_score = home_score_max - home_score_min,
    away_score = away_score_max - away_score_min,
    .groups = "drop"
  ) %>%
  mutate(
    game_home_margin_period = home_score - away_score,
    game_total_period = home_score + away_score
  ) %>%
  group_by(game_id) %>%
  mutate(
    # Use actual time if available, otherwise approximate from game_id
    game_date = coalesce(as.Date(min(period_start_time, na.rm = TRUE)), first(game_date_from_id)),
    game_start_time = coalesce(
      min(period_start_time, na.rm = TRUE),
      as.POSIXct(paste(first(game_date_from_id), "12:00:00"), tz = "UTC")
    ),
    home_final_score = sum(home_score, na.rm = TRUE),
    away_final_score = sum(away_score, na.rm = TRUE),
    total_final_score = home_final_score + away_final_score,
    home_margin = home_final_score - away_final_score,
    home_team = unname(team_lookup[home_team]),
    away_team = unname(team_lookup[away_team]),
    home_winner = ifelse(home_final_score == away_final_score, NA, ifelse(home_final_score > away_final_score, 1, 0))
  ) %>%
  ungroup() %>%
  select(-home_score, -away_score, -period_start_time, -away_score_max, -away_score_min,
         -home_score_max, -home_score_min, -game_date_from_id) %>%
  pivot_wider(names_from = qtr, values_from = c(game_home_margin_period, game_total_period))

game_by_half <- pbp_all %>%
  group_by(game_id, game_half) %>%
  summarize(
    # Try time_of_day first, fall back if missing
    period_start_time = if_else(
      all(is.na(time_of_day)),
      NA_POSIXct_,
      ymd_hms(min(time_of_day, na.rm = TRUE), tz = "UTC")
    ),
    # Extract date from game_id as fallback
    game_date_from_id = as.Date(paste0(substr(game_id[1], 1, 4), "-01-01")) + 
      (as.numeric(substr(game_id[1], 6, 7)) - 1) * 7,
    home_score_max = max(total_home_score, na.rm = TRUE),
    away_score_max = max(total_away_score, na.rm = TRUE),
    home_score_min = min(total_home_score, na.rm = TRUE),
    away_score_min = min(total_away_score, na.rm = TRUE),
    home_score = home_score_max - home_score_min,
    away_score = away_score_max - away_score_min,
    .groups = "drop"
  ) %>%
  mutate(
    game_home_margin_period = home_score - away_score,
    game_total_period = home_score + away_score
  ) %>% 
  ungroup() %>%
  select(-home_score, -away_score, -period_start_time, -away_score_max, -away_score_min,
         -home_score_max, -home_score_min, -game_date_from_id) %>%
  pivot_wider(names_from = game_half, values_from = c(game_home_margin_period, game_total_period))

game_by_half_period <- game_by_period %>% inner_join(game_by_half, by = "game_id")

joined_dt <- join_pbp_odds(game_by_half_period, nfl_betting_history)

nfl_weights <- spread_results$weights %>% inner_join(total_results$weights, by = "bookmaker_key")

con <- dbConnect(duckdb(), dbdir = "pbp.duckdb")
dbWriteTable(con, "nfl_betting_pbp", joined_dt, overwrite = TRUE)
dbWriteTable(con, "nfl_weights", nfl_weights, overwrite = TRUE)
dbDisconnect(con)
