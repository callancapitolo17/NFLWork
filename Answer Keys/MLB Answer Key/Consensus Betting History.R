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

#think about timing of betting lines

# 1) Load Data ----
con <- dbConnect(duckdb(), dbdir = "pbp.duckdb", read_only = TRUE)
on.exit(tryCatch(dbDisconnect(con), error = function(e) NULL), add = TRUE)
pbp_all <- dbGetQuery(con, "
  SELECT game_pk, year(CAST(game_date AS DATE)) AS season, \"about.inning\",
         game_date, home_team, away_team, \"about.startTime\",
         \"result.homeScore\", \"result.awayScore\"
  FROM mlb_pbp_all
  WHERE year(CAST(game_date AS DATE)) >= 2020
")
# mlb_betting_history contains closing odds (T-15min snapshot, captured by
# get_event_odds_by_id in Acquire New MLB Data.R)
flat_odds <- dbGetQuery(con, "SELECT * FROM mlb_betting_history WHERE commence_time > '2020-01-01'")
dbDisconnect(con)
on.exit(NULL)

# 2) Prep Box Scores ----
game_box_score <- pbp_all %>%
  group_by(game_pk,game_date) %>%
  summarize(
    game_date = first(game_date),
    home_team = first(home_team),
    start_time = ymd_hms(min(`about.startTime`, na.rm = TRUE), tz = "UTC"),
    away_team = first(away_team),
    home_score = max(`result.homeScore`),
    away_score = max(`result.awayScore`)
  ) %>%
  mutate(
    total = home_score + away_score,
    home_winner = ifelse(home_score == away_score, NA, ifelse(home_score > away_score, 1, 0)),
    home_margin = home_score - away_score,
    game_date = as.Date(game_date)
  )

# 3) Prep Odds ----
clean_odds <- flat_odds %>%
  mutate(
    devig_home_odds = devig_american(ml_home_odds, ml_away_odds)[[1]],
    devig_away_odds = devig_american(ml_home_odds, ml_away_odds)[[2]],
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

# 5) Moneyline Consensus ----
ml_results <- build_market_consensus(
  eval_dt = eval_dt,
  market_type = "moneyline",
  game_id_col = "game_pk",
  outcome_col = "home_winner",
  prob_col = "devig_home_odds"
)

# 6) Totals Consensus ----
total_results <- build_market_consensus(
  eval_dt = eval_dt %>% mutate(over_hit = ifelse(total > total_line, 1, 0)),
  market_type = "totals",
  game_id_col = "game_pk",
  outcome_col = "over_hit",
  prob_col = "devig_over_odds"
)

# 7) Build Final Consensus Lines ----
ml_consensus_df <- moneyline_consensus(
  df = eval_dt %>% inner_join(ml_results$weights, by = "bookmaker_key"),
  game_id_col = "game_pk",
  weight_col = "moneyline_weight",
  date_col = "game_date",
  time_col = "game_start_time",
  market1 = "devig_home_odds",
  market2 = "devig_away_odds"
)

totals_consensus_df <- pick_consensus_line(
  df = eval_dt %>% inner_join(total_results$weights, by = "bookmaker_key"),
  game_id_col = "game_pk",
  line_col = "total_line",
  weight_col = "totals_weight",
  date_col = "game_date",
  time_col = "game_start_time",
  market1 = "devig_over_odds",
  market2 = "devig_under_odds"
)

mlb_betting_history <- ml_consensus_df %>%
  inner_join(totals_consensus_df, by = c("game_pk", "game_date", "game_start_time", "home_team", "away_team"))

# 8) Add PBP Periods & Save ----
game_by_inning <- pbp_all %>%
  group_by(game_pk, game_date,`about.inning`) %>%
  summarize(
    game_date = first(game_date),
    home_team = first(home_team),
    away_team = first(away_team),
    inning_start_time = min(`about.startTime`, na.rm = TRUE),
    home_score = max(`result.homeScore`, na.rm = TRUE),
    away_score = max(`result.awayScore`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    game_home_margin_inning = home_score - away_score,
    game_total_inning = home_score + away_score
  ) %>%
  group_by(game_pk,game_date) %>%
  mutate(game_start_time = min(inning_start_time)) %>%
  mutate(game_start_time = ymd_hms(game_start_time, tz = "UTC")) %>%
  mutate(
    home_final_score = max(home_score),
    away_final_score = max(away_score),
    total_final_score = home_final_score + away_final_score,
    home_margin = home_final_score - away_final_score,
    home_winner = ifelse(home_final_score == away_final_score, NA, ifelse(home_final_score > away_final_score, 1, 0))
  ) %>%
  ungroup() %>%
  select(-home_score, -away_score, -inning_start_time) %>%
  pivot_wider(
    names_from = `about.inning`, 
    values_from = c(game_home_margin_inning, game_total_inning),
    names_prefix = "inning_"
  )

joined_dt <- join_pbp_odds(game_by_inning, mlb_betting_history, 
                           join_cols = c("game_pk", "game_date", "home_team", "away_team", "game_start_time"))

mlb_weights <- ml_results$weights %>% inner_join(total_results$weights, by = "bookmaker_key")

con <- dbConnect(duckdb(), dbdir = "pbp.duckdb")
on.exit(tryCatch(dbDisconnect(con), error = function(e) NULL), add = TRUE)
dbWriteTable(con, "mlb_betting_pbp", joined_dt, overwrite = TRUE)
dbWriteTable(con, "mlb_weights", mlb_weights, overwrite = TRUE)
message(sprintf("Wrote %d rows to mlb_betting_pbp.", nrow(joined_dt)))
dbDisconnect(con)
on.exit(NULL)
