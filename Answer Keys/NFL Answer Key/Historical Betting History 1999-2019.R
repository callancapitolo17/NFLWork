setwd("~/NFLWork/Answer Keys")
library(data.table)
library(duckdb)
library(dplyr)
library(tidyr)
library(lubridate)
library(DBI)
library(nflfastR)
source("Tools.R")

# =============================================================================
# Historical NFL Closing Odds (1999-2019)
# Uploads historical betting lines to nfl_closing_odds table
# Schema matches: id, commence_time, home_team, away_team, bookmaker_key,
#                 bookmaker_title, bookmaker_update, ml_home_odds, ml_away_odds,
#                 home_spread, away_spread, spread_home_odds, spread_away_odds,
#                 total_line, tot_over_odds, tot_under_odds
# =============================================================================

# 1) Load Historical Schedules ----
cat("Loading historical schedules (1999-2019)...\n")
historical_schedules <- fast_scraper_schedules(1999:2019) %>%
  filter(!is.na(home_score), !is.na(away_score)) %>%  # Completed games only
  filter(!is.na(spread_line), !is.na(total_line)) %>%  # Must have betting lines
  mutate(
    # Handle missing odds - default to -110
    home_spread_odds = ifelse(is.na(home_spread_odds), -110L, as.integer(home_spread_odds)),
    away_spread_odds = ifelse(is.na(away_spread_odds), -110L, as.integer(away_spread_odds)),
    over_odds = ifelse(is.na(over_odds), -110L, as.integer(over_odds)),
    under_odds = ifelse(is.na(under_odds), -110L, as.integer(under_odds)),
    # Moneylines - use from schedule if available, otherwise NA
    home_moneyline = as.integer(home_moneyline),
    away_moneyline = as.integer(away_moneyline)
  )

cat(sprintf("Loaded %d historical games\n", nrow(historical_schedules)))

# 2) Format to match nfl_closing_odds schema ----
cat("Formatting to match nfl_closing_odds schema...\n")

# Map team abbreviations to full names (matching odds API format)
data("teams_colors_logos", package = "nflfastR")
team_lookup <- setNames(teams_colors_logos$team_name, teams_colors_logos$team_abbr)

# Handle historical team relocations
historical_lookup <- c(
  team_lookup,
  "SD" = "Los Angeles Chargers",
  "STL" = "Los Angeles Rams",
  "OAK" = "Las Vegas Raiders",
  "LAR" = "Los Angeles Rams"
)

historical_closing_odds <- historical_schedules %>%
  mutate(
    # Match nfl_closing_odds schema exactly
    id = game_id,  # Use game_id as the event id
    commence_time = as.POSIXct(paste(gameday, ifelse(is.na(gametime), "13:00", gametime)), 
                                tz = "America/New_York") %>% with_tz("UTC"),
    home_team = unname(historical_lookup[home_team]),
    away_team = unname(historical_lookup[away_team]),
    bookmaker_key = "consensus_espn",  # Identify as historical ESPN consensus
    bookmaker_title = "ESPN Consensus (Historical)",
    bookmaker_update = commence_time,  # Use game time as update time
    
    # Moneylines
    ml_home_odds = home_moneyline,
    ml_away_odds = away_moneyline,
    
    # Spreads - nflfastR spread_line is from home perspective (negative = home favored)
    home_spread = -spread_line,  # Convert: spread_line is "away spread", so negate for home
    away_spread = spread_line,
    spread_home_odds = home_spread_odds,
    spread_away_odds = away_spread_odds,
    
    # Totals
    total_line = total_line,
    tot_over_odds = over_odds,
    tot_under_odds = under_odds
  ) %>%
  select(
    id, commence_time, home_team, away_team,
    bookmaker_key, bookmaker_title, bookmaker_update,
    ml_home_odds, ml_away_odds,
    home_spread, away_spread, spread_home_odds, spread_away_odds,
    total_line, tot_over_odds, tot_under_odds
  )

# 3) Summary Stats ----
cat("\n=== SUMMARY ===\n")
cat(sprintf("Total records: %d\n", nrow(historical_closing_odds)))
cat(sprintf("Years: %d - %d\n", 
            year(min(historical_closing_odds$commence_time)),
            year(max(historical_closing_odds$commence_time))))
cat(sprintf("Date range: %s to %s\n",
            min(as.Date(historical_closing_odds$commence_time)),
            max(as.Date(historical_closing_odds$commence_time))))

cat("\nSpread distribution:\n")
print(summary(historical_closing_odds$home_spread))

cat("\nTotal line distribution:\n")
print(summary(historical_closing_odds$total_line))

cat("\nMoneyline availability:\n")
cat(sprintf("  Home ML available: %d (%.1f%%)\n", 
            sum(!is.na(historical_closing_odds$ml_home_odds)),
            100 * mean(!is.na(historical_closing_odds$ml_home_odds))))

# 4) Save to Database ----
cat("\n=== SAVING TO DATABASE ===\n")
cat("NOTE: Make sure to close any RStudio connections to pbp.duckdb first!\n")
cat("Run: DBI::dbDisconnect(con) in your RStudio session\n\n")

# Uncomment these lines when ready to save:
# con <- dbConnect(duckdb(), dbdir = "pbp.duckdb")
# dbWriteTable(con, "nfl_odds_pre_20", historical_closing_odds, append = TRUE)
# dbDisconnect(con, shutdown = TRUE)
# cat("Appended to nfl_closing_odds!\n")
con <- dbConnect(duckdb(), dbdir = "pbp.duckdb")
pbp_all <- dbGetQuery(con, "SELECT * FROM nfl_pbp where season < 2020")
flat_odds <- dbGetQuery(con, "SELECT * FROM nfl_odds_pre_20")
dbDisconnect(con, shutdown = TRUE)
data("teams_colors_logos", package = "nflfastR")
team_lookup <- setNames(teams_colors_logos$team_name, teams_colors_logos$team_abbr)
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

# Transform flat_odds to match nfl_betting_history format from Consensus Betting History
# Since historical data has only one "bookmaker" (ESPN consensus), create consensus columns directly
nfl_betting_history_pre20 <- flat_odds %>%
  mutate(
    game_id = id,
    game_date = as.Date(commence_time),
    game_start_time = ymd_hms(commence_time, tz = "UTC"),
    # Devig the spread odds
    devig_home_odds = devig_american(spread_home_odds, spread_away_odds)[[1]],
    devig_away_odds = devig_american(spread_home_odds, spread_away_odds)[[2]],
    # Devig the total odds  
    devig_over_odds = devig_american(tot_over_odds, tot_under_odds)[[1]],
    devig_under_odds = devig_american(tot_over_odds, tot_under_odds)[[2]],
    # Create consensus columns (same as devig since single source)
    consensus_devig_home_odds = devig_home_odds,
    consensus_devig_away_odds = devig_away_odds,
    consensus_devig_over_odds = devig_over_odds,
    consensus_devig_under_odds = devig_under_odds,
    # Rename spread to match expected column name
    spread = home_spread
  ) %>%
  select(
    game_id, game_date, game_start_time, home_team, away_team,
    spread, total_line,
    consensus_devig_home_odds, consensus_devig_away_odds,
    consensus_devig_over_odds, consensus_devig_under_odds
  )

# Join with PBP data using the same function as Consensus Betting History
pre_20_betting_history <- join_pbp_odds(game_by_half_period, nfl_betting_history_pre20)

con <- dbConnect(duckdb(), dbdir = "pbp.duckdb")
dbWriteTable(con, "nfl_pre_20_betting_history", pre_20_betting_history, overwrite = TRUE)
dbDisconnect(con)
