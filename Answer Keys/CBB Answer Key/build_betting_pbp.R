setwd("~/NFLWork/Answer Keys")
library(duckdb)
library(DBI)
library(dplyr)

# Connect to DB
con <- dbConnect(duckdb(), dbdir = "cbb.duckdb")

# Load data
message("Loading data from database...")

# Load team name mapping (OddsAPI -> ESPN)
team_mapping <- dbGetQuery(con, "SELECT * FROM team_name_mapping")
message(sprintf("Team mapping: %d entries", nrow(team_mapping)))

# Get odds data and map to ESPN names
odds_data <- dbGetQuery(con, "
  SELECT DISTINCT
    id as odds_id,
    DATE(commence_time) as game_date,
    home_team as odds_home,
    away_team as odds_away
  FROM cbb_closing_odds
") %>%
  left_join(team_mapping, by = c("odds_home" = "odds_name")) %>%
  rename(home_espn = espn_name) %>%
  left_join(team_mapping, by = c("odds_away" = "odds_name")) %>%
  rename(away_espn = espn_name)

# Check which race-to-X columns exist (may not be backfilled yet)
race_cols <- c("first_to_10_h1", "first_to_10_fg", "first_to_20_fg", "first_to_40_fg")
available_race_cols <- character()
for (col in race_cols) {
  exists <- tryCatch({
    dbGetQuery(con, sprintf("SELECT %s FROM cbb_pbp_v2 LIMIT 0", col))
    TRUE
  }, error = function(e) FALSE)
  if (exists) available_race_cols <- c(available_race_cols, col)
}
race_col_sql <- if (length(available_race_cols) > 0) paste0(", ", paste(available_race_cols, collapse = ", ")) else ""

# Get PBP data (CBBpy uses ESPN names)
pbp_data <- dbGetQuery(con, sprintf("
  SELECT
    game_id,
    DATE(game_date) as game_date,
    home_team,
    away_team,
    -- Half 1
    home_h1_score,
    away_h1_score,
    game_home_margin_h1,
    game_total_h1,
    -- Half 2 (regulation only)
    home_h2_score,
    away_h2_score,
    game_home_margin_h2,
    game_total_h2,
    -- Overtime
    home_ot_score,
    away_ot_score,
    game_home_margin_ot,
    game_total_ot,
    -- Full game
    home_final_score,
    away_final_score,
    game_home_margin_fg,
    game_total_fg,
    home_winner,
    went_to_ot
    %s
  FROM cbb_pbp_v2
  WHERE game_home_margin_h2 IS NOT NULL
", race_col_sql))

# Ensure all race-to-X columns exist in dataframe (NA if not backfilled)
for (col in c("first_to_10_h1", "first_to_10_fg", "first_to_20_fg", "first_to_40_fg")) {
  if (!(col %in% names(pbp_data))) {
    pbp_data[[col]] <- NA_integer_
  }
}

message(sprintf("Odds: %d unique games", nrow(odds_data)))
message(sprintf("PBP: %d complete games", nrow(pbp_data)))
message(sprintf("Odds with ESPN home mapping: %d (%.1f%%)",
                sum(!is.na(odds_data$home_espn)),
                100 * sum(!is.na(odds_data$home_espn)) / nrow(odds_data)))

# Join on date + ESPN home team name
# OddsAPI uses UTC, so late-night Pacific games appear on the next day
# Try matching on both original date and date-1
message("\nJoining on date + ESPN team name (with timezone adjustment)...")

# First try original date match
joined_same_day <- odds_data %>%
  filter(!is.na(home_espn)) %>%
  inner_join(
    pbp_data,
    by = c("game_date" = "game_date", "home_espn" = "home_team"),
    suffix = c("_odds", "_pbp")
  ) %>%
  mutate(match_type = "same_day")

# For unmatched, try date-1 (timezone offset)
unmatched_odds <- odds_data %>%
  filter(!is.na(home_espn)) %>%
  anti_join(joined_same_day, by = "odds_id")

joined_day_before <- unmatched_odds %>%
  mutate(game_date_adj = game_date - 1) %>%
  inner_join(
    pbp_data,
    by = c("game_date_adj" = "game_date", "home_espn" = "home_team"),
    suffix = c("_odds", "_pbp")
  ) %>%
  select(-game_date_adj) %>%
  mutate(match_type = "day_before")

# Combine both matches
joined <- bind_rows(joined_same_day, joined_day_before)

message(sprintf("  Same-day matches: %d", nrow(joined_same_day)))
message(sprintf("  Day-before matches: %d", nrow(joined_day_before)))

message(sprintf("Matched: %d games", nrow(joined)))
message(sprintf("Match rate: %.1f%%", 100 * nrow(joined) / nrow(odds_data)))

# Check match rate by season
message("\nMatch rate by season:")
for (yr in 2021:2026) {
  start_date <- as.Date(paste0(yr - 1, "-11-01"))
  end_date <- as.Date(paste0(yr, "-04-15"))
  season_odds <- odds_data %>% filter(game_date >= start_date & game_date <= end_date)
  season_joined <- joined %>% filter(game_date >= start_date & game_date <= end_date)
  if (nrow(season_odds) > 0) {
    message(sprintf("  %d-%d: %d/%d matched (%.1f%%)",
                    yr - 1, yr, nrow(season_joined), nrow(season_odds),
                    100 * nrow(season_joined) / nrow(season_odds)))
  }
}

# Create betting_pbp table (scores now come directly from PBP, including OT tracking)
betting_pbp <- joined %>%
  select(
    odds_id,
    game_id,
    game_date,
    home_team = odds_home,
    away_team = odds_away,
    # Half 1
    home_h1_score,
    away_h1_score,
    game_home_margin_h1,
    game_total_h1,
    # Half 2 (regulation only)
    home_h2_score,
    away_h2_score,
    game_home_margin_h2,
    game_total_h2,
    # Overtime
    home_ot_score,
    away_ot_score,
    game_home_margin_ot,
    game_total_ot,
    # Full game
    home_final_score,
    away_final_score,
    game_home_margin_fg,
    game_total_fg,
    home_winner,
    went_to_ot,
    first_to_10_h1,
    first_to_10_fg,
    first_to_20_fg,
    first_to_40_fg
  ) %>%
  distinct()

message(sprintf("\nFinal betting_pbp: %d rows", nrow(betting_pbp)))

# Save to database
dbWriteTable(con, "cbb_betting_pbp", betting_pbp, overwrite = TRUE)
message("Saved to cbb_betting_pbp table")

# Summary by season
summary_stats <- betting_pbp %>%
  mutate(season = ifelse(as.numeric(format(game_date, "%m")) >= 11,
                         as.numeric(format(game_date, "%Y")) + 1,
                         as.numeric(format(game_date, "%Y")))) %>%
  group_by(season) %>%
  summarize(games = n(), .groups = "drop")
message("\nGames by season in cbb_betting_pbp:")
print(summary_stats)

dbDisconnect(con, shutdown = TRUE)
message("Done!")
