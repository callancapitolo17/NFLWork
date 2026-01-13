library(baseballr)
library(dplyr)
library(purrr)
library(tibble)
library(progress)
library(tidyr)
library(lubridate)

# Function to safely pull PBP for one season
get_season_pbp <- function(season_year) {
  message("Pulling schedule for ", season_year, "...")
  
  sched <- baseballr::mlb_schedule(season = season_year)
  
  game_pks <- sched %>%
    filter(
      !series_description %in% c("Exhibition", "Spring Training"),
      status_coded_game_state == "F"
    ) %>%
    pull(game_pk)
  
  pb <- progress_bar$new(
    format = paste0("  ", season_year, " [:bar] :percent eta: :eta"),
    total = length(game_pks), clear = FALSE, width = 60
  )
  
  safely_get <- possibly(~ {
    pb$tick()
    get_pbp_mlb(.x)
  }, otherwise = NULL)
  
  season_pbp <- map_dfr(game_pks, safely_get, .id = "game_pk")
  
  return(season_pbp)
}

# Loop over all seasons since 2007
# seasons <- 2015:2024
seasons <- 2025

for (season in seasons) {
  pbp <- get_season_pbp(season)
  saveRDS(pbp, file = paste0("pbp_", season, ".rds"))
  message("Saved season: ", season)
}
df_pbp_all <- pbp_all %>% select(-reviewDetails.additionalReviews)

seasons <- 2015:2024

# convert to .parquet once
library(arrow)
library(purrr)
library(dplyr)

# Helper: Drop all list-columns
drop_list_cols <- function(df) {
  df %>%
    select(where(~ !is.list(.)))
}

# Convert and save
for (season in seasons) {
  df <- readRDS(paste0("pbp_", season, ".rds"))
  df_clean <- drop_list_cols(df)
  write_parquet(df_clean, paste0("pbp_", season, ".parquet"))
}

parquet_files <- list.files(pattern = "pbp_\\d{4}\\.parquet$")

# Optional: extract the season year from filename
get_season <- function(file) {
  gsub("pbp_(\\d{4})\\.parquet", "\\1", file)
}

# Read and combine
pbp_all <- map_dfr(parquet_files[6:length(parquet_files)], function(file) { #get just pbp for games with betting data
  df <- read_parquet(file)
  df$season <- get_season(file)
  df
})

library(DBI)
library(duckdb)
library(glue)


# files <- parquet_files
# files_sql <- paste0("['", paste(files, collapse = "','"), "']")
# 
# # Create a table by scanning the parquet files. filename=true adds a filename column.
# qry <- glue("
#   CREATE OR REPLACE TABLE pbp_all AS
#   SELECT
#     t.*,
#     CAST(regexp_extract(filename, '([0-9]{{4}})', 1) AS INTEGER) AS season
#   FROM read_parquet({files_sql}, filename=true) AS t
# ")
# dbExecute(con, qry)
con <- dbConnect(duckdb(), dbdir = "pbp.duckdb")


pbp_all <- dbGetQuery(con, "
  SELECT game_pk,season,\"about.inning\", game_date, home_team, away_team,\"about.startTime\",\"result.homeScore\",\"result.awayScore\",
  FROM pbp_all
  WHERE SEASON >= 2020
")
DBI::dbDisconnect(con, shutdown = TRUE)

#Analysis ----
game_by_inning <- pbp_all %>% 
  # filter(!is.na(about.inning)) %>%
  group_by(game_pk,season, about.inning) %>% 
  summarize(game_date = first(game_date),home_team = first(home_team), away_team = first(away_team),inning_start_time = min(about.startTime), home_score = max(result.homeScore),away_score = max(result.awayScore)) %>% 
  mutate(game_home_margin_in = home_score-away_score) %>% 
  mutate(game_total_in = home_score+away_score) %>% 
  ungroup() %>% 
  group_by(season,game_pk) %>% 
  mutate(game_start_time = min(inning_start_time)) %>% 
  mutate(game_start_time = ymd_hms(game_start_time, tz = "UTC")) %>% 
  mutate(home_final_score = max(home_score), away_final_score = max(away_score), 
         total_final_score = home_final_score + away_final_score,
         home_winner = ifelse(home_final_score == away_final_score, NA,ifelse(home_final_score > away_final_score, 1,0))) %>% 
  ungroup() %>% 
  select(-home_score,-away_score,-inning_start_time) %>% #allow a proper pivt - score isn't need because of game score
  pivot_wider(names_from = about.inning, values_from = c(game_home_margin,full_game_total_inning)) %>% 
  mutate(game_date = as.Date(game_start_time)) #ensures dates are in the same format
