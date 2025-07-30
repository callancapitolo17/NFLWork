

library(baseballr)
library(dplyr)
library(purrr)
library(tibble)
library(progress)
library(tidyr)

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
seasons <- 2015:2024

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
pbp_all <- map_dfr(parquet_files, function(file) {
  df <- read_parquet(file)
  df$season <- get_season(file)
  df
})

#Analysis ----
game_by_inning <- pbp_all %>% 
  # filter(!is.na(about.inning)) %>% 
  group_by(game_pk,season, about.inning) %>% 
  summarize(game_date = first(game_date),home_team = first(home_team), away_team = first(away_team),inning_start_time = min(about.startTime), home_score = max(result.homeScore),away_score = max(result.awayScore)) %>% 
  mutate(game_home_ml_inning = ifelse(home_score > away_score,1,ifelse(home_score== away_score,NA,0))) %>% 
  mutate(full_game_total_inning = home_score+away_score) %>% 
  ungroup() %>% 
  group_by(season,game_pk) %>% 
  mutate(game_start_time = min(inning_start_time)) %>% 
  ungroup() %>% 
  select(-home_score,-away_score,-inning_start_time) %>% #allow a proper pivt - score isn't need because of game score
  pivot_wider(names_from = about.inning, values_from = c(game_home_ml_inning,full_game_total_inning))


game_total_by_inning <- pbp_all %>% 
  filter(!is.na(about.inning)) %>% 
  group_by(game_pk,season,about.inning) %>% 
  summarize(first(game_date),first(home_team), first(away_team),home_score = max(result.homeScore),away_score = max(result.awayScore),
            start_time = min(about.startTime)) %>% 
  mutate(total = home_score+away_score) %>% 
  ungroup() %>% 
  select(-home_score,-away_score) %>% #allow a proper pivt - score isn't need because of game score
  pivot_wider(names_from = about.inning, values_from = total,names_prefix = "Total")
