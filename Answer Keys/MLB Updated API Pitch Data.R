# install + load
if (!requireNamespace("baseballr", quietly=TRUE)) install.packages("baseballr")
library(baseballr)
library(dplyr)
library(purrr)

# 1) Pick your season
season_year <- 2007:2025

# 2) Fetch the full MLB schedule for that season
sched <- baseballr::mlb_schedule(season = season_year)

# 3) Extract all game_pk values
game_pks <- pull(sched %>% 
  filter(!series_description %in% c("Exhibition","Spring Training") & status_coded_game_state == "F") %>% 
    select(game_pk))

# 4) Pull play‐by‐play for each game_pk and row‐bind
#    (this will take a few minutes for ~2,430 games)
pbp_season <- map_dfr(
  game_pks,
  ~ get_pbp_mlb(game_pk = .x),
  .id = "game_pk"
)

# 5) Quick sanity check
glimpse(pbp_season)

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

pbp_all <- map_dfr(seasons, get_season_pbp, .id = "season")
df_pbp_all <- pbp_all %>% select(-reviewDetails.additionalReviews)
write.csv(df_pbp_all,"Pitch_Level_Data_15-24.csv")

first6 <- head(df_pbp_all)

games_by_inning <- df_pbp_all %>% 
  group_by(game_pk,season,about.inning) %>% 
  summarize(first(game_date),first(home_team), first(away_team),home_score = max(result.homeScore),away_score = max(result.awayScore)) %>% 
  mutate(home_winner_using_game_score = ifelse(home_score > away_score,1,0)) %>% 
  pivot_wider(names_from = about.inning, values_from = home_winner_using_game_score)
