# install + load
if (!requireNamespace("baseballr", quietly=TRUE)) install.packages("baseballr")
library(baseballr)
library(dplyr)
library(purrr)

# 1) Pick your season
season_year <- 2023

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