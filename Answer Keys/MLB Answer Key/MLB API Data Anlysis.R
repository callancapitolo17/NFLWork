if (!requireNamespace("baseballr", quietly=TRUE)) install.packages("baseballr")
library(baseballr)
library(dplyr)
library(purrr)
library(tidyr)
df_pbp_all <- read.csv("Pitch_Level_Data_15-24.csv")

first6 <- head(df_pbp_all)

game_ml_by_inning <- df_pbp_all %>% 
  # filter(!is.na(about.inning)) %>% 
  group_by(game_pk,season,about.inning) %>% 
  summarize(first(game_date),first(home_team), first(away_team),home_score = max(result.homeScore),away_score = max(result.awayScore)) %>% 
  mutate(home_winner_using_game_score = ifelse(home_score > away_score,1,ifelse(home_score== away_score,NA,0))) %>% 
  ungroup() %>% 
  select(-home_score,-away_score) %>% #allow a proper pivt - score isn't need because of game score
  pivot_wider(names_from = about.inning, values_from = home_winner_using_game_score, names_prefix = "HomeML")


game_total_by_inning <- df_pbp_all %>% 
  filter(!is.na(about.inning)) %>% 
  group_by(game_pk,season,about.inning) %>% 
  summarize(first(game_date),first(home_team), first(away_team),home_score = max(result.homeScore),away_score = max(result.awayScore)) %>% 
  mutate(total = home_score+away_score) %>% 
  ungroup() %>% 
  select(-home_score,-away_score) %>% #allow a proper pivt - score isn't need because of game score
  pivot_wider(names_from = about.inning, values_from = total,names_prefix = "Total")
