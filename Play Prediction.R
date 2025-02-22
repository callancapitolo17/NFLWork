library(nflfastR)
library(ggplot2)
library(tidyverse)
library(ggimage)
library(ggthemes)
library(dplyr)
library(ggrepel)
library(nflreadr)
library(gt)
library(ggrepel)
library(nflplotR)
library(gtExtras)
library(lubridate)

nfl99all <- load_pbp(1999:2024)
nfl99 <- nfl99all %>%
  filter(pass == 1 | rush == 1 | qb_kneel == 1 | qb_spike == 1) %>% #kneels are not included
  mutate(explosive = ifelse((yards_gained>20 & pass_attempt == 1) | (yards_gained >12 & (qb_scramble == 1 | rush == 1)),1,0),
         negative = ifelse(yards_gained < 0, 1,0))

nfl_play_data <- nfl99 %>% #no huddle, drive length, 3rd down/4th down
  filter(two_point_attempt!=1) %>%
  filter(season > 2005) %>% 
  filter(play_type_nfl != "FIELD_GOAL", play_type_nfl != "PENALTY" ) %>% 
  group_by(game_id,posteam,week,season) %>% 
  summarize(offensive_plays = n(),total_pass_oe = sum(pass_oe,na.rm = T), proe_eligible_plays = sum(!is.na(pass_oe),na.rm = T), proe = total_pass_oe/proe_eligible_plays) %>% 
  filter(!is.na(posteam)) %>% 
  arrange(game_id,posteam,week) %>% 
  ungroup() %>% 
  group_by(posteam,season) %>% 
  mutate(plays_per_game = cummean(offensive_plays))

#drive stats
drive <- nfl99 %>% #no huddle, drive length, 3rd down/4th down
  mutate(
    # Approximate play end time using next play's start time
    game_seconds_end = lead(game_seconds_remaining),
    
    # Handle the last play of the game (set to 0 if it's missing)
    game_seconds_end = ifelse(is.na(game_seconds_end), 0, game_seconds_end),
    
    # Time elapsed on the play (play start - next play start)
    play_duration = game_seconds_remaining - game_seconds_end
  ) %>%
  filter(two_point_attempt!=1) %>%
  filter(season > 2005) %>% 
  filter(play_type_nfl != "FIELD_GOAL", play_type_nfl != "PENALTY" ) %>%
  arrange(game_id, play_id) %>%
  mutate(
    # Detect when a new drive should start
    new_drive_flag = (posteam != lag(posteam, default = first(posteam))) |  # Change in possession
      (lag(qtr, default = first(qtr)) == 2 & qtr == 3) |      # Halftime reset
      (lag(play_type, default = "NA") %in% c("punt", "field_goal", "interception", "fumble")),  # Turnover plays
    
    # Assign new drive numbers
    corrected_drive = cumsum(coalesce(new_drive_flag, 0)) + 1
  ) %>% 
  mutate(drive_id = paste(game_id,corrected_drive),
         drive_top = as.numeric(ms(drive_time_of_possession))) %>% 
  filter(game_id %in% c("2019_15_CLE_ARI")) %>% 
  select(desc,drive_id, drive_play_count,drive_time_of_possession,fixed_drive,drive,drive_game_clock_start,drive_game_clock_end,corrected_drive,game_seconds_remaining,time)
  group_by(game_id,drive_id,posteam,week) %>% 
  summarize(plays_drive = n(),top = max(drive_top), time_per_play = top/plays_drive, max_time = max(game_seconds_remaining),
            min_time = min(game_seconds_remaining), drive_time = max_time - min_time, time_per_play = drive_time/plays_drive)
  group_by(game_id,posteam,week,season) %>% 
  summarize(offensive_plays = n(),total_pass_oe = sum(pass_oe,na.rm = T), proe_eligible_plays = sum(!is.na(pass_oe),na.rm = T), proe = total_pass_oe/proe_eligible_plays) %>% 
  filter(!is.na(posteam)) %>% 
  arrange(game_id,posteam,week) %>% 
  ungroup() %>% 
  group_by(posteam,season) %>% 
  mutate(plays_per_game = cummean(offensive_plays))
  
all_schedules <- load_schedules(1999:2024)

american_to_decimal_odds <- function(odds){ifelse(odds < 0, 1 - (100/odds), 1+(odds/100))}
no_vig_odds <- function(odds_1,odds_2) {
  # Ensure odds are numeric
  prob_odds_1 <- 1/american_to_decimal_odds(odds_1)
  prob_odds_2 <- 1/american_to_decimal_odds(odds_2)
  
  # Remove the vig (normalize so that the sum is 1)
  no_vig_odds_1 <- prob_odds_1 / (prob_odds_1+prob_odds_2)
  no_vig_odds_2 <- prob_odds_2 / (prob_odds_1+prob_odds_2)
  
  return(c(no_vig_odds_1,no_vig_odds_2))
}

#Explore model types with different variables,spread, moneyline, ml+spread
clean_schedule <- all_schedules %>% 
  select(game_id,spread_line,wind,temp,surface,roof, location,total_line,home_moneyline,away_moneyline,home_team,away_team) %>% 
  mutate(
    home_prob = map2_dbl(home_moneyline, away_moneyline, ~ no_vig_odds(.x, .y)[1]),
    away_prob = map2_dbl(home_moneyline, away_moneyline, ~ no_vig_odds(.x, .y)[2])
  ) %>% 
  select(-home_moneyline,-away_moneyline) %>% 
  pivot_longer(cols = c(home_team,away_team), names_to = "home_or_away", values_to = "team") %>% 
  mutate(clean_spread = ifelse(home_or_away == "home_team" , spread_line*-1,spread_line),
         win_prob = ifelse(home_or_away == "home_team", home_prob,away_prob)) %>% 
  left_join(all_schedules %>% 
              select(game_id,home_team,away_team), by = c("game_id")) %>% 
  mutate(opponent = ifelse(home_team == team, away_team,home_team)) %>% 
  select(-away_team,-home_team)
joined_play_data <- nfl_play_data %>% 
  left_join(clean_schedule, by = c("game_id", "posteam" = "team"))
