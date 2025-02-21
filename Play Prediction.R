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

test <- nfl99 %>% #no huddle, drive length, 3rd down/4th down
  filter(two_point_attempt!=1) %>%
  filter(season > 2005) %>% 
  filter(play_type_nfl != "FIELD_GOAL", play_type_nfl != "PENALTY" ) %>% 
  mutate(drive_id = paste(game_id,drive),
         drive_top = as.numeric(ms(drive_time_of_possession))) %>% 
  group_by(drive_id) %>% 
  summarize(plays_drive = n(),top = max(drive_top))
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



# Load necessary libraries
library(nflfastR)
library(dplyr)

# Load play-by-play data (example for 2023 season)
pbp <- load_pbp(2023)

# Compute time differences and infer clock stoppages
time_between_plays <- nfl99 %>%
  filter(qb_kneel != 1) %>% 
  arrange(game_id, posteam, desc(game_seconds_remaining)) %>%
  mutate(
    next_play_time = lead(game_seconds_remaining),  # Game clock at next play
    time_between_plays = game_seconds_remaining - next_play_time,  # Time diff
    
    # Clock stopped if the game clock doesn't move OR an excessive gap exists
    clock_stopped = ifelse(
      is.na(next_play_time) | next_play_time == game_seconds_remaining | time_between_plays > 40,
      1, 0
    )
  ) %>% 
  filter(clock_stopped == 0, time_between_plays >0) %>%
  select(game_seconds_remaining,next_play_time,time_between_plays,desc)
  group_by(game_id, posteam) %>%
  summarise(
    avg_time_between_plays = mean(time_between_plays, na.rm = TRUE),
    median_time_between_plays = median(time_between_plays, na.rm = TRUE),
    n_plays = n()
  ) %>%
  ungroup()

# View results
head(pbp_filtered)
