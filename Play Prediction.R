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

nfl99all <- load_pbp(1999:2024)
nfl99 <- nfl99all %>%
  filter(pass == 1 | rush == 1 | qb_kneel == 1 | qb_spike == 1) %>% #kneels are not included
  mutate(explosive = ifelse((yards_gained>20 & pass_attempt == 1) | (yards_gained >12 & (qb_scramble == 1 | rush == 1)),1,0),
         negative = ifelse(yards_gained < 0, 1,0))

nfl_play_data <- nfl99 %>%
  filter(two_point_attempt!=1) %>%
  filter(play_type_nfl != "FIELD_GOAL", play_type_nfl != "PENALTY" ) %>% 
  group_by(game_id,posteam) %>% 
  summarize(offensive_plays = n(),total_pass_oe = sum(pass_oe,na.rm = T), proe = total_pass_oe/offensive_plays, PROE = mean(pass_oe)) %>% 
  filter(!is.na(posteam)) %>% 
  mutate(check = proe -  PROE)
  
all_schedules <- load_schedules(1999:2024)
check <- all_schedules %>% 
  # select(spread_line,wind,temp,surface,location,total_line,home_moneyline,away_moneyline)
  select(spread_line,wind,temp,surface,location,total_line)


x<- pbp_rp %>% 
  filter(is.na(pass_oe)) %>%
  filter(two_point_attempt!=1) %>% 
  select(desc,pass,pass_oe, two_point_attempt,epa)

help <- nfl99 %>% 
  filter(play_type_nfl == "PENALTY") %>% 
  select(desc,pass_oe,pass,play_type, play_type_nfl)
