#Part 1 ----
library(dplyr)
library(stringr)
race_info <- read.csv("race_info_2024.csv") #overall race info
race_sessions <- read.csv("race_sessions_2024.csv") #information about each session in race
session_competitors <- read.csv("session_competitors_2024.csv") #competitor results
start_list <- read.csv("start_list.csv")

# keeps the first occurrence of each ID, with all other columns
race_info_unique <- race_info %>% 
  distinct(ID, .keep_all = TRUE)

race_sessions_unique <- race_sessions %>% 
  distinct(ID, .keep_all = TRUE)

race_competitors_unique <- session_competitors %>% 
  distinct(ID, .keep_all = TRUE)

clean_racing_df <- race_competitors_unique %>% 
  rename("race_competitors_ID" = ID) %>% 
  left_join(race_sessions_unique %>% rename("race_sessions_ID" = ID) , by = c("SessionID" = "race_sessions_ID", "RaceID")) %>% 
  left_join(race_info_unique %>%  rename("race_info_ID" = ID), by = c("RaceID" = "race_info_ID"))

#Question 1----

regulars <- clean_racing_df %>% 
  filter(season == 2024) %>% 
  mutate(total_races = n_distinct(RaceID)) %>% 
  group_by(fullName) %>% 
  summarize(races = n_distinct(RaceID), total_races = max(total_races), participation_rate = races/total_races) %>% 
  filter(participation_rate > 0.75) %>% 
  mutate(regular = "Yes")

#Question 2----
clean_racing_df %>% 
  filter(season == 2024) %>% 
  left_join(regulars %>% select(fullName,regular), by = c("fullName")) %>% 
  filter(regular == "Yes") %>% 
  filter(str_detect(CategoryString, "A Feature"))  

test <- race_sessions_unique %>% 
  filter(season == 2024) %>% 
  
