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
  distinct(ID, .keep_all = TRUE) %>% 
  rename("Sessions_Name" = Name)

race_competitors_unique <- session_competitors %>% 
  distinct(ID, .keep_all = TRUE)

clean_racing_df <- race_competitors_unique %>% 
  rename("race_competitors_ID" = ID) %>% 
  left_join(race_sessions_unique %>% rename("race_sessions_ID" = ID) , by = c("SessionID" = "race_sessions_ID", "RaceID")) %>% 
  left_join(race_info_unique %>%  rename("race_info_ID" = ID), by = c("RaceID" = "race_info_ID")) %>% 
  mutate(Position = as.numeric(Position))

#Question 1----

regulars <- clean_racing_df %>% 
  filter(season == 2024) %>% 
  mutate(total_races = n_distinct(RaceID)) %>% 
  group_by(fullName) %>% 
  summarize(races = n_distinct(RaceID), total_races = max(total_races), participation_rate = races/total_races) %>% 
  filter(participation_rate > 0.75) %>% 
  mutate(regular = "Yes")

#Question 2----
regular_racing_df <- clean_racing_df %>% 
  filter(season == 2024) %>% 
  left_join(regulars %>% select(fullName,regular), by = c("fullName"))

race_results <- regular_racing_df %>% 
  filter(regular == "Yes") %>% 
  filter(str_detect(Sessions_Name,"A Feature")) %>% #check this
  group_by(fullName) %>% 
  summarize(wins = mean(Position == 1), top_3 = mean(Position <=3), top_5 = mean(Position <= 5), top_10 = mean(Position <= 10))


#Question 3----
#ChatGPT assistance
h2h_all <- regular_racing_df %>%
  filter(regular == "Yes") %>% 
  filter(str_detect(Sessions_Name,"A Feature")) %>%
  select(RaceID,SessionID,fullName,Position) %>% 
  inner_join(regular_racing_df %>%   filter(regular == "Yes") %>% 
               filter(str_detect(Sessions_Name,"A Feature")) %>% select(RaceID,SessionID,fullName,Position), by = c("RaceID","SessionID"), suffix = (c(".A",".B"))) %>% 
  filter(fullName.A != fullName.B) %>% 
  mutate(
    A_beats_B = case_when(
      Position.A < Position.B              ~ 1,
      Position.A > Position.B              ~ 0,
      Position.A == Position.B             ~ NA_real_,  # tie
      TRUE                                  ~ NA_real_
    )
  )

h2h_matrix <- h2h_all %>%
  group_by(fullName.A, fullName.B) %>%
  summarise(
    races_together = n(),
    wins           = sum(A_beats_B, na.rm = TRUE),
    losses         = sum(1 - A_beats_B, na.rm = TRUE),
    win_pct        = wins / races_together,
    .groups        = "drop"
  )

#Question 4----
non_regular_results <- regular_racing_df %>% 
  filter(is.na(regular)) %>% 
  filter(str_detect(Sessions_Name,"A Feature")) %>% #check this
  group_by(fullName) %>% 
  summarize(wins = mean(Position == 1), top_3 = mean(Position <=3), top_5 = mean(Position <= 5), top_10 = mean(Position <= 10), 
            races_count = n_distinct(RaceID)) %>% 
  filter(fullName != "")

#Question 5
all_results <- regular_racing_df %>% 
  filter(str_detect(Sessions_Name,"A Feature")) %>% #check this
  group_by(fullName) %>% 
  summarize(wins = mean(Position == 1), top_3 = mean(Position <=3), top_5 = mean(Position <= 5), top_10 = mean(Position <= 10), 
            races_count = n_distinct(RaceID)) %>% 
  filter(fullName != "")

start_list_results <- start_list %>% 
  left_join(all_results, by = "fullName") %>% 
  mutate(composite_metric = wins*10+top_3*7+top_5*5+top_10)
#Based on last year's perfromance top contenders are Brad Sweet, Kyle Larson, and Tyler Courtney as the major favorites, then Rico Abreu and Corey Day

#Question 6
regular25 <- start_list %>% filter(isRegular == 1) %>% 
  left_join(regulars, by = "fullName") %>% 
  filter(is.na(regular)) %>% 
  left_join(all_results, by = "fullName")
#Aaron Reutzel has experience in Ridge & Sons Racing