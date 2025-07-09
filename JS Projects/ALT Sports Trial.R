#Part 1 ----
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(gt)
library(gtExtras)
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

regulars %>% 
  select(fullName,participation_rate) %>% 
  arrange(-participation_rate) %>% 
  gt() %>% 
  cols_align(align = "center") %>%
  cols_label(fullName = "Driver",participation_rate = "Participation Rate") %>% 
  gt_theme_538(quiet = TRUE) %>% 
  tab_header(
    title = md("Regular Drivers")
  ) %>% 
  tab_options(
    heading.align = "center"    # centers both title and subtitle
  )

#Question 2----
regular_racing_df <- clean_racing_df %>% 
  filter(season == 2024) %>% 
  left_join(regulars %>% select(fullName,regular), by = c("fullName"))

race_results <- regular_racing_df %>% 
  filter(regular == "Yes") %>% 
  filter(str_detect(Sessions_Name,"A Feature")) %>% #check this
  group_by(fullName) %>% 
  summarize(wins = mean(Position == 1), top_3 = mean(Position <=3), top_5 = mean(Position <= 5), top_10 = mean(Position <= 10))

race_results %>% 
  gt() %>% 
  fmt_number(
    columns = c(wins, top_3, top_5, top_10),
    decimals = 2
  ) %>% 
  cols_align(align = "center") %>%
  cols_label(fullName = "Driver",wins = "Win Rate", top_3 = "Top 3 Rate", top_5 = "Top 5 Rate", top_10 = "Top 10 Rate") %>% 
  gt_theme_538(quiet = TRUE) %>% 
  gt_hulk_col_numeric(columns = -fullName) %>% 
  tab_header(
    title = md("Regular Drivers Race Results")
  ) %>% 
  tab_options(
    heading.align = "center"    # centers both title and subtitle
  )


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

#Part 2----
f1_2024 <- read.csv("f1_2024.csv") %>% mutate(position = as.numeric(position))#results
f1_2025 <- read.csv("f1_2025.csv") %>% mutate(position = as.numeric(position)) #results
f1_odds_2024 <- read.csv("f1_odds_24.csv") #betting odds


#Question 1 ----
clean_odds <- f1_odds_2024 %>% mutate(probability = 1/EU.Odds) %>% group_by(Race,Market) %>% mutate(total_market_probability = sum(probability)) %>% ungroup() %>% 
  mutate(no_vig_prob = ifelse(Market == "Top 6", probability * (6/total_market_probability),
                              ifelse(Market == "Top 10",probability * (10/total_market_probability),
                              ifelse(Market == "Event Podium", probability * (3/total_market_probability),
                                     probability/total_market_probability))))

results_vs_prob <- clean_odds %>%   filter(Market %in% c("Event Winner","Event Podium","Top 6","Top 10")) %>%  filter(Name.Result != "") %>% 
  left_join(f1_2024 %>% filter(sessionName == "Race") %>% select(season, eventName,eventNum,driverName,position) %>% mutate(driverName = ifelse(driverName == "Alexander Albon", "Alex Albon", 
                                                                                                                                                ifelse(driverName == "Guanyu Zhou", "Zhou Guanyu",driverName))) %>% 
              group_by(driverName) %>% mutate(races_completed = sum(position >0,na.rm = T)) %>%  ungroup(), 
            by = c("Season" = "season", "RaceStop" = "eventNum", "Name.Result" = "driverName")) %>% 
  mutate(result = ifelse(Market == "Event Winner" & position == 1, 1, ifelse(Market == "Event Podium" & position <= 3, 1, ifelse(Market == "Top 6" & position <= 6, 1,
                                                                                                                                ifelse(Market == "Top 10" & position <= 10, 1,0))))) %>% 
  filter(!is.na(result),!is.na(no_vig_prob)) %>% 
  arrange(RaceStop) %>% 
  group_by(Name.Result,Market) %>% 
  mutate(total_prob = cumsum(no_vig_prob), actual_result = cumsum(result), races = max(races_completed, na.rm = T)) %>% 
  filter(races>=18)

results_vs_prob %>% pivot_longer(c(total_prob,actual_result), names_to = "Type", values_to = "Values") %>% 
  filter(Name.Result == "Lewis Hamilton") %>% 
  ggplot(aes(x = RaceStop, y = Values, color = Market, linetype = Type))+
  geom_line(lwd = 4)+
  labs(x = "Race Number", y = "Race Results", title = "Lewis Hamilton Race Results vs No Vig Probability by Race", 
       caption = "Lines de-vigged through proportional scaling approach")+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold", size = 20),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 26),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 22),
        plot.caption = element_text(colour = "white", size = 14),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(face = "bold", colour = "white",size = 22),
        axis.title = element_text(color = "white", size = 24),
        panel.border = element_rect(colour = "white", fill = NA, size = 1))+
  scale_linetype_manual(name   = "Series",
                        values = c(total_prob   = "dashed",
                                   actual_result = "solid"),
                        labels = c("Implied No Vig Probability",
                                   "Cumulative Results")) +
  scale_color_manual( values = c("Event Podium"   = "purple",  # cyan
                       "Event Winner"   = "red",  # magenta
                       "Top 10"  = "orange",  # yellow
                       "Top 6"   = "darkgrey"   # lime green
  ))+
  guides(
    linetype = guide_legend(
      override.aes = list(
        size = 5,           # make legend lines thicker
        color = "white",
        lwd = 1
      )),
    color = guide_legend(
      override.aes = list(
        size = 5,           # make legend lines thicker
        lwd = 1
      ))
  )

#Double check de-vig to make sure sums properly
#Question 2 ----
max_odds <- clean_odds %>%   filter(Market %in% c("Event Winner","Event Podium","Top 6","Top 10")) %>%  filter(Name.Result != "") %>% 
  left_join(f1_2024 %>% filter(sessionName == "Race") %>% select(season, eventName,eventNum,driverName,position) %>% mutate(driverName = ifelse(driverName == "Alexander Albon", "Alex Albon", 
                                                                                                                                                ifelse(driverName == "Guanyu Zhou", "Zhou Guanyu",driverName))), 
              by = c("Season" = "season", "RaceStop" = "eventNum", "Name.Result" = "driverName")) %>% 
  mutate(result = ifelse(Market == "Event Winner" & position == 1, 1, ifelse(Market == "Event Podium" & position <= 3, 1, ifelse(Market == "Top 6" & position <= 6, 1,
                                                                                                                                 ifelse(Market == "Top 10" & position <= 10, 1,0))))) %>% 
  group_by(Market) %>% 
  summarize(max_odds = max(EU.Odds[result == 1],na.rm = T))
#Question 3----
results_25 <- f1_2025 %>% 
  filter(sessionName == "Race") %>% 
  group_by(driverName) %>% 
  summarize(wins = mean(position == 1,na.rm =T), top_3 = mean(position <=3,na.rm =T), 
            top_5 = mean(position <= 5,na.rm =T), top_10 = mean(position <= 10,na.rm =T),
            races = n_distinct(eventKey))
#Max pedigree, firing of director --> lean into variance
#Charles Leclerc similar top 3, top 5, top 10 rate to heavy favorites but way off rest pricing