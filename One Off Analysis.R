library(ggimage)
library(gt)
library(nflfastR)
library(tidyverse)
library(gtExtras)
library(nflplotR)

pbp <- load_pbp(2010:2023)

nfl99all <- load_pbp(1999:2024)
nfl99 <- nfl99all %>% 
  filter(pass == 1 | rush == 1) %>% 
  mutate(explosive = ifelse((yards_gained>20 & pass_attempt == 1) | (yards_gained >12 & (qb_scramble == 1 | rush == 1)),1,0),
         negative = ifelse(yards_gained < 0, 1,0))


pbp_rp <- pbp %>%
  filter(pass == 1 | rush == 1) %>%
  filter(!is.na(epa))

pass_efficiency <- pbp_rp %>%
  filter(pass ==1) %>%
  group_by(posteam,season) %>%
  summarize(passes = n(), pass_epa = mean(epa))

rush_efficiency <- pbp_rp %>%
  filter(rush ==1) %>%
  group_by(posteam, season) %>%
  summarize(rushes = n(), rush_epa= mean(epa))

pass_eff_2010_2023 <- pass_efficiency %>%
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr")) 

pass_eff_2010_2023 <- pass_eff_2010_2023 %>%
  mutate(combined_year_name = ifelse(season == 2023, paste("*",paste(team_nick,season, sep = " "), sep = ""),paste(team_nick,season, sep = " ")))

top_n_pass_eff_2010_2023 <- pass_eff_2010_2023 %>%
  arrange(-pass_epa) %>%
  head(25)

top_n_pass_eff_2010_2023 %>%
  ggplot(aes(x = pass_epa, y = fct_reorder(combined_year_name,pass_epa)))+
  geom_bar(aes(fill = team_color, color = team_color2), stat = "identity", alpha = 0.9)+
  scale_color_identity(aesthetics = c("fill","color"))+
  theme_bw()+
  geom_image(aes(image = team_logo_espn, x = ifelse(pass_epa>0, pass_epa+0.01, pass_epa-0.01)),
             asp = 16/9, size = 0.035)+
  labs(x = "EPA Per Pass Play", y = "Team & Season", title = "Top 25 Passing Offenses Based On EPA Per Pass Since 2010", subtitle = "The 49ers are off to a historic passing start" 
       ,caption = "*Current Season                                                                              By Callan Capitolo | @CapAnalytics7")+
  theme(panel.grid.major = element_blank())+
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 16, hjust = 0.5))
ggsave("49ersPassingOffense.png", width = 14, height =10, dpi = "retina")

#Bears Analysis

pbp_rp_2023 <- pbp_2023 %>%
  filter(pass == 1 | rush == 1) %>%
  filter(!is.na(epa))

pbp3_p_2023_chi <- pbp_2023 %>%
  filter(pass == 1) %>%
  filter(!is.na(epa)) %>%
  filter(posteam =="CHI") %>%
  filter(week %in% c(1,2,3)) %>%
  group_by(shotgun) %>%
  summarize(air_yards = mean(air_yards, na.rm = T), YAC = mean(yards_after_catch, na.rm = T), sack = mean(sack), count = n(), epa_play = mean(epa),qb_hit = mean(qb_hit))

pbp5_p_2023_chi <- pbp_2023 %>%
  filter(pass == 1) %>%
  filter(!is.na(epa)) %>%
  filter(posteam =="CHI") %>%
  filter(week %in% c(4,5)) %>%
  group_by(shotgun) %>%
  summarize(air_yards = mean(air_yards, na.rm = T), YAC = mean(yards_after_catch, na.rm = T), sack = mean(sack), count = n(), epa_play = mean(epa),qb_hit = mean(qb_hit))


bears_opp = pbp_2023 %>%
  filter(pass == 1) %>%
  filter(!is.na(epa)) %>%
  filter(shotgun == 1) %>%
  group_by(defteam) %>%
  summarize(air_yards = mean(air_yards, na.rm = T), YAC = mean(yards_after_catch, na.rm = T), sack = mean(sack), count = n(), epa_play = mean(epa),qb_hit = mean(qb_hit))%>%
  arrange(epa_play)




chi_pbp_2023<- pbp_rp_2023 %>%
  filter(posteam =="CHI")

chi_pass_efficiency <- chi_pbp_2023 %>%
  filter(pass ==1) %>%
  group_by(week,defteam,posteam) %>%
  summarize(passes = n(), pass_epa = mean(epa), cpoe = mean(cpoe, na.rm = TRUE))

chi3_pass_efficiency <- chi_pbp_2023 %>%
  filter(pass ==1, week %in% c(1,2,3)) %>%
  group_by(posteam) %>%
  summarize(passes = n()/3, pass_epa = mean(epa), cpoe = mean(cpoe, na.rm = TRUE))

chi5_pass_efficiency <- chi_pbp_2023 %>%
  filter(pass ==1, week >= 4) %>%
  group_by(posteam) %>%
  summarize(passes = n()/2, pass_epa = mean(epa), cpoe = mean(cpoe, na.rm = TRUE))

chi_rush_efficiency <- chi_pbp_2023 %>%
  filter(rush ==1) %>%
  group_by(week,defteam) %>%
  summarize(rushes = n(), rush_epa= mean(epa))

chi_efficiency <- chi_pass_efficiency %>%
  inner_join(chi_rush_efficiency)

chi_efficiency <- chi_efficiency %>%
  left_join(teams_colors_logos, by = c("defteam" = "team_abbr")) 

chi_efficiency %>%
  ggplot(aes(x = pass_epa, y = rush_epa, label = week))+
  geom_text()+
  theme_bw()+
  labs(x = "EPA Per Pass",
       y = "EPA Per Rush", title = "Bears Offensive Efficiency by Week", subtitle = "Bears Have Seen a Significant Improvement the Past 2 Weeks",
       caption = "Callan Capitolo | @CapAnalytics7")+
  theme_bw()

#Better Way to show efficency by week
bears_gt<-ungroup(chi_efficiency) %>%
  arrange(week) %>%
  select(week,pass_epa, rush_epa, team_wordmark)%>%
  mutate(pass_epa = round(pass_epa,2), rush_epa = round(rush_epa,2))%>%
  gt() %>%
  cols_align(align = "center") %>%
  gtExtras::gt_img_rows(team_wordmark) %>%
  cols_label(week = "Week",
              pass_epa = "EPA Per Pass",
              team_wordmark = "Opponent",
              rush_epa = "EPA Per Rush") %>%
  gtExtras::gt_theme_538() 
gtsave(bears_gt, "Bears_week_efficiency.png") 

qb_epa_play %>%
  arrange(-epa_play) %>%
  mutate(rank = row_number()) %>%
  dplyr::select(rank,name, team_wordmark,pass_attempts,pass_rate, epa_play)%>%
  mutate(pass_rate = 100*round(pass_rate,3), epa_play = round(epa_play,2)) %>%
  gt() %>%
  cols_align(align = "center") %>%
  gtExtras::gt_img_rows(team_wordmark) %>%
  cols_label( rank = "Rank",
              name = "Quarterback",
              team_wordmark = "",
              pass_attempts = "Pass Attempts",
              pass_rate = "Pass Rate %",
              epa_play = "EPA Per Play") %>%
  gtExtras::gt_theme_538() %>%
  gtExtras::gt_hulk_col_numeric(epa_play)
gtsave(qb_gt, "qb_efficiency_week5.png") 
  



#Brock vs Jimmy ----
nfl99 %>% 
  filter(posteam == "SF") %>% 
  mutate(wp = round(wp,2)) %>% 
  group_by(posteam,wp,id) %>% 
  summarize(name = first(name), epa = mean(epa,na.rm = T)) %>% 
  filter(name %in% c("J.Garoppolo", "B.Purdy")) %>% 
  ggplot(aes(x = wp, y = epa,color = name)) + 
  geom_smooth(se = FALSE, span =0.6, method = "loess", lwd = 3)+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 12),
        plot.caption = element_text(colour = "white", size = 10),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(face = "bold", colour = "white",size = 12),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1))+
  labs(x = "Win Probability", y = "EPA/Play (Pass or Rush by QB)", title = "How do Purdy and Garoppolo Perform at Different Win Probabilities?",
       caption = "@CapAnalytics7 | nflfastR", subtitle = "Purdy Performs Better in High Win Probability Situations than Garoppolo")
ggsave("BrockvPurdy.png", width = 14, height =10, dpi = "retina")

nfl99 %>% 
  filter(posteam == "SF") %>% 
  mutate(wp = round(wp,2)) %>%
  left_join(nfl99 %>% 
              mutate(pur_g = ifelse(name %in% c("B.Purdy","J.Garoppolo"),name, 0)) %>% 
              group_by(game_id) %>% 
              summarize(qb = max(pur_g)), by = c("game_id")) %>% 
  filter(qb != 0) %>% 
  # group_by(qb,wp) %>% 
  group_by(qb,wp, season) %>%
  filter(qb == "B.Purdy") %>% 
  summarize(pass = mean(pass_oe, na.rm = T)) %>% 
  # ggplot(aes(x= wp, y = pass, color = as.factor(qb)))+
  ggplot(aes(x= wp, y = pass, color = as.factor(season)))+
  geom_smooth(se = F, method = "loess", span = .6, lwd = 3)+
  # scale_color_brewer(palette = "Set2")+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 12),
        plot.caption = element_text(colour = "white", size = 10),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(face = "bold", colour = "white",size = 12),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1))+
  labs(x = "Win Probability", y = "Pass Rate Over Expected", title = "Has Shanahan Shown More Trust in Purdy as He Has Developed?",
       caption = "@CapAnalytics7 | nflfastR", subtitle = "In high win probability situations, Shanahan has let Purdy air the ball out more in years past")
ggsave("BrockvPurdyRate.png", width = 14, height =10, dpi = "retina")


#Running Efficinecy by Win Prob and Exp Pass---
nfl99 %>% 
  filter(rush == 1, qb_kneel!= 1) %>% 
  mutate(wp = round(wp,2)) %>%
  group_by(wp) %>% 
  summarize(rush_epa = mean(epa,na.rm =T)) %>% 
  ggplot(aes(x = wp, y = rush_epa)) + 
  geom_smooth(se = FALSE, method = "loess" ,span = 0.6, lwd = 3)+
  geom_point(color = "white", alpha = 0.2)+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 12),
        plot.caption = element_text(colour = "white", size = 10),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(face = "bold", colour = "white",size = 12),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1))+
  labs(x = "Win Probability", y = "EPA/Rush", title = "Does Rushing Efficiency Change With Win Probability?",
       caption = "@CapAnalytics7 | nflfastR", subtitle = "Rushing is most efficient when the offense has a low win probablity")
ggsave("RushvsWP.png", width = 14, height =10, dpi = "retina")

nfl99 %>% 
  filter(rush == 1, qb_kneel!= 1) %>% 
  mutate(xpass = round(xpass,2)) %>%
  group_by(xpass) %>% 
  summarize(rush_epa = mean(epa,na.rm =T)) %>% 
  ggplot(aes(x = xpass, y = rush_epa)) + 
  geom_smooth(se = FALSE, method = "loess" ,span = 0.6, lwd = 3)+
  geom_point(color = "white", alpha = 0.2)+
  ylim(-0.4,0.05)+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 12),
        plot.caption = element_text(colour = "white", size = 10),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(face = "bold", colour = "white",size = 12),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1))+
  labs(x = "Pass Play Probability", y = "EPA/Rush", title = "Does Rushing Efficiency Change With the Likelihood of a Pass Play?",
       caption = "@CapAnalytics7 | nflfastR", subtitle = "Rushing is most efficient when the offense has low win probablity")
ggsave("RushXP.png", width = 14, height =10, dpi = "retina")

nfl99 %>% 
  filter(rush == 1, qb_kneel!= 1) %>% 
  mutate(xpass = round(xpass,2)) %>%
  filter(posteam == "SF", season>=2017) %>% 
  group_by(xpass) %>% 
  summarize(rush_epa = mean(epa,na.rm =T)) %>% 
  ggplot(aes(x = xpass, y = rush_epa)) + 
  geom_smooth(se = FALSE, method = "loess" ,span = 0.6, lwd = 3)+
  geom_point(color = "white", alpha = 0.2)+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 12),
        plot.caption = element_text(colour = "white", size = 10),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(face = "bold", colour = "white",size = 12),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1))+
  labs(x = "Pass Play Probability", y = "EPA/Rush", title = "Does Rushing Efficiency Change With the Likelihood of a Pass Play?",
       caption = "@CapAnalytics7 | nflfastR", subtitle = "Rushing is most efficient when the offense has low win probablity")
nfl99 %>% 
  filter(rush == 1, qb_kneel!= 1) %>% 
  mutate(wp = round(wp,2)) %>%
  filter(posteam == "SF", season>=2017) %>% 
  group_by(wp) %>% 
  summarize(rush_epa = mean(epa,na.rm =T)) %>% 
  ggplot(aes(x = wp, y = rush_epa)) + 
  geom_smooth(se = FALSE, method = "loess" ,span = 0.6, lwd = 3)+
  geom_point(color = "white", alpha = 0.2)+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 12),
        plot.caption = element_text(colour = "white", size = 10),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(face = "bold", colour = "white",size = 12),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1))+
  labs(x = "Win Probability", y = "EPA/Rush", title = "Does Shanahan's Offense Rushing Efficiency Change With Win Probability?",
       caption = "@CapAnalytics7 | nflfastR", subtitle = "Shanahan's offense's rushing efficiency stays fairly constant in respect to their win probablity")
ggsave("ShanWP.png", width = 14, height =10, dpi = "retina")



#EPA Decay Rate
# Load necessary libraries
library(nflfastR)
library(dplyr)
library(purrr)
library(Metrics) # For RMSE calculation

# Load NFL data
pbp_data <- load_pbp(1999:2023)# Load data from 1999 to 2023

pbp_data <- pbp_data %>% 
  filter(season>2010)
calculate_adjusted_epa <- function(data, decay_rate) {
  data %>%
    group_by(season, posteam, week) %>%
    summarize(
      weekly_epa = mean(epa, na.rm = TRUE), .groups = 'drop'
    ) %>%
    arrange(season, week) %>%
    group_by(posteam,season,week) %>%
    mutate(
      decay_factor = decay_rate ^ (max(week) - week),
      weighted_epa = weekly_epa * decay_factor
    ) 
  # %>%
  #   summarize(
  #     adjusted_epa = sum(weighted_epa) / sum(decay_factor)
  #   )
}

# Define decay rates to test
decay_rates <- seq(0.9, 1, by = 0.01)

# Cross-validation function
# Adjusted cross-validation function
cross_validate_decay <- function(data, decay_rate) {
  # Calculate adjusted EPA for the chosen decay rate
  adjusted_data <- calculate_adjusted_epa(data, decay_rate)
  
  # Shift the adjusted EPA forward by one week for validation
  adjusted_data <- adjusted_data %>%
    mutate(week = week +1)
  
  # Join with actual EPA data for the next week
  validation_data <- data %>%
    group_by(season, posteam, week) %>%
    summarize(
      actual_epa = mean(epa, na.rm = TRUE), .groups = 'drop'
    ) %>%
    left_join(adjusted_data, by = c("season", "posteam", "week"))
  # Filter out any rows where `next_week_adjusted_epa` is NA
  validation_data <- validation_data %>% 
  filter(!is.na(adjusted_epa))
  
  # Calculate RMSE between adjusted EPA (from prior week) and actual EPA
  rmse(validation_data$adjusted_epa, validation_data$actual_epa)
}

# Apply cross-validation for each decay rate
decay_results <- map_df(decay_rates, ~ tibble(
  decay_rate = .x,
  rmse = cross_validate_decay(pbp_data, .x)
))

# Find the decay rate with the lowest RMSE
best_decay_rate <- decay_results %>%
  arrange(rmse) %>%
  slice(1)

print(best_decay_rate)


#Method of Obtaining Ball----
nfl99 %>% 
  filter(!is.na(drive_start_transition)) %>% 
  mutate(drive_start_transition = ifelse(grepl("blocked",tolower(drive_start_transition)),"Blocked Kick/Punt",
                                         ifelse(grepl("safety",tolower(drive_start_transition)), "Safety", 
                                         ifelse(grepl("missed",tolower(drive_start_transition)),"Missed FG", 
                                         ifelse(grepl("interception|fumble|muffed",tolower(drive_start_transition)),"Turnover",drive_start_transition))))) %>% 
  group_by(toupper(drive_start_transition)) %>% 
  summarize(count = n(),epa_play = mean(epa,na.rm = T), success_rate = mean(success,na.rm =T)) %>% 
  filter(count>100) %>% 
  arrange(-epa_play) %>% 
  ggplot(aes(x = reorder(`toupper(drive_start_transition)`,-epa_play), y = epa_play))+
  geom_bar(stat = "identity")+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 12),
        plot.caption = element_text(colour = "white", size = 10),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(face = "bold", colour = "white",size = 11),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1))+
  labs(x = "How Offense Acquired Ball", y = "EPA/Play on Drive", title = "Does the Way an Offense Acquire the Ball Influence Drive Success?",
       caption = "@CapAnalytics7 | nflfastR", subtitle = "The way an Offense Acquires the Ball does Have an Influence on Offensive Efficiency")
ggsave("EPABallAcq.png", width = 16, height =10, dpi = "retina")

nfl99 %>% 
  filter(!is.na(drive_start_transition)) %>% 
  mutate(drive_start_transition = ifelse(grepl("blocked",tolower(drive_start_transition)),"Blocked Kick/Punt",
                                         ifelse(grepl("safety",tolower(drive_start_transition)), "Safety",
                                                ifelse(grepl("missed",tolower(drive_start_transition)),"Missed FG",
                                                       ifelse(grepl("interception|fumble|muffed",tolower(drive_start_transition)),"Turnover",drive_start_transition))))) %>%
  group_by(toupper(drive_start_transition)) %>% 
  summarize(count = n(),epa_play = mean(epa,na.rm = T), success_rate = mean(success,na.rm =T)) %>% 
  filter(count>100) %>% 
  arrange(-epa_play) %>% 
  ggplot(aes(x = reorder(`toupper(drive_start_transition)`,-success_rate), y = success_rate))+
  geom_bar(stat = "identity")+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 12),
        plot.caption = element_text(colour = "white", size = 10),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(face = "bold", colour = "white",size = 11),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1))+
  labs(x = "How Offense Acquired Ball", y = "Success Rate on Drive", title = "Does the Way an Offense Acquire the Ball Influence Drive Success?",
       caption = "@CapAnalytics7 | nflfastR", subtitle = "The way an Offense Acquires the Ball does Have an Influence on Offensive Efficiency")
ggsave("SuccessRateBallAcq.png", width = 16, height =10, dpi = "retina")

mutate(drive_start = ifelse(yrdln == drive_start_yard_line, yardline_100,NA)) %>% 

comp_data<-nfl99 %>% 
  filter(!is.na(drive_start_transition)) %>% 
  mutate(drive_start_transition = tolower(ifelse(grepl("blocked",tolower(drive_start_transition)),"Blocked Kick/Punt",
                                         ifelse(grepl("safety",tolower(drive_start_transition)), "Safety", 
                                                ifelse(grepl("missed",tolower(drive_start_transition)),"Missed FG", 
                                                       ifelse(grepl("interception|fumble|muffed",tolower(drive_start_transition)),"Turnover",drive_start_transition)))))) %>% 
  filter(!is.na(epa))
anova_result <- (aov(success~drive_start_transition, data =comp_data))

tukey_result <- TukeyHSD(anova_result, "drive_start_transition")
summary(anova_result)

# View the results
print(tukey_result)

#Jets----
Jets<- pbp_rp %>% 
  filter(posteam == "NYJ") %>% 
  group_by(week) %>% 
  summarize( `Motion Rate` = mean(is_motion), `EPA/Play With Motion` = mean(epa[is_motion == T], na.rm = T), `EPA/Play With No Motion` = mean(epa[is_motion == F], na.rm = T), `Difference Between EPAs` = `EPA/Play With Motion` - `EPA/Play With No Motion`) %>% 
  mutate_if(is.numeric, ~round(.,2)) %>% 
  gt() %>% 
  cols_align(align = "center") %>%
  gtExtras::gt_theme_538() %>%
  gtExtras::gt_hulk_col_numeric(`Difference Between EPAs`) %>%
  tab_header(
    title = md("How Has the Jets Motion Usage Changed?"),
    subtitle = md("The Jets offense increased their motion rate under Downing and saw an improvement in their efficiency with motion")
  )
gtsave(Jets, "Jets.png")

#49ers Radar Plot
pbp23 <- load_pbp(2023:2024)
ftn_data <- nflreadr::load_ftn_charting(2023:2024) %>%
  select(-week, -season)
pbp23 <- pbp23 %>%
  left_join(ftn_data, by = c("game_id" = "nflverse_game_id",
                             "play_id" = "nflverse_play_id")) 

nfl23 <- pbp23 %>% 
  filter(pass == 1 | rush == 1) %>% 
  filter(qb_kneel== 0, qb_spike ==0) %>% 
  mutate(explosive = ifelse((yards_gained>20 & pass_attempt == 1) | (yards_gained >12 & (qb_scramble == 1 | rush == 1)),1,0),
         negative = ifelse(yards_gained < 0, 1,0))

team_stats <- nfl23 %>% 
  mutate(turnover = ifelse(interception == 0 & fumble_lost == 0, 0,1)) %>% 
  group_by(season,posteam) %>% 
  summarize(`EPA/Play`= mean(epa,na.rm =T), `Success Rate` = mean(success,na.rm =T) ,`EPA/Dropback` = mean(epa[pass == 1],na.rm = T), ADoT = mean(air_yards,na.rm = T), `EPA/Rush` = mean(epa[rush == 1], na.rm = T), CPOE = mean(cpoe,na.rm = T), `Explosive Play Rate` = mean(explosive, na.rm = T),
            `% of Yards From YAC` = sum(yards_after_catch,na.rm = T)/sum(yards_gained[complete_pass == 1],na.rm = T),`Negative Play Rate` = mean(negative, na.rm = T), 
            `Motion Rate`= mean(is_motion, na.rm = T), `Turnover Rate` = mean(turnover,na.rm = T),
            # `Catchable Ball Rate` = mean(is_catchable_ball[!is.na(air_yards) & is_throw_away == 0], na.rm = T),
            # `Contested Throw Rate` = mean(is_contested_ball[!is.na(air_yards)]), 
            `YAC Over Expected/Catch` = mean(yards_after_catch-xyac_mean_yardage,na.rm = T)
            ) %>% 
  mutate(name = paste (posteam,season, sep = " ")) %>% 
  ungroup %>% 
  select(-season,-posteam)

replace_with_values_and_ranks <- function(column) {
  values <- column
  ranks <- rank(column*-1,ties.method = "max")
}
team_rank <- as.data.frame(cbind(team_stats$name,(apply(team_stats %>% 
                                                                    select(-name,) %>% 
                                                                    mutate(`Negative Play Rate` = `Negative Play Rate` * -1, , `Turnover Rate` = `Turnover Rate` * -1), 2, replace_with_values_and_ranks)))) %>% 
  pivot_longer(cols = -`V1`, names_to = "statistic", values_to = "rank") %>% 
  mutate(rank = as.numeric(rank))

team_all_comb<- team_stats %>% pivot_longer(cols = c(-name), names_to = "statistic", values_to = "value") %>% 
  inner_join(team_rank, by = c("name" = "V1", "statistic"))

team1 <- "SF 2023" #Rework to be more flexible and incorporate more players
team2 <- "SF 2024"
#Production Stats ----


temp <- (360/n_distinct(team_all_comb$statistic))/2

myAng <- seq(-temp, -360+temp, length.out = n_distinct(team_all_comb$statistic))

ang <- ifelse(myAng < -90, myAng+180,myAng)

ang <- ifelse(ang < -90, ang+180, ang)


pizza_prod<- team_all_comb %>% 
  mutate(rank = max(rank) + 1 - rank) %>% 
  filter(name %in% c(team1,team2)) %>% 
  mutate(name = factor(name, levels = c(team1, team2))) %>% 
  arrange(name) %>% 
  ggplot(aes(x = statistic, y = rank, fill = name, label = value))+
  geom_bar(stat = "identity", position = position_identity(), alpha = 0.6)+
  # geom_label(color = "white", size=2.5, fontface="bold", show.legend = FALSE, position = position_jitterdodge())+
  coord_polar()+
  geom_bar(aes(y = max(team_all_comb$rank)/n_distinct(name)),stat = "identity", width =1, alpha = 0.1, fill = "grey")+
  geom_hline(yintercept = seq(1, max(team_all_comb$rank), by = max(team_all_comb$rank)),
             color = "white",
             size = 1)+
  geom_vline(xintercept = seq(.5, n_distinct(team_all_comb$statistic), by = 1),
             color = "white",
             size = .5)+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 8),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(face = "bold", size = 8, colour = "white"),
        axis.title = element_blank(),
        axis.text.x = element_text(face = "bold", size = 12, angle = ang)) +
  labs(x = NULL, y = NULL)+
  scale_fill_brewer(palette = "Set1")+
  annotate("text", x = (pi/12) * 2, y = seq(10,n_distinct(team_all_comb$name), by = 10), label = seq(10*n_distinct(team_all_comb$name)%/%10, 10 ,by = -10), hjust = 1.15, 
           color = "White", size = 5)


tab_prod <- team_all_comb %>%
  filter(name %in% c(team1,team2)) %>%
  mutate(value = round(value,3)) %>% 
  mutate(name = factor(name, levels = c(team1, team2))) %>% 
  pivot_wider(names_from = name, values_from = c(value,rank)) %>% 
  rename("Value1" = paste("value_",team1,sep = ""), "Rank1" = paste("rank_",team1,sep = ""), "Value2" = paste("value_",team2,sep = ""),"Rank2" = paste("rank_",team2,sep = "")) %>%
  select(statistic,Value1,Rank1,Value2,Rank2) %>% 
  arrange(statistic) %>% 
  gt() %>% 
  cols_align(align = "center") %>% 
  tab_spanner(label = team1, columns = c(Value1,Rank1)) %>% 
  tab_spanner(label = team2, columns = c(Value2,Rank2)) %>% 
  cols_label(Value1 = "Value", Rank1 = "Rank",Value2 = "Value",Rank2 = "Rank") %>% 
  gtExtras::gt_theme_538() %>% 
  tab_options(
    table.background.color = "black",        # Set the entire table background to black
    heading.background.color = "black",      # Set the header background to black
    column_labels.background.color = "black", # Set the column label background to black
    row_group.background.color = "black",    # Set row group background to black (if any)
    summary_row.background.color = "black",  # Set summary row background to black (if any)
    grand_summary_row.background.color = "black", # Set grand summary row background to black (if any)
    footnotes.background.color = "black",    # Set footnotes background to black
    source_notes.background.color = "black", # Set source notes background to black
    table.border.top.color = "black",        # Set table top border to black
    table.border.bottom.color = "black",     # Set table bottom border to black
    heading.border.bottom.color = "black",   # Set header bottom border to black
    column_labels.border.top.color = "black",# Set column label top border to black
    column_labels.border.bottom.color = "black" # Set column label bottom border to black
  ) %>% 
  gt_hulk_col_numeric(columns = c(Rank1,Rank2),reverse = TRUE) %>%
  tab_style(
    style = cell_text(size = px(16), weight = "bold", color = "white"),  # Change font and size for column labels
    locations = cells_column_labels(columns = everything())
  ) %>% 
  tab_style(
    style = cell_text(color = "white", size = px(16)),  # Change font and size for the body text
    locations = cells_body(columns = c(statistic, Value1, Value2))
  )

if (exists("f") && inherits(f, "ChromoteSession")) {
  try(f$shutdown(), silent = TRUE)
}

# Start a new session
f <- ChromoteSession$new()

gtsave(tab_prod, "prod_temp_table.png")

table_image <- image_read("prod_temp_table.png")
table_image_transparent <- image_transparent(table_image, "white")
table_grob <- rasterGrob(table_image, interpolate = TRUE)
spacer <- ggplot() + theme_void() + theme(panel.background = element_rect(fill = "black"))

pizza_prod + table_grob+spacer + plot_layout(ncol = 3, widths = c(6,3,.1))& 
  theme(
    plot.background = element_rect(fill = "black", color = "black"),
    panel.background = element_rect(fill = "black", color = "black"),
    plot.caption = element_text(size = 10, hjust = 0.5, face = "italic", color = "white"),
    plot.title = element_text(size =20,hjust = 0.5, face = "italic", color = "white"),
    plot.subtitle = element_text(hjust = 0.5, face = "bold", color = "white"))&
  plot_annotation(
    caption = glue("Compared to {n_distinct(qb_all_stats$name)} QB seasons with 60+ dropbacks in 2022-2023"),
    title = "2022-2023 QB Production Comparison",
    subtitle =  "@CapAnalytics7 | nflfastR"
  )
ggsave("QB_Comp_Prod.png", bg = "black", ,width = 14, height =10, dpi = "retina") 


#Fix this matchup----
currentweek <- load_schedules(2024) %>% 
  filter(week ==10)

matchups <- currentweek %>% 
  left_join(total_first_half, by = c("home_team" = "posteam"))

matchups <- matchups %>% 
  left_join(total_first_half, by = c("away_team" = "posteam"))

epa_1h_matchups <- matchups %>% 
  select(team_wordmark.x,def_epa.x,off_epa.x, team_wordmark.y,def_epa.y,off_epa.y) %>% 
  mutate(home_team_off_diff = off_epa.x - def_epa.y*-1) %>% 
  mutate(home_team_def_diff = def_epa.x*-1 - off_epa.y) %>% 
  mutate_if(is.numeric,round,2) %>% 
  mutate(EPAdifference = home_team_def_diff+home_team_off_diff) %>% 
  mutate_if(is.numeric,round,2)

first_half <- epa_1h_matchups %>%
  select(-home_team_off_diff,-home_team_def_diff, -EPAdifference) %>% 
  # arrange(-EPAdifference) %>%
  gt() %>%
  cols_align(align = "center") %>%
  gtExtras::gt_img_rows(team_wordmark.x) %>%
  gtExtras::gt_img_rows(team_wordmark.y) %>%
  cols_label(team_wordmark.x = "Home Team",
             team_wordmark.y = "Away Team",
             def_epa.x = "Home 1H Def EPA/Play",
             off_epa.x = "Home 1H Off EPA/Play",
             def_epa.y = "Away 1H Def EPA/Play",
             off_epa.y = "Away 1H Off EPA/Play") %>% 
  # EPAdifference = "Difference in Total 1H EPA/Play") %>%
  gtExtras::gt_theme_538() %>%
  # gtExtras::gt_hulk_col_numeric(EPAdifference) %>% 
  tab_header(
    title = md("1st Half Efficiency by Week 12 Matchups"),
    subtitle = md("Negative Defensive EPA is Good, Positive Offensive EPA is Good")
  )
gtsave(first_half, "FirstHalfEfficiency.png") 


#Passing Charts----
nfl99all %>%
  mutate(air_yards_bins = cut(air_yards,
                              breaks = c(-Inf, 0, 10, 20, Inf),
                              labels = c("<=0", "1-10", "11-20", "21+"))) %>%
  filter(posteam == "SF") %>% 
  filter(passer_player_name %in% c("B.Purdy", "J.Garoppolo") & !is.na(pass_location)) %>%
  group_by(passer_player_name, pass_location, air_yards_bins) %>%
  summarize(count_air = n(), epa_pass = mean(epa), success_rate = mean(success),
            YPA = mean(yards_gained),
            cmp = sum(complete_pass) / n()) %>%
  ungroup() %>%
  group_by(passer_player_name) %>%
  mutate(pct_throws = count_air / sum(count_air)) %>%
  ggplot(aes(x = pass_location, y = air_yards_bins, fill = pct_throws)) +
  geom_tile() +
  facet_wrap(~ passer_player_name) +
  geom_text(aes(label = paste("% Throws:", round(pct_throws, 2), "\nEPA:", round(epa_pass, 2),
                              "\nSuccess:", round(success_rate, 2), "\nYPA:", round(YPA, 2),
                              "\nComp %:", round(cmp, 2))),
            size = 3) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "Pass Location", y = "Air Yards", title = "Brock Purdy vs Jimmy Garropolo (Super Bowl Seasons Only)") +
  theme_minimal()
ggsave("JimmyGvsPurdy.png", width = 14, height =10, dpi = "retina")


#Rolling average----
pbp_rp %>% 
  filter(week == 22, posteam == "SF") %>%
  mutate(roll_epa = rollapply(epa, width = 15, align = "right", FUN = mean,partial = TRUE),
         play_num = row_number()) %>%
  ggplot(aes(x = play_num, y = roll_epa))+
  # geom_point()
  geom_line(linewidth =1 )+
  geom_vline(xintercept = 5,linetype = "dashed", color = 'red')+
  # geom_vline(xintercept = 45,linetype = "dashed", color = 'red')+
  labs(x = "Play Number", y = "49ers Offense Rolling EPA/Play for Previous 15 Plays", title = "49ers Rolling EPA/Play by Play Number for Super Bowl")+
  annotate("text", x = 9, y = 0.1, label = "CMC Fumble", color = "red",size =4)
# annotate("text", x = 48, y = -0.05, label = "Muffed Punt", color = "red",size =4)
ggsave("SuperBowlBreakdown.png", width = 14, height =10, dpi = "retina")

#Recovering from blowout?----

bet_data <- nfl99 %>%
  group_by(game_id) %>% 
  summarize(year = max(season), week = max(week),homescore = max(home_score), away_score = max(away_score),spread = max(spread_line), result = max(result), bet_total = max(total_line), result_total = max(total),location = max(location),
            home = max(home_team), away = max(away_team)) %>% 
  mutate(favorite = ifelse(spread >= 0, "home", "away")) %>% 
  mutate(underdog = ifelse(spread >= 0, "away", "home")) %>% 
  mutate(favorite = ifelse(spread >= 0, "home", "away"), favorite_cover = ifelse((favorite == "home" & result > spread) |(favorite == "away" & result < spread),1,0 )) %>% 
  mutate(underdog = ifelse(spread >= 0, "away", "home"), underdog_cover = ifelse((underdog == "home" & result > spread) |(underdog == "away" & result < spread),1,0 )) %>% 
  pivot_longer(cols = c("favorite_cover", "underdog_cover"), values_to = "cover", names_to = "side") %>% 
  mutate(team_bet = ifelse((side == "favorite_cover" & favorite == "away") | (side == "underdog_cover" & underdog == "away"),away,home)) %>% 
  # mutate(coverby = ifelse(side == "favorite_cover" & favorite == "home"),result - spread) %>%
  mutate(coverby = result - spread) %>% 
  mutate(coverby = ifelse((cover == 0 & coverby >0) |(cover == 1 & coverby <0) , coverby*-1,coverby))



following_week <- bet_data %>% 
  group_by(year,week,team_bet) %>% summarize(next_week_cover = sum(cover)) %>% 
  mutate(week = week-1)

bet_data_comb <- bet_data %>% 
  left_join(following_week,by = c("year","week","team_bet"))

corr <- bet_data_comb %>% 
  group_by(coverby) %>% 
  summarize(following_week_cover = mean(next_week_cover,na.rm = T), count = n()) %>% 
  filter(count > 20) %>% 
  ggplot(aes(x = coverby, y = following_week_cover))+
  geom_point(aes(size =count),color = "white")+
  geom_smooth(se = FALSE, linewidth = 3)+
  labs(x = "Points Covered By", y = "Following Week Cover Rate", title= "Do NFL Teams Who Get Blown Out Recover the Following Week ATS?",
       subtitle = "Teams on either side of a blow out tend to perform better the following week than teams with closer contested spreads", 
       caption = "Dotted line represents break even win rate for -110 spread; data since 1999                      @CapAnalytics7 | nflfastR")+
  theme(legend.position = "none",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 12),
        plot.caption = element_text(colour = "white", size = 10),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(face = "bold", colour = "white",size = 12),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1))+
  geom_hline(yintercept = 0.5236, color = "white",linetype = "dashed")
ggsave("Blowouts.png", width = 14, height =10, dpi = "retina")

#Success by Air Yards----
nfl99 %>% 
  filter(!is.na(air_yards)) %>% 
  group_by(air_yards) %>% 
  summarize(epa_pass = mean(epa,na.rm = T), count = n()) %>% 
  filter(count > 200) %>% 
  ggplot(aes(x = air_yards, y = epa_pass))+
  geom_point(color = "white")+
  geom_smooth(se = FALSE, linewidth = 3)+
  labs(x = "Air Yards", y = "EPA/Pass", title= "How Does Pass Efficiency Change By Air Yards?",
       subtitle = "Between 11-35 Air Yards Pass Efficiency Reaches a Plateau", 
       caption = "Minimum 200 passes,data since 2006                      @CapAnalytics7 | nflfastR")+
  theme(legend.position = "none",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 12),
        plot.caption = element_text(colour = "white", size = 10),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(face = "bold", colour = "white",size = 12),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1))
ggsave("AirYardsEff.png", width = 14, height =10, dpi = "retina")

#Lucky----
lucky <- pbp_rp %>% 
  filter(qb_kneel == 0, qb_spike == 0,pass == 1 |rush  == 1) %>% 
  filter(season == 2024) %>%
  group_by(game_id) %>% 
  summarize(homescore = max(home_score),homeepa = mean(epa[home_team == posteam], na.rm = T), hometeam = max(home_team),
            awayscore = max(away_score),awayepa = mean(epa[away_team == posteam], na.rm = T), awayteam = max(away_team)) %>% 
  mutate(winningteam = ifelse(homescore == awayscore, NA,ifelse(homescore>awayscore,hometeam,awayteam)),
         losingteam = ifelse(homescore == awayscore, NA,ifelse(homescore>awayscore,awayteam,hometeam)),
         winning_epa = ifelse(hometeam == winningteam, homeepa, awayepa),
         losing_epa = ifelse(hometeam == losingteam, homeepa,awayepa)) %>% 
  # mutate(luckywin = ifelse(homescore == awayscore,NA,ifelse((homescore>awayscore & homeepa<awayepa)|(homescore<awayscore & homeepa>awayepa),1,0)))
  mutate(luckyepawin = ifelse(homescore == awayscore, NA,ifelse(winning_epa<losing_epa,1,0)))
unlucky_teams <- lucky %>% 
  group_by(losingteam) %>% 
  summarize(pct_unlucky_losses = mean(luckyepawin,na.rm = T),
            total_unlucky_win = sum(luckyepawin,na.rm = T))
lucky_teams <- lucky %>% 
  group_by(winningteam) %>% 
  summarize(pct_lucky_losses = mean(luckyepawin,na.rm = T),
            total_lucky_win = sum(luckyepawin,na.rm = T))

check <- pbp %>% 
  filter(!(play_type %in% c("run","pass"))) %>% 
  select(desc,epa,posteam,defteam,play_type,pass, rush)

lucky_data <- pbp %>% 
  group_by(game_id,posteam) %>% 
  summarize(total_epa_off = sum(epa,na.rm = T), off_plays = n(), homescore = max(home_score), hometeam = max(home_team),
            awayscore = max(away_score), awayteam = max(away_team),off_total_success = sum(success,na.rm = T)) %>% 
  mutate(winningteam = ifelse(homescore == awayscore, NA,ifelse(homescore>awayscore,hometeam,awayteam)),
                              losingteam = ifelse(homescore == awayscore, NA,ifelse(homescore>awayscore,awayteam,hometeam))) %>% 
  filter(!is.na(posteam)) %>% 
  left_join(pbp %>% 
              group_by(game_id,defteam) %>% 
              summarize(total_epa_def = sum(epa,na.rm = T)*-1, def_plays = n(), total_success = sum(success,na.rm = T)) %>%
              mutate(def_total_success = def_plays-total_success) %>% 
              filter(!is.na(defteam)), by = c("game_id", "posteam" = "defteam")) %>% 
  mutate(total_epa = total_epa_off+total_epa_def,
         all_total_success = def_total_success+off_total_success) %>% 
  mutate(result = ifelse(is.na(winningteam), NA, ifelse(winningteam == posteam, "win","loss"))) %>% 
  select(game_id, winningteam,losingteam,total_epa,result,all_total_success) %>% 
  pivot_wider(names_from = result, values_from = c(total_epa,all_total_success)) %>% 
  mutate(epa_lucky_win = ifelse(total_epa_loss>total_epa_win,1,0),
         success_lucky_win = ifelse(all_total_success_loss>all_total_success_win,1,0))

lucky_data %>% 
  group_by(winningteam) %>% 
  summarize(`Lucky Wins By EPA` = sum(epa_lucky_win), `Lucky Wins by Success Rate` = sum(success_lucky_win)) %>% 
  left_join(lucky_data %>% 
              group_by(losingteam) %>% 
              summarize(`Unlucky Losses by EPA` = sum(epa_lucky_win), `Unlucky Losses by Success Rate` = sum(success_lucky_win)),
            by = c("winningteam" = "losingteam")) %>% 
  mutate(across(where(is.numeric), ~if_else(is.na(.), 0, .))) %>% 
  mutate(`Net Lucky Wins By Success Rate` = `Lucky Wins by Success Rate`- `Unlucky Losses by Success Rate`) %>% 
  mutate(`Net Lucky Wins By EPA` = `Lucky Wins By EPA` - `Unlucky Losses by EPA`) %>% 
  rename(" " = "winningteam") %>% 
  select(" ", `Net Lucky Wins By EPA`,`Lucky Wins By EPA`, `Unlucky Losses by EPA`, `Net Lucky Wins By Success Rate`,`Lucky Wins by Success Rate`, `Unlucky Losses by Success Rate`) %>% 
  arrange(`Net Lucky Wins By Success Rate`) %>% 
  gt() %>% 
  gt_nfl_wordmarks(columns = c(" ")) %>% 
  cols_align(align = "center") %>% 
  # gt_hulk_col_numeric(columns = c("Net Lucky Wins By Success Rate", "Net Lucky Wins By EPA")) %>% 
  gt_hulk_col_numeric(columns = -" ") %>% 
  gt_theme_538() %>% 
  tab_header(
    title = md("Which Teams Have Been Luckiest this Season?"),
    subtitle = md("Unlucky/Lucky Win = Losing Team Outperformed Winning Team in Metric")
  ) %>% 
  tab_footnote(footnote = md("@CapAnalytics7|nflfastr"))
gtsave(check, "my_table.png")


#Purdy Game vs Bucs----
test <- nfl99 %>% 
  filter(name == "B.Purdy") %>% 
  group_by(game_id,id) %>% 
  summarize(`EPA/Play` = mean(epa,na.rm = T), success_rate = mean(success,na.rm = T)) %>% 
  arrange(-`EPA/Play`)

#Playground----  
summarize(total_epa_off = sum(epa,na.rm = T), off_plays = n(), homescore = max(home_score), hometeam = max(home_team),
          awayscore = max(away_score), awayteam = max(away_team),off_total_success = sum(success,na.rm = T),
          off_success_rate = mean(success, na.rm = T), off_epa_play = mean(epa,na.rm = T), off_fumbles = sum(fumble,na.rm = T), 
          off_fumbles_lost = sum(fumble_lost,na.rm =T), off_int = sum(interception,na.rm = T), off_int_worth = sum(is_interception_worthy,na.rm = T),
          off_expected_fumbles = recovery_rate * off_fumbles, off_expected_ints = interception_rate*off_int_worth,
          off_xp_made = sum(extra_result,na.rm = T), off_xp_attempts = sum(extra_point_attempt,na.rm = T), off_exp_xp = off_xp_attempts * pa_rate,
          off_2p_att = sum(two_point_attempt), off_2p_conv = sum(two_point_conv_result,na.rm = T), off_exp_2p = off_2p_conv*two_point_rate,
          off_fg_made = sum(field_goal_result,na.rm = T), off_fg_att = sum(field_goal_attempt,na.rm = T), off_fg_exp = sum(fg_prob[field_goal_attempt == 1],na.rm = T)) %>%
  
  
  
  lucky_data <- pbp22 %>% 
  group_by(game_id,posteam_type) %>% 
  summarize(total_epa_off = sum(epa,na.rm = T), off_plays = n(), homescore = max(home_score), hometeam = max(home_team),
            awayscore = max(away_score), awayteam = max(away_team),off_total_success = sum(success,na.rm = T)) %>% 
  mutate(winningteam = ifelse(homescore == awayscore, NA,ifelse(homescore>awayscore,hometeam,awayteam)),
         losingteam = ifelse(homescore == awayscore, NA,ifelse(homescore>awayscore,awayteam,hometeam))) %>% 
  filter(!is.na(posteam)) %>% 
  left_join(pbp22 %>% 
              group_by(game_id,defteam) %>% 
              summarize(total_epa_def = sum(epa,na.rm = T)*-1, def_plays = n(), total_success = sum(success,na.rm = T)) %>%
              mutate(def_total_success = def_plays-total_success) %>% 
              filter(!is.na(defteam)), by = c("game_id", "posteam" = "defteam")) %>% 
  mutate(total_epa = total_epa_off+total_epa_def,
         all_total_success = def_total_success+off_total_success) %>% 
  mutate(result = ifelse(is.na(winningteam), NA, ifelse(winningteam == posteam, "win","loss")))

avg_epa_int <- mean(nfl99$epa[nfl99$interception==1],na.rm = T)
avg_epa_fumble <- mean(nfl99$epa[nfl99$fumble_lost==1],na.rm = T)

tests <- pbp_rp %>% 
  filter(!is.na(air_yards)) %>% 
  mutate(unlucky_epa = ifelse(is_interception_worthy == FALSE & interception == 1,xyac_epa+air_epa-epa,
                              ifelse(is_interception_worthy == TRUE & interception == 0,avg_epa_int - epa,
                                     ifelse(fumble_lost == 0 & fumble == 1,avg_epa_fumble-epa,NA)))) %>% 
  select(desc,fumble,fumble_lost, epa, yac_epa,air_epa, xyac_epa,unlucky_epa)
table(pbp$extra_point_result)

median(pbp$epa[pbp$extra_point_result != "good"],na.rm = T)


#Expected Wins----
nfl99all <- load_pbp(1999:2024)
nfl99 <- nfl99all %>%
  filter(pass == 1 | rush == 1) %>%
  mutate(explosive = ifelse((yards_gained>20 & pass_attempt == 1) | (yards_gained >12 & (qb_scramble == 1 | rush == 1)),1,0),
         negative = ifelse(yards_gained < 0, 1,0))
ftn_data <- nflreadr::load_ftn_charting(2022:2024) %>%
  select(-week, -season)
participation <- load_participation(2022:2024) %>% 
  select(-old_game_id)
pbp22 <- load_pbp(2022:2024)
pbp22 <- pbp22 %>%
  left_join(ftn_data, by = c("game_id" = "nflverse_game_id",
                             "play_id" = "nflverse_play_id")) %>% 
  left_join(participation,by = c("game_id" = "nflverse_game_id",
                                 "play_id" = "play_id"))

library(tidyr)
recovery_rate <- mean(nfl99all$fumble_lost[nfl99all$fumble == 1],na.rm = T)
interception_rate <- mean(pbp22$interception[pbp22$is_interception_worthy],na.rm = T)
pa_rate <- nfl99all %>% mutate(extra_result = ifelse(extra_point_result == "good",1,0)) %>% 
  filter(season > 2014) %>%
  summarize(avg_conv = mean(extra_result[extra_point_attempt == 1],na.rm = T)) %>% 
  pull(avg_conv)
two_point_rate <- nfl99all %>% 
  mutate(two_point_conv_result = ifelse(two_point_conv_result == "success",1,0)) %>% 
  summarize(avg_conv = mean(two_point_conv_result[two_point_attempt == 1],na.rm = T)) %>% 
  pull(avg_conv)
lucky_data <- pbp22 %>% #add extra point totals, simulation vs model
  mutate(extra_result = ifelse(extra_point_result == "good",1,0)) %>% 
  mutate(two_point_conv_result = ifelse(two_point_conv_result == "success",1,0)) %>%
  mutate(field_goal_result = ifelse(field_goal_result == "made",1,0)) %>% 
  group_by(game_id,posteam) %>% 
  summarize(total_epa_off = sum(epa,na.rm = T), off_plays = n(), homescore = max(home_score), hometeam = max(home_team),
            awayscore = max(away_score), awayteam = max(away_team),off_total_success = sum(success,na.rm = T),
            off_success_rate = mean(success, na.rm = T), off_epa_play = mean(epa,na.rm = T),
            off_expected_ints_oe = sum(interception,na.rm = T)-(interception_rate* sum(is_interception_worthy,na.rm = T)),
            off_expected_fumbles_oe = sum(fumble_lost,na.rm =T)-(recovery_rate * sum(fumble,na.rm = T)),
            off_xp_oe = sum(extra_result,na.rm = T) - sum(extra_point_attempt,na.rm = T) * pa_rate,
            off_2p_oe = sum(two_point_conv_result,na.rm = T) - sum(two_point_attempt,na.rm = T)*two_point_rate,
            off_fg_oe = sum(field_goal_result,na.rm = T) -  sum(fg_prob[field_goal_attempt == 1],na.rm = T)) %>%
  filter(!is.na(posteam)) %>% 
  mutate(location = ifelse(hometeam==posteam, "hometeam", "awayteam"), home_team_win = ifelse(homescore>awayscore,1,0)) %>% 
  # mutate(winningteam = ifelse(homescore == awayscore, NA,ifelse(homescore>awayscore,hometeam,awayteam)),
  #        losingteam = ifelse(homescore == awayscore, NA,ifelse(homescore>awayscore,awayteam,hometeam))) %>% 
  left_join(pbp22 %>% 
              mutate(extra_result = ifelse(extra_point_result == "good",1,0)) %>% 
              mutate(two_point_conv_result = ifelse(two_point_conv_result == "success",1,0)) %>%
              mutate(field_goal_result = ifelse(field_goal_result == "made",1,0)) %>% 
              group_by(game_id,defteam) %>% 
              summarize(total_epa_def = sum(epa,na.rm = T)*-1, def_plays = n(), off_def_success = sum(success,na.rm = T), 
                        def_success_rate_allowed = mean(success,na.rm = T), def_epa_play = mean(epa,na.rm = T),
                        def_expected_ints_oe = sum(interception,na.rm = T)-(interception_rate* sum(is_interception_worthy,na.rm = T)),
                        def_expected_fumbles_oe = sum(fumble_lost,na.rm =T)-(recovery_rate * sum(fumble,na.rm = T)),
                        def_xp_oe = sum(extra_result,na.rm = T) - sum(extra_point_attempt,na.rm = T) * pa_rate,
                        def_2p_oe = sum(two_point_conv_result,na.rm = T) - sum(two_point_attempt,na.rm = T)*two_point_rate,
                        def_fg_oe = sum(field_goal_result,na.rm = T) -  sum(fg_prob[field_goal_attempt == 1],na.rm = T)) %>%
              mutate(def_total_success = def_plays-off_def_success) %>% 
              filter(!is.na(defteam)), by = c("game_id", "posteam" = "defteam")) %>% 
  mutate(total_epa = total_epa_off+total_epa_def,
         all_total_success = def_total_success+off_total_success, net_success_rate = all_total_success/(off_plays+def_plays),
         net_epa_play = total_epa/(off_plays+def_plays)) %>%
  ungroup() %>% 
  filter(location == "hometeam") %>% 
  select(-homescore,-awayscore, -all_total_success, -def_total_success,-posteam,-hometeam, -awayteam,-off_total_success,
         -off_def_success, -total_epa_off, -total_epa_def, -location)
  # mutate(result = ifelse(is.na(winningteam), NA, ifelse(winningteam == posteam, "win","loss")))
  # select(game_id, winningteam,losingteam,total_epa,result,all_total_success) %>% 
  # pivot_wider(names_from = location, values_from = c("off_plays", "off_success_rate", "off_epa_play", 
  #                                                    "def_success_rate_allowed", "def_epa_play", 
  #                                                    "total_epa", "net_success_rate", "net_epa_play"))

# Load necessary libraries
library(caret)
library(glmnet)
library(dplyr)

# Prepare the data: ensure no NA values and set predictors and response
lucky_data <- lucky_data %>%
  filter(!is.na(home_team_win)) # Ensure no missing values in the response variable

# Define predictor variables and the response
predictors <- lucky_data %>% 
  select(
  total_epa, net_success_rate, net_epa_play,
  off_success_rate,
                                    off_epa_play,
  off_expected_ints_oe, off_expected_fumbles_oe, off_xp_oe,
                                    off_2p_oe, off_fg_oe, def_success_rate_allowed, 
  def_epa_play,
                                    def_expected_ints_oe, def_expected_fumbles_oe, def_xp_oe, def_2p_oe, def_fg_oe)
response <- lucky_data$home_team_win

# Convert predictors to matrix format as required by glmnet
X <- as.matrix(predictors)
y <- factor(response, levels = c(0, 1), labels = c("loss", "win"))
# Set up cross-validation with the caret package
set.seed(123)  # for reproducibility
train_control <- trainControl(
  method = "cv",         # Use cross-validation
  number = 10,           # 10-fold cross-validation
  summaryFunction = twoClassSummary,  # For AUC metric
  classProbs = TRUE,     # Needed for AUC
  savePredictions = "final"
)

# Set up a tuning grid for regularization parameters
tune_grid <- expand.grid(
  alpha = c(0, 0.5, 1),   # Elastic net mixing parameter: 0 = Ridge, 1 = Lasso
  lambda = 10^seq(-4, 1, length = 10)  # Regularization parameter
)

# Fit the logistic regression model with cross-validation and parameter tuning
model <- train(
  X, y,
  method = "glmnet",
  trControl = train_control,
  tuneGrid = tune_grid,
  metric = "ROC"  # Use Area Under the ROC Curve to select best model
)

# Print best tuning parameters and model summary
print(model$bestTune)
print(model)

# Make predictions on the training data
lucky_data$predicted_home_win_prob <- predict(model, X, type = "prob")[, "win"]  # Probabilities for "1"
lucky_data$predicted_home_win <- ifelse(lucky_data$predicted_home_win_prob > 0.5, "win", "loss")

 #Evaluate the model’s performance
 conf_matrix <- confusionMatrix(factor(lucky_data$predicted_home_win), factor(y))
 print(conf_matrix)

 # Extract and print overall model metrics
 cat("Accuracy:", conf_matrix$overall["Accuracy"], "\n")
cat("AUC:", max(model$results$ROC), "\n")

variable_importance <- varImp(model, scale = FALSE)

# Print the variable importance
print(variable_importance)

#Lucky Data Approach 2----
lucky_data_all <- pbp22 %>% #look to incorporate epa with luck elements?
  mutate(extra_result = ifelse(extra_point_result == "good",1,0)) %>% 
  mutate(two_point_conv_result = ifelse(two_point_conv_result == "success",1,0)) %>%
  mutate(field_goal_result = ifelse(field_goal_result == "made",1,0)) %>% 
  group_by(game_id,posteam) %>% 
  summarize(total_epa_off = sum(epa,na.rm = T), off_plays = n(), homescore = max(home_score), hometeam = max(home_team),
            awayscore = max(away_score), awayteam = max(away_team),off_total_success = sum(success,na.rm = T),
            off_success_rate = mean(success, na.rm = T), off_epa_play = mean(epa,na.rm = T),
            off_expected_ints_oe = sum(interception,na.rm = T)-(interception_rate* sum(is_interception_worthy,na.rm = T)),
            off_expected_fumbles_oe = sum(fumble_lost,na.rm =T)-(recovery_rate * sum(fumble,na.rm = T)),
            off_xp_oe = sum(extra_result,na.rm = T) - sum(extra_point_attempt,na.rm = T) * pa_rate,
            off_2p_oe = sum(two_point_conv_result,na.rm = T) - sum(two_point_attempt,na.rm = T)*two_point_rate,
            off_fg_oe = sum(field_goal_result,na.rm = T) -  sum(fg_prob[field_goal_attempt == 1],na.rm = T)) %>%
  filter(!is.na(posteam)) %>% 
  mutate(location = ifelse(hometeam==posteam, "hometeam", "awayteam")) %>% 
  # mutate(winningteam = ifelse(homescore == awayscore, NA,ifelse(homescore>awayscore,hometeam,awayteam)),
  #        losingteam = ifelse(homescore == awayscore, NA,ifelse(homescore>awayscore,awayteam,hometeam))) %>% 
  left_join(pbp22 %>% 
              mutate(extra_result = ifelse(extra_point_result == "good",1,0)) %>% 
              mutate(two_point_conv_result = ifelse(two_point_conv_result == "success",1,0)) %>%
              mutate(field_goal_result = ifelse(field_goal_result == "made",1,0)) %>% 
              group_by(game_id,defteam) %>% 
              summarize(total_epa_def = sum(epa,na.rm = T)*-1, def_plays = n(), off_def_success = sum(success,na.rm = T), 
                        def_success_rate_allowed = mean(success,na.rm = T), def_epa_play = mean(epa,na.rm = T),
                        def_expected_ints_oe = sum(interception,na.rm = T)-(interception_rate* sum(is_interception_worthy,na.rm = T)),
                        def_expected_fumbles_oe = sum(fumble_lost,na.rm =T)-(recovery_rate * sum(fumble,na.rm = T)),
                        def_xp_oe = sum(extra_result,na.rm = T) - sum(extra_point_attempt,na.rm = T) * pa_rate,
                        def_2p_oe = sum(two_point_conv_result,na.rm = T) - sum(two_point_attempt,na.rm = T)*two_point_rate,
                        def_fg_oe = sum(field_goal_result,na.rm = T) -  sum(fg_prob[field_goal_attempt == 1],na.rm = T),
                        season = max(season), week = max(week)) %>%
              mutate(def_total_success = def_plays-off_def_success) %>% 
              filter(!is.na(defteam)), by = c("game_id", "posteam" = "defteam")) %>% 
  mutate(total_epa = total_epa_off+total_epa_def,
         all_total_success = def_total_success+off_total_success, net_success_rate = all_total_success/(off_plays+def_plays),
         net_epa_play = total_epa/(off_plays+def_plays), home_team_deserved_win = ifelse(net_epa_play>=0,1,0)) %>%
  ungroup() %>% 
  filter(location == "hometeam") %>% 
  select(-homescore,-awayscore, -all_total_success, -def_total_success,-posteam,-hometeam, -awayteam,-off_total_success,
         -off_def_success, -total_epa_off, -total_epa_def, -location)
# mutate(result = ifelse(is.na(winningteam), NA, ifelse(winningteam == posteam, "win","loss")))
# select(game_id, winningteam,losingteam,total_epa,result,all_total_success) %>% 
# pivot_wider(names_from = location, values_from = c("off_plays", "off_success_rate", "off_epa_play", 
#                                                    "def_success_rate_allowed", "def_epa_play", 
#                                                    "total_epa", "net_success_rate", "net_epa_play"))

# Load necessary libraries
library(caret)
library(glmnet)
library(dplyr)

# Prepare the data: ensure no NA values and set predictors and response
lucky_data <- lucky_data_all %>%
  filter(season <2024) %>% 
  filter(!is.na(home_team_deserved_win)) # Ensure no missing values in the response variable

# Define predictor variables and the response
predictors <- lucky_data %>% 
  select(-total_epa, net_success_rate, -net_epa_play,
    off_success_rate,
    -off_epa_play,-week,-season,
    off_expected_ints_oe, off_expected_fumbles_oe, off_xp_oe,
    off_2p_oe, off_fg_oe, def_success_rate_allowed, 
    -def_epa_play,
    def_expected_ints_oe, def_expected_fumbles_oe, def_xp_oe, def_2p_oe, def_fg_oe, -home_team_deserved_win, -game_id,-def_plays,-off_plays)
response <- lucky_data$home_team_deserved_win

# Convert predictors to matrix format as required by glmnet
X <- as.matrix(predictors)
y <- factor(response, levels = c(0, 1), labels = c("loss", "win"))
# Set up cross-validation with the caret package
set.seed(123)  # for reproducibility
train_control <- trainControl(
  method = "cv",         # Use cross-validation
  number = 10,           # 10-fold cross-validation
  summaryFunction = twoClassSummary,  # For AUC metric
  classProbs = TRUE,     # Needed for AUC
  savePredictions = "final"
)

# Set up a tuning grid for regularization parameters
tune_grid <- expand.grid(
  alpha = c(0, 0.5, 1),   # Elastic net mixing parameter: 0 = Ridge, 1 = Lasso
  lambda = 10^seq(-4, 1, length = 10)  # Regularization parameter
)

# Fit the logistic regression model with cross-validation and parameter tuning
model <- train(
  X, y,
  method = "glmnet",
  trControl = train_control,
  tuneGrid = tune_grid,
  metric = "ROC"  # Use Area Under the ROC Curve to select best model
)

# Print best tuning parameters and model summary
print(model$bestTune)
print(model)
lucky_data_24 <- lucky_data_all %>% 
  filter(season == 2024)

# Make predictions on the training data
lucky_data_24$predicted_home_win_prob<- predict(model,lucky_data_24 %>% select(-home_team_deserved_win), type = "prob")[, "win"]  # Probabilities for "1"
lucky_data_24$predicted_home_win <- ifelse(lucky_data_24$predicted_home_win_prob > 0.5, "win", "loss")

lucky_data$predicted_home_win_prob <- predict(model, X, type = "prob")[, "win"]  # Probabilities for "1"
lucky_data$predicted_home_win <- ifelse(lucky_data$predicted_home_win_prob > 0.5, "win", "loss")
#Evaluate the model’s performance
conf_matrix <- confusionMatrix(factor(lucky_data$predicted_home_win), factor(y))
print(conf_matrix)

# Extract and print overall model metrics
cat("Accuracy:", conf_matrix$overall["Accuracy"], "\n")
cat("AUC:", max(model$results$ROC), "\n")

variable_importance <- varImp(model, scale = FALSE)

# Print the variable importance
print(variable_importance)


test <- lucky_data_24 %>% mutate(away_win_prob = 1 - predicted_home_win_prob) %>% 
  left_join(pbp22 %>%
              filter(season == 2024) %>% 
              group_by(game_id) %>% 
              summarize(homescore = max(home_score), hometeam = max(home_team),
                        awayscore = max(away_score), awayteam = max(away_team),
                        season = max(season), week = max(week)),
             by = c("game_id")) %>% 
  select(hometeam,awayteam, predicted_home_win_prob, away_win_prob) %>% 
  pivot_longer(
    cols = c(hometeam, awayteam),
    names_to = "team_type",
    values_to = "team"
  ) %>% 
    mutate(
      win_prob = ifelse(team_type == "hometeam", predicted_home_win_prob, away_win_prob)
    ) %>% 
  group_by(team) %>% 
  summarize(total_win = sum(win_prob))
  
#Approach #3----
lucky_data_all_3 <- pbp22 %>% #look to incorporate epa with luck elements?
  mutate(extra_result = ifelse(extra_point_result == "good",1,0)) %>% 
  mutate(two_point_conv_result = ifelse(two_point_conv_result == "success",1,0)) %>%
  mutate(field_goal_result = ifelse(field_goal_result == "made",1,0)) %>%
  mutate(epa_noiseless = ifelse(fumble == 1 | (interception == 1 & is_interception_worthy == 0), NA, epa)) %>% 
  group_by(game_id,posteam) %>% 
  summarize(total_epa_off = sum(epa,na.rm = T), off_plays = n(), homescore = max(home_score), hometeam = max(home_team),
            awayscore = max(away_score), awayteam = max(away_team),off_total_success = sum(success,na.rm = T),
            off_success_rate = mean(success, na.rm = T), off_epa_play = mean(epa_noiseless,na.rm = T),
            off_expected_ints_oe = sum(interception,na.rm = T)-(interception_rate* sum(is_interception_worthy,na.rm = T)),
            off_expected_fumbles_oe = sum(fumble_lost,na.rm =T)-(recovery_rate * sum(fumble,na.rm = T)),
            off_xp_oe = sum(extra_result,na.rm = T) - sum(extra_point_attempt,na.rm = T) * pa_rate,
            off_2p_oe = sum(two_point_conv_result,na.rm = T) - sum(two_point_attempt,na.rm = T)*two_point_rate,
            off_fg_oe = sum(field_goal_result,na.rm = T) -  sum(fg_prob[field_goal_attempt == 1],na.rm = T)) %>%
  filter(!is.na(posteam)) %>% 
  mutate(location = ifelse(hometeam==posteam, "hometeam", "awayteam")) %>% 
  # mutate(winningteam = ifelse(homescore == awayscore, NA,ifelse(homescore>awayscore,hometeam,awayteam)),
  #        losingteam = ifelse(homescore == awayscore, NA,ifelse(homescore>awayscore,awayteam,hometeam))) %>% 
  left_join(pbp22 %>% 
              mutate(extra_result = ifelse(extra_point_result == "good",1,0)) %>% 
              mutate(two_point_conv_result = ifelse(two_point_conv_result == "success",1,0)) %>%
              mutate(field_goal_result = ifelse(field_goal_result == "made",1,0)) %>% 
              mutate(epa_noiseless = ifelse(fumble == 1 | (interception == 1 & is_interception_worthy == 0), NA, epa)) %>% 
              group_by(game_id,defteam) %>% 
              summarize(total_epa_def = sum(epa,na.rm = T)*-1, def_plays = n(), off_def_success = sum(success,na.rm = T), 
                        def_success_rate_allowed = mean(success,na.rm = T), def_epa_play = mean(epa_noiseless,na.rm = T),
                        def_expected_ints_oe = sum(interception,na.rm = T)-(interception_rate* sum(is_interception_worthy,na.rm = T)),
                        def_expected_fumbles_oe = sum(fumble_lost,na.rm =T)-(recovery_rate * sum(fumble,na.rm = T)),
                        def_xp_oe = sum(extra_result,na.rm = T) - sum(extra_point_attempt,na.rm = T) * pa_rate,
                        def_2p_oe = sum(two_point_conv_result,na.rm = T) - sum(two_point_attempt,na.rm = T)*two_point_rate,
                        def_fg_oe = sum(field_goal_result,na.rm = T) -  sum(fg_prob[field_goal_attempt == 1],na.rm = T),
                        season = max(season), week = max(week)) %>%
              mutate(def_total_success = def_plays-off_def_success) %>% 
              filter(!is.na(defteam)), by = c("game_id", "posteam" = "defteam")) %>% 
  mutate(total_epa = total_epa_off+total_epa_def,
         all_total_success = def_total_success+off_total_success, net_success_rate = all_total_success/(off_plays+def_plays),
         net_epa_play = total_epa/(off_plays+def_plays), home_team_deserved_win = ifelse(net_epa_play>=0,1,0)) %>%
  ungroup() %>% 
  filter(location == "hometeam") %>% 
  select(-homescore,-awayscore, -all_total_success, -def_total_success,-posteam,-hometeam, -awayteam,-off_total_success,
         -off_def_success, -total_epa_off, -total_epa_def, -location)
# mutate(result = ifelse(is.na(winningteam), NA, ifelse(winningteam == posteam, "win","loss")))
# select(game_id, winningteam,losingteam,total_epa,result,all_total_success) %>% 
# pivot_wider(names_from = location, values_from = c("off_plays", "off_success_rate", "off_epa_play", 
#                                                    "def_success_rate_allowed", "def_epa_play", 
#                                                    "total_epa", "net_success_rate", "net_epa_play"))

# Load necessary libraries
library(caret)
library(glmnet)
library(dplyr)

# Prepare the data: ensure no NA values and set predictors and response
lucky_data <- lucky_data_all_3 %>%
  filter(season <2024) %>% 
  filter(!is.na(home_team_deserved_win)) # Ensure no missing values in the response variable

# Define predictor variables and the response
predictors <- lucky_data %>% 
  select(-total_epa, net_success_rate, -net_epa_play,
         off_success_rate,
         off_epa_play,-week,-season,
         off_expected_ints_oe, off_expected_fumbles_oe, off_xp_oe,
         off_2p_oe, off_fg_oe, def_success_rate_allowed, 
         def_epa_play,
         def_expected_ints_oe, def_expected_fumbles_oe, def_xp_oe, def_2p_oe, def_fg_oe, -home_team_deserved_win, -game_id,-def_plays,-off_plays)
response <- lucky_data$home_team_deserved_win

# Convert predictors to matrix format as required by glmnet
X <- as.matrix(predictors)
y <- factor(response, levels = c(0, 1), labels = c("loss", "win"))
# Set up cross-validation with the caret package
set.seed(123)  # for reproducibility
train_control <- trainControl(
  method = "cv",         # Use cross-validation
  number = 10,           # 10-fold cross-validation
  summaryFunction = twoClassSummary,  # For AUC metric
  classProbs = TRUE,     # Needed for AUC
  savePredictions = "final"
)

# Set up a tuning grid for regularization parameters
tune_grid <- expand.grid(
  alpha = c(0, 0.5, 1),   # Elastic net mixing parameter: 0 = Ridge, 1 = Lasso
  lambda = 10^seq(-4, 1, length = 10)  # Regularization parameter
)

# Fit the logistic regression model with cross-validation and parameter tuning
model <- train(
  X, y,
  method = "glmnet",
  trControl = train_control,
  tuneGrid = tune_grid,
  metric = "ROC"  # Use Area Under the ROC Curve to select best model
)

# Print best tuning parameters and model summary
print(model$bestTune)
print(model)
lucky_data_24 <- lucky_data_all %>% 
  filter(season == 2024)

# Make predictions on the training data
lucky_data_24$predicted_home_win_prob<- predict(model,lucky_data_24 %>% select(-home_team_deserved_win), type = "prob")[, "win"]  # Probabilities for "1"
lucky_data_24$predicted_home_win <- ifelse(lucky_data_24$predicted_home_win_prob > 0.5, "win", "loss")

lucky_data$predicted_home_win_prob <- predict(model, X, type = "prob")[, "win"]  # Probabilities for "1"
lucky_data$predicted_home_win <- ifelse(lucky_data$predicted_home_win_prob > 0.5, "win", "loss")
#Evaluate the model’s performance
conf_matrix <- confusionMatrix(factor(lucky_data$predicted_home_win), factor(y))
print(conf_matrix)

# Extract and print overall model metrics
cat("Accuracy:", conf_matrix$overall["Accuracy"], "\n")
cat("AUC:", max(model$results$ROC), "\n")

variable_importance <- varImp(model, scale = FALSE)

# Print the variable importance
print(variable_importance)


x_wins <- lucky_data_24 %>% mutate(away_win_prob = 1 - predicted_home_win_prob) %>% 
  left_join(pbp22 %>%
              filter(season == 2024) %>% 
              group_by(game_id) %>% 
              summarize(homescore = max(home_score), hometeam = max(home_team),
                        awayscore = max(away_score), awayteam = max(away_team),
                        season = max(season), week = max(week)),
            by = c("game_id")) %>% 
  select(hometeam,awayteam, predicted_home_win_prob, away_win_prob,homescore,awayscore) %>% 
  pivot_longer(
    cols = c(hometeam, awayteam),
    names_to = "team_type",
    values_to = "team"
  ) %>% 
  mutate(
    win_prob = ifelse(team_type == "hometeam", predicted_home_win_prob, away_win_prob),
    win = ifelse((team_type == "awayteam" & awayscore>homescore )|(team_type == "hometeam" & awayscore<homescore ) ,1,0)
  ) %>% 
  group_by(team) %>% 
  summarize(exp_total_win = round(sum(win_prob),2),total_wins = sum(win)) %>% 
  mutate(wins_oe = total_wins-exp_total_win)
xw_tab <- x_wins %>% 
  arrange(-exp_total_win) %>% 
  gt() %>% 
  cols_align(align = "center") %>% 
  cols_label(team = "",
             exp_total_win = "Expected Wins",
             total_wins = "Total Wins",
             wins_oe = "Wins Over Expectation") %>%
  gtExtras::gt_theme_538() %>% 
  gt_nfl_wordmarks(columns = "team") %>% 
  gt_hulk_col_numeric(columns = wins_oe) %>% 
  tab_header(
    title = md("How Lucky Are We?"),
    subtitle = md("Expected Wins vs Actual Wins; Green Represents Overperformance, Purple Represents Underperformance")
  ) %>% 
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_title()
  )
gtsave(xw_tab, "XW.png")
