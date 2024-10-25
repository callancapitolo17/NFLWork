library(ggimage)
library(gt)
library(nflfastR)
library(tidyverse)

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
