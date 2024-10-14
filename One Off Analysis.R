library(ggimage)
library(gt)
library(nflfastR)
library(tidyverse)

pbp <- load_pbp(2010:2023)

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


