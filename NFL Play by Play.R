library(ggimage)
library(gt)
library(nflfastR)
library(tidyverse)

pbp <- load_pbp(2023)

pbp_rp <- pbp %>%
  filter(pass == 1 | rush == 1) %>%
  filter(!is.na(epa))


#Offensive Efficiency Plot
pass_efficiency <- pbp_rp %>%
  filter(season == 2023, pass ==1) %>%
  group_by(posteam) %>%
  summarize(passes = n(), pass_epa = mean(epa))

rush_efficiency <- pbp_rp %>%
  filter(season == 2023, rush ==1) %>%
  group_by(posteam) %>%
  summarize(rushes = n(), rush_epa= mean(epa))

success_rate <- pbp_rp %>% 
  group_by(posteam) %>% 
  summarize(total_epa = mean(epa), success = mean(success))

success_rate <- success_rate %>% 
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr"))
success_rate %>%
  ggplot(aes(x = success, y = total_epa))+
  geom_image(aes(image = team_logo_espn), size = 0.05, asp = 16/9)+
  theme_bw()+
  labs(y = "EPA Per Play",
       x = "Success Rate", title = "NFL Offensive Efficency Landscape Since Week 13",
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")
ggsave("OffensiveLandscape.png", width = 14, height =10, dpi = "retina")

total_eff <- pass_efficiency %>% 
  left_join(rush_efficiency, by = "posteam")

total_eff <- total_eff %>%
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr"))

total_eff %>%
  ggplot(aes(x = pass_epa, y = rush_epa))+
  geom_hline(yintercept = mean(total_eff$rush_epa), linetype = "dashed")+
  geom_vline(xintercept = mean(total_eff$pass_epa), linetype = "dashed")+
  geom_image(aes(image = team_logo_espn), size = 0.05, asp = 16/9)+
  theme_bw()+
  labs(x = "EPA Per Pass",
       y = "EPA Per Rush", title = "NFL Offensive Pass vs Rush Efficiency", subtitle = "Dotted Lines Represent League Average",
       caption = "Callan Capitolo | @CapAnalytics7")
ggsave("Week9OffensiveEfficiency.png", width = 14, height =10, dpi = "retina")

#Red zone efficiency
pbp_rp_red <- pbp %>%
  filter(pass == 1 | rush == 1) %>%
  filter(!is.na(epa))%>%
  filter(yardline_100 <=20)

pass_efficiency_red <- pbp_rp_red %>%
  filter(season == 2023, pass ==1) %>%
  group_by(posteam) %>%
  summarize(passes = n(), pass_epa = mean(epa))

rush_efficiency_red <- pbp_rp_red %>%
  filter(season == 2023, rush ==1) %>%
  group_by(posteam) %>%
  summarize(rushes = n(), rush_epa= mean(epa))

total_eff_red <- pass_efficiency_red %>% 
  left_join(rush_efficiency_red, by = "posteam")

total_eff_red <- total_eff_red %>%
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr"))

total_eff_red <- total_eff_red %>%
  mutate(plays = passes+rushes)

total_eff_red %>%
  ggplot(aes(x = pass_epa, y = rush_epa))+
  geom_text(aes(label = posteam, size = plays, color = team_color))+
  scale_color_identity(aesthetics = "color")+
  labs(x = "EPA Per Pass",
       y = "EPA Per Rush", title = "NFL Offensive Red Zone Efficiency Following Week 17", subtitle = "Size of Text Represents Number of Red Zone Plays",
       caption = "Callan Capitolo | @CapAnalytics7", size = "Red Zone Plays")+
  theme_bw()+
  #scale_color_manual(values = total_eff_red$team_color) +
  geom_vline(xintercept = 0, color = "black", size = 0.5) +  # Vertical line at x=0
  geom_hline(yintercept = 0, color = "black", size = 0.5)# Horizontal line at y=0
ggsave("Week6RedZoneEfficiency.png", width = 14, height =10, dpi = "retina")


#Offense/ Defense Efficiency

total_offensive_efficiency <- pbp_rp %>%
  filter(season == 2023) %>%
  group_by(posteam) %>%
  summarize(offensive_epa = mean(epa))

total_defensive_efficiency <- pbp_rp %>%
  filter(season == 2023) %>%
  group_by(defteam) %>%
  summarize(defensive_epa = mean(epa))

total_efficiency_both <- total_offensive_efficiency %>%
  left_join(total_defensive_efficiency, by = c("posteam" = "defteam"))

total_efficiency_both <- total_efficiency_both %>%
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr"))

total_efficiency_both %>%
  ggplot(aes(x = defensive_epa, y = offensive_epa))+
  geom_image(aes(image = team_logo_espn), size = 0.05, asp = 16/9)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  geom_image(aes(image = team_logo_espn), size = 0.05, asp = 16/9)+
  theme_minimal()+
  scale_x_reverse()+
  labs(x = "Defensive Efficiency",
       y = "Offensive Efficiency", title = "NFL Offensive & Defensive Efficiency Following Week 17", subtitle = "Efficiency Represented by EPA/Play",
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")+
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())+
  geom_text(aes(x = 0.05, y = 0.25, label = "Bad Defense, Good Offense"), size = 3) +
  geom_text(aes(x = -0.15, y = 0.25, label = "Contenders"), size = 3) +
  geom_text(aes(x = 0.05, y = -0.25, label = "Bad Defense, Bad Offense"), size = 3) +
  geom_text(aes(x = -0.15, y = -0.25, label = "Bad Offense, Good Defense"), size = 3)
ggsave("CombinedEfficiency.png", width = 14, height =10, dpi = "retina")


##Previous week efficiency
prior_week <- pbp_rp %>% 
  filter(week %in% c(1:max(week)-1))

prior_offensive_efficiency <- prior_week %>%
  filter(season == 2023) %>%
  group_by(posteam) %>%
  summarize(offensive_epa = mean(epa), week = max(week))

prior_defensive_efficiency <- prior_week %>%
  filter(season == 2023) %>%
  group_by(defteam) %>%
  summarize(defensive_epa = mean(epa))

prior_efficiency_both <- prior_offensive_efficiency %>%
  left_join(prior_defensive_efficiency, by = c("posteam" = "defteam"))

prior_efficiency_both <- prior_efficiency_both %>%
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr"))



# Comparison to Previous Week
prior_week = function(x){
  pass_efficiency_fun <- x %>%
    filter(season == 2023, pass ==1, week %in% c(1:(max(week)-1))) %>%
    group_by(posteam) %>%
    summarize(passes = n(), pass_epa = mean(epa))
  
  rush_efficiency_fun <- x %>%
    filter(season == 2023, rush ==1, week %in% c(1:(max(week)-1))) %>% 
    group_by(posteam) %>%
    summarize(rushes = n(), rush_epa= mean(epa))
  
  total_eff_fun <- pass_efficiency_fun %>% 
    left_join(rush_efficiency_fun, by = "posteam")
  
  total_eff_fun <- total_eff_fun %>%
    left_join(teams_colors_logos, by = c("posteam" = "team_abbr"))
  return(total_eff_fun)
}

current_week = function(x){
  pass_efficiency_fun <- x %>%
    filter(season == 2023, pass ==1, week == max(week)) %>%
    group_by(posteam) %>%
    summarize(passes = n(), pass_epa = mean(epa))
  
  rush_efficiency_fun <- x %>%
    filter(season == 2023, rush ==1, week == max(week)) %>% 
    group_by(posteam) %>%
    summarize(rushes = n(), rush_epa= mean(epa))
  
  total_eff_fun <- pass_efficiency_fun %>% 
    left_join(rush_efficiency_fun, by = "posteam")
  
  total_eff_fun <- total_eff_fun %>%
    left_join(teams_colors_logos, by = c("posteam" = "team_abbr"))
  
  return(total_eff_fun)
}

prior_week_data = prior_week(pbp_rp)
current_week_data = current_week(pbp_rp)



change_eff = current_week_data %>%
  mutate(change_pass_epa = (pass_epa - prior_week_data$pass_epa)/abs(prior_week_data$pass_epa), change_rush_epa = (rush_epa - prior_week_data$rush_epa)/abs(prior_week_data$rush_epa))


change_eff %>%
  ggplot(aes(x = change_pass_epa*100, y = change_rush_epa*100))+
  scale_x_continuous(labels = percent_format(scale = 1)) +
  scale_y_continuous(labels = percent_format(scale = 1)) +
  geom_image(aes(image = team_logo_espn), size = 0.05, asp = 16/9)+
  theme_bw()+
  labs(x = "% Change in EPA Per Pass",
       y = "% Change in EPA Per Rush", title = "Week 4 % Change in Offensive Efficiency Compared to Season Averages", 
       subtitle  = "Eagles Pass Offense and 49ers Rush Offensive Signficantly Improve in Week 4", caption = "By Callan Capitolo")+
  theme(plot.title = element_text(size = 12))
ggsave("Week4ChangeinEfficiency.png", width = 14, height =10, dpi = "retina")

