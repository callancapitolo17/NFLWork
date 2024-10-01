library(nflfastR)
library(ggplot2)
library(tidyverse)
library(ggimage)
library(ggthemes)
library(dplyr)
library(ggrepel)
library(nflreadr)
library(gt)

pbp <- load_pbp(2023)

pbp_rp <- pbp %>%
  filter(pass == 1 | rush == 1) %>%
  filter(!is.na(epa))


#Home vs Road NFL Efficiency
home <- pbp_rp %>% 
  filter(posteam == home_team) %>% 
  group_by(posteam) %>% 
  summarize(home_epa = mean(epa))

away <- pbp_rp %>% 
  filter(posteam == away_team) %>% 
  group_by(posteam) %>% 
  summarize(away_epa = mean(epa))

combined <- home %>% 
  left_join(away, "posteam") %>% 
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr"))

combined %>% 
  ggplot(aes(x = away_epa, y = home_epa))+   
  geom_hline(yintercept = mean(combined$away_epa), linetype = "dashed")+
  geom_vline(xintercept = mean(combined$home_epa), linetype = "dashed")+
  geom_image(aes(image = team_logo_espn), size = 0.05, asp = 16/9)+
  theme_bw()+
  labs(x = "Away EPA/Play",
       y = "Home EPA/Play", title = "Home vs Away Offensive Efficiency Following Week 10",
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR", subtitle = "Dashed Lines Represent League Average")
ggsave("HomevsAwayOffense.png", width = 14, height =10, dpi = "retina")


#Running Back Efficiency
rushing_player<-pbp_rp %>% 
  filter(rush == 1, !is.na(rusher_player_name)) %>% 
  group_by(rusher_player_id,rusher_player_name) %>% 
  summarize(rushes = n(), rush_epa = mean(epa), success_rate = mean(success)) %>% 
  filter(rushes>=30) %>% 
  mutate(rusher_player_name = ifelse(rusher_player_name == "T.Hill", "Taysom Hill", rusher_player_name))
    


#EPA by late vs early down by team
rushing_player %>%   
  ggplot(aes(x = success_rate, y = rush_epa))+
  geom_point()+
  labs(x = "Success Rate", y = "EPA/Rush", title = "Rushing Efficiency by Rusher Following Week 10", subtitle = "Minimum 30 Rushes",
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")+
  theme_bw()+
  # geom_text(
  #   data = subset(rushing_player, rush_epa >0.35 | rush_epa < -0.35),
  #   aes(label = rusher_player_name),
  #   nudge_x = 0.01, # adjust this value for proper positioning of labels
  #   nudge_y = 0.01  # adjust this value for proper positioning of labels
  # )
  geom_text_repel(
    data = subset(rushing_player, rush_epa >0.13 | rush_epa < -0.3 | success_rate > 0.5 | success_rate < 0.26),
    aes(label = rusher_player_name),
    box.padding = 0.05,  # adjust this value for padding around the labels
    point.padding = 0.01,  # adjust this value for padding around the points
    segment.color = "grey",
    segment.size = 0.2
  )
ggsave("RunningEfficiency.png", width = 14, height =10, dpi = "retina")

#Explosive Pass Play Rate vs aDot

explosive_pass <- pbp %>% 
  filter(pass == 1, !is.na(passer_player_name), !is.na(pass_length)) %>% 
  mutate(explosive = ifelse(yards_gained >= 20, 1,0)) %>% 
  group_by(passer_id,passer_player_name, posteam) %>% 
  summarize(passes = n(), adot = mean(air_yards,na.rm = T), explosive_rate = mean(explosive)) %>% 
  filter(passes>=50)

explosive_pass <- explosive_pass %>% 
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr"))

explosive_pass %>% 
  ggplot(aes(x = adot, y = explosive_rate))+
  geom_image(aes(image = team_logo_espn), size = 0.02, asp = 16/9)+
  geom_text_repel(
    aes(label = passer_player_name),
    box.padding = 0.05,  # adjust this value for padding around the labels
    point.padding = 0.1,  # adjust this value for padding around the points
    segment.color = "grey",
    segment.size = 0.2
  )+
  theme_bw()+
  labs(x = "Average Depth of Target",
       y = "Explosive Pass Rate", title = "Explosive Pass Rate vs Average Depth of Target Following Week 10",
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR", subtitle = "Minimum 50 Passes")
ggsave("ExplosivePassRate.png", width = 14, height =10, dpi = "retina")

#Lead vs Trailing
leading <- pbp_rp %>% 
  filter(posteam_score >= defteam_score) %>% 
  group_by(posteam) %>% 
  summarize(plays = n(), epa_leading_play = mean(epa))

trailing <- pbp_rp %>% 
  filter(posteam_score < defteam_score) %>% 
  group_by(posteam) %>% 
  summarize(plays = n(), epa_trailing_play = mean(epa))

trailvslead <-leading %>% 
  left_join(trailing, by = "posteam") %>% 
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr"))

trailvslead %>% 
  ggplot(aes(x = epa_trailing_play, y = epa_leading_play)) +
  geom_image(aes(image = team_logo_espn), size = 0.1, asp = 16/19)+
  theme_bw()+
  labs(x = "EPA Per Play When Trailing", y = "EPA Per Playing When Tied or Leading", title = "Offensive Efficiency When Leading/Tied vs Trailing",
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")
ggsave("LeadingvsTrailing.png", width = 14, height =10, dpi = "retina")

#Early down vs Late Down Efficiency
early_down <- pbp_rp %>% 
  filter(down <= 2) %>% 
  group_by(posteam) %>% 
  summarize(early_down_epa = mean(epa))

late_down <- pbp_rp %>% 
  filter(down > 2) %>% 
  group_by(posteam) %>% 
  summarize(late_down_epa = mean(epa))

earlyvslate <-early_down %>% 
  left_join(late_down, by = "posteam") %>% 
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr"))

earlyvslate %>% 
  ggplot(aes(x = early_down_epa, y = late_down_epa)) +
  geom_image(aes(image = team_logo_espn), size = 0.1, asp = 16/19)+
  theme_bw()+
  labs(x = "EPA Per Play on Early Downs (1st & 2nd down)", y = "EPA Per Play on Late Down (3rd & 4th down)", title = "Offensive Efficiency Late Down vs Early Down Entering Week 11",
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")
ggsave("EarlyvsLateOffensiveEfficiency.png", width = 14, height =10, dpi = "retina")
#1st Half Total Efficiency
first_half_off <- pbp_rp %>% 
  filter(qtr %in% c(1,2)) %>% 
  group_by(posteam)%>% 
  summarize(off_epa = mean(epa))

first_half_def <- pbp_rp %>% 
  filter(qtr %in% c(1,2)) %>% 
  group_by(defteam)%>% 
  summarize(def_epa = mean(epa))

total_first_half <- first_half_off %>% 
  left_join(first_half_def, by = c("posteam" = "defteam")) %>% 
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr"))

total_first_half <- total_first_half %>% 
  mutate(total_epa_play = ifelse(def_epa<0,abs(def_epa),def_epa) + off_epa)


total_first_half %>% 
  ggplot(aes(x = def_epa, y = off_epa)) +
  geom_image(aes(image = team_logo_espn), size = 0.1, asp = 16/19)+
  theme_bw()+
  scale_x_reverse()+
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())+
  labs(x = "Defensive EPA Per Play", y = "Offesnive EPA Per Play", title = "1st Half Offensive and Defensive Efficiency Entering Week 11",
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")
ggsave("First Half Performance", width = 14, height =10, dpi = "retina")

week11 <- load_schedules(2023) %>% 
  filter(week ==11)

matchups <- week11 %>% 
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

qb_gt <- epa_1h_matchups %>%
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
    title = md("1st Half Efficiency by Week 11 Matchups"),
    subtitle = md("Negative Defensive EPA is Good, Positive Offensive EPA is Good")
  )
gtsave(qb_gt, "FirstHalfEfficiency.png") 

home_wp <- pbp %>% 
  group_by(game_seconds_remaining,home_team) %>% 
  summarize(avg_home_wp = mean(home_wp), home_count = n())

away <- pbp %>% 
  group_by(game_seconds_remaining,away_team) %>% 
  summarize(avg_away_wp = mean(away_wp), away_count = n())

wp_total <- home_wp %>% 
  full_join(away, by = c("home_team" = "away_team","game_seconds_remaining" = "game_seconds_remaining")) %>%
  mutate(avg_home_wp = ifelse(is.na(avg_home_wp),0,avg_home_wp), avg_away_wp = ifelse(is.na(avg_away_wp),0,avg_away_wp), 
         away_count = ifelse(is.na(away_count),0,away_count), home_count = ifelse(is.na(home_count),0,home_count)) %>% 
  mutate(avg_total_wp = (home_count*avg_home_wp+away_count*avg_away_wp)/(home_count+away_count)) %>% 
  left_join(teams_colors_logos, by = c("home_team" = "team_abbr"))

wp_total %>% 
  ggplot(aes(x = game_seconds_remaining, y = avg_total_wp, group = home_team))+
  stat_smooth(se = FALSE, show.legend = FALSE, aes(color =team_color))+
  xlim(0,3600)+
  scale_x_reverse()+
  scale_color_identity()+
  geom_image(
    data = subset(wp_total, game_seconds_remaining == 0),
    aes(image = team_logo_espn), size = 0.02, asp = 16/9, nudge_y = 0) +
  # geom_label_repel(aes(label = home_team), box.padding = 0.5, point.padding = 0.5) +
  labs(x = "Time Remaining in Game (Seconds)",
       y = "Average Win Probability", title = "Average Win Probability vs Time Remaining")

