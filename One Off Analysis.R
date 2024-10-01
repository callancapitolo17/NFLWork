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
  





