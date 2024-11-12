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
year <- 2024

pfr_stats_pass <- nflreadr::load_pfr_advstats(seasons = 2023,
                                              stat_type = "pass",
                                              summary_level = "week")
ngs_data_passing <- nflreadr::load_nextgen_stats(seasons = 2023,
                                                 stat_type = "passing")

ftn_data <- nflreadr::load_ftn_charting(2022:2024) %>%
  select(-week, -season)
participation <- load_participation(2022:2024) %>% 
  select(-old_game_id)
pbp <- load_pbp(2024)
pbp <- pbp %>%
  left_join(ftn_data, by = c("game_id" = "nflverse_game_id",
                             "play_id" = "nflverse_play_id")) %>% 
  left_join(participation,by = c("game_id" = "nflverse_game_id",
                                 "play_id" = "play_id"))

pbp_rp <- pbp %>%
  filter(pass == 1 | rush == 1) %>%
  filter(qb_kneel == 0,qb_spike == 0) %>% 
  filter(!is.na(epa)) %>% 
  mutate(blitz = ifelse(n_pass_rushers> 4,1,0),
         lightbox = ifelse(n_defense_box<=6,1,0),
         heavybox = ifelse(n_defense_box>=8,1,0)) %>% 
  mutate(explosive = ifelse((yards_gained>20 & pass_attempt == 1) | (yards_gained >12 & (qb_scramble == 1 | rush == 1)),1,0),
        negative = ifelse(yards_gained < 0, 1,0))

# nfl99all <- load_pbp(1999:2024)
# nfl99 <- nfl99all %>% 
#   filter(pass == 1 | rush == 1) %>% 
#   mutate(explosive = ifelse((yards_gained>20 & pass_attempt == 1) | (yards_gained >12 & (qb_scramble == 1 | rush == 1)),1,0),
#          negative = ifelse(yards_gained < 0, 1,0))

#Off vs Defensive Efficiency ----
total_offensive_efficiency <- pbp_rp %>%
  filter(season == year) %>%
  group_by(posteam) %>%
  summarize(offensive_epa = mean(epa))

total_defensive_efficiency <- pbp_rp %>%
  filter(season == year) %>%
  group_by(defteam) %>%
  summarize(defensive_epa = mean(epa))

total_efficiency_both <- total_offensive_efficiency %>%
  left_join(total_defensive_efficiency, by = c("posteam" = "defteam"))

total_efficiency_both <- total_efficiency_both %>%
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr"))

total_efficiency_both %>% 
  ggplot(aes(x = defensive_epa, y = offensive_epa)) +
  geom_nfl_logos(aes(team_abbr = posteam), width = 0.05, alpha = 0.95)+
  scale_x_reverse()+
  labs(x = "Defensive EPA/Play", y = "Offensive EPA/Play", title = "Efficiency Landscape",
       subtitle = "Dotted Lines Represent League Average",
       caption = "@CapAnalytics7 | nflfastR")+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 12),
        plot.caption = element_text(colour = "white", size = 8),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(face = "bold", colour = "white",size = 12),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1))+
  geom_hline(yintercept = mean(total_efficiency_both$offensive_epa), linetype = "dashed",color = "white")+
  geom_vline(xintercept = mean(total_efficiency_both$defensive_epa), linetype = "dashed", color = "white")
ggsave("EffLandscape.png", width = 14, height =10, dpi = "retina")



#EPA vs Success----
pbp_rp %>% 
  # group_by(posteam) %>%
  group_by(defteam) %>%
  summarize(success_rate = mean(success,na.rm = T), epa_play = mean(epa,na.rm = T)) %>% 
  ggplot(aes(x = success_rate, y = epa_play )) +
  # geom_nfl_logos(aes(team_abbr = posteam), width = 0.045, alpha = 0.95)+
  geom_nfl_logos(aes(team_abbr = defteam), width = 0.045)+
  # labs(x = "Success Rate", y = "EPA/Play", title = "Offensive Success Rate vs EPA/Play",caption = "@CapAnalytics7 | nflfastR")+
  labs(x = "Success Rate Allowed", y = "EPA/Play",subtite = "Dotted lines represent league average", title = "Defensive Success Rate vs EPA/Play",caption = "@CapAnalytics7 | nflfastR")+
  scale_x_reverse()+
  scale_y_reverse()+
  theme(legend.position = "top",
          legend.direction = "horizontal",
          legend.background = element_rect(fill = "white", color="white"),
          legend.title = element_blank(),
          legend.text = element_text(colour = "black", face = "bold"),
          plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
          plot.subtitle = element_text(hjust = .5, colour = "white", size = 12),
          plot.caption = element_text(colour = "white", size = 8),
          panel.grid = element_blank(),
          plot.background = element_rect(fill = "black", color="black"),
          panel.background = element_rect(fill = "black", color="black"),
          axis.ticks = element_line(color = "white"),
          axis.text = element_text(face = "bold", colour = "white",size = 12),
          axis.title = element_text(color = "white", size = 14),
          panel.border = element_rect(colour = "white", fill = NA, size = 1))+
  geom_smooth(se = F, method = "lm", color = "white")
ggsave("SuccessLandscape.png", width = 14, height =10, dpi = "retina")
#Side of Ball Breakdown----
pbp_rp %>% 
  filter(season == year) %>% 
  # group_by(posteam) %>%
  group_by(defteam) %>%
  summarize(epa_db = mean(epa[pass == 1],na.rm = T), epa_rush = mean(epa[rush == 1], na.rm = T)) %>% 
  ggplot(aes(x = epa_rush, y = epa_db))+
  # geom_nfl_logos(aes(team_abbr = posteam), width = 0.045, alpha = 0.95)+
  geom_nfl_logos(aes(team_abbr = defteam), width = 0.045, alpha = 0.955)+
  scale_x_reverse()+
  scale_y_reverse()+
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
  labs(x = "EPA/Rush", y = "EPA/Dropback", title = "Defensive Efficiency Landscape",subtitle = "Dotted lines represent league average", caption = "@CapAnalytics7 | nflfastR")+
  # labs(x = "EPA/Rush", y = "EPA/Dropback", title = "Offensive Efficiency Landscape",subtitle = "Dotted lines represent league average", caption = "@CapAnalytics7 | nflfastR")+
  geom_hline(yintercept = mean(pbp_rp$epa[pbp_rp$season == 2024 & pbp_rp$pass == 1],na.rm = T), linetype = "dashed",color = "white")+
  geom_vline(xintercept = mean(pbp_rp$epa[pbp_rp$season == 2024 & pbp_rp$rush == 1],na.rm = T), linetype = "dashed",color = "white")
ggsave("DefOffBreakdown.png", width = 14, height =10, dpi = "retina")



#NFL Offense Breakout----
pbp_rp %>%
  filter(rush == 1 | pass == 1, qb_kneel == 0, qb_spike == 0) %>%
  mutate(detailed_play_type = case_when(
    penalty == 1 ~ "Penalty",
    interception == 1 | fumble == 1 ~ "Turnover",
    rush == 1 ~ "Designed Run",
    qb_scramble == 1 ~ "QB Scramble",
    sack == 1 ~ "Sack",
    air_yards <= 0 ~ "At/Behind LOS Pass",
    air_yards > 0 & air_yards <= 10 ~ "0-10 Air Yard Pass",
    air_yards > 10 & air_yards <= 20 ~ "10-20 Air Yard Pass",
    air_yards > 20 ~ "20+ Air Yard Pass",
    TRUE ~ "Other"  # This acts as the catch-all for anything not matched
  )) %>%
  group_by(posteam,detailed_play_type) %>%
  summarize(total_epa = sum(epa,na.rm = T), epa_play = mean(epa,na.rm = T), count = n()) %>%
  mutate(sum_epa = sum(total_epa), total_count = sum(count), epa_tot_play = total_epa/total_count) %>% #epa_tot_play sums to epa/play
  ggplot(aes(x  = epa_tot_play, y =reorder(posteam,epa_tot_play), fill = detailed_play_type))+
  geom_bar(stat = "identity")+
  scale_fill_brewer(palette = "Set3") +
  # geom_nfl_logos(aes(team_abbr = max(posteam)), width = 0.05, alpha = 0.7)+
  labs(y = "Offense", x = "EPA/Play Components", title = "Where are Offenses Generating Success From?", subtitle = "Sections represent different components of an offense's EPA/Play",
       caption = "@CapAnalytics7 | nflfastR")+
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
        panel.border = element_rect(colour = "white", fill = NA, size = 1),
        axis.title.y = element_blank(),
        axis.text.y = element_nfl_logo(size = 0.9)
        )
ggsave("OffBreakout.png", width = 14, height =10, dpi = "retina")
 

#Early down vs Late Down Efficiency----
pbp_rp %>%
  filter(season == year) %>% 
  # group_by(posteam) %>%
  group_by(defteam) %>%
  summarize(early_down_epa = mean(epa[down<=2],na.rm = T), late_down_epa = mean(epa[down>2],na.rm = T)) %>% 
  ggplot(aes(x = early_down_epa, y = late_down_epa)) +
  # geom_nfl_logos(aes(team_abbr = posteam), width = 0.05)+
  geom_nfl_logos(aes(team_abbr = defteam), width = 0.045, alpha = 0.95)+
  scale_x_reverse()+
  scale_y_reverse()+
  labs(x = "EPA/Early Down (1st & 2nd down)", y = "EPA/Late Down (3rd & 4th down)", title = "Defensive Efficiency Late Down vs Early Down",
       subtitle = "Dotted lines represent average",
       caption = "@CapAnalytics7 | nflfastR")+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 10),
        plot.caption = element_text(colour = "white", size = 10),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(face = "bold", colour = "white",size = 12),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1),
        panel.grid = element_blank())+
  geom_hline(yintercept = mean(pbp_rp$epa[pbp_rp$down > 2],na.rm = T), linetype = "dashed",color = "white")+
  geom_vline(xintercept = mean(pbp_rp$epa[pbp_rp$down <= 2], na.rm = T), linetype = "dashed", color = "white")
ggsave("EarlyvsLateEfficiency.png", width = 14, height =10, dpi = "retina")


#Explosive vs Negative ----
pbp_rp %>% 
  filter(season == year) %>% 
  group_by(posteam) %>%
  # group_by(defteam) %>%
  summarize(negative_rate = mean(negative,na.rm = T),explosive_rate = mean(explosive,na.rm = T)) %>% 
  ggplot(aes(x = negative_rate, y = explosive_rate))+
  geom_point()+
  scale_x_reverse()+
  # scale_y_reverse()+
  # geom_nfl_logos(aes(team_abbr = defteam), width = 0.05)+
  geom_nfl_logos(aes(team_abbr = posteam), width = 0.05)+
  labs(x = "Negative Play Rate", y = "Explosive Play Rate*", title = "Which Offenses Create Explosive Plays and Prevent Negatives?", caption = "*Passes that gained greater than 20 yards or runs that gained greater than 12 yards                         @CapAnalytics7 | nflfastR", subtitle = "Dotted Lines Represent League Average")+
  # labs(x = "Negative Play Rate", y = "Explosive Play Rate*", title = "Which Defenses Prevent Explosives and Create Negatives?", caption = "*Passes that gained greater than 20 yards or runs that gained greater than 12 yards                         @CapAnalytics7 | nflfastR", subtitle = "Dotted Lines Represent League Average")+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 10),
        plot.caption = element_text(colour = "white", size = 10),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(face = "bold", colour = "white",size = 12),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1))+
  geom_hline(yintercept = mean(pbp_rp$explosive, na.rm = TRUE), linetype = "dashed",color = "white")+
  geom_vline(xintercept = mean(pbp_rp$negative, na.rm = TRUE), linetype = "dashed", color = "white")

ggsave("ExpvsNeg.png", width = 14, height =10, dpi = "retina")

#Ability to recover from negative plays----

negatives <- pbp %>% 
  filter(!is.na(series)) %>% 
  mutate(
    negative_play = ifelse((!is.na(penalty_team) & !is.na(penalty_yards) & penalty_team == posteam & penalty_yards >0) | (yards_gained < 0 & down <4), 1, 0),
    unique_series_identifier = paste(game_id, series)) %>%
  filter(qb_kneel !=1, !is.na(negative_play)) %>%
  # select(negative_play, desc,penalty,pass, rush)#check drive_points
  group_by(unique_series_identifier) %>% 
  summarize(posteam = first(posteam),negative_plays = max(negative_play),series_suc = max(series_success),
            total_neg = sum(negative_play)) %>% 
  group_by(negative_plays,posteam) %>% 
  summarize(success_series = mean(series_suc,na.rm = T), count = n(), total_negatives = sum(total_neg), average_negative_plays_neg_drive = mean(total_neg)) %>% 
  filter(!is.na(posteam)) %>% 
  mutate(negative_plays = ifelse(negative_plays == 1, "Yes", "No")) %>% 
  pivot_wider(values_from = c(count,success_series, total_negatives, average_negative_plays_neg_drive), names_from = negative_plays) %>% 
  mutate(avg_neg = total_negatives_Yes/(count_Yes+count_No))

negatives %>% 
  ggplot(aes(x = success_series_No, y = success_series_Yes))+
  geom_nfl_logos(aes(team_abbr = posteam), width = 0.05)+
  labs(x = "Series Success With No Negative Plays", y = "Series Success With 1+ Negative Plays",
       title = "Which Offenses are Best at Recovering from Negative Plays?", subtitle = "Series Success Defined as Series Results in First Down or Touchdown",
       caption = "@CapAnalytics7 | nflfastR" )+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 10),
        plot.caption = element_text(colour = "white", size = 10),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(face = "bold", colour = "white",size = 12),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1),
        panel.grid = element_blank())+
  geom_smooth(method = "lm", se = FALSE, color = "white", linetype = "dashed")
ggsave("NegativePlays.png", width = 14, height =10, dpi = "retina")

#Negative Series Frequency----
off_neg <- negatives %>% 
  mutate(neg_series_rate = (count_Yes/(count_No+count_Yes))) %>% 
  arrange(neg_series_rate) %>% 
  select(posteam, neg_series_rate,average_negative_plays_neg_drive_Yes) %>% 
  mutate_if(is.numeric, ~round(.,3)) %>% 
  gt() %>% 
  cols_align(align = "center") %>% 
  gt_nfl_wordmarks(columns = "posteam") %>% 
  cols_label(posteam = "Team",
             neg_series_rate = "Series W/Negative Play Rate",
             average_negative_plays_neg_drive_Yes = "Average Negative Plays Per Series W/Negative Play") %>% 
  gtExtras::gt_theme_538() %>%
  gtExtras::gt_hulk_col_numeric(c(neg_series_rate,average_negative_plays_neg_drive_Yes)) %>% 
  tab_header(
    title = md("Which Offenses Avoid Negative Plays?")
  )
gtsave(off_neg, "NegativeOffPerformance.png")
#Negative Plays on Early Downs, issue with negative plays per series with grouping



#1st Half Total Efficiency----
first_half_off <- pbp_rp %>% 
  filter(season == year) %>% 
  group_by(posteam)%>% 
  summarize(off_1h_epa = mean(epa[qtr %in% c(1,2)]), off_2h_epa = mean(epa[qtr > 2])) %>% 
  mutate(off_2h_impr = off_2h_epa-off_1h_epa)

first_half_def <- pbp_rp %>% 
  filter(season == year) %>% 
  group_by(defteam)%>% 
  summarize(def_1h_epa = mean(epa[qtr %in% c(1,2)]), def_2h_epa = mean(epa[qtr > 2],na.rm = T)) %>% 
  mutate(def_2h_impr = def_2h_epa-def_1h_epa)

total_first_half <- first_half_off %>% 
  left_join(first_half_def, by = c("posteam" = "defteam")) %>% 
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr"))



total_first_half %>% 
  ggplot(aes(x = def_2h_impr, y = off_2h_impr)) +
  geom_nfl_logos(aes(team_abbr = posteam), width = 0.055)+
  theme_bw()+
  scale_x_reverse()+
  labs(x = "2nd Half Defensive EPA/Play Improvement", y = "2nd Half Offensive EPA/Play Improvement", title = "How Do Teams Performances Differ From Between Halves?",
       subtitle = "Dotted Lines Represents League Average",
       caption = "@CapAnalytics7 | nflfastR")+
  geom_mean_lines(aes(x0 = def_2h_impr, y0 = off_2h_impr), linetype = "dashed", color = "white")+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 12),
        plot.caption = element_text(colour = "white", size = 8),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(face = "bold", colour = "white",size = 12),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1))
ggsave("First Half Performance.png", width = 14, height =10, dpi = "retina")



#QB Explosive Pass Play Rate vs aDot----

explosive_pass <- pbp %>% 
  filter(pass == 1, !is.na(passer_player_name), !is.na(air_yards)) %>% 
  mutate(explosive = ifelse(yards_gained >= 20, 1,0)) %>% 
  group_by(passer_id,passer_player_name, posteam) %>% 
  summarize(passes = n(), adot = mean(air_yards,na.rm = T), explosive_rate = mean(explosive)) %>% 
  filter(passes>=50)

explosive_pass %>% 
  ggplot(aes(x = adot, y = explosive_rate))+
  geom_nfl_logos(aes(team_abbr = posteam), width = 0.025, alpha = 0.7)+
  geom_text_repel(
    aes(label = passer_player_name),
    box.padding = 0.05,  # adjust this value for padding around the labels
    point.padding = 0.1,  # adjust this value for padding around the points
    segment.color = "grey",
    segment.size = 0.2,
    color = "white"
  )+
  labs(x = "Average Depth of Target",
       y = "Explosive Pass Rate*", title = "Explosive Pass Rate vs Average Depth of Target",
       caption = "*Passes that gain more than 20 yards      @CapAnalytics7 | nflfastR", subtitle = "Minimum 50 Dropbacks")+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 10),
        plot.caption = element_text(colour = "white", size = 10),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(face = "bold", colour = "white",size = 12),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1),
        panel.grid = element_blank())+
  geom_hline(yintercept = mean(pbp_rp$explosive[pbp_rp$pass == 1 & !is.na(pbp_rp$air_yards)],na.rm = T), linetype = "dashed",color = "white")+
  geom_vline(xintercept = mean(pbp_rp$air_yards,na.rm = T), linetype = "dashed", color = "white")
ggsave("ExplosivePassRate.png", width = 14, height =10, dpi = "retina")






#Pass Over Exp vs EPA/Exp ----
pbp_rp %>% 
  filter(season == year) %>% 
  group_by(posteam) %>% 
  summarize(xpassoe = mean(pass_oe,na.rm = T), epa_xpass = mean(epa[xpass>=0.9 & pass == 1],na.rm = T)) %>% 
  ggplot(aes(x = xpassoe, y = epa_xpass))+
  geom_nfl_logos(aes(team_abbr = posteam), width = 0.05)+
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
  labs(x = "Pass Rate Over Expectation", y = "EPA/Expected Pass Dropback*", title = "Which Teams Are Passing More than Expected? How do They Perform in Expected Passing Situations?",
       subtitle = "Dotted lines represent league average", 
       caption = "*Expected pass situation is down with xpass >0.9      @CapAnalytics7 | nflfastR")+
  geom_hline(yintercept = mean(pbp_rp$epa[pbp_rp$season == 2024 & pbp_rp$pass == 1 & pbp_rp$xpass>=0.9],na.rm = T), linetype = "dashed",color = "white")+
  geom_vline(xintercept = mean(pbp_rp$pass_oe[pbp_rp$season == 2024],na.rm = T), linetype = "dashed",color = "white")
ggsave("xPass.png", width = 14, height =10, dpi = "retina")

#Non Red vs Red----

test <- pbp_rp %>% 
  mutate(unique_series_identifier = paste(game_id, series)) %>%
  group_by(unique_series_identifier) %>% 
  summarize(posteam = first(posteam),drive_20 = max(drive_inside20),drive_td = max(drive_end_transition)) %>% 
  mutate(touchdown = ifelse(drive_td == "TOUCHDOWN",1,0)) %>% 
  group_by(posteam) %>% 
  summarize(red_td_rate = mean(touchdown[drive_20 == 1], na.rm = T))

pbp_rp %>% 
  group_by(defteam) %>%
  # group_by(posteam) %>%
  summarize(epa_red = mean(epa[yardline_100<= 20],na.rm = T), epa_non_red = mean(epa[yardline_100> 20],na.rm = T)) %>% 
  left_join( pbp_rp %>% 
               mutate(unique_series_identifier = paste(game_id, series)) %>%
               group_by(unique_series_identifier) %>% 
               summarize(posteam = first(posteam),defteam = first(defteam),drive_20 = max(drive_inside20),drive_td = max(drive_end_transition)) %>% 
               mutate(touchdown = ifelse(drive_td == "TOUCHDOWN",1,0)) %>% 
               # group_by(posteam) %>% 
               group_by(defteam) %>%
               summarize(red_td_rate = mean(touchdown[drive_20 == 1], na.rm = T)), 
             # by = c("posteam")
             by = c("defteam")
             ) %>% 
  ggplot(aes(x = red_td_rate, y = epa_non_red))+
  # geom_nfl_logos(aes(team_abbr = posteam), width = 0.04, alpha = 0.95)+
  geom_mean_lines(aes(x0 = red_td_rate, y0 = epa_non_red), color = "white", linetype = "dashed")+
  geom_nfl_logos(aes(team_abbr = defteam), width = 0.04, alpha = 0.95)+
  scale_x_reverse()+
  scale_y_reverse()+
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
  # labs(x = "Redzone TD Rate", y = "Non Redzone EPA/Play", title = "Offense Efficiency Inside vs Outside Red Zone", subtitle = "Dotted lines represent league average", caption = "@CapAnalytics7 | nflfastR")
  labs(x = "Redzone TD Rate Allowed", y = "Non Redzone EPA/Play", title = "Defensive Efficiency Inside vs Outside Red Zone", subtitle = "Dotted lines represent league average", caption = "@CapAnalytics7 | nflfastR")
ggsave("RedBreakOut.png", width = 14, height =10, dpi = "retina")

#Biggest Plays ----
big_play <- pbp_rp %>% 
  filter(season == 2024) %>% 
  filter(week ==9) %>%
  select(desc,wpa) %>% 
  mutate(wpa = abs(wpa)) %>% arrange(-wpa)

#Win Probability----

home_wp <- pbp %>% 
  filter(season == year) %>%
  group_by(game_seconds_remaining,home_team) %>% 
  summarize(avg_home_wp = mean(home_wp), home_count = n())

away <- pbp %>% 
  group_by(game_seconds_remaining,away_team) %>% 
  filter(season == year) %>%
  summarize(avg_away_wp = mean(away_wp), away_count = n())

wp_total <- home_wp %>% 
  full_join(away, by = c("home_team" = "away_team","game_seconds_remaining" = "game_seconds_remaining")) %>%
  mutate(avg_home_wp = ifelse(is.na(avg_home_wp),0,avg_home_wp), avg_away_wp = ifelse(is.na(avg_away_wp),0,avg_away_wp), 
         away_count = ifelse(is.na(away_count),0,away_count), home_count = ifelse(is.na(home_count),0,home_count)) %>% 
  mutate(avg_total_wp = (home_count*avg_home_wp+away_count*avg_away_wp)/(home_count+away_count)) %>% 
  left_join(teams_colors_logos, by = c("home_team" = "team_abbr"))

wp_total %>%
  ggplot(aes(x = game_seconds_remaining, y = avg_total_wp, group = home_team)) +
  stat_smooth(se = FALSE, show.legend = FALSE, aes(color = team_color2, fill = team_color), span = 0.25, method = "loess", lwd = 2) +
  xlim(0, 3600) +
  scale_x_reverse() +
  scale_color_identity() +
  geom_vline(xintercept = 2700,color = "white") +
  geom_vline(xintercept = 1800 ,color = "white") +
  geom_vline(xintercept = 900,color = "white") +# Vertical line at x = 1800+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 12),
        plot.caption = element_text(colour = "white", size = 10),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = alpha("black",0.4), color="black") ,
        panel.background = element_rect(fill = alpha("black",0.4), color="black"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(face = "bold", colour = "white",size = 12),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1),
        strip.text = nflplotR::element_nfl_wordmark(size = 1))+
  facet_wrap(~ home_team, ncol = 8, nrow =4) +  # Create separate facets for each home_team
  labs(x = "Time Remaining in Game",
       y = "Average Win Probability",
       title = "How Do Games Play Out?", subtitle = "Vertical Lines Denote End of Quarters",
       caption = "@CapAnalytics7 | nflfastR")
ggsave("WinProbvsTime.png", width = 14, height =10, dpi = "retina")


#Home vs Road NFL Efficiency----
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
       y = "Home EPA/Play", title = "Home vs Away Offensive Efficiency Following Week 17",
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR", subtitle = "Dashed Lines Represent League Average")
ggsave("HomevsAwayOffense.png", width = 14, height =10, dpi = "retina")


#Running Back Efficiency----
rushing_player<-pbp_rp %>% 
  filter(rush == 1, !is.na(rusher_player_name)) %>% 
  filter(season == year) %>% 
  group_by(rusher_player_id,rusher_player_name) %>% 
  summarize(rusher_player_name = max(rusher_player_name),posteam = last(posteam),rushes = n(), rush_epa = mean(epa), success_rate = mean(success)) %>% 
  filter(rushes>=50) %>% 
  mutate(rusher_player_name = ifelse(rusher_player_name == "T.Hill", "Taysom Hill", rusher_player_name)) %>% 
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr"))
    



rushing_player %>%   
  ggplot(aes(x = success_rate, y = rush_epa))+
  geom_nfl_logos(aes(team_abbr = posteam), width = 0.02)+
  labs(x = "Success Rate", y = "EPA/Rush", title = "Rushing Efficiency by Rusher Following Week 3", subtitle = "Minimum 50 Rushes",
       caption = "@CapAnalytics7 | nflfastR")+
  theme_bw()+
  geom_text_repel(
    aes(label = rusher_player_name),
    box.padding = 0.05,  # adjust this value for padding around the labels
    point.padding = 0.01,  # adjust this value for padding around the points
    segment.color = "grey",
    segment.size = 0.2,
    color = "white",
    size = 5
  )+
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
  geom_hline(yintercept = mean(rushing_player$rush_epa,na.rm = T), linetype = "dashed",color = "white")+
  geom_vline(xintercept = mean(rushing_player$success_rate,na.rm = T), linetype = "dashed",color = "white")
ggsave("RunningEfficiency.png", width = 14, height =10, dpi = "retina")


#Motion Usage----
pbp_rp %>% 
  group_by(posteam) %>% 
  summarize(motion_rate = mean(is_motion,na.rm = T), epa_motion = mean(epa[is_motion == 1],na.rm = T), epa_no_motion = mean(epa[is_motion == 0],na.rm = T)) %>% 
  mutate(epa_change = ((epa_motion - epa_no_motion)/abs(epa_no_motion))*100) %>% 
  ggplot(aes(x = motion_rate, y = epa_change))+
  geom_nfl_logos(aes(team_abbr = posteam), width = 0.05)+
  geom_mean_lines(aes(x0 = motion_rate, y0 = epa_change), color = "white", linetype = "dotted")+
  labs(x = "Motion Rate", y = "% Change in EPA/Play With Motion", title = "Which Teams Should Increase/Decrease Their Motion Usage?",
       subtitle = "Dotted Lines Represent League Average", caption = "@CapAnalytics7 | nflfastR")+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 12),
        plot.caption = element_text(colour = "white", size = 8),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(face = "bold", colour = "white",size = 12),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1))
ggsave("MotionRate.png", width = 14, height =10, dpi = "retina")


#Play Action Rate---- 
pbp_rp %>% 
  filter(season == year) %>% 
  group_by(posteam) %>% 
  summarize(pa_rate = mean(is_play_action[pass == 1],na.rm = T), pa24 = mean(epa[is_play_action == 1],na.rm = T),
            no_pa = mean(epa[is_play_action == 0 & pass == 1], na.rm = T)) %>% 
  mutate(pa_improve = (pa24 - no_pa)/abs(no_pa)*100) %>% 
  ggplot(aes(x = pa_rate, y = pa_improve))+
  geom_nfl_logos(aes(team_abbr = posteam), width = 0.045,alpha = 0.95)+
  geom_mean_lines(aes(x0 = pa_rate, y0 = pa_improve),color = "white", linetype = "dashed")+
  labs(x = "Play Action Rate", y = "% Change in EPA/Play When Using Play Action", title = "Which Teams Should Increase/Decrease Play Action Usage?", caption ="@CapAnalytics7 | nflfastR",
       subtitle = "Dotted Lines Represent League Average")+
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
        panel.border = element_rect(colour = "white", fill = NA, size = 1))
ggsave("PARate.png", width = 14, height =10, dpi = "retina")


#Lead vs Trailing----
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
  labs(x = "EPA Per Play When Trailing", y = "EPA Per Playing When Tied or Leading", 
       title = "Offensive Efficiency When Leading/Tied vs Trailing Following Week 17",
       caption = "@CapAnalytics7 | nflfastR")+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 10),
        plot.caption = element_text(colour = "white", size = 10),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(face = "bold", colour = "white",size = 12),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1),
        panel.grid = element_blank())
ggsave("LeadingvsTrailing.png", width = 14, height =10, dpi = "retina")





# Yards After Catch----
yac_passing <- pbp_rp %>% 
  filter(pass == 1) %>% 
  filter(season == year) %>%
  group_by(passer_player_id,passer_player_name,posteam) %>% 
  summarise(count = n(), total_passing_yards = sum(yards_gained, na.rm =T), yac = sum(yards_after_catch, na.rm = T), epa_pass = mean(epa),
            adot = mean(air_yards,na.rm =T)) %>% 
  mutate(yac_pct = yac/total_passing_yards) %>% 
  filter(count >= 50 & !is.na(passer_player_id))

yac_passing <- yac_passing %>% 
  left_join(teams_colors_logos ,by = c("posteam" = "team_abbr"))

yac_passing %>% 
  ggplot(aes(x=yac_pct, y = epa_pass))+
  geom_nfl_logos(aes(team_abbr = posteam), width = 0.05)+
  labs(x = "% of Passing Yards from YAC", y = "EPA/Dropback", title = "EPA Per Pass Attempt vs % of Passing Yards from YAC", subtitle = "Minimum 50 Pass Attempts",
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")+
  # geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  geom_text_repel(
    aes(label = passer_player_name),
    box.padding = 0.05,  # adjust this value for padding around the labels
    point.padding = 0.01,  # adjust this value for padding around the points
    segment.color = "grey",
    segment.size = 0.2,
    color = "white"
  )+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 10),
        plot.caption = element_text(colour = "white", size = 10),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(face = "bold", colour = "white",size = 12),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1),
        panel.grid = element_blank())+
  geom_smooth(method = "lm", se = FALSE, color = "white", linetype = "dashed")
ggsave("YAC.png", width = 14, height =10, dpi = "retina")

#Time to throw vs aDoT----
pbp_rp %>% 
  filter(season == year) %>% 
  group_by(posteam) %>% 
  summarize(aDot = mean(air_yards,na.rm = T), avg_time_throw = mean(time_to_throw,na.rm = T)) %>% 
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr")) %>% 
  ggplot(aes(x = aDot, y = avg_time_throw))+
  geom_point()+
  scale_x_reverse()+
  geom_image(aes(image = team_logo_espn), size = 0.05, asp = 16/9)+
  labs(x = "Negative Play Rate", y = "Explosive Play Rate", title = "How has Play Action Rate Rate Changed for Teams from Prior Season?")+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 10),
        plot.caption = element_text(colour = "white", size = 8),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(face = "bold", colour = "white",size = 12),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1))




#Middle 8 ----

#Decay Efficiency Landscape----
decay_rate <- 0.95
weight_epa <- pbp_rp %>%
  filter(season == year) %>%
  mutate(decay_factor = decay_rate^ ((max(week) - week)),
         weighted_epa = epa * decay_factor) %>% 
  group_by(posteam) %>%
  summarize(offensive_epa = sum(weighted_epa)/sum(decay_factor), off_nw_epa = mean(epa,na.rm = T)) %>% 
  left_join(pbp_rp %>%
  filter(season == year) %>%
  mutate(decay_factor = decay_rate^ ((max(week) - week)),
         weighted_epa = epa * decay_factor) %>% 
  group_by(defteam) %>%
  summarize(defensive_epa = sum(weighted_epa)/sum(decay_factor),def_nw_epa = mean(epa,na.rm = T)), by = c("posteam" = "defteam"))
weight_epa %>% 
  ggplot(aes(x = defensive_epa, y = offensive_epa)) +
  geom_nfl_logos(aes(team_abbr = posteam), width = 0.04,alpha = 0.8)+
  scale_x_reverse()+
  labs(x = "Weighted Defensive EPA/Play", y = "Weighted Offensive EPA/Play", title = "What Have You Done Lately? Time Weighted Efficiency Landscane",
       subtitle = "Logos represent weighted EPA/Play, dots represent season long EPA/Play",
       caption = "Weighted EPA/Play factors performance in recent weeks more; Decay Rate is 0.95                   @CapAnalytics7 | nflfastR")+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 12),
        plot.caption = element_text(colour = "white", size = 8),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(face = "bold", colour = "white",size = 12),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1))+
  geom_point(aes(x = def_nw_epa, y = off_nw_epa), color = "white",fill = "blue", shape = 1, size = 8)+
  scale_color_nfl(type = "primary")+
  geom_segment(aes(x = def_nw_epa, y = off_nw_epa, 
                   xend = defensive_epa, yend = offensive_epa), 
               color = "lightgrey", alpha = 0.8, lwd =1.5)# Add points for regular EPA
ggsave("WeightedLandscape.png", width = 14, height =10, dpi = "retina")

#Field Position----
pbp %>% 
  mutate(drive_start = ifelse(yrdln == drive_start_yard_line, yardline_100, -1)) %>% 
  mutate(unique_drive = paste(game_id,drive)) %>% 
  group_by(unique_drive) %>%summarize(strat_field = max(drive_start,na.rm = T),posteam = max(posteam, na.rm = T)) %>% 
filter(strat_field>0) %>% 
  ungroup() %>% 
  group_by(posteam) %>% 
  filter(!is.na(posteam)) %>% 
  summarize(off_start = mean(strat_field)) %>% 
  left_join(pbp %>% 
              mutate(drive_start = ifelse(yrdln == drive_start_yard_line, yardline_100,-1)) %>% 
              mutate(unique_drive = paste(game_id,drive)) %>% 
              group_by(unique_drive) %>% 
              summarize(strat_field = max(drive_start,na.rm = T),defteam = max(defteam)) %>% 
              filter(strat_field>0) %>% 
              ungroup() %>% 
              group_by(defteam) %>% 
              filter(!is.na(defteam)) %>% 
            summarize(def_start = mean(strat_field)), by = c("posteam" = "defteam")) %>% 
  ggplot(aes(x = def_start, y = off_start))+
  geom_nfl_logos(aes(team_abbr = posteam), width = 0.05)+
  scale_y_reverse()+
  geom_mean_lines(aes(x0 = def_start, y0 = off_start),color = "white", linetype = "dashed")+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 10),
        plot.caption = element_text(colour = "white", size = 10),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(face = "bold", colour = "white",size = 12),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1),
        panel.grid = element_blank())+
  labs(x = "Defensive Starting Field Position (Yards to Goalline)", y = "Offensive Starting Field Position (Yards to Goaline)",
       title = "Which Teams Have Benefitted the Most From Field Position?", subtitle = "Dotted Lines Represent League Averages")
ggsave("FieldPosition.png", width = 14, height =10, dpi = "retina")
#Maybe add expected ponts at start of drive?

#Stadium Location Performance ----
pbp_rp %>% 
  mutate(outdoor = ifelse(roof == "outdoors", "outdoor","indoor")) %>% 
  group_by(posteam,outdoor) %>% 
  summarize(epa_play = mean(epa)) %>% 
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr")) %>% 
  pivot_wider(names_from = outdoor, values_from = epa_play) %>% 
  ggplot(aes(x = indoor, y = outdoor))+
  geom_image(aes(image = team_logo_espn), size = 0.05, asp = 16/9)+
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") 

#Offense Predictability----
pbp_rp %>% 
  group_by(posteam) %>% 
  summarize(given_no_uc_pass_rate = mean(pass_oe[qb_location %in% c("S","U")],na.rm = T),
            given_under_center_rush_rate = mean(pass_oe[qb_location == "U"],na.rm = T)) %>% 
  summarize(given_no_uc_pass_rate = mean(pass[qb_location %in% c("S","U")],na.rm = T),
            given_under_center_rush_rate = mean(rush[qb_location == "U"],na.rm = T)) %>% 
  ggplot(aes(x = given_no_uc_pass_rate, y = given_under_center_rush_rate))+
  geom_nfl_logos(aes(team_abbr = posteam), width = 0.04,alpha = 0.95)+
  geom_mean_lines(aes(x0 = given_no_uc_pass_rate, y0 = given_under_center_rush_rate), linetype = "dashed", color = "white")+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 10),
        plot.caption = element_text(colour = "white", size = 10),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_line(color = "white"),
        axis.text = element_text(face = "bold", colour = "white",size = 12),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1),
        panel.grid = element_blank())+
  labs(title = "How Predictable is The Offensive Playcalling?", subtitle = "Dotted Lines Represent League Average", 
       caption =  "@CapAnalytics7 | nflfastR", x = "Non-Under Center Pass Rate", y = "Under Center Rush Rate")

#Strength of Schedule----
strength_efficiency <- pbp_rp %>%
  filter(season == year) %>%
  group_by(posteam) %>%
  summarize(offensive_epa = mean(epa)) %>% 
  left_join(pbp_rp %>%
              filter(season == year) %>%
              group_by(defteam) %>%
              summarize(defensive_epa = mean(epa)), by = c("posteam" = "defteam"))
replace_with_ranks<- function(column){
  values <- column
  ranks <- as.numeric(rank(column*-1,ties.method = "max"))
}
strength_ranks <- apply(strength_efficiency %>% 
        mutate(defensive_epa = defensive_epa*-1) %>% select(-posteam), 2, replace_with_ranks)
strength_teams <- as.data.frame(cbind(strength_efficiency$posteam,strength_ranks)) %>% 
  rename("team" = "V1") %>% 
  mutate(offensive_epa = as.numeric(offensive_epa), defensive_epa = as.numeric(defensive_epa))
schedule_strength <- pbp_rp %>% 
  group_by(game_id) %>% 
  summarize(home_team = max(home_team), away_team = max(away_team)) %>% 
  pivot_longer(cols = c(home_team,away_team),names_to = "type", values_to = "group_team") %>% 
  left_join(pbp_rp %>% 
              group_by(game_id) %>% 
              summarize(home_team = max(home_team), away_team = max(away_team)), by = c("game_id")) %>% 
  mutate(opponent = ifelse(home_team == group_team, away_team, home_team)) %>% 
  select(-home_team,-away_team) %>% 
  left_join(strength_teams, by = c("opponent" = "team")) %>% 
  group_by(group_team) %>% 
  summarize(past_off_rank = mean(offensive_epa), past_def_rank = mean(defensive_epa))

schedule_strength %>% 
  ggplot(aes(x = past_def_rank, y = past_off_rank)) +
  geom_nfl_logos(aes(team_abbr = group_team),width = 0.045, alpha = 0.95) +
  geom_mean_lines(aes(x0 = def_rank, y0 = off_rank), color = "white", linetype = "dotted")+
  scale_x_reverse()+
  scale_y_reverse()+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 10),
        plot.caption = element_text(colour = "white", size = 10),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1),
        panel.grid = element_blank())+
  labs(y = "Offenses Played Strength", x = "Defenses Played Strength", title = "How Hard is a Team's Schedule So Far?", caption = "@CapAnalytics7 | nflfastR", subtitle = "Dotted Lines Represent League Average")+
  annotate("text",label = "Played Good Offenses, Played Bad Defenses", y = quantile(schedule_strength$off_rank,0.01), x = quantile(schedule_strength$def_rank,0.95),color = "white")+
  annotate("text",label = "Played Bad Offenses, Played Bad Defenses", y = quantile(schedule_strength$off_rank,0.95), x = quantile(schedule_strength$def_rank,0.95),color = "white")+
  annotate("text",label = "Played Good Offenses, Played Good Defenses", y = quantile(schedule_strength$off_rank,0.01), x = quantile(schedule_strength$def_rank,0.05),color = "white")+
  annotate("text",label = "Played Bad Offenses, Played Good Defenses", y = quantile(schedule_strength$off_rank,0.97), x = quantile(schedule_strength$def_rank,0.03),color = "white")
ggsave("SOS.png", width = 14, height =10, dpi = "retina")


#Forward Looking----
sched24 <- load_schedules(2024)
upcoming_sos <- sched24 %>% 
  filter(week > max(pbp_rp$week)) %>% 
  pivot_longer(cols = c(home_team,away_team),names_to = "type", values_to = "group_team") %>% 
  left_join(sched24, by = c("game_id")) %>% 
  mutate(opponent = ifelse(home_team == group_team, away_team, home_team)) %>%
  select(game_id, group_team, opponent,home_team,away_team) %>% 
  left_join(strength_teams, by = c("opponent" = "team")) %>% 
  group_by(group_team) %>% 
  summarize(up_off_rank = mean(offensive_epa), up_def_rank = mean(defensive_epa))

upcoming_sos %>% 
  ggplot(aes(x = up_def_rank, y = up_off_rank)) +
  geom_nfl_logos(aes(team_abbr = group_team),width = 0.045, alpha = 0.95) +
  geom_mean_lines(aes(x0 = def_rank, y0 = off_rank), color = "white", linetype = "dotted")+
  scale_x_reverse()+
  scale_y_reverse()+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 10),
        plot.caption = element_text(colour = "white", size = 10),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1),
        panel.grid = element_blank())+
  labs(y = "Upcoming Offensive Strength", x = "Upcoming Defensive Strength", title = "How Difficult is a Team's Upcoming Schedule?", caption = "@CapAnalytics7 | nflfastR", subtitle = "Dotted Lines Represent League Average")+
  annotate("text",label = "Upcoming Good Offense\nUpcoming Bad Defenses", y = quantile(upcoming_sos$off_rank,0.02), x = quantile(upcoming_sos$def_rank,0.98),color = "white")+
  annotate("text",label = "Upcoming Bad Offenses\nUpcoming Bad Defenses", y = quantile(upcoming_sos$off_rank,0.97), x = quantile(upcoming_sos$def_rank,0.98),color = "white")+
  annotate("text",label = "Upcoming Good Offenses\nUpcoming Good Defenses", y = quantile(upcoming_sos$off_rank,0.02), x = quantile(upcoming_sos$def_rank,0.02),color = "white")+
  annotate("text",label = "Upcoming Bad Offenses\nUpcoming Good Defenses", y = quantile(upcoming_sos$off_rank,0.97), x = quantile(upcoming_sos$def_rank,0.02),color = "white")
ggsave("UpcomingSOS.png", width = 14, height =10, dpi = "retina")

#SOS Difference----
sos_diff <- upcoming_sos %>% 
  left_join(schedule_strength, by = c("group_team")) %>% 
  mutate(schedule_diff_def = up_def_rank-past_def_rank,schedule_diff_off = up_off_rank-past_def_rank)
sos_diff %>% 
  ggplot(aes(x = schedule_diff_def, y = schedule_diff_off))+
  geom_nfl_logos(aes(team_abbr = group_team),width = 0.045, alpha = 0.95) +
  geom_mean_lines(aes(x0 = schedule_diff_def, y0 = schedule_diff_off), color = "white", linetype = "dotted")+
  scale_x_reverse()+
  scale_y_reverse()+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 10),
        plot.caption = element_text(colour = "white", size = 10),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(color = "white", size = 14),
        panel.border = element_rect(colour = "white", fill = NA, size = 1),
        panel.grid = element_blank())+
  labs(y = "Strength Difference Between Upcoming vs Played Offenses", x = "Strength Difference Between Upcoming vs Played Offenses", title = "Which Teams See the Biggest Shift in Schedule Difficulty?", caption = "@CapAnalytics7 | nflfastR", subtitle = "Dotted Lines Represent League Average")+
  annotate("text",label = "Better Offenses\nWorse Defenses", y = quantile(sos_diff$schedule_diff_off,0.02), x = quantile(sos_diff$schedule_diff_def,0.98),color = "white")+
  annotate("text",label = "Worse Offenses\nWorse Defenses", y = quantile(sos_diff$schedule_diff_off,0.97), x = quantile(sos_diff$schedule_diff_def,0.98),color = "white")+
  annotate("text",label = "Better Offenses\nBetter Defenses", y = quantile(sos_diff$schedule_diff_off,0.02), x = quantile(sos_diff$schedule_diff_def,0.02),color = "white")+
  annotate("text",label = "Worse Offenses\nBetter Defenses", y = quantile(sos_diff$schedule_diff_off,0.97), x = quantile(sos_diff$schedule_diff_def,0.02),color = "white")
ggsave("DiffSOS.png", width = 14, height =10, dpi = "retina")


#Adjusted Indexing----
replace_with_index<- function(column){
  index <- column/mean(column)
}
strength_efficiency_ind <- pbp_rp %>%
  filter(season == year) %>%
  group_by(posteam) %>%
  summarize(offensive_epa = mean(epa)) %>% 
  left_join(pbp_rp %>%
              filter(season == year) %>%
              group_by(defteam) %>%
              summarize(defensive_epa = mean(epa)), by = c("posteam" = "defteam")) %>% 
  mutate(def_index = scale(defensive_epa), off_index = offensive_epa/mean(offensive_epa), def_check_mean = scale(defensive_epa)*-1, def_check = defensive_epa)


#Lucky----
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

lucky_tab <- lucky_data %>% 
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
  # select(" ", `Net Lucky Wins By EPA`,`Lucky Wins By EPA`, `Unlucky Losses by EPA`, `Net Lucky Wins By Success Rate`,`Lucky Wins by Success Rate`, `Unlucky Losses by Success Rate`) %>% 
  select(" ", `Net Lucky Wins By Success Rate`,`Lucky Wins by Success Rate`, `Unlucky Losses by Success Rate`) %>% 
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
gtsave(lucky_tab, "lucky_table.png")
