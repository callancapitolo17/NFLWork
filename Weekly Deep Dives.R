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

nfl99all <- load_pbp(1999:2023)
nfl99 <- nfl99all %>% 
  filter(pass == 1 | rush == 1) %>% 
  mutate(explosive = ifelse((yards_gained>20 & pass_attempt == 1) | (yards_gained >12 & (qb_scramble == 1 | rush == 1)),1,0),
         negative = ifelse(yards_gained < 0, 1,0))

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
  geom_image(aes(image = team_logo_espn), size = 0.05, asp = 16/9)+
  theme_bw()+
  scale_x_reverse()+
  labs(x = "Defensive EPA/Play", y = "Offensive EPA/Play", title = "Efficiency Landscape Following Week 4",
       subtitle = "Dotted Lines Represent League Average",
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")+
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



#Side of Ball Breakdown----
pbp_rp %>% 
  filter(season == year) %>% 
  # group_by(posteam) %>%
  group_by(defteam) %>%
  summarize(epa_db = mean(epa[pass == 1],na.rm = T), epa_rush = mean(epa[rush == 1], na.rm = T)) %>% 
  # left_join(teams_colors_logos ,by = c("posteam" = "team_abbr")) %>%
  left_join(teams_colors_logos ,by = c("defteam" = "team_abbr")) %>%
  ggplot(aes(x = epa_rush, y = epa_db))+
  scale_x_reverse()+
  scale_y_reverse()+
  geom_image(aes(image = team_logo_espn), size = 0.05, asp = 16/9)+
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
  labs(x = "EPA/Rush", y = "EPA/Dropback", title = "Defensive Efficiency Landscape",
       subtitle = "Dotted lines represent league average", 
       caption = "@CapAnalytics7 | nflfastR")+
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
  
  summarize(total_epa = sum(epa,na.rm = T)) %>%
  mutate(sum_epa = sum(total_epa)) %>%
  ggplot(aes(x  = total_epa, y =reorder(posteam,sum_epa), fill = detailed_play_type))+
  geom_bar(stat = "identity")+
  scale_fill_brewer(palette = "Set3") +
  # geom_nfl_logos(aes(team_abbr = max(posteam)), width = 0.05, alpha = 0.7)+
  labs(y = "Offense", x = "Total EPA", title = "Where are Offenses Generating Success From?",
       
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
  group_by(posteam) %>%
  # group_by(defteam) %>%
  summarize(early_down_epa = mean(epa[down<=2],na.rm = T), late_down_epa = mean(epa[down>2],na.rm = T)) %>% 
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr")) %>%
  # left_join(teams_colors_logos, by = c("defteam" = "team_abbr")) %>% 
  ggplot(aes(x = early_down_epa, y = late_down_epa)) +
  geom_image(aes(image = team_logo_espn), size = 0.1, asp = 16/19)+
  theme_bw()+
  # scale_x_reverse()+
  # scale_y_reverse()+
  labs(x = "EPA/Early Down (1st & 2nd down)", y = "EPA/Late Down (3rd & 4th down)", title = "Defensive Efficiency Late Down vs Early Down Following Week 3",
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
  # geom_nfl_logos(aes(team_abbr = defteam), width = 0.065)+
  geom_nfl_logos(aes(team_abbr = posteam), width = 0.065)+
  labs(x = "Negative Play Rate", y = "Explosive Play Rate*", title = "Which Offenses Avoid Negative Plays and Create Explosives?",
       caption = "*Passes that gained greater than 20 yards or runs that gained greater than 12 yards                         @CapAnalytics7 | nflfastR",
       subtitle = "Dotted Lines Represent League Average")+
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

#1st Half Total Efficiency----
first_half_off <- pbp_rp %>% 
  filter(qtr %in% c(1,2)) %>% 
  filter(season == year) %>% 
  group_by(posteam)%>% 
  summarize(off_epa = mean(epa))

first_half_def <- pbp_rp %>% 
  filter(qtr %in% c(1,2)) %>% 
  filter(season == year) %>% 
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
  labs(x = "First Half Defensive EPA/Play", y = "First Half Offensive EPA/Play", title = "1st Half Offensive and Defensive Efficiency Following Week 1",
       subtitle = "Dotted Lines Represents League Average",
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
  geom_hline(yintercept = mean(total_first_half$off_epa), linetype = "dashed",color = "white")+
  geom_vline(xintercept = mean(total_first_half$def_epa), linetype = "dashed", color = "white")
ggsave("First Half Performance.png", width = 14, height =10, dpi = "retina")



#QB Explosive Pass Play Rate vs aDot----

explosive_pass <- pbp %>% 
  filter(pass == 1, !is.na(passer_player_name), !is.na(air_yards)) %>% 
  mutate(explosive = ifelse(yards_gained >= 20, 1,0)) %>% 
  group_by(passer_id,passer_player_name, posteam) %>% 
  summarize(passes = n(), adot = mean(air_yards,na.rm = T), explosive_rate = mean(explosive)) %>% 
  filter(passes>=50)

explosive_pass <- explosive_pass %>% 
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr"))

explosive_pass %>% 
  ggplot(aes(x = adot, y = explosive_rate))+
  geom_image(aes(image = team_logo_espn), size = 0.03, asp = 16/9)+
  geom_text_repel(
    aes(label = passer_player_name),
    box.padding = 0.05,  # adjust this value for padding around the labels
    point.padding = 0.1,  # adjust this value for padding around the points
    segment.color = "grey",
    segment.size = 0.2,
    color = "white"
  )+
  labs(x = "Average Depth of Target",
       y = "Explosive Pass Rate", title = "Explosive Pass Rate vs Average Depth of Target",
       caption = "@CapAnalytics7 | nflfastR", subtitle = "Minimum 50 Passes")+
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
  left_join(teams_colors_logos ,by = c("posteam" = "team_abbr")) %>% 
  ggplot(aes(x = xpassoe, y = epa_xpass))+
  geom_image(aes(image = team_logo_espn), size = 0.05, asp = 16/9)+
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
  labs(x = "Pass Rate Over Expectation", y = "EPA/Expected Pass Dropback*", title = "Which Teams Are Passing More than Expected How do They Perform in Expected Passing Situations?",
       subtitle = "Dotted lines represent league average", 
       caption = "*Expected pass situation is down with xpass >0.9      @CapAnalytics7 | nflfastR")+
  geom_hline(yintercept = mean(pbp_rp$epa[pbp_rp$season == 2024 & pbp_rp$pass == 1 & pbp_rp$xpass>=0.9],na.rm = T), linetype = "dashed",color = "white")+
  geom_vline(xintercept = mean(pbp_rp$pass_oe[pbp_rp$season == 2024],na.rm = T), linetype = "dashed",color = "white")
ggsave("xPass.png", width = 14, height =10, dpi = "retina")

#Non Red vs Red----
pbp_rp %>% 
  # group_by(defteam) %>% 
  group_by(posteam) %>%
  summarize(epa_red = mean(epa[yardline_100<= 20],na.rm = T), epa_non_red = mean(epa[yardline_100> 20],na.rm = T)) %>% 
  left_join(teams_colors_logos ,by = c("posteam" = "team_abbr")) %>%
  # left_join(teams_colors_logos ,by = c("defteam" = "team_abbr")) %>%
  ggplot(aes(x = epa_red, y = epa_non_red))+
  geom_image(aes(image = team_logo_espn), size = 0.05, asp = 16/9)+
  # scale_x_reverse()+
  # scale_y_reverse()+
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
  labs(x = "EPA/Redzone", y = "EPA/Outside Redzone", title = "Offensive Efficiency Inside vs Outside Red Zone",
       subtitle = "Dotted lines represent league average", 
       caption = "@CapAnalytics7 | nflfastR")+
  geom_hline(yintercept = mean(pbp_rp$epa[pbp_rp$yardline_100>20],na.rm = T), linetype = "dashed",color = "white")+
  geom_vline(xintercept = mean(pbp_rp$epa[pbp_rp$yardline_100<=20],na.rm = T), linetype = "dashed",color = "white")
ggsave("RedBreakOut.png", width = 14, height =10, dpi = "retina")

#Biggest Plays ----
big_play <- pbp_rp %>% 
  filter(season == 2024) %>% 
  filter(week == 2) %>% 
  select(desc,wpa) %>% 
  mutate(wpa = abs(wpa))
arrange(-wpa)

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
       title = "Average Win Probability vs Time Remaining", subtitle = "Vertical Lines Denote End of Quarters",
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
  filter(rushes>=20) %>% 
  mutate(rusher_player_name = ifelse(rusher_player_name == "T.Hill", "Taysom Hill", rusher_player_name)) %>% 
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr"))
    



rushing_player %>%   
  ggplot(aes(x = success_rate, y = rush_epa))+
  geom_image(aes(image = team_logo_espn), size = 0.03, asp = 16/9)+
  labs(x = "Success Rate", y = "EPA/Rush", title = "Rushing Efficiency by Rusher Following Week 3", subtitle = "Minimum 20 Rushes",
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")+
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
        panel.border = element_rect(colour = "white", fill = NA, size = 1))
ggsave("RunningEfficiency.png", width = 14, height =10, dpi = "retina")


#Motion Usage----
pbp_rp %>% 
  group_by(posteam) %>% 
  summarize(motion_rate = mean(is_motion,na.rm = T), epa_motion = mean(epa[is_motion == 1],na.rm = T), epa_no_motion = mean(epa[is_motion == 0],na.rm = T)) %>% 
  mutate(epa_change = epa_motion - epa_no_motion) %>% 
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr")) %>% 
  ggplot(aes(x = motion_rate, y = epa_change))+
  geom_point()+
  geom_image(aes(image = team_logo_espn),size = 0.05, asp = 16/9)+
  labs(x = "Motion Rate", y = "Change in EPA/Play With Motion (EPA/Motion Play - EPA/No Motion Play)", title = "Which Teams Should Increase/Decrease Their Motion Usage?",
       subtitle = "Dotted Lines Represent League Average")+
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
  geom_hline(yintercept = mean(pbp_rp$epa[pbp_rp$is_motion == 1], na.rm = TRUE) - mean(pbp_rp$epa[pbp_rp$is_motion == 0], na.rm = TRUE), linetype = "dashed",color = "white")+
  geom_vline(xintercept = mean(pbp_rp$is_motion, na.rm = TRUE), linetype = "dashed", color = "white")
ggsave("MotionRate.png", width = 14, height =10, dpi = "retina")

#Play Action Rate---- 
pbp_rp %>% 
  filter(season == year) %>% 
  group_by(posteam) %>% 
  summarize(pa_rate = mean(is_play_action[pass == 1],na.rm = T), pa24 = mean(epa[is_play_action == 1],na.rm = T),
            no_pa = mean(epa[is_play_action == 0 & pass == 1], na.rm = T)) %>% 
  mutate(pa_improve = pa24 - no_pa) %>% 
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr")) %>% 
  ggplot(aes(x = pa_rate, y = pa_improve))+
  # geom_image(aes(image = team_logo_espn), size = 0.05, asp = 16/9)+
  geom_nfl_logos(aes(team_abbr = posteam), width = 0.065)+
  labs(x = "Play Action Rate", y = "EPA Improvement With Play Action (EPA/DB w/PA - EPA/DB no PA)", title = "Which Teams Should Increase/Decrease Play Action Usage?", caption ="@CapAnalytics7 | nflfastR",
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
        panel.border = element_rect(colour = "white", fill = NA, size = 1))+
  geom_hline(yintercept = mean(pbp_rp$epa[pbp_rp$is_play_action==1], na.rm = TRUE), linetype = "dashed",color = "white")+
  geom_vline(xintercept = mean(pbp_rp$is_play_action[pbp_rp$pass == 1], na.rm = TRUE), linetype = "dashed", color = "white")
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
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")+
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
  geom_hline(yintercept = mean(earlyvslate$late_down_epa), linetype = "dashed",color = "white")+
  geom_vline(xintercept = mean(earlyvslate$early_down_epa), linetype = "dashed", color = "white")
ggsave("LeadingvsTrailing.png", width = 14, height =10, dpi = "retina")


currentweek <- load_schedules(2023) %>% 
  filter(week ==12)

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
  geom_image(aes(image = team_logo_espn), size = 0.02, asp = 16/9)+
  labs(x = "% of Passing Yards from YAC", y = "EPA Per Dropback", title = "EPA Per Pass Attempt vs % of Passing Yards from YAC", subtitle = "Minimum 50 Pass Attempts",
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")+
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
  geom_text_repel(
    aes(label = passer_player_name),
    box.padding = 0.05,  # adjust this value for padding around the labels
    point.padding = 0.01,  # adjust this value for padding around the points
    segment.color = "grey",
    segment.size = 0.2,
  )+
  theme_bw()
ggsave("YAC.png", width = 14, height =10, dpi = "retina")


#Ability to recover from negative plays----

negative_plays <- pbp_rp %>% 
  mutate(
    negative_play = ifelse(yards_gained < 0 | penalty_team == posteam, 1, 0),
    unique_drive_identifier = paste(game_id, drive),
    drive_points = ifelse(fixed_drive_result == "Touchdown", 6,
                          ifelse(fixed_drive_result == "Field goal", 3, 0))) %>% #check drive_points
  filter(!grepl("punt", desc, ignore.case = TRUE)) %>% 
  group_by(unique_drive_identifier) %>% 
  summarize(posteam = first(posteam),negative_plays = sum(negative_play, na.rm = T),points_drive = max(drive_points)) %>% 
  group_by(negative_plays,posteam) %>% 
  summarize(avg_points_drive = mean(points_drive), count = n()) %>% 
  filter(!is.na(posteam)) %>% 
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr"))


negative_plays %>% 
  filter(count>2) %>% 
  ggplot(aes(x = negative_plays, y = avg_points_drive))+
  geom_image(aes(image = team_logo_espn), size = 0.02, asp = 16/9)+
  labs(x = "Number of Negative Plays on Drive (Includes Penalties)", y = "Average Points Per Drive",
       title = "Points Per Drive vs Negative Plays Following Week 17 TNF", subtitle = "Minimum 3 Drives",
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR" )+
  theme_bw()
ggsave("NegativePlays.png", width = 14, height =10, dpi = "retina")

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

#Biggest Plays ----
big_play <- pbp_rp %>% 
  filter(season == 2024) %>% 
  filter(week == 2) %>% 
  select(desc,wpa) %>% 
  mutate(wpa = abs(wpa))
  arrange(-wpa)


#Success by Air Yards
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

#Middle 8 ----