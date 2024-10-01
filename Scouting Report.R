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
library(gridExtra)
library(patchwork)
library(zoo)

# add explosive and negative play rate ----
pfr_stats_pass <- nflreadr::load_pfr_advstats(seasons = 2023,
                                              stat_type = "pass",
                                              summary_level = "week")
ngs_data_passing <- nflreadr::load_nextgen_stats(seasons = 2023,
                                                 stat_type = "passing")

ftn_data <- nflreadr::load_ftn_charting(2024) %>%
  select(-week, -season)
pbp <- load_pbp(2024)
pbp <- pbp %>%
  left_join(ftn_data, by = c("game_id" = "nflverse_game_id",
                             "play_id" = "nflverse_play_id"))


pbp_rp <- pbp %>%
  filter(pass == 1 | rush == 1) %>%
  filter(!is.na(epa)) %>% 
  mutate(blitz = ifelse(n_pass_rushers<= 4 | is.na(n_pass_rushers),0,1),
         lightbox = ifelse(n_defense_box<=6,1,0),
         heavybox = ifelse(n_defense_box>=8,1,0),
         short_throw = ifelse(air_yards<=10,1,0),
         medium_throw = ifelse(air_yards>10&air_yards<=20,1,0),
         long_throw = ifelse(air_yards>20,1,0)) %>% 
  mutate(explosive = ifelse((yards_gained>20 & pass_attempt == 1) | (yards_gained >12 & (qb_scramble == 1 | rush == 1)),1,0),
         negative = ifelse(yards_gained < 0, 1,0))


replace_with_values_and_ranks <- function(column) {
  values <- column
  ranks <- rank(column*-1,ties.method = "max")
  # reversed_ranks <- max(ranks) + 1 - ranks
  formatted <- paste(values, "[", ranks, "]", sep = "")
  formatted
}

def_replace_with_values_and_ranks <- function(column) {
  values <- column
  ranks <- 32-rank(column*-1,ties.method = "max")+1
  formatted <- paste(values, "[", ranks, "]", sep = "")
  formatted
}

#----Scouting Report
offensive_scouting <- pbp_rp %>% 
  # filter(week == 22,play_id !=150) %>% 
  group_by(posteam) %>% 
  summarize(
    off_epa = mean(epa),
    off_pass_epa = mean(epa[pass == 1]),
    off_rush_epa = mean(epa[pass == 0]),
    off_explosive = mean(explosive,na.rm = T),
    off_explosive_rush = mean(explosive[rush == 1],na.rm = T),
    off_explosive_pass = mean(explosive[pass == 1],na.rm = T),
    off_negative = mean(yards_gained<0,na.rm = T),
    off_early_down_epa = mean(epa[down %in% c(1,2)]),
    off_early_down_pass_epa = mean(epa[down %in% c(1,2) & pass ==1],na.rm =T),
    off_early_down_rush_epa = mean(epa[down %in% c(1,2) & rush ==1],na.rm =T),
    off_early_down_pass_rate = mean(pass[down == 1 | down == 2],na.rm =T),
    off_early_down_1st_pct = mean(first_down[down == 1 |down == 2],na.rm =T),
    off_3rd_down_1st_pct = mean(first_down[down == 3],na.rm =T),
    off_third_down_dist = mean(ydstogo[down ==3],na.rm =T),
    off_third_down_epa = mean(epa[down == 3], na.rm =T),
    off_red_zone_td_pct = mean(touchdown[yardline_100<=20], na.rm =T),
    off_red_zone_td_epa = mean(epa[yardline_100<=20], na.rm =T),
    off_yac_pct = sum(yards_after_catch,na.rm =T)/sum(yards_gained[pass==1],na.rm =T),
    off_scramble_rate = mean(qb_scramble[pass == 1]),
    off_scramble_epa = mean(epa[pass == 1 & qb_scramble == 1]),
    short_pass_rate = mean(short_throw[pass == 1],na.rm =T),
    medium_pass_rate = mean(medium_throw[pass == 1],na.rm =T),
    long_pass_rate = mean(long_throw[pass == 1],na.rm =T),
    off_short_pass_epa = mean(epa[short_throw==1],na.rm =T),
    off_medium_pass_epa = mean(epa[medium_throw ==1],na.rm =T),
    off_long_pass_epa = mean(epa[long_throw == 1],na.rm =T),
    pass_rate = mean(pass),
    passoe = mean(pass_oe,na.rm = T),
    off_epa_xpass_sit = mean(epa[pass_oe>=0.9],na.rm = T),
    pass_right_rate = mean(pass_location == "right",na.rm = T),
    off_pass_right_epa = mean(epa[pass_location == "right"],na.rm = T),
    pass_left_rate = mean(pass_location == "left",na.rm = T),
    off_pass_left_epa = mean(epa[pass_location == "left"],na.rm = T),
    pass_middle_rate = mean(pass_location == "middle",na.rm = T),
    off_pass_middle_epa = mean(epa[pass_location == "middle"],na.rm = T),
    aDot = mean(air_yards, na.rm = TRUE),
    off_1h_epa = mean(epa[qtr<=21],na.rm =T),
    off_1h_pass_epa = mean(epa[qtr<=2&pass ==1], na.rm =T),
    off_1h_rush_epa = mean(epa[qtr<=2&rush ==1], na.rm =T),
    off_2h_epa = mean(epa[qtr>2],na.rm =T),
    off_2h_pass_epa = mean(epa[qtr>2&pass ==1], na.rm =T),
    off_2h_rush_epa = mean(epa[qtr>2&rush ==1], na.rm =T),
    off_no_blitz_epa = mean(epa[blitz == 0 & pass == 1],na.rm = T),
    off_blitz_epa = mean(epa[blitz == 1 & pass == 1],na.rm = T),
    off_blitz_scramble_epa = mean(epa[blitz == 1 & pass == 1 & qb_scramble == 1],na.rm = T),
    off_blitz_scramble_ypc = mean(yards_gained[blitz == 1 & pass == 1&qb_scramble == 1],na.rm = T),
    off_adot_blitz = mean(air_yards[blitz == 1 & pass ==1],na.rm =T),
    blitz_sens = off_blitz_epa - off_no_blitz_epa,
    motion_rate = mean(is_motion),
    pass_motion_rate = mean(is_motion[pass == 1]),
    rush_motion_rate = mean(is_motion[rush == 1]),
    off_no_motion_epa = mean(epa[is_motion == FALSE]),
    off_pass_no_motion_epa = mean(epa[is_motion == FALSE & pass ==1 ]),
    off_rush_no_motion_epa = mean(epa[is_motion == FALSE & rush == 1 ]),
    off_motion_epa = mean(epa[is_motion == TRUE]),
    off_pass_motion_epa = mean(epa[is_motion == TRUE & pass ==1 ]),
    off_rush_motion_epa = mean(epa[is_motion == TRUE & rush == 1 ]),
    off_motion_impr = off_motion_epa-off_no_motion_epa,
    off_pass_motion_impr = off_pass_motion_epa - off_pass_no_motion_epa,
    off_rush_motion_impr = off_rush_motion_epa - off_rush_no_motion_epa,
    rush_right_rate = mean(run_location == "right",na.rm = T),
    off_rush_right_epa = mean(epa[run_location == "right"],na.rm = T),
    rush_left_rate = mean(run_location == "left",na.rm = T),
    off_rush_left_epa = mean(epa[run_location == "left"],na.rm = T),
    rush_middle_rate = mean(run_location == "middle",na.rm = T),
    off_rush_middle_epa = mean(epa[run_location == "middle"],na.rm = T),
    play_action_rate = mean(is_play_action[pass == 1]),
    off_play_action_epa = mean(epa[pass == 1 & is_play_action == TRUE]),
    off_no_play_action_epa = mean(epa[pass == 1 & is_play_action == FALSE]),
    off_shotgun_rate = mean(shotgun),
    given_pass_pct_from_shotgun = mean(shotgun[pass==1]),
    given_rush_pct_from_shotgun = mean(shotgun[rush==1]),
    given_shotgun_pass_rate = mean(pass[shotgun==1]),
    given_under_center_pass_rate = mean(pass[shotgun==0]),
    off_shotgun_epa = mean(epa[shotgun == 1]),
    off_shotgun_pass_epa = mean(epa[shotgun == 1 & pass == 1]),
    off_shotgun_rush_epa = mean(epa[shotgun == 1 & rush == 1]),
    off_no_shotgun_epa = mean(epa[shotgun == 0]),
    off_no_shotgun_pass_epa = mean(epa[shotgun == 0 & pass == 1]),
    off_no_shotgun_rush_epa = mean(epa[shotgun == 0 & rush == 1]),
    off_leave_pocket_rate = mean(is_qb_out_of_pocket[pass == 1]),
    off_epa_in_pocket = mean(epa[pass == 1 & is_qb_out_of_pocket == 0]),
    off_epa_out_pocket  = mean(epa[pass == 1 & is_qb_out_of_pocket == 1]),
    intereception_worthy_rate = mean(is_interception_worthy[pass==1]),
    intereception_worthy_epa = mean(epa[pass==1&is_interception_worthy ==1]),
    rpo_rate = mean(is_rpo),
    off_rpo_epa = mean(epa[is_rpo == 1]),
    screen_rate = mean(is_screen_pass[pass == 1]),
    off_screen_epa = mean(epa[is_screen_pass == 1]),
    off_rush_in_heavy_box_rate = mean(rush[heavybox == 1]),
    off_rush_in_heavy_box_epa = mean(epa[heavybox == 1 & rush == 1]),
    off_pass_in_heavy_box_epa = mean(epa[heavybox == 1 & pass == 1]),
    off_rush_in_light_box_rate = mean(rush[lightbox == 1]),
    off_rush_in_light_box_epa = mean(epa[lightbox == 1 & rush == 1]),
    off_pass_in_light_box_epa = mean(epa[lightbox == 1 & pass == 1]),
    off_no_huddle_rate = mean(is_no_huddle),
    off_no_huddle_epa = mean(epa[is_no_huddle == 1]),
    home_off_epa = mean(epa[posteam == home_team]),
    home_pass_off_epa = mean(epa[posteam == home_team&pass == 1]),
    home_rush_off_epa = mean(epa[posteam == home_team&pass == 0]),
    away_off_epa = mean(epa[posteam == away_team]),
    away_pass_off_epa = mean(epa[posteam == away_team&pass == 1]),
    away_rush_off_epa = mean(epa[posteam == away_team&pass == 0])
  )

defensive_scouting <- pbp_rp %>%
  group_by(defteam) %>% 
  summarize(
    def_epa = mean(epa),
    def_rush_epa = mean(epa[pass == 0]),
    def_pass_epa = mean(epa[pass == 1]),
    def_epa_xpass_sit = mean(epa[pass_oe>=0.9],na.rm = T),
    def_explosive = mean(explosive,na.rm = T),
    def_explosive_rush = mean(explosive[rush == 1],na.rm = T),
    def_explosive_pass = mean(explosive[pass == 1],na.rm = T),
    def_negative = mean(yards_gained<0,na.rm = T),
    def_early_down_epa = mean(epa[down %in% c(1,2)]),
    def_early_down_pass_epa = mean(epa[down %in% c(1,2) & pass ==1],na.rm =T),
    def_early_down_rush_epa = mean(epa[down %in% c(1,2) & rush ==1],na.rm =T),
    def_early_down_1st_pct = mean(first_down[down == 1 |down == 2],na.rm =T),
    def_3rd_down_1st_pct = mean(first_down[down == 3],na.rm =T),
    def_third_down_dist = mean(ydstogo[down ==3],na.rm =T),
    def_third_down_epa = mean(epa[down == 3], na.rm =T),
    def_red_zone_td_pct = mean(touchdown[yardline_100<=20], na.rm =T),
    def_red_zone_td_epa = mean(epa[yardline_100<=20], na.rm =T),
    def_short_pass_epa = mean(epa[short_throw == 1],na.rm =T),
    def_medium_pass_epa = mean(epa[medium_throw == 1],na.rm =T),
    def_long_pass_epa = mean(epa[long_throw == 1],na.rm =T),
    def_yac_pct = sum(yards_after_catch,na.rm =T)/sum(yards_gained[pass==1],na.rm =T),
    def_1h_epa = mean(epa[qtr<=2],na.rm =T),
    def_1h_pass_epa = mean(epa[qtr<=2&pass ==1], na.rm =T),
    def_1h_rush_epa = mean(epa[qtr<=2&rush ==1], na.rm =T),
    def_2h_epa = mean(epa[qtr>2],na.rm =T),
    def_2h_pass_epa = mean(epa[qtr>2&pass ==1], na.rm =T),
    def_2h_rush_epa = mean(epa[qtr>2&rush ==1], na.rm =T),
    def_scrabmle_epa = mean(epa[pass == 1 & qb_scramble == 1]),
    def_blitz_rate = mean(blitz),
    def_no_blitz_epa = mean(epa[blitz == 0 & pass == 1]),
    def_blitz_epa = mean(epa[blitz == 1 & pass == 1],na.rm = T),
    def_blitz_scramble_epa = mean(epa[blitz == 1 & pass == 1&qb_scramble == 1],na.rm = T),
    def_blitz_scramble_ypc = mean(yards_gained[blitz == 1 & pass == 1&qb_scramble == 1],na.rm = T),
    def_no_motion_epa = mean(epa[is_motion == FALSE]),
    def_pass_no_motion_epa = mean(epa[is_motion == FALSE & pass ==1 ]),
    def_rush_no_motion_epa = mean(epa[is_motion == FALSE & rush == 1 ]),
    def_motion_epa = mean(epa[is_motion == TRUE]),
    def_pass_motion_epa = mean(epa[is_motion == TRUE & pass ==1 ]),
    def_rush_motion_epa = mean(epa[is_motion == TRUE & rush == 1 ]),
    def_motion_sens = def_motion_epa-def_no_motion_epa,
    def_pass_motion_sens = def_pass_motion_epa - def_pass_no_motion_epa,
    def_rush_motion_sens = def_rush_motion_epa - def_rush_no_motion_epa,
    def_play_action_epa = mean(epa[pass == 1 & is_play_action == TRUE]),
    def_no_play_action_epa = mean(epa[pass == 1 & is_play_action == FALSE]),
    def_shotgun_epa = mean(epa[shotgun == 1]),
    def_shotgun_pass_epa = mean(epa[shotgun == 1 & pass == 1]),
    def_shotgun_rush_epa = mean(epa[shotgun == 1 & rush == 1]),
    def_no_shotgun_epa = mean(epa[shotgun == 0]),
    def_no_shotgun_pass_epa = mean(epa[shotgun == 0 & pass == 1]),
    def_no_shotgun_rush_epa = mean(epa[shotgun == 0 & rush == 1]),
    def_pass_right_epa = mean(epa[pass_location == "right"],na.rm = T),
    def_pass_left_epa = mean(epa[pass_location == "left"],na.rm = T),
    def_pass_middle_epa = mean(epa[pass_location == "middle"],na.rm = T),
    def_rush_right_epa = mean(epa[run_location == "right"],na.rm = T),
    def_rush_left_epa = mean(epa[run_location == "left"],na.rm = T),
    def_rush_middle_epa = mean(epa[run_location == "middle"],na.rm = T),
    def_leave_pocket_rate = mean(is_qb_out_of_pocket[pass == 1]),
    def_epa_in_pocket = mean(epa[pass == 1 & is_qb_out_of_pocket == 0]),
    def_epa_out_pocket  = mean(epa[pass == 1 & is_qb_out_of_pocket == 1]),
    def_rpo_epa = mean(epa[is_rpo == 1]),
    def_screen_epa = mean(epa[is_screen_pass == 1]),
    def_heavy_box_rate = mean(heavybox),
    def_rush_in_heavy_box_epa = mean(epa[heavybox == 1 & rush == 1]),
    def_pass_in_heavy_box_epa = mean(epa[heavybox == 1 & pass == 1]),
    def_light_box_rate = mean(lightbox),
    def_rush_in_light_box_epa = mean(epa[lightbox == 1 & rush == 1]),
    def_pass_in_light_box_epa = mean(epa[lightbox == 1 & pass == 1]),
    def_no_huddle_epa = mean(epa[is_no_huddle == 1]),
    home_def_epa = mean(epa[defteam == home_team]),
    home_pass_def_epa = mean(epa[defteam == home_team&pass == 1]),
    home_rush_def_epa = mean(epa[defteam == home_team&pass == 0]),
    away_def_epa = mean(epa[defteam == away_team]),
    away_pass_def_epa = mean(epa[defteam == away_team&pass == 1]),
    away_rush_def_epa = mean(epa[defteam == away_team&pass == 0])
  )

scouting<- offensive_scouting %>% 
  left_join(defensive_scouting, by = c("posteam" = "defteam"))


team_stats_numeric <-
  scouting %>% 
  select(-posteam) %>% 
  mutate_all(~round(.,3))

data_with_values_and_ranks <- apply(team_stats_numeric, 2, replace_with_values_and_ranks)

data_with_ranks <- cbind(scouting$posteam,data.frame(data_with_values_and_ranks)) %>%
  rename("team" = "scouting$posteam")

teams <- c("DAL","NYG")

scout_post <- data_with_ranks %>%
  filter(team %in% teams) %>%
  t() %>%
  as.data.frame()


# mahomes_pass_map<-pbp_rp %>%
#   mutate(air_yards_bins = cut(air_yards,
#                               breaks = c(-Inf, 0, 10, 20, Inf),
#                               labels = c("<=0", "1-10", "11-20", "21+"))) %>%
#   filter(posteam == "KC") %>% 
#   filter(passer_player_name %in% c("P.Mahomes") & !is.na(pass_location)) %>%
#   group_by(passer_player_name, pass_location, air_yards_bins) %>%
#   summarize(count_air = n(), epa_pass = mean(epa), success_rate = mean(success),
#             YPA = mean(yards_gained),
#             cmp = sum(complete_pass) / n()) %>%
#   ungroup() %>%
#   group_by(passer_player_name) %>%
#   mutate(pct_throws = count_air / sum(count_air)) %>%
#   ggplot(aes(x = pass_location, y = air_yards_bins, fill = epa_pass)) +
#   geom_tile(show.legend = FALSE)+
#   geom_text(aes(label = paste("Throw Rate:", round(pct_throws, 2), "\nEPA:", round(epa_pass, 2),
#                               "\nSuccess:", round(success_rate, 2), "\nYPA:", round(YPA, 2),
#                               "\nComp %:", round(cmp, 2))),
#             size = 3) +
#   scale_fill_gradient(low = "red", high = "green", limits = c(-0.5,0.7), na.value = "green") +
#   labs(x = "Pass Location", y = "Air Yards", title = "Mahomes Pass Map") +
#   theme(legend.position = "none")+
#   theme_minimal()
# 
# sf_pass_map<- pbp_rp %>%
#   mutate(air_yards_bins = cut(air_yards,
#                               breaks = c(-Inf, 0, 10, 20, Inf),
#                               labels = c("<=0", "1-10", "11-20", "21+"))) %>%
#   filter(defteam == "SF", !is.na(pass_location)) %>%
#   group_by(defteam, pass_location, air_yards_bins) %>%
#   summarize(count_air = n(), epa_pass = mean(epa), success_rate = mean(success),
#             YPA = mean(yards_gained),
#             cmp = sum(complete_pass) / n()) %>%
#   ungroup() %>%
#   group_by(defteam) %>%
#   mutate(pct_throws = count_air / sum(count_air)) %>%
#   ggplot(aes(x = pass_location, y = air_yards_bins, fill = epa_pass)) +
#   geom_tile() +
#   geom_text(aes(label = paste("Throw Rate:", round(pct_throws, 2), "\nEPA:", round(epa_pass, 2),
#                               "\nSuccess:", round(success_rate, 2), "\nYPA:", round(YPA, 2),
#                               "\nComp %:", round(cmp, 2))),
#             size = 3) +
#   scale_fill_gradient(low = "red", high = "green", limits = c(-0.5,0.5), na.value = "red") +
#   labs(x = "Pass Location", y = NULL, title = "49ers Defense Pass Map") +
#   theme_minimal()
# 
# # Assuming purdy_pass_map and chiefs_pass_map are ggplot objects
# 
# # Arrange the plots
# combined_plot <- mahomes_pass_map + sf_pass_map
# combined_plot
# # Save the combined plot
# ggsave("combined_plot.png", combined_plot, width = 10, height = 5)
# 
# pbp_rp %>% 
#   filter(defteam == "SF") %>% 
#   group_by(week) %>%
#     summarize(pass = mean(epa[pass ==1]), rush = mean(epa[rush == 1])) %>% 
#   pivot_longer(c(pass,rush),names_to = "play_type",values_to = "epa_play") %>% 
#   ggplot(aes(x = week,y = epa_play))+
#   geom_point(aes(color = play_type))+
#   geom_line(aes(color = play_type))+
#   labs(x = "Week Number", y = "EPA/Play", title = "49ers Defense EPA/Play by Week and Play Type",
#        caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")
# ggsave("49ersbyWeek.png", width = 14, height =10, dpi = "retina")


test<-pbp_rp %>% 
  filter(week == 22,play_id !=150,posteam == "KC") %>% 
  group_by(game_half) %>% 
  summarize(
    off_epa = mean(epa),
    off_pass_epa = mean(epa[pass == 1]),
    off_rush_epa = mean(epa[pass == 0]),
    off_early_down_epa = mean(epa[down %in% c(1,2)]),
    off_early_down_pass_epa = mean(epa[down %in% c(1,2) & pass ==1],na.rm =T),
    off_early_down_rush_epa = mean(epa[down %in% c(1,2) & rush ==1],na.rm =T),
    off_early_down_pass_rate = mean(pass[down == 1 | down == 2],na.rm =T),
    off_early_down_1st_pct = mean(first_down[down == 1 |down == 2],na.rm =T),
    off_3rd_down_1st_pct = mean(first_down[down == 3],na.rm =T),
    off_third_down_dist = mean(ydstogo[down ==3],na.rm =T),
    off_third_down_epa = mean(epa[down == 3], na.rm =T),
    off_red_zone_td_pct = mean(touchdown[yardline_100<=20], na.rm =T),
    off_red_zone_td_epa = mean(epa[yardline_100<=20], na.rm =T),
    off_yac_pct = sum(yards_after_catch,na.rm =T)/sum(yards_gained[pass==1],na.rm =T),
    off_1h_epa = mean(epa[qtr<=2],na.rm =T),
    off_1h_pass_epa = mean(epa[qtr<=2&pass ==1], na.rm =T),
    off_1h_rush_epa = mean(epa[qtr<=2&rush ==1], na.rm =T),
    off_2h_epa = mean(epa[qtr>2],na.rm =T),
    off_2h_pass_epa = mean(epa[qtr>2&pass ==1], na.rm =T),
    off_2h_rush_epa = mean(epa[qtr>2&rush ==1], na.rm =T),
    off_scramble_rate = mean(qb_scramble[pass == 1]),
    off_scramble_epa = mean(epa[pass == 1 & qb_scramble == 1]),
    off_no_blitz_epa = mean(epa[blitz == 0 & pass == 1],na.rm = T),
    off_blitz_epa = mean(epa[blitz == 1 & pass == 1],na.rm = T),
    off_blitz_scramble_epa = mean(epa[blitz == 1 & pass == 1 & qb_scramble == 1],na.rm = T),
    off_blitz_scramble_ypc = mean(yards_gained[blitz == 1 & pass == 1&qb_scramble == 1],na.rm = T),
    off_adot_blitz = mean(air_yards[blitz == 1 & pass ==1],na.rm =T),
    blitz_sens = off_blitz_epa - off_no_blitz_epa,
    motion_rate = mean(is_motion),
    pass_motion_rate = mean(is_motion[pass == 1]),
    short_pass_rate = mean(short_throw[pass == 1],na.rm =T),
    medium_pass_rate = mean(medium_throw[pass == 1],na.rm =T),
    long_pass_rate = mean(long_throw[pass == 1],na.rm =T),
    off_short_pass_epa = mean(epa[short_throw==1],na.rm =T),
    off_medium_pass_epa = mean(epa[medium_throw ==1],na.rm =T),
    off_long_pass_epa = mean(epa[long_throw == 1],na.rm =T),
    rush_motion_rate = mean(is_motion[rush == 1]),
    pass_rate = mean(pass),
    pass_right_rate = mean(pass_location == "right",na.rm = T),
    off_pass_right_epa = mean(epa[pass_location == "right"],na.rm = T),
    pass_left_rate = mean(pass_location == "left",na.rm = T),
    off_pass_left_epa = mean(epa[pass_location == "left"],na.rm = T),
    pass_middle_rate = mean(pass_location == "middle",na.rm = T),
    off_pass_middle_epa = mean(epa[pass_location == "middle"],na.rm = T),
    aDot = mean(air_yards, na.rm = TRUE),
    rush_right_rate = mean(run_location == "right",na.rm = T),
    off_rush_right_epa = mean(epa[run_location == "right"],na.rm = T),
    rush_left_rate = mean(run_location == "left",na.rm = T),
    off_rush_left_epa = mean(epa[run_location == "left"],na.rm = T),
    rush_middle_rate = mean(run_location == "middle",na.rm = T),
    off_rush_middle_epa = mean(epa[run_location == "middle"],na.rm = T),
    off_no_motion_epa = mean(epa[is_motion == FALSE]),
    off_pass_no_motion_epa = mean(epa[is_motion == FALSE & pass ==1 ]),
    off_rush_no_motion_epa = mean(epa[is_motion == FALSE & rush == 1 ]),
    off_motion_epa = mean(epa[is_motion == TRUE]),
    off_pass_motion_epa = mean(epa[is_motion == TRUE & pass ==1 ]),
    off_rush_motion_epa = mean(epa[is_motion == TRUE & rush == 1 ]),
    off_motion_impr = off_motion_epa-off_no_motion_epa,
    off_pass_motion_impr = off_pass_motion_epa - off_pass_no_motion_epa,
    off_rush_motion_impr = off_rush_motion_epa - off_rush_no_motion_epa,
    play_action_rate = mean(is_play_action[pass == 1]),
    off_play_action_epa = mean(epa[pass == 1 & is_play_action == TRUE]),
    off_no_play_action_epa = mean(epa[pass == 1 & is_play_action == FALSE]),
    off_shotgun_rate = mean(shotgun),
    given_pass_pct_from_shotgun = mean(shotgun[pass==1]),
    given_rush_pct_from_shotgun = mean(shotgun[rush==1]),
    given_shotgun_pass_rate = mean(pass[shotgun==1]),
    given_under_center_pass_rate = mean(pass[shotgun==0]),
    off_shotgun_epa = mean(epa[shotgun == 1]),
    off_shotgun_pass_epa = mean(epa[shotgun == 1 & pass == 1]),
    off_shotgun_rush_epa = mean(epa[shotgun == 1 & rush == 1]),
    off_no_shotgun_epa = mean(epa[shotgun == 0]),
    off_no_shotgun_pass_epa = mean(epa[shotgun == 0 & pass == 1]),
    off_no_shotgun_rush_epa = mean(epa[shotgun == 0 & rush == 1]),
    off_leave_pocket_rate = mean(is_qb_out_of_pocket[pass == 1]),
    off_epa_in_pocket = mean(epa[pass == 1 & is_qb_out_of_pocket == 0]),
    off_epa_out_pocket  = mean(epa[pass == 1 & is_qb_out_of_pocket == 1]),
    intereception_worthy_rate = mean(is_interception_worthy[pass==1]),
    intereception_worthy_epa = mean(epa[pass==1&is_interception_worthy ==1]),
    rpo_rate = mean(is_rpo),
    off_rpo_epa = mean(epa[is_rpo == 1]),
    screen_rate = mean(is_screen_pass[pass == 1]),
    off_screen_epa = mean(epa[is_screen_pass == 1]),
    off_rush_in_heavy_box_rate = mean(rush[heavybox == 1]),
    off_rush_in_heavy_box_epa = mean(epa[heavybox == 1 & rush == 1]),
    off_pass_in_heavy_box_epa = mean(epa[heavybox == 1 & pass == 1]),
    off_rush_in_light_box_rate = mean(rush[lightbox == 1]),
    off_rush_in_light_box_epa = mean(epa[lightbox == 1 & rush == 1]),
    off_pass_in_light_box_epa = mean(epa[lightbox == 1 & pass == 1]),
    off_no_huddle_rate = mean(is_no_huddle),
    off_no_huddle_epa = mean(epa[is_no_huddle == 1]),
    home_off_epa = mean(epa[posteam == home_team]),
    home_pass_off_epa = mean(epa[posteam == home_team&pass == 1]),
    home_rush_off_epa = mean(epa[posteam == home_team&pass == 0]),
    away_off_epa = mean(epa[posteam == away_team]),
    away_pass_off_epa = mean(epa[posteam == away_team&pass == 1]),
    away_rush_off_epa = mean(epa[posteam == away_team&pass == 0])
  ) %>% 
  t()
