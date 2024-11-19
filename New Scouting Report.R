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
library(tidyr)
library(gt)
library(gtExtras)

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

replace_with_ranks<- function(column){
  values <- column
  ranks <- rank(column*-1,ties.method = "max")
}

#----Scouting Report (add blitz and fix ranking for yards and stuff, int rate)
offensive_scouting <- pbp_rp %>% 
  group_by(posteam) %>% 
  summarize(
    off_epa = mean(epa,na.rm = T),
    `EPA Last 4 Weeks` = mean(epa[week> max(week)-4],na.rm = T),
    `Success Rate` = mean(success, na.rm = T),
    off_pass_epa = mean(epa[pass == 1],na.rm = T),
    `Dropback Success Rate` = mean(success[pass==1],na.rm = T),
    off_rush_epa = mean(epa[pass == 0],na.rm = T),
    `Rush Success Rate` = mean(success[rush==1],na.rm = T),
    off_explosive = mean(explosive,na.rm = T),
    off_explosive_rush = mean(explosive[rush == 1],na.rm = T),
    off_explosive_pass = mean(explosive[pass == 1],na.rm = T),
    off_negative = mean(yards_gained<0,na.rm = T),
    off_negative_pass = mean(yards_gained[pass == 1]<0,na.rm = T),
    off_negative_rush = mean(yards_gained[rush == 1]<0,na.rm = T),
    turnover_rate = mean(interception+fumble_lost,na.rm = T),
    fumble_rate = mean(fumble, na.rm = T),
    fumble_lost_rate =  mean(fumble_lost, na.rm = T),
    off_early_down_epa = mean(epa[down %in% c(1,2)],na.rm = T),
    off_early_down_pass_epa = mean(epa[down %in% c(1,2) & pass ==1],na.rm =T),
    off_early_down_rush_epa = mean(epa[down %in% c(1,2) & rush ==1],na.rm =T),
    off_early_down_pass_rate = mean(pass[down == 1 | down == 2],na.rm =T),
    early_down_pass_oe = mean(pass_oe[down == 1 | down == 2],na.rm =T),
    off_early_down_1st_pct = mean(first_down[down == 1 |down == 2],na.rm =T),
    off_3rd_down_1st_pct = mean(first_down[down == 3],na.rm =T),
    off_third_down_dist = mean(ydstogo[down ==3],na.rm =T),
    off_third_down_epa = mean(epa[down == 3], na.rm =T),
    off_red_zone_td_pct = mean(touchdown[yardline_100<=20], na.rm =T),
    off_red_zone_td_epa = mean(epa[yardline_100<=20], na.rm =T),
    off_yac_pct = sum(yards_after_catch,na.rm =T)/sum(yards_gained[pass==1],na.rm =T),
    sack_rate = mean(sack[pass == 1], na.rm = T),
    off_scramble_rate = mean(qb_scramble[pass == 1],na.rm = T),
    off_scramble_epa = mean(epa[pass == 1 & qb_scramble == 1],na.rm = T),
    short_pass_rate = mean(short_throw[pass == 1],na.rm =T),
    medium_pass_rate = mean(medium_throw[pass == 1],na.rm =T),
    long_pass_rate = mean(long_throw[pass == 1],na.rm =T),
    off_short_pass_epa = mean(epa[short_throw==1],na.rm =T),
    off_medium_pass_epa = mean(epa[medium_throw ==1],na.rm =T),
    off_long_pass_epa = mean(epa[long_throw == 1],na.rm =T),
    pass_rate = mean(pass,na.rm = T),
    passoe = mean(pass_oe,na.rm = T),
    off_epa_xpass_sit = mean(epa[pass_oe>=0.9],na.rm = T),
    pass_right_rate = mean(pass_location == "right",na.rm = T),
    off_pass_right_epa = mean(epa[pass_location == "right"],na.rm = T),
    pass_left_rate = mean(pass_location == "left",na.rm = T),
    off_pass_left_epa = mean(epa[pass_location == "left"],na.rm = T),
    pass_middle_rate = mean(pass_location == "middle",na.rm = T),
    off_pass_middle_epa = mean(epa[pass_location == "middle"],na.rm = T),
    aDot = mean(air_yards, na.rm = TRUE),
    play_action_rate = mean(is_play_action[pass == 1],na.rm = T),
    off_play_action_epa = mean(epa[pass == 1 & is_play_action == TRUE],na.rm = T),
    off_no_play_action_epa = mean(epa[pass == 1 & is_play_action == FALSE],na.rm = T),
    off_leave_pocket_rate = mean(is_qb_out_of_pocket[pass == 1],na.rm = T),
    off_epa_in_pocket = mean(epa[pass == 1 & is_qb_out_of_pocket == 0],na.rm = T),
    off_epa_out_pocket  = mean(epa[pass == 1 & is_qb_out_of_pocket == 1],na.rm = T),
    interception_rate = mean(interception[pass == 1], na.rm = T),
    interception_worthy_rate = mean(is_interception_worthy[pass==1],na.rm = T),
    intereception_worthy_epa = mean(epa[pass==1&is_interception_worthy ==1],na.rm = T),
    screen_rate = mean(is_screen_pass[pass == 1],na.rm = T),
    off_screen_epa = mean(epa[is_screen_pass == 1],na.rm = T),
    rush_right_rate = mean(run_location == "right",na.rm = T),
    off_rush_right_epa = mean(epa[run_location == "right"],na.rm = T),
    rush_left_rate = mean(run_location == "left",na.rm = T),
    off_rush_left_epa = mean(epa[run_location == "left"],na.rm = T),
    rush_middle_rate = mean(run_location == "middle",na.rm = T),
    off_rush_middle_epa = mean(epa[run_location == "middle"],na.rm = T),
    off_shotgun_rate = mean(shotgun,na.rm = T),
    given_pass_pct_from_shotgun = mean(shotgun[pass==1],na.rm = T),
    given_rush_pct_from_shotgun = mean(shotgun[rush==1],na.rm = T),
    given_shotgun_pass_rate = mean(pass[shotgun==1],na.rm = T),
    given_under_center_pass_rate = mean(pass[shotgun==0],na.rm = T),
    off_shotgun_epa = mean(epa[shotgun == 1],na.rm = T),
    off_shotgun_pass_epa = mean(epa[shotgun == 1 & pass == 1],na.rm = T),
    off_shotgun_rush_epa = mean(epa[shotgun == 1 & rush == 1],na.rm = T),
    off_no_shotgun_epa = mean(epa[shotgun == 0],na.rm = T),
    off_no_shotgun_pass_epa = mean(epa[shotgun == 0 & pass == 1],na.rm = T),
    off_no_shotgun_rush_epa = mean(epa[shotgun == 0 & rush == 1],na.rm = T),
    off_1h_epa = mean(epa[qtr<=21],na.rm =T),
    off_1h_pass_epa = mean(epa[qtr<=2&pass ==1], na.rm =T),
    off_1h_rush_epa = mean(epa[qtr<=2&rush ==1], na.rm =T),
    off_2h_epa = mean(epa[qtr>2],na.rm =T),
    off_2h_pass_epa = mean(epa[qtr>2&pass ==1], na.rm =T),
    off_2h_rush_epa = mean(epa[qtr>2&rush ==1], na.rm =T),
    blitz_rate = mean(blitz,na.rm = T),
    off_no_blitz_epa = mean(epa[blitz == 0 & pass == 1],na.rm = T),
    off_blitz_epa = mean(epa[blitz == 1 & pass == 1],na.rm = T),
    off_blitz_scramble_epa = mean(epa[blitz == 1 & pass == 1 & qb_scramble == 1],na.rm = T),
    off_blitz_scramble_ypc = mean(yards_gained[blitz == 1 & pass == 1&qb_scramble == 1],na.rm = T),
    off_adot_blitz = mean(air_yards[blitz == 1 & pass ==1],na.rm =T),
    blitz_sens = off_blitz_epa - off_no_blitz_epa,
    motion_rate = mean(is_motion,na.rm = T),
    pass_motion_rate = mean(is_motion[pass == 1],na.rm = T),
    rush_motion_rate = mean(is_motion[rush == 1],na.rm = T),
    off_no_motion_epa = mean(epa[is_motion == FALSE],na.rm = T),
    off_motion_epa = mean(epa[is_motion == TRUE],na.rm = T),
    off_pass_no_motion_epa = mean(epa[is_motion == FALSE & pass ==1 ],na.rm = T),
    off_rush_no_motion_epa = mean(epa[is_motion == FALSE & rush == 1 ],na.rm = T),
    off_pass_motion_epa = mean(epa[is_motion == TRUE & pass ==1 ],na.rm = T),
    off_rush_motion_epa = mean(epa[is_motion == TRUE & rush == 1 ],na.rm = T),
    off_motion_impr = off_motion_epa-off_no_motion_epa,
    off_pass_motion_impr = off_pass_motion_epa - off_pass_no_motion_epa,
    off_rush_motion_impr = off_rush_motion_epa - off_rush_no_motion_epa,
    rpo_rate = mean(is_rpo,na.rm = T),
    off_rpo_epa = mean(epa[is_rpo == 1],na.rm = T),
    off_rush_in_heavy_box_rate = mean(rush[heavybox == 1],na.rm = T),
    off_rush_in_heavy_box_epa = mean(epa[heavybox == 1 & rush == 1],na.rm = T),
    off_pass_in_heavy_box_epa = mean(epa[heavybox == 1 & pass == 1],na.rm = T),
    off_rush_in_light_box_rate = mean(rush[lightbox == 1],na.rm = T),
    off_rush_in_light_box_epa = mean(epa[lightbox == 1 & rush == 1],na.rm = T),
    off_pass_in_light_box_epa = mean(epa[lightbox == 1 & pass == 1],na.rm = T),
    off_no_huddle_rate = mean(is_no_huddle,na.rm = T),
    off_no_huddle_epa = mean(epa[is_no_huddle == 1],na.rm = T),
    home_off_epa = mean(epa[posteam == home_team],na.rm = T),
    home_pass_off_epa = mean(epa[posteam == home_team&pass == 1],na.rm = T),
    home_rush_off_epa = mean(epa[posteam == home_team&pass == 0],na.rm = T),
    away_off_epa = mean(epa[posteam == away_team],na.rm = T),
    away_pass_off_epa = mean(epa[posteam == away_team&pass == 1],na.rm = T),
    away_rush_off_epa = mean(epa[posteam == away_team&pass == 0],na.rm = T)
  )

#Copy to defense finish checking rankings shit----
defensive_scouting <- pbp_rp %>% 
  group_by(defteam) %>% 
  summarize(
    off_epa = mean(epa,na.rm = T),
    `EPA Last 4 Weeks` = mean(epa[week> max(week)-4]),
    `Success Rate` = mean(success, na.rm = T),
    off_pass_epa = mean(epa[pass == 1],na.rm = T),
    `Dropback Success Rate` = mean(success[pass==1],na.rm = T),
    off_rush_epa = mean(epa[pass == 0],na.rm = T),
    `Rush Success Rate` = mean(success[rush==1],na.rm = T),
    off_explosive = mean(explosive,na.rm = T),
    off_explosive_rush = mean(explosive[rush == 1],na.rm = T),
    off_explosive_pass = mean(explosive[pass == 1],na.rm = T),
    off_negative = mean(yards_gained<0,na.rm = T),
    off_negative_pass = mean(yards_gained[pass == 1]<0,na.rm = T),
    off_negative_rush = mean(yards_gained[rush == 1]<0,na.rm = T),
    turnover_rate = mean(interception+fumble_lost,na.rm = T),
    fumble_rate = mean(fumble, na.rm = T),
    fumble_lost_rate =  mean(fumble_lost, na.rm = T),
    off_early_down_epa = mean(epa[down %in% c(1,2)],na.rm = T),
    off_early_down_pass_epa = mean(epa[down %in% c(1,2) & pass ==1],na.rm =T),
    off_early_down_rush_epa = mean(epa[down %in% c(1,2) & rush ==1],na.rm =T),
    off_early_down_pass_rate = mean(pass[down == 1 | down == 2],na.rm =T),
    early_down_pass_oe = mean(pass_oe[down == 1 | down == 2],na.rm =T),
    off_early_down_1st_pct = mean(first_down[down == 1 |down == 2],na.rm =T),
    off_3rd_down_1st_pct = mean(first_down[down == 3],na.rm =T),
    off_third_down_dist = mean(ydstogo[down ==3],na.rm =T),
    off_third_down_epa = mean(epa[down == 3], na.rm =T),
    off_red_zone_td_pct = mean(touchdown[yardline_100<=20], na.rm =T),
    off_red_zone_td_epa = mean(epa[yardline_100<=20], na.rm =T),
    off_yac_pct = sum(yards_after_catch,na.rm =T)/sum(yards_gained[pass==1],na.rm =T),
    sack_rate = mean(sack[pass == 1], na.rm = T),
    off_scramble_rate = mean(qb_scramble[pass == 1],na.rm = T),
    off_scramble_epa = mean(epa[pass == 1 & qb_scramble == 1],na.rm = T),
    short_pass_rate = mean(short_throw[pass == 1],na.rm =T),
    medium_pass_rate = mean(medium_throw[pass == 1],na.rm =T),
    long_pass_rate = mean(long_throw[pass == 1],na.rm =T),
    off_short_pass_epa = mean(epa[short_throw==1],na.rm =T),
    off_medium_pass_epa = mean(epa[medium_throw ==1],na.rm =T),
    off_long_pass_epa = mean(epa[long_throw == 1],na.rm =T),
    pass_rate = mean(pass,na.rm = T),
    passoe = mean(pass_oe,na.rm = T),
    off_epa_xpass_sit = mean(epa[pass_oe>=0.9],na.rm = T),
    pass_right_rate = mean(pass_location == "right",na.rm = T),
    off_pass_right_epa = mean(epa[pass_location == "right"],na.rm = T),
    pass_left_rate = mean(pass_location == "left",na.rm = T),
    off_pass_left_epa = mean(epa[pass_location == "left"],na.rm = T),
    pass_middle_rate = mean(pass_location == "middle",na.rm = T),
    off_pass_middle_epa = mean(epa[pass_location == "middle"],na.rm = T),
    aDot = mean(air_yards, na.rm = TRUE),
    play_action_rate = mean(is_play_action[pass == 1],na.rm = T),
    off_play_action_epa = mean(epa[pass == 1 & is_play_action == TRUE],na.rm = T),
    off_no_play_action_epa = mean(epa[pass == 1 & is_play_action == FALSE],na.rm = T),
    off_leave_pocket_rate = mean(is_qb_out_of_pocket[pass == 1],na.rm = T),
    off_epa_in_pocket = mean(epa[pass == 1 & is_qb_out_of_pocket == 0],na.rm = T),
    off_epa_out_pocket  = mean(epa[pass == 1 & is_qb_out_of_pocket == 1],na.rm = T),
    interception_rate = mean(interception[pass == 1], na.rm = T),
    interception_worthy_rate = mean(is_interception_worthy[pass==1],na.rm = T),
    intereception_worthy_epa = mean(epa[pass==1&is_interception_worthy ==1],na.rm = T),
    screen_rate = mean(is_screen_pass[pass == 1],na.rm = T),
    off_screen_epa = mean(epa[is_screen_pass == 1],na.rm = T),
    rush_right_rate = mean(run_location == "right",na.rm = T),
    off_rush_right_epa = mean(epa[run_location == "right"],na.rm = T),
    rush_left_rate = mean(run_location == "left",na.rm = T),
    off_rush_left_epa = mean(epa[run_location == "left"],na.rm = T),
    rush_middle_rate = mean(run_location == "middle",na.rm = T),
    off_rush_middle_epa = mean(epa[run_location == "middle"],na.rm = T),
    off_shotgun_rate = mean(shotgun,na.rm = T),
    given_pass_pct_from_shotgun = mean(shotgun[pass==1],na.rm = T),
    given_rush_pct_from_shotgun = mean(shotgun[rush==1],na.rm = T),
    given_shotgun_pass_rate = mean(pass[shotgun==1],na.rm = T),
    given_under_center_pass_rate = mean(pass[shotgun==0],na.rm = T),
    off_shotgun_epa = mean(epa[shotgun == 1],na.rm = T),
    off_shotgun_pass_epa = mean(epa[shotgun == 1 & pass == 1],na.rm = T),
    off_shotgun_rush_epa = mean(epa[shotgun == 1 & rush == 1],na.rm = T),
    off_no_shotgun_epa = mean(epa[shotgun == 0],na.rm = T),
    off_no_shotgun_pass_epa = mean(epa[shotgun == 0 & pass == 1],na.rm = T),
    off_no_shotgun_rush_epa = mean(epa[shotgun == 0 & rush == 1],na.rm = T),
    off_1h_epa = mean(epa[qtr<=21],na.rm =T),
    off_1h_pass_epa = mean(epa[qtr<=2&pass ==1], na.rm =T),
    off_1h_rush_epa = mean(epa[qtr<=2&rush ==1], na.rm =T),
    off_2h_epa = mean(epa[qtr>2],na.rm =T),
    off_2h_pass_epa = mean(epa[qtr>2&pass ==1], na.rm =T),
    off_2h_rush_epa = mean(epa[qtr>2&rush ==1], na.rm =T),
    blitz_rate = mean(blitz,na.rm = T),
    off_no_blitz_epa = mean(epa[blitz == 0 & pass == 1],na.rm = T),
    off_blitz_epa = mean(epa[blitz == 1 & pass == 1],na.rm = T),
    off_blitz_scramble_epa = mean(epa[blitz == 1 & pass == 1 & qb_scramble == 1],na.rm = T),
    off_blitz_scramble_ypc = mean(yards_gained[blitz == 1 & pass == 1&qb_scramble == 1],na.rm = T),
    off_adot_blitz = mean(air_yards[blitz == 1 & pass ==1],na.rm =T),
    blitz_sens = off_blitz_epa - off_no_blitz_epa,
    motion_rate = mean(is_motion,na.rm = T),
    pass_motion_rate = mean(is_motion[pass == 1],na.rm = T),
    rush_motion_rate = mean(is_motion[rush == 1],na.rm = T),
    off_no_motion_epa = mean(epa[is_motion == FALSE],na.rm = T),
    off_motion_epa = mean(epa[is_motion == TRUE],na.rm = T),
    off_pass_no_motion_epa = mean(epa[is_motion == FALSE & pass ==1 ],na.rm = T),
    off_rush_no_motion_epa = mean(epa[is_motion == FALSE & rush == 1 ],na.rm = T),
    off_pass_motion_epa = mean(epa[is_motion == TRUE & pass ==1 ],na.rm = T),
    off_rush_motion_epa = mean(epa[is_motion == TRUE & rush == 1 ],na.rm = T),
    off_motion_impr = off_motion_epa-off_no_motion_epa,
    off_pass_motion_impr = off_pass_motion_epa - off_pass_no_motion_epa,
    off_rush_motion_impr = off_rush_motion_epa - off_rush_no_motion_epa,
    rpo_rate = mean(is_rpo,na.rm = T),
    off_rpo_epa = mean(epa[is_rpo == 1],na.rm = T),
    off_rush_in_heavy_box_rate = mean(rush[heavybox == 1],na.rm = T),
    off_rush_in_heavy_box_epa = mean(epa[heavybox == 1 & rush == 1],na.rm = T),
    off_pass_in_heavy_box_epa = mean(epa[heavybox == 1 & pass == 1],na.rm = T),
    off_rush_in_light_box_rate = mean(rush[lightbox == 1],na.rm = T),
    off_rush_in_light_box_epa = mean(epa[lightbox == 1 & rush == 1],na.rm = T),
    off_pass_in_light_box_epa = mean(epa[lightbox == 1 & pass == 1],na.rm = T),
    off_no_huddle_rate = mean(is_no_huddle,na.rm = T),
    off_no_huddle_epa = mean(epa[is_no_huddle == 1],na.rm = T),
    home_off_epa = mean(epa[posteam == home_team],na.rm = T),
    home_pass_off_epa = mean(epa[posteam == home_team&pass == 1],na.rm = T),
    home_rush_off_epa = mean(epa[posteam == home_team&pass == 0],na.rm = T),
    away_off_epa = mean(epa[posteam == away_team],na.rm = T),
    away_pass_off_epa = mean(epa[posteam == away_team&pass == 1],na.rm = T),
    away_rush_off_epa = mean(epa[posteam == away_team&pass == 0],na.rm = T)
  )


offense <- "HOU"
defense <- "DAL"





team_stats_numeric_off <- offensive_scouting %>% 
  select(-posteam) %>% 
  mutate_all(~round(.,3))

team_stats_numeric_def <-defensive_scouting %>% 
  select(-defteam) %>% 
  mutate_all(~round(.,3))

data_ranks_off <- apply(team_stats_numeric_off %>% 
                          mutate(off_negative  = off_negative*-1, off_negative_pass  = off_negative_pass*-1, off_negative_rush = off_negative_rush*-1,off_third_down_dist = off_third_down_dist*-1,
                        sack_rate = sack_rate*-1, interception_worthy_rate = interception_worthy_rate*-1, interception_rate = interception_rate*-1,turnover_rate = turnover_rate*-1, fumble_rate = fumble_rate*-1, fumble_lost_rate = fumble_lost_rate*-1), 2, replace_with_ranks) %>% 
  as.data.frame(.) %>% 
  cbind(offensive_scouting$posteam,.) %>% 
  rename("posteam" = "offensive_scouting$posteam") %>% 
  pivot_longer(cols = -posteam,names_to = "Stat", values_to = "Offense Rank")

data_with_values_and_ranks_off <- apply(team_stats_numeric_off %>% 
                                          mutate(off_negative  = off_negative*-1, off_negative_pass  = off_negative_pass*-1, off_negative_rush = off_negative_rush*-1,off_third_down_dist = off_third_down_dist*-1,
                                                 sack_rate = sack_rate*-1, interception_worthy_rate = interception_worthy_rate*-1, interception_rate = interception_rate*-1,
                                                 turnover_rate = turnover_rate*-1, fumble_rate = fumble_rate*-1, fumble_lost_rate = fumble_lost_rate*-1)
                                        , 2, replace_with_values_and_ranks) %>%
  as.data.frame(.) %>% 
  cbind(offensive_scouting$posteam,.) %>% 
  rename("posteam" = "offensive_scouting$posteam") %>% 
  pivot_longer(cols = -posteam,names_to = "Stat", values_to = "Offense Value") %>% 
  left_join(data_ranks_off, by = c("posteam","Stat"))

data_ranks_def <- apply(team_stats_numeric_def %>% 
                          mutate_all(~.*-1) %>% 
                          mutate(off_negative  = off_negative*-1, off_negative_pass  = off_negative_pass*-1, off_negative_rush = off_negative_rush*-1,off_third_down_dist = off_third_down_dist*-1, 
                                 turnover_rate = turnover_rate*-1, fumble_rate = fumble_rate*-1, fumble_lost_rate = fumble_lost_rate*-1, sack_rate = sack_rate*-1, 
                                 interception_rate =  interception_rate * -1, interception_worthy_rate = interception_worthy_rate*-1, blitz_rate = blitz_rate*-1,
                                 ), 2, replace_with_ranks) %>% 
  as.data.frame(.) %>% 
  cbind(defensive_scouting$defteam,.) %>% 
  rename("defteam" = "defensive_scouting$defteam") %>% 
  pivot_longer(cols = -defteam,names_to = "Stat", values_to = "Defense Rank")


data_with_values_and_ranks_def <- apply(team_stats_numeric_def %>% 
                                          mutate_all(~.*-1) %>% 
                                          mutate(off_negative  = off_negative*-1, off_negative_pass  = off_negative_pass*-1, off_negative_rush = off_negative_rush*-1,off_third_down_dist = off_third_down_dist*-1, 
                                                 turnover_rate = turnover_rate*-1, fumble_rate = fumble_rate*-1, fumble_lost_rate = fumble_lost_rate*-1, sack_rate = sack_rate*-1, 
                                                 interception_rate =  interception_rate * -1, interception_worthy_rate = interception_worthy_rate*-1, blitz_rate = blitz_rate*-1,
                                                 ), 2, replace_with_values_and_ranks) %>% 
  as.data.frame(.) %>% 
  cbind(defensive_scouting$defteam,.) %>% 
  rename("defteam" = "defensive_scouting$defteam") %>% 
  pivot_longer(cols = -defteam,names_to = "Stat", values_to = "Defense Value") %>% 
  left_join(data_ranks_def, by = c("defteam","Stat"))


scouting_report <- data_with_values_and_ranks_off %>% 
  filter(posteam == offense) %>% 
  left_join(data_with_values_and_ranks_def %>% filter(defteam == defense),
            by = c("Stat")) %>% 
  select(-posteam,-defteam)

scouting_report %>% 
  select(Stat,`Offense Value`,`Defense Value`,`Offense Rank`,`Defense Rank`) %>% 
  mutate(rank_diff = `Offense Rank`-`Defense Rank`) %>% 
  gt() %>% 
  cols_align(align = "center") %>%
  cols_label(`Stat` = "Stat",`Offense Value` = paste(offense,"Offense Value", sep = " "), `Defense Value` = paste(defense,"Defense Value", sep = " "),
             `Offense Rank` = paste(offense,"Offense Rank", sep = " "), `Defense Rank` = paste(defense,"Defense Rank", sep = " "),
             rank_diff = "Rank Difference") %>% 
  gtExtras::gt_theme_538() %>% 
  gt_hulk_col_numeric(columns = c(`Offense Rank`,`Defense Rank`,rank_diff)) %>% 
  tab_header(title = md("Scouting Report"), subtitle = md("Purple Represents Offense Advantage, Green Represents Defense Advantage"))


