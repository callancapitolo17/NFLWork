library(dplyr)
library(tidyverse)
library(nflfastR)
pbp <- load_pbp(c(2006:2023)) %>% 
  filter(pass == 1 | rush == 1) %>%
  filter(!is.na(epa),!is.na(posteam),!is.na(defteam))%>% 
  mutate(short_throw = ifelse(air_yards<=10,1,0),
         medium_throw = ifelse(air_yards>10&air_yards<=20,1,0),
         long_throw = ifelse(air_yards>20,1,0))

offensive_scouting <- pbp %>% 
  group_by(posteam,season) %>% 
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
    off_1h_epa = mean(epa[qtr<=21],na.rm =T),
    off_1h_pass_epa = mean(epa[qtr<=2&pass ==1], na.rm =T),
    off_1h_rush_epa = mean(epa[qtr<=2&rush ==1], na.rm =T),
    off_2h_epa = mean(epa[qtr>2],na.rm =T),
    off_2h_pass_epa = mean(epa[qtr>2&pass ==1], na.rm =T),
    off_2h_rush_epa = mean(epa[qtr>2&rush ==1], na.rm =T),
    off_scramble_rate = mean(qb_scramble[pass == 1]),
    off_scramble_epa = mean(epa[pass == 1 & qb_scramble == 1]),
    # off_no_blitz_epa = mean(epa[blitz == 0 & pass == 1],na.rm = T),
    # off_blitz_epa = mean(epa[blitz == 1 & pass == 1],na.rm = T),
    # off_blitz_scramble_epa = mean(epa[blitz == 1 & pass == 1 & qb_scramble == 1],na.rm = T),
    # off_blitz_scramble_ypc = mean(yards_gained[blitz == 1 & pass == 1&qb_scramble == 1],na.rm = T),
    # off_adot_blitz = mean(air_yards[blitz == 1 & pass ==1],na.rm =T),
    # blitz_sens = off_blitz_epa - off_no_blitz_epa,
    # motion_rate = mean(is_motion),
    # pass_motion_rate = mean(is_motion[pass == 1]),
    short_pass_rate = mean(short_throw[pass == 1],na.rm =T),
    medium_pass_rate = mean(medium_throw[pass == 1],na.rm =T),
    long_pass_rate = mean(long_throw[pass == 1],na.rm =T),
    off_short_pass_epa = mean(epa[short_throw==1],na.rm =T),
    off_medium_pass_epa = mean(epa[medium_throw ==1],na.rm =T),
    off_long_pass_epa = mean(epa[long_throw == 1],na.rm =T),
    # rush_motion_rate = mean(is_motion[rush == 1]),
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
    # off_no_motion_epa = mean(epa[is_motion == FALSE]),
    # off_pass_no_motion_epa = mean(epa[is_motion == FALSE & pass ==1 ]),
    # off_rush_no_motion_epa = mean(epa[is_motion == FALSE & rush == 1 ]),
    # off_motion_epa = mean(epa[is_motion == TRUE]),
    # off_pass_motion_epa = mean(epa[is_motion == TRUE & pass ==1 ]),
    # off_rush_motion_epa = mean(epa[is_motion == TRUE & rush == 1 ]),
    # off_motion_impr = off_motion_epa-off_no_motion_epa,
    # off_pass_motion_impr = off_pass_motion_epa - off_pass_no_motion_epa,
    # off_rush_motion_impr = off_rush_motion_epa - off_rush_no_motion_epa,
    # play_action_rate = mean(is_play_action[pass == 1]),
    # off_play_action_epa = mean(epa[pass == 1 & is_play_action == TRUE]),
    # off_no_play_action_epa = mean(epa[pass == 1 & is_play_action == FALSE]),
    # off_shotgun_rate = mean(shotgun),
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
    # off_leave_pocket_rate = mean(is_qb_out_of_pocket[pass == 1]),
    # off_epa_in_pocket = mean(epa[pass == 1 & is_qb_out_of_pocket == 0]),
    # off_epa_out_pocket  = mean(epa[pass == 1 & is_qb_out_of_pocket == 1]),
    # intereception_worthy_rate = mean(is_interception_worthy[pass==1]),
    # intereception_worthy_epa = mean(epa[pass==1&is_interception_worthy ==1]),
    # rpo_rate = mean(is_rpo),
    # off_rpo_epa = mean(epa[is_rpo == 1]),
    # screen_rate = mean(is_screen_pass[pass == 1]),
    # off_screen_epa = mean(epa[is_screen_pass == 1]),
    # off_rush_in_heavy_box_rate = mean(rush[heavybox == 1]),
    # off_rush_in_heavy_box_epa = mean(epa[heavybox == 1 & rush == 1]),
    # off_pass_in_heavy_box_epa = mean(epa[heavybox == 1 & pass == 1]),
    # off_rush_in_light_box_rate = mean(rush[lightbox == 1]),
    # off_rush_in_light_box_epa = mean(epa[lightbox == 1 & rush == 1]),
    # off_pass_in_light_box_epa = mean(epa[lightbox == 1 & pass == 1]),
    # off_no_huddle_rate = mean(is_no_huddle),
    # off_no_huddle_epa = mean(epa[is_no_huddle == 1]),
    home_off_epa = mean(epa[posteam == home_team]),
    home_pass_off_epa = mean(epa[posteam == home_team&pass == 1]),
    home_rush_off_epa = mean(epa[posteam == home_team&pass == 0]),
    away_off_epa = mean(epa[posteam == away_team]),
    away_pass_off_epa = mean(epa[posteam == away_team&pass == 1]),
    away_rush_off_epa = mean(epa[posteam == away_team&pass == 0])
  ) %>% 
  mutate(team = paste(posteam,season)) %>% 
  select(-season)

defensive_scouting <- pbp %>%
  group_by(defteam,season) %>% 
  summarize(
    def_epa = mean(epa),
    def_rush_epa = mean(epa[pass == 0]),
    def_pass_epa = mean(epa[pass == 1]),
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
    def_2h_epa = mean(epa[qtr<=2],na.rm =T),
    def_2h_pass_epa = mean(epa[qtr<=2&pass ==1], na.rm =T),
    def_2h_rush_epa = mean(epa[qtr<=2&rush ==1], na.rm =T),
    def_scrabmle_epa = mean(epa[pass == 1 & qb_scramble == 1]),
    # def_blitz_rate = mean(blitz),
    # def_no_blitz_epa = mean(epa[blitz == 0 & pass == 1]),
    # def_blitz_epa = mean(epa[blitz == 1 & pass == 1],na.rm = T),
    # def_blitz_scramble_epa = mean(epa[blitz == 1 & pass == 1&qb_scramble == 1],na.rm = T),
    # def_blitz_scramble_ypc = mean(yards_gained[blitz == 1 & pass == 1&qb_scramble == 1],na.rm = T),
    # def_no_motion_epa = mean(epa[is_motion == FALSE]),
    # def_pass_no_motion_epa = mean(epa[is_motion == FALSE & pass ==1 ]),
    # def_rush_no_motion_epa = mean(epa[is_motion == FALSE & rush == 1 ]),
    # def_motion_epa = mean(epa[is_motion == TRUE]),
    # def_pass_motion_epa = mean(epa[is_motion == TRUE & pass ==1 ]),
    # def_rush_motion_epa = mean(epa[is_motion == TRUE & rush == 1 ]),
    # def_motion_sens = def_motion_epa-def_no_motion_epa,
    # def_pass_motion_sens = def_pass_motion_epa - def_pass_no_motion_epa,
    # def_rush_motion_sens = def_rush_motion_epa - def_rush_no_motion_epa,
    # def_play_action_epa = mean(epa[pass == 1 & is_play_action == TRUE]),
    # def_no_play_action_epa = mean(epa[pass == 1 & is_play_action == FALSE]),
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
    # def_leave_pocket_rate = mean(is_qb_out_of_pocket[pass == 1]),
    # def_epa_in_pocket = mean(epa[pass == 1 & is_qb_out_of_pocket == 0]),
    # def_epa_out_pocket  = mean(epa[pass == 1 & is_qb_out_of_pocket == 1]),
    # def_rpo_epa = mean(epa[is_rpo == 1]),
    # def_screen_epa = mean(epa[is_screen_pass == 1]),
    # def_heavy_box_rate = mean(heavybox),
    # def_rush_in_heavy_box_epa = mean(epa[heavybox == 1 & rush == 1]),
    # def_pass_in_heavy_box_epa = mean(epa[heavybox == 1 & pass == 1]),
    # def_light_box_rate = mean(lightbox),
    # def_rush_in_light_box_epa = mean(epa[lightbox == 1 & rush == 1]),
    # def_pass_in_light_box_epa = mean(epa[lightbox == 1 & pass == 1]),
    # def_no_huddle_epa = mean(epa[is_no_huddle == 1]),
    home_def_epa = mean(epa[defteam == home_team]),
    home_pass_def_epa = mean(epa[defteam == home_team&pass == 1]),
    home_rush_def_epa = mean(epa[defteam == home_team&pass == 0]),
    away_def_epa = mean(epa[defteam == away_team]),
    away_pass_def_epa = mean(epa[defteam == away_team&pass == 1]),
    away_rush_def_epa = mean(epa[defteam == away_team&pass == 0])
  ) %>% 
  mutate(team = paste(defteam,season))

scouting<- offensive_scouting %>% 
  left_join(defensive_scouting, by = c("team")) %>% 
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr"))

euclidean_distance <- function(x, y) {
  sqrt(sum((x - y)^2))
}



manhattan_distance <- function(point1, point2) {
  sum(abs(point1 - point2))
}

cosine_distance <- function(vector1, vector2) {
  dot_product <- sum(vector1 * vector2)
  magnitude1 <- sqrt(sum(vector1^2))
  magnitude2 <- sqrt(sum(vector2^2))
  
  return(1 - dot_product / (magnitude1 * magnitude2))
}


# Function to calculate Euclidean distance between one row and all rows in the dataframe
calculate_similarity <- function(row_index, df) { #how can I look at feature selection
  row <- as.numeric(df[row_index, ])
  eu_distances <- apply(df, 1, function(x) euclidean_distance(row, as.numeric(x)))
  man_distances<- apply(df, 1, function(x) manhattan_distance(row, as.numeric(x)))
  cos_distances<- apply(df, 1, function(x) cosine_distance(row, as.numeric(x)))
  
  distances <- scale(cbind(eu_distances,man_distances,cos_distances)) %>% 
    as.data.frame(.) %>% 
    mutate(avg_dist =rowMeans(.))
  return(distances)
}

library(dbscan)
library(zoo)
test<-kNN(scale(Filter(is.numeric,scouting) %>% select(-season) %>% 
                         mutate_all(~na.aggregate(.))),k = 575)
scouting[test$id[which(scouting$team == "SF 2023"),c(1:10)],"team"]
1- test$dist[which(scouting$team == "SF 2023"),c(1:10)]/max(test$dist[which(scouting$team == "SF 2023"),])

top_10_indices <- test$id[which(scouting$team == "SF 2023"),c(1:10)]
similarity <- c(1,1- test$dist[which(scouting$team == "SF 2023"),c(1:10)]/max(test$dist[which(scouting$team == "SF 2023"),]))

scouting[c(which(scouting$team == "SF 2023"),top_10_indices),] %>% 
  ungroup() %>% 
  mutate(rank = ifelse(row_number()==1,"",row_number()-1),
         sim_score = ifelse(similarity < 1,paste(round(similarity,4)*100,"%"),"")) %>%
  select(rank,sim_score,team_wordmark,season,off_epa, off_pass_epa, off_rush_epa, pass_rate, def_epa, def_pass_epa, def_rush_epa)
# Row index for which similarity is to be calculated
hold<-test$dist
# Calculate similarity between the specified row and all rows in the dataframe
similarities <- cbind(calculate_similarity(which(scouting$team== "SF 2023"), scale(Filter(is.numeric,scouting %>% select(-season)))),scouting$team)
sorted_indices <- order(similarities$avg_dist)
top_11_indices <- sorted_indices[1:11]
scouting[top_11_indices,"team"]
