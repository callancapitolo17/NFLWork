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
  mutate(explosive = ifelse((yardline_100 <  20 & pass_attempt == 1) | (yardline_100<12 & (qb_scramble ==1 |rush ==1)), NA, ifelse((yards_gained>20 & pass_attempt == 1) | (yards_gained >12 & (qb_scramble == 1 | rush == 1)),1,0)),
         negative = ifelse(yards_gained < 0, 1,0))

replace_with_ranks<- function(column){
  values <- column
  ranks <- rank(column*-1,ties.method = "max")
}

summarize_pbp_data <- function(data,group_col) {
  data %>% 
    group_by(across(all_of(group_col))) %>% 
    summarize(
      epa_avg = mean(epa, na.rm = TRUE),
      success_rate = mean(success, na.rm = TRUE),
      pass_epa = mean(epa[pass == 1], na.rm = TRUE),
      dropback_success_rate = mean(success[pass == 1], na.rm = TRUE),
      rush_epa = mean(epa[pass == 0], na.rm = TRUE),
      rush_success_rate = mean(success[rush == 1], na.rm = TRUE),
      explosive_rate = mean(explosive, na.rm = TRUE),
      explosive_rush_rate = mean(explosive[rush == 1], na.rm = TRUE),
      explosive_pass_rate = mean(explosive[pass == 1], na.rm = TRUE),
      negative_rate = mean(yards_gained < 0, na.rm = TRUE),
      negative_pass_rate = mean(yards_gained[pass == 1] < 0, na.rm = TRUE),
      negative_rush_rate = mean(yards_gained[rush == 1] < 0, na.rm = TRUE),
      turnover_rate = mean(interception + fumble_lost, na.rm = TRUE),
      fumble_rate = mean(fumble, na.rm = TRUE),
      fumble_lost_rate = mean(fumble_lost, na.rm = TRUE),
      early_down_epa = mean(epa[down %in% c(1, 2)], na.rm = TRUE),
      early_down_pass_epa = mean(epa[down %in% c(1, 2) & pass == 1], na.rm = TRUE),
      early_down_rush_epa = mean(epa[down %in% c(1, 2) & rush == 1], na.rm = TRUE),
      early_down_pass_rate = mean(pass[down == 1 | down == 2], na.rm = TRUE),
      early_down_pass_oe = mean(pass_oe[down == 1 | down == 2], na.rm = TRUE),
      early_down_first_pct = mean(first_down[down == 1 | down == 2], na.rm = TRUE),
      third_down_first_pct = mean(first_down[down == 3], na.rm = TRUE),
      third_down_distance_avg = mean(ydstogo[down == 3], na.rm = TRUE),
      third_down_epa = mean(epa[down == 3], na.rm = TRUE),
      red_zone_td_pct = mean(touchdown[yardline_100 <= 20], na.rm = TRUE),
      red_zone_td_epa = mean(epa[yardline_100 <= 20], na.rm = TRUE),
      yac_pct = sum(yards_after_catch, na.rm = TRUE) / sum(yards_gained[pass == 1], na.rm = TRUE),
      sack_rate = mean(sack[pass == 1], na.rm = TRUE),
      scramble_rate = mean(qb_scramble[pass == 1], na.rm = TRUE),
      scramble_epa = mean(epa[pass == 1 & qb_scramble == 1], na.rm = TRUE),
      short_pass_rate = mean(short_throw[pass == 1], na.rm = TRUE),
      medium_pass_rate = mean(medium_throw[pass == 1], na.rm = TRUE),
      long_pass_rate = mean(long_throw[pass == 1], na.rm = TRUE),
      short_pass_epa = mean(epa[short_throw == 1], na.rm = TRUE),
      medium_pass_epa = mean(epa[medium_throw == 1], na.rm = TRUE),
      long_pass_epa = mean(epa[long_throw == 1], na.rm = TRUE),
      pass_rate = mean(pass, na.rm = TRUE),
      pass_oe_avg = mean(pass_oe, na.rm = TRUE),
      epa_xpass_sit = mean(epa[pass_oe >= 0.9], na.rm = TRUE),
      pass_right_rate = mean(pass_location == "right", na.rm = TRUE),
      pass_right_epa = mean(epa[pass_location == "right"], na.rm = TRUE),
      pass_left_rate = mean(pass_location == "left", na.rm = TRUE),
      pass_left_epa = mean(epa[pass_location == "left"], na.rm = TRUE),
      pass_middle_rate = mean(pass_location == "middle", na.rm = TRUE),
      pass_middle_epa = mean(epa[pass_location == "middle"], na.rm = TRUE),
      aDot = mean(air_yards, na.rm = TRUE),
      play_action_rate = mean(is_play_action[pass == 1], na.rm = TRUE),
      play_action_epa = mean(epa[pass == 1 & is_play_action == TRUE], na.rm = TRUE),
      no_play_action_epa = mean(epa[pass == 1 & is_play_action == FALSE], na.rm = TRUE),
      leave_pocket_rate = mean(is_qb_out_of_pocket[pass == 1], na.rm = TRUE),
      epa_in_pocket = mean(epa[pass == 1 & is_qb_out_of_pocket == 0], na.rm = TRUE),
      epa_out_pocket = mean(epa[pass == 1 & is_qb_out_of_pocket == 1], na.rm = TRUE),
      interception_rate = mean(interception[pass == 1], na.rm = TRUE),
      interception_worthy_rate = mean(is_interception_worthy[pass == 1], na.rm = TRUE),
      interception_worthy_epa = mean(epa[pass == 1 & is_interception_worthy == 1], na.rm = TRUE),
      screen_rate = mean(is_screen_pass[pass == 1], na.rm = TRUE),
      screen_epa = mean(epa[is_screen_pass == 1], na.rm = TRUE),
      rush_right_rate = mean(run_location == "right", na.rm = TRUE),
      rush_right_epa = mean(epa[run_location == "right"], na.rm = TRUE),
      rush_left_rate = mean(run_location == "left", na.rm = TRUE),
      rush_left_epa = mean(epa[run_location == "left"], na.rm = TRUE),
      rush_middle_rate = mean(run_location == "middle", na.rm = TRUE),
      rush_middle_epa = mean(epa[run_location == "middle"], na.rm = TRUE),
      shotgun_rate = mean(shotgun, na.rm = TRUE),
      pass_pct_from_shotgun = mean(shotgun[pass == 1], na.rm = TRUE),
      rush_pct_from_shotgun = mean(shotgun[rush == 1], na.rm = TRUE),
      shotgun_pass_rate = mean(pass[shotgun == 1], na.rm = TRUE),
      under_center_pass_rate = mean(pass[shotgun == 0], na.rm = TRUE),
      shotgun_epa = mean(epa[shotgun == 1], na.rm = TRUE),
      shotgun_pass_epa = mean(epa[shotgun == 1 & pass == 1], na.rm = TRUE),
      shotgun_rush_epa = mean(epa[shotgun == 1 & rush == 1], na.rm = TRUE),
      no_shotgun_epa = mean(epa[shotgun == 0], na.rm = TRUE),
      no_shotgun_pass_epa = mean(epa[shotgun == 0 & pass == 1], na.rm = TRUE),
      no_shotgun_rush_epa = mean(epa[shotgun == 0 & rush == 1], na.rm = TRUE),
      first_half_epa = mean(epa[qtr <= 2], na.rm = TRUE),
      first_half_pass_epa = mean(epa[qtr <= 2 & pass == 1], na.rm = TRUE),
      first_half_rush_epa = mean(epa[qtr <= 2 & rush == 1], na.rm = TRUE),
      second_half_epa = mean(epa[qtr > 2], na.rm = TRUE),
      second_half_pass_epa = mean(epa[qtr > 2 & pass == 1], na.rm = TRUE),
      second_half_rush_epa = mean(epa[qtr > 2 & rush == 1], na.rm = TRUE),
      blitz_rate = mean(blitz, na.rm = TRUE),
      no_blitz_epa = mean(epa[blitz == 0 & pass == 1], na.rm = TRUE),
      blitz_epa = mean(epa[blitz == 1 & pass == 1], na.rm = TRUE),
      blitz_scramble_epa = mean(epa[blitz == 1 & pass == 1 & qb_scramble == 1], na.rm = TRUE),
      blitz_scramble_ypc = mean(yards_gained[blitz == 1 & pass == 1 & qb_scramble == 1], na.rm = TRUE),
      aDot_blitz = mean(air_yards[blitz == 1 & pass == 1], na.rm = TRUE),
      blitz_sensitivity = blitz_epa - no_blitz_epa,
      motion_rate = mean(is_motion, na.rm = TRUE),
      pass_motion_rate = mean(is_motion[pass == 1], na.rm = TRUE),
      rush_motion_rate = mean(is_motion[rush == 1], na.rm = TRUE),
      no_motion_epa = mean(epa[is_motion == FALSE], na.rm = TRUE),
      motion_epa = mean(epa[is_motion == TRUE], na.rm = TRUE),
      pass_no_motion_epa = mean(epa[is_motion == FALSE & pass == 1], na.rm = TRUE),
      rush_no_motion_epa = mean(epa[is_motion == FALSE & rush == 1], na.rm = TRUE),
      pass_motion_epa = mean(epa[is_motion == TRUE & pass == 1], na.rm = TRUE),
      rush_motion_epa = mean(epa[is_motion == TRUE & rush == 1], na.rm = TRUE),
      motion_improvement = motion_epa - no_motion_epa,
      pass_motion_improvement = pass_motion_epa - pass_no_motion_epa,
      rush_motion_improvement = rush_motion_epa - rush_no_motion_epa,
      rpo_rate = mean(is_rpo, na.rm = TRUE),
      rpo_epa = mean(epa[is_rpo == 1], na.rm = TRUE),
      rush_heavy_box_rate = mean(rush[heavybox == 1], na.rm = TRUE),
      rush_heavy_box_epa = mean(epa[heavybox == 1 & rush == 1], na.rm = TRUE),
      pass_heavy_box_epa = mean(epa[heavybox == 1 & pass == 1], na.rm = TRUE),
      rush_light_box_rate = mean(rush[lightbox == 1], na.rm = TRUE),
      rush_light_box_epa = mean(epa[lightbox == 1 & rush == 1], na.rm = TRUE),
      pass_light_box_epa = mean(epa[lightbox == 1 & pass == 1], na.rm = TRUE),
      no_huddle_rate = mean(is_no_huddle, na.rm = TRUE),
      no_huddle_epa = mean(epa[is_no_huddle == 1], na.rm = TRUE),
      home_epa = mean(epa[posteam == home_team], na.rm = TRUE),
      home_pass_epa = mean(epa[posteam == home_team & pass == 1], na.rm = TRUE),
      home_rush_epa = mean(epa[posteam == home_team & pass == 0], na.rm = TRUE),
      away_epa = mean(epa[posteam == away_team], na.rm = TRUE),
      away_pass_epa = mean(epa[posteam == away_team & pass == 1], na.rm = TRUE),
      away_rush_epa = mean(epa[posteam == away_team & pass == 0], na.rm = TRUE)
    )
}

offensive_scouting <- summarize_pbp_data(data = pbp_rp,group_col = "posteam")
defensive_scouting <- summarize_pbp_data(data = pbp_rp,group_col = "defteam")
rec_offensive_scouting <- summarize_pbp_data(data = pbp_rp %>% filter(week >= (max(week)-4)),group_col = "posteam")
rec_defensive_scouting <- summarize_pbp_data(data = pbp_rp %>% filter(week >= (max(week)-4)),group_col = "defteam")



# Helper function for preprocessing and ranking
prepare_rank_data <- function(data, id_col, negate_cols = NULL, rank_function, name) {
  data_numeric <- data %>% 
    select(-all_of(id_col)) %>% 
    mutate(across(everything(), ~ round(., 3)))
  
  if (!is.null(negate_cols)) {
    data_numeric <- data_numeric %>% 
      mutate(across(all_of(negate_cols), ~ . * -1))
  }
  
  ranked_data <- data_numeric %>% 
    apply(2, rank_function) %>% 
    as.data.frame() %>% 
    cbind(id = data[[id_col]], .) %>% 
    pivot_longer(cols = -id, names_to = "Stat", values_to = paste(name,"Rank", sep =" "))
  
  ranked_data
}

# Define columns to negate
negate_columns_off <- c("negative_rate", "negative_pass_rate", "negative_rush_rate", "third_down_distance_avg", 
                        "sack_rate", "interception_worthy_rate", "interception_rate", "turnover_rate", 
                        "fumble_rate", "fumble_lost_rate")

negate_columns_def <- c("negative_rate", "negative_pass_rate", "negative_rush_rate", "third_down_distance_avg", 
                        "turnover_rate", "fumble_rate", "fumble_lost_rate", "sack_rate", 
                        "interception_rate", "interception_worthy_rate", "blitz_rate")

# Process offense and defense data
data_ranks_off <- prepare_rank_data(
  offensive_scouting, 
  id_col = "posteam", 
  negate_cols = negate_columns_off, 
  rank_function = replace_with_ranks,
  name = "Off"
)

data_ranks_def <- prepare_rank_data(
  defensive_scouting %>% mutate_if(is.numeric,~.*-1), 
  id_col = "defteam", 
  negate_cols = negate_columns_def, 
  rank_function = replace_with_ranks,
  name = "Def"
)

rec_data_ranks_off <- prepare_rank_data(
  rec_offensive_scouting, 
  id_col = "posteam", 
  negate_cols = negate_columns_off, 
  rank_function = replace_with_ranks,
  name = "Rec Off"
)

rec_data_ranks_def <- prepare_rank_data(
  rec_defensive_scouting %>% mutate_if(is.numeric,~.*-1), 
  id_col = "defteam", 
  negate_cols = negate_columns_def, 
  rank_function = replace_with_ranks,
  name = "Rec Def"
)

joined_offense <- data_ranks_off %>% left_join(rec_data_ranks_off, by = c("id", "Stat"))
joined_defense <- data_ranks_def %>% left_join(rec_data_ranks_def, by = c("id", "Stat"))

scouting_guide <- function(){
offense <- readline(prompt = "Enter Offense:")
defense <- readline(prompt = "Enter Defense:")

scouting_report <- joined_offense %>% 
  filter(id == offense) %>% rename("offense" = id) %>%  
  left_join(joined_defense %>% filter(id == defense) %>% rename("defense" = id),
            by = c("Stat")) %>% 
  select(-offense,-defense)

scouting_report %>% 
  mutate(`Rank Diff` = `Off Rank`-`Def Rank`, `Rec Rank Diff` = `Rec Off Rank`-`Rec Def Rank`) %>% 
  select(Stat,`Off Rank`,`Def Rank`,`Rank Diff`,`Rec Off Rank`,`Rec Def Rank`,`Rec Rank Diff`) %>% 
  gt() %>% 
  cols_align(align = "center") %>%
  cols_label(`Stat` = "Stat",`Off Rank` = paste(offense,"Offense Rank", sep = " "), `Def Rank` = paste(defense,"Defense Rank", sep = " "),
             `Rank Diff` = "Rank Difference",
             `Rec Off Rank` = paste(offense,"Rec Offense Rank", sep = " "), `Rec Def Rank` = paste(defense,"Rec Defense Rank", sep = " "),
             `Rec Rank Diff` = "Rec Rank Difference") %>% 
  gtExtras::gt_theme_538() %>% 
  gt_hulk_col_numeric(columns = c(`Rank Diff`, `Rec Rank Diff`)) %>% 
  tab_header(title = md("Scouting Report"), subtitle = md("Purple Represents Offense Advantage, Green Represents Defense Advantage"))
}
