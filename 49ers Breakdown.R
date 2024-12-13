#49ers This Year vs Last----
shanahan_in<- nfl99 %>% filter(pass == 1, season >2022, posteam == "SF", name == "B.Purdy") %>% 
  left_join(ftn_data, by = c("game_id" = "nflverse_game_id",
                             "play_id" = "nflverse_play_id"))
shanahan_in %>% 
  group_by(season) %>% 
  summarize(`EPA/PA Play` = mean(epa[is_play_action],na.rm = T),
            `EPA/RPO` = mean(epa[is_rpo == 1],na.rm = T),
            `EPA/Screen` = mean(epa[is_screen_pass == 1],na.rm = T),
            `Non Help DB` = mean(epa[is_screen_pass == 0 & is_rpo == 0 & is_play_action == 0],na.rm = T),
            `YAC Over Expected` = mean(yards_after_catch-xyac_mean_yardage,na.rm = T)
            ) %>% 
  mutate_if(is.numeric,~round(.,2)) %>% 
  gt() %>% 
  cols_align(align = "center") %>% 
  gtExtras::gt_theme_538() %>% 
  tab_header(title = md("Where has it Gone Wrong for the 49ers Offense?"),
            subtitle = md("The 49ers are just not as good at what made them great last year with the big issue coming from a lack of YAC")) %>% 
  tab_footnote(footnote = md("Only Purdy Dropbacks @CapAnalytics7|nflfastr | FTN"))
shanahan_in %>% 
    group_by(season) %>% 
    summarize(`EPA/PA Play` = mean(comp_yac_epa[is_play_action],na.rm = T),
              `EPA/RPO` = mean(comp_yac_epa[is_rpo == 1],na.rm = T),
              `EPA/Screen` = mean(comp_yac_epa[is_screen_pass == 1],na.rm = T),
              `Non Help DB` = mean(comp_yac_epa[is_screen_pass == 0 & is_rpo == 0 & is_play_action == 0],na.rm = T),
              `YAC EPA/All Passes` = mean(comp_yac_epa, na.rm = T)
              ) %>% 
  mutate_if(is.numeric,~round(.,2)) %>% 
  gt() %>% 
  cols_align(align = "center") %>% 
  gtExtras::gt_theme_538() %>% 
  tab_header(title = md("49ers Have seen a Major Regression in YAC Across All Pass Categories")) %>% 
  tab_footnote(footnote = md("Only Purdy Dropbacks and Completed Passes @CapAnalytics7|nflfastr | FTN"))

shanahan_in %>% 
  group_by(season) %>% 
  summarize(`Air Yards EPA/PA Play` = mean(comp_air_epa[is_play_action],na.rm = T),
            `Air Yards EPA/RPO` = mean(comp_air_epa[is_rpo == 1],na.rm = T),
            `Air Yards EPA/Screen` = mean(comp_air_epa[is_screen_pass == 1],na.rm = T),
            `Air Yards Non Help DB` = mean(comp_air_epa[is_screen_pass == 0 & is_rpo == 0 & is_play_action == 0],na.rm = T),
            `Air Yards EPA/All Passes` = mean(comp_air_epa, na.rm = T), 
            `ADoT` = mean(air_yards,na.rm = T)) %>% 
    mutate_if(is.numeric,~round(.,2)) %>% 
    gt() %>% 
    cols_align(align = "center") %>% 
    gtExtras::gt_theme_538() %>% 
    tab_header(title = md("When Removing YAC EPA the Niners Offense Has Actually Improved Compared to Last Year")) %>% 
    tab_footnote(footnote = md("Only Purdy Dropbacks and Completed Passes @CapAnalytics7|nflfastr | FTN"))
