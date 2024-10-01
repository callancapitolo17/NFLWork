library(nflreadr)
library(dplyr)
library(gt)
library(ggplot2)
library(gtExtras)

teams <- load_teams() 
draft_data <- load_draft_picks() %>%
  mutate(clean_team = ifelse(team %in% c("OAK","LVR","RAI"),"LV",ifelse(team %in% c("STL", "LAR","RAM"),"LA",
                ifelse(team == "GNB","GB", ifelse(team == "NOR", "NO",ifelse(team == "SFO", "SF",
                ifelse(team == "KAN", "KC",ifelse(team == "NWE","NE",ifelse(team == "TAM","TB",
                ifelse(team %in% c("LAC","SDG"), "LAC", ifelse(team == "PHO","ARI",team)))))))))))

draft_data %>% 
  group_by(category,season) %>% 
  summarize(earliest_selection = min(pick), first_round = sum(round == 1)) %>% 
  filter(!(category %in% c("P","K","KR","LS")), category %in% c("RB","WR","TE"), season >= 2000) %>% 
  ggplot(aes(x = season, y = first_round, color = category))+
  geom_line(lwd = 2)

age_draft <- draft_data %>% 
  group_by(season,round) %>% 
  filter(round<=7) %>% 
  # summarize(draft_avg_age = mean(age,na.rm = T), round_1_age = mean(age[round == 1])) %>% 
  summarize(draft_avg_age = mean(age,na.rm = T))

draft_data %>% 
  group_by(season) %>% 
  summarize(avg_age = mean(age,na.rm = T)) %>% 
  ggplot(aes(x = season, y = avg_age))+
  geom_line(lwd = 1.5) +
  labs(title = "How has Average Age of Draftees Evolved?", 
       subtitle = "Draftees Have Gotten Younger Since the Mid 2000s, but we Have Seen a Sharp Increase in Age the Past 2 Drafts",
       x = "Draft Year", y = "Average Age", caption = "Callan Capitolo | @CapAnalytics7 | nflverse")+
  theme_classic()
ggsave("DraftAgeChange.png", width = 14, height =10, dpi = "retina")

age_draft%>% 
  ggplot(aes(x = season, y = draft_avg_age, color = as.factor(round)))+
  geom_line(lwd = 1.5,aes(alpha = ifelse(round == 1, 1, 0.99)))+
  labs(title = "How has the Average Age of Draftees Changed Based on Round Selected?",
       subtitle = "First Round Selections Have Always Been the Youngest, but Since 2014 they are Getting Even Younger", x = "Draft Year", y = "Average Age",
       color = "Round", caption = "Callan Capitolo | @CapAnalytics7 | nflverse")+
  geom_vline(xintercept = 2014, linetype = "dashed", color = "black")+
  theme_classic()+guides(alpha = FALSE)
ggsave("FirstRoundAge.png", width = 14, height =10, dpi = "retina")

test<-draft_data %>% 
  filter(round <= 7) %>% 
  group_by(round) %>% 
  summarize(onwards14 = mean(age[season >= 2014],na.rm = T),
            pre14 = mean(age[season < 2014],na.rm = T))
#which organizations drove these changes post 2014 team by first round draft age with number of draft picks

team_first_round <- draft_data %>% 
  filter(round == 1) %>% 
  group_by(clean_team) %>% 
  summarize(post14selections = sum(season>=2014), onwards14 = round(mean(age[season >= 2014],na.rm = T),2),
            numpre14 = sum(season<2014),pre14 = round(mean(age[season < 2014],na.rm = T),2)) %>% 
  mutate(diff = onwards14-pre14) %>% 
  left_join(teams, by = c("clean_team" = "team_abbr")) %>% 
  arrange(onwards14)

teamfirstbehavior<- team_first_round %>% 
  select(team_wordmark,onwards14,pre14,diff) %>% 
  gt() %>% 
  cols_align(align = "center") %>%
  gtExtras::gt_img_rows(team_wordmark) %>% 
  cols_label(team_wordmark = "Team",
             onwards14 = "1st Rounders Avg Age Since 2014",
             pre14 = "1st Rounders Avg Age Before 2014",
             diff = "Avg Age Difference") %>% 
  gt_theme_538() %>% 
  gt_hulk_col_numeric(columns = is.numeric) %>% 
  tab_header(
    title = md("How has Team Drafting Behavior Surrounding 1st Rounders' Age Changed Since 2014?"),
    subtitle = md("The Chiefs and Bills Have Led the Change in Draft Strategy Surrounding 1st Rounders' Age")
  )
gtsave(teamfirstbehavior, "teamfirstbehavior.png")
#break age down by position too and maybe age by position and round


# QB by Age ----

allqbs<- draft_data %>% 
  filter(position == "QB",season >=2000) %>% 
  group_by(age) %>% 
  summarize(count = n(), avg_pick_num = round(mean(pick),0), avg_w_av = mean(w_av,na.rm = T),
            avg_szn = mean(seasons_started, na.rm = T),
            avg_all_pro= mean(allpro,na.rm = T),
            w_av_season = mean(w_av/ifelse(seasons_started==0,1,seasons_started),na.rm = T),
            pct_first_rounds = mean(round == 1,na.rm = T)) %>% 
  filter(count >20) %>% 
  mutate_if(is.numeric,~round(.,2)) %>% 
  gt() %>% 
  cols_align(align = "center") %>% 
  cols_label(age = "Age", count = "# of Selections", avg_pick_num = "Avg Pick Number",
             avg_w_av = "Avg Career Weighted Approximate Value", avg_szn = "Avg Seasons as Starter",
             avg_all_pro = "Avg Career All Pros", w_av_season = "Avg Weighted Approximate Value per Season as Starter",
             pct_first_rounds = "% Selections in 1st Round") %>% 
  gt_theme_538() %>% 
  gt_hulk_col_numeric(columns = -age) %>% 
  tab_header(
    title = md("How do QBs Compare Based on their Draft Age Since 2000?"),
    subtitle = md("It would appear that the younger the QB the more sucessful career, but that does not tell the full story as younger QBs are normally better prospects")
  )
gtsave(allqbs, "allqbs.png")
  

firstqbs <- draft_data %>% 
  filter(position == "QB",season >=2000,round==1) %>% 
  group_by(age) %>% 
  summarize(count = n(), avg_pick_num = round(mean(pick),0), avg_w_av = mean(w_av,na.rm = T),
            avg_szn = mean(seasons_started, na.rm = T),
            avg_all_pro= mean(allpro,na.rm = T),
            w_av_season = mean(w_av/ifelse(seasons_started==0,1,seasons_started),na.rm = T)) %>% 
  filter(age <25) %>% 
  mutate_if(is.numeric,~round(.,2)) %>% 
  gt() %>% 
  cols_align(align = "center") %>% 
  cols_label(age = "Age", count = "# of Selections", avg_pick_num = "Avg Pick Number",
             avg_w_av = "Avg Career Weighted Approximate Value", avg_szn = "Avg Seasons as Starter",
             avg_all_pro = "Avg Career All Pros", w_av_season = "Avg Weighted Approximate Value per Season as Starter") %>% 
  gt_theme_538() %>% 
  gt_hulk_col_numeric(columns = -age) %>% 
  tab_header(
    title = md("How do 1st Round QBs Compare Based on their Draft Age Since 2000?"),
    subtitle = md("The picture becomes less clear when looking at just first rounders and actually suggests that age might not be that important when drafting a 1st round QB")
  )
gtsave(firstqbs, "firstqbs.png")

draft_data %>%
  ggplot(aes(x = probowls, y = w_av))+
  geom_point()
#younger players have better careers but also related to quality of player coming out as shown by average pick number
#option to look at by round to account for quality of player
#avg draft positon by season by age
draft_data %>% 
  filter(position == "QB",season >=2000,round == 1, age == 24)


# GM Draft Traits ----
combine <- load_combine()

lynch_avg <- combine %>% 
  filter(season >= 2017 & draft_team == "San Francisco 49ers") %>% 
  mutate(inches = as.numeric(substr(ht,1,1))*12+ as.numeric(sub(".*-(\\d+)", "\\1", ht))) %>% 
  left_join(draft_data, by = c("pfr_id" = "pfr_player_id")) %>% 
  group_by(position) %>% 
  summarize(count = n(), avg_age = mean(age,na.rm = T), first_round = sum(round ==1),avg_wt = mean(wt,na.rm = T),
            avg_forty = mean(forty,na.rm = T), avg_bench = mean(bench,na.rm = T),
            avg_vert = mean(vertical,na.rm = T), avg_broad = mean(broad_jump,na.rm = T),
            avg_cone = mean(cone,na.rm = T), avg_shuttle = mean(shuttle,na.rm = T),
            avg_height = mean(inches,na.rm = T)) %>% 
  # filter(count >1) %>% 
  mutate_if(is.numeric, ~round(.,2)) %>% 
  gt() %>% 
  cols_align(align = "center") %>% 
  cols_label(count = "# of Draft Picks", avg_age = "Avg Age", first_round = "# in 1st round", avg_wt = "Avg Weight",
             avg_forty = "Avg 40",avg_bench = "Avg Bench", avg_vert = "Avg Vertical", avg_broad = "Avg Broad",
             avg_cone = "Avg Cone", avg_shuttle = "Avg Shuttle", avg_height = "Avg Height") %>% 
  gtExtras::gt_theme_538() %>% 
  gt_hulk_col_numeric(columns = is.numeric) %>% 
  tab_header(title = md("What are John Lynch's draft tendencies based on combine metrics since coming to the 49ers?"))
gtsave(lynch_avg,"lynchdraft.png")

combine17 <- combine %>% 
  mutate(inches = as.numeric(substr(ht,1,1))*12+ as.numeric(sub(".*-(\\d+)", "\\1", ht))) 
league_avg <- draft_data %>%  
  filter(season>=2017) %>% 
  left_join(combine17, by = c("pfr_player_id" = "pfr_id")) %>% 
  group_by(position) %>% 
  summarize(count = n(), avg_age = mean(age,na.rm = T), first_round = sum(round ==1),avg_wt = mean(wt,na.rm = T),
            avg_forty = mean(forty,na.rm = T), avg_bench = mean(bench,na.rm = T),
            avg_vert = mean(vertical,na.rm = T), avg_broad = mean(broad_jump,na.rm = T),
            avg_cone = mean(cone,na.rm = T), avg_shuttle = mean(shuttle,na.rm = T),
            avg_height = mean(inches,na.rm = T)) %>% 
  filter(count >20) %>% 
  mutate_if(is.numeric, ~round(.,2)) %>% 
  gt() %>% 
  cols_align(align = "center") %>% 
  cols_label(count = "# of Draft Picks", avg_age = "Avg Age", first_round = "# in 1st round", avg_wt = "Avg Weight",
             avg_forty = "Avg 40",avg_bench = "Avg Bench", avg_vert = "Avg Vertical", avg_broad = "Avg Broad",
             avg_cone = "Avg Cone", avg_shuttle = "Avg Shuttle", avg_height = "Avg Height") %>% 
  gtExtras::gt_theme_538() %>% 
  gt_hulk_col_numeric(columns = is.numeric) %>% 
  tab_header(title = md("Combine Measurables for Draftees Since 2017"))
gtsave(league_avg,"leaguedraft.png")

comparison <- cbind(lynch_avg$position,lynch_avg %>% select(-position) - league_avg %>% filter(position != "C") %>%  select(-position))
lynch_comp<-comparison %>% 
  select(-count,-first_round) %>%
  filter(`lynch_avg$position` != "K", `lynch_avg$position` != "P") %>% 
  gt() %>% 
  cols_align(align = "center") %>% 
  cols_label(`lynch_avg$position` = "Position Category", avg_age = "Avg Age", avg_wt = "Avg Weight",
             avg_forty = "Avg 40",avg_bench = "Avg Bench", avg_vert = "Avg Vertical", avg_broad = "Avg Broad",
             avg_cone = "Avg Cone", avg_shuttle = "Avg Shuttle", avg_height = "Avg Height") %>% 
  gtExtras::gt_theme_538() %>% 
  gt_hulk_col_numeric(columns = is.numeric) %>% 
  tab_header(title = md("How do John Lynch's Draft Tendencies Pre 2024 Draft Compare to League Average"),subtitle = md("Each value represents 49ers Average Under Lynch - League Average Since 2017"))
gtsave(lynch_comp,"lynchcomp.png")

lynch_comp

# Patriots Offensive Line
test<-draft_data %>%  
  filter(season >= 2000, !is.na(pfr_player_id), team == "NWE", category == "OL") %>% 
  left_join(combine17, by = c("pfr_player_id" = "pfr_id"))
