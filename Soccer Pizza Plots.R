library(worldfootballR)
library(ggshakeR)
library(dplyr)




data1 <- fb_player_scouting_report(pull(player_dictionary_mapping() %>% filter(PlayerFBref  == "Rafael Leão") %>% select(UrlFBref)), pos_versus = "primary", time_pause =  10)
data2 <- fb_player_scouting_report(pull(player_dictionary_mapping() %>% filter(PlayerFBref  == "Nico Williams") %>% select(UrlFBref)), pos_versus = "primary",time_pause =  10)

strikers <- c(#"Goals","Assists","Non-Penalty Goals",
              "xG: Expected Goals", "npxG: Non-Penalty xG","Progressive Carries","Progressive Passes",
              "npxG/Shot","xA: Expected Assists","Shot-Creating Actions","Goal-Creating Actions",
              "Tackles (Att 3rd)","Interceptions","Touches (Att Pen)","Successful Take-Ons",
              "Successful Take-On %","Ball Recoveries","Aerials Won","% of Aerials Won")

statgroup_rankings <- rbind(data1, data2) %>%
  count(StatGroup, name = "n") %>%
  arrange(n) %>%
  mutate(stat_priority = row_number())

# 3. Join back to the main data
data <- rbind(data1, data2) %>%
  filter(Statistic %in% strikers) %>% 
  filter(scouting_period == "Last 365 Days Men's Big 5 Leagues, UCL, UEL") %>% 
  left_join(statgroup_rankings, by = "StatGroup") %>%
  group_by(Player, Statistic) %>%
  arrange(stat_priority, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>% 
  mutate(Statistic = factor(Statistic, levels = strikers)) %>% 
  arrange(Player, Statistic)
  


plot_pizza(data = data, type = "comparison", template = "custom",
                         player_1 = first(data$Player), player_2 = last(data$Player),
                         season_player_1 = "Last 365 Days Men's Big 5 Leagues, UCL, UEL", 
                         season_player_2 = "Last 365 Days Men's Big 5 Leagues, UCL, UEL",
                         color_compare = "red", theme = "black")
