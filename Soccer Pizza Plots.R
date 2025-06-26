library(worldfootballR)
library(ggshakeR)
library(dplyr)
library(stringi)
library(ggplot2)
library(ggtext)
player_comps<- function(player1,player2,stat_type = "strikers",league1 = "Last 365 Days Men's Big 5 Leagues, UCL, UEL",league2 = "Last 365 Days Men's Big 5 Leagues, UCL, UEL"){

data1 <- fb_player_scouting_report(pull(player_dictionary_mapping() %>% filter(PlayerFBref  == player1) %>% select(UrlFBref)%>% slice(1)), pos_versus = "primary", time_pause =  5)
data2 <- fb_player_scouting_report(pull(player_dictionary_mapping() %>% filter(PlayerFBref  == player2) %>% select(UrlFBref) %>% slice(1)), pos_versus = "primary",time_pause =  5)

strikers <- c(#"Goals","Assists","Non-Penalty Goals",
              "xG: Expected Goals", "npxG: Non-Penalty xG","npxG/Shot","xA: Expected Assists","Shot-Creating Actions","Goal-Creating Actions",
              "Progressive Passes", "Progressive Carries","Successful Take-Ons", "Successful Take-On %","Touches (Att Pen)","Aerials Won","% of Aerials Won",
              "Tackles (Att 3rd)","Interceptions","Ball Recoveries")
midfielders <- c("xG: Expected Goals", "npxG: Non-Penalty xG","xA: Expected Assists", "Progressive Passes", 
                 "Progressive Passes Rec","Passes into Final Third", "Shot-Creating Actions",
                 "Progressive Carries","Successful Take-Ons", "Successful Take-On %", "Carries into Final Third",
                 "Tackles Won", "Interceptions",
                 "Ball Recoveries","% of Aerials Won")

# 3. Join back to the main data
data <- rbind(data1, data2) %>%
  filter(scouting_period %in% c(league1,league2)) %>% 
  filter(Statistic %in% get(stat_type)) %>% 
  group_by(Player, Statistic) %>%
  arrange(StatGroup, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>% 
  mutate(Statistic = factor(Statistic, levels = get(stat_type))) %>% 
  arrange(Player, Statistic)
  


comp_chart <- plot_pizza(data = data, type = "comparison", template = "custom",
                         player_1 = stri_trans_general(player1, "Latin-ASCII"), player_2 = stri_trans_general(player2, "Latin-ASCII"),
                         # season_player_1 = "Last 365 Days Men's Big 5 Leagues, UCL, UEL", 
                         # season_player_2 = "Last 365 Days Men's Big 5 Leagues, UCL, UEL",
                         color_compare = "red", theme = "white")+
  labs(
    title = paste0(
      "<span style='color:red'>",player1,"</span> vs ",
      "<span style='color:black'>",player2,"</span>"
    ),
    subtitle = "@CapAnalytics7")+
  theme(
    # switch title to element_markdown so those spans get rendered
    plot.title = element_markdown(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_markdown(hjust = 0.5, size = 12)
  )+
  scale_x_discrete(limits = get(stat_type),labels = function(x) stringr::str_wrap(x, width = 8))
comp_chart
ggsave("image.png", width = 4200, height = 2800, units = "px") 
print(comp_chart)
}