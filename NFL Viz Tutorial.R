library(tidyverse)
library(nflfastR)
library(ggimage)
library(gt)
library(ggthemes)
library(gtExtras)

pbp_2023 <- load_pbp(2023)

qb_epa_play <- pbp_2023 %>%
  filter(pass == 1 | rush ==1) %>%
  filter(!is.na(epa)) %>%
  #filter(week %in% c(1,2,3,4)) %>%
  group_by(id) %>%
  summarize(name = first(name), team = last(posteam), plays = n(), epa_play = mean(epa),
            pass_attempts = sum(incomplete_pass+ complete_pass, na.rm = T)) %>%
  filter(pass_attempts >=50) %>%
  mutate(pass_rate = pass_attempts / plays) %>%
  left_join(teams_colors_logos, by = c("team" = "team_abbr"))
#scatter plot
qb_epa_play %>%
  ggplot(aes(x = pass_rate, y = epa_play)) +
  geom_point(aes(fill = team_color, color = team_color2, size = plays), shape = 21, alpha = 0.9)+
  scale_color_identity(aesthetics =  c("fill","color"))+
  ggrepel::geom_text_repel(aes(label = name))+
  theme_bw()+
  geom_hline(yintercept = mean(qb_epa_play$epa_play), linetype = "dashed") +
  geom_vline(xintercept = mean(qb_epa_play$pass_rate), linetype = "dashed") +
  labs(x = "Pass Rate", y = "EPA/Play", title = "EPA/Play and Pass Rate Following Week 5",
       subtitle = "Minimum 50 Pass Attempts", caption = "Dotted Line Represents NFL Average                                       By Callan Capitolo | @CapAnalytics7")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 8))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 8))+
  theme(plot.title = element_text(size = 22, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 16, hjust = 0.5))
ggsave("Week5QBEfficiency23.png", width = 14, height =10, dpi = "retina")
#bar graph
qb_epa_play %>%
  ggplot(aes(x = epa_play, y = fct_reorder(name,epa_play)))+
  geom_bar(aes(fill = team_color, color = team_color2), stat = "identity", alpha = 0.9)+
  scale_color_identity(aesthetics = c("fill","color"))+
  theme_bw()+
  geom_image(aes(image = team_logo_espn, x = ifelse(epa_play>0, epa_play+0.01, epa_play-0.01)),
             asp = 16/9, size = 0.035)+
  labs(x = "EPA/Play", y = "Quarteback", title = "Quarterback's EPA/Play", caption = "Minimum 50 Attempts")+
  theme(panel.grid.major = element_blank())+
  theme(plot.title = element_text(size = 22, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 16, hjust = 0.5))
#table
qb_gt <- qb_epa_play %>%
  arrange(-epa_play) %>%
  mutate(rank = row_number()) %>%
  dplyr::select(rank,name, team_wordmark,pass_attempts,pass_rate, epa_play)%>%
  mutate(pass_rate = 100*round(pass_rate,3), epa_play = round(epa_play,2)) %>%
  gt() %>%
  cols_align(align = "center") %>%
  gtExtras::gt_img_rows(team_wordmark) %>%
  cols_label( rank = "Rank",
              name = "Quarterback",
              team_wordmark = "",
              pass_attempts = "Pass Attempts",
              pass_rate = "Pass Rate %",
              epa_play = "EPA Per Pass") %>%
  gtExtras::gt_theme_538() %>%
  gtExtras::gt_hulk_col_numeric(epa_play) 
 gtsave(qb_gt, "qb_efficiency_week12.png") 
  
