library(plotly)
library(ggimage)
library(gt)
library(nflfastR)
library(tidyverse)

pbp <- load_pbp(2023)

pbp_rp <- pbp %>%
  filter(pass == 1 | rush == 1) %>%
  filter(!is.na(epa))

off_result_df <- data.frame()  # Initialize an empty list to store results

for (i in c(1:max(pbp_rp$week))) {
  result <- pbp_rp %>%
    filter(week %in% c(1:i)) %>% 
    group_by(posteam) %>%
    summarize(week = max(week), cum_off_epa = round(mean(epa,na.rm =T),3)) %>% 
    mutate(week = i) %>% 
  do(data.frame(.))
  off_result_df <- bind_rows(off_result_df, result)
}

def_result_df <- data.frame()  # Initialize an empty list to store results

for (i in c(1:max(pbp_rp$week))) {
  result <- pbp_rp %>%
    filter(week %in% c(1:i)) %>% 
    group_by(defteam) %>%
    summarize(cum_def_epa = round(mean(epa,na.rm =T),3)) %>% 
    mutate(week = i) %>% 
    do(data.frame(.))
  def_result_df <- bind_rows(def_result_df, result)
}
total_efficiency_both <- off_result_df %>%
  inner_join(def_result_df, by = c("posteam" = "defteam","week"))
total_efficiency_both <- total_efficiency_both %>%
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr"))

plot_ly(data = total_efficiency_both, x = ~cum_def_epa, y = ~cum_off_epa, frame = ~week,
        showlegend = F) %>% 
  add_text(
    text = ~posteam,  # Specify the text labels
    x = ~cum_def_epa, y = ~cum_off_epa,  # Coordinates for text annotations
    textposition = "top center",  # Adjust text position as needed
    showlegend = FALSE,  # Hide legend for text annotations
    textfont = list(color = ~team_color)
  ) %>%
  layout(xaxis = list(title = "Cumulative Defensive EPA/Play", autorange = "reversed", range = c(-0.3,0.3)), yaxis = list(title = "Cumulative Offensive EPA", range = c(-0.3,0.3)),
         title = list(text = "NFL Efficiency Landscape"),
         margin = list(l = 50, r = 10, b = 50, t = 50))  %>% 
  animation_opts(
    1250, easing = "cubic-in-out", redraw = FALSE
  )
  


pbp <- load_pbp(as.numeric(input$year))
pbp_rp <- pbp %>%
  filter(pass == 1 | rush == 1) %>%
  filter(!is.na(epa))

xside = "Defense"
yside = "Offense"
x_side <- ifelse(xside == "Defense","defteam","posteam")
y_side <- ifelse(yside == "Defense","defteam","posteam")
xstat = "epa"
ystat = "epa"
x_df <- data.frame()  # Initialize an empty list to store results

for (i in c(1:max(pbp_rp$week))) {
  result <- pbp_rp %>%
    filter(week %in% c(1:i)) %>% 
    group_by(!!sym(x_side)) %>%
    summarize(week = max(week), cum_off_epa = round(mean(!!sym(xstat),na.rm =T),3)) %>% 
    mutate(week = i) %>% 
    do(data.frame(.))
  x_df <- bind_rows(x_df, result)
}
x_df <- x_df %>% 
  rename("team" = x_side)

y_df <- data.frame()  # Initialize an empty list to store results

for (i in c(1:max(pbp_rp$week))) {
  result <- pbp_rp %>%
    filter(week %in% c(1:i)) %>% 
    group_by(!!sym(y_side)) %>%
    summarize(cum_def_epa = round(mean(!!sym(ystat),na.rm =T),3)) %>% 
    mutate(week = i) %>% 
    do(data.frame(.))
  y_df <- bind_rows(y_df, result)
}

y_df <- y_df %>% 
  rename("team" = y_side)
total_efficiency_both <- x_df %>%
  inner_join(y_df, by = c("team", "week"))


total_efficiency_both <- total_efficiency_both %>%
  left_join(teams_colors_logos, by = c("team" = "team_abbr"))

plot_ly(data = total_efficiency_both, x = ~cum_def_epa, y = ~cum_off_epa, frame = ~week,
        showlegend = F) %>% 
  add_text(
    text = ~posteam,  # Specify the text labels
    x = ~cum_def_epa, y = ~cum_off_epa,  # Coordinates for text annotations
    textposition = "top center",  # Adjust text position as needed
    showlegend = FALSE,  # Hide legend for text annotations
    textfont = list(color = ~team_color)
  ) %>%
  layout(xaxis = list(title = "Cumulative Defensive EPA/Play", autorange = "reversed"), 
         yaxis = list(title = "Cumulative Offensive EPA/Play"),
         title = list(text = "NFL Efficiency Landscape"),
         margin = list(l = 50, r = 10, b = 50, t = 50))  %>% 
  animation_opts(
    1250, easing = "cubic-in-out", redraw = FALSE
)
  