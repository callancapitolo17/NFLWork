if (!requireNamespace("plotly", quietly = TRUE)) {install.packages("plotly")}
if (!requireNamespace("nflfastR", quietly = TRUE)) {install.packages("nflfastR")}
if (!requireNamespace("tidyverse", quietly = TRUE)) {install.packages("tidyverse")}
if (!requireNamespace("gt", quietly = TRUE)) {install.packages("gt")}
if (!requireNamespace("dbscan", quietly = TRUE)) {install.packages("dbscan")}
if (!requireNamespace("zoo", quietly = TRUE)) {install.packages("zoo")}
if (!requireNamespace("gtExtras", quietly = TRUE)) {install.packages("gtExtras")}
library(plotly)
library(nflfastR)
library(tidyverse)
library(gt)
library(dbscan)
library(zoo)
library(gtExtras)
#loads play by play data
pbp <- load_pbp(c(2006:2023)) %>% 
  filter(pass == 1 | rush == 1) %>%
  filter(!is.na(epa),!is.na(posteam),!is.na(defteam))%>% 
  mutate(short_throw = ifelse(air_yards<=10,1,0),
         medium_throw = ifelse(air_yards>10&air_yards<=20,1,0),
         long_throw = ifelse(air_yards>20,1,0))
#creates offensive scouting report
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
    off_1h_epa = mean(epa[qtr<=2],na.rm =T),
    off_1h_pass_epa = mean(epa[qtr<=2&pass ==1], na.rm =T),
    off_1h_rush_epa = mean(epa[qtr<=2&rush ==1], na.rm =T),
    off_2h_epa = mean(epa[qtr>2],na.rm =T),
    off_2h_pass_epa = mean(epa[qtr>2&pass ==1], na.rm =T),
    off_2h_rush_epa = mean(epa[qtr>2&rush ==1], na.rm =T),
    off_scramble_rate = mean(qb_scramble[pass == 1]),
    off_scramble_epa = mean(epa[pass == 1 & qb_scramble == 1]),
    short_pass_rate = mean(short_throw[pass == 1],na.rm =T),
    medium_pass_rate = mean(medium_throw[pass == 1],na.rm =T),
    long_pass_rate = mean(long_throw[pass == 1],na.rm =T),
    off_short_pass_epa = mean(epa[short_throw==1],na.rm =T),
    off_medium_pass_epa = mean(epa[medium_throw ==1],na.rm =T),
    off_long_pass_epa = mean(epa[long_throw == 1],na.rm =T),
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
    home_off_epa = mean(epa[posteam == home_team]),
    home_pass_off_epa = mean(epa[posteam == home_team&pass == 1]),
    home_rush_off_epa = mean(epa[posteam == home_team&pass == 0]),
    away_off_epa = mean(epa[posteam == away_team]),
    away_pass_off_epa = mean(epa[posteam == away_team&pass == 1]),
    away_rush_off_epa = mean(epa[posteam == away_team&pass == 0])
  ) %>% 
  mutate(team = paste(posteam,season)) %>% 
  select(-season)

#creates defensive scouting report
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
    def_scramble_epa = mean(epa[pass == 1 & qb_scramble == 1]),
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
    home_def_epa = mean(epa[defteam == home_team]),
    home_pass_def_epa = mean(epa[defteam == home_team&pass == 1]),
    home_rush_def_epa = mean(epa[defteam == home_team&pass == 0]),
    away_def_epa = mean(epa[defteam == away_team]),
    away_pass_def_epa = mean(epa[defteam == away_team&pass == 1]),
    away_rush_def_epa = mean(epa[defteam == away_team&pass == 0])
  ) %>% 
  mutate(team = paste(defteam,season))


#joins scouting reports for both sides of the ball
scouting<- offensive_scouting %>% 
  left_join(defensive_scouting, by = c("team")) %>% 
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr"))



#creates the k-d tree for nearest neighbor search
similarity_scores<-kNN(scale(Filter(is.numeric,scouting) %>% select(-season) %>% 
                               mutate_all(~na.aggregate(.))),k = nrow(scouting)-1,sort = TRUE)


ui <- fluidPage(
  
  # Application title
  titlePanel("NFL Team Similarity Tool"),
  
  #team selection
  sidebarLayout(
    sidebarPanel(
      selectInput("teamcomp", "Select a Team Since the 2006 NFL Season to Compare", choices = scouting$team, selected = "SF 2023")),
    
    #create table
    mainPanel(
      gt_output("plot")
    )
  )
)


server <- function(input, output) {
  
  output$plot <- render_gt({
    top_10_indices <- similarity_scores$id[which(scouting$team == input$teamcomp),c(1:10)] #find 10 most similar teams indexes
    #create similarity score vector
    similarity <- c(1,1- similarity_scores$dist[which(scouting$team == input$teamcomp),c(1:10)]/max(similarity_scores$dist[which(scouting$team == input$teamcomp),]))
    #create table
    scouting[c(which(scouting$team == input$teamcomp),top_10_indices),] %>% 
      ungroup() %>% 
      mutate(rank = ifelse(row_number()==1,"",row_number()-1),
             sim_score = ifelse(similarity < 1,paste(round(similarity,4)*100,"%"),"")) %>%
      select(rank,sim_score,team_wordmark,season,off_epa, off_pass_epa, off_rush_epa, pass_rate, def_epa, def_pass_epa, def_rush_epa) %>% 
      mutate_if(is.numeric, ~ round(., 2)) %>% 
      gt() %>%
      tab_style(
        style = cell_fill(color = "yellow"),
        locations = cells_body(rows = 1)  # Specify the row to be styled
      ) %>% 
      gtExtras::gt_img_rows(team_wordmark) %>%
      cols_align(align = "center") %>%
      cols_label(rank = "Similarity Ranking",
                 sim_score = "Similarity Score",
                 team_wordmark = "Team",
                 season = "Season",
                 off_epa = "Off EPA/Play",
                 off_pass_epa = "Off EPA/Dropback",
                 off_rush_epa = "Off EPA/Rush",
                 pass_rate = "Pass Rate",
                 def_epa = "Def EPA/Play",
                 def_pass_epa = "Def EPA/Dropback",
                 def_rush_epa = "Def EPA/Rush") %>% 
      gtExtras::gt_theme_538() %>%
      tab_header(
        title = md(paste(c("10 Most Similar Teams to ", input$teamcomp)))
      )
    
      })
}

# Run the application 
shinyApp(ui = ui, server = server)
