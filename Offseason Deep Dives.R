library(ggimage)
library(gt)
library(nflfastR)
library(tidyverse)
library(png)
library(nflreadr)
library(gtExtras)

# Kyler Murray ----
card <- load_pbp(2019:2023) %>% 
  filter(posteam == "ARI" | defteam == "ARI")

card_pbp <- card %>%
  filter(pass == 1 | rush == 1) %>%
  filter(!is.na(epa))

cards_results<-card %>% 
  group_by(game_id) %>% 
  summarize(hometeam = max(home_team),homescore = max(home_score),awayteam = max(away_team),awayscore = max(away_score)) %>% 
  mutate(game_result = ifelse((hometeam == "ARI"&homescore>awayscore) | (awayteam == "ARI" & awayscore>homescore), "W",ifelse(homescore == awayscore,"T","L")))

test_data <- card %>% 
  left_join(cards_results, by = "game_id")

murray_eff<-test_data %>%
  filter(qb_kneel == 0) %>% 
  group_by(game_result) %>%
  summarize(rushes = mean((qb_scramble == 1 | rush == 1) & id == "00-0035228") ,avg_rush_epa = mean(epa[(qb_scramble == 1 | rush == 1) & id == "00-0035228"],na.rm = T), avg_pass_epa =  mean(epa[(qb_scramble == 0 & pass == 1) & passer_id == "00-0035228"],na.rm = T)) %>%
  filter(game_result != "T") %>% 
  mutate_if(is.numeric, ~round(.,2)) %>% 
  gt() %>% 
  cols_align(align = "center") %>%
  cols_label(game_result = "Game Result",rushes = "Rush Rate", avg_rush_epa = "EPA/Rush", avg_pass_epa = "EPA/Pass") %>% 
  gtExtras::gt_theme_538() %>%
  tab_header(title = md("Kyler Murray Efficiency by Game Result"), subtitle = md("Murray has drastic splits in EPA/Pass in wins vs losses"))
gtsave(murray_eff, "murraybreakdowwn.png")


#kyler actually has better rush performance in losses, but the splits between his pass_epa win vs loss is insane
#additionally and interestingly kyler's rush rate increases in wins despite less efficiency

murray_rush<- test_data %>%
  filter(qb_kneel == 0) %>% 
  group_by(game_id) %>%
  summarize(rushes = mean((qb_scramble == 1 | rush == 1) & id == "00-0035228") ,avg_rush_epa = mean(epa[(qb_scramble == 1 | rush == 1) & id == "00-0035228"],na.rm = T), avg_pass_epa =  mean(epa[(qb_scramble == 0 & pass == 1) & passer_id == "00-0035228"],na.rm = T))
rush_results <- left_join(murray_rush,cards_results, by = c("game_id"))


rush_results %>% 
  filter(avg_rush_epa> -3.5) %>% 
  ggplot(aes(x = avg_rush_epa, y = avg_pass_epa, label = game_result, color = game_result))+
  geom_text(size = 5)+
  theme(legend.position = "none")+
  labs(x = "EPA/Rush (Includes QB Scrambles)", y = "EPA/Pass", title = "EPA/Pass vs EPA Rush EPA by Game Result for Kyler Murray",
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR", subtitle = "Kyler Murray's Passing Effiency is More Impactful on Results Than His Rushing Efficiency")
ggsave("MurrayScatter.png", width = 14, height =10, dpi = "retina")

#Eagles Rushing Look ----
eagles <- load_pbp(2023) %>% 
  filter(rush == 1)

eagles %>% 
  filter(rusher_player_name != "J.Hurts") %>% 
  group_by(posteam) %>% 
  summarize(epa_rush = mean(epa), success_rate = mean(success)) %>% 
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr")) %>% 
  ggplot(aes(x = success_rate, y = epa_rush))+
  geom_image(aes(image = team_logo_espn), size = 0.05, asp = 16/9)+
  labs(x = "Success Rate", y = "EPA/Rush",title = "EPA/Rush vs Success Rate (No Hurts Rushes)")

# Home Field Advantage Breakdown----

football <- load_pbp(c(1999:2023))

game_by_game <- football %>% 
  # filter(season == 2023) %>% 
  group_by(game_id) %>% 
  summarise(homescore = max(home_score),hometeam = max(home_team),awayscore = max(away_score),awayteam = max(away_team),
            year = max(season), field = max(game_stadium), field_id = max(stadium_id), loc = max(location)) %>% 
  filter(loc != "Neutral") %>% 
  mutate(total_points = homescore+awayscore)


home_scoring <- game_by_game %>% 
  group_by(hometeam,field_id) %>% 
  summarize(avg_home_game_points = mean(total_points),lastyear = max(year),name = max(field), num_games = n())  # filter(num_games>10)

away_scoring <- game_by_game %>% 
  left_join(home_scoring %>% select(hometeam, lastyear), by = c("awayteam" = "hometeam")) %>% 
  mutate(correct = ifelse(year <= lastyear, "yes","no")) %>%
  filter(correct == "yes") %>% 
  group_by(game_id) %>% 
  summarize(min_year=min(lastyear),awayteam = max(awayteam), awayscore = max(awayscore), 
            hometeam = max(hometeam), homescore = max(homescore)) %>% 
  mutate(totalpoints = homescore+awayscore) %>% 
  group_by(awayteam,min_year) %>% 
  summarize(avg_away_score = mean(totalpoints))

parkfactor <- home_scoring %>% 
  left_join(away_scoring, by = c("hometeam" = "awayteam", "lastyear" = "min_year")) %>% 
  filter(num_games>=10) %>% 
  mutate(stadium_factor = avg_home_game_points/avg_away_score, stadium_id = paste(hometeam,name,lastyear, sep = " ")) %>% 
  left_join(teams_colors_logos, by = c("hometeam" = "team_abbr"))

parkfactor %>%
  ggplot(aes(x = reorder(stadium_id, stadium_factor), y = stadium_factor, fill = team_color))+
  geom_bar(stat = "identity")+
  geom_image(aes(image = team_logo_espn), size = 0.03, asp = 16/9)+
  labs(x = "Stadium (Team, Stadium Name, Last Year With Games)", y = "Stadium Factor (Total Points Per Home Game/Total Points Per Away Game)", title = "What Stadium is it Hardest to Score Points at?",
       subtitle = "Minimum of 10 Games Played at the Stadium Since 1999", caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")+
  theme(axis.text.x = element_text(angle = 65, hjust = 1))+
  scale_fill_identity()
ggsave("StadiumFactor.png", width = 14, height =10, dpi = "retina")


#below is probably useless

factorbyteam<- game_by_game %>% 
  group_by(awayteam) %>% 
  summarize(avg_away_game_points = mean(total_points)) %>% 
  left_join(home_scoring,by = c("awayteam" = "hometeam")) %>% 
  mutate(stadium_factor = avg_home_game_points/avg_away_game_points)
#look to input last year of stadium use into away data, then group_by(away_team,year)
#join on stadium id and the only other thing be select(id,year)

factorbystadium<- game_by_game %>% 
  group_by(awayteam,field) %>% 
  summarize(avg_away_game_points = mean(total_points)) %>% 
  left_join(home_scoring,by = c("awayteam" = "hometeam", "field")) %>% 
  mutate(stadium_factor = avg_home_game_points/avg_away_game_points)

## 2023 Home Field ----
game_by_game <- football %>% 
  filter(season == 2023) %>% 
  group_by(game_id) %>% 
  summarise(homescore = max(home_score),hometeam = max(home_team),awayscore = max(away_score),awayteam = max(away_team),
            year = max(season), field = max(game_stadium), field_id = max(stadium_id), loc = max(location)) %>% 
  filter(loc != "Neutral") %>% 
  mutate(total_points = homescore+awayscore)


home_scoring <- game_by_game %>% 
  group_by(hometeam,field_id) %>% 
  summarize(avg_home_game_points = mean(total_points),lastyear = max(year),name = max(field), num_games = n())  # filter(num_games>10)

away_scoring <- game_by_game %>% 
  left_join(home_scoring %>% select(hometeam, lastyear), by = c("awayteam" = "hometeam")) %>% 
  mutate(correct = ifelse(year <= lastyear, "yes","no")) %>%
  filter(correct == "yes") %>% 
  group_by(game_id) %>% 
  summarize(min_year=min(lastyear),awayteam = max(awayteam), awayscore = max(awayscore), 
            hometeam = max(hometeam), homescore = max(homescore)) %>% 
  mutate(totalpoints = homescore+awayscore) %>% 
  group_by(awayteam,min_year) %>% 
  summarize(avg_away_score = mean(totalpoints))

parkfactor <- home_scoring %>% 
  left_join(away_scoring, by = c("hometeam" = "awayteam", "lastyear" = "min_year")) %>% 
  mutate(stadium_factor = avg_home_game_points/avg_away_score, stadium_id = paste(hometeam,name,lastyear, sep = " ")) %>% 
  left_join(teams_colors_logos, by = c("hometeam" = "team_abbr"))

parkfactor %>%
  ggplot(aes(x = reorder(stadium_id, stadium_factor), y = stadium_factor, fill = team_color))+
  geom_bar(stat = "identity")+
  geom_image(aes(image = team_logo_espn), size = 0.03, asp = 16/9)+
  labs(x = "Stadium (Team, Stadium Name)", y = "Stadium Factor (Total Points Per Home Game/Total Points Per Away Game)", 
       title = "What Stadium was it Hardest to Score Points at in 2023?",
       subtitle = "FirstEnergy was the Hardest Stadium to Score at, while Lincoln Financial was the Easiest According to Stadium Factor", caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")+
  theme(axis.text.x = element_text(angle = 65, hjust = 1))+
  scale_fill_identity()
ggsave("StadiumFactor23.png", width = 14, height =10, dpi = "retina")


## Conditions Factor ----
library(stringr)
library(rsample)      # data splitting 
library(gbm)          # basic implementation
library(xgboost)      # a faster implementation of gbm
library(caret)        # an aggregator package for performing many machine learning models
library(h2o)          # a java-based platform
library(pdp)          # model visualization
library(lime)         # model visualization
pbp15 <- load_pbp(c(2015:2023))


categories <- c(
  snow = c("snow", "snowy"),
  rainy = c("rain", "rainy", "showers", "thunderstorms"),
  clear = c("clear", "sunny", "fair", "warm", "hot", "no", "0%"),
  cloudy = c("cloudy", "overcast", "partly cloudy", "mostly cloudy")
)

# Function to classify values into categories
classify_weather <- function(value) {
  for (category in names(categories)) {
    if (any(sapply(categories[[category]], grepl, value, ignore.case = TRUE))) {
      category_cleaned <- gsub("\\d+$", "", category)
      return(category)
    }
  }
  return("other")  # Default category if not matched
}


weather_data <- pbp15 %>% 
  filter(roof == "outdoors", season < 2023) %>% 
  group_by(game_id) %>% 
  summarize(weather = max(weather), dome = max(roof), temperature = max(temp), windmph = max(wind), bet_total = max(total_line), surface = max(surface),
            game_total = max(total), pass_epa = mean(epa[qb_scramble == 0 & pass_attempt == 1], na.rm = T), rush_epa = mean(epa[qb_scramble==1 | rush == 1],na.rm =T)) %>% 
  mutate(weather_desc2 = tolower(str_extract(weather, ".*(?=T)")),
         clean_weather = as.factor(gsub("\\d+$", "",sapply(weather_desc2, classify_weather))),
         # wind_direction = ifelse(grepl("from", tolower(weather)),str_extract(weather, "(?<=Wind: From )[A-Z]+"),str_extract(weather, "(?<=Wind: )[A-Z]+")),
         # clean_wind_direction = ifelse(grepl("EAST",wind_direction), "E", ifelse(grepl("WE", wind_direction),"W",wind_direction)),
         humidity = as.numeric(str_match(weather, "Humidity: (\\d+)%")[, 2]),
         total_diff = game_total - bet_total,
         clean_surface = as.factor(ifelse(grepl("turf",surface),"turf","grass"))) %>% 
  ungroup() %>% 
  select(-weather,-game_id,-dome,-weather_desc2,-surface)

pairs(select_if(weather_data,is.numeric))
cor(na.omit(select_if(weather_data,is.numeric)))

test2<-weather_data %>% 
  group_by(clean_wind_direction) %>% 
  summarize(bet = mean(bet_total), game = mean(game_total), pass = mean(pass_epa),
            rush = mean(rush_epa), diff = mean(total_diff))

weather_data %>% 
  filter(windmph<60) %>% 
  ggplot(aes(x = bet_total, y = rush_epa, size= windmph))+
  geom_point()+
  geom_smooth(method = "lm")

table(test$clean_weather)

passing_data <- weather_data %>% 
  select(-rush_epa, -game_total, -total_diff)

pairs(select_if(passing_data,is.numeric))
cor(na.omit(select_if(passing_data,is.numeric)))

## modeling

passing_model <- lm(pass_epa ~ ., data = passing_data)
summary(passing_model)

set.seed(123)
pass_split <- initial_split(na.omit(passing_data), prop = .7)
pass_train <- training(pass_split)
pass_test  <- testing(pass_split)

features <- setdiff(names(pass_train),"pass_epa")

treatplan <- vtreat::designTreatmentsZ(pass_train, features, verbose = FALSE)

new_vars <- treatplan %>%
  magrittr::use_series(scoreFrame) %>%        
  dplyr::filter(code %in% c("clean", "lev")) %>% 
  magrittr::use_series(varName)     

features_train <- vtreat::prepare(treatplan, pass_train, varRestriction = new_vars) %>% as.matrix()
response_train <- pass_train$pass_epa

features_test <- vtreat::prepare(treatplan, pass_test, varRestriction = new_vars) %>% as.matrix()
response_test <- pass_test$pass_epa

set.seed(123)



hyper_grid <- expand.grid(
  eta = c(.01, .05, .1, .3),
  max_depth = c(1, 3, 5, 7),
  min_child_weight = c(1, 3, 5, 7),
  subsample = c(.65, .8, 1), 
  colsample_bytree = c(.8, .9, 1),
  optimal_trees = 0,               # a place to dump results
  min_RMSE = 0                     # a place to dump results
)

for(i in 1:nrow(hyper_grid)) {
  
  # create parameter list
  params <- list(
    eta = hyper_grid$eta[i],
    max_depth = hyper_grid$max_depth[i],
    min_child_weight = hyper_grid$min_child_weight[i],
    subsample = hyper_grid$subsample[i],
    colsample_bytree = hyper_grid$colsample_bytree[i]
  )
  
  # reproducibility
  set.seed(123)
  
  # train model
  xgb.tune <- xgb.cv(
    params = params,
    data = features_train,
    label = response_train,
    nrounds = 5000,
    nfold = 5,
    objective = "reg:linear",  # for regression models
    verbose = 0,               # silent,
    early_stopping_rounds = 10 # stop if no improvement for 10 consecutive trees
  )
  
  # add min training error and trees to grid
  hyper_grid$optimal_trees[i] <- which.min(xgb.tune$evaluation_log$test_rmse_mean)
  hyper_grid$min_RMSE[i] <- min(xgb.tune$evaluation_log$test_rmse_mean)
}


hyper_grid %>%
  dplyr::arrange(min_RMSE) %>%
  head(10)

params <- list(
  eta = 0.1,
  max_depth = 1,
  min_child_weight = 7,
  subsample = 0.65,
  colsample_bytree = 0.8
)

xgb.final <- xgb.cv(
  params = params,
  data = features_train,
  label = response_train,
  nrounds = 50,
  nfold = 5,
  objective = "reg:linear",  # for regression models
  verbose = 0               # silent
  
)

weather_data <- pbp15 %>% 
  filter(roof == "outdoors", season == 2023) %>% 
  group_by(game_id) %>% 
  summarize(weather = max(weather), dome = max(roof), temperature = max(temp), windmph = max(wind), bet_total = max(total_line), surface = max(surface),
            game_total = max(total), pass_epa = mean(epa[qb_scramble == 0 & pass_attempt == 1], na.rm = T), rush_epa = mean(epa[qb_scramble==1 | rush == 1],na.rm =T)) %>% 
  mutate(weather_desc2 = tolower(str_extract(weather, ".*(?=T)")),
         clean_weather = as.factor(gsub("\\d+$", "",sapply(weather_desc2, classify_weather))),
         # wind_direction = ifelse(grepl("from", tolower(weather)),str_extract(weather, "(?<=Wind: From )[A-Z]+"),str_extract(weather, "(?<=Wind: )[A-Z]+")),
         # clean_wind_direction = ifelse(grepl("EAST",wind_direction), "E", ifelse(grepl("WE", wind_direction),"W",wind_direction)),
         humidity = as.numeric(str_match(weather, "Humidity: (\\d+)%")[, 2]),
         total_diff = game_total - bet_total,
         clean_surface = as.factor(ifelse(grepl("turf",surface),"turf","grass"))) %>% 
  ungroup() %>% 
  select(-weather,-game_id,-dome,-weather_desc2,-surface)


## Important Metrics ----
nfl99all <- load_pbp(1999:2023)
nfl99 <- nfl99all %>% 
  filter(pass == 1 | rush == 1) %>% 
  mutate(explosive = ifelse((yards_gained>20 & pass_attempt == 1) | (yards_gained >12 & (qb_scramble == 1 | rush == 1)),1,0),
         negative = ifelse(yards_gained < 0, 1,0))
metrics <- nfl99 %>% 
  filter(qb_kneel == 0, qb_spike == 0, pass == 1 | rush == 1) %>% 
  group_by(game_id) %>% 
  summarize(homescore = max(home_score),homeepa = mean(epa[home_team == posteam], na.rm = T), homepassepa = mean(epa[home_team == posteam & pass == 1], na.rm = T), homerushepa = mean(epa[home_team == posteam & rush == 1], na.rm = T),
            homesuccess = mean(success[home_team == posteam], na.rm = T), homepasssuccess = mean(success[home_team == posteam & pass == 1], na.rm = T), homerushsuccess = mean(success[home_team == posteam & rush == 1], na.rm = T), 
            homeyards = mean(yards_gained[home_team == posteam], na.rm = T), homepassyards = mean(yards_gained[home_team == posteam & pass == 1], na.rm = T), homerushyards = mean(yards_gained[home_team == posteam & rush == 1], na.rm = T),
            awayscore = max(away_score),awayepa = mean(epa[away_team == posteam], na.rm = T), awaypassepa = mean(epa[away_team == posteam & pass == 1], na.rm = T), awayrushepa = mean(epa[away_team == posteam & rush == 1], na.rm = T),
            awaysuccess = mean(success[away_team == posteam], na.rm = T), awaypasssuccess = mean(success[away_team == posteam & pass == 1], na.rm = T), awayrushsuccess = mean(success[away_team == posteam & rush == 1], na.rm = T), 
            awayyards = mean(yards_gained[away_team == posteam], na.rm = T), awaypassyards = mean(yards_gained[away_team == posteam & pass == 1], na.rm = T), awayrushyards = mean(yards_gained[away_team == posteam & rush == 1], na.rm = T),
            year = max(season)) %>% 
  mutate(epahigherwin = ifelse(homescore == awayscore,NA,ifelse((homescore>awayscore & homeepa>awayepa)|(homescore<awayscore & homeepa<awayepa),1,0)),
         epapasshigherwin = ifelse(homescore == awayscore,NA,ifelse((homescore>awayscore & homepassepa>awaypassepa)|(homescore<awayscore & homepassepa<awaypassepa),1,0)),
         eparushhigherwin = ifelse(homescore == awayscore,NA,ifelse((homescore>awayscore & homerushepa>awayrushepa)|(homescore<awayscore & homerushepa<awayrushepa),1,0)),
         successhigherwin = ifelse(homescore == awayscore,NA,ifelse((homescore>awayscore & homesuccess>awaysuccess)|(homescore<awayscore & homesuccess<awaysuccess),1,0)),
         successpasshigherwin = ifelse(homescore == awayscore,NA,ifelse((homescore>awayscore & homepasssuccess>awaypasssuccess)|(homescore<awayscore & homepasssuccess<awaypasssuccess),1,0)),
         successrushhigherwin = ifelse(homescore == awayscore,NA,ifelse((homescore>awayscore & homerushsuccess>awayrushsuccess)|(homescore<awayscore & homerushsuccess<awayrushsuccess),1,0)),
         yardshigherwin = ifelse(homescore == awayscore,NA,ifelse((homescore>awayscore & homeyards>awayyards)|(homescore<awayscore & homeyards<awayyards),1,0)),
         yardspasshigherwin = ifelse(homescore == awayscore,NA,ifelse((homescore>awayscore & homepassyards>awaypassyards)|(homescore<awayscore & homepassyards<awaypassyards),1,0)),
         yardsrushhigherwin = ifelse(homescore == awayscore,NA,ifelse((homescore>awayscore & homerushyards>awayrushyards)|(homescore<awayscore & homerushyards<awayrushyards),1,0))) %>% 
  group_by(year) %>% 
  summarize(pct_epa_high_win = mean(epahigherwin,na.rm = T),
            pct_epa_pass_high_win = mean(epapasshigherwin,na.rm = T),
            pct_epa_rush_high_win = mean(eparushhigherwin,na.rm = T),
            pct_success_high_win = mean(successhigherwin,na.rm = T),
            pct_success_pass_high_win = mean(successpasshigherwin,na.rm = T),
            pct_success_rush_high_win = mean(successrushhigherwin,na.rm = T),
            pct_yards_high_win = mean(yardshigherwin,na.rm = T),
            pct_yards_pass_high_win = mean(yardspasshigherwin,na.rm = T),
            pct_yards_rush_high_win = mean(yardsrushhigherwin,na.rm = T)) %>% 
  pivot_longer(cols = starts_with("pct") ,names_to = "type",values_to = "pct_winners")

metrics %>% 
  filter(type %in% c("pct_epa_high_win","pct_success_high_win","pct_yards_high_win")) %>% 
  ggplot(aes(x = year, y = pct_winners, color = type)) +
  geom_line(lwd = 1.5)+
  labs(x = "Season", y = "Rate of Game Winners that Outperformed Opponent in Given Metric",
       title = "Which Offensive Metric Gives the Best Information About who the Game Winner Will be?",
       color = "Offensive Metric", subtitle = "Teams that outperformed their opponent in EPA/Play have won over 90% of their games in every season since 1999")+
  scale_color_manual(values = c("pct_epa_high_win" = "red", "pct_success_high_win"= "blue",
                                "pct_yards_high_win" = "black"),
                     labels = c("pct_epa_high_win" = "EPA/Play", "pct_success_high_win"= "Success Rate",
                                "pct_yards_high_win" = "Yards/Play"))
ggsave("Metrics.png", width = 14, height =10, dpi = "retina")


metrics %>% 
  filter(type %in% c("pct_epa_pass_high_win","pct_success_pass_high_win","pct_yards_pass_high_win")) %>% 
  ggplot(aes(x = year, y = pct_winners, color = type)) +
  geom_line(lwd = 1.5)+
  labs(x = "Season", y = "Rate of Game Winners that Outperformed Opponent in Given Metric",
       title = "Which Offensive Metric Gives the Best Information About who the Game Winner Will be?",
       color = "Offensive Metric", subtitle = "Teams that outperformed their opponent in EPA/Dropback have won over 80% of their games in 22 out of the last 24 seasons")+
  scale_color_manual(values = c("pct_epa_pass_high_win" = "red", "pct_success_pass_high_win"= "blue",
                                "pct_yards_pass_high_win" = "black"),
                     labels = c("pct_epa_pass_high_win" = "EPA/Dropback", "pct_success_pass_high_win"= "Dropback Success Rate",
                                "pct_yards_pass_high_win" = "Yards/Dropback"))
ggsave("PassMetrics.png", width = 14, height =10, dpi = "retina")


metrics %>% 
  filter(type %in% c("pct_epa_rush_high_win","pct_success_rush_high_win","pct_yards_rush_high_win")) %>% 
  ggplot(aes(x = year, y = pct_winners, color = type)) +
  geom_line(lwd = 1.5)+
  labs(x = "Season", y = "Rate of Game Winners that Outperformed Opponent in Given Metric",
       title = "Which Offensive Metric Gives the Best Information About who the Game Winner Will be?",
       color = "Offensive Metric", subtitle = "EPA/Rush only has a slight advantage over the other metrics in determining game winner")+
  scale_color_manual(values = c("pct_epa_rush_high_win" = "red", "pct_success_rush_high_win"= "blue",
                                "pct_yards_rush_high_win" = "black"),
                     labels = c("pct_epa_rush_high_win" = "EPA/Rush", "pct_success_rush_high_win"= "Rush Success Rate",
                                "pct_yards_rush_high_win" = "Yards/Rush"))
ggsave("RushMetrics.png", width = 14, height =10, dpi = "retina")
#from this lens who are the luckiest and unluckiest teams

metrics %>% 
  filter(type %in% c("pct_epa_high_win","pct_epa_pass_high_win","pct_epa_rush_high_win")) %>% 
  ggplot(aes(x = year, y = pct_winners, color = type)) +
  geom_line(lwd = 1.5)+
  labs(x = "Season", y = "Rate of Game Winners that Outperformed Opponent in Given Metric",
       title = "Which Offensive EPA Metric Gives the Best Information About who the Game Winner Will be?",
       color = "Offensive Metric", subtitle = "EPA/Play is by far the best metric in determing winners")+
  scale_color_manual(values = c("pct_epa_high_win" = "red", "pct_epa_pass_high_win"= "blue",
                                "pct_epa_rush_high_win" = "black"),
                     labels = c("pct_epa_high_win" = "EPA/Play", "pct_epa_pass_high_win"= "EPA/Dropback",
                                "pct_epa_rush_high_win" = "EPA/Rush"))
ggsave("EPA.png", width = 14, height =10, dpi = "retina")

metrics %>% 
  filter(type %in% c("pct_success_high_win","pct_success_pass_high_win","pct_success_rush_high_win")) %>% 
  ggplot(aes(x = year, y = pct_winners, color = type)) +
  geom_line(lwd = 1.5)+
  labs(x = "Season", y = "Rate of Game Winners that Outperformed Opponent in Given Metric",
       title = "Which Offensive Success Rate Metric Gives the Best Information About who the Game Winner Will be?",
       color = "Offensive Metric", subtitle = "Dropback success rate and overall success rate are by far the biggest factors in determining game winner")+
  scale_color_manual(values = c("pct_success_high_win" = "red", "pct_success_pass_high_win"= "blue",
                                "pct_success_rush_high_win" = "black"),
                     labels = c("pct_success_high_win" = "Success Rate", "pct_success_pass_high_win"= "Pass Success Rate",
                                "pct_success_rush_high_win" = "Rush Success Rate"))
ggsave("Success.png", width = 14, height =10, dpi = "retina")

metrics %>% 
  filter(type %in% c("pct_yards_high_win","pct_yards_pass_high_win","pct_yards_rush_high_win")) %>% 
  ggplot(aes(x = year, y = pct_winners, color = type)) +
  geom_line(lwd = 1.5)+
  labs(x = "Season", y = "Rate of Game Winners that Outperformed Opponent in Given Metric",
       title = "Which Offensive Yards Metric Gives the Best Information About who the Game Winner Will be?",
       color = "Offensive Metric", subtitle = "Yards/dropback is more important factor in determining game winners than yards/rush and yards/play")+
  scale_color_manual(values = c("pct_yards_high_win" = "red", "pct_yards_pass_high_win"= "blue",
                                "pct_yards_rush_high_win" = "black"),
                     labels = c("pct_yards_high_win" = "Yards/Play", "pct_yards_pass_high_win"= "Yards/Dropback",
                                "pct_yards_rush_high_win" = "Yards/Rush"))
ggsave("Yards.png", width = 14, height =10, dpi = "retina")

# Luckiness ----
test <- nfl99 %>% 
  filter(qb_kneel == 0, qb_spike == 0,pass == 1 |rush  == 1) %>% 
  group_by(game_id) %>% 
  summarize(homescore = max(home_score),homeepa = mean(epa[home_team == posteam], na.rm = T), hometeam = max(home_team),
            awayscore = max(away_score),awayepa = mean(epa[away_team == posteam], na.rm = T), awayteam = max(away_team), year = max(season)) %>% 
  mutate(winningteam = ifelse(homescore == awayscore, NA,ifelse(homescore>awayscore,hometeam,awayteam)),
         losingteam = ifelse(homescore == awayscore, NA,ifelse(homescore>awayscore,awayteam,hometeam))) %>% 
  mutate(luckywin = ifelse(homescore == awayscore,NA,ifelse((homescore>awayscore & homeepa<awayepa)|(homescore<awayscore & homeepa>awayepa),1,0))) %>% 
  group_by(losingteam,year) %>% 
  summarize(pct_lucky_win = mean(luckywin,na.rm = T),
            total_lucky_win = sum(luckywin,na.rm = T))


# Biggest Change in offensive EPA----

epa99<- nfl99 %>% 
  filter(rush == 1 | pass == 1) %>% 
  group_by(posteam,season) %>%
  summarize(epa_play = mean(epa,na.rm = T), epa_pass = mean(epa[pass == 1],na.rm = T))

teamlist <- na.omit(unique(epa99$posteam))
epa99$epachange <- NA
for(team in teamlist){
  for(x in c((min(epa99[which(epa99$posteam == team),"season"])+1):max(epa99[which(epa99$posteam == team),"season"]))){
    index<- which(epa99$posteam == team & epa99$season == x)
    priorindex<- which(epa99$posteam == team & epa99$season == x-1)
    epa99[index,"epachange"]<- epa99[[index,"epa_play"]]-epa99[[ifelse(length(priorindex) == 0,index,priorindex),"epa_play"]]
  }
}

top10change<-epa99 %>% 
  arrange(epachange) %>% 
  head(.,10) %>% 
  ungroup() %>% 
  mutate(changes = c("Hired Ron Rivera, Drafted Cam Newton","Hired Kliff Kingsbury, Drafted Kyler Murray","Hired Sean McVay","Hired Bruce Arians, Traded for Carson Palmer",
                     "Matt Ryan MVP Season","Hired Mike Smith, Drafted Matt Ryan","Replaced OC with Mike Martz, Mike Nolan Fired Midseason","Drafted Dak Prescott and Ezekiel Elliot","Josh Freeman full time starter",
                     "Hired Mike McCoy")) %>% 
  left_join(teams_colors_logos,by = c("posteam" = "team_abbr"))


top10change %>% 
  select(team_wordmark,season,epa_play,epachange,changes) %>% 
  mutate_if(is.numeric,~round(.,2)) %>% 
  gt() %>% 
  cols_align(align = "center") %>% 
  gtExtras::gt_img_rows(team_wordmark) %>% 
  cols_label(team_wordmark = "Team", season = "Season", epa_play = "EPA/Play",
             epachange = "EPA/Play Improvement", changes = "Organization Changes") %>% 
  gtExtras::gt_theme_538() %>% 
  tab_header(title = md("Top 10 Improvement in EPA/Play from Prior Season Since 1999"))
gtsave(offchange,"epachange.png")

logos <- epa99 %>% 
  filter(!is.na(posteam)) %>% 
  left_join(teams_colors_logos,by = c("posteam" = "team_abbr")) 

image_url <- logos$team_logo_espn[which(logos$posteam == "CHI" & logos$season == 2003)]

# Download the image
download.file(image_url, destfile = "logo.png", mode = "wb")

# Read the downloaded image
logo <- readPNG("logo.png")

epa99 %>% 
  filter(!is.na(posteam)) %>% 
  left_join(teams_colors_logos,by = c("posteam" = "team_abbr")) %>% 
  ggplot(aes(x = season, y = epa_pass, color = team_color))+
  geom_line(aes(alpha = ifelse(posteam == "CHI", 1, 0.99)), lwd = 1.5)+
  scale_color_identity()+guides(alpha = FALSE)+
  # geom_vline(xintercept = 2002, linetype = "dashed")+
  # annotate("text", x = 2003, y = 0.2, label = "Mariucci \n Final Season", color = "black",size =4)+
  # geom_vline(xintercept = 2005, linetype = "dashed")+
  # annotate("text", x = 2006, y = 0.2, label = "Nolan \n 1st Season", color = "black",size =4)+
  # geom_vline(xintercept = 2008, linetype = "dashed")+
  # annotate("text", x = 2009.25, y = 0.2, label = "Nolan \n Fired Midseason", color = "black",size =4)+
  # geom_vline(xintercept = 2011, linetype = "dashed")+
  # annotate("text", x = 2012, y = 0.2, label = "Harbaugh \n 1st Season", color = "black",size =4)+
  # geom_vline(xintercept = 2014, linetype = "dashed")+
  # annotate("text", x = 2015, y = 0.2, label = "Harbaugh \n Final Season", color = "black",size =4)+
  # geom_vline(xintercept = 2017, linetype = "dashed")+
  # annotate("text", x = 2018, y = 0.2, label = "Shanahan \n First Season", color = "black",size =4)+
  labs(y = "Offensive EPA/Dropback", x = "Season", title = "Will Caleb Williams Finally be the Answer at QB the Bears Have Been Looking for to Overcome Their Historical Struggles in the Passing Game?",
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")+
  annotation_custom(grid::rasterGrob(logo), ymin = -0.04, ymax = -0.1, xmin = 2022.5, xmax = 2024) +
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())
ggsave("BearsOffense.png", width = 14, height =10, dpi = "retina")

epa99 %>% filter(posteam == "NE")

#Clutch ----
clutch23 <- nfl99 %>% 
  mutate(pointdiff = abs(total_home_score - total_away_score)) %>% 
  filter(season >=2020, pass == 1 | rush == 1, (qtr == 4 & quarter_seconds_remaining < 300)| game_seconds_remaining < 300, pointdiff <= 8)

off_clutch23 <- clutch23 %>% 
  group_by(posteam) %>% 
  summarize(epa_off = mean(epa,na.rm = T), off_plays = n())

def_clutch23 <- clutch23 %>% 
  group_by(defteam) %>% 
  summarize(epa_def = mean(epa,na.rm = T), def_plays = n())
off_clutch23 %>% 
  left_join(def_clutch23, by = c("posteam" = "defteam")) %>% 
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr")) %>% 
  ggplot(aes(x = epa_def, y = epa_off))+
  geom_point()+
  geom_image(aes(image = team_logo_espn), size = 0.05, asp = 16/9)+
  labs(x = "Clutch Defensive EPA/Play", y = "Clutch Offensive EPA/Play", title = "How Did Teams Perform in the Clutch Since 2020?",
       subtitle = "The Lions were the clutch kings, the Ravens defense seems to implode late game", caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")+ 
  scale_x_reverse()+
  theme_minimal()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
ggsave("Clutch23.png", width = 14, height =10, dpi = "retina")

coach_clutch <- nfl99 %>% 
  mutate(pointdiff = abs(total_home_score - total_away_score)) %>% 
  filter(pass == 1 | rush == 1, (qtr == 4 & quarter_seconds_remaining < 300)| game_seconds_remaining < 300, pointdiff <= 8, qb_kneel == 0, qb_spike ==0) %>% 
  mutate(poscoach = ifelse(home_team == posteam, home_coach,away_coach),
         defcoach = ifelse(home_team == defteam, home_coach,away_coach))
off_coach <- coach_clutch %>%   
  group_by(poscoach) %>% 
  summarize(off_plays = n(),off_epa = mean(epa,na.rm = T))
comb_clutch <- coach_clutch %>%   
  group_by(defcoach) %>% 
  summarize(def_plays = n(),def_epa = mean(epa,na.rm = T)) %>% 
  left_join(off_coach, by = c("defcoach" = "poscoach")) %>% 
  mutate(total_plays = off_plays+def_plays)

library(ggrepel)

comb_clutch %>% 
  filter(total_plays>=200) %>% 
  ggplot(aes(x = def_epa, y = off_epa)) +
  # geom_text(aes(label = defcoach))+
  geom_text_repel(
    aes(label = defcoach),
    box.padding = 0.01,  # adjust this value for padding around the labels
    point.padding = 0.01,  # adjust this value for padding around the points
    segment.color = "grey",
    segment.size = 0.1
  )+
  scale_x_reverse()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(x = "Clutch Defensive EPA/Play", y = "Clutch Offensive EPA/Play", title = "Who are the Coaches who have Performed the Best in the Clutch Since 1999?",
       subtitle = "Clutch Defined as 1 Score Game with Less than 5 Minutes Left, Minimum 200 Clutch Plays", caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")
ggsave("CoachClutch.png", width = 14, height =10, dpi = "retina")

# Eberflus Evaluation ----
eberflus <- nfl99 %>% 
  mutate(pointdiff = abs(total_home_score - total_away_score)) %>% 
  filter(home_coach == "Matt Eberflus" | away_coach == "Matt Eberflus",qb_kneel == 0, qb_spike ==0,pass == 1 | rush == 1) %>% 
  mutate(isclutch = ifelse((qtr == 4 & quarter_seconds_remaining < 300 & pointdiff <= 8) | (game_seconds_remaining < 300 & pointdiff <= 8), "Clutch","Non-Clutch")) %>% 
  mutate(poscoach = ifelse(home_team == posteam, home_coach,away_coach),
         defcoach = ifelse(home_team == defteam, home_coach,away_coach)) %>% 
  filter(poscoach == "Matt Eberflus") %>% 
  mutate(gamestate = ifelse((total_home_score >total_away_score & home_team == posteam) | (total_home_score < total_away_score & away_team == posteam),"Winning","Tied/Losing")) %>% 
  mutate(short_throw = ifelse(air_yards<=10,1,0),
         medium_throw = ifelse(air_yards>10&air_yards<=20,1,0),
         long_throw = ifelse(air_yards>20,1,0)) %>% 
  group_by(isclutch) %>% 
  summarize(
    plays = n(),
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
  )%>% 
  mutate_if(is.numeric,~round(.,2))
eberflus_tab <- eberflus %>% 
  select(isclutch,off_epa,off_pass_epa,pass_rate, aDot, off_rush_epa,off_early_down_epa,off_third_down_dist,off_3rd_down_1st_pct) %>% 
  gt() %>% 
  cols_align(align = "center") %>% 
  cols_label(isclutch = "Game State", off_epa = "EPA/Play", off_pass_epa = "EPA/Dropback",
             off_rush_epa = "EPA/Rush", off_early_down_epa = "Early Down EPA", off_third_down_dist = "Avg 3rd Down Distance",
             pass_rate = "Pass Rate", off_3rd_down_1st_pct = "3rd Down Conversion Rate") %>% 
  gtExtras::gt_theme_538() %>% 
  tab_header(title = md("Where does it go wrong for Eberflus's offense?"), subtitle = "Eberflus's offense sees regression in just about every category in the clutch")
gtsave(eberflus_tab,"eberflus.png")

clutch_vs_nonclutch <- nfl99 %>% 
  mutate(pointdiff = abs(total_home_score - total_away_score)) %>% 
  filter(qb_kneel == 0, qb_spike ==0,pass == 1 | rush == 1) %>% 
  mutate(isclutch = ifelse((qtr == 4 & quarter_seconds_remaining < 300 & pointdiff <= 8) | (game_seconds_remaining < 300 & pointdiff <= 8), "Clutch","Non-Clutch")) %>% 
  mutate(poscoach = ifelse(home_team == posteam, home_coach,away_coach))

clutch_off <- clutch_vs_nonclutch %>% 
  group_by(posteam,season,isclutch) %>% 
  summarize(epa_play = mean(epa,na.rm = T)) %>% 
  filter(!is.na(posteam)) %>% 
  pivot_wider(names_from = isclutch, values_from = epa_play) %>% 
  filter(!is.na(Clutch),!is.na(`Non-Clutch`))
clutch_off %>% 
  ggplot(aes(x = `Non-Clutch`, y = Clutch))+
  geom_point()+
  labs(x = "Non-Clutch Offensive EPA/Play",y = "Clutch Offensive EPA/Play",
       title = "Is Offensive Performance Outside of Clutch Situations Correlated to Clutch Offensive Performance?",
       subtitle = "Having a good offense in non-clutch situations does not translate to an efficient clutch offense",
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")+
  geom_smooth(method = "lm", se = F)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())
cor(clutch_off$Clutch,clutch_off$`Non-Clutch`)
ggsave("OverallClutch.png", width = 14, height =10, dpi = "retina")

clutch_pass <- clutch_vs_nonclutch %>% 
  filter(pass == 1) %>%
  group_by(posteam,season,isclutch) %>%
  summarize(epa_play = mean(epa,na.rm = T)) %>%
  filter(!is.na(posteam)) %>%
  pivot_wider(names_from = isclutch, values_from = epa_play) %>%
  filter(!is.na(Clutch),!is.na(`Non-Clutch`))
clutch_pass %>% 
  ggplot(aes(x = `Non-Clutch`, y = Clutch))+
  geom_point()+
  labs(x = "Non-Clutch Offensive EPA/Dropback",y = "Clutch Offensive EPA/Dropback",
       title = "What is the Relationship Between Clutch Passing Efficiency and Non-Clutch Passing Efficiency?",
       subtitle = "There is a weak relationship between clutch & non-clutch passing efficiency",
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")+
  geom_smooth(method = "lm", se = F)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())
cor(clutch_pass$Clutch,clutch_pass$`Non-Clutch`)
ggsave("ClutchPassing.png", width = 14, height =10, dpi = "retina")

clutch_rush <- clutch_vs_nonclutch %>% 
  filter(rush == 1) %>% 
  group_by(posteam,season,isclutch) %>% 
  summarize(epa_play = mean(epa,na.rm = T)) %>% 
  filter(!is.na(posteam)) %>% 
  pivot_wider(names_from = isclutch, values_from = epa_play) %>% 
  filter(!is.na(Clutch),!is.na(`Non-Clutch`))
clutch_rush %>% 
  ggplot(aes(x = `Non-Clutch`, y = Clutch))+
  geom_point()+
  labs(x = "Non-Clutch Offensive EPA/Rush",y = "Clutch Offensive EPA/Rush",
       title = "Does an efficient run offense translate to an efficient clutch run offense?",
       subtitle = "There is a weak correlation between clutch rush efficiency and non-clutch rush efficiency",
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")+
  geom_smooth(method = "lm", se = F)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())
cor(clutch_rush$Clutch,clutch_rush$`Non-Clutch`)
ggsave("RushClutch.png", width = 14, height =10, dpi = "retina")


xpass_performance <- clutch_vs_nonclutch %>% 
  filter(season >= 2006) %>% 
  mutate(game_state = ifelse(score_differential_post >= 0, "winning/tied","losing")) %>% 
  group_by(posteam,season) %>% 
  summarize(epa_play = mean(epa[xpass>=0.8 & isclutch == "Non-Clutch"],na.rm = T), clutch_play = mean(epa[isclutch == "Clutch"],na.rm = T)) %>% 
  # filter(!is.na(posteam)) %>% 
  # pivot_wider(names_from = isclutch, values_from = epa_play) %>% 
  filter(!is.na(epa_play),!is.na(clutch_play))
xpass_performance %>% 
  ggplot(aes(x = epa_play, y = clutch_play))+
  geom_point()+
  labs(x = "EPA/Play in Non-Clutch Expected Pass Situations",y = "Clutch Offensive EPA/Play",
       title = "Does Offensive Performance Outside in Expected Pass Situtions Correlate to Clutch Offensive Performance?",
       subtitle = "Despite expected pass situations mirroring potential clutch scenarios, there is a weak relationship between performance in expected pass situations and clutch offensive efficiency",
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")+
  geom_smooth(method = "lm", se = F)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())
ggsave("Xpass.png", width = 14, height =10, dpi = "retina")
cor(xpass_performance$epa_play,xpass_performance$clutch_play)
#EPA on Expected Pass Situations

shanahan <- clutch_vs_nonclutch %>% 
  filter(poscoach == "Kyle Shanahan") %>% 
  group_by(season) %>% 
  summarize(epa_play = round(mean(epa[isclutch == "Non-Clutch"],na.rm = T),2), clutch_play = round(mean(epa[isclutch == "Clutch"],na.rm = T),2)) %>% 
  gt() %>% 
  cols_align(align = "center") %>% 
  cols_label(season = "Season", epa_play = "Non-Clutch EPA/Play", clutch_play = "Clutch EPA/Play") %>% 
  gtExtras::gt_theme_538() %>% 
  tab_header(title = md("How Does Kyle Shanahan Fair in the Clutch?"))
gtsave(shanahan,"shanahan.png")

clutch_time <- nfl99 %>% 
  mutate(pointdiff = abs(total_home_score - total_away_score)) %>% 
  filter(pass == 1 | rush == 1, (qtr == 4 & quarter_seconds_remaining < 300)| game_seconds_remaining < 300, pointdiff <= 8)
clutch_time %>% 
  filter(qb_kneel == 0, qb_spike == 0) %>% 
  group_by(id) %>% 
  summarize(epa_play = mean(epa, na.rm = T), plays = n(), player_name = max(name), passing_plays = sum(pass), team = names(sort(-table(posteam)))[1]) %>% 
  filter(plays > 200,passing_plays > 0) %>% 
  left_join(teams_colors_logos, by = c("team" = "team_abbr")) %>% 
  ggplot(aes(x = reorder(player_name,-epa_play),y = epa_play, fill = team_color, color = team_color2))+
  geom_bar(stat = "identity")+
  scale_fill_identity()+
  scale_color_identity()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 35, hjust = 1))+
  labs(x = "Player", y = "Clutch EPA/Play",title = "Who are the Best QBs Since 1999 in the Clutch Minimum 200 Clutch Plays?", subtitle = "Mahomes is the King of the Clutch", 
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")
ggsave("BarClutch.png", width = 14, height =10, dpi = "retina")

clutch_time %>% 
  filter(qb_kneel == 0, qb_spike == 0) %>% 
  group_by(id) %>% 
  summarize(epa_pass = mean(epa[pass_attempt == 1],na.rm = T), epa_rush  = mean(epa[rush == 1 | qb_scramble == 1],na.rm = T), 
            rush_plays = sum(rush == 1 | qb_scramble == 1), plays = n(), player_name = max(name), passing_plays = sum(pass), team = names(sort(-table(posteam)))[1]) %>% 
  filter(plays > 200,passing_plays > 0) %>% 
  left_join(teams_colors_logos, by = c("team" = "team_abbr")) %>% 
  ggplot(aes(x = epa_rush,y = epa_pass, color = team_color))+
  # geom_text(aes(label = player_name))+
  geom_text_repel(
    aes(label = player_name),
    box.padding = 0.01,  # adjust this value for padding around the labels
    point.padding = 0.01,  # adjust this value for padding around the points
    # segment.color = "grey",
    segment.size = 0.1
  )+
  scale_color_identity()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 35, hjust = 1))+
  labs(x = "Clutch EPA/Rush Includes Scrambles", y = "Clutch EPA/Pass",title = "Who are the Best QBs Since 1999 in the Clutch at Rushing and Passing in the Clutch?", subtitle = "Mahomes Dominates on the Ground and Through the Air", 
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")
ggsave("ClutchPlot.png", width = 14, height =10, dpi = "retina")


clutch_time %>% 
  filter(qb_kneel == 0, qb_spike == 0) %>% 
  group_by(id) %>% 
  summarize(epa_play = mean(epa, na.rm = T), plays = n(), player_name = max(name), passing_plays = sum(pass), team = names(sort(-table(posteam)))[1],
            last_season = max(season)) %>% 
  filter(plays <= 200,last_season == 2023, passing_plays > 20) %>% 
  left_join(teams_colors_logos, by = c("team" = "team_abbr")) %>% 
  ggplot(aes(x = reorder(player_name,-epa_play),y = epa_play, fill = team_color, color = team_color2))+
  geom_bar(stat = "identity")+
  scale_fill_identity()+
  scale_color_identity()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())+
  theme(axis.text.x = element_text(angle = 35, hjust = 1))+
  labs(x = "Player", y = "Clutch EPA/Play",title = "Who are the Best QBs Who Played in 2023, but have less than 200 Career Clutch Plays?", subtitle = "What an incredible clutch performance season for Stroud in his rookie year", 
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")
ggsave("BarClutch.png", width = 14, height =10, dpi = "retina")

clutch_time %>% 
  filter(qb_kneel == 0, qb_spike == 0) %>% 
  group_by(id) %>% 
  summarize(epa_pass = mean(epa[pass_attempt == 1],na.rm = T), epa_rush  = mean(epa[rush == 1 | qb_scramble == 1],na.rm = T), 
            rush_plays = sum(rush == 1 | qb_scramble == 1), plays = n(), player_name = max(name), passing_plays = sum(pass), team = names(sort(-table(posteam)))[1],
            last_season = max(season)) %>% 
  filter(plays <= 200,last_season == 2023, passing_plays > 20) %>% 
  left_join(teams_colors_logos, by = c("team" = "team_abbr")) %>% 
  ggplot(aes(x = epa_rush,y = epa_pass, color = team_color))+
  # geom_text(aes(label = player_name))+
  geom_text_repel(
    aes(label = player_name),
    box.padding = 0.01,  # adjust this value for padding around the labels
    point.padding = 0.01,  # adjust this value for padding around the points
    # segment.color = "grey",
    segment.size = 0.1
  )+
  scale_color_identity()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())+
  labs(x = "Clutch EPA/Rush Includes Scrambles", y = "Clutch EPA/Pass",title = "Stroud Dominates Clutch Situations Through the Air", 
       subtitle = "QB Must Have Played in 2023 and Had Less than 200 Clutch Career Plays", 
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")
ggsave("ClutchPlot.png", width = 14, height =10, dpi = "retina")

qb_cor <- nfl99 %>% 
  filter(pass ==1 | rush == 1, qb_kneel == 0, qb_spike == 0) %>%
  mutate(pointdiff = abs(total_home_score - total_away_score)) %>% 
  mutate(isclutch = ifelse((qtr == 4 & quarter_seconds_remaining < 300 & pointdiff <= 8) | (game_seconds_remaining < 300 & pointdiff <= 8), "Clutch","Non-Clutch")) %>% 
  group_by(id) %>% 
  summarize(epa_clutch = mean(epa[isclutch == "Clutch"], na.rm = T), epa_n_clutch = mean(epa[isclutch == "Non-Clutch"],na.rm = T),plays = sum(isclutch == "Clutch"), 
            passing_plays = sum(pass), epa_clutch_pa = mean(epa[isclutch == "Clutch" & pass_attempt == 1], na.rm = T), 
            epa_clutch_db = mean(epa[isclutch == "Clutch" & pass == 1], na.rm = T), epa_nclutch_pa = mean(epa[isclutch == "Non-Clutch" & pass_attempt == 1], na.rm = T), 
            epa_nclutch_db = mean(epa[isclutch == "Non-Clutch" & pass == 1], na.rm = T)) %>% 
  filter(plays > 200,passing_plays > 0)
qb_cor%>% 
  ggplot(aes(x = epa_n_clutch,y = epa_clutch))+
  geom_point()+
  geom_smooth(method = "lm",se=F)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())+
  labs(x = " Non-Clutch EPA/ Play", y = "Clutch EPA/Play",title = "Does Performance in Non-Clutch Situations Translate to Clutch Performance?", 
       subtitle = "QBs who perform best in non-clutch situations tend to perform the best in the clutch", 
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")
ggsave("Clutchcor.png", width = 14, height =10, dpi = "retina")

test <- nfl99 %>%
  filter(pass == 1,qb_kneel == 0, qb_spike == 0) %>% 
  mutate(throw_type = ifelse(air_yards<=10,"short",ifelse(air_yards>10&air_yards<=20,"medium",ifelse(air_yards>20,"long","Behind LOS")))) %>% 
  group_by(throw_type,pass_location, passer_player_id) %>% 
  summarize(name = max(passer_player_name),attempts = n(), epa_pass = mean(epa,na.rm = T)) %>% 
  filter(attempts > 75) %>% 
  filter(!is.na(passer_player_id),!is.na(pass_location))

key_num_freq <- nfl99 %>% 
  group_by(game_id) %>% 
  summarize(key_num = abs(max(result)), year = max(season)) %>% 
  group_by(year, key_num) %>%
  summarise(count = n()) %>%
  mutate(relative_freq = count / sum(count)) %>%
  ungroup() %>% 
  select(-count)

library(RColorBrewer)

key_num_freq %>% 
  #filter(key_num == 3 | key_num ==  7) %>% 
  filter(key_num <= 11 | key_num == 14) %>%
  filter(key_num>1) %>% 
  ggplot(aes(x = year, y = relative_freq, color = as.factor(key_num))) +
  #geom_point()+
  geom_smooth(se = FALSE, lwd = 1.5) +
  #scale_color_manual(values = c("3" = "red", "7" = "grey"))+
  scale_color_manual(values = c("3" = "red", "7" = "grey", "11" = "black", "2" = "darkblue", "4" = "lightblue", "5" = "darkgreen", "6" = "purple", "8" = "brown" , "9" = "lightgreen", "10" = "pink", "14"= "orange"))+
  labs(title = "Have the Traditional Key Numbers Become Less Important in Handicapping NFL Sides?",
       subtitle = "3 & 7 are still the most important key numbers, but 4 & 6 have seen a rise in their importance",
       x = "Season",
       y = "Frequency",
       color = "Key Number",
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())
ggsave("Clutchcor.png", width = 14, height =10, dpi = "retina")

rates <- nfl99 %>% 
  mutate(attempted4th = ifelse(fourth_down_converted == 1 | fourth_down_failed == 1,1,0),
         pointaftertry = ifelse(extra_point_attempt == 1 | two_point_attempt == 1, 1, 0),
         extra_point_result_bin = ifelse(extra_point_result == "good",1,0)) %>% 
  #filter(pointaftertry == 1) %>% 
  group_by(season) %>% 
  summarize(go_for_2_rate = mean(two_point_attempt[pointaftertry == 1], na.rm = T),
            go_for_4th_down_rate = mean(attempted4th[down == 4], na.rm = T),
            pat_success_rate = mean(extra_point_result_bin[extra_point_attempt==1],na.rm = T)) %>% 
  pivot_longer(cols = c("go_for_2_rate","go_for_4th_down_rate","pat_success_rate"), values_to = "relative_freq",names_to = "key_num") %>% 
  rename("year" = season)

test <- rbind(key_num_freq,rates) %>% 
  #filter(key_num %in% c(1,2,3,5,7,4,6,8,"go_for_2_rate", "go_for_4th_down_rate")) %>% 
  filter(key_num %in% c(3,7,4,6,8,"go_for_2_rate", "go_for_4th_down_rate")) %>% 
  ggplot(aes(x = year, y = relative_freq, color = as.factor(key_num), linetype = as.factor(key_num))) +
  geom_line(se = FALSE, lwd = 1.5) +
  scale_color_manual(values = c("3" = "red", "7" = "grey", "4" = "lightblue", "6" = "purple","8" = "brown", "go_for_4th_down_rate" = "black", "go_for_2_rate" = "pink"))+
  scale_linetype_manual(values = c("3" = "solid", "7" = "solid", "4" = "solid", "6" = "solid","8" = "solid", "go_for_4th_down_rate" = "dotdash", "go_for_2_rate" = "dotdash"), guide = FALSE)+
  labs(title = "Is There a Relationship Between Changes in NFL In-Game Decision Making and Key Numbers?",
       subtitle = "The resurgence in going for 2 coincided with the decline of importance of 7 as a key number and a rise in importance of 6 & 8",
       x = "Season",
       y = "Frequency",
       color = "Key Number",
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())+
  guides(
    color = guide_legend(override.aes = list(linetype = c("solid", "solid", "solid", "solid","solid", "dotdash", "dashed")))
  )
ggsave("Clutchcor.png", width = 14, height =10, dpi = "retina")

# Wong Teaser Research ----


teaser_data <- nfl99 %>%
  group_by(game_id) %>% 
  summarize(year = max(season), homescore = max(home_score), away_score = max(away_score),spread = max(spread_line), result = max(result), bet_total = max(total_line), result_total = max(total),
            go_for_2_rate = sum(two_point_attempt, na.rm = T),date = min(game_date)) %>% 
  mutate(favorite = ifelse(spread >= 0, "home", "away"), favorite_cover = ifelse((favorite == "home" & result > spread-6) |(favorite == "away" & result < spread+6),1,0 )) %>% 
  mutate(underdog = ifelse(spread >= 0, "away", "home"), underdog_cover = ifelse((underdog == "home" & result > spread-6) |(underdog == "away" & result < spread+6),1,0 ))


teaser_legs <- teaser_data %>% 
  select(-favorite,-underdog) %>% 
  pivot_longer(cols = c("favorite_cover", "underdog_cover"), values_to = "cover", names_to = "spread_type") %>% 
  mutate(spread = ifelse(spread_type == "underdog_cover" & spread < 0,spread*-1, ifelse(spread_type == "favorite_cover" & spread > 0, spread*-1,spread))) %>% 
  group_by(spread) %>% 
  summarize(cover_rate = round(mean(cover),3), count = n())

teaser_tab <- teaser_legs %>% 
  filter(count >= 10) %>%
  arrange(desc(cover_rate)) %>% 
  gt() %>% 
  cols_align(align = "center") %>%
  cols_label(spread = "Game Spread",cover_rate = "Cover Rate in a 6 Point Teaser Leg", count = "Number of Occurences") %>% 
  tab_style(
    style = list(
      cell_fill(color = "green")
    ),
    locations = cells_body(
      columns = c(spread,cover_rate,count),
      rows = (cover_rate) > 0.7386))%>%
  gtExtras::gt_theme_538() %>%
  tab_header(title = md("Which Spreads are the Best to Put into a 6 Point Teaser Since 1999?"), subtitle = md("+1.5, -7.5, +2.5, -8.5 are the most frequent occurring profitable winners in a 6-point. (Cells in green surpass 73.86% cover rate based on averages odds for a 6-point teaser across sportsbooks)"))
gtsave(teaser_tab,"Teaser.png")
# ggplot(aes(x = spread, y = cover_rate))+
# geom_bar(stat = "identity")
#mutate(teaser_line = spread_value+6, covered = ifelse(result > 0 & spread > 0,1 ifelse()))


teaser_years <- teaser_data %>% 
  pivot_longer(cols = c("favorite_cover", "underdog_cover"), values_to = c("cover"), names_to = c("spread_type")) %>% 
  mutate(push = ifelse(((favorite == "home" & result == spread-6 & spread_type == "favorite_cover" ) |(favorite == "away" & spread_type == "favorite_cover" & result ==  spread+6)) | ((underdog == "home" & spread_type == "underdog_cover"  & result == spread-6) |(underdog == "away" & spread_type == "underdog_cover"  & result == spread+6)) ,1,0 )) %>% 
  mutate(spread = ifelse(spread_type == "underdog_cover" & spread < 0,spread*-1, ifelse(spread_type == "favorite_cover" & spread > 0, spread*-1,spread))) %>% 
  mutate(pre_2015 = ifelse(year>= 2015, "2015", "2014"), teaser_leg = ifelse(spread_type == "underdog_cover", spread+6,spread-6)) %>% 
  group_by(spread,pre_2015) %>% 
  summarize(cover_rate = round(mean(cover),3), count = n(), push_rate = round(mean(push),3), pt2 = mean(go_for_2_rate,na.rm = T)) %>% 
  pivot_wider( values_from = c(cover_rate,count,push_rate, pt2), names_from = pre_2015) %>% 
  filter(count_2015 > 30, count_2014 > 30) %>% 
  mutate(diff = cover_rate_2015 - cover_rate_2014, pt2diff = pt2_2015 - pt2_2014)

rates %>% 
  filter(key_num %in% c("go_for_2_rate")) %>% 
  ggplot(aes(x = year, y = relative_freq)) +
  geom_line(lwd = 1.5) +
  geom_point()+
  labs(title = "How Has the 2 Point Attempt Rate Changed Since 2015?",
       subtitle = "2 Point Attempt Rate Has Doubled Since 2014 With a Major Shift Occurring in 2015 (The Season When the PAT was Moved Back)",
       x = "Season",
       y = "2 Point Attempt Rate",
       color = "Key Number",
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())
ggsave("2PTsuccess.png", width = 14, height =10, dpi = "retina")

rates %>% 
  filter(key_num %in% c("pat_success_rate")) %>% 
  ggplot(aes(x = year, y = relative_freq)) +
  geom_line(lwd = 1.5) +
  geom_point()+
  labs(title = "How Has the PAT Rate Changed Since 2015?",
       subtitle = "PAT success rate has seen a decrease since the PAT was moved back in 2015",
       x = "Season",
       y = "PAT Succcess Rate",
       color = "Key Number",
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())
ggsave("PATsuccess.png", width = 14, height =10, dpi = "retina")

negative_shift <- teaser_years %>% 
  ungroup %>% 
  arrange(diff) %>% 
  slice_head(n = 10) %>% 
  select(spread,cover_rate_2014,cover_rate_2015,diff) %>% 
  gt() %>% 
  cols_align(align = "center") %>% 
  cols_label(spread = "Game Spread", cover_rate_2015 = "Cover Rate Starting in 2015", cover_rate_2014 = "Cover Rate Before 2015",
             diff = "Cover Rate Difference") %>% 
  tab_style(
    style = list(
      cell_fill(color = "green")
    ),
    locations = cells_body(
      columns = c(cover_rate_2014),
      rows = (cover_rate_2014) > 0.7386))%>%
  tab_style(
    style = list(
      cell_fill(color = "green")
    ),
    locations = cells_body(
      columns = c(cover_rate_2015),
      rows = (cover_rate_2015) > 0.7386))%>%
  gtExtras::gt_theme_538() %>% 
  tab_header(title = md("Which Spreads Have Seen the Biggest Negative Shift in Cover Rate in 6 Point Teasers Since 2015?"), subtitle = md("Since 2015, -4 & -2 are no longer profitable plays with changes in their cover rates by more than 9%")) %>% 
  tab_footnote(footnote = md("Cells in green surpass 73.86% cover rate based on averages odds for a 6-point teaser across sportsbooks. Spreads must have occurred more than 30 times before and after 2015. Data dates back to 1999.")) %>% 
  tab_source_note("Callan Capitolo | @CapAnalytics7 | nflfastR")
gtsave(negative_shift,"negshift.png")

pos_shift <- teaser_years %>% 
  ungroup %>% 
  arrange(-diff) %>% 
  slice_head(n = 10) %>% 
  select(spread,cover_rate_2014,cover_rate_2015,diff) %>% 
  gt() %>% 
  cols_align(align = "center") %>% 
  cols_label(spread = "Game Spread", cover_rate_2015 = "Cover Rate Starting in 2015", cover_rate_2014 = "Cover Rate Before 2015",
             diff = "Cover Rate Difference") %>% 
  tab_style(
    style = list(
      cell_fill(color = "green")
    ),
    locations = cells_body(
      columns = c(cover_rate_2014),
      rows = (cover_rate_2014) > 0.7386))%>%
  tab_style(
    style = list(
      cell_fill(color = "green")
    ),
    locations = cells_body(
      columns = c(cover_rate_2015),
      rows = (cover_rate_2015) > 0.7386))%>%
  gtExtras::gt_theme_538() %>% 
  tab_header(title = md("Which Spreads Have Seen the Biggest Positive Shift in Cover Rate in 6 Point Teasers Since 2015?"), subtitle = md("+5.5, +4, and +10 have seen their cover rates improve by more than 10% in the past 9 seasons")) %>% 
  tab_footnote(footnote = md("Cells in green surpass 73.86% cover rate based on averages odds for a 6-point teaser across sportsbooks. Spreads must have occurred more than 30 times before and after 2015. Data dates back to 1999.")) %>% 
  tab_source_note("Callan Capitolo | @CapAnalytics7 | nflfastR")
gtsave(pos_shift,"posshift.png")
#investiate push rates for numbers like 2 & - 8
#Fanduel: individual teaser legs are -279 for a 4 leg teaser payout of +240, 2 leg teaser is -311 with payout of -134, 3 leg teaser -295 payout of +140, -294 5 leg payout of 333, -285 for 6 leg payout of 500

cover_trends <- teaser_years %>% 
  ungroup %>% 
  filter(as.numeric(cover_rate_2015) > 0.7386) %>% 
  arrange(-cover_rate_2015) %>% 
  select(spread,cover_rate_2015, count_2015)
cover_trends %>% 
  gt() %>% 
  cols_align(align = "center") %>% 
  cols_label(spread = "Game Spread", cover_rate_2015 = "Cover Rate Starting in 2015", count_2015 = "# of Occurrences") %>% 
  gtExtras::gt_theme_538() %>% 
  tab_header(title = md("What Spreads Have Been Most Profitable in the Past 9 Seasons?"), subtitle = md("+5.5, +4, and +10 have seen their cover rates improve by more than 10% in the past 9 seasons")) %>% 
  tab_footnote(footnote = md("Cells in green surpass 73.86% cover rate based on averages odds for a 6-point teaser across sportsbooks. Spreads must have occurred more than 30 after the 2014 season.")) %>% 
  tab_source_note("Callan Capitolo | @CapAnalytics7 | nflfastR") %>% 
  gtExtras::gt_hulk_col_numeric(columns = count_2015)

library(zoo)
teaser_data %>% 
  select(-favorite,-underdog) %>% 
  pivot_longer(cols = c("favorite_cover", "underdog_cover"), values_to = "cover", names_to = "spread_type") %>% 
  mutate(spread = ifelse(spread_type == "underdog_cover" & spread < 0,spread*-1, ifelse(spread_type == "favorite_cover" & spread > 0, spread*-1,spread))) %>% 
  right_join(cover_trends %>% select(spread), by = "spread") %>% 
  arrange(spread,game_id) %>% 
  group_by(spread) %>% 
  mutate(game_number = row_number()) %>% 
  mutate(roll_cover = rollapply(cover, width = 70, align = "right", FUN = mean,partial = TRUE))%>% 
  ggplot(aes(x = game_number, y = roll_cover, color = as.factor(spread))) +
  geom_smooth(se = F, span = 0.25, lwd = 1.5)+
  geom_hline(yintercept = 0.7386, linetype = "dashed")+
  #geom_point()+
  labs(y = "Rolling 70 Game Cover Rate", x = "Game Number",legend = "spread", color = "Spread",
       title = "How Have Cover Rates in 6 Point Teasers Since 1999?", subtitle = "+4 stands out as the spread trending in the right direction in recent seasons",
       caption = "*Spreads must have a cover rate above 0.7386 since 2015                 Callan Capitolo | @CapAnalytics7 | nflfastR")+
  scale_color_manual(values = c("-8.5" = "red", "-8" = "grey", "-7.5" = "black", "-6.5" = "purple", "1.5" = "lightblue", "2" = "darkgreen", "2.5" = "blue", "4" = "brown" , "4.5" = "lightgreen",
                                "5.5" = "orange", "10" = "pink"))+
  xlim(20,500)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())
ggsave("covertrend.png", width = 14, height =10, dpi = "retina")


teaser_data %>% 
  pivot_longer(cols = c("favorite_cover", "underdog_cover"), values_to = "cover", names_to = "spread_type") %>% 
  mutate(prop = abs(round(spread/bet_total,2))) %>% 
  group_by(bet_total) %>% 
  summarize(cover = mean(cover), count = n()/2) %>% 
  filter(count>20) %>% 
  ggplot(aes(x = bet_total, y = cover)) +
  geom_hline(yintercept = 0.7386, linetype = "dashed")+
  geom_point()+
  geom_line()+
  guides(size = "none")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())+
  annotate("text", x = 36, y = 0.7386, label = "Break Even Cover Rate", vjust = -1, hjust = 1.1)+
  labs(x = "Total", y = "Spread Cover Rate in 6 Point Teaser", title = "Do Side Teaser Legs Perform Better With Low Totals?",
       subtitle = "Surpisingly there is no clear difference in spread cover rates in low totals vs high totals",
       caption = "Callan Capitolo | @CapAnalytics7 | nflfastR" )
ggsave("covertotal.png", width = 14, height =10, dpi = "retina")

teaser_data %>% 
  pivot_longer(cols = c("favorite_cover", "underdog_cover"), values_to = "cover", names_to = "spread_type") %>% 
  mutate(prop = abs(round(spread/bet_total,2))) %>% 
  group_by(bet_total,spread_type) %>% 
  summarize(cover = mean(cover), count = n()) %>% 
  filter(count>20) %>% 
  ggplot(aes(x = bet_total, y = cover, color = spread_type)) +
  geom_hline(yintercept = 0.7386, linetype = "dashed")+
  geom_point(aes(size = count))+
  geom_line()+
  guides(size = "none")+
  scale_color_manual(values = c("favorite_cover" = "grey", "underdog_cover" = "red"),
                     labels = c("favorite_cover" = "Favorite", "underdog_cover" = "Underdog"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())+
  annotate("text", x = 36, y = 0.7386, label = "Break Even Cover Rate", vjust = -1, hjust = 1.1)+
  labs(x = "Bet Total", y = "Cover Rate", title = "How do Cover Rates for Underdogs and Favorites in 6 Point Teasers Change with Game Totals Since 1999?",
       subtitle = "Underdogs perform much better with very low totals than favorites",color = "Side",caption = "@CapAnalytics7 | nflfastR" )
ggsave("underdogtoal.png", width = 14, height =10, dpi = "retina")

teaser_data %>% 
  pivot_longer(cols = c("favorite_cover", "underdog_cover"), values_to = "cover", names_to = "spread_type") %>% 
  mutate(prop = abs(round(spread/bet_total,2))) %>% 
  group_by(bet_total,spread_type,year) %>% 
  summarize(cover = mean(cover), count = n()) %>% 
  filter(count > 5) %>% 
  ggplot(aes(x = year, y = bet_total, fill = cover)) +
  geom_tile()+
  labs(x = "Bet Total", y = "Cover Rate", title = "How do Cover Rates for Underdogs and Favorites in 6 Point Teasers Change with Game Totals Since 1999?",
       subtitle = "Underdogs perform much better with very low totals than favorites",color = "Side",caption = "@CapAnalytics7 | nflfastR" )


side_tot <- teaser_data %>% 
  pivot_longer(cols = c("favorite_cover", "underdog_cover"), values_to = "cover", names_to = "spread_type") %>% 
  mutate(spread = ifelse(spread_type == "underdog_cover" & spread < 0,spread*-1, ifelse(spread_type == "favorite_cover" & spread > 0, spread*-1,spread))) %>% 
  group_by(bet_total,spread) %>%
  summarize(cover = round(mean(cover),4), count = n()) %>% 
  filter(count>20) %>% 
  ungroup %>% 
  arrange(-cover) %>% 
  slice_head(n = 10) %>% 
  gt() %>% 
  cols_align(align = "center") %>% 
  cols_label(spread = "Game Spread", count = "# of Occurrences", bet_total = "Total", cover = "6 Point Teaser Side Cover Rate") %>% 
  gtExtras::gt_theme_538() %>% 
  tab_header(title = md("Which Spread/Total Combinations Have the Highest Side Cover Rates in a 6 Point Teaser Since 1999?"), subtitle = md("46 and +2.5 have an astronomical side cover rate as a combination")) %>% 
  tab_footnote(footnote = md("Top 10 spread/total combinations based on side cover rate with a miminum of 20 occurrences")) %>% 
  tab_source_note("Callan Capitolo | @CapAnalytics7 | nflfastR") %>% 
  gtExtras::gt_theme_538()
gtsave(side_tot,"side_tot.png")  




# Teaser Modeling ----
library(caret)
library(xgboost)
library(lightgbm)


model_data <- teaser_data %>% 
  filter(year>=2015) %>% 
  pivot_longer(cols = c("favorite_cover", "underdog_cover"), values_to = "cover", names_to = "spread_type") %>% 
  mutate(spread = ifelse(spread_type == "underdog_cover" & spread < 0,spread*-1, ifelse(spread_type == "favorite_cover" & spread > 0, spread*-1,spread))) %>% 
  select(spread,bet_total,cover, year)
X<- model_data[,c('spread','bet_total','year')]
y <- ifelse(model_data$cover == 1, "yes", "no")

# Define the current year for weighting
set.seed(42)
train_indices <- createDataPartition(y, p = 0.8, list = FALSE)
X_train <- X[train_indices, ]
y_train <- y[train_indices]
X_test <- X[-train_indices, ]
y_test <- y[-train_indices]

train_control <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  verboseIter = TRUE
)

# Define parameter grid for XGBoost
xgb_grid <- expand.grid(
  nrounds = c(100, 200),
  max_depth = c(3, 6, 9),
  eta = c(0.01, 0.1, 0.3),
  gamma = c(0, 1),
  colsample_bytree = c(0.5, 0.7),
  min_child_weight = c(1, 5),
  subsample = c(0.5, 0.7)
)

# Train XGBoost model
set.seed(42)
xgb_model <- train(
  x = X_train,
  y = as.factor(y_train),
  method = "xgbTree",
  trControl = train_control,
  tuneGrid = xgb_grid,
  metric = "ROC"
)


predictions <- predict(xgb_model, newdata = X_test, type = "prob")

test <- cbind(X_test,predictions$yes)



# Negative vs Explosive Plays ----

nfl99 %>% 
  filter(pass == 1 | rush == 1, season == 2023) %>% 
  group_by(defteam) %>% 
  summarize(exp_rate = mean(explosive,na.rm = T), neg_rate = mean(negative,na.rm = T)) %>% 
  left_join(teams_colors_logos, by = c("defteam" = "team_abbr")) %>% 
  ggplot(aes(x = neg_rate, y = exp_rate))+
  geom_image(aes(image = team_logo_espn), size = 0.05, asp = 16/9)+
  labs(x = "Negative Play Rate", y = "Explosive Play Rate*", title = "Which Defenses Created Negative Plays and Limited Explosive Plays in 2023?",
       subtitle = "The Browns Created a Negative Play Once Every 8 Plays",
       caption = "*Passes that gained greater than 20 yards or runs that gained greater than 12 yards                         Callan Capitolo | @CapAnalytics7 | nflfastR")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())+
  scale_y_reverse()
ggsave("DefExpvsNeg.png", width = 14, height =10, dpi = "retina")

nfl99 %>% 
  filter(pass == 1 | rush == 1, season == 2023) %>% 
  group_by(posteam) %>% 
  summarize(exp_rate = mean(explosive,na.rm = T), neg_rate = mean(negative,na.rm = T)) %>% 
  left_join(teams_colors_logos, by = c("posteam" = "team_abbr")) %>% 
  ggplot(aes(x = neg_rate, y = exp_rate))+
  geom_image(aes(image = team_logo_espn), size = 0.05, asp = 16/9)+
  labs(x = "Negative Play Rate", y = "Explosive Play Rate*", title = "Which Offenses Created Explosive Plays and Limited Negative Plays in 2023?",
       subtitle = "The 49ers Created the Most Explosive Plays While the Giants Were Hurt by an Extraordinarily High Negative Play Rate",
       caption = "*Passes that gained greater than 20 yards or runs that gained greater than 12 yards                         Callan Capitolo | @CapAnalytics7 | nflfastR")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())+
  scale_x_reverse()
ggsave("OffExpvsNeg.png", width = 14, height =10, dpi = "retina")

#Drive Efficiency with Negative Play vs Explosive vs none
#Running on 1st down by Yard

nfl99 %>% 
  filter(pass == 1 | rush == 1) %>% 
  mutate(drive_points = ifelse(fixed_drive_result == "Touchdown", 6,
                               ifelse(fixed_drive_result == "Field goal", 3, 0))) %>% 
  group_by(game_id,fixed_drive) %>% 
  summarize(points = max(drive_points), explosive = sum(explosive,na.rm = T)) %>% 
  group_by(explosive) %>% 
  summarize(avg_points = mean(points,na.rm = T), count = n()) %>%
  filter(count>23) %>% 
  ggplot(aes(x = explosive, y = avg_points, fill = "blue"))+
  geom_bar(stat = "identity")+
  geom_text(aes(label = round(avg_points, 2)), vjust = -0.5, color = "black")+
  scale_fill_identity()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())+
  labs(x = "Number of Explosive Plays on Drive*", y = "Average Points Scored on Drive",
       title = "Why Do Explosive Plays Matter?", subtitle = "As the # of explosive plays increase, the average points per drive increases. The biggest jump in average points per drive comes from 0 to 1 explosive plays",
       caption = "*Passes that gained greater than 20 yards or runs that gained greater than 12 yards; data since 1999; TD counts as 6 points (PAT and 2 point attempt excluded)                        @CapAnalytics7 | nflfastR")
ggsave("ExpDrive.png", width = 14, height =10, dpi = "retina")
#explosiveness by season
nfl99 %>% 
  mutate(drive_points = ifelse(fixed_drive_result == "Touchdown", 6,
                               ifelse(fixed_drive_result == "Field goal", 3, 0))) %>% 
  group_by(game_id,fixed_drive) %>% 
  summarize(points = max(drive_points), negative = sum(negative,na.rm = T)) %>% 
  group_by(negative) %>% 
  summarize(avg_points = mean(points,na.rm = T), count = n())

nfl99 %>% 
  group_by(season) %>% 
  summarize(explosive_rate = mean(explosive, na.rm = T)) %>% 
  ggplot(aes(x = season, y = explosive_rate))+
  geom_line(lwd = 2)+
  ylim(0.06,0.08)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())+
  labs(y = "Explosive Play Rate*", x = "Season", title = "How Have Explosive Play Rates Trended?",
       subtitle = "Explosive play rate has seen a major jump since the early 2000s",
       caption = "*Passes that gained greater than 20 yards or runs that gained greater than 12 yards; data since 1999; TD counts as 6 points (PAT and 2 point attempt excluded)                        @CapAnalytics7 | nflfastR")
ggsave("ExpDrive.png", width = 14, height =10, dpi = "retina")


nfl99 %>% 
  filter(play_type %in% c("pass","run")) %>% 
  group_by(season,as.factor(play_type)) %>% 
  summarize(explosive_rate = mean(explosive, na.rm = T)) %>% 
  ggplot(aes(x = season, y = explosive_rate, color = `as.factor(play_type)`))+
  geom_line(lwd = 2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())+
  scale_color_manual(values = c("pass" = "red", "run" = "black"))+
  labs(y = "Explosive Play Rate*", x = "Season", title = "How Have Explosive Play Rates Trended by Play Type?",
       subtitle = "Both runs and passes have a steep increase in the late 2000s, but pass plays see a much sharper increase in explosiveness", color = "Play Type",
       caption = "*Passes that gained greater than 20 yards or runs that gained greater than 12 yards; data since 1999; TD counts as 6 points (PAT and 2 point attempt excluded)                        @CapAnalytics7 | nflfastR") 
ggsave("ExpPlayType.png", width = 14, height =10, dpi = "retina")  

nfl99 %>% 
  group_by(season) %>% 
  summarize(aDot = mean(air_yards, na.rm = T)) %>% 
  ggplot(aes(x = season, y = aDot))+
  geom_line(lwd = 2)+
  ylim(7,9)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())+
  labs(x = "Season", y = "Average Depth of Target", title = "How Has Average Depth of Target Trended Since 2006?",
       subtitle = "Depth of target has been decreasing since the early 2010s", caption = "@CapAnalytics7 | nflfastR")+
  xlim(2006,2023)
ggsave("ADotTrend.png", width = 14, height =10, dpi = "retina") 

nfl99 %>% 
  filter(air_yards<=0) %>% 
  group_by(season) %>% 
  summarize(explosive = mean(explosive, na.rm = T), count = n()) %>% 
  ggplot(aes(x = season, y = explosive, fill = "blue"))+
  scale_fill_identity()+
  geom_bar(stat = "identity")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())+
  labs(x = "Season", y = "Explosive Play Rate*", title = "How has explosiveness on passes at or behind the line of scrimmage changed?",
       subtitle = "Explosive pass rate on passes at or behind the line of scrimmage have dipped in recent seasons", caption = "*Passes that gained greater than 20 yards or runs that gained greater than 12 yards;                      @CapAnalytics7 | nflfastR")+
  xlim(2006,2023)
ggsave("LOSExplosivePassRate.png", width = 14, height =10, dpi = "retina") 

nfl99 %>% 
  mutate(depth = ifelse(air_yards<=0,"At Or Behind LOS",ifelse(air_yards>0 & air_yards<=10,"0-10 Yards",ifelse(air_yards > 10 & air_yards <=20,"11-20 yards", ">20 yards")))) %>% 
  group_by(season,depth) %>% 
  filter(!is.na(depth)) %>% 
  summarize(explosive = mean(explosive, na.rm = T), count = n()) %>% 
  ggplot(aes(x = season, y = explosive, color = depth))+
  geom_line(lwd = 2)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank())+
  labs(x = "Season", y = "Explosive Play Rate*", title = "How has Explosiveness by Depth of Target Changed?",
       subtitle = "There has been an increase in explosiveness on deeper targets in recent years", caption = "*Passes that gained greater than 20 yards or runs that gained greater than 12 yards;                      @CapAnalytics7 | nflfastR",
       color = "Depth of Target")+
  xlim(2006,2023)
ggsave("Depth_Exp.png", width = 14, height =10, dpi = "retina") 



nfl99 %>% 
  filter(complete_pass == 1,season > 2013) %>% 
  group_by(air_yards) %>%
  summarize(avg_yards_gained = mean(yards_gained,na.rm = T), count = n()) %>% filter(count >600) %>% 
  ggplot(aes(x = air_yards, y = avg_yards_gained))+
  geom_line()+labs(x = "Air Yards")

explosive_wins <- nfl99 %>% 
  group_by(game_id) %>% 
  summarize(home = sum(explosive[posteam_type == "home"], na.rm = T), away = sum(explosive[posteam_type == "away"], na.rm = T),
            result = max(result), spread = max(spread_line), more_exp = ifelse(away == home, NA, ifelse(away > home, "away", "home")) ,exp_won = ifelse((result < 0 & away>home) | (result>0 & home>away),1,ifelse(home == away | result == 0,NA,0)),
            exp_cover = ifelse((more_exp == "home" & result > spread) | (more_exp == "away" & result < spread) ,1,0)) %>% 
  summarize(mean(exp_won, na.rm = T), mean(exp_cover, na.rm = T))



negative_wins <- nfl99 %>% 
  group_by(game_id) %>% 
  summarize(home = sum(negative[posteam_type == "home"], na.rm = T), away = sum(negative[posteam_type == "away"], na.rm = T),
            result = max(result), spread = max(spread_line), less_neg = ifelse(away == home, NA, ifelse(away < home, "away", "home")) ,neg_won = ifelse((result < 0 & away<home) | (result>0 & home<away),1,ifelse(home == away | result == 0,NA,0)),
            neg_cover = ifelse((less_neg == "home" & result > spread) | (less_neg == "away" & result < spread) ,1,0)) %>% 
  summarize(mean(neg_won, na.rm = T), mean(neg_cover, na.rm = T), mean(interception[pass_attempt == 1]))


#Thoughts about iterating through the explosive play threshold (plays > 10,11,12,13 yard etc)

#Radar Plots ----

qbs <- nfl99 %>% filter(season >=2022) %>% 
  group_by(id,season) %>% 
  summarize(name = max(name),`EPA/Pass` = mean(epa[pass == 1 & qb_scramble!= 1],na.rm = T), aDoT = mean(air_yards,na.rm = T), `EPA/Rush` = mean(epa[rush == 1 | qb_scramble == 1]), `Success Rate/DB` = mean(success[pass == 1]), CPOE = mean(cpoe,na.rm = T),dropbacks = sum(pass),
            `Sack Rate` = mean(sack[pass==1], na.rm = T), `Explosive DB Rate` = mean(explosive[pass==1], na.rm = T),
            `Int Rate` = mean(interception[!is.na(air_yards)], na.rm = T), `% of Yards From YAC` = sum(yards_after_catch,na.rm = T)/sum(yards_gained[complete_pass == 1],na.rm = T),
             `Comp Rate` = sum(complete_pass,na.rm = T)/sum(!is.na(air_yards)), `Negative DB Rate` = mean(negative[pass == 1], na.rm = T), `Early Down EPA/DB` = mean(epa[pass == 1 & down %in% c(1,2)],na.rm = T)) %>% 
  filter(dropbacks>=100) %>% 
  mutate(name = paste (name,season, sep = " ")) %>% 
  ungroup %>% 
  select(-id,-dropbacks,-season)

replace_with_values_and_ranks <- function(column) {
  values <- column
  ranks <- rank(column*-1,ties.method = "max")
  # reversed_ranks <- max(ranks) + 1 - ranks
  # reversed_ranks
}


qb_rank <- as.data.frame(cbind(qbs$name,(apply(qbs %>% 
                                                 select(-name) %>% 
                                                 mutate(`Sack Rate` = `Sack Rate`*-1,`Int Rate` = `Int Rate`*-1,
                                                        `Negative DB Rate` = `Negative DB Rate`*-1), 2, replace_with_values_and_ranks)))) %>% 
  pivot_longer(cols = -`V1`, names_to = "statistic", values_to = "rank") %>% 
  mutate(rank = as.numeric(rank))
qbs_comb<- qbs %>% pivot_longer(cols = -name, names_to = "statistic", values_to = "value") %>% 
  inner_join(qb_rank, by = c("name" = "V1", "statistic"))


library(ggrepel)
library(cowplot)

temp <- (360/n_distinct(qbs_comb$statistic))/2

myAng <- seq(-temp, -360+temp, length.out = n_distinct(qbs_comb$statistic))

ang <- ifelse(myAng < -90, myAng+180,myAng)

ang <- ifelse(ang < -90, ang+180, ang)

player1 <- "T.Lawrence 2022"
player2 <- "T.Lawrence 2023"

pizza_comp <- qbs_comb %>% 
  filter(statistic != "dropbacks") %>%
  mutate(value = ifelse(statistic == "Int Rate", round(value,3),round(value,2)),
         rank = max(rank) + 1 - rank) %>% 
  filter(name %in% c(player1,player2)) %>% 
  mutate(name = factor(name, levels = c(player1, player2))) %>% 
  ggplot(aes(x = reorder(statistic,value), y = rank, fill = name, label = value))+
  geom_bar(stat = "identity", position = position_identity(), alpha = 0.6)+
  # geom_label(color = "white", size=2.5, fontface="bold", show.legend = FALSE, position = position_jitterdodge())+
  coord_polar()+
  geom_bar(aes(y = max(qbs_comb$rank)/n_distinct(name)),stat = "identity", width =1, alpha = 0.1, fill = "grey")+
  geom_hline(yintercept = seq(1, max(qbs_comb$rank), by = max(qbs_comb$rank)),
             color = "white",
             size = 1)+
  geom_vline(xintercept = seq(.5, n_distinct(qbs_comb$statistic), by = 1),
             color = "white",
             size = .5)+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.background = element_rect(fill = "white", color="white"),
        legend.title = element_blank(),
        legend.text = element_text(colour = "black", face = "bold"),
        plot.title = element_text(hjust = .5, colour = "white", face = "bold", size = 16),
        plot.subtitle = element_text(hjust = .5, colour = "white", size = 8),
        plot.background = element_rect(fill = "black", color="black"),
        panel.background = element_rect(fill = "black", color="black"),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(face = "bold", size = 8, colour = "white"),
        axis.title = element_blank(),
        axis.text.x = element_text(face = "bold", size = 14, angle = ang)) +
  # labs(title = "2023 QB Comparison",
  #      subtitle = "@CapAnalytics7 | nflfastR", x = NULL, y = NULL)+
  labs(x = NULL, y = NULL)+
  scale_fill_brewer(palette = "Set1")
ggsave("pizzachart.png",width = 10, height = 10, dpi ="retina")

tab_comp <- qbs_comb %>%
  filter(name %in% c(player1,player2)) %>% 
  mutate(value = ifelse(statistic == "Int Rate", round(value,3),round(value,3))) %>% 
  mutate(name = factor(name, levels = c(player1, player2))) %>% 
  pivot_wider(names_from = name, values_from = c(value,rank)) %>% 
  rename("Value1" = paste("value_",player1,sep = ""), "Rank1" = paste("rank_",player1,sep = ""), "Value2" = paste("value_",player2,sep = ""),"Rank2" = paste("rank_",player2,sep = "")) %>%
  select(statistic,Value1,Rank1,Value2,Rank2) %>% 
  gt() %>% 
  cols_align(align = "center") %>% 
  tab_spanner(label = player1, columns = c(Value1,Rank1)) %>% 
  tab_spanner(label = player2, columns = c(Value2,Rank2)) %>% 
  cols_label(Value1 = "Value", Rank1 = "Rank",Value2 = "Value",Rank2 = "Rank") %>% 
  gtExtras::gt_theme_538() %>% 
  tab_options(
    table.background.color = "black",        # Set the entire table background to black
    heading.background.color = "black",      # Set the header background to black
    column_labels.background.color = "black", # Set the column label background to black
    row_group.background.color = "black",    # Set row group background to black (if any)
    summary_row.background.color = "black",  # Set summary row background to black (if any)
    grand_summary_row.background.color = "black", # Set grand summary row background to black (if any)
    footnotes.background.color = "black",    # Set footnotes background to black
    source_notes.background.color = "black", # Set source notes background to black
    table.border.top.color = "black",        # Set table top border to black
    table.border.bottom.color = "black",     # Set table bottom border to black
    heading.border.bottom.color = "black",   # Set header bottom border to black
    column_labels.border.top.color = "black",# Set column label top border to black
    column_labels.border.bottom.color = "black" # Set column label bottom border to black
  ) %>% 
  gt_hulk_col_numeric(columns = c(Rank1,Rank2),reverse = TRUE) %>%
  tab_style(
    style = cell_text(size = px(16), weight = "bold", color = "white"),  # Change font and size for column labels
    locations = cells_column_labels(columns = everything())
  ) %>% 
  tab_style(
    style = cell_text(color = "white", size = px(16)),  # Change font and size for the body text
    locations = cells_body(columns = c(statistic, Value1, Value2))
  )
library(chromote)
if (exists("f") && inherits(f, "ChromoteSession")) {
  try(f$shutdown(), silent = TRUE)
}

# Start a new session
f <- ChromoteSession$new()

gtsave(tab_comp, "temp_table.png")

table_image <- image_read("temp_table.png")
table_image_transparent <- image_transparent(table_image, "white")
table_grob <- rasterGrob(table_image, interpolate = TRUE)

library(glue)
spacer <- ggplot() + theme_void() + theme(panel.background = element_rect(fill = "black"))

pizza_comp + table_grob+spacer + plot_layout(ncol = 3, widths = c(6,3,.1))& 
  theme(
    plot.background = element_rect(fill = "black", color = "black"),
    panel.background = element_rect(fill = "black", color = "black"),
    plot.caption = element_text(size = 10, hjust = 0.5, face = "italic", color = "white"),
    plot.title = element_text(size =20,hjust = 0.5, face = "italic", color = "white"),
    plot.subtitle = element_text(hjust = 0.5, face = "bold", color = "white"))&
  # plot_annotation(
  #   caption = glue("Compared to {n_distinct(qbs_comb$name)} QBs with 100+ dropbacks in 2023"),
  #   title = "2023 QB Comparison",
  #   subtitle =  "@CapAnalytics7 | nflfastR"
  # )
  plot_annotation(
    caption = glue("Compared to {n_distinct(qbs_comb$name)} QBs seasons with 100+ dropbacks in 2022-2023"),
    title = "2022-2023 QB Comparison",
    subtitle =  "@CapAnalytics7 | nflfastR"
  )
ggsave("QB_Comp2.png", bg = "black", ,width = 14, height =10, dpi = "retina") 
advanced_stats <- nfl99 %>% filter(season == 2023) %>% 
  group_by(id) %>% 
  summarize(name = max(name), dropbacks = sum(pass), `EPA Outside Num Throws` = mean(epa[pass_location %in% c("left","right")],na.rm = T),
            `EPA Btw Num Throws` = mean(epa[pass_location == "middle"], na.rm = T), `Scramble Rate` = mean(qb_scramble[pass == 1]),`EPA/Exp Pass DB` = mean(epa[xpass>0.85 & pass==1], na.rm = T),
            `EPA/Late Down DB` = mean(epa[down %in% c(3,4)], na.rm = T), `EPA/4Q DB` = mean(epa[qtr>=4 & pass == 1]), epa_trailing = mean(epa[posteam_score<defteam_score],na.rm = T),
            `EPA/Garbage Time DB` = mean(epa[(def_wp>0.9 | def_wp < 0.1) & pass == 1],na.rm = T), `EPA/Non Garbage Time DB` = mean(epa[(def_wp <= 0.9 | def_wp >= 0.1) & pass == 1],na.rm = T), 
            `EPA/Shotgun DB` = mean(epa[pass == 1 & shotgun == 1], na.rm = T), `EPA/Under Center DB` = mean(epa[pass == 1 & shotgun == 0],na.rm =T)) %>% 
  filter(dropbacks > 100)


