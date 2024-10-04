library(ggimage)
library(gt)
library(nflfastR)
library(tidyverse)
library(png)
library(nflreadr)
library(tidyverse)   # Data manipulation
library(caret)       # Model training and tuning
library(randomForest) # Random Forest
library(xgboost)     
library(dplyr)
library(caret)

sched <- load_schedules(2010:2024)

model_update_data <- sched %>% 
  mutate(favorite = ifelse(spread_line >= 0, "home", "away"), favorite_cover = ifelse((favorite == "home" & result > spread_line-6) |(favorite == "away" & result < spread_line+6),1,0 )) %>% 
  mutate(underdog = ifelse(spread_line >= 0, "away", "home"), underdog_cover = ifelse((underdog == "home" & result > spread_line-6) |(underdog == "away" & result < spread_line+6),1,0 )) %>% 
  pivot_longer(cols = c("favorite_cover", "underdog_cover"), values_to = "cover", names_to = "spread_type") %>% 
  mutate(teaser_under = total_line+6, teaser_over = total_line-6) %>% 
  pivot_longer(cols = c("teaser_under", "teaser_over"), values_to = "teaser_total", names_to = "total_type") %>% 
  mutate(place = ifelse(location == "Neutral",location,ifelse(spread_type == "favorite_cover",favorite, underdog))) %>% 
  mutate(spread_line = ifelse(spread_type == "underdog_cover" & spread_line < 0,spread_line*-1, ifelse(spread_type == "favorite_cover" & spread_line > 0, spread_line*-1,spread_line))) %>% 
  mutate(total_cover = ifelse((total_type == "teaser_over" & teaser_total < total) | (total_type == "teaser_under" & teaser_total > total),1,0),
         combo_cover = ifelse(cover == 1 & total_cover == 1,1,0)) %>% 
  select(total_line,total_type,spread_line,place,away_spread_odds,home_spread_odds,under_odds,over_odds,favorite,underdog,combo_cover,season,week) %>% 
  mutate(away_prob = 1/ifelse(away_spread_odds <0,(100/abs(away_spread_odds)) +1, (abs(away_spread_odds)/100) + 1),
  home_prob = 1/ifelse(home_spread_odds <0,(100/abs(home_spread_odds)) +1, (abs(home_spread_odds)/100) + 1),
  away_vig_free = away_prob/(home_prob+away_prob),
  home_vig_free = home_prob/(home_prob+away_prob),
  under_prob = 1/ifelse(under_odds <0,(100/abs(under_odds)) +1, (abs(under_odds)/100) + 1),
  over_prob = 1/ifelse(over_odds <0,(100/abs(over_odds)) +1, (abs(over_odds)/100) + 1),,
  under_vig_free = under_prob/(over_prob+under_prob),
  over_vig_free = over_prob/(over_prob+under_prob)) %>% 
  mutate(spread_odds = ifelse(spread_line < 0 & favorite == "away",away_vig_free, 
                              ifelse(spread_line < 0 & favorite == "home",home_vig_free,ifelse(spread_line > 0 & underdog == "away",away_vig_free,home_vig_free))),
         total_odds = ifelse(total_type == "teaser_under",under_vig_free,over_vig_free)) %>% 
  select(total_line,total_type,spread_line,place,spread_odds,total_odds, combo_cover,season,week) %>% 
  drop_na()

#Model Testing----
seasons <- c(2015:2015)
units_year <- numeric(length(seasons))
threshold <- numeric(length(seasons))

for (szn in seasons) {
  
  edge_model_data <- model_update_data %>% 
    filter(season >= szn) %>%
    mutate(combo_cover = ifelse(combo_cover == 1, "yes","no"))
  
  # Gradient Boosting (XGBoost)
  
  # Set a random seed for reproducibility
  set.seed(123)
  
  # Assume `edge_model_data` has been preprocessed as per previous steps
  
  # Convert categorical variables to factors
  edge_model_data$total_type <- as.factor(edge_model_data$total_type)
  edge_model_data$combo_cover <- as.factor(edge_model_data$combo_cover)
  
  # Split data into training and testing sets
  train_index <- createDataPartition(edge_model_data$combo_cover, p = 0.8, list = FALSE)
  train_data <- edge_model_data %>% filter(season < 2023) 
  # %>% select(-season)
  test_data <- edge_model_data %>% filter(season==2023) 
  # %>% select(-season)
  
  # Define cross-validation method with more detailed summary metrics
  cv_control <- trainControl(method = "cv", number = 5, classProbs = TRUE,
                             summaryFunction = twoClassSummary, savePredictions = "final")
xgb_grid_expanded <- expand.grid(
  nrounds = c(25,50, 100, 150), 
  max_depth = c(3,4,5,6), #depth of tree
  eta = c(0.01, 0.1, 0.3),#learning rate
  gamma = c(0, 1),
  colsample_bytree = c(0.5, 0.75, 1),
  min_child_weight = c(1, 5, 10),
  subsample = c(0.5, 0.75, 1)
)
xgb_model_expanded <- train(combo_cover ~ ., data = train_data, method = "xgbTree",
                            trControl = cv_control, tuneGrid = xgb_grid_expanded, metric = "ROC")
  
  # Best Model Performance on Test Set
  best_model_expanded <- xgb_model_expanded  # Assuming XGBoost performs best
  predictions_expanded <- predict(best_model_expanded, test_data)
  

  prob_pred <- predict(best_model_expanded, test_data, type = "prob")$yes
  prob_threshold <- seq(0.45,0.85,0.005)
  winnings <- numeric(length(prob_threshold))
  for (i in c(1:length(prob_threshold))) {
    temp_units<-cbind(test_data,prob_pred) %>% 
      mutate(pred_bin = ifelse(prob_pred >prob_threshold[i],"yes","no")) %>% 
      filter(pred_bin == "yes") %>%
      mutate(correct = ifelse(pred_bin == combo_cover,1,0)) %>% 
      mutate(units = ifelse(correct == 1, 0.91,-1))
    winnings[i] <- sum(temp_units$units)
  }
  units_year[(szn-min(seasons)+1)] = max(winnings)
  threshold[(szn-min(seasons)+1)] <- prob_threshold[which.max(winnings)]
  # as.data.frame(cbind(prob_threshold,winnings)) %>% 
  #   ggplot(aes(x = prob_threshold, y = winnings, fill = "lightblue"))+
  #   geom_bar(stat = "identity")+
  #   geom_text(label = round(winnings,2))+
  #   labs(x = "Probability Threshold", y = "Units +/-")+
  #   theme_bw()+
  #   theme(panel.grid.major = element_blank(),  # Remove major gridlines
  #         panel.grid.minor = element_blank())+
  #   scale_fill_identity()
  # Print confusion matrix
  # print(conf_matrix_expanded)
  # variable_importances2<- varImp(best_model_expanded)[[1]]
  # variable_importances_exp2 <- cbind(Variables = rownames(variable_importances),variable_importances)
  # var_imp_comb2 <- if(szn == min(seasons)) {variable_importances_exp} else{var_imp_comb %>%
  #     left_join(variable_importances_exp, by = "Variables")}
  # 
  # colnames(var_imp_comb2)[(szn-min(seasons)+2)] <- toString(szn)
}
units_winning2<- cbind(seasons,threshold, units_year)

#Finalize Model ----

best_params <- best_model_expanded$bestTune

best_params_manual_upd <- expand.grid(
  nrounds = 150,         # Number of boosting rounds
  max_depth = 6,        # Maximum tree depth
  eta = 0.3,            # Learning rate
  gamma = 0,            # Minimum loss reduction
  colsample_bytree = 0.75, # Subsample ratio of columns when constructing each tree
  min_child_weight = 1, # Minimum sum of instance weight (hessian) needed in a child
  subsample = 1      # Subsample ratio of the training instance
)
cv_control <- trainControl(method = "cv", number = 5, classProbs = TRUE,
                           summaryFunction = twoClassSummary, savePredictions = "final")

edge_model_data <- model_update_data %>% 
  filter(season >= 2015) %>%
  mutate(combo_cover = ifelse(combo_cover == 1, "yes","no")) %>% 
  select(total_line,total_type,spread_line,place,spread_odds,total_odds, combo_cover,season,week) %>% 
  drop_na()

# Convert categorical variables to factors
edge_model_data$total_type <- as.factor(edge_model_data$total_type)
edge_model_data$combo_cover <- as.factor(edge_model_data$combo_cover)

set.seed(123)
xgb_model_prod <- train(combo_cover ~ ., data = edge_model_data, method = "xgbTree",
                        trControl = cv_control, tuneGrid = best_params_manual_upd, metric = "ROC")

sched24 <- load_schedules(2024)
# 
# teaser_sched <- sched24 %>% 
#   mutate(favorite = ifelse(spread_line >= 0, "home", "away"), favorite_cover = ifelse((favorite == "home" & result > spread_line-6) |(favorite == "away" & result < spread_line+6),1,0 )) %>% 
#   mutate(underdog = ifelse(spread_line >= 0, "away", "home"), underdog_cover = ifelse((underdog == "home" & result > spread_line-6) |(underdog == "away" & result < spread_line+6),1,0 )) %>% 
#   pivot_longer(cols = c("favorite_cover", "underdog_cover"), values_to = "cover", names_to = "spread_type") %>% 
#   mutate(teaser_under = total_line+6, teaser_over = total_line-6) %>% 
#   pivot_longer(cols = c("teaser_under", "teaser_over"), values_to = "teaser_total", names_to = "total_type") %>% 
#   mutate(place = ifelse(location == "Neutral",location,ifelse(spread_type == "favorite_cover",favorite, underdog))) %>% 
#   mutate(spread_line = ifelse(spread_type == "underdog_cover" & spread_line < 0,spread_line*-1, ifelse(spread_type == "favorite_cover" & spread_line > 0, spread_line*-1,spread_line))) %>% 
#   filter(week == 2) 

teaser_sched <- sched24 %>% 
  mutate(favorite = ifelse(spread_line >= 0, "home", "away"), favorite_cover = ifelse((favorite == "home" & result > spread_line-6) |(favorite == "away" & result < spread_line+6),1,0 )) %>% 
  mutate(underdog = ifelse(spread_line >= 0, "away", "home"), underdog_cover = ifelse((underdog == "home" & result > spread_line-6) |(underdog == "away" & result < spread_line+6),1,0 )) %>% 
  pivot_longer(cols = c("favorite_cover", "underdog_cover"), values_to = "cover", names_to = "spread_type") %>% 
  mutate(teaser_under = total_line+6, teaser_over = total_line-6) %>% 
  pivot_longer(cols = c("teaser_under", "teaser_over"), values_to = "teaser_total", names_to = "total_type") %>% 
  mutate(place = ifelse(location == "Neutral",location,ifelse(spread_type == "favorite_cover",favorite, underdog))) %>% 
  mutate(spread_line = ifelse(spread_type == "underdog_cover" & spread_line < 0,spread_line*-1, ifelse(spread_type == "favorite_cover" & spread_line > 0, spread_line*-1,spread_line))) %>% 
  mutate(away_prob = 1/ifelse(away_spread_odds <0,(100/abs(away_spread_odds)) +1, (abs(away_spread_odds)/100) + 1),
         home_prob = 1/ifelse(home_spread_odds <0,(100/abs(home_spread_odds)) +1, (abs(home_spread_odds)/100) + 1),
         away_vig_free = away_prob/(home_prob+away_prob),
         home_vig_free = home_prob/(home_prob+away_prob),
         under_prob = 1/ifelse(under_odds <0,(100/abs(under_odds)) +1, (abs(under_odds)/100) + 1),
         over_prob = 1/ifelse(over_odds <0,(100/abs(over_odds)) +1, (abs(over_odds)/100) + 1),,
         under_vig_free = under_prob/(over_prob+under_prob),
         over_vig_free = over_prob/(over_prob+under_prob)) %>% 
  mutate(spread_odds = ifelse(spread_line < 0 & favorite == "away",away_vig_free, 
                              ifelse(spread_line < 0 & favorite == "home",home_vig_free,ifelse(spread_line > 0 & underdog == "away",away_vig_free,home_vig_free))),
         total_odds = ifelse(total_type == "teaser_under",under_vig_free,over_vig_free)) %>% 
  filter(week == 5)

prediction_data <- teaser_sched %>% 
  select(total_line,total_type,spread_line,place,spread_odds,total_odds,season,week)
  # rename("bet_total" = "total_line","spread" = "spread_line")

pred_1 <- cbind(teaser_sched,predict(xgb_model_prod,prediction_data,type = "prob")$yes) %>% 
  rename("prob" = "predict(xgb_model_prod, prediction_data, type = \"prob\")$yes") %>% 
  rename("side" = "place") %>% 
  select(away_team,home_team, total_line,total_type,spread_line,side,prob) %>% 
  mutate(place = ifelse(prob > 0.62,"bet","no"))


