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
nfl99all <- load_pbp(1999:2024)
nfl99 <- nfl99all %>% 
  filter(pass == 1 | rush == 1) %>% 
  mutate(explosive = ifelse((yards_gained>20 & pass_attempt == 1) | (yards_gained >12 & (qb_scramble == 1 | rush == 1)),1,0),
         negative = ifelse(yards_gained < 0, 1,0))

teaser_data <- nfl99 %>%
  group_by(game_id) %>% 
  summarize(year = max(season), homescore = max(home_score), away_score = max(away_score),spread = max(spread_line), result = max(result), bet_total = max(total_line), result_total = max(total),location = max(location),
            max(home_team), max(away_team), week = max(week)) %>% 
  mutate(favorite = ifelse(spread >= 0, "home", "away"), favorite_cover = ifelse((favorite == "home" & result > spread-6) |(favorite == "away" & result < spread+6),1,0 )) %>% 
  mutate(underdog = ifelse(spread >= 0, "away", "home"), underdog_cover = ifelse((underdog == "home" & result > spread-6) |(underdog == "away" & result < spread+6),1,0 ))

# SGP Teaser ----
edge_data <- teaser_data %>% 
  pivot_longer(cols = c("favorite_cover", "underdog_cover"), values_to = "cover", names_to = "spread_type") %>% 
  mutate(teaser_under = bet_total+6, teaser_over = bet_total-6) %>% 
  pivot_longer(cols = c("teaser_under", "teaser_over"), values_to = "teaser_total", names_to = "total_type") %>% 
  mutate(team = ifelse(location == "Neutral",location,ifelse(spread_type == "favorite_cover",favorite, underdog))) %>% 
  mutate(total_cover = ifelse((total_type == "teaser_over" & teaser_total < result_total) | (total_type == "teaser_under" & teaser_total > result_total),1,0),
         combo_cover = ifelse(cover == 1 & total_cover == 1,1,0)) %>% 
  mutate(spread = ifelse(spread_type == "underdog_cover" & spread < 0,spread*-1, ifelse(spread_type == "favorite_cover" & spread > 0, spread*-1,spread)))

edge_sum <- edge_data %>% 
 group_by(total_type,bet_total,spread) %>% 
  summarize(combination_cover_rate = mean(combo_cover,na.rm = T), occurrences = n()) %>%
  filter(occurrences > 20) %>% 
  mutate(total_type = ifelse(total_type == "teaser_under","under","over")) %>% 
  mutate(combination_cover_rate = round(combination_cover_rate,2)) %>% 
  filter(combination_cover_rate > 0.56) %>% 
  arrange(-combination_cover_rate) %>% 
  ungroup() %>% 
  gt() %>% 
  cols_align(align = "center") %>% 
  cols_label( total_type = "Total Type",
    bet_total = "Bet Total",
              spread = "Spread",
              combination_cover_rate = "Cover Rate",
              occurrences = "Occurrences") %>%
  tab_header(
    title = md("Top 20 Spread/Total Combination Cover Rates"),
    subtitle = md("Breakeven win rate is 54.55%; data since 2015")) %>% 
  gtExtras::gt_theme_538()
gtsave(edge_sum,"edge_sum.png")

seasons <- c(2015:2022)
units_year <- numeric(length(seasons))
threshold <- numeric(length(seasons))

for (season in seasons) {

edge_model_data <- edge_data %>% 
  filter(year >= 2022) %>%
  mutate(combo_cover = ifelse(combo_cover == 1, "yes","no")) %>% 
  select(bet_total,total_type,spread,combo_cover,team,year,week)

# Gradient Boosting (XGBoost)

# Set a random seed for reproducibility
set.seed(123)

# Assume `edge_model_data` has been preprocessed as per previous steps

# Convert categorical variables to factors
edge_model_data$total_type <- as.factor(edge_model_data$total_type)
edge_model_data$combo_cover <- as.factor(edge_model_data$combo_cover)

# Split data into training and testing sets
train_index <- createDataPartition(edge_model_data$combo_cover, p = 0.8, list = FALSE)
train_data <- edge_model_data %>% filter(year < 2023) %>% select(-year)
test_data <- edge_model_data %>% filter(year==2023) %>% select(-year)

# Define cross-validation method with more detailed summary metrics
cv_control <- trainControl(method = "cv", number = 5, classProbs = TRUE, 
                       summaryFunction = twoClassSummary, savePredictions = "final")

# Random Forest Model with expanded grid
# rf_grid_expanded <- expand.grid(mtry = c(2,3) 
#                                 # ntrees = c(500,1000),
#                                 # nodesize = c(5,10)
#                                 )
# rf_model_expanded <- train(combo_cover ~ ., data = train_data, method = "rf",
#                            trControl = cv_control, tuneGrid = rf_grid_expanded, metric = "ROC")

# XGBoost Model with expanded grid
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

# Compare model accuracies
# results_expanded <- resamples(list(RandomForest = rf_model_expanded, XGBoost = xgb_model_expanded))
# summary(results_expanded)

# Plotting model comparison
# bwplot(results_expanded)

# Best Model Performance on Test Set
best_model_expanded <- xgb_model_expanded  # Assuming XGBoost performs best
predictions_expanded <- predict(best_model_expanded, test_data)
# conf_matrix_expanded <- confusionMatrix(predictions_expanded, test_data$combo_cover)
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
  units_year[(season-min(seasons)+1)] = max(winnings)
  threshold[(season-min(seasons)+1)] <- prob_threshold[which.max(winnings)]
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
 variable_importances<- varImp(best_model_expanded)[[1]]
 variable_importances_exp <- cbind(Variables = rownames(variable_importances),variable_importances)
 var_imp_comb <- if(season == min(seasons)) {variable_importances_exp} else{var_imp_comb %>%
     left_join(variable_importances_exp, by = "Variables")}

 colnames(var_imp_comb)[(season-min(seasons)+2)] <- toString(season)
}
units_winning<- cbind(seasons,threshold, units_year)

best_params <- best_model_expanded$bestTune

best_params_manual <- expand.grid(
  nrounds = 25,         # Number of boosting rounds
  max_depth = 3,        # Maximum tree depth
  eta = 0.3,            # Learning rate
  gamma = 0,            # Minimum loss reduction
  colsample_bytree = 0.5, # Subsample ratio of columns when constructing each tree
  min_child_weight = 1, # Minimum sum of instance weight (hessian) needed in a child
  subsample = 0.5      # Subsample ratio of the training instance
)

edge_model_data <- edge_data %>% 
  filter(year >= 2022) %>%
  mutate(combo_cover = ifelse(combo_cover == 1, "yes","no")) %>% 
  select(bet_total,total_type,spread,combo_cover,team)

# Convert categorical variables to factors
edge_model_data$total_type <- as.factor(edge_model_data$total_type)
edge_model_data$combo_cover <- as.factor(edge_model_data$combo_cover)

set.seed(123)
xgb_model_prod <- train(combo_cover ~ ., data = edge_model_data, method = "xgbTree",
                            trControl = cv_control, tuneGrid = best_params_manual, metric = "ROC")

sched24 <- load_schedules(2024)

teaser_sched <- sched24 %>% 
  mutate(favorite = ifelse(spread_line >= 0, "home", "away"), favorite_cover = ifelse((favorite == "home" & result > spread_line-6) |(favorite == "away" & result < spread_line+6),1,0 )) %>% 
  mutate(underdog = ifelse(spread_line >= 0, "away", "home"), underdog_cover = ifelse((underdog == "home" & result > spread_line-6) |(underdog == "away" & result < spread_line+6),1,0 )) %>% 
  pivot_longer(cols = c("favorite_cover", "underdog_cover"), values_to = "cover", names_to = "spread_type") %>% 
  mutate(teaser_under = total_line+6, teaser_over = total_line-6) %>% 
  pivot_longer(cols = c("teaser_under", "teaser_over"), values_to = "teaser_total", names_to = "total_type") %>% 
  mutate(team = ifelse(location == "Neutral",location,ifelse(spread_type == "favorite_cover",favorite, underdog))) %>% 
  mutate(spread_line = ifelse(spread_type == "underdog_cover" & spread_line < 0,spread_line*-1, ifelse(spread_type == "favorite_cover" & spread_line > 0, spread_line*-1,spread_line))) %>% 
  filter(week == 2) 

prediction_data <- teaser_sched %>% 
  select(total_line,total_type,spread_line,team) %>% 
  rename("bet_total" = "total_line","spread" = "spread_line")

pred_1 <- cbind(teaser_sched,predict(xgb_model_prod,prediction_data,type = "prob")$yes) %>% 
  rename("prob" = "predict(xgb_model_prod, prediction_data, type = \"prob\")$yes") %>% 
  rename("side" = "team") %>% 
  select(away_team,home_team, total_line,total_type,spread_line,side,prob) %>% 
  mutate(place = ifelse(prob > 0.575,"bet","no"))

predict(xgb_model_expanded,c())

