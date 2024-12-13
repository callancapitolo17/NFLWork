library(caret)
library(tidyverse)
library(data.table)
library(randomForest)
library(xgboost)
library(Metrics)
library(goftest)
library(fitdistrplus)
library(nflfastR)
library(ggplot2)
library(tidyverse)
library(ggimage)
library(ggthemes)
library(ggrepel)
library(nflreadr)
library(gt)
library(ggrepel)
library(gtExtras)
library(dplyr)
#Data Prep----
ftn_data <- nflreadr::load_ftn_charting(2022:2024) %>%
  select(-week, -season)
wr_long_data <- load_pbp(2022:2024) %>% left_join(ftn_data, by = c("game_id" = "nflverse_game_id",
                                                                              "play_id" = "nflverse_play_id")) 
game_wr <- wr_long_data %>% 
  filter(!is.na(air_yards)) %>% 
  group_by(receiver_player_id,game_id) %>% 
  summarize(name = first(receiver_player_name),posteam = max(posteam),longest_rec = max(yards_gained[!is.na(air_yards)],na.rm = T), week = max(week), targets = n(),
            total_air_yards = sum(air_yards,na.rm = T), total_yards_gained = sum(yards_gained), YAC = sum(yards_after_catch,na.rm = T),
            epa_target = mean(epa,na.rm = T), total_epa = sum(epa,na.rm = T), success_rate = mean(success,na.rm = T), spread = max(spread_line),
            bet_total = max(total_line), adot = mean(air_yards,na.rm = T), year = max(season), defteam = max(defteam),
            catchable_targ = sum(is_catchable_ball,na.rm = T), contested_ball = sum(is_contested_ball),
            created_rec = sum(is_created_reception,na.rm = T)) %>% 
  mutate(longest_rec = ifelse(longest_rec< -100,0,longest_rec)) %>% 
  filter(!is.na(name)) %>% 
  group_by(game_id,posteam) %>% 
  mutate(pct_air_yards = total_air_yards/sum(total_air_yards), yards_gained_share = total_yards_gained/sum(total_yards_gained), target_share = targets/sum(targets))


season_wr <- game_wr %>%
  group_by(receiver_player_id, posteam) %>%
  arrange(game_id,week) %>% # Ensure the data is sorted by week and game
  mutate(
    targets_per_game = cummean(targets),
    air_yards_per_game = cummean(total_air_yards),
    yards_game_per_game = cummean(total_yards_gained),
    yac__per_game = cummean(YAC),
    epa_per_game = cummean(total_epa),
    success_rate__per_game = cummean(success_rate), 
    longest_rec__per_game = cummean(longest_rec),
    pct_air_yards__per_game = cummean(pct_air_yards), # Weighted average try?
    target_share__per_game = cummean(target_share),
    adot_per_game = cummean(adot),
    cum_targets = cumsum(targets),
    cum_target_rate = cumsum(catchable_targ)/cumsum(targets),
    cum_contested_rate = cumsum(contested_ball)/cumsum(targets),
    cum_created_reception_rate = cumsum(created_rec)/cumsum(targets)
    
  ) %>%
  ungroup() %>% 
  group_by(receiver_player_id) %>%
  mutate(next_week = lead(week)) %>% 
  filter(next_week >1) %>% 
  filter(!is.na(week)) %>% 
  select(week = next_week,receiver_player_id,year ,posteam,ends_with("per_game"),starts_with("cum"))



off_scout <- wr_long_data %>%
  filter(!is.na(air_yards)) %>% 
  group_by(posteam, week,game_id,season) %>%
  summarize(
    plays_in_week = n(),                # Total plays in the week
    total_epa_in_week = sum(epa, na.rm = TRUE), # Total EPA in the week
    total_dropbacks = sum(pass,na.rm = T),
    total_db_epa = sum(epa[pass == 1],na.rm =T),
    total_db_success = sum(success[pass==1], na.rm = T),
    total_success = sum(success,na.rm = T),
    total_proe = sum(pass_oe,na.rm = T),
    off_game_longest_rec = max(yards_gained[!is.na(air_yards)]),
    .groups = "drop"
  ) %>%
  group_by(posteam,season) %>%
  arrange(game_id, week) %>%
  mutate(
    # Cumulative EPA
    off_cum_avg_epa_per_play = cumsum(total_epa_in_week) /cumsum(plays_in_week), # Cumulative avg EPA/play
    off_cum_avg_epa_per_db = cumsum(total_db_epa) /cumsum(total_dropbacks),
    off_cum_avg_success_per_db = cumsum(total_db_success) /cumsum(total_dropbacks),
    off_cum_avg_success_rate = cumsum(total_success) /cumsum(plays_in_week),
    off_cum_proe = cumsum(total_proe) /cumsum(plays_in_week),
    off_cum_longest_rec_pg = cummean(off_game_longest_rec)
  ) %>% 
  group_by(posteam) %>%
  mutate(next_week = lead(week)) %>% 
  select(posteam, week = next_week,starts_with("off_cum"),year = season) %>% 
  filter(!is.na(week))

def_scout <- wr_long_data %>%
  filter(!is.na(air_yards)) %>% 
  group_by(defteam, week,game_id,season) %>%
  summarize(
    plays_in_week = n(),                # Total plays in the week
    total_epa_in_week = sum(epa, na.rm = TRUE), # Total EPA in the week
    total_dropbacks = sum(pass,na.rm = T),
    total_db_epa = sum(epa[pass == 1],na.rm =T),
    total_db_success = sum(success[pass==1], na.rm = T),
    total_success = sum(success,na.rm = T),
    total_proe = sum(pass_oe,na.rm = T),
    def_game_longest_rec = max(yards_gained[!is.na(air_yards)]),
    .groups = "drop"
  ) %>%
  group_by(defteam,season) %>%
  arrange(game_id, week) %>%
  mutate(
    # Cumulative EPA
    def_cum_avg_epa_per_play = cumsum(total_epa_in_week) /cumsum(plays_in_week), # Cumulative avg EPA/play
    def_cum_avg_epa_per_db = cumsum(total_db_epa) /cumsum(total_dropbacks),
    def_cum_avg_success_per_db = cumsum(total_db_success) /cumsum(total_dropbacks),
    def_cum_avg_success_rate = cumsum(total_success) /cumsum(plays_in_week),
    def_cum_proe = cumsum(total_proe) /cumsum(plays_in_week),
    def_cum_longest_rec_pg = cummean(def_game_longest_rec)
  ) %>% 
  group_by(defteam) %>%
  mutate(next_week = lead(week)) %>% 
  select(defteam, week = next_week,starts_with("def_cum"),year = season) %>% 
  filter(!is.na(week))



longest_rec_data <- game_wr %>% 
  filter(longest_rec >0 & targets <=3) %>% 
  select(receiver_player_id,name,game_id,week,posteam,spread,bet_total,longest_rec,year,defteam)

longest_rec_joined_data <- longest_rec_data %>% 
  left_join(season_wr, by = c("receiver_player_id","week","year","posteam")) %>% 
  left_join(off_scout, by = c("week","posteam","year")) %>% 
  left_join(def_scout, by = c("week","defteam","year")) %>% 
  ungroup() %>% 
  select(-posteam,-defteam) %>% 
  filter(week>1)
   #Filter out 0

checks <- longest_rec_joined_data %>% 
  filter(longest_rec == 0)

test <- longest_rec_joined_data %>%
  group_by(longest_rec) %>% 
  summarize(count = n())
  

longest_rec_joined_data %>% 
  ggplot(aes(x =sqrt(longest_rec)))+
  geom_density()

tree_predictions_vector <- sqrt(longest_rec_joined_data$longest_rec)
# tree_predictions_vector <- ifelse(tree_predictions_vector==0,0.0000000000000001,tree_predictions_vector)

library(goftest)

# Extract the individual tree predictions

# Load the fitdistrplus package
library(fitdistrplus)

fit_normal <- fitdist(tree_predictions_vector, "norm")
fit_lognorm <- fitdist(tree_predictions_vector,"lnorm")
fit_gamma <- fitdist(tree_predictions_vector, "gamma")
fit_weibull <- fitdist(tree_predictions_vector, "weibull")
# fit_beta <- fitdist(tree_predictions_vector / max(tree_predictions_vector), "beta")
fit_exp = fitdist(tree_predictions_vector, "exp")
# fit_t <- fitdist(tree_predictions_vector, "t")

# Compare AIC and BIC for all fitted distributions
aic_values <- c(
  normal = fit_normal$aic,
  lognormal = fit_lognorm$aic,
  gamma = fit_gamma$aic,
  weibull = fit_weibull$aic,
  # beta = fit_beta$aic,
  exp = fit_exp$aic
  # t = fit_t$aic
)

bic_values <- c(
  normal = BIC(fit_normal),
  lognormal = BIC(fit_lognorm),
  gamma = BIC(fit_gamma),
  weibull = BIC(fit_weibull),
  # beta = BIC(fit_beta),
  exp = BIC(fit_exp)
  # t = BIC(fit_t)
)

aic_values
bic_values

plot.legend <- c("Normal", "Log-Normal", "Gamma", "Weibull", "Exp")
denscomp(list(fit_normal, fit_lognorm, fit_gamma, fit_weibull), legendtext = plot.legend)

qqcomp(list(fit_normal, fit_lognorm, fit_gamma, fit_weibull), legendtext = plot.legend)
ppcomp(list(fit_normal, fit_lognorm, fit_gamma, fit_weibull), legendtext = plot.legend)
# Perform KS test for the gamma fit
set.seed(123)  # For reproducibility

ks.test(
  tree_predictions_vector,
  "pgamma",
  rate = fit_gamma$estimate["rate"],
  shape = fit_gamma$estimate["shape"]
)

# Perform AD test for the normal fit
ad.test(
  tree_predictions_vector,
  pgamma,
  rate = fit_gamma$estimate["rate"],
  shape = fit_gamma$estimate["shape"]
)

# Define bins (intervals for the test)
bins <- hist(tree_predictions_vector, breaks = 10, plot = FALSE)$breaks

# Calculate expected frequencies for normal distribution
expected <- diff(pnorm(bins, mean = fit_normal$estimate["mean"], sd = fit_normal$estimate["sd"])) * length(tree_predictions_vector)

# Calculate observed frequencies
observed <- hist(tree_predictions_vector, breaks = bins, plot = FALSE)$counts

# Perform chi-squared test
chisq.test(observed, p = expected / sum(expected))

# Generate the fitted gamma density
gamma_density <- data.frame(
  x = seq(min(tree_predictions_vector), max(tree_predictions_vector), length.out = 100),
  y = dgamma(seq(min(tree_predictions_vector), max(tree_predictions_vector), length.out = 100),
             shape = fit_gamma$estimate["shape"],
             rate = fit_gamma$estimate["rate"])
)

# Plot the histogram and fitted gamma density
ggplot(data.frame(tree_predictions_vector), aes(x = tree_predictions_vector)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "blue", alpha = 0.4) +
  geom_line(data = gamma_density, aes(x = x, y = y), color = "red", size = 1) +
  labs(title = "Gamma Distribution Fit",
       x = "Predicted Values",
       y = "Density")

# Generate Q-Q plot for gamma
qqcomp(list(fit_gamma), legendtext = c("Gamma"))

# Create an empirical CDF and compare with gamma CDF
ecdf_data <- ecdf(tree_predictions_vector)
plot(ecdf_data, main = "CDF Comparison", xlab = "Predicted Values", ylab = "Cumulative Probability")
curve(pgamma(x, shape = fit_gamma$estimate["shape"], rate = fit_gamma$estimate["rate"]),
      add = TRUE, col = "red", lwd = 2)
legend("bottomright", legend = c("Empirical CDF", "Theoretical Gamma CDF"), col = c("black", "red"), lty = 1)

# Residual Analysis using gamma fit
observed <- tree_predictions_vector
expected <- qgamma(
  p = ecdf(tree_predictions_vector)(tree_predictions_vector),
  shape = fit_gamma$estimate["shape"],
  rate = fit_gamma$estimate["rate"]
)
shape <- fit_gamma$estimate["shape"]
rate <- fit_gamma$estimate["rate"]

# Calculate fitted (expected) values as the mean of the gamma distribution
fitted <- qgamma(
  p = pgamma(tree_predictions_vector, shape = shape, rate = rate), 
  shape = shape, 
  rate = rate
)

residuals <- observed - fitted

# Plot 1: Residuals vs Observed Values
ggplot(data.frame(observed, residuals), aes(x = observed, y = residuals)) +
  geom_point(alpha = 0.5, color = "blue") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Residuals vs Observed Values (Gamma)",
       x = "Observed Values",
       y = "Residuals") +
  theme_minimal()

# Plot 2: Histogram of Residuals
ggplot(data.frame(residuals), aes(x = residuals)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "blue", alpha = 0.4) +
  geom_density(color = "red", size = 1) +
  labs(title = "Histogram of Residuals (Gamma)",
       x = "Residuals",
       y = "Density") +
  theme_minimal()


#Longest Completion Transformed Model----

# Load required libraries
library(caret)
library(tidyverse)
library(data.table)
library(randomForest)
library(xgboost)
library(Metrics)

# Data Preparation
set.seed(42) # Set seed for reproducibility

# Clean and preprocess data
sqrt_data <- longest_rec_joined_data %>%
  drop_na() %>%  # Remove rows with missing values
  mutate_if(is.character, as.factor) %>%  # Convert character columns to factors
  select(-receiver_player_id, -game_id, -name) %>%  # Drop non-predictive or irrelevant columns
  mutate(longest_rec = sqrt(ifelse(longest_rec<0,0,longest_rec)))
# Split data into training and testing sets
sqrt_trainIndex <- createDataPartition(sqrt_data$longest_rec, p = 0.8, list = FALSE)
sqrt_trainData <- sqrt_data[sqrt_trainIndex, ]
sqrt_testData <- sqrt_data[-sqrt_trainIndex, ]

# Separate predictors and response
sqrt_x_train <- sqrt_trainData %>% select(-longest_rec)
sqrt_y_train <- sqrt_trainData$longest_rec
sqrt_x_test <- sqrt_testData %>% select(-longest_rec)
sqrt_y_test <- sqrt_testData$longest_rec

# Feature Scaling
sqrt_preProc <- preProcess(sqrt_x_train, method = c("center", "scale"))
sqrt_x_train <- predict(sqrt_preProc, sqrt_x_train)
sqrt_x_test <- predict(sqrt_preProc, sqrt_x_test)

# Model Training and Cross-Validation

# Linear Regression
sqrt_lm_model <- train(longest_rec ~ ., data = sqrt_trainData, method = "lm",
                       trControl = trainControl(method = "cv", number = 5))

# Random Forest
sqrt_rf_grid <- expand.grid(
  mtry = seq(2, ncol(sqrt_x_train), by = 2)  # Test mtry values for the transformed model
)

sqrt_rf_model <- train(
  longest_rec ~ ., 
  data = sqrt_trainData, 
  method = "rf",
  trControl = trainControl(method = "cv", number = 5),  # 5-fold cross-validation
  tuneGrid = sqrt_rf_grid,
  ntree = 500
)

# Gradient Boosting (XGBoost)
sqrt_xgb_grid <- expand.grid(
  nrounds = c(100, 200),
  max_depth = c(3, 5, 7),
  eta = c(0.01, 0.1),
  gamma = 0,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.8
)

sqrt_xgb_model <- train(longest_rec ~ ., data = sqrt_trainData, method = "xgbTree",
                        trControl = trainControl(method = "cv", number = 5),
                        tuneGrid = sqrt_xgb_grid)

# Evaluate Models
sqrt_models <- list(lm = sqrt_lm_model, rf = sqrt_rf_model, xgb = sqrt_xgb_model)
sqrt_results <- resamples(sqrt_models)

# Summarize results
summary(sqrt_results)

# Compare RMSE on the test set
sqrt_model_performance <- sapply(sqrt_models, function(model) {
  predictions <- predict(model, newdata = sqrt_x_test)
  rmse(sqrt_y_test, predictions)
})

# Select the best model
sqrt_best_model_name <- names(sqrt_model_performance)[which.min(sqrt_model_performance)]
sqrt_best_model <- sqrt_models[[sqrt_best_model_name]]

cat("Best model:", sqrt_best_model_name, "with RMSE:", min(sqrt_model_performance), "\n")

# Feature Importance (for Random Forest or XGBoost)
if (sqrt_best_model_name %in% c("rf", "xgb")) {
  varImp(sqrt_best_model) %>% plot()
}

sqrt_final_predictions <- as.data.frame(predict(sqrt_best_model, newdata = sqrt_x_test))
sqrt_final_check<- as.data.frame(predict(sqrt_best_model, newdata = sqrt_x_test, interval = "prediction", level = 0.5))

sqrt_residuals <- sqrt_y_test - sqrt_final_predictions[[1]]

# Plot residuals
plot(sqrt_final_predictions[[1]], sqrt_residuals, main = "Residuals vs Predictions", xlab = "Predicted Values", ylab = "Residuals")
abline(h = 0, col = "red")



#No Transform----
# Load required libraries

# Data Preparation
set.seed(42) # Set seed for reproducibility

# Clean and preprocess data
data <- longest_rec_joined_data %>%
  drop_na() %>%  # Remove rows with missing values
  mutate_if(is.character, as.factor) %>%  # Convert character columns to factors
  select(-receiver_player_id, -game_id, -name) # Drop non-predictive or irrelevant columns

# Split data into training and testing sets
trainIndex <- createDataPartition(data$longest_rec, p = 0.8, list = FALSE)
trainData <- data[trainIndex, ]
testData <- data[-trainIndex, ]

# Separate predictors and response
x_train <- trainData %>% select(-longest_rec)
y_train <- trainData$longest_rec
x_test <- testData %>% select(-longest_rec)
y_test <- testData$longest_rec

# Feature Scaling
preProc <- preProcess(x_train, method = c("center", "scale"))
x_train <- predict(preProc, x_train)
x_test <- predict(preProc, x_test)

# Model Training and Cross-Validation

# Linear Regression
lm_model <- train(longest_rec ~ ., data = trainData, method = "lm",
                  trControl = trainControl(method = "cv", number = 5))

# Random Forest
rf_grid <- expand.grid(
  mtry = seq(2, ncol(x_train), by = 2)  # Test mtry values from 2 to the number of predictors
)

# Retraining Original RF Model
rf_model <- train(
  longest_rec ~ ., 
  data = trainData, 
  method = "rf",
  trControl = trainControl(method = "cv", number = 5),  # 5-fold cross-validation
  tuneGrid = rf_grid,
  ntree = 500
)

# Gradient Boosting (XGBoost)
xgb_grid <- expand.grid(
  nrounds = c(100, 200),
  max_depth = c(3, 5, 7),
  eta = c(0.01, 0.1),
  gamma = 0,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.8
)

xgb_model <- train(longest_rec ~ ., data = trainData, method = "xgbTree",
                   trControl = trainControl(method = "cv", number = 5),
                   tuneGrid = xgb_grid)

# Evaluate Models
models <- list(lm = lm_model, rf = rf_model, xgb = xgb_model)
results <- resamples(models)

# Summarize results
summary(results)

# Compare RMSE on the test set
model_performance <- sapply(models, function(model) {
  predictions <- predict(model, newdata = x_test)
  rmse(y_test, predictions)
})

# Select the best model
best_model_name <- names(model_performance)[which.min(model_performance)]
best_model <- models[[best_model_name]]

cat("Best model:", best_model_name, "with RMSE:", min(model_performance), "\n")

# Feature Importance (for Random Forest or XGBoost)
if (best_model_name %in% c("rf", "xgb")) {
  varImp(best_model) %>% plot()
}
final_predictions <- as.data.frame(predict(best_model, newdata = x_test))
residuals <- y_test - final_predictions[[1]]

# Plot residuals
plot(final_predictions[[1]], residuals, main = "Residuals vs Predictions", xlab = "Predicted Values", ylab = "Residuals")
abline(h = 0, col = "red")

# Final predictions
final_check<- as.data.frame(predict(best_model, newdata = x_test, interval = "prediction", level = 0.5))

test_identifiers <- testData %>% select(receiver_player_id, game_id, name,longest_rec) # Drop non-predictive or irrelevant columns

# Make predictions using the best model

# Combine predictions with identifiers
predictions_with_identifiers <- test_identifiers %>%
  mutate(predicted_longest_rec = final_predictions)

#Model Comparison----
# Reverse the square root transformation
sqrt_final_predictions_original_scale <- sqrt_final_predictions[[1]]^2
sqrt_residuals_original_scale <- sqrt_y_test^2 - sqrt_final_predictions_original_scale

# RMSE and MAE for the original model
original_rmse <- rmse(y_test, final_predictions[[1]])
original_mae <- mae(y_test, final_predictions[[1]])

# RMSE and MAE for the transformed model (on original scale)
sqrt_rmse <- rmse(y_test, sqrt_final_predictions_original_scale)
sqrt_mae <- mae(y_test, sqrt_final_predictions_original_scale)

cat("Original Model RMSE:", original_rmse, "\n")
cat("Original Model MAE:", original_mae, "\n")
cat("Square Root Model RMSE (Original Scale):", sqrt_rmse, "\n")
cat("Square Root Model MAE (Original Scale):", sqrt_mae, "\n")

# Residuals plot for the original model
original_residuals <- y_test - final_predictions[[1]]
plot(final_predictions[[1]], original_residuals, main = "Residuals vs Predictions (Original Model)", xlab = "Predicted Values", ylab = "Residuals")
abline(h = 0, col = "red")

# Residuals plot for the square root-transformed model (on original scale)
plot(sqrt_final_predictions_original_scale, sqrt_residuals_original_scale, main = "Residuals vs Predictions (Sqrt Model)", xlab = "Predicted Values", ylab = "Residuals")
abline(h = 0, col = "red")

# Paired t-test on residuals (original scale)
t.test(abs(original_residuals), abs(sqrt_residuals_original_scale), paired = TRUE)

# Check the best model name
cat("Best model:", sqrt_best_model_name, "\n")

# Extract the best tuning parameters for the model
if (sqrt_best_model_name == "xgb") {
  best_params <- sqrt_xgb_model$bestTune
} else if (sqrt_best_model_name == "rf") {
  best_params <- sqrt_rf_model$bestTune
} else if (sqrt_best_model_name == "lm") {
  best_params <- "No tuning parameters for linear regression"
}

print(best_params)



#Longest Completion Forward Data Prep----

wr_pred_data <- pbp_rp

game_wr <- wr_pred_data %>% 
  filter(!is.na(air_yards)) %>% 
  group_by(receiver_player_id,game_id) %>% 
  summarize(name = first(receiver_player_name),posteam = max(posteam),longest_rec = max(yards_gained[!is.na(air_yards)],na.rm = T), week = max(week), targets = n(),
            total_air_yards = sum(air_yards,na.rm = T), total_yards_gained = sum(yards_gained), YAC = sum(yards_after_catch,na.rm = T),
            epa_target = mean(epa,na.rm = T), total_epa = sum(epa,na.rm = T), success_rate = mean(success,na.rm = T), spread = max(spread_line),
            bet_total = max(total_line), adot = mean(air_yards,na.rm = T), year = max(season), defteam = max(defteam),
            catchable_targ = sum(is_catchable_ball,na.rm = T), contested_ball = sum(is_contested_ball),
            created_rec = sum(is_created_reception,na.rm = T)) %>% 
  mutate(longest_rec = ifelse(longest_rec< -100,0,longest_rec)) %>% 
  filter(!is.na(name)) %>% 
  group_by(game_id,posteam) %>% 
  mutate(pct_air_yards = total_air_yards/sum(total_air_yards), yards_gained_share = total_yards_gained/sum(total_yards_gained), target_share = targets/sum(targets))


season_wr <- game_wr %>%
  group_by(receiver_player_id, posteam,name) %>%
  summarize(
    targets_per_game = mean(targets),
    air_yards_per_game = mean(total_air_yards),
    yards_game_per_game = mean(total_yards_gained),
    yac__per_game = mean(YAC),
    epa_per_game = mean(total_epa),
    success_rate__per_game = mean(success_rate), 
    longest_rec__per_game = mean(longest_rec),
    pct_air_yards__per_game = mean(pct_air_yards), # Weighted average try?
    target_share__per_game = mean(target_share),
    adot_per_game = mean(adot),
    cum_targets = sum(targets),
    cum_target_rate = sum(catchable_targ)/sum(targets),
    cum_contested_rate = sum(contested_ball)/sum(targets),
    cum_created_reception_rate = sum(created_rec)/sum(targets)
    
  ) %>%
  ungroup() %>% 
  select(receiver_player_id,posteam,name,ends_with("per_game"),starts_with("cum"))



off_scout <- wr_pred_data %>%
  group_by(posteam, week,game_id,season) %>%
  summarize(
    plays_in_week = n(),                # Total plays in the week
    total_epa_in_week = sum(epa, na.rm = TRUE), # Total EPA in the week
    total_dropbacks = sum(pass,na.rm = T),
    total_db_epa = sum(epa[pass == 1],na.rm =T),
    total_db_success = sum(success[pass==1], na.rm = T),
    total_success = sum(success,na.rm = T),
    total_proe = sum(pass_oe,na.rm = T),
    off_game_longest_rec = max(yards_gained[!is.na(air_yards)]),
    .groups = "drop"
  ) %>%
  group_by(posteam,season) %>%
  summarize(
    # Cumulative EPA
    off_cum_avg_epa_per_play = sum(total_epa_in_week) /sum(plays_in_week), # Cumulative avg EPA/play
    off_cum_avg_epa_per_db = sum(total_db_epa) /sum(total_dropbacks),
    off_cum_avg_success_per_db = sum(total_db_success) /sum(total_dropbacks),
    off_cum_avg_success_rate = sum(total_success) /sum(plays_in_week),
    off_cum_proe = sum(total_proe) /sum(plays_in_week),
    off_cum_longest_rec_pg = mean(off_game_longest_rec)
  )

def_scout <- wr_pred_data %>%
  group_by(defteam, week,game_id,season) %>%
  summarize(
    plays_in_week = n(),                # Total plays in the week
    total_epa_in_week = sum(epa, na.rm = TRUE), # Total EPA in the week
    total_dropbacks = sum(pass,na.rm = T),
    total_db_epa = sum(epa[pass == 1],na.rm =T),
    total_db_success = sum(success[pass==1], na.rm = T),
    total_success = sum(success,na.rm = T),
    total_proe = sum(pass_oe,na.rm = T),
    def_game_longest_rec = max(yards_gained[!is.na(air_yards)]),
    .groups = "drop"
  ) %>%
  group_by(defteam,season) %>%
  arrange(week) %>%
  summarize(
    # Cumulative EPA
    def_cum_avg_epa_per_play = sum(total_epa_in_week) /sum(plays_in_week), # Cumulative avg EPA/play
    def_cum_avg_epa_per_db = sum(total_db_epa) /sum(total_dropbacks),
    def_cum_avg_success_per_db = sum(total_db_success) /sum(total_dropbacks),
    def_cum_avg_success_rate = sum(total_success) /sum(plays_in_week),
    def_cum_proe = sum(total_proe) /sum(plays_in_week),
    def_cum_longest_rec_pg = mean(def_game_longest_rec)
  ) %>% 
  select(defteam, starts_with("def_cum"))# Remove rows year = seasonwhere the next week doesn't exist


schedules <- load_schedules(2024)

wr_sched <- schedules %>% 
  # filter(is.na(away_score), week == pbp_rp$week +1)
  filter(is.na(away_score), week == 16) %>% 
  select(week, home_team ,spread = spread_line,bet_total = total_line,year = season,away_team) %>% 
  mutate(posteam_home = "yes", posteam_away = "no") %>% 
  pivot_longer(cols = c(posteam_home,posteam_away), values_to = "home_pos_team", names_to = "posteam") %>% 
  mutate(posteam = ifelse(home_pos_team == "yes", home_team, away_team),
         defteam = ifelse(home_pos_team == "yes", away_team, home_team)) %>% 
  select(-home_team,-away_team,-home_pos_team)

longest_rec_joined_data_pred <- wr_sched %>% 
  right_join(season_wr, by = c("posteam")) %>% 
  left_join(off_scout, by = c("posteam")) %>% 
  left_join(def_scout, by = c("defteam")) %>% 
  ungroup() 



#Final Model----

# Combine training and testing data for final retraining
final_data <- rbind(trainData, testData)

# Separate predictors and target variable
final_x <- final_data %>% select(-longest_rec)
final_y <- final_data$longest_rec

# Step 1: Preprocess data (center and scale)
final_preProc <- preProcess(final_x, method = c("center", "scale"))
final_x <- predict(final_preProc, final_x)

# Step 2: Retrain the best model (assume Random Forest is the best)
final_rf_model <- randomForest(
  sqrt(longest_rec) ~ ., 
  data = cbind(final_x, longest_rec = final_y), 
  mtry = 2,  # Best mtry value identified earlier
  ntree = 500,      # Optimal ntree value
  importance = TRUE  # Enable variable importance
)

# Save preprocessing steps and model
saveRDS(final_preProc, file = "final_preprocessing.rds")
saveRDS(final_rf_model, file = "final_rf_model.rds")

# Evaluate the retrained model on an independent test set
processed_test_x <- predict(final_preProc, newdata = x_test)

# Make predictions on the test set
test_predictions <- predict(final_rf_model, newdata = processed_test_x)

# Reverse transformation for predictions (if necessary)
test_predictions_original_scale <- test_predictions^2

# Calculate evaluation metrics
final_rmse <- rmse(y_test, test_predictions_original_scale)
final_mae <- mae(y_test, test_predictions_original_scale)

cat("Final Model RMSE:", final_rmse, "\n")
cat("Final Model MAE:", final_mae, "\n")


#Predictions ----

# Load the saved preprocessing object and final model
final_preProc <- readRDS("final_preprocessing.rds")
final_rf_model <- readRDS("final_rf_model.rds")

# Preprocess new data
processed_new_data <- predict(final_preProc, newdata = longest_rec_joined_data_pred)

# Make predictions
predictions_df <- longest_rec_joined_data_pred %>%  mutate(longest_reception_predctions = predict(final_rf_model, newdata = processed_new_data))

sqrt_longest_rec <- sqrt(longest_rec_joined_data$longest_rec)
fit_gamma <- fitdist(sqrt_longest_rec, "gamma")

# Extract parameters
shape <- fit_gamma$estimate["shape"]
rate <- fit_gamma$estimate["rate"]

market_line <- 15.5

sqrt_market_line <- sqrt(market_line)


# Calculate probabilities for each prediction
predictions_df <- longest_rec_joined_data_pred %>%
  mutate(
    longest_reception_predictions = predict(final_rf_model, newdata = processed_new_data),
    adjusted_rate = shape / longest_reception_predictions,  # rate = shape / mean
    
    # Calculate the probability of exceeding the market line
    prob_exceed_threshold = 1 - pgamma(
      sqrt_market_line,       # Transformed market line
      shape = shape,          # Shape remains constant
      rate = adjusted_rate    # Player-specific rate
    )
  )

check <- predictions_df %>% select(posteam,name, longest_reception_predictions, prob_exceed_threshold)
american_to_decimal_odds <- function(odds){ifelse(odds < 0, 1 - (100/odds), 1+(odds/100))}
prob_to_decimal <- function(prob){ifelse(prob > 0, 1 / prob, NA)}
decimal_to_american = function(fair_decimal_odds)ifelse(
  fair_decimal_odds > 2.0,
  (fair_decimal_odds - 1) * 100,          # Positive American odds
  -100 / (fair_decimal_odds - 1)          # Negative American odds
)
ev <- function(predicted_prob,decimal_odds){(predicted_prob * (decimal_odds-1)) - (1 - predicted_prob)}
longest_reception_prediction_output <- function(){
  player <- readline(prompt = "Enter Player Name:")
  betting_line <- as.numeric(readline(prompt = "Enter Line:"))
  over_odds <- as.numeric(readline(promp = "Enter Over Odds:"))
  under_odds <- as.numeric(readline(promp = "Enter Under Odds:"))
  over_dec <- american_to_decimal_odds(over_odds)
  under_dec <- american_to_decimal_odds(under_odds)
  over_prob <- 1/over_dec
  under_prob <- 1/under_dec
  results <- predictions_df %>% 
    filter(name == player) %>% 
    mutate(betting_line = betting_line,
      adjusted_rate = shape / longest_reception_predictions,
      predicted_over_prob = 1 - pgamma(
      sqrt(betting_line),       # Transformed market line
      shape = shape,          # Shape remains constant
      rate = adjusted_rate    # Player-specific rate
    ),
    predicted_under_prob = 1-predicted_over_prob,
    pred_over_dec_odds = prob_to_decimal(predicted_over_prob),
    pred_over_american_odds = decimal_to_american(pred_over_dec_odds),
    `Over ROI %` = ev(predicted_over_prob,over_dec)*100,
    pred_under_dec_odds = prob_to_decimal(predicted_under_prob),
    pred_under_american_odds = decimal_to_american(pred_under_dec_odds),
    `Under ROI %` = ev(predicted_under_prob,under_dec)*100,
    pred_longest_rec = longest_reception_predictions^2) %>% 
    select(posteam,name, pred_longest_rec,betting_line, predicted_over_prob,`Over ROI %`, `Under ROI %`)
  return(results)
}
longest_reception_prediction_output()
