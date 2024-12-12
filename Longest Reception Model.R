#Data Prep----
wr_long_data <- nfl99 %>% filter(season>=2022) %>% left_join(ftn_data, by = c("game_id" = "nflverse_game_id",
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
  select(receiver_player_id,name,game_id,week, posteam,spread,bet_total,longest_rec,year,defteam)

longest_rec_joined_data <- longest_rec_data %>% 
  left_join(season_wr, by = c("receiver_player_id","week","year","posteam")) %>% 
  left_join(off_scout, by = c("week","posteam","year")) %>% 
  left_join(def_scout, by = c("week","defteam","year")) %>% 
  ungroup() %>% 
  select(-posteam,-defteam) %>% 
  filter(week>1)
#No Transform----
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
rf_model <- train(longest_rec ~ ., data = trainData, method = "rf",
                  trControl = trainControl(method = "cv", number = 5),
                  tuneLength = 5)

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

# Load necessary libraries


tree_predictions_vector <- predict(best_model, x_test, predict.all = TRUE)

# Extract the individual tree predictions

# Load the fitdistrplus package
library(fitdistrplus)

# Fit various distributions
# Perform KS test for the Weibull fit
# Replace gamma-specific parts with log-normal equivalents

# Perform KS test for the log-normal fit
ks.test(
  tree_predictions_vector,
  "plnorm",
  meanlog = fit_lognorm$estimate["meanlog"],
  sdlog = fit_lognorm$estimate["sdlog"]
)

# Perform AD test for the log-normal fit
ad.test(
  tree_predictions_vector,
  plnorm,
  meanlog = fit_lognorm$estimate["meanlog"],
  sdlog = fit_lognorm$estimate["sdlog"]
)

# Define bins (intervals for the test)
bins <- hist(tree_predictions_vector, breaks = 10, plot = FALSE)$breaks

# Calculate expected frequencies for log-normal
expected <- diff(plnorm(bins, meanlog = fit_lognorm$estimate["meanlog"], sdlog = fit_lognorm$estimate["sdlog"])) * length(tree_predictions_vector)

# Calculate observed frequencies
observed <- hist(tree_predictions_vector, breaks = bins, plot = FALSE)$counts

# Perform chi-squared test
chisq.test(observed, p = expected / sum(expected))

# Generate the fitted log-normal density
lognorm_density <- data.frame(
  x = seq(min(tree_predictions_vector), max(tree_predictions_vector), length.out = 100),
  y = dlnorm(seq(min(tree_predictions_vector), max(tree_predictions_vector), length.out = 100),
             meanlog = fit_lognorm$estimate["meanlog"],
             sdlog = fit_lognorm$estimate["sdlog"])
)

# Plot the histogram and fitted density
ggplot(data.frame(tree_predictions_vector), aes(x = tree_predictions_vector)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "blue", alpha = 0.4) +
  geom_line(data = lognorm_density, aes(x = x, y = y), color = "red", size = 1) +
  labs(title = "Log-Normal Distribution Fit",
       x = "Predicted Values",
       y = "Density")

# Generate Q-Q plot for log-normal
qqcomp(list(fit_lognorm), legendtext = c("Log-Normal"))

# Create an empirical CDF and compare with log-normal CDF
ecdf_data <- ecdf(tree_predictions_vector)

# Plot the empirical and theoretical CDFs
plot(ecdf_data, main = "CDF Comparison", xlab = "Predicted Values", ylab = "Cumulative Probability")
curve(plnorm(x, meanlog = fit_lognorm$estimate["meanlog"], sdlog = fit_lognorm$estimate["sdlog"]),
      add = TRUE, col = "red")
legend("bottomright", legend = c("Empirical CDF", "Theoretical Log-Normal CDF"), col = c("black", "red"), lty = 1)

# Residual Analysis using log-normal fit
observed <- tree_predictions_vector
expected <- qlnorm(
  p = ecdf(tree_predictions_vector)(tree_predictions_vector),
  meanlog = fit_lognorm$estimate["meanlog"],
  sdlog = fit_lognorm$estimate["sdlog"]
)

residuals <- observed - expected

# Plot 1: Residuals vs. Observed Values
ggplot(data.frame(observed, residuals), aes(x = observed, y = residuals)) +
  geom_point(alpha = 0.5, color = "blue") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Residuals vs Observed Values (Log-Normal)",
       x = "Observed Values",
       y = "Residuals") +
  theme_minimal()

# Plot 2: Histogram of Residuals
ggplot(data.frame(residuals), aes(x = residuals)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "blue", alpha = 0.4) +
  geom_density(color = "red", size = 1) +
  labs(title = "Histogram of Residuals (Log-Normal)",
       x = "Residuals",
       y = "Density") +
  theme_minimal()

# Plot 3: Residuals vs Fitted Values
fitted <- expected
ggplot(data.frame(fitted, residuals), aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.5, color = "blue") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Residuals vs Fitted Values (Log-Normal)",
       x = "Fitted Values",
       y = "Residuals") +
  theme_minimal()

#Issues with residuals



#Longest Completion Transformed Model----


log_transformed <- data %>% 
  mutate(longest_rec = ifelse(longest_rec <=0,0.1,longest_rec))
sqrt_transformed <- data %>% 
  mutate(longest_rec = ifelse(longest_rec <=0,0,longest_rec))
par(mfrow = c(1, 3))
qqnorm(data$longest_rec, main = "Original")
qqline(data$longest_rec, col = "red")

qqnorm(log(log_transformed$longest_rec), main = "Log-Transformed")
qqline(log(log_transformed$longest_rec), col = "blue")

qqnorm(sqrt(sqrt_transformed$longest_rec), main = "Square Root-Transformed")
qqline(sqrt(sqrt_transformed$longest_rec), col = "green")

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
data <- longest_rec_joined_data %>%
  drop_na() %>%  # Remove rows with missing values
  mutate_if(is.character, as.factor) %>%  # Convert character columns to factors
  select(-receiver_player_id, -game_id, -name) %>%  # Drop non-predictive or irrelevant columns
  mutate(longest_rec = sqrt(ifelse(longest_rec<0,0,longest_rec)))
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
rf_model <- train(longest_rec ~ ., data = trainData, method = "rf",
                  trControl = trainControl(method = "cv", number = 5),
                  tuneLength = 5)

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

residuals <- y_test - final_predictions[[1]]

# Plot residuals
plot(final_predictions[[1]], residuals, main = "Residuals vs Predictions", xlab = "Predicted Values", ylab = "Residuals")
abline(h = 0, col = "red")

# Final predictions
final_predictions <- as.data.frame(predict(best_model, newdata = x_test))
final_check<- as.data.frame(predict(best_model, newdata = x_test, interval = "prediction", level = 0.5))

test_identifiers <- testData %>% select(receiver_player_id, game_id, name,longest_rec) # Drop non-predictive or irrelevant columns

# Make predictions using the best model

# Combine predictions with identifiers
predictions_with_identifiers <- test_identifiers %>%
  mutate(predicted_longest_rec = final_predictions)

# Load necessary libraries


tree_predictions_vector <- predict(rf_model, x_test, predict.all = TRUE)

# Extract the individual tree predictions

# Load the fitdistrplus package
library(fitdistrplus)

# Fit various distributions
# Perform KS test for the Weibull fit
# Replace gamma-specific parts with log-normal equivalents

# Perform KS test for the log-normal fit
ks.test(
  tree_predictions_vector,
  "plnorm",
  meanlog = fit_lognorm$estimate["meanlog"],
  sdlog = fit_lognorm$estimate["sdlog"]
)

# Perform AD test for the log-normal fit
ad.test(
  tree_predictions_vector,
  plnorm,
  meanlog = fit_lognorm$estimate["meanlog"],
  sdlog = fit_lognorm$estimate["sdlog"]
)

# Define bins (intervals for the test)
bins <- hist(tree_predictions_vector, breaks = 10, plot = FALSE)$breaks

# Calculate expected frequencies for log-normal
expected <- diff(plnorm(bins, meanlog = fit_lognorm$estimate["meanlog"], sdlog = fit_lognorm$estimate["sdlog"])) * length(tree_predictions_vector)

# Calculate observed frequencies
observed <- hist(tree_predictions_vector, breaks = bins, plot = FALSE)$counts

# Perform chi-squared test
chisq.test(observed, p = expected / sum(expected))

# Generate the fitted log-normal density
lognorm_density <- data.frame(
  x = seq(min(tree_predictions_vector), max(tree_predictions_vector), length.out = 100),
  y = dlnorm(seq(min(tree_predictions_vector), max(tree_predictions_vector), length.out = 100),
             meanlog = fit_lognorm$estimate["meanlog"],
             sdlog = fit_lognorm$estimate["sdlog"])
)

# Plot the histogram and fitted density
ggplot(data.frame(tree_predictions_vector), aes(x = tree_predictions_vector)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "blue", alpha = 0.4) +
  geom_line(data = lognorm_density, aes(x = x, y = y), color = "red", size = 1) +
  labs(title = "Log-Normal Distribution Fit",
       x = "Predicted Values",
       y = "Density")

# Generate Q-Q plot for log-normal
qqcomp(list(fit_lognorm), legendtext = c("Log-Normal"))

# Create an empirical CDF and compare with log-normal CDF
ecdf_data <- ecdf(tree_predictions_vector)

# Plot the empirical and theoretical CDFs
plot(ecdf_data, main = "CDF Comparison", xlab = "Predicted Values", ylab = "Cumulative Probability")
curve(plnorm(x, meanlog = fit_lognorm$estimate["meanlog"], sdlog = fit_lognorm$estimate["sdlog"]),
      add = TRUE, col = "red")
legend("bottomright", legend = c("Empirical CDF", "Theoretical Log-Normal CDF"), col = c("black", "red"), lty = 1)

# Residual Analysis using log-normal fit
observed <- tree_predictions_vector
expected <- qlnorm(
  p = ecdf(tree_predictions_vector)(tree_predictions_vector),
  meanlog = fit_lognorm$estimate["meanlog"],
  sdlog = fit_lognorm$estimate["sdlog"]
)

residuals <- log(observed) - log(expected)

# Plot 1: Residuals vs. Observed Values
ggplot(data.frame(log(observed), residuals), aes(x = observed, y = residuals)) +
  geom_point(alpha = 0.5, color = "blue") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Residuals vs Observed Values (Log-Normal)",
       x = "Observed Values",
       y = "Residuals") +
  theme_minimal()

# Plot 2: Histogram of Residuals
ggplot(data.frame(residuals), aes(x = residuals)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "blue", alpha = 0.4) +
  geom_density(color = "red", size = 1) +
  labs(title = "Histogram of Residuals (Log-Normal)",
       x = "Residuals",
       y = "Density") +
  theme_minimal()

# Plot 3: Residuals vs Fitted Values
fitted <- expected
ggplot(data.frame(fitted, residuals), aes(x = fitted, y = residuals)) +
  geom_point(alpha = 0.5, color = "blue") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Residuals vs Fitted Values (Log-Normal)",
       x = "Fitted Values",
       y = "Residuals") +
  theme_minimal()

#Issues with residuals






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
  filter(is.na(away_score), week == 15) %>% 
  select(week, home_team ,spread = spread_line,bet_total = total_line,year = season,away_team) %>% 
  mutate(posteam_home = "yes", posteam_away = "no") %>% 
  pivot_longer(cols = c(posteam_home,posteam_away), values_to = "home_pos_team", names_to = "posteam") %>% 
  mutate(posteam = ifelse(home_pos_team == "yes", home_team, away_team),
         defteam = ifelse(home_pos_team == "yes", away_team, home_team)) %>% 
  select(-home_team,-away_team,-home_pos_team)

longest_rec_joined_data <- wr_sched %>% 
  right_join(season_wr, by = c("posteam")) %>% 
  left_join(off_scout, by = c("posteam")) %>% 
  left_join(def_scout, by = c("defteam")) %>% 
  ungroup() 

week_predictions <- as.data.frame(predict(best_model, newdata = longest_rec_joined_data))[[1]]

longest_rec_predictions <- longest_rec_joined_data %>% 
  filter(!is.na(spread)) %>% 
  select(name) %>% 
  mutate(pred_longest_rec = week_predictions)

