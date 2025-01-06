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
library(kernlab)
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

#Distribution Fitting----

longest_rec_joined_data %>% 
  ggplot(aes(x =(longest_rec)))+
  geom_density()

tree_predictions_vector <- (longest_rec_joined_data$longest_rec)
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

cross_validate_distribution <- function(data, target_var, dist_name, k_folds = 5) {
  # Create folds
  folds <- createFolds(data[[target_var]], k = k_folds, list = TRUE)
  
  # Initialize metrics storage
  log_likelihoods <- numeric(k_folds)
  
  # Cross-validation loop
  for (i in seq_along(folds)) {
    # Split data
    train_data <- data[-folds[[i]], ]
    test_data <- data[folds[[i]], ]
    
    # Fit distribution to training data
    fit <- fitdist(train_data[[target_var]], dist_name)
    
    # Predict probabilities for test data
    test_values <- test_data[[target_var]]
    
    if (dist_name == "gamma") {
      predicted_probs <- dgamma(test_values, 
                                shape = fit$estimate["shape"], 
                                rate = fit$estimate["rate"])
    } else if (dist_name == "norm") {
      predicted_probs <- dnorm(test_values, 
                               mean = fit$estimate["mean"], 
                               sd = fit$estimate["sd"])
    } else if (dist_name == "lnorm") {
      predicted_probs <- dlnorm(test_values, 
                                meanlog = fit$estimate["meanlog"], 
                                sdlog = fit$estimate["sdlog"])
    } else {
      stop("Unsupported distribution type.")
    }
    
    # Compute log-likelihood for the test set
    log_likelihoods[i] <- sum(log(predicted_probs + 1e-9))  # Adding a small value to avoid log(0)
  }
  
  # Return average log-likelihood across folds
  mean(log_likelihoods)
}



# Example: Assuming your target variable is sqrt(longest_rec)
data <- longest_rec_joined_data %>% 
  drop_na()

# Compare distributions using 5-fold CV
gamma_cv_ll <- cross_validate_distribution(data, "longest_rec", "gamma", k_folds = 5)
normal_cv_ll <- cross_validate_distribution(data, "longest_rec", "norm", k_folds = 5)

cat("Gamma Distribution Log-Likelihood (CV):", gamma_cv_ll, "\n")
cat("Normal Distribution Log-Likelihood (CV):", normal_cv_ll, "\n")
lognormal_cv_ll <- cross_validate_distribution(data, "longest_rec", "lnorm", k_folds = 5)
cat("Log-Normal Distribution Log-Likelihood (CV):", lognormal_cv_ll, "\n")

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

#Weighted Average----
season_weighted <- game_wr %>%
  group_by(receiver_player_id, posteam) %>%
  arrange(game_id, week) %>% # Ensure the data is sorted by week
  mutate(
    max_week = max(week), # Find the most recent week
    week_weight = (week - min(week) + 1) / (max_week - min(week) + 1), # Scale weights to give higher value to recent weeks
    
    # Weighted averages for per-game metrics
    weighted_targets_per_game = cumsum(targets * week_weight) / cumsum(week_weight),
    weighted_air_yards_per_game = cumsum(total_air_yards * week_weight) / cumsum(week_weight),
    weighted_yards_game_per_game = cumsum(total_yards_gained * week_weight) / cumsum(week_weight),
    weighted_yac_per_game = cumsum(YAC * week_weight) / cumsum(week_weight),
    weighted_epa_per_game = cumsum(total_epa * week_weight) / cumsum(week_weight),
    weighted_success_rate_per_game = cumsum(success_rate * week_weight) / cumsum(week_weight),
    weighted_longest_rec_per_game = cumsum(longest_rec * week_weight) / cumsum(week_weight),
    weighted_pct_air_yards_per_game = cumsum(pct_air_yards * week_weight) / cumsum(week_weight),
    weighted_target_share_per_game = cumsum(target_share * week_weight) / cumsum(week_weight),
    weighted_adot_per_game = cumsum(adot * week_weight) / cumsum(week_weight),
    
    # Cumulative metrics
    cum_targets = cumsum(targets),
    cum_target_rate = cumsum(catchable_targ) / cumsum(targets),
    cum_contested_rate = cumsum(contested_ball) / cumsum(targets),
    cum_created_reception_rate = cumsum(created_rec) / cumsum(targets)
  ) %>%
  ungroup() %>% 
  group_by(receiver_player_id) %>%
  mutate(next_week = lead(week)) %>% 
  filter(next_week > 1) %>% 
  filter(!is.na(week)) %>% 
  select(
    week = next_week, receiver_player_id, year, posteam,
    ends_with("per_game"), starts_with("cum")
  )
longest_rec_joined_weighted_data <- longest_rec_data %>% 
  left_join(season_weighted, by = c("receiver_player_id","week","year","posteam")) %>% 
  left_join(off_scout, by = c("week","posteam","year")) %>% 
  left_join(def_scout, by = c("week","defteam","year")) %>% 
  ungroup() %>% 
  select(-posteam,-defteam) %>% 
  filter(week>1)

weighted_data <- longest_rec_joined_weighted_data %>%
  drop_na() %>%  # Remove rows with missing values
  mutate_if(is.character, as.factor) %>%  # Convert character columns to factors
  select(-receiver_player_id, -game_id, -name) %>%  # Drop non-predictive or irrelevant columns
  mutate(longest_rec = sqrt(ifelse(longest_rec<0,0,longest_rec)))
# Split data into training and testing sets
weighted_trainIndex <- createDataPartition(weighted_data$longest_rec, p = 0.8, list = FALSE)
weighted_trainData <- weighted_data[weighted_trainIndex, ]
weighted_testData <- weighted_data[-weighted_trainIndex, ]

# Separate predictors and response
weighted_x_train <- weighted_trainData %>% select(-longest_rec)
weighted_y_train <- weighted_trainData$longest_rec
weighted_x_test <- weighted_testData %>% select(-longest_rec)
weighted_y_test <- weighted_testData$longest_rec

# Feature Scaling
weighted_preProc <- preProcess(weighted_x_train, method = c("center", "scale"))
weighted_x_train <- predict(weighted_preProc, weighted_x_train)
weighted_x_test <- predict(weighted_preProc, weighted_x_test)

# Model Training and Cross-Validation

# Linear Regression
weighted_lm_model <- train(longest_rec ~ ., data = weighted_trainData, method = "lm",
                       trControl = trainControl(method = "cv", number = 5))

# Random Forest
weighted_rf_grid <- expand.grid(
  mtry = seq(2, ncol(weighted_x_train), by = 2)  # Test mtry values for the transformed model
)

weighted_rf_model <- train(
  longest_rec ~ ., 
  data = weighted_trainData, 
  method = "rf",
  trControl = trainControl(method = "cv", number = 5),  # 5-fold cross-validation
  tuneGrid = weighted_rf_grid,
  ntree = 500
)

# Gradient Boosting (XGBoost)
weighted_xgb_grid <- expand.grid(
  nrounds = c(100, 200),
  max_depth = c(3, 5, 7),
  eta = c(0.01, 0.1),
  gamma = 0,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.8
)

weighted_xgb_model <- train(longest_rec ~ ., data = weighted_trainData, method = "xgbTree",
                        trControl = trainControl(method = "cv", number = 5),
                        tuneGrid = weighted_xgb_grid)

# Evaluate Models
weighted_models <- list(lm = weighted_lm_model, rf = weighted_rf_model, xgb = weighted_xgb_model)
weighted_results <- resamples(weighted_models)

# Summarize results
summary(weighted_results)

# Compare RMSE on the test set
weighted_model_performance <- sapply(weighted_models, function(model) {
  predictions <- predict(model, newdata = weighted_x_test)
  rmse(weighted_y_test, predictions)
})

# Select the best model
weighted_best_model_name <- names(weighted_model_performance)[which.min(weighted_model_performance)]
weighted_best_model <- weighted_models[[weighted_best_model_name]]

cat("Best model:", weighted_best_model_name, "with RMSE:", min(weighted_model_performance), "\n")

# Feature Importance (for Random Forest or XGBoost)
if (weighted_best_model_name %in% c("rf", "xgb")) {
  varImp(weighted_best_model) %>% plot()
}

weighted_final_predictions <- as.data.frame(predict(weighted_best_model, newdata = weighted_x_test))
weighted_final_check<- as.data.frame(predict(weighted_best_model, newdata = weighted_x_test, interval = "prediction", level = 0.5))

weighted_residuals <- weighted_y_test - weighted_final_predictions[[1]]

# Plot residuals
plot(weighted_final_predictions[[1]], weighted_residuals, main = "Residuals vs Predictions", xlab = "Predicted Values", ylab = "Residuals")
abline(h = 0, col = "red")
#Season and Weighted Model ----
season_weighted_combo <- game_wr %>%
  group_by(receiver_player_id, posteam) %>%
  arrange(game_id, week) %>% # Ensure the data is sorted by week
  mutate(
    max_week = max(week), # Find the most recent week
    week_weight = (week - min(week) + 1) / (max_week - min(week) + 1), # Scale weights to give higher value to recent weeks
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
    # Weighted averages for per-game metrics
    weighted_targets_per_game = cumsum(targets * week_weight) / cumsum(week_weight),
    weighted_air_yards_per_game = cumsum(total_air_yards * week_weight) / cumsum(week_weight),
    weighted_yards_game_per_game = cumsum(total_yards_gained * week_weight) / cumsum(week_weight),
    weighted_yac_per_game = cumsum(YAC * week_weight) / cumsum(week_weight),
    weighted_epa_per_game = cumsum(total_epa * week_weight) / cumsum(week_weight),
    weighted_success_rate_per_game = cumsum(success_rate * week_weight) / cumsum(week_weight),
    weighted_longest_rec_per_game = cumsum(longest_rec * week_weight) / cumsum(week_weight),
    weighted_pct_air_yards_per_game = cumsum(pct_air_yards * week_weight) / cumsum(week_weight),
    weighted_target_share_per_game = cumsum(target_share * week_weight) / cumsum(week_weight),
    weighted_adot_per_game = cumsum(adot * week_weight) / cumsum(week_weight),
    
    # Cumulative metrics
    cum_targets = cumsum(targets),
    cum_target_rate = cumsum(catchable_targ) / cumsum(targets),
    cum_contested_rate = cumsum(contested_ball) / cumsum(targets),
    cum_created_reception_rate = cumsum(created_rec) / cumsum(targets)
  ) %>%
  ungroup() %>% 
  group_by(receiver_player_id) %>%
  mutate(next_week = lead(week)) %>% 
  filter(next_week > 1) %>% 
  filter(!is.na(week)) %>% 
  select(
    week = next_week, receiver_player_id, year, posteam,
    ends_with("per_game"), starts_with("cum")
  )
longest_rec_joined_combo_data <- longest_rec_data %>% 
  left_join(season_weighted_combo, by = c("receiver_player_id","week","year","posteam")) %>% 
  left_join(off_scout, by = c("week","posteam","year")) %>%
  left_join(def_scout, by = c("week","defteam","year")) %>% 
  ungroup() %>% 
  select(-posteam,-defteam) %>% 
  filter(week>1)
combo_data <- longest_rec_joined_combo_data %>%
  drop_na() %>%  # Remove rows with missing values
  mutate_if(is.character, as.factor) %>%  # Convert character columns to factors
  select(-receiver_player_id, -game_id, -name) %>%  # Drop non-predictive or irrelevant columns
  mutate(longest_rec = sqrt(ifelse(longest_rec<0,0,longest_rec)))
# Split data into training and testing sets
combo_trainIndex <- createDataPartition(combo_data$longest_rec, p = 0.8, list = FALSE)
combo_trainData <- combo_data[combo_trainIndex, ]
combo_testData <- combo_data[-combo_trainIndex, ]

# Separate predictors and response
combo_x_train <- combo_trainData %>% select(-longest_rec)
combo_y_train <- combo_trainData$longest_rec
combo_x_test <- combo_testData %>% select(-longest_rec)
combo_y_test <- combo_testData$longest_rec

# Feature Scaling
combo_preProc <- preProcess(combo_x_train, method = c("center", "scale"))
combo_x_train <- predict(combo_preProc, combo_x_train)
combo_x_test <- predict(combo_preProc, combo_x_test)

# Model Training and Cross-Validation

# Linear Regression
combo_lm_model <- train(longest_rec ~ ., data = combo_trainData, method = "lm",
                           trControl = trainControl(method = "cv", number = 5))

# Random Forest
combo_rf_grid <- expand.grid(
  mtry = seq(2, ncol(combo_x_train), by = 2)  # Test mtry values for the transformed model
)

combo_rf_model <- train(
  longest_rec ~ ., 
  data = combo_trainData, 
  method = "rf",
  trControl = trainControl(method = "cv", number = 5),  # 5-fold cross-validation
  tuneGrid = combo_rf_grid,
  ntree = 500
)

# Gradient Boosting (XGBoost)
combo_xgb_grid <- expand.grid(
  nrounds = c(100, 200),
  max_depth = c(3, 5, 7),
  eta = c(0.01, 0.1),
  gamma = 0,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.8
)

combo_xgb_model <- train(longest_rec ~ ., data = combo_trainData, method = "xgbTree",
                            trControl = trainControl(method = "cv", number = 5),
                            tuneGrid = combo_xgb_grid)

# Evaluate Models
combo_models <- list(lm = combo_lm_model, rf = combo_rf_model, xgb = combo_xgb_model)
combo_results <- resamples(combo_models)

# Summarize results
summary(combo_results)

# Compare RMSE on the test set
combo_model_performance <- sapply(combo_models, function(model) {
  predictions <- predict(model, newdata = combo_x_test)
  rmse(combo_y_test, predictions)
})

# Select the best model
combo_best_model_name <- names(combo_model_performance)[which.min(combo_model_performance)]
combo_best_model <- combo_models[[combo_best_model_name]]

cat("Best model:", combo_best_model_name, "with RMSE:", min(combo_model_performance), "\n")

# Feature Importance (for Random Forest or XGBoost)
if (combo_best_model_name %in% c("rf", "xgb")) {
  varImp(combo_best_model) %>% plot()
}

combo_final_predictions <- as.data.frame(predict(combo_best_model, newdata = combo_x_test))
combo_final_check<- as.data.frame(predict(combo_best_model, newdata = combo_x_test, interval = "prediction", level = 0.5))

combo_residuals <- combo_y_test - combo_final_predictions[[1]]
#PCA----
# PCA Preprocessing
combo_pca_preProc <- preProcess(combo_x_train, method = "pca", pcaComp = 10) # Keep top 10 principal components
combo_x_train_pca <- predict(combo_pca_preProc, combo_x_train)
combo_x_test_pca <- predict(combo_pca_preProc, combo_x_test)

# Random Forest with PCA
combo_pca_rf_model <- train(
  longest_rec ~ ., 
  data = data.frame(longest_rec = combo_y_train, combo_x_train_pca),
  method = "rf",
  trControl = trainControl(method = "cv", number = 5),
  tuneGrid = expand.grid(mtry = seq(2, ncol(combo_x_train_pca), by = 2)),
  ntree = 500
)

# SVM with PCA
combo_pca_svm_model <- train(
  longest_rec ~ ., 
  data = data.frame(longest_rec = combo_y_train, combo_x_train_pca),
  method = "svmRadial",
  trControl = trainControl(method = "cv", number = 5),
  tuneGrid = expand.grid(
    sigma = 0.01, # SVM RBF kernel hyperparameter
    C = seq(1, 10, by = 1) # Regularization strength
  )
)

# Gradient Boosting (XGBoost) with PCA
combo_pca_xgb_model <- train(
  longest_rec ~ ., 
  data = data.frame(longest_rec = combo_y_train, combo_x_train_pca),
  method = "xgbTree",
  trControl = trainControl(method = "cv", number = 5),
  tuneGrid = expand.grid(
    nrounds = c(100, 200),
    max_depth = c(3, 5, 7),
    eta = c(0.01, 0.1),
    gamma = 0,
    colsample_bytree = 0.8,
    min_child_weight = 1,
    subsample = 0.8
  )
)

# KNN with PCA
combo_pca_knn_model <- train(
  longest_rec ~ ., 
  data = data.frame(longest_rec = combo_y_train, combo_x_train_pca),
  method = "knn",
  trControl = trainControl(method = "cv", number = 5),
  tuneGrid = expand.grid(k = seq(1, 20, by = 2)) # Number of neighbors
)

# Combine all PCA models for comparison
combo_pca_models <- list(
  rf = combo_pca_rf_model,
  svm = combo_pca_svm_model,
  xgb = combo_pca_xgb_model,
  knn = combo_pca_knn_model
)

# Evaluate PCA Models
combo_pca_results <- resamples(combo_pca_models)
summary(combo_pca_results)

# Compare RMSE on the test set
combo_pca_model_performance <- sapply(combo_pca_models, function(model) {
  predictions <- predict(model, newdata = combo_x_test_pca)
  rmse(combo_y_test, predictions)
})

# Select the Best PCA Model
combo_best_pca_model_name <- names(combo_pca_model_performance)[which.min(combo_pca_model_performance)]
combo_best_pca_model <- combo_pca_models[[combo_best_pca_model_name]]

cat("Best PCA model:", combo_best_pca_model_name, "with RMSE:", min(combo_pca_model_performance), "\n")

# Feature Importance (if applicable)
if (combo_best_pca_model_name %in% c("rf", "xgb")) {
  varImp(combo_best_pca_model) %>% plot()
}

# Final Predictions with Best PCA Model
combo_final_predictions_pca <- as.data.frame(predict(combo_best_pca_model, newdata = combo_x_test_pca))

# Residual Analysis
combo_residuals_pca <- combo_y_test - combo_final_predictions_pca[[1]]

# Output RMSE Comparison
cat("PCA Model RMSEs:\n")
print(combo_pca_model_performance)

#Exponential -----
# Integrate exponential decay into feature engineering
library(zoo)

# Create a range of decay factors to tune
decay_factors <- seq(0.7, 1, by = 0.05)  # Adjust range as needed

tune_results <- list()

for (decay_factor in decay_factors) {
  # Exponentially decay weights for key metrics
  season_weighted_exp <- game_wr %>%
    group_by(receiver_player_id, posteam) %>%
    arrange(game_id, week) %>% # Ensure the data is sorted by week
    mutate(
      # Exponential decay for metrics
      exp_weighted_targets_per_game = zoo::rollapply(targets, width = 10, 
                                                     FUN = function(x) sum(x * decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
      exp_weighted_air_yards_per_game = zoo::rollapply(total_air_yards, width = 10, 
                                                       FUN = function(x) sum(x * decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
      exp_weighted_yards_game_per_game = zoo::rollapply(total_yards_gained, width = 10, 
                                                        FUN = function(x) sum(x * decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
      exp_weighted_longest_rec_per_game = zoo::rollapply(longest_rec, width = 10, 
                                                         FUN = function(x) sum(x * decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
      
      # Cumulative metrics (no change from prior implementation)
      cum_targets = cumsum(targets),
      cum_target_rate = cumsum(catchable_targ) / cumsum(targets),
      cum_contested_rate = cumsum(contested_ball) / cumsum(targets),
      cum_created_reception_rate = cumsum(created_rec) / cumsum(targets)
    ) %>%
    ungroup() %>% 
    group_by(receiver_player_id) %>%
    mutate(next_week = lead(week)) %>% 
    filter(next_week > 1) %>% 
    filter(!is.na(week)) %>% 
    select(
      week = next_week, receiver_player_id, year, posteam,
      starts_with("exp_weighted"), starts_with("cum")
    )
  
  # Merge with additional data sources
  longest_rec_joined_pexp_data <- longest_rec_data %>% 
    left_join(season_weighted_pexp, by = c("receiver_player_id","week","year","posteam")) %>% 
    left_join(off_scout, by = c("week","posteam","year")) %>%
    left_join(def_scout, by = c("week","defteam","year")) %>% 
    ungroup() %>% 
    select(-posteam,-defteam) %>% 
    filter(week > 1)
  
  # Data preparation for modeling
  exp_data <- longest_rec_joined_pexp_data %>%
    drop_na() %>%  # Remove rows with missing values
    mutate_if(is.character, as.factor) %>%  # Convert character columns to factors
    select(-receiver_player_id, -game_id, -name) %>%  # Drop non-predictive or irrelevant columns
    mutate(longest_rec = sqrt(ifelse(longest_rec < 0, 0, longest_rec)))
  
  # Split data into training and testing sets
  exp_trainIndex <- createDataPartition(exp_data$longest_rec, p = 0.8, list = FALSE)
  exp_trainData <- exp_data[exp_trainIndex, ]
  exp_testData <- exp_data[-exp_trainIndex, ]
  
  # Separate predictors and response
  exp_x_train <- exp_trainData %>% select(-longest_rec)
  exp_y_train <- exp_trainData$longest_rec
  exp_x_test <- exp_testData %>% select(-longest_rec)
  exp_y_test <- exp_testData$longest_rec
  
  # Feature Scaling
  exp_preProc <- preProcess(exp_x_train, method = c("center", "scale"))
  exp_x_train <- predict(exp_preProc, exp_x_train)
  exp_x_test <- predict(exp_preProc, exp_x_test)
  
  # Model Training and Cross-Validation
  exp_rf_grid <- expand.grid(
    mtry = seq(2, ncol(exp_x_train), by = 2)  # Test mtry values for the transformed model
  )
  
  exp_rf_model <- train(
    longest_rec ~ ., 
    data = exp_trainData, 
    method = "rf",
    trControl = trainControl(method = "cv", number = 5),  # 5-fold cross-validation
    tuneGrid = exp_rf_grid,
    ntree = 500
  )
  
  # Evaluate Model Performance
  predictions <- predict(exp_rf_model, newdata = exp_x_test)
  rmse_value <- rmse(exp_y_test, predictions)
  
  # Store results
  tune_results[[as.character(decay_factor)]] <- list(
    decay_factor = decay_factor,
    rmse = rmse_value
  )
}

# Find the best decay factor
best_tune <- do.call(rbind, tune_results) %>% as.data.frame()
best_decay_factor <- best_tune[which.min(best_tune$rmse), "decay_factor"][[1]] #0.85

# cat("Best Decay Factor:", best_decay_factor, "with RMSE:", min(best_tune$rmse), "\n")

# Re-run the process with the best decay factor
# Exponentially decay weights for key metrics with the best decay factor
season_weighted_pexp <- game_wr %>%
  group_by(receiver_player_id, posteam) %>%
  arrange(game_id, week) %>% # Ensure the data is sorted by week
  mutate(
    exp_weighted_targets_per_game = zoo::rollapply(targets, width = 10, 
                                                   FUN = function(x) sum(x * best_decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
    exp_weighted_air_yards_per_game = zoo::rollapply(total_air_yards, width = 10, 
                                                     FUN = function(x) sum(x * best_decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
    exp_weighted_yards_game_per_game = zoo::rollapply(total_yards_gained, width = 10, 
                                                      FUN = function(x) sum(x * best_decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
    exp_weighted_longest_rec_per_game = zoo::rollapply(longest_rec, width = 10, 
                                                       FUN = function(x) sum(x * best_decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
    
    cum_targets = cumsum(targets),
    cum_target_rate = cumsum(catchable_targ) / cumsum(targets),
    cum_contested_rate = cumsum(contested_ball) / cumsum(targets),
    cum_created_reception_rate = cumsum(created_rec) / cumsum(targets)
  ) %>%
  ungroup() %>% 
  group_by(receiver_player_id) %>%
  mutate(next_week = lead(week)) %>% 
  filter(next_week > 1) %>% 
  filter(!is.na(week)) %>% 
  select(
    week = next_week, receiver_player_id, year, posteam,
    starts_with("exp_weighted"), starts_with("cum")
  )

# Merge with additional data sources
longest_rec_joined_pexp_data <- longest_rec_data %>% 
  left_join(season_weighted_pexp, by = c("receiver_player_id","week","year","posteam")) %>% 
  left_join(off_scout, by = c("week","posteam","year")) %>%
  left_join(def_scout, by = c("week","defteam","year")) %>% 
  ungroup() %>% 
  select(-posteam,-defteam) %>% 
  filter(week > 1)

# Data preparation for modeling
exp_data <- longest_rec_joined_pexp_data %>%
  drop_na() %>%  # Remove rows with missing values
  mutate_if(is.character, as.factor) %>%  # Convert character columns to factors
  select(-receiver_player_id, -game_id, -name) %>%  # Drop non-predictive or irrelevant columns
  mutate(longest_rec = sqrt(ifelse(longest_rec < 0, 0, longest_rec)))

# Split data into training and testing sets
exp_trainIndex <- createDataPartition(exp_data$longest_rec, p = 0.8, list = FALSE)
exp_trainData <- exp_data[exp_trainIndex, ]
exp_testData <- exp_data[-exp_trainIndex, ]

# Separate predictors and response
exp_x_train <- exp_trainData %>% select(-longest_rec)
exp_y_train <- exp_trainData$longest_rec
exp_x_test <- exp_testData %>% select(-longest_rec)
exp_y_test <- exp_testData$longest_rec

# Feature Scaling
exp_preProc <- preProcess(exp_x_train, method = c("center", "scale"))
exp_x_train <- predict(exp_preProc, exp_x_train)
exp_x_test <- predict(exp_preProc, exp_x_test)

# Re-train all models with best decay factor

# Linear Regression
exp_lm_model <- train(longest_rec ~ ., data = exp_trainData, method = "lm",
                      trControl = trainControl(method = "cv", number = 5))

# Random Forest
exp_rf_grid <- expand.grid(
  mtry = seq(2, ncol(exp_x_train), by = 2)  # Test mtry values for the transformed model
)

exp_rf_model <- train(
  longest_rec ~ ., 
  data = exp_trainData, 
  method = "rf",
  trControl = trainControl(method = "cv", number = 5),  # 5-fold cross-validation
  tuneGrid = exp_rf_grid,
  ntree = 500
)

# Gradient Boosting (XGBoost)
exp_xgb_grid <- expand.grid(
  nrounds = c(100, 200),
  max_depth = c(3, 5, 7),
  eta = c(0.01, 0.1),
  gamma = 0,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.8
)

exp_xgb_model <- train(longest_rec ~ ., data = exp_trainData, method = "xgbTree",
                       trControl = trainControl(method = "cv", number = 5),
                       tuneGrid = exp_xgb_grid)

# Evaluate Models
exp_models <- list(lm = exp_lm_model, rf = exp_rf_model, xgb = exp_xgb_model)
exp_results <- resamples(exp_models)

# Summarize results
summary(exp_results)

# Compare RMSE on the test set
exp_model_performance <- sapply(exp_models, function(model) {
  predictions <- predict(model, newdata = exp_x_test)
  rmse(exp_y_test, predictions)
})

# Select the best model
exp_best_model_name <- names(exp_model_performance)[which.min(exp_model_performance)]
exp_best_model <- exp_models[[exp_best_model_name]]

cat("Best model:", exp_best_model_name, "with RMSE:", min(exp_model_performance), "
")
exp_final_predictions <- as.data.frame(predict(exp_best_model, newdata = exp_x_test))




#Exponential No Transform -----
# Integrate nt_exponential decay into feature engineering
library(zoo)

# Create a range of decay factors to tune
decay_factors <- seq(0.7, 1, by = 0.05)  # Adjust range as needed

tune_results <- list()

for (decay_factor in decay_factors) {
  # nt_exponentially decay weights for key metrics
  season_weighted_nt_exp <- game_wr %>%
    group_by(receiver_player_id, posteam) %>%
    arrange(game_id, week) %>% # Ensure the data is sorted by week
    mutate(
      # nt_exponential decay for metrics
      nt_exp_weighted_targets_per_game = zoo::rollapply(targets, width = 10, 
                                                     FUN = function(x) sum(x * decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
      nt_exp_weighted_air_yards_per_game = zoo::rollapply(total_air_yards, width = 10, 
                                                       FUN = function(x) sum(x * decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
      nt_exp_weighted_yards_game_per_game = zoo::rollapply(total_yards_gained, width = 10, 
                                                        FUN = function(x) sum(x * decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
      nt_exp_weighted_longest_rec_per_game = zoo::rollapply(longest_rec, width = 10, 
                                                         FUN = function(x) sum(x * decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
      
      # Cumulative metrics (no change from prior implementation)
      cum_targets = cumsum(targets),
      cum_target_rate = cumsum(catchable_targ) / cumsum(targets),
      cum_contested_rate = cumsum(contested_ball) / cumsum(targets),
      cum_created_reception_rate = cumsum(created_rec) / cumsum(targets)
    ) %>%
    ungroup() %>% 
    group_by(receiver_player_id) %>%
    mutate(next_week = lead(week)) %>% 
    filter(next_week > 1) %>% 
    filter(!is.na(week)) %>% 
    select(
      week = next_week, receiver_player_id, year, posteam,
      starts_with("nt_exp_weighted"), starts_with("cum")
    )
  
  # Merge with additional data sources
  longest_rec_joined_pnt_exp_data <- longest_rec_data %>% 
    left_join(season_weighted_nt_exp, by = c("receiver_player_id","week","year","posteam")) %>% 
    left_join(off_scout, by = c("week","posteam","year")) %>%
    left_join(def_scout, by = c("week","defteam","year")) %>% 
    ungroup() %>% 
    select(-posteam,-defteam) %>% 
    filter(week > 1)
  
  # Data preparation for modeling
  nt_exp_data <- longest_rec_joined_pnt_exp_data %>%
    drop_na() %>%  # Remove rows with missing values
    mutate_if(is.character, as.factor) %>%  # Convert character columns to factors
    select(-receiver_player_id, -game_id, -name)
  
  # Split data into training and testing sets
  nt_exp_trainIndex <- createDataPartition(nt_exp_data$longest_rec, p = 0.8, list = FALSE)
  nt_exp_trainData <- nt_exp_data[nt_exp_trainIndex, ]
  nt_exp_testData <- nt_exp_data[-nt_exp_trainIndex, ]
  
  # Separate predictors and response
  nt_exp_x_train <- nt_exp_trainData %>% select(-longest_rec)
  nt_exp_y_train <- nt_exp_trainData$longest_rec
  nt_exp_x_test <- nt_exp_testData %>% select(-longest_rec)
  nt_exp_y_test <- nt_exp_testData$longest_rec
  
  # Feature Scaling
  nt_exp_preProc <- preProcess(nt_exp_x_train, method = c("center", "scale"))
  nt_exp_x_train <- predict(nt_exp_preProc, nt_exp_x_train)
  nt_exp_x_test <- predict(nt_exp_preProc, nt_exp_x_test)
  
  # Model Training and Cross-Validation
  nt_exp_rf_grid <- expand.grid(
    mtry = seq(2, ncol(nt_exp_x_train), by = 2)  # Test mtry values for the transformed model
  )
  
  nt_exp_rf_model <- train(
    longest_rec ~ ., 
    data = nt_exp_trainData, 
    method = "rf",
    trControl = trainControl(method = "cv", number = 5),  # 5-fold cross-validation
    tuneGrid = nt_exp_rf_grid,
    ntree = 500
  )
  
  # Evaluate Model Performance
  predictions <- predict(nt_exp_rf_model, newdata = nt_exp_x_test)
  rmse_value <- rmse(nt_exp_y_test, predictions)
  
  # Store results
  tune_results[[as.character(decay_factor)]] <- list(
    decay_factor = decay_factor,
    rmse = rmse_value
  )
}

# Find the best decay factor
best_tune <- do.call(rbind, tune_results) %>% as.data.frame()
best_decay_factor <- best_tune[which.min(best_tune$rmse), "decay_factor"][[1]] 

# cat("Best Decay Factor:", best_decay_factor, "with RMSE:", min(best_tune$rmse), "\n")

# Re-run the process with the best decay factor
# nt_exponentially decay weights for key metrics with the best decay factor
season_weighted_pnt_exp <- game_wr %>%
  group_by(receiver_player_id, posteam) %>%
  arrange(game_id, week) %>% # Ensure the data is sorted by week
  mutate(
    nt_exp_weighted_targets_per_game = zoo::rollapply(targets, width = 10, 
                                                   FUN = function(x) sum(x * best_decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
    nt_exp_weighted_air_yards_per_game = zoo::rollapply(total_air_yards, width = 10, 
                                                     FUN = function(x) sum(x * best_decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
    nt_exp_weighted_yards_game_per_game = zoo::rollapply(total_yards_gained, width = 10, 
                                                      FUN = function(x) sum(x * best_decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
    nt_exp_weighted_longest_rec_per_game = zoo::rollapply(longest_rec, width = 10, 
                                                       FUN = function(x) sum(x * best_decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
    
    cum_targets = cumsum(targets),
    cum_target_rate = cumsum(catchable_targ) / cumsum(targets),
    cum_contested_rate = cumsum(contested_ball) / cumsum(targets),
    cum_created_reception_rate = cumsum(created_rec) / cumsum(targets)
  ) %>%
  ungroup() %>% 
  group_by(receiver_player_id) %>%
  mutate(next_week = lead(week)) %>% 
  filter(next_week > 1) %>% 
  filter(!is.na(week)) %>% 
  select(
    week = next_week, receiver_player_id, year, posteam,
    starts_with("nt_exp_weighted"), starts_with("cum")
  )

# Merge with additional data sources
longest_rec_joined_pnt_exp_data <- longest_rec_data %>% 
  left_join(season_weighted_pnt_exp, by = c("receiver_player_id","week","year","posteam")) %>% 
  left_join(off_scout, by = c("week","posteam","year")) %>%
  left_join(def_scout, by = c("week","defteam","year")) %>% 
  ungroup() %>% 
  select(-posteam,-defteam) %>% 
  filter(week > 1)

# Data preparation for modeling
nt_exp_data <- longest_rec_joined_pnt_exp_data %>%
  drop_na() %>%  # Remove rows with missing values
  mutate_if(is.character, as.factor) %>%  # Convert character columns to factors
  select(-receiver_player_id, -game_id, -name) 

# Split data into training and testing sets
nt_exp_trainIndex <- createDataPartition(nt_exp_data$longest_rec, p = 0.8, list = FALSE)
nt_exp_trainData <- nt_exp_data[nt_exp_trainIndex, ]
nt_exp_testData <- nt_exp_data[-nt_exp_trainIndex, ]

# Separate predictors and response
nt_exp_x_train <- nt_exp_trainData %>% select(-longest_rec)
nt_exp_y_train <- nt_exp_trainData$longest_rec
nt_exp_x_test <- nt_exp_testData %>% select(-longest_rec)
nt_exp_y_test <- nt_exp_testData$longest_rec

# Feature Scaling
nt_exp_preProc <- preProcess(nt_exp_x_train, method = c("center", "scale"))
nt_exp_x_train <- predict(nt_exp_preProc, nt_exp_x_train)
nt_exp_x_test <- predict(nt_exp_preProc, nt_exp_x_test)

# Re-train all models with best decay factor

# Linear Regression
nt_exp_lm_model <- train(longest_rec ~ ., data = nt_exp_trainData, method = "lm",
                      trControl = trainControl(method = "cv", number = 5))

# Random Forest
nt_exp_rf_grid <- expand.grid(
  mtry = seq(2, ncol(nt_exp_x_train), by = 2)  # Test mtry values for the transformed model
)

nt_exp_rf_model <- train(
  longest_rec ~ ., 
  data = nt_exp_trainData, 
  method = "rf",
  trControl = trainControl(method = "cv", number = 5),  # 5-fold cross-validation
  tuneGrid = nt_exp_rf_grid,
  ntree = 500
)

# Gradient Boosting (XGBoost)
nt_exp_xgb_grid <- expand.grid(
  nrounds = c(100, 200),
  max_depth = c(3, 5, 7),
  eta = c(0.01, 0.1),
  gamma = 0,
  colsample_bytree = 0.8,
  min_child_weight = 1,
  subsample = 0.8
)

nt_exp_xgb_model <- train(longest_rec ~ ., data = nt_exp_trainData, method = "xgbTree",
                       trControl = trainControl(method = "cv", number = 5),
                       tuneGrid = nt_exp_xgb_grid)

# Evaluate Models
nt_exp_models <- list(lm = nt_exp_lm_model, rf = nt_exp_rf_model, xgb = nt_exp_xgb_model)
nt_exp_results <- resamples(nt_exp_models)

# Summarize results
summary(nt_exp_results)

# Compare RMSE on the test set
nt_exp_model_performance <- sapply(nt_exp_models, function(model) {
  predictions <- predict(model, newdata = nt_exp_x_test)
  rmse(nt_exp_y_test, predictions)
})

# Select the best model
nt_exp_best_model_name <- names(nt_exp_model_performance)[which.min(nt_exp_model_performance)]
nt_exp_best_model <- nt_exp_models[[nt_exp_best_model_name]]

cat("Best model:", nt_exp_best_model_name, "with RMSE:", min(nt_exp_model_performance), "
")
nt_exp_final_predictions <- as.data.frame(predict(nt_exp_best_model, newdata = nt_exp_x_test))

#No Transform Roll Exp PCA----
decay_factor <- 0.7 #Could be useful to further tune
season_nt_pexp <- game_wr %>%
  group_by(receiver_player_id, posteam) %>%
  arrange(game_id, week) %>% # Ensure the data is sorted by week
  mutate(
    exp_weighted_targets_per_game = zoo::rollapply(targets, width = 10, 
                                                   FUN = function(x) sum(x * decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
    exp_weighted_air_yards_per_game = zoo::rollapply(total_air_yards, width = 10, 
                                                     FUN = function(x) sum(x * decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
    exp_weighted_yards_game_per_game = zoo::rollapply(total_yards_gained, width = 10, 
                                                      FUN = function(x) sum(x * decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
    exp_weighted_longest_rec_per_game = zoo::rollapply(longest_rec, width = 10, 
                                                       FUN = function(x) sum(x * decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
    max_week = max(week), # Find the most recent week
    week_weight = (week - min(week) + 1) / (max_week - min(week) + 1), # Scale weights to give higher value to recent weeks
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
    # Weighted averages for per-game metrics
    weighted_targets_per_game = cumsum(targets * week_weight) / cumsum(week_weight),
    weighted_air_yards_per_game = cumsum(total_air_yards * week_weight) / cumsum(week_weight),
    weighted_yards_game_per_game = cumsum(total_yards_gained * week_weight) / cumsum(week_weight),
    weighted_yac_per_game = cumsum(YAC * week_weight) / cumsum(week_weight),
    weighted_epa_per_game = cumsum(total_epa * week_weight) / cumsum(week_weight),
    weighted_success_rate_per_game = cumsum(success_rate * week_weight) / cumsum(week_weight),
    weighted_longest_rec_per_game = cumsum(longest_rec * week_weight) / cumsum(week_weight),
    weighted_pct_air_yards_per_game = cumsum(pct_air_yards * week_weight) / cumsum(week_weight),
    weighted_target_share_per_game = cumsum(target_share * week_weight) / cumsum(week_weight),
    weighted_adot_per_game = cumsum(adot * week_weight) / cumsum(week_weight),
    
    # Cumulative metrics
    cum_targets = cumsum(targets),
    cum_target_rate = cumsum(catchable_targ) / cumsum(targets),
    cum_contested_rate = cumsum(contested_ball) / cumsum(targets),
    cum_created_reception_rate = cumsum(created_rec) / cumsum(targets)
  ) %>%
  ungroup() %>% 
  group_by(receiver_player_id) %>%
  mutate(next_week = lead(week)) %>% 
  filter(next_week > 1) %>% 
  filter(!is.na(week)) %>% 
  select(
    week = next_week, receiver_player_id, year, posteam,
    ends_with("per_game"), starts_with("cum")
  )
longest_rec_joined_nt_pexp_data <- longest_rec_data %>% 
  left_join(season_nt_pexp, by = c("receiver_player_id","week","year","posteam")) %>% 
  left_join(off_scout, by = c("week","posteam","year")) %>%
  left_join(def_scout, by = c("week","defteam","year")) %>% 
  ungroup() %>% 
  select(-posteam,-defteam) %>% 
  filter(week>1)
nt_pexp_data <- longest_rec_joined_nt_pexp_data %>%
  drop_na() %>%  # Remove rows with missing values
  mutate_if(is.character, as.factor) %>%  # Convert character columns to factors
  select(-receiver_player_id, -game_id, -name)
# Split data into training and testing sets
nt_pexp_trainIndex <- createDataPartition(nt_pexp_data$longest_rec, p = 0.8, list = FALSE)
nt_pexp_trainData <- nt_pexp_data[nt_pexp_trainIndex, ]
nt_pexp_testData <- nt_pexp_data[-nt_pexp_trainIndex, ]

# Separate predictors and response
nt_pexp_x_train <- nt_pexp_trainData %>% select(-longest_rec)
nt_pexp_y_train <- nt_pexp_trainData$longest_rec
nt_pexp_x_test <- nt_pexp_testData %>% select(-longest_rec)
nt_pexp_y_test <- nt_pexp_testData$longest_rec


nt_pexp_pca_preProc <- preProcess(nt_pexp_x_train, method = "pca", pcaComp = 10) # Keep top 10 principal components
nt_pexp_x_train_pca <- predict(nt_pexp_pca_preProc, nt_pexp_x_train)
nt_pexp_x_test_pca <- predict(nt_pexp_pca_preProc, nt_pexp_x_test)

# Random Forest with PCA
nt_pexp_pca_rf_model <- train(
  longest_rec ~ ., 
  data = data.frame(longest_rec = nt_pexp_y_train, nt_pexp_x_train_pca),
  method = "rf",
  trControl = trainControl(method = "cv", number = 5),
  tuneGrid = expand.grid(mtry = seq(2, ncol(nt_pexp_x_train_pca), by = 2)),
  ntree = 500
)

# SVM with PCA
nt_pexp_pca_svm_model <- train(
  longest_rec ~ ., 
  data = data.frame(longest_rec = nt_pexp_y_train, nt_pexp_x_train_pca),
  method = "svmRadial",
  trControl = trainControl(method = "cv", number = 5),
  tuneGrid = expand.grid(
    sigma = 0.01, # SVM RBF kernel hyperparameter
    C = seq(1, 10, by = 1) # Regularization strength
  )
)

# Gradient Boosting (XGBoost) with PCA
nt_pexp_pca_xgb_model <- train(
  longest_rec ~ ., 
  data = data.frame(longest_rec = nt_pexp_y_train, nt_pexp_x_train_pca),
  method = "xgbTree",
  trControl = trainControl(method = "cv", number = 5),
  tuneGrid = expand.grid(
    nrounds = c(100, 200, 300),
    max_depth = c(3, 5, 7),
    eta = c(0.01, 0.1, 0.2),
    gamma = c(0, 1, 5),
    colsample_bytree = c(0.6, 0.8),
    min_child_weight = c(1, 5, 10),
    subsample = c(0.7, 0.8, 0.9)
  )
)

library(xgboost)

# Prepare Data
dtrain <- xgb.DMatrix(data = as.matrix(nt_pexp_x_train_pca), label = nt_pexp_y_train)
dtest <- xgb.DMatrix(data = as.matrix(nt_pexp_x_test_pca), label = nt_pexp_y_test)

# Grid of parameters to tune
param_grid <- expand.grid(
  max_depth = c(3, 5, 7),
  eta = c(0.01, 0.1, 0.2),
  gamma = c(0, 1, 5),
  min_child_weight = c(1, 5, 10),
  subsample = c(0.7, 0.8),
  colsample_bytree = c(0.6, 0.8)
)

# Storage for results
tuning_results <- data.frame()

# Loop over parameter combinations
for (i in 1:nrow(param_grid)) {
  params <- list(
    objective = "reg:pseudohubererror",  # Use the correct objective function
    max_depth = param_grid$max_depth[i],
    eta = param_grid$eta[i],
    gamma = param_grid$gamma[i],
    min_child_weight = param_grid$min_child_weight[i],
    subsample = param_grid$subsample[i],
    colsample_bytree = param_grid$colsample_bytree[i]
  )
  
  # Cross-validation
  cv_results <- xgb.cv(
    params = params,
    data = dtrain,
    nrounds = 200,
    nfold = 5,
    metrics = list("rmse"),
    early_stopping_rounds = 10,
    verbose = FALSE
  )
  
  # Save results
  tuning_results <- rbind(
    tuning_results,
    cbind(param_grid[i, ], best_rmse = min(cv_results$evaluation_log$test_rmse_mean))
  )
}

# Find the best parameters
best_params <- tuning_results[which.min(tuning_results$best_rmse), ]
print(best_params)

# Convert best_params to a list for final training
final_params <- list(
  objective = "reg:pseudohubererror",  # Use Huber-like loss
  max_depth = best_params$max_depth,
  eta = best_params$eta,
  gamma = best_params$gamma,
  min_child_weight = best_params$min_child_weight,
  subsample = best_params$subsample,
  colsample_bytree = best_params$colsample_bytree
)

# Train the final model
final_xgb_model <- xgb.train(
  params = final_params,
  data = dtrain,
  nrounds = 200,  # Use the same maximum rounds as during cross-validation
  watchlist = list(train = dtrain),
  verbose = TRUE
)

# Make predictions on the test set
test_predictions <- predict(final_xgb_model, dtest)

# Calculate RMSE
final_rmse <- sqrt(mean((test_predictions - nt_pexp_y_test)^2))
cat("Final RMSE on Test Set:", final_rmse, "\n")

# Plot residuals
residuals <- nt_pexp_y_test - test_predictions

# Residual plot
plot(test_predictions, residuals,
     main = "Residuals vs Predictions (Final Model)",
     xlab = "Predicted Values",
     ylab = "Residuals",
     pch = 20)
abline(h = 0, col = "red")



# KNN with PCA
nt_pexp_pca_knn_model <- train(
  longest_rec ~ ., 
  data = data.frame(longest_rec = nt_pexp_y_train, nt_pexp_x_train_pca),
  method = "knn",
  trControl = trainControl(method = "cv", number = 5),
  tuneGrid = expand.grid(k = seq(1, 20, by = 2)) # Number of neighbors
)

# Combine all PCA models for comparison
nt_pexp_pca_models <- list(
  rf = nt_pexp_pca_rf_model,
  svm = nt_pexp_pca_svm_model,
  xgb = nt_pexp_pca_xgb_model,
  knn = nt_pexp_pca_knn_model
)

# Evaluate PCA Models
nt_pexp_pca_results <- resamples(nt_pexp_pca_models)
summary(nt_pexp_pca_results)

# Compare RMSE on the test set
nt_pexp_pca_model_performance <- sapply(nt_pexp_pca_models, function(model) {
  predictions <- predict(model, newdata = nt_pexp_x_test_pca)
  rmse(nt_pexp_y_test, predictions)
})

# Select the Best PCA Model
nt_pexp_best_pca_model_name <- names(nt_pexp_pca_model_performance)[which.min(nt_pexp_pca_model_performance)]
nt_pexp_best_pca_model <- nt_pexp_pca_models[[nt_pexp_best_pca_model_name]]

cat("Best PCA model:", nt_pexp_best_pca_model_name, "with RMSE:", min(nt_pexp_pca_model_performance), "\n")

# Feature Importance (if applicable)
if (nt_pexp_best_pca_model_name %in% c("rf", "xgb")) {
  varImp(nt_pexp_best_pca_model) %>% plot()
}

# Final Predictions with Best PCA Model
nt_pexp_final_predictions_pca <- as.data.frame(predict(nt_pexp_best_pca_model, newdata = nt_pexp_x_test_pca))

# Residual Analysis
nt_pexp_residuals_pca <- nt_pexp_y_test - nt_pexp_final_predictions_pca[[1]]

# Output RMSE Comparison
cat("PCA Model RMSEs:\n")
print(nt_pexp_pca_model_performance)

#Roll Exp PCA----
decay_factor <- 0.85
season_pexp <- game_wr %>%
  group_by(receiver_player_id, posteam) %>%
  arrange(game_id, week) %>% # Ensure the data is sorted by week
  mutate(
    exp_weighted_targets_per_game = zoo::rollapply(targets, width = 10, 
                                                   FUN = function(x) sum(x * decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
    exp_weighted_air_yards_per_game = zoo::rollapply(total_air_yards, width = 10, 
                                                     FUN = function(x) sum(x * decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
    exp_weighted_yards_game_per_game = zoo::rollapply(total_yards_gained, width = 10, 
                                                      FUN = function(x) sum(x * decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
    exp_weighted_longest_rec_per_game = zoo::rollapply(longest_rec, width = 10, 
                                                       FUN = function(x) sum(x * decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
    max_week = max(week), # Find the most recent week
    week_weight = (week - min(week) + 1) / (max_week - min(week) + 1), # Scale weights to give higher value to recent weeks
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
    # Weighted averages for per-game metrics
    weighted_targets_per_game = cumsum(targets * week_weight) / cumsum(week_weight),
    weighted_air_yards_per_game = cumsum(total_air_yards * week_weight) / cumsum(week_weight),
    weighted_yards_game_per_game = cumsum(total_yards_gained * week_weight) / cumsum(week_weight),
    weighted_yac_per_game = cumsum(YAC * week_weight) / cumsum(week_weight),
    weighted_epa_per_game = cumsum(total_epa * week_weight) / cumsum(week_weight),
    weighted_success_rate_per_game = cumsum(success_rate * week_weight) / cumsum(week_weight),
    weighted_longest_rec_per_game = cumsum(longest_rec * week_weight) / cumsum(week_weight),
    weighted_pct_air_yards_per_game = cumsum(pct_air_yards * week_weight) / cumsum(week_weight),
    weighted_target_share_per_game = cumsum(target_share * week_weight) / cumsum(week_weight),
    weighted_adot_per_game = cumsum(adot * week_weight) / cumsum(week_weight),
    
    # Cumulative metrics
    cum_targets = cumsum(targets),
    cum_target_rate = cumsum(catchable_targ) / cumsum(targets),
    cum_contested_rate = cumsum(contested_ball) / cumsum(targets),
    cum_created_reception_rate = cumsum(created_rec) / cumsum(targets)
  ) %>%
  ungroup() %>% 
  group_by(receiver_player_id) %>%
  mutate(next_week = lead(week)) %>% 
  filter(next_week > 1) %>% 
  filter(!is.na(week)) %>% 
  select(
    week = next_week, receiver_player_id, year, posteam,
    ends_with("per_game"), starts_with("cum")
  )
longest_rec_joined_pexp_data <- longest_rec_data %>% 
  left_join(season_pexp, by = c("receiver_player_id","week","year","posteam")) %>% 
  left_join(off_scout, by = c("week","posteam","year")) %>%
  left_join(def_scout, by = c("week","defteam","year")) %>% 
  ungroup() %>% 
  select(-posteam,-defteam) %>% 
  filter(week>1)
pexp_data <- longest_rec_joined_pexp_data %>%
  drop_na() %>%  # Remove rows with missing values
  mutate_if(is.character, as.factor) %>%  # Convert character columns to factors
  select(-receiver_player_id, -game_id, -name) %>%  # Drop non-predictive or irrelevant columns
  mutate(longest_rec = sqrt(ifelse(longest_rec<0,0,longest_rec)))
# Split data into training and testing sets
pexp_trainIndex <- createDataPartition(pexp_data$longest_rec, p = 0.8, list = FALSE)
pexp_trainData <- pexp_data[pexp_trainIndex, ]
pexp_testData <- pexp_data[-pexp_trainIndex, ]

# Separate predictors and response
pexp_x_train <- pexp_trainData %>% select(-longest_rec)
pexp_y_train <- pexp_trainData$longest_rec
pexp_x_test <- pexp_testData %>% select(-longest_rec)
pexp_y_test <- pexp_testData$longest_rec


pexp_pca_preProc <- preProcess(pexp_x_train, method = "pca", pcaComp = 10) # Keep top 10 principal components
pexp_x_train_pca <- predict(pexp_pca_preProc, pexp_x_train)
pexp_x_test_pca <- predict(pexp_pca_preProc, pexp_x_test)

# Random Forest with PCA
pexp_pca_rf_model <- train(
  longest_rec ~ ., 
  data = data.frame(longest_rec = pexp_y_train, pexp_x_train_pca),
  method = "rf",
  trControl = trainControl(method = "cv", number = 5),
  tuneGrid = expand.grid(mtry = seq(2, ncol(pexp_x_train_pca), by = 2)),
  ntree = 500
)

# SVM with PCA
pexp_pca_svm_model <- train(
  longest_rec ~ ., 
  data = data.frame(longest_rec = pexp_y_train, pexp_x_train_pca),
  method = "svmRadial",
  trControl = trainControl(method = "cv", number = 5),
  tuneGrid = expand.grid(
    sigma = 0.01, # SVM RBF kernel hyperparameter
    C = seq(1, 10, by = 1) # Regularization strength
  )
)

# Gradient Boosting (XGBoost) with PCA
pexp_pca_xgb_model <- train(
  longest_rec ~ ., 
  data = data.frame(longest_rec = pexp_y_train, pexp_x_train_pca),
  method = "xgbTree",
  trControl = trainControl(method = "cv", number = 5),
  tuneGrid = expand.grid(
    nrounds = c(100, 200),
    max_depth = c(3, 5, 7),
    eta = c(0.01, 0.1),
    gamma = 0,
    colsample_bytree = 0.8,
    min_child_weight = 1,
    subsample = 0.8
  )
)

# KNN with PCA
pexp_pca_knn_model <- train(
  longest_rec ~ ., 
  data = data.frame(longest_rec = pexp_y_train, pexp_x_train_pca),
  method = "knn",
  trControl = trainControl(method = "cv", number = 5),
  tuneGrid = expand.grid(k = seq(1, 20, by = 2)) # Number of neighbors
)

# Combine all PCA models for comparison
pexp_pca_models <- list(
  rf = pexp_pca_rf_model,
  svm = pexp_pca_svm_model,
  xgb = pexp_pca_xgb_model,
  knn = pexp_pca_knn_model
)

# Evaluate PCA Models
pexp_pca_results <- resamples(pexp_pca_models)
summary(pexp_pca_results)

# Compare RMSE on the test set
pexp_pca_model_performance <- sapply(pexp_pca_models, function(model) {
  predictions <- predict(model, newdata = pexp_x_test_pca)
  rmse(pexp_y_test, predictions)
})

# Select the Best PCA Model
pexp_best_pca_model_name <- names(pexp_pca_model_performance)[which.min(pexp_pca_model_performance)]
pexp_best_pca_model <- pexp_pca_models[[pexp_best_pca_model_name]]

cat("Best PCA model:", pexp_best_pca_model_name, "with RMSE:", min(pexp_pca_model_performance), "\n")

# Feature Importance (if applicable)
if (pexp_best_pca_model_name %in% c("rf", "xgb")) {
  varImp(pexp_best_pca_model) %>% plot()
}

# Final Predictions with Best PCA Model
pexp_final_predictions_pca <- as.data.frame(predict(pexp_best_pca_model, newdata = pexp_x_test_pca))

# Residual Analysis
pexp_residuals_pca <- pexp_y_test - pexp_final_predictions_pca[[1]]

# Output RMSE Comparison
cat("PCA Model RMSEs:\n")
print(pexp_pca_model_performance)


#Exp PCA Needs to Be Adjusted----
decay_factor <- 0.85
season_pexp <- game_wr %>%
  group_by(receiver_player_id, posteam) %>%
  arrange(game_id, week) %>% # Ensure the data is sorted by week
  mutate(
    exp_weighted_targets_per_game = zoo::rollapply(targets, width = 10, 
                                                   FUN = function(x) sum(x * decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
    exp_weighted_air_yards_per_game = zoo::rollapply(total_air_yards, width = 10, 
                                                     FUN = function(x) sum(x * decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
    exp_weighted_yards_game_per_game = zoo::rollapply(total_yards_gained, width = 10, 
                                                      FUN = function(x) sum(x * decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
    exp_weighted_longest_rec_per_game = zoo::rollapply(longest_rec, width = 10, 
                                                       FUN = function(x) sum(x * decay_factor^(rev(seq_along(x)) - 1)), fill = NA, align = "right"),
    max_week = max(week), # Find the most recent week
    week_weight = (week - min(week) + 1) / (max_week - min(week) + 1), # Scale weights to give higher value to recent weeks
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
    # Weighted averages for per-game metrics
    weighted_targets_per_game = cumsum(targets * week_weight) / cumsum(week_weight),
    weighted_air_yards_per_game = cumsum(total_air_yards * week_weight) / cumsum(week_weight),
    weighted_yards_game_per_game = cumsum(total_yards_gained * week_weight) / cumsum(week_weight),
    weighted_yac_per_game = cumsum(YAC * week_weight) / cumsum(week_weight),
    weighted_epa_per_game = cumsum(total_epa * week_weight) / cumsum(week_weight),
    weighted_success_rate_per_game = cumsum(success_rate * week_weight) / cumsum(week_weight),
    weighted_longest_rec_per_game = cumsum(longest_rec * week_weight) / cumsum(week_weight),
    weighted_pct_air_yards_per_game = cumsum(pct_air_yards * week_weight) / cumsum(week_weight),
    weighted_target_share_per_game = cumsum(target_share * week_weight) / cumsum(week_weight),
    weighted_adot_per_game = cumsum(adot * week_weight) / cumsum(week_weight),
    
    # Cumulative metrics
    cum_targets = cumsum(targets),
    cum_target_rate = cumsum(catchable_targ) / cumsum(targets),
    cum_contested_rate = cumsum(contested_ball) / cumsum(targets),
    cum_created_reception_rate = cumsum(created_rec) / cumsum(targets)
  ) %>%
  ungroup() %>% 
  group_by(receiver_player_id) %>%
  mutate(next_week = lead(week)) %>% 
  filter(next_week > 1) %>% 
  filter(!is.na(week)) %>% 
  select(
    week = next_week, receiver_player_id, year, posteam,
    ends_with("per_game"), starts_with("cum")
  )
longest_rec_joined_pexp_data <- longest_rec_data %>% 
  left_join(season_pexp, by = c("receiver_player_id","week","year","posteam")) %>% 
  left_join(off_scout, by = c("week","posteam","year")) %>%
  left_join(def_scout, by = c("week","defteam","year")) %>% 
  ungroup() %>% 
  select(-posteam,-defteam) %>% 
  filter(week>1)
pexp_data <- longest_rec_joined_pexp_data %>%
  drop_na() %>%  # Remove rows with missing values
  mutate_if(is.character, as.factor) %>%  # Convert character columns to factors
  select(-receiver_player_id, -game_id, -name) %>%  # Drop non-predictive or irrelevant columns
  mutate(longest_rec = sqrt(ifelse(longest_rec<0,0,longest_rec)))
# Split data into training and testing sets
pexp_trainIndex <- createDataPartition(pexp_data$longest_rec, p = 0.8, list = FALSE)
pexp_trainData <- pexp_data[pexp_trainIndex, ]
pexp_testData <- pexp_data[-pexp_trainIndex, ]

# Separate predictors and response
pexp_x_train <- pexp_trainData %>% select(-longest_rec)
pexp_y_train <- pexp_trainData$longest_rec
pexp_x_test <- pexp_testData %>% select(-longest_rec)
pexp_y_test <- pexp_testData$longest_rec


pexp_pca_preProc <- preProcess(pexp_x_train, method = "pca", pcaComp = 10) # Keep top 10 principal components
pexp_x_train_pca <- predict(pexp_pca_preProc, pexp_x_train)
pexp_x_test_pca <- predict(pexp_pca_preProc, pexp_x_test)

# Random Forest with PCA
pexp_pca_rf_model <- train(
  longest_rec ~ ., 
  data = data.frame(longest_rec = pexp_y_train, pexp_x_train_pca),
  method = "rf",
  trControl = trainControl(method = "cv", number = 5),
  tuneGrid = expand.grid(mtry = seq(2, ncol(pexp_x_train_pca), by = 2)),
  ntree = 500
)

# SVM with PCA
pexp_pca_svm_model <- train(
  longest_rec ~ ., 
  data = data.frame(longest_rec = pexp_y_train, pexp_x_train_pca),
  method = "svmRadial",
  trControl = trainControl(method = "cv", number = 5),
  tuneGrid = expand.grid(
    sigma = 0.01, # SVM RBF kernel hyperparameter
    C = seq(1, 10, by = 1) # Regularization strength
  )
)

# Gradient Boosting (XGBoost) with PCA
pexp_pca_xgb_model <- train(
  longest_rec ~ ., 
  data = data.frame(longest_rec = pexp_y_train, pexp_x_train_pca),
  method = "xgbTree",
  trControl = trainControl(method = "cv", number = 5),
  tuneGrid = expand.grid(
    nrounds = c(100, 200),
    max_depth = c(3, 5, 7),
    eta = c(0.01, 0.1),
    gamma = 0,
    colsample_bytree = 0.8,
    min_child_weight = 1,
    subsample = 0.8
  )
)

# KNN with PCA
pexp_pca_knn_model <- train(
  longest_rec ~ ., 
  data = data.frame(longest_rec = pexp_y_train, pexp_x_train_pca),
  method = "knn",
  trControl = trainControl(method = "cv", number = 5),
  tuneGrid = expand.grid(k = seq(1, 20, by = 2)) # Number of neighbors
)

# Combine all PCA models for comparison
pexp_pca_models <- list(
  rf = pexp_pca_rf_model,
  svm = pexp_pca_svm_model,
  xgb = pexp_pca_xgb_model,
  knn = pexp_pca_knn_model
)

# Evaluate PCA Models
pexp_pca_results <- resamples(pexp_pca_models)
summary(pexp_pca_results)

# Compare RMSE on the test set
pexp_pca_model_performance <- sapply(pexp_pca_models, function(model) {
  predictions <- predict(model, newdata = pexp_x_test_pca)
  rmse(pexp_y_test, predictions)
})

# Select the Best PCA Model
pexp_best_pca_model_name <- names(pexp_pca_model_performance)[which.min(pexp_pca_model_performance)]
pexp_best_pca_model <- pexp_pca_models[[pexp_best_pca_model_name]]

cat("Best PCA model:", pexp_best_pca_model_name, "with RMSE:", min(pexp_pca_model_performance), "\n")

# Feature Importance (if applicable)
if (pexp_best_pca_model_name %in% c("rf", "xgb")) {
  varImp(pexp_best_pca_model) %>% plot()
}

# Final Predictions with Best PCA Model
pexp_final_predictions_pca <- as.data.frame(predict(pexp_best_pca_model, newdata = pexp_x_test_pca))

# Residual Analysis
pexp_residuals_pca <- pexp_y_test - pexp_final_predictions_pca[[1]]

# Output RMSE Comparison
cat("PCA Model RMSEs:\n")
print(pexp_pca_model_performance)


#Model Comparison----
# Reverse the square root transformation
sqrt_final_predictions_original_scale <- sqrt_final_predictions[[1]]^2
sqrt_residuals_original_scale <- sqrt_y_test^2 - sqrt_final_predictions_original_scale

pca_predictions_original_scale <- pexp_final_predictions_pca[[1]]^2
pca_residuals_original_scale <- pexp_y_test^2 - pca_predictions_original_scale
exp_predictions_original_scale <- exp_final_predictions[[1]]^2
exp_residuals_original_scale <- exp_y_test^2 - exp_predictions_original_scale
pexp_predictions_original_scale <- pexp_final_predictions[[1]]^2
pexp_residuals_original_scale <- pexp_y_test^2 - pexp_predictions_original_scale

# RMSE and MAE for the original model
original_rmse <- rmse(y_test, final_predictions[[1]])
original_mae <- mae(y_test, final_predictions[[1]])

# RMSE and MAE for the transformed model (on original scale)

sqrt_rmse <- rmse(sqrt_y_test^2, sqrt_final_predictions_original_scale)
sqrt_mae <- mae(sqrt_y_test^2, sqrt_final_predictions_original_scale)

pexp_rmse <- rmse(pexp_y_test^2, pexp_final_predictions_original_scale)
pexp_mae <- mae(pexp_y_test^2, pexp_final_predictions_original_scale)

weighted_rmse <- rmse(weighted_y_test^2, weighted_final_predictions_original_scale)
weighted_mae <- mae(weighted_y_test^2, weighted_final_predictions_original_scale)

pca_rmse <- rmse(pexp_y_test^2,pca_predictions_original_scale)
pca_mae <- mae(pexp_y_test^2,pca_predictions_original_scale)

nt_pca_rmse <- rmse(nt_pexp_y_test,nt_pexp_final_predictions_pca[[1]])
nt_pca_mae <- mae(nt_pexp_y_test,nt_pexp_final_predictions_pca[[1]])

exp_rmse <- rmse(exp_y_test^2,exp_predictions_original_scale)
exp_mae <- mae(exp_y_test^2,exp_predictions_original_scale)

#Look to explore predicting untransformed longest_rec
cat("Original Model RMSE:", original_rmse, "\n")
cat("Original Model MAE:", original_mae, "\n")
cat("Square Root Model RMSE (Original Scale):", sqrt_rmse, "\n")
cat("Square Root Model MAE (Original Scale):", sqrt_mae, "\n")
cat("Combo Root Model RMSE (Original Scale):", combo_rmse, "\n")
cat("Combo Root Model MAE (Original Scale):", combo_mae, "\n")
# cat("Weighted Model RMSE (Original Scale):", weighted_rmse, "\n")
# cat("Weighted Model MAE (Original Scale):", weighted_mae, "\n")
cat("PCA Model RMSE (Original Scale):", pca_rmse, "\n")
cat("PCA Model MAE (Original Scale):", pca_mae, "\n")
cat("Exp Model RMSE (Original Scale):", exp_rmse, "\n")
cat("Exp Model MAE (Original Scale):", exp_mae, "\n")
cat("Pexp Model RMSE (Original Scale):", pexp_rmse, "\n")
cat("Pexp Model MAE (Original Scale):", pexp_mae, "\n")
cat("Nt_Pexp Model RMSE (Original Scale):", nt_pca_rmse, "\n")
cat("Nt_Pexp Model MAE (Original Scale):", nt_pca_mae, "\n")


# Residuals plot for the original model
nt_pexp_original_residuals <- nt_pexp_y_test-nt_pexp_final_predictions_pca[[1]]
original_residuals <- y_test - final_predictions[[1]]
plot(nt_pexp_final_predictions_pca[[1]], nt_pexp_original_residuals, main = "Residuals vs Predictions (Original Model)", xlab = "Predicted Values", ylab = "Residuals")
abline(h = 0, col = "red")
# Residuals plot for the square root-transformed model (on original scale)
plot(pca_predictions_original_scale, pca_residuals_original_scale, main = "Residuals vs Predictions (Sqrt Model)", xlab = "Predicted Values", ylab = "Residuals")
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
pbp <- load_pbp(2024)
pbp <- pbp %>%
  left_join(ftn_data, by = c("game_id" = "nflverse_game_id",
                             "play_id" = "nflverse_play_id")) %>% 
  left_join(participation,by = c("game_id" = "nflverse_game_id",
                                 "play_id" = "play_id"))

pbp_rp <- pbp %>%
  filter(pass == 1 | rush == 1) %>%
  filter(qb_kneel == 0,qb_spike == 0) %>% 
  filter(!is.na(epa)) %>% 
  mutate(blitz = ifelse(n_pass_rushers> 4,1,0),
         lightbox = ifelse(n_defense_box<=6,1,0),
         heavybox = ifelse(n_defense_box>=8,1,0)) %>% 
  mutate(explosive = ifelse((yardline_100 <  20 & pass_attempt == 1) | (yardline_100<12 & (qb_scramble ==1 |rush ==1)), NA, ifelse((yards_gained>20 & pass_attempt == 1) | (yards_gained >12 & (qb_scramble == 1 | rush == 1)),1,0)),
         negative = ifelse(yards_gained < 0, 1,0))

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
decay_factor <- 0.7 #change if necessary
season_wr <- game_wr %>%
  group_by(receiver_player_id, posteam, name) %>% 
  mutate(max_week = max(week),
  week_weight = (week - min(week) + 1) / (max_week - min(week) + 1)) %>% 
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
    # Weighted averages for per-game metrics
    weighted_targets_per_game = sum(targets * week_weight) / sum(week_weight),
    weighted_air_yards_per_game = sum(total_air_yards * week_weight) / sum(week_weight),
    weighted_yards_game_per_game = sum(total_yards_gained * week_weight) / sum(week_weight),
    weighted_yac_per_game = sum(YAC * week_weight) / sum(week_weight),
    weighted_epa_per_game = sum(total_epa * week_weight) / sum(week_weight),
    weighted_success_rate_per_game = sum(success_rate * week_weight) / sum(week_weight),
    weighted_longest_rec_per_game = sum(longest_rec * week_weight) / sum(week_weight),
    weighted_pct_air_yards_per_game = sum(pct_air_yards * week_weight) / sum(week_weight),
    weighted_target_share_per_game = sum(target_share * week_weight) / sum(week_weight),
    weighted_adot_per_game = sum(adot * week_weight) / sum(week_weight),
    #Exponential weights
    exp_weighted_targets_per_game = sum(targets * decay_factor^(max_week - week)), # Exponential weights
    exp_weighted_air_yards_per_game = sum(total_air_yards * decay_factor^(max_week - week)),
    exp_weighted_yards_game_per_game = sum(total_yards_gained * decay_factor^(max_week - week)),
    exp_weighted_longest_rec_per_game = sum(longest_rec * decay_factor^(max_week - week)),
    
    # cumulative metrics
    cum_targets = sum(targets),
    cum_target_rate = sum(catchable_targ) / sum(targets),
    cum_contested_rate = sum(contested_ball) / sum(targets),
    cum_created_reception_rate = sum(created_rec) / sum(targets)
  ) %>%
  ungroup() %>%
  select(
    receiver_player_id, posteam, name, 
    ends_with("per_game"), starts_with("cum"),starts_with("weighted"),starts_with("exp")
  )

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
    off_cum_avg_epa_per_play = sum(total_epa_in_week,na.rm = T) /sum(plays_in_week,na.rm = T), # Cumulative avg EPA/play
    off_cum_avg_epa_per_db = sum(total_db_epa,na.rm = T) /sum(total_dropbacks,na.rm = T),
    off_cum_avg_success_per_db = sum(total_db_success,na.rm = T) /sum(total_dropbacks,na.rm = T),
    off_cum_avg_success_rate = sum(total_success,na.rm = T) /sum(plays_in_week,na.rm = T),
    off_cum_proe = sum(total_proe,na.rm = T) /sum(plays_in_week,na.rm = T),
    off_cum_longest_rec_pg = mean(off_game_longest_rec,na.rm = T)
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
    def_cum_avg_epa_per_play = sum(total_epa_in_week,na.rm = T) /sum(plays_in_week,na.rm = T), # Cumulative avg EPA/play
    def_cum_avg_epa_per_db = sum(total_db_epa,na.rm = T) /sum(total_dropbacks,na.rm = T),
    def_cum_avg_success_per_db = sum(total_db_success,na.rm = T) /sum(total_dropbacks,na.rm = T),
    def_cum_avg_success_rate = sum(total_success,na.rm = T) /sum(plays_in_week,na.rm = T),
    def_cum_proe = sum(total_proe,na.rm = T) /sum(plays_in_week,na.rm = T),
    def_cum_longest_rec_pg = mean(def_game_longest_rec,na.rm = T)
  ) %>% 
  select(defteam, starts_with("def_cum"))# Remove rows year = seasonwhere the next week doesn't exist


schedules <- load_schedules(2024)

wr_sched <- schedules %>% 
  # filter(is.na(away_score), week == pbp_rp$week +1)
  filter(is.na(away_score)) %>% 
  mutate(min_week = min(week)) %>% 
  filter(week == min_week) %>% 
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
# Step 1: Apply PCA to Full Dataset
nt_pexp_full_x <- nt_pexp_data %>% select(-longest_rec)
nt_pexp_full_y <- nt_pexp_data$longest_rec

# Apply the same PCA preprocessing fitted earlier to the full dataset
nt_pexp_full_x_pca <- predict(nt_pexp_pca_preProc, nt_pexp_full_x)

# Combine the transformed predictors with the target variable
nt_pexp_full_data <- data.frame(longest_rec = nt_pexp_full_y, nt_pexp_full_x_pca)

# Step 2: Train the Best PCA Model on Full Dataset
# Assuming the best model is the PCA-based Random Forest
best_xgb_params <- nt_pexp_pca_xgb_model$bestTune

# Train the final XGBoost model using the best parameters
final_model <- train(
  longest_rec ~ ., 
  data = nt_pexp_full_data,
  method = "xgbTree",
  trControl = trainControl(method = "none"), # No cross-validation for the final model
  tuneGrid = best_xgb_params # Use the best parameters
)
# Step 3: Save Preprocessing and Model
saveRDS(nt_pexp_pca_preProc, file = "final_pca_preprocessing.rds")
saveRDS(final_model, file = "final_model.rds")

# Step 4: Validate Model on Independent Test Data
# Apply PCA to test data using the saved preprocessing steps
test_x_pca <- predict(nt_pexp_pca_preProc, nt_pexp_x_test)

# Make predictions using the final trained model
final_test_predictions <- predict(final_model, newdata = test_x_pca)

# Reverse transformations if needed (e.g., square-root)
final_test_predictions_original_scale <- final_test_predictions

# Calculate RMSE and MAE for the final model
final_test_rmse <- rmse(nt_pexp_y_test, final_test_predictions_original_scale)
final_test_mae <- mae(nt_pexp_y_test, final_test_predictions_original_scale)

cat("Final PCA Model RMSE on Test Data:", final_test_rmse, "\n")
cat("Final PCA Model MAE on Test Data:", final_test_mae, "\n")

#Predictions ----

# Load the saved preprocessing object and final model
final_preProc <- readRDS("final_pca_preprocessing.rds")
final_model <- readRDS("final_model.rds")

# Preprocess new data
processed_new_data <- predict(final_preProc, newdata = longest_rec_joined_data_pred %>% drop_na())

# Make predictions
predictions_df <- longest_rec_joined_data_pred %>%filter(complete.cases(.)) %>%   
  mutate(longest_reception_predictions = predict(final_model, newdata = processed_new_data %>% select(starts_with("PC"))))

fit_gamma <- fitdist(longest_rec_joined_data$longest_rec, "gamma")

# Extract parameters
shape <- fit_gamma$estimate["shape"]
rate <- fit_gamma$estimate["rate"]

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
      betting_line,       # Transformed market line
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
    pred_longest_rec = longest_reception_predictions) %>% 
    select(posteam,name, pred_longest_rec,betting_line, predicted_over_prob,`Over ROI %`, `Under ROI %`)
  return(results)
}
longest_reception_prediction_output()


