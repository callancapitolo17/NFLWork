library(nflfastR)
library(ggplot2)
library(tidyverse)
library(ggimage)
library(ggthemes)
library(dplyr)
library(ggrepel)
library(nflreadr)
library(gt)
library(ggrepel)
library(nflplotR)
library(gtExtras)

prob_to_decimal <- function(prob){ifelse(prob > 0, 1 / prob, NA)}
decimal_to_american = function(fair_decimal_odds){
  ifelse(
  fair_decimal_odds > 2.0,
  (fair_decimal_odds - 1) * 100,          # Positive American odds
  -100 / (fair_decimal_odds - 1)          # Negative American odds
)}
prob_to_american <- function(prob){
  fair_decimal_odds <- ifelse(prob > 0, 1 / prob, NA)
  ifelse(
    fair_decimal_odds > 2.0,
    (fair_decimal_odds - 1) * 100,          # Positive American odds
    -100 / (fair_decimal_odds - 1)) # Negative American odds
}
american_to_decimal_odds <- function(odds){ifelse(odds < 0, 1 - (100/odds), 1+(odds/100))}
no_vig_odds <- function(odds_1,odds_2) {
  # Ensure odds are numeric
  prob_odds_1 <- 1/american_to_decimal_odds(odds_1)
  prob_odds_2 <- 1/american_to_decimal_odds(odds_2)
  
  # Remove the vig (normalize so that the sum is 1)
  no_vig_odds_1 <- prob_odds_1 / (prob_odds_1+prob_odds_2)
  no_vig_odds_2 <- prob_odds_2 / (prob_odds_1+prob_odds_2)
  
  return(c(no_vig_odds_1,no_vig_odds_2))
}



nfl99all <- load_pbp(1999:2024)
nfl99 <- nfl99all %>% 
  filter(pass == 1 | rush == 1) %>% 
  mutate(explosive = ifelse((yardline_100 <  20 & pass_attempt == 1) | (yardline_100<12 & (qb_scramble ==1 |rush ==1)), NA, ifelse((yards_gained>20 & pass_attempt == 1) | (yards_gained >12 & (qb_scramble == 1 | rush == 1)),1,0)),
         negative = ifelse(yards_gained < 0, 1,0))

# First TD Jersey Number----
pbp_rp %>% 
  filter(td_team %in% c("KC", "PHI")) %>% 
  group_by(game_id) %>% 
  mutate(td_count = cumsum(touchdown)) %>% 
  mutate(jersey = case_when(
    passer_id == td_player_id ~ passer_jersey_number,
    rusher_id == td_player_id ~ rusher_jersey_number,
    receiver_id == td_player_id ~ receiver_jersey_number,
    TRUE ~ NA_real_
  )) %>% 
  ungroup() %>% 
  # group_by(td_team) %>% 
  summarize(mean(jersey[td_count == 1],na.rm = T), median(jersey[td_count == 1],na.rm = T), sd(jersey[td_count == 1],na.rm = T))


pbp_rp %>% 
  filter(td_team %in% c("KC", "PHI"), touchdown == 1) %>% 
  group_by(game_id) %>% 
  # slice_tail(n = 1) %>% 
  mutate(jersey = case_when(
    passer_id == td_player_id   ~ passer_jersey_number,
    rusher_id == td_player_id   ~ rusher_jersey_number,
    receiver_id == td_player_id ~ receiver_jersey_number,
    TRUE ~ NA_real_
  )) %>% 
  ungroup() %>% 
  summarize(mean_jersey = mean(jersey, na.rm = TRUE),
            median_jersey = median(jersey, na.rm = TRUE),
            sd_jersey = sd(jersey, na.rm = TRUE))

pbp_rp %>% 
  filter(td_team %in% c("KC", "PHI"), touchdown == 1) %>% 
  group_by(game_id) %>% 
  slice_tail(n = 1) %>%
  mutate(jersey = case_when(
    passer_id == td_player_id   ~ passer_jersey_number,
    rusher_id == td_player_id   ~ rusher_jersey_number,
    receiver_id == td_player_id ~ receiver_jersey_number,
    TRUE ~ NA_real_
  )) %>% 
  ungroup() %>% 
  summarize(mean_jersey = mean(jersey, na.rm = TRUE),
            median_jersey = median(jersey, na.rm = TRUE),
            sd_jersey = sd(jersey, na.rm = TRUE))





#Timeouts----
pbp %>% 
  filter(timeout_team %in% c("PHI","KC")) %>% 
  filter(game_half == "Half1") %>% 
  group_by(game_id) %>%
  mutate(tos = cumsum(timeout)) %>% 
  ungroup() %>% 
  group_by(timeout_team) %>% 
  summarize(first_to_mean = mean(half_seconds_remaining[tos == 1]), first_to_median = median(half_seconds_remaining[tos == 1]))

nfl99all %>%
  filter(season>=2021) %>% 
  filter(timeout_team %in% c("PHI","KC")) %>% 
  filter(game_half == "Half1") %>% 
  group_by(game_id) %>%
  mutate(tos = cumsum(timeout)) %>% 
  ungroup() %>% 
  group_by(timeout_team) %>% 
  summarize(first_to_mean = mean(half_seconds_remaining[tos == 1]),first_to_median = median(half_seconds_remaining[tos == 1]))

library(dplyr)

timeout_stats <- pbp %>%
  filter(timeout_team %in% c("PHI", "KC"),
         game_half == "Half1") %>%
  group_by(game_id) %>%
  mutate(tos = cumsum(timeout)) %>% 
  ungroup() %>% 
  group_by(game_id, timeout_team) %>%
  filter(tos == 1) %>%  # First timeout for that team in each game
  summarize(timeout_time = half_seconds_remaining, .groups = "drop") %>%
  group_by(timeout_team) %>%
  summarize(mean_time = mean(timeout_time),
            sd_time = sd(timeout_time))
print(timeout_stats)

mu_KC <- timeout_stats %>% filter(timeout_team == "KC") %>% pull(mean_time)
mu_PHI <- timeout_stats %>% filter(timeout_team == "PHI") %>% pull(mean_time)
sd_KC <- timeout_stats %>% filter(timeout_team == "KC") %>% pull(sd_time)
sd_PHI <- timeout_stats %>% filter(timeout_team == "PHI") %>% pull(sd_time)

# Calculate the difference in means and the combined standard deviation
delta <- mu_KC - mu_PHI
combined_sd <- sqrt(sd_KC^2 + sd_PHI^2)

# Calculate the probability that the Chiefs call their timeout earlier
p_chiefs <- pnorm(delta / combined_sd)


timeout_stats <- nfl99all %>%
  filter(season>=2021) %>% 
  filter(timeout_team %in% c("PHI", "KC"),
         game_half == "Half1") %>%
  group_by(game_id) %>%
  mutate(tos = cumsum(timeout)) %>% 
  ungroup() %>% 
  group_by(game_id, timeout_team) %>%
  filter(tos == 1) %>%  # First timeout for that team in each game
  summarize(timeout_time = half_seconds_remaining, .groups = "drop") %>%
  group_by(timeout_team) %>%
  summarize(mean_time = mean(timeout_time,na.rm = T),
            sd_time = sd(timeout_time, na.rm = T))
print(timeout_stats)

mu_KC <- timeout_stats %>% filter(timeout_team == "KC") %>% pull(mean_time)
mu_PHI <- timeout_stats %>% filter(timeout_team == "PHI") %>% pull(mean_time)
sd_KC <- timeout_stats %>% filter(timeout_team == "KC") %>% pull(sd_time)
sd_PHI <- timeout_stats %>% filter(timeout_team == "PHI") %>% pull(sd_time)

# Calculate the difference in means and the combined standard deviation
delta <- mu_KC - mu_PHI
combined_sd <- sqrt(sd_KC^2 + sd_PHI^2)

# Calculate the probability that the Chiefs call their timeout earlier
p_chiefs <- pnorm(delta / combined_sd)

# Load required package
library(dplyr)

# -----------------------------
# 1. Get Each Team's First Timeout in the First Half
# -----------------------------
# Assuming 'pbp' is your play-by-play dataset.

first_timeouts <- nfl99all   %>%
  filter(season>=2021) %>% 
  filter(timeout_team %in% c("PHI", "KC"),
         game_half == "Half1") %>% 
  group_by(game_id, timeout_team) %>%
  mutate(tos = cumsum(timeout)) %>%
  filter(tos == 1) %>% 
  ungroup()

# -----------------------------
# 2. Calculate Summary Statistics for Each Team
# -----------------------------
timeout_stats <- first_timeouts %>%
  group_by(timeout_team) %>%
  summarize(
    mean_time = mean(half_seconds_remaining, na.rm = TRUE),
    sd_time = sd(half_seconds_remaining, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

# Extract summary statistics
mu_KC <- timeout_stats %>% filter(timeout_team == "KC") %>% pull(mean_time)
mu_PHI <- timeout_stats %>% filter(timeout_team == "PHI") %>% pull(mean_time)
sd_KC <- timeout_stats %>% filter(timeout_team == "KC") %>% pull(sd_time)
sd_PHI <- timeout_stats %>% filter(timeout_team == "PHI") %>% pull(sd_time)

# -----------------------------
# 3. Monte Carlo Simulation
# -----------------------------

n_simulations <- 100000  # Number of simulated games

# Simulate first timeout times for each team
sim_KC <- rnorm(n_simulations, mean = mu_KC, sd = sd_KC)
sim_PHI <- rnorm(n_simulations, mean = mu_PHI, sd = sd_PHI)

# Determine in how many simulations the Chiefs call the timeout first
# (i.e., KC has more seconds remaining)
chiefs_first_count <- sum(sim_KC > sim_PHI)

# Estimate probability
p_chiefs_simulated <- chiefs_first_count / n_simulations

# -----------------------------
# 4. Compare to Market Odds & Calculate Edge
# -----------------------------
implied_prob <- 120 / (120 + 100)  # Implied probability from -120 odds (~54.55%)
edge <- p_chiefs_simulated - implied_prob

# -----------------------------
# 5. Confidence Intervals for the Simulated Probability
# -----------------------------
# Using the normal approximation for the binomial proportion
se <- sqrt(p_chiefs_simulated * (1 - p_chiefs_simulated) / n_simulations)
ci_lower <- p_chiefs_simulated - 1.96 * se
ci_upper <- p_chiefs_simulated + 1.96 * se

# -----------------------------
# 6. Print the Results
# -----------------------------
cat("Chiefs' Mean First Timeout Time (seconds remaining):", mu_KC, "\n")
cat("Eagles' Mean First Timeout Time (seconds remaining):", mu_PHI, "\n")
cat("Estimated Probability Chiefs Timeout First:", round(p_chiefs_simulated, 4), "\n")
cat("95% Confidence Interval:", round(ci_lower, 4), "-", round(ci_upper, 4), "\n")
cat("Market Implied Probability (from -120 odds):", round(implied_prob, 4), "\n")
cat("Edge (Estimated - Implied):", round(edge, 4), "\n")

# Load required packages

# Load required packages
library(dplyr)
library(rstan)

# Calculate mean and variance for each team
team_stats <-  pbp %>% 
  filter(timeout_team %in% c("PHI", "KC"),
         game_half == "Half1") %>% 
  group_by(game_id, timeout_team) %>%
  mutate(tos = cumsum(timeout)) %>%
  filter(tos == 1) %>% 
  ungroup() %>% 
  summarize(
    mean_time = mean(half_seconds_remaining + 1),  # Add 1 to avoid zero
    var_time = var(half_seconds_remaining + 1)
  )

# Calculate alpha and beta using method of moments
team_stats <- team_stats %>%
  mutate(
    alpha = (mean_time^2) / var_time,
    beta = mean_time / var_time
  )

# Display the results
print(team_stats)


# -----------------------------
# 1. Prepare the Data
# -----------------------------
first_timeouts <- pbp %>% 
  filter(timeout_team %in% c("PHI", "KC"),
         game_half == "Half1") %>% 
  group_by(game_id, timeout_team) %>%
  mutate(tos = cumsum(timeout)) %>%
  filter(tos == 1) %>% 
  ungroup()

# Recode teams: KC = 1, PHI = 2
stan_data <- first_timeouts %>%
  mutate(team_code = ifelse(timeout_team == "KC", 1, 2)) %>%
  select(half_seconds_remaining, team_code)

# Prepare data for Stan
N <- nrow(stan_data)
y <- stan_data$half_seconds_remaining + 1  # Add 1 to avoid issues with zeros
team <- stan_data$team_code

data_list <- list(
  N = N,
  y = y,
  team = team
)

# -----------------------------
# 2. Stan Model (Gamma Distribution)
# -----------------------------
stan_model_code <- "
data {
  int<lower=0> N;                    // number of observations
  vector<lower=0>[N] y;              // timeout times (seconds remaining)
  int<lower=1,upper=2> team[N];      // team indicator: 1 for KC, 2 for PHI
}
parameters {
  vector<lower=0>[2] alpha;          // shape parameters for each team
  vector<lower=0>[2] beta;           // rate parameters for each team
}
model {
  // Priors (adjust these based on your season data)
  alpha ~ normal(0.57, 0.1);              // Weakly informative prior for shape
  beta ~ normal(0.0017, 0.005);        // Prior for rate (controls skewness)

  // Likelihood: Gamma distribution for each team
  for (n in 1:N) {
    y[n] ~ gamma(alpha[team[n]], beta[team[n]]);
  }
}
generated quantities {
  real mean_KC = alpha[1] / beta[1];     // Mean timeout time for KC
  real mean_PHI = alpha[2] / beta[2];    // Mean timeout time for PHI
  real diff = mean_KC - mean_PHI;        // Difference in means
}
"

# -----------------------------
# 3. Compile and Run the Stan Model
# -----------------------------
fit <- stan(model_code = stan_model_code, data = data_list,
            iter = 4000, chains = 4, seed = 123)

# -----------------------------
# 4. Extract Posterior Samples and Compute Probabilities
# -----------------------------
post <- rstan::extract(fit)
prob_KC_earlier <- mean(post$diff > 0)

# -----------------------------
# 5. Compare to Market Odds and Calculate the Edge
# -----------------------------
implied_prob <- 120 / (120 + 100)  # ~54.55%
edge <- prob_KC_earlier - implied_prob

# -----------------------------
# 6. Print Results
# -----------------------------
cat("Posterior probability that Chiefs call timeout earlier:", round(prob_KC_earlier, 4), "\n")
cat("Market implied probability (from -120 odds):", round(implied_prob, 4), "\n")
cat("Edge (Posterior probability - Market implied):", round(edge, 4), "\n")


#Safety----
safety_prob <- nfl99all %>% 
  filter(season>=2024) %>% 
  group_by(game_id) %>% 
  summarize(safeties = sum(safety,na.rm = T)) %>% 
  mutate(bin_safety = ifelse(safeties ==0,0,1)) %>% 
  ungroup() %>% 
  summarize(mean_safety = mean(bin_safety)) %>% 
  pull(mean_safety)
prob_to_american(safety_prob)

#Overtime----
overtime_data <- load_schedules(1999:2024) %>% 
  filter(!is.na(away_score), !is.na(home_moneyline)) %>% 
  mutate(
    home_prob = map2_dbl(home_moneyline, away_moneyline, ~ no_vig_odds(.x, .y)[1]),
    away_prob = map2_dbl(home_moneyline, away_moneyline, ~ no_vig_odds(.x, .y)[2])
  ) %>% 
  select(overtime,total_line,home_prob,spread_line) %>% 
  mutate(overtime = ifelse(overtime == 1, "Yes","No")) 

library(caret)
library(pROC)

set.seed(123)  # For reproducibility
train_index <- createDataPartition(overtime_data$overtime, p = 0.8, list = FALSE)
train_data  <- overtime_data[train_index, ]
test_data   <- overtime_data[-train_index, ]

# Set up cross-validation parameters (10-fold CV)
cv_control <- trainControl(method = "cv", 
                           number = 10, 
                           summaryFunction = twoClassSummary,
                           classProbs = TRUE)  # required for ROC summary

# Fit the logistic regression model using caret's train() function.
# We use the predictors spread, total, and location.
set.seed(123)  # For reproducibility of CV
logistic_model <- train(overtime ~ .+spread_line:home_prob,
                        data = train_data,
                        method = "glmnet",
                        family = "binomial",
                        trControl = cv_control,
                        preProcess = c("center", "scale"),
                        metric = "ROC")  # Optimize for AUC

# Retrieve the best lambda from the tuning results
best_lambda <- logistic_model$bestTune$lambda

# Extract the coefficients at that lambda
coef_matrix <- coef(logistic_model$finalModel, s = best_lambda)
print(coef_matrix)

var_imp <- varImp(logistic_model)
print(var_imp)
plot(var_imp, main = "Variable Importance - glmnet Model")


# Obtain predicted probabilities using the final glmnet model:
# Here we assume that test_data has only the predictors (without the outcome) or that you subset appropriately.
test_data$pred_prob <- predict(logistic_model, newdata = test_data, type = "prob")[[2]]

library(caret)
library(ggplot2)

# Assume test_data has the true outcome in "overtime" and predicted probability in "pred_prob"
# Convert overtime to numeric for computing observed frequencies (Yes -> 1, No -> 0)
test_data$overtime_numeric <- ifelse(test_data$overtime == "Yes", 1, 0)

# Create bins for the predicted probabilities
test_data$prob_bin <- cut(test_data$pred_prob, breaks = seq(round(quantile(test_data$pred_prob,.005)[[1]],2), round(quantile(test_data$pred_prob,.995)[[1]],2), by = 0.005), include.lowest = TRUE)

# Compute observed frequency for each bin
calibration_df <- test_data %>%
  group_by(prob_bin) %>%
  summarize(mean_pred = mean(pred_prob),
            obs_freq = mean(overtime_numeric),
            count = n())

# Plot the calibration curve
ggplot(calibration_df, aes(x = mean_pred, y = obs_freq)) +
  geom_point(size = 3) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Mean Predicted Probability", y = "Observed Frequency",
       title = "Calibration Plot") +
  theme_minimal()


# Calculate the Brier score
brier_score <- mean((test_data$pred_prob - test_data$overtime_numeric)^2)
cat("Brier Score:", brier_score, "\n")

basic_brier_score <- mean((mean(test_data$overtime_numeric)-test_data$overtime_numeric)^2)
1-(brier_score/basic_brier_score)

final_model <- train(
  overtime ~ . + spread_line:home_prob,   # include the interaction if desired
  data = overtime_data,
  method = "glmnet",
  tuneGrid = logistic_model$bestTune,       # use the best tuning parameters from earlier
  preProcess = c("center", "scale"),
  family = "binomial"
)
# Calculate the Brier score.
brier_score <- mean((test_data$pred_prob - test_data$overtime_numeric)^2)
cat("Final Model Brier Score:", brier_score, "\n")

# Optionally, compute a baseline Brier score (e.g., always predicting the overall overtime rate).
baseline_prob <- mean(test_data$overtime_numeric)
baseline_brier <- mean((baseline_prob - test_data$overtime_numeric)^2)
cat("Baseline Brier Score:", baseline_brier, "\n")

# Calculate the Brier Skill Score (BSS), where positive values indicate improvement over baseline.
brier_skill_score <- 1 - (brier_score / baseline_brier)
cat("Final Model Brier Skill Score:", brier_skill_score, "\n")

#Super bowl Prediction
prediction_data <- load_schedules(2024) %>% filter(is.na(away_score)) %>% 
  mutate(
    home_prob = map2_dbl(home_moneyline, away_moneyline, ~ no_vig_odds(.x, .y)[1]),
    away_prob = map2_dbl(home_moneyline, away_moneyline, ~ no_vig_odds(.x, .y)[2])
  )
  
prob_to_american(predict(final_model, newdata = prediction_data, type = "prob")[[2]])
prob_to_american(mean(overtime_data$overtime == "Yes",na.rm = T))



library(caret)
library(ranger)

set.seed(123)  # For reproducibility

tune_grid <- expand.grid(
  mtry = c(1, 2, 3),        # try 1, 2, and 3 variables at each split
  splitrule = "gini",       # for classification, "gini" is typical
  min.node.size = c(1, 5)    # try two different minimum node sizes
)

set.seed(123)
rf_model <- train(
  overtime ~ .,
  data = train_data,
  method = "ranger",
  trControl = cv_control,
  metric = "ROC",
  preProcess = c("center", "scale"),
  num.trees = 500,
  importance = "impurity",
  tuneGrid = tune_grid      # Use our custom grid
)

# Print model summary and variable importance
print(rf_model)
rf_var_imp <- varImp(rf_model)
print(rf_var_imp)
plot(rf_var_imp, main = "Variable Importance - Random Forest Model")
test_data$rf_pred_prob <- predict(rf_model, newdata = test_data, type = "prob")[[2]]

# If not already done, convert the outcome to numeric for calibration and Brier score (Yes = 1, No = 0)
test_data$overtime_numeric <- ifelse(test_data$overtime == "Yes", 1, 0)

# Create bins for the predicted probabilities
test_data$prob_bin_rf <- cut(
  test_data$rf_pred_prob,
  breaks = seq(
    round(quantile(test_data$rf_pred_prob, .005), 2),
    round(quantile(test_data$rf_pred_prob, .995), 2),
    by = 0.005
  ),
  include.lowest = TRUE
)

# Compute observed frequency for each bin
calibration_df_rf <- test_data %>%
  group_by(prob_bin_rf) %>%
  summarize(
    mean_pred = mean(rf_pred_prob),
    obs_freq = mean(overtime_numeric),
    count = n()
  )

# Plot the calibration curve
ggplot(calibration_df_rf, aes(x = mean_pred, y = obs_freq)) +
  geom_point(size = 3) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Mean Predicted Probability", y = "Observed Frequency",
       title = "Calibration Plot - Random Forest Model") +
  theme_minimal()

# Calculate the Brier score for the random forest model
brier_score_rf <- mean((test_data$rf_pred_prob - test_data$overtime_numeric)^2)
cat("Random Forest Brier Score:", brier_score_rf, "\n")

# For comparison, calculate the baseline (climatological) Brier score
basic_brier_score <- mean((mean(test_data$overtime_numeric) - test_data$overtime_numeric)^2)
brier_skill_score_rf <- 1 - (brier_score_rf / basic_brier_score)
cat("Random Forest Brier Skill Score:", brier_skill_score_rf, "\n")

# Calculate and plot the ROC curve and AUC for the random forest model
roc_rf <- roc(test_data$overtime, test_data$rf_pred_prob)
auc_rf <- auc(roc_rf)
cat("Random Forest ROC AUC:", auc_rf, "\n")
plot(roc_rf, main = paste("ROC Curve - Random Forest (AUC =", round(auc_rf, 2), ")"))



