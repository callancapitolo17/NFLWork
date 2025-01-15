### Bayesian Sports Models in R: Chapter 13 - Bayesian ELO| Andrew Mack | @Gingfacekillah

# Load libraries | Install required packages prior to loading!
library(ggplot2)        # ggplot: plotting functions
library(ggridges)       # density ridges plot add-on
library(viridis)        # viridis color palette for plots
library(bayesplot)      # Plot mcmc results
library(rstan)          # R interface for Stan programming language
library(tidyverse)      # data wrangling functions
library(lubridate)      # time & date functions
library(nflfastR)
library(nflreadr)
all <- load_schedules(1999:2024) %>% 
  filter(!is.na(away_score)) %>% 
  filter(location == "Home") %>% 
  rename(away_team_score = away_score, home_team_score = home_score) %>%
  mutate(
    outcome = ifelse(home_team_score == away_team_score,0.5, ifelse(home_team_score > away_team_score, 1, 0)))

home_win_rate <- mean(all$outcome)  # Outcome: 1 if home wins, 0 otherwise
hfa_points <- 400 * log10(home_win_rate / (1 - home_win_rate))

#### 1. Load the data ####
matches <- load_schedules(2024) %>% 
  filter(!is.na(away_score)) %>% 
  select(date = gameday,away_team,away_score,home_team,home_score,location)
head(matches)

#### 2. Data wrangling ####
matches <- matches %>%
  rename(away_team_score = away_score, home_team_score = home_score) %>%
  mutate(
    outcome = ifelse(home_team_score == away_team_score,0.5, ifelse(home_team_score > away_team_score, 1, 0)))

# Print data frame
print(head(matches))

# Create unique IDs for each team
teams <- unique(c(matches$home_team, matches$away_team))
team_ids <- setNames(1:length(teams), teams)

# Map team names to their corresponding IDs
matches <- matches %>%
  mutate(
    home_team_id = team_ids[home_team],
    away_team_id = team_ids[away_team])

#### 3. Using win totals to estimate starting ELO ratings ####
# Gather win totals from reasonably sharp sportsbook and enter them into this vector:, how should I account for different odds?
win_totals <- c(
  "ARI" = 6.5,
  "ATL" = 9.5,
  "BAL" = 11.5,
  "BUF" = 10.5,
  "CAR" = 5.5,
  "CHI" = 8.5,
  "CIN" = 10.5,
  "CLE" = 9.5,
  "DAL" = 10.5,
  "DEN" = 5.5,
  "DET" = 10.5,
  "GB" = 9.5,
  "HOU" = 9.5,
  "IND" = 8.5,
  "JAX" = 8.5,
  "KC" = 11.5,
  "LV" = 6.5,
  "LAC" = 8.5,
  "LAR" = 8.5,
  "MIA" = 9.5,
  "MIN" = 6.5,
  "NE" = 5.5,
  "NO" = 7.5,
  "NYG" = 6.5,
  "NYJ" = 9.5,
  "PHI" = 10.5,
  "PIT" = 7.5,
  "SF" = 11.5,
  "SEA" = 7.5,
  "TB" = 7.5,
  "TEN" = 6.5,
  "WAS" = 6.5
)

# Normalize win totals
total_wins <- sum(win_totals)
normalized_win_totals <- win_totals / total_wins * 265.5 #<- total wins
win_probabilities <- normalized_win_totals / 17 #<- total games in NBA season

# Function to convert win probability to starting ELO rating
elo_from_win_prob <- function(win_prob, avg_opponent_elo = 1500, home_field_adv = hfa_points) { #Where does 73 come from?
  team_elo <- avg_opponent_elo + home_field_adv - 400 * log10((1 / win_prob) - 1)
  return(team_elo)
}

# Calculate starting ELO ratings for each team
starting_elo_ratings <- sapply(win_probabilities, elo_from_win_prob)
team_starting_elo <- data.frame(
  Team = names(starting_elo_ratings),
  Starting_Elo_Rating = starting_elo_ratings)

# Print ELO ratings
print(team_starting_elo)

# Prepare the starting ELO ratings to match the team IDs
starting_elo <- team_starting_elo$Starting_Elo_Rating[match(names(team_ids), team_starting_elo$Team)]

# Prepare data for Stan (using win total priors)
data_list <- list(
  N = nrow(matches),
  T = length(teams),
  home_team = matches$home_team_id,
  away_team = matches$away_team_id,
  outcome = matches$outcome,
  starting_elo = starting_elo)

# Updated Stan code (using win total priors)
stan_model_code_new <- "
data {
    int<lower=0> N;                     // Number of matches
    int<lower=0> T;                     // Number of teams
    int<lower=1, upper=T> home_team[N]; // Home team for each match
    int<lower=1, upper=T> away_team[N]; // Away team for each match
    int<lower=0, upper=1> outcome[N];   // Outcome of each match
    vector[T] starting_elo;             // Starting ELO ratings
}
parameters {
    vector[T] rating;                   // Ratings for each team
    real<lower=0> home_adv;             // Home field advantage)
  real<lower=0, upper=100> K;           // ELO sensitivity constant
}
model {
  // Priors
  rating ~ normal(starting_elo, 50);  // Team rating prior using starting ELO ratings
  home_adv ~ normal(50, 25);          // Home advantage prior
  K ~ normal(50, 25);                 // K parameter prior

  for (n in 1:N) {
    // Likelihood
    real prob = 1 / (1 + pow(10, (rating[away_team[n]] - (rating[home_team[n]] + home_adv)) / 400));
    outcome[n] ~ bernoulli(prob);
  }
}
generated quantities {
  vector[T] new_rating = rating;
  vector[N] rating_change;

  for (n in 1:N) {
    // Calculate the probability of home win with ELo formula
    real prob = 1 / (1 + pow(10, (new_rating[away_team[n]] - (new_rating[home_team[n]] + home_adv)) / 400));

    // Calculate the rating change
    rating_change[n] = K * (outcome[n] - prob);

    // Apply the rating changes incrementally
    new_rating[home_team[n]] += rating_change[n];
    new_rating[away_team[n]] -= rating_change[n];
  }
}
"

# fit Stan model with mcmc
fit <- stan(model_code = stan_model_code_new, data = data_list,
            iter = 2000, warmup = 500, chains = 4, seed = 123)

# Print Stan fit
print(fit)

# Print selected parameter trace plots
traceplot(fit, pars = c("K", "home_adv", "rating[1]", "rating[2]"))

# Print team ELO ratings + home_adv & K parameter estimates
fit_summary <- summary(fit, pars = c("rating", "home_adv", "K"))$summary
rating_estimates <- fit_summary[, "mean"]
team_ratings <- data.frame(
  Team = names(c(team_ids, "home_adv", "K")),
  ELO_Rating = rating_estimates)

# Print ELO ratings
print(team_ratings)


# Print Stan fit
print(fit)

# Print selected parameter trace plots
traceplot(fit, pars = c("K", "home_adv", "rating[1]", "rating[2]"))

# Print team ELO ratings + home_adv & K parameter estimates
fit_summary <- summary(fit, pars = c("rating", "home_adv", "K"))$summary
rating_estimates <- fit_summary[, "mean"]
team_ratings <- data.frame(
  Team = names(c(team_ids, "home_adv", "K")),
  ELO_Rating = rating_estimates)

# Print ELO ratings
print(team_ratings)


#### 4. Plot estimated posterior team strength ####
# Extract parameters from the fitted model
fit_summary <- summary(fit, pars = "rating")$summary
rating_estimates <- fit_summary[, "mean"]
team_ratings <- data.frame(
  Team = names(team_ids),
  ELO_Rating = rating_estimates)
team_ratings <- team_ratings %>% arrange(desc(ELO_Rating))
posterior <- rstan::extract(fit)
rating_samples <- as.data.frame(posterior$rating)
colnames(rating_samples) <- teams
rating_samples_long <- tidyr::gather(rating_samples, key = "team", value = "rating")
rating_samples_long$team <- factor(rating_samples_long$team, levels = rev(team_ratings$Team))

# Plot team ELO ratings highest to lowest
ggplot(rating_samples_long, aes(x = rating, y = team, fill = after_stat(x))) +
  geom_density_ridges_gradient() +
  scale_fill_viridis_c(name = "Rating", option = "C") + # G or D also good color schemes
  labs(title = "Estimated Team Ratings",
       x = "ELO Posterior Rating",
       y = "Team") +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_vline(xintercept = 1500, linetype = "dashed", color = viridis::viridis(1))+
  geom_vline(xintercept = 1300, linetype = "dashed", color = viridis::viridis(1))+
  geom_vline(xintercept = 1700, linetype = "dashed", color = viridis::viridis(1))

# Calculate maximum a posteriori (MAP) for each team's rating
team_map <- rating_samples_long %>%
  group_by(team) %>%
  summarize(map = mean(rating))

# Print MAP
team_map


#### 5. Simulating future games ####
predict_match_outcome <- function(home_team_name, away_team_name, n_simulations = 10000) {
  home_id <- team_ids[home_team_name]
  away_id <- team_ids[away_team_name]
  
  sim_outcomes <- replicate(n_simulations, {
    rating_home <- sample(rating_samples[, home_id], 1)
    rating_away <- sample(rating_samples[, away_id], 1)
    home_adv <- sample(posterior$home_adv, 1)
    
    win_prob <- 1 / (1 + 10^((rating_away - (rating_home + home_adv)) / 400))
    rbinom(1, 1, win_prob)
  })
  
  win_prob_home <- mean(sim_outcomes)
  return(win_prob_home)
}

# Example matchup prediction: Detroit at Boston
home_team <- "Boston Celtics"
away_team <- "Detroit Pistons"
win_prob_home <- predict_match_outcome(home_team, away_team)
print(paste("Win probability for", home_team, "against", away_team, ":", win_prob_home))

# Sanity check: Boston vs Boston - should show us the home advantage
home_team <- "Boston Celtics"
away_team <- "Boston Celtics"
win_prob_home <- predict_match_outcome(home_team, away_team)
print(paste("Win probability for", home_team, "against", away_team, ":", win_prob_home))

# Convert win probability to decimal odds
exp_home_odds <- 1/win_prob_home
exp_home_odds

# Convert win probability to spread | NBA value = 0.16 | NFL value  = 0.143
exp_home_spread <- (log((1/win_prob_home)-1)/0.16)
exp_home_spread