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
win_total_data <- read.csv("Clean Win Totals - Sheet1.csv")

# Function to convert American odds to implied probability
american_odds_to_prob <- function(odds) {
  return(ifelse(odds > 0,
         100 / (odds + 100),  # For positive odds
         -odds / (-odds + 100)))
}

# Function to de-vig American odds
devig_american_odds <- function(over_odds, under_odds) {
  # Convert odds to implied probabilities
  over_prob = american_odds_to_prob(over_odds)
  under_prob = american_odds_to_prob(under_odds)
  
  # Calculate the total implied probability
  total_prob = over_prob + under_prob
  
  # Remove the vig by normalizing probabilities
  over_prob_no_vig = over_prob / total_prob
  under_prob_no_vig = under_prob / total_prob
  
  # Return de-vig probabilities
  return(list(over_prob_no_vig = over_prob_no_vig, under_prob_no_vig = under_prob_no_vig))
}


FD_nv_odds <- devig_american_odds(win_total_data$Fanduel.Over,win_total_data$Fanduel.Under)
DK_nv_odds <- devig_american_odds(win_total_data$Draftkings.Over,win_total_data$Draftkings.Under)

adjusted_win <- win_total_data %>% 
  mutate(FD_adjusted_total = (Fanduel.Total+0.5)*(FD_nv_odds$over_prob_no_vig) + (Fanduel.Total-0.5)*(FD_nv_odds$under_prob_no_vig),
         DK_adjusted_total = (DraftKings.Total+0.5)*(DK_nv_odds$over_prob_no_vig) + (DraftKings.Total-0.5)*(DK_nv_odds$under_prob_no_vig),
         )
fd_win_data_list <- list(
  N = nrow(win_total_data),
  total = win_total_data$Fanduel.Total,
  prob_over = FD_nv_odds$over_prob_no_vig,
  prob_under = FD_nv_odds$under_prob_no_vig
)

dk_win_data_list <- list(
  N = nrow(win_total_data),
  total = win_total_data$DraftKings.Total,
  prob_over = DK_nv_odds$over_prob_no_vig,
  prob_under = DK_nv_odds$under_prob_no_vig
)
stan_wins <- "data {
  int<lower=0> N;                     // Number of teams
  real<lower=0, upper=17> total[N];             // Listed win totals from sportsbook
  real<lower=0, upper=1> prob_over[N];         // Implied probability for over
  real<lower=0, upper=1> prob_under[N];        // Implied probability for under
}

parameters {
  real<lower=0, upper=17> team_wins[N];         // Latent win totals for each team
  real<lower=0> alpha;                // Dispersion parameter (now inferred)
}

model {
  // Priors for alpha (e.g., weakly informative prior)
  alpha ~ normal(2, 0.5);  // Assume alpha ~ N(2, 1) but constrain to positive
  
  // Priors for win totals shaped by sportsbook probabilities
  for (i in 1:N) {
    team_wins[i] ~ normal(total[i] + (prob_over[i] - prob_under[i]), 2);
  }
}

generated quantities {
  real simulated_wins[N];     // Simulated win totals
  for (i in 1:N) {
    simulated_wins[i] = normal_rng(team_wins[i], alpha);
  }
}
"
fd_win_fit <- stan(model_code = stan_wins, data = fd_win_data_list,
            iter = 2000, warmup = 500, chains = 4, seed = 123) 
print(fd_win_fit)

dk_win_fit <- stan(model_code = stan_wins, data = dk_win_data_list,
                   iter = 2000, warmup = 500, chains = 4, seed = 123)

traceplot(dk_win_fit, pars = c("team_wins[1]", "team_wins[2]"))
print(dk_win_fit)


# Combine results with the original data
results <- betting_data %>%
  mutate(expected_wins = expected_wins)

# Print the results
print(results)

# Normalize win totals
total_wins <- sum(win_totals)
normalized_win_totals <- win_totals / total_wins * 272 #<- total wins
win_probabilities <- normalized_win_totals / 17 #<- total games in Nfl season

# Function to convert win probability to starting ELO rating
elo_from_win_prob <- function(win_prob, avg_opponent_elo = 1500, home_field_adv = 73) { #Where does 73 come from?
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
print(fit) #mean = estimate, se_mean = uncertainty of estimate, sd = variability in posterior, n_eff = number of independent samples, rhat = diagnose convergence, 1 means chains have converged

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
home_team <- "KC"
away_team <- "BUF"
win_prob_home <- predict_match_outcome(home_team, away_team)
print(paste("Win probability for", home_team, "against", away_team, ":", win_prob_home))

# Convert win probability to decimal odds
exp_home_odds <- 1/win_prob_home

# Convert win probability to spread | NBA value = 0.16 | NFL value  = 0.143
exp_home_spread <- (log((1/win_prob_home)-1)/0.143)
exp_home_spread