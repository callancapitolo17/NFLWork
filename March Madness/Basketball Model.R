# ============================
# NCAA Tournament Simulator — Pre-Tournament (Full Bracket)
# ============================
# Simulates the entire tournament from scratch using Monte Carlo.
# Run BEFORE the tournament starts (or to model the full bracket ignoring results).

source("shared.R")
source("espn_bracket.R")

# --- 1. Load bracket + ratings ---
bracket_result <- fetch_espn_bracket()
teams_std <- get_teams_std()

# Use all 68 teams (ignoring elimination status from actual results)
final_bracket <- bracket_result$bracket %>% select(team, seed, region, play_in)
bracket_with_ratings <- fetch_bracket_with_ratings(final_bracket, teams_std)

# --- 2. Monte Carlo Simulation (vectorized) ---
# Pre-resolve First Four once (not 10K times)
games_df <- bracket_result$games
bracket_64 <- as.data.frame(resolve_first_four(bracket_with_ratings, games_df))
region_order <- get_region_order(bracket_64)
cat(sprintf("Bracket after First Four: %d teams\n", nrow(bracket_64)))

n_simulations <- 10000
cat(sprintf("Running %s simulations...\n", format(n_simulations, big.mark = ",")))
t0 <- Sys.time()
sim_results <- map_dfr(1:n_simulations, ~ simulate_tournament_fast(bracket_64, region_order))
cat(sprintf("Done in %.0f seconds.\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# --- 4. Results ---
team_results <- sim_results %>%
  group_by(team, seed) %>%
  summarise(
    `Round 32`   = mean(Round_32),
    `Sweet 16`   = mean(Sweet_16),
    `Elite 8`    = mean(Elite_8),
    `Final 4`    = mean(Final_4),
    `Title Game` = mean(Title_Game),
    Champion     = mean(Champion),
    .groups = "drop"
  ) %>%
  arrange(desc(Champion))

# Display as American odds
team_results_odds <- team_results %>%
  mutate(across(-c(team, seed), prob_to_american)) %>%
  mutate(across(where(is.numeric), ~ round(., 0)))

bets <- team_results_odds %>%
  gt() %>%
  tab_header(title = sprintf("March Madness Betting Guide (%s Simulations)", format(n_simulations, big.mark = ",")))
bets
