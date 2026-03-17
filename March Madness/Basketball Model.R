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

# --- 2. Simulation Functions ---

simulate_round <- function(teams, game_number = 1, region_order_auto = NULL) {
  teams <- get_bracket_matchups(teams, region_order_auto)
  team_pairs <- split(teams, rep(1:(nrow(teams) / 2), each = 2))
  winners <- map(team_pairs, ~ {
    if (nrow(.x) < 2) return(.x[1, ])
    simulate_game(.x[1, ], .x[2, ], game_number)$winner
  })
  bind_rows(winners)
}

simulate_tournament <- function(bracket) {
  round_num <- 1
  teams_round <- bracket
  region_order_auto <- get_region_order(bracket)

  team_progress <- teams_round %>%
    select(team, seed) %>%
    mutate(Round_32 = 0, Sweet_16 = 0, Elite_8 = 0, Final_4 = 0, Title_Game = 0, Champion = 0)

  while (nrow(teams_round) > 1) {
    teams_round <- simulate_round(teams_round, game_number = round_num, region_order_auto = region_order_auto)
    team_progress <- team_progress %>%
      mutate(
        Round_32   = ifelse(nrow(teams_round) < 64 & team %in% teams_round$team, 1, Round_32),
        Sweet_16   = ifelse(nrow(teams_round) < 32 & team %in% teams_round$team, 1, Sweet_16),
        Elite_8    = ifelse(nrow(teams_round) < 16 & team %in% teams_round$team, 1, Elite_8),
        Final_4    = ifelse(nrow(teams_round) <  8 & team %in% teams_round$team, 1, Final_4),
        Title_Game = ifelse(nrow(teams_round) <  4 & team %in% teams_round$team, 1, Title_Game),
        Champion   = ifelse(nrow(teams_round) == 1 & team %in% teams_round$team, 1, Champion)
      )
    round_num <- round_num + 1
  }
  team_progress
}

# --- 3. Monte Carlo Simulation ---
n_simulations <- 10000
sim_results <- map_dfr(1:n_simulations, ~ simulate_tournament(bracket_with_ratings))

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
