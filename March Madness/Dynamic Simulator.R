# ============================
# NCAA Tournament Simulator — Dynamic (Mid-Tournament)
# ============================
# Simulates remaining games from the current tournament state.
# Handles partially completed rounds where some teams have already advanced.
# Includes Survivor Value metric for survivor pool strategy.

source("shared.R")
source("espn_bracket.R")

# --- 1. Load bracket + ratings + current state ---
bracket_result <- fetch_espn_bracket()
teams_std <- get_teams_std()

final_bracket <- bracket_result$bracket %>% select(team, seed, region, play_in)
bracket_with_ratings <- fetch_bracket_with_ratings(final_bracket, teams_std)

# Build current_bracket: only surviving teams
tourney_state <- bracket_result$tournament_state
games_played  <- bracket_result$games

current_bracket <- bracket_with_ratings %>%
  filter(!team %in% {
    games_played %>%
      filter(status == "final", !is.na(winner)) %>%
      mutate(loser = ifelse(winner == team1, team2, team1)) %>%
      pull(loser)
  }) %>%
  mutate(status = "pending")

cat(sprintf("Tournament: %s | Round: %s | Teams remaining: %d\n",
            tourney_state$state, tourney_state$current_round %||% "N/A", nrow(current_bracket)))

# --- 2. Simulation Functions ---

# Dynamic round: if one team already advanced, carry it forward
simulate_round_dynamic <- function(teams, game_number = 1, region_order_auto = NULL) {
  teams_ordered <- get_bracket_matchups(teams, region_order_auto)
  team_pairs <- split(teams_ordered, rep(1:(nrow(teams_ordered)/2), each = 2))

  winners <- map(team_pairs, function(pair) {
    if (all(pair$status == "advanced")) {
      pair[1, ]
    } else if (any(pair$status == "advanced") && any(pair$status == "pending")) {
      pair %>% filter(status == "advanced") %>% slice(1)
    } else {
      simulate_game(pair[1, ], pair[2, ], game_number)$winner
    }
  })
  bind_rows(winners)
}

simulate_remaining_tournament_dynamic <- function(bracket, region_order_auto = NULL) {
  if (is.null(region_order_auto)) region_order_auto <- get_region_order(bracket)

  round_names <- get_remaining_rounds(nrow(bracket))
  team_progress <- bracket %>% select(team, seed)
  for (r in round_names) team_progress[[r]] <- 0

  round_num <- 1
  teams_round <- bracket
  while (nrow(teams_round) > 1) {
    teams_round <- simulate_round_dynamic(teams_round, game_number = round_num, region_order_auto = region_order_auto)
    if (round_num <= length(round_names)) {
      rn <- round_names[round_num]
      team_progress <- team_progress %>%
        mutate("{rn}" := ifelse(team %in% teams_round$team, 1, .data[[rn]]))
    }
    round_num <- round_num + 1
  }
  team_progress
}

# --- 3. Monte Carlo Simulation ---
set.seed(12)
n_simulations <- 10000
sim_results <- map_dfr(1:n_simulations, ~ simulate_remaining_tournament_dynamic(current_bracket))

# --- 4. Results with Survivor Value ---
team_results <- sim_results %>%
  group_by(team, seed) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  {
    sim_rounds <- setdiff(names(.), c("team", "seed"))
    if (length(sim_rounds) >= 1) {
      p_current_col <- sim_rounds[1]
      future_rounds <- sim_rounds[-1]
      mutate(.,
        p_current = .data[[p_current_col]],
        f_future = if (length(future_rounds) > 0) rowMeans(select(., all_of(future_rounds))) else 0,
        `Survivor Value` = p_current * (1 - f_future)
      )
    } else .
  } %>%
  arrange(desc(`Survivor Value`))

# Display
team_results_display <- team_results %>%
  select(-p_current, -f_future) %>%
  arrange(desc(Champion)) %>%
  mutate(across(-c(team, seed, `Survivor Value`), prob_to_american))

bets <- team_results_display %>%
  gt() %>%
  tab_header(title = "March Madness — Dynamic Simulator with Survivor Values") %>%
  fmt_number(columns = where(is.numeric) & !matches("Survivor"), decimals = 0) %>%
  fmt_number(columns = `Survivor Value`, decimals = 3)
bets
