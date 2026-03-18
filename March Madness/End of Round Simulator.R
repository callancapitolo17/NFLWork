# ============================
# NCAA Tournament Simulator — End of Round (Between Rounds)
# ============================
# Simulates remaining tournament from the end of a completed round.
# Unlike the Dynamic Simulator, assumes the current round is fully complete.
# Includes Survivor Value Rank for survivor pool picks.

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

# Resolve First Four (use actual results where available, simulate the rest)
bracket_64 <- resolve_first_four(bracket_with_ratings, games_played)

# Then filter out any additional losers from later rounds
current_bracket <- bracket_64 %>%
  filter(!team %in% {
    games_played %>%
      filter(status == "final", round != "First Four", !is.na(winner)) %>%
      mutate(loser = ifelse(winner == team1, team2, team1)) %>%
      pull(loser)
  })

cat(sprintf("Tournament: %s | Round: %s | Teams remaining: %d\n",
            tourney_state$state, tourney_state$current_round %||% "N/A", nrow(current_bracket)))

# --- 2. Simulation Functions ---

simulate_round <- function(teams, game_number = 1, region_order_auto = NULL) {
  teams <- get_bracket_matchups(teams, region_order_auto)
  team_pairs <- split(teams, rep(1:(nrow(teams)/2), each = 2))
  winners <- map(team_pairs, function(pair) {
    if (nrow(pair) < 2) return(pair[1, ])
    simulate_game(pair[1, ], pair[2, ], game_number)$winner
  })
  bind_rows(winners)
}

simulate_remaining_tournament <- function(bracket, region_order_auto = NULL) {
  if (is.null(region_order_auto)) region_order_auto <- get_region_order(bracket)

  round_names <- get_remaining_rounds(nrow(bracket))
  team_progress <- bracket %>% select(team, seed)
  for (r in round_names) team_progress[[r]] <- 0

  round_num <- 1
  teams_round <- bracket
  while (nrow(teams_round) > 1) {
    teams_round <- simulate_round(teams_round, game_number = round_num, region_order_auto = region_order_auto)
    if (round_num <= length(round_names)) {
      rn <- round_names[round_num]
      team_progress <- team_progress %>%
        mutate("{rn}" := ifelse(team %in% teams_round$team, 1, .data[[rn]]))
    }
    round_num <- round_num + 1
  }
  team_progress
}

# --- 3. Monte Carlo Simulation (parallelized in batches) ---
library(furrr)
n_workers <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = n_workers)
n_simulations <- 10000
batch_size <- ceiling(n_simulations / n_workers)
sim_results <- future_map_dfr(1:n_workers, function(w) {
  map_dfr(1:batch_size, ~ simulate_remaining_tournament(current_bracket))
}, .options = furrr_options(seed = TRUE))

# --- 4. Results with Survivor Value Rank ---
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

# Display with Survivor Value Rank
team_results_display <- team_results %>%
  select(-p_current, -f_future) %>%
  arrange(desc(Champion)) %>%
  mutate(`Survivor Value Rank` = min_rank(desc(`Survivor Value`))) %>%
  mutate(across(-c(team, seed, `Survivor Value`, `Survivor Value Rank`), prob_to_american)) %>%
  select(-`Survivor Value`)

bets <- team_results_display %>%
  gt() %>%
  tab_header(title = "March Madness — End of Round with Survivor Recommendations") %>%
  fmt_number(columns = where(is.numeric) & !matches("Rank"), decimals = 0) %>%
  gtExtras::gt_hulk_col_numeric(columns = `Survivor Value Rank`)
bets
