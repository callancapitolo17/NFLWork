# =============================================================================
# March Madness Props Analysis
# =============================================================================
# Derives prop bet probabilities from simulation output and compares to
# Kalshi/sportsbook markets. Based on methodology from:
# https://unabated.com/articles/solving-march-madness-prop-bets
#
# Props covered:
#   - Upset props (will a X-seed beat a Y-seed?)
#   - Seed advancement (how many 1-seeds in F4?)
#   - Conference win totals
#   - Cumulative seed sum props
#   - Team advancement (make Sweet 16, Elite 8, etc.)


source("shared.R")
source("espn_bracket.R")
library(hoopR)

# =============================================================================
# 1. Run simulation and keep raw per-sim data
# =============================================================================

cat("Loading bracket and ratings...\n")
br <- fetch_espn_bracket()
ts <- get_teams_std()
fb <- br$bracket %>% select(team, seed, region, play_in)
bwr <- fetch_bracket_with_ratings(fb, ts)
bracket_64 <- as.data.frame(resolve_first_four(bwr, br$games))
region_order <- get_region_order(bracket_64)

n_sims <- 10000
cat(sprintf("Running %s simulations...\n", format(n_sims, big.mark = ",")))
t0 <- Sys.time()
raw_sims <- map_dfr(1:n_sims, function(i) {
  result <- simulate_tournament_fast(bracket_64, region_order)
  result$sim_id <- i
  result
})
cat(sprintf("Done in %.1f seconds.\n\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# Add conference info (from cbbdata BPI which has conf column)
conf_map <- tryCatch({
  cbd_bpi_ratings() %>%
    transmute(team_short = team, conference = conf) %>%
    mutate(team = map_chr(team_short, ~ get_standard_team(.x, teams_std = ts))) %>%
    select(team, conference) %>%
    distinct(team, .keep_all = TRUE)
}, error = function(e) NULL)
if (!is.null(conf_map)) {
  raw_sims <- raw_sims %>% left_join(conf_map, by = "team")
}

# =============================================================================
# 2. Prop Calculations
# =============================================================================

prob_to_american <- function(prob) {
  ifelse(prob >= 0.5, -100 * prob / (1 - prob), 100 / prob - 100)
}

cat("=== CHAMPIONSHIP ODDS ===\n\n")
champ <- raw_sims %>%
  group_by(team, seed) %>%
  summarise(prob = mean(Champion), .groups = "drop") %>%
  arrange(desc(prob)) %>%
  head(15)
cat(sprintf("%-30s %4s %7s %8s\n", "Team", "Seed", "Prob", "Odds"))
for (i in 1:nrow(champ)) {
  r <- champ[i, ]
  cat(sprintf("%-30s %4d %6.2f%% %+8.0f\n", r$team, r$seed, r$prob * 100, prob_to_american(r$prob)))
}

# --- Upset Props ---
cat("\n=== UPSET PROPS (Round of 64) ===\n\n")

# For each seed matchup, probability of upset
seed_matchups <- list(
  c(1, 16), c(2, 15), c(3, 14), c(4, 13),
  c(5, 12), c(6, 11), c(7, 10), c(8, 9)
)

cat(sprintf("%-12s %8s %8s %8s\n", "Matchup", "P(Fav)", "P(Upset)", "P(1+ Upset)"))
for (m in seed_matchups) {
  fav_seed <- m[1]
  dog_seed <- m[2]

  # Each sim has 4 games of this type (one per region)
  # A favorite survives R32 if they advance
  sim_by_sim <- raw_sims %>%
    filter(seed == fav_seed) %>%
    group_by(sim_id) %>%
    summarise(fav_wins = sum(Round_32), .groups = "drop")

  # P(individual favorite wins) = mean across all 4 games
  p_fav_individual <- mean(raw_sims$Round_32[raw_sims$seed == fav_seed])
  p_upset_individual <- 1 - p_fav_individual

  # P(at least one upset in 4 games)
  p_all_fav_win <- mean(sim_by_sim$fav_wins == 4)
  p_at_least_one_upset <- 1 - p_all_fav_win

  cat(sprintf("#%d vs #%d     %7.1f%% %7.1f%% %7.1f%%\n",
              fav_seed, dog_seed, p_fav_individual * 100, p_upset_individual * 100, p_at_least_one_upset * 100))
}

# --- Specific upset counts ---
cat("\n=== 12-SEED UPSET COUNT (how many 12-seeds beat 5-seeds?) ===\n\n")
twelve_wins <- raw_sims %>%
  filter(seed == 12) %>%
  group_by(sim_id) %>%
  summarise(n_wins = sum(Round_32), .groups = "drop")

for (k in 0:4) {
  p <- mean(twelve_wins$n_wins == k)
  cat(sprintf("Exactly %d: %6.2f%% (%+.0f)\n", k, p * 100, prob_to_american(p)))
}
cat(sprintf("At least 1: %6.2f%% (%+.0f)\n",
            mean(twelve_wins$n_wins >= 1) * 100, prob_to_american(mean(twelve_wins$n_wins >= 1))))

# --- 1-Seed Props ---
cat("\n=== 1-SEED ADVANCEMENT ===\n\n")
one_seeds <- raw_sims %>% filter(seed == 1)
one_seed_f4 <- one_seeds %>%
  group_by(sim_id) %>%
  summarise(in_f4 = sum(Final_4), champ = sum(Champion), .groups = "drop")

cat("1-seeds in Final Four:\n")
for (k in 0:4) {
  p <- mean(one_seed_f4$in_f4 == k)
  cat(sprintf("  Exactly %d: %6.2f%% (%+.0f)\n", k, p * 100, prob_to_american(p)))
}

cat(sprintf("\nA 1-seed wins championship: %6.2f%% (%+.0f)\n",
            mean(one_seed_f4$champ >= 1) * 100, prob_to_american(mean(one_seed_f4$champ >= 1))))

# --- Seed Sum Props ---
cat("\n=== FINAL FOUR SEED SUM ===\n\n")
f4_seeds <- raw_sims %>%
  filter(Final_4 == 1) %>%
  group_by(sim_id) %>%
  summarise(seed_sum = sum(seed), lowest_seed = max(seed), .groups = "drop")

cat(sprintf("Mean seed sum: %.1f\n", mean(f4_seeds$seed_sum)))
cat(sprintf("Median seed sum: %.0f\n", median(f4_seeds$seed_sum)))

# Common over/under lines
for (line in c(8, 10, 12, 15)) {
  p_over <- mean(f4_seeds$seed_sum > line)
  cat(sprintf("Over %d: %6.2f%% (%+.0f)\n", line, p_over * 100, prob_to_american(p_over)))
}

cat("\nLowest seed in Final Four:\n")
for (s in sort(unique(f4_seeds$lowest_seed))) {
  p <- mean(f4_seeds$lowest_seed == s)
  if (p > 0.005) cat(sprintf("  %2d-seed: %6.2f%%\n", s, p * 100))
}

# --- Conference Props ---
if ("conference" %in% names(raw_sims)) {
  cat("\n=== CONFERENCE WIN TOTALS ===\n\n")
  conf_wins <- raw_sims %>%
    filter(!is.na(conference)) %>%
    group_by(sim_id, conference) %>%
    summarise(
      total_wins = sum(Round_32 + Sweet_16 + Elite_8 + Final_4 + Title_Game + Champion),
      champ = max(Champion),
      .groups = "drop"
    )

  conf_summary <- conf_wins %>%
    group_by(conference) %>%
    summarise(
      mean_wins = mean(total_wins),
      p_champ = mean(champ >= 1),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_wins))

  cat(sprintf("%-8s %8s %10s\n", "Conf", "Avg Wins", "P(Champ)"))
  for (i in 1:min(15, nrow(conf_summary))) {
    r <- conf_summary[i, ]
    cat(sprintf("%-8s %8.1f %9.1f%%\n", r$conference, r$mean_wins, r$p_champ * 100))
  }
}

# --- Team Advancement Props ---
cat("\n=== TEAM ADVANCEMENT PROPS (Sweet 16, Elite 8, Final 4) ===\n\n")
team_adv <- raw_sims %>%
  group_by(team, seed) %>%
  summarise(
    R32 = mean(Round_32), S16 = mean(Sweet_16), E8 = mean(Elite_8),
    F4 = mean(Final_4), TG = mean(Title_Game), Champ = mean(Champion),
    .groups = "drop"
  ) %>%
  arrange(desc(Champ))

cat(sprintf("%-30s %4s %6s %6s %6s %6s %6s %7s\n",
            "Team", "Seed", "R32", "S16", "E8", "F4", "Title", "Champ"))
cat(paste(rep("-", 95), collapse = ""), "\n")
for (i in 1:min(25, nrow(team_adv))) {
  r <- team_adv[i, ]
  cat(sprintf("%-30s %4d %5.1f%% %5.1f%% %5.1f%% %5.1f%% %5.1f%% %6.2f%%\n",
              r$team, r$seed, r$R32*100, r$S16*100, r$E8*100, r$F4*100, r$TG*100, r$Champ*100))
}

cat("\n=== DONE ===\n")
