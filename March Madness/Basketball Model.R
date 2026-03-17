# ============================
# NCAA Tournament Simulator
# ============================
# This script:
# 1. Acquires bracket data (team names, seeds, regions, and play-in flag).
# 2. Acquires power ratings from
# 3. Combines these ratings into a composite rating (weighted average).
# 4. Merges the composite ratings with the bracket data.
# 5. Runs a Monte Carlo simulation of the tournament,
#    including simulating the First Four (play-in) games.
#
# Load required libraries:
library(httr)
library(jsonlite)
library(dplyr)
library(purrr)
library(readr)
library(tidyr)
library(googlesheets4)
library(stringr)
library(hoopR)
library(cbbdata)
library(gt)
# gs4_auth()

# --- Configuration ---
# Season year: Nov-Apr academic year (e.g., 2025-26 season = 2026)
current_year <- if (as.integer(format(Sys.Date(), "%m")) >= 11) {
  as.integer(format(Sys.Date(), "%Y")) + 1L
} else {
  as.integer(format(Sys.Date(), "%Y"))
}

# ------------------------------
# 1. Acquire Bracket Data (via ESPN API)
# ------------------------------
source("espn_bracket.R")
bracket_result <- fetch_espn_bracket()
# Basketball Model always simulates the full tournament from scratch.
# Use all 68 teams (ignoring elimination status from actual results).
final_bracket <- bracket_result$bracket %>%
  select(team, seed, region, play_in)

# ------------------------------
# 2. Acquire Power Ratings
# ------------------------------

teams_std <- espn_mbb_teams(current_year) %>%
  mutate(team = ifelse(team == "McNeese", "McNeese State", team))

# Helper for team name matching
clean_text <- function(x) {
  sub("St\\.$", "State", x, ignore.case = TRUE)
}

get_standard_team <- function(team, teams_std) {
  if (is.na(team)) return(NA_character_)
  team_clean <- clean_text(team)
  for (i in 1:nrow(teams_std)) {
    variants <- teams_std[i, c("abbreviation", "display_name", "short_name", "mascot", "nickname", "team")]
    for (variant in variants) {
      if (!is.na(variant) && team_clean == clean_text(variant)) {
        return(teams_std$display_name[i])
      }
    }
  }
  team
}

# a) ESPN BPI (via cbbdata)
cat("Fetching BPI ratings...\n")
clean_bpi_data <- tryCatch({
  cbd_bpi_ratings() %>%
    transmute(team = team, bpi = bpi_value, standard_team = team)
}, error = function(e) {
  cat(sprintf("Warning: BPI fetch failed: %s\n", e$message))
  tibble(team = character(), bpi = numeric(), standard_team = character())
})
cat(sprintf("BPI: %d teams\n", nrow(clean_bpi_data)))

# b) KenPom Ratings (from Google Sheets)
cat("Fetching KenPom ratings...\n")
sheet_url <- "https://docs.google.com/spreadsheets/d/10o9dwZeyREliM8iOIevpHIhNeCHYeVGE-y50N3KHujQ/edit?gid=0#gid=0"
clean_kenpom_data <- tryCatch({
  kenpom_data <- read_sheet(sheet_url) %>%
    mutate(Team = str_replace(Team, "\\s\\d+$", "")) %>%
    as.data.frame() %>%
    mutate(Team = ifelse(Team == "Connecticut", "UConn",
                  ifelse(Team == "Mississippi", "Ole Miss",
                  ifelse(Team == "Nebraska Omaha", "Omaha", Team))))
  kenpom_data %>%
    mutate(standard_team = map_chr(Team, ~ get_standard_team(.x, teams_std = teams_std))) %>%
    select(standard_team, NetRtg) %>%
    mutate(NetRtg = map_dbl(NetRtg, ~ ifelse(is.null(.x), NA_real_, as.numeric(.x))))
}, error = function(e) {
  cat(sprintf("Warning: KenPom fetch failed: %s\n", e$message))
  tibble(standard_team = character(), NetRtg = numeric())
})
cat(sprintf("KenPom: %d teams\n", nrow(clean_kenpom_data)))

# c) Torvik Ratings (via cbbdata, fallback to direct API)
cat("Fetching Torvik ratings...\n")
clean_torvik_data <- tryCatch({
  torvik <- cbd_torvik_ratings(year = current_year)
  if (nrow(torvik) == 0) stop("No data for current year")
  torvik %>%
    mutate(TorvikMargin = ((adj_o - adj_d) / 100) * 70) %>%
    transmute(standard_team = team, TorvikMargin)
}, error = function(e) {
  cat(sprintf("  cbbdata Torvik failed (%s), trying direct API...\n", e$message))
  tryCatch({
    resp <- GET(
      sprintf("https://barttorvik.com/trank.php?year=%d&t=0&json=1", current_year),
      add_headers("User-Agent" = "Mozilla/5.0")
    )
    txt <- content(resp, "text", encoding = "UTF-8")
    if (!startsWith(trimws(txt), "[")) stop("Cloudflare blocked")
    torvik_data <- as.data.frame(fromJSON(txt))
    colnames(torvik_data)[c(1, 2, 3)] <- c("Team", "OffEff", "DefEff")
    torvik_data %>%
      mutate(OffEff = as.numeric(OffEff), DefEff = as.numeric(DefEff),
             TorvikMargin = ((OffEff - DefEff) / 100) * 70,
             standard_team = map_chr(Team, ~ get_standard_team(.x, teams_std = teams_std))) %>%
      select(standard_team, TorvikMargin)
  }, error = function(e2) {
    cat(sprintf("  Warning: Torvik unavailable: %s\n", e2$message))
    tibble(standard_team = character(), TorvikMargin = numeric())
  })
})
cat(sprintf("Torvik: %d teams\n", nrow(clean_torvik_data)))

# d) EvanMiya Ratings (from Google Sheets - fragile, may fail)
cat("Fetching EvanMiya ratings...\n")
clean_evan_miya_wide <- tryCatch({
  evan_miya <- read_sheet("https://docs.google.com/spreadsheets/u/0/d/1zaDMGo6bRMis_VD7yH4uZY9aFXbp0QUPp3HRoE8x3EY/edit?usp=drive_web&pli=1&authuser=0")
  group_size <- 21
  evan_miya$group <- rep(1:(nrow(evan_miya)/group_size), each = group_size)[1:nrow(evan_miya)]
  evan_miya <- evan_miya %>%
    group_by(group) %>% mutate(row_number = row_number()) %>% ungroup() %>% as.data.frame()
  evan_miya_wide <- evan_miya %>%
    pivot_wider(names_from = row_number, values_from = 1) %>%
    select(-group) %>%
    mutate(`2` = `2` %>% str_replace_all("[^A-Za-z0-9 ]", "") %>% str_squish())
  colnames(evan_miya_wide) <- c(
    "Relative Ranking", "Team", "O-Rate", "D-Rate", "Relative Rating",
    "Opponent Adjust", "Pace Adjust", "Off Rank", "Def Rank", "True Tempo",
    "Tempo Rank", "Injury Rank", "Home Rank", "Roster Rank",
    "Kill Shots Per Game", "Kill Shots Conceded Per Game", "Kill Shots Margin Per Game",
    "Total Kill Shots", "Total Kill Shots Conceded", "D1 Wins", "D1 Losses"
  )
  evan_miya_wide %>%
    mutate(Team = ifelse(Team == "Connecticut", "UConn",
                  ifelse(Team == "Mississippi", "Ole Miss", Team))) %>%
    mutate(standard_team = map_chr(Team, ~ get_standard_team(.x, teams_std = teams_std))) %>%
    select(standard_team, `Relative Rating`) %>%
    mutate(`Relative Rating` = map_dbl(`Relative Rating`, ~ ifelse(is.null(.x), NA_real_, as.numeric(.x))))
}, error = function(e) {
  cat(sprintf("Warning: EvanMiya fetch failed: %s\n", e$message))
  tibble(standard_team = character(), `Relative Rating` = numeric())
})
cat(sprintf("EvanMiya: %d teams\n", nrow(clean_evan_miya_wide)))

# Combine all power ratings
power_ratings <- clean_bpi_data %>%
  left_join(clean_kenpom_data %>% rename(KenPomRating = NetRtg), by = "standard_team") %>%
  left_join(clean_torvik_data, by = "standard_team") %>%
  left_join(clean_evan_miya_wide %>% rename(EvanMiyaRating = `Relative Rating`), by = "standard_team") %>%
  mutate(
    EvanMiyaRating = as.numeric(unlist(EvanMiyaRating)),
    KenPomMargin = (as.numeric(KenPomRating) / 100) * 70,
    EvanMiyaMargin = (EvanMiyaRating / 100) * 70,
    BPIMargin = bpi
  ) %>%
  select(standard_team, KenPomMargin, BPIMargin, EvanMiyaMargin, TorvikMargin)

# ------------------------------
# 3. Merge Bracket Data with Power Ratings
# ------------------------------
# ESPN bracket names match BPI names directly (both from ESPN)
# For other sources, use get_standard_team() to resolve
final_bracket <- final_bracket %>%
  mutate(standard_team = map_chr(team, ~ get_standard_team(.x, teams_std = teams_std))) %>%
  mutate(seed = as.numeric(seed))

bracket_with_ratings <- left_join(final_bracket, power_ratings, by = "standard_team") %>%
  rowwise() %>%
  mutate(composite_rating = mean(c_across(KenPomMargin:TorvikMargin), na.rm = TRUE)) %>%
  ungroup()

n_matched <- sum(!is.na(bracket_with_ratings$composite_rating))
cat(sprintf("Bracket: %d teams, %d with ratings\n", nrow(bracket_with_ratings), n_matched))

# ------------------------------
# 4. Simulation Functions
# ------------------------------

# Function to simulate a game between two teams using composite ratings.
# ============================
# NCAA Tournament Simulator
# ============================

# Load required libraries
library(dplyr)
library(purrr)
library(gt)

# ------------------------------
# 1. Helper Functions
# ------------------------------

# Get correct NCAA Tournament matchups
get_bracket_matchups <- function(teams) {
  match_order <- tibble(
    seed = c(1, 16, 8, 9, 5, 12, 4, 13, 6, 11, 3, 14, 7, 10, 2, 15),
    matchup_order = 1:16
  )
  
  teams <- teams %>%
    left_join(match_order, by = "seed") %>%  # Ensure correct pairing order
    arrange(region, matchup_order) %>%       # Order by region and seeding
    select(-matchup_order)                   # Drop helper column
  
  return(teams)
}

# Simulate a single game
simulate_game <- function(team1, team2, game_number = 1, beta1 = 0.1, beta2 = 0.05, sd_margin = 11.2) {
  expected_diff <- team1$composite_rating - team2$composite_rating
  win_prob_team1 <- 1 / (1 + exp(-0.16 * expected_diff))
  actual_margin <- rnorm(1, mean = expected_diff, sd = sd_margin)
  
  winner <- if (actual_margin > 0) team1 else team2
  loser <- if (actual_margin > 0) team2 else team1
  
  rating_change <- beta1 * (actual_margin - expected_diff) +
    beta2 * (actual_margin - expected_diff) * log(game_number + 1)
  
  winner$composite_rating <- winner$composite_rating + rating_change
  loser$composite_rating <- loser$composite_rating - rating_change
  
  return(list(winner = winner))
}

# Simulate a single round of the tournament
simulate_round <- function(teams, game_number = 1) {
  teams <- get_bracket_matchups(teams)
  team_pairs <- split(teams, rep(1:(nrow(teams) / 2), each = 2))
  
  winners <- map(team_pairs, ~ {
    team1 <- .x[1, ]
    team2 <- .x[2, ]
    if (nrow(.x) < 2) return(team1) # Auto-advance if no opponent
    simulate_game(team1, team2, game_number)$winner
  })
  
  return(bind_rows(winners))
}

# Simulate the full NCAA tournament
simulate_tournament <- function(bracket) {
  round_num <- 1
  teams_round <- bracket
  
  team_progress <- teams_round %>%
    select(team, seed) %>%
    mutate(Round_32 = 0, Sweet_16 = 0, Elite_8 = 0, Final_4 = 0, Title_Game = 0, Champion = 0)
  
  while (nrow(teams_round) > 1) {
    teams_round <- simulate_round(teams_round, game_number = round_num)
    
    team_progress <- team_progress %>%
      mutate(
        Round_32 = ifelse(nrow(teams_round) < 64 & team %in% teams_round$team, 1, Round_32),
        Sweet_16 = ifelse(nrow(teams_round) < 32 & team %in% teams_round$team, 1, Sweet_16),
        Elite_8 = ifelse(nrow(teams_round) < 16 & team %in% teams_round$team, 1, Elite_8),
        Final_4 = ifelse(nrow(teams_round) < 8 & team %in% teams_round$team, 1, Final_4),
        Title_Game = ifelse(nrow(teams_round) < 4 & team %in% teams_round$team, 1, Title_Game),
        Champion = ifelse(nrow(teams_round) == 1 & team %in% teams_round$team, 1, Champion)
      )
    
    round_num <- round_num + 1
  }
  
  return(team_progress)
}

# ------------------------------
# 2. Monte Carlo Tournament Simulation
# ------------------------------

n_simulations <- 10000

team_results <- bracket_with_ratings %>%
  select(team, seed) %>%
  mutate(Round_32 = 0, Sweet_16 = 0, Elite_8 = 0, Final_4 = 0, Title_Game = 0, Champion = 0)

sim_results <- map_dfr(1:n_simulations, ~ simulate_tournament(bracket_with_ratings))

prob_to_american <- function(prob) {
  ifelse(prob >= 0.5, 
         -100 * prob / (1 - prob),  # Favorite (negative odds)
         100 / prob - 100)          # Underdog (positive odds)
}
team_results <- sim_results %>%
  group_by(team, seed) %>%
  summarise(
    `Round 32` = prob_to_american(mean(Round_32)),
    `Sweet 16` = prob_to_american(mean(Sweet_16)),
    `Elite 8` = prob_to_american(mean(Elite_8)),
    `Final 4` = prob_to_american(mean(Final_4)),
    `Title Game` = prob_to_american(mean(Title_Game)),
    Champion = prob_to_american(mean(Champion)),
    .groups = "drop"
  ) %>%
  arrange((Champion)) %>% 
  mutate_if(is.numeric, ~round(.,0)) #%>%
  # left_join(bracket_with_ratings %>% select(team,region), by = "team") #too lazy to fix can remove

team_results <- sim_results %>%
  group_by(team, seed) %>%
  summarise(
    `Round 32` = (mean(Round_32)),
    `Sweet 16` = (mean(Sweet_16)),
    `Elite 8` = (mean(Elite_8)),
    `Final 4` = (mean(Final_4)),
    `TitleGame` = (mean(Title_Game)),
    Champion = (mean(Champion)),
    .groups = "drop"
  ) %>%
  arrange((Champion)) 



# ------------------------------
# 3. Display Results
# ------------------------------

bets <- team_results %>%
  gt() %>%
  tab_header(title = "March Madness Betting Guide Based on 10,000 Simulations") %>%
  fmt_number(columns = 3:8, decimals = 0)
gtsave(bets,"bets.png")
prob_to_american(
sim_results %>% 
  left_join(bracket_with_ratings %>% select(team,standard_team), by = c("team")) %>% 
  left_join(teams_std %>% select(display_name,conference_short_name), by = c("standard_team" = "display_name")) %>% 
  mutate(simulation = rep(1:n_simulations, each = 64)) %>% 
  filter(conference_short_name == "SEC") %>% 
  group_by(simulation) %>% 
  summarize(total_wins = sum(c_across(c(`Round_32`, `Sweet_16`, `Elite_8`, Final_4, Title_Game, Champion)), na.rm = TRUE)) %>%
  summarize(mean(total_wins >= 21)) %>% 
  pull())

prob_to_american(sim_results %>% 
  mutate(simulation = rep(1:n_simulations, each = 64)) %>% 
  filter(seed %in% c(15)) %>%
    # Only keep 12-seeds
  group_by(simulation) %>%                  # Group by simulation
  summarise(num_12_wins = sum(Round_32, na.rm = TRUE)) %>%  # Count 12-seed wins per sim
  summarise(prob = mean(num_12_wins >= 1)) %>%   # Calculate probability of 2+ wins
  pull())

prob_to_american(
sim_results %>% 
  mutate(simulation = rep(1:n_simulations, each = 64)) %>% 
  filter(Final_4 == 1) %>% 
  group_by(simulation) %>% 
  summarize(total_wins = sum(as.numeric(seed), na.rm = TRUE)) %>%
  summarize(mean(total_wins > 11)) %>% 
  pull())
#Worst Seed Final 4----
  sim_results %>% 
    mutate(simulation = rep(1:n_simulations, each = 64)) %>% 
    filter(Elite_8 == 1) %>% 
    group_by(simulation) %>% 
    summarize(low_seed = max(as.numeric(seed), na.rm = TRUE)) %>%
    group_by(low_seed) %>% 
    summarize(prob_to_american(n()/n_simulations))

prob_to_american(team_results %>% 
  filter(seed == 1) %>% 
  summarize(sum(Champion)) %>% 
  pull())

prob_to_american(sim_results %>% 
  left_join(bracket_with_ratings %>% select(team,standard_team), by = c("team")) %>% 
  left_join(teams_std %>% select(display_name,conference_short_name), by = c("standard_team" = "display_name")) %>% 
  mutate(simulation = rep(1:n_simulations, each = 64)) %>% 
  filter(conference_short_name == "SEC") %>% 
  group_by(simulation) %>% 
  summarize(total_wins = sum(Champion, na.rm = TRUE)) %>%
  summarize(mean(total_wins)) %>% 
  pull())

library(tidyverse)

# Load your simulated data


# Load sportsbook odds
sportsbook_odds <- read_csv("ncaa_elite8_odds.csv")

# Function to convert probability to American odds
prob_to_american <- function(prob) {
  ifelse(prob >= 0.5, 
         -100 * prob / (1 - prob),  # Favorite (negative odds)
         100 / prob - 100)          # Underdog (positive odds)
}

# Function to convert American odds to implied probability
american_to_prob <- function(odds) {
  ifelse(odds < 0, abs(odds) / (abs(odds) + 100), 100 / (odds + 100))
}

# Function to calculate expected value (EV) for a $100 bet
calculate_ev <- function(prob, odds) {
  implied_prob <- american_to_prob(odds)
  payout <- ifelse(odds > 0, (odds / 100) * 100, (100 / abs(odds)) * 100)  # Payout per $100 bet
  (prob * payout) - ((1 - prob) * 100)  # EV formula
}

# Convert simulated probabilities to American odds
team_results <- team_results %>%
  mutate(Odds = prob_to_american(`Elite 8`))  # Assuming the column is "Elite 8"

# Merge sportsbook and simulated data
comparison <- sportsbook_odds %>%
  left_join(team_results %>% select(team, Odds), by = c("Team"="team")) %>%
  rename(Simulated_Odds = Odds.y, Sportsbook_Odds = Odds.x) %>% 
  mutate(Simulated_Odds = ifelse(Bet == "NO", Simulated_Odds*-1,Simulated_Odds))

# Calculate expected value (EV)
comparison <- comparison %>%
  mutate(
    Simulated_Prob = american_to_prob(Simulated_Odds),
    Sportsbook_Prob = american_to_prob(Sportsbook_Odds),
    EV = calculate_ev(Simulated_Prob, Sportsbook_Odds),
    Value_Bet = case_when(
      Bet == "YES" & Simulated_Prob > Sportsbook_Prob ~ "Yes Bet (Value)",
      Bet == "NO" & Simulated_Prob < Sportsbook_Prob ~ "No Bet (Value)",
      TRUE ~ "No Value"
    )
  )

# Show best value bets sorted by EV
best_bets <- comparison %>%
  filter(Value_Bet != "No Value") %>%
  arrange(desc(EV))

# View results
print(best_bets)

# Save the results as a CSV for further analysis


