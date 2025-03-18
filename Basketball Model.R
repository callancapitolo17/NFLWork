# ============================
# NCAA Tournament Simulator
# ============================
# This script:
# 1. Acquires bracket data (team names, seeds, regions, and play-in flag).
# 2. Acquires power ratings from three sources:
#    - ESPN’s BPI (scraped live)
#    - KenPom ratings (from a local CSV file)
#    - PoolGenius ratings (from a local CSV file)
# 3. Combines these ratings into a composite rating (weighted average).
# 4. Merges the composite ratings with the bracket data.
# 5. Runs a Monte Carlo simulation of the tournament,
#    including simulating the First Four (play-in) games.
#
# Load required libraries:
library(rvest)
library(httr)
library(jsonlite)
library(dplyr)
library(purrr)
library(readr)
library(tidyr)

# ------------------------------
# 1. Acquire Bracket Data
# ------------------------------
# (Scrape main bracket teams and play-in teams separately, then join.)

# --- Main Bracket Data ---
bracket_url <- "https://www.ncaa.com/march-madness-live/bracket"  # update if needed
bracket_page <- read_html(bracket_url)

# Get region containers (4 regions in this example)
region_nodes <- bracket_page %>% html_nodes("div.bracket-container div.region")

# Initialize an empty data frame to hold all main bracket teams
bracket_data <- data.frame(
  team = character(),
  seed = character(),
  region = character(),
  stringsAsFactors = FALSE
)

# Loop over each region container for the main bracket
for (region in region_nodes) {
  # Extract region name from the <span> element inside the region container
  region_name <- region %>% 
    html_node("span.subtitle") %>% 
    html_text(trim = TRUE)
  
  # Extract team names from <p class="body_2"> within this region
  teams <- region %>% 
    html_nodes("p.body_2") %>% 
    html_text(trim = TRUE)
  teams <- teams[nzchar(teams)]
  
  # Extract seeds from the team containers; adjust the selector if needed.
  seeds <- region %>% 
    html_nodes("div.team > span") %>% 
    html_text(trim = TRUE)
  seeds <- seeds[nzchar(seeds)]
  
  # Optional: warn if teams and seeds don't match
  if (length(teams) != length(seeds)) {
    cat("Warning: For region", region_name, "found", length(teams), 
        "teams and", length(seeds), "seeds.\n")
  }
  
  # Create temporary data frame for this region (mark play_in = 0 for main bracket)
  temp_df <- data.frame(
    team = teams,
    seed = seeds,
    region = rep(region_name, length(teams)),
    play_in = 0,
    stringsAsFactors = FALSE
  )
  
  # Append to overall bracket_data
  bracket_data <- bind_rows(bracket_data, temp_df)
}

# --- Play-In Data ---
# Scrape play-in game pods.
game_nodes <- bracket_page %>% html_nodes("div.game-pod")

# Loop over each play-in game pod and extract details
play_in_games <- lapply(game_nodes, function(node) {
  
  # Extract the game ID from the <a> element (if needed)
  game_id <- node %>% html_node("a") %>% html_attr("id")
  
  # Extract top team details:
  team_top_seed <- node %>% 
    html_node("div.team.team-top span.overline") %>% 
    html_text(trim = TRUE)
  team_top_name <- node %>% 
    html_node("div.team.team-top p.body") %>% 
    html_text(trim = TRUE)
  
  # Extract bottom team details:
  team_bottom_seed <- node %>% 
    html_node("div.team.team-bottom span.overline") %>% 
    html_text(trim = TRUE)
  team_bottom_name <- node %>% 
    html_node("div.team.team-bottom p.body") %>% 
    html_text(trim = TRUE)
  
  # Extract the region from the <span class="subtitle"> after the <a> element.
  region <- node %>% 
    html_node("span.subtitle") %>% 
    html_text(trim = TRUE)
  
  data.frame(
    game_id = game_id,
    team_top = team_top_name,
    seed_top = team_top_seed,
    team_bottom = team_bottom_name,
    seed_bottom = team_bottom_seed,
    region = region,
    stringsAsFactors = FALSE
  )
})

# Combine all play-in game data into one data frame
play_in_games_df <- bind_rows(play_in_games)
# Reshape play-in data to match main bracket format (one row per team)
play_in_bracket <- play_in_games_df %>%
  transmute(team = team_top, seed = seed_top, region = region) %>%
  bind_rows(
    play_in_games_df %>% transmute(team = team_bottom, seed = seed_bottom, region = region)
  ) %>%
  mutate(play_in = 1)

# --- Combine Main Bracket and Play-In Data ---
final_bracket <- bind_rows(bracket_data, play_in_bracket) %>% 
  # Optionally, sort by region and seed (note: seeds are still character type here)
  arrange(region, seed) %>%
  # Recode region abbreviations from play-in if necessary:
  mutate(region = recode(region,
                         "E"  = "East",
                         "MW" = "Midwest",
                         "S"  = "South"))

# View final bracket data
print(final_bracket)

# ------------------------------
# 2. Acquire Power Ratings
# ------------------------------
# a) ESPN's BPI (scraped live)
library(httr)
library(jsonlite)
library(dplyr)
library(purrr)
library(tibble)

# Base URL
base_url <- "https://site.web.api.espn.com/apis/fitt/v3/sports/basketball/mens-college-basketball/powerindex?region=us&lang=en&groups=50&limit=50&page="

# Initialize empty list to store results
all_teams <- list()

# Loop through pages (assuming a max of 8 pages based on pagination info)
for (page in 1:8) {
  url <- paste0(base_url, page)  # Construct page-specific URL
  response <- GET(url)
  
  if (status_code(response) != 200) {
    cat("Error fetching page", page, "\n")
    next
  }
  
  content_text <- content(response, "text", encoding = "UTF-8")
  data <- fromJSON(content_text)
  
  # Extract team names
  team_names <- data$teams$team$displayName
  
  # Extract BPI values
  bpi_values <- map_dbl(data$teams$categories, function(cat) {
    if (is.data.frame(cat) && "name" %in% names(cat)) {
      bpi_row <- cat %>% filter(name == "bpi")
      if (nrow(bpi_row) > 0 && !is.null(bpi_row$values[[1]][1])) {
        return(bpi_row$values[[1]][1])
      }
    }
    return(NA_real_)
  })
  
  # Store results as a tibble
  team_data <- tibble(
    team = team_names,
    bpi = bpi_values
  )
  
  all_teams[[page]] <- team_data
}

# Combine all pages into a single data frame
final_bpi_data <- bind_rows(all_teams)

# View final results
print(final_bpi_data)

test <- read.csv("export (1).csv")


# b) KenPom Ratings (local CSV file)
library(googlesheets4)
gs4_auth()
# Replace with your actual Google Sheets URL
sheet_url <- "https://docs.google.com/spreadsheets/d/10o9dwZeyREliM8iOIevpHIhNeCHYeVGE-y50N3KHujQ/edit?gid=0#gid=0"

# Read the sheet (first sheet by default)
kenpom_data <- read_sheet(sheet_url) %>% 
  mutate(Team = str_replace(Team, "\\s\\d+$", ""))

#torvik

torvik_url <- "https://barttorvik.com/trank.php?year=2025&t=0&json=1"
torvik_data <- as.data.frame(fromJSON(torvik_url))
colnames(torvik_data)[c(1, 2, 3,4)] <- c("Team", "OffEff", "DefEff", "TorvikPower")

# Convert columns to numeric
torvik_data <- torvik_data %>%
  mutate(
    OffEff = as.numeric(OffEff),
    DefEff = as.numeric(DefEff),
    TorvikPower = as.numeric(TorvikPower),
    TorvikMargin = ((OffEff-DefEff)/100)*70
  ) %>% 
  select(Team,TorvikMargin)

sagarin_url <- "https://sagarin.com/sports/cbsend.htm"
sagarin_page <- read_html(sagarin_url)
sagarin_table <- html_table(sagarin_page, fill = TRUE)[[1]]
sagarin_ratings <- sagarin_table %>% select(team = Team, sagarin = Predictor)

evanmiya_url <- "https://evanmiya.com/?team_ratings"
evanmiya_page <- read_html(evanmiya_url)
library(tidyr)
library(tidyverse)
library(googlesheets4)

# Read data from Google Sheets
library(tidyverse)
library(googlesheets4)

# Read in EvanMiya data
evan_miya <- read_sheet(sheet_url2)

# Determine the grouping factor (every 22 rows should be one row)
group_size <- 21  # Set the number of rows per team
evan_miya$group <- rep(1:(nrow(evan_miya)/group_size), each = group_size)[1:nrow(evan_miya)]  # Assign group numbers

# Reshape data: Convert every 22 rows into a single row
evan_miya <- evan_miya %>%
  group_by(group) %>%
  mutate(row_number = row_number()) %>%
  ungroup()

# Pivot into wide format (ensuring teams align correctly)
evan_miya_wide <- evan_miya %>%
  pivot_wider(names_from = row_number, values_from = 1) %>%  
  select(-group) %>% 
  mutate(`2` = gsub("[^[:alnum:][:space:]]", "", evan_miya_wide$Team))

colnames(evan_miya_wide) <- c(
  "Relative Ranking", "Team", "O-Rate", "D-Rate", "Relative Rating",
  "Opponent Adjust", "Pace Adjust", "Off Rank", "Def Rank", "True Tempo",
  "Tempo Rank", "Injury Rank", "Home Rank", "Roster Rank",
  "Kill Shots Per Game", "Kill Shots Conceded Per Game", "Kill Shots Margin Per Game",
  "Total Kill Shots", "Total Kill Shots Conceded", "D1 Wins", "D1 Losses"
)
# Merge power ratings from all sources.
power_ratings <- espn_bpi_ratings %>% 
  full_join(kenpom_ratings, by = "team") %>% 
  full_join(poolgenius_ratings, by = "team") %>%
  na.omit() %>%
  mutate(composite_rating = 0.4 * espn_bpi + 0.35 * kenpom_rating + 0.25 * poolgenius_rating)

# ------------------------------
# 3. Merge Bracket Data with Power Ratings
# ------------------------------
# Standardize team names and convert seeds if possible.
final_bracket <- final_bracket %>%
  mutate(team = str_trim(team))
power_ratings <- power_ratings %>%
  mutate(team = str_trim(team))

# If seeds are purely numeric strings, convert them; otherwise, keep as character.
final_bracket <- final_bracket %>%
  mutate(seed = as.numeric(seed))

# Merge using inner_join to keep only teams for which we have ratings.
bracket_with_ratings <- inner_join(final_bracket, power_ratings, by = "team")
cat("Teams with power ratings available:", nrow(bracket_with_ratings), "\n\n")

# ------------------------------
# 4. Simulation Functions
# ------------------------------

# Function to simulate a game between two teams using composite ratings.
simulate_game <- function(team1, team2, game_number = 1, beta1 = 0.1, beta2 = 0.05, sd_margin = 11.2) {
  expected_diff <- team1$composite_rating - team2$composite_rating
  win_prob_team1 <- 1 / (1 + exp(-0.163 * expected_diff))
  actual_margin <- rnorm(1, mean = expected_diff, sd = sd_margin)
  
  if (actual_margin > 0) {
    winner <- team1
    loser <- team2
  } else {
    winner <- team2
    loser <- team1
  }
  
  rating_change <- beta1 * (actual_margin - expected_diff) +
    beta2 * (actual_margin - expected_diff) * log(game_number + 1)
  
  updated_team1 <- team1
  updated_team2 <- team2
  if (winner$team == team1$team) {
    updated_team1$composite_rating <- team1$composite_rating + rating_change
    updated_team2$composite_rating <- team2$composite_rating - rating_change
  } else {
    updated_team1$composite_rating <- team1$composite_rating - rating_change
    updated_team2$composite_rating <- team2$composite_rating + rating_change
  }
  
  return(list(winner = winner,
              updated_team1 = updated_team1,
              updated_team2 = updated_team2))
}

# Function to simulate one round (pair teams two-by-two).
simulate_round <- function(teams, game_number = 1) {
  n <- nrow(teams)
  winners <- data.frame()
  
  # Sort teams by seed for pairing (assuming lower seed numbers are higher-ranked)
  teams <- teams %>% arrange(seed)
  
  for (i in seq(1, n, by = 2)) {
    team1 <- teams[i, ]
    team2 <- teams[i + 1, ]
    game_result <- simulate_game(team1, team2, game_number)
    winners <- bind_rows(winners, game_result$winner)
  }
  
  return(winners)
}

# Function to simulate the full tournament.
simulate_tournament <- function(bracket) {
  round_num <- 1
  teams_round <- bracket
  
  # Process play-in games first based on play_in flag
  if (any(teams_round$play_in == 1)) {
    play_in_subset <- teams_round %>% filter(play_in == 1)
    play_in_winners <- simulate_round(play_in_subset, game_number = round_num)
    
    # Remove play-in teams and add winners back into the main bracket
    teams_round <- teams_round %>% filter(play_in == 0)
    teams_round <- bind_rows(teams_round, play_in_winners)
  }
  
  # Simulate subsequent rounds until one champion remains.
  while (nrow(teams_round) > 1) {
    teams_round <- simulate_round(teams_round, game_number = round_num)
    round_num <- round_num + 1
  }
  
  return(teams_round$team)
}

# ------------------------------
# 5. Monte Carlo Tournament Simulation
# ------------------------------
n_simulations <- 1000
champion_counts <- setNames(rep(0, nrow(bracket_with_ratings)), bracket_with_ratings$team)

for (i in 1:n_simulations) {
  champ <- simulate_tournament(bracket_with_ratings)
  champion_counts[champ] <- champion_counts[champ] + 1
}

champion_probabilities <- champion_counts / n_simulations
champion_probabilities <- sort(champion_probabilities, decreasing = TRUE)

cat("Champion Probabilities (after", n_simulations, "simulations):\n")
print(champion_probabilities)
