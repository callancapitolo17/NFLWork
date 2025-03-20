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
library(rvest)
library(httr)
library(jsonlite)
library(dplyr)
library(purrr)
library(readr)
library(tidyr)
library(googlesheets4)
library(stringr)
# gs4_auth()

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
final_bracket <- bracket_data
# # --- Play-In Data ---
# # Scrape play-in game pods.
# game_nodes <- bracket_page %>% html_nodes("div.game-pod")
# 
# # Loop over each play-in game pod and extract details
# play_in_games <- lapply(game_nodes, function(node) {
#   
#   # Extract the game ID from the <a> element (if needed)
#   game_id <- node %>% html_node("a") %>% html_attr("id")
#   
#   # Extract top team details:
#   team_top_seed <- node %>% 
#     html_node("div.team.team-top span.overline") %>% 
#     html_text(trim = TRUE)
#   team_top_name <- node %>% 
#     html_node("div.team.team-top p.body") %>% 
#     html_text(trim = TRUE)
#   
#   # Extract bottom team details:
#   team_bottom_seed <- node %>% 
#     html_node("div.team.team-bottom span.overline") %>% 
#     html_text(trim = TRUE)
#   team_bottom_name <- node %>% 
#     html_node("div.team.team-bottom p.body") %>% 
#     html_text(trim = TRUE)
#   
#   # Extract the region from the <span class="subtitle"> after the <a> element.
#   region <- node %>% 
#     html_node("span.subtitle") %>% 
#     html_text(trim = TRUE)
#   
#   data.frame(
#     game_id = game_id,
#     team_top = team_top_name,
#     seed_top = team_top_seed,
#     team_bottom = team_bottom_name,
#     seed_bottom = team_bottom_seed,
#     region = region,
#     stringsAsFactors = FALSE
#   )
# })
# 
# # Combine all play-in game data into one data frame
# play_in_games_df <- bind_rows(play_in_games)
# # Reshape play-in data to match main bracket format (one row per team)
# play_in_bracket <- play_in_games_df %>%
#   transmute(team = team_top, seed = seed_top, region = region) %>%
#   bind_rows(
#     play_in_games_df %>% transmute(team = team_bottom, seed = seed_bottom, region = region)
#   ) %>%
#   mutate(play_in = 1)
# 
# # --- Combine Main Bracket and Play-In Data ---
# final_bracket <- bind_rows(bracket_data, play_in_bracket) %>% 
#   # Optionally, sort by region and seed (note: seeds are still character type here)
#   arrange(region, seed) %>%
#   # Recode region abbreviations from play-in if necessary:
#   mutate(region = recode(region,
#                          "E"  = "East",
#                          "MW" = "Midwest",
#                          "S"  = "South"))
# 

# ------------------------------
# 2. Acquire Power Ratings
# ------------------------------
# a) ESPN's BPI (scraped live)
library(httr)
library(jsonlite)
library(dplyr)
library(purrr)
library(tibble)

library(hoopR)
teams_std <- espn_mbb_teams(2025) %>% 
  mutate(team = ifelse(team == "McNeese", "McNeese State",team))

# Define a helper function to clean text (lowercase, remove non-alphanumerics)
clean_text <- function(x) {
  x <- sub("St\\.$", "State", x, ignore.case = TRUE)
  # Now remove any non-alphanumeric characters.
  return(x)
}

get_standard_team <- function(team, teams_std) {
  # If the input team is NA, return NA immediately.
  if (is.na(team)) return(NA_character_)
  
  team_clean <- clean_text(team)
  
  # Loop through each row in the hoopR table.
  for(i in 1:nrow(teams_std)) {
    # Check each variant column.
    variants <- teams_std[i, c("abbreviation", "display_name", "short_name", "mascot", "nickname","team")]
    for(variant in variants) {
      # Only compare if the variant is not NA.
      if (!is.na(variant)) {
        # If both cleaned texts are not NA and are equal, return the canonical name.
        if (!is.na(variant) && (team_clean == clean_text(variant))) {
          return(teams_std$display_name[i])
        }
      }
    }
  }
  # If no match is found, return the original team name.
  return(team)
}
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

clean_bpi_data <- final_bpi_data %>% 
  mutate(standard_team = map_chr(team, ~ get_standard_team(.x, teams_std = teams_std)))


# b) KenPom Ratings (local CSV file)

# Replace with your actual Google Sheets URL
sheet_url <- "https://docs.google.com/spreadsheets/d/10o9dwZeyREliM8iOIevpHIhNeCHYeVGE-y50N3KHujQ/edit?gid=0#gid=0"

# Read the sheet (first sheet by default)
kenpom_data <- read_sheet(sheet_url) %>% 
  mutate(Team = str_replace(Team, "\\s\\d+$", "")) %>% 
  as.data.frame() %>% 
  mutate(Team = ifelse(Team == "Connecticut", "UConn",ifelse(Team == "Mississippi", "Ole Miss",
                                                             ifelse(Team == "Nebraska Omaha", "Omaha", Team)))) 

clean_kenpom_data <- kenpom_data %>%
  mutate(standard_team = map_chr(Team, ~ get_standard_team(.x, teams_std = teams_std))) %>% 
  select(standard_team, NetRtg) %>% 
  mutate(NetRtg = map_dbl(NetRtg, ~ ifelse(is.null(.x), NA_real_, as.numeric(.x))))
#torvik

library(httr)
library(jsonlite)

torvik_url <- "https://barttorvik.com/trank.php?year=2025&t=0&json=1"

# Set user-agent to mimic a browser request
response <- GET(torvik_url, add_headers("User-Agent" = "Mozilla/5.0"))

# Check if request was successful
if (status_code(response) == 200) {
  torvik_data <- as.data.frame(fromJSON(content(response, "text", encoding = "UTF-8")))
  
  colnames(torvik_data)[c(1, 2, 3, 4)] <- c("Team", "OffEff", "DefEff", "TorvikPower")
  
  # Convert columns to numeric
  torvik_data <- torvik_data %>%
    mutate(
      OffEff = as.numeric(OffEff),
      DefEff = as.numeric(DefEff),
      TorvikPower = as.numeric(TorvikPower),
      TorvikMargin = ((OffEff - DefEff) / 100) * 70
    ) %>% 
    select(Team, TorvikMargin) 
} else {
  stop("Failed to retrieve Torvik data. HTTP status:", status_code(response))
}

clean_torvik_data <- torvik_data %>% 
  mutate(Team = ifelse(Team == "Connecticut", "UConn",ifelse(Team == "Mississippi", "Ole Miss",
                                                             ifelse(Team == "Nebraska Omaha", "Omaha", Team)))) %>% 
  mutate(standard_team = map_chr(Team, ~ get_standard_team(.x, teams_std = teams_std))) %>% 
  select(standard_team, TorvikMargin)

# Read in EvanMiya data
evan_miya <- read_sheet("https://docs.google.com/spreadsheets/u/0/d/1zaDMGo6bRMis_VD7yH4uZY9aFXbp0QUPp3HRoE8x3EY/edit?usp=drive_web&pli=1&authuser=0")

# Determine the grouping factor (every 22 rows should be one row)
group_size <- 21  # Set the number of rows per team
evan_miya$group <- rep(1:(nrow(evan_miya)/group_size), each = group_size)[1:nrow(evan_miya)]  # Assign group numbers

# Reshape data: Convert every 22 rows into a single row
evan_miya <- evan_miya %>%
  group_by(group) %>%
  mutate(row_number = row_number()) %>%
  ungroup() %>% 
  as.data.frame()

# Pivot into wide format (ensuring teams align correctly)
evan_miya_wide <- evan_miya %>%
  pivot_wider(names_from = row_number, values_from = 1) %>%  
  select(-group) %>% 
  mutate(`2` = `2` %>% str_replace_all("[^A-Za-z0-9 ]", "") %>% 
           str_squish()) 

colnames(evan_miya_wide) <- c(
  "Relative Ranking", "Team", "O-Rate", "D-Rate", "Relative Rating",
  "Opponent Adjust", "Pace Adjust", "Off Rank", "Def Rank", "True Tempo",
  "Tempo Rank", "Injury Rank", "Home Rank", "Roster Rank",
  "Kill Shots Per Game", "Kill Shots Conceded Per Game", "Kill Shots Margin Per Game",
  "Total Kill Shots", "Total Kill Shots Conceded", "D1 Wins", "D1 Losses"
)

clean_evan_miya_wide <- evan_miya_wide %>% 
  mutate(Team = ifelse(Team == "Connecticut", "UConn",ifelse(Team == "Mississippi", "Ole Miss",Team))) %>% 
  mutate(Team = sub("AM","A&M",Team)) %>% 
  mutate(standard_team = map_chr(Team, ~ get_standard_team(.x, teams_std = teams_std %>% 
                                                             mutate(nickname = sub("'", "", nickname),
                                                                    nickname = sub("St.", "Saint ",nickname)
                                                                    )))) %>% 
  select(standard_team, `Relative Rating`) %>% 
  mutate(`Relative Rating` = map_dbl(`Relative Rating`, ~ ifelse(is.null(.x), NA_real_, as.numeric(.x))))

 
# Now power_ratings$standard_team will contain the standardized team names.

power_ratings <- clean_bpi_data %>%
  left_join(clean_kenpom_data %>% rename(KenPomRating = NetRtg), by = c("standard_team")) %>%
  left_join(clean_torvik_data, by = c("standard_team")) %>%
  left_join(clean_evan_miya_wide %>% rename(EvanMiyaRating = `Relative Rating`), by = c("standard_team")) %>%  #fix evan miya
  mutate(,
    EvanMiyaRating = as.numeric(unlist(EvanMiyaRating)),
    # Normalize KenPom & EvanMiya (convert from per 100 possessions to per 70 possessions)
    KenPomMargin = (as.numeric(KenPomRating) / 100) * 70,
    EvanMiyaMargin = (EvanMiyaRating / 100) * 70,
    
    # Sagarin needs empirical scaling
    
    # BPI & Torvik are already correct
    BPIMargin = bpi
  ) %>%
  select(team, KenPomMargin, BPIMargin, EvanMiyaMargin, TorvikMargin)

# ------------------------------
# 3. Merge Bracket Data with Power Ratings
# ------------------------------
# Standardize team names and convert seeds if possible.
final_bracket <- final_bracket %>%
  mutate(standard_team = map_chr(team, ~ get_standard_team(.x, teams_std = teams_std)))

# If seeds are purely numeric strings, convert them; otherwise, keep as character.
final_bracket <- final_bracket %>%
  mutate(seed = as.numeric(seed))

# Merge using inner_join to keep only teams for which we have ratings.
bracket_with_ratings <- left_join(final_bracket, power_ratings, by = c("standard_team" = "team")) %>% 
  rowwise() %>%
  mutate(composite_rating = mean(c_across(KenPomMargin:TorvikMargin), na.rm = TRUE)) %>%
  ungroup()

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


