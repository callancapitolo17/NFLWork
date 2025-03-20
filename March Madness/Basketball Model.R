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
library(googlesheets4)
library(stringr)
gs4_auth()

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
  mutate(Team = ifelse(Team == "Connecticut", "UConn",Team))

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
  mutate(Team = ifelse(Team == "Connecticut", "UConn",ifelse(Team == "Mississippi", "Ole Miss",Team))) %>% 
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
bracket_with_ratings <- left_join(final_bracket, power_ratings, by = c("standard_team" = "team"))
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
  teams <- teams %>% arrange(seed)
  matchups <- teams %>%
    mutate(pairing = rep(1:(n() / 2), each = 2)) %>%
    group_by(pairing) %>%
    summarise(winner = list(simulate_game(cur_data()[1,], cur_data()[2,], game_number)$winner))
  return(bind_rows(matchups$winner))
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
