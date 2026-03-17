library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(hoopR)
library(googlesheets4)
library(jsonlite)
library(readr)
library(gt)
library(httr)
library(cbbdata)
# gs4_auth()

# --- Configuration ---
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
final_bracket <- bracket_result$bracket %>%
  select(team, seed, region, play_in)

# ------------------------------
# 2. Acquire Power Ratings (shared logic with Basketball Model.R)
# ------------------------------
teams_std <- espn_mbb_teams(current_year) %>%
  mutate(team = ifelse(team == "McNeese", "McNeese State", team))

clean_text <- function(x) sub("St\\.$", "State", x, ignore.case = TRUE)

get_standard_team <- function(team, teams_std) {
  if (is.na(team)) return(NA_character_)
  team_clean <- clean_text(team)
  for (i in 1:nrow(teams_std)) {
    variants <- teams_std[i, c("abbreviation", "display_name", "short_name", "mascot", "nickname", "team")]
    for (variant in variants) {
      if (!is.na(variant) && team_clean == clean_text(variant)) return(teams_std$display_name[i])
    }
  }
  team
}

# a) BPI via cbbdata
clean_bpi_data <- tryCatch({
  cbd_bpi_ratings() %>% transmute(team = team, bpi = bpi_value, standard_team = team)
}, error = function(e) {
  cat(sprintf("Warning: BPI fetch failed: %s\n", e$message))
  tibble(team = character(), bpi = numeric(), standard_team = character())
})

# b) KenPom from Google Sheets
sheet_url <- "https://docs.google.com/spreadsheets/d/10o9dwZeyREliM8iOIevpHIhNeCHYeVGE-y50N3KHujQ/edit?gid=0#gid=0"
clean_kenpom_data <- tryCatch({
  read_sheet(sheet_url) %>%
    mutate(Team = str_replace(Team, "\\s\\d+$", "")) %>% as.data.frame() %>%
    mutate(Team = ifelse(Team == "Connecticut", "UConn",
                  ifelse(Team == "Mississippi", "Ole Miss",
                  ifelse(Team == "Nebraska Omaha", "Omaha", Team)))) %>%
    mutate(standard_team = map_chr(Team, ~ get_standard_team(.x, teams_std = teams_std))) %>%
    select(standard_team, NetRtg) %>%
    mutate(NetRtg = map_dbl(NetRtg, ~ ifelse(is.null(.x), NA_real_, as.numeric(.x))))
}, error = function(e) {
  cat(sprintf("Warning: KenPom fetch failed: %s\n", e$message))
  tibble(standard_team = character(), NetRtg = numeric())
})

# c) Torvik via cbbdata (fallback to direct API)
clean_torvik_data <- tryCatch({
  torvik <- cbd_torvik_ratings(year = current_year)
  if (nrow(torvik) == 0) stop("No data for current year")
  torvik %>% mutate(TorvikMargin = ((adj_o - adj_d) / 100) * 70) %>%
    transmute(standard_team = team, TorvikMargin)
}, error = function(e) {
  tryCatch({
    resp <- GET(sprintf("https://barttorvik.com/trank.php?year=%d&t=0&json=1", current_year),
                add_headers("User-Agent" = "Mozilla/5.0"))
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
    cat(sprintf("Warning: Torvik unavailable: %s\n", e2$message))
    tibble(standard_team = character(), TorvikMargin = numeric())
  })
})

# d) EvanMiya from Google Sheets
clean_evan_miya_wide <- tryCatch({
  evan_miya <- read_sheet("https://docs.google.com/spreadsheets/u/0/d/1zaDMGo6bRMis_VD7yH4uZY9aFXbp0QUPp3HRoE8x3EY/edit?usp=drive_web&pli=1&authuser=0")
  group_size <- 21
  evan_miya$group <- rep(1:(nrow(evan_miya)/group_size), each = group_size)[1:nrow(evan_miya)]
  evan_miya <- evan_miya %>% group_by(group) %>% mutate(row_number = row_number()) %>% ungroup() %>% as.data.frame()
  evan_miya_wide <- evan_miya %>% pivot_wider(names_from = row_number, values_from = 1) %>% select(-group) %>%
    mutate(`2` = `2` %>% str_replace_all("[^A-Za-z0-9 ]", "") %>% str_squish())
  colnames(evan_miya_wide) <- c(
    "Relative Ranking", "Team", "O-Rate", "D-Rate", "Relative Rating",
    "Opponent Adjust", "Pace Adjust", "Off Rank", "Def Rank", "True Tempo",
    "Tempo Rank", "Injury Rank", "Home Rank", "Roster Rank",
    "Kill Shots Per Game", "Kill Shots Conceded Per Game", "Kill Shots Margin Per Game",
    "Total Kill Shots", "Total Kill Shots Conceded", "D1 Wins", "D1 Losses")
  if (any(grepl("^[0-9.]+$", evan_miya_wide$Team))) {
    stop("EvanMiya pivot failed: Team column contains numeric values. Sheet format may have changed.")
  }
  evan_miya_wide %>%
    mutate(Team = case_when(
      Team == "Connecticut" ~ "UConn",
      Team == "Mississippi" ~ "Ole Miss",
      grepl("^Texas AM$", Team) ~ "Texas A&M",
      grepl("^Prairie View AM$", Team) ~ "Prairie View A&M",
      grepl("^Florida AM$", Team) ~ "Florida A&M",
      TRUE ~ Team
    )) %>%
    mutate(standard_team = map_chr(Team, ~ get_standard_team(.x, teams_std = teams_std))) %>%
    select(standard_team, `Relative Rating`) %>%
    mutate(`Relative Rating` = map_dbl(`Relative Rating`, ~ ifelse(is.null(.x), NA_real_, as.numeric(.x))))
}, error = function(e) {
  cat(sprintf("Warning: EvanMiya fetch failed: %s\n", e$message))
  tibble(standard_team = character(), `Relative Rating` = numeric())
})

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
# 3. Merge Bracket Data with Power Ratings & Define Current State
# ------------------------------
final_bracket <- final_bracket %>%
  mutate(standard_team = map_chr(team, ~ get_standard_team(.x, teams_std = teams_std))) %>%
  mutate(seed = as.numeric(seed))
bracket_with_ratings <- left_join(final_bracket, power_ratings, by = "standard_team") %>%
  rowwise() %>%
  mutate(composite_rating = mean(c_across(KenPomMargin:TorvikMargin), na.rm = TRUE)) %>%
  ungroup()

# --- Extract the Current Bracket State (from ESPN API) ---
tourney_state <- bracket_result$tournament_state
games_played  <- bracket_result$games

# Build current_bracket: only surviving teams
current_bracket <- bracket_with_ratings %>%
  filter(!team %in% {
    games_played %>%
      filter(status == "final", !is.na(winner)) %>%
      mutate(loser = ifelse(winner == team1, team2, team1)) %>%
      pull(loser)
  })

# ------------------------------
# 4. Simulation Functions
# ------------------------------

# (A) Helper to generate matchups.
get_region_order <- function(final_bracket) {
  # Extract the unique regions in the order they appear.
  final_bracket %>% distinct(region) %>% pull(region)
}

# Updated get_bracket_matchups for simulation.
get_bracket_matchups <- function(teams, region_order_auto = NULL) {
  n <- nrow(teams)
  if(n == 4) {
    if(is.null(region_order_auto)) {
      region_order_auto <- c("South", "East", "West", "Midwest")
    }
    teams <- teams %>% mutate(region_index = match(region, region_order_auto))
    teams <- teams %>% arrange(region_index)
    # For Final Four, we want the first vs. third and second vs. fourth.
    teams <- teams[c(1, 3, 2, 4), ]
    teams <- teams %>% select(-region_index)
    return(teams)
  } else {
    match_order <- tibble(
      seed = c(1,16,8,9,5,12,4,13,6,11,3,14,7,10,2,15),
      matchup_order = 1:16
    )
    teams %>%
      left_join(match_order, by = "seed") %>%
      arrange(region, matchup_order) %>%
      select(-matchup_order)
  }
}

simulate_game <- function(team1, team2, game_number = 1, beta1 = 0.1, beta2 = 0.05,
                          sd_margin = 11.2, sd_rating_shift = 0.5) {
  if(is.na(team1$composite_rating)) team1$composite_rating <- 0
  if(is.na(team2$composite_rating)) team2$composite_rating <- 0
  expected_diff <- team1$composite_rating - team2$composite_rating
  actual_margin <- rnorm(1, mean = expected_diff, sd = sd_margin)
  winner <- if(actual_margin > 0) team1 else team2
  loser <- if(actual_margin > 0) team2 else team1
  # Rating update with residual sampling (per Unabated article)
  rating_change <- beta1 * (actual_margin - expected_diff) +
    beta2 * (actual_margin - expected_diff) * log(game_number + 1) +
    rnorm(1, 0, sd_rating_shift)
  winner$composite_rating <- winner$composite_rating + rating_change
  loser$composite_rating <- loser$composite_rating - rating_change
  list(winner = winner)
}

simulate_round <- function(teams, game_number = 1, region_order_auto = NULL) {
  teams <- get_bracket_matchups(teams, region_order_auto)
  team_pairs <- split(teams, rep(1:(nrow(teams)/2), each = 2))
  winners <- map(team_pairs, function(pair) {
    team1 <- pair[1, ]
    team2 <- pair[2, ]
    if(nrow(pair) < 2) return(team1)
    simulate_game(team1, team2, game_number)$winner
  })
  bind_rows(winners)
}

get_remaining_rounds <- function(current_N) {
  if(current_N == 64) {
    list(round_names = c("Round 32", "Sweet 16", "Elite 8", "Final 4", "Title Game", "Champion"))
  } else if(current_N == 32) {
    list(round_names = c("Sweet 16", "Elite 8", "Final 4", "Title Game", "Champion"))
  } else if(current_N == 16) {
    list(round_names = c("Elite 8", "Final 4", "Title Game", "Champion"))
  } else if(current_N == 8) {
    list(round_names = c("Final 4", "Title Game", "Champion"))
  } else if(current_N == 4) {
    list(round_names = c("Title Game", "Champion"))
  } else if(current_N == 2) {
    list(round_names = c("Champion"))
  } else {
    stop("Bracket size must be 2, 4, 8, 16, 32, or 64.")
  }
}

simulate_remaining_tournament <- function(bracket, region_order_auto = NULL) {
  if(is.null(region_order_auto)) {
    region_order_auto <- get_region_order(bracket)
  }
  current_N <- nrow(bracket)
  rounds_info <- get_remaining_rounds(current_N)
  round_names <- rounds_info$round_names
  team_progress <- bracket %>% select(team, seed)
  for(r in round_names) {
    team_progress[[r]] <- 0
  }
  round_num <- 1
  teams_round <- bracket
  while(nrow(teams_round) > 1) {
    teams_round <- simulate_round(teams_round, game_number = round_num, region_order_auto = region_order_auto)
    if(round_num <= length(round_names)) {
      current_round_name <- round_names[round_num]
      team_progress <- team_progress %>%
        mutate("{current_round_name}" := ifelse(team %in% teams_round$team, 1, .data[[current_round_name]]))
    }
    round_num <- round_num + 1
  }
  team_progress
}

# ------------------------------
# 5. Monte Carlo Tournament Simulation
# ------------------------------

n_simulations <- 10000
sim_results <- map_dfr(1:n_simulations, ~ simulate_remaining_tournament(current_bracket))

# Define prob_to_american (must be defined before using in mutate).
prob_to_american <- function(prob) {
  ifelse(prob >= 0.5, -100 * prob / (1 - prob), 100 / prob - 100)
}

team_results <- sim_results %>%
  group_by(team, seed) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  { 
    # Identify simulation round columns (all columns except team and seed)
    sim_rounds <- setdiff(names(.), c("team", "seed"))
    
    # If there is at least one simulated round, set the first as p_current
    if (length(sim_rounds) >= 1) {
      p_current_col <- sim_rounds[1]
      # All remaining rounds are used to compute a cumulative future win probability.
      future_rounds <- sim_rounds[-1]
      
      # Compute the option value:
      # p_current = probability of advancing from the current round (first column)
      # f_future = average win probability in later rounds (can be adjusted to different weights)
      mutate(., 
             p_current = .data[[p_current_col]],
             f_future = if(length(future_rounds) > 0) rowMeans(select(., all_of(future_rounds))) else 0,
             `Survivor Value` = p_current * (1 - f_future)
      )
    } else {
      .
    }
  } %>%
  arrange(desc(`Survivor Value`))

# Optionally, convert other probabilities to American odds (leave sim_option_value as a probability)
team_results <- team_results %>%
  select(-p_current,-f_future) %>% 
  arrange(desc(`Champion`)) %>% 
  mutate(`Survivor Value Rank` = min_rank(desc(`Survivor Value`))) %>% 
  mutate(across(-c(team, seed, `Survivor Value`,`Survivor Value Rank`), prob_to_american)) %>%  #Convert Survivor Value to Rank?
  select(-`Survivor Value`)
# Display the results with gt (or any method of your choice)
bets <- team_results %>%
  gt() %>%
  tab_header(title = "March Madness Betting Guide with Survivor Recommendations") %>%
  fmt_number(columns = c(3:ncol(team_results)), decimals = 0) %>% 
  gtExtras::gt_hulk_col_numeric(columns = `Survivor Value Rank`)
bets
gtsave(bets,"MarchMadness.png")
