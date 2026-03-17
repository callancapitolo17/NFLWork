# =============================================================================
# March Madness Simulator - Shared Utilities
# =============================================================================
# Common code used by Basketball Model.R, Dynamic Simulator.R, and
# End of Round Simulator.R. Source this file instead of duplicating.
#
# Provides:
#   - Library loading + current_year config
#   - Team name matching (clean_text, get_standard_team)
#   - Power rating fetchers (BPI, KenPom, Torvik, EvanMiya)
#   - fetch_power_ratings() — main entry point, returns merged ratings
#   - fetch_bracket_with_ratings() — bracket + ratings + composite
#   - simulate_game() — single game simulation with residual sampling
#   - get_bracket_matchups() — NCAA bracket seeding order
#   - get_remaining_rounds() — round names by bracket size
#   - prob_to_american() — probability to American odds conversion
# =============================================================================

suppressPackageStartupMessages({
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
})

# --- Configuration ---
current_year <- if (as.integer(format(Sys.Date(), "%m")) >= 11) {
  as.integer(format(Sys.Date(), "%Y")) + 1L
} else {
  as.integer(format(Sys.Date(), "%Y"))
}

# =============================================================================
# ESPN Team Dictionary
# =============================================================================

#' Load ESPN team dictionary for the current season
get_teams_std <- function() {
  espn_mbb_teams(current_year) %>%
    mutate(team = ifelse(team == "McNeese", "McNeese State", team))
}

# =============================================================================
# Team Name Matching
# =============================================================================

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

# =============================================================================
# Power Rating Fetchers
# =============================================================================

fetch_bpi <- function(teams_std) {
  cat("Fetching BPI ratings...\n")
  # cbbdata BPI uses short names that don't always match hoopR display_name
  bpi_name_fixes <- c(
    "Connecticut"  = "UConn",
    "Miami FL"     = "Miami",
    "Hawaii"       = "Hawai'i",
    "LIU Brooklyn" = "LIU",
    "Long Island University" = "LIU",
    "Cal Baptist"  = "California Baptist",
    "Queens"       = "Queens University",
    "NC State"     = "NC State",
    "North Carolina St." = "NC State"
  )
  tryCatch({
    bpi_raw <- cbd_bpi_ratings() %>%
      transmute(team = team, bpi = bpi_value) %>%
      mutate(team = ifelse(team %in% names(bpi_name_fixes), bpi_name_fixes[team], team))
    # Resolve to full display names ("Duke" -> "Duke Blue Devils")
    result <- bpi_raw %>%
      mutate(standard_team = map_chr(team, ~ get_standard_team(.x, teams_std = teams_std)))
    cat(sprintf("BPI: %d teams\n", nrow(result)))
    result
  }, error = function(e) {
    cat(sprintf("Warning: BPI fetch failed: %s\n", e$message))
    tibble(team = character(), bpi = numeric(), standard_team = character())
  })
}

fetch_kenpom <- function(teams_std) {
  cat("Fetching KenPom ratings...\n")
  sheet_url <- "https://docs.google.com/spreadsheets/d/10o9dwZeyREliM8iOIevpHIhNeCHYeVGE-y50N3KHujQ/edit?gid=0#gid=0"
  tryCatch({
    kenpom_data <- read_sheet(sheet_url) %>%
      mutate(Team = str_replace(Team, "\\s\\d+$", "")) %>%
      as.data.frame() %>%
      mutate(Team = ifelse(Team == "Connecticut", "UConn",
                    ifelse(Team == "Mississippi", "Ole Miss",
                    ifelse(Team == "Nebraska Omaha", "Omaha", Team))))
    result <- kenpom_data %>%
      mutate(standard_team = map_chr(Team, ~ get_standard_team(.x, teams_std = teams_std))) %>%
      select(standard_team, NetRtg) %>%
      mutate(NetRtg = map_dbl(NetRtg, ~ ifelse(is.null(.x), NA_real_, as.numeric(.x))))
    cat(sprintf("KenPom: %d teams\n", nrow(result)))
    result
  }, error = function(e) {
    cat(sprintf("Warning: KenPom fetch failed: %s\n", e$message))
    tibble(standard_team = character(), NetRtg = numeric())
  })
}

fetch_torvik <- function(teams_std) {
  cat("Fetching Torvik ratings...\n")
  tryCatch({
    torvik <- cbd_torvik_ratings(year = current_year)
    if (nrow(torvik) == 0) stop("No data for current year")
    result <- torvik %>%
      mutate(TorvikMargin = ((adj_o - adj_d) / 100) * 70) %>%
      transmute(standard_team = team, TorvikMargin)
    cat(sprintf("Torvik: %d teams\n", nrow(result)))
    result
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
      result <- torvik_data %>%
        mutate(OffEff = as.numeric(OffEff), DefEff = as.numeric(DefEff),
               TorvikMargin = ((OffEff - DefEff) / 100) * 70,
               standard_team = map_chr(Team, ~ get_standard_team(.x, teams_std = teams_std))) %>%
        select(standard_team, TorvikMargin)
      cat(sprintf("Torvik: %d teams\n", nrow(result)))
      result
    }, error = function(e2) {
      cat(sprintf("  Warning: Torvik unavailable: %s\n", e2$message))
      tibble(standard_team = character(), TorvikMargin = numeric())
    })
  })
}

fetch_evan_miya <- function(teams_std) {
  cat("Fetching EvanMiya ratings...\n")
  tryCatch({
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
    # Validate pivot didn't silently scramble data
    if (any(grepl("^[0-9.]+$", evan_miya_wide$Team))) {
      stop("EvanMiya pivot failed: Team column contains numeric values. Sheet format may have changed.")
    }
    result <- evan_miya_wide %>%
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
    cat(sprintf("EvanMiya: %d teams\n", nrow(result)))
    result
  }, error = function(e) {
    cat(sprintf("Warning: EvanMiya fetch failed: %s\n", e$message))
    tibble(standard_team = character(), `Relative Rating` = numeric())
  })
}

# =============================================================================
# Combined Power Ratings + Bracket Merge
# =============================================================================

#' Fetch all power ratings and merge into a single table
fetch_power_ratings <- function(teams_std) {
  clean_bpi_data <- fetch_bpi(teams_std)
  clean_kenpom_data <- fetch_kenpom(teams_std)
  clean_torvik_data <- fetch_torvik(teams_std)
  clean_evan_miya_data <- fetch_evan_miya(teams_std)

  clean_bpi_data %>%
    left_join(clean_kenpom_data %>% rename(KenPomRating = NetRtg), by = "standard_team") %>%
    left_join(clean_torvik_data, by = "standard_team") %>%
    left_join(clean_evan_miya_data %>% rename(EvanMiyaRating = `Relative Rating`), by = "standard_team") %>%
    mutate(
      EvanMiyaRating = as.numeric(unlist(EvanMiyaRating)),
      KenPomMargin = (as.numeric(KenPomRating) / 100) * 70,
      EvanMiyaMargin = (EvanMiyaRating / 100) * 70,
      BPIMargin = bpi
    ) %>%
    select(standard_team, KenPomMargin, BPIMargin, EvanMiyaMargin, TorvikMargin)
}

#' Fetch bracket from ESPN, merge with power ratings, compute composite
fetch_bracket_with_ratings <- function(bracket_df, teams_std) {
  power_ratings <- fetch_power_ratings(teams_std)

  bracket_merged <- bracket_df %>%
    mutate(
      standard_team = map_chr(team, ~ get_standard_team(.x, teams_std = teams_std)),
      seed = as.numeric(seed)
    ) %>%
    left_join(power_ratings, by = "standard_team") %>%
    rowwise() %>%
    mutate(composite_rating = mean(c_across(KenPomMargin:TorvikMargin), na.rm = TRUE)) %>%
    ungroup()

  n_matched <- sum(!is.na(bracket_merged$composite_rating))
  cat(sprintf("Bracket: %d teams, %d with ratings\n", nrow(bracket_merged), n_matched))
  bracket_merged
}

# =============================================================================
# Simulation Functions
# =============================================================================

#' Simulate a single game between two teams
simulate_game <- function(team1, team2, game_number = 1, beta1 = 0.1, beta2 = 0.05,
                          sd_margin = 11.2, sd_rating_shift = 0.5) {
  if (is.na(team1$composite_rating)) team1$composite_rating <- 0
  if (is.na(team2$composite_rating)) team2$composite_rating <- 0
  expected_diff <- team1$composite_rating - team2$composite_rating
  actual_margin <- rnorm(1, mean = expected_diff, sd = sd_margin)

  winner <- if (actual_margin > 0) team1 else team2
  loser  <- if (actual_margin > 0) team2 else team1

  # Rating update with residual sampling (per Unabated article: capture
  # unmeasured factors like shooting variance, foul trouble, injuries)
  rating_change <- beta1 * (actual_margin - expected_diff) +
    beta2 * (actual_margin - expected_diff) * log(game_number + 1) +
    rnorm(1, 0, sd_rating_shift)

  winner$composite_rating <- winner$composite_rating + rating_change
  loser$composite_rating  <- loser$composite_rating - rating_change

  list(winner = winner)
}

#' Order teams by NCAA bracket seeding for correct matchups
#' Handles both regional rounds (seed-based) and Final Four (region-based)
get_bracket_matchups <- function(teams, region_order_auto = NULL) {
  n <- nrow(teams)
  if (n == 4) {
    # Final Four: pair by region bracket position
    if (is.null(region_order_auto)) {
      region_order_auto <- c("South", "East", "West", "Midwest")
    }
    teams <- teams %>% mutate(region_index = match(region, region_order_auto))
    teams <- teams %>% arrange(region_index)
    teams <- teams[c(1, 3, 2, 4), ]
    teams <- teams %>% select(-region_index)
    return(teams)
  } else {
    match_order <- tibble(
      seed = c(1, 16, 8, 9, 5, 12, 4, 13, 6, 11, 3, 14, 7, 10, 2, 15),
      matchup_order = 1:16
    )
    teams %>%
      left_join(match_order, by = "seed") %>%
      arrange(region, matchup_order) %>%
      select(-matchup_order)
  }
}

#' Extract region order from bracket data
get_region_order <- function(bracket) {
  bracket %>% distinct(region) %>% pull(region)
}

#' Get remaining round names based on current bracket size
get_remaining_rounds <- function(current_N) {
  switch(as.character(current_N),
    "64" = c("Round 32", "Sweet 16", "Elite 8", "Final 4", "Title Game", "Champion"),
    "32" = c("Sweet 16", "Elite 8", "Final 4", "Title Game", "Champion"),
    "16" = c("Elite 8", "Final 4", "Title Game", "Champion"),
    "8"  = c("Final 4", "Title Game", "Champion"),
    "4"  = c("Title Game", "Champion"),
    "2"  = c("Champion"),
    stop("Bracket size must be 2, 4, 8, 16, 32, or 64.")
  )
}

#' Convert probability to American odds
prob_to_american <- function(prob) {
  ifelse(prob >= 0.5,
         -100 * prob / (1 - prob),
         100 / prob - 100)
}
