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

# Centralized name fixes: maps short/alternate names to the canonical form
# that hoopR's get_standard_team() can resolve. Add entries here when a
# rating source uses a name that doesn't match any hoopR variant.
TEAM_NAME_FIXES <- c(
  # BPI / cbbdata short names
  "Connecticut"          = "UConn",
  "Miami FL"             = "Miami",
  "Hawaii"               = "Hawai'i",
  "LIU Brooklyn"         = "LIU",
  "Long Island University" = "LIU",
  "Cal Baptist"          = "California Baptist",
  "NC State"             = "NC State",
  "North Carolina St."   = "NC State",
  # TeamRankings short names
  "Ohio St"              = "Ohio State",
  "Iowa St"              = "Iowa State",
  "Utah St"              = "Utah State",
  "South Fla"            = "South Florida",
  "USF"                  = "South Florida",
  "N Iowa"               = "Northern Iowa",
  "UNI"                  = "Northern Iowa",
  # EvanMiya names (after punctuation stripping)
  "Texas AM"             = "Texas A&M",
  "Prairie View AM"      = "Prairie View A&M",
  "Florida AM"           = "Florida A&M",
  "Saint Marys"          = "Saint Mary's",
  "St Johns"             = "St. John's",
  "Miami Ohio"           = "Miami (OH)",
  "Miami Fla"            = "Miami",
  # TeamRankings names
  "S Florida"            = "South Florida",
  # KenPom names
  "N.C. State"           = "NC State",
  "Mississippi"          = "Ole Miss",
  "Nebraska Omaha"       = "Omaha",
  # Queens (new D1 program, not in hoopR — map to bracket name)
  "Queens"               = "Queens University Royals",
  "Queens University"    = "Queens University Royals"
)

clean_text <- function(x) {
  sub("St\\.$", "State", x, ignore.case = TRUE)
}

#' Normalize a team name: apply TEAM_NAME_FIXES first, then match against hoopR
normalize_team_name <- function(team) {
  if (team %in% names(TEAM_NAME_FIXES)) {
    return(TEAM_NAME_FIXES[[team]])
  }
  team
}

get_standard_team <- function(team, teams_std) {
  if (is.na(team)) return(NA_character_)
  team <- normalize_team_name(team)
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
  cat("Fetching BPI ratings (direct ESPN API)...\n")
  base_url <- "https://site.web.api.espn.com/apis/fitt/v3/sports/basketball/mens-college-basketball/powerindex?region=us&lang=en&groups=50&limit=50&page="
  tryCatch({
    all_teams <- list()
    for (page in 1:8) {
      resp <- GET(paste0(base_url, page))
      if (status_code(resp) != 200) next
      data <- fromJSON(content(resp, "text", encoding = "UTF-8"))
      team_names <- data$teams$team$displayName
      bpi_values <- map_dbl(data$teams$categories, function(cat) {
        if (is.data.frame(cat) && "name" %in% names(cat)) {
          bpi_row <- cat %>% filter(name == "bpi")
          if (nrow(bpi_row) > 0 && !is.null(bpi_row$values[[1]][1])) {
            return(bpi_row$values[[1]][1])
          }
        }
        NA_real_
      })
      all_teams[[page]] <- tibble(team = team_names, bpi = bpi_values)
    }
    bpi_raw <- bind_rows(all_teams)
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
      as.data.frame()
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
  # Try 1: Google Sheets (manual paste from barttorvik.com)
  torvik_sheet <- "https://docs.google.com/spreadsheets/d/1Vt6AXl0rACgTAF1oBaO55mFOLtYsnbMMp0WKbXZFQug/edit?gid=0#gid=0"
  result <- tryCatch({
    raw <- read_sheet(torvik_sheet)
    # Skip header row (row 1 = column labels), data starts at row 2
    # Col 2 = Team (with seed/info on newline), Col 6 = AdjOE (with rank on newline), Col 7 = AdjDE (with rank on newline)
    torvik_df <- tibble(
      Team = as.character(raw[[2]][-1]),
      AdjOE = as.character(raw[[6]][-1]),
      AdjDE = as.character(raw[[7]][-1])
    ) %>%
      # Strip newline suffixes (rank numbers, seed info)
      mutate(
        Team = str_replace(Team, "\n.*$", ""),
        AdjOE = as.numeric(str_replace(AdjOE, "\n.*$", "")),
        AdjDE = as.numeric(str_replace(AdjDE, "\n.*$", ""))
      ) %>%
      filter(!is.na(AdjOE), !is.na(AdjDE), Team != "Team") %>%
      distinct(Team, .keep_all = TRUE) %>%  # Remove duplicate pastes
      mutate(
        TorvikMargin = ((AdjOE - AdjDE) / 100) * 70,
        standard_team = map_chr(Team, ~ get_standard_team(.x, teams_std = teams_std))
      ) %>%
      select(standard_team, TorvikMargin)
    # Validate: should have 300+ unique teams, reject if too few
    n_unique <- n_distinct(torvik_df$standard_team)
    if (n_unique < 300) {
      cat(sprintf("  Warning: Torvik sheet only has %d unique teams (expected 350+). Data may be corrupt.\n", n_unique))
      stop("Torvik sheet data quality too low")
    }
    cat(sprintf("Torvik: %d teams\n", nrow(torvik_df)))
    torvik_df
  }, error = function(e) {
    cat(sprintf("  Google Sheet failed (%s), trying cbbdata...\n", e$message))
    # Try 2: cbbdata package (may have older season data)
    tryCatch({
      torvik <- cbd_torvik_ratings(year = current_year)
      if (nrow(torvik) == 0) stop("No data for current year")
      torvik %>%
        mutate(TorvikMargin = ((adj_o - adj_d) / 100) * 70) %>%
        transmute(standard_team = team, TorvikMargin)
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
      mutate(Team = ifelse(Team == "Mississippi", "Ole Miss", Team)) %>%
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

fetch_teamrankings <- function(teams_std) {
  cat("Fetching TeamRankings predictive ratings...\n")
  tryCatch({
    page <- rvest::read_html("https://www.teamrankings.com/ncaa-basketball/ranking/predictive-by-other")
    df <- (page %>% rvest::html_table())[[1]]
    # Strip record from team name: "Duke (32-2)" -> "Duke"
    result <- df %>%
      transmute(
        team = str_replace(Team, "\\s*\\(.*\\)$", ""),
        TeamRankingsMargin = Rating
      ) %>%
      mutate(standard_team = map_chr(team, ~ get_standard_team(.x, teams_std = teams_std)))
    cat(sprintf("TeamRankings: %d teams\n", nrow(result)))
    result
  }, error = function(e) {
    cat(sprintf("Warning: TeamRankings fetch failed: %s\n", e$message))
    tibble(team = character(), TeamRankingsMargin = numeric(), standard_team = character())
  })
}

fetch_powerrank <- function(teams_std) {
  cat("Fetching Power Rank ratings...\n")
  tryCatch({
    page <- rvest::read_html("https://thepowerrank.com/college-basketball-rankings/")
    df <- (page %>% rvest::html_table())[[1]]
    result <- df %>%
      transmute(
        team = str_replace(Team, "\\s*\\(.*\\)$", ""),
        PowerRankMargin = Rating
      ) %>%
      mutate(standard_team = map_chr(team, ~ get_standard_team(.x, teams_std = teams_std)))
    cat(sprintf("PowerRank: %d teams\n", nrow(result)))
    result
  }, error = function(e) {
    cat(sprintf("Warning: PowerRank fetch failed: %s\n", e$message))
    tibble(team = character(), PowerRankMargin = numeric(), standard_team = character())
  })
}

# =============================================================================
# Combined Power Ratings + Bracket Merge
# =============================================================================

#' Fetch all power ratings and merge into a single table
#' Uses parallel fetching to speed up (6 sources simultaneously)
fetch_power_ratings <- function(teams_std) {
  # Parallel fetch: all 6 sources at once
  if (requireNamespace("future", quietly = TRUE) && requireNamespace("future.apply", quietly = TRUE)) {
    cat("Fetching ratings in parallel...\n")
    future::plan(future::multisession, workers = 6)
    fetchers <- list(
      bpi = function() fetch_bpi(teams_std),
      kenpom = function() fetch_kenpom(teams_std),
      torvik = function() fetch_torvik(teams_std),
      evanmiya = function() fetch_evan_miya(teams_std),
      teamrankings = function() fetch_teamrankings(teams_std),
      powerrank = function() fetch_powerrank(teams_std)
    )
    results <- future.apply::future_lapply(fetchers, function(f) f(), future.seed = TRUE)
    future::plan(future::sequential)  # Reset
    clean_bpi_data <- results$bpi
    clean_kenpom_data <- results$kenpom
    clean_torvik_data <- results$torvik
    clean_evan_miya_data <- results$evanmiya
    clean_teamrankings_data <- results$teamrankings
    clean_powerrank_data <- results$powerrank
  } else {
    # Fallback: sequential
    clean_bpi_data <- fetch_bpi(teams_std)
    clean_kenpom_data <- fetch_kenpom(teams_std)
    clean_torvik_data <- fetch_torvik(teams_std)
    clean_evan_miya_data <- fetch_evan_miya(teams_std)
    clean_teamrankings_data <- fetch_teamrankings(teams_std)
    clean_powerrank_data <- fetch_powerrank(teams_std)
  }

  # Deduplicate each source before joining (prevents many-to-many)
  clean_kenpom_data <- clean_kenpom_data %>% filter(!is.na(standard_team), standard_team != "Team") %>% distinct(standard_team, .keep_all = TRUE)
  clean_evan_miya_data <- clean_evan_miya_data %>% filter(!is.na(standard_team)) %>% distinct(standard_team, .keep_all = TRUE)

  clean_bpi_data %>%
    left_join(clean_kenpom_data %>% rename(KenPomRating = NetRtg), by = "standard_team") %>%
    left_join(clean_torvik_data, by = "standard_team") %>%
    left_join(clean_evan_miya_data %>% rename(EvanMiyaRating = `Relative Rating`), by = "standard_team") %>%
    left_join(clean_teamrankings_data %>% select(standard_team, TeamRankingsMargin), by = "standard_team") %>%
    left_join(clean_powerrank_data %>% select(standard_team, PowerRankMargin), by = "standard_team") %>%
    mutate(
      EvanMiyaRating = as.numeric(unlist(EvanMiyaRating)),
      KenPomMargin = (as.numeric(KenPomRating) / 100) * 70,
      EvanMiyaMargin = (EvanMiyaRating / 100) * 70,
      BPIMargin = bpi
    ) %>%
    select(standard_team, KenPomMargin, BPIMargin, EvanMiyaMargin, TorvikMargin, TeamRankingsMargin, PowerRankMargin)
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
    mutate(composite_rating = median(c_across(KenPomMargin:PowerRankMargin), na.rm = TRUE)) %>%
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

#' Fast vectorized tournament simulation
#' Pre-resolves First Four, uses vectorized rnorm for all games per round,
#' tracks progress with base R matrix (no dplyr in hot loop).
#'
#' @param bracket_64 A 64-team data.frame (First Four already resolved)
#' @param region_order Character vector of region names for F4 pairing
simulate_tournament_fast <- function(bracket_64, region_order = NULL) {
  # Convert to base R data.frame and extract key vectors
  df <- as.data.frame(bracket_64)
  if (is.null(region_order)) region_order <- unique(df$region)

  all_teams <- df$team
  all_seeds <- df$seed
  n_all <- length(all_teams)
  progress <- matrix(0L, nrow = n_all, ncol = 6)

  # Seed-based matchup order (1v16, 8v9, 5v12, 4v13, 6v11, 3v14, 7v10, 2v15)
  seed_order <- c(1, 16, 8, 9, 5, 12, 4, 13, 6, 11, 3, 14, 7, 10, 2, 15)

  # Working copies
  teams <- df$team
  seeds <- df$seed
  regions <- df$region
  ratings <- ifelse(is.na(df$composite_rating), 0, df$composite_rating)
  n <- length(teams)
  round_num <- 1

  while (n > 1) {
    if (n == 4) {
      # Final Four: pair by region order (1v3, 2v4)
      reg_idx <- match(regions, region_order)
      ord <- order(reg_idx)
      ord <- ord[c(1, 3, 2, 4)]  # Cross-bracket pairing
    } else if (n == 2) {
      ord <- 1:2
    } else {
      # Regional rounds: order by region then seed bracket position
      seed_pos <- match(seeds, seed_order)
      ord <- order(match(regions, region_order), seed_pos)
    }

    # Reorder
    teams <- teams[ord]
    seeds <- seeds[ord]
    regions <- regions[ord]
    ratings <- ratings[ord]

    # Vectorized game sim: odd indices vs even indices
    n_games <- n %/% 2
    i1 <- seq(1, n, by = 2)
    i2 <- seq(2, n, by = 2)

    expected_diff <- ratings[i1] - ratings[i2]
    actual_margin <- rnorm(n_games, mean = expected_diff, sd = 11.2)
    w <- ifelse(actual_margin > 0, i1, i2)

    # Rating update
    residual <- actual_margin - expected_diff
    rating_change <- 0.1 * residual + 0.05 * residual * log(round_num + 1) + rnorm(n_games, 0, 0.5)

    # Keep winners
    teams <- teams[w]
    seeds <- seeds[w]
    regions <- regions[w]
    ratings <- ratings[w] + abs(rating_change)
    n <- length(teams)

    # Track progress
    col <- min(round_num, 6)
    progress[match(teams, all_teams), col] <- 1L
    round_num <- round_num + 1
  }

  data.frame(
    team = all_teams, seed = all_seeds,
    Round_32 = progress[, 1], Sweet_16 = progress[, 2], Elite_8 = progress[, 3],
    Final_4 = progress[, 4], Title_Game = progress[, 5], Champion = progress[, 6],
    stringsAsFactors = FALSE
  )
}

#' Resolve First Four: reduce 68-team bracket to 64 teams.
#' If a First Four game has been played (in games_df), use the actual winner.
#' If not yet played, simulate the game. Returns a 64-team bracket.
resolve_first_four <- function(bracket_with_ratings, games_df = NULL) {
  play_in_teams <- bracket_with_ratings %>% filter(play_in == 1)
  main_teams    <- bracket_with_ratings %>% filter(play_in == 0)

  if (nrow(play_in_teams) == 0) return(bracket_with_ratings)

  matchups <- play_in_teams %>%
    group_by(region, seed) %>%
    group_split()

  first_four_winners <- map_dfr(matchups, function(pair) {
    if (nrow(pair) != 2) return(pair[1, ])

    if (!is.null(games_df) && nrow(games_df) > 0) {
      played <- games_df %>%
        filter(round == "First Four", status == "final",
               (team1 == pair$team[1] & team2 == pair$team[2]) |
               (team1 == pair$team[2] & team2 == pair$team[1]))
      if (nrow(played) > 0) {
        winner_name <- played$winner[1]
        return(pair %>% filter(team == winner_name))
      }
    }

    simulate_game(pair[1, ], pair[2, ])$winner
  })

  bind_rows(main_teams, first_four_winners)
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
