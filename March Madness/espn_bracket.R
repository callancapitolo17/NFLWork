# =============================================================================
# ESPN Bracket Fetcher for March Madness Simulator
# =============================================================================
# Replaces fragile NCAA.com HTML scraping with ESPN's public scoreboard API.
# Returns bracket data (teams, seeds, regions), game results, and tournament state.
#
# Usage:
#   source("espn_bracket.R")
#   result <- fetch_espn_bracket()       # current year
#   result <- fetch_espn_bracket(2025)   # specific year
#
#   final_bracket   <- result$bracket          # data.frame: team, seed, region, play_in
#   games           <- result$games            # data.frame: game-level results
#   tourney_state   <- result$tournament_state # list: state, current_round, etc.
# =============================================================================

suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(tidyr)
})

# Optionally source Tools.R for team resolution utilities (suppress Rcpp warnings)
tryCatch(
  suppressWarnings(suppressMessages(source("~/NFLWork/Answer Keys/Tools.R"))),
  error = function(e) message("Tools.R not available: ", e$message)
)

# --- Constants ---
ESPN_SCOREBOARD_URL <- "https://site.api.espn.com/apis/site/v2/sports/basketball/mens-college-basketball/scoreboard"
NCAA_TOURNAMENT_ID <- "22"

ROUND_ORDER <- c(
  "First Four", "1st Round", "2nd Round", "Sweet 16",
  "Elite 8", "Final Four", "National Championship"
)

# =============================================================================
# Core Functions
# =============================================================================

#' Generate date range to scan for tournament games
get_tournament_dates <- function(year) {
  start <- as.Date(paste0(year, "-03-14"))
  end   <- as.Date(paste0(year, "-04-10"))
  format(seq(start, end, by = "day"), "%Y%m%d")
}

#' Fetch one day of scoreboard data, return only tournament games (raw JSON lists)
fetch_tournament_day <- function(date_str) {
  url <- sprintf("%s?dates=%s&groups=100&limit=100", ESPN_SCOREBOARD_URL, date_str)

  resp <- tryCatch(
    GET(url, timeout(15)),
    error = function(e) {
      warning(sprintf("ESPN API request failed for %s: %s", date_str, e$message))
      return(NULL)
    }
  )

  if (is.null(resp) || status_code(resp) != 200) {
    return(list())
  }

  data <- fromJSON(content(resp, "text", encoding = "UTF-8"), simplifyVector = FALSE)
  events <- data$events %||% list()

  # Filter to NCAA tournament games only
  Filter(function(ev) {
    comps <- ev$competitions
    if (is.null(comps) || length(comps) == 0) return(FALSE)
    tid <- comps[[1]]$tournamentId
    !is.null(tid) && as.character(tid) == NCAA_TOURNAMENT_ID
  }, events)
}

#' Parse ESPN notes headline into region + round
#' Handles patterns like:
#'   "Men's Basketball Championship - Midwest Region - 1st Round"
#'   "Men's Basketball Championship - Final Four"
#'   "Men's Basketball Championship - National Championship"
parse_headline <- function(headline) {
  if (is.null(headline) || is.na(headline) || headline == "") {
    return(list(region = NA_character_, round = NA_character_))
  }

  parts <- strsplit(headline, " - ")[[1]]

  if (length(parts) == 3) {
    region_str <- trimws(parts[2])
    round_str  <- trimws(parts[3])
    # Strip " Region" suffix from region name
    region_str <- sub(" Region$", "", region_str)
    return(list(region = region_str, round = round_str))
  } else if (length(parts) == 2) {
    round_str <- trimws(parts[2])
    return(list(region = NA_character_, round = round_str))
  }

  list(region = NA_character_, round = NA_character_)
}

#' Parse a single ESPN event into a one-row game tibble
parse_event <- function(ev) {
  comp <- ev$competitions[[1]]
  competitors <- comp$competitors

  # Get headline from notes
  headline <- tryCatch({
    comp$notes[[1]]$headline
  }, error = function(e) NA_character_)
  if (is.null(headline)) headline <- NA_character_

  parsed <- parse_headline(headline)

  # Extract competitor details (ESPN always has exactly 2 competitors)
  get_competitor <- function(idx) {
    c <- competitors[[idx]]
    list(
      team = c$team$displayName %||% NA_character_,
      seed = as.integer(c$curatedRank$current %||% NA_integer_),
      score = as.integer(c$score %||% NA_integer_),
      winner = isTRUE(c$winner)
    )
  }

  c1 <- get_competitor(1)
  c2 <- get_competitor(2)

  # Game status
  status_raw <- comp$status$type$name %||% "STATUS_SCHEDULED"
  status <- case_when(
    status_raw == "STATUS_FINAL" ~ "final",
    status_raw == "STATUS_IN_PROGRESS" ~ "in_progress",
    TRUE ~ "scheduled"
  )

  # Determine winner
  winner_name <- NA_character_
  if (status == "final") {
    if (c1$winner) winner_name <- c1$team
    else if (c2$winner) winner_name <- c2$team
  }

  tibble(
    game_id = as.character(ev$id %||% NA_character_),
    round   = parsed$round,
    region  = parsed$region,
    team1   = c1$team,
    seed1   = c1$seed,
    score1  = c1$score,
    team2   = c2$team,
    seed2   = c2$seed,
    score2  = c2$score,
    winner  = winner_name,
    status  = status
  )
}

#' Detect tournament state from games data
detect_tournament_state <- function(games_df) {
  if (nrow(games_df) == 0) {
    return(list(
      state = "pre_tournament",
      current_round = NA_character_,
      completed_rounds = character(0),
      teams_remaining = 68L
    ))
  }

  round_status <- games_df %>%
    group_by(round) %>%
    summarise(
      total = n(),
      completed = sum(status == "final"),
      in_progress = sum(status == "in_progress"),
      .groups = "drop"
    ) %>%
    mutate(round_idx = match(round, ROUND_ORDER)) %>%
    filter(!is.na(round_idx)) %>%
    arrange(round_idx)

  completed_rounds <- round_status %>%
    filter(completed == total) %>%
    pull(round)

  active_round <- round_status %>%
    filter(completed < total) %>%
    arrange(round_idx) %>%
    slice(1) %>%
    pull(round)

  if (length(active_round) == 0 && "National Championship" %in% completed_rounds) {
    state <- "complete"
    current_round <- "National Championship"
  } else if (length(active_round) > 0) {
    state <- "in_progress"
    current_round <- active_round
  } else {
    state <- "pre_tournament"
    current_round <- NA_character_
  }

  # Map current round to remaining teams
  teams_remaining <- switch(current_round,
    "First Four"              = 68L,
    "1st Round"               = 64L,
    "2nd Round"               = 32L,
    "Sweet 16"                = 16L,
    "Elite 8"                 = 8L,
    "Final Four"              = 4L,
    "National Championship"   = 2L,
    68L  # default
  )

  list(
    state = state,
    current_round = current_round,
    completed_rounds = completed_rounds,
    teams_remaining = teams_remaining
  )
}

#' Build the 64-team bracket from games data
#' For pre-tournament: all teams from 1st Round games + First Four teams
#' For in-progress: tracks who has been eliminated
build_bracket <- function(games_df) {
  if (nrow(games_df) == 0) return(tibble())

  # Get all unique teams with their seeds and regions from 1st Round
  first_round <- games_df %>%
    filter(round == "1st Round")

  first_four <- games_df %>%
    filter(round == "First Four")

  # Teams from 1st Round (both sides of each matchup)
  bracket_teams <- bind_rows(
    first_round %>% transmute(team = team1, seed = seed1, region = region),
    first_round %>% transmute(team = team2, seed = seed2, region = region)
  ) %>%
    distinct(team, .keep_all = TRUE)

  # Add First Four teams that didn't make it into the 1st Round
  if (nrow(first_four) > 0) {
    first_four_teams <- bind_rows(
      first_four %>% transmute(team = team1, seed = seed1, region = region),
      first_four %>% transmute(team = team2, seed = seed2, region = region)
    ) %>%
      filter(!team %in% bracket_teams$team) %>%
      distinct(team, .keep_all = TRUE)

    bracket_teams <- bind_rows(bracket_teams, first_four_teams)
  }

  # Remove TBD/placeholder entries (ESPN pre-populates future rounds with "TBD")
  bracket_teams <- bracket_teams %>%
    filter(!is.na(team), team != "", team != "TBD")

  # Mark play-in teams
  first_four_team_names <- c(first_four$team1, first_four$team2)
  bracket_teams <- bracket_teams %>%
    mutate(play_in = as.integer(team %in% first_four_team_names))

  # Determine eliminated teams (losers of completed games)
  losers <- games_df %>%
    filter(status == "final", !is.na(winner)) %>%
    mutate(loser = ifelse(winner == team1, team2, team1)) %>%
    pull(loser)

  bracket_teams <- bracket_teams %>%
    mutate(eliminated = team %in% losers)

  bracket_teams
}

# =============================================================================
# Main Entry Point
# =============================================================================

#' Fetch the NCAA Tournament bracket from ESPN's scoreboard API
#'
#' @param year Tournament year (defaults to current year)
#' @return list with $bracket (team-level), $games (game-level), $tournament_state
fetch_espn_bracket <- function(year = as.integer(format(Sys.Date(), "%Y"))) {
  cat(sprintf("Fetching %d NCAA Tournament data from ESPN API...\n", year))

  dates <- get_tournament_dates(year)
  all_events <- list()

  for (date_str in dates) {
    events <- fetch_tournament_day(date_str)
    if (length(events) > 0) {
      all_events <- c(all_events, events)
      cat(sprintf("  %s: %d tournament games\n", date_str, length(events)))
    }
  }

  cat(sprintf("Found %d total tournament games.\n", length(all_events)))

  if (length(all_events) == 0) {
    cat("No tournament games found. Tournament may not have started yet.\n")
    return(list(
      bracket = tibble(
        team = character(), seed = integer(), region = character(),
        play_in = integer(), eliminated = logical()
      ),
      games = tibble(
        game_id = character(), round = character(), region = character(),
        team1 = character(), seed1 = integer(), score1 = integer(),
        team2 = character(), seed2 = integer(), score2 = integer(),
        winner = character(), status = character()
      ),
      tournament_state = list(
        state = "pre_tournament", current_round = NA_character_,
        completed_rounds = character(0), teams_remaining = 68L
      )
    ))
  }

  # Parse all events into games data frame
  games_df <- map_dfr(all_events, parse_event) %>%
    distinct(game_id, .keep_all = TRUE)

  # Build bracket and detect state
  bracket_df <- build_bracket(games_df)
  tourney_state <- detect_tournament_state(games_df)

  cat(sprintf("Bracket: %d teams | State: %s | Round: %s\n",
              nrow(bracket_df), tourney_state$state,
              tourney_state$current_round %||% "N/A"))

  list(
    bracket = bracket_df,
    games = games_df,
    tournament_state = tourney_state
  )
}
