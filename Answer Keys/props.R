#!/usr/bin/env Rscript
# Props Fair Odds Calculator
# Interactive tool for pricing arbitrary props from sample data
#
# Usage:
#   source("props.R")
#   prop_odds("touchdown", "> 4.5")
#   prop_odds("pass_touchdown", "> 2.5", period = "Half1")
#   prop_odds("field_goal_attempt", "> 3.5", team = "home")

suppressPackageStartupMessages({
  library(duckdb)
  library(dplyr)
  library(purrr)
})

setwd("~/NFLWork/Answer Keys")
source("Tools.R")

# =============================================================================
# CONFIGURATION
# =============================================================================

VALID_PERIODS <- c("1", "2", "3", "4", "Half1", "Half2", "Full")
VALID_TEAMS <- c("home", "away")
VALID_OPS <- c(">", ">=", "<", "<=", "==")

# Common prop columns (for help display)
COMMON_COLUMNS <- c(
  "touchdown", "pass_touchdown", "rush_touchdown", "return_touchdown",
  "field_goal_attempt", "safety",
  "extra_point_attempt", "two_point_attempt"
)

# =============================================================================
# SAMPLE FRESHNESS CHECK
# =============================================================================

check_samples_fresh <- function(con, max_age_minutes = 5) {
  tryCatch({
    meta <- dbGetQuery(con, "SELECT generated_at FROM nfl_samples_meta")
    if (nrow(meta) == 0) return(FALSE)
    generated_at <- as.POSIXct(meta$generated_at, tz = "UTC")
    age_minutes <- as.numeric(difftime(Sys.time(), generated_at, units = "mins"))
    return(age_minutes <= max_age_minutes)
  }, error = function(e) {
    return(FALSE)
  })
}

generate_samples <- function() {
  cat("Generating fresh samples...\n")
  result <- system("Rscript 'NFL Answer Key/NFLPrepare.R'", intern = FALSE)
  if (result != 0) {
    stop("Error generating samples.")
  }
}

ensure_fresh_samples <- function() {
  con <- dbConnect(duckdb(), dbdir = "pbp.duckdb")
  is_fresh <- check_samples_fresh(con, max_age_minutes = 5)
  dbDisconnect(con)

  if (!is_fresh) {
    generate_samples()
  }
}

# =============================================================================
# COLUMN VALIDATION
# =============================================================================

get_valid_columns <- function(con) {
  dbGetQuery(con, "SELECT column_name FROM information_schema.columns WHERE table_name = 'nfl_pbp'")$column_name
}

validate_column <- function(column, valid_columns) {
  if (!column %in% valid_columns) {
    # Find similar columns for suggestion
    similar <- valid_columns[grepl(substr(column, 1, 4), valid_columns, ignore.case = TRUE)]
    if (length(similar) > 0) {
      stop(sprintf("Invalid column '%s'. Did you mean: %s?", column, paste(similar[1:min(5, length(similar))], collapse = ", ")))
    } else {
      stop(sprintf("Invalid column '%s'. Run show_columns() to see available columns.", column))
    }
  }
}

# =============================================================================
# PARSE CONDITION
# =============================================================================

parse_condition <- function(condition) {
  # Parse condition string like "> 4.5" into operator and value
  condition <- trimws(condition)

  # Try to match operator at start
  op_match <- regmatches(condition, regexpr("^(>=|<=|==|>|<)", condition))

  if (length(op_match) == 0 || op_match == "") {
    stop("Invalid condition. Must start with an operator: >, >=, <, <=, ==")
  }

  op <- op_match
  value_str <- trimws(sub("^(>=|<=|==|>|<)", "", condition))
  value <- as.numeric(value_str)

  if (is.na(value)) {
    stop(sprintf("Invalid value '%s'. Must be a number.", value_str))
  }

  list(op = op, value = value)
}

# =============================================================================
# PERIOD FILTER
# =============================================================================

get_period_filter <- function(period) {
  if (is.null(period) || period == "Full") {
    return(NULL)  # No filter, all quarters
  }

  switch(period,
    "1" = "qtr = 1",
    "2" = "qtr = 2",
    "3" = "qtr = 3",
    "4" = "qtr = 4",
    "Half1" = "qtr IN (1, 2)",
    "Half2" = "qtr IN (3, 4)",
    stop(sprintf("Invalid period '%s'. Valid: %s", period, paste(VALID_PERIODS, collapse = ", ")))
  )
}

# =============================================================================
# MAIN FUNCTION
# =============================================================================

#' Calculate fair odds for a prop bet
#'
#' @param column Column from nfl_pbp to aggregate (e.g., "touchdown")
#' @param condition Operator and line (e.g., "> 4.5")
#' @param period Optional period filter: "1", "2", "3", "4", "Half1", "Half2", "Full"
#' @param team Optional team filter: "home", "away"
#' @param agg Aggregation function: "SUM" (default), "COUNT", "MAX"
#' @return List with probability and fair odds
#'
#' @examples
#' prop_odds("touchdown", "> 4.5")
#' prop_odds("pass_touchdown", "> 2.5", period = "Half1")
#' prop_odds("field_goal_attempt", "> 3.5", team = "home")
prop_odds <- function(column = NULL, condition = NULL, period = NULL, team = NULL, agg = "SUM") {

  # Show help if no args
  if (is.null(column) || is.null(condition)) {
    show_help()
    return(invisible(NULL))
  }

  # Ensure fresh samples
  ensure_fresh_samples()

  # Connect to DB
  con <- dbConnect(duckdb(), dbdir = "pbp.duckdb")
  on.exit(dbDisconnect(con))

  # Validate column
  valid_columns <- get_valid_columns(con)
  validate_column(column, valid_columns)

  # Parse condition
  cond <- parse_condition(condition)

  # Validate period
  if (!is.null(period) && !period %in% VALID_PERIODS) {
    stop(sprintf("Invalid period '%s'. Valid: %s", period, paste(VALID_PERIODS, collapse = ", ")))
  }

  # Validate team
  if (!is.null(team) && !team %in% VALID_TEAMS) {
    stop(sprintf("Invalid team '%s'. Valid: %s", team, paste(VALID_TEAMS, collapse = ", ")))
  }

  # Get sample game IDs by joining to betting history tables (both pre-2020 and post-2020)
  samples <- dbGetQuery(con, "
    SELECT DISTINCT
      s.game_date, s.home_team, s.away_team,
      b.game_id as pbp_game_id
    FROM nfl_samples_temp s
    LEFT JOIN (
      SELECT game_id, game_date, home_team, away_team FROM nfl_betting_pbp
      UNION ALL
      SELECT game_id, game_date, home_team, away_team FROM nfl_pre_20_betting_history
    ) b
      ON s.game_date = b.game_date
      AND s.home_team = b.home_team
      AND s.away_team = b.away_team
    WHERE b.game_id IS NOT NULL
  ")

  if (nrow(samples) == 0) {
    stop("No samples found or could not match to PBP data.")
  }

  game_ids <- samples$pbp_game_id

  # Build the query
  period_filter <- get_period_filter(period)

  # Team filter - need to handle based on posteam or specific columns
  team_filter <- NULL
  if (!is.null(team)) {
    # For most stats, filter by posteam
    # But for things like td_team, we might need different logic
    team_filter <- sprintf("posteam = %s_team", team)
  }

  # Construct WHERE clause
  where_clauses <- c(sprintf("game_id IN ('%s')", paste(game_ids, collapse = "','")))
  if (!is.null(period_filter)) where_clauses <- c(where_clauses, period_filter)

  where_sql <- paste(where_clauses, collapse = " AND ")

  # Query to aggregate stat per game
  if (!is.null(team)) {
    # nfl_pbp has home_team/away_team as abbreviations, posteam is also abbreviation
    # Use a subquery to get the team abbreviation for each game, then filter
    query <- sprintf("
      SELECT
        p.game_id,
        %s(CASE WHEN p.posteam = p.%s_team THEN p.%s ELSE 0 END) as stat_value
      FROM nfl_pbp p
      WHERE %s
      GROUP BY p.game_id
    ", agg, team, column, where_sql)
  } else {
    query <- sprintf("
      SELECT
        game_id,
        %s(%s) as stat_value
      FROM nfl_pbp
      WHERE %s
      GROUP BY game_id
    ", agg, column, where_sql)
  }

  # Execute query
  results <- tryCatch({
    dbGetQuery(con, query)
  }, error = function(e) {
    stop(sprintf("Query error: %s", e$message))
  })

  # Handle games with no plays matching filter (stat_value = 0)
  # Left join with all sample game_ids to include games with 0
  all_games <- data.frame(game_id = game_ids)
  results <- all_games %>%
    left_join(results, by = "game_id")
  results$stat_value[is.na(results$stat_value)] <- 0

  # Evaluate condition
  results <- results %>%
    mutate(condition_met = case_when(
      cond$op == ">" ~ stat_value > cond$value,
      cond$op == ">=" ~ stat_value >= cond$value,
      cond$op == "<" ~ stat_value < cond$value,
      cond$op == "<=" ~ stat_value <= cond$value,
      cond$op == "==" ~ stat_value == cond$value
    ))

  # Calculate probability
  n_games <- nrow(results)
  n_hits <- sum(results$condition_met, na.rm = TRUE)
  prob <- n_hits / n_games

  # Convert to odds
  fair_american <- prob_to_american(prob)
  fair_decimal <- 1 / prob

  # Build result
  result <- list(
    column = column,
    condition = condition,
    period = if (is.null(period)) "Full" else period,
    team = if (is.null(team)) "Both" else team,
    agg = agg,
    prob = prob,
    fair_american_odds = fair_american,
    fair_decimal_odds = round(fair_decimal, 2),
    n_games = n_games,
    n_hits = n_hits,
    stat_distribution = summary(results$stat_value),
    raw_values = results$stat_value  # For histogram generation
  )

  # Print results
  print_prop_result(result)

  invisible(result)
}

# =============================================================================
# DISPLAY FUNCTIONS
# =============================================================================

print_prop_result <- function(result) {
  cat("\n=== PROP FAIR ODDS ===\n\n")

  # Build prop description
  prop_desc <- sprintf("%s(%s) %s", result$agg, result$column, result$condition)
  if (result$team != "Both") {
    prop_desc <- paste0(result$team, " ", prop_desc)
  }
  if (result$period != "Full") {
    prop_desc <- paste0(prop_desc, " (", result$period, ")")
  } else {
    prop_desc <- paste0(prop_desc, " (Full Game)")
  }

  cat(sprintf("Prop: %s\n\n", prop_desc))
  cat(sprintf("Fair probability: %.1f%%\n", result$prob * 100))
  if (is.finite(result$fair_american_odds)) {
    cat(sprintf("Fair American odds: %+d\n", as.integer(result$fair_american_odds)))
  } else {
    cat("Fair American odds: N/A (probability too extreme)\n")
  }
  if (is.finite(result$fair_decimal_odds)) {
    cat(sprintf("Fair decimal odds: %.2f\n", result$fair_decimal_odds))
  } else {
    cat("Fair decimal odds: N/A\n")
  }
  cat(sprintf("\nGames evaluated: %d\n", result$n_games))
  cat(sprintf("Condition met: %d (%.1f%%)\n", result$n_hits, result$prob * 100))

  cat("\nStat distribution:\n")
  print(result$stat_distribution)
}

show_help <- function() {
  cat("\n=== PROPS CALCULATOR ===\n\n")
  cat("Usage:\n")
  cat("  prop_odds(column, condition, period = NULL, team = NULL)\n\n")
  cat("Arguments:\n")
  cat("  column    - Column from nfl_pbp (e.g., 'touchdown')\n")
  cat("  condition - Operator and line (e.g., '> 4.5')\n")
  cat("  period    - Optional: '1', '2', '3', '4', 'Half1', 'Half2', 'Full'\n")
  cat("  team      - Optional: 'home', 'away'\n\n")
  cat("Examples:\n")
  cat("  prop_odds('touchdown', '> 4.5')                    # Total TDs over 4.5\n")
  cat("  prop_odds('pass_touchdown', '> 2.5', period='Half1') # 1H passing TDs\n")
  cat("  prop_odds('field_goal_attempt', '> 3.5', team='home') # Home FG attempts\n")
  cat("  prop_odds('safety', '>= 1')                        # Any safety\n\n")
  cat("Common columns:\n")
  cat(sprintf("  %s\n", paste(COMMON_COLUMNS, collapse = ", ")))
  cat("\nRun show_columns() to see all available columns.\n")
}

show_columns <- function(pattern = NULL) {
  con <- dbConnect(duckdb(), dbdir = "pbp.duckdb")
  cols <- get_valid_columns(con)
  dbDisconnect(con)

  if (!is.null(pattern)) {
    cols <- cols[grepl(pattern, cols, ignore.case = TRUE)]
  }

  cat(sprintf("\nAvailable columns (%d):\n", length(cols)))
  cat(paste(cols, collapse = "\n"))
  cat("\n")

  invisible(cols)
}

# Show help on source
cat("Props calculator loaded. Run prop_odds() for help.\n")
