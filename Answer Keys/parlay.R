#!/usr/bin/env Rscript
# Parlay Fair Odds Calculator
# Usage: Rscript parlay.R "1H home -3" "1H under 22.5" "1Q home ML"
#
# Leg format: "[period] [side] [line/market]"
#   Periods: 1Q, 2Q, 3Q, 4Q, 1H, 2H, FG (full game)
#
#   Spreads:    "1H home -3" or "1H away +3"
#   Totals:     "1H over 22.5" or "1H under 22.5"
#   Team Total: "1H home over 10.5" or "1H away under 7.5"
#   Moneyline:  "1H home ML" (2-way) or "1H home ML3" (3-way)
#   Tie:        "1Q tie" (3-way only)

suppressPackageStartupMessages({
  library(duckdb)
  library(dplyr)
  library(purrr)
})

setwd("~/NFLWork/Answer Keys")
source("Tools.R")

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  cat("Usage: Rscript parlay.R \"1H home -3\" \"1H under 22.5\" ...\n\n")
  cat("Leg formats:\n")
  cat("  Spread:     \"1H home -3\" or \"1H away +3\"\n")
  cat("  Total:      \"1H over 22.5\" or \"1H under 22.5\"\n")
  cat("  Team Total: \"1H home over 10.5\" or \"1H away under 7.5\"\n")
  cat("  Moneyline:  \"1H home ML\" (2-way) or \"1H home ML3\" (3-way)\n")
  cat("  Tie:        \"1Q tie\" (3-way)\n\n")
  cat("Periods: 1Q, 2Q, 3Q, 4Q, 1H, 2H, FG\n")
  quit(status = 0)
}

# Parse period string to internal format
parse_period <- function(p) {
  switch(toupper(p),
    "1Q" = "1",
    "2Q" = "2",
    "3Q" = "3",
    "4Q" = "4",
    "1H" = "Half1",
    "2H" = "Half2",
    "FG" = "Full",
    "FULL" = "Full",
    p  # fallback
 )
}

# Parse a leg string into a leg specification
parse_leg <- function(leg_str) {
  parts <- strsplit(trimws(leg_str), "\\s+")[[1]]

  if (length(parts) < 2) {
    stop(sprintf("Invalid leg format: '%s'", leg_str))
  }

  period <- parse_period(parts[1])

  # Check for tie (3-way ML)
  if (tolower(parts[2]) == "tie") {
    return(list(market = "moneyline_3way", period = period, side = "tie", line = NA))
  }

  # Check for over/under (game total)
  if (tolower(parts[2]) %in% c("over", "under")) {
    side <- tolower(parts[2])
    line <- as.numeric(parts[3])
    return(list(market = "total", period = period, side = side, line = line))
  }

  # Must be home/away
  side_raw <- tolower(parts[2])
  if (!side_raw %in% c("home", "away")) {
    stop(sprintf("Invalid side '%s' in leg: '%s'", side_raw, leg_str))
  }

  # Check what comes next
  if (length(parts) < 3) {
    stop(sprintf("Incomplete leg: '%s'", leg_str))
  }

  third <- parts[3]

  # Moneyline?
  if (toupper(third) == "ML") {
    return(list(market = "moneyline", period = period, side = side_raw, line = NA))
  }
  if (toupper(third) == "ML3") {
    return(list(market = "moneyline_3way", period = period, side = side_raw, line = NA))
  }

  # Team total? (home over 10.5, away under 7.5)
  if (tolower(third) %in% c("over", "under")) {
    direction <- tolower(third)
    line <- as.numeric(parts[4])
    team_side <- paste0(side_raw, "_", direction)
    return(list(market = "team_total", period = period, side = team_side, line = line))
  }

  # Must be spread (home -3, away +3)
  line <- as.numeric(third)
  return(list(market = "spread", period = period, side = side_raw, line = line))
}

# Parse all legs
legs <- tryCatch({
  lapply(args, parse_leg)
}, error = function(e) {
  cat("Error parsing legs:", e$message, "\n")
  quit(status = 1)
})

# =============================================================================
# CHECK SAMPLE FRESHNESS AND REGENERATE IF NEEDED
# =============================================================================

check_samples_fresh <- function(con, max_age_minutes = 5) {
  # Check if samples exist and are fresh (generated within max_age_minutes)
  tryCatch({
    meta <- dbGetQuery(con, "SELECT generated_at FROM nfl_samples_meta")
    if (nrow(meta) == 0) return(FALSE)

    generated_at <- as.POSIXct(meta$generated_at, tz = "UTC")
    age_minutes <- as.numeric(difftime(Sys.time(), generated_at, units = "mins"))

    return(age_minutes <= max_age_minutes)
  }, error = function(e) {
    return(FALSE)  # Table doesn't exist or other error
  })
}

generate_samples <- function() {
  cat("Generating fresh samples...\n")
  result <- system("Rscript 'NFL Answer Key/NFLPrepare.R'", intern = FALSE)
  if (result != 0) {
    cat("Error generating samples.\n")
    quit(status = 1)
  }
}

# Check freshness
con <- dbConnect(duckdb(), dbdir = "pbp.duckdb")
is_fresh <- check_samples_fresh(con, max_age_minutes = 5)
dbDisconnect(con)

if (!is_fresh) {
  generate_samples()
}

# Load samples
con <- dbConnect(duckdb(), dbdir = "pbp.duckdb")
samples <- dbGetQuery(con, "SELECT * FROM nfl_samples_temp")
dbDisconnect(con)

if (nrow(samples) == 0) {
  cat("No samples found after generation. Check NFLPrepare.R for errors.\n")
  quit(status = 1)
}

# Compute fair odds
result <- compute_parlay_fair_odds(samples, legs)

# Print results
print_parlay_result(result, legs)
