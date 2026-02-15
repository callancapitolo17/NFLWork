# Run CBB Derivative Odds Acquisition
# Starting from Nov 2023 when derivative data became available

setwd("~/NFLWork/Answer Keys/CBB Answer Key")
source("Acquire CBB Derivative Odds.R")

cat("\n\n========================================\n")
cat("STARTING CBB DERIVATIVE ODDS ACQUISITION\n")
cat("========================================\n\n")

# Start from Nov 2023 (derivative data available from May 2023)
start_date <- as.Date("2023-11-01")
end_date <- Sys.Date()

cat("Date range:", as.character(start_date), "to", as.character(end_date), "\n\n")

# Run acquisition with main markets only (no alternates to save quota)
result <- acquire_derivative_odds(
  markets = MAIN_MARKETS,
  start_date = start_date,
  end_date = end_date,
  max_games = NULL,  # All games
  append = FALSE
)

cat("\nAcquisition complete!\n")
