# NFL Combine - Part 2 of parallel pipeline
# Runs after scrapers and NFLPrepare.R complete
# Loads samples from DuckDB, merges with scraped odds, finds edge

setwd("~/NFLWork/Answer Keys")
library(data.table)
library(oddsapiR)
library(duckdb)
library(dplyr)
library(purrr)
library(lubridate)
library(DBI)
library(httr)
library(jsonlite)
library(tidyverse)
source("Tools.R")

cat("=== NFL COMBINE: Loading data and finding edge ===\n")

# =============================================================================
# LOAD PREPARED DATA FROM DUCKDB
# =============================================================================

con <- dbConnect(duckdb(), dbdir = "pbp.duckdb")

samples_df <- dbGetQuery(con, "SELECT * FROM nfl_samples_temp")
nfl_odds <- dbGetQuery(con, "SELECT * FROM nfl_odds_temp")
events <- dbGetQuery(con, "SELECT * FROM nfl_events_temp")
params <- dbGetQuery(con, "SELECT * FROM nfl_params_temp")

dbDisconnect(con)

# Extract betting parameters
bankroll <- params %>% filter(param == "bankroll") %>% pull(value)
kelly_mult <- params %>% filter(param == "kelly_mult") %>% pull(value)

cat(sprintf("Loaded samples for %d games.\n", n_distinct(samples_df$game_id)))

# Reconstruct samples list from dataframe
# The original samples structure is: samples[[game_id]]$sample = dataframe
# Group by game_id, remove game_id col, wrap in list with $sample
game_ids <- unique(samples_df$game_id)
samples <- map(game_ids, function(gid) {
  df <- samples_df %>% filter(game_id == gid) %>% select(-game_id)
  list(sample = df)
}) %>%
  set_names(game_ids)

# =============================================================================
# LOAD SCRAPED BOOKMAKER ODDS
# =============================================================================

wagerzon_odds <- get_wagerzon_odds("nfl")
cat(sprintf("Loaded %d Wagerzon records.\n", nrow(wagerzon_odds)))

hoop88_odds <- get_hoop88_odds("nfl")
cat(sprintf("Loaded %d Hoop88 records.\n", nrow(hoop88_odds)))

bfa_odds <- get_bfa_odds("nfl")
cat(sprintf("Loaded %d BFA records.\n", nrow(bfa_odds)))

# =============================================================================
# GENERATE PREDICTIONS FOR ALL DERIVATIVE MARKETS
# =============================================================================

# Base periods (Q1-Q4 + H1 + H2)
periods_base <- c("1", "2", "3", "4", "Half1", "Half2")

# --- MONEYLINES (Q1-Q4 + H1 + H2) ---
cat("Building moneyline predictions...\n")
ml_results <- build_moneylines_from_samples(
  samples         = samples,
  consensus_odds  = nfl_odds,
  events          = events,
  periods         = periods_base,
  markets         = c("h2h_q1", "h2h_q2", "h2h_q3", "h2h_q4", "h2h_h1", "h2h_h2"),
  sport_key       = "americanfootball_nfl",
  bankroll        = bankroll,
  kelly_mult      = kelly_mult,
  margin_col      = "game_home_margin_period"
)

ml_bets <- ml_results$bets

# --- 3-WAY MONEYLINES (Q1-Q4 + H1 + H2) ---
cat("Building 3-way moneyline predictions...\n")
ml_3way_results <- build_3way_from_samples(
  samples         = samples,
  consensus_odds  = nfl_odds,
  events          = events,
  periods         = periods_base,
  markets         = c("h2h_3_way_q1", "h2h_3_way_q2", "h2h_3_way_q3", "h2h_3_way_q4", "h2h_3_way_h1", "h2h_3_way_h2"),
  sport_key       = "americanfootball_nfl",
  bankroll        = bankroll,
  kelly_mult      = kelly_mult,
  margin_col      = "game_home_margin_period"
)

ml_3way_bets <- ml_3way_results$bets

# --- TOTALS (Q1-Q4 + H1 + H2) + ALTERNATE TOTALS ---
cat("Building totals predictions...\n")
totals_markets <- c(
  "totals_q1", "totals_q2", "totals_q3", "totals_q4", "totals_h1", "totals_h2",
  "alternate_totals_q1", "alternate_totals_q2", "alternate_totals_q3", "alternate_totals_q4", "alternate_totals_h1", "alternate_totals_h2"
)
totals_periods <- c(
  "1", "2", "3", "4", "Half1", "Half2",
  "1", "2", "3", "4", "Half1", "Half2"
)

total_results <- build_totals_from_samples(
  samples         = samples,
  consensus_odds  = nfl_odds,
  events          = events,
  periods         = totals_periods,
  markets         = totals_markets,
  sport_key       = "americanfootball_nfl",
  bankroll        = bankroll,
  kelly_mult      = kelly_mult
)

total_bets <- total_results$bets

# --- SPREADS (Q1-Q4 + H1 + H2) + ALTERNATE SPREADS ---
cat("Building spread predictions...\n")
spreads_markets <- c(
  "spreads_q1", "spreads_q2", "spreads_q3", "spreads_q4", "spreads_h1", "spreads_h2",
  "alternate_spreads_q1", "alternate_spreads_q2", "alternate_spreads_q3", "alternate_spreads_q4", "alternate_spreads_h1", "alternate_spreads_h2"
)
spreads_periods <- c(
  "1", "2", "3", "4", "Half1", "Half2",
  "1", "2", "3", "4", "Half1", "Half2"
)

spread_results <- build_spreads_from_samples(
  samples         = samples,
  consensus_odds  = nfl_odds,
  events          = events,
  periods         = spreads_periods,
  markets         = spreads_markets,
  sport_key       = "americanfootball_nfl",
  bankroll        = bankroll,
  kelly_mult      = kelly_mult
)

spread_bets <- spread_results$bets

# --- TEAM TOTALS (Q1-Q4 + H1 + H2) + ALTERNATE TEAM TOTALS ---
cat("Building team totals predictions...\n")
team_totals_markets <- c(
  "team_totals_q1", "team_totals_q2", "team_totals_q3", "team_totals_q4",
  "team_totals_h1", "team_totals_h2",
  "alternate_team_totals_q1", "alternate_team_totals_q2", "alternate_team_totals_q3", "alternate_team_totals_q4",
  "alternate_team_totals_h1", "alternate_team_totals_h2"
)
team_totals_periods <- c(
  "1", "2", "3", "4", "Half1", "Half2",
  "1", "2", "3", "4", "Half1", "Half2"
)

team_totals_results <- build_team_totals_from_samples(
  samples         = samples,
  consensus_odds  = nfl_odds,
  events          = events,
  periods         = team_totals_periods,
  markets         = team_totals_markets,
  sport_key       = "americanfootball_nfl",
  bankroll        = bankroll,
  kelly_mult      = kelly_mult
)

team_totals_bets <- team_totals_results$bets

# =============================================================================
# COMBINE ALL API BETS
# =============================================================================

all_bets <- bind_rows(
  ml_bets %>% mutate(market_type = "moneyline"),
  ml_3way_bets %>% mutate(market_type = "moneyline_3way"),
  total_bets %>% mutate(market_type = "totals"),
  spread_bets %>% mutate(market_type = "spreads"),
  team_totals_bets %>% mutate(market_type = "team_totals")
) %>%
  arrange(desc(ev))

cat("\n=== API BETS SUMMARY ===\n")
all_bets %>%
  group_by(market_type) %>%
  summarise(
    n_bets = n(),
    total_stake = sum(bet_size),
    avg_ev = mean(ev),
    .groups = "drop"
  ) %>%
  print()

# =============================================================================
# ADD WAGERZON AS BOOKMAKER TO PREDICTIONS
# =============================================================================

# Wagerzon bets are generated through the same prediction logic,
# just compared against Wagerzon's odds instead of API bookmakers
if (nrow(wagerzon_odds) > 0) {
  # Silently compare to Wagerzon (suppress warnings and messages)
  invisible(capture.output(suppressWarnings({
    wz_spread_bets <- compare_spreads_to_wagerzon(
      spread_results, wagerzon_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
    wz_total_bets <- compare_totals_to_wagerzon(
      total_results, wagerzon_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
    wz_ml_bets <- compare_moneylines_to_wagerzon(
      ml_results, wagerzon_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
  })))

  # Combine Wagerzon bets
  wagerzon_bets <- bind_rows(
    wz_spread_bets$bets %>% mutate(market_type = "spreads"),
    wz_total_bets$bets %>% mutate(market_type = "totals"),
    wz_ml_bets$bets %>% mutate(market_type = "moneyline")
  )

  cat(sprintf("Added %d Wagerzon bets to predictions.\n", nrow(wagerzon_bets)))
} else {
  wagerzon_bets <- tibble()
}

# =============================================================================
# ADD HOOP88 AS BOOKMAKER TO PREDICTIONS
# =============================================================================

if (nrow(hoop88_odds) > 0) {
  invisible(capture.output(suppressWarnings({
    h88_spread_bets <- compare_spreads_to_wagerzon(
      spread_results, hoop88_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
    h88_total_bets <- compare_totals_to_wagerzon(
      total_results, hoop88_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
    h88_ml_bets <- compare_moneylines_to_wagerzon(
      ml_results, hoop88_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
  })))

  hoop88_bets <- bind_rows(
    h88_spread_bets$bets %>% mutate(market_type = "spreads"),
    h88_total_bets$bets %>% mutate(market_type = "totals"),
    h88_ml_bets$bets %>% mutate(market_type = "moneyline")
  )

  cat(sprintf("Added %d Hoop88 bets to predictions.\n", nrow(hoop88_bets)))
} else {
  hoop88_bets <- tibble()
}

# =============================================================================
# ADD BFA AS BOOKMAKER TO PREDICTIONS
# =============================================================================

if (nrow(bfa_odds) > 0) {
  invisible(capture.output(suppressWarnings({
    bfa_spread_bets <- compare_spreads_to_wagerzon(
      spread_results, bfa_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
    bfa_total_bets <- compare_totals_to_wagerzon(
      total_results, bfa_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
    bfa_ml_bets <- compare_moneylines_to_wagerzon(
      ml_results, bfa_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
  })))

  bfa_bets <- bind_rows(
    bfa_spread_bets$bets %>% mutate(market_type = "spreads"),
    bfa_total_bets$bets %>% mutate(market_type = "totals"),
    bfa_ml_bets$bets %>% mutate(market_type = "moneyline")
  )

  cat(sprintf("Added %d BFA bets to predictions.\n", nrow(bfa_bets)))
} else {
  bfa_bets <- tibble()
}

# =============================================================================
# COMBINE ALL BETS (API + OFFSHORE as unified output)
# =============================================================================

all_bets_combined <- bind_rows(
  all_bets,
  wagerzon_bets,
  hoop88_bets,
  bfa_bets
) %>%
  arrange(desc(ev))

# Final summary by market type (not by source)
cat("\n=== BETTING SUMMARY ===\n")
all_bets_combined %>%
  group_by(market_type) %>%
  summarise(
    n_bets = n(),
    total_stake = sum(bet_size),
    avg_ev = mean(ev),
    .groups = "drop"
  ) %>%
  arrange(desc(total_stake)) %>%
  print()

# Top bets across all bookmakers
cat("\n=== TOP 20 BETS ===\n")
print(all_bets_combined %>% head(20))

# Open interactive table in browser (temp file auto-cleaned by OS)
if (nrow(all_bets_combined) > 0) {
  library(DT)

  # Format for display
  display_df <- all_bets_combined %>%
    select(bookmaker_key, market, bet_on, line, ev, prob, odds, bet_size, home_team, away_team) %>%
    mutate(
      ev = round(ev * 100, 1),  # Convert to percentage
      prob = round(prob * 100, 1),
      bet_size = round(bet_size, 2)
    ) %>%
    rename(
      Book = bookmaker_key,
      Market = market,
      Bet = bet_on,
      Line = line,
      `EV%` = ev,
      `Prob%` = prob,
      Odds = odds,
      `Bet$` = bet_size,
      Home = home_team,
      Away = away_team
    )

  # Create interactive datatable
  dt <- datatable(
    display_df,
    caption = paste("NFL Bets -", Sys.Date()),
    filter = "top",
    options = list(
      pageLength = 50,
      order = list(list(4, "desc")),  # Sort by EV% descending
      dom = "Bfrtip"
    )
  ) %>%
    formatStyle("EV%", backgroundColor = styleInterval(c(5, 10), c("white", "lightgreen", "green")))

  # Save to temp and open in browser
  tmp <- tempfile(fileext = ".html")
  htmlwidgets::saveWidget(dt, tmp, selfcontained = TRUE)

  # Normalize path and open with file:// URL scheme
  tmp_path <- normalizePath(tmp, mustWork = FALSE)
  if (file.exists(tmp_path)) {
    browseURL(paste0("file://", tmp_path))
    Sys.sleep(2)  # Wait for browser to load before R exits
    cat("\nOpened bets in browser.\n")
  } else {
    cat("\nWarning: Could not create HTML file.\n")
  }
}

# =============================================================================
# SAVE TO DUCKDB FOR DASHBOARD
# =============================================================================

con <- dbConnect(duckdb(), dbdir = "pbp.duckdb")

# Save all_bets_combined for dashboard to read
dbExecute(con, "DROP TABLE IF EXISTS nfl_bets_combined")
dbWriteTable(con, "nfl_bets_combined", all_bets_combined)

dbDisconnect(con)

cat(sprintf("Saved %d bets to nfl_bets_combined table.\n", nrow(all_bets_combined)))
cat("\n=== NFL COMBINE: Complete ===\n")
