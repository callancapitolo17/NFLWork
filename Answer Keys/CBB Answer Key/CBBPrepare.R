# CBB Prepare - Part 1 of parallel pipeline
# Runs in parallel with scrapers
# Generates samples and saves to DuckDB for CBBCombine.R

.t_script_start <- Sys.time()
setwd("~/NFLWork/Answer Keys")
suppressPackageStartupMessages({
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
  library(hoopR)
})
source("Tools.R")
timer <- pipeline_timer()
# Backfill startup time (packages + source before timer existed)
startup_secs <- as.numeric(difftime(Sys.time(), .t_script_start, units = "secs"))
timer$mark(sprintf("r_startup (%.1fs total)", startup_secs))

cat("=== CBB PREPARE: Starting sample generation ===\n")

# =============================================================================
# LOAD HISTORICAL DATA
# =============================================================================

con <- dbConnect(duckdb(), dbdir = "cbb.duckdb")

# Load betting PBP with game results
betting_pbp <- dbGetQuery(con, "SELECT * FROM cbb_betting_pbp")

# Load closing odds and build consensus spread/total for historical games
# Group by game (odds_id) and compute weighted average across books
closing_odds <- dbGetQuery(con, "
  SELECT
    id as odds_id,
    home_spread,
    total_line,
    spread_home_odds,
    spread_away_odds,
    tot_over_odds,
    tot_under_odds,
    bookmaker_key
  FROM cbb_closing_odds
  WHERE home_spread IS NOT NULL
    AND total_line IS NOT NULL
")

# For historical data, compute simple consensus (median line, devigged probs)
historical_consensus <- closing_odds %>%
  group_by(odds_id) %>%
  summarize(
    home_spread = median(home_spread, na.rm = TRUE),
    total_line = median(total_line, na.rm = TRUE),
    # Average devigged probabilities across books
    .groups = "drop"
  ) %>%
  # Compute devigged probabilities from median odds

  mutate(
    consensus_devig_home_odds = 0.5,  # Spread markets ~50/50
    consensus_devig_away_odds = 0.5,
    consensus_devig_over_odds = 0.5,
    consensus_devig_under_odds = 0.5
  )

# Join PBP with consensus lines
DT <- betting_pbp %>%
  inner_join(historical_consensus, by = "odds_id") %>%
  # Rename columns to match Tools.R expectations
  rename(
    home_margin = game_home_margin_fg,
    total_final_score = game_total_fg,
    # Period columns - Tools.R expects game_home_margin_period_Half1 format
    game_home_margin_period_Half1 = game_home_margin_h1,
    game_home_margin_period_Half2 = game_home_margin_h2,
    game_total_period_Half1 = game_total_h1,
    game_total_period_Half2 = game_total_h2,
    home_score_period_Half1 = home_h1_score,
    home_score_period_Half2 = home_h2_score,
    away_score_period_Half1 = away_h1_score,
    away_score_period_Half2 = away_h2_score,
    home_spread_odds = consensus_devig_home_odds,
    away_spread_odds = consensus_devig_away_odds,
    over_odds = consensus_devig_over_odds,
    under_odds = consensus_devig_under_odds
  ) %>%
  mutate(
    actual_over = ifelse(total_final_score > total_line, 1, 0),
    actual_cover = ifelse(home_margin > -home_spread, 1, 0)
  ) %>%
  as.data.table()

cat(sprintf("Historical data loaded: %d games\n", nrow(DT)))

# Compute dispersion for distance weighting
disp <- compute_dispersion(DT, moneyline = FALSE)
ss <- disp$ss
st <- disp$st

# Create book weights table if it doesn't exist
# Use equal weights for now (can optimize later based on CLV analysis)
book_weights <- dbGetQuery(con, "
  SELECT DISTINCT bookmaker_key
  FROM cbb_closing_odds
") %>%
  mutate(
    spread_weight = 1.0,
    totals_weight = 1.0
  )

# Check if weights table exists, if not create it
tables <- dbListTables(con)
if (!"cbb_weights" %in% tables) {
  dbWriteTable(con, "cbb_weights", book_weights, overwrite = TRUE)
  cat("Created cbb_weights table with equal weights.\n")
} else {
  book_weights <- dbGetQuery(con, "SELECT * FROM cbb_weights")
}

# Close DB before forking (DuckDB doesn't support concurrent access)
dbDisconnect(con)
cat("Historical data loaded.\n")
timer$mark("historical_load")

# =============================================================================
# GET CURRENT GAME ODDS & BUILD CONSENSUS
# =============================================================================

cat("Fetching Odds API...\n")

game_odds <- tryCatch({
  toa_sports_odds(
    sport_key = "basketball_ncaab",
    regions = "us,us2,eu",
    markets = "spreads,totals",
    odds_format = "american",
    date_format = "iso"
  )
}, error = function(e) {
  cat(sprintf("Warning: Could not fetch odds: %s\n", e$message))
  return(data.frame())
})

cat(sprintf("API returned %d rows (%d unique games) before filtering.\n",
            nrow(game_odds), n_distinct(game_odds$id)))

# Parse commence_time from ISO 8601 character to proper POSIXct
# (API returns "2026-02-15T17:00:00Z" which R can't auto-coerce correctly)
if (is.character(game_odds$commence_time)) {
  game_odds$commence_time <- ymd_hms(game_odds$commence_time, tz = "UTC")
}

# Filter to only games that haven't started yet (no in-progress games)
n_before <- n_distinct(game_odds$id)
game_odds <- game_odds %>%
  filter(commence_time > Sys.time())
n_after <- n_distinct(game_odds$id)
cat(sprintf("Filtered in-progress: %d -> %d games (removed %d).\n", n_before, n_after, n_before - n_after))

if (nrow(game_odds) == 0) {
  cat("No upcoming games found from Odds API (filtered out in-progress). Exiting.\n")
  dbDisconnect(con)
  quit(status = 0)
}

# Build consensus spread
consensus_spread <- prepare_two_way_odds(
  game_odds    = game_odds,
  mkt_key      = "spreads",
  book_weights = book_weights,
  prob_fun     = devig_american,
  prob_names   = c("prob_home", "prob_away"),
  odds_names   = c("outcomes_price_home", "outcomes_price_away")
) %>%
  select(-outcomes_point_away) %>%
  rename(spread = outcomes_point_home) %>%
  pick_consensus_line(
    game_id_col = "id",
    line_col    = "spread",
    weight_col  = "spread_weight",
    date_col    = "date",
    time_col    = "commence_time",
    market1     = "prob_home",
    market2     = "prob_away"
  )
cat(sprintf("Consensus spreads: %d games.\n", n_distinct(consensus_spread$id)))

# Build consensus total
consensus_total <- prepare_two_way_odds(
  game_odds    = game_odds,
  mkt_key      = "totals",
  book_weights = book_weights,
  prob_fun     = devig_american,
  prob_names   = c("prob_over", "prob_under"),
  odds_names   = c("outcomes_price_Over", "outcomes_price_Under")
) %>%
  select(-outcomes_point_Under) %>%
  rename(total_line = outcomes_point_Over) %>%
  pick_consensus_line(
    game_id_col = "id",
    line_col    = "total_line",
    weight_col  = "totals_weight",
    date_col    = "date",
    time_col    = "commence_time",
    market1     = "prob_over",
    market2     = "prob_under"
  )
cat(sprintf("Consensus totals: %d games.\n", n_distinct(consensus_total$id)))

# Combine into single odds table
cbb_odds <- consensus_spread %>%
  inner_join(
    consensus_total %>% ungroup() %>% select(-home_team, -away_team, -date, -commence_time),
    by = "id"
  )
cat(sprintf("After inner_join (spread+total): %d games.\n", nrow(cbb_odds)))

cbb_odds <- cbb_odds %>%
  filter(if_all(everything(), ~ !is.na(.)))
cat(sprintf("After NA filter: %d games.\n", nrow(cbb_odds)))
timer$mark("consensus")

# =============================================================================
# BUILD TEAM NAME DICTIONARY (ESPN + Odds API)
# =============================================================================

cat("Building team name dictionary from ESPN...\n")
espn_teams <- tryCatch({
  espn_mbb_teams() %>%
    select(team_id, abbreviation, display_name, short_name, nickname, mascot)
}, error = function(e) {
  cat(sprintf("Warning: Could not fetch ESPN teams: %s\n", e$message))
  data.frame()
})

if (nrow(espn_teams) > 0) {
  odds_api_names <- unique(c(cbb_odds$home_team, cbb_odds$away_team))

  team_dict <- espn_teams %>%
    mutate(odds_api_name = NA_character_)

  for (oa_name in odds_api_names) {
    # Exact match on display_name
    idx <- which(tolower(team_dict$display_name) == tolower(oa_name))
    if (length(idx) == 0) {
      # nickname + mascot match
      idx <- which(sapply(seq_len(nrow(team_dict)), function(i) {
        grepl(team_dict$nickname[i], oa_name, ignore.case = TRUE) &&
          grepl(team_dict$mascot[i], oa_name, ignore.case = TRUE)
      }))
    }
    if (length(idx) == 1) {
      team_dict$odds_api_name[idx] <- oa_name
    }
  }

  matched <- sum(!is.na(team_dict$odds_api_name))
  cat(sprintf("Mapped %d/%d Odds API teams to ESPN dictionary.\n", matched, length(odds_api_names)))

  con_dict <- dbConnect(duckdb(), dbdir = "cbb.duckdb")
  dbWriteTable(con_dict, "cbb_team_dict", team_dict, overwrite = TRUE)
  dbDisconnect(con_dict)
}
timer$mark("espn_teams")

# =============================================================================
# SETUP PREDICTION PARAMETERS
# =============================================================================

# Targets for Answer Key algorithm
targets <- cbb_odds %>%
  transmute(
    id,
    parent_spread = spread,
    parent_total  = total_line,
    target_cover  = consensus_prob_home,
    target_over   = consensus_prob_over
  )

# Betting parameters (read from dashboard DB if saved, otherwise defaults)
bankroll   <- 100
kelly_mult <- 0.25
dash_db <- file.path(getwd(), "CBB Dashboard", "cbb_dashboard.duckdb")
if (file.exists(dash_db)) {
  tryCatch({
    dash_con <- dbConnect(duckdb(), dbdir = dash_db, read_only = TRUE)
    saved <- dbGetQuery(dash_con, "SELECT param, value FROM sizing_settings")
    dbDisconnect(dash_con)
    if ("bankroll" %in% saved$param) bankroll <- saved$value[saved$param == "bankroll"]
    if ("kelly_mult" %in% saved$param) kelly_mult <- saved$value[saved$param == "kelly_mult"]
  }, error = function(e) NULL)
}
cat(sprintf("Using bankroll=$%.0f, kelly=%.2f\n", bankroll, kelly_mult))
N <- round(nrow(DT) * 0.02, 0)  # 2% sample (optimized from parameter sweep)

# Get upcoming CBB events
events <- tryCatch({
  get_events("basketball_ncaab", regions = "us")
}, error = function(e) {
  cat(sprintf("Warning: Could not fetch events: %s\n", e$message))
  return(data.frame())
})

# =============================================================================
# GENERATE SAMPLES + FETCH DERIVATIVE ODDS (IN PARALLEL)
# =============================================================================

# All derivative markets to pre-fetch for CBBCombine
all_deriv_markets <- c(
  "h2h_h1", "h2h_h2",
  "spreads_h1", "spreads_h2", "alternate_spreads_h1", "alternate_spreads_h2",
  "totals_h1", "totals_h2", "alternate_totals_h1", "alternate_totals_h2",
  "team_totals_h1", "team_totals_h2", "alternate_team_totals_h1", "alternate_team_totals_h2"
)

# Fork child process to fetch derivative odds SIMULTANEOUSLY with sample generation
odds_job <- NULL
if (nrow(events) > 0) {
  n_combos <- nrow(events) * length(all_deriv_markets)
  cat(sprintf("Forking child process to fetch %d derivative odds (parallel with samples)...\n", n_combos))
  odds_job <- parallel::mcparallel({
    fetch_odds_bulk(events$id, all_deriv_markets, "basketball_ncaab")
  })
}

# Main process: generate samples (this takes ~5 min)
cat("Generating samples for all games...\n")
samples <- generate_all_samples(
  targets         = targets,
  DT              = DT,
  ss              = ss,
  st              = st,
  N               = N,
  use_spread_line = TRUE
)
cat(sprintf("Generated %d samples.\n", length(samples)))
timer$mark("sample_gen")

# Collect derivative odds from child process
prefetched_raw <- NULL
if (!is.null(odds_job)) {
  cat("Collecting pre-fetched derivative odds from child process...\n")
  prefetched_raw <- parallel::mccollect(odds_job)[[1]]
  prefetched_raw <- prefetched_raw[!is.na(prefetched_raw$json_response), ]
  cat(sprintf("Pre-fetched %d API responses (%d events x %d markets).\n",
              nrow(prefetched_raw), n_distinct(events$id), length(all_deriv_markets)))
}
timer$mark("prefetch_odds")

# =============================================================================
# SAVE TO DUCKDB (shared state, no temp files)
# =============================================================================

cat("Saving to DuckDB...\n")
con <- dbConnect(duckdb(), dbdir = "cbb.duckdb")

# Flatten samples list to dataframe for storage
samples_df <- samples %>%
  imap_dfr(~ .x$sample %>% mutate(game_id = .y))

# Save to DuckDB tables
dbWriteTable(con, "cbb_samples_temp", samples_df, overwrite = TRUE)
dbWriteTable(con, "cbb_odds_temp", cbb_odds, overwrite = TRUE)
if (nrow(events) > 0) {
  dbWriteTable(con, "cbb_events_temp", events, overwrite = TRUE)
}
if (!is.null(prefetched_raw) && nrow(prefetched_raw) > 0) {
  dbWriteTable(con, "cbb_prefetched_odds_temp", prefetched_raw, overwrite = TRUE)
  cat(sprintf("Saved %d pre-fetched odds to DuckDB.\n", nrow(prefetched_raw)))
}

# Also save betting parameters
params_df <- tibble(
  param = c("bankroll", "kelly_mult"),
  value = c(bankroll, kelly_mult)
)
dbWriteTable(con, "cbb_params_temp", params_df, overwrite = TRUE)

# Save generation timestamp for freshness checking
dbExecute(con, "CREATE OR REPLACE TABLE cbb_samples_meta AS SELECT CURRENT_TIMESTAMP as generated_at")

timer$mark("duckdb_save")

dbDisconnect(con)

cat("=== CBB PREPARE: Complete. Samples saved to DuckDB. ===\n")
