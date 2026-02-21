# CBB Combine - Part 2 of parallel pipeline
# Runs after scrapers and CBBPrepare.R complete
# Loads samples from DuckDB, merges with scraped odds, finds edge

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
})
source("Tools.R")

cat("=== CBB COMBINE: Loading data and finding edge ===\n")

# =============================================================================
# LOAD PREPARED DATA FROM DUCKDB
# =============================================================================

con <- dbConnect(duckdb(), dbdir = "cbb.duckdb")

samples_df <- dbGetQuery(con, "SELECT * FROM cbb_samples_temp")
cbb_odds <- dbGetQuery(con, "SELECT * FROM cbb_odds_temp")

# Events table might not exist if API failed
events <- tryCatch({
  dbGetQuery(con, "SELECT * FROM cbb_events_temp")
}, error = function(e) {
  data.frame()
})

params <- dbGetQuery(con, "SELECT * FROM cbb_params_temp")

# Load pre-fetched derivative odds (cached JSON from parallel fetch in CBBPrepare)
prefetched_odds <- tryCatch({
  dbGetQuery(con, "SELECT * FROM cbb_prefetched_odds_temp")
}, error = function(e) {
  cat("No pre-fetched odds found, will fetch live from API.\n")
  NULL
})

dbDisconnect(con)

# Extract betting parameters
bankroll <- params %>% filter(param == "bankroll") %>% pull(value)
kelly_mult <- params %>% filter(param == "kelly_mult") %>% pull(value)

# Read enabled books from dashboard settings (controls API books only)
dashboard_db <- file.path(getwd(), "CBB Dashboard", "cbb_dashboard.duckdb")
if (file.exists(dashboard_db)) {
  dash_con <- dbConnect(duckdb(), dbdir = dashboard_db)
  # Seed offshore books as enabled by default
  for (book in c("wagerzon", "hoop88", "bfa")) {
    dbExecute(dash_con, "INSERT INTO book_settings (bookmaker_key, enabled)
      VALUES (?, TRUE) ON CONFLICT (bookmaker_key) DO NOTHING", list(book))
  }
  enabled_books <- tryCatch({
    dbGetQuery(dash_con, "SELECT bookmaker_key FROM book_settings WHERE enabled = TRUE")$bookmaker_key
  }, error = function(e) NULL)
  dbDisconnect(dash_con)
  if (length(enabled_books) == 0) enabled_books <- NULL
  cat(sprintf("Enabled books: %s\n", if (is.null(enabled_books)) "ALL (none configured)" else paste(enabled_books, collapse = ", ")))
} else {
  enabled_books <- NULL
  cat("No dashboard DB found, using all books.\n")
}

cat(sprintf("Loaded samples for %d games.\n", n_distinct(samples_df$game_id)))

# Reconstruct samples list from dataframe
game_ids <- unique(samples_df$game_id)
samples <- map(game_ids, function(gid) {
  df <- samples_df %>% filter(game_id == gid) %>% select(-game_id)
  list(sample = df)
}) %>%
  set_names(game_ids)

# =============================================================================
# LOAD SCRAPED BOOKMAKER ODDS
# =============================================================================

wagerzon_odds <- get_wagerzon_odds("cbb")
cat(sprintf("Loaded %d Wagerzon records.\n", nrow(wagerzon_odds)))

hoop88_odds <- get_hoop88_odds("cbb")
cat(sprintf("Loaded %d Hoop88 records.\n", nrow(hoop88_odds)))

bfa_odds <- get_bfa_odds("cbb")
cat(sprintf("Loaded %d BFA records.\n", nrow(bfa_odds)))

# =============================================================================
# GENERATE PREDICTIONS FOR CBB MARKETS
# =============================================================================

# CBB has 2 halves (not 4 quarters)
periods_base <- c("Half1", "Half2")

# --- MONEYLINES (H1 + H2) ---
cat("Building moneyline predictions...\n")
ml_results <- build_moneylines_from_samples(
  samples          = samples,
  consensus_odds   = cbb_odds,
  events           = events,
  periods          = periods_base,
  markets          = c("h2h_h1", "h2h_h2"),
  sport_key        = "basketball_ncaab",
  bankroll         = bankroll,
  kelly_mult       = kelly_mult,
  margin_col       = "game_home_margin_period",
  pre_fetched_odds = prefetched_odds,
  books            = enabled_books
)

ml_bets <- ml_results$bets

# --- TOTALS (H1 + H2) + ALTERNATE TOTALS ---
cat("Building totals predictions...\n")
totals_markets <- c(
  "totals_h1", "totals_h2",
  "alternate_totals_h1", "alternate_totals_h2"
)
totals_periods <- c(
  "Half1", "Half2",
  "Half1", "Half2"
)

total_results <- build_totals_from_samples(
  samples          = samples,
  consensus_odds   = cbb_odds,
  events           = events,
  periods          = totals_periods,
  markets          = totals_markets,
  sport_key        = "basketball_ncaab",
  bankroll         = bankroll,
  kelly_mult       = kelly_mult,
  pre_fetched_odds = prefetched_odds,
  books            = enabled_books
)

total_bets <- total_results$bets

# --- SPREADS (H1 + H2) + ALTERNATE SPREADS ---
cat("Building spread predictions...\n")
spreads_markets <- c(
  "spreads_h1", "spreads_h2",
  "alternate_spreads_h1", "alternate_spreads_h2"
)
spreads_periods <- c(
  "Half1", "Half2",
  "Half1", "Half2"
)

spread_results <- build_spreads_from_samples(
  samples          = samples,
  consensus_odds   = cbb_odds,
  events           = events,
  periods          = spreads_periods,
  markets          = spreads_markets,
  sport_key        = "basketball_ncaab",
  bankroll         = bankroll,
  kelly_mult       = kelly_mult,
  pre_fetched_odds = prefetched_odds,
  books            = enabled_books
)

spread_bets <- spread_results$bets

# --- TEAM TOTALS (H1 + H2) + ALTERNATE TEAM TOTALS ---
cat("Building team totals predictions...\n")
team_totals_markets <- c(
  "team_totals_h1", "team_totals_h2",
  "alternate_team_totals_h1", "alternate_team_totals_h2"
)
team_totals_periods <- c(
  "Half1", "Half2",
  "Half1", "Half2"
)

team_totals_results <- build_team_totals_from_samples(
  samples          = samples,
  consensus_odds   = cbb_odds,
  events           = events,
  periods          = team_totals_periods,
  markets          = team_totals_markets,
  sport_key        = "basketball_ncaab",
  bankroll         = bankroll,
  kelly_mult       = kelly_mult,
  pre_fetched_odds = prefetched_odds,
  books            = enabled_books
)

team_totals_bets <- team_totals_results$bets

# =============================================================================
# COMBINE ALL API BETS
# =============================================================================

EV_THRESHOLD <- 0.05  # 5% minimum EV (optimized from parameter sweep)

all_bets <- bind_rows(
  ml_bets %>% mutate(market_type = "moneyline"),
  total_bets %>% mutate(market_type = "totals"),
  spread_bets %>% mutate(market_type = "spreads"),
  team_totals_bets %>% mutate(market_type = "team_totals")
) %>%
  filter(ev >= EV_THRESHOLD) %>%  # Filter to 5%+ EV bets only
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

if (nrow(wagerzon_odds) > 0) {
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

  # Compute EV for BFA alt lines directly from samples (no API dependency)
  bfa_alt_bets <- compare_alts_to_samples(
    samples = samples,
    offshore_odds = bfa_odds,
    consensus_odds = cbb_odds,
    bankroll = bankroll,
    kelly_mult = kelly_mult,
    ev_threshold = 0.05
  )
  if (nrow(bfa_alt_bets) > 0) {
    bfa_alt_bets <- bfa_alt_bets %>%
      mutate(market_type = ifelse(grepl("spread", market), "spreads", "totals"))
  }
  cat(sprintf("Added %d BFA alt line bets from samples.\n", nrow(bfa_alt_bets)))
} else {
  bfa_bets <- tibble()
  bfa_alt_bets <- tibble()
}

# =============================================================================
# COMBINE ALL BETS (API + OFFSHORE)
# =============================================================================

all_bets_combined <- bind_rows(
  all_bets,
  wagerzon_bets,
  hoop88_bets,
  bfa_bets,
  bfa_alt_bets
) %>%
  # Filter out games that have already started
  filter(is.na(pt_start_time) | pt_start_time > Sys.time()) %>%
  # Collapse main + alt into one base market (e.g., alternate_spreads_h1 -> spreads_h1)
  mutate(base_market = gsub("^alternate_", "", market)) %>%
  # Global dedup: keep best EV per (game, base_market, side) across all books/lines
  group_by(id, base_market, bet_on) %>%
  filter(ev == max(ev)) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(-base_market) %>%
  arrange(desc(ev))

# Final summary by market type
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

# =============================================================================
# SAVE TO DUCKDB FOR DASHBOARD
# =============================================================================

con <- dbConnect(duckdb(), dbdir = "cbb.duckdb")

dbExecute(con, "DROP TABLE IF EXISTS cbb_bets_combined")
dbWriteTable(con, "cbb_bets_combined", all_bets_combined)

dbDisconnect(con)

cat(sprintf("Saved %d bets to cbb_bets_combined table.\n", nrow(all_bets_combined)))
cat("\n=== CBB COMBINE: Complete ===\n")
