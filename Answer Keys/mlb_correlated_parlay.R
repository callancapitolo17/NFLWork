#!/usr/bin/env Rscript
# MLB Correlated Parlay Edge Finder
#
# Finds +EV spread+total parlays on Wagerzon where correlation creates edge.
# Wagerzon prices these independently (multiplies odds), but spread and total
# are correlated:
#   - Favorite + Over: positively correlated
#   - Underdog + Under: positively correlated
#
# Uses the same sampling engine and compute_parlay_fair_odds() from Tools.R
# that powers the NFL parlay calculator.
#
# Usage: Rscript mlb_correlated_parlay.R

suppressPackageStartupMessages({
  library(duckdb)
  library(dplyr)
  library(purrr)
})

setwd("~/NFLWork/Answer Keys")
source("Tools.R")

# =============================================================================
# CONFIG
# =============================================================================

EV_THRESHOLD <- 0.02   # 2% minimum edge to flag
MLB_DB <- "mlb.duckdb"

# Load sizing from dashboard if available
bankroll   <- 100
kelly_mult <- 0.25
dash_db <- "MLB Dashboard/mlb_dashboard.duckdb"
if (file.exists(dash_db)) {
  tryCatch({
    dash_con <- dbConnect(duckdb(), dbdir = dash_db, read_only = TRUE)
    saved <- dbGetQuery(dash_con, "SELECT param, value FROM sizing_settings")
    dbDisconnect(dash_con)
    if ("bankroll" %in% saved$param) bankroll <- saved$value[saved$param == "bankroll"]
    if ("kelly_mult" %in% saved$param) kelly_mult <- saved$value[saved$param == "kelly_mult"]
  }, error = function(e) NULL)
}

# Helper: American odds → decimal odds
american_to_decimal <- function(odds) {
  ifelse(odds > 0, 1 + odds / 100, 1 + 100 / abs(odds))
}

# =============================================================================
# CHECK SAMPLE FRESHNESS (same pattern as parlay.R)
# =============================================================================

check_mlb_samples_fresh <- function(max_age_minutes = 5) {
  tryCatch({
    con <- dbConnect(duckdb(), dbdir = MLB_DB, read_only = TRUE)
    on.exit(dbDisconnect(con))
    meta <- dbGetQuery(con, "SELECT generated_at FROM mlb_samples_meta")
    if (nrow(meta) == 0) return(FALSE)
    generated_at <- as.POSIXct(meta$generated_at, tz = "UTC")
    age_minutes <- as.numeric(difftime(Sys.time(), generated_at, units = "mins"))
    return(age_minutes <= max_age_minutes)
  }, error = function(e) {
    return(FALSE)
  })
}

if (!check_mlb_samples_fresh(max_age_minutes = 5)) {
  cat("Samples stale or missing. Regenerating via MLB.R...\n")
  result <- system("Rscript 'MLB Answer Key/MLB.R'", intern = FALSE)
  if (result != 0) {
    cat("Error running MLB.R pipeline.\n")
    quit(status = 1)
  }
}

# =============================================================================
# LOAD SAMPLES (same pattern as parlay.R:155-162)
# =============================================================================

con <- dbConnect(duckdb(), dbdir = MLB_DB, read_only = TRUE)
on.exit(dbDisconnect(con))

samples_df <- dbGetQuery(con, "SELECT * FROM mlb_game_samples")
consensus  <- dbGetQuery(con, "SELECT * FROM mlb_consensus_temp")
dbDisconnect(con)
on.exit(NULL)

if (nrow(samples_df) == 0) {
  cat("No samples found in mlb_game_samples. Run MLB.R first.\n")
  quit(status = 0)
}

# Split into per-game list of data frames
samples_list <- split(samples_df, samples_df$game_id)
cat(sprintf("Loaded %d samples across %d games.\n",
            nrow(samples_df), length(samples_list)))

# =============================================================================
# LOAD WAGERZON FG ODDS
# =============================================================================

wz <- get_wagerzon_odds("mlb")

if (nrow(wz) == 0) {
  cat("No Wagerzon MLB odds found. Run scraper first.\n")
  quit(status = 0)
}

# Filter to full game MAIN line only (exclude alternates and H1)
wz_spreads <- wz %>% filter(period == "fg", market == "spreads")
wz_totals  <- wz %>% filter(period == "fg", market == "totals")

if (nrow(wz_spreads) == 0 || nrow(wz_totals) == 0) {
  cat(sprintf("Wagerzon FG data: %d spreads, %d totals. Need both.\n",
              nrow(wz_spreads), nrow(wz_totals)))
  quit(status = 0)
}

# Join spreads + totals by (home_team, away_team) into one row per game
wz_combined <- wz_spreads %>%
  select(home_team, away_team,
         home_spread, home_spread_price = odds_home,
         away_spread, away_spread_price = odds_away) %>%
  inner_join(
    wz_totals %>% select(home_team, away_team,
                          total_line = line,
                          over_price = odds_over,
                          under_price = odds_under),
    by = c("home_team", "away_team")
  )

cat(sprintf("Wagerzon FG games with spread + total: %d\n", nrow(wz_combined)))

# =============================================================================
# MATCH WAGERZON GAMES TO SAMPLES
# =============================================================================

wz_matched <- wz_combined %>%
  inner_join(
    consensus %>% select(id, home_team, away_team, consensus_prob_home),
    by = c("home_team", "away_team")
  )

if (nrow(wz_matched) == 0) {
  cat("No Wagerzon games matched to consensus. Check team name resolution.\n")
  cat("Wagerzon teams:\n")
  for (i in seq_len(nrow(wz_combined))) {
    cat(sprintf("  %s @ %s\n", wz_combined$away_team[i], wz_combined$home_team[i]))
  }
  cat("Consensus teams:\n")
  for (i in seq_len(nrow(consensus))) {
    cat(sprintf("  %s @ %s\n", consensus$away_team[i], consensus$home_team[i]))
  }
  quit(status = 0)
}

cat(sprintf("Matched %d games between Wagerzon and samples.\n", nrow(wz_matched)))

# =============================================================================
# COMPUTE PARLAY FAIR ODDS FOR EACH GAME
# =============================================================================

results <- list()

for (i in seq_len(nrow(wz_matched))) {
  row <- wz_matched[i, ]
  game_id <- row$id

  samp <- samples_list[[game_id]]
  if (is.null(samp)) {
    cat(sprintf("  Skipping %s @ %s: no samples for game_id %s\n",
                row$away_team, row$home_team, game_id))
    next
  }

  # Build 4 parlay combos: each spread side × each total side
  combos <- list(
    list(
      name = "Home Spread + Over",
      legs = list(
        list(market = "spread", period = "Full", side = "home", line = row$home_spread),
        list(market = "total",  period = "Full", side = "over", line = row$total_line)
      ),
      wz_spread_price = row$home_spread_price,
      wz_total_price  = row$over_price
    ),
    list(
      name = "Home Spread + Under",
      legs = list(
        list(market = "spread", period = "Full", side = "home", line = row$home_spread),
        list(market = "total",  period = "Full", side = "under", line = row$total_line)
      ),
      wz_spread_price = row$home_spread_price,
      wz_total_price  = row$under_price
    ),
    list(
      name = "Away Spread + Over",
      legs = list(
        list(market = "spread", period = "Full", side = "away", line = row$away_spread),
        list(market = "total",  period = "Full", side = "over", line = row$total_line)
      ),
      wz_spread_price = row$away_spread_price,
      wz_total_price  = row$over_price
    ),
    list(
      name = "Away Spread + Under",
      legs = list(
        list(market = "spread", period = "Full", side = "away", line = row$away_spread),
        list(market = "total",  period = "Full", side = "under", line = row$total_line)
      ),
      wz_spread_price = row$away_spread_price,
      wz_total_price  = row$under_price
    )
  )

  for (combo in combos) {
    # Skip if either price is missing
    if (is.na(combo$wz_spread_price) || is.na(combo$wz_total_price)) next

    # Our fair price (correlation-adjusted via historical samples)
    fair <- compute_parlay_fair_odds(samp, combo$legs)

    # Wagerzon's price (independent: multiply the two leg decimals)
    wz_dec <- american_to_decimal(combo$wz_spread_price) *
              american_to_decimal(combo$wz_total_price)
    wz_american <- prob_to_american(1 / wz_dec)

    # Edge: positive means Wagerzon overpays relative to fair
    fair_dec <- fair$fair_decimal_odds
    edge_pct <- (wz_dec - fair_dec) / fair_dec * 100

    # Kelly sizing
    p <- fair$joint_prob
    b <- wz_dec - 1
    kelly_full <- (b * p - (1 - p)) / b
    kelly_bet  <- max(0, round(kelly_full * kelly_mult * bankroll, 2))

    combo_name <- combo$name
    combo_spread <- if (grepl("Home", combo_name)) row$home_spread else row$away_spread

    results[[length(results) + 1]] <- tibble(
      game        = sprintf("%s @ %s", row$away_team, row$home_team),
      combo       = combo_name,
      spread_line = combo_spread,
      total_line  = row$total_line,
      fair_odds   = fair$fair_american_odds,
      wz_odds     = wz_american,
      fair_dec    = round(fair_dec, 3),
      wz_dec      = round(wz_dec, 3),
      corr_factor = fair$correlation_factor,
      edge_pct    = round(edge_pct, 1),
      kelly_bet   = kelly_bet,
      joint_prob  = round(fair$joint_prob * 100, 1),
      n_samples   = fair$n_samples_resolved
    )
  }
}

# =============================================================================
# OUTPUT
# =============================================================================

if (length(results) == 0) {
  cat("No parlay combos could be computed. Check data availability.\n")
  quit(status = 0)
}

all_results <- bind_rows(results) %>% arrange(desc(edge_pct))

cat("\n")
cat("=================================================================\n")
cat("  MLB CORRELATED PARLAY ANALYSIS\n")
cat(sprintf("  Bankroll: $%.0f | Kelly: %.0f%% | EV Threshold: %.0f%%\n",
            bankroll, kelly_mult * 100, EV_THRESHOLD * 100))
cat("=================================================================\n\n")

# Show all combos grouped by game
for (game_name in unique(all_results$game)) {
  game_rows <- all_results %>% filter(game == game_name)
  cat(sprintf("--- %s ---\n", game_name))
  cat(sprintf("  %-22s  %6s  %6s  %5s  %6s  %6s\n",
              "Combo", "Fair", "WZ", "Corr", "Edge%", "Kelly$"))

  for (j in seq_len(nrow(game_rows))) {
    r <- game_rows[j, ]
    edge_flag <- if (r$edge_pct >= EV_THRESHOLD * 100) " ***" else ""
    cat(sprintf("  %-22s  %+6d  %+6d  %5.3f  %+5.1f%%  $%5.2f%s\n",
                r$combo, r$fair_odds, r$wz_odds, r$corr_factor,
                r$edge_pct, r$kelly_bet, edge_flag))
  }
  cat("\n")
}

# Summary of edges
edges <- all_results %>% filter(edge_pct >= EV_THRESHOLD * 100)
if (nrow(edges) > 0) {
  cat(sprintf("=== %d EDGES FOUND (>= %.0f%% EV) ===\n\n", nrow(edges), EV_THRESHOLD * 100))
  for (j in seq_len(nrow(edges))) {
    e <- edges[j, ]
    cat(sprintf("  %s | %s\n", e$game, e$combo))
    cat(sprintf("    Spread: %+.1f | Total: %.1f\n", e$spread_line, e$total_line))
    cat(sprintf("    Fair: %+d (%.1f%%) | WZ: %+d | Correlation: %.3f\n",
                e$fair_odds, e$joint_prob, e$wz_odds, e$corr_factor))
    cat(sprintf("    Edge: %+.1f%% | Kelly bet: $%.2f\n\n", e$edge_pct, e$kelly_bet))
  }
} else {
  cat("No edges found above threshold.\n")
  cat("(Correlation exists but may not overcome Wagerzon's vig on these lines.)\n")
}
