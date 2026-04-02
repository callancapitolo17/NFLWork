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
  library(digest)
  library(lubridate)
})

setwd("~/NFLWork/Answer Keys")
source("Tools.R")

# =============================================================================
# CONFIG
# =============================================================================

EV_THRESHOLD <- 0.02   # 2% minimum edge to flag (default, overridden by dashboard)
KELLY_EDGE_MIN <- 0.03 # 3% minimum edge to enter conditional Kelly sizing
WZ_PARLAY_SHAVE <- 0.989  # Wagerzon takes ~1.1% off independent multiply
MLB_DB <- "mlb.duckdb"

# Load sizing from dashboard if available (parlay-specific settings, falls back to main)
bankroll   <- 4000
kelly_mult <- 0.25
dash_db <- "MLB Dashboard/mlb_dashboard.duckdb"
if (file.exists(dash_db)) {
  tryCatch({
    dash_con <- dbConnect(duckdb(), dbdir = dash_db, read_only = TRUE)
    saved <- dbGetQuery(dash_con, "SELECT param, value FROM sizing_settings")
    dbDisconnect(dash_con)
    # Prefer parlay-specific settings, fall back to main bankroll/kelly
    if ("parlay_bankroll" %in% saved$param) {
      bankroll <- saved$value[saved$param == "parlay_bankroll"]
    } else if ("bankroll" %in% saved$param) {
      bankroll <- saved$value[saved$param == "bankroll"]
    }
    if ("parlay_kelly_mult" %in% saved$param) {
      kelly_mult <- saved$value[saved$param == "parlay_kelly_mult"]
    } else if ("kelly_mult" %in% saved$param) {
      kelly_mult <- saved$value[saved$param == "kelly_mult"]
    }
    if ("parlay_min_edge" %in% saved$param) {
      EV_THRESHOLD <- saved$value[saved$param == "parlay_min_edge"] / 100
    }
  }, error = function(e) NULL)
}

# Helper: American odds → decimal odds
american_to_decimal <- function(odds) {
  ifelse(odds > 0, 1 + odds / 100, 1 + 100 / abs(odds))
}

# =============================================================================
# CONDITIONAL KELLY FOR CORRELATED PARLAYS
# =============================================================================
#
# Each parlay is treated as a single "bet" that wins when both legs hit.
# We build a binary outcome matrix (n_sims x n_parlays) from game samples,
# compute the covariance of returns across parlays, and solve the multivariate
# Kelly: f* = Sigma^{-1} * mu. This adjusts sizes downward when parlays on
# the same game are positively correlated (e.g., Fav+Over and Dog+Under tend
# to win/lose in the same game states).
#
# Same math as multivariate_kelly() in Tools.R, but operates on parlay-level
# outcomes (both legs ANDed) rather than single-leg outcomes.

parlay_multivariate_kelly <- function(parlay_group, samples, bankroll, kelly_mult) {
  # parlay_group: list of lists, each with $legs, $joint_prob, $wz_dec
  # samples: game sample data frame (one game)
  # Returns: numeric vector of dollar wager amounts (one per parlay)

  n <- length(parlay_group)

  # --- Build parlay outcome matrix: 1 = both legs hit, 0 = miss, NA = push ---
  outcome_matrix <- matrix(NA_real_, nrow = nrow(samples), ncol = n)
  for (i in seq_len(n)) {
    legs <- parlay_group[[i]]$legs
    # Evaluate each leg, AND them together
    leg_results <- map(legs, ~evaluate_leg(samples, .x))
    # Parlay wins only when ALL legs resolve AND all hit
    all_resolved <- map(leg_results, ~!is.na(.x)) %>% reduce(`&`)
    all_hit <- map(leg_results, ~.x == TRUE) %>% reduce(`&`)
    # Mark unresolved (push on any leg) as NA
    outcomes <- ifelse(all_resolved, as.numeric(all_hit), NA_real_)
    outcome_matrix[, i] <- outcomes
  }

  # Drop rows with any NA (push on any parlay in the group)
  complete_rows <- complete.cases(outcome_matrix)
  if (sum(complete_rows) < 30) {
    # Too few resolved samples — fall back to independent Kelly
    return(independent_kelly(parlay_group, bankroll, kelly_mult))
  }
  outcome_matrix <- outcome_matrix[complete_rows, , drop = FALSE]

  # --- Correlation matrix from binary outcomes ---
  R <- cor(outcome_matrix)
  if (any(is.na(R))) {
    return(independent_kelly(parlay_group, bankroll, kelly_mult))
  }

  # --- Covariance of returns ---
  # sigma_i = sqrt(p_i * (1 - p_i)) * (b_i + 1)
  # where p_i = fair joint prob, b_i = wz decimal odds - 1
  probs <- sapply(parlay_group, `[[`, "joint_prob")
  b <- sapply(parlay_group, `[[`, "wz_dec") - 1
  sigmas <- sqrt(probs * (1 - probs)) * (b + 1)

  Sigma <- R * outer(sigmas, sigmas)
  # Ridge regularization to prevent singularity
  Sigma <- Sigma + diag(n) * 0.01

  if (kappa(Sigma) > 100) {
    return(fallback_rho_scaling(parlay_group, R, bankroll, kelly_mult))
  }

  # --- EV vector: expected return per dollar wagered ---
  mu <- probs * (b + 1) - 1  # p * wz_dec - 1

  # --- Solve f* = Sigma^{-1} * mu ---
  f_star <- tryCatch(solve(Sigma, mu), error = function(e) NULL)
  if (is.null(f_star)) {
    return(fallback_rho_scaling(parlay_group, R, bankroll, kelly_mult))
  }

  # Clamp negatives to 0 (Kelly says don't bet these)
  f_star <- pmax(f_star, 0)

  # Convert fractions to dollar wagers
  wagers <- f_star * kelly_mult * bankroll
  round(pmax(wagers, 0))
}


# Fallback: scale independent Kelly by average correlation
fallback_rho_scaling <- function(parlay_group, R, bankroll, kelly_mult) {
  n <- length(parlay_group)
  wagers <- numeric(n)
  for (i in seq_len(n)) {
    p <- parlay_group[[i]]$joint_prob
    b <- parlay_group[[i]]$wz_dec - 1
    kelly_full <- (b * p - (1 - p)) / b
    if (kelly_full <= 0) next

    # Average absolute correlation with other parlays in the group
    rho_vals <- abs(R[i, -i])
    avg_rho <- mean(rho_vals)
    # Scale down: more correlated = smaller bet
    scale <- 1 / sqrt(1 + (n - 1) * avg_rho)
    wagers[i] <- round(kelly_full * scale * kelly_mult * bankroll)
  }
  pmax(wagers, 0)
}


# Independent Kelly: no correlation adjustment (single parlay or fallback)
independent_kelly <- function(parlay_group, bankroll, kelly_mult) {
  wagers <- numeric(length(parlay_group))
  for (i in seq_along(parlay_group)) {
    p <- parlay_group[[i]]$joint_prob
    b <- parlay_group[[i]]$wz_dec - 1
    kelly_full <- (b * p - (1 - p)) / b
    wagers[i] <- max(0, round(kelly_full * kelly_mult * bankroll))
  }
  wagers
}


# Nudge wager so Wagerzon "to win" rounds UP (0.50 rounds DOWN on WZ)
nudge_wager_rounds_up <- function(wager, b) {
  if (wager <= 0) return(wager)
  rounds_up <- function(w) (w * b) %% 1 > 0.50
  if (!rounds_up(wager)) {
    if (rounds_up(wager + 1)) {
      wager <- wager + 1
    } else if (wager > 1 && rounds_up(wager - 1)) {
      wager <- wager - 1
    }
  }
  wager
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
    # MLB.R writes local time via format(Sys.time()), but DuckDB TIMESTAMP
    # assumes UTC. Re-parse as local to compare correctly against Sys.time().
    generated_at <- as.POSIXct(format(meta$generated_at), tz = Sys.timezone())
    age_minutes <- as.numeric(difftime(Sys.time(), generated_at, units = "mins"))
    return(age_minutes <= max_age_minutes)
  }, error = function(e) {
    return(FALSE)
  })
}

if (!check_mlb_samples_fresh(max_age_minutes = 10)) {
  cat("Samples stale or missing. Run the MLB pipeline first (run.py mlb).\n")
  quit(status = 0)
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

# Join spreads + totals by (home_team, away_team, game_date) into one row per game.
# game_date is "MM/DD" from Wagerzon's own API field — used later to match the
# correct consensus game when teams play on consecutive days (series).
wz_combined <- wz_spreads %>%
  select(home_team, away_team, game_date,
         home_spread, home_spread_price = odds_home,
         away_spread, away_spread_price = odds_away) %>%
  inner_join(
    wz_totals %>% select(home_team, away_team, game_date,
                          total_line = line,
                          over_price = odds_over,
                          under_price = odds_under),
    by = c("home_team", "away_team", "game_date")
  )

cat(sprintf("Wagerzon FG games with spread + total: %d\n", nrow(wz_combined)))

# Try to load exact parlay prices from ConfirmWagerHelper API
wz_db <- normalizePath(path.expand("~/NFLWork/wagerzon_odds/wagerzon.duckdb"), mustWork = FALSE)
exact_prices <- tryCatch({
  con_wz <- dbConnect(duckdb(), dbdir = wz_db, read_only = TRUE)
  on.exit(dbDisconnect(con_wz, shutdown = TRUE))
  pp <- dbGetQuery(con_wz, "SELECT * FROM mlb_parlay_prices")
  dbDisconnect(con_wz, shutdown = TRUE)
  on.exit(NULL)
  pp
}, error = function(e) data.frame())

use_exact <- nrow(exact_prices) > 0
if (use_exact) {
  cat(sprintf("Loaded %d exact parlay prices from ConfirmWagerHelper API.\n", nrow(exact_prices)))
} else {
  cat("No exact parlay prices found. Using shave+round approximation.\n")
}

# =============================================================================
# MATCH WAGERZON GAMES TO SAMPLES
# =============================================================================

# Add a game_date column to consensus ("MM/DD" in Eastern time) so we can
# join on date and avoid fan-out when teams play on consecutive days in a series.
consensus_dated <- consensus %>%
  select(id, home_team, away_team, consensus_prob_home, any_of("commence_time")) %>%
  mutate(game_date = if ("commence_time" %in% names(.)) {
    format(with_tz(as.POSIXct(commence_time, tz = "UTC"), "America/New_York"), "%m/%d")
  } else NA_character_)

wz_matched <- wz_combined %>%
  inner_join(consensus_dated, by = c("home_team", "away_team", "game_date"))

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
legs_store <- list()  # stash legs for conditional Kelly pass, keyed by game_id|combo

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

    # Wagerzon's parlay price — use exact API price if available, else shave+round
    combo_name <- combo$name
    exact_row <- if (use_exact) {
      exact_prices %>% filter(home_team == row$home_team, combo == combo_name)
    } else {
      data.frame()
    }

    if (nrow(exact_row) > 0) {
      wz_dec <- exact_row$wz_decimal[1]
      wz_american <- exact_row$wz_american[1]
    } else {
      # Fallback: multiply legs independently, apply 1.1% shave, round to whole American
      wz_dec_raw <- american_to_decimal(combo$wz_spread_price) *
                    american_to_decimal(combo$wz_total_price)
      wz_dec_shaved <- wz_dec_raw * WZ_PARLAY_SHAVE
      wz_american <- round(prob_to_american(1 / wz_dec_shaved))
      wz_dec <- american_to_decimal(wz_american)
    }

    # Edge: positive means Wagerzon overpays relative to fair
    fair_dec <- fair$fair_decimal_odds
    edge_pct <- (wz_dec - fair_dec) / fair_dec * 100

    combo_spread       <- if (grepl("Home", combo_name)) row$home_spread else row$away_spread
    combo_spread_price <- combo$wz_spread_price
    combo_total_price  <- combo$wz_total_price

    # Store WITHOUT sizing — conditional Kelly sizes in pass 2
    results[[length(results) + 1]] <- tibble(
      game_id     = game_id,
      game        = sprintf("%s @ %s", row$away_team, row$home_team),
      home_team   = row$home_team,
      away_team   = row$away_team,
      game_time   = if ("commence_time" %in% names(row) && !is.na(row$commence_time)) {
        format(as.POSIXct(row$commence_time, tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
      } else NA_character_,
      combo        = combo_name,
      spread_line  = combo_spread,
      total_line   = row$total_line,
      spread_price = combo_spread_price,
      total_price  = combo_total_price,
      fair_odds   = fair$fair_american_odds,
      wz_odds     = wz_american,
      fair_dec    = round(fair_dec, 3),
      wz_dec      = round(wz_dec, 3),
      corr_factor = fair$correlation_factor,
      edge_pct    = round(edge_pct, 1),
      joint_prob_raw = fair$joint_prob,  # unrounded, needed for Kelly
      joint_prob  = round(fair$joint_prob * 100, 1),
      n_samples   = fair$n_samples_resolved
    )

    # Stash legs for conditional Kelly pass
    legs_store[[paste0(game_id, "|", combo_name)]] <- combo$legs
  }
}

# =============================================================================
# PASS 2: CONDITIONAL KELLY SIZING
# =============================================================================
#
# Group parlays by game. For each game:
#   - Filter to parlays above KELLY_EDGE_MIN
#   - 0 qualify → $0 wager
#   - 1 qualifies → standard single-bet Kelly
#   - 2+ qualify → multivariate Kelly accounting for parlay-level correlation

if (length(results) == 0) {
  cat("No parlay combos could be computed. Check data availability.\n")
  quit(status = 0)
}

all_results <- bind_rows(results)
all_results$kelly_bet <- 0  # initialize, fill below

for (gid in unique(all_results$game_id)) {
  game_mask <- all_results$game_id == gid
  game_rows <- all_results[game_mask, ]

  # Filter to parlays above the Kelly edge threshold
  qualifies <- game_rows$edge_pct >= KELLY_EDGE_MIN * 100
  if (sum(qualifies) == 0) next  # no edges worth sizing

  qual_rows <- game_rows[qualifies, ]
  samp <- samples_list[[gid]]

  # Build the parlay group list that parlay_multivariate_kelly() expects
  parlay_group <- list()
  for (j in seq_len(nrow(qual_rows))) {
    r <- qual_rows[j, ]
    leg_key <- paste0(gid, "|", r$combo)
    parlay_group[[j]] <- list(
      legs = legs_store[[leg_key]],
      joint_prob = r$joint_prob_raw,
      wz_dec = r$wz_dec
    )
  }

  # Size with conditional Kelly (multivariate if 2+, independent if 1)
  if (length(parlay_group) >= 2) {
    wagers <- parlay_multivariate_kelly(parlay_group, samp, bankroll, kelly_mult)
  } else {
    wagers <- independent_kelly(parlay_group, bankroll, kelly_mult)
  }

  # Apply Wagerzon "rounds up" nudge to each wager
  for (j in seq_len(length(wagers))) {
    b <- parlay_group[[j]]$wz_dec - 1
    wagers[j] <- nudge_wager_rounds_up(wagers[j], b)
  }

  # Write sized wagers back into the qualifying rows
  qual_indices <- which(game_mask)[qualifies]
  all_results$kelly_bet[qual_indices] <- wagers
}

all_results <- all_results %>% arrange(desc(edge_pct))

# =============================================================================
# OUTPUT
# =============================================================================

cat("\n")
cat("=================================================================\n")
cat("  MLB CORRELATED PARLAY ANALYSIS (Conditional Kelly)\n")
cat(sprintf("  Bankroll: $%.0f | Kelly: %.0f%% | Edge Min: %.0f%% | Flag: %.0f%%\n",
            bankroll, kelly_mult * 100, KELLY_EDGE_MIN * 100, EV_THRESHOLD * 100))
cat("=================================================================\n\n")

# Show all combos grouped by game
for (game_name in unique(all_results$game)) {
  game_rows <- all_results %>% filter(game == game_name)
  cat(sprintf("--- %s ---\n", game_name))
  cat(sprintf("  %-22s  %6s  %6s  %5s  %6s  %6s  %6s\n",
              "Combo", "Fair", "WZ", "Corr", "Edge%", "Wager", "ToWin"))

  for (j in seq_len(nrow(game_rows))) {
    r <- game_rows[j, ]
    to_win <- round(r$kelly_bet * (r$wz_dec - 1))
    edge_flag <- if (r$edge_pct >= KELLY_EDGE_MIN * 100) " ***" else ""
    cat(sprintf("  %-22s  %+6d  %+6d  %5.3f  %+5.1f%%  $%4d  $%4d%s\n",
                r$combo, r$fair_odds, r$wz_odds, r$corr_factor,
                r$edge_pct, round(r$kelly_bet), to_win, edge_flag))
  }
  cat("\n")
}

# Summary of edges (only sized parlays — above KELLY_EDGE_MIN)
edges <- all_results %>% filter(kelly_bet > 0)
if (nrow(edges) > 0) {
  total_wager <- sum(round(edges$kelly_bet))
  total_to_win <- sum(round(edges$kelly_bet * (edges$wz_dec - 1)))
  cat(sprintf("=== %d EDGES SIZED (>= %.0f%% EV) | Total Wager: $%d | Total To Win: $%d ===\n\n",
              nrow(edges), KELLY_EDGE_MIN * 100, total_wager, total_to_win))
  for (j in seq_len(nrow(edges))) {
    e <- edges[j, ]
    wager <- round(e$kelly_bet)
    to_win <- round(wager * (e$wz_dec - 1))
    cat(sprintf("  %s | %s\n", e$game, e$combo))
    cat(sprintf("    Spread: %+.1f | Total: %.1f\n", e$spread_line, e$total_line))
    cat(sprintf("    Fair: %+d (%.1f%%) | WZ: %+d | Correlation: %.3f\n",
                e$fair_odds, e$joint_prob, e$wz_odds, e$corr_factor))
    cat(sprintf("    Edge: %+.1f%% | Wager: $%d | To Win: $%d\n\n",
                e$edge_pct, wager, to_win))
  }
} else {
  cat("No edges found above threshold.\n")
  cat("(Correlation exists but may not overcome Wagerzon's vig on these lines.)\n")
}

# =============================================================================
# WRITE TO DUCKDB (for dashboard consumption)
# =============================================================================

# Add parlay hash for dedup in dashboard
all_results <- all_results %>%
  mutate(parlay_hash = pmap_chr(list(game_id, combo), function(gid, cmb) {
    digest(paste(gid, cmb, sep = "|"), algo = "sha256", serialize = FALSE)
  }))

write_con <- NULL
tryCatch({
  write_con <- dbConnect(duckdb(), dbdir = MLB_DB)
  dbExecute(write_con, "DROP TABLE IF EXISTS mlb_parlay_opportunities")
  dbWriteTable(write_con, "mlb_parlay_opportunities", all_results)
  cat(sprintf("Wrote %d parlay opportunities to %s.\n", nrow(all_results), MLB_DB))
}, error = function(e) {
  cat(sprintf("Warning: Failed to write parlays to DB: %s\n", e$message))
})
if (!is.null(write_con)) dbDisconnect(write_con)
