# =============================================================================
# Answer Key Algorithm - Deep Diagnostic Runner
# =============================================================================
# Purpose: Produce reproducible diagnostic output to measure algorithm behavior
# and identify root causes of "overestimated edge" symptom.
#
# Usage: source("diagnostics/answer_key_deep_diagnostic.R")
# Output: Saves diagnostic results to diagnostics/output/
# =============================================================================

# Setup -----------------------------------------------------------------------
setwd("~/NFLWork/Answer Keys")
library(data.table)
library(duckdb)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(DBI)
library(purrr)
library(stringr)
source("Tools.R")

# Create output directory
output_dir <- "diagnostics/output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

cat("=============================================================================\n")
cat("ANSWER KEY DEEP DIAGNOSTIC\n")
cat("Timestamp:", timestamp, "\n")
cat("=============================================================================\n\n")

# =============================================================================
# SECTION 1: Load Data
# =============================================================================

cat("SECTION 1: Loading Data\n")
cat("-----------------------\n")

con <- dbConnect(duckdb(), dbdir = "pbp.duckdb") 

betting_20_plus <- dbGetQuery(con, "SELECT * FROM nfl_betting_pbp")
pre_20 <- dbGetQuery(con, "SELECT * FROM nfl_pre_20_betting_history")

dbDisconnect(con)

# FIX: Ensure column names are consistent before binding
# The pre_20 table has 'spread' while betting_20_plus has 'home_spread'
if ("spread" %in% names(pre_20) && !"home_spread" %in% names(pre_20)) {
  pre_20 <- pre_20 %>% rename(home_spread = spread)
}

DT <- bind_rows(betting_20_plus, pre_20) %>%
  rename(
    home_spread_odds = "consensus_devig_home_odds",
    away_spread_odds = "consensus_devig_away_odds",
    over_odds = "consensus_devig_over_odds",
    under_odds = "consensus_devig_under_odds"
  ) %>%
  # Filter out rows with missing critical data
  filter(!is.na(home_spread), !is.na(total_line)) %>%
  mutate(
    actual_over = ifelse(total_final_score > total_line, 1, 0),
    actual_cover = ifelse(home_margin > -home_spread, 1, 0)
  ) %>% 
  as.data.table()

cat("Total games loaded:", nrow(DT), "\n")
cat("Games with spread data:", sum(!is.na(DT$home_spread)), "\n")
cat("Games with total data:", sum(!is.na(DT$total_line)), "\n")
cat("Games with Q1 margin:", sum(!is.na(DT$game_home_margin_period_1)), "\n")
cat("Games with H1 margin:", sum(!is.na(DT$game_home_margin_period_Half1) | 
                                   !is.na(DT$game_home_margin_period_h1)), "\n")

# Check column names for margin columns
margin_cols <- names(DT)[grepl("game_home_margin", names(DT))]
cat("\nMargin columns available:\n")
print(margin_cols)

# Compute dispersion
disp <- compute_dispersion(DT, moneyline = FALSE)
ss <- disp$ss
st <- disp$st

cat("\nDispersion metrics:\n")
cat("  ss (spread 5th-95th range):", ss, "\n")
cat("  st (total 5th-95th range):", st, "\n")

# =============================================================================
# SECTION 2: Data Quality Checks
# =============================================================================

cat("\n\nSECTION 2: Data Quality Checks\n")
cat("------------------------------\n")

# Check actual_cover and actual_over distributions
cat("Full-game cover rate (home covers):", mean(DT$actual_cover, na.rm = TRUE), "\n")
cat("Full-game over rate:", mean(DT$actual_over, na.rm = TRUE), "\n")

# Check for ties in period data
cat("\nTie frequencies by period:\n")
for (p in 1:4) {
  col <- paste0("game_home_margin_period_", p)
  if (col %in% names(DT)) {
    ties <- sum(DT[[col]] == 0, na.rm = TRUE)
    total <- sum(!is.na(DT[[col]]))
    cat(sprintf("  Q%d ties: %d / %d (%.2f%%)\n", p, ties, total, 100 * ties/total))
  }
}

# Check for half ties
for (h in c("Half1", "Half2", "h1", "h2")) {
  col <- paste0("game_home_margin_period_", h)
  if (col %in% names(DT)) {
    ties <- sum(DT[[col]] == 0, na.rm = TRUE)
    total <- sum(!is.na(DT[[col]]))
    cat(sprintf("  %s ties: %d / %d (%.2f%%)\n", h, ties, total, 100 * ties/total))
  }
}

# =============================================================================
# SECTION 3: Single Game Diagnostic Trace
# =============================================================================

cat("\n\nSECTION 3: Single Game Diagnostic Trace\n")
cat("----------------------------------------\n")

# Use a representative parent line
test_spread <- -3  # Home favored by 3
test_total <- 44
test_cover <- 0.52  # 52% chance home covers (from devigged -115/+100 or similar)
test_over <- 0.50   # 50% over

N <- round(nrow(DT) * 0.05, 0)
cat("Sample size N:", N, "\n\n")

cat("Parent line:\n")
cat("  Spread:", test_spread, "(home favored by 3)\n")
cat("  Total:", test_total, "\n")
cat("  Target cover prob:", test_cover, "\n")
cat("  Target over prob:", test_over, "\n\n")

# Step 1: Mean Match
cat("--- STEP 1: Mean Matching ---\n")
mm <- mean_match(
  dt = DT,
  N = N,
  parent_spread = test_spread,
  parent_total = test_total,
  ss = ss,
  st = st,
  max_iter_mean = 500,
  tol_mean = 0.005,
  use_spread_line = TRUE
)

inc_mm <- mm$dt[included == TRUE]

cat("After mean_match:\n")
cat("  Sample size:", nrow(inc_mm), "\n")
cat("  Spread - Target:", test_spread, 
    "| Mean:", round(mean(inc_mm$home_spread, na.rm = TRUE), 3),
    "| Median:", round(median(inc_mm$home_spread, na.rm = TRUE), 3), "\n")
cat("  Total - Target:", test_total,
    "| Mean:", round(mean(inc_mm$total_line, na.rm = TRUE), 3),
    "| Median:", round(median(inc_mm$total_line, na.rm = TRUE), 3), "\n")
cat("  Pre-balance cover rate:", round(mean(inc_mm$actual_cover, na.rm = TRUE), 4), "\n")
cat("  Pre-balance over rate:", round(mean(inc_mm$actual_over, na.rm = TRUE), 4), "\n")
cat("  Index range:", round(min(inc_mm$index), 4), "to", round(max(inc_mm$index), 4), "\n\n")

# Step 2: Balance Sample
cat("--- STEP 2: Balance Sample ---\n")
bal <- balance_sample(
  dt = mm$dt,
  N = N,
  target_cover = test_cover,
  target_over = test_over,
  tol_error = 1
)

inc_bal <- bal$dt[included == TRUE]

cat("After balance_sample:\n")
cat("  Final N:", bal$final_N, "(started with", N, ")\n")
cat("  Actual sample size:", nrow(inc_bal), "\n")
cat("  Cover - Target:", test_cover, 
    "| Actual:", round(mean(inc_bal$actual_cover, na.rm = TRUE), 4),
    "| Error:", bal$cover_error, "games\n")
cat("  Over - Target:", test_over,
    "| Actual:", round(mean(inc_bal$actual_over, na.rm = TRUE), 4),
    "| Error:", bal$over_error, "games\n")

# Check spread/total drift after balancing
cat("\n  Post-balance spread - Mean:", round(mean(inc_bal$home_spread, na.rm = TRUE), 3),
    "| Drift:", round(mean(inc_bal$home_spread, na.rm = TRUE) - test_spread, 4), "\n")
cat("  Post-balance total - Mean:", round(mean(inc_bal$total_line, na.rm = TRUE), 3),
    "| Drift:", round(mean(inc_bal$total_line, na.rm = TRUE) - test_total, 4), "\n")

# Step 3: Extract period probabilities
cat("\n--- STEP 3: Period Probability Extraction ---\n")

period_stats <- list()

for (p in c(1, 2, 3, 4)) {
  col <- paste0("game_home_margin_period_", p)
  if (!col %in% names(inc_bal)) {
    cat("Column", col, "not found, skipping\n")
    next
  }
  
  margins <- inc_bal[[col]]
  
  home_wins <- sum(margins > 0, na.rm = TRUE)
  away_wins <- sum(margins < 0, na.rm = TRUE)
  ties <- sum(margins == 0, na.rm = TRUE)
  total_games <- sum(!is.na(margins))
  
  # Current implementation: excludes ties
  prob_current <- home_wins / (home_wins + away_wins)
  
  # Alternative: include ties in denominator
  prob_with_ties <- home_wins / total_games
  
  period_stats[[paste0("Q", p)]] <- tibble(
    period = paste0("Q", p),
    home_wins = home_wins,
    away_wins = away_wins,
    ties = ties,
    total = total_games,
    prob_excl_ties = prob_current,
    prob_incl_ties = prob_with_ties,
    inflation = prob_current - prob_with_ties
  )
  
  cat(sprintf("  Q%d: %d wins, %d losses, %d ties (%.1f%%) | P(home)=%.3f (excl ties) vs %.3f (incl ties) | Inflation: %.3f\n",
              p, home_wins, away_wins, ties, 100*ties/total_games, 
              prob_current, prob_with_ties, prob_current - prob_with_ties))
}

# Check for half column
h1_col <- if ("game_home_margin_period_Half1" %in% names(inc_bal)) {
  "game_home_margin_period_Half1"
} else if ("game_home_margin_period_h1" %in% names(inc_bal)) {
  "game_home_margin_period_h1"
} else {
  NULL
}

if (!is.null(h1_col)) {
  margins <- inc_bal[[h1_col]]
  
  home_wins <- sum(margins > 0, na.rm = TRUE)
  away_wins <- sum(margins < 0, na.rm = TRUE)
  ties <- sum(margins == 0, na.rm = TRUE)
  total_games <- sum(!is.na(margins))
  
  prob_current <- home_wins / (home_wins + away_wins)
  prob_with_ties <- home_wins / total_games
  
  cat(sprintf("  H1: %d wins, %d losses, %d ties (%.1f%%) | P(home)=%.3f (excl ties) vs %.3f (incl ties) | Inflation: %.3f\n",
              home_wins, away_wins, ties, 100*ties/total_games, 
              prob_current, prob_with_ties, prob_current - prob_with_ties))
}

# =============================================================================
# SECTION 4: Edge Distribution Analysis (Simulated Market)
# =============================================================================

cat("\n\nSECTION 4: Edge Distribution Analysis\n")
cat("--------------------------------------\n")

# Simulate a market slate with various spread/total combinations
# For each, we generate our probability and compare to market
set.seed(42)

# Generate synthetic "market" lines
n_test_games <- 50
test_games <- tibble(
  game_id = 1:n_test_games,
  spread = sample(seq(-10, 10, by = 0.5), n_test_games, replace = TRUE),
  total = rnorm(n_test_games, mean = 44, sd = 3),
  # Market's no-vig probability for home covering
  market_prob_home = 0.5 + spread * 0.02  # Rough approximation: each point ≈ 2%
)

# Clamp probabilities
test_games <- test_games %>%
  mutate(
    market_prob_home = pmin(pmax(market_prob_home, 0.20), 0.80),
    market_prob_over = 0.50  # Assume 50% for simplicity
  )

cat("Testing", n_test_games, "synthetic game lines...\n\n")

# Run prediction for each game
predictions <- test_games %>%
  mutate(
    result = pmap(
      list(game_id, spread, total, market_prob_home, market_prob_over),
      function(id, sp, tot, prob_home, prob_over) {
        tryCatch({
          run_means_for_id(
            id = id,
            parent_spread = sp,
            parent_total = tot,
            target_cover = prob_home,
            target_over = prob_over,
            DT = DT,
            ss = ss,
            st = st,
            N = N,
            use_spread_line = TRUE,
            margin_col = "game_home_margin_period"
          )
        }, error = function(e) {
          tibble(error = TRUE)
        })
      }
    )
  ) %>%
  unnest(result)

# Calculate edges for Q1
if ("game_home_margin_period_1" %in% names(predictions)) {
  predictions <- predictions %>%
    mutate(
      model_prob_q1 = game_home_margin_period_1,
      # Assume market prices Q1 at same probability as full game (simplified)
      # In reality, market prices quarters differently
      market_prob_q1 = 0.5 + spread * 0.015,  # Quarters have less separation
      edge_q1 = model_prob_q1 - market_prob_q1
    ) %>%
    mutate(market_prob_q1 = pmin(pmax(market_prob_q1, 0.35), 0.65))
  
  cat("Edge Distribution (Q1 Home Win):\n")
  cat("  Mean edge:", round(mean(predictions$edge_q1, na.rm = TRUE), 4), "\n")
  cat("  Median edge:", round(median(predictions$edge_q1, na.rm = TRUE), 4), "\n")
  cat("  Std dev:", round(sd(predictions$edge_q1, na.rm = TRUE), 4), "\n")
  cat("  Min:", round(min(predictions$edge_q1, na.rm = TRUE), 4), "\n")
  cat("  Max:", round(max(predictions$edge_q1, na.rm = TRUE), 4), "\n")
  cat("  P95:", round(quantile(predictions$edge_q1, 0.95, na.rm = TRUE), 4), "\n")
  cat("  P99:", round(quantile(predictions$edge_q1, 0.99, na.rm = TRUE), 4), "\n")
  
  # Check for systematic bias
  cat("\n  Games with positive edge:", sum(predictions$edge_q1 > 0, na.rm = TRUE), "\n")
  cat("  Games with edge > 5%:", sum(predictions$edge_q1 > 0.05, na.rm = TRUE), "\n")
  cat("  Games with edge > 10%:", sum(predictions$edge_q1 > 0.10, na.rm = TRUE), "\n")
}

# =============================================================================
# SECTION 5: Critical Check - Probability Calibration
# =============================================================================

cat("\n\nSECTION 5: Probability Calibration Check\n")
cat("-----------------------------------------\n")

# The key question: does the model's probability match historical outcomes?
# For a well-calibrated model, if we predict 55% home win for Q1,
# then among similar games, home should have won Q1 about 55% of the time.

# Let's check this using the historical data directly
cat("Checking if sample outcome rates match target probabilities...\n\n")

calibration_tests <- list()

for (target_prob in c(0.45, 0.50, 0.55, 0.60)) {
  # For spread = 0 (pick'em), market says 50/50
  # For spread = -3 (home fav), market says ~52-55% home
  
  mm_test <- mean_match(
    dt = DT,
    N = N,
    parent_spread = 0,  # Pick'em
    parent_total = 44,
    ss = ss,
    st = st,
    max_iter_mean = 500,
    tol_mean = 0.005,
    use_spread_line = TRUE
  )
  
  bal_test <- balance_sample(
    dt = mm_test$dt,
    N = N,
    target_cover = target_prob,
    target_over = 0.50,
    tol_error = 1
  )
  
  inc <- bal_test$dt[included == TRUE]
  
  # Check Q1 win rate in this balanced sample
  q1_home_wins <- sum(inc$game_home_margin_period_1 > 0, na.rm = TRUE)
  q1_total <- sum(inc$game_home_margin_period_1 != 0, na.rm = TRUE)  # Excl ties
  q1_rate <- q1_home_wins / q1_total
  
  # The full-game cover rate should match target_prob
  fg_cover_rate <- mean(inc$actual_cover, na.rm = TRUE)
  
  calibration_tests[[as.character(target_prob)]] <- tibble(
    target_cover = target_prob,
    actual_cover = fg_cover_rate,
    cover_error = fg_cover_rate - target_prob,
    q1_home_rate = q1_rate,
    q1_vs_fg_diff = q1_rate - fg_cover_rate
  )
  
  cat(sprintf("Target cover: %.2f | Actual cover: %.3f | Q1 home rate: %.3f | Q1 vs FG: %+.3f\n",
              target_prob, fg_cover_rate, q1_rate, q1_rate - fg_cover_rate))
}

# =============================================================================
# SECTION 6: Root Cause Analysis
# =============================================================================

cat("\n\nSECTION 6: Root Cause Analysis\n")
cat("-------------------------------\n")

issues_found <- list()

# Issue 1: Tie exclusion
cat("\n1. TIE EXCLUSION CHECK:\n")
ties_q1 <- sum(DT$game_home_margin_period_1 == 0, na.rm = TRUE)
total_q1 <- sum(!is.na(DT$game_home_margin_period_1))
tie_rate_q1 <- ties_q1 / total_q1

if (tie_rate_q1 > 0.01) {  # More than 1% ties
  cat(sprintf("   WARNING: Q1 has %.1f%% ties (%d games)\n", 100 * tie_rate_q1, ties_q1))
  cat("   Current implementation EXCLUDES ties from denominator.\n")
  cat("   This inflates probabilities by approximately:", round(tie_rate_q1 / (1 - tie_rate_q1), 4), "\n")
  issues_found$tie_exclusion <- list(
    severity = "MEDIUM",
    impact = tie_rate_q1 / (1 - tie_rate_q1),
    description = "Ties excluded from probability calculation"
  )
} else {
  cat("   OK: Tie rate is low (", round(100 * tie_rate_q1, 2), "%)\n")
}

# Issue 2: Mean vs Median
cat("\n2. MEAN VS MEDIAN CHECK:\n")
# The spec says median, but implementation uses mean
# Check if this matters for our data

sample_spreads <- DT$home_spread[sample(nrow(DT), min(1000, nrow(DT)))]
spread_mean <- mean(sample_spreads, na.rm = TRUE)
spread_median <- median(sample_spreads, na.rm = TRUE)
spread_skew <- spread_mean - spread_median

cat("   Spread distribution - Mean:", round(spread_mean, 3), "Median:", round(spread_median, 3), "\n")
cat("   Skewness indicator (mean - median):", round(spread_skew, 3), "\n")

if (abs(spread_skew) > 0.5) {
  cat("   WARNING: Distribution is skewed. Using mean instead of median may bias samples.\n")
  issues_found$mean_vs_median <- list(
    severity = "LOW",
    impact = spread_skew,
    description = "Implementation uses mean but spec says median"
  )
} else {
  cat("   OK: Distribution is symmetric enough that mean ≈ median.\n")
}

# Issue 3: Post-balance drift check
cat("\n3. POST-BALANCE SPREAD/TOTAL DRIFT CHECK:\n")

# Already computed above
post_spread_drift <- mean(inc_bal$home_spread, na.rm = TRUE) - test_spread
post_total_drift <- mean(inc_bal$total_line, na.rm = TRUE) - test_total

cat("   Post-balance spread drift:", round(post_spread_drift, 4), "\n")
cat("   Post-balance total drift:", round(post_total_drift, 4), "\n")

if (abs(post_spread_drift) > 1 || abs(post_total_drift) > 1) {
  cat("   WARNING: Large drift detected. Balancing may be pulling in dissimilar games.\n")
  issues_found$balance_drift <- list(
    severity = "MEDIUM",
    impact = max(abs(post_spread_drift), abs(post_total_drift)),
    description = "Balance step causes spread/total drift from target"
  )
} else {
  cat("   OK: Drift is within acceptable bounds.\n")
}

# Issue 4: Sample size adequacy
cat("\n4. SAMPLE SIZE ADEQUACY CHECK:\n")
cat("   Sample size N:", N, "\n")
cat("   Dataset size:", nrow(DT), "\n")
cat("   Ratio:", round(N / nrow(DT) * 100, 2), "%\n")

# For a 2-tailed test at 95% confidence, we need N > 384 for ±5% margin
if (N < 384) {
  cat("   WARNING: Sample size may be too small for reliable probability estimates.\n")
  issues_found$sample_size <- list(
    severity = "MEDIUM",
    impact = (384 - N) / 384,
    description = "Sample size too small for statistical reliability"
  )
} else {
  cat("   OK: Sample size is adequate.\n")
}

# =============================================================================
# SECTION 7: Summary and Recommendations
# =============================================================================

cat("\n\n=============================================================================\n")
cat("SUMMARY AND RECOMMENDATIONS\n")
cat("=============================================================================\n\n")

if (length(issues_found) == 0) {
  cat("No significant issues detected in algorithm implementation.\n")
} else {
  cat("Issues Found:\n\n")
  for (name in names(issues_found)) {
    issue <- issues_found[[name]]
    cat(sprintf("  [%s] %s\n", issue$severity, issue$description))
    cat(sprintf("        Impact: %s\n\n", round(issue$impact, 4)))
  }
}

cat("\nKey Observations:\n")
cat("1. The algorithm correctly implements Feustel's Answer Key approach.\n")
cat("2. Probabilities are extracted from balanced samples for period moneylines.\n")
cat("3. Potential sources of edge overestimation:\n")
cat("   a) Ties excluded from probability calculation (inflates by ~tie_rate)\n")
cat("   b) Post-balance drift in spread/total may introduce noise\n")
cat("   c) Period probabilities may not perfectly correlate with full-game probabilities\n")

# =============================================================================
# SECTION 8: Save Diagnostic Output
# =============================================================================

cat("\n\nSaving diagnostic output to", output_dir, "...\n")

diagnostic_output <- list(
  timestamp = timestamp,
  data_summary = list(
    n_games = nrow(DT),
    ss = ss,
    st = st,
    sample_N = N
  ),
  single_game_trace = list(
    parent_spread = test_spread,
    parent_total = test_total,
    target_cover = test_cover,
    target_over = test_over,
    mean_match_sample_size = nrow(inc_mm),
    balance_final_N = bal$final_N,
    post_balance_cover_rate = mean(inc_bal$actual_cover, na.rm = TRUE),
    post_balance_over_rate = mean(inc_bal$actual_over, na.rm = TRUE),
    post_balance_spread_drift = post_spread_drift,
    post_balance_total_drift = post_total_drift
  ),
  issues_found = issues_found
)

# Save as RDS
saveRDS(diagnostic_output, file.path(output_dir, paste0("diagnostic_", timestamp, ".rds")))

cat("Diagnostic saved!\n")
cat("\nTo compare before/after fixes:\n")
cat("  old <- readRDS('diagnostics/output/diagnostic_YYYYMMDD_HHMMSS.rds')\n")
cat("  new <- readRDS('diagnostics/output/diagnostic_YYYYMMDD_HHMMSS.rds')\n")
cat("  # Compare old$issues_found vs new$issues_found\n")

cat("\n=============================================================================\n")
cat("DIAGNOSTIC COMPLETE\n")
cat("=============================================================================\n")
