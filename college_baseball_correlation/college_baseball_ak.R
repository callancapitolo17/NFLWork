# College Baseball Answer Key Backtest
# Leave-one-out validation: verify correlation factors match expectations
# Medium totals (12-13.5) should show CFs of 1.11-1.17 for positively-correlated combos

setwd("~/NFLWork/Answer Keys")
suppressPackageStartupMessages({
  library(data.table)
  library(duckdb)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(lubridate)
  library(DBI)
})
source("Tools.R")

cat("=== COLLEGE BASEBALL ANSWER KEY BACKTEST ===\n")

# =============================================================================
# LOAD HISTORICAL DATA
# =============================================================================

db_path <- "~/NFLWork/college_baseball_correlation/college_baseball.duckdb"
if (!file.exists(db_path)) {
  stop("Database not found. Run: python validate.py --pull-scores 2024 2025 --pull-odds 2024 2025")
}

con <- dbConnect(duckdb(), dbdir = db_path, read_only = TRUE)

# Join games to consensus odds
matched <- dbGetQuery(con, "
  WITH consensus AS (
    SELECT
      home_team, away_team,
      CAST(commence_time AS DATE) as game_date,
      MEDIAN(home_ml) as home_ml,
      MEDIAN(away_ml) as away_ml,
      MEDIAN(total_line) as total_line
    FROM odds
    WHERE home_ml IS NOT NULL AND away_ml IS NOT NULL AND total_line IS NOT NULL
    GROUP BY home_team, away_team, CAST(commence_time AS DATE)
  )
  SELECT
    g.game_id,
    g.home_team, g.away_team,
    g.home_score, g.away_score,
    g.total_runs, g.home_margin,
    CASE WHEN g.home_won THEN 1 ELSE 0 END as home_won,
    c.home_ml, c.away_ml, c.total_line
  FROM games g
  INNER JOIN consensus c
    ON g.home_team = c.home_team AND g.away_team = c.away_team
    AND g.date = c.game_date
  WHERE NOT g.went_extras AND NOT g.mercy_rule
    AND c.home_ml IS NOT NULL AND c.away_ml IS NOT NULL
")

dbDisconnect(con)

cat(sprintf("Loaded %d matched regulation games\n", nrow(matched)))

if (nrow(matched) < 100) {
  stop("Insufficient matched games for backtest. Need 100+, got ", nrow(matched))
}

# =============================================================================
# BUILD DT (data.table for answer key)
# =============================================================================

# Baseball has no spread — use logit(P(home)) as implied spread
# P(home) = devigged ML probability
devig_ml <- function(home_ml, away_ml) {
  # Convert American to implied probability
  imp_home <- ifelse(home_ml > 0, 100 / (home_ml + 100), abs(home_ml) / (abs(home_ml) + 100))
  imp_away <- ifelse(away_ml > 0, 100 / (away_ml + 100), abs(away_ml) / (abs(away_ml) + 100))
  total_imp <- imp_home + imp_away
  # Devig
  list(p_home = imp_home / total_imp, p_away = imp_away / total_imp)
}

probs <- devig_ml(matched$home_ml, matched$away_ml)

DT <- data.table(
  id = matched$game_id,
  home_team = matched$home_team,
  away_team = matched$away_team,
  home_spread = qlogis(probs$p_home),  # logit of P(home) as pseudo-spread
  total_line = matched$total_line,
  home_margin = matched$home_margin,
  total_final_score = matched$total_runs,
  actual_cover = matched$home_won,
  actual_over = as.integer(matched$total_runs > matched$total_line),
  home_ml_odds = matched$home_ml,
  away_ml_odds = matched$away_ml
)

cat(sprintf("DT: %d rows, spread range [%.2f, %.2f], total range [%.1f, %.1f]\n",
            nrow(DT), min(DT$home_spread), max(DT$home_spread),
            min(DT$total_line), max(DT$total_line)))

# Dispersion
disp <- compute_dispersion(DT, moneyline = TRUE,
                           odds_cols = c("home_ml_odds", "away_ml_odds"))
ss <- disp$ss
st <- disp$st
cat(sprintf("Dispersion: ss=%.2f, st=%.2f\n", ss, st))

N <- max(50, round(nrow(DT) * 0.05))
cat(sprintf("Sample size N=%d (5%% of pool)\n", N))

# =============================================================================
# LEAVE-ONE-OUT BACKTEST
# =============================================================================

cat("\nRunning leave-one-out backtest...\n")

results <- list()
n_total <- nrow(DT)
t0 <- Sys.time()

for (i in seq_len(n_total)) {
  # Leave this game out
  game <- DT[i, ]
  pool <- DT[-i, ]

  # Generate sample for this game's line
  sample_result <- tryCatch({
    run_answer_key_sample(
      id = game$id,
      parent_spread = game$home_spread,
      parent_total = game$total_line,
      target_cover = 0.5,  # Neutral (no prior on who wins)
      target_over = 0.5,   # Neutral
      DT = pool,
      ss = ss, st = st, N = N,
      use_spread_line = FALSE,
      max_iter_mean = 300,
      tol_mean = 0.01,
      tol_error = 2,
      shrink_factor = 0.9,
      min_N = 30
    )
  }, error = function(e) NULL)

  if (is.null(sample_result) || is.null(sample_result$sample) || nrow(sample_result$sample) < 30) {
    next
  }

  samp <- sample_result$sample

  # Compute parlay fair odds for 4 combos
  combos <- list(
    home_over  = list(
      list(market = "moneyline", period = "Full", side = "home", line = NA),
      list(market = "total", period = "Full", side = "over", line = game$total_line)
    ),
    home_under = list(
      list(market = "moneyline", period = "Full", side = "home", line = NA),
      list(market = "total", period = "Full", side = "under", line = game$total_line)
    ),
    away_over  = list(
      list(market = "moneyline", period = "Full", side = "away", line = NA),
      list(market = "total", period = "Full", side = "over", line = game$total_line)
    ),
    away_under = list(
      list(market = "moneyline", period = "Full", side = "away", line = NA),
      list(market = "total", period = "Full", side = "under", line = game$total_line)
    )
  )

  for (combo_name in names(combos)) {
    parlay <- tryCatch({
      compute_parlay_fair_odds(samp, combos[[combo_name]])
    }, error = function(e) NULL)

    if (is.null(parlay)) next

    results[[length(results) + 1]] <- data.frame(
      game_id = game$id,
      total_line = game$total_line,
      combo = combo_name,
      joint_prob = parlay$joint_prob,
      independent_prob = parlay$independent_prob,
      correlation_factor = parlay$correlation_factor,
      n_samples = parlay$n_samples_resolved,
      # Actual outcome
      home_won = game$actual_cover,
      went_over = game$actual_over,
      stringsAsFactors = FALSE
    )
  }

  if (i %% 100 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    rate <- i / elapsed
    eta <- (n_total - i) / rate
    cat(sprintf("  %d/%d (%.0f/s, ETA %.0fs)\n", i, n_total, rate, eta))
  }
}

elapsed_total <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("Backtest complete: %d results in %.1f seconds\n",
            length(results), elapsed_total))

# =============================================================================
# ANALYZE RESULTS
# =============================================================================

if (length(results) == 0) {
  stop("No results generated. Check data quality.")
}

res_df <- bind_rows(results)

# Total line buckets
res_df <- res_df %>%
  mutate(total_bucket = case_when(
    total_line < 12 ~ "low (<12)",
    total_line <= 13.5 ~ "medium (12-13.5)",
    TRUE ~ "high (>13.5)"
  ))

cat("\n=== CORRELATION FACTORS BY TOTAL BUCKET ===\n")
summary_df <- res_df %>%
  group_by(total_bucket, combo) %>%
  summarise(
    n = n(),
    mean_cf = mean(correlation_factor, na.rm = TRUE),
    median_cf = median(correlation_factor, na.rm = TRUE),
    sd_cf = sd(correlation_factor, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(total_bucket, combo)

print(summary_df, n = 50)

# Summary by combo across all buckets
cat("\n=== OVERALL CORRELATION FACTORS ===\n")
overall <- res_df %>%
  group_by(combo) %>%
  summarise(
    n = n(),
    mean_cf = mean(correlation_factor, na.rm = TRUE),
    median_cf = median(correlation_factor, na.rm = TRUE),
    .groups = "drop"
  )
print(overall)

# Calibration: does AK fair prob predict actual outcomes?
cat("\n=== CALIBRATION CHECK ===\n")
cal <- res_df %>%
  mutate(
    actual_hit = case_when(
      combo == "home_over"  ~ (home_won == 1 & went_over == 1),
      combo == "home_under" ~ (home_won == 1 & went_over == 0),
      combo == "away_over"  ~ (home_won == 0 & went_over == 1),
      combo == "away_under" ~ (home_won == 0 & went_over == 0),
    ),
    prob_bucket = round(joint_prob * 20) / 20  # 5% buckets
  ) %>%
  group_by(prob_bucket) %>%
  summarise(
    n = n(),
    predicted_rate = mean(joint_prob),
    actual_rate = mean(actual_hit),
    .groups = "drop"
  ) %>%
  filter(n >= 20)

print(cal, n = 30)

cat("\n=== BACKTEST COMPLETE ===\n")
cat("Key finding: Look for CFs > 1.10 in medium-total bucket for Away+Over and Home+Under.\n")
cat("These are the positively-correlated combos that books underprice.\n")
