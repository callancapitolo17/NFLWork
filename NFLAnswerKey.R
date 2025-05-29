# Answer Key Builder: NFL Data with nflfastR (Indexed & Efficient Median Refinement)
# Language: R

# 1. Install and load required packages
# -------------------------------------
# install.packages(c("nflfastR", "tidyverse", "data.table"))
library(nflfastR)
library(tidyverse)
library(data.table)

# 2. Fetch historical odds + moneyline data
# -------------------------------------------
seasons <- 1999:2024
# grab schedule (lines & scores) and moneyline odds
dt <- map_df(seasons, ~ nflfastR::fast_scraper_schedules(.x)) %>%
  select(game_id, season, week,
         home_score, away_score,
         spread_line, total_line,
         home_moneyline, away_moneyline) %>%
  filter(!is.na(spread_line) & !is.na(total_line) &
           !is.na(home_moneyline) & !is.na(away_moneyline)) %>%
  as.data.table()

# compute implied de-vig probabilities for underdog cover
# --------------------------------------------------------
dt[, `:=`(
  # implied win probabilities from moneyline (raw)
  p_home_raw = fifelse(home_moneyline > 0,
                       100 / (home_moneyline + 100),
                       -home_moneyline / (-home_moneyline + 100)),
  p_away_raw = fifelse(away_moneyline > 0,
                       100 / (away_moneyline + 100),
                       -away_moneyline / (-away_moneyline + 100))
)]
# normalize to remove vig
dt[, `:=`(
  p_home = p_home_raw / (p_home_raw + p_away_raw),
  p_away = p_away_raw / (p_home_raw + p_away_raw)
)]
# flag which side is underdog for spread (moneyline-based)  
# dt[, underdog_ml := ifelse(spread_line > 0, p_home, p_away)]

# compute implied de-vig probabilities for spread odds (if available)
if (all(c("home_spread_odds","away_spread_odds") %in% names(dt))) {
  dt[, `:=`(
    p_home_spread_raw = fifelse(home_spread_odds > 0,
                                100 / (home_spread_odds + 100),
                                -home_spread_odds / (-home_spread_odds + 100)),
    p_away_spread_raw = fifelse(away_spread_odds > 0,
                                100 / (away_spread_odds + 100),
                                -away_spread_odds / (-away_spread_odds + 100))
  )]
  # normalize to remove vig for spread
  dt[, `:=`(
    p_home_spread = p_home_spread_raw / (p_home_spread_raw + p_away_spread_raw),
    p_away_spread = p_away_spread_raw / (p_home_spread_raw + p_away_spread_raw)
  )]
  # set spread target probability based on underdog side for cover
  dt[, spread_target := fifelse(spread_line > 0, p_home_spread, p_away_spread)]
} else {
  # fallback to moneyline implied probabilities
  dt[, spread_target := underdog_ml]
}

# compute implied de-vig over probability if odds available
# ---------------------------------------------------------
# assume columns over_moneyline and under_moneyline exist
if (all(c("over_moneyline","under_moneyline") %in% names(dt))) {
  dt[, `:=`(
    p_over_raw = fifelse(over_moneyline > 0,
                         100 / (over_moneyline + 100),
                         -over_moneyline / (-under_moneyline + 100)),
    p_under_raw = fifelse(under_moneyline > 0,
                          100 / (under_moneyline + 100),
                          -under_moneyline / (-under_moneyline + 100))
  )]
  dt[, `:=`(
    p_over = p_over_raw / (p_over_raw + p_under_raw),
    p_under = p_under_raw / (p_over_raw + p_under_raw)
  )]
}

# 3. Compute actual cover + over outcomes
# ----------------------------------------
dt[, actual_cover := as.integer((home_score - away_score + spread_line) > 0)]
dt[, actual_over  := as.integer((home_score + away_score) > total_line)]

# 4. Distance index helper (compute only index)
# ---------------------------------------------
distance_index <- function(DT, ps, pt, ss, st) {
  DT[, index := ((spread_line - ps)/ss)^2 + ((total_line - pt)/st)^2]
}

# 5. Answer Key sampler with de-vig targets
# ------------------------------------------
answer_key_sample <- function(DT,
                              parent_spread, parent_total,
                              N = ceiling(0.025 * nrow(DT)),
                              tol_rate = 0.01,
                              max_iter = 1e5) {
  
  dt <- copy(DT)
  
  # a) compute typical variation (5th–95th percentile)
  qs <- dt[, quantile(spread_line, probs = c(0.05, 0.95), na.rm = TRUE)]
  ss <- qs[2] - qs[1]
  qt <- dt[, quantile(total_line,  probs = c(0.05, 0.95), na.rm = TRUE)]
  st <- qt[2] - qt[1]
  
  # b) sort by similarity and initial sample
  distance_index(dt, parent_spread, parent_total, ss, st)
  setorder(dt, index)
  dt[, included := FALSE]
  dt[1:N, included := TRUE]
  
  # c) compute de-vigged target rates
  target_cover <- mean(dt[spread_line == parent_spread]$underdog_ml)
  if ("p_over" %in% names(dt)) {
    target_over <- mean(dt[total_line == parent_total]$p_over)
  } else {
    target_over <- 0.5
  }
  
  # d) one-step adjustment function (flag_col, target_rate)
  one_pass <- function(flag_col, target_rate) {
    current_rate <- mean(dt[included == TRUE][[flag_col]])
    err <- current_rate - target_rate
    if (abs(err) <= tol_rate) return(FALSE)
    
    if (err > 0) {
      rem_bit <- 1L; add_bit <- 0L
    } else {
      rem_bit <- 0L; add_bit <- 1L
    }
    # remove worst-fit included of rem_bit
    rem_idx <- dt[included == TRUE & get(flag_col) == rem_bit, .I[.N]]
    if (length(rem_idx) == 0) return(FALSE)
    dt[rem_idx, included := FALSE]
    
    # add best-fit excluded of add_bit
    add_idx <- dt[included == FALSE & get(flag_col) == add_bit, .I[1]]
    if (length(add_idx) == 0) return(FALSE)
    dt[add_idx, included := TRUE]
    
    TRUE
  }
  
  # e) iterative adjustment: covers then overs
  iter <- 0
  while (iter < max_iter && one_pass("actual_cover", target_cover)) {
    iter <- iter + 1
  }
  while (iter < max_iter && one_pass("actual_over",  target_over)) {
    iter <- iter + 1
  }
  
  dt_final <- dt[included == TRUE]
  dt_final
}

# 6. Example end-to-end
# ---------------------
# Set your target lines
parent_spread <- 3.5   # target spread for underdog cover
parent_total  <- 46.0  # target total for over/under

# Run the Answer Key sampler
final_sample <- answer_key_sample(dt, parent_spread, parent_total)

# Inspect results
cat("Final sample size:", nrow(final_sample), "games
")

# Compare sample rates to de-vigged targets
# De-vigged targets:
target_cover <- mean(dt[spread_line == parent_spread]$underdog_ml)
if ("p_over" %in% names(dt)) {
  target_over <- mean(dt[total_line == parent_total]$p_over)
} else {
  target_over <- 0.5
}
cat("Target cover (de-vig):", round(target_cover,3), " | Sample cover:", round(mean(final_sample$actual_cover),3), "
")
cat("Target over  (de-vig):", round(target_over,3),  " | Sample over: ", round(mean(final_sample$actual_over),3), "
")

# View top 10 games in final sample by similarity
final_sample[order(index)][1:10, .(
  game_id, season, week,
  spread_line, total_line,
  underdog_ml, actual_cover, actual_over, index
)]
