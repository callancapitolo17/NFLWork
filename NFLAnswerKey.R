# Answer Key Builder: NFL Data with nflfastR (Exact Stepwise Implementation)
# Language: R

# 1. Install and load required packages
# -------------------------------------
# install.packages(c("nflfastR", "data.table", "dplyr"))
library(nflfastR)
library(data.table)
library(dplyr)

# Helper: de-vig two American-style odds, return data.frame
devig_american <- function(odd1, odd2) {
  p1_raw <- ifelse(odd1 > 0,
                   100 / (odd1 + 100),
                   -odd1 / (-odd1 + 100))
  p2_raw <- ifelse(odd2 > 0,
                   100 / (odd2 + 100),
                   -odd2 / (-odd2 + 100))
  total_raw <- p1_raw + p2_raw
  data.frame(p1 = p1_raw/total_raw, p2 = p2_raw/total_raw)
}

# 2. Fetch historical odds and compute de-vig probabilities
DT <- nflfastR::fast_scraper_schedules(1999:2024) %>%
  select(game_id, season, week,
         home_score, away_score,
         spread_line, total_line,
         home_spread_odds, away_spread_odds,
         over_odds, under_odds) %>%
  filter(!is.na(spread_line), !is.na(total_line),
         !is.na(home_spread_odds), !is.na(away_spread_odds),
         !is.na(over_odds), !is.na(under_odds)) %>%
  as.data.table()

# Compute de-vigged probabilities for spread and over
spr <- devig_american(DT$home_spread_odds, DT$away_spread_odds)
DT[, `:=`(home_spread_prob = spr$p1,
          away_spread_prob = spr$p2)]
ou <- devig_american(DT$over_odds, DT$under_odds)
DT[, `:=`(over_prob = ou$p1,
          under_prob = ou$p2)]

# 4. Distance index helper (Pythagorean)
distance_index <- function(dt, ps, pt, ss, st) {
  dt[, index := ((spread_line - ps)/ss)^2 + ((total_line - pt)/st)^2]
}

# 5. Answer Key sampler: mean-match then stepwise refine
answer_key_sample <- function(dt,
                              parent_spread, parent_total,
                              home_sp_odds, away_sp_odds,
                              over_odds,    under_odds,
                              N = 500,
                              tol_mean = 0.01,
                              tol_refine = 1,
                              max_iter_mean = 20,
                              max_iter_refine = 10000) {
  dt <- copy(dt)
  # Compute de-vig market targets
  target_cover <- devig_american(home_sp_odds, away_sp_odds)$p1
  target_over  <- devig_american(over_odds,    under_odds)$p1
  # Compute actual outcomes based on the PARENT lines, not game lines
  dt[, actual_cover := as.integer((home_score - away_score + parent_spread) > 0)]
  dt[, actual_over  := as.integer((home_score + away_score) > parent_total)]
  # Pre-calc scales for indexing
  qs <- dt[, quantile(spread_line, probs = c(0.05, 0.95))]
  ss <- qs[2] - qs[1]
  qt <- dt[, quantile(total_line,  probs = c(0.05, 0.95))]
  st <- qt[2] - qt[1]
  
  # --- Phase 1: Mean-matching loop ---
  for (i in seq_len(max_iter_mean)) {
    distance_index(dt, parent_spread, parent_total, ss, st)
    setorder(dt, index)
    dt[, included := FALSE]
    dt[1:N, included := TRUE]
    mean_s <- dt[included == TRUE, mean(spread_line)]
    mean_t <- dt[included == TRUE, mean(total_line)]
    err_s  <- mean_s - parent_spread
    err_t  <- mean_t - parent_total
    if (abs(err_s) < tol_mean && abs(err_t) < tol_mean) break
    parent_spread <- parent_spread + err_s
    parent_total  <- parent_total  + err_t
    # update outcomes after shifting parent lines
    dt[, actual_cover := as.integer((home_score - away_score + parent_spread) > 0)]
    dt[, actual_over  := as.integer((home_score + away_score) > parent_total)]
  }
  
  # --- Phase 2 & 3: Interleaved stepwise refine (cover & over) ---
  cover_error <- dt[included == TRUE, sum(actual_cover)] - round(target_cover * N)
  over_error  <- dt[included == TRUE, sum(actual_over)]  - round(target_over  * N)
  iter <- 0
  while ((cover_error != 0 || over_error != 0) && iter < max_iter_refine) {
    if (abs(cover_error) >= abs(over_error)) {
      # refine cover
      bit_rem <- as.integer(cover_error > 0)
      idx_rem <- dt[included == TRUE & actual_cover == bit_rem, .I[.N]]
      if (length(idx_rem) == 0) break
      dt[idx_rem, included := FALSE]
      cover_error <- cover_error - dt[idx_rem, actual_cover]
      bit_add <- as.integer(cover_error < 0)
      idx_add <- dt[included == FALSE & actual_cover == bit_add, .I[1]]
      if (length(idx_add) == 0) break
      dt[idx_add, included := TRUE]
      cover_error <- cover_error + dt[idx_add, actual_cover]
    } else {
      # refine over
      bit_rem <- as.integer(over_error > 0)
      idx_rem <- dt[included == TRUE & actual_over == bit_rem, .I[.N]]
      if (length(idx_rem) == 0) break
      dt[idx_rem, included := FALSE]
      over_error <- over_error - dt[idx_rem, actual_over]
      bit_add <- as.integer(over_error < 0)
      idx_add <- dt[included == FALSE & actual_over == bit_add, .I[1]]
      if (length(idx_add) == 0) break
      dt[idx_add, included := TRUE]
      over_error <- over_error + dt[idx_add, actual_over]
    }
    iter <- iter + 1
  }
  
  # Return final sample + targets
  list(
    sample = dt[included == TRUE],
    target_cover = target_cover,
    target_over  = target_over
  )
}

# 6. Example end-to-end
parent_spread <- 3.5
parent_total  <- 46.0
home_sp_odds <- +110; away_sp_odds <- -130
over_odds    <- -110; under_odds    <- +100

result <- answer_key_sample(
  DT, parent_spread, parent_total,
  home_sp_odds, away_sp_odds,
  over_odds,    under_odds,
  N = 500
)
final_sample <- result$sample
target_cover <- result$target_cover
target_over  <- result$target_over

cat("Sample size: ", nrow(final_sample), "games\n")
cat("Cover rate: ", round(mean(final_sample$actual_cover),3), "| Target:", round(target_cover,3), "\n")
cat("Over rate:  ", round(mean(final_sample$actual_over),3),   "| Target:", round(target_over,3),   "\n")
final_sample[1:10, .(game_id, spread_line, total_line, actual_cover, actual_over, index)]
