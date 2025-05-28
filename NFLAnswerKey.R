# Answer Key Builder: NFL Data with nflfastR (Indexed & Efficient Median Refinement)
# Language: R

# 1. Install and load required packages
# -------------------------------------
# install.packages(c("nflfastR", "tidyverse", "data.table"))
library(nflfastR)
library(tidyverse)
library(data.table)

# 2. Fetch historical odds data
# -----------------------------
seasons <- 1999:2024
odds_raw <- map_df(seasons, ~ nflfastR::fast_scraper_schedules(.x) %>%
                     select(game_id, season, week, spread_line, total_line)) %>%
  filter(!is.na(spread_line) & !is.na(total_line))

# 3. Distance index (Pythagorean) for spread and total
# -----------------------------------------------------
distance_index <- function(data, parent_spread, parent_total, scale_spread, scale_total) {
  factor <- scale_spread / scale_total
  data %>% mutate(
    dx    = (spread_line - parent_spread) / scale_spread,
    dy    = (total_line  - parent_total)  / scale_total,
    index = dx^2 + dy^2
  )
}

# 4. Efficient mean-matching sample builder (data.table + dynamic N)
# ----------------------------------------------------------------------
mean_matching_sample <- function(data, parent_spread, parent_total,
                                 base_pct = 0.025,  # fraction of data for base sample
                                 min_N = 50,        # minimum allowed sample size
                                 tol = 0.01, max_iter = 20) {
  dt <- as.data.table(data)
  # Compute base_N solely from base_pct
  base_N <- ceiling(base_pct * nrow(dt))
  
  # use the 5th to 95th percentile range for typical variations
  qs <- dt[, quantile(spread_line, probs = c(0.05, 0.95), na.rm = TRUE)]
  scale_spread <- qs[2] - qs[1]
  qt <- dt[, quantile(total_line, probs = c(0.05, 0.95), na.rm = TRUE)]
  scale_total  <- qt[2] - qt[1]
  
  # ---------------
  # Mean-matching loop to adjust parent lines
  # ---------------
  iter <- 0
  samp_dt <- NULL
  repeat {
    iter <- iter + 1
    # recalc distances
    dt <- distance_index(dt, parent_spread, parent_total, scale_spread, scale_total)
    idx_vec <- as.numeric(dt$index)
    
    # dynamic sample size based on error magnitude
    if (is.null(samp_dt) || iter == 1) {
      N_current <- base_N
    } else {
      err_spread <- mean(samp_dt$spread_line) - parent_spread
      err_total  <- mean(samp_dt$total_line)   - parent_total
      error_ratio <- max(abs(err_spread), abs(err_total)) / tol
      N_current <- min(nrow(dt), max(min_N, ceiling(base_N * (1 + error_ratio))))
    }
    
    # select top N_current by distance
    thr <- sort(idx_vec, partial = N_current)[N_current]
    samp_dt <- dt[idx_vec <= thr]
    if (nrow(samp_dt) > N_current) samp_dt <- samp_dt[order(index)][1:N_current]
    
    # compute mean errors
    err_spread <- mean(samp_dt$spread_line) - parent_spread
    err_total  <- mean(samp_dt$total_line)   - parent_total
    # check convergence
    if ((abs(err_spread) < tol && abs(err_total) < tol) || iter >= max_iter) break
    
    # adjust parent lines
    parent_spread <- parent_spread + err_spread
    parent_total  <- parent_total  + err_total
  }
  
  # ---------------
  # Maximize sample size under constraints
  # ---------------
  sorted_dt <- copy(dt)[order(index)]
  sorted_dt[, cum_spread := cumsum(spread_line) / seq_len(.N)]
  sorted_dt[, cum_total  := cumsum(total_line)  / seq_len(.N)]
  
  valid_idx <- which(abs(sorted_dt$cum_spread - parent_spread) <= tol &
                       abs(sorted_dt$cum_total  - parent_total ) <= tol)
  if (length(valid_idx) > 0) {
    N_max <- max(valid_idx)
    samp_dt <- sorted_dt[1:N_max]
  }
  
  samp_dt[, in_sample := TRUE]
  pool_dt <- sorted_dt[!seq_len(.N) %in% seq_len(ifelse(exists("N_max"), N_max, 0))]
  pool_dt[, in_sample := FALSE]
  
  list(sample        = samp_dt,
       pool          = pool_dt,
       parent_spread = parent_spread,
       parent_total  = parent_total,
       iterations    = iter)
}

# 5. Directional 2D median refinement (data.table). Directional 2D median refinement (data.table)
# -----------------------------------------------------
refine_sample <- function(initial, parent_spread, parent_total,
                          tol_med = 0.005, max_iter = 1000) {
  samp_dt <- copy(initial$sample)
  pool_dt <- copy(initial$pool)
  iter <- 0
  
  repeat {
    iter <- iter + 1
    # compute current medians and errors
    med_s <- samp_dt[, median(spread_line)]
    med_t <- samp_dt[, median(total_line)]
    err_s <- med_s - parent_spread
    err_t <- med_t - parent_total
    
    # stop if both within tolerance
    if (abs(err_s) < tol_med && abs(err_t) < tol_med) break
    
    # choose dimension with larger relative error
    if (abs(err_s) >= abs(err_t)) {
      if (err_s > 0) {
        to_remove <- samp_dt[spread_line > parent_spread]
        to_add    <- pool_dt[spread_line < parent_spread]
      } else {
        to_remove <- samp_dt[spread_line < parent_spread]
        to_add    <- pool_dt[spread_line > parent_spread]
      }
    } else {
      if (err_t > 0) {
        to_remove <- samp_dt[total_line > parent_total]
        to_add    <- pool_dt[total_line < parent_total]
      } else {
        to_remove <- samp_dt[total_line < parent_total]
        to_add    <- pool_dt[total_line > parent_total]
      }
    }
    
    # if no candidates to remove or add, break
    if (nrow(to_remove) == 0 || nrow(to_add) == 0) break
    
    # remove the worst by 2D index
    worst_idx <- to_remove[, which.max(index)]
    worst     <- to_remove[worst_idx]
    samp_dt   <- samp_dt[game_id != worst$game_id]
    pool_dt   <- rbindlist(list(pool_dt, worst), use.names = TRUE)
    
    # add the best by 2D index
    add_candidates <- to_add
    best_idx <- add_candidates[, which.min(index)]
    best     <- add_candidates[best_idx]
    samp_dt  <- rbindlist(list(samp_dt, best), use.names = TRUE)
    pool_dt  <- pool_dt[game_id != best$game_id]
    
    if (iter >= max_iter) break
  }
  
  samp_dt
}

# 6. Example end-to-end
# ---------------------
parent_spread <- 3   # initial target spread
parent_total  <- 44  # initial target total
mm_init      <- mean_matching_sample(odds_raw, parent_spread, parent_total)
final_sample <- refine_sample(mm_init, mm_init$parent_spread, mm_init$parent_total)

# 7. Inspect final sample
# -----------------------
cat("Mean-match iterations:", mm_init$iterations, "\n")
cat("Adjusted parent line: Spread =", round(mm_init$parent_spread,2),
    ", Total =", round(mm_init$parent_total,2), "\n")
cat("Final sample size:", nrow(final_sample), "games\n")
cat("Median spread:", median(final_sample$spread_line),
    "| Median total:", median(final_sample$total_line), "\n")

final_sample[order(index), .(game_id, season, week, spread_line, total_line, index)]

# 8. Next steps: integrate probability model for Answer Key EV computation
# -------------------------------------------------------------------------
