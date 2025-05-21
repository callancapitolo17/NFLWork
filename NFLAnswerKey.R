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

# 4. Efficient mean-matching sample builder (data.table + partial sort fix)
# -------------------------------------------------------------------------
mean_matching_sample <- function(data, parent_spread, parent_total,
                                 N = max(ceiling(0.025 * nrow(data)), 300),
                                 tol = 0.01, max_iter = 20) {
  dt <- as.data.table(data)
  scale_spread <- dt[, max(spread_line, na.rm = TRUE) - min(spread_line, na.rm = TRUE)]
  scale_total  <- dt[, max(total_line,  na.rm = TRUE) - min(total_line,  na.rm = TRUE)]
  iter <- 0
  
  repeat {
    iter <- iter + 1
    dt <- distance_index(dt, parent_spread, parent_total, scale_spread, scale_total)
    idx_vec <- as.numeric(dt$index)
    thr     <- sort(idx_vec, partial = N)[N]
    samp_dt <- dt[idx_vec <= thr]
    if (nrow(samp_dt) > N) samp_dt <- samp_dt[order(index)][1:N]
    err_spread <- mean(samp_dt$spread_line) - parent_spread
    err_total  <- mean(samp_dt$total_line)   - parent_total
    if ((abs(err_spread) < tol && abs(err_total) < tol) || iter >= max_iter) break
    parent_spread <- parent_spread + err_spread
    parent_total  <- parent_total  + err_total
  }
  samp_dt$in_sample <- TRUE
  pool_dt <- dt[!dt$game_id %in% samp_dt$game_id]
  pool_dt$in_sample <- FALSE
  list(sample        = samp_dt,
       pool          = pool_dt,
       parent_spread = parent_spread,
       parent_total  = parent_total,
       iterations    = iter)
}

# 5. Directional median-based refinement (data.table)
# ----------------------------------------------------------------
refine_sample <- function(initial, parent_spread, parent_total,
                          tol_med = 0.005, max_iter = 1000) {
  samp_dt <- copy(initial$sample)
  pool_dt <- copy(initial$pool)
  iter <- 0
  
  repeat {
    iter <- iter + 1
    # compute current median and error for spread
    med_s <- samp_dt[, median(spread_line)]
    err_s <- med_s - parent_spread
    # convergence check on spread median
    if (abs(err_s) < tol_med) break
    
    # 1) Remove the game contributing furthest in the direction of the error
    if (err_s > 0) {
      # median too high: remove largest > parent
      candidates <- samp_dt[spread_line > parent_spread]
    } else {
      # median too low: remove smallest < parent
      candidates <- samp_dt[spread_line < parent_spread]
    }
    if (nrow(candidates) == 0) break
    # pick worst by index
    worst_row <- candidates[which.max(candidates$index)]
    samp_dt <- samp_dt[game_id != worst_row$game_id]
    # re-add removed to pool for potential later re-entry
    pool_dt <- rbindlist(list(pool_dt, worst_row), use.names = TRUE)
    
    # 2) Add a game that pulls median toward parent
    if (err_s > 0) {
      # want a lower spread: pick pool game below parent
      pool_cand <- pool_dt[spread_line < parent_spread]
    } else {
      pool_cand <- pool_dt[spread_line > parent_spread]
    }
    if (nrow(pool_cand) == 0) next
    # pick candidate closest to parent on spread
    cand_idx <- pool_cand[, which.min(abs(spread_line - parent_spread))]
    add_row <- pool_cand[cand_idx]
    samp_dt <- rbindlist(list(samp_dt, add_row), use.names = TRUE)
    pool_dt <- pool_dt[game_id != add_row$game_id]
    
    if (iter >= max_iter) break
  }
  
  samp_dt
}

# 6. Example end-to-end. Example end-to-end. Example end-to-end Example end-to-end Example end-to-end
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
