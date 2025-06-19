# Answer Key Builder: NFL Data with nflfastR (Elihu-Style Greedy Error Minimization, Robust Swap Logic)
# Language: R

library(nflreadr)
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
DT <- load_schedules(1999:2024) %>%
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

distance_index <- function(dt, ps, pt, ss, st) {
  dt[, index := ((spread_line - ps)/ss)^2 + ((total_line - pt)/st)^2]
}

# Elihu-Style Add/Remove Refinement: improved logic
eliu_refine_addremove <- function(dt, N, target_cover, target_over, tol_error = 1, max_iter = 2000) {
  iter <- 0
  dt <- copy(dt)
  errors <- c(
    cover_error = sum(dt[included == TRUE, actual_cover]) - round(target_cover * N),
    over_error  = sum(dt[included == TRUE, actual_over ]) - round(target_over  * N)
  )
  while (any(abs(errors) > tol_error) && iter < max_iter) {
    iter <- iter + 1
    improved <- FALSE
    inc <- dt[included == TRUE]
    exc <- dt[included == FALSE]
    # Remove logic: worst-included (highest index), one at a time
    if (abs(errors['cover_error']) >= abs(errors['over_error'])) {
      # Cover error biggest
      bit <- as.integer(errors['cover_error'] > 0) # Remove 1 if too many covers, 0 if too few
      worst_inc <- inc[actual_cover == bit]
      if (nrow(worst_inc) > 0) {
        remove_idx <- worst_inc[which.max(index), which = TRUE] # the highest index
        test_incl <- dt$included
        test_incl[remove_idx] <- FALSE
        new_cover_error <- sum(dt$actual_cover[test_incl]) - round(target_cover * N)
        if (abs(new_cover_error) < abs(errors['cover_error'])) {
          dt$included[inc$.I[remove_idx]] <- FALSE
          errors['cover_error'] <- new_cover_error
          improved <- TRUE
        }
      }
      # Add logic: best-excluded (lowest index)
      if (!improved) {
        bit_add <- as.integer(errors['cover_error'] < 0)
        best_exc <- exc[actual_cover == bit_add]
        if (nrow(best_exc) > 0) {
          add_idx <- best_exc[which.min(index), which = TRUE] # the lowest index
          test_incl <- dt$included
          test_incl[exc$.I[add_idx]] <- TRUE
          new_cover_error <- sum(dt$actual_cover[test_incl]) - round(target_cover * N)
          if (abs(new_cover_error) < abs(errors['cover_error'])) {
            dt$included[exc$.I[add_idx]] <- TRUE
            errors['cover_error'] <- new_cover_error
            improved <- TRUE
          }
        }
      }
    } else {
      # Over error biggest
      bit <- as.integer(errors['over_error'] > 0) # Remove 1 if too many overs, 0 if too few
      worst_inc <- inc[actual_over == bit]
      if (nrow(worst_inc) > 0) {
        remove_idx <- worst_inc[which.max(index), which = TRUE]
        test_incl <- dt$included
        test_incl[inc$.I[remove_idx]] <- FALSE
        new_over_error <- sum(dt$actual_over[test_incl]) - round(target_over * N)
        if (abs(new_over_error) < abs(errors['over_error'])) {
          dt$included[inc$.I[remove_idx]] <- FALSE
          errors['over_error'] <- new_over_error
          improved <- TRUE
        }
      }
      if (!improved) {
        bit_add <- as.integer(errors['over_error'] < 0)
        best_exc <- exc[actual_over == bit_add]
        if (nrow(best_exc) > 0) {
          add_idx <- best_exc[which.min(index), which = TRUE]
          test_incl <- dt$included
          test_incl[exc$.I[add_idx]] <- TRUE
          new_over_error <- sum(dt$actual_over[test_incl]) - round(target_over * N)
          if (abs(new_over_error) < abs(errors['over_error'])) {
            dt$included[exc$.I[add_idx]] <- TRUE
            errors['over_error'] <- new_over_error
            improved <- TRUE
          }
        }
      }
    }
    if (!improved) break
  }
  list(
    dt = dt,
    final_error = errors,
    iterations = iter
  )
}

answer_key_sample <- function(dt,
                              parent_spread, parent_total,
                              home_sp_odds, away_sp_odds,
                              over_odds,    under_odds,
                              N = 500,
                              tol_mean = 0.01,
                              max_iter_mean = 20,
                              max_iter_refine = 2000,
                              tol_error = 1) {
  dt <- copy(dt)
  target_cover <- devig_american(home_sp_odds, away_sp_odds)$p1
  target_over  <- devig_american(over_odds,    under_odds)$p1
  dt[, actual_cover := as.integer((home_score - away_score) > parent_spread)]
  dt[, actual_over  := as.integer((home_score + away_score) > parent_total)]
  qs <- dt[, quantile(spread_line, probs = c(0.05, 0.95))]
  ss <- qs[2] - qs[1]
  qt <- dt[, quantile(total_line,  probs = c(0.05, 0.95))]
  st <- qt[2] - qt[1]
  # 1. Mean-match
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
  }
  # 2. Greedy swap refinement
  swap_result <- eliu_refine_addremove(dt, N, target_cover, target_over, tol_error, max_iter_refine)
  dt <- swap_result$dt
  final_error <- swap_result$final_error
  iterations  <- swap_result$iterations
  list(
    sample = dt[included == TRUE],
    target_cover = target_cover,
    target_over  = target_over,
    iterations   = iterations,
    final_error  = final_error
  )
}

# Example end-to-end run
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
