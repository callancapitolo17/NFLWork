# Assumptions (already defined):
#   dt               : data.table with columns:
#                       - spread_line, total_line
#                       - actual_cover (0/1)
#                       - actual_over  (0/1)
#                       - included     (logical)
#   N                : integer size of the original sample (e.g. 500)
#   parent_spread    : numeric
#   parent_total     : numeric
#   ss, st           : extra args for distance_index()
#   max_iter_mean    : integer
#   tol_mean         : numeric
#   target_cover     : numeric (e.g. 0.532)
#   target_over      : numeric (e.g. 0.49)
#   tol_error        : integer (e.g. 1)

library(data.table)

# 1) “Mean‐matching” loop (unchanged)
for (iter in seq_len(max_iter_mean)) {
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

# 2) Initialize “raw” cover_error and over_error from the original sample
#    (always use the fixed N=500; do NOT adjust N when removing/adding)
set_cover_error <- dt[included == TRUE, sum(actual_cover)] -
  round(target_cover * N)
set_over_error  <- dt[included == TRUE, sum(actual_over )] -
  round(target_over  * N)

# 3) Copy into working vars
cover_error <- set_cover_error
over_error  <- set_over_error

# 4) Repeat removal/addition until within ±tol_error or no improvement
repeat {
  cover_adjustments <- abs(cover_error)
  over_adjustments  <- abs(over_error)
  
  # If both metrics are within tolerance, stop.
  if (cover_adjustments <= tol_error && over_adjustments <= tol_error) {
    break
  }
  
  #
  # ----- 4a) REMOVE one game from rows 1..N, preferring “improves both” -----
  #
  removal_failed <- TRUE
  
  # First try to find a game i in 1..N where removing helps BOTH cover & over.
  for (i in seq(N, 1)) {
    if (!dt[i, included]) next
    
    # Would removing this game help cover?
    improves_cover <- (cover_error > 0  && dt[i, actual_cover] == 1) ||
      (cover_error < 0  && dt[i, actual_cover] == 0)
    # Would removing this game help over?
    improves_over  <- (over_error > 0   && dt[i, actual_over] == 1) ||
      (over_error < 0   && dt[i, actual_over] == 0)
    
    if (improves_cover && improves_over) {
      # Remove the game
      dt[i, included := FALSE]
      removal_failed <- FALSE
      
      # Recompute BOTH errors from scratch (using fixed N=500)
      cover_error <- dt[included == TRUE, sum(actual_cover)] -
        round(target_cover * N)
      over_error  <- dt[included == TRUE, sum(actual_over )] -
        round(target_over  * N)
      break
    }
  }
  
  # If no single game improved both, try to find one that improves at least one:
  if (removal_failed) {
    for (i in seq(N, 1)) {
      if (!dt[i, included]) next
      
      # Removing helps cover?
      improves_cover <- (cover_error > 0  && dt[i, actual_cover] == 1) ||
        (cover_error < 0  && dt[i, actual_cover] == 0)
      # Removing hurts cover?
      hurts_cover    <- (!improves_cover) &&
        ((cover_error > 0 && dt[i, actual_cover] == 0) ||
           (cover_error < 0 && dt[i, actual_cover] == 1))
      
      # Removing helps over?
      improves_over  <- (over_error > 0   && dt[i, actual_over] == 1) ||
        (over_error < 0   && dt[i, actual_over] == 0)
      # Removing hurts over?
      hurts_over     <- (!improves_over) &&
        ((over_error > 0  && dt[i, actual_over] == 0) ||
           (over_error < 0  && dt[i, actual_over] == 1))
      
      if (improves_cover || improves_over) {
        dt[i, included := FALSE]
        removal_failed <- FALSE
        
        # Recompute BOTH errors from scratch (still using N=500)
        cover_error <- dt[included == TRUE, sum(actual_cover)] -
          round(target_cover * N)
        over_error  <- dt[included == TRUE, sum(actual_over )] -
          round(target_over  * N)
        break
      }
    }
  }
  
  #
  # ----- 4b) ADD one game from rows (N+1)..nrow(dt), preferring “improves both” -----
  #
  addition_failed <- TRUE
  
  # First try to find an out‐of‐sample game that helps BOTH cover & over.
  for (i in seq(N + 1, nrow(dt))) {
    if (dt[i, included]) next
    
    # Adding helps cover?
    improves_cover <- (cover_error > 0  && dt[i, actual_cover] == 0) ||
      (cover_error < 0  && dt[i, actual_cover] == 1)
    # Adding helps over?
    improves_over  <- (over_error > 0   && dt[i, actual_over] == 0) ||
      (over_error < 0   && dt[i, actual_over] == 1)
    
    if (improves_cover && improves_over) {
      dt[i, included := TRUE]
      addition_failed <- FALSE
      
      # Recompute BOTH errors from scratch with fixed N=500
      cover_error <- dt[included == TRUE, sum(actual_cover)] -
        round(target_cover * N)
      over_error  <- dt[included == TRUE, sum(actual_over )] -
        round(target_over  * N)
      break
    }
  }
  
  # If no “both” found, second pass: find one that helps at least one
  if (addition_failed) {
    for (i in seq(N + 1, nrow(dt))) {
      if (dt[i, included]) next
      
      # Adding helps cover?
      improves_cover <- (cover_error > 0  && dt[i, actual_cover] == 0) ||
        (cover_error < 0  && dt[i, actual_cover] == 1)
      # Adding hurts cover?
      hurts_cover    <- (!improves_cover) &&
        ((cover_error > 0 && dt[i, actual_cover] == 1) ||
           (cover_error < 0 && dt[i, actual_cover] == 0))
      
      # Adding helps over?
      improves_over  <- (over_error > 0   && dt[i, actual_over] == 0) ||
        (over_error < 0   && dt[i, actual_over] == 1)
      # Adding hurts over?
      hurts_over     <- (!improves_over) &&
        ((over_error > 0  && dt[i, actual_over] == 1) ||
           (over_error < 0  && dt[i, actual_over] == 0))
      
      if (improves_cover || improves_over) {
        dt[i, included := TRUE]
        addition_failed <- FALSE
        
        # Recompute BOTH errors from scratch with fixed N=500
        cover_error <- dt[included == TRUE, sum(actual_cover)] -
          round(target_cover * N)
        over_error  <- dt[included == TRUE, sum(actual_over )] -
          round(target_over  * N)
        break
      }
    }
  }
  
  # 4c) If neither removal nor addition helped, exit
  if (removal_failed && addition_failed) {
    break
  }
}

# 5) Final diagnostics
final_N      <- sum(dt[, included])
final_covers <- dt[included == TRUE, sum(actual_cover)]
final_overs  <- dt[included == TRUE, sum(actual_over )]

cat("Final sample size:   ", final_N, "\n")
cat("Final cover rate:    ", final_covers, " / ", final_N, " = ",
    round(final_covers / final_N, 4), "\n")
cat("Final over rate:     ", final_overs,  " / ", final_N, " = ",
    round(final_overs / final_N, 4), "\n")
cat("Residual cover_error:", cover_error, "\n")
cat("Residual over_error:  ", over_error,  "\n")
