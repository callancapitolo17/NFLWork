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

WRite R for this please using the base code I have below:
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

cover_error = sum(dt[included == TRUE, actual_cover]) - round(target_cover * N)
over_error  = sum(dt[included == TRUE, actual_over ]) - round(target_over  * N)

repeat{
  cover_adjustments <- abs(cover_error)
  over_adjustments <- abs(over_error)
  
  if(cover_adjustments == 0 & over_adjustments == 0 ){
    break
  }
  removed <- FALSE
  if(abs(cover_error) > tol_error | abs(over_error) > tol_error){
    for(i in seq(N, 1)){
      if(dt[i,"included"] == TRUE){
        improved_cover <- as.integer(cover_error > 0) == dt[i,actual_cover]
        improved_over <- as.integer(over_error > 0) == dt[i,actual_over] 
        print(c(improved_cover,improved_over))
      }
      if(improved_cover | improved_over){
        dt[i, "included"] <- FALSE
        removed <- TRUE
      }
      cover_adjustments <- ifelse(improved_cover,cover_adjustments - 1, cover_adjustments +1)
      over_adjustments <- ifelse(improved_cover,over_adjustments - 1, over_adjustments +1)
      cover_error  = sum(dt[included == TRUE, actual_cover ]) - round(target_cover  * nrow(dt[included == TRUE,]))
      over_error  = sum(dt[included == TRUE, actual_over ]) - round(target_over  * nrow(dt[included == TRUE,]))
    }
    if(i == 1){break}
    
  }}
  
  if (!is.na(best_idx) && best_improvement > 0) {
    dt[best_idx, included := FALSE]
  } else {
    # No improvement possible, exit
    break
  }
}

  
#   added <- FALSE
#   if(abs(cover_error) > tol_error | abs(over_error) > tol_error){
#     for(i in seq(N+1, nrow(dt))){
#       if(dt[i,"included"] == FALSE){
#         improved_cover <- as.integer(cover_error > 0) != dt[i,actual_cover]
#         improved_over <- as.integer(over_error > 0) != dt[i,actual_over] }
#       if(improved_cover | improved_over){
#         dt[i, "included"] <- TRUE
#         added <- TRUE
#         break
#       }
#     }
#   }
#   if (added) {
#     cover_error = sum(dt[included == TRUE, actual_cover]) - round(target_cover * N)
#     over_error  = sum(dt[included == TRUE, actual_over ]) - round(target_over  * N)
#     next
#   }
#   break
# }