
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
cover_adjustments <- abs(cover_error)
over_adjustments <- abs(over_error)
for (i in seq(N,1)) {
  if(over_adjustments == 0 & cover_adjustments == 0){break}
  if(dt[i,"included"] == TRUE){
    improved_cover <- (as.integer(cover_error > 0) == dt[i, actual_cover])
    improved_over  <- (as.integer(over_error  > 0) == dt[i, actual_over ])
    if(improved_over | improved_cover){
      dt[i,"included"] <- FALSE
      over_adjustments <- ifelse(improved_over, over_adjustments -1, over_adjustments+1)
      cover_adjustments <- ifelse(improved_cover, cover_adjustments -1, cover_adjustments+1)
      cover_error <- sum(dt[included == TRUE, actual_cover]) - round(target_cover * sum(dt[, included]) )
      over_error  <- sum(dt[included == TRUE, actual_over])  - round(target_over  * sum(dt[, included]) )
      print(c(over_adjustments,cover_adjustments,i))
    }
  }
}

#Look into how to keep track of adjustments, use spreadsheet