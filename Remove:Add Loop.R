improved_cover <- FALSE
improved_over <- FALSE
repeat{
  if(abs(cover_error) < tol_error & abs(over_error) < tol_error){
    break
  }
  removed <- FALSE
  if(abs(cover_error) > tol_error | abs(over_error) > tol_error){
    for(i in seq(N, 1)){
      if(dt[i,"included"] == TRUE){
        improved_cover <- as.integer(cover_error > 0) == dt[i,actual_cover]
        improved_over <- as.integer(over_error > 0) == dt[i,actual_over] 
        }
      if(improved_cover | improved_over){
        dt[i, "included"] <- FALSE
        removed <- TRUE
        break
      }
    }
  }
  if (removed) {
    cover_error = sum(dt[included == TRUE, actual_cover]) - round(target_cover * N)
    over_error  = sum(dt[included == TRUE, actual_over ]) - round(target_over  * N)
    next
  } 
  
  added <- FALSE
  if(abs(cover_error) > tol_error | abs(over_error) > tol_error){
    for(i in seq(N+1, nrow(dt))){
      if(dt[i,"included"] == FALSE){
        improved_cover <- as.integer(cover_error > 0) != dt[i,actual_cover]
        improved_over <- as.integer(over_error > 0) != dt[i,actual_over] }
      if(improved_cover | improved_over){
        dt[i, "included"] <- TRUE
        added <- TRUE
        break
      }
    }
  }
  if (added) {
    cover_error = sum(dt[included == TRUE, actual_cover]) - round(target_cover * N)
    over_error  = sum(dt[included == TRUE, actual_over ]) - round(target_over  * N)
    next
  }
  break
}