# Answer Keys/triple_play_helpers.R

#' Determine if the home team scored first in the game.
#'
#' Baseball half-inning order: away bats top (first), home bats bottom (second).
#' So if the away team scored any runs in the first inning where scoring
#' occurred, away scored first — even if home also scored in the bottom half.
#'
#' Inputs are cumulative inning margins (home-minus-away) and totals after
#' innings 1 through 5.
#'
#' @return 1L if home scored first, 0L if away scored first, NA_integer_ if
#'   no scoring through inning 5 (rare — ~5% of games).
determine_home_scored_first <- function(m1, t1, m2, t2, m3, t3, m4, t4, m5, t5) {
  m <- c(m1, m2, m3, m4, m5)
  t <- c(t1, t2, t3, t4, t5)
  prev_m <- 0
  prev_t <- 0
  for (i in seq_len(5)) {
    if (is.na(t[i]) || is.na(m[i])) return(NA_integer_)
    if (t[i] > prev_t) {
      dt <- t[i] - prev_t
      dm <- m[i] - prev_m
      # Defensive: impossible-input guard. Total runs in an inning must be >=
      # absolute margin change. Corrupt inputs return NA rather than a
      # plausible-looking wrong answer.
      if (abs(dm) > dt) return(NA_integer_)
      away_runs <- (dt - dm) / 2
      home_runs <- (dt + dm) / 2
      if (away_runs > 0) return(0L)
      if (home_runs > 0) return(1L)
      return(NA_integer_)
    }
    prev_m <- m[i]
    prev_t <- t[i]
  }
  NA_integer_
}

#' Vectorized wrapper for use inside data.table/dplyr mutates.
determine_home_scored_first_vec <- function(m1, t1, m2, t2, m3, t3, m4, t4, m5, t5) {
  mapply(determine_home_scored_first,
         m1, t1, m2, t2, m3, t3, m4, t4, m5, t5,
         USE.NAMES = FALSE)
}
