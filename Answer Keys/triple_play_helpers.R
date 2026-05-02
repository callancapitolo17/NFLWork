# Answer Keys/triple_play_helpers.R

#' Determine whether each team scored ≥1 run in the 1st inning.
#'
#' Wagerzon's "SCR 1ST" prop pays out when the named team scored at least one
#' run during the 1st inning specifically — independent of the opposing team.
#' Both teams can win their respective SCR 1ST bets in the same game.
#'
#' Computed from cumulative-after-inning-1 fields:
#'   t1 = total runs scored in inning 1 (away_runs_1 + home_runs_1)
#'   m1 = home margin after inning 1 (home_runs_1 - away_runs_1)
#' which gives:
#'   away_runs_1 = (t1 - m1) / 2
#'   home_runs_1 = (t1 + m1) / 2
#'
#' @return list(home = 0L|1L|NA_integer_, away = 0L|1L|NA_integer_).
#'   NA when either input is NA, or when |m1| > t1 (impossible inputs —
#'   corrupt PBP data, return NA rather than a plausible-looking wrong answer).
determine_inning_1_scoring <- function(m1, t1) {
  if (is.na(m1) || is.na(t1) || abs(m1) > t1) {
    return(list(home = NA_integer_, away = NA_integer_))
  }
  away_runs_1 <- (t1 - m1) / 2
  home_runs_1 <- (t1 + m1) / 2
  list(
    home = as.integer(home_runs_1 > 0),
    away = as.integer(away_runs_1 > 0)
  )
}

#' Vectorized wrapper. Returns list(home = integer vector, away = integer vector).
determine_inning_1_scoring_vec <- function(m1, t1) {
  results <- mapply(determine_inning_1_scoring, m1, t1,
                    SIMPLIFY = FALSE, USE.NAMES = FALSE)
  list(
    home = vapply(results, `[[`, integer(1), "home"),
    away = vapply(results, `[[`, integer(1), "away")
  )
}
