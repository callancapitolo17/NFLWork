#!/usr/bin/env Rscript
# MLB Triple-Play Pricer
#
# Prices "SCR 1ST, 1H & GM" props: team scores first in the game AND wins the
# 1st half (F5, strict lead) AND wins the game. Uses the dispersion-matched
# samples from mlb_game_samples — same framework as mlb_correlated_parlay.R.
#
# Usage: Rscript mlb_triple_play.R

suppressPackageStartupMessages({
  library(duckdb)
  library(dplyr)
  library(tibble)
  library(DBI)
})

# =============================================================================
# CORE PRICER FUNCTIONS (pure — unit-tested)
# =============================================================================

#' Fair probability for a team's triple-play, computed empirically on sample rows.
#'
#' @param samples data.frame with columns home_margin, home_margin_f5,
#'   home_scored_first (0/1/NA)
#' @param side "home" or "away"
#' @return scalar probability in [0, 1], or NA if no valid rows
compute_triple_play_fair <- function(samples, side = c("home", "away")) {
  side <- match.arg(side)
  samples <- samples[!is.na(samples$home_scored_first), ]
  if (nrow(samples) == 0) return(NA_real_)
  if (side == "home") {
    mean(samples$home_scored_first == 1L &
         samples$home_margin_f5   > 0 &
         samples$home_margin      > 0)
  } else {
    mean(samples$home_scored_first == 0L &
         samples$home_margin_f5   < 0 &
         samples$home_margin      < 0)
  }
}

#' Probability → American odds (integer). Returns NA for p <= 0 or p >= 1.
prob_to_american <- function(p) {
  if (is.na(p) || p <= 0 || p >= 1) return(NA_integer_)
  if (p >= 0.5) return(as.integer(round(-100 * p / (1 - p))))
  as.integer(round(100 * (1 - p) / p))
}

#' American odds → implied probability (includes vig). Returns NA_real_ for
#' NA input or 0 (0 is not a valid American odds value).
american_to_prob <- function(o) {
  if (is.na(o) || o == 0) return(NA_real_)
  if (o > 0) return(100 / (o + 100))
  (-o) / ((-o) + 100)
}

# =============================================================================
# MAIN (guarded — filled in by Task 4)
# =============================================================================
if (!interactive() && sys.nframe() == 0L) {
  invisible(NULL)  # Task 4 replaces this block with the actual pricer invocation
}
