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
# MAIN
# =============================================================================

if (!interactive() && sys.nframe() == 0L) {

  setwd("~/NFLWork/Answer Keys")
  source("parse_legs.R")
  MLB_DB <- "mlb.duckdb"

  # Today's book lines — edit this tribble whenever new props post.
  # home_team / away_team must match Odds API canonical names in mlb_consensus_temp.
  # description is the Wagerzon label verbatim; parse_legs() derives leg logic from it.
  todays_lines <- tribble(
    ~home_team,              ~away_team,             ~target_team, ~side,   ~book_odds, ~description,
    "Colorado Rockies",      "San Diego Padres",     "Rockies",    "home",  +530,       "ROCKIES TRIPLE-PLAY (SCR 1ST, 1H & GM)",
    "Colorado Rockies",      "San Diego Padres",     "Padres",     "away",  +190,       "PADRES TRIPLE-PLAY (SCR 1ST, 1H & GM)",
    "San Francisco Giants",  "Los Angeles Dodgers",  "Giants",     "home",  +750,       "GIANTS TRIPLE-PLAY (SCR 1ST, 1H & GM)",
    "San Francisco Giants",  "Los Angeles Dodgers",  "Dodgers",    "away",  +155,       "DODGERS TRIPLE-PLAY (SCR 1ST, 1H & GM)",
    "Seattle Mariners",      "Athletics",            "Mariners",   "home",  +215,       "MARINERS TRIPLE-PLAY (SCR 1ST, 1H & GM)",
    "Seattle Mariners",      "Athletics",            "Athletics",  "away",  +455,       "ATHLETICS TRIPLE-PLAY (SCR 1ST, 1H & GM)",
    "Arizona Diamondbacks",  "Chicago White Sox",    "DBacks",     "home",  +240,       "DBACKS TRIPLE-PLAY (SCR 1ST, 1H & GM)",
    "Arizona Diamondbacks",  "Chicago White Sox",    "White Sox",  "away",  +415,       "WHITE SOX TRIPLE-PLAY (SCR 1ST, 1H & GM)"
  )

  con <- dbConnect(duckdb(), dbdir = MLB_DB, read_only = TRUE)
  on.exit(tryCatch(dbDisconnect(con), error = function(e) NULL), add = TRUE)

  samples_df <- dbGetQuery(con,
    "SELECT game_id, home_margin, total_final_score,
            home_margin_f3, home_margin_f5, home_margin_f7,
            home_scored_first
     FROM mlb_game_samples")
  consensus  <- dbGetQuery(con,
    "SELECT id, home_team, away_team, commence_time FROM mlb_consensus_temp")

  if (!"home_scored_first" %in% names(samples_df)) {
    stop("mlb_game_samples is missing home_scored_first. Re-run MLB.R to regenerate.")
  }

  # Restrict consensus to games starting within 12 hours. mlb_consensus_temp
  # now carries multiple days of slate; without this filter, doubleheaders +
  # next-day games cause each tribble row to match >1 game_id.
  consensus <- consensus %>%
    filter(commence_time > Sys.time() &
           commence_time < Sys.time() + 12 * 3600) %>%
    select(-commence_time)

  # Join today's lines to their game_id via team names
  matched <- todays_lines %>%
    inner_join(consensus, by = c("home_team", "away_team"))

  n_matched <- nrow(matched)
  n_posted  <- nrow(todays_lines)
  if (n_matched < n_posted) {
    dropped <- anti_join(todays_lines, consensus, by = c("home_team", "away_team"))
    warning(sprintf("Dropped %d/%d lines with no matching consensus game:\n%s",
                    n_posted - n_matched, n_posted,
                    paste(capture.output(print(dropped)), collapse = "\n")))
  }

  # Price each line
  priced <- matched %>%
    rowwise() %>%
    mutate(
      game_samples = list(samples_df[samples_df$game_id == id, ]),
      n_samples    = nrow(game_samples),
      legs         = list(parse_legs(description)),
      fair_prob    = compute_prop_fair(game_samples, side, legs),
      fair_odds    = prob_to_american(fair_prob),
      book_prob    = american_to_prob(book_odds),
      edge_pct     = (fair_prob / book_prob - 1) * 100
    ) %>%
    ungroup() %>%
    select(target_team, side, n_samples, fair_prob, fair_odds, book_odds, edge_pct) %>%
    arrange(desc(edge_pct))

  cat("\n=== MLB Triple-Play Fair Prices (SCR 1ST + F5 + GM) ===\n")
  cat(sprintf("Matched %d / %d posted lines.\n\n", n_matched, n_posted))

  display <- priced %>%
    mutate(
      fair_prob = sprintf("%.3f", fair_prob),
      fair_odds = ifelse(fair_odds > 0,
                         paste0("+", fair_odds), as.character(fair_odds)),
      book_odds = ifelse(book_odds > 0,
                         paste0("+", book_odds), as.character(book_odds)),
      edge_pct  = sprintf("%+.1f%%", edge_pct)
    )
  print(as.data.frame(display), row.names = FALSE)
}
