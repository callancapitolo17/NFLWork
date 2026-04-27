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

  # Read today's posted specials from the wagerzon_specials scraper output.
  # Uses the most recent snapshot. Pricer is robust to empty / off-day cases:
  # if zero rows are posted, prints a clear message and exits.
  WZ_DB <- "~/NFLWork/wagerzon_odds/wagerzon.duckdb"
  wz_con <- dbConnect(duckdb(), dbdir = path.expand(WZ_DB), read_only = TRUE)
  on.exit(tryCatch(dbDisconnect(wz_con), error = function(e) NULL), add = TRUE)
  specials <- dbGetQuery(wz_con, "
    SELECT team, prop_type, description, odds AS book_odds
    FROM wagerzon_specials
    WHERE sport = 'mlb'
      AND prop_type IN ('TRIPLE-PLAY', 'GRAND-SLAM')
      AND scraped_at = (SELECT MAX(scraped_at) FROM wagerzon_specials WHERE sport = 'mlb')
  ")
  dbDisconnect(wz_con)

  if (nrow(specials) == 0) {
    cat("No priceable specials found in wagerzon_specials. Run scraper_specials.py first.\n")
    quit(status = 0)
  }

  # Translate Wagerzon team names (UPPER) to Odds API canonical names so the
  # consensus join works. Map covers the 30 MLB teams.
  WZ_TO_CANONICAL <- c(
    "ANGELS"        = "Los Angeles Angels",
    "ASTROS"        = "Houston Astros",
    "ATHLETICS"     = "Athletics",
    "BLUE JAYS"     = "Toronto Blue Jays",
    "BRAVES"        = "Atlanta Braves",
    "BREWERS"       = "Milwaukee Brewers",
    "CARDINALS"     = "St. Louis Cardinals",
    "CUBS"          = "Chicago Cubs",
    "DBACKS"        = "Arizona Diamondbacks",
    "DIAMONDBACKS"  = "Arizona Diamondbacks",
    "DODGERS"       = "Los Angeles Dodgers",
    "GIANTS"        = "San Francisco Giants",
    "GUARDIANS"     = "Cleveland Guardians",
    "MARINERS"      = "Seattle Mariners",
    "MARLINS"       = "Miami Marlins",
    "METS"          = "New York Mets",
    "NATIONALS"     = "Washington Nationals",
    "ORIOLES"       = "Baltimore Orioles",
    "PADRES"        = "San Diego Padres",
    "PHILLIES"      = "Philadelphia Phillies",
    "PIRATES"       = "Pittsburgh Pirates",
    "RANGERS"       = "Texas Rangers",
    "RAYS"          = "Tampa Bay Rays",
    "RED SOX"       = "Boston Red Sox",
    "REDS"          = "Cincinnati Reds",
    "ROCKIES"       = "Colorado Rockies",
    "ROYALS"        = "Kansas City Royals",
    "TIGERS"        = "Detroit Tigers",
    "TWINS"         = "Minnesota Twins",
    "WHITE SOX"     = "Chicago White Sox",
    "YANKEES"       = "New York Yankees"
  )

  specials$canonical_team <- WZ_TO_CANONICAL[specials$team]
  unmapped <- specials[is.na(specials$canonical_team), ]
  if (nrow(unmapped) > 0) {
    warning(sprintf("Dropped %d specials with unmapped team names: %s",
                    nrow(unmapped),
                    paste(unique(unmapped$team), collapse = ", ")))
    specials <- specials[!is.na(specials$canonical_team), ]
  }

  # Resolve home/away by joining to mlb_consensus_temp (same approach used
  # for the prior tribble path). For each canonical_team, find the consensus
  # row where it's home OR away.
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

  # 12-hour filter (same as before): mlb_consensus_temp carries multiple days
  consensus <- consensus %>%
    filter(commence_time > Sys.time() &
           commence_time < Sys.time() + 12 * 3600) %>%
    select(-commence_time)

  todays_lines <- specials %>%
    inner_join(consensus, by = c("canonical_team" = "home_team"),
               relationship = "many-to-many", keep = TRUE) %>%
    mutate(side = "home", target_team = team,
           home_team = canonical_team) %>%
    select(home_team, away_team, target_team, side, book_odds, description, id) %>%
    bind_rows(
      specials %>%
        inner_join(consensus, by = c("canonical_team" = "away_team"),
                   relationship = "many-to-many", keep = TRUE) %>%
        mutate(side = "away", target_team = team,
               away_team = canonical_team) %>%
        select(home_team, away_team, target_team, side, book_odds, description, id)
    )

  matched <- todays_lines  # consensus join already done; keep var name compatible
  n_matched <- nrow(matched)
  n_posted  <- nrow(specials)
  if (n_matched < n_posted) {
    warning(sprintf("Matched %d/%d posted specials. Some teams may have no game tonight.",
                    n_matched, n_posted))
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
      edge_pct     = (fair_prob / book_prob - 1) * 100,
      prop_type    = {
        # Anchor on known prop-type tokens so multi-word team names work
        # ("WHITE SOX TRIPLE-PLAY ..." and "GIANTS GRAND-SLAM ..." both parse).
        # Extend this pattern as new prop types are added to TOKEN_REGISTRY.
        known_props <- "(TRIPLE-PLAY|GRAND-SLAM)"
        m <- regmatches(description,
                        regexec(paste0("\\b", known_props, "\\b"), description))[[1]]
        if (length(m) >= 2) m[[2]] else NA_character_
      }
    ) %>%
    ungroup() %>%
    select(target_team, prop_type, side, n_samples,
           fair_prob, fair_odds, book_odds, edge_pct) %>%
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
