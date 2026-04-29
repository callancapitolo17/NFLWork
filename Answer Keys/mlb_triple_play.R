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
  library(jsonlite)
  library(digest)   # for trifecta_hash
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

  # Ensure DK trifecta SGP table exists (idempotent; created lazily on first
  # pricer run). The Plan #2 scraper writes here; Plan #1 reads it (empty).
  con_mig <- dbConnect(duckdb(), dbdir = MLB_DB)
  dbExecute(con_mig, "
    CREATE TABLE IF NOT EXISTS mlb_trifecta_sgp_odds (
      fetch_time     TIMESTAMP,
      game_id        VARCHAR,
      prop_type      VARCHAR,
      side           VARCHAR,
      legs_json      VARCHAR,
      selection_ids  VARCHAR,
      sgp_decimal    DOUBLE,
      sgp_american   INTEGER,
      source         VARCHAR
    )
  ")
  dbDisconnect(con_mig)

  # Read today's posted specials from the wagerzon_specials scraper output.
  # Uses the most recent snapshot AND filters to upcoming games (game_time > NOW()).
  # The game_time filter is the real safeguard against stale data: if WZ has no
  # specials posted today, MAX(scraped_at) returns yesterday's snapshot, but every
  # row in it has a past game_time and gets filtered out. Off-day → zero rows →
  # empty Trifectas tab (correct behavior).
  WZ_DB <- "~/NFLWork/wagerzon_odds/wagerzon.duckdb"
  wz_con <- dbConnect(duckdb(), dbdir = path.expand(WZ_DB), read_only = TRUE)
  on.exit(tryCatch(dbDisconnect(wz_con), error = function(e) NULL), add = TRUE)
  specials <- dbGetQuery(wz_con, "
    SELECT team, prop_type, description, odds AS book_odds
    FROM wagerzon_specials
    WHERE sport = 'mlb'
      AND prop_type IN ('TRIPLE-PLAY', 'GRAND-SLAM')
      AND scraped_at = (SELECT MAX(scraped_at) FROM wagerzon_specials WHERE sport = 'mlb')
      AND game_time > NOW()
  ")
  dbDisconnect(wz_con)

  if (nrow(specials) == 0) {
    cat("No priceable specials found in wagerzon_specials. Run scraper_specials.py first.\n")
    # Clear mlb_trifecta_opportunities so the dashboard reflects "nothing
    # priceable today" instead of leftover rows from a previous run. Without
    # this, a prior run's rows linger with future game_time values and the
    # dashboard's load filter (game_time > NOW()) can't distinguish them
    # from genuinely fresh pricing.
    cleanup_con <- tryCatch(
      dbConnect(duckdb(), dbdir = MLB_DB),
      error = function(e) NULL
    )
    if (!is.null(cleanup_con)) {
      tryCatch({
        dbExecute(cleanup_con, "DROP TABLE IF EXISTS mlb_trifecta_opportunities")
      }, error = function(e) {
        cat(sprintf("Warning: failed to clear mlb_trifecta_opportunities: %s\n", e$message))
      })
      duckdb::dbDisconnect(cleanup_con, shutdown = TRUE)
    }
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
           commence_time < Sys.time() + 12 * 3600)

  todays_lines <- specials %>%
    inner_join(consensus, by = c("canonical_team" = "home_team"),
               relationship = "many-to-many", keep = TRUE) %>%
    mutate(side = "home", target_team = team,
           home_team = canonical_team) %>%
    select(home_team, away_team, target_team, side, book_odds, description, prop_type, id, commence_time) %>%
    bind_rows(
      specials %>%
        inner_join(consensus, by = c("canonical_team" = "away_team"),
                   relationship = "many-to-many", keep = TRUE) %>%
        mutate(side = "away", target_team = team,
               away_team = canonical_team) %>%
        select(home_team, away_team, target_team, side, book_odds, description, prop_type, id, commence_time)
    )

  matched <- todays_lines  # consensus join already done; keep var name compatible

  # Derive prop_type from the description so the scraper invocation can
  # send it (the rowwise mutate below previously did this — moved up so
  # the scraper request includes it).
  matched <- matched %>%
    mutate(
      prop_type = vapply(description, function(d) {
        m <- regmatches(d, regexec("\\b(TRIPLE-PLAY|GRAND-SLAM)\\b", d))[[1]]
        if (length(m) >= 2) m[[2]] else NA_character_
      }, character(1))
    )

  n_matched <- nrow(matched)
  n_posted  <- nrow(specials)
  if (n_matched < n_posted) {
    warning(sprintf("Matched %d/%d posted specials. Some teams may have no game tonight.",
                    n_matched, n_posted))
  }

  # ===== DK trifecta SGP scraper invocation (Plan #2) =====
  # Build a JSON request file with one entry per matched line, then call the
  # Python scraper via system2(). Scraper writes rows to mlb_trifecta_sgp_odds;
  # we read them back below in the dk_sgp <- dbGetQuery(...) block.
  #
  # If the scraper exits non-zero, we log a warning and continue — the blend
  # gracefully degrades to model-only when the table has no fresh rows.
  SGP_VENV_PYTHON <- "/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python"
  SGP_DIR <- "/Users/callancapitolo/NFLWork/mlb_sgp"

  # Build request list: keep legs as a native R list so write_json serializes
  # it as a JSON array (not a double-encoded string).
  trifecta_list <- lapply(seq_len(nrow(matched)), function(i) {
    row <- matched[i, ]
    list(
      game_id    = row$id,
      home_team  = row$home_team,
      away_team  = row$away_team,
      prop_type  = row$prop_type,
      side       = row$side,
      legs       = parse_legs(row$description)
    )
  })
  request_path <- tempfile(fileext = ".json")
  jsonlite::write_json(trifecta_list, request_path, auto_unbox = TRUE)

  # Close mlb.duckdb connection before invoking the Python scraper.
  # DuckDB doesn't allow concurrent read-only + write access on the same
  # file; the scraper needs exclusive write access while it runs.
  dbDisconnect(con)

  dk_status <- tryCatch({
    system2(SGP_VENV_PYTHON,
            args = c(file.path(SGP_DIR, "scraper_draftkings_trifecta.py"),
                     "--input", request_path, "--db", MLB_DB),
            stdout = TRUE, stderr = TRUE)
    0
  }, warning = function(w) {
    cat(sprintf("DK trifecta scraper warning: %s\n", conditionMessage(w)))
    1
  }, error = function(e) {
    cat(sprintf("DK trifecta scraper failed: %s\n", conditionMessage(e)))
    1
  })
  if (dk_status != 0) {
    cat("DK trifecta blend disabled for this run; using model-only fair odds.\n")
  }
  unlink(request_path)
  # ========================================================

  # Reopen mlb.duckdb (read-only) now that the scraper has finished writing.
  # query sees the rows just written by the Python process.
  con <- dbConnect(duckdb(), dbdir = MLB_DB, read_only = TRUE)
  on.exit(tryCatch(dbDisconnect(con), error = function(e) NULL), add = TRUE)
  dk_sgp <- tryCatch({
    dbGetQuery(con, "
      SELECT game_id, prop_type, side, sgp_decimal
      FROM mlb_trifecta_sgp_odds
      WHERE source = 'draftkings_direct'
        AND fetch_time = (SELECT MAX(fetch_time) FROM mlb_trifecta_sgp_odds)
    ")
  }, error = function(e) data.frame(
    game_id     = character(0),
    prop_type   = character(0),
    side        = character(0),
    sgp_decimal = double(0)
  ))

  # DK SGP vig fallback for trifectas. Higher than mlb_correlated_parlay.R's
  # 1.10 because trifectas are 3-4 legs (TRIPLE-PLAY = 3, GRAND-SLAM = 4) and
  # DK shaves more aggressively as legs grow. With only 2 observations per
  # game (home + away) we cannot fit per-game vig empirically, so we use a
  # conservative fixed estimate. Revisit once Plan #2 collects real DK data.
  DK_SGP_VIG_DEFAULT <- 1.25

  # Price each line
  priced <- matched %>%
    rowwise() %>%
    mutate(
      game_samples    = list(samples_df[samples_df$game_id == id, ]),
      n_samples       = nrow(game_samples),
      legs            = list(parse_legs(description)),
      model_fair_prob = compute_prop_fair(game_samples, side, legs),

      dk_match        = list(dk_sgp[dk_sgp$game_id   == id        &
                                    dk_sgp$prop_type == prop_type &
                                    dk_sgp$side      == side, ]),
      dk_decimal      = if (nrow(dk_match) > 0) dk_match$sgp_decimal[1] else NA_real_,
      dk_fair_prob    = if (!is.na(dk_decimal) && dk_decimal > 0) {
                          (1 / dk_decimal) / DK_SGP_VIG_DEFAULT
                        } else NA_real_,

      fair_prob       = blend_dk_with_model(model_fair_prob, dk_decimal,
                                            DK_SGP_VIG_DEFAULT),

      model_odds      = prob_to_american(model_fair_prob),
      dk_odds         = prob_to_american(dk_fair_prob),
      fair_odds       = prob_to_american(fair_prob),
      book_prob       = american_to_prob(book_odds),
      edge_pct        = (fair_prob / book_prob - 1) * 100
    ) %>%
    ungroup() %>%
    select(id, target_team, prop_type, side, description,
           home_team, away_team, commence_time, n_samples,
           model_odds, dk_odds, fair_odds, book_odds, edge_pct) %>%
    arrange(desc(edge_pct))

  cat("\n=== MLB Triple-Play Fair Prices (SCR 1ST + F5 + GM) ===\n")
  cat(sprintf("Matched %d / %d posted lines.\n\n", n_matched, n_posted))

  display <- priced %>%
    mutate(
      model_odds = ifelse(is.na(model_odds), NA_character_,
                          ifelse(model_odds > 0,
                                 paste0("+", model_odds), as.character(model_odds))),
      dk_odds    = ifelse(is.na(dk_odds), NA_character_,
                          ifelse(dk_odds > 0,
                                 paste0("+", dk_odds), as.character(dk_odds))),
      fair_odds  = ifelse(is.na(fair_odds), NA_character_,
                          ifelse(fair_odds > 0,
                                 paste0("+", fair_odds), as.character(fair_odds))),
      book_odds  = ifelse(book_odds > 0,
                          paste0("+", book_odds), as.character(book_odds)),
      edge_pct   = sprintf("%+.1f%%", edge_pct)
    )
  print(as.data.frame(display), row.names = FALSE)

  # =============================================================================
  # WRITE TO DUCKDB (for dashboard consumption)
  # =============================================================================

  # Read trifecta sizing settings from the dashboard DB (same pattern as
  # mlb_correlated_parlay.R reads parlay sizing). Falls back to safe defaults
  # if the dashboard DB or rows are missing.
  trifecta_bankroll   <- 100
  trifecta_kelly_mult <- 0.10
  trifecta_min_edge   <- 0.05
  dash_db_path <- path.expand("~/NFLWork/Answer Keys/MLB Dashboard/mlb_dashboard.duckdb")
  if (file.exists(dash_db_path)) {
    dash_con <- tryCatch(
      dbConnect(duckdb(), dbdir = dash_db_path, read_only = TRUE),
      error = function(e) NULL
    )
    if (!is.null(dash_con)) {
      saved <- tryCatch(
        dbGetQuery(dash_con, "SELECT param, value FROM sizing_settings"),
        error = function(e) data.frame(param = character(0), value = numeric(0))
      )
      if ("trifecta_bankroll"   %in% saved$param) trifecta_bankroll   <- saved$value[saved$param == "trifecta_bankroll"]
      if ("trifecta_kelly_mult" %in% saved$param) trifecta_kelly_mult <- saved$value[saved$param == "trifecta_kelly_mult"]
      if ("trifecta_min_edge"   %in% saved$param) trifecta_min_edge   <- saved$value[saved$param == "trifecta_min_edge"]
      dbDisconnect(dash_con)
    }
  }

  # Compute Kelly per row. Trifectas are single-ticket bets from the
  # operator's POV (one ticket per row), so independent Kelly is correct.
  # Filter via min_edge: rows below the threshold get kelly_bet = 0.
  priced <- priced %>%
    rowwise() %>%
    mutate(
      win_prob   = if (!is.na(fair_odds)) american_to_prob(fair_odds) else NA_real_,
      dec_odds   = if (!is.na(book_odds)) {
                     if (book_odds > 0) 1 + book_odds / 100 else 1 + 100 / abs(book_odds)
                   } else NA_real_,
      kelly_frac = if (!is.na(win_prob) && !is.na(dec_odds) && dec_odds > 1) {
                     b <- dec_odds - 1
                     p <- win_prob
                     q <- 1 - p
                     max(0, (b * p - q) / b)
                   } else 0,
      kelly_bet  = if (!is.na(edge_pct) && edge_pct >= trifecta_min_edge * 100) {
                     trifecta_bankroll * trifecta_kelly_mult * kelly_frac
                   } else 0
    ) %>%
    ungroup() %>%
    select(-win_prob, -dec_odds, -kelly_frac)

  # Add display columns + dedup hash
  priced <- priced %>%
    mutate(
      game          = sprintf("%s @ %s", away_team, home_team),
      game_time     = commence_time,
      game_id       = as.character(id),
      trifecta_hash = sapply(seq_len(n()), function(i) {
        digest::digest(
          paste(game_id[i], target_team[i], prop_type[i], side[i], sep = "|"),
          algo = "sha256", serialize = FALSE
        )
      })
    ) %>%
    select(trifecta_hash, game_id, game, game_time, target_team, prop_type, side,
           description, n_samples, model_odds, dk_odds, fair_odds, book_odds,
           edge_pct, kelly_bet)

  # Close the read-only mlb.duckdb connection before opening it for writes.
  # DuckDB does not allow a simultaneous read-only + read-write connection to
  # the same database file. The on.exit() handlers registered earlier will
  # attempt dbDisconnect(con) at script exit — tryCatch there silences the
  # "already closed" error harmlessly.
  tryCatch(dbDisconnect(con), error = function(e) NULL)

  # Drop + rewrite (same pattern as mlb_correlated_parlay.R)
  write_con <- NULL
  tryCatch({
    write_con <- dbConnect(duckdb(), dbdir = MLB_DB)
    on.exit(if (!is.null(write_con)) duckdb::dbDisconnect(write_con, shutdown = TRUE), add = TRUE)
    dbExecute(write_con, "DROP TABLE IF EXISTS mlb_trifecta_opportunities")
    dbWriteTable(write_con, "mlb_trifecta_opportunities", priced)
    cat(sprintf("Wrote %d trifecta opportunities to %s.\n", nrow(priced), MLB_DB))
  }, error = function(e) {
    cat(sprintf("Warning: Failed to write trifectas to DB: %s\n", e$message))
  })

}
