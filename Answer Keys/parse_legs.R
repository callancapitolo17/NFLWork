# Answer Keys/parse_legs.R
# Wagerzon prop-description parser and generic leg-based joint probability pricer.

#' Convert a Wagerzon fraction string to a numeric.
#'
#' Handles three formats Wagerzon uses interchangeably:
#'   "3"   → 3
#'   "2.5" → 2.5
#'   "2½"  → 2.5 (using the U+00BD unicode fraction)
#' Also handles ¼ (U+00BC → 0.25) and ¾ (U+00BE → 0.75) for future-proofing.
#'
#' Returns NA_real_ on unparseable input.
parse_fraction <- function(s) {
  if (is.null(s) || is.na(s) || !nzchar(s)) return(NA_real_)
  frac_map <- c("¼" = 0.25, "½" = 0.5, "¾" = 0.75)
  # Try pure numeric first (handles "3", "2.5")
  n <- suppressWarnings(as.numeric(s))
  if (!is.na(n)) return(n)
  # Check for unicode fraction suffix
  last_char <- substr(s, nchar(s), nchar(s))
  if (last_char %in% names(frac_map)) {
    int_part <- suppressWarnings(as.numeric(substr(s, 1L, nchar(s) - 1L)))
    if (is.na(int_part)) int_part <- 0
    return(int_part + frac_map[[last_char]])
  }
  NA_real_
}

#' Token registry: regex → function(match) returning a leg spec.
#'
#' Each entry has:
#'   pattern: POSIX-extended regex (anchored with ^ and $) matched against the
#'            stripped, uppercased token string.
#'   spec:    function(character vector of capture groups) returning a named list.
#'
#' Add new tokens here to teach the parser new prop types.
TOKEN_REGISTRY <- list(
  list(pattern = "^F3$",
       spec    = function(m) list(type = "wins_period", period = "F3")),
  list(pattern = "^1H$",
       spec    = function(m) list(type = "wins_period", period = "F5")),
  list(pattern = "^F7$",
       spec    = function(m) list(type = "wins_period", period = "F7")),
  list(pattern = "^GM$",
       spec    = function(m) list(type = "wins_period", period = "FG")),
  list(pattern = "^SCR 1ST$",
       spec    = function(m) list(type = "scores_first")),
  # For patterns with a capture group, m[[1]] is the full match and
  # m[[2]] is the first capture. parse_fraction needs the capture only.
  list(pattern = "^SCR U(.+)$",
       spec    = function(m) list(type = "team_total_under",
                                  line = parse_fraction(m[[2]]))),
  list(pattern = "^SCR O(.+)$",
       spec    = function(m) list(type = "team_total_over",
                                  line = parse_fraction(m[[2]])))
)

#' Match a single token against TOKEN_REGISTRY. Returns leg spec or NULL.
.match_token <- function(token) {
  for (entry in TOKEN_REGISTRY) {
    m <- regmatches(token, regexec(entry$pattern, token, perl = TRUE))[[1]]
    if (length(m) > 0) return(entry$spec(m))
  }
  NULL
}

#' Parse a Wagerzon description string into a list of leg specs.
#'
#' "GIANTS TRIPLE-PLAY (SCR 1ST, 1H & GM)" → 3 leg specs.
#' Unknown tokens trigger a warning and return NULL (pricer treats NULL as
#' "cannot price" → fair_prob = NA).
#' Descriptions with no parenthetical also return NULL.
parse_legs <- function(description) {
  if (is.null(description) || length(description) == 0 ||
      is.na(description) || !nzchar(description)) {
    return(NULL)
  }
  # Extract the parenthetical content (e.g. "SCR 1ST, 1H & GM")
  inside <- regmatches(description,
                       regexec("\\((.+)\\)", description, perl = TRUE))[[1]]
  if (length(inside) < 2) return(NULL)
  raw <- inside[[2]]
  # Tokenize: split on "," and "&", trim whitespace
  tokens <- trimws(unlist(strsplit(raw, "[,&]")))
  tokens <- tokens[nzchar(tokens)]
  # Map each token via registry
  legs <- vector("list", length(tokens))
  for (i in seq_along(tokens)) {
    spec <- .match_token(tokens[[i]])
    if (is.null(spec)) {
      warning(sprintf("Unknown leg token: '%s' in description '%s'",
                      tokens[[i]], description))
      return(NULL)
    }
    legs[[i]] <- spec
  }
  legs
}

#' Evaluate a single leg spec against a samples frame. Returns a logical vector
#' of length nrow(samples).
#'
#' @param leg Named list returned by parse_legs() (e.g. list(type="scores_first")).
#' @param samples data.frame with columns matching what the leg type needs.
#'   scores_first → home_scored_in_1st (home side) / away_scored_in_1st (away side)
#'   wins_period  → home_margin, home_margin_f3, home_margin_f5, or home_margin_f7
#'   team_total_* → nothing from samples directly (uses team_runs arg)
#'   opp_total_*  → uses opp_runs arg
#' @param side "home" or "away" — direction of inequalities and indicator value.
#' @param team_runs numeric vector of length nrow(samples), runs scored by the
#'   prop's subject team. Computed once per game by caller from home_margin +
#'   total_final_score.
#' @param opp_runs  numeric vector, runs scored by the opposing team.
eval_leg <- function(leg, samples, side, team_runs, opp_runs) {
  switch(leg$type,
    scores_first = {
      col <- if (side == "home") samples$home_scored_in_1st
             else                samples$away_scored_in_1st
      col == 1L
    },
    wins_period = {
      col <- switch(leg$period,
                    F3 = samples$home_margin_f3,
                    F5 = samples$home_margin_f5,
                    F7 = samples$home_margin_f7,
                    FG = samples$home_margin,
                    stop("Unknown period: ", leg$period))
      if (side == "home") col > 0 else col < 0
    },
    team_total_under = team_runs <  leg$line,
    team_total_over  = team_runs >  leg$line,
    opp_total_under  = opp_runs  <  leg$line,
    opp_total_over   = opp_runs  >  leg$line,
    stop("Unknown leg type: ", leg$type)
  )
}

#' Generic joint-probability pricer for any multi-leg prop.
#'
#' Drops legs' NA evaluator results after AND-reducing across samples;
#' derives team/opp runs from margin + total.
#'
#' @param samples data.frame with at minimum: home_margin, home_margin_f5,
#'   home_scored_in_1st, away_scored_in_1st, total_final_score. wins_period
#'   legs for F3/F7 also need home_margin_f3 / home_margin_f7 respectively.
#' @param side "home" or "away" — which team is the prop's subject.
#' @param legs list of leg specs as returned by parse_legs().
compute_prop_fair <- function(samples, side, legs) {
  if (is.null(legs) || length(legs) == 0) return(NA_real_)
  if (nrow(samples) == 0) return(NA_real_)

  # Derive team-specific run totals once from margin + total.
  home_runs <- (samples$total_final_score + samples$home_margin) / 2
  away_runs <- (samples$total_final_score - samples$home_margin) / 2
  team_runs <- if (side == "home") home_runs else away_runs
  opp_runs  <- if (side == "home") away_runs else home_runs

  hits <- rep(TRUE, nrow(samples))
  for (leg in legs) {
    hits <- hits & eval_leg(leg, samples, side, team_runs, opp_runs)
  }
  # NA filter is per-leg now — eval_leg returns NA for rows where the
  # side-specific column is NA. Drop those after AND-reduce so legs that
  # don't depend on inning-1 scoring at all aren't affected by inning-1 NAs.
  hits <- hits[!is.na(hits)]
  if (length(hits) == 0) return(NA_real_)
  mean(hits)
}

#' Blend a model fair probability with a DK SGP devigged probability.
#'
#' If `dk_decimal` is NA / 0 / negative, treats DK as missing and returns
#' `model_prob` alone.
#' If `model_prob` is NA but DK is valid, returns DK's devigged probability
#' alone.
#' If both are missing, returns NA_real_.
#' Otherwise returns the mean of the two devigged probabilities — the
#' blending strategy used by mlb_correlated_parlay.R for DK + FD SGPs.
#'
#' @param model_prob numeric, the historical fair probability from
#'   compute_prop_fair (may be NA for unpriceable props).
#' @param dk_decimal numeric, DK SGP decimal odds (>1.0 if valid; NA / 0
#'   / negative when DK couldn't price the prop).
#' @param vig numeric, the SGP vig fallback (default convention is 1.10,
#'   matching DK_SGP_VIG_DEFAULT in mlb_correlated_parlay.R).
#' @return numeric scalar — the blended probability — or NA_real_ if no
#'   inputs are available.
blend_dk_with_model <- function(model_prob, dk_decimal, vig) {
  dk_valid <- !is.na(dk_decimal) && dk_decimal > 0
  dk_prob  <- if (dk_valid) (1 / dk_decimal) / vig else NA_real_
  probs    <- c(model_prob, dk_prob)
  probs    <- probs[!is.na(probs)]
  if (length(probs) == 0) return(NA_real_)
  mean(probs)
}
