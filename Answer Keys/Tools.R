# =============================================================================
# TOOLS.R - Answer Key Utility Functions
# =============================================================================
#
# PERFORMANCE NOTES:
# -----------------
# The main bottleneck is run_answer_key_sample() which runs mean_match +
# balance_sample (~5-6 sec/game with 21k historical games).
#
# Rcpp acceleration: mean_match_cpp and balance_sample_cpp in src/sampling.cpp
# provide 10-20x speedup via nth_element and incremental bookkeeping.
# Falls back to pure R if Rcpp is unavailable.
#
# =============================================================================

# --- Rcpp acceleration (optional, falls back to R if unavailable) ---
.use_rcpp <- tryCatch({
  suppressMessages(library(Rcpp))
  sourceCpp("src/sampling.cpp")
  TRUE
}, error = function(e) {
  message("Rcpp acceleration unavailable, using pure R: ", e$message)
  FALSE
})

# --- Pipeline timing utility ---
pipeline_timer <- function() {
  t0 <- Sys.time()
  list(
    mark = function(label) {
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      cat(sprintf("  [%5.1fs] %s\n", elapsed, label))
      t0 <<- Sys.time()
      invisible(elapsed)
    }
  )
}

odds_to_prob <- function(odds) {
  ifelse(odds > 0, 100 / (odds + 100), -odds / (-odds + 100))
}

devig_american <- function(odd1, odd2) {
  p1_raw <- ifelse(odd1 > 0,
                   100 / (odd1 + 100),
                   -odd1 / (-odd1 + 100)
  )
  p2_raw <- ifelse(odd2 > 0,
                   100 / (odd2 + 100),
                   -odd2 / (-odd2 + 100)
  )
  total_raw <- p1_raw + p2_raw
  data.frame(p1 = p1_raw / total_raw, p2 = p2_raw / total_raw)
}

#' Devig 3-way American odds (home/away/tie)
#' @param odd_home Home win odds (American format)
#' @param odd_away Away win odds (American format)
#' @param odd_tie Tie/Draw odds (American format)
#' @return data.frame with p_home, p_away, p_tie (devigged probabilities summing to 1)
devig_american_3way <- function(odd_home, odd_away, odd_tie) {
  p_home_raw <- ifelse(odd_home > 0, 100 / (odd_home + 100), -odd_home / (-odd_home + 100))
  p_away_raw <- ifelse(odd_away > 0, 100 / (odd_away + 100), -odd_away / (-odd_away + 100))
  p_tie_raw <- ifelse(odd_tie > 0, 100 / (odd_tie + 100), -odd_tie / (-odd_tie + 100))

  total_raw <- p_home_raw + p_away_raw + p_tie_raw
  data.frame(
    p_home = p_home_raw / total_raw,
    p_away = p_away_raw / total_raw,
    p_tie = p_tie_raw / total_raw
  )
}
american_prob <- function(odd1, odd2) {
  data.frame(
    p1 = ifelse(odd1 > 0,
                100 / (odd1 + 100),
                -odd1 / (-odd1 + 100)
    ),
    p2 = ifelse(odd2 > 0,
                100 / (odd2 + 100),
                -odd2 / (-odd2 + 100)
    )
  )
}
# Sharp books for live consensus — weight > 0 included, 0 excluded.
# Pinnacle + Bookmaker at 1.1 = tiebreaker preference over 1.0-weight sharps.
SHARP_BOOKS <- list(
  pinnacle    = 1.1,
  bookmaker   = 1.1,
  lowvig      = 1.0,
  circasports = 1.0,
  bet105      = 1.0
)

logit_ <- function(p) log(p / (1 - p))
invlogit_ <- function(x) 1 / (1 + exp(-x))
logloss_ <- function(p, y) -(y * log(p) + (1 - y) * log(1 - p))

pick_consensus_line <- function(df,
                                game_id_col = "game_pk",
                                line_col = "total_line",
                                weight_col = "weight",
                                date_col = "game_date",
                                market1 = "devig_over_odds",
                                market2 = "devig_under_odds",
                                time_col = "game_start_time",
                                home = "home_team",
                                away = "away_team") {
  # Purpose: For each game, pick the line (spread or total) that has the
  # highest total weight across all sportsbooks posting it, then calculate
  # weighted average probability across all books posting that line.
  #
  # Args:
  #   df: data frame containing at least (game_id_col, date_col, line_col, weight_col)
  #   game_id_col: column name for unique game ID (string)
  #   line_col: column name for the line (spread or total)
  #   weight_col: column name for book weights
  #   date_col: column name for game date
  #
  # Returns:
  #   data frame with one row per game_id containing the consensus line

  df %>%
    # Step 1: Calculate total weight for each line option
    group_by(.data[[game_id_col]], .data[[date_col]], .data[[line_col]]) %>%
    mutate(total_weight = sum(.data[[weight_col]], na.rm = TRUE)) %>%
    ungroup() %>%
    # Step 2: Keep only rows where line has max weight for that game
    group_by(.data[[game_id_col]], .data[[date_col]]) %>%
    filter(total_weight == max(total_weight, na.rm = TRUE)) %>%
    ungroup() %>%
    # Step 3: Calculate weighted probabilities (keeping ALL books posting winning line)
    mutate(market1_weighted_prob = .data[[market1]] * .data[[weight_col]],
           market2_weighted_prob = .data[[market2]] * .data[[weight_col]]) %>%
    # Step 4: Aggregate to get weighted average across all books
    group_by(.data[[game_id_col]], .data[[date_col]], .data[[line_col]], .data[[time_col]], .data[[home]], .data[[away]]) %>%
    summarize(consensus_market1 = sum(market1_weighted_prob, na.rm = TRUE) / sum(.data[[weight_col]], na.rm = TRUE),
              consensus_market2 = sum(market2_weighted_prob, na.rm = TRUE) / sum(.data[[weight_col]], na.rm = TRUE),
              .groups = "drop") %>%
    # Step 5: Handle ties (if two lines had same weight, pick one)
    group_by(.data[[game_id_col]], .data[[date_col]]) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    rename(!!paste0("consensus_", market1) := consensus_market1,
           !!paste0("consensus_", market2) := consensus_market2)
}

moneyline_consensus <- function(df, 
                                game_id_col = "game_pk", 
                                weight_col = "weight", 
                                date_col = "game_date",
                                market1 = "devig_home_odds",
                                market2 = "devig_away_odds",
                                time_col = "game_start_time",
                                home = "home_team",
                                away = "away_team"){
  df %>% 
    group_by(.data[[game_id_col]], .data[[date_col]], .data[[time_col]], .data[[home]],.data[[away]]) %>%
    mutate(market1_weighted_prob = .data[[market1]] * .data[[weight_col]],
           market2_weighted_prob = .data[[market2]] * .data[[weight_col]]) %>% 
    summarize(consensus_market1 = sum(market1_weighted_prob, na.rm = T) / sum(.data[[weight_col]], na.rm = T),
              consensus_market2 = sum(market2_weighted_prob, na.rm = T) / sum(.data[[weight_col]], na.rm = T)) %>% 
    rename(!!paste0("consensus_",market1,sep = "") := consensus_market1,
           !!paste0("consensus_",market2,sep = "") := consensus_market2)
  
}

scraper_to_odds_api_format <- function(offshore_odds, game_odds) {
  # Convert scraper output (wide format) to Odds API long format so it can

  # be fed into prepare_two_way_odds() -> pick_consensus_line() for consensus.
  # Joins to game_odds on team names to get `id` and `commence_time`.
  # Returns data frame with Odds API columns, or empty data frame on failure.

  if (nrow(offshore_odds) == 0 || nrow(game_odds) == 0) return(data.frame())

  # Build lookup: canonical home/away -> game id + commence_time (one row per game)
  game_lookup <- game_odds %>%
    group_by(id, home_team, away_team) %>%
    summarize(commence_time = first(commence_time), .groups = "drop")

  # Join scraper rows to canonical games (full-game main lines only)
  matched <- offshore_odds %>%
    filter(market %in% c("spreads", "totals")) %>%
    inner_join(game_lookup, by = c("home_team", "away_team"))

  if (nrow(matched) == 0) return(data.frame())

  # Build spread rows (long format: one row per outcome)
  spread_rows <- matched %>%
    filter(market_type == "spreads", !is.na(line)) %>%
    rowwise() %>%
    reframe(
      id = id, commence_time = commence_time,
      home_team = home_team, away_team = away_team,
      bookmaker_key = bookmaker_key, market_key = "spreads",
      outcomes_name  = c(home_team, away_team),
      outcomes_price = c(odds_home, odds_away),
      outcomes_point = c(line, -line)
    )

  # Build total rows
  total_rows <- matched %>%
    filter(market_type == "totals", !is.na(line)) %>%
    rowwise() %>%
    reframe(
      id = id, commence_time = commence_time,
      home_team = home_team, away_team = away_team,
      bookmaker_key = bookmaker_key, market_key = "totals",
      outcomes_name  = c("Over", "Under"),
      outcomes_price = c(odds_over, odds_under),
      outcomes_point = c(line, line)
    )

  result <- bind_rows(spread_rows, total_rows)

  if (nrow(result) > 0) {
    bk <- unique(result$bookmaker_key)[1]
    cat(sprintf("Added %d %s rows to consensus pool (spreads=%d, totals=%d)\n",
                nrow(result), bk,
                sum(result$market_key == "spreads"),
                sum(result$market_key == "totals")))
  }

  result
}

get_event_odds_by_id <- function(event_id, commence_time,
                                 url = "https://api.the-odds-api.com/v4/historical/sports/baseball_mlb/odds") {
  snapshot <- format(commence_time - minutes(15), "%Y-%m-%dT%H:%M:%SZ")
  
  res <- GET(
    url = url,
    query = list(
      apiKey = Sys.getenv("ODDS_API_KEY"),
      date = snapshot,
      eventIds = event_id,
      regions = "us,eu,us2",
      markets = "h2h,totals,spreads",
      oddsFormat = "american",
      dateFormat = "iso"
    )
  )
  
  if (http_status(res)$category != "Success") {
    warning(paste("Request failed for", event_id, "at snapshot:", snapshot))
    return(tibble())
  }
  
  parsed <- fromJSON(content(res, as = "text"), flatten = TRUE)
  
  if (length(parsed$data) == 0) {
    return(tibble())
  }
  as_tibble(parsed$data)
}

distance_index <- function(dt, ps, pt, ss, st) {
  dt[, index := ((home_ml_odds - ps) / ss)^2 + ((total_line - pt) / st)^2]
}

# Function 1: mean‐matching to select initial N games
# Now supports EITHER spread lines OR moneyline probabilities
mean_match <- function(dt, N, parent_spread, parent_total,
                       ss, st, max_iter_mean, tol_mean,
                       use_spread_line = FALSE) {
  dt <- copy(dt)

  # Defensive: ensure N is valid integer
  N <- as.integer(N)
  if (is.na(N) || N <= 0) {
    stop(paste("Invalid N in mean_match. N =", N))
  }
  if (N > nrow(dt)) {
    warning(paste("N (", N, ") > nrow(dt) (", nrow(dt), "). Using nrow(dt) instead."))
    N <- nrow(dt)
  }

  # Validate parent_spread and parent_total are single numerics
  if (!is.numeric(parent_spread) || length(parent_spread) != 1 || is.na(parent_spread)) {
    stop(paste("Invalid parent_spread in mean_match. Value:", parent_spread,
               "Type:", class(parent_spread), "Length:", length(parent_spread)))
  }
  if (!is.numeric(parent_total) || length(parent_total) != 1 || is.na(parent_total)) {
    stop(paste("Invalid parent_total in mean_match. Value:", parent_total,
               "Type:", class(parent_total), "Length:", length(parent_total)))
  }

  # Determine which column to match against
  spread_col <- if (use_spread_line) "home_spread" else "home_ml_odds"

  # Check column exists
  if (!spread_col %in% names(dt)) {
    stop(paste("Column", spread_col, "not found in dt. Available columns:",
               paste(head(names(dt), 20), collapse = ", ")))
  }

  # Check total_line exists
  if (!"total_line" %in% names(dt)) {
    stop("Column 'total_line' not found in dt")
  }

  # --- C++ fast path ---
  if (.use_rcpp) {
    result <- mean_match_cpp(
      spread_vals    = dt[[spread_col]],
      total_vals     = dt[["total_line"]],
      game_date_num  = as.numeric(dt[["game_date"]]),
      N              = N,
      parent_spread  = parent_spread,
      parent_total   = parent_total,
      ss             = ss,
      st             = st,
      max_iter_mean  = as.integer(max_iter_mean),
      tol_mean       = tol_mean
    )

    # Reorder dt to match C++ distance-sorted output
    dt <- dt[result$order, ]
    set(dt, j = "included", value = FALSE)
    set(dt, i = 1L:result$N, j = "included", value = TRUE)
    # Recompute distance index on the reordered table (for downstream use)
    distance_index_generic(dt, parent_spread, parent_total, ss, st, spread_col)

    return(list(
      dt = dt,
      parent_spread = parent_spread,
      parent_total = parent_total,
      converged = result$converged
    ))
  }

  # --- Pure R fallback ---
  adj_spread <- parent_spread
  adj_total <- parent_total
  mm_converged <- FALSE

  for (iter in seq_len(max_iter_mean)) {
    distance_index_generic(dt, adj_spread, adj_total, ss, st, spread_col)
    setorder(dt, index, -game_date)

    set(dt, j = "included", value = FALSE)
    set(dt, i = 1L:N, j = "included", value = TRUE)

    mean_s <- mean(dt[included == TRUE][[spread_col]], na.rm = TRUE)
    mean_t <- mean(dt[included == TRUE][["total_line"]], na.rm = TRUE)
    err_s <- mean_s - parent_spread
    err_t <- mean_t - parent_total

    if (abs(err_s) < tol_mean && abs(err_t) < tol_mean) {
      mm_converged <- TRUE
      break
    }

    adj_spread <- adj_spread - err_s
    adj_total <- adj_total - err_t
  }
  list(
    dt = dt,
    parent_spread = parent_spread,
    parent_total = parent_total,
    converged = mm_converged
  )
}

# New generic distance calculator
distance_index_generic <- function(dt, ps, pt, ss, st, spread_col = "home_ml_odds") {
  # Validate inputs
  if (!spread_col %in% names(dt)) {
    stop(paste("Column", spread_col, "not found in distance_index_generic"))
  }
  if (!is.numeric(ps) || length(ps) != 1) {
    stop(paste("ps must be single numeric. Value:", ps, "Type:", class(ps), "Length:", length(ps)))
  }
  if (!is.numeric(pt) || length(pt) != 1) {
    stop(paste("pt must be single numeric. Value:", pt, "Type:", class(pt), "Length:", length(pt)))
  }
  if (!is.numeric(ss) || length(ss) != 1 || ss == 0) {
    stop(paste("ss must be single non-zero numeric. Value:", ss, "Type:", class(ss)))
  }
  if (!is.numeric(st) || length(st) != 1 || st == 0) {
    stop(paste("st must be single non-zero numeric. Value:", st, "Type:", class(st)))
  }
  
  # Use set() for direct column assignment to avoid data.table scoping issues
  index_vals <- ((dt[[spread_col]] - ps) / ss)^2 + ((dt[["total_line"]] - pt) / st)^2
  set(dt, j = "index", value = index_vals)
  invisible(dt)
}

# Function 2: balance sample by greedy remove/add
# Per Feustel spec: if both add and remove fail, return converged=FALSE
# and let caller restart with smaller N
balance_sample <- function(dt, N, target_cover, target_over, tol_error) {

  # --- C++ fast path ---
  if (.use_rcpp) {
    result <- balance_sample_cpp(
      order_indices = seq_len(nrow(dt)),  # dt is already sorted by distance from mean_match
      actual_cover  = as.integer(dt[["actual_cover"]]),
      actual_over   = as.integer(dt[["actual_over"]]),
      N             = as.integer(N),
      target_cover  = target_cover,
      target_over   = target_over,
      tol_error     = as.integer(tol_error)
    )

    dt <- copy(dt)
    set(dt, j = "included", value = FALSE)
    set(dt, i = result$included_indices, j = "included", value = TRUE)

    return(list(
      dt = dt,
      final_N = length(result$included_indices),
      cover_error = result$cover_error,
      over_error = result$over_error,
      converged = result$converged
    ))
  }

  # --- Pure R fallback ---
  dt <- copy(dt)
  n_sample <- N

  cover_error <- dt[included == TRUE, sum(actual_cover, na.rm = T)] -
    round(target_cover * n_sample)
  over_error <- dt[included == TRUE, sum(actual_over, na.rm = T)] -
    round(target_over * n_sample)
  M <- nrow(dt)

  repeat {
    if ((abs(cover_error) <= tol_error) && (abs(over_error) <= tol_error)) {
      return(list(
        dt = dt,
        final_N = n_sample,
        cover_error = cover_error,
        over_error = over_error,
        converged = TRUE
      ))
    }

    removal_failed <- TRUE
    addition_failed <- TRUE

    # ---- remove one game from 1..n_sample ----
    in1toN <- dt[seq_len(n_sample), included]
    cov_help1 <- (cover_error > 0 & dt[seq_len(n_sample), actual_cover] == 1) |
      (cover_error < 0 & dt[seq_len(n_sample), actual_cover] == 0)
    ov_help1 <- (over_error > 0 & dt[seq_len(n_sample), actual_over] == 1) |
      (over_error < 0 & dt[seq_len(n_sample), actual_over] == 0)

    both1 <- which(in1toN & cov_help1 & ov_help1)
    either1 <- which(in1toN & (cov_help1 | ov_help1))

    if (length(both1) > 0) {
      i <- max(both1)
      set(dt, i = as.integer(i), j = "included", value = FALSE)
      removal_failed <- FALSE
    } else if (length(either1) > 0) {
      i <- max(either1)
      set(dt, i = as.integer(i), j = "included", value = FALSE)
      removal_failed <- FALSE
    }

    if (!removal_failed) {
      cover_error <- dt[included == TRUE, sum(actual_cover, na.rm = T)] -
        round(target_cover * n_sample)
      over_error <- dt[included == TRUE, sum(actual_over, na.rm = T)] -
        round(target_over * n_sample)
    }

    # ---- add one game from (n_sample+1)..M ----
    idx_range <- seq(n_sample + 1, M)
    inNplus <- dt[idx_range, included]
    cov_help2 <- (cover_error > 0 & dt[idx_range, actual_cover] == 0) |
      (cover_error < 0 & dt[idx_range, actual_cover] == 1)
    ov_help2 <- (over_error > 0 & dt[idx_range, actual_over] == 0) |
      (over_error < 0 & dt[idx_range, actual_over] == 1)

    both2 <- which(!inNplus & cov_help2 & ov_help2)
    either2 <- which(!inNplus & (cov_help2 | ov_help2))

    if (length(both2) > 0) {
      j <- min(both2) + n_sample
      set(dt, i = as.integer(j), j = "included", value = TRUE)
      addition_failed <- FALSE
    } else if (length(either2) > 0) {
      j <- min(either2) + n_sample
      set(dt, i = as.integer(j), j = "included", value = TRUE)
      addition_failed <- FALSE
    }

    if (!addition_failed) {
      cover_error <- dt[included == TRUE, sum(actual_cover, na.rm = T)] -
        round(target_cover * n_sample)
      over_error <- dt[included == TRUE, sum(actual_over, na.rm = T)] -
        round(target_over * n_sample)
    }

    # Per Feustel: "restart the procedure looking for a smaller sample"
    if (removal_failed && addition_failed) {
      return(list(
        dt = dt,
        final_N = n_sample,
        cover_error = cover_error,
        over_error = over_error,
        converged = FALSE
      ))
    }
  }
}

# =============================================================================
# ANSWER KEY CORE FUNCTIONS (New Architecture)
# =============================================================================

#' Run Answer Key sampling algorithm once per game
#'
#' This is the core function that runs mean_match + balance_sample to select
#' a balanced historical sample. The returned sample can then be used to
#' generate predictions for multiple market types (moneylines, spreads, totals).
#'
#' @param id Game identifier
#' @param parent_spread Target spread (or moneyline probability if use_spread_line=FALSE)
#' @param parent_total Target total line
#' @param target_cover Target cover rate (0.0-1.0)
#' @param target_over Target over rate (0.0-1.0)
#' @param DT Historical data table with game outcomes
#' @param ss Spread/probability standard deviation for distance weighting
#' @param st Total line standard deviation for distance weighting
#' @param N Target sample size
#' @param use_spread_line If TRUE, match on spread line; if FALSE, match on moneyline probability
#' @return List containing: sample (data frame), final_N, cover_error, over_error, converged, metadata
run_answer_key_sample <- function(
    id,
    parent_spread,
    parent_total,
    target_cover,
    target_over,
    DT, ss, st, N,
    max_iter_mean = 500,
    tol_mean = 0.005,
    tol_mean_gate = 0.1,
    tol_error = 1,
    use_spread_line = FALSE,
    shrink_factor = 0.9,
    min_N = 50
) {
  # Defensive: ensure N is valid
  N <- as.integer(N)
  if (is.na(N) || N <= 0) {
    stop(paste("Invalid N in run_answer_key_sample. N =", N))
  }

  spread_col <- if (use_spread_line) "home_spread" else "home_ml_odds"
  current_N <- N

  repeat {
    # Step 1: mean_match to get sample with correct spread/total means
    mm <- mean_match(DT, current_N, parent_spread, parent_total, ss, st,
                     max_iter_mean, tol_mean, use_spread_line)

    # Per Feustel: if mean is meaningfully off, restart with smaller sample.
    # Use tol_mean_gate (looser) for the shrink decision — tol_mean (tight) is
    # only the iteration target inside mean_match.
    inc <- mm$dt[mm$dt$included == TRUE, ]
    mean_err_s <- abs(mean(inc[[spread_col]], na.rm = TRUE) - parent_spread)
    mean_err_t <- abs(mean(inc[["total_line"]], na.rm = TRUE) - parent_total)
    mm_acceptable <- mean_err_s < tol_mean_gate && mean_err_t < tol_mean_gate

    if (!mm_acceptable) {
      current_N <- as.integer(current_N * shrink_factor)
      if (current_N < min_N) {
        warning(paste("mean_match could not converge at min_N =", min_N,
                      "for spread =", parent_spread, "total =", parent_total,
                      "(err_s =", round(mean_err_s, 3), "err_t =", round(mean_err_t, 3), ")"))
        return(NULL)
      }
      next
    }

    # Step 2: balance_sample to match cover/over rates
    bal <- balance_sample(mm$dt, current_N, target_cover, target_over, tol_error)

    if (bal$converged) {
      # Success - use this sample
      break
    }

    # Per Feustel: "restart the procedure looking for a smaller sample"
    current_N <- as.integer(current_N * shrink_factor)

    if (current_N < min_N) {
      # Cannot shrink further - use current sample (accept as-is per spec)
      warning(paste("balance_sample could not converge. Using sample with N =",
                    bal$final_N, "cover_error =", bal$cover_error,
                    "over_error =", bal$over_error))
      break
    }
  }

  # Return the balanced sample with metadata
  list(
    sample = bal$dt[bal$dt$included == TRUE, ],
    final_N = sum(bal$dt$included),
    cover_error = bal$cover_error,
    over_error = bal$over_error,
    converged = bal$converged,
    metadata = list(
      id = id,
      parent_spread = parent_spread,
      parent_total = parent_total,
      target_cover = target_cover,
      target_over = target_over
    )
  )
}

#' Generate moneyline predictions from a balanced sample
#'
#' @param sample Data frame of selected games from run_answer_key_sample()
#' @param margin_col Column name prefix for margin columns (e.g., "game_home_margin_period")
#' @return Tibble with 2-way and 3-way probabilities for each period
predict_moneyline_from_sample <- function(
    sample,
    margin_col = "game_home_margin_period"
) {
  n_games <- nrow(sample)

  # 2-way probabilities (excludes ties): home_wins / (home_wins + away_wins)
  probs_2way <- sample %>%
    summarise(across(starts_with(margin_col),
                     ~ sum(.x > 0, na.rm = TRUE) / sum(.x != 0, na.rm = TRUE))) %>%
    ungroup()

  # 3-way probabilities (includes ties): each outcome as fraction of total
  probs_3way_home <- sample %>%
    summarise(across(starts_with(margin_col), ~ sum(.x > 0, na.rm = TRUE) / n_games)) %>%
    rename_with(~ paste0(.x, "_3way_home"), starts_with(margin_col))

  probs_3way_away <- sample %>%
    summarise(across(starts_with(margin_col), ~ sum(.x < 0, na.rm = TRUE) / n_games)) %>%
    rename_with(~ paste0(.x, "_3way_away"), starts_with(margin_col))

  probs_3way_tie <- sample %>%
    summarise(across(starts_with(margin_col), ~ sum(.x == 0, na.rm = TRUE) / n_games)) %>%
    rename_with(~ paste0(.x, "_3way_tie"), starts_with(margin_col))

  bind_cols(probs_2way, probs_3way_home, probs_3way_away, probs_3way_tie)
}

#' Generate spread cover predictions from a balanced sample
#'
#' @param sample Data frame of selected games from run_answer_key_sample()
#' @param spreads Vector of spread values to calculate cover probability for
#' @param margin_col Column name prefix for margin columns
#' @param period Period suffix (e.g., "1" for Q1, "Half1" for first half)
#' @return Tibble with cover probability for each spread value
predict_spreads_from_sample <- function(
    sample,
    spreads,
    margin_col = "game_home_margin_period",
    period
) {
  col_name <- paste0(margin_col, "_", period)

  if (!(col_name %in% names(sample))) {
    warning(paste("Column", col_name, "not found in sample"))
    return(tibble())
  }

  margins <- sample[[col_name]]

  # For each spread, calculate probability of home covering
  # Home covers when margin > -spread (e.g., if spread is -7, home covers when margin > 7)
  pct_cols <- set_names(
    map(spreads, ~ sum(margins > -.x, na.rm = TRUE) / sum(margins != -.x, na.rm = TRUE)),
    paste0("pct_home_cover_", gsub("-", "neg", gsub("\\.", "_", as.character(spreads))))
  )

  tibble(!!!pct_cols)
}

#' Generate total over predictions from a balanced sample
#'
#' @param sample Data frame of selected games from run_answer_key_sample()
#' @param totals Vector of total values to calculate over probability for
#' @param total_col Column name prefix for total columns
#' @param period Period suffix (e.g., "1" for Q1, "Half1" for first half)
#' @return Tibble with over probability for each total value
predict_totals_from_sample <- function(
    sample,
    totals,
    total_col = "game_total_period",
    period
) {
  col_name <- paste0(total_col, "_", period)

  if (!(col_name %in% names(sample))) {
    warning(paste("Column", col_name, "not found in sample"))
    return(tibble())
  }

  total_vals <- sample[[col_name]]

  # For each total, calculate probability of over
  pct_cols <- set_names(
    map(totals, ~ sum(total_vals > .x, na.rm = TRUE) / sum(total_vals != .x, na.rm = TRUE)),
    paste0("pct_over_", gsub("\\.", "_", as.character(totals)))
  )

  tibble(!!!pct_cols)
}

#' Generate samples for all games upfront
#'
#' This function runs run_answer_key_sample() once per game and returns
#' a named list of samples that can be reused across multiple market types.
#'
#' @param targets Data frame with id, parent_spread, parent_total, target_cover, target_over
#' @param DT Historical data table
#' @param ss Spread/probability standard deviation
#' @param st Total line standard deviation
#' @param N Target sample size
#' @param use_spread_line If TRUE, match on spread line; if FALSE, match on moneyline probability
#' @return Named list of sample results, keyed by game id
generate_all_samples <- function(
    targets,
    DT,
    ss, st, N,
    use_spread_line = TRUE,
    max_iter_mean = 500,
    tol_mean = 0.005,
    tol_error = 1,
    shrink_factor = 0.9,
    min_N = 50
) {
  n_games <- nrow(targets)
  n_cores <- min(parallel::detectCores() - 1L, n_games, na.rm = TRUE)
  n_cores <- max(n_cores, 1L)
  cat(sprintf("Generating samples for %d games using %d cores...\n", n_games, n_cores))

  N <- as.integer(N)

  # Split targets into a list of single-row lists for mclapply
  target_list <- lapply(seq_len(n_games), function(i) {
    list(
      id = targets$id[i],
      parent_spread = targets$parent_spread[i],
      parent_total = targets$parent_total[i],
      target_cover = targets$target_cover[i],
      target_over = targets$target_over[i]
    )
  })

  samples <- parallel::mclapply(target_list, function(tgt) {
    run_answer_key_sample(
      id = tgt$id,
      parent_spread = tgt$parent_spread,
      parent_total = tgt$parent_total,
      target_cover = tgt$target_cover,
      target_over = tgt$target_over,
      DT = DT, ss = ss, st = st, N = N,
      max_iter_mean = max_iter_mean,
      tol_mean = tol_mean,
      tol_error = tol_error,
      use_spread_line = use_spread_line,
      shrink_factor = shrink_factor,
      min_N = min_N
    )
  }, mc.cores = n_cores)

  # Name the list by game id for easy lookup
  names(samples) <- targets$id
  # Drop games where mean_match could not converge (returned NULL)
  samples <- samples[!vapply(samples, is.null, logical(1))]

  cat(sprintf("Generated %d samples (%d skipped)\n",
              length(samples), nrow(targets) - length(samples)))
  samples
}

# =============================================================================
# END ANSWER KEY CORE FUNCTIONS
# =============================================================================

fetch_event_odds <- function(event_id, market, sport_key = "baseball_mlb") {
  res <- GET(
    paste0("https://api.the-odds-api.com/v4/sports/", sport_key, "/events/", event_id, "/odds"),
    query = list(
      apiKey     = Sys.getenv("ODDS_API_KEY"),
      regions    = "us,us2,us_ex",
      markets    = market,
      oddsFormat = "american",
      dateFormat = "iso"
    )
  )
  
  if (http_error(res)) {
    warning(paste("Failed for event:", event_id))
    return(tibble())
  }
  
  odds_data <- fromJSON(content(res, "text"), flatten = TRUE)
  if (is.null(odds_data) || length(odds_data) == 0) {
    return(tibble())
  }
  as_tibble(odds_data)
}

fetch_odds_bulk <- function(event_ids, markets, sport_key = "basketball_ncaab",
                            batch_size = 5) {
  api_key <- Sys.getenv("ODDS_API_KEY")
  markets_param <- paste(markets, collapse = ",")
  n <- length(event_ids)
  # Use environment for reference semantics (callbacks can write from any scope)
  res_env <- new.env(parent = emptyenv())

  # Process in batches to avoid overwhelming the API with concurrent connections
  for (batch_start in seq(1, n, by = batch_size)) {
    batch_end <- min(batch_start + batch_size - 1, n)
    pool <- curl::new_pool()

    for (i in batch_start:batch_end) {
      local({
        idx <- i
        eid <- event_ids[i]
        url <- paste0(
          "https://api.the-odds-api.com/v4/sports/", sport_key, "/events/",
          eid, "/odds?",
          "apiKey=", utils::URLencode(api_key, reserved = TRUE),
          "&regions=us,us2,us_ex",
          "&markets=", markets_param,
          "&oddsFormat=american",
          "&dateFormat=iso"
        )
        curl::curl_fetch_multi(url, done = function(resp) {
          res_env[[as.character(idx)]] <- data.frame(
            event_id = eid,
            json_response = if (resp$status_code == 200) rawToChar(resp$content) else NA_character_,
            stringsAsFactors = FALSE
          )
        }, fail = function(msg) {
          res_env[[as.character(idx)]] <- data.frame(
            event_id = eid,
            json_response = NA_character_,
            stringsAsFactors = FALSE
          )
        }, pool = pool)
      })
    }

    curl::multi_run(pool = pool)
  }

  do.call(rbind, as.list(res_env))
}

prob_to_american <- function(prob) {
  ifelse(prob >= 0.5,
         # For favorites
         -(prob / (1 - prob)) * 100,
         # For underdogs
         ((1 - prob) / prob) * 100
  ) %>% round(0) # round to nearest integer
}

#understand spread of line
compute_dispersion <- function(
    DT,
    odds_cols = c("home_ml_odds", "away_ml_odds"),
    total_col = "total_line",
    spread_col = "home_spread",
    lower = 0.05,
    upper = 0.95,
    moneyline = TRUE
) {
  # combine all price/probability columns into one vector
  if(moneyline){
    odds_values <- DT %>%
      select(all_of(odds_cols)) %>%
      pivot_longer(everything()) %>%
      pull(value)
    qs <- quantile(odds_values, probs = c(lower, upper), na.rm = TRUE)
  }
  else{
    qs <- quantile(DT[[spread_col]], probs = c(lower, upper), na.rm = TRUE)
  }
  
  qt <- quantile(DT[[total_col]], probs = c(lower, upper), na.rm = TRUE)
  
  list(
    ss = as.numeric(qs[[2]] - qs[[1]]),
    st = as.numeric(qt[[2]] - qt[[1]])
  )
}

prepare_two_way_odds <- function(
    game_odds,
    mkt_key,                    # "h2h", "totals", "spreads", etc
    book_weights,
    prob_fun,                   # american_prob or devig_american
    prob_names = c("prob_1", "prob_2"),
    odds_names
) { #prepare odds for consensus
  game_odds %>%
    filter(market_key == mkt_key) %>%
    mutate(
      date      = as.Date(commence_time),
      odds_type = ifelse(
        home_team == outcomes_name, "home",
        ifelse(away_team == outcomes_name, "away", outcomes_name)
      )
    ) %>%
    pivot_wider(
      id_cols      = c(id, commence_time, home_team, away_team, bookmaker_key, date),
      names_from   = odds_type,
      values_from  = c(outcomes_price,outcomes_point)
    ) %>% 
    mutate(as_tibble(prob_fun(.data[[odds_names[1]]], .data[[odds_names[2]]]))) %>%
    rename(
      !!prob_names[1] := p1,
      !!prob_names[2] := p2
    ) %>%
    left_join(book_weights, by = "bookmaker_key")
}

flatten_event_odds <- function(raw_odds) {
  raw_odds %>%
    # 1. one row per bookmaker
    unnest_longer(bookmakers) %>%
    unnest_wider(bookmakers, names_sep = "_") %>%
    # 2. one row per market (handles multi-market JSON responses)
    mutate(bookmakers_markets = map(bookmakers_markets, ~ {
      if (is.list(.x) && !is.data.frame(.x)) bind_rows(.x) else .x
    })) %>%
    unnest(bookmakers_markets, names_sep = "_") %>%
    # 3. extract outcome columns as list columns (one row per bookmaker-market)
    mutate(
      bookmakers_markets_outcomes_name = map(bookmakers_markets_outcomes, "name"),
      bookmakers_markets_outcomes_price = map(bookmakers_markets_outcomes, "price"),
      bookmakers_markets_outcomes_point = map(bookmakers_markets_outcomes, ~ {
        if ("point" %in% names(.x)) .x$point else NULL
      })
    ) %>%
    select(-bookmakers_markets_outcomes) %>%
    # 4. parse datetimes + standardize column names
    mutate(
      commence_time = ymd_hms(commence_time, tz = "UTC"),
      market_update = ymd_hms(bookmakers_markets_last_update, tz = "UTC")
    ) %>%
    rename(
      bookmaker_key   = bookmakers_key,
      bookmaker_title = bookmakers_title,
      market_key      = bookmakers_markets_key,
      outcome_name    = bookmakers_markets_outcomes_name,
      closing_odds    = bookmakers_markets_outcomes_price
    )
}

compute_ev <- function(pred_prob, book_prob) {
  pred_prob * ((1 / book_prob) - 1) - (1 - pred_prob)
}

#' Compute EV for 3-way market (home/away/tie)
#' For 3-way, EV = pred_prob * (decimal_odds - 1) - (1 - pred_prob)
#' This is the same formula but we calculate for each of the three outcomes independently
compute_ev_3way <- function(pred_prob, book_odds) {
  # Convert American to decimal

  decimal_odds <- ifelse(book_odds > 0, 1 + book_odds / 100, 1 + 100 / abs(book_odds))
  # EV calculation
  pred_prob * (decimal_odds - 1) - (1 - pred_prob)
}

kelly_stake <- function(ev, book_prob, bankroll, kelly_mult) {
  edge_fraction <- ev / ((1 / book_prob) - 1)
  stake <- edge_fraction * kelly_mult * bankroll
  round(stake, 2)
}

get_events <- function(
    sport_key,
    regions = "us",
    tz_local = "America/Los_Angeles"
) {
  res <- httr::GET(
    paste0("https://api.the-odds-api.com/v4/sports/", sport_key, "/events"),
    query = list(
      apiKey = Sys.getenv("ODDS_API_KEY"),
      regions = regions
    )
  )
  httr::stop_for_status(res)
  
  now_local   <- lubridate::with_tz(Sys.time(), tz_local)
  
  fromJSON(content(res, "text"), flatten = TRUE) %>%
    as_tibble() %>%
    mutate(
      pt_start_time = lubridate::with_tz(
        lubridate::ymd_hms(commence_time, tz = "UTC"),
        tz_local
      )
    ) %>%
    filter(
      pt_start_time > now_local
    )
}

#standardize bets table
format_bets_table <- function(
    df,
    pred1, pred2,          # names of predicted prob columns (e.g. "home_predicted_prob", "away_predicted_prob")
    ev1, ev2,              # names of EV columns           (e.g. "home_ev", "away_ev")
    size1, size2,          # names of bet size columns     (e.g. "home_bet_size", "away_bet_size")
    odds1, odds2,          # names of odds columns         (e.g. "book_home_market", "book_away_market")
    pred3 = NULL, ev3 = NULL, size3 = NULL, odds3 = NULL,  # optional 3rd side for 3-way markets
    cents1 = NULL, cents2 = NULL, cents3 = NULL,  # optional: effective cents (Kalshi precision)
    line_col = NULL,       # optional: single line column (e.g. "book_total_line" for totals)
    line_col_1 = NULL,     # optional: line for side 1 (e.g. "book_home_spread" for spreads)
    line_col_2 = NULL,     # optional: line for side 2 (e.g. "book_away_spread" for spreads)
    books   = NULL,
    time_col = "commence_time",
    tz_out   = "America/Los_Angeles",
    ev_threshold = 0.05    # minimum EV to include bet (default 5%)
) {
  if (is.null(books)) books <- unique(df$bookmaker_key)

  s1 <- str_extract(size1, "home|away|over|under")
  s2 <- str_extract(size2, "home|away|over|under")

  # Map columns to standard names: metric_side (e.g. prob_home, size_over)
  cols_map <- c(
    set_names(pred1, paste0("prob_", s1)), set_names(pred2, paste0("prob_", s2)),
    set_names(ev1,   paste0("ev_", s1)),   set_names(ev2,   paste0("ev_", s2)),
    set_names(size1, paste0("size_", s1)), set_names(size2, paste0("size_", s2)),
    set_names(odds1, paste0("odds_", s1)), set_names(odds2, paste0("odds_", s2))
  )

  # Add cents columns if provided (for Kalshi effective price precision)
  if (!is.null(cents1) && cents1 %in% names(df)) {
    cols_map <- c(cols_map,
      set_names(cents1, paste0("cents_", s1)),
      set_names(cents2, paste0("cents_", s2))
    )
  }

  # Add 3rd side if provided (for 3-way markets)
  if (!is.null(size3)) {
    s3 <- str_extract(size3, "home|away|over|under|tie")
    cols_map <- c(cols_map,
      set_names(pred3, paste0("prob_", s3)),
      set_names(ev3,   paste0("ev_", s3)),
      set_names(size3, paste0("size_", s3)),
      set_names(odds3, paste0("odds_", s3))
    )
    if (!is.null(cents3) && cents3 %in% names(df)) {
      cols_map <- c(cols_map, set_names(cents3, paste0("cents_", s3)))
    }
  }

  # Handle line columns - support both single line_col and separate line_col_1/line_col_2
  has_single_line <- !is.null(line_col) && line_col %in% names(df)
  has_separate_lines <- !is.null(line_col_1) && !is.null(line_col_2) &&
                        line_col_1 %in% names(df) && line_col_2 %in% names(df)

  result <- df %>%
    filter(bookmaker_key %in% books) %>%
    mutate(pt_start_time = lubridate::with_tz(as.POSIXct(.data[[time_col]], tz = "UTC"), tzone = tz_out))

  if (has_separate_lines) {
    # For spreads: separate line columns for each side (included in pivot)
    cols_map <- c(cols_map,
      set_names(line_col_1, paste0("line_", s1)),
      set_names(line_col_2, paste0("line_", s2))
    )
    result <- result %>%
      select(id, home_team, away_team, pt_start_time, bookmaker_key, market, all_of(cols_map))
  } else if (has_single_line) {
    # For totals/moneylines: single line column (excluded from pivot)
    result <- result %>%
      mutate(line = .data[[line_col]]) %>%
      select(id, home_team, away_team, pt_start_time, bookmaker_key, market, line, all_of(cols_map))
  } else {
    result <- result %>%
      mutate(line = NA_real_) %>%
      select(id, home_team, away_team, pt_start_time, bookmaker_key, market, line, all_of(cols_map))
  }

  # Pivot based on whether we have separate lines
  if (has_separate_lines) {
    result <- result %>%
      pivot_longer(
        cols = -c(id, home_team, away_team, pt_start_time, bookmaker_key, market),
        names_to = c(".value", "bet_on"),
        names_sep = "_"
      )
  } else {
    result <- result %>%
      pivot_longer(
        cols = -c(id, home_team, away_team, pt_start_time, bookmaker_key, market, line),
        names_to = c(".value", "bet_on"),
        names_sep = "_"
      )
  }

  result %>%
    # Map 'home'/'away' to team names; 'tie' to 'Tie'; keep 'Over'/'Under' as is
    mutate(bet_on = case_when(
      bet_on == "home" ~ home_team,
      bet_on == "away" ~ away_team,
      bet_on == "tie" ~ "Tie",
      TRUE ~ str_to_title(bet_on)
    )) %>%
    # Filter by EV threshold
    filter(ev >= ev_threshold) %>%
    arrange(desc(size)) %>%
    select(id, home_team, away_team, pt_start_time, bookmaker_key, market, bet_on, line, bet_size = size, ev, odds, prob, any_of("cents"))
}


# Entire Moneyline Process ----
build_moneyline_market <- function(
    DT,
    consensus_odds,
    ss, st, N,
    period,
    events,
    market,
    bankroll   = 200,
    kelly_mult = 0.25,
    targets = NULL,
    use_spread_line = FALSE,  # NEW: set TRUE for NFL spreads
    sport_key,
    margin_col
) {
  
  # If targets not provided, build them from consensus_odds
  if (is.null(targets)) {
    if (use_spread_line) {
      targets <- consensus_odds %>%
        transmute(
          id,
          parent_spread = spread,                    # actual spread line
          parent_total  = total_line,
          target_cover  = consensus_devig_home_odds, # probability of covering
          target_over   = consensus_devig_over_odds
        )
    } else {
      targets <- consensus_odds %>%
        transmute(
          id,
          parent_spread = consensus_devig_home_odds, # moneyline probability
          parent_total  = total_line,
          target_cover  = consensus_devig_home_odds,
          target_over   = consensus_devig_over_odds
        )
    }
  }
  # 2) Elihu engine for each id - using new architecture
  final_preds <- targets %>%
    mutate(res = pmap(
      list(id, parent_spread, parent_total, target_cover, target_over),
      function(id, ps, pt, tc, to) {
        # Step 1: Get balanced sample
        sample_result <- run_answer_key_sample(
          id = id, parent_spread = ps, parent_total = pt,
          target_cover = tc, target_over = to,
          DT = DT, ss = ss, st = st, N = N,
          use_spread_line = use_spread_line,
          max_iter_mean = 500
        )
        if (is.null(sample_result)) return(NULL)
        # Step 2: Generate moneyline predictions from sample
        predict_moneyline_from_sample(sample_result$sample, margin_col = margin_col)
      }
    )) %>%
    unnest(res) %>%
    ungroup() %>%
    inner_join(
      consensus_odds %>%
        ungroup() %>%
        select(id, home_team, away_team, commence_time),
      by = c("id")
    ) %>%
    select(home_team, away_team, commence_time, everything())

  predictions <- final_preds %>%
    mutate(
      across(
        starts_with(margin_col),
        ~ prob_to_american(.x),
        .names = "{.col}_american"
      )
    ) %>%
    relocate(
      ends_with("_american"),
      .before = starts_with(margin_col)
    )
  
  # 3) Get current odds for the chosen market
  all_odds <- events$id %>%
    map(~ fetch_event_odds(.x, market, sport_key)) %>%
    keep(~ nrow(.x) > 0) %>%
    bind_rows()
  
  flat_betting_odds <- all_odds %>%
    flatten_event_odds() %>%
    unnest_wider(outcome_name, names_sep = "_") %>%
    unnest_wider(closing_odds, names_sep = "_") %>%
    mutate(
      home_odds = ifelse(home_team == outcome_name_1, closing_odds_1, closing_odds_2),
      away_odds = ifelse(away_team == outcome_name_1, closing_odds_1, closing_odds_2)
    ) %>%
    group_by(
      id, commence_time, home_team, away_team,
      bookmaker_key, bookmaker_title
    ) %>%
    summarise(
      book_home_market = first(home_odds),
      book_away_market = first(away_odds),
      .groups = "drop"
    )
  
  # 4) Join predictions + consensus + market odds, compute EV + Kelly
  prediction_set <- flat_betting_odds %>%
    left_join(
      predictions %>% ungroup() %>% select(id, contains(paste0("_", period))),
      by = "id"
    ) %>%
    left_join(
      consensus_odds %>% ungroup() %>%  select(-home_team, -away_team,-commence_time),
      by = "id"
    ) %>%
    rename(
      book_full_game_home_prob  = consensus_prob_home,
      book_full_game_away_prob  = consensus_prob_away,
      book_full_game_over_prob  = consensus_prob_over,
      book_full_game_under_prob = consensus_prob_under,
      book_full_game_total_line = total_line
    ) %>%
    rename(
      home_predicted_prob          = paste0(margin_col, "_", period),
      home_predicted_american_odds = paste0(margin_col,"_", period, "_american")
    ) %>%
    mutate(
      away_predicted_prob = 1 - home_predicted_prob
    ) %>%
    mutate(as_tibble(american_prob(book_home_market, book_away_market))) %>%
    rename(
      book_market_prob_home = p1,
      book_market_prob_away = p2
    ) %>%
    filter(!is.na(spread)) %>% 
    mutate(
      home_ev       = compute_ev(home_predicted_prob, book_market_prob_home),
      away_ev       = compute_ev(away_predicted_prob, book_market_prob_away),
      home_bet_size = kelly_stake(home_ev, book_market_prob_home, bankroll, kelly_mult),
      away_bet_size = kelly_stake(away_ev, book_market_prob_away, bankroll, kelly_mult),
      market        = market
    )
  bets <- prediction_set %>%
    format_bets_table(
      pred1 = "home_predicted_prob",
      pred2 = "away_predicted_prob",
      ev1   = "home_ev",
      ev2   = "away_ev",
      size1 = "home_bet_size",
      size2 = "away_bet_size",
      odds1 = "book_home_market",
      odds2 = "book_away_market"
    )
  list(
    predictions    = predictions,
    prediction_set = prediction_set,
    bets           = bets
  )
}


#Entire Total Process ----

build_totals_market <- function(
    DT,
    consensus_odds,
    ss, st, N,
    period,             # inning / quarter / half etc
    events,
    market,             # e.g. "alternate_totals_1st_5_innings"
    sport_key,          # e.g. "baseball_mlb", "americanfootball_nfl"
    bankroll   = 200,
    kelly_mult = 0.25
) {
  # 1) Targets from consensus
  targets <- consensus_odds %>%
    transmute(
      id,
      parent_spread = consensus_prob_home,
      parent_total  = total_line,
      target_cover  = consensus_prob_home,
      target_over   = consensus_over
    )
  
  # 2) Current totals odds
  all_odds <- map_dfr(events$id, ~ fetch_event_odds(.x, market, sport_key))
  
  flat_betting_odds <- all_odds %>%
    flatten_event_odds() %>%
    unnest_wider(outcome_name, names_sep = "_") %>%
    unnest_wider(closing_odds, names_sep = "_") %>%
    unnest_wider(bookmakers_markets_outcomes_point, names_sep = "_") %>%
    group_by(
      id, commence_time, home_team, away_team,
      bookmaker_key, bookmaker_title
    ) %>%
    summarise(
      book_over_market  = first(closing_odds_1),
      book_under_market = first(closing_odds_2),
      book_market_total = first(bookmakers_markets_outcomes_point_1),
      .groups = "drop"
    )
  
  # 3) Elihu engine for totals at all posted numbers - using new architecture
  totals_list <- unique(flat_betting_odds$book_market_total)

  total_final_preds <- targets %>%
    mutate(res = pmap(
      list(id, parent_spread, parent_total, target_cover, target_over),
      function(id, ps, pt, tc, to) {
        # Step 1: Get balanced sample
        sample_result <- run_answer_key_sample(
          id = id, parent_spread = ps, parent_total = pt,
          target_cover = tc, target_over = to,
          DT = DT, ss = ss, st = st, N = N
        )
        if (is.null(sample_result)) return(NULL)
        # Step 2: Generate total predictions from sample
        predict_totals_from_sample(sample_result$sample, totals = totals_list, period = period)
      }
    )) %>%
    unnest(res) %>%
    inner_join(
      consensus_odds %>%
        select(id, home_team, away_team, commence_time),
      by = "id"
    ) %>%
    select(home_team, away_team, commence_time, everything())
  
  total_predictions <- total_final_preds %>%
    pivot_longer(
      starts_with("pct_over"),
      names_to  = "market_total_line",
      values_to = "over_prediction"
    ) %>%
    mutate(
      market_total_line = as.numeric(gsub("_", ".", sub("^pct_over_", "", market_total_line)))
    ) %>%
    mutate(
      across(
        starts_with("pct_over"),
        ~ prob_to_american(.x),
        .names = "{.col}_american"
      )
    )
  
  # 4) Join book totals + Elihu probs + consensus, EV + Kelly
  prediction_set <- flat_betting_odds %>%
    left_join(
      total_predictions %>%
        select(id, market_total_line, over_prediction),
      by = c("id", "book_market_total" = "market_total_line")
    ) %>%
    left_join(
      consensus_odds %>% select(-home_team, -away_team),
      by = "id"
    ) %>%
    rename(
      book_full_game_home_prob  = consensus_prob_home,
      book_full_game_away_prob  = consensus_prob_away,
      book_full_game_over_prob  = consensus_over,
      book_full_game_under_prob = consensus_under,
      book_full_game_total_line = total_line
    ) %>%
    mutate(
      under_prediction = 1 - over_prediction
    ) %>%
    mutate(as_tibble(american_prob(book_over_market, book_under_market))) %>%
    rename(
      book_market_prob_over  = p1,
      book_market_prob_under = p2
    ) %>%
    mutate(
      over_ev  = compute_ev(over_prediction,  book_market_prob_over),
      under_ev = compute_ev(under_prediction, book_market_prob_under),
      over_bet_size  = kelly_stake(over_ev,  book_market_prob_over,  bankroll, kelly_mult),
      under_bet_size = kelly_stake(under_ev, book_market_prob_under, bankroll, kelly_mult),
      market         = market
    )
  
  bets <- prediction_set %>%
    format_bets_table(
      pred1 = "over_prediction",
      pred2 = "under_prediction",
      ev1   = "over_ev",
      ev2   = "under_ev",
      size1 = "over_bet_size",
      size2 = "under_bet_size",
      odds1 = "book_over_market",
      odds2 = "book_under_market"
    )
  
  list(
    predictions    = total_predictions,
    prediction_set = prediction_set,
    bets           = bets
  )
}

#Entire Spread Process
build_spread_market <- function(
    DT,
    consensus_odds,
    ss, st, N,
    period,             # inning / quarter / half etc
    events,
    market,             # e.g. "spreads_1st_5_innings"
    sport_key,          # e.g. "baseball_mlb", "americanfootball_nfl"
    bankroll   = 200,
    kelly_mult = 0.25
) {
  # 1) Targets from consensus
  targets <- consensus_odds %>%
    transmute(
      id,
      parent_spread = consensus_prob_home,
      parent_total  = total_line,
      target_cover  = consensus_prob_home,
      target_over   = consensus_over
    )
  
  # 2) Current spread odds
  all_odds <- map_dfr(events$id, ~ fetch_event_odds(.x, market, sport_key))

  flat_betting_odds <- all_odds %>%
    flatten_event_odds() %>%
    # Each bookmaker row has LISTS of outcomes - unnest to get 1 row per outcome
    unnest_longer(c(outcome_name, closing_odds, bookmakers_markets_outcomes_point)) %>%
    mutate(
      side = if_else(outcome_name == home_team, "home", "away"),
      spread = round(bookmakers_markets_outcomes_point, 1),
      odds = closing_odds
    ) %>%
    # Group by spread_line (absolute value) to pair home/away
    mutate(spread_line = abs(spread)) %>%
    group_by(id, commence_time, home_team, away_team, bookmaker_key, bookmaker_title, spread_line) %>%
    summarise(
      book_home_spread = spread[side == "home"][1],
      book_away_spread = spread[side == "away"][1],
      book_home_odds   = odds[side == "home"][1],
      book_away_odds   = odds[side == "away"][1],
      .groups = "drop"
    ) %>%
    select(-spread_line) %>%
    filter(!is.na(book_home_spread), !is.na(book_away_spread))
  
  # 3) Elihu engine for spreads at all posted numbers - using new architecture
  spreads_list <- unique(flat_betting_odds$book_home_spread)

  spread_final_preds <- targets %>%
    mutate(res = pmap(
      list(id, parent_spread, parent_total, target_cover, target_over),
      function(id, ps, pt, tc, to) {
        # Step 1: Get balanced sample
        sample_result <- run_answer_key_sample(
          id = id, parent_spread = ps, parent_total = pt,
          target_cover = tc, target_over = to,
          DT = DT, ss = ss, st = st, N = N
        )
        if (is.null(sample_result)) return(NULL)
        # Step 2: Generate spread predictions from sample
        predict_spreads_from_sample(sample_result$sample, spreads = spreads_list, period = period)
      }
    )) %>%
    unnest(res) %>%
    inner_join(
      consensus_odds %>%
        select(id, home_team, away_team, commence_time),
      by = "id"
    ) %>%
    select(home_team, away_team, commence_time, everything())
  
  spread_predictions <- spread_final_preds %>%
    pivot_longer(
      starts_with("pct_home_cover"),
      names_to  = "market_spread_line",
      values_to = "home_spread_prediction"
    ) %>%
    mutate(
      market_spread = as.numeric(gsub("_", ".", sub("^pct_home_cover.", "", market_spread_line)))
    ) %>%
    mutate(
      across(
        starts_with("home_spread_prediction"),
        ~ prob_to_american(.x),
        .names = "{.col}_american"
      )
    )
  
  # 4) Join spreads + predictions + consensus, EV + Kelly
  prediction_set <- flat_betting_odds %>%
    left_join(
      spread_predictions %>%
        select(id, market_spread, home_spread_prediction),
      by = c("id", "book_home_spread" = "market_spread")
    ) %>%
    left_join(
      consensus_odds %>% select(-home_team, -away_team),
      by = "id"
    ) %>%
    rename(
      book_full_game_home_prob  = consensus_prob_home,
      book_full_game_away_prob  = consensus_prob_away,
      book_full_game_over_prob  = consensus_over,
      book_full_game_under_prob = consensus_under,
      book_full_game_total_line = total_line
    ) %>%
    mutate(
      away_spread_prediction = 1 - home_spread_prediction
    ) %>%
    mutate(as_tibble(american_prob(book_home_odds, book_away_odds))) %>%
    rename(
      book_market_prob_home = p1,
      book_market_prob_away = p2
    ) %>%
    mutate(
      home_ev       = compute_ev(home_spread_prediction, book_market_prob_home),
      away_ev       = compute_ev(away_spread_prediction, book_market_prob_away),
      home_bet_size = kelly_stake(home_ev, book_market_prob_home, bankroll, kelly_mult),
      away_bet_size = kelly_stake(away_ev, book_market_prob_away, bankroll, kelly_mult),
      market        = market
    )
  
  bets <- prediction_set %>%
    format_bets_table(
      pred1 = "home_spread_prediction",
      pred2 = "away_spread_prediction",
      ev1   = "home_ev",
      ev2   = "away_ev",
      size1 = "home_bet_size",
      size2 = "away_bet_size",
      odds1 = "book_home_odds",
      odds2 = "book_away_odds",
      line_col_1 = "book_home_spread",  # Home bet shows home spread
      line_col_2 = "book_away_spread"   # Away bet shows away spread
    ) %>%
    arrange(desc(home_bet_size))
  
  list(
    predictions    = spread_predictions,
    prediction_set = prediction_set,
    bets           = bets
  )
}

# =============================================================================
# BUILD FROM SAMPLES FUNCTIONS (use pre-computed samples)
# =============================================================================

#' Build moneyline predictions from pre-computed samples
#'
#' @param samples Named list from generate_all_samples()
#' @param consensus_odds Consensus odds data frame
#' @param events Events data frame
#' @param periods Vector of periods
#' @param markets Vector of market names
#' @param sport_key Sport key for API
#' @param bankroll Bankroll for Kelly sizing
#' @param kelly_mult Kelly multiplier
#' @param margin_col Column name prefix for margins
build_moneylines_from_samples <- function(
    samples,
    consensus_odds,
    events,
    periods,
    markets,
    sport_key,
    bankroll = 200,
    kelly_mult = 0.25,
    margin_col = "game_home_margin_period",
    pre_fetched_odds = NULL,
    books = NULL
) {
  # Validate inputs
  if (length(periods) != length(markets)) {
    stop("periods and markets must be same length and correspond to each other")
  }

  # 1) Get odds: either from pre-fetched cache or live API
  if (!is.null(pre_fetched_odds)) {
    cat(sprintf("Using %d cached event responses for %d markets\n", nrow(pre_fetched_odds), length(markets)))
    all_odds <- pre_fetched_odds %>%
      mutate(data = map(json_response, ~ {
        parsed <- tryCatch(fromJSON(.x, flatten = TRUE), error = function(e) NULL)
        if (!is.null(parsed) && length(parsed) > 0 && "bookmakers" %in% names(parsed)) {
          as_tibble(parsed)
        } else {
          tibble()
        }
      })) %>%
      filter(map_lgl(data, ~ nrow(.x) > 0)) %>%
      unnest(data) %>%
      select(-json_response, -event_id)
  } else {
    cat(sprintf("Fetching %d markets for %d events (bulk)...\n", length(markets), length(events$id)))
    cached <- fetch_odds_bulk(events$id, markets, sport_key)
    all_odds <- cached %>%
      filter(!is.na(json_response)) %>%
      mutate(data = map(json_response, ~ {
        parsed <- tryCatch(fromJSON(.x, flatten = TRUE), error = function(e) NULL)
        if (!is.null(parsed) && length(parsed) > 0 && "bookmakers" %in% names(parsed)) {
          as_tibble(parsed)
        } else {
          tibble()
        }
      })) %>%
      filter(map_lgl(data, ~ nrow(.x) > 0)) %>%
      unnest(data) %>%
      select(-json_response, -event_id)
  }

  if (nrow(all_odds) == 0) {
    cat("No moneyline odds data returned from API\n")
    return(list(predictions = tibble(), prediction_set = tibble(), bets = tibble(),
                markets_summary = tibble(market = character(), n_bets = integer(),
                                         total_stake = numeric(), avg_ev = numeric(), max_ev = numeric())))
  }

  # 2) Flatten odds - each row has a 2-element list (home/away outcomes)
  # unnest_wider expands these into _1 and _2 columns
  flat_odds <- all_odds %>%
    flatten_event_odds() %>%
    filter(market_key %in% markets) %>%
    rename(market = market_key) %>%
    unnest_wider(outcome_name, names_sep = "_") %>%
    unnest_wider(closing_odds, names_sep = "_") %>%
    mutate(
      home_odds = if_else(home_team == outcome_name_1, closing_odds_1, closing_odds_2),
      away_odds = if_else(away_team == outcome_name_1, closing_odds_1, closing_odds_2)
    ) %>%
    group_by(market, id, commence_time, home_team, away_team, bookmaker_key, bookmaker_title) %>%
    summarise(
      book_home_market = first(home_odds),
      book_away_market = first(away_odds),
      .groups = "drop"
    )

  # 3) Generate predictions from pre-computed samples
  cat(sprintf("Generating predictions from %d pre-computed samples...\n", length(samples)))

  predictions_raw <- map_dfr(names(samples), function(game_id) {
    sample_result <- samples[[game_id]]
    preds <- predict_moneyline_from_sample(sample_result$sample, margin_col = margin_col)
    preds$id <- game_id
    preds
  })

  # Join with consensus to get team names and time
  consensus_info <- consensus_odds %>%
    ungroup() %>%
    select(id, home_team, away_team, commence_time)

  if (is.character(consensus_info$commence_time)) {
    consensus_info <- consensus_info %>%
      mutate(commence_time = ymd_hms(commence_time, tz = "UTC"))
  }

  predictions <- predictions_raw %>%
    inner_join(consensus_info, by = "id")

  # 4) Create market-to-period mapping and pivot predictions
  market_period_map <- tibble(market = markets, period = as.character(periods))

  margin_cols <- names(predictions)[startsWith(names(predictions), margin_col)]
  if (length(margin_cols) == 0) {
    stop(paste("No columns starting with", margin_col, "found in predictions"))
  }

  predictions_long <- predictions %>%
    select(id, home_team, away_team, commence_time, all_of(margin_cols)) %>%
    pivot_longer(
      cols = starts_with(margin_col),
      names_to = "period_col",
      values_to = "home_win_prob"
    ) %>%
    mutate(
      period = str_extract(period_col, "(?<=_)(\\d+|Half\\d+)$")
    ) %>%
    select(-period_col) %>%
    inner_join(market_period_map, by = "period") %>%
    mutate(away_win_prob = 1 - home_win_prob)

  # 5) Join predictions to markets and compute EV
  prediction_set <- flat_odds %>%
    left_join(predictions_long %>% select(-commence_time), by = c("market", "id", "home_team", "away_team")) %>%
    left_join(
      consensus_odds %>%
        ungroup() %>%
        select(id, spread, total_line, starts_with("consensus")),
      by = "id"
    ) %>%
    mutate(as_tibble(american_prob(book_home_market, book_away_market))) %>%
    rename(book_market_prob_home = p1, book_market_prob_away = p2) %>%
    filter(!is.na(spread), !is.na(home_win_prob)) %>%
    mutate(
      home_ev = compute_ev(home_win_prob, book_market_prob_home),
      away_ev = compute_ev(away_win_prob, book_market_prob_away),
      home_bet_size = kelly_stake(home_ev, book_market_prob_home, bankroll, kelly_mult),
      away_bet_size = kelly_stake(away_ev, book_market_prob_away, bankroll, kelly_mult)
    )

  # 6) Format bets
  bets <- prediction_set %>%
    format_bets_table(
      pred1 = "home_win_prob", pred2 = "away_win_prob",
      ev1 = "home_ev", ev2 = "away_ev",
      size1 = "home_bet_size", size2 = "away_bet_size",
      odds1 = "book_home_market", odds2 = "book_away_market",
      books = books
    )

  # 7) Market summary
  summary <- bets %>%
    group_by(market) %>%
    summarise(
      n_bets = n(),
      total_stake = sum(bet_size),
      avg_ev = mean(ev),
      max_ev = max(ev),
      .groups = "drop"
    ) %>%
    arrange(desc(total_stake))

  cat(sprintf("Complete! Generated %d bets across %d markets\n", nrow(bets), length(markets)))

  list(
    predictions = predictions_long,
    prediction_set = prediction_set,
    bets = bets,
    markets_summary = summary
  )
}

#' Build 3-way moneyline predictions from pre-computed samples
#' Handles markets with home/away/tie outcomes (e.g., h2h_3_way_q1)
#'
#' @param samples Named list of sample results from generate_all_samples()
#' @param consensus_odds Consensus odds data frame
#' @param events Events data frame from get_events()
#' @param periods Vector of periods (e.g., c("1", "2", "3", "4", "Half1"))
#' @param markets Vector of 3-way market names (e.g., c("h2h_3_way_q1", ...))
#' @param sport_key Sport key for API
#' @param bankroll Bankroll for Kelly sizing
#' @param kelly_mult Kelly multiplier
#' @param margin_col Column name prefix for margins
#' @return List with predictions, prediction_set, bets, and markets_summary
build_3way_from_samples <- function(
    samples,
    consensus_odds,
    events,
    periods,
    markets,
    sport_key,
    bankroll = 200,
    kelly_mult = 0.25,
    margin_col = "game_home_margin_period"
) {
  # Validate inputs
  if (length(periods) != length(markets)) {
    stop("periods and markets must be same length and correspond to each other")
  }

  cat(sprintf("Fetching %d 3-way markets for %d events (bulk)...\n", length(markets), length(events$id)))

  # 1) Batch fetch all markets for all events
  cached <- fetch_odds_bulk(events$id, markets, sport_key)
  all_odds <- cached %>%
    filter(!is.na(json_response)) %>%
    mutate(data = map(json_response, ~ {
      parsed <- tryCatch(fromJSON(.x, flatten = TRUE), error = function(e) NULL)
      if (!is.null(parsed) && length(parsed) > 0 && "bookmakers" %in% names(parsed)) {
        as_tibble(parsed)
      } else {
        tibble()
      }
    })) %>%
    filter(map_lgl(data, ~ nrow(.x) > 0)) %>%
    unnest(data) %>%
    select(-json_response, -event_id)

  if (nrow(all_odds) == 0) {
    cat("No 3-way odds data returned\n")
    return(list(predictions = tibble(), prediction_set = tibble(), bets = tibble(),
                markets_summary = tibble(market = character(), n_bets = integer(),
                                         total_stake = numeric(), avg_ev = numeric(), max_ev = numeric())))
  }

  # 2) Flatten 3-way odds - each row has 3 outcomes (home/away/tie)
  flat_odds <- all_odds %>%
    flatten_event_odds() %>%
    rename(market = market_key) %>%
    unnest_longer(c(outcome_name, closing_odds)) %>%
    mutate(
      side = case_when(
        outcome_name == home_team ~ "home",
        outcome_name == away_team ~ "away",
        outcome_name %in% c("Draw", "Tie") ~ "tie",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(side)) %>%
    group_by(market, id, commence_time, home_team, away_team, bookmaker_key, bookmaker_title) %>%
    summarise(
      book_home_odds = closing_odds[side == "home"][1],
      book_away_odds = closing_odds[side == "away"][1],
      book_tie_odds = closing_odds[side == "tie"][1],
      .groups = "drop"
    ) %>%
    filter(!is.na(book_home_odds), !is.na(book_away_odds), !is.na(book_tie_odds))

  if (nrow(flat_odds) == 0) {
    cat("No valid 3-way odds after flattening\n")
    return(list(predictions = tibble(), prediction_set = tibble(), bets = tibble(),
                markets_summary = tibble(market = character(), n_bets = integer(),
                                         total_stake = numeric(), avg_ev = numeric(), max_ev = numeric())))
  }

  # 3) Generate predictions from pre-computed samples
  cat(sprintf("Generating 3-way predictions from %d pre-computed samples...\n", length(samples)))

  predictions_raw <- map_dfr(names(samples), function(game_id) {
    sample_result <- samples[[game_id]]
    preds <- predict_moneyline_from_sample(sample_result$sample, margin_col = margin_col)
    preds$id <- game_id
    preds
  })

  # Join with consensus to get team names and time
  consensus_info <- consensus_odds %>%
    ungroup() %>%
    select(id, home_team, away_team, commence_time)

  if (is.character(consensus_info$commence_time)) {
    consensus_info <- consensus_info %>%
      mutate(commence_time = ymd_hms(commence_time, tz = "UTC"))
  }

  predictions <- predictions_raw %>%
    inner_join(consensus_info, by = "id")

  # 4) Create market-to-period mapping and pivot 3-way predictions
  market_period_map <- tibble(market = markets, period = as.character(periods))

  # Get the 3-way prediction columns
  home_3way_cols <- names(predictions)[grepl(paste0(margin_col, ".*_3way_home$"), names(predictions))]
  away_3way_cols <- names(predictions)[grepl(paste0(margin_col, ".*_3way_away$"), names(predictions))]
  tie_3way_cols <- names(predictions)[grepl(paste0(margin_col, ".*_3way_tie$"), names(predictions))]

  if (length(home_3way_cols) == 0) {
    stop("No 3-way prediction columns found")
  }

  # Pivot home probs
  predictions_home <- predictions %>%
    select(id, home_team, away_team, commence_time, all_of(home_3way_cols)) %>%
    pivot_longer(
      cols = all_of(home_3way_cols),
      names_to = "period_col",
      values_to = "home_prob"
    ) %>%
    mutate(period = str_extract(period_col, "(?<=_)(\\d+|Half\\d+)(?=_3way)")) %>%
    select(-period_col)

  # Pivot away probs
  predictions_away <- predictions %>%
    select(id, all_of(away_3way_cols)) %>%
    pivot_longer(
      cols = all_of(away_3way_cols),
      names_to = "period_col",
      values_to = "away_prob"
    ) %>%
    mutate(period = str_extract(period_col, "(?<=_)(\\d+|Half\\d+)(?=_3way)")) %>%
    select(id, period, away_prob)

  # Pivot tie probs
  predictions_tie <- predictions %>%
    select(id, all_of(tie_3way_cols)) %>%
    pivot_longer(
      cols = all_of(tie_3way_cols),
      names_to = "period_col",
      values_to = "tie_prob"
    ) %>%
    mutate(period = str_extract(period_col, "(?<=_)(\\d+|Half\\d+)(?=_3way)")) %>%
    select(id, period, tie_prob)

  # Combine all 3 probs
  predictions_long <- predictions_home %>%
    inner_join(predictions_away, by = c("id", "period")) %>%
    inner_join(predictions_tie, by = c("id", "period")) %>%
    inner_join(market_period_map, by = "period")

  # 5) Join predictions to markets and compute EV
  prediction_set <- flat_odds %>%
    left_join(predictions_long %>% select(-commence_time), by = c("market", "id", "home_team", "away_team")) %>%
    left_join(
      consensus_odds %>%
        ungroup() %>%
        select(id, spread, total_line, starts_with("consensus")),
      by = "id"
    ) %>%
    # Raw implied probs (with vig) — matches 2-way pattern
    mutate(
      book_prob_home = odds_to_prob(book_home_odds),
      book_prob_away = odds_to_prob(book_away_odds),
      book_prob_tie  = odds_to_prob(book_tie_odds)
    ) %>%
    filter(!is.na(spread), !is.na(home_prob)) %>%
    mutate(
      home_ev = compute_ev(home_prob, book_prob_home),
      away_ev = compute_ev(away_prob, book_prob_away),
      tie_ev = compute_ev(tie_prob, book_prob_tie),
      home_bet_size = kelly_stake(home_ev, book_prob_home, bankroll, kelly_mult),
      away_bet_size = kelly_stake(away_ev, book_prob_away, bankroll, kelly_mult),
      tie_bet_size = kelly_stake(tie_ev, book_prob_tie, bankroll, kelly_mult)
    )

  # 6) Format bets (using extended format_bets_table with 3rd side)
  bets <- prediction_set %>%
    format_bets_table(
      pred1 = "home_prob", pred2 = "away_prob",
      ev1 = "home_ev", ev2 = "away_ev",
      size1 = "home_bet_size", size2 = "away_bet_size",
      odds1 = "book_home_odds", odds2 = "book_away_odds",
      pred3 = "tie_prob", ev3 = "tie_ev", size3 = "tie_bet_size", odds3 = "book_tie_odds",
      books = books
    )

  # 7) Market summary
  summary <- bets %>%
    group_by(market) %>%
    summarise(
      n_bets = n(),
      total_stake = sum(bet_size),
      avg_ev = mean(ev),
      max_ev = max(ev),
      .groups = "drop"
    ) %>%
    arrange(desc(total_stake))

  cat(sprintf("Complete! Generated %d 3-way bets across %d markets\n", nrow(bets), length(unique(bets$market))))

  list(
    predictions = predictions_long,
    prediction_set = prediction_set,
    bets = bets,
    markets_summary = summary
  )
}

#' Build totals predictions from pre-computed samples
build_totals_from_samples <- function(
    samples,
    consensus_odds,
    events,
    periods,
    markets,
    sport_key,
    bankroll = 200,
    kelly_mult = 0.25,
    total_col = "game_total_period",
    pre_fetched_odds = NULL,
    books = NULL
) {
  # Validate inputs
  if (length(periods) != length(markets)) {
    stop("periods and markets must be same length and correspond to each other")
  }

  # 1) Get odds: either from pre-fetched cache or live API
  if (!is.null(pre_fetched_odds)) {
    cat(sprintf("Using %d cached event responses for %d totals markets\n", nrow(pre_fetched_odds), length(markets)))
    all_odds <- pre_fetched_odds %>%
      mutate(data = map(json_response, ~ {
        parsed <- tryCatch(fromJSON(.x, flatten = TRUE), error = function(e) NULL)
        if (!is.null(parsed) && length(parsed) > 0 && "bookmakers" %in% names(parsed)) {
          as_tibble(parsed)
        } else {
          tibble()
        }
      })) %>%
      filter(map_lgl(data, ~ nrow(.x) > 0)) %>%
      unnest(data) %>%
      select(-json_response, -event_id)
  } else {
    cat(sprintf("Fetching %d totals markets for %d events (bulk)...\n", length(markets), length(events$id)))
    cached <- fetch_odds_bulk(events$id, markets, sport_key)
    all_odds <- cached %>%
      filter(!is.na(json_response)) %>%
      mutate(data = map(json_response, ~ {
        parsed <- tryCatch(fromJSON(.x, flatten = TRUE), error = function(e) NULL)
        if (!is.null(parsed) && length(parsed) > 0 && "bookmakers" %in% names(parsed)) {
          as_tibble(parsed)
        } else {
          tibble()
        }
      })) %>%
      filter(map_lgl(data, ~ nrow(.x) > 0)) %>%
      unnest(data) %>%
      select(-json_response, -event_id)
  }

  if (nrow(all_odds) == 0) {
    cat("No totals odds data returned from API\n")
    return(list(predictions = tibble(), prediction_set = tibble(), bets = tibble(),
                markets_summary = tibble(market = character(), n_bets = integer(),
                                         total_stake = numeric(), avg_ev = numeric(), max_ev = numeric())))
  }

  # 2) Flatten odds - each row has a 2-element list (over/under outcomes)
  # unnest_wider expands these into _1 and _2 columns
  flat_odds <- all_odds %>%
    flatten_event_odds() %>%
    filter(market_key %in% markets) %>%
    rename(market = market_key) %>%
    unnest_wider(outcome_name, names_sep = "_") %>%
    unnest_wider(closing_odds, names_sep = "_") %>%
    unnest_wider(bookmakers_markets_outcomes_point, names_sep = "_") %>%
    mutate(
      over_odds = if_else(outcome_name_1 == "Over", closing_odds_1, closing_odds_2),
      under_odds = if_else(outcome_name_1 == "Under", closing_odds_1, closing_odds_2),
      total_line_market = if_else(outcome_name_1 == "Over",
                                   bookmakers_markets_outcomes_point_1,
                                   bookmakers_markets_outcomes_point_2)
    ) %>%
    # Group by line value to support alternate lines (multiple lines per bookmaker)
    group_by(market, id, commence_time, home_team, away_team, bookmaker_key, bookmaker_title, total_line_market) %>%
    summarise(
      book_over_market = first(over_odds),
      book_under_market = first(under_odds),
      .groups = "drop"
    ) %>%
    rename(book_total_line = total_line_market)

  # Create market-to-period mapping
  market_period_map <- tibble(market = markets, period = as.character(periods))

  # Get unique total lines per market/period
  totals_by_market <- flat_odds %>%
    left_join(market_period_map, by = "market") %>%
    group_by(market, period) %>%
    summarise(total_lines = list(unique(book_total_line)), .groups = "drop")

  # 3) Generate predictions from pre-computed samples
  cat(sprintf("Generating predictions from %d pre-computed samples...\n", length(samples)))

  predictions_list <- map_dfr(names(samples), function(game_id) {
    sample_result <- samples[[game_id]]

    preds <- market_period_map %>%
      pmap_dfr(function(market, period) {
        totals_for_market <- totals_by_market %>%
          filter(market == !!market) %>%
          pull(total_lines) %>%
          unlist()

        if (length(totals_for_market) == 0) return(tibble())

        pred <- predict_totals_from_sample(
          sample_result$sample,
          totals = totals_for_market,
          total_col = total_col,
          period = period
        )

        if (nrow(pred) == 0) return(tibble())

        pred %>%
          pivot_longer(everything(), names_to = "total_col_name", values_to = "over_prob") %>%
          mutate(
            book_total_line = as.numeric(gsub("_", ".", sub("^pct_over_", "", total_col_name))),
            market = market,
            period = period
          ) %>%
          select(-total_col_name)
      })

    preds$id <- game_id
    preds
  })

  # Join with consensus
  consensus_info <- consensus_odds %>%
    ungroup() %>%
    select(id, home_team, away_team, commence_time)

  if (is.character(consensus_info$commence_time)) {
    consensus_info <- consensus_info %>%
      mutate(commence_time = ymd_hms(commence_time, tz = "UTC"))
  }

  predictions_long <- predictions_list %>%
    inner_join(consensus_info, by = "id") %>%
    mutate(under_prob = 1 - over_prob)

  # 4) Join and compute EV
  prediction_set <- flat_odds %>%
    left_join(market_period_map, by = "market") %>%
    left_join(
      predictions_long %>% select(id, market, book_total_line, over_prob, under_prob),
      by = c("market", "id", "book_total_line")
    ) %>%
    left_join(
      consensus_odds %>% ungroup() %>% select(id, spread, total_line, starts_with("consensus")),
      by = "id"
    ) %>%
    mutate(as_tibble(american_prob(book_over_market, book_under_market))) %>%
    rename(book_market_prob_over = p1, book_market_prob_under = p2) %>%
    filter(!is.na(over_prob)) %>%
    mutate(
      over_ev = compute_ev(over_prob, book_market_prob_over),
      under_ev = compute_ev(under_prob, book_market_prob_under),
      over_bet_size = kelly_stake(over_ev, book_market_prob_over, bankroll, kelly_mult),
      under_bet_size = kelly_stake(under_ev, book_market_prob_under, bankroll, kelly_mult)
    )

  # 5) Format bets
  bets <- prediction_set %>%
    format_bets_table(
      pred1 = "over_prob", pred2 = "under_prob",
      ev1 = "over_ev", ev2 = "under_ev",
      size1 = "over_bet_size", size2 = "under_bet_size",
      odds1 = "book_over_market", odds2 = "book_under_market",
      line_col = "book_total_line",
      books = books
    )

  # 6) Market summary
  summary <- bets %>%
    group_by(market) %>%
    summarise(n_bets = n(), total_stake = sum(bet_size), avg_ev = mean(ev), max_ev = max(ev), .groups = "drop") %>%
    arrange(desc(total_stake))

  cat(sprintf("Complete! Generated %d totals bets across %d markets\n", nrow(bets), length(markets)))

  list(predictions = predictions_long, prediction_set = prediction_set, bets = bets, markets_summary = summary)
}

#' Build spreads predictions from pre-computed samples
build_spreads_from_samples <- function(
    samples,
    consensus_odds,
    events,
    periods,
    markets,
    sport_key,
    bankroll = 200,
    kelly_mult = 0.25,
    margin_col = "game_home_margin_period",
    pre_fetched_odds = NULL,
    books = NULL
) {
  # Validate inputs
  if (length(periods) != length(markets)) {
    stop("periods and markets must be same length and correspond to each other")
  }

  # 1) Get odds: either from pre-fetched cache or live API
  if (!is.null(pre_fetched_odds)) {
    cat(sprintf("Using %d cached event responses for %d spreads markets\n", nrow(pre_fetched_odds), length(markets)))
    all_odds <- pre_fetched_odds %>%
      mutate(data = map(json_response, ~ {
        parsed <- tryCatch(fromJSON(.x, flatten = TRUE), error = function(e) NULL)
        if (!is.null(parsed) && length(parsed) > 0 && "bookmakers" %in% names(parsed)) {
          as_tibble(parsed)
        } else {
          tibble()
        }
      })) %>%
      filter(map_lgl(data, ~ nrow(.x) > 0)) %>%
      unnest(data) %>%
      select(-json_response, -event_id)
  } else {
    cat(sprintf("Fetching %d spreads markets for %d events (bulk)...\n", length(markets), length(events$id)))
    cached <- fetch_odds_bulk(events$id, markets, sport_key)
    all_odds <- cached %>%
      filter(!is.na(json_response)) %>%
      mutate(data = map(json_response, ~ {
        parsed <- tryCatch(fromJSON(.x, flatten = TRUE), error = function(e) NULL)
        if (!is.null(parsed) && length(parsed) > 0 && "bookmakers" %in% names(parsed)) {
          as_tibble(parsed)
        } else {
          tibble()
        }
      })) %>%
      filter(map_lgl(data, ~ nrow(.x) > 0)) %>%
      unnest(data) %>%
      select(-json_response, -event_id)
  }

  if (nrow(all_odds) == 0) {
    cat("No spreads odds data returned from API\n")
    return(list(predictions = tibble(), prediction_set = tibble(), bets = tibble(),
                markets_summary = tibble(market = character(), n_bets = integer(),
                                         total_stake = numeric(), avg_ev = numeric(), max_ev = numeric())))
  }

  # 2) Flatten odds - each bookmaker row has LISTS of outcomes (multiple spread lines)
  # First unnest_longer to get one row per outcome, then group to pair home/away
  flat_odds <- all_odds %>%
    flatten_event_odds() %>%
    filter(market_key %in% markets) %>%
    rename(market = market_key) %>%
    # Each row has lists of 24+ outcomes - unnest to get 1 row per outcome
    unnest_longer(c(outcome_name, closing_odds, bookmakers_markets_outcomes_point)) %>%
    # Assign home/away side based on team name
    mutate(
      side = if_else(outcome_name == home_team, "home", "away"),
      spread = round(bookmakers_markets_outcomes_point, 1),
      odds = closing_odds
    ) %>%
    # Group by spread_line (absolute value) to pair home/away for each line
    mutate(spread_line = abs(spread)) %>%
    group_by(market, id, commence_time, home_team, away_team, bookmaker_key, bookmaker_title, spread_line) %>%
    summarise(
      book_home_spread = spread[side == "home"][1],
      book_away_spread = spread[side == "away"][1],
      book_home_market = odds[side == "home"][1],
      book_away_market = odds[side == "away"][1],
      .groups = "drop"
    ) %>%
    select(-spread_line) %>%
    # Remove rows where we don't have both sides
    filter(!is.na(book_home_spread), !is.na(book_away_spread))

  # Create market-to-period mapping
  market_period_map <- tibble(market = markets, period = as.character(periods))

  # Get unique spreads per market/period
  spreads_by_market <- flat_odds %>%
    left_join(market_period_map, by = "market") %>%
    group_by(market, period) %>%
    summarise(spread_lines = list(unique(book_home_spread)), .groups = "drop")

  # 3) Generate predictions from pre-computed samples
  cat(sprintf("Generating predictions from %d pre-computed samples...\n", length(samples)))

  predictions_list <- map_dfr(names(samples), function(game_id) {
    sample_result <- samples[[game_id]]

    preds <- market_period_map %>%
      pmap_dfr(function(market, period) {
        spreads_for_market <- spreads_by_market %>%
          filter(market == !!market) %>%
          pull(spread_lines) %>%
          unlist()

        if (length(spreads_for_market) == 0) return(tibble())

        pred <- predict_spreads_from_sample(
          sample_result$sample,
          spreads = spreads_for_market,
          margin_col = margin_col,
          period = period
        )

        if (nrow(pred) == 0) return(tibble())

        pred %>%
          pivot_longer(everything(), names_to = "spread_col_name", values_to = "home_cover_prob") %>%
          mutate(
            spread_str = sub("^pct_home_cover_", "", spread_col_name),
            spread_str = gsub("neg", "-", spread_str),
            spread_str = gsub("_", ".", spread_str),
            # Round to match the rounding in flat_odds
            book_home_spread = round(as.numeric(spread_str), 1),
            market = market,
            period = period
          ) %>%
          select(-spread_col_name, -spread_str)
      })

    preds$id <- game_id
    preds
  })

  # Join with consensus
  consensus_info <- consensus_odds %>%
    ungroup() %>%
    select(id, home_team, away_team, commence_time)

  if (is.character(consensus_info$commence_time)) {
    consensus_info <- consensus_info %>%
      mutate(commence_time = ymd_hms(commence_time, tz = "UTC"))
  }

  predictions_long <- predictions_list %>%
    inner_join(consensus_info, by = "id") %>%
    mutate(away_cover_prob = 1 - home_cover_prob)

  # 4) Join and compute EV
  prediction_set <- flat_odds %>%
    left_join(market_period_map, by = "market") %>%
    left_join(
      predictions_long %>% select(id, market, period, book_home_spread, home_cover_prob, away_cover_prob),
      by = c("market", "id", "period", "book_home_spread")
    ) %>%
    left_join(
      consensus_odds %>% ungroup() %>% select(id, spread, total_line, starts_with("consensus")),
      by = "id"
    ) %>%
    mutate(as_tibble(american_prob(book_home_market, book_away_market))) %>%
    rename(book_market_prob_home = p1, book_market_prob_away = p2) %>%
    filter(!is.na(home_cover_prob)) %>%
    mutate(
      home_ev = compute_ev(home_cover_prob, book_market_prob_home),
      away_ev = compute_ev(away_cover_prob, book_market_prob_away),
      home_bet_size = kelly_stake(home_ev, book_market_prob_home, bankroll, kelly_mult),
      away_bet_size = kelly_stake(away_ev, book_market_prob_away, bankroll, kelly_mult)
    )

  # 5) Format bets
  bets <- prediction_set %>%
    format_bets_table(
      pred1 = "home_cover_prob", pred2 = "away_cover_prob",
      ev1 = "home_ev", ev2 = "away_ev",
      size1 = "home_bet_size", size2 = "away_bet_size",
      odds1 = "book_home_market", odds2 = "book_away_market",
      line_col_1 = "book_home_spread",  # Home bet shows home spread (e.g., +1.5)
      line_col_2 = "book_away_spread",  # Away bet shows away spread (e.g., -1.5)
      books = books
    )

  # 6) Market summary
  summary <- bets %>%
    group_by(market) %>%
    summarise(n_bets = n(), total_stake = sum(bet_size), avg_ev = mean(ev), max_ev = max(ev), .groups = "drop") %>%
    arrange(desc(total_stake))

  cat(sprintf("Complete! Generated %d spreads bets across %d markets\n", nrow(bets), length(markets)))

  list(predictions = predictions_long, prediction_set = prediction_set, bets = bets, markets_summary = summary)
}

# =============================================================================
# CONVENIENCE WRAPPERS (generate samples internally)
# =============================================================================

# Efficient Multi-Market Moneyline Builder ----
# NOTE: This is now a convenience wrapper around the new explicit sampling architecture.
# For better efficiency when using multiple market types, use generate_all_samples() once
# and then call build_moneylines_from_samples(), build_totals_from_samples(),
# build_spreads_from_samples() with the pre-computed samples.
build_multi_moneyline_markets <- function(
    DT,
    consensus_odds,
    ss, st, N,
    periods,            # Vector of periods: c("1", "2", "3", "4", "h1", "h2")
    events,
    markets,            # Vector of markets: c("h2h_q1", "h2h_q2", "h2h_q3", "h2h_q4", "h2h_h1", "h2h_h2")
    sport_key,
    bankroll   = 200,
    kelly_mult = 0.25,
    targets = NULL,
    use_spread_line = FALSE,
    margin_col = "game_home_margin_period"
) {

  # Build targets if not provided
  if (is.null(targets)) {
    targets <- consensus_odds %>%
      ungroup() %>%
      transmute(
        id,
        parent_spread = if_else(use_spread_line, spread, consensus_prob_home),
        parent_total  = total_line,
        target_cover  = consensus_prob_home,
        target_over   = consensus_prob_over
      )
  }

  # Generate samples (the expensive part)
  samples <- generate_all_samples(
    targets         = targets,
    DT              = DT,
    ss              = ss,
    st              = st,
    N               = N,
    use_spread_line = use_spread_line
  )

  # Build moneylines from samples
  build_moneylines_from_samples(
    samples        = samples,
    consensus_odds = consensus_odds,
    events         = events,
    periods        = periods,
    markets        = markets,
    sport_key      = sport_key,
    bankroll       = bankroll,
    kelly_mult     = kelly_mult,
    margin_col     = margin_col
  )
}

# Efficient Multi-Market Totals Builder ----
# NOTE: This is now a convenience wrapper around the new explicit sampling architecture.
# For better efficiency when using multiple market types, use generate_all_samples() once
# and then call build_moneylines_from_samples(), build_totals_from_samples(),
# build_spreads_from_samples() with the pre-computed samples.
build_multi_totals_markets <- function(
    DT,
    consensus_odds,
    ss, st, N,
    periods,            # Vector of periods: c("1", "2", "3", "4", "Half1")
    events,
    markets,            # Vector of markets: c("totals_q1", "totals_q2", "totals_q3", "totals_q4", "totals_h1")
    sport_key,
    bankroll   = 200,
    kelly_mult = 0.25,
    targets = NULL,
    use_spread_line = TRUE,
    total_col = "game_total_period"
) {

  # Build targets if not provided
  if (is.null(targets)) {
    targets <- consensus_odds %>%
      ungroup() %>%
      transmute(
        id,
        parent_spread = if_else(use_spread_line, spread, consensus_prob_home),
        parent_total  = total_line,
        target_cover  = consensus_prob_home,
        target_over   = consensus_prob_over
      )
  }

  # Generate samples (the expensive part)
  samples <- generate_all_samples(
    targets         = targets,
    DT              = DT,
    ss              = ss,
    st              = st,
    N               = N,
    use_spread_line = use_spread_line
  )

  # Build totals from samples
  build_totals_from_samples(
    samples        = samples,
    consensus_odds = consensus_odds,
    events         = events,
    periods        = periods,
    markets        = markets,
    sport_key      = sport_key,
    bankroll       = bankroll,
    kelly_mult     = kelly_mult,
    total_col      = total_col
  )
}

# Efficient Multi-Market Spreads Builder ----
# NOTE: This is now a convenience wrapper around the new explicit sampling architecture.
# For better efficiency when using multiple market types, use generate_all_samples() once
# and then call build_moneylines_from_samples(), build_totals_from_samples(),
# build_spreads_from_samples() with the pre-computed samples.
build_multi_spreads_markets <- function(
    DT,
    consensus_odds,
    ss, st, N,
    periods,            # Vector of periods: c("1", "2", "3", "4", "Half1")
    events,
    markets,            # Vector of markets: c("spreads_q1", "spreads_q2", "spreads_q3", "spreads_q4", "spreads_h1")
    sport_key,
    bankroll   = 200,
    kelly_mult = 0.25,
    targets = NULL,
    use_spread_line = TRUE,
    margin_col = "game_home_margin_period"
) {

  # Build targets if not provided
  if (is.null(targets)) {
    targets <- consensus_odds %>%
      ungroup() %>%
      transmute(
        id,
        parent_spread = if_else(use_spread_line, spread, consensus_prob_home),
        parent_total  = total_line,
        target_cover  = consensus_prob_home,
        target_over   = consensus_prob_over
      )
  }

  # Generate samples (the expensive part)
  samples <- generate_all_samples(
    targets         = targets,
    DT              = DT,
    ss              = ss,
    st              = st,
    N               = N,
    use_spread_line = use_spread_line
  )

  # Build spreads from samples
  build_spreads_from_samples(
    samples        = samples,
    consensus_odds = consensus_odds,
    events         = events,
    periods        = periods,
    markets        = markets,
    sport_key      = sport_key,
    bankroll       = bankroll,
    kelly_mult     = kelly_mult,
    margin_col     = margin_col
  )
}

# Generic Consensus Builder ----
build_market_consensus <- function(
    eval_dt,
    market_type,           # "moneyline" or "totals"
    game_id_col = "game_id",
    date_col = "game_date",
    outcome_col,           # "home_winner" for ML, "over_hit" for totals
    prob_col,              # "devig_home_odds" for ML, "devig_over_odds" for totals
    min_count_1yr = 50
) {
  # 1) Filter valid rows
  df <- eval_dt %>%
    filter(!is.na(.data[[prob_col]]), !is.na(.data[[outcome_col]]))
  
  if (market_type != "moneyline") {
    df <- df %>% filter(.data[[prob_col]] > 0.4, .data[[prob_col]] < 0.6)
  }
  
  # 2) Compute logit
  df <- df %>% mutate(logit = logit_(.data[[prob_col]]))
  
  # 3) Score books
  score_tbl <- df %>%
    mutate(
      logloss = logloss_(.data[[prob_col]], .data[[outcome_col]]),
      age_days = as.numeric(Sys.Date() - as.Date(.data[[date_col]]))
    ) %>%
    group_by(bookmaker_key) %>%
    summarise(
      avg_logloss = mean(logloss),
      avg_logloss_1yr = mean(logloss[age_days <= 365], na.rm = TRUE),
      count_1yr = sum(age_days <= 365),
      .groups = "drop"
    ) %>%
    arrange(avg_logloss)
  
  # 4) Weights
  weight_col_name <- paste0(market_type, "_weight")
  weights <- score_tbl %>%
    mutate(!!weight_col_name := ifelse(count_1yr > min_count_1yr, 1 / avg_logloss_1yr, 0)) %>%
    filter(.data[[weight_col_name]] > 0) %>%
    select(bookmaker_key, !!weight_col_name)
  
  # 5) Build consensus
  consensus <- df %>%
    inner_join(weights, by = "bookmaker_key") %>%
    group_by(.data[[game_id_col]], .data[[date_col]]) %>%
    summarise(
      consensus_logit = sum(.data[[weight_col_name]] * logit) / sum(.data[[weight_col_name]]),
      consensus_p = invlogit_(consensus_logit),
      .groups = "drop"
    )
  
  # 6) Evaluate
  eval_tbl <- df %>%
    inner_join(consensus, by = c(game_id_col, date_col)) %>%
    mutate(
      logloss_book = logloss_(.data[[prob_col]], .data[[outcome_col]]),
      logloss_cons = logloss_(consensus_p, .data[[outcome_col]])
    )
  
  comparison <- tibble(
    model = c("Consensus", "Average Book"),
    logloss = c(mean(eval_tbl$logloss_cons, na.rm = TRUE),
                mean(eval_tbl$logloss_book, na.rm = TRUE))
  )
  
  list(
    score_tbl = score_tbl,
    weights = weights,
    consensus = consensus,
    eval_tbl = eval_tbl,
    comparison = comparison
  )
}

# ============================================================================
# TEAM TOTALS FUNCTIONS
# ============================================================================

#' Predict team totals from a sample
#'
#' Calculates probability of team scoring over/under various lines
#'
#' @param sample Sample data frame containing team score columns
#' @param team_totals Vector of total lines to evaluate
#' @param team "home" or "away"
#' @param period Period identifier (e.g., "1", "2", "Half1")
#' @return Tibble with probability columns for each total line
predict_team_totals_from_sample <- function(
    sample,
    team_totals,
    team = c("home", "away"),
    period
) {
  team <- match.arg(team)
  col_name <- paste0(team, "_score_period_", period)

  if (!(col_name %in% names(sample))) {
    warning(paste("Column", col_name, "not found in sample"))
    return(tibble())
  }

  score_vals <- sample[[col_name]]

  # For each total, calculate probability of over
  pct_cols <- set_names(
    map(team_totals, ~ sum(score_vals > .x, na.rm = TRUE) / sum(score_vals != .x, na.rm = TRUE)),
    paste0("pct_over_", team, "_", gsub("\\.", "_", as.character(team_totals)))
  )

  tibble(!!!pct_cols)
}

#' Build team totals predictions from pre-computed samples
#'
#' @param samples Named list of sample results from generate_all_samples()
#' @param consensus_odds Consensus odds data frame
#' @param events Events data frame from get_events()
#' @param periods Vector of periods (e.g., c("1", "2", "3", "4", "Half1"))
#' @param markets Vector of team totals market names (e.g., c("team_totals_q1", ...))
#' @param sport_key Sport key for API
#' @param bankroll Bankroll for Kelly sizing
#' @param kelly_mult Kelly multiplier
#' @return List with predictions, prediction_set, bets, and markets_summary
build_team_totals_from_samples <- function(
    samples,
    consensus_odds,
    events,
    periods,
    markets,
    sport_key,
    bankroll = 200,
    kelly_mult = 0.25,
    ev_threshold = 0.05,
    pre_fetched_odds = NULL,
    books = NULL
) {
  # Validate inputs
  if (length(periods) != length(markets)) {
    stop("periods and markets must have same length")
  }

  # Create market-to-period mapping
  market_period_map <- tibble(market = markets, period = as.character(periods))

  # 1) Get odds: either from pre-fetched cache or live API
  if (!is.null(pre_fetched_odds)) {
    cat(sprintf("Using %d cached event responses for %d team totals markets\n", nrow(pre_fetched_odds), length(markets)))
    raw_odds <- pre_fetched_odds %>%
      mutate(data = map(json_response, ~ {
        parsed <- tryCatch(fromJSON(.x, flatten = FALSE), error = function(e) NULL)
        if (!is.null(parsed) && length(parsed$bookmakers) > 0) {
          as_tibble(parsed)
        } else {
          tibble()
        }
      })) %>%
      filter(map_lgl(data, ~ nrow(.x) > 0)) %>%
      unnest(data) %>%
      select(-json_response, -event_id)
  } else {
    cat(sprintf("Fetching %d team totals markets for %d events (bulk)...\n", length(markets), nrow(events)))
    cached <- fetch_odds_bulk(events$id, markets, sport_key)
    raw_odds <- cached %>%
      filter(!is.na(json_response)) %>%
      mutate(data = map(json_response, ~ {
        parsed <- tryCatch(fromJSON(.x, flatten = FALSE), error = function(e) NULL)
        if (!is.null(parsed) && length(parsed$bookmakers) > 0) {
          as_tibble(parsed)
        } else {
          tibble()
        }
      })) %>%
      filter(map_lgl(data, ~ nrow(.x) > 0)) %>%
      unnest(data) %>%
      select(-json_response, -event_id)
  }

  if (nrow(raw_odds) == 0) {
    cat("No team totals odds returned from API\n")
    return(list(predictions = tibble(), prediction_set = tibble(), bets = tibble(),
                markets_summary = tibble(market = character(), n_bets = integer(),
                                         total_stake = numeric(), avg_ev = numeric(), max_ev = numeric())))
  }

  # 2) Flatten the odds - team totals have description field for team name
  flat_odds <- raw_odds %>%
    unnest_longer(bookmakers) %>%
    unnest_wider(bookmakers, names_sep = "_") %>%
    # Properly handle multi-market data frames
    mutate(bookmakers_markets = map(bookmakers_markets, ~ {
      if (is.list(.x) && !is.data.frame(.x)) bind_rows(.x) else .x
    })) %>%
    unnest(bookmakers_markets, names_sep = "_") %>%
    filter(bookmakers_markets_key %in% markets) %>%
    # outcomes is a list of lists - unnest to individual rows
    unnest_longer(bookmakers_markets_outcomes) %>%
    unnest_wider(bookmakers_markets_outcomes, names_sep = "_") %>%
    # These columns may still be lists if multiple values - unnest them
    unnest_longer(c(bookmakers_markets_outcomes_name,
                    bookmakers_markets_outcomes_description,
                    bookmakers_markets_outcomes_price,
                    bookmakers_markets_outcomes_point)) %>%
    mutate(
      commence_time = ymd_hms(commence_time, tz = "UTC")
    ) %>%
    rename(
      bookmaker_key = bookmakers_key,
      bookmaker_title = bookmakers_title,
      market = bookmakers_markets_key,
      outcome_name = bookmakers_markets_outcomes_name,  # "Over" or "Under"
      closing_odds = bookmakers_markets_outcomes_price,
      team_total_line = bookmakers_markets_outcomes_point,
      team_name = bookmakers_markets_outcomes_description  # Team name
    ) %>%
    filter(!is.na(team_name), !is.na(team_total_line)) %>%
    # Determine if this is home or away team
    mutate(
      team_side = if_else(team_name == home_team, "home", "away")
    ) %>%
    # Pivot to have over/under odds per team per line
    group_by(market, id, commence_time, home_team, away_team, bookmaker_key, bookmaker_title, team_name, team_side, team_total_line) %>%
    summarise(
      book_over_odds = closing_odds[outcome_name == "Over"][1],
      book_under_odds = closing_odds[outcome_name == "Under"][1],
      .groups = "drop"
    ) %>%
    filter(!is.na(book_over_odds), !is.na(book_under_odds))

  # Filter to enabled books
  if (!is.null(books)) {
    flat_odds <- flat_odds %>% filter(bookmaker_key %in% books)
  }

  if (nrow(flat_odds) == 0) {
    cat("No valid team totals odds after flattening\n")
    return(list(predictions = tibble(), prediction_set = tibble(), bets = tibble(),
                markets_summary = tibble(market = character(), n_bets = integer(),
                                         total_stake = numeric(), avg_ev = numeric(), max_ev = numeric())))
  }

  # Get unique team totals per market/period/team
  totals_by_market <- flat_odds %>%
    left_join(market_period_map, by = "market") %>%
    group_by(market, period, team_side) %>%
    summarise(total_lines = list(unique(team_total_line)), .groups = "drop")

  # 3) Generate predictions from pre-computed samples
  cat(sprintf("Generating team totals predictions from %d pre-computed samples...\n", length(samples)))

  predictions_list <- map_dfr(names(samples), function(game_id) {
    sample_result <- samples[[game_id]]

    # For each market/period/team combination
    preds <- totals_by_market %>%
      pmap_dfr(function(market, period, team_side, total_lines) {
        totals_for_market <- unlist(total_lines)
        if (length(totals_for_market) == 0) return(tibble())

        pred <- predict_team_totals_from_sample(
          sample_result$sample,
          team_totals = totals_for_market,
          team = team_side,
          period = period
        )

        if (nrow(pred) == 0) return(tibble())

        pred %>%
          pivot_longer(everything(), names_to = "total_col_name", values_to = "over_prob") %>%
          mutate(
            # Extract team and line from column name: pct_over_home_10_5 -> home, 10.5
            parts = str_match(total_col_name, "pct_over_(home|away)_(.+)"),
            team_side = parts[,2],
            team_total_line = as.numeric(gsub("_", ".", parts[,3])),
            market = market,
            period = period
          ) %>%
          select(-total_col_name, -parts)
      })

    preds$id <- game_id
    preds
  })

  # Join with consensus
  consensus_info <- consensus_odds %>%
    ungroup() %>%
    select(id, home_team, away_team, commence_time)

  if (is.character(consensus_info$commence_time)) {
    consensus_info <- consensus_info %>%
      mutate(commence_time = ymd_hms(commence_time, tz = "UTC"))
  }

  predictions_long <- predictions_list %>%
    inner_join(consensus_info, by = "id") %>%
    mutate(under_prob = 1 - over_prob)

  # 4) Join and compute EV
  prediction_set <- flat_odds %>%
    left_join(market_period_map, by = "market") %>%
    left_join(
      predictions_long %>% select(id, market, team_side, team_total_line, over_prob, under_prob),
      by = c("market", "id", "team_side", "team_total_line")
    ) %>%
    left_join(
      consensus_odds %>% ungroup() %>% select(id, spread, total_line, starts_with("consensus")),
      by = "id"
    ) %>%
    mutate(as_tibble(american_prob(book_over_odds, book_under_odds))) %>%
    rename(book_market_prob_over = p1, book_market_prob_under = p2) %>%
    filter(!is.na(over_prob)) %>%
    mutate(
      over_ev = compute_ev(over_prob, book_market_prob_over),
      under_ev = compute_ev(under_prob, book_market_prob_under),
      over_bet_size = kelly_stake(over_ev, book_market_prob_over, bankroll, kelly_mult),
      under_bet_size = kelly_stake(under_ev, book_market_prob_under, bankroll, kelly_mult)
    )

  # 5) Format bets - need custom formatting since bet_on should be team name
  bets_over <- prediction_set %>%
    filter(over_bet_size > 0, over_ev >= ev_threshold) %>%
    transmute(
      id,
      market,
      home_team,
      away_team,
      commence_time,
      bet_on = paste0(team_name, " Over"),
      line = team_total_line,
      prob = over_prob,
      ev = over_ev,
      bet_size = over_bet_size,
      odds = book_over_odds,
      bookmaker_key
    )

  bets_under <- prediction_set %>%
    filter(under_bet_size > 0, under_ev >= ev_threshold) %>%
    transmute(
      id,
      market,
      home_team,
      away_team,
      commence_time,
      bet_on = paste0(team_name, " Under"),
      line = team_total_line,
      prob = under_prob,
      ev = under_ev,
      bet_size = under_bet_size,
      odds = book_under_odds,
      bookmaker_key
    )

  bets <- bind_rows(bets_over, bets_under) %>%
    # Convert commence_time to pt_start_time for consistency with other bet tables
    mutate(
      pt_start_time = with_tz(commence_time, "America/Los_Angeles")
    ) %>%
    select(-commence_time) %>%
    arrange(desc(ev))

  # 6) Market summary
  summary <- bets %>%
    group_by(market) %>%
    summarise(n_bets = n(), total_stake = sum(bet_size), avg_ev = mean(ev), max_ev = max(ev), .groups = "drop") %>%
    arrange(desc(total_stake))

  cat(sprintf("Complete! Generated %d team totals bets across %d markets\n", nrow(bets), length(markets)))

  list(
    predictions = predictions_long,
    prediction_set = prediction_set,
    bets = bets,
    markets_summary = summary
  )
}

# Generic PBP + Odds Join ----
join_pbp_odds <- function(
    pbp_dt,
    odds_dt,
    join_cols = c("home_team", "away_team", "game_date", "game_id", "game_start_time")
) {
  pbp_dt <- as.data.table(pbp_dt)
  odds_dt <- as.data.table(odds_dt)

  pbp_dt[, game_start_time := ymd_hms(game_start_time, tz = "UTC")]
  pbp_dt[, game_date := as.Date(game_date)]
  odds_dt[, game_start_time := ymd_hms(game_start_time, tz = "UTC")]
  odds_dt[, game_date := as.Date(game_date)]

  setkeyv(pbp_dt, join_cols)
  setkeyv(odds_dt, join_cols)

  return(pbp_dt[odds_dt, on = join_cols, roll = "nearest", nomatch = 0L])
}

# =============================================================================
# Wagerzon Offshore Odds Integration
# =============================================================================

#' Run Wagerzon scraper to fetch current odds
#'
#' @param sport Sport to scrape ("nfl", "cbb", "nba")
#' @param scraper_dir Path to the wagerzon_odds directory
#' @param venv_path Path to Python virtual environment
#' @return TRUE if successful, FALSE otherwise
run_wagerzon_scraper <- function(
    sport = "nfl",
    scraper_dir = "~/NFLWork/wagerzon_odds",
    venv_path = "venv"
) {
  scraper_dir <- normalizePath(path.expand(scraper_dir), mustWork = TRUE)

  # Build the command to run the scraper
  cmd <- sprintf(
    "cd '%s' && source %s/bin/activate && python3 scraper_v2.py %s 2>&1",
    scraper_dir,
    venv_path,
    sport
  )

  cat(sprintf("Running Wagerzon %s scraper...\n", toupper(sport)))

  # Execute the scraper
  result <- system(cmd, intern = TRUE)

  # Print output
  cat(paste(result, collapse = "\n"), "\n")

  # Check for success (look for "Saved" in output)
  success <- any(grepl("Saved|records", result))

  if (success) {
    cat("Wagerzon scraper completed successfully.\n")
  } else {
    warning("Wagerzon scraper may have failed. Check output above.")
  }

  return(success)
}


#' Get Wagerzon odds from database
#'
#' @param sport Sport to retrieve ("nfl", "cbb", "nba")
#' @param db_path Path to wagerzon.duckdb
#' @return data.frame of odds in standardized format
get_wagerzon_odds <- function(
    sport = "nfl",
    db_path = "~/NFLWork/wagerzon_odds/wagerzon.duckdb"
) {
  db_path <- normalizePath(path.expand(db_path), mustWork = TRUE)

  table_name <- paste0(sport, "_odds")

  con <- dbConnect(duckdb(), dbdir = db_path, read_only = TRUE)
  on.exit(dbDisconnect(con, shutdown = TRUE))

  # Check if table exists
  tables <- dbListTables(con)
  if (!table_name %in% tables) {
    warning(sprintf("Table '%s' not found in Wagerzon database. Run scraper first.", table_name))
    return(data.frame())
  }

  # Read raw data
  raw_odds <- dbGetQuery(con, sprintf("SELECT * FROM %s", table_name))

  if (nrow(raw_odds) == 0) {
    warning("No odds found in Wagerzon database.")
    return(data.frame())
  }

  # Transform to match The Odds API format
  # Create separate rows for spreads, totals, and moneylines

  result_list <- list()

  for (i in seq_len(nrow(raw_odds))) {
    row <- raw_odds[i, ]

    # Base info for all record types
    base <- list(
      bookmaker_key = "wagerzon",
      sport_key = row$sport_key,
      home_team = row$home_team,
      away_team = row$away_team,
      game_date = row$game_date,
      game_time = row$game_time,
      wagerzon_game_id = row$game_id,
      fetch_time = row$fetch_time,
      period = row$period
    )

    # Spreads record
    if (!is.na(row$away_spread)) {
      market_name <- row$market
      spread_rec <- c(base, list(
        market = market_name,
        market_type = "spreads",
        line = row$home_spread,
        odds_away = row$away_spread_price,
        odds_home = row$home_spread_price,
        odds_over = NA_integer_,
        odds_under = NA_integer_,
        away_spread = row$away_spread,
        home_spread = row$home_spread
      ))
      result_list[[length(result_list) + 1]] <- spread_rec
    }

    # Totals record
    if (!is.na(row$total)) {
      totals_market <- gsub("spreads", "totals", row$market)
      totals_rec <- c(base, list(
        market = totals_market,
        market_type = "totals",
        line = row$total,
        odds_away = NA_integer_,
        odds_home = NA_integer_,
        odds_over = row$over_price,
        odds_under = row$under_price
      ))
      result_list[[length(result_list) + 1]] <- totals_rec
    }

    # Moneyline record
    if (!is.na(row$away_ml)) {
      ml_market <- gsub("spreads", "h2h", row$market)
      ml_rec <- c(base, list(
        market = ml_market,
        market_type = "h2h",
        line = NA_real_,
        odds_away = row$away_ml,
        odds_home = row$home_ml,
        odds_over = NA_integer_,
        odds_under = NA_integer_
      ))
      result_list[[length(result_list) + 1]] <- ml_rec
    }
  }

  # Convert list to data frame
  if (length(result_list) == 0) {
    return(data.frame())
  }

  result <- bind_rows(lapply(result_list, as.data.frame))

  cat(sprintf("Loaded %d Wagerzon odds records (%d spreads, %d totals, %d ML)\n",
              nrow(result),
              sum(result$market_type == "spreads"),
              sum(result$market_type == "totals"),
              sum(result$market_type == "h2h")))

  result <- resolve_offshore_teams(result, sport = sport)
  return(result)
}


#' Load Hoop88 odds from DuckDB
#'
#' @param sport Sport key (e.g., "nfl", "ncaaf")
#' @param db_path Path to the Hoop88 DuckDB database
#' @return Data frame with Hoop88 odds in standardized format
get_hoop88_odds <- function(
    sport = "nfl",
    db_path = "~/NFLWork/hoop88_odds/hoop88.duckdb"
) {
  db_path <- normalizePath(path.expand(db_path), mustWork = FALSE)

  if (!file.exists(db_path)) {
    warning(sprintf("Hoop88 database not found at %s. Run scraper first.", db_path))
    return(data.frame())
  }

  table_name <- paste0(sport, "_odds")

  con <- dbConnect(duckdb(), dbdir = db_path, read_only = TRUE)
  on.exit(dbDisconnect(con, shutdown = TRUE))

  tables <- dbListTables(con)
  if (!table_name %in% tables) {
    warning(sprintf("Table '%s' not found in Hoop88 database. Run scraper first.", table_name))
    return(data.frame())
  }

  raw_odds <- dbGetQuery(con, sprintf("SELECT * FROM %s", table_name))

  if (nrow(raw_odds) == 0) {
    warning("No odds found in Hoop88 database.")
    return(data.frame())
  }

  result_list <- list()

  for (i in seq_len(nrow(raw_odds))) {
    row <- raw_odds[i, ]

    base <- list(
      bookmaker_key = "hoop88",
      sport_key = row$sport_key,
      home_team = row$home_team,
      away_team = row$away_team,
      game_date = row$game_date,
      game_time = row$game_time,
      hoop88_game_id = row$game_id,
      fetch_time = row$fetch_time,
      period = row$period
    )

    # Spreads record
    if (!is.na(row$away_spread)) {
      market_name <- row$market
      spread_rec <- c(base, list(
        market = market_name,
        market_type = "spreads",
        line = row$home_spread,
        odds_away = row$away_spread_price,
        odds_home = row$home_spread_price,
        odds_over = NA_integer_,
        odds_under = NA_integer_,
        away_spread = row$away_spread,
        home_spread = row$home_spread
      ))
      result_list[[length(result_list) + 1]] <- spread_rec
    }

    # Totals record
    if (!is.na(row$total)) {
      totals_market <- gsub("spreads", "totals", row$market)
      totals_rec <- c(base, list(
        market = totals_market,
        market_type = "totals",
        line = row$total,
        odds_away = NA_integer_,
        odds_home = NA_integer_,
        odds_over = row$over_price,
        odds_under = row$under_price
      ))
      result_list[[length(result_list) + 1]] <- totals_rec
    }

    # Moneyline record
    if (!is.na(row$away_ml)) {
      ml_market <- gsub("spreads", "h2h", row$market)
      ml_rec <- c(base, list(
        market = ml_market,
        market_type = "h2h",
        line = NA_real_,
        odds_away = row$away_ml,
        odds_home = row$home_ml,
        odds_over = NA_integer_,
        odds_under = NA_integer_
      ))
      result_list[[length(result_list) + 1]] <- ml_rec
    }
  }

  if (length(result_list) == 0) {
    return(data.frame())
  }

  result <- bind_rows(lapply(result_list, as.data.frame))

  # Filter to correct sport only (scraper may include wrong sport data)
  expected_key <- switch(sport,
    "nfl" = "americanfootball_nfl",
    "ncaaf" = "americanfootball_ncaaf",
    "cbb" = "basketball_ncaab",
    "nba" = "basketball_nba",
    "college_baseball" = "baseball_ncaa",
    sport)
  if ("sport_key" %in% names(result)) {
    result <- result %>% filter(sport_key == expected_key)
  }

  cat(sprintf("Loaded %d Hoop88 odds records (%d spreads, %d totals, %d ML)\n",
              nrow(result),
              sum(result$market_type == "spreads"),
              sum(result$market_type == "totals"),
              sum(result$market_type == "h2h")))

  result <- resolve_offshore_teams(result, sport = sport)
  return(result)
}


#' Load BFA Gaming odds from DuckDB
#'
#' Loads scraped BFA Gaming derivative odds and transforms them into a standardized
#' format compatible with the Wagerzon comparison functions.
#'
#' @param sport Sport key (e.g., "nfl", "nba")
#' @param db_path Path to the BFA DuckDB database
#' @return Data frame with BFA odds in standardized format
get_bfa_odds <- function(
    sport = "nfl",
    db_path = "~/NFLWork/bfa_odds/bfa.duckdb"
) {
  db_path <- normalizePath(path.expand(db_path), mustWork = FALSE)

  if (!file.exists(db_path)) {
    warning(sprintf("BFA database not found at %s. Run scraper first.", db_path))
    return(data.frame())
  }

  table_name <- paste0(sport, "_odds")

  con <- dbConnect(duckdb(), dbdir = db_path, read_only = TRUE)
  on.exit(dbDisconnect(con, shutdown = TRUE))

  tables <- dbListTables(con)
  if (!table_name %in% tables) {
    warning(sprintf("Table '%s' not found in BFA database. Run scraper first.", table_name))
    return(data.frame())
  }

  raw_odds <- dbGetQuery(con, sprintf("SELECT * FROM %s", table_name))

  if (nrow(raw_odds) == 0) {
    warning("No odds found in BFA database.")
    return(data.frame())
  }

  result_list <- list()

  for (i in seq_len(nrow(raw_odds))) {
    row <- raw_odds[i, ]

    base <- list(
      bookmaker_key = "bfa",
      sport_key = row$sport_key,
      home_team = row$home_team,
      away_team = row$away_team,
      game_date = row$game_date,
      game_time = row$game_time,
      bfa_game_id = row$game_id,
      fetch_time = row$fetch_time,
      period = row$period
    )

    # Spreads record
    if (!is.na(row$away_spread)) {
      market_name <- row$market
      spread_rec <- c(base, list(
        market = market_name,
        market_type = "spreads",
        line = row$home_spread,
        odds_away = row$away_spread_price,
        odds_home = row$home_spread_price,
        odds_over = NA_integer_,
        odds_under = NA_integer_,
        away_spread = row$away_spread,
        home_spread = row$home_spread
      ))
      result_list[[length(result_list) + 1]] <- spread_rec
    }

    # Totals record — skip team totals (handled separately below)
    if (!is.na(row$total) && !grepl("team_totals_", row$market)) {
      totals_market <- gsub("spreads", "totals", row$market)
      totals_rec <- c(base, list(
        market = totals_market,
        market_type = "totals",
        line = row$total,
        odds_away = NA_integer_,
        odds_home = NA_integer_,
        odds_over = row$over_price,
        odds_under = row$under_price
      ))
      result_list[[length(result_list) + 1]] <- totals_rec
    }

    # Moneyline record
    if (!is.na(row$away_ml)) {
      ml_market <- gsub("spreads", "h2h", row$market)
      ml_rec <- c(base, list(
        market = ml_market,
        market_type = "h2h",
        line = NA_real_,
        odds_away = row$away_ml,
        odds_home = row$home_ml,
        odds_over = NA_integer_,
        odds_under = NA_integer_
      ))
      result_list[[length(result_list) + 1]] <- ml_rec
    }

    # Team totals record
    if (grepl("team_totals_", row$market) && !is.na(row$total)) {
      tt_rec <- c(base, list(
        market = row$market,
        market_type = "team_totals",
        line = row$total,
        odds_away = NA_integer_,
        odds_home = NA_integer_,
        odds_over = row$over_price,
        odds_under = row$under_price
      ))
      result_list[[length(result_list) + 1]] <- tt_rec
    }
  }

  if (length(result_list) == 0) {
    return(data.frame())
  }

  result <- bind_rows(lapply(result_list, as.data.frame))

  cat(sprintf("Loaded %d BFA odds records (%d spreads, %d totals, %d ML, %d team_totals)\n",
              nrow(result),
              sum(result$market_type == "spreads"),
              sum(result$market_type == "totals"),
              sum(result$market_type == "h2h"),
              sum(result$market_type == "team_totals")))

  result <- resolve_offshore_teams(result, sport = sport)
  return(result)
}


get_bookmaker_odds <- function(
    sport = "cbb",
    db_path = "~/NFLWork/bookmaker_odds/bookmaker.duckdb"
) {
  db_path <- normalizePath(path.expand(db_path), mustWork = FALSE)

  if (!file.exists(db_path)) {
    warning(sprintf("Bookmaker database not found at %s. Run scraper first.", db_path))
    return(data.frame())
  }

  table_name <- paste0(sport, "_odds")

  con <- dbConnect(duckdb(), dbdir = db_path, read_only = TRUE)
  on.exit(dbDisconnect(con, shutdown = TRUE))

  tables <- dbListTables(con)
  if (!table_name %in% tables) {
    warning(sprintf("Table '%s' not found in Bookmaker database. Run scraper first.", table_name))
    return(data.frame())
  }

  raw_odds <- dbGetQuery(con, sprintf("SELECT * FROM %s", table_name))

  if (nrow(raw_odds) == 0) {
    warning("No odds found in Bookmaker database.")
    return(data.frame())
  }

  result_list <- list()

  for (i in seq_len(nrow(raw_odds))) {
    row <- raw_odds[i, ]

    base <- list(
      bookmaker_key = "bookmaker",
      sport_key = row$sport_key,
      home_team = row$home_team,
      away_team = row$away_team,
      game_date = row$game_date,
      game_time = row$game_time,
      bfa_game_id = row$game_id,
      fetch_time = row$fetch_time,
      period = row$period
    )

    # Spreads record
    if (!is.na(row$away_spread)) {
      market_name <- row$market
      spread_rec <- c(base, list(
        market = market_name,
        market_type = "spreads",
        line = row$home_spread,
        odds_away = row$away_spread_price,
        odds_home = row$home_spread_price,
        odds_over = NA_integer_,
        odds_under = NA_integer_,
        away_spread = row$away_spread,
        home_spread = row$home_spread
      ))
      result_list[[length(result_list) + 1]] <- spread_rec
    }

    # Totals record
    if (!is.na(row$total)) {
      totals_market <- gsub("spreads", "totals", row$market)
      totals_rec <- c(base, list(
        market = totals_market,
        market_type = "totals",
        line = row$total,
        odds_away = NA_integer_,
        odds_home = NA_integer_,
        odds_over = row$over_price,
        odds_under = row$under_price
      ))
      result_list[[length(result_list) + 1]] <- totals_rec
    }

    # Moneyline record
    if (!is.na(row$away_ml)) {
      ml_market <- gsub("spreads", "h2h", row$market)
      ml_rec <- c(base, list(
        market = ml_market,
        market_type = "h2h",
        line = NA_real_,
        odds_away = row$away_ml,
        odds_home = row$home_ml,
        odds_over = NA_integer_,
        odds_under = NA_integer_
      ))
      result_list[[length(result_list) + 1]] <- ml_rec
    }
  }

  if (length(result_list) == 0) {
    return(data.frame())
  }

  result <- bind_rows(lapply(result_list, as.data.frame))

  cat(sprintf("Loaded %d Bookmaker odds records (%d spreads, %d totals, %d ML)\n",
              nrow(result),
              sum(result$market_type == "spreads"),
              sum(result$market_type == "totals"),
              sum(result$market_type == "h2h")))

  result <- resolve_offshore_teams(result, sport = sport)
  return(result)
}

get_bet105_odds <- function(
    sport = "cbb",
    db_path = "~/NFLWork/bet105_odds/bet105.duckdb"
) {
  db_path <- normalizePath(path.expand(db_path), mustWork = FALSE)

  if (!file.exists(db_path)) {
    warning(sprintf("Bet105 database not found at %s. Run scraper first.", db_path))
    return(data.frame())
  }

  table_name <- paste0(sport, "_odds")

  con <- dbConnect(duckdb(), dbdir = db_path, read_only = TRUE)
  on.exit(dbDisconnect(con, shutdown = TRUE))

  tables <- dbListTables(con)
  if (!table_name %in% tables) {
    warning(sprintf("Table '%s' not found in Bet105 database. Run scraper first.", table_name))
    return(data.frame())
  }

  raw_odds <- dbGetQuery(con, sprintf("SELECT * FROM %s", table_name))

  if (nrow(raw_odds) == 0) {
    warning("No odds found in Bet105 database.")
    return(data.frame())
  }

  result_list <- list()

  for (i in seq_len(nrow(raw_odds))) {
    row <- raw_odds[i, ]

    base <- list(
      bookmaker_key = "bet105",
      sport_key = row$sport_key,
      home_team = row$home_team,
      away_team = row$away_team,
      game_date = row$game_date,
      game_time = row$game_time,
      bfa_game_id = row$game_id,
      fetch_time = row$fetch_time,
      period = row$period
    )

    # Spreads record
    if (!is.na(row$away_spread)) {
      market_name <- row$market
      spread_rec <- c(base, list(
        market = market_name,
        market_type = "spreads",
        line = row$home_spread,
        odds_away = row$away_spread_price,
        odds_home = row$home_spread_price,
        odds_over = NA_integer_,
        odds_under = NA_integer_,
        away_spread = row$away_spread,
        home_spread = row$home_spread
      ))
      result_list[[length(result_list) + 1]] <- spread_rec
    }

    # Totals record (exclude team totals which have their own block)
    if (!is.na(row$total) && !grepl("team_totals_", row$market)) {
      totals_market <- gsub("spreads", "totals", row$market)
      totals_rec <- c(base, list(
        market = totals_market,
        market_type = "totals",
        line = row$total,
        odds_away = NA_integer_,
        odds_home = NA_integer_,
        odds_over = row$over_price,
        odds_under = row$under_price
      ))
      result_list[[length(result_list) + 1]] <- totals_rec
    }

    # Moneyline record
    if (!is.na(row$away_ml)) {
      ml_market <- gsub("spreads", "h2h", row$market)
      ml_rec <- c(base, list(
        market = ml_market,
        market_type = "h2h",
        line = NA_real_,
        odds_away = row$away_ml,
        odds_home = row$home_ml,
        odds_over = NA_integer_,
        odds_under = NA_integer_
      ))
      result_list[[length(result_list) + 1]] <- ml_rec
    }

    # Team totals record
    if (grepl("team_totals_", row$market) && !is.na(row$total)) {
      tt_rec <- c(base, list(
        market = row$market,
        market_type = "team_totals",
        line = row$total,
        odds_away = NA_integer_,
        odds_home = NA_integer_,
        odds_over = row$over_price,
        odds_under = row$under_price
      ))
      result_list[[length(result_list) + 1]] <- tt_rec
    }
  }

  if (length(result_list) == 0) {
    return(data.frame())
  }

  result <- bind_rows(lapply(result_list, as.data.frame))

  cat(sprintf("Loaded %d Bet105 odds records (%d spreads, %d totals, %d ML, %d team_totals)\n",
              nrow(result),
              sum(result$market_type == "spreads"),
              sum(result$market_type == "totals"),
              sum(result$market_type == "h2h"),
              sum(result$market_type == "team_totals")))

  result <- resolve_offshore_teams(result, sport = sport)
  return(result)
}


get_kalshi_odds <- function(
    sport = "cbb",
    db_path = "~/NFLWork/kalshi_odds/kalshi.duckdb"
) {
  db_path <- normalizePath(path.expand(db_path), mustWork = FALSE)

  if (!file.exists(db_path)) {
    warning(sprintf("Kalshi database not found at %s. Run scraper first.", db_path))
    return(data.frame())
  }

  table_name <- paste0(sport, "_odds")

  con <- dbConnect(duckdb(), dbdir = db_path, read_only = TRUE)
  on.exit(dbDisconnect(con, shutdown = TRUE))

  tables <- dbListTables(con)
  if (!table_name %in% tables) {
    warning(sprintf("Table '%s' not found in Kalshi database. Run scraper first.", table_name))
    return(data.frame())
  }

  raw_odds <- dbGetQuery(con, sprintf("SELECT * FROM %s", table_name))

  if (nrow(raw_odds) == 0) {
    warning("No odds found in Kalshi database.")
    return(data.frame())
  }

  result_list <- list()

  for (i in seq_len(nrow(raw_odds))) {
    row <- raw_odds[i, ]

    base <- list(
      bookmaker_key = "kalshi",
      sport_key = row$sport_key,
      home_team = row$home_team,
      away_team = row$away_team,
      game_date = row$game_date,
      game_time = row$game_time,
      bfa_game_id = row$game_id,
      fetch_time = row$fetch_time,
      period = row$period
    )

    # Spreads record
    if (!is.na(row$away_spread)) {
      spread_rec <- c(base, list(
        market = row$market,
        market_type = "spreads",
        line = row$home_spread,
        odds_away = row$away_spread_price,
        odds_home = row$home_spread_price,
        cents_away = if ("away_spread_cents" %in% names(row)) row$away_spread_cents else NA_real_,
        cents_home = if ("home_spread_cents" %in% names(row)) row$home_spread_cents else NA_real_,
        odds_over = NA_integer_,
        odds_under = NA_integer_,
        away_spread = row$away_spread,
        home_spread = row$home_spread
      ))
      result_list[[length(result_list) + 1]] <- spread_rec
    }

    # Totals record (Kalshi writes separate rows per market, so market is already correct)
    if (!is.na(row$total) && !grepl("team_totals_", row$market)) {
      totals_rec <- c(base, list(
        market = row$market,
        market_type = "totals",
        line = row$total,
        odds_away = NA_integer_,
        odds_home = NA_integer_,
        odds_over = row$over_price,
        odds_under = row$under_price,
        cents_over = if ("over_cents" %in% names(row)) row$over_cents else NA_real_,
        cents_under = if ("under_cents" %in% names(row)) row$under_cents else NA_real_
      ))
      result_list[[length(result_list) + 1]] <- totals_rec
    }

    # Moneyline record
    if (!is.na(row$away_ml)) {
      has_tie <- "tie_ml" %in% names(row) && !is.na(row$tie_ml)
      ml_rec <- c(base, list(
        market = row$market,
        market_type = if (has_tie) "h2h_3way" else "h2h",
        line = NA_real_,
        odds_away = row$away_ml,
        odds_home = row$home_ml,
        odds_tie = if (has_tie) row$tie_ml else NA_integer_,
        cents_away = if ("away_ml_cents" %in% names(row)) row$away_ml_cents else NA_real_,
        cents_home = if ("home_ml_cents" %in% names(row)) row$home_ml_cents else NA_real_,
        cents_tie = if (has_tie && "tie_ml_cents" %in% names(row)) row$tie_ml_cents else NA_real_,
        odds_over = NA_integer_,
        odds_under = NA_integer_
      ))
      result_list[[length(result_list) + 1]] <- ml_rec
    }
  }

  if (length(result_list) == 0) {
    return(data.frame())
  }

  result <- bind_rows(lapply(result_list, as.data.frame))

  cat(sprintf("Loaded %d Kalshi odds records (%d spreads, %d totals, %d ML 2-way, %d ML 3-way)\n",
              nrow(result),
              sum(result$market_type == "spreads"),
              sum(result$market_type == "totals"),
              sum(result$market_type == "h2h"),
              sum(result$market_type == "h2h_3way")))

  result <- resolve_offshore_teams(result, sport = sport)
  return(result)
}


#' Resolve offshore scraped team names to Odds API format
#'
#' Two-layer approach matching canonical_match.py:
#' Layer 1: cbb_team_dict lookup (short_name/nickname/abbreviation -> odds_api_name)
#' Layer 2: Game-level matching against canonical games from cbb_odds_temp
#'
#' @param odds_df Data frame with home_team and away_team columns
#' @param sport Sport key ("cbb", "nfl", etc.)
#' @return odds_df with resolved team names
resolve_offshore_teams <- function(odds_df, sport = "cbb") {
  if (nrow(odds_df) == 0) return(odds_df)

  dict_db <- sprintf("~/NFLWork/Answer Keys/%s.duckdb", sport)
  dict_db <- normalizePath(path.expand(dict_db), mustWork = FALSE)
  if (!file.exists(dict_db)) return(odds_df)

  con <- dbConnect(duckdb(), dbdir = dict_db, read_only = TRUE)
  on.exit(dbDisconnect(con, shutdown = TRUE))

  # Layer 1: team dict lookup
  team_dict <- tryCatch(
    dbGetQuery(con, paste0(sport, "_team_dict") %>%
      {sprintf("SELECT short_name, nickname, abbreviation, odds_api_name FROM %s WHERE odds_api_name IS NOT NULL", .)}),
    error = function(e) data.frame()
  )

  # Layer 2: canonical games from API consensus
  canonical <- tryCatch(
    dbGetQuery(con, paste0(sport, "_odds_temp") %>%
      {sprintf("SELECT DISTINCT home_team, away_team FROM %s", .)}),
    error = function(e) data.frame()
  )

  if (nrow(team_dict) == 0 && nrow(canonical) == 0) return(odds_df)

  # Strip quotes + punctuation (keep spaces)
  strip_chars <- function(x) {
    x <- gsub("\u2018|\u2019|\u201C|\u201D|'|`|\"", "", x)
    gsub("[.\\-]", "", tolower(trimws(x)))
  }

  # Strip everything non-alphanumeric (for space-vs-dash matching)
  strip_all <- function(x) gsub("[^a-z0-9]", "", tolower(x))

  # Build lookup: lowercased name variants -> odds_api_name
  lookup <- character(0)
  if (nrow(team_dict) > 0) {
    for (i in seq_len(nrow(team_dict))) {
      api_name <- team_dict$odds_api_name[i]
      for (col in c("short_name", "nickname", "abbreviation")) {
        val <- team_dict[[col]][i]
        if (!is.na(val) && nchar(trimws(val)) > 0) {
          lookup[tolower(trimws(val))] <- api_name
          lookup[strip_chars(val)] <- api_name
          lookup[strip_all(val)] <- api_name
        }
      }
    }
  }

  # Also index canonical names (API format) and their no-mascot variants
  all_canonical <- character(0)
  if (nrow(canonical) > 0) {
    all_canonical <- unique(c(canonical$home_team, canonical$away_team))
    for (name in all_canonical) {
      lookup[tolower(name)] <- name
      lookup[strip_chars(name)] <- name
      lookup[strip_all(name)] <- name
      no_mascot <- sub(" [^ ]+$", "", name)
      if (nchar(no_mascot) > 0) {
        lookup[tolower(no_mascot)] <- name
        lookup[strip_chars(no_mascot)] <- name
        lookup[strip_all(no_mascot)] <- name
      }
    }
  }

  # Resolution function for a single name (Layer 1 + substring matching)
  resolve_one <- function(raw_name) {
    low <- tolower(trimws(raw_name))
    stripped <- strip_chars(raw_name)
    squished <- strip_all(raw_name)

    # Direct lookup (try low, stripped, then squished/no-space variants)
    hit <- lookup[low]
    if (!is.na(hit)) return(unname(hit))
    hit <- lookup[stripped]
    if (!is.na(hit)) return(unname(hit))
    hit <- lookup[squished]
    if (!is.na(hit)) return(unname(hit))

    # Substring match against canonical names (both sides stripped)
    if (length(all_canonical) > 0) {
      for (cn in all_canonical) {
        cn_stripped <- strip_chars(cn)
        cn_no_mascot <- strip_chars(sub(" [^ ]+$", "", cn))
        cn_squished <- strip_all(sub(" [^ ]+$", "", cn))
        if (nchar(cn_no_mascot) >= 4 &&
            (stripped == cn_no_mascot ||
             squished == cn_squished ||
             grepl(stripped, cn_stripped, fixed = TRUE) ||
             grepl(cn_no_mascot, stripped, fixed = TRUE))) {
          return(cn)
        }
      }
    }

    NULL  # unresolved
  }

  # Game-level matching (Layer 2): match scraped game pair to canonical game pair
  # Uses word-overlap scoring to find best canonical game match
  resolve_game <- function(scraped_away, scraped_home) {
    if (nrow(canonical) == 0) return(list(away = scraped_away, home = scraped_home))

    sa_words <- unlist(strsplit(strip_chars(scraped_away), "\\s+"))
    sh_words <- unlist(strsplit(strip_chars(scraped_home), "\\s+"))
    sa_words <- sa_words[nchar(sa_words) >= 3]
    sh_words <- sh_words[nchar(sh_words) >= 3]

    best_score <- 0
    best_away <- scraped_away
    best_home <- scraped_home

    for (j in seq_len(nrow(canonical))) {
      ca <- canonical$away_team[j]
      ch <- canonical$home_team[j]
      ca_words <- unlist(strsplit(strip_chars(ca), "\\s+"))
      ch_words <- unlist(strsplit(strip_chars(ch), "\\s+"))
      ca_words <- ca_words[nchar(ca_words) >= 3]
      ch_words <- ch_words[nchar(ch_words) >= 3]

      # Score: how many scraped words appear in canonical name
      away_score <- sum(sa_words %in% ca_words) + sum(ca_words %in% sa_words)
      home_score <- sum(sh_words %in% ch_words) + sum(ch_words %in% sh_words)

      # Both sides must have some match
      if (away_score >= 2 && home_score >= 2) {
        total <- away_score + home_score
        if (total > best_score) {
          best_score <- total
          best_away <- ca
          best_home <- ch
        }
      }
    }

    list(away = best_away, home = best_home)
  }

  # Step 1: Try resolve_one for each team individually
  # Step 2: For unresolved games, try game-level matching
  resolved <- 0
  for (i in seq_len(nrow(odds_df))) {
    old_away <- odds_df$away_team[i]
    old_home <- odds_df$home_team[i]

    r_away <- resolve_one(old_away)
    r_home <- resolve_one(old_home)

    if (!is.null(r_away) && !is.null(r_home)) {
      odds_df$away_team[i] <- r_away
      odds_df$home_team[i] <- r_home
    } else {
      # Game-level matching fallback
      game_match <- resolve_game(
        if (!is.null(r_away)) r_away else old_away,
        if (!is.null(r_home)) r_home else old_home
      )
      odds_df$away_team[i] <- game_match$away
      odds_df$home_team[i] <- game_match$home
    }

    if (odds_df$away_team[i] != old_away || odds_df$home_team[i] != old_home) {
      resolved <- resolved + 1
    }
  }

  n_games <- n_distinct(paste(odds_df$away_team, odds_df$home_team))
  cat(sprintf("Resolved %d records (%d unique games) via team name matching.\n", resolved, n_games))
  odds_df
}


#' Merge Wagerzon odds with The Odds API consensus odds
#'
#' @param api_odds Odds from The Odds API (data.frame)
#' @param wagerzon_odds Odds from Wagerzon (data.frame from get_wagerzon_odds)
#' @param match_cols Columns to match games on
#' @return Combined data.frame with both sources
merge_odds_sources <- function(
    api_odds,
    wagerzon_odds,
    match_cols = c("home_team", "away_team")
) {
  if (nrow(wagerzon_odds) == 0) {
    cat("No Wagerzon odds to merge.\n")
    return(api_odds)
  }

  # Ensure api_odds has bookmaker_key column
  if (!"bookmaker_key" %in% names(api_odds)) {
    api_odds$bookmaker_key <- "consensus"
  }

  # Find common columns between the two data sources
  common_cols <- intersect(names(api_odds), names(wagerzon_odds))

  # Identify games that match between sources
  api_games <- api_odds %>%
    distinct(across(all_of(match_cols))) %>%
    mutate(in_api = TRUE)

  wagerzon_games <- wagerzon_odds %>%
    distinct(across(all_of(match_cols))) %>%
    mutate(in_wagerzon = TRUE)

  matched_games <- inner_join(api_games, wagerzon_games, by = match_cols)

  cat(sprintf("Game matching: %d API games, %d Wagerzon games, %d matched\n",
              nrow(api_games), nrow(wagerzon_games), nrow(matched_games)))

  if (nrow(matched_games) == 0) {
    warning("No games matched between API and Wagerzon. Check team name normalization.")
    return(api_odds)
  }

  # Filter Wagerzon odds to only matched games
  wagerzon_matched <- wagerzon_odds %>%
    semi_join(matched_games, by = match_cols)

  # Select common columns and combine
  api_subset <- api_odds %>% select(any_of(common_cols))
  wagerzon_subset <- wagerzon_matched %>% select(any_of(common_cols))

  # Combine the data
  combined <- bind_rows(
    api_subset,
    wagerzon_subset
  )

  cat(sprintf("Merged odds: %d total records (%d from API, %d from Wagerzon)\n",
              nrow(combined),
              nrow(api_subset),
              nrow(wagerzon_subset)))

  return(combined)
}


#' Get Wagerzon odds for specific markets, formatted for betting
#'
#' @param sport Sport ("nfl", "cbb", "nba")
#' @param markets Vector of market types to include (e.g., c("spreads", "totals"))
#' @param periods Vector of periods to include (e.g., c("fg", "h1", "q1"))
#' @param db_path Path to wagerzon.duckdb
#' @return data.frame of filtered odds
get_wagerzon_betting_odds <- function(
    sport = "nfl",
    markets = c("spreads", "totals"),
    periods = c("fg", "h1", "q1", "q2", "q3", "q4"),
    db_path = "~/NFLWork/wagerzon_odds/wagerzon.duckdb"
) {
  odds <- get_wagerzon_odds(sport = sport, db_path = db_path)

  if (nrow(odds) == 0) {
    return(odds)
  }

  # Filter by market type and period
  filtered <- odds %>%
    filter(market_type %in% markets) %>%
    filter(period %in% periods)

  cat(sprintf("Filtered to %d records for markets: %s, periods: %s\n",
              nrow(filtered),
              paste(markets, collapse = ", "),
              paste(periods, collapse = ", ")))

  return(filtered)
}


#' Compare model predictions to Wagerzon odds for spreads
#'
#' @param spread_results Output from build_spreads_from_samples()
#' @param wagerzon_odds Wagerzon odds from get_wagerzon_odds()
#' @param bankroll Bankroll for Kelly sizing
#' @param kelly_mult Kelly multiplier (fraction of Kelly)
#' @param ev_threshold Minimum EV to include bet
#' @return List with predictions, prediction_set, bets, and markets_summary
compare_spreads_to_wagerzon <- function(
    spread_results,
    wagerzon_odds,
    bankroll = 100,
    kelly_mult = 0.25,
    ev_threshold = 0.05
) {
  predictions <- spread_results$predictions

  if (is.null(predictions) || nrow(predictions) == 0) {
    warning("No predictions found in spread_results")
    return(list(bets = data.frame()))
  }

  # Filter Wagerzon odds to spreads only
  wz_spreads <- wagerzon_odds %>%
    filter(market_type == "spreads")

  if (nrow(wz_spreads) == 0) {
    warning("No spreads found in Wagerzon odds")
    return(list(bets = data.frame()))
  }

  # Join predictions with Wagerzon odds
  # Match on: home_team, away_team, market, spread line
  joined <- wz_spreads %>%
    inner_join(
      predictions %>%
        select(id, home_team, away_team, market, book_home_spread, home_cover_prob, away_cover_prob, commence_time),
      by = c("home_team", "away_team", "market"),
      relationship = "many-to-many"
    ) %>%
    filter(abs(home_spread - book_home_spread) < 0.1)

  if (nrow(joined) == 0) {
    cat("No matches found between predictions and Wagerzon odds.\n")
    cat("Prediction markets:", unique(predictions$market), "\n")
    cat("Wagerzon markets:", unique(wz_spreads$market), "\n")
    return(list(bets = data.frame()))
  }

  cat(sprintf("Found %d matches between predictions and Wagerzon odds\n", nrow(joined)))

  # Compute EV against Wagerzon prices
  prediction_set <- joined %>%
    mutate(as_tibble(american_prob(odds_away, odds_home))) %>%
    rename(wz_prob_away = p1, wz_prob_home = p2) %>%
    mutate(
      home_ev = compute_ev(home_cover_prob, wz_prob_home),
      away_ev = compute_ev(away_cover_prob, wz_prob_away),
      home_bet_size = kelly_stake(home_ev, wz_prob_home, bankroll, kelly_mult),
      away_bet_size = kelly_stake(away_ev, wz_prob_away, bankroll, kelly_mult)
    )

  # Get bookmaker key from input data (works for wagerzon, hoop88, etc.)
  book_key <- unique(wz_spreads$bookmaker_key)[1]

  bets <- prediction_set %>%
    format_bets_table(
      pred1 = "home_cover_prob", pred2 = "away_cover_prob",
      ev1 = "home_ev", ev2 = "away_ev",
      size1 = "home_bet_size", size2 = "away_bet_size",
      odds1 = "odds_home", odds2 = "odds_away",
      cents1 = "cents_home", cents2 = "cents_away",
      line_col_1 = "home_spread",
      line_col_2 = "away_spread",
      books = book_key,
      ev_threshold = ev_threshold
    )

  summary <- bets %>%
    group_by(market) %>%
    summarise(n_bets = n(), total_stake = sum(bet_size), avg_ev = mean(ev), max_ev = max(ev), .groups = "drop") %>%
    arrange(desc(total_stake))

  cat(sprintf("Generated %d %s spread bets\n", nrow(bets), book_key))

  list(predictions = predictions, prediction_set = prediction_set, bets = bets, markets_summary = summary)
}


#' Compare model predictions to Wagerzon odds for totals
#'
#' @param totals_results Output from build_totals_from_samples()
#' @param wagerzon_odds Wagerzon odds from get_wagerzon_odds()
#' @param bankroll Bankroll for Kelly sizing
#' @param kelly_mult Kelly multiplier (fraction of Kelly)
#' @param ev_threshold Minimum EV to include bet
#' @return List with bets and markets_summary
compare_totals_to_wagerzon <- function(
    totals_results,
    wagerzon_odds,
    bankroll = 100,
    kelly_mult = 0.25,
    ev_threshold = 0.05
) {
  predictions <- totals_results$predictions

  if (is.null(predictions) || nrow(predictions) == 0) {
    warning("No predictions found in totals_results")
    return(list(bets = data.frame()))
  }

  wz_totals <- wagerzon_odds %>%
    filter(market_type == "totals")

  if (nrow(wz_totals) == 0) {
    warning("No totals found in Wagerzon odds")
    return(list(bets = data.frame()))
  }

  joined <- wz_totals %>%
    inner_join(
      predictions %>%
        select(id, home_team, away_team, market, book_total_line, over_prob, under_prob, commence_time),
      by = c("home_team", "away_team", "market"),
      relationship = "many-to-many"
    ) %>%
    filter(abs(line - book_total_line) < 0.1)

  if (nrow(joined) == 0) {
    cat("No matches found between predictions and Wagerzon totals.\n")
    return(list(bets = data.frame()))
  }

  cat(sprintf("Found %d matches between predictions and Wagerzon totals\n", nrow(joined)))

  prediction_set <- joined %>%
    mutate(as_tibble(american_prob(odds_over, odds_under))) %>%
    rename(wz_prob_over = p1, wz_prob_under = p2) %>%
    mutate(
      over_ev = compute_ev(over_prob, wz_prob_over),
      under_ev = compute_ev(under_prob, wz_prob_under),
      over_bet_size = kelly_stake(over_ev, wz_prob_over, bankroll, kelly_mult),
      under_bet_size = kelly_stake(under_ev, wz_prob_under, bankroll, kelly_mult)
    )

  # Get bookmaker key from input data
  book_key <- unique(wz_totals$bookmaker_key)[1]

  bets <- prediction_set %>%
    format_bets_table(
      pred1 = "over_prob", pred2 = "under_prob",
      ev1 = "over_ev", ev2 = "under_ev",
      size1 = "over_bet_size", size2 = "under_bet_size",
      odds1 = "odds_over", odds2 = "odds_under",
      cents1 = "cents_over", cents2 = "cents_under",
      line_col_1 = "line", line_col_2 = "line",
      books = book_key,
      ev_threshold = ev_threshold
    )

  summary <- bets %>%
    group_by(market) %>%
    summarise(n_bets = n(), total_stake = sum(bet_size), avg_ev = mean(ev), max_ev = max(ev), .groups = "drop") %>%
    arrange(desc(total_stake))

  cat(sprintf("Generated %d %s totals bets\n", nrow(bets), book_key))

  list(predictions = predictions, prediction_set = prediction_set, bets = bets, markets_summary = summary)
}


#' Compare model predictions to Wagerzon odds for moneylines
#'
#' @param ml_results Output from build_moneylines_from_samples()
#' @param wagerzon_odds Wagerzon odds from get_wagerzon_odds()
#' @param bankroll Bankroll for Kelly sizing
#' @param kelly_mult Kelly multiplier (fraction of Kelly)
#' @param ev_threshold Minimum EV to include bet
#' @return List with bets and markets_summary
compare_moneylines_to_wagerzon <- function(
    ml_results,
    wagerzon_odds,
    bankroll = 100,
    kelly_mult = 0.25,
    ev_threshold = 0.05
) {
  predictions <- ml_results$predictions

  if (is.null(predictions) || nrow(predictions) == 0) {
    warning("No predictions found in ml_results")
    return(list(bets = data.frame()))
  }

  wz_ml <- wagerzon_odds %>%
    filter(market_type == "h2h")

  if (nrow(wz_ml) == 0) {
    warning("No moneylines found in Wagerzon odds")
    return(list(bets = data.frame()))
  }

  joined <- wz_ml %>%
    inner_join(
      predictions %>%
        select(id, home_team, away_team, market, home_win_prob, away_win_prob, commence_time),
      by = c("home_team", "away_team", "market")
    )

  if (nrow(joined) == 0) {
    cat("No matches found between predictions and Wagerzon moneylines.\n")
    return(list(bets = data.frame()))
  }

  cat(sprintf("Found %d matches between predictions and Wagerzon moneylines\n", nrow(joined)))

  prediction_set <- joined %>%
    mutate(as_tibble(american_prob(odds_away, odds_home))) %>%
    rename(wz_prob_away = p1, wz_prob_home = p2) %>%
    mutate(
      home_ev = compute_ev(home_win_prob, wz_prob_home),
      away_ev = compute_ev(away_win_prob, wz_prob_away),
      home_bet_size = kelly_stake(home_ev, wz_prob_home, bankroll, kelly_mult),
      away_bet_size = kelly_stake(away_ev, wz_prob_away, bankroll, kelly_mult)
    )

  # Get bookmaker key from input data
  book_key <- unique(wz_ml$bookmaker_key)[1]

  bets <- prediction_set %>%
    format_bets_table(
      pred1 = "home_win_prob", pred2 = "away_win_prob",
      ev1 = "home_ev", ev2 = "away_ev",
      size1 = "home_bet_size", size2 = "away_bet_size",
      odds1 = "odds_home", odds2 = "odds_away",
      cents1 = "cents_home", cents2 = "cents_away",
      books = book_key,
      ev_threshold = ev_threshold
    )

  summary <- bets %>%
    group_by(market) %>%
    summarise(n_bets = n(), total_stake = sum(bet_size), avg_ev = mean(ev), max_ev = max(ev), .groups = "drop") %>%
    arrange(desc(total_stake))

  cat(sprintf("Generated %d %s moneyline bets\n", nrow(bets), book_key))

  list(predictions = predictions, prediction_set = prediction_set, bets = bets, markets_summary = summary)
}


#' Compare Kalshi 3-way moneylines (home/away/tie) using 3-way predictions
#'
#' Kalshi 1H winner markets are 3-way: tie is a distinct losing outcome (not a push).
#' Uses devig_american_3way() to properly normalize all 3 outcomes, and
#' predict_moneyline_from_sample() which already computes _3way_home/_3way_away/_3way_tie.
compare_moneylines_3way_to_kalshi <- function(
    ml_results,
    kalshi_odds,
    samples,
    consensus_odds,
    bankroll = 100,
    kelly_mult = 0.25,
    ev_threshold = 0.05,
    margin_col = "game_home_margin_period"
) {
  # Filter to 3-way Kalshi records only
  kal_3way <- kalshi_odds %>% filter(market_type == "h2h_3way")

  if (nrow(kal_3way) == 0) {
    return(list(bets = tibble()))
  }

  # Build consensus lookup for team name → game_id mapping
  consensus_info <- consensus_odds %>%
    ungroup() %>%
    select(id, home_team, away_team, commence_time)
  if (is.character(consensus_info$commence_time)) {
    consensus_info <- consensus_info %>%
      mutate(commence_time = ymd_hms(commence_time, tz = "UTC"))
  }

  # Only generate predictions for games that have Kalshi 3-way odds
  needed_ids <- consensus_info %>%
    inner_join(kal_3way, by = c("home_team", "away_team")) %>%
    pull(id) %>%
    unique()

  needed_samples <- samples[intersect(names(samples), needed_ids)]
  if (length(needed_samples) == 0) {
    cat("No matching samples for Kalshi 3-way games.\n")
    return(list(bets = tibble()))
  }

  # Generate 3-way predictions only for matched games
  predictions_raw <- map_dfr(names(needed_samples), function(game_id) {
    sample_result <- needed_samples[[game_id]]
    preds <- predict_moneyline_from_sample(sample_result$sample, margin_col = margin_col)
    preds$id <- game_id
    preds
  })

  predictions <- predictions_raw %>%
    inner_join(consensus_info, by = "id")

  # Extract Half1 3-way predictions
  home_col <- paste0(margin_col, "_Half1_3way_home")
  away_col <- paste0(margin_col, "_Half1_3way_away")
  tie_col  <- paste0(margin_col, "_Half1_3way_tie")

  if (!home_col %in% names(predictions)) {
    cat("No Half1 3-way prediction columns found\n")
    return(list(bets = tibble()))
  }

  half1_preds <- predictions %>%
    transmute(
      id, home_team, away_team, commence_time,
      home_prob = .data[[home_col]],
      away_prob = .data[[away_col]],
      tie_prob  = .data[[tie_col]]
    )

  # Join with Kalshi odds
  joined <- kal_3way %>%
    inner_join(half1_preds, by = c("home_team", "away_team"))

  if (nrow(joined) == 0) {
    cat("No matches between 3-way predictions and Kalshi moneylines.\n")
    return(list(bets = tibble()))
  }

  cat(sprintf("Found %d matches between 3-way predictions and Kalshi moneylines\n", nrow(joined)))

  # Raw implied probs (with vig) — matches 2-way pattern
  prediction_set <- joined %>%
    mutate(bookmaker_key = "kalshi") %>%
    mutate(
      book_prob_home = odds_to_prob(odds_home),
      book_prob_away = odds_to_prob(odds_away),
      book_prob_tie  = odds_to_prob(odds_tie)
    ) %>%
    mutate(
      home_ev = compute_ev(home_prob, book_prob_home),
      away_ev = compute_ev(away_prob, book_prob_away),
      tie_ev  = compute_ev(tie_prob, book_prob_tie),
      home_bet_size = kelly_stake(home_ev, book_prob_home, bankroll, kelly_mult),
      away_bet_size = kelly_stake(away_ev, book_prob_away, bankroll, kelly_mult),
      tie_bet_size  = kelly_stake(tie_ev, book_prob_tie, bankroll, kelly_mult)
    )

  # Format using format_bets_table with 3rd side
  bets <- prediction_set %>%
    format_bets_table(
      pred1 = "home_prob", pred2 = "away_prob",
      ev1 = "home_ev", ev2 = "away_ev",
      size1 = "home_bet_size", size2 = "away_bet_size",
      odds1 = "odds_home", odds2 = "odds_away",
      cents1 = "cents_home", cents2 = "cents_away",
      pred3 = "tie_prob", ev3 = "tie_ev", size3 = "tie_bet_size", odds3 = "odds_tie",
      cents3 = "cents_tie",
      books = "kalshi",
      ev_threshold = ev_threshold
    )

  cat(sprintf("Generated %d Kalshi 3-way moneyline bets\n", nrow(bets)))
  list(bets = bets)
}


#' Compare offshore alt lines directly against samples (no API dependency)
#'
#' For alt markets not covered by the Odds API (e.g., CBB alternate_spreads_h1),
#' compute probabilities from the sample distribution and calculate EV against
#' the offshore book's prices.
#'
#' @param samples Named list of samples (keyed by game_id from Odds API)
#' @param offshore_odds Data frame from get_bfa_odds() / get_wagerzon_odds() etc.
#' @param consensus_odds API consensus odds (needs id, home_team, away_team, commence_time)
#' @param bankroll Bankroll for Kelly sizing
#' @param kelly_mult Kelly multiplier (fraction of Kelly)
#' @param ev_threshold Minimum EV to include bet
#' @param margin_col Column prefix for margin data in samples
#' @param total_col Column prefix for total data in samples
#' @return Data frame of bets in standard format (same as format_bets_table output)
compare_alts_to_samples <- function(
    samples,
    offshore_odds,
    consensus_odds,
    bankroll = 100,
    kelly_mult = 0.25,
    ev_threshold = 0.05,
    margin_col = "game_home_margin_period",
    total_col = "game_total_period"
) {
  # Market suffix → sample period column suffix
  period_map <- c(h1 = "Half1", h2 = "Half2", q1 = "1", q2 = "2", q3 = "3", q4 = "4")

  # Filter to alt markets and team totals
  alt_odds <- offshore_odds %>%
    filter(grepl("^alternate_", market) | grepl("^team_totals_", market))

  if (nrow(alt_odds) == 0) return(tibble())

  # Match offshore games to API game_ids via team names
  game_lookup <- consensus_odds %>%
    select(id, home_team, away_team, commence_time) %>%
    distinct(home_team, away_team, .keep_all = TRUE)

  alt_odds <- alt_odds %>%
    inner_join(game_lookup, by = c("home_team", "away_team"))

  if (nrow(alt_odds) == 0) {
    cat("No alt line games matched to API games.\n")
    return(tibble())
  }

  book_key <- unique(alt_odds$bookmaker_key)[1]
  all_bets <- list()

  for (i in seq_len(nrow(alt_odds))) {
    row <- alt_odds[i, ]
    game_id <- row$id

    # Look up sample for this game
    if (!game_id %in% names(samples)) next
    sample_df <- samples[[game_id]]$sample
    if (is.null(sample_df) || nrow(sample_df) == 0) next

    # Extract period from market name (e.g., "alternate_spreads_h1" → "h1" → "Half1")
    suffix <- sub(".*_", "", row$market)
    period <- period_map[suffix]
    if (is.na(period)) next

    pt_start_time <- tryCatch(
      lubridate::with_tz(as.POSIXct(row$commence_time, tz = "UTC"), tzone = "America/Los_Angeles"),
      error = function(e) as.POSIXct(NA)
    )

    if (grepl("spread", row$market) && !is.na(row$home_spread)) {
      # Alt spread: compute P(home covers) from sample margins
      col_name <- paste0(margin_col, "_", period)
      if (!col_name %in% names(sample_df)) next
      margins <- sample_df[[col_name]]
      margins <- margins[!is.na(margins)]
      if (length(margins) == 0) next

      home_spread <- row$home_spread
      non_push <- margins[margins != -home_spread]
      if (length(non_push) == 0) next
      p_home_cover <- sum(non_push > -home_spread) / length(non_push)
      p_away_cover <- 1 - p_home_cover

      # Convert BFA odds to implied probabilities (devigged)
      probs <- american_prob(row$odds_away, row$odds_home)
      if (any(is.na(probs)) || any(probs == 0)) next

      home_ev <- compute_ev(p_home_cover, probs$p2)
      away_ev <- compute_ev(p_away_cover, probs$p1)
      home_size <- kelly_stake(home_ev, probs$p2, bankroll, kelly_mult)
      away_size <- kelly_stake(away_ev, probs$p1, bankroll, kelly_mult)

      if (home_ev >= ev_threshold) {
        all_bets[[length(all_bets) + 1]] <- tibble(
          id = game_id, home_team = row$home_team, away_team = row$away_team,
          pt_start_time = pt_start_time, bookmaker_key = book_key,
          market = row$market, bet_on = row$home_team,
          line = home_spread, bet_size = home_size, ev = home_ev,
          odds = row$odds_home, prob = p_home_cover
        )
      }
      if (away_ev >= ev_threshold) {
        all_bets[[length(all_bets) + 1]] <- tibble(
          id = game_id, home_team = row$home_team, away_team = row$away_team,
          pt_start_time = pt_start_time, bookmaker_key = book_key,
          market = row$market, bet_on = row$away_team,
          line = row$away_spread, bet_size = away_size, ev = away_ev,
          odds = row$odds_away, prob = p_away_cover
        )
      }

    } else if (grepl("total", row$market) && !grepl("team_totals_", row$market) && !is.na(row$line)) {
      # Alt total: compute P(over) from sample totals (skip team totals, handled below)
      col_name <- paste0(total_col, "_", period)
      if (!col_name %in% names(sample_df)) next
      total_vals <- sample_df[[col_name]]
      total_vals <- total_vals[!is.na(total_vals)]
      if (length(total_vals) == 0) next

      total_line <- row$line
      non_push <- total_vals[total_vals != total_line]
      if (length(non_push) == 0) next
      p_over <- sum(non_push > total_line) / length(non_push)
      p_under <- 1 - p_over

      probs <- american_prob(row$odds_over, row$odds_under)
      if (any(is.na(probs)) || any(probs == 0)) next

      over_ev <- compute_ev(p_over, probs$p1)
      under_ev <- compute_ev(p_under, probs$p2)
      over_size <- kelly_stake(over_ev, probs$p1, bankroll, kelly_mult)
      under_size <- kelly_stake(under_ev, probs$p2, bankroll, kelly_mult)

      if (over_ev >= ev_threshold) {
        all_bets[[length(all_bets) + 1]] <- tibble(
          id = game_id, home_team = row$home_team, away_team = row$away_team,
          pt_start_time = pt_start_time, bookmaker_key = book_key,
          market = row$market, bet_on = "Over",
          line = total_line, bet_size = over_size, ev = over_ev,
          odds = row$odds_over, prob = p_over
        )
      }
      if (under_ev >= ev_threshold) {
        all_bets[[length(all_bets) + 1]] <- tibble(
          id = game_id, home_team = row$home_team, away_team = row$away_team,
          pt_start_time = pt_start_time, bookmaker_key = book_key,
          market = row$market, bet_on = "Under",
          line = total_line, bet_size = under_size, ev = under_ev,
          odds = row$odds_under, prob = p_under
        )
      }

    } else if (grepl("team_totals_", row$market) && !is.na(row$line)) {
      # Team total: compute P(over) from team-specific score samples
      team_side <- ifelse(grepl("_home_", row$market), "home", "away")
      score_col <- paste0(team_side, "_score_period_", period)
      if (!score_col %in% names(sample_df)) next
      scores <- sample_df[[score_col]]
      scores <- scores[!is.na(scores)]
      if (length(scores) == 0) next

      total_line <- row$line
      non_push <- scores[scores != total_line]
      if (length(non_push) == 0) next
      p_over <- sum(non_push > total_line) / length(non_push)
      p_under <- 1 - p_over

      probs <- american_prob(row$odds_over, row$odds_under)
      if (any(is.na(probs)) || any(probs == 0)) next

      team_name <- ifelse(team_side == "home", row$home_team, row$away_team)

      over_ev <- compute_ev(p_over, probs$p1)
      under_ev <- compute_ev(p_under, probs$p2)
      over_size <- kelly_stake(over_ev, probs$p1, bankroll, kelly_mult)
      under_size <- kelly_stake(under_ev, probs$p2, bankroll, kelly_mult)

      if (over_ev >= ev_threshold) {
        all_bets[[length(all_bets) + 1]] <- tibble(
          id = game_id, home_team = row$home_team, away_team = row$away_team,
          pt_start_time = pt_start_time, bookmaker_key = book_key,
          market = row$market, bet_on = paste(team_name, "Over"),
          line = total_line, bet_size = over_size, ev = over_ev,
          odds = row$odds_over, prob = p_over
        )
      }
      if (under_ev >= ev_threshold) {
        all_bets[[length(all_bets) + 1]] <- tibble(
          id = game_id, home_team = row$home_team, away_team = row$away_team,
          pt_start_time = pt_start_time, bookmaker_key = book_key,
          market = row$market, bet_on = paste(team_name, "Under"),
          line = total_line, bet_size = under_size, ev = under_ev,
          odds = row$odds_under, prob = p_under
        )
      }
    }
  }

  if (length(all_bets) == 0) return(tibble())
  result <- bind_rows(all_bets) %>% arrange(desc(ev))
  cat(sprintf("Generated %d %s alt/team-total bets from samples\n", nrow(result), book_key))
  result
}


#' Run full Wagerzon comparison pipeline
#'
#' Scrapes Wagerzon, then compares predictions to Wagerzon odds
#'
#' @param spread_results Output from build_spreads_from_samples()
#' @param totals_results Output from build_totals_from_samples()
#' @param ml_results Output from build_moneylines_from_samples()
#' @param sport Sport to scrape ("nfl")
#' @param bankroll Bankroll for Kelly sizing
#' @param kelly_mult Kelly multiplier
#' @param ev_threshold Minimum EV to include
#' @param scrape_first Whether to run the scraper before comparing
#' @return List with combined bets and summaries
run_wagerzon_comparison <- function(
    spread_results = NULL,
    totals_results = NULL,
    ml_results = NULL,
    sport = "nfl",
    bankroll = 100,
    kelly_mult = 0.25,
    ev_threshold = 0.05,
    scrape_first = TRUE
) {
  if (scrape_first) {
    run_wagerzon_scraper(sport = sport)
  }

  wz_odds <- get_wagerzon_odds(sport = sport)

  if (nrow(wz_odds) == 0) {
    warning("No Wagerzon odds available")
    return(list(bets = data.frame()))
  }

  all_bets <- list()

  if (!is.null(spread_results)) {
    wz_spreads <- compare_spreads_to_wagerzon(spread_results, wz_odds, bankroll, kelly_mult, ev_threshold)
    if (nrow(wz_spreads$bets) > 0) {
      all_bets$spreads <- wz_spreads$bets %>% mutate(market_type = "spreads")
    }
  }

  if (!is.null(totals_results)) {
    wz_totals <- compare_totals_to_wagerzon(totals_results, wz_odds, bankroll, kelly_mult, ev_threshold)
    if (nrow(wz_totals$bets) > 0) {
      all_bets$totals <- wz_totals$bets %>% mutate(market_type = "totals")
    }
  }

  if (!is.null(ml_results)) {
    wz_ml <- compare_moneylines_to_wagerzon(ml_results, wz_odds, bankroll, kelly_mult, ev_threshold)
    if (nrow(wz_ml$bets) > 0) {
      all_bets$moneylines <- wz_ml$bets %>% mutate(market_type = "moneyline")
    }
  }

  combined_bets <- bind_rows(all_bets) %>%
    arrange(desc(ev))

  cat(sprintf("\n=== WAGERZON BETTING SUMMARY ===\n"))
  cat(sprintf("Total bets: %d\n", nrow(combined_bets)))

  if (nrow(combined_bets) > 0) {
    summary <- combined_bets %>%
      group_by(market_type) %>%
      summarise(n_bets = n(), total_stake = sum(bet_size), avg_ev = mean(ev), .groups = "drop")
    print(summary)
  }

  list(bets = combined_bets, wagerzon_odds = wz_odds)
}

# =============================================================================
# PARLAY FAIR ODDS CALCULATOR
# =============================================================================

#' Evaluate a single leg of a parlay against sample data
#'
#' @param samples Dataframe of sample rows (historical similar games)
#' @param leg List with: market, period, side, line
#' @return Logical vector: TRUE (win), FALSE (lose), NA (push)
#'
#' @examples
#' # 1H home -3
#' evaluate_leg(samples, list(market = "spread", period = "Half1", side = "home", line = -3))
#' # 1H under 22.5
#' evaluate_leg(samples, list(market = "total", period = "Half1", side = "under", line = 22.5))
evaluate_leg <- function(samples, leg) {

  period <- leg$period

  # Get the right columns based on period

if (period == "Full") {
    margin_col <- "home_margin"
    total_col <- "total_final_score"
  } else {
    margin_col <- paste0("game_home_margin_period_", period)
    total_col <- paste0("game_total_period_", period)
  }

  margin <- samples[[margin_col]]
  total <- samples[[total_col]]

  result <- switch(leg$market,

    "spread" = {
      if (leg$side == "home") {
        # Home -3 means line = -3, home covers if margin > 3
        threshold <- -leg$line
        case_when(
          margin > threshold ~ TRUE,
          margin < threshold ~ FALSE,
          TRUE ~ NA
        )
      } else {
        # Away +3 means line = 3, away covers if margin < 3
        threshold <- leg$line
        case_when(
          margin < threshold ~ TRUE,
          margin > threshold ~ FALSE,
          TRUE ~ NA
        )
      }
    },

    "total" = {
      if (leg$side == "over") {
        case_when(
          total > leg$line ~ TRUE,
          total < leg$line ~ FALSE,
          TRUE ~ NA
        )
      } else {
        case_when(
          total < leg$line ~ TRUE,
          total > leg$line ~ FALSE,
          TRUE ~ NA
        )
      }
    },

    "team_total" = {
      home_score <- (total + margin) / 2
      away_score <- (total - margin) / 2

      switch(leg$side,
        "home_over" = case_when(
          home_score > leg$line ~ TRUE,
          home_score < leg$line ~ FALSE,
          TRUE ~ NA
        ),
        "home_under" = case_when(
          home_score < leg$line ~ TRUE,
          home_score > leg$line ~ FALSE,
          TRUE ~ NA
        ),
        "away_over" = case_when(
          away_score > leg$line ~ TRUE,
          away_score < leg$line ~ FALSE,
          TRUE ~ NA
        ),
        "away_under" = case_when(
          away_score < leg$line ~ TRUE,
          away_score > leg$line ~ FALSE,
          TRUE ~ NA
        )
      )
    },

    "moneyline" = {
      # 2-way moneyline: tie = push
      if (leg$side == "home") {
        case_when(
          margin > 0 ~ TRUE,
          margin < 0 ~ FALSE,
          TRUE ~ NA  # tie = push
        )
      } else {
        case_when(
          margin < 0 ~ TRUE,
          margin > 0 ~ FALSE,
          TRUE ~ NA  # tie = push
        )
      }
    },

    "moneyline_3way" = {
      # 3-way moneyline: can bet on tie (no pushes)
      switch(leg$side,
        "home" = case_when(
          margin > 0 ~ TRUE,
          TRUE ~ FALSE  # lose if tie or away wins
        ),
        "away" = case_when(
          margin < 0 ~ TRUE,
          TRUE ~ FALSE  # lose if tie or home wins
        ),
        "tie" = case_when(
          margin == 0 ~ TRUE,
          TRUE ~ FALSE  # lose if either team wins
        )
      )
    }
  )

  return(result)
}


#' Compute fair parlay odds from sample data
#'
#' @param samples Dataframe of sample rows for a single game
#' @param legs List of leg specifications, each with: market, period, side, line
#' @return List with joint_prob, fair_american_odds, leg_probs, correlation_factor, etc.
#'
#' @examples
#' legs <- list(
#'   list(market = "spread", period = "Half1", side = "home", line = -3),
#'   list(market = "total", period = "Half1", side = "under", line = 22.5)
#' )
#' result <- compute_parlay_fair_odds(samples, legs)
compute_parlay_fair_odds <- function(samples, legs) {

  # Evaluate each leg
  leg_results <- map(legs, ~evaluate_leg(samples, .x))

  # Find samples where ALL legs resolved (no NAs)
  all_resolved <- map(leg_results, ~!is.na(.x)) %>%
    reduce(`&`)

  # Find samples where ALL legs hit
  all_hit <- map(leg_results, ~.x == TRUE) %>%
    reduce(`&`)
  all_hit[is.na(all_hit)] <- FALSE

  # Count resolved and hits
  n_samples_total <- nrow(samples)
  n_samples_resolved <- sum(all_resolved)
  n_hits <- sum(all_hit & all_resolved)
  n_pushes <- n_samples_total - n_samples_resolved

  # Joint probability
  joint_prob <- n_hits / n_samples_resolved

  # Individual leg probabilities (for correlation factor)
  leg_probs <- map_dbl(leg_results, function(res) {
    resolved <- !is.na(res)
    sum(res[resolved] == TRUE) / sum(resolved)
  })

  # Correlation factor: actual joint prob vs independent assumption
  independent_prob <- prod(leg_probs)
  correlation_factor <- joint_prob / independent_prob

  # Convert to odds
  fair_american_odds <- prob_to_american(joint_prob)
  fair_decimal_odds <- 1 / joint_prob

  list(
    joint_prob = joint_prob,
    fair_american_odds = fair_american_odds,
    fair_decimal_odds = round(fair_decimal_odds, 2),
    n_samples_resolved = n_samples_resolved,
    n_samples_total = n_samples_total,
    n_hits = n_hits,
    n_pushes = n_pushes,
    leg_probs = leg_probs,
    independent_prob = independent_prob,
    correlation_factor = round(correlation_factor, 3)
  )
}


#' Format a leg specification as a readable string
#'
#' @param leg List with: market, period, side, line
#' @return Character string describing the leg
format_leg <- function(leg) {
  period_display <- switch(leg$period,
    "1" = "1Q",
    "2" = "2Q",
    "3" = "3Q",
    "4" = "4Q",
    "Half1" = "1H",
    "Half2" = "2H",
    "Full" = "FG",
    leg$period
  )

  switch(leg$market,
    "spread" = sprintf("%s %s %+.1f", period_display, leg$side, leg$line),
    "total" = sprintf("%s %s %.1f", period_display, leg$side, leg$line),
    "team_total" = sprintf("%s %s %.1f", period_display, gsub("_", " ", leg$side), leg$line),
    "moneyline" = sprintf("%s %s ML", period_display, leg$side),
    "moneyline_3way" = sprintf("%s %s ML (3-way)", period_display, leg$side),
    sprintf("%s %s", leg$market, leg$side)
  )
}


#' Print parlay fair odds result in a readable format
#'
#' @param result Output from compute_parlay_fair_odds
#' @param legs The legs used (for display)
print_parlay_result <- function(result, legs) {
  cat("\n=== PARLAY FAIR ODDS ===\n\n")

  cat("Legs:\n")
  for (i in seq_along(legs)) {
    cat(sprintf("  %d. %s (%.1f%% individual)\n",
        i, format_leg(legs[[i]]), result$leg_probs[i] * 100))
  }

  cat(sprintf("\nJoint probability: %.2f%%\n", result$joint_prob * 100))
  cat(sprintf("Fair American odds: %+d\n", result$fair_american_odds))
  cat(sprintf("Fair decimal odds: %.2f\n", result$fair_decimal_odds))

  cat(sprintf("\nCorrelation factor: %.3f", result$correlation_factor))
  if (result$correlation_factor > 1) {
    cat(" (positive correlation - legs help each other)\n")
  } else if (result$correlation_factor < 1) {
    cat(" (negative correlation - legs hurt each other)\n")
  } else {
    cat(" (no correlation)\n")
  }

  cat(sprintf("\nSamples: %d resolved, %d pushes excluded\n",
      result$n_samples_resolved, result$n_pushes))
}


# =============================================================================
# CORRELATION-ADJUSTED KELLY SIZING
# =============================================================================

#' Convert a bet table row into an evaluate_leg() spec
#'
#' @param bet_row A single row from the bets table (as a list or 1-row data frame)
#' @return List with market, period, side, line (compatible with evaluate_leg())
bet_to_leg <- function(bet_row) {
  market_raw <- bet_row$market

  # Strip alternate_ prefix
  clean_market <- gsub("^alternate_", "", market_raw)

  # Parse period from suffix
  period <- if (grepl("_h1$", clean_market)) {
    "Half1"
  } else if (grepl("_h2$", clean_market)) {
    "Half2"
  } else {
    "Full"
  }

  # Parse market type
  market_base <- gsub("_(h1|h2)$", "", clean_market)
  leg_market <- switch(market_base,
    "h2h"               = "moneyline",
    "spreads"            = "spread",
    "totals"             = "total",
    "team_totals"        = "team_total",
    "team_totals_home"   = "team_total",
    "team_totals_away"   = "team_total",
    stop(paste("Unknown market base:", market_base))
  )

  # Parse side from bet_on
  bet_on <- bet_row$bet_on
  home_team <- bet_row$home_team
  away_team <- bet_row$away_team

  if (leg_market == "team_total") {
    # bet_on is like "Duke Over" or "Duke Under"
    is_over <- grepl(" Over$", bet_on)
    team_name <- sub(" (Over|Under)$", "", bet_on)
    team_side <- if (team_name == home_team) "home" else "away"
    side <- paste0(team_side, "_", ifelse(is_over, "over", "under"))
  } else if (leg_market %in% c("spread", "moneyline")) {
    side <- if (bet_on == home_team) "home" else "away"
  } else {
    # totals: bet_on is "Over" or "Under"
    side <- tolower(bet_on)
  }

  list(
    market = leg_market,
    period = period,
    side   = side,
    line   = bet_row$line
  )
}


#' Compute portfolio-optimal Kelly fractions for correlated bets on one game
#'
#' Uses the multivariate Kelly approximation: f* = Σ⁻¹ · μ
#' where Σ is the covariance matrix of bet returns and μ is the EV vector.
#'
#' @param bets_group Data frame of bets on the same game (2+ rows)
#' @param sample Data frame of historical sample games for this game
#' @return Named list with f_star (optimal fractions) and cor_matrix, or NULL if ill-conditioned
multivariate_kelly <- function(bets_group, sample) {
  n <- nrow(bets_group)

  # Convert each bet to a leg and evaluate on the sample
  outcome_matrix <- matrix(NA_real_, nrow = nrow(sample), ncol = n)
  for (i in seq_len(n)) {
    leg <- tryCatch(bet_to_leg(bets_group[i, ]), error = function(e) NULL)
    if (is.null(leg)) return(NULL)
    outcomes <- tryCatch(evaluate_leg(sample, leg), error = function(e) NULL)
    if (is.null(outcomes) || length(outcomes) == 0) return(NULL)
    outcome_matrix[, i] <- as.numeric(outcomes)
  }

  # Remove rows with any NA (pushes on any leg)
  complete_rows <- complete.cases(outcome_matrix)
  if (sum(complete_rows) < 30) return(NULL)  # too few resolved samples
  outcome_matrix <- outcome_matrix[complete_rows, , drop = FALSE]

  # Compute Pearson correlation matrix of binary outcomes
  R <- cor(outcome_matrix)
  if (any(is.na(R))) return(NULL)

  # Build return covariance matrix
  # For binary bet: σ_i = sqrt(p_i * (1-p_i)) * (b_i + 1)
  # where b_i = decimal_odds - 1, p_i = predicted probability
  pred_probs <- bets_group$prob
  odds_american <- bets_group$odds
  decimal_odds <- ifelse(odds_american > 0, 1 + odds_american / 100, 1 + 100 / abs(odds_american))
  b <- decimal_odds - 1

  sigmas <- sqrt(pred_probs * (1 - pred_probs)) * (b + 1)

  # Covariance matrix: Σ_ij = R_ij * σ_i * σ_j
  Sigma <- R * outer(sigmas, sigmas)

  # Ridge regularization: prevents singularity from near-duplicate bets
  # (e.g., Under +75.5 placed vs Under +75 pipeline have ρ ≈ 0.99)
  Sigma <- Sigma + diag(ncol(Sigma)) * 0.01

  # Check condition number (rarely triggers after regularization)
  if (kappa(Sigma) > 100) return(NULL)

  # EV vector (already computed in the bet table)
  mu <- bets_group$ev

  # Solve: f* = Σ⁻¹ · μ
  f_star <- tryCatch(
    solve(Sigma, mu),
    error = function(e) NULL
  )
  if (is.null(f_star)) return(NULL)

  # Clamp negatives to 0
  f_star <- pmax(f_star, 0)

  list(
    f_star     = f_star,
    cor_matrix = R,
    cov_matrix = Sigma
  )
}


#' Adjust Kelly bet sizes for within-game correlation
#'
#' Groups bets by game, measures empirical correlations from the game's sample,
#' and adjusts sizes using multivariate Kelly (primary) or per-bet average ρ (fallback).
#'
#' @param bets_df Data frame of bets (output of Phase 7 dedup)
#' @param samples Named list of sample results, keyed by game id
#' @param bankroll Current bankroll
#' @param kelly_mult Fractional Kelly multiplier
#' @return bets_df with adjusted bet_size and new correlation_adj column
adjust_kelly_for_correlation <- function(bets_df, samples, bankroll, kelly_mult,
                                         placed_bets = NULL) {
  if (nrow(bets_df) == 0) return(mutate(bets_df, correlation_adj = numeric(0)))

  bets_df$correlation_adj <- 1.0

  # Normalize placed_bets to match bets_df column names
  has_placed <- !is.null(placed_bets) && nrow(placed_bets) > 0
  if (has_placed) {
    # Rename game_id -> id to match bets_df
    if ("game_id" %in% names(placed_bets) && !"id" %in% names(placed_bets)) {
      placed_bets$id <- placed_bets$game_id
    }
    placed_bets$.is_placed <- TRUE
    n_placed_games <- length(unique(placed_bets$id[placed_bets$id %in% unique(bets_df$id)]))
    cat(sprintf("Placed bets loaded: %d bets on %d games overlapping with new bets\n",
                nrow(placed_bets), n_placed_games))
  }

  game_ids <- unique(bets_df$id)
  n_mv <- 0L
  n_fb <- 0L
  n_skip <- 0L
  n_placed_used <- 0L

  for (gid in game_ids) {
    idx <- which(bets_df$id == gid)
    new_bets <- bets_df[idx, ]
    new_bets$.is_placed <- FALSE

    # Check for already-placed bets on this game
    placed_game <- NULL
    if (has_placed) {
      placed_game <- placed_bets[placed_bets$id == gid, , drop = FALSE]
      if (nrow(placed_game) == 0) placed_game <- NULL
    }

    # When a placed bet matches a pipeline bet (same market/bet_on/line),
    # keep the placed bet as a fixed position and remove the matching
    # pipeline bet (already acted on). This ensures other correlated
    # pipeline bets are penalized for actual placed exposure.
    if (!is.null(placed_game)) {
      placed_keys <- paste(placed_game$market, placed_game$bet_on,
                           ifelse(is.na(placed_game$line), "ML", placed_game$line))
      new_keys <- paste(new_bets$market, new_bets$bet_on,
                        ifelse(is.na(new_bets$line), "ML", new_bets$line))

      # Remove pipeline bets that have already been placed
      pipeline_dup_mask <- new_keys %in% placed_keys
      if (any(pipeline_dup_mask)) {
        # Zero out matched pipeline bets in bets_df so they don't appear in output
        removed_idx <- idx[pipeline_dup_mask]
        bets_df$bet_size[removed_idx] <- 0
        bets_df$correlation_adj[removed_idx] <- 0

        idx <- idx[!pipeline_dup_mask]
        new_bets <- new_bets[!pipeline_dup_mask, , drop = FALSE]
      }

      # No new bets left to size — all were already placed
      if (nrow(new_bets) == 0) next
    }

    # Build full correlation group: new bets + placed bets
    n_total <- nrow(new_bets) + ifelse(is.null(placed_game), 0L, nrow(placed_game))
    if (n_total < 2) next  # need 2+ bets in the group for correlation

    # Look up game sample
    if (is.null(samples[[gid]])) {
      n_skip <- n_skip + 1L
      next
    }
    game_sample <- samples[[gid]]$sample

    # Combine placed + new for correlation computation
    if (!is.null(placed_game)) {
      # Ensure placed_game has the columns multivariate_kelly needs
      bets_group <- bind_rows(placed_game, new_bets)
      n_placed_used <- n_placed_used + 1L
    } else {
      bets_group <- new_bets
    }

    n_new <- nrow(new_bets)
    n_placed_in_group <- nrow(bets_group) - n_new
    new_positions <- (n_placed_in_group + 1):nrow(bets_group)  # indices of new bets in bets_group

    mv_result <- multivariate_kelly(bets_group, game_sample)

    kelly_applied <- FALSE

    if (!is.null(mv_result)) {
      original_sizes <- new_bets$bet_size

      if (n_placed_in_group > 0) {
        # Conditional Kelly: account for existing placed exposure
        # f_new* = Σ_nn⁻¹ × (μ_new − Σ_np × f_placed)
        # Penalty term reduces recommendations proportionally to how much
        # is already wagered on correlated positions
        p_idx <- 1:n_placed_in_group
        n_idx <- new_positions
        Sigma <- mv_result$cov_matrix
        Sigma_nn <- Sigma[n_idx, n_idx, drop = FALSE]
        Sigma_np <- Sigma[n_idx, p_idx, drop = FALSE]

        f_placed <- bets_group$bet_size[p_idx] / (kelly_mult * bankroll)
        mu_new <- bets_group$ev[n_idx]
        adjusted_mu <- as.numeric(mu_new - Sigma_np %*% f_placed)

        f_new_star <- tryCatch(solve(Sigma_nn, adjusted_mu), error = function(e) NULL)
        if (!is.null(f_new_star)) {
          mv_sizes <- pmax(as.numeric(f_new_star), 0) * kelly_mult * bankroll
          mv_sizes <- round(mv_sizes, 2)
          mv_sizes <- pmin(mv_sizes, original_sizes * 1.5)

          bets_df$bet_size[idx] <- mv_sizes
          bets_df$correlation_adj[idx] <- ifelse(original_sizes > 0, mv_sizes / original_sizes, 1.0)
          n_mv <- n_mv + 1L
          kelly_applied <- TRUE
        }
      } else {
        # Standard multivariate Kelly (no placed bets in group)
        mv_fracs <- mv_result$f_star[new_positions]
        mv_sizes <- mv_fracs * kelly_mult * bankroll
        mv_sizes <- round(mv_sizes, 2)
        mv_sizes <- pmin(mv_sizes, original_sizes * 1.5)
        mv_sizes <- pmax(mv_sizes, 0)

        bets_df$bet_size[idx] <- mv_sizes
        bets_df$correlation_adj[idx] <- ifelse(original_sizes > 0, mv_sizes / original_sizes, 1.0)
        n_mv <- n_mv + 1L
        kelly_applied <- TRUE
      }
    }

    if (!kelly_applied) {
      # Fallback: per-bet average ρ on the full group
      n_all <- nrow(bets_group)
      outcome_matrix <- matrix(NA_real_, nrow = nrow(game_sample), ncol = n_all)
      leg_failed <- FALSE
      for (i in seq_len(n_all)) {
        leg <- tryCatch(bet_to_leg(bets_group[i, ]), error = function(e) NULL)
        if (is.null(leg)) { leg_failed <- TRUE; break }
        outcomes <- tryCatch(evaluate_leg(game_sample, leg), error = function(e) NULL)
        if (is.null(outcomes) || length(outcomes) == 0) { leg_failed <- TRUE; break }
        outcome_matrix[, i] <- as.numeric(outcomes)
      }
      if (leg_failed) { n_skip <- n_skip + 1L; next }

      complete_rows <- complete.cases(outcome_matrix)
      if (sum(complete_rows) < 30) {
        n_skip <- n_skip + 1L
        next
      }
      outcome_matrix <- outcome_matrix[complete_rows, , drop = FALSE]
      R <- cor(outcome_matrix)
      if (any(is.na(R))) {
        n_skip <- n_skip + 1L
        next
      }

      # Only adjust new bets (positions after placed bets in the group)
      p_idx <- if (n_placed_in_group > 0) 1:n_placed_in_group else integer(0)
      for (j in seq_len(n_new)) {
        i <- new_positions[j]  # position in full group

        # If any placed bet is highly correlated (ρ > 0.90), set to $0
        # (fallback can't properly compute conditional Kelly subtraction)
        if (n_placed_in_group > 0 && max(R[i, p_idx]) > 0.90) {
          bets_df$bet_size[idx[j]] <- 0
          bets_df$correlation_adj[idx[j]] <- 0
          next
        }

        avg_rho_i <- mean(R[i, -i])
        denom <- 1 + (n_all - 1) * avg_rho_i
        scale_i <- if (denom <= 0) 1.5 else min(1.5, 1 / sqrt(denom))
        scale_i <- max(scale_i, 0)
        bets_df$bet_size[idx[j]] <- round(bets_df$bet_size[idx[j]] * scale_i, 2)
        bets_df$correlation_adj[idx[j]] <- scale_i
      }
      n_fb <- n_fb + 1L
    }
  }

  cat(sprintf("\n=== CORRELATION ADJUSTMENT ===\n"))
  cat(sprintf("Games with 2+ bets: %d (multivariate: %d, fallback: %d, skipped: %d)\n",
              n_mv + n_fb + n_skip, n_mv, n_fb, n_skip))
  if (n_placed_used > 0) {
    cat(sprintf("Games using placed bets for correlation: %d\n", n_placed_used))
  }
  adjusted <- bets_df$correlation_adj[bets_df$correlation_adj != 1.0]
  if (length(adjusted) > 0) {
    cat(sprintf("Adjusted bets: %d, avg scale: %.3f, range: [%.3f, %.3f]\n",
                length(adjusted), mean(adjusted), min(adjusted), max(adjusted)))
  }

  bets_df
}


# =============================================================================
# GENERALIZED DERIVATIVE BACKTEST FUNCTION
# =============================================================================

#' Run a generalized derivative backtest for any sport
#'
#' This function runs a leave-one-out backtest on derivative markets (H1, H2, etc.)
#' comparing answer key predictions to book implied probabilities and actual outcomes.
#'
#' @param pbp_data Data frame with play-by-play outcomes (margins, totals by period)
#' @param sport_config List with sport-specific configuration:
#'   - sport: "nfl", "cbb", etc.
#'   - margin_col_prefix: e.g., "game_home_margin_"
#'   - total_col_prefix: e.g., "game_total_"
#'   - periods: list with period names and their column suffixes
#'   - parent_spread_col: column name for FG spread
#'   - parent_total_col: column name for FG total
#'   - consensus_home_odds_col: column for consensus home cover prob
#'   - consensus_over_odds_col: column for consensus over prob
#'   - home_score_prefix: e.g., "home_" for home_h1_score
#'   - away_score_prefix: e.g., "away_" for away_h1_score
#' @param derivative_odds Optional data frame with derivative odds (for log-loss vs book)
#' @param kelly_mult Kelly multiplier for ROI simulation (default 0.25)
#' @param bankroll Starting bankroll for ROI simulation (default 1000)
#' @param sample_pct Percentage of historical data to use for sample size (default 0.10)
#' @param test_game_ids Optional vector of game IDs to test (NULL = test all games)
#' @param verbose Print progress messages (default TRUE)
#' @return List with: predictions, by_market, by_period, overall, plots
run_derivative_backtest <- function(
    pbp_data,
    sport_config,
    derivative_odds = NULL,
    kelly_mult = 0.25,
    bankroll = 1000,
    sample_pct = 0.10,
    test_game_ids = NULL,
    verbose = TRUE
) {

  # ============================================================================
  # HELPER FUNCTIONS
  # ============================================================================

  american_to_prob <- function(odds) {
    ifelse(odds > 0, 100 / (odds + 100), abs(odds) / (abs(odds) + 100))
  }

  devig_american_pair <- function(odds1, odds2) {
    p1 <- american_to_prob(odds1)
    p2 <- american_to_prob(odds2)
    total <- p1 + p2
    list(prob1 = p1 / total, prob2 = p2 / total)
  }

  calc_ev <- function(pred_prob, book_odds) {
    decimal_odds <- ifelse(book_odds > 0, 1 + book_odds / 100, 1 + 100 / abs(book_odds))
    pred_prob * decimal_odds - 1
  }

  kelly_stake <- function(ev, book_odds, bankroll, kelly_mult = 0.25) {
    if (ev <= 0) return(0)
    decimal_odds <- ifelse(book_odds > 0, 1 + book_odds / 100, 1 + 100 / abs(book_odds))
    b <- decimal_odds - 1
    p <- (ev + 1) / decimal_odds
    q <- 1 - p
    kelly_full <- max(0, (b * p - q) / b)
    bankroll * kelly_mult * kelly_full
  }

  # ============================================================================
  # SETUP
  # ============================================================================

  if (verbose) cat("\n===========================================\n")
  if (verbose) cat("DERIVATIVE BACKTEST:", toupper(sport_config$sport), "\n")
  if (verbose) cat("===========================================\n\n")

  # Convert to data.table for efficiency
  DT <- as.data.table(pbp_data)

  # Calculate dispersion for sample weighting
  disp <- compute_dispersion(DT, moneyline = FALSE,
                             spread_col = sport_config$parent_spread_col,
                             total_col = sport_config$parent_total_col)
  ss <- disp$ss
  st <- disp$st
  N <- round(nrow(DT) * sample_pct, 0)

  if (verbose) cat("Total games:", nrow(DT), "\n")
  if (verbose) cat("Sample size N:", N, "\n")
  if (verbose) cat("Dispersion - ss:", round(ss, 2), "st:", round(st, 2), "\n\n")

  # Get list of games to test (use subset if provided)
  if (!is.null(test_game_ids)) {
    all_game_ids <- test_game_ids[test_game_ids %in% unique(DT$game_id)]
    if (verbose) cat("Testing subset of", length(all_game_ids), "games\n\n")
  } else {
    all_game_ids <- unique(DT$game_id)
  }

  # ============================================================================
  # DEFINE MARKETS TO TEST
  # ============================================================================

  # Build market definitions based on sport config
  markets <- list()

  for (period_name in names(sport_config$periods)) {
    period_info <- sport_config$periods[[period_name]]
    suffix <- period_info$suffix

    # Spread market
    margin_col <- paste0(sport_config$margin_col_prefix, suffix)
    if (margin_col %in% names(DT)) {
      markets[[paste0("spread_", period_name)]] <- list(
        type = "spread",
        period = period_name,
        margin_col = margin_col,
        display_name = paste0(period_info$display, " Spread")
      )
    }

    # Total market
    total_col <- paste0(sport_config$total_col_prefix, suffix)
    if (total_col %in% names(DT)) {
      markets[[paste0("total_", period_name)]] <- list(
        type = "total",
        period = period_name,
        total_col = total_col,
        display_name = paste0(period_info$display, " Total")
      )
    }

    # Moneyline market (uses margin column)
    if (margin_col %in% names(DT)) {
      markets[[paste0("ml_", period_name)]] <- list(
        type = "moneyline",
        period = period_name,
        margin_col = margin_col,
        display_name = paste0(period_info$display, " Moneyline")
      )
    }

    # Team totals (home and away)
    home_score_col <- paste0(sport_config$home_score_prefix, suffix, "_score")
    away_score_col <- paste0(sport_config$away_score_prefix, suffix, "_score")

    if (home_score_col %in% names(DT)) {
      markets[[paste0("home_tt_", period_name)]] <- list(
        type = "team_total",
        period = period_name,
        score_col = home_score_col,
        side = "home",
        display_name = paste0(period_info$display, " Home Team Total")
      )
    }

    if (away_score_col %in% names(DT)) {
      markets[[paste0("away_tt_", period_name)]] <- list(
        type = "team_total",
        period = period_name,
        score_col = away_score_col,
        side = "away",
        display_name = paste0(period_info$display, " Away Team Total")
      )
    }
  }

  if (verbose) cat("Markets to test:", length(markets), "\n")
  if (verbose) cat("Market names:", paste(names(markets), collapse = ", "), "\n\n")

  # ============================================================================
  # RUN PREDICTIONS
  # ============================================================================

  if (verbose) cat("Running predictions...\n")

  all_predictions <- list()

  for (i in seq_along(all_game_ids)) {
    if (verbose && i %% 100 == 0) cat("  ", i, "/", length(all_game_ids), "\n")

    test_game_id <- all_game_ids[i]
    test_game <- DT[DT$game_id == test_game_id, ][1, ]

    # Get parent lines for this game
    parent_spread <- test_game[[sport_config$parent_spread_col]]
    parent_total <- test_game[[sport_config$parent_total_col]]

    # Get consensus probabilities for target cover/over rates
    target_cover <- test_game[[sport_config$consensus_home_odds_col]]
    target_over <- test_game[[sport_config$consensus_over_odds_col]]

    # Skip if missing parent data
    if (is.na(parent_spread) || is.na(parent_total)) next
    if (is.na(target_cover) || is.na(target_over)) next

    # Leave-one-out training data
    DT_train <- DT[DT$game_id != test_game_id, ]

    # Run answer key sampling
    sample_result <- tryCatch({
      run_answer_key_sample(
        id = test_game_id,
        parent_spread = parent_spread,
        parent_total = parent_total,
        target_cover = target_cover,
        target_over = target_over,
        DT = DT_train,
        ss = ss, st = st, N = N,
        use_spread_line = TRUE
      )
    }, error = function(e) NULL)

    if (is.null(sample_result)) next

    sample <- sample_result$sample

    # Generate predictions for each market
    for (market_name in names(markets)) {
      market <- markets[[market_name]]

      if (market$type == "spread") {
        # Spread: predict P(home covers) at line 0 (since we're testing the prediction accuracy)
        margin_col <- market$margin_col
        if (!(margin_col %in% names(sample))) next

        margins <- sample[[margin_col]]
        actual_margin <- test_game[[margin_col]]
        if (is.na(actual_margin)) next

        # Predict home cover prob (margin > 0)
        non_push <- margins != 0
        if (sum(non_push) == 0) next

        home_win_prob <- mean(margins[non_push] > 0)
        away_win_prob <- 1 - home_win_prob

        # Actual outcome (at 0 line = moneyline result for this period)
        if (actual_margin == 0) next
        actual_home_win <- ifelse(actual_margin > 0, 1, 0)

        # Home side
        all_predictions[[length(all_predictions) + 1]] <- list(
          game_id = test_game_id,
          market = market_name,
          market_type = "spread",
          period = market$period,
          display_name = market$display_name,
          bet_side = "home",
          algo_prob = home_win_prob,
          actual_outcome = actual_home_win,
          has_book_odds = FALSE,
          book_prob = NA_real_,
          book_odds = NA_real_,
          ev = NA_real_
        )

        # Away side
        all_predictions[[length(all_predictions) + 1]] <- list(
          game_id = test_game_id,
          market = market_name,
          market_type = "spread",
          period = market$period,
          display_name = market$display_name,
          bet_side = "away",
          algo_prob = away_win_prob,
          actual_outcome = 1 - actual_home_win,
          has_book_odds = FALSE,
          book_prob = NA_real_,
          book_odds = NA_real_,
          ev = NA_real_
        )

      } else if (market$type == "total") {
        # Total: predict P(over) at the average total for this period
        total_col <- market$total_col
        if (!(total_col %in% names(sample))) next

        totals <- sample[[total_col]]
        actual_total <- test_game[[total_col]]
        if (is.na(actual_total)) next

        # Use median as reference line
        ref_line <- median(totals, na.rm = TRUE)

        non_push <- totals != ref_line
        if (sum(non_push) == 0) next

        over_prob <- mean(totals[non_push] > ref_line)
        under_prob <- 1 - over_prob

        # Actual outcome
        if (actual_total == ref_line) next
        actual_over <- ifelse(actual_total > ref_line, 1, 0)

        # Over side
        all_predictions[[length(all_predictions) + 1]] <- list(
          game_id = test_game_id,
          market = market_name,
          market_type = "total",
          period = market$period,
          display_name = market$display_name,
          bet_side = "over",
          algo_prob = over_prob,
          actual_outcome = actual_over,
          has_book_odds = FALSE,
          book_prob = NA_real_,
          book_odds = NA_real_,
          ev = NA_real_
        )

        # Under side
        all_predictions[[length(all_predictions) + 1]] <- list(
          game_id = test_game_id,
          market = market_name,
          market_type = "total",
          period = market$period,
          display_name = market$display_name,
          bet_side = "under",
          algo_prob = under_prob,
          actual_outcome = 1 - actual_over,
          has_book_odds = FALSE,
          book_prob = NA_real_,
          book_odds = NA_real_,
          ev = NA_real_
        )

      } else if (market$type == "moneyline") {
        # Moneyline: predict P(home wins)
        margin_col <- market$margin_col
        if (!(margin_col %in% names(sample))) next

        margins <- sample[[margin_col]]
        actual_margin <- test_game[[margin_col]]
        if (is.na(actual_margin)) next

        # 2-way moneyline (excluding ties)
        non_tie <- margins != 0
        if (sum(non_tie) == 0) next

        home_win_prob <- mean(margins[non_tie] > 0)
        away_win_prob <- 1 - home_win_prob

        # Skip ties in actual outcome
        if (actual_margin == 0) next
        actual_home_win <- ifelse(actual_margin > 0, 1, 0)

        # Home side
        all_predictions[[length(all_predictions) + 1]] <- list(
          game_id = test_game_id,
          market = market_name,
          market_type = "moneyline",
          period = market$period,
          display_name = market$display_name,
          bet_side = "home",
          algo_prob = home_win_prob,
          actual_outcome = actual_home_win,
          has_book_odds = FALSE,
          book_prob = NA_real_,
          book_odds = NA_real_,
          ev = NA_real_
        )

        # Away side
        all_predictions[[length(all_predictions) + 1]] <- list(
          game_id = test_game_id,
          market = market_name,
          market_type = "moneyline",
          period = market$period,
          display_name = market$display_name,
          bet_side = "away",
          algo_prob = away_win_prob,
          actual_outcome = 1 - actual_home_win,
          has_book_odds = FALSE,
          book_prob = NA_real_,
          book_odds = NA_real_,
          ev = NA_real_
        )

      } else if (market$type == "team_total") {
        # Team total: predict P(team scores > median)
        score_col <- market$score_col
        if (!(score_col %in% names(sample))) next

        scores <- sample[[score_col]]
        actual_score <- test_game[[score_col]]
        if (is.na(actual_score)) next

        # Use median as reference line
        ref_line <- median(scores, na.rm = TRUE)

        non_push <- scores != ref_line
        if (sum(non_push) == 0) next

        over_prob <- mean(scores[non_push] > ref_line)
        under_prob <- 1 - over_prob

        # Actual outcome
        if (actual_score == ref_line) next
        actual_over <- ifelse(actual_score > ref_line, 1, 0)

        # Over side
        all_predictions[[length(all_predictions) + 1]] <- list(
          game_id = test_game_id,
          market = market_name,
          market_type = "team_total",
          period = market$period,
          display_name = market$display_name,
          bet_side = paste0(market$side, "_over"),
          algo_prob = over_prob,
          actual_outcome = actual_over,
          has_book_odds = FALSE,
          book_prob = NA_real_,
          book_odds = NA_real_,
          ev = NA_real_
        )

        # Under side
        all_predictions[[length(all_predictions) + 1]] <- list(
          game_id = test_game_id,
          market = market_name,
          market_type = "team_total",
          period = market$period,
          display_name = market$display_name,
          bet_side = paste0(market$side, "_under"),
          algo_prob = under_prob,
          actual_outcome = 1 - actual_over,
          has_book_odds = FALSE,
          book_prob = NA_real_,
          book_odds = NA_real_,
          ev = NA_real_
        )
      }
    }
  }

  # Convert to data frame
  results <- map_dfr(all_predictions, as_tibble)

  if (nrow(results) == 0) {
    warning("No predictions generated!")
    return(NULL)
  }

  if (verbose) cat("\nTotal predictions:", nrow(results), "\n")
  if (verbose) cat("Unique games:", n_distinct(results$game_id), "\n\n")

  # ============================================================================
  # CALCULATE METRICS
  # ============================================================================

  # Log-loss for algorithm predictions
  results <- results %>%
    mutate(
      algo_logloss = -(actual_outcome * log(pmax(algo_prob, 1e-10)) +
                       (1 - actual_outcome) * log(pmax(1 - algo_prob, 1e-10)))
    )

  # ============================================================================
  # RESULTS BY MARKET
  # ============================================================================

  by_market <- results %>%
    group_by(market, market_type, period, display_name) %>%
    summarize(
      n_predictions = n(),
      n_games = n_distinct(game_id),
      algo_logloss = mean(algo_logloss, na.rm = TRUE),
      avg_algo_prob = mean(algo_prob),
      avg_actual = mean(actual_outcome),
      calibration_error = mean(algo_prob) - mean(actual_outcome),
      abs_cal_error = abs(mean(algo_prob) - mean(actual_outcome)),
      .groups = "drop"
    ) %>%
    arrange(market_type, period)

  # ============================================================================
  # RESULTS BY PERIOD
  # ============================================================================

  by_period <- results %>%
    group_by(period) %>%
    summarize(
      n_predictions = n(),
      n_games = n_distinct(game_id),
      algo_logloss = mean(algo_logloss, na.rm = TRUE),
      avg_algo_prob = mean(algo_prob),
      avg_actual = mean(actual_outcome),
      calibration_error = mean(algo_prob) - mean(actual_outcome),
      .groups = "drop"
    )

  # ============================================================================
  # OVERALL SUMMARY
  # ============================================================================

  overall <- results %>%
    summarize(
      total_predictions = n(),
      total_games = n_distinct(game_id),
      algo_logloss = mean(algo_logloss, na.rm = TRUE),
      avg_algo_prob = mean(algo_prob),
      avg_actual = mean(actual_outcome),
      calibration_error = mean(algo_prob) - mean(actual_outcome)
    )

  # ============================================================================
  # CALIBRATION ANALYSIS
  # ============================================================================

  # Bin predictions into deciles for calibration plot
  calibration <- results %>%
    mutate(prob_bin = cut(algo_prob, breaks = seq(0, 1, 0.1), include.lowest = TRUE)) %>%
    group_by(prob_bin) %>%
    summarize(
      n = n(),
      avg_predicted = mean(algo_prob),
      avg_actual = mean(actual_outcome),
      .groups = "drop"
    ) %>%
    filter(!is.na(prob_bin))

  # ============================================================================
  # PRINT RESULTS
  # ============================================================================

  if (verbose) {
    cat("===========================================\n")
    cat("RESULTS BY MARKET\n")
    cat("===========================================\n\n")
    print(by_market %>% select(display_name, n_predictions, n_games, algo_logloss, calibration_error))

    cat("\n===========================================\n")
    cat("RESULTS BY PERIOD\n")
    cat("===========================================\n\n")
    print(by_period)

    cat("\n===========================================\n")
    cat("OVERALL SUMMARY\n")
    cat("===========================================\n\n")
    cat("Total predictions:", overall$total_predictions, "\n")
    cat("Total games:", overall$total_games, "\n")
    cat("Algorithm log-loss:", round(overall$algo_logloss, 4), "\n")
    cat("Calibration error:", round(overall$calibration_error, 4), "\n")

    cat("\n===========================================\n")
    cat("CALIBRATION TABLE\n")
    cat("===========================================\n\n")
    print(calibration)
  }

  # ============================================================================
  # RETURN RESULTS
  # ============================================================================

  list(
    predictions = results,
    by_market = by_market,
    by_period = by_period,
    overall = overall,
    calibration = calibration,
    config = sport_config
  )
}


#' Create sport-specific configuration for derivative backtest
#'
#' @param sport One of "nfl", "cbb"
#' @return List with sport-specific column mappings and period definitions
get_sport_backtest_config <- function(sport) {

  if (sport == "cbb") {
    list(
      sport = "cbb",
      margin_col_prefix = "game_home_margin_",
      total_col_prefix = "game_total_",
      home_score_prefix = "home_",
      away_score_prefix = "away_",
      parent_spread_col = "home_spread",
      parent_total_col = "total_line",
      consensus_home_odds_col = "consensus_home_cover_prob",
      consensus_over_odds_col = "consensus_over_prob",
      periods = list(
        h1 = list(suffix = "h1", display = "H1"),
        h2 = list(suffix = "h2", display = "H2")
      )
    )

  } else if (sport == "nfl") {
    list(
      sport = "nfl",
      margin_col_prefix = "game_home_margin_period_",
      total_col_prefix = "game_total_period_",
      home_score_prefix = "home_score_period_",
      away_score_prefix = "away_score_period_",
      parent_spread_col = "home_spread",
      parent_total_col = "total_line",
      consensus_home_odds_col = "consensus_devig_home_odds",
      consensus_over_odds_col = "consensus_devig_over_odds",
      periods = list(
        q1 = list(suffix = "1", display = "Q1"),
        q2 = list(suffix = "2", display = "Q2"),
        q3 = list(suffix = "3", display = "Q3"),
        q4 = list(suffix = "4", display = "Q4"),
        h1 = list(suffix = "Half1", display = "H1"),
        h2 = list(suffix = "Half2", display = "H2")
      )
    )

  } else {
    stop(paste("Unknown sport:", sport))
  }
}

