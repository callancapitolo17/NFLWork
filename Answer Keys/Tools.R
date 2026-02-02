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
                       use_spread_line = FALSE) {  # NEW parameter
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
  
  adj_spread <- parent_spread
  adj_total <- parent_total
  
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
  
  for (iter in seq_len(max_iter_mean)) {
    distance_index_generic(dt, adj_spread, adj_total, ss, st, spread_col)  # NEW helper
    setorder(dt, index)
    
    # Use set() to avoid data.table scoping issues
    set(dt, j = "included", value = FALSE)
    set(dt, i = 1L:N, j = "included", value = TRUE)
    
    # compute current means and errors - use direct column access
    mean_s <- mean(dt[included == TRUE][[spread_col]], na.rm = TRUE)
    mean_t <- mean(dt[included == TRUE][["total_line"]], na.rm = TRUE)
    err_s <- mean_s - parent_spread
    err_t <- mean_t - parent_total
    
    # stop if within tolerance
    if (abs(err_s) < tol_mean && abs(err_t) < tol_mean) break
    
    # adjust the parent targets
    adj_spread <- adj_spread - err_s
    adj_total <- adj_total - err_t
  }
  # return updated dt and targets
  list(
    dt = dt,
    parent_spread = parent_spread,
    parent_total = parent_total
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
  dt <- copy(dt)
  n_sample <- N

  # initialize errors
  cover_error <- dt[included == TRUE, sum(actual_cover, na.rm = T)] -
    round(target_cover * n_sample)
  over_error <- dt[included == TRUE, sum(actual_over, na.rm = T)] -
    round(target_over * n_sample)
  M <- nrow(dt)

  repeat {
    if ((abs(cover_error) <= tol_error) && (abs(over_error) <= tol_error)) {
      # Successfully converged
      return(list(
        dt = dt,
        final_N = n_sample,
        cover_error = cover_error,
        over_error = over_error,
        converged = TRUE
      ))
    }

    # mark failures
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

    # ---- if neither helped, signal non-convergence ----
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

  current_N <- N

  repeat {
    # Step 1: mean_match to get sample with correct spread/total means
    mm <- mean_match(DT, current_N, parent_spread, parent_total, ss, st,
                     max_iter_mean, tol_mean, use_spread_line)

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
  cat(sprintf("Generating samples for %d games...\n", nrow(targets)))

  N <- as.integer(N)

  samples <- targets %>%
    pmap(function(id, parent_spread, parent_total, target_cover, target_over) {
      run_answer_key_sample(
        id = id,
        parent_spread = parent_spread,
        parent_total = parent_total,
        target_cover = target_cover,
        target_over = target_over,
        DT = DT, ss = ss, st = st, N = N,
        max_iter_mean = max_iter_mean,
        tol_mean = tol_mean,
        tol_error = tol_error,
        use_spread_line = use_spread_line,
        shrink_factor = shrink_factor,
        min_N = min_N
      )
    })

  # Name the list by game id for easy lookup
  names(samples) <- targets$id

  cat(sprintf("Generated %d samples\n", length(samples)))
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
    # 2. one row per market (h2h, totals, spreads, etc.)
    unnest_longer(bookmakers_markets) %>%
    unnest_wider(bookmakers_markets, names_sep = "_") %>%
    # 3. one row per single outcome (e.g. "Detroit Tigers" @ -113)
    unnest_longer(bookmakers_markets_outcomes) %>%
    unnest_wider(bookmakers_markets_outcomes, names_sep = "_") %>%
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
    line_col = NULL,       # optional: single line column (e.g. "book_total_line" for totals)
    line_col_1 = NULL,     # optional: line for side 1 (e.g. "book_home_spread" for spreads)
    line_col_2 = NULL,     # optional: line for side 2 (e.g. "book_away_spread" for spreads)
    books   = c("betonlineag","kalshi","draftkings"),
    time_col = "commence_time",
    tz_out   = "America/Los_Angeles",
    ev_threshold = 0.05    # minimum EV to include bet (default 5%)
) {
  s1 <- str_extract(size1, "home|away|over|under")
  s2 <- str_extract(size2, "home|away|over|under")

  # Map columns to standard names: metric_side (e.g. prob_home, size_over)
  cols_map <- c(
    set_names(pred1, paste0("prob_", s1)), set_names(pred2, paste0("prob_", s2)),
    set_names(ev1,   paste0("ev_", s1)),   set_names(ev2,   paste0("ev_", s2)),
    set_names(size1, paste0("size_", s1)), set_names(size2, paste0("size_", s2)),
    set_names(odds1, paste0("odds_", s1)), set_names(odds2, paste0("odds_", s2))
  )

  # Add 3rd side if provided (for 3-way markets)
  if (!is.null(size3)) {
    s3 <- str_extract(size3, "home|away|over|under|tie")
    cols_map <- c(cols_map,
      set_names(pred3, paste0("prob_", s3)),
      set_names(ev3,   paste0("ev_", s3)),
      set_names(size3, paste0("size_", s3)),
      set_names(odds3, paste0("odds_", s3))
    )
  }

  # Handle line columns - support both single line_col and separate line_col_1/line_col_2
  has_single_line <- !is.null(line_col) && line_col %in% names(df)
  has_separate_lines <- !is.null(line_col_1) && !is.null(line_col_2) &&
                        line_col_1 %in% names(df) && line_col_2 %in% names(df)

  result <- df %>%
    filter(bookmaker_key %in% books) %>%
    mutate(pt_start_time = lubridate::with_tz(lubridate::ymd_hms(.data[[time_col]], tz = "UTC"), tzone = tz_out))

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
    # Keep only highest bet size per game per market per side (best book for each side)
    group_by(id, market, bet_on) %>%
    filter(size == max(size)) %>%
    ungroup() %>%
    arrange(desc(size)) %>%
    select(home_team, away_team, pt_start_time, bookmaker_key, market, bet_on, line, bet_size = size, ev, odds, prob)
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
    margin_col = "game_home_margin_period"
) {
  # Validate inputs
  if (length(periods) != length(markets)) {
    stop("periods and markets must be same length and correspond to each other")
  }

  cat(sprintf("Fetching %d markets for %d events...\n", length(markets), length(events$id)))

  # 1) Batch fetch all markets for all events
  all_odds <- expand_grid(event_id = events$id, market = markets) %>%
    mutate(
      data = map2(event_id, market, possibly(
        ~ {
          result <- fetch_event_odds(.x, .y, sport_key)
          if (nrow(result) > 0 && "bookmakers" %in% names(result)) {
            result
          } else {
            tibble()
          }
        },
        otherwise = tibble()
      ))
    ) %>%
    filter(map_lgl(data, ~ nrow(.x) > 0)) %>%
    unnest(data) %>%
    select(-event_id)

  if (nrow(all_odds) == 0) stop("No odds data returned")

  # 2) Flatten odds - each row has a 2-element list (home/away outcomes)
  # unnest_wider expands these into _1 and _2 columns
  flat_odds <- all_odds %>%
    flatten_event_odds() %>%
    select(-market_key) %>%
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
      odds1 = "book_home_market", odds2 = "book_away_market"
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

  cat(sprintf("Fetching %d 3-way markets for %d events...\n", length(markets), length(events$id)))

  # 1) Batch fetch all markets for all events
  all_odds <- expand_grid(event_id = events$id, market = markets) %>%
    mutate(
      data = map2(event_id, market, possibly(
        ~ {
          result <- fetch_event_odds(.x, .y, sport_key)
          if (nrow(result) > 0 && "bookmakers" %in% names(result)) {
            result
          } else {
            tibble()
          }
        },
        otherwise = tibble()
      ))
    ) %>%
    filter(map_lgl(data, ~ nrow(.x) > 0)) %>%
    unnest(data) %>%
    select(-event_id)

  if (nrow(all_odds) == 0) {
    cat("No 3-way odds data returned\n")
    return(list(predictions = tibble(), prediction_set = tibble(), bets = tibble(),
                markets_summary = tibble(market = character(), n_bets = integer(),
                                         total_stake = numeric(), avg_ev = numeric(), max_ev = numeric())))
  }

  # 2) Flatten 3-way odds - each row has 3 outcomes (home/away/tie)
  flat_odds <- all_odds %>%
    flatten_event_odds() %>%
    select(-market_key) %>%
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
    # Devig 3-way odds
    mutate(devigged = pmap(list(book_home_odds, book_away_odds, book_tie_odds),
                           ~ devig_american_3way(..1, ..2, ..3))) %>%
    unnest_wider(devigged) %>%
    rename(book_prob_home = p_home, book_prob_away = p_away, book_prob_tie = p_tie) %>%
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
      pred3 = "tie_prob", ev3 = "tie_ev", size3 = "tie_bet_size", odds3 = "book_tie_odds"
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
    total_col = "game_total_period"
) {
  # Validate inputs
  if (length(periods) != length(markets)) {
    stop("periods and markets must be same length and correspond to each other")
  }

  cat(sprintf("Fetching %d totals markets for %d events...\n", length(markets), length(events$id)))

  # 1) Batch fetch all markets
  all_odds <- expand_grid(event_id = events$id, market = markets) %>%
    mutate(
      data = map2(event_id, market, possibly(
        ~ {
          result <- fetch_event_odds(.x, .y, sport_key)
          if (nrow(result) > 0 && "bookmakers" %in% names(result)) {
            result
          } else {
            tibble()
          }
        },
        otherwise = tibble()
      ))
    ) %>%
    filter(map_lgl(data, ~ nrow(.x) > 0)) %>%
    unnest(data) %>%
    select(-event_id)

  if (nrow(all_odds) == 0) stop("No totals odds data returned")

  # 2) Flatten odds - each row has a 2-element list (over/under outcomes)
  # unnest_wider expands these into _1 and _2 columns
  flat_odds <- all_odds %>%
    flatten_event_odds() %>%
    select(-market_key) %>%
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
      line_col = "book_total_line"
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
    margin_col = "game_home_margin_period"
) {
  # Validate inputs
  if (length(periods) != length(markets)) {
    stop("periods and markets must be same length and correspond to each other")
  }

  cat(sprintf("Fetching %d spreads markets for %d events...\n", length(markets), length(events$id)))

  # 1) Batch fetch all markets
  all_odds <- expand_grid(event_id = events$id, market = markets) %>%
    mutate(
      data = map2(event_id, market, possibly(
        ~ {
          result <- fetch_event_odds(.x, .y, sport_key)
          if (nrow(result) > 0 && "bookmakers" %in% names(result)) {
            result
          } else {
            tibble()
          }
        },
        otherwise = tibble()
      ))
    ) %>%
    filter(map_lgl(data, ~ nrow(.x) > 0)) %>%
    unnest(data) %>%
    select(-event_id)

  if (nrow(all_odds) == 0) stop("No spreads odds data returned")

  # 2) Flatten odds - each bookmaker row has LISTS of outcomes (multiple spread lines)
  # First unnest_longer to get one row per outcome, then group to pair home/away
  flat_odds <- all_odds %>%
    flatten_event_odds() %>%
    select(-market_key) %>%
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
      line_col_2 = "book_away_spread"   # Away bet shows away spread (e.g., -1.5)
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
    ev_threshold = 0.05
) {
  # Validate inputs
  if (length(periods) != length(markets)) {
    stop("periods and markets must have same length")
  }

  # Create market-to-period mapping
  market_period_map <- tibble(market = markets, period = as.character(periods))

  # 1) Fetch odds for all markets
  cat(sprintf("Fetching %d team totals markets for %d events...\n", length(markets), nrow(events)))

  raw_odds <- map_dfr(markets, function(m) {
    event_ids <- events$id
    map_dfr(event_ids, function(eid) {
      tryCatch({
        res <- httr::GET(
          paste0("https://api.the-odds-api.com/v4/sports/", sport_key, "/events/", eid, "/odds"),
          query = list(
            apiKey = Sys.getenv("ODDS_API_KEY"),
            regions = "us,us2,eu",
            markets = m,
            oddsFormat = "american"
          )
        )
        httr::stop_for_status(res)

        out <- fromJSON(content(res, "text"), flatten = FALSE)
        if (length(out$bookmakers) == 0) return(tibble())
        out$market <- m
        as_tibble(out)
      }, error = function(e) tibble())
    })
  })

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
    unnest_longer(bookmakers_markets) %>%
    unnest_wider(bookmakers_markets, names_sep = "_") %>%
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
    # Filter to best book per bet
    group_by(id, market, bet_on, line) %>%
    filter(bet_size == max(bet_size)) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
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

  return(result)
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
        select(id, home_team, away_team, market, book_home_spread, home_cover_prob, away_cover_prob),
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
      away_bet_size = kelly_stake(away_ev, wz_prob_away, bankroll, kelly_mult),
      # Add commence_time for format_bets_table (use fetch_time as proxy)
      commence_time = as.POSIXct(fetch_time, tz = "UTC")
    )

  # Format bets - pass books = "wagerzon" to include Wagerzon in filter
  bets <- prediction_set %>%
    format_bets_table(
      pred1 = "home_cover_prob", pred2 = "away_cover_prob",
      ev1 = "home_ev", ev2 = "away_ev",
      size1 = "home_bet_size", size2 = "away_bet_size",
      odds1 = "odds_home", odds2 = "odds_away",
      line_col_1 = "home_spread",
      line_col_2 = "away_spread",
      books = "wagerzon",
      ev_threshold = ev_threshold
    )

  summary <- bets %>%
    group_by(market) %>%
    summarise(n_bets = n(), total_stake = sum(bet_size), avg_ev = mean(ev), max_ev = max(ev), .groups = "drop") %>%
    arrange(desc(total_stake))

  cat(sprintf("Generated %d Wagerzon spread bets\n", nrow(bets)))

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
        select(id, home_team, away_team, market, book_total_line, over_prob, under_prob),
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
      under_bet_size = kelly_stake(under_ev, wz_prob_under, bankroll, kelly_mult),
      commence_time = as.POSIXct(fetch_time, tz = "UTC")
    )

  bets <- prediction_set %>%
    format_bets_table(
      pred1 = "over_prob", pred2 = "under_prob",
      ev1 = "over_ev", ev2 = "under_ev",
      size1 = "over_bet_size", size2 = "under_bet_size",
      odds1 = "odds_over", odds2 = "odds_under",
      line_col_1 = "line", line_col_2 = "line",
      books = "wagerzon",
      ev_threshold = ev_threshold
    )

  summary <- bets %>%
    group_by(market) %>%
    summarise(n_bets = n(), total_stake = sum(bet_size), avg_ev = mean(ev), max_ev = max(ev), .groups = "drop") %>%
    arrange(desc(total_stake))

  cat(sprintf("Generated %d Wagerzon totals bets\n", nrow(bets)))

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
        select(id, home_team, away_team, market, home_win_prob, away_win_prob),
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
      away_bet_size = kelly_stake(away_ev, wz_prob_away, bankroll, kelly_mult),
      commence_time = as.POSIXct(fetch_time, tz = "UTC")
    )

  bets <- prediction_set %>%
    format_bets_table(
      pred1 = "home_win_prob", pred2 = "away_win_prob",
      ev1 = "home_ev", ev2 = "away_ev",
      size1 = "home_bet_size", size2 = "away_bet_size",
      odds1 = "odds_home", odds2 = "odds_away",
      books = "wagerzon",
      ev_threshold = ev_threshold
    )

  summary <- bets %>%
    group_by(market) %>%
    summarise(n_bets = n(), total_stake = sum(bet_size), avg_ev = mean(ev), max_ev = max(ev), .groups = "drop") %>%
    arrange(desc(total_stake))

  cat(sprintf("Generated %d Wagerzon moneyline bets\n", nrow(bets)))

  list(predictions = predictions, prediction_set = prediction_set, bets = bets, markets_summary = summary)
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

