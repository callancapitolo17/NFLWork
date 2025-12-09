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
  # highest total weight across all sportsbooks posting it.
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
    group_by(.data[[game_id_col]], .data[[date_col]], .data[[line_col]]) %>%
    mutate(total_weight = sum(.data[[weight_col]],na.rm =T)) %>% 
    ungroup() %>% 
    group_by(.data[[game_id_col]],.data[[date_col]]) %>% 
    filter(total_weight == max(total_weight, na.rm = TRUE)) %>%
    slice_head(n = 1) %>%   # handle ties safely
    ungroup() %>% 
    mutate(market1_weighted_prob = .data[[market1]] * .data[[weight_col]],
           market2_weighted_prob = .data[[market2]] * .data[[weight_col]]) %>% 
    group_by(.data[[game_id_col]], .data[[date_col]], .data[[line_col]],.data[[time_col]], .data[[home]],.data[[away]]) %>%
    summarize(consensus_market1 = sum(market1_weighted_prob, na.rm = T) / sum(.data[[weight_col]], na.rm = T),
              consensus_market2 = sum(market2_weighted_prob, na.rm = T) / sum(.data[[weight_col]], na.rm = T),
              .groups = "drop") %>% 
    rename(!!paste0("consensus_",market1,sep = "") := consensus_market1,
           !!paste0("consensus_",market2,sep = "") := consensus_market2)
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
  adj_spread <- parent_spread
  adj_total <- parent_total
  
  # Determine which column to match against
  spread_col <- if (use_spread_line) "home_spread" else "home_ml_odds"
  
  for (iter in seq_len(max_iter_mean)) {
    distance_index_generic(dt, adj_spread, adj_total, ss, st, spread_col)  # NEW helper
    setorder(dt, index)
    dt[, included := FALSE]
    dt[1:N, included := TRUE]
    
    # compute current means and errors
    mean_s <- dt[included == TRUE, mean(get(spread_col))]
    mean_t <- dt[included == TRUE, mean(total_line)]
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
  dt[, index := ((get(spread_col) - ps) / ss)^2 + ((total_line - pt) / st)^2]
}

# Function 2: balance sample by greedy remove/add + shrink‐on‐stall
balance_sample <- function(dt, N, target_cover, target_over, tol_error) {
  dt <- copy(dt)
  n_sample <- N  # CHANGED: use different variable name to avoid confusion
  
  # initialize errors
  cover_error <- dt[included == TRUE, sum(actual_cover, na.rm =T)] -
    round(target_cover * n_sample)
  over_error <- dt[included == TRUE, sum(actual_over,na.rm = T)] -
    round(target_over * n_sample)
  M <- nrow(dt)
  
  repeat {
    if ((abs(cover_error) <= tol_error) && (abs(over_error) <= tol_error)) {break}
    
    # mark failures
    removal_failed <- TRUE
    addition_failed <- TRUE
    
    # ---- remove one game from 1..n_sample ----
    in1toN <- dt[seq_len(n_sample), included]  # CHANGED: use seq_len
    cov_help1 <- (cover_error > 0 & dt[seq_len(n_sample), actual_cover] == 1) |
      (cover_error < 0 & dt[seq_len(n_sample), actual_cover] == 0)
    ov_help1 <- (over_error > 0 & dt[seq_len(n_sample), actual_over] == 1) |
      (over_error < 0 & dt[seq_len(n_sample), actual_over] == 0)
    
    both1 <- which(in1toN & cov_help1 & ov_help1)
    either1 <- which(in1toN & (cov_help1 | ov_help1))
    
    if (length(both1) > 0) {
      i <- max(both1)
      dt[i, included := FALSE]
      removal_failed <- FALSE
    } else if (length(either1) > 0) {
      i <- max(either1)
      dt[i, included := FALSE]
      removal_failed <- FALSE
    }
    
    if (!removal_failed) {
      cover_error <- dt[included == TRUE, sum(actual_cover,na.rm = T)] -
        round(target_cover * n_sample)
      over_error <- dt[included == TRUE, sum(actual_over, na.rm = T)] -
        round(target_over * n_sample)
    }
    
    # ---- add one game from (n_sample+1)..M ----
    idx_range <- seq(n_sample + 1, M)  # CHANGED: explicit range
    inNplus <- dt[idx_range, included]
    cov_help2 <- (cover_error > 0 & dt[idx_range, actual_cover] == 0) |
      (cover_error < 0 & dt[idx_range, actual_cover] == 1)
    ov_help2 <- (over_error > 0 & dt[idx_range, actual_over] == 0) |
      (over_error < 0 & dt[idx_range, actual_over] == 1)
    
    both2 <- which(!inNplus & cov_help2 & ov_help2)
    either2 <- which(!inNplus & (cov_help2 | ov_help2))
    
    if (length(both2) > 0) {
      j <- max(both2) + n_sample
      dt[j, included := TRUE]
      addition_failed <- FALSE
    } else if (length(either2) > 0) {
      j <- max(either2) + n_sample
      dt[j, included := TRUE]
      addition_failed <- FALSE
    }
    
    if (!addition_failed) {
      cover_error <- dt[included == TRUE, sum(actual_cover, na.rm = T)] -
        round(target_cover * n_sample)
      over_error <- dt[included == TRUE, sum(actual_over, na.rm = T)] -
        round(target_over * n_sample)
    }
    
    # ---- if neither helped, shrink sample size ----
    if (removal_failed && addition_failed) {
      worst_i <- dt[seq_len(n_sample), which.max(index)]  # CHANGED
      dt[worst_i, included := FALSE]
      n_sample <- n_sample - 1  # CHANGED
      cover_error <- dt[included == TRUE, sum(actual_cover, na.rm = T)] -
        round(target_cover * n_sample)
      over_error <- dt[included == TRUE, sum(actual_over, na.rm = T)] -
        round(target_over * n_sample)
      next
    }
  }
  
  # return final dt and stats
  list(
    dt = dt,
    final_N = n_sample,  # CHANGED
    cover_error = cover_error,
    over_error = over_error
  )
}

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

run_means_for_id_total <- function(id, parent_spread, parent_total, target_cover, target_over, totals, inning,
                                   DT, ss, st, N,
                                   max_iter_mean = 500, tol_mean = 0.005, tol_error = 3) {
  mm <- mean_match(DT, N, parent_spread, parent_total, ss, st, max_iter_mean, tol_mean)
  bal <- balance_sample(mm$dt, N, target_cover, target_over, tol_error)
  
  inc_df <- bal$dt %>% filter(included == TRUE)
  
  # summarise all game_home_ml* columns
  vals <- inc_df %>%
    select("inning_total" = !!sym(paste0("full_game_total_inning_", inning, sep = ""))) %>%
    pull()
  pct_cols <- set_names(
    map(totals, ~ sum(vals > .x, na.rm = T) / sum(vals != .x, na.rm = T)),
    paste0("pct_over_", gsub("\\.", "_", totals))
  )
  
  tibble(!!!pct_cols)
}

run_means_for_id_spread <- function(id, parent_spread, parent_total, target_cover, target_over, spreads, inning,
                                    DT, ss, st, N,
                                    max_iter_mean = 500, tol_mean = 0.005, tol_error = 2) {
  mm <- mean_match(DT, N, parent_spread, parent_total, ss, st, max_iter_mean, tol_mean)
  bal <- balance_sample(mm$dt, N, target_cover, target_over, tol_error)

  inc_df <- bal$dt %>% filter(included == TRUE)

  spreads <- spreads %>%
    filter(`id` == id) %>%
    pull(book_home_spread) %>%
    unique()

  # summarise all spread_colums
  vals <- inc_df %>%
    select("spread" = !!sym(paste0("game_home_margin_in_", inning, sep = ""))) %>% # select the correct column based on the inning
    pull()
  pct_cols <- set_names(
    map(spreads, ~ sum(vals > -.x, na.rm = T) / sum(vals != -.x, na.rm = T)), # flip sign of spread - to cover -0.5 margin needs to be greater than 0.5
    paste0("pct_home_cover_", gsub("\\.", "_", spreads))
  )

  tibble(!!!pct_cols)
}

#generates bets for moneylines
run_means_for_id <- function(id, parent_spread, parent_total, target_cover, target_over,
                             DT, ss, st, N,
                             max_iter_mean = 500, tol_mean = 0.005, tol_error = 1,
                             use_spread_line = FALSE,
                             margin_col = "game_home_margin_in") {  
  mm <- mean_match(DT, N, parent_spread, parent_total, ss, st, 
                   max_iter_mean, tol_mean, use_spread_line)
  bal <- balance_sample(mm$dt, N, target_cover, target_over, tol_error)
  inc_df <- bal$dt %>% filter(included == TRUE)
  bets_summary <- inc_df %>%
    summarise(across(starts_with(margin_col), ~ sum(.x > 0, na.rm = TRUE) / sum(.x != 0, na.rm = T))) %>% 
    ungroup()
  bets_summary
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
    books   = c("fliff", "rebet", "novig", "prophetx"),
    time_col = "commence_time",
    tz_out   = "America/Los_Angeles"
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
  
  df %>%
    filter(bookmaker_key %in% books) %>%
    mutate(pt_start_time = lubridate::with_tz(lubridate::ymd_hms(.data[[time_col]], tz = "UTC"), tzone = tz_out)) %>%
    # Select and rename dynamically in one step
    select(home_team, away_team, pt_start_time, bookmaker_key, market, all_of(cols_map)) %>%
    pivot_longer(
      cols = -c(home_team, away_team, pt_start_time, bookmaker_key, market),
      names_to = c(".value", "bet_on"),
      names_sep = "_"
    ) %>%
    # Map 'home'/'away' to team names; keep 'Over'/'Under' as is
    mutate(bet_on = case_when(
      bet_on == "home" ~ home_team,
      bet_on == "away" ~ away_team,
      TRUE ~ str_to_title(bet_on)
    )) %>%
    filter(size > 0) %>% 
    arrange(desc(size)) %>%
    select(home_team, away_team, pt_start_time, bookmaker_key, market, bet_on, bet_size = size, ev, odds, prob)
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
  # 2) Elihu engine for each id
  final_preds <- targets %>%
    mutate(res = pmap(
      list(id, parent_spread, parent_total, target_cover, target_over),
      ~ run_means_for_id(..1, ..2, ..3, ..4, ..5,
                         DT = DT, ss = ss, st = st, N = N,
                         use_spread_line = use_spread_line,
                         max_iter_mean = 500,
                         margin_col = margin_col))) %>%
    unnest(res) %>%
    ungroup() %>% 
    inner_join(
      consensus_odds %>%
        ungroup() %>% 
        select(id, home_team, away_team, commence_time),
      by = c("id", "home_team", "commence_time")
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
  
  # 3) Elihu engine for totals at all posted numbers
  total_final_preds <- targets %>%
    mutate(res = pmap(
      list(id, parent_spread, parent_total, target_cover, target_over),
      ~ run_means_for_id_total(..1, ..2, ..3, ..4, ..5,
                               totals = unique(flat_betting_odds$book_market_total),
                               period = period,        # NOTE: you'll need to rename arg in run_means_for_id_total if you change from inning
                               DT = DT, ss = ss, st = st, N = N)
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
    unnest_wider(outcome_name, names_sep = "_") %>%
    unnest_wider(closing_odds, names_sep = "_") %>%
    unnest_wider(bookmakers_markets_outcomes_point, names_sep = "_") %>%
    mutate(
      home_odds   = ifelse(home_team == outcome_name_1, closing_odds_1, closing_odds_2),
      away_odds   = ifelse(away_team == outcome_name_1, closing_odds_1, closing_odds_2),
      home_spread = ifelse(home_team == outcome_name_1, bookmakers_markets_outcomes_point_1, bookmakers_markets_outcomes_point_2),
      away_spread = ifelse(away_team == outcome_name_1, bookmakers_markets_outcomes_point_1, bookmakers_markets_outcomes_point_2)
    ) %>%
    group_by(
      id, commence_time, home_team, away_team,
      bookmaker_key, bookmaker_title
    ) %>%
    summarise(
      book_home_odds   = first(home_odds),
      book_home_spread = first(home_spread),
      book_away_odds   = first(away_odds),
      book_away_spread = first(away_spread),
      .groups = "drop"
    )
  
  # 3) Elihu engine for spreads at all posted numbers
  spread_final_preds <- targets %>%
    mutate(res = pmap(
      list(id, parent_spread, parent_total, target_cover, target_over),
      ~ run_means_for_id_spread(..1, ..2, ..3, ..4, ..5,
                                spreads = flat_betting_odds %>% ungroup() %>% select(id, book_home_spread),
                                period  = period,    # rename in run_means_for_id_spread if needed
                                DT = DT, ss = ss, st = st, N = N)
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
      odds2 = "book_away_odds"
    ) %>%
    arrange(desc(home_bet_size))
  
  list(
    predictions    = spread_predictions,
    prediction_set = prediction_set,
    bets           = bets
  )
}

# Efficient Multi-Market Moneyline Builder ----
build_multi_moneyline_markets <- function(
    DT,
    consensus_odds,
    ss, st, N,
    periods,            # Vector of periods: c(1, 2, 3, 4) for quarters
    events,
    markets,            # Vector of markets: c("h2h_q1", "h2h_q2", "h2h_q3", "h2h_q4")
    sport_key,
    bankroll   = 200,
    kelly_mult = 0.25,
    targets = NULL,
    use_spread_line = FALSE,
    margin_col = "game_home_margin_period"
) {
  
  # Validate inputs
  if (length(periods) != length(markets)) {
    stop("periods and markets must be same length and correspond to each other")
  }
  
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
  
  cat(sprintf("Fetching %d markets for %d events (%d API calls)...\n", 
              length(markets), length(events$id), length(markets) * length(events$id)))
  
  # 1) Batch fetch all markets for all events
  all_odds <- expand_grid(event_id = events$id, market = markets) %>%
    mutate(
      data = map2(event_id, market, possibly(
        ~ {
          result <- fetch_event_odds(.x, .y, sport_key)
          if (nrow(result) > 0 && "bookmakers" %in% names(result)) {
            result  # Don't add market here, it's already in the parent tibble
          } else {
            tibble()
          }
        },
        otherwise = tibble()
      ))
    ) %>%
    filter(map_lgl(data, ~ nrow(.x) > 0)) %>%  # Filter out empty results before unnesting
    unnest(data) %>%
    select(-event_id)
  
  if (nrow(all_odds) == 0) stop("No odds data returned")
  
  # 2) Flatten all odds
  flat_odds <- all_odds %>%
    flatten_event_odds() %>%
    select(-market_key) %>%  # Remove market_key from API response, we already have 'market' column
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
  
  # 3) Run predictions ONCE - generates probabilities for all periods
  cat(sprintf("Running prediction engine for %d events across %d periods...\n", 
              nrow(targets), length(periods)))
  
  predictions_raw <- targets %>%
    ungroup() %>%  # Remove any grouping from targets
    select(id, parent_spread, parent_total, target_cover, target_over) %>%  # Only keep needed columns
    mutate(res = pmap(
      list(id, parent_spread, parent_total, target_cover, target_over),
      ~ run_means_for_id(..1, ..2, ..3, ..4, ..5, DT, ss, st, N,
                         use_spread_line = use_spread_line, 
                         margin_col = margin_col)
    )) %>%
    unnest(res)
  
  # Join with consensus to get team names and time
  consensus_info <- consensus_odds %>% 
    ungroup()
  
  # Debug: check what columns exist
  if (!all(c("id", "home_team", "away_team", "commence_time") %in% names(consensus_info))) {
    stop(paste("Missing columns in consensus_odds. Available columns:", 
               paste(names(consensus_info), collapse = ", ")))
  }
  
  consensus_info <- consensus_info %>%
    select(id, home_team, away_team, commence_time)
  
  # Convert commence_time if it's character
  if (is.character(consensus_info$commence_time)) {
    consensus_info <- consensus_info %>%
      mutate(commence_time = ymd_hms(commence_time, tz = "UTC"))
  }
  
  predictions <- predictions_raw %>%
    inner_join(consensus_info, by = "id")
  
  # Debug: verify predictions has the needed columns after join
  if (nrow(predictions) == 0) {
    stop("No predictions after joining with consensus_info. Check that IDs match between targets and consensus_odds.")
  }
  if (!"home_team" %in% names(predictions)) {
    stop(paste("home_team missing from predictions after join. predictions columns:", 
               paste(names(predictions), collapse = ", "),
               "\nconsensus_info columns:", paste(names(consensus_info), collapse = ", ")))
  }
  
  # 4) Create market-to-period mapping and pivot predictions
  market_period_map <- tibble(market = markets, period = periods)
  
  # Select only columns we need - use any_of to handle missing gracefully
  margin_cols <- names(predictions)[startsWith(names(predictions), margin_col)]
  if (length(margin_cols) == 0) {
    stop(paste("No columns starting with", margin_col, "found in predictions. Available columns:", 
               paste(names(predictions), collapse = ", ")))
  }
  
  predictions_long <- predictions %>%
    select(id, home_team, away_team, commence_time, all_of(margin_cols)) %>%
    pivot_longer(
      cols = starts_with(margin_col),
      names_to = "period_col",
      values_to = "home_win_prob"
    ) %>%
    mutate(
      period = as.numeric(str_extract(period_col, "\\d+"))
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


