# =============================================================================
# Kalshi March Madness Edge Finder
# =============================================================================
# Fetches tournament prop markets from Kalshi, runs bracket sim,
# maps each market to a sim probability, computes EV after taker fees.
#
# Returns a data frame of edges sorted by best EV%.

suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(tidyr)
})

KALSHI_BASE <- "https://api.elections.kalshi.com/trade-api/v2"
TAKER_FEE_RATE <- 0.07

# =============================================================================
# KALSHI API
# =============================================================================

kalshi_request <- function(path) {
  url <- paste0(KALSHI_BASE, path)
  tryCatch({
    fromJSON(url, simplifyDataFrame = FALSE)
  }, error = function(e) {
    warning(sprintf("Kalshi API error: %s", e$message))
    NULL
  })
}

fetch_kalshi_markets <- function(series_ticker = NULL, event_ticker = NULL, status = "open") {
  markets <- list()
  cursor <- NULL
  repeat {
    path <- "/markets?"
    if (!is.null(series_ticker)) path <- paste0(path, "series_ticker=", series_ticker, "&")
    if (!is.null(event_ticker)) path <- paste0(path, "event_ticker=", event_ticker, "&")
    path <- paste0(path, "status=", status, "&limit=200")
    if (!is.null(cursor) && nchar(cursor) > 0) path <- paste0(path, "&cursor=", cursor)

    data <- kalshi_request(path)
    if (is.null(data)) break

    batch <- data$markets
    if (length(batch) == 0) break
    markets <- c(markets, batch)

    cursor <- data$cursor
    if (is.null(cursor) || nchar(cursor) == 0) break
    Sys.sleep(0.1)
  }
  markets
}

# =============================================================================
# EV CALCULATION (matches taker.py exactly)
# =============================================================================
# Named kalshi_compute_ev to avoid collision with cbbdata::compute_ev

kalshi_fee_cents <- function(price_cents) {
  p <- price_cents / 100
  TAKER_FEE_RATE * p * (1 - p) * 100
}

kalshi_compute_ev <- function(sim_prob, yes_ask_cents, yes_bid_cents) {
  fair_cents <- sim_prob * 100
  no_ask_cents <- 100 - yes_bid_cents

  # YES side
  yes_fee <- kalshi_fee_cents(yes_ask_cents)
  yes_edge <- fair_cents - yes_ask_cents - yes_fee
  yes_ev_pct <- if (yes_ask_cents > 0) yes_edge / yes_ask_cents else -Inf

  # NO side
  no_fee <- kalshi_fee_cents(no_ask_cents)
  no_edge <- (100 - fair_cents) - no_ask_cents - no_fee
  no_ev_pct <- if (no_ask_cents > 0) no_edge / no_ask_cents else -Inf

  list(
    yes_ev_pct = yes_ev_pct,
    no_ev_pct = no_ev_pct,
    yes_fee = yes_fee,
    no_fee = no_fee,
    yes_edge = yes_edge,
    no_edge = no_edge,
    best_side = if (yes_ev_pct >= no_ev_pct) "YES" else "NO",
    best_ev_pct = max(yes_ev_pct, no_ev_pct)
  )
}

# =============================================================================
# BRACKET SIMULATION
# =============================================================================

run_bracket_sim <- function(nflwork_root, n_sims = 10000) {
  mm_dir <- file.path(nflwork_root, "March Madness")

  old_wd <- getwd()
  setwd(mm_dir)
  on.exit(setwd(old_wd))
  source("shared.R", local = TRUE)
  source("espn_bracket.R", local = TRUE)

  bracket_result <- fetch_espn_bracket()
  teams_std <- get_teams_std()
  final_bracket <- bracket_result$bracket %>% select(team, seed, region, play_in)
  bracket_with_ratings <- fetch_bracket_with_ratings(final_bracket, teams_std)

  games_played <- bracket_result$games
  bracket_64 <- resolve_first_four(bracket_with_ratings, games_played)

  # Filter out eliminated teams for accurate remaining-team probs
  eliminated <- games_played %>%
    filter(status == "final", round != "First Four", !is.na(winner)) %>%
    mutate(loser = ifelse(winner == team1, team2, team1)) %>%
    pull(loser)

  region_order <- get_region_order(bracket_64)
  n_remaining <- nrow(bracket_64) - length(eliminated)

  cat(sprintf("  Running %s bracket sims (%d teams remaining)...\n",
              format(n_sims, big.mark = ","), n_remaining))
  t0 <- Sys.time()

  sim_results <- map_dfr(1:n_sims, ~ simulate_tournament_fast(bracket_64, region_order))

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("  Done in %.0f seconds.\n", elapsed))

  # Team-level probabilities
  team_probs <- sim_results %>%
    group_by(team, seed) %>%
    summarise(
      Round_32 = mean(Round_32),
      Sweet_16 = mean(Sweet_16),
      Elite_8 = mean(Elite_8),
      Final_4 = mean(Final_4),
      Title_Game = mean(Title_Game),
      Champion = mean(Champion),
      .groups = "drop"
    )

  list(raw_sims = sim_results, team_probs = team_probs, n_sims = n_sims)
}

# =============================================================================
# TEAM NAME MATCHING
# =============================================================================

parse_kalshi_team <- function(title) {
  m <- str_match(title, "^Will (.+?) qualify for")
  if (!is.na(m[1, 2])) return(m[1, 2])
  NULL
}

match_team <- function(kalshi_name, team_probs) {
  if (is.null(kalshi_name)) return(NA_character_)

  normalize <- function(x) {
    x <- tolower(x)
    x <- str_replace_all(x, "['\\.\\-]", "")
    x <- str_replace_all(x, "\\bst\\.?\\b", "st")
    x <- str_replace_all(x, "\\bstate\\b", "st")
    x <- str_replace_all(x, "\\buniversity\\b", "")
    str_squish(x)
  }

  kn <- normalize(kalshi_name)
  candidates <- team_probs$team
  cn <- normalize(candidates)

  # Exact
  idx <- which(cn == kn)
  if (length(idx) == 1) return(candidates[idx])

  # Substring
  idx <- which(str_detect(cn, fixed(kn)) | str_detect(kn, cn))
  if (length(idx) == 1) return(candidates[idx])

  # Word overlap
  kn_words <- str_split(kn, "\\s+")[[1]]
  scores <- sapply(cn, function(c) {
    c_words <- str_split(c, "\\s+")[[1]]
    sum(kn_words %in% c_words) / max(length(kn_words), 1)
  })
  best <- which.max(scores)
  if (scores[best] >= 0.5) return(candidates[best])

  NA_character_
}

# =============================================================================
# MARKET-TO-SIM MAPPING
# =============================================================================

map_round_advancement <- function(markets, team_probs, round_col) {
  map_dfr(markets, function(m) {
    title <- m$title %||% ""
    kalshi_team <- parse_kalshi_team(title)
    sim_team <- match_team(kalshi_team, team_probs)
    if (is.na(sim_team)) return(tibble())

    sim_prob <- team_probs[[round_col]][team_probs$team == sim_team]
    if (length(sim_prob) == 0) return(tibble())

    yes_ask <- round(as.numeric(m$yes_ask_dollars %||% 0) * 100)
    yes_bid <- round(as.numeric(m$yes_bid_dollars %||% 0) * 100)
    ev <- kalshi_compute_ev(sim_prob, yes_ask, yes_bid)

    tibble(
      ticker = m$ticker %||% "",
      category = paste0("Round: ", gsub("_", " ", round_col)),
      title = title,
      team = sim_team,
      sim_prob = sim_prob,
      yes_ask = yes_ask, yes_bid = yes_bid,
      yes_fee = ev$yes_fee, no_fee = ev$no_fee,
      yes_ev_pct = ev$yes_ev_pct, no_ev_pct = ev$no_ev_pct,
      yes_edge = ev$yes_edge, no_edge = ev$no_edge,
      best_side = ev$best_side, best_ev_pct = ev$best_ev_pct
    )
  })
}

map_seed_props <- function(markets, raw_sims, category_name, compute_fn) {
  # Pre-compute sim_id once
  # Detect rows-per-sim by finding the second occurrence of the first team
  first_team <- raw_sims$team[1]
  second_occurrence <- which(raw_sims$team == first_team)[2]
  n_per_sim <- second_occurrence - 1L
  n_sims <- nrow(raw_sims) / n_per_sim
  raw_sims$sim_id <- rep(1:n_sims, each = n_per_sim)

  map_dfr(markets, function(m) {
    title <- m$title %||% ""
    floor_strike <- as.numeric(m$floor_strike %||% NA)
    cap_strike <- as.numeric(m$cap_strike %||% NA)
    strike_type <- m$strike_type %||% "greater_or_equal"
    event_ticker <- m$event_ticker %||% ""

    sim_prob <- compute_fn(raw_sims, floor_strike, cap_strike, strike_type, event_ticker, title, n_sims)
    if (is.na(sim_prob)) return(tibble())

    yes_ask <- round(as.numeric(m$yes_ask_dollars %||% 0) * 100)
    yes_bid <- round(as.numeric(m$yes_bid_dollars %||% 0) * 100)
    ev <- kalshi_compute_ev(sim_prob, yes_ask, yes_bid)

    tibble(
      ticker = m$ticker %||% "",
      category = category_name,
      title = title,
      team = NA_character_,
      sim_prob = sim_prob,
      yes_ask = yes_ask, yes_bid = yes_bid,
      yes_fee = ev$yes_fee, no_fee = ev$no_fee,
      yes_ev_pct = ev$yes_ev_pct, no_ev_pct = ev$no_ev_pct,
      yes_edge = ev$yes_edge, no_edge = ev$no_edge,
      best_side = ev$best_side, best_ev_pct = ev$best_ev_pct
    )
  })
}

# =============================================================================
# SEED PROP COMPUTE FUNCTIONS
# =============================================================================

compute_seed_sum <- function(raw_sims, floor_strike, cap_strike, strike_type, event_ticker, title, n_sims) {
  round_col <- if (grepl("F4", event_ticker)) "Final_4"
               else if (grepl("T2", event_ticker)) "Title_Game"
               else return(NA_real_)

  advancing <- raw_sims %>%
    filter(.data[[round_col]] == 1) %>%
    group_by(sim_id) %>%
    summarise(seed_sum = sum(seed), .groups = "drop")

  if (nrow(advancing) == 0) return(NA_real_)

  all_sims <- tibble(sim_id = 1:n_sims)
  advancing <- all_sims %>%
    left_join(advancing, by = "sim_id") %>%
    mutate(seed_sum = coalesce(seed_sum, 0L))

  if (strike_type == "between" && !is.na(cap_strike)) {
    mean(advancing$seed_sum >= floor_strike & advancing$seed_sum <= cap_strike)
  } else {
    mean(advancing$seed_sum >= floor_strike)
  }
}

compute_seed_count_in_round <- function(raw_sims, floor_strike, cap_strike, strike_type, event_ticker, title, n_sims) {
  seed_match <- str_match(event_ticker, "S(\\d+)")
  if (is.na(seed_match[1, 2])) return(NA_real_)
  target_seed <- as.integer(seed_match[1, 2])

  round_col <- case_when(
    grepl("R8", event_ticker) ~ "Elite_8",
    grepl("R16", event_ticker) ~ "Sweet_16",
    grepl("R32", event_ticker) ~ "Round_32",
    grepl("F4", event_ticker) ~ "Final_4",
    TRUE ~ NA_character_
  )
  if (is.na(round_col)) return(NA_real_)

  seed_counts <- raw_sims %>%
    filter(seed == target_seed, .data[[round_col]] == 1) %>%
    group_by(sim_id) %>%
    summarise(n = n(), .groups = "drop")

  all_sims <- tibble(sim_id = 1:n_sims)
  seed_counts <- all_sims %>%
    left_join(seed_counts, by = "sim_id") %>%
    mutate(n = coalesce(n, 0L))

  if (grepl("exactly", title, ignore.case = TRUE) || strike_type == "between") {
    mean(seed_counts$n == floor_strike)
  } else {
    mean(seed_counts$n >= floor_strike)
  }
}

compute_highest_seed <- function(raw_sims, floor_strike, cap_strike, strike_type, event_ticker, title, n_sims) {
  round_col <- case_when(
    grepl("-26R32$", event_ticker) ~ "Round_32",
    grepl("-26R16$", event_ticker) ~ "Sweet_16",
    grepl("-26R8$", event_ticker) ~ "Elite_8",
    grepl("-26F4$", event_ticker) ~ "Final_4",
    grepl("-26$", event_ticker) ~ "Champion",
    TRUE ~ NA_character_
  )
  if (is.na(round_col)) return(NA_real_)

  max_seeds <- raw_sims %>%
    filter(.data[[round_col]] == 1) %>%
    group_by(sim_id) %>%
    summarise(max_seed = max(seed), .groups = "drop")

  all_sims <- tibble(sim_id = 1:n_sims)
  max_seeds <- all_sims %>%
    left_join(max_seeds, by = "sim_id") %>%
    mutate(max_seed = coalesce(max_seed, 0L))

  if (strike_type == "between" && !is.na(cap_strike)) {
    mean(max_seeds$max_seed >= floor_strike & max_seeds$max_seed <= cap_strike)
  } else if (grepl("exactly", title, ignore.case = TRUE)) {
    mean(max_seeds$max_seed == floor_strike)
  } else {
    mean(max_seeds$max_seed >= floor_strike)
  }
}

compute_upset_count <- function(raw_sims, floor_strike, cap_strike, strike_type, event_ticker, title, n_sims) {
  # Map round code to sim column and define which seeds are "underdogs"
  # R64 upsets: higher seed (9-16) advancing to R32
  # R32 upsets: in each R32 matchup, the higher seed advancing
  # etc.
  round_col <- case_when(
    grepl("R64", event_ticker) ~ "Round_32",
    grepl("R32", event_ticker) ~ "Sweet_16",
    grepl("R16", event_ticker) ~ "Elite_8",
    grepl("R8", event_ticker) ~ "Final_4",
    TRUE ~ NA_character_
  )
  if (is.na(round_col)) return(NA_real_)

  # In R64: standard matchups are 1v16,8v9,5v12,4v13,6v11,3v14,7v10,2v15
  # An upset = the higher-seeded team (underdog) wins
  # Seeds 9-16 are always underdogs in R64
  if (round_col == "Round_32") {
    upset_seeds <- 9:16
  } else if (round_col == "Sweet_16") {
    # R32 upsets: harder to define without matchup reconstruction
    # Heuristic: seeds 5-16 winning in R32 are generally upsets
    upset_seeds <- 5:16
  } else if (round_col == "Elite_8") {
    upset_seeds <- 3:16
  } else {
    upset_seeds <- 2:16
  }

  upset_counts <- raw_sims %>%
    filter(.data[[round_col]] == 1, seed %in% upset_seeds) %>%
    group_by(sim_id) %>%
    summarise(n_upsets = n(), .groups = "drop")

  all_sims <- tibble(sim_id = 1:n_sims)
  upset_counts <- all_sims %>%
    left_join(upset_counts, by = "sim_id") %>%
    mutate(n_upsets = coalesce(n_upsets, 0L))

  if (grepl("exactly", title, ignore.case = TRUE)) {
    mean(upset_counts$n_upsets == floor_strike)
  } else {
    mean(upset_counts$n_upsets >= floor_strike)
  }
}

compute_seed_win <- function(raw_sims, floor_strike, cap_strike, strike_type, event_ticker, title, n_sims) {
  # All seed win props are about R64 (advancing to Round_32)
  round_col <- "Round_32"

  # Parse target seeds from event ticker
  seeds_match <- str_match(event_ticker, "S(\\d+)$")
  if (!is.na(seeds_match[1, 2])) {
    target_seeds <- as.integer(seeds_match[1, 2])
  } else {
    # Multi-seed: extract numbers from ticker, exclude year prefix
    nums <- as.integer(str_extract_all(event_ticker, "\\d+")[[1]])
    nums <- nums[nums >= 10 & nums <= 16]
    if (length(nums) == 0) return(NA_real_)
    target_seeds <- nums
  }

  wins <- raw_sims %>%
    filter(seed %in% target_seeds, .data[[round_col]] == 1) %>%
    group_by(sim_id) %>%
    summarise(n_wins = n(), .groups = "drop")

  all_sims <- tibble(sim_id = 1:n_sims)
  wins <- all_sims %>%
    left_join(wins, by = "sim_id") %>%
    mutate(n_wins = coalesce(n_wins, 0L))

  if (grepl("exactly", title, ignore.case = TRUE)) {
    mean(wins$n_wins == floor_strike)
  } else {
    mean(wins$n_wins >= floor_strike)
  }
}

# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

fetch_kalshi_edges <- function(nflwork_root, n_sims = 10000) {
  cat("=== Kalshi Edge Finder ===\n")

  # 1. Run bracket simulation
  cat("Running bracket simulation...\n")
  sim_data <- run_bracket_sim(nflwork_root, n_sims)
  team_probs <- sim_data$team_probs
  raw_sims <- sim_data$raw_sims

  # 2. Fetch all Kalshi tournament prop markets
  cat("Fetching Kalshi tournament props...\n")

  all_edges <- tibble()

  # --- Round Advancement (team-level) ---
  round_events <- list(
    list(event = "KXMARMADROUND-26RO32", col = "Round_32"),
    list(event = "KXMARMADROUND-26S16", col = "Sweet_16"),
    list(event = "KXMARMADROUND-26E8", col = "Elite_8"),
    list(event = "KXMARMADROUND-26F4", col = "Final_4"),
    list(event = "KXMARMADROUND-26T2", col = "Title_Game")
  )

  for (re in round_events) {
    markets <- fetch_kalshi_markets(event_ticker = re$event)
    cat(sprintf("  %s: %d markets\n", re$event, length(markets)))
    if (length(markets) > 0) {
      edges <- map_round_advancement(markets, team_probs, re$col)
      if (nrow(edges) > 0) all_edges <- bind_rows(all_edges, edges)
    }
  }

  # --- Seed Sum ---
  seed_sum_markets <- fetch_kalshi_markets(series_ticker = "KXMARMADSEEDSUM")
  cat(sprintf("  KXMARMADSEEDSUM: %d markets\n", length(seed_sum_markets)))
  if (length(seed_sum_markets) > 0) {
    edges <- map_seed_props(seed_sum_markets, raw_sims, "Seed Sum", compute_seed_sum)
    if (nrow(edges) > 0) all_edges <- bind_rows(all_edges, edges)
  }

  # --- Seed Count Per Round ---
  seed_round_markets <- fetch_kalshi_markets(series_ticker = "KXMARMADSEEDROUND")
  cat(sprintf("  KXMARMADSEEDROUND: %d markets\n", length(seed_round_markets)))
  if (length(seed_round_markets) > 0) {
    edges <- map_seed_props(seed_round_markets, raw_sims, "Seed Count", compute_seed_count_in_round)
    if (nrow(edges) > 0) all_edges <- bind_rows(all_edges, edges)
  }

  # --- Highest Seed ---
  highest_seed_markets <- fetch_kalshi_markets(series_ticker = "KXMARMADSEED")
  cat(sprintf("  KXMARMADSEED: %d markets\n", length(highest_seed_markets)))
  if (length(highest_seed_markets) > 0) {
    edges <- map_seed_props(highest_seed_markets, raw_sims, "Highest Seed", compute_highest_seed)
    if (nrow(edges) > 0) all_edges <- bind_rows(all_edges, edges)
  }

  # --- Upset Count ---
  upset_markets <- fetch_kalshi_markets(series_ticker = "KXMARMADUPSET")
  cat(sprintf("  KXMARMADUPSET: %d markets\n", length(upset_markets)))
  if (length(upset_markets) > 0) {
    edges <- map_seed_props(upset_markets, raw_sims, "Upsets", compute_upset_count)
    if (nrow(edges) > 0) all_edges <- bind_rows(all_edges, edges)
  }

  # --- Seed Win Props ---
  seed_win_markets <- fetch_kalshi_markets(series_ticker = "KXMARMADSEEDWIN")
  cat(sprintf("  KXMARMADSEEDWIN: %d markets\n", length(seed_win_markets)))
  if (length(seed_win_markets) > 0) {
    edges <- map_seed_props(seed_win_markets, raw_sims, "Seed Wins", compute_seed_win)
    if (nrow(edges) > 0) all_edges <- bind_rows(all_edges, edges)
  }

  if (nrow(all_edges) == 0) {
    cat("No edges found.\n")
    return(all_edges)
  }

  all_edges <- all_edges %>% arrange(desc(best_ev_pct))

  cat(sprintf("\n%d total markets analyzed, %d with +EV\n",
              nrow(all_edges), sum(all_edges$best_ev_pct > 0)))

  all_edges
}
