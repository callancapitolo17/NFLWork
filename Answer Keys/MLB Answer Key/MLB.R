# MLB Answer Key - Merged pipeline
# Single R process: prepare samples, build F5 predictions, compare to offshore
# Runs in parallel with scrapers via run.py; waits for sentinel before loading scraper data

.t_script_start <- Sys.time()
setwd("~/NFLWork/Answer Keys")
suppressPackageStartupMessages({
  library(data.table)
  library(oddsapiR)
  library(duckdb)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(lubridate)
  library(DBI)
  library(httr)
  library(jsonlite)
})
source("Tools.R")
timer <- pipeline_timer()
startup_secs <- as.numeric(difftime(Sys.time(), .t_script_start, units = "secs"))
timer$mark(sprintf("r_startup (%.1fs total)", startup_secs))

cat("=== MLB ANSWER KEY ===\n")

# =============================================================================
# PHASE 1: LOAD HISTORICAL DATA
# =============================================================================

con <- dbConnect(duckdb(), dbdir = "pbp.duckdb", read_only = TRUE)
on.exit(tryCatch(dbDisconnect(con), error = function(e) NULL), add = TRUE)

DT <- dbGetQuery(con, "
  SELECT *
  FROM mlb_betting_pbp
  WHERE game_home_margin_inning_inning_5 IS NOT NULL
") %>%
  rename(
    # Map inning cumulative columns to Tools.R period convention
    game_home_margin_period_F5 = game_home_margin_inning_inning_5,
    game_total_period_F5       = game_total_inning_inning_5,
    # Full game columns
    game_home_margin_period_FG = home_margin,
    game_total_period_FG       = total_final_score,
    # Consensus columns → standard names for sampling
    home_ml_odds  = consensus_devig_home_odds,
    away_ml_odds  = consensus_devig_away_odds,
    over_odds     = consensus_devig_over_odds,
    under_odds    = consensus_devig_under_odds,
    actual_cover  = home_winner
  ) %>%
  mutate(
    actual_over = ifelse(game_total_period_FG > total_line, 1, 0)
  ) %>%
  as.data.table()

dbDisconnect(con)
on.exit(NULL)  # clear the on.exit since we disconnected

cat(sprintf("Historical data loaded: %d games\n", nrow(DT)))

disp <- compute_dispersion(DT, moneyline = TRUE)
ss <- disp$ss
st <- disp$st
timer$mark("historical_load")

# =============================================================================
# PHASE 2: GET CURRENT GAME ODDS & BUILD CONSENSUS
# =============================================================================

cat("Fetching Odds API...\n")

game_odds <- tryCatch({
  toa_sports_odds(
    sport_key = "baseball_mlb",
    regions = "us,us2,eu",
    markets = "h2h,totals",
    odds_format = "american",
    date_format = "iso"
  )
}, error = function(e) {
  cat(sprintf("Warning: Could not fetch odds: %s\n", e$message))
  return(data.frame())
})

if (nrow(game_odds) == 0) {
  cat("No MLB games found from Odds API. Exiting.\n")
  quit(status = 0)
}

cat(sprintf("API returned %d rows (%d unique games) before filtering.\n",
            nrow(game_odds), n_distinct(game_odds$id)))

if (is.character(game_odds$commence_time)) {
  game_odds$commence_time <- ymd_hms(game_odds$commence_time, tz = "UTC")
}

n_before <- n_distinct(game_odds$id)
game_odds <- game_odds %>%
  filter(commence_time > Sys.time())
n_after <- n_distinct(game_odds$id)
cat(sprintf("Filtered in-progress: %d -> %d games (removed %d).\n",
            n_before, n_after, n_before - n_after))

if (nrow(game_odds) == 0) {
  cat("No upcoming MLB games found. Exiting.\n")
  quit(status = 0)
}

# --- Load sharp scraper data into consensus pool ---
bookmaker_rows <- tryCatch({
  scraper_to_odds_api_format(get_bookmaker_odds("mlb"), game_odds)
}, error = function(e) { cat(sprintf("Bookmaker scraper skip: %s\n", e$message)); data.frame() })

bet105_rows <- tryCatch({
  scraper_to_odds_api_format(get_bet105_odds("mlb"), game_odds)
}, error = function(e) { cat(sprintf("Bet105 scraper skip: %s\n", e$message)); data.frame() })

if (nrow(bookmaker_rows) > 0 || nrow(bet105_rows) > 0) {
  game_odds <- bind_rows(game_odds, bookmaker_rows, bet105_rows)
  cat(sprintf("Consensus pool: %d total rows after adding scrapers.\n", nrow(game_odds)))
}

# --- Build sharp-only weights ---
sharp_names <- names(SHARP_BOOKS)
sharp_weights <- vapply(sharp_names, function(bk) SHARP_BOOKS[[bk]], numeric(1))

all_books <- unique(game_odds$bookmaker_key)
book_weights <- data.frame(bookmaker_key = all_books, stringsAsFactors = FALSE) %>%
  mutate(
    ml_weight     = ifelse(bookmaker_key %in% sharp_names,
                           sharp_weights[match(bookmaker_key, sharp_names)], 0.0),
    totals_weight = ml_weight
  )

# Fallback: if zero sharp books found, use all books equally
if (all(book_weights$ml_weight == 0)) {
  warning("No sharp books in odds data. Falling back to all-book consensus.")
  book_weights$ml_weight <- 1.0
  book_weights$totals_weight <- 1.0
} else {
  # Drop games where NO sharp book posts a line
  sharp_book_keys <- book_weights$bookmaker_key[book_weights$ml_weight > 0]
  games_with_sharp <- game_odds %>%
    filter(bookmaker_key %in% sharp_book_keys) %>%
    pull(id) %>% unique()
  n_before_sharp <- n_distinct(game_odds$id)
  game_odds <- game_odds %>% filter(id %in% games_with_sharp)
  n_dropped <- n_before_sharp - n_distinct(game_odds$id)
  if (n_dropped > 0) {
    cat(sprintf("Dropped %d games with no sharp book coverage.\n", n_dropped))
  }
}

sharp_books_found <- book_weights %>% filter(ml_weight > 0)
cat(sprintf("Sharp consensus: %d sharp books (%s), %d rec books zeroed out.\n",
            nrow(sharp_books_found),
            paste(sharp_books_found$bookmaker_key, collapse = ", "),
            sum(book_weights$ml_weight == 0)))

# MLB consensus: moneyline (not spread) + totals
ml_consensus <- prepare_two_way_odds(
  game_odds    = game_odds,
  mkt_key      = "h2h",
  book_weights = book_weights,
  prob_fun     = devig_american,
  prob_names   = c("prob_home", "prob_away"),
  odds_names   = c("outcomes_price_home", "outcomes_price_away")
) %>%
  moneyline_consensus(
    game_id_col = "id",
    weight_col  = "ml_weight",
    date_col    = "date",
    time_col    = "commence_time",
    market1     = "prob_home",
    market2     = "prob_away"
  )
cat(sprintf("Consensus moneyline: %d games.\n", n_distinct(ml_consensus$id)))

consensus_total <- prepare_two_way_odds(
  game_odds    = game_odds,
  mkt_key      = "totals",
  book_weights = book_weights,
  prob_fun     = devig_american,
  prob_names   = c("prob_over", "prob_under"),
  odds_names   = c("outcomes_price_Over", "outcomes_price_Under")
) %>%
  select(-outcomes_point_Under) %>%
  rename(total_line = outcomes_point_Over) %>%
  pick_consensus_line(
    game_id_col = "id",
    line_col    = "total_line",
    weight_col  = "totals_weight",
    date_col    = "date",
    time_col    = "commence_time",
    market1     = "prob_over",
    market2     = "prob_under"
  )
cat(sprintf("Consensus totals: %d games.\n", n_distinct(consensus_total$id)))

mlb_odds <- ml_consensus %>%
  inner_join(
    consensus_total %>% ungroup() %>% select(-home_team, -away_team, -date, -commence_time),
    by = "id"
  )
cat(sprintf("After inner_join (ML+total): %d games.\n", nrow(mlb_odds)))

mlb_odds <- mlb_odds %>%
  filter(if_all(everything(), ~ !is.na(.)))
cat(sprintf("After NA filter: %d games.\n", nrow(mlb_odds)))
timer$mark("consensus")

# =============================================================================
# PHASE 3: TEAM NAME DICTIONARY
# =============================================================================

# MLB has 30 fixed teams; build dict from Odds API canonical names
cat("Building team name dictionary from Odds API...\n")
odds_api_names <- unique(c(mlb_odds$home_team, mlb_odds$away_team))
team_dict <- data.frame(
  odds_api_name = odds_api_names,
  stringsAsFactors = FALSE
)

con_mlb <- duckdb_connect_retry("mlb.duckdb")
on.exit(tryCatch(dbDisconnect(con_mlb), error = function(e) NULL), add = TRUE)
dbExecute(con_mlb, "DROP TABLE IF EXISTS mlb_team_dict")
dbWriteTable(con_mlb, "mlb_team_dict", team_dict)
cat(sprintf("Stored %d MLB teams in dictionary.\n", nrow(team_dict)))
timer$mark("team_dict")

# =============================================================================
# PHASE 4: GENERATE SAMPLES & FETCH DERIVATIVE ODDS
# =============================================================================

# MLB targets: use ML probability as parent_spread (use_spread_line = FALSE)
targets <- mlb_odds %>%
  transmute(
    id,
    parent_spread = consensus_prob_home,
    parent_total  = total_line,
    target_cover  = consensus_prob_home,
    target_over   = consensus_prob_over
  )

bankroll   <- 100
kelly_mult <- 0.25
dash_db <- file.path(getwd(), "MLB Dashboard", "mlb_dashboard.duckdb")
if (file.exists(dash_db)) {
  tryCatch({
    dash_con <- dbConnect(duckdb(), dbdir = dash_db, read_only = TRUE)
    saved <- dbGetQuery(dash_con, "SELECT param, value FROM sizing_settings")
    dbDisconnect(dash_con)
    if ("bankroll" %in% saved$param) bankroll <- saved$value[saved$param == "bankroll"]
    if ("kelly_mult" %in% saved$param) kelly_mult <- saved$value[saved$param == "kelly_mult"]
  }, error = function(e) NULL)
}
cat(sprintf("Using bankroll=$%.0f, kelly=%.2f\n", bankroll, kelly_mult))
N <- round(nrow(DT) * 0.10, 0)  # 10% sample — validated via parameter sweep

events <- tryCatch({
  get_events("baseball_mlb", regions = "us")
}, error = function(e) {
  cat(sprintf("Warning: Could not fetch events: %s\n", e$message))
  return(data.frame())
})

# F5 derivative markets
all_deriv_markets <- c(
  "h2h_1st_5_innings",
  "spreads_1st_5_innings",
  "totals_1st_5_innings",
  "alternate_totals_1st_5_innings"
)

cat("Generating samples for all games...\n")
samples <- generate_all_samples(
  targets         = targets,
  DT              = DT,
  ss              = ss,
  st              = st,
  N               = N,
  use_spread_line = FALSE   # MLB uses ML probability, not spread line
)
cat(sprintf("Generated %d samples.\n", length(samples)))
timer$mark("sample_gen")

prefetched_odds <- NULL
if (nrow(events) > 0) {
  cat(sprintf("Fetching derivative odds for %d events (%d markets per call)...\n",
              nrow(events), length(all_deriv_markets)))
  prefetched_odds <- fetch_odds_bulk(events$id, all_deriv_markets, "baseball_mlb")
  prefetched_odds <- prefetched_odds[!is.na(prefetched_odds$json_response), ]
  cat(sprintf("Pre-fetched %d/%d event responses.\n", nrow(prefetched_odds), nrow(events)))
}
timer$mark("prefetch_odds")

# =============================================================================
# PHASE 5: BUILD PREDICTIONS (F5 MARKETS)
# =============================================================================

# Read enabled books from dashboard settings
dashboard_db <- file.path(getwd(), "MLB Dashboard", "mlb_dashboard.duckdb")
if (file.exists(dashboard_db)) {
  dash_con <- dbConnect(duckdb(), dbdir = dashboard_db)
  for (book in c("wagerzon", "hoop88", "bfa", "bookmaker", "bet105", "kalshi")) {
    dbExecute(dash_con, "INSERT INTO book_settings (bookmaker_key, enabled)
      VALUES (?, TRUE) ON CONFLICT (bookmaker_key) DO NOTHING", list(book))
  }
  enabled_books <- tryCatch({
    dbGetQuery(dash_con, "SELECT bookmaker_key FROM book_settings WHERE enabled = TRUE")$bookmaker_key
  }, error = function(e) NULL)
  dbDisconnect(dash_con)
  if (length(enabled_books) == 0) enabled_books <- NULL
  cat(sprintf("Enabled books: %s\n",
              if (is.null(enabled_books)) "ALL (none configured)" else paste(enabled_books, collapse = ", ")))
} else {
  enabled_books <- NULL
  cat("No dashboard DB found, using all books.\n")
}

# --- F5 MONEYLINES ---
cat("Building F5 moneyline predictions...\n")
ml_results <- build_moneylines_from_samples(
  samples          = samples,
  consensus_odds   = mlb_odds,
  events           = events,
  periods          = "F5",
  markets          = "h2h_1st_5_innings",
  sport_key        = "baseball_mlb",
  bankroll         = bankroll,
  kelly_mult       = kelly_mult,
  margin_col       = "game_home_margin_period",
  pre_fetched_odds = prefetched_odds,
  books            = enabled_books
)
ml_bets <- ml_results$bets
timer$mark("build_moneylines")

# --- F5 TOTALS ---
cat("Building F5 totals predictions...\n")
total_results <- build_totals_from_samples(
  samples          = samples,
  consensus_odds   = mlb_odds,
  events           = events,
  periods          = c("F5", "F5"),
  markets          = c("totals_1st_5_innings", "alternate_totals_1st_5_innings"),
  sport_key        = "baseball_mlb",
  bankroll         = bankroll,
  kelly_mult       = kelly_mult,
  pre_fetched_odds = prefetched_odds,
  books            = enabled_books
)
total_bets <- total_results$bets
timer$mark("build_totals")

# --- F5 SPREADS ---
cat("Building F5 spread predictions...\n")
spread_results <- build_spreads_from_samples(
  samples          = samples,
  consensus_odds   = mlb_odds,
  events           = events,
  periods          = "F5",
  markets          = "spreads_1st_5_innings",
  sport_key        = "baseball_mlb",
  bankroll         = bankroll,
  kelly_mult       = kelly_mult,
  pre_fetched_odds = prefetched_odds,
  books            = enabled_books
)
spread_bets <- spread_results$bets
timer$mark("build_spreads")

# Combine API bets
EV_THRESHOLD <- 0.02  # 2% — maximizes total profit per parameter sweep

all_bets <- bind_rows(
  ml_bets %>% mutate(market_type = "moneyline"),
  total_bets %>% mutate(market_type = "totals"),
  spread_bets %>% mutate(market_type = "spreads")
) %>%
  filter(ev >= EV_THRESHOLD) %>%
  arrange(desc(ev))

cat("\n=== API BETS SUMMARY ===\n")
all_bets %>%
  group_by(market_type) %>%
  summarise(
    n_bets = n(),
    total_stake = sum(bet_size),
    avg_ev = mean(ev),
    .groups = "drop"
  ) %>%
  print()

# =============================================================================
# PHASE 6: WAIT FOR SCRAPERS & COMPARE TO OFFSHORE
# =============================================================================

sentinel <- file.path(getwd(), ".scrapers_done_mlb")
waited <- 0
while (!file.exists(sentinel) && waited < 120) {
  Sys.sleep(0.5)
  waited <- waited + 0.5
}
if (file.exists(sentinel)) {
  cat(sprintf("Scrapers done (waited %.1fs for sentinel).\n", waited))
} else {
  cat("Warning: Scraper sentinel not found after 120s, proceeding anyway.\n")
}

# Load scraped odds
wagerzon_odds <- tryCatch(get_wagerzon_odds("mlb"), error = function(e) { cat(sprintf("Wagerzon skip: %s\n", e$message)); tibble() })
cat(sprintf("Loaded %d Wagerzon records.\n", nrow(wagerzon_odds)))

hoop88_odds <- tryCatch(get_hoop88_odds("mlb"), error = function(e) { cat(sprintf("Hoop88 skip: %s\n", e$message)); tibble() })
cat(sprintf("Loaded %d Hoop88 records.\n", nrow(hoop88_odds)))

bfa_odds <- tryCatch(get_bfa_odds("mlb"), error = function(e) { cat(sprintf("BFA skip: %s\n", e$message)); tibble() })
cat(sprintf("Loaded %d BFA records.\n", nrow(bfa_odds)))

bookmaker_odds <- tryCatch(get_bookmaker_odds("mlb"), error = function(e) { cat(sprintf("Bookmaker skip: %s\n", e$message)); tibble() })
cat(sprintf("Loaded %d Bookmaker records.\n", nrow(bookmaker_odds)))

bet105_odds <- tryCatch(get_bet105_odds("mlb"), error = function(e) { cat(sprintf("Bet105 skip: %s\n", e$message)); tibble() })
cat(sprintf("Loaded %d Bet105 records.\n", nrow(bet105_odds)))

kalshi_odds <- tryCatch(get_kalshi_odds("mlb"), error = function(e) { cat(sprintf("Kalshi skip: %s\n", e$message)); tibble() })
cat(sprintf("Loaded %d Kalshi records.\n", nrow(kalshi_odds)))
timer$mark("load_scrapers")

# --- WAGERZON ---
if (nrow(wagerzon_odds) > 0) {
  suppressWarnings({
    wz_spread_bets <- compare_spreads_to_wagerzon(
      spread_results, wagerzon_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
    wz_total_bets <- compare_totals_to_wagerzon(
      total_results, wagerzon_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
    wz_ml_bets <- compare_moneylines_to_wagerzon(
      ml_results, wagerzon_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
  })
  wagerzon_bets <- bind_rows(
    wz_spread_bets$bets %>% mutate(market_type = "spreads"),
    wz_total_bets$bets %>% mutate(market_type = "totals"),
    wz_ml_bets$bets %>% mutate(market_type = "moneyline")
  )
  cat(sprintf("Added %d Wagerzon bets.\n", nrow(wagerzon_bets)))
} else {
  wagerzon_bets <- tibble()
}

# --- HOOP88 ---
if (nrow(hoop88_odds) > 0) {
  suppressWarnings({
    h88_spread_bets <- compare_spreads_to_wagerzon(
      spread_results, hoop88_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
    h88_total_bets <- compare_totals_to_wagerzon(
      total_results, hoop88_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
    h88_ml_bets <- compare_moneylines_to_wagerzon(
      ml_results, hoop88_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
  })
  hoop88_bets <- bind_rows(
    h88_spread_bets$bets %>% mutate(market_type = "spreads"),
    h88_total_bets$bets %>% mutate(market_type = "totals"),
    h88_ml_bets$bets %>% mutate(market_type = "moneyline")
  )
  cat(sprintf("Added %d Hoop88 bets.\n", nrow(hoop88_bets)))
} else {
  hoop88_bets <- tibble()
}

# --- BFA ---
if (nrow(bfa_odds) > 0) {
  suppressWarnings({
    bfa_spread_bets <- compare_spreads_to_wagerzon(
      spread_results, bfa_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
    bfa_total_bets <- compare_totals_to_wagerzon(
      total_results, bfa_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
    bfa_ml_bets <- compare_moneylines_to_wagerzon(
      ml_results, bfa_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
  })
  bfa_bets <- bind_rows(
    bfa_spread_bets$bets %>% mutate(market_type = "spreads"),
    bfa_total_bets$bets %>% mutate(market_type = "totals"),
    bfa_ml_bets$bets %>% mutate(market_type = "moneyline")
  )
  cat(sprintf("Added %d BFA bets.\n", nrow(bfa_bets)))

  bfa_alt_bets <- compare_alts_to_samples(
    samples = samples, offshore_odds = bfa_odds,
    consensus_odds = mlb_odds,
    bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
  )
  if (nrow(bfa_alt_bets) > 0) {
    bfa_alt_bets <- bfa_alt_bets %>%
      mutate(market_type = ifelse(grepl("spread", market), "spreads", "totals"))
  }
  cat(sprintf("Added %d BFA alt line bets.\n", nrow(bfa_alt_bets)))
} else {
  bfa_bets <- tibble()
  bfa_alt_bets <- tibble()
}

# --- BOOKMAKER ---
if (nrow(bookmaker_odds) > 0) {
  suppressWarnings({
    bkm_spread_bets <- compare_spreads_to_wagerzon(
      spread_results, bookmaker_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
    bkm_total_bets <- compare_totals_to_wagerzon(
      total_results, bookmaker_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
    bkm_ml_bets <- compare_moneylines_to_wagerzon(
      ml_results, bookmaker_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
  })
  bookmaker_bets <- bind_rows(
    bkm_spread_bets$bets %>% mutate(market_type = "spreads"),
    bkm_total_bets$bets %>% mutate(market_type = "totals"),
    bkm_ml_bets$bets %>% mutate(market_type = "moneyline")
  )
  cat(sprintf("Added %d Bookmaker bets.\n", nrow(bookmaker_bets)))
} else {
  bookmaker_bets <- tibble()
}

# --- BET105 ---
if (nrow(bet105_odds) > 0) {
  suppressWarnings({
    b105_spread_bets <- compare_spreads_to_wagerzon(
      spread_results, bet105_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
    b105_total_bets <- compare_totals_to_wagerzon(
      total_results, bet105_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
    b105_ml_bets <- compare_moneylines_to_wagerzon(
      ml_results, bet105_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
  })
  bet105_bets <- bind_rows(
    b105_spread_bets$bets %>% mutate(market_type = "spreads"),
    b105_total_bets$bets %>% mutate(market_type = "totals"),
    b105_ml_bets$bets %>% mutate(market_type = "moneyline")
  )
  cat(sprintf("Added %d Bet105 bets.\n", nrow(bet105_bets)))

  bet105_alt_bets <- compare_alts_to_samples(
    samples = samples, offshore_odds = bet105_odds,
    consensus_odds = mlb_odds,
    bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
  )
  if (nrow(bet105_alt_bets) > 0) {
    bet105_alt_bets <- bet105_alt_bets %>%
      mutate(market_type = ifelse(grepl("spread", market), "spreads", "totals"))
  }
  cat(sprintf("Added %d Bet105 alt line bets.\n", nrow(bet105_alt_bets)))
} else {
  bet105_bets <- tibble()
  bet105_alt_bets <- tibble()
}

# --- KALSHI ---
if (nrow(kalshi_odds) > 0) {
  suppressWarnings({
    kal_spread_bets <- compare_spreads_to_wagerzon(
      spread_results, kalshi_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
    kal_total_bets <- compare_totals_to_wagerzon(
      total_results, kalshi_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
    kal_ml_bets <- compare_moneylines_to_wagerzon(
      ml_results, kalshi_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = 0.05
    )
  })
  kalshi_bets <- bind_rows(
    kal_spread_bets$bets %>% mutate(market_type = "spreads"),
    kal_total_bets$bets %>% mutate(market_type = "totals"),
    kal_ml_bets$bets %>% mutate(market_type = "moneyline")
  )
  cat(sprintf("Added %d Kalshi bets.\n", nrow(kalshi_bets)))
} else {
  kalshi_bets <- tibble()
}
timer$mark("compare_offshore")

# =============================================================================
# PHASE 7: COMBINE ALL BETS & CORRELATION-ADJUSTED KELLY
# =============================================================================

all_bets_combined <- bind_rows(
  all_bets,
  wagerzon_bets,
  hoop88_bets,
  bfa_bets,
  bfa_alt_bets,
  bookmaker_bets,
  bet105_bets,
  bet105_alt_bets,
  kalshi_bets
) %>%
  { if (!is.null(enabled_books)) filter(., bookmaker_key %in% enabled_books) else . } %>%
  filter(is.na(pt_start_time) | pt_start_time > Sys.time()) %>%
  mutate(base_market = gsub("^alternate_", "", market)) %>%
  group_by(id, base_market, bet_on) %>%
  filter(ev == max(ev)) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(-base_market) %>%
  arrange(desc(ev))

# Load placed bets for correlation awareness
placed_bets <- NULL
if (file.exists(dash_db)) {
  tryCatch({
    dash_con <- dbConnect(duckdb(), dbdir = dash_db, read_only = TRUE)
    now_pacific <- format(Sys.time(), tz = "America/Los_Angeles", usetz = FALSE)
    placed_bets <- dbGetQuery(dash_con, sprintf("
      SELECT game_id, home_team, away_team, game_time, market, bet_on,
             line, model_prob AS prob, model_ev AS ev,
             COALESCE(actual_size, recommended_size) AS bet_size, odds
      FROM placed_bets
      WHERE (game_time > TIMESTAMP '%s' OR game_time IS NULL)
        AND status = 'pending'
    ", now_pacific))
    dbDisconnect(dash_con)
    if (nrow(placed_bets) == 0) placed_bets <- NULL
  }, error = function(e) { placed_bets <<- NULL })
}

all_bets_combined <- adjust_kelly_for_correlation(
  all_bets_combined, samples, bankroll, kelly_mult, placed_bets = placed_bets
)
timer$mark("correlation_adj")

cat("\n=== BETTING SUMMARY ===\n")
all_bets_combined %>%
  group_by(market_type) %>%
  summarise(
    n_bets = n(),
    total_stake = sum(bet_size),
    avg_ev = mean(ev),
    .groups = "drop"
  ) %>%
  arrange(desc(total_stake)) %>%
  print()

cat("\n=== TOP 20 BETS ===\n")
print(all_bets_combined %>% head(20))

# =============================================================================
# PHASE 8: SAVE TO DUCKDB
# =============================================================================

dbExecute(con_mlb, "DROP TABLE IF EXISTS mlb_bets_combined")
dbWriteTable(con_mlb, "mlb_bets_combined", all_bets_combined)
dbDisconnect(con_mlb)
on.exit(NULL)

timer$mark("save_bets")
cat(sprintf("Saved %d bets to mlb_bets_combined table.\n", nrow(all_bets_combined)))
cat("\n=== MLB ANSWER KEY: Complete ===\n")
