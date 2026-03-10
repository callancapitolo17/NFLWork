# CBB Answer Key - Merged pipeline
# Single R process: prepare samples, build predictions, compare to offshore
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
  library(hoopR)
})
source("Tools.R")
timer <- pipeline_timer()
startup_secs <- as.numeric(difftime(Sys.time(), .t_script_start, units = "secs"))
timer$mark(sprintf("r_startup (%.1fs total)", startup_secs))

cat("=== CBB ANSWER KEY ===\n")

# =============================================================================
# PHASE 1: LOAD HISTORICAL DATA
# =============================================================================

con <- dbConnect(duckdb(), dbdir = "cbb.duckdb")

betting_pbp <- dbGetQuery(con, "SELECT * FROM cbb_betting_pbp")

closing_odds <- dbGetQuery(con, "
  SELECT
    id as odds_id,
    home_spread,
    total_line,
    spread_home_odds,
    spread_away_odds,
    tot_over_odds,
    tot_under_odds,
    bookmaker_key
  FROM cbb_closing_odds
  WHERE home_spread IS NOT NULL
    AND total_line IS NOT NULL
")

historical_consensus <- closing_odds %>%
  group_by(odds_id) %>%
  summarize(
    home_spread = median(home_spread, na.rm = TRUE),
    total_line = median(total_line, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    consensus_devig_home_odds = 0.5,
    consensus_devig_away_odds = 0.5,
    consensus_devig_over_odds = 0.5,
    consensus_devig_under_odds = 0.5
  )

DT <- betting_pbp %>%
  inner_join(historical_consensus, by = "odds_id") %>%
  rename(
    home_margin = game_home_margin_fg,
    total_final_score = game_total_fg,
    game_home_margin_period_Half1 = game_home_margin_h1,
    game_home_margin_period_Half2 = game_home_margin_h2,
    game_total_period_Half1 = game_total_h1,
    game_total_period_Half2 = game_total_h2,
    home_score_period_Half1 = home_h1_score,
    home_score_period_Half2 = home_h2_score,
    away_score_period_Half1 = away_h1_score,
    away_score_period_Half2 = away_h2_score,
    home_spread_odds = consensus_devig_home_odds,
    away_spread_odds = consensus_devig_away_odds,
    over_odds = consensus_devig_over_odds,
    under_odds = consensus_devig_under_odds
  ) %>%
  mutate(
    actual_over = ifelse(total_final_score > total_line, 1, 0),
    actual_cover = ifelse(home_margin > -home_spread, 1, 0)
  ) %>%
  as.data.table()

cat(sprintf("Historical data loaded: %d games\n", nrow(DT)))

disp <- compute_dispersion(DT, moneyline = FALSE)
ss <- disp$ss
st <- disp$st

book_weights <- dbGetQuery(con, "
  SELECT DISTINCT bookmaker_key
  FROM cbb_closing_odds
") %>%
  mutate(
    spread_weight = 1.0,
    totals_weight = 1.0
  )

tables <- dbListTables(con)
if (!"cbb_weights" %in% tables) {
  dbWriteTable(con, "cbb_weights", book_weights, overwrite = TRUE)
  cat("Created cbb_weights table with equal weights.\n")
} else {
  book_weights <- dbGetQuery(con, "SELECT * FROM cbb_weights")
}

dbDisconnect(con)
cat("Historical data loaded.\n")
timer$mark("historical_load")

# =============================================================================
# PHASE 2: GET CURRENT GAME ODDS & BUILD CONSENSUS
# =============================================================================

cat("Fetching Odds API...\n")

game_odds <- tryCatch({
  toa_sports_odds(
    sport_key = "basketball_ncaab",
    regions = "us,us2,eu",
    markets = "spreads,totals",
    odds_format = "american",
    date_format = "iso"
  )
}, error = function(e) {
  cat(sprintf("Warning: Could not fetch odds: %s\n", e$message))
  return(data.frame())
})

cat(sprintf("API returned %d rows (%d unique games) before filtering.\n",
            nrow(game_odds), n_distinct(game_odds$id)))

if (is.character(game_odds$commence_time)) {
  game_odds$commence_time <- ymd_hms(game_odds$commence_time, tz = "UTC")
}

n_before <- n_distinct(game_odds$id)
game_odds <- game_odds %>%
  filter(commence_time > Sys.time())
n_after <- n_distinct(game_odds$id)
cat(sprintf("Filtered in-progress: %d -> %d games (removed %d).\n", n_before, n_after, n_before - n_after))

if (nrow(game_odds) == 0) {
  cat("No upcoming games found from Odds API. Exiting.\n")
  quit(status = 0)
}

consensus_spread <- prepare_two_way_odds(
  game_odds    = game_odds,
  mkt_key      = "spreads",
  book_weights = book_weights,
  prob_fun     = devig_american,
  prob_names   = c("prob_home", "prob_away"),
  odds_names   = c("outcomes_price_home", "outcomes_price_away")
) %>%
  select(-outcomes_point_away) %>%
  rename(spread = outcomes_point_home) %>%
  pick_consensus_line(
    game_id_col = "id",
    line_col    = "spread",
    weight_col  = "spread_weight",
    date_col    = "date",
    time_col    = "commence_time",
    market1     = "prob_home",
    market2     = "prob_away"
  )
cat(sprintf("Consensus spreads: %d games.\n", n_distinct(consensus_spread$id)))

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

cbb_odds <- consensus_spread %>%
  inner_join(
    consensus_total %>% ungroup() %>% select(-home_team, -away_team, -date, -commence_time),
    by = "id"
  )
cat(sprintf("After inner_join (spread+total): %d games.\n", nrow(cbb_odds)))

cbb_odds <- cbb_odds %>%
  filter(if_all(everything(), ~ !is.na(.)))
cat(sprintf("After NA filter: %d games.\n", nrow(cbb_odds)))
timer$mark("consensus")

# =============================================================================
# PHASE 3: BUILD TEAM NAME DICTIONARY
# =============================================================================

cat("Building team name dictionary from ESPN...\n")
espn_teams <- tryCatch({
  espn_mbb_teams() %>%
    select(team_id, abbreviation, display_name, short_name, nickname, mascot)
}, error = function(e) {
  cat(sprintf("Warning: Could not fetch ESPN teams: %s\n", e$message))
  data.frame()
})

if (nrow(espn_teams) > 0) {
  odds_api_names <- unique(c(cbb_odds$home_team, cbb_odds$away_team))

  team_dict <- espn_teams %>%
    mutate(odds_api_name = NA_character_)

  for (oa_name in odds_api_names) {
    idx <- which(tolower(team_dict$display_name) == tolower(oa_name))
    if (length(idx) == 0) {
      idx <- which(sapply(seq_len(nrow(team_dict)), function(i) {
        grepl(team_dict$nickname[i], oa_name, ignore.case = TRUE) &&
          grepl(team_dict$mascot[i], oa_name, ignore.case = TRUE)
      }))
    }
    if (length(idx) == 1) {
      team_dict$odds_api_name[idx] <- oa_name
    }
  }

  matched <- sum(!is.na(team_dict$odds_api_name))
  cat(sprintf("Mapped %d/%d Odds API teams to ESPN dictionary.\n", matched, length(odds_api_names)))

  con_dict <- dbConnect(duckdb(), dbdir = "cbb.duckdb")
  dbWriteTable(con_dict, "cbb_team_dict", team_dict, overwrite = TRUE)
  dbDisconnect(con_dict)
}
timer$mark("espn_teams")

# =============================================================================
# PHASE 4: GENERATE SAMPLES & FETCH DERIVATIVE ODDS
# =============================================================================

targets <- cbb_odds %>%
  transmute(
    id,
    parent_spread = spread,
    parent_total  = total_line,
    target_cover  = consensus_prob_home,
    target_over   = consensus_prob_over
  )

bankroll   <- 100
kelly_mult <- 0.25
dash_db <- file.path(getwd(), "CBB Dashboard", "cbb_dashboard.duckdb")
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
N <- round(nrow(DT) * 0.02, 0)

events <- tryCatch({
  get_events("basketball_ncaab", regions = "us")
}, error = function(e) {
  cat(sprintf("Warning: Could not fetch events: %s\n", e$message))
  return(data.frame())
})

all_deriv_markets <- c(
  "h2h_h1", "h2h_h2",
  "spreads_h1", "spreads_h2", "alternate_spreads_h1", "alternate_spreads_h2",
  "totals_h1", "totals_h2", "alternate_totals_h1", "alternate_totals_h2",
  "team_totals_h1", "team_totals_h2", "alternate_team_totals_h1", "alternate_team_totals_h2"
)

cat("Generating samples for all games...\n")
samples <- generate_all_samples(
  targets         = targets,
  DT              = DT,
  ss              = ss,
  st              = st,
  N               = N,
  use_spread_line = TRUE
)
cat(sprintf("Generated %d samples.\n", length(samples)))
timer$mark("sample_gen")

prefetched_odds <- NULL
if (nrow(events) > 0) {
  cat(sprintf("Fetching derivative odds for %d events (%d markets per call)...\n",
              nrow(events), length(all_deriv_markets)))
  prefetched_odds <- fetch_odds_bulk(events$id, all_deriv_markets, "basketball_ncaab")
  prefetched_odds <- prefetched_odds[!is.na(prefetched_odds$json_response), ]
  cat(sprintf("Pre-fetched %d/%d event responses.\n", nrow(prefetched_odds), nrow(events)))
}
timer$mark("prefetch_odds")

# =============================================================================
# PHASE 5: BUILD PREDICTIONS (API BETS)
# =============================================================================

# Read enabled books from dashboard settings
dashboard_db <- file.path(getwd(), "CBB Dashboard", "cbb_dashboard.duckdb")
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
  cat(sprintf("Enabled books: %s\n", if (is.null(enabled_books)) "ALL (none configured)" else paste(enabled_books, collapse = ", ")))
} else {
  enabled_books <- NULL
  cat("No dashboard DB found, using all books.\n")
}

periods_base <- c("Half1", "Half2")

# --- MONEYLINES ---
cat("Building moneyline predictions...\n")
ml_results <- build_moneylines_from_samples(
  samples          = samples,
  consensus_odds   = cbb_odds,
  events           = events,
  periods          = periods_base,
  markets          = c("h2h_h1", "h2h_h2"),
  sport_key        = "basketball_ncaab",
  bankroll         = bankroll,
  kelly_mult       = kelly_mult,
  margin_col       = "game_home_margin_period",
  pre_fetched_odds = prefetched_odds,
  books            = enabled_books
)
ml_bets <- ml_results$bets
timer$mark("build_moneylines")

# --- TOTALS + ALTERNATE TOTALS ---
cat("Building totals predictions...\n")
totals_markets <- c(
  "totals_h1", "totals_h2",
  "alternate_totals_h1", "alternate_totals_h2"
)
totals_periods <- c("Half1", "Half2", "Half1", "Half2")

total_results <- build_totals_from_samples(
  samples          = samples,
  consensus_odds   = cbb_odds,
  events           = events,
  periods          = totals_periods,
  markets          = totals_markets,
  sport_key        = "basketball_ncaab",
  bankroll         = bankroll,
  kelly_mult       = kelly_mult,
  pre_fetched_odds = prefetched_odds,
  books            = enabled_books
)
total_bets <- total_results$bets
timer$mark("build_totals")

# --- SPREADS + ALTERNATE SPREADS ---
cat("Building spread predictions...\n")
spreads_markets <- c(
  "spreads_h1", "spreads_h2",
  "alternate_spreads_h1", "alternate_spreads_h2"
)
spreads_periods <- c("Half1", "Half2", "Half1", "Half2")

spread_results <- build_spreads_from_samples(
  samples          = samples,
  consensus_odds   = cbb_odds,
  events           = events,
  periods          = spreads_periods,
  markets          = spreads_markets,
  sport_key        = "basketball_ncaab",
  bankroll         = bankroll,
  kelly_mult       = kelly_mult,
  pre_fetched_odds = prefetched_odds,
  books            = enabled_books
)
spread_bets <- spread_results$bets
timer$mark("build_spreads")

# --- TEAM TOTALS + ALTERNATE TEAM TOTALS ---
cat("Building team totals predictions...\n")
team_totals_markets <- c(
  "team_totals_h1", "team_totals_h2",
  "alternate_team_totals_h1", "alternate_team_totals_h2"
)
team_totals_periods <- c("Half1", "Half2", "Half1", "Half2")

team_totals_results <- build_team_totals_from_samples(
  samples          = samples,
  consensus_odds   = cbb_odds,
  events           = events,
  periods          = team_totals_periods,
  markets          = team_totals_markets,
  sport_key        = "basketball_ncaab",
  bankroll         = bankroll,
  kelly_mult       = kelly_mult,
  pre_fetched_odds = prefetched_odds,
  books            = enabled_books
)
team_totals_bets <- team_totals_results$bets
timer$mark("build_team_totals")

# Combine API bets
EV_THRESHOLD <- 0.05

all_bets <- bind_rows(
  ml_bets %>% mutate(market_type = "moneyline"),
  total_bets %>% mutate(market_type = "totals"),
  spread_bets %>% mutate(market_type = "spreads"),
  team_totals_bets %>% mutate(market_type = "team_totals")
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

# Wait for run.py to signal that all scrapers are complete
sentinel <- file.path(getwd(), ".scrapers_done_cbb")
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
wagerzon_odds <- get_wagerzon_odds("cbb")
cat(sprintf("Loaded %d Wagerzon records.\n", nrow(wagerzon_odds)))

hoop88_odds <- get_hoop88_odds("cbb")
cat(sprintf("Loaded %d Hoop88 records.\n", nrow(hoop88_odds)))

bfa_odds <- get_bfa_odds("cbb")
cat(sprintf("Loaded %d BFA records.\n", nrow(bfa_odds)))

bookmaker_odds <- get_bookmaker_odds("cbb")
cat(sprintf("Loaded %d Bookmaker records.\n", nrow(bookmaker_odds)))

bet105_odds <- get_bet105_odds("cbb")
cat(sprintf("Loaded %d Bet105 records.\n", nrow(bet105_odds)))

kalshi_odds <- get_kalshi_odds("cbb")
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
  cat(sprintf("Added %d Wagerzon bets to predictions.\n", nrow(wagerzon_bets)))
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
  cat(sprintf("Added %d Hoop88 bets to predictions.\n", nrow(hoop88_bets)))
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
  cat(sprintf("Added %d BFA bets to predictions.\n", nrow(bfa_bets)))

  bfa_alt_bets <- compare_alts_to_samples(
    samples = samples,
    offshore_odds = bfa_odds,
    consensus_odds = cbb_odds,
    bankroll = bankroll,
    kelly_mult = kelly_mult,
    ev_threshold = 0.05
  )
  if (nrow(bfa_alt_bets) > 0) {
    bfa_alt_bets <- bfa_alt_bets %>%
      mutate(market_type = ifelse(grepl("spread", market), "spreads", "totals"))
  }
  cat(sprintf("Added %d BFA alt line bets from samples.\n", nrow(bfa_alt_bets)))
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
  cat(sprintf("Added %d Bookmaker bets to predictions.\n", nrow(bookmaker_bets)))
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
  cat(sprintf("Added %d Bet105 bets to predictions.\n", nrow(bet105_bets)))

  bet105_alt_bets <- compare_alts_to_samples(
    samples = samples,
    offshore_odds = bet105_odds,
    consensus_odds = cbb_odds,
    bankroll = bankroll,
    kelly_mult = kelly_mult,
    ev_threshold = 0.05
  )
  if (nrow(bet105_alt_bets) > 0) {
    bet105_alt_bets <- bet105_alt_bets %>%
      mutate(market_type = ifelse(grepl("spread", market), "spreads", "totals"))
  }
  cat(sprintf("Added %d Bet105 alt line bets from samples.\n", nrow(bet105_alt_bets)))
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
  cat(sprintf("Added %d Kalshi bets to predictions.\n", nrow(kalshi_bets)))
} else {
  kalshi_bets <- tibble()
}
timer$mark("compare_offshore")

# =============================================================================
# PHASE 7: COMBINE ALL BETS & SAVE
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
  filter(is.na(pt_start_time) | pt_start_time > Sys.time()) %>%
  mutate(base_market = gsub("^alternate_", "", market)) %>%
  group_by(id, base_market, bet_on) %>%
  filter(ev == max(ev)) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(-base_market) %>%
  arrange(desc(ev))

# =============================================================================
# PHASE 7.5: CORRELATION-ADJUSTED KELLY SIZING
# =============================================================================

# Load already-placed bets for correlation awareness (Cases 2-5)
placed_bets <- NULL
if (file.exists(dash_db)) {
  tryCatch({
    dash_con <- dbConnect(duckdb(), dbdir = dash_db, read_only = TRUE)
    # game_time is stored as Pacific time (naive timestamp), so compare against Pacific now
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

con <- dbConnect(duckdb(), dbdir = "cbb.duckdb")
dbExecute(con, "DROP TABLE IF EXISTS cbb_bets_combined")
dbWriteTable(con, "cbb_bets_combined", all_bets_combined)
timer$mark("save_bets")
dbDisconnect(con)

cat(sprintf("Saved %d bets to cbb_bets_combined table.\n", nrow(all_bets_combined)))
cat("\n=== CBB ANSWER KEY: Complete ===\n")
