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
  library(digest)
})
source("Tools.R")
source("triple_play_helpers.R")
source("MLB Answer Key/odds_screen.R")
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
    game_home_margin_period_F3 = game_home_margin_inning_inning_3,
    game_total_period_F3       = game_total_inning_inning_3,
    game_home_margin_period_F5 = game_home_margin_inning_inning_5,
    game_total_period_F5       = game_total_inning_inning_5,
    game_home_margin_period_F7 = game_home_margin_inning_inning_7,
    game_total_period_F7       = game_total_inning_inning_7,
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

# Compute per-team 1st-inning scoring indicators for SCR 1ST prop pricing.
# Carried through sampling by preserving as extra columns on DT.
inning_1 <- determine_inning_1_scoring_vec(
  m1 = DT$game_home_margin_inning_inning_1,
  t1 = DT$game_total_inning_inning_1
)
DT[, home_scored_in_1st := inning_1$home]
DT[, away_scored_in_1st := inning_1$away]
cat(sprintf("inning-1 scoring coverage: %.1f%% of games have determinable inning 1.\n",
            100 * mean(!is.na(DT$home_scored_in_1st))))

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
# Filter to FG-only markets: Phase 2 consensus uses full-game h2h + totals.
# Including F5 markets (h2h_f5, totals_f5) creates duplicate rows in the
# pivot that crash devig_american().
bookmaker_rows <- tryCatch({
  raw <- get_bookmaker_odds("mlb")
  if (nrow(raw) > 0) raw <- raw %>% filter(period %in% c("fg", "Full", "full"))
  scraper_to_odds_api_format(raw, game_odds)
}, error = function(e) { cat(sprintf("Bookmaker scraper skip: %s\n", e$message)); data.frame() })

bet105_rows <- tryCatch({
  raw <- get_bet105_odds("mlb")
  if (nrow(raw) > 0) raw <- raw %>% filter(period %in% c("fg", "Full", "full"))
  scraper_to_odds_api_format(raw, game_odds)
}, error = function(e) { cat(sprintf("Bet105 scraper skip: %s\n", e$message)); data.frame() })

if (nrow(bookmaker_rows) > 0 || nrow(bet105_rows) > 0) {
  game_odds <- bind_rows(game_odds, bookmaker_rows, bet105_rows) %>%
    # Deduplicate: doubleheaders cause many-to-many joins in scraper_to_odds_api_format,
    # creating duplicate (id, bookmaker, market, outcome) rows that crash the pivot.
    distinct(id, bookmaker_key, market_key, outcomes_name, .keep_all = TRUE)
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

# Tools.R build_*_from_samples() functions expect a 'spread' column.
# MLB doesn't have a FG spread — use consensus_prob_home as the equivalent
# (it's used for joining/filtering, not for spread evaluation).
mlb_odds <- mlb_odds %>%
  mutate(spread = consensus_prob_home) %>%
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

# Derivative markets fetched per-event via fetch_odds_bulk().
# Covers FG mains + alts, F3 mains, F5 mains + alts, F7 mains.
# NOTE: do NOT add these FG mains to the Phase 2 toa_sports_odds() call above —
# mixing FG and F-period markets in that pivot crashes devig_american().
all_deriv_markets <- c(
  # Full game
  "h2h", "totals", "spreads",
  "alternate_totals", "alternate_spreads",
  # First 3 innings
  "h2h_1st_3_innings", "totals_1st_3_innings", "spreads_1st_3_innings",
  # First 5 innings
  "h2h_1st_5_innings", "totals_1st_5_innings", "spreads_1st_5_innings",
  "alternate_totals_1st_5_innings", "alternate_spreads_1st_5_innings",
  # First 7 innings
  "h2h_1st_7_innings", "totals_1st_7_innings", "spreads_1st_7_innings"
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

# --- Store samples for downstream tools (parlay edge finder) ---
cat("Storing samples to DuckDB for parlay analysis...\n")
sample_rows <- imap_dfr(samples, function(s, game_id) {
  samp <- s$sample
  tibble(
    game_id            = game_id,
    sim_idx            = seq_len(nrow(samp)),
    home_margin        = samp$game_home_margin_period_FG,
    total_final_score  = samp$game_total_period_FG,
    home_margin_f3     = samp$game_home_margin_period_F3,
    total_f3           = samp$game_total_period_F3,
    home_margin_f5     = samp$game_home_margin_period_F5,
    total_f5           = samp$game_total_period_F5,
    home_margin_f7     = samp$game_home_margin_period_F7,
    total_f7           = samp$game_total_period_F7,
    home_scored_in_1st = samp$home_scored_in_1st,
    away_scored_in_1st = samp$away_scored_in_1st
  )
})

# --- Export MM tables to separate DB (avoids lock contention with RFQ bot) ---
# The RFQ bot reads mlb_mm.duckdb on a tight cycle; using a separate file
# means the pipeline's write lock on mlb.duckdb never blocks the bot.
con_mm <- duckdb_connect_retry("mlb_mm.duckdb")
on.exit(tryCatch(dbDisconnect(con_mm), error = function(e) NULL), add = TRUE)
dbExecute(con_mm, "DROP TABLE IF EXISTS mlb_game_samples")
dbWriteTable(con_mm, "mlb_game_samples", sample_rows)

dbExecute(con_mm, "DROP TABLE IF EXISTS mlb_samples_meta")
dbExecute(con_mm, "CREATE TABLE mlb_samples_meta (generated_at TIMESTAMP)")
dbExecute(con_mm, sprintf(
  "INSERT INTO mlb_samples_meta VALUES ('%s')",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
))
dbDisconnect(con_mm)

# Store consensus for game matching (id ↔ team names)
# Also stored as mlb_odds_temp so resolve_offshore_teams() Layer 2 can
# resolve Wagerzon team names against canonical Odds API names.
consensus_export <- mlb_odds %>%
  select(id, home_team, away_team, total_line, consensus_prob_home, commence_time)
dbExecute(con_mlb, "DROP TABLE IF EXISTS mlb_consensus_temp")
dbWriteTable(con_mlb, "mlb_consensus_temp", consensus_export)
dbExecute(con_mlb, "DROP TABLE IF EXISTS mlb_odds_temp")
dbWriteTable(con_mlb, "mlb_odds_temp", consensus_export)

cat(sprintf("Stored %d samples (%d games) + consensus for parlay tool.\n",
            nrow(sample_rows), length(unique(sample_rows$game_id))))

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

# Log which Odds API books returned F5 data
api_f5_books <- unique(c(
  if (nrow(ml_bets) > 0) ml_bets$bookmaker_key else character(0),
  if (nrow(total_bets) > 0) total_bets$bookmaker_key else character(0),
  if (nrow(spread_bets) > 0) spread_bets$bookmaker_key else character(0)
))
if (length(api_f5_books) > 0) {
  cat(sprintf("Odds API F5 books: %s\n", paste(api_f5_books, collapse = ", ")))
} else {
  cat("Warning: No Odds API books returned F5 data.\n")
}

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

# Map scraper market names to standardized period convention.
# Scrapers use short names (spreads, h2h) + period column;
# predictions use Odds API names (spreads_1st_5_innings, h2h_1st_5_innings).
# compare_alts_to_samples uses suffixes like _f3, _f5, _f7 for period lookup.
map_scraper_markets_mlb <- function(odds) {
  if (nrow(odds) == 0) return(odds)
  odds %>% mutate(market = case_when(
    # --- F3 (first 3 innings) ---
    period == "F3" & market == "spreads"  ~ "spreads_1st_3_innings",
    period == "F3" & market == "totals"   ~ "totals_1st_3_innings",
    period == "F3" & market == "h2h"      ~ "h2h_1st_3_innings",
    market == "spreads_f3"                ~ "spreads_1st_3_innings",
    market == "totals_f3"                 ~ "totals_1st_3_innings",
    market == "h2h_f3"                    ~ "h2h_1st_3_innings",
    # --- F5 (first 5 innings) ---
    period == "F5" & market == "spreads"  ~ "spreads_1st_5_innings",
    period == "F5" & market == "totals"   ~ "totals_1st_5_innings",
    period == "F5" & market == "h2h"      ~ "h2h_1st_5_innings",
    market == "spreads_f5"                ~ "spreads_1st_5_innings",
    market == "totals_f5"                 ~ "totals_1st_5_innings",
    market == "h2h_f5"                    ~ "h2h_1st_5_innings",
    # _h1 suffix — in baseball "1st half" = first 5 innings (Wagerzon, BFA)
    market == "spreads_h1"                ~ "spreads_1st_5_innings",
    market == "totals_h1"                 ~ "totals_1st_5_innings",
    market == "h2h_h1"                    ~ "h2h_1st_5_innings",
    # --- F7 (first 7 innings) ---
    market == "spreads_f7"                ~ "spreads_1st_7_innings",
    market == "totals_f7"                 ~ "totals_1st_7_innings",
    market == "h2h_f7"                    ~ "h2h_1st_7_innings",
    # --- Alt lines: remap _h1 → _f5 so compare_alts_to_samples period lookup works ---
    # (BFA/Bet105 use "alternate_totals_h1" but MLB samples use "F5" not "Half1")
    market == "alternate_totals_h1"       ~ "alternate_totals_f5",
    market == "alternate_spreads_h1"      ~ "alternate_spreads_f5",
    TRUE ~ market
  ))
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

# Map all scraper markets to Odds API F5 convention
wagerzon_odds  <- map_scraper_markets_mlb(wagerzon_odds)
hoop88_odds    <- map_scraper_markets_mlb(hoop88_odds)
bfa_odds       <- map_scraper_markets_mlb(bfa_odds)
bookmaker_odds <- map_scraper_markets_mlb(bookmaker_odds)
bet105_odds    <- map_scraper_markets_mlb(bet105_odds)
kalshi_odds    <- map_scraper_markets_mlb(kalshi_odds)
timer$mark("load_scrapers")

# --- WAGERZON ---
if (nrow(wagerzon_odds) > 0) {
  suppressWarnings({
    wz_spread_bets <- compare_spreads_to_wagerzon(
      spread_results, wagerzon_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = EV_THRESHOLD
    )
    wz_total_bets <- compare_totals_to_wagerzon(
      total_results, wagerzon_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = EV_THRESHOLD
    )
    wz_ml_bets <- compare_moneylines_to_wagerzon(
      ml_results, wagerzon_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = EV_THRESHOLD
    )
  })
  wagerzon_bets <- bind_rows(
    wz_spread_bets$bets %>% mutate(market_type = "spreads"),
    wz_total_bets$bets %>% mutate(market_type = "totals"),
    wz_ml_bets$bets %>% mutate(market_type = "moneyline")
  )
  cat(sprintf("Added %d Wagerzon bets.\n", nrow(wagerzon_bets)))

  wz_alt_bets <- compare_alts_to_samples(
    samples = samples, offshore_odds = wagerzon_odds,
    consensus_odds = mlb_odds,
    bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = EV_THRESHOLD
  )
  if (nrow(wz_alt_bets) > 0) {
    wz_alt_bets <- wz_alt_bets %>%
      mutate(market_type = ifelse(grepl("spread", market), "spreads", "totals"))
  }
  cat(sprintf("Added %d Wagerzon derivative/alt bets.\n", nrow(wz_alt_bets)))
} else {
  wagerzon_bets <- tibble()
  wz_alt_bets <- tibble()
}

# --- HOOP88 ---
if (nrow(hoop88_odds) > 0) {
  suppressWarnings({
    h88_spread_bets <- compare_spreads_to_wagerzon(
      spread_results, hoop88_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = EV_THRESHOLD
    )
    h88_total_bets <- compare_totals_to_wagerzon(
      total_results, hoop88_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = EV_THRESHOLD
    )
    h88_ml_bets <- compare_moneylines_to_wagerzon(
      ml_results, hoop88_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = EV_THRESHOLD
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
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = EV_THRESHOLD
    )
    bfa_total_bets <- compare_totals_to_wagerzon(
      total_results, bfa_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = EV_THRESHOLD
    )
    bfa_ml_bets <- compare_moneylines_to_wagerzon(
      ml_results, bfa_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = EV_THRESHOLD
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
    bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = EV_THRESHOLD
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
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = EV_THRESHOLD
    )
    bkm_total_bets <- compare_totals_to_wagerzon(
      total_results, bookmaker_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = EV_THRESHOLD
    )
    bkm_ml_bets <- compare_moneylines_to_wagerzon(
      ml_results, bookmaker_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = EV_THRESHOLD
    )
  })
  bookmaker_bets <- bind_rows(
    bkm_spread_bets$bets %>% mutate(market_type = "spreads"),
    bkm_total_bets$bets %>% mutate(market_type = "totals"),
    bkm_ml_bets$bets %>% mutate(market_type = "moneyline")
  )
  cat(sprintf("Added %d Bookmaker bets.\n", nrow(bookmaker_bets)))

  bkm_alt_bets <- compare_alts_to_samples(
    samples = samples, offshore_odds = bookmaker_odds,
    consensus_odds = mlb_odds,
    bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = EV_THRESHOLD
  )
  if (nrow(bkm_alt_bets) > 0) {
    bkm_alt_bets <- bkm_alt_bets %>%
      mutate(market_type = ifelse(grepl("spread", market), "spreads", "totals"))
  }
  cat(sprintf("Added %d Bookmaker derivative/alt bets.\n", nrow(bkm_alt_bets)))
} else {
  bookmaker_bets <- tibble()
  bkm_alt_bets <- tibble()
}

# --- BET105 ---
if (nrow(bet105_odds) > 0) {
  suppressWarnings({
    b105_spread_bets <- compare_spreads_to_wagerzon(
      spread_results, bet105_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = EV_THRESHOLD
    )
    b105_total_bets <- compare_totals_to_wagerzon(
      total_results, bet105_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = EV_THRESHOLD
    )
    b105_ml_bets <- compare_moneylines_to_wagerzon(
      ml_results, bet105_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = EV_THRESHOLD
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
    bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = EV_THRESHOLD
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
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = EV_THRESHOLD
    )
    kal_total_bets <- compare_totals_to_wagerzon(
      total_results, kalshi_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = EV_THRESHOLD
    )
    kal_ml_bets <- compare_moneylines_to_wagerzon(
      ml_results, kalshi_odds,
      bankroll = bankroll, kelly_mult = kelly_mult, ev_threshold = EV_THRESHOLD
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
  wz_alt_bets,
  hoop88_bets,
  bfa_bets,
  bfa_alt_bets,
  bookmaker_bets,
  bkm_alt_bets,
  bet105_bets,
  bet105_alt_bets,
  kalshi_bets
) %>%
  { if (!is.null(enabled_books)) filter(., bookmaker_key %in% enabled_books) else . } %>%
  filter(is.na(pt_start_time) | pt_start_time > Sys.time()) %>%
  # bet_row_id: stable hash of (game_id, market, line, bet_on) — shared key
  # across mlb_bets_combined and mlb_bets_book_prices. Excludes bookmaker so
  # all books for the same bet line/side share an id.
  mutate(bet_row_id = vapply(
    paste(id, market, ifelse(is.na(line), "", as.character(line)), bet_on, sep = "|"),
    function(s) digest::digest(s, algo = "md5"),
    character(1)
  )) %>%
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
# PHASE 7b: BUILD mlb_bets_book_prices (per-book pill odds)
# =============================================================================

# Game id lookup: join scraper wide frames to Odds API game IDs via home/away team.
.game_id_lookup <- mlb_odds %>%
  select(id, home_team, away_team) %>%
  distinct()

# Convert a wide-format scraper frame (wagerzon/hoop88/bfa/bookmaker/bet105)
# to canonical shape for normalize_book_odds_frame:
#   game_id, market_name, bet_on, line, american_odds, fetch_time.
# The scraper frames are wide (odds_home, odds_away, odds_over, odds_under).
# We join to .game_id_lookup on home_team/away_team to get the Odds API id,
# then pivot to long format.
.scraper_to_canonical <- function(raw) {
  if (is.null(raw) || nrow(raw) == 0) return(NULL)
  if (!"fetch_time" %in% names(raw)) raw$fetch_time <- as.POSIXct(NA, tz = "UTC")

  # Join to Odds API game IDs. Scraper frames have home_team + away_team
  # already resolved to canonical names by resolve_offshore_teams().
  joined <- raw %>%
    inner_join(.game_id_lookup, by = c("home_team", "away_team")) %>%
    rename(game_id = id)

  if (nrow(joined) == 0) return(NULL)

  rows <- list()

  for (i in seq_len(nrow(joined))) {
    row <- joined[i, , drop = FALSE]
    mkt <- row$market
    ft  <- row$fetch_time
    gid <- row$game_id

    # Spread side (home / away)
    if (!is.na(row$odds_home) && !is.na(row$line)) {
      rows[[length(rows) + 1]] <- tibble(
        game_id = gid, market_name = mkt,
        bet_on = row$home_team, line = as.numeric(row$line),
        american_odds = as.integer(row$odds_home), fetch_time = ft
      )
    }
    if (!is.na(row$odds_away) && !is.na(row$line)) {
      rows[[length(rows) + 1]] <- tibble(
        game_id = gid, market_name = mkt,
        bet_on = row$away_team,
        line = -as.numeric(row$line),   # away spread is negated home spread
        american_odds = as.integer(row$odds_away), fetch_time = ft
      )
    }
    # Totals (Over / Under) — scraper rows are one-market-per-row, so the
    # row's `market` is already a totals_* string when odds_over/under populate.
    # (The previous gsub("spreads","totals", mkt) was a misleading no-op.)
    if (!is.na(row$odds_over) && !is.na(row$line)) {
      rows[[length(rows) + 1]] <- tibble(
        game_id = gid, market_name = mkt,
        bet_on = "Over", line = as.numeric(row$line),
        american_odds = as.integer(row$odds_over), fetch_time = ft
      )
    }
    if (!is.na(row$odds_under) && !is.na(row$line)) {
      rows[[length(rows) + 1]] <- tibble(
        game_id = gid, market_name = mkt,
        bet_on = "Under", line = as.numeric(row$line),
        american_odds = as.integer(row$odds_under), fetch_time = ft
      )
    }
    # Moneyline (h2h) — odds_home/odds_away with no line. Row's `market`
    # is already an h2h_* string by the time we hit this branch.
    if (!is.na(row$odds_home) && is.na(row$line)) {
      rows[[length(rows) + 1]] <- tibble(
        game_id = gid, market_name = mkt,
        bet_on = row$home_team, line = NA_real_,
        american_odds = as.integer(row$odds_home), fetch_time = ft
      )
    }
    if (!is.na(row$odds_away) && is.na(row$line)) {
      rows[[length(rows) + 1]] <- tibble(
        game_id = gid, market_name = mkt,
        bet_on = row$away_team, line = NA_real_,
        american_odds = as.integer(row$odds_away), fetch_time = ft
      )
    }
  }

  if (length(rows) == 0) return(NULL)
  result <- bind_rows(rows)
  normalize_book_odds_frame(result)
}

# Convert an Odds API long-format frame (already one row per outcome) to
# canonical shape. Expected columns: id, bookmaker_key, market_key,
# outcomes_name, outcomes_price, outcomes_point. Optionally honors a
# `fetch_time` column if the caller has already derived one (e.g., from the
# bookmaker `last_update` field); falls back to Sys.time() when absent so
# downstream staleness chips don't render "always stale."
.odds_api_to_canonical <- function(raw) {
  if (is.null(raw) || nrow(raw) == 0) return(NULL)
  ft <- if ("fetch_time" %in% names(raw)) raw$fetch_time else rep(Sys.time(), nrow(raw))
  # Replace any NA fetch_times with now (NA breaks staleness chip rendering)
  ft <- as.POSIXct(ft, tz = "UTC")
  ft[is.na(ft)] <- Sys.time()

  result <- raw %>%
    transmute(
      game_id      = id,
      market_name  = market_key,
      bet_on       = outcomes_name,
      line         = as.numeric(outcomes_point),
      american_odds = as.integer(outcomes_price),
      fetch_time   = ft
    )
  normalize_book_odds_frame(result)
}

# Small NULL-coalesce helper (matches rlang's %||% behaviour without the dep).
# Defined before the parser so its body resolves cleanly even on first call.
`%||%` <- function(x, y) if (is.null(x)) y else x

# Parse the prefetched_odds frame (event_id + raw JSON response from the
# Odds API /events/{id}/odds endpoint) into long format matching
# .odds_api_to_canonical's expected column shape.
#
# This is the F-period and alt-market coverage path for DK/FD/Pinnacle.
# game_odds (the canonical h2h+totals frame consumed by Phase 2 devigging)
# is INTENTIONALLY narrow — widening it crashes the devig pivot. So we
# parse the same DK/FD/Pinn data out of prefetched_odds (which holds all
# 16 derivative markets per event) and emit a long frame here.
#
# Each response JSON has top-level fields {id, sport_key, commence_time,
# home_team, away_team, bookmakers[]} where each bookmaker has
# {key, title, last_update, markets[]} and each market has
# {key, last_update, outcomes[]} where each outcome is {name, price, point}.
#
# @param prefetched_odds tibble with columns event_id, json_response (char)
# @param bookmaker_keys character vector of book keys to keep
# @return tibble with columns id, bookmaker_key, market_key, outcomes_name,
#         outcomes_price, outcomes_point, fetch_time
.parse_prefetched_to_long <- function(prefetched_odds, bookmaker_keys) {
  if (is.null(prefetched_odds) || nrow(prefetched_odds) == 0) {
    return(tibble(
      id = character(), bookmaker_key = character(), market_key = character(),
      outcomes_name = character(), outcomes_price = integer(),
      outcomes_point = numeric(), fetch_time = as.POSIXct(character(), tz = "UTC")
    ))
  }

  out <- list()
  for (i in seq_len(nrow(prefetched_odds))) {
    js <- prefetched_odds$json_response[i]
    if (is.na(js) || !nzchar(js)) next
    parsed <- tryCatch(jsonlite::fromJSON(js, simplifyVector = FALSE),
                       error = function(e) NULL)
    if (is.null(parsed) || is.null(parsed$bookmakers)) next
    eid <- parsed$id %||% prefetched_odds$event_id[i]

    for (bm in parsed$bookmakers) {
      bk <- bm$key
      if (is.null(bk) || !(bk %in% bookmaker_keys)) next
      bm_lu <- bm$last_update  # ISO8601 string, or NULL

      for (mk in bm$markets %||% list()) {
        mkey <- mk$key
        if (is.null(mkey)) next
        # Prefer the market-level last_update (more granular), fall back to bookmaker-level.
        ft_str <- mk$last_update %||% bm_lu
        ft_val <- if (!is.null(ft_str)) {
          tryCatch(ymd_hms(ft_str, tz = "UTC", quiet = TRUE),
                   error = function(e) Sys.time())
        } else Sys.time()
        if (is.na(ft_val)) ft_val <- Sys.time()

        for (oc in mk$outcomes %||% list()) {
          out[[length(out) + 1]] <- tibble(
            id             = as.character(eid),
            bookmaker_key  = bk,
            market_key     = mkey,
            outcomes_name  = oc$name %||% NA_character_,
            outcomes_price = if (!is.null(oc$price)) as.integer(oc$price) else NA_integer_,
            outcomes_point = if (!is.null(oc$point)) as.numeric(oc$point) else NA_real_,
            fetch_time     = ft_val
          )
        }
      }
    }
  }

  if (length(out) == 0) {
    return(tibble(
      id = character(), bookmaker_key = character(), market_key = character(),
      outcomes_name = character(), outcomes_price = integer(),
      outcomes_point = numeric(), fetch_time = as.POSIXct(character(), tz = "UTC")
    ))
  }
  bind_rows(out)
}

# Parse the prefetched_odds JSON cache into long format so DK/FD/Pinnacle
# have F-period (F3/F5/F7) and alt-market coverage in book_odds_by_book.
# Without this, the pill grid only shows DK/FD/Pinn odds for FG h2h + FG
# totals (the only markets in game_odds) and is empty for the F-period bets
# that dominate all_bets_combined.
prefetched_long <- .parse_prefetched_to_long(
  prefetched_odds,
  bookmaker_keys = c("draftkings", "fanduel", "pinnacle")
)

# Kalshi is intentionally excluded from book_odds_by_book — it's a
# peer-to-peer venue with its own dashboard tab/flow, not a sportsbook
# to compare against for the odds-screen view.
book_odds_by_book <- list(
  wagerzon  = .scraper_to_canonical(wagerzon_odds),
  hoop88    = .scraper_to_canonical(hoop88_odds),
  bfa       = .scraper_to_canonical(bfa_odds),
  bookmaker = .scraper_to_canonical(bookmaker_odds),
  bet105    = .scraper_to_canonical(bet105_odds),
  draftkings = .odds_api_to_canonical(
                 prefetched_long %>% filter(bookmaker_key == "draftkings")),
  fanduel   = .odds_api_to_canonical(
                 prefetched_long %>% filter(bookmaker_key == "fanduel")),
  pinnacle  = .odds_api_to_canonical(
                 prefetched_long %>% filter(bookmaker_key == "pinnacle"))
)

# Strip NULL entries so expand_bets_to_book_prices doesn't iterate empty lists.
book_odds_by_book <- Filter(Negate(is.null), book_odds_by_book)

book_prices_long <- expand_bets_to_book_prices(all_bets_combined,
                                               book_odds_by_book)

# Data-bug guard: pick book should always be on the exact line.
pick_book_mismatches <- book_prices_long %>%
  inner_join(all_bets_combined %>%
               select(bet_row_id, pick_book = bookmaker_key),
             by = "bet_row_id") %>%
  filter(side == "pick" & bookmaker == pick_book & !is_exact_line)
if (nrow(pick_book_mismatches) > 0) {
  warning(sprintf(
    "[odds_screen] %d pick-book rows on mismatched line — investigate (data bug)",
    nrow(pick_book_mismatches)))
}

cat(sprintf("Built %d book_prices rows for %d bets across %d books.\n",
            nrow(book_prices_long), nrow(all_bets_combined),
            length(unique(book_prices_long$bookmaker))))

# =============================================================================
# PHASE 8: SAVE TO DUCKDB
# =============================================================================

# Disconnect the long-held mlb.duckdb writer FIRST so it cannot block the
# brief mlb_mm.duckdb write below.
dbDisconnect(con_mlb)

# Write mlb_bets_combined to mlb_mm.duckdb (moved from mlb.duckdb on
# 2026-05-05 to free the dashboard's bet loader from the pipeline's long
# write lock on mlb.duckdb — finishes the migration started 2026-04-30).
con_bets <- duckdb_connect_retry("mlb_mm.duckdb")
on.exit(tryCatch(dbDisconnect(con_bets), error = function(e) NULL), add = TRUE)
dbExecute(con_bets, "DROP TABLE IF EXISTS mlb_bets_combined")
dbWriteTable(con_bets, "mlb_bets_combined", all_bets_combined)
dbExecute(con_bets, "DROP TABLE IF EXISTS mlb_bets_book_prices")
dbWriteTable(con_bets, "mlb_bets_book_prices", as.data.frame(book_prices_long))
dbDisconnect(con_bets)
on.exit(NULL)

timer$mark("save_bets")
cat(sprintf("Saved %d bets to mlb_bets_combined table.\n", nrow(all_bets_combined)))
cat("\n=== MLB ANSWER KEY: Complete ===\n")
