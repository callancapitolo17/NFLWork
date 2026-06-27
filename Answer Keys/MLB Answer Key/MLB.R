# MLB Answer Key - Merged pipeline
# Single R process: prepare samples, build F5 predictions, compare to offshore
# Runs in parallel with scrapers via run.py; waits for sentinel before loading scraper data

.t_script_start <- Sys.time()
# Resolve cwd from the script location so a worktree copy stays self-contained.
# MLB.R lives at <root>/Answer Keys/MLB Answer Key/MLB.R, so its Answer Keys dir
# is the parent of the script dir. On main this is identical to ~/NFLWork/Answer Keys.
.args <- commandArgs(trailingOnly = FALSE)
.file_arg <- grep("^--file=", .args, value = TRUE)
.script_dir <- if (length(.file_arg) > 0) {
  .raw <- gsub("~\\+~", " ", sub("^--file=", "", .file_arg[1]))
  normalizePath(dirname(.raw), mustWork = FALSE)
} else {
  normalizePath("~/NFLWork/Answer Keys/MLB Answer Key", mustWork = FALSE)
}
setwd(dirname(.script_dir))  # <root>/Answer Keys
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
set_nflwork_root(derive_repo_root())  # so get_*_odds read this repo's DBs
source("triple_play_helpers.R")
source("MLB Answer Key/odds_screen.R")
source("MLB Answer Key/market_edge.R")
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

# Per-book staleness cutoff — drop rows older than this before canonicalization.
# Prevents stale book quotes from rendering as if they were fresh on the dashboard.
BOOK_STALENESS_CUTOFF_MIN <- 30

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

# DK + FD single-leg odds are scraped to their own DuckDB files; used only
# for the bets-tab odds-screen pills (book_odds_by_book below). Pinnacle
# still flows through the Odds API prefetched cache. Helpers in Tools.R
# mirror get_bookmaker_odds(): wide 18-col schema -> canonical long format.
dk_odds <- tryCatch(get_dk_odds("mlb"), error = function(e) { cat(sprintf("DraftKings skip: %s\n", e$message)); tibble() })
cat(sprintf("Loaded %d DraftKings records.\n", nrow(dk_odds)))

fd_odds <- tryCatch(get_fd_odds("mlb"), error = function(e) { cat(sprintf("FanDuel skip: %s\n", e$message)); tibble() })
cat(sprintf("Loaded %d FanDuel records.\n", nrow(fd_odds)))

kalshi_odds <- tryCatch(get_kalshi_odds("mlb"), error = function(e) { cat(sprintf("Kalshi skip: %s\n", e$message)); tibble() })
cat(sprintf("Loaded %d Kalshi records.\n", nrow(kalshi_odds)))

# Helper — drop scraper rows older than BOOK_STALENESS_CUTOFF_MIN minutes by
# fetch_time. Defensive: a stuck scraper would otherwise render hours-old
# quotes on the bets tab as if they were live. NA fetch_times are kept (they
# get replaced with Sys.time() downstream in scraper_to_canonical).
.drop_stale_book_rows <- function(df, book_label) {
  if (is.null(df) || nrow(df) == 0) return(df)
  if (!"fetch_time" %in% names(df)) return(df)
  cutoff <- Sys.time() - lubridate::minutes(BOOK_STALENESS_CUTOFF_MIN)
  fresh <- df %>% filter(is.na(fetch_time) | fetch_time >= cutoff)
  dropped_n <- nrow(df) - nrow(fresh)
  if (dropped_n > 0) {
    message(sprintf("[%s] dropped %d row(s) older than %d min (fetch_time gate)",
                    book_label, dropped_n, BOOK_STALENESS_CUTOFF_MIN))
  }
  fresh
}

wagerzon_odds  <- .drop_stale_book_rows(wagerzon_odds,  "wagerzon")
hoop88_odds    <- .drop_stale_book_rows(hoop88_odds,    "hoop88")
bfa_odds       <- .drop_stale_book_rows(bfa_odds,       "bfa")
bookmaker_odds <- .drop_stale_book_rows(bookmaker_odds, "bookmaker")
bet105_odds    <- .drop_stale_book_rows(bet105_odds,    "bet105")
dk_odds        <- .drop_stale_book_rows(dk_odds,        "dk")
fd_odds        <- .drop_stale_book_rows(fd_odds,        "fd")
kalshi_odds    <- .drop_stale_book_rows(kalshi_odds,    "kalshi")

# Map all scraper markets to Odds API F5 convention
wagerzon_odds  <- map_scraper_markets_mlb(wagerzon_odds)
hoop88_odds    <- map_scraper_markets_mlb(hoop88_odds)
bfa_odds       <- map_scraper_markets_mlb(bfa_odds)
bookmaker_odds <- map_scraper_markets_mlb(bookmaker_odds)
bet105_odds    <- map_scraper_markets_mlb(bet105_odds)
dk_odds        <- map_scraper_markets_mlb(dk_odds)
fd_odds        <- map_scraper_markets_mlb(fd_odds)
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
  # bet_row_id: canonical identity hash shared with the market-edge path so
  # find_market_edges() output joins to model bets for the same wager. Recipe:
  # (id, base_market_type, period, normalized_line, bet_on) — see
  # odds_screen.R::compute_bet_row_id. Book-agnostic (all books for the same
  # line/side share an id) AND market-string-agnostic (the canonical identity
  # is reproducible from canonical per-book odds, unlike the verbose market).
  mutate(bet_row_id = compute_bet_row_id(
    id, .derive_market_type(market), .derive_period(market), line, bet_on
  )) %>%
  mutate(edge_source = "model", model_ev = ev, market_ev = NA_real_) %>%
  mutate(base_market = gsub("^alternate_", "", market)) %>%
  group_by(id, base_market, bet_on) %>%
  filter(ev == max(ev)) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(-base_market) %>%
  arrange(desc(ev))

# Game id lookup: join scraper wide frames to Odds API game IDs via home/away team.
# Passed as `lookup` argument to scraper_to_canonical() (defined in odds_screen.R).
.game_id_lookup <- mlb_odds %>%
  select(id, home_team, away_team) %>%
  distinct()

# Parse the prefetched_odds JSON cache into long format for Pinnacle pill
# coverage. DK and FD come from their per-book DuckDBs via get_dk_odds /
# get_fd_odds (see Tools.R); the Odds API path is used only for Pinnacle,
# which has no public REST API.
prefetched_long <- parse_prefetched_to_long(
  prefetched_odds,
  bookmaker_keys = "pinnacle"
)

# Kalshi is intentionally excluded from book_odds_by_book — it's a
# peer-to-peer venue with its own dashboard tab/flow, not a sportsbook
# to compare against for the odds-screen view.
book_odds_by_book <- list(
  wagerzon  = scraper_to_canonical(wagerzon_odds,  .game_id_lookup, book_name = "wagerzon"),
  hoop88    = scraper_to_canonical(hoop88_odds,    .game_id_lookup, book_name = "hoop88"),
  bfa       = scraper_to_canonical(bfa_odds,       .game_id_lookup, book_name = "bfa"),
  bookmaker = scraper_to_canonical(bookmaker_odds, .game_id_lookup, book_name = "bookmaker"),
  bet105    = scraper_to_canonical(bet105_odds,    .game_id_lookup, book_name = "bet105"),
  draftkings = scraper_to_canonical(dk_odds, .game_id_lookup, book_name = "dk"),
  fanduel    = scraper_to_canonical(fd_odds, .game_id_lookup, book_name = "fd"),
  pinnacle  = odds_api_to_canonical(
                 prefetched_long %>% filter(bookmaker_key == "pinnacle"))
)

# Strip NULL entries so expand_bets_to_book_prices doesn't iterate empty lists.
book_odds_by_book <- Filter(Negate(is.null), book_odds_by_book)

# ---- Market-consensus edges (model OR market flagging) ----
# Second, independent edge signal: flag a bet when a book's price is out of line
# with the leave-one-out devigged consensus of the OTHER books (>= EV_THRESHOLD),
# even when the model is neutral. See market_edge.R::find_market_edges.
# game_info gives market-only cards their teams + tipoff and the future-game gate.
# ungroup() is REQUIRED: mlb_odds arrives grouped (by id) from the consensus
# pipeline, and dplyr transmute() silently RETAINS grouping columns. Without
# ungroup, game_info keeps an `id` column alongside game_id; that `id` then
# collides with find_market_edges()'s own `id` in its left_join(by="game_id"),
# renaming it `id.x`. The market rows lose their `id`, MLB.R's merge NA-fills it,
# and expand_bets_to_book_prices() (which joins books on the game id) matches
# nothing — so every MARKET-edge card renders all em-dashes.
.market_game_info <- mlb_odds %>%
  ungroup() %>%
  transmute(game_id = id, home_team, away_team, pt_start_time = commence_time) %>%
  distinct(game_id, .keep_all = TRUE)

market_bets <- find_market_edges(
  book_odds_by_book,
  game_info     = .market_game_info,
  threshold     = EV_THRESHOLD,
  min_others    = 1,
  staleness_min = BOOK_STALENESS_CUTOFF_MIN,
  now           = Sys.time(),
  bankroll      = bankroll,
  kelly_mult    = kelly_mult
) %>%
  filter(is.na(pt_start_time) | pt_start_time > Sys.time())

cat(sprintf("find_market_edges: %d market-edge candidate(s).\n", nrow(market_bets)))

# Merge model + market on the shared canonical bet_row_id.
.model_cols <- names(all_bets_combined)
market_new  <- market_bets %>% filter(!(bet_row_id %in% all_bets_combined$bet_row_id))
market_both <- market_bets %>% filter(bet_row_id %in% all_bets_combined$bet_row_id) %>%
  select(bet_row_id, market_ev_join = market_ev)

# (a) Promote model rows that ALSO have a market edge to "both".
all_bets_combined <- all_bets_combined %>%
  left_join(market_both, by = "bet_row_id") %>%
  mutate(
    # as.character() guards the empty-frame case: base ifelse() seeds its
    # result from the logical condition, so on 0 rows (no model bets today)
    # it returns logical(0), which then can't bind_rows() with the character
    # market-only rows below. Coercing keeps edge_source character-stable.
    edge_source = as.character(ifelse(!is.na(market_ev_join), "both", edge_source)),
    market_ev   = as.numeric(ifelse(!is.na(market_ev_join), market_ev_join, market_ev))
  ) %>%
  select(-market_ev_join)

# (b) Append market-only rows, aligning columns to the model schema.
if (nrow(market_new) > 0) {
  miss <- setdiff(.model_cols, names(market_new))
  for (cc in miss) market_new[[cc]] <- NA
  market_new <- market_new %>% select(all_of(.model_cols))
  all_bets_combined <- bind_rows(all_bets_combined, market_new)
}

# Rank by the better of the two signals; drop any exact-id duplicate.
all_bets_combined <- all_bets_combined %>%
  mutate(ev = pmax(model_ev, market_ev, na.rm = TRUE)) %>%
  distinct(bet_row_id, .keep_all = TRUE) %>%
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

# NOTE: book_odds_by_book (and .game_id_lookup / prefetched_long) are now built
# earlier, at the end of PHASE 7, so the market-consensus merge can run before
# correlation-adjusted Kelly. The expand call below reuses that same object.

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
# Pre-declare schema so fetch_time keeps its UTC timezone (otherwise
# dbWriteTable downgrades POSIXct to naive TIMESTAMP and breaks
# downstream NOW()-fetch_time math).
dbExecute(con_bets, "DROP TABLE IF EXISTS mlb_bets_book_prices")
dbExecute(con_bets, "
  CREATE TABLE mlb_bets_book_prices (
    bet_row_id      VARCHAR,
    game_id         VARCHAR,
    market          VARCHAR,
    period          VARCHAR,
    side            VARCHAR,
    bookmaker       VARCHAR,
    line            DOUBLE,
    line_quoted     DOUBLE,
    is_exact_line   BOOLEAN,
    american_odds   INTEGER,
    derived_fair_odds DOUBLE,
    fetch_time      TIMESTAMPTZ,
    game_start_time TIMESTAMPTZ
  )
")
# Ensure POSIXct columns are in UTC so DuckDB sees TZ-aware values.
book_prices_long$fetch_time      <- as.POSIXct(book_prices_long$fetch_time,      tz = "UTC")
book_prices_long$game_start_time <- as.POSIXct(book_prices_long$game_start_time, tz = "UTC")
dbAppendTable(con_bets, "mlb_bets_book_prices", book_prices_long)
dbDisconnect(con_bets)
on.exit(NULL)

timer$mark("save_bets")
cat(sprintf("Saved %d bets to mlb_bets_combined table.\n", nrow(all_bets_combined)))
cat("\n=== MLB ANSWER KEY: Complete ===\n")
