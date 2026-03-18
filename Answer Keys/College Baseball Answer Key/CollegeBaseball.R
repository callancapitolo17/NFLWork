# College Baseball Correlated Parlay Answer Key
# Production pipeline: finds +EV ML+Total parlays at offshore books
#
# The edge: Books multiply decimal odds independently (no correlation adjustment).
# Answer key captures real ML+Total correlation → fair price differs.
# Positively-correlated combos (Away ML + Over, Home ML + Under) are underpriced.
#
# Follows CBB.R merged single-script pattern (7 phases).

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

cat("=== COLLEGE BASEBALL CORRELATED PARLAY ANSWER KEY ===\n")

# =============================================================================
# PHASE 1: LOAD HISTORICAL DATA
# =============================================================================

hist_db <- "~/NFLWork/college_baseball_correlation/college_baseball.duckdb"
if (!file.exists(hist_db)) {
  cat("Error: Historical database not found. Run validate.py first.\n")
  quit(status = 1)
}

con <- dbConnect(duckdb(), dbdir = hist_db, read_only = TRUE)

matched <- dbGetQuery(con, "
  WITH consensus AS (
    SELECT
      home_team, away_team,
      CAST(commence_time AS DATE) as game_date,
      MEDIAN(home_ml) as home_ml,
      MEDIAN(away_ml) as away_ml,
      MEDIAN(total_line) as total_line
    FROM odds
    WHERE home_ml IS NOT NULL AND away_ml IS NOT NULL AND total_line IS NOT NULL
    GROUP BY home_team, away_team, CAST(commence_time AS DATE)
  )
  SELECT
    g.game_id,
    g.home_team, g.away_team,
    g.home_score, g.away_score,
    g.total_runs, g.home_margin,
    CASE WHEN g.home_won THEN 1 ELSE 0 END as home_won,
    c.home_ml, c.away_ml, c.total_line
  FROM games g
  INNER JOIN consensus c
    ON g.home_team = c.home_team AND g.away_team = c.away_team
    AND g.date = c.game_date
  WHERE NOT g.went_extras AND NOT g.mercy_rule
    AND c.home_ml IS NOT NULL AND c.away_ml IS NOT NULL
")
dbDisconnect(con)

if (nrow(matched) < 100) {
  cat(sprintf("Insufficient historical data: %d games (need 100+). Exiting.\n", nrow(matched)))
  quit(status = 0)
}

# Devig moneylines to get P(home)
imp_home <- ifelse(matched$home_ml > 0,
                   100 / (matched$home_ml + 100),
                   abs(matched$home_ml) / (abs(matched$home_ml) + 100))
imp_away <- ifelse(matched$away_ml > 0,
                   100 / (matched$away_ml + 100),
                   abs(matched$away_ml) / (abs(matched$away_ml) + 100))
total_imp <- imp_home + imp_away
p_home <- imp_home / total_imp

DT <- data.table(
  id = matched$game_id,
  home_team = matched$home_team,
  away_team = matched$away_team,
  home_spread = qlogis(p_home),  # logit(P(home)) as pseudo-spread
  total_line = matched$total_line,
  home_margin = matched$home_margin,
  total_final_score = matched$total_runs,
  actual_cover = matched$home_won,
  actual_over = as.integer(matched$total_runs > matched$total_line),
  home_ml_odds = matched$home_ml,
  away_ml_odds = matched$away_ml
)

cat(sprintf("Historical data loaded: %d games\n", nrow(DT)))

disp <- compute_dispersion(DT, moneyline = TRUE,
                           odds_cols = c("home_ml_odds", "away_ml_odds"))
ss <- disp$ss
st <- disp$st

N <- max(50, round(nrow(DT) * 0.05))
cat(sprintf("Dispersion: ss=%.2f, st=%.2f, N=%d\n", ss, st, N))
timer$mark("historical_load")

# =============================================================================
# PHASE 2: GET TODAY'S ODDS & BUILD CONSENSUS
# =============================================================================

cat("Fetching Odds API for college baseball...\n")

game_odds <- tryCatch({
  toa_sports_odds(
    sport_key = "baseball_ncaa",
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
  cat("No games found from Odds API. Exiting.\n")
  quit(status = 0)
}

cat(sprintf("API returned %d rows (%d unique games) before filtering.\n",
            nrow(game_odds), n_distinct(game_odds$id)))

if (is.character(game_odds$commence_time)) {
  game_odds$commence_time <- ymd_hms(game_odds$commence_time, tz = "UTC")
}

# Filter to future games only
n_before <- n_distinct(game_odds$id)
game_odds <- game_odds %>% filter(commence_time > Sys.time())
n_after <- n_distinct(game_odds$id)
cat(sprintf("Filtered in-progress: %d -> %d games.\n", n_before, n_after))

if (nrow(game_odds) == 0) {
  cat("No upcoming games. Exiting.\n")
  quit(status = 0)
}

# Build consensus ML
ml_odds <- game_odds %>%
  filter(market_key == "h2h") %>%
  mutate(
    date = as.Date(commence_time),
    odds_type = ifelse(
      home_team == outcomes_name, "home",
      ifelse(away_team == outcomes_name, "away", outcomes_name)
    )
  ) %>%
  pivot_wider(
    id_cols = c(id, commence_time, home_team, away_team, bookmaker_key, date),
    names_from = odds_type,
    values_from = c(outcomes_price, outcomes_point)
  )

# Devig each book's ML
if (!"outcomes_price_home" %in% names(ml_odds) || !"outcomes_price_away" %in% names(ml_odds)) {
  cat("No moneyline odds found. Exiting.\n")
  quit(status = 0)
}

ml_devigged <- ml_odds %>%
  filter(!is.na(outcomes_price_home) & !is.na(outcomes_price_away)) %>%
  mutate(
    imp_home = ifelse(outcomes_price_home > 0,
                      100 / (outcomes_price_home + 100),
                      abs(outcomes_price_home) / (abs(outcomes_price_home) + 100)),
    imp_away = ifelse(outcomes_price_away > 0,
                      100 / (outcomes_price_away + 100),
                      abs(outcomes_price_away) / (abs(outcomes_price_away) + 100)),
    total_imp = imp_home + imp_away,
    p_home = imp_home / total_imp,
    p_away = imp_away / total_imp
  )

# Consensus ML per game (median across books)
consensus_ml <- ml_devigged %>%
  group_by(id, home_team, away_team, date, commence_time) %>%
  summarise(
    consensus_p_home = median(p_home),
    consensus_p_away = median(p_away),
    n_books_ml = n(),
    .groups = "drop"
  )

# Build consensus total
total_odds <- game_odds %>%
  filter(market_key == "totals") %>%
  mutate(
    date = as.Date(commence_time),
    odds_type = ifelse(
      outcomes_name == "Over", "Over",
      ifelse(outcomes_name == "Under", "Under", outcomes_name)
    )
  ) %>%
  pivot_wider(
    id_cols = c(id, commence_time, home_team, away_team, bookmaker_key, date),
    names_from = odds_type,
    values_from = c(outcomes_price, outcomes_point)
  )

if (!"outcomes_point_Over" %in% names(total_odds)) {
  cat("No totals found. Exiting.\n")
  quit(status = 0)
}

consensus_total <- total_odds %>%
  filter(!is.na(outcomes_point_Over)) %>%
  group_by(id) %>%
  summarise(
    total_line = median(outcomes_point_Over, na.rm = TRUE),
    .groups = "drop"
  )

# Merge ML + totals
today_odds <- consensus_ml %>%
  inner_join(consensus_total, by = "id") %>%
  mutate(
    implied_spread = qlogis(consensus_p_home)  # logit transform
  )

cat(sprintf("Today's games with ML + total: %d\n", nrow(today_odds)))
timer$mark("consensus")

if (nrow(today_odds) == 0) {
  cat("No games with both ML and total odds. Exiting.\n")
  quit(status = 0)
}

# =============================================================================
# PHASE 3: GENERATE SAMPLES
# =============================================================================

cat("Generating samples for all games...\n")

targets <- today_odds %>%
  transmute(
    id,
    parent_spread = implied_spread,
    parent_total = total_line,
    target_cover = consensus_p_home,
    target_over = 0.5
  )

samples <- generate_all_samples(
  targets         = targets,
  DT              = DT,
  ss              = ss,
  st              = st,
  N               = N,
  use_spread_line = FALSE
)
cat(sprintf("Generated %d samples.\n", length(samples)))
timer$mark("sample_gen")

# =============================================================================
# PHASE 4: COMPUTE FAIR PARLAY ODDS FOR EACH GAME
# =============================================================================

cat("Computing fair parlay odds (4 combos per game)...\n")

parlay_results <- list()

for (game_id in names(samples)) {
  samp <- samples[[game_id]]$sample
  if (is.null(samp) || nrow(samp) < 30) next

  game_info <- today_odds %>% filter(id == game_id)
  if (nrow(game_info) == 0) next
  tl <- game_info$total_line[1]

  combos <- list(
    home_over  = list(
      list(market = "moneyline", period = "Full", side = "home", line = NA),
      list(market = "total", period = "Full", side = "over", line = tl)
    ),
    home_under = list(
      list(market = "moneyline", period = "Full", side = "home", line = NA),
      list(market = "total", period = "Full", side = "under", line = tl)
    ),
    away_over  = list(
      list(market = "moneyline", period = "Full", side = "away", line = NA),
      list(market = "total", period = "Full", side = "over", line = tl)
    ),
    away_under = list(
      list(market = "moneyline", period = "Full", side = "away", line = NA),
      list(market = "total", period = "Full", side = "under", line = tl)
    )
  )

  for (combo_name in names(combos)) {
    parlay <- tryCatch(
      compute_parlay_fair_odds(samp, combos[[combo_name]]),
      error = function(e) NULL
    )
    if (is.null(parlay)) next

    parlay_results[[length(parlay_results) + 1]] <- data.frame(
      id = game_id,
      home_team = game_info$home_team[1],
      away_team = game_info$away_team[1],
      commence_time = game_info$commence_time[1],
      total_line = tl,
      combo = combo_name,
      fair_joint_prob = parlay$joint_prob,
      fair_american_odds = parlay$fair_american_odds,
      fair_decimal_odds = parlay$fair_decimal_odds,
      correlation_factor = parlay$correlation_factor,
      independent_prob = parlay$independent_prob,
      n_samples = parlay$n_samples_resolved,
      leg1_prob = parlay$leg_probs[1],
      leg2_prob = parlay$leg_probs[2],
      stringsAsFactors = FALSE
    )
  }
}

if (length(parlay_results) == 0) {
  cat("No parlay results generated. Check sample quality.\n")
  quit(status = 0)
}

parlays <- bind_rows(parlay_results)
cat(sprintf("Computed %d parlay fair prices across %d games.\n",
            nrow(parlays), n_distinct(parlays$id)))
timer$mark("parlay_fair_odds")

# =============================================================================
# PHASE 5: WAIT FOR SCRAPERS, LOAD BOOK ODDS
# =============================================================================

sentinel <- file.path(getwd(), ".scrapers_done_college_baseball")
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
wagerzon_odds <- tryCatch(get_wagerzon_odds("college_baseball"), error = function(e) data.frame())
cat(sprintf("Loaded %d Wagerzon records.\n", nrow(wagerzon_odds)))

hoop88_odds <- tryCatch(get_hoop88_odds("college_baseball"), error = function(e) data.frame())
cat(sprintf("Loaded %d Hoop88 records.\n", nrow(hoop88_odds)))
timer$mark("load_scrapers")

# =============================================================================
# PHASE 6: COMPARE PARLAYS — THE KEY STEP
# =============================================================================

# Helper: Convert American odds to decimal
american_to_decimal <- function(odds) {
  ifelse(odds > 0, odds / 100 + 1, 100 / abs(odds) + 1)
}

# Helper: Convert decimal to American
decimal_to_american <- function(dec) {
  ifelse(dec >= 2,
         round((dec - 1) * 100),
         round(-100 / (dec - 1)))
}

# Helper: Apply Hoop88 shave (round toward house)
apply_shave <- function(odds) {
  ifelse(odds > 0,
         floor(odds / 5) * 5,      # +264 -> +260
         -ceiling(abs(odds) / 5) * 5)  # -245 -> -250 (note: reversed for R)
}

# Sizing
bankroll <- 100
kelly_mult <- 0.25

compare_parlays <- function(parlays, offshore_odds, book_name) {
  if (nrow(offshore_odds) == 0) return(tibble())

  # Get ML and total odds per game
  ml_odds <- offshore_odds %>%
    filter(market_type == "h2h", period == "Full") %>%
    select(home_team, away_team, odds_home, odds_away)

  total_offshore <- offshore_odds %>%
    filter(market_type == "totals", period == "Full") %>%
    select(home_team, away_team, line, odds_over, odds_under)

  if (nrow(ml_odds) == 0 || nrow(total_offshore) == 0) return(tibble())

  edges <- list()

  for (i in seq_len(nrow(parlays))) {
    p <- parlays[i, ]

    # Find matching game in offshore odds
    ml_match <- ml_odds %>%
      filter(home_team == p$home_team & away_team == p$away_team)
    tot_match <- total_offshore %>%
      filter(home_team == p$home_team & away_team == p$away_team &
               abs(line - p$total_line) < 0.5)

    if (nrow(ml_match) == 0 || nrow(tot_match) == 0) next

    # Get the right leg odds based on combo
    ml_side <- ifelse(grepl("^home", p$combo), "home", "away")
    total_side <- ifelse(grepl("over$", p$combo), "over", "under")

    book_ml_american <- if (ml_side == "home") ml_match$odds_home[1] else ml_match$odds_away[1]
    book_total_american <- if (total_side == "over") tot_match$odds_over[1] else tot_match$odds_under[1]

    if (is.na(book_ml_american) || is.na(book_total_american)) next

    # Book prices parlay independently (no correlation adjustment)
    book_ml_dec <- american_to_decimal(book_ml_american)
    book_total_dec <- american_to_decimal(book_total_american)
    book_parlay_dec <- book_ml_dec * book_total_dec
    book_parlay_american <- decimal_to_american(book_parlay_dec)

    # Apply shave for Hoop88
    if (book_name == "hoop88") {
      book_parlay_american <- apply_shave(book_parlay_american)
      book_parlay_dec <- american_to_decimal(book_parlay_american)
    }

    # Fair price from answer key (correlation-adjusted)
    fair_dec <- p$fair_decimal_odds

    # Edge
    edge_pct <- (book_parlay_dec / fair_dec - 1) * 100

    # Book implied prob
    book_implied_prob <- 1 / book_parlay_dec

    # Kelly sizing
    if (edge_pct > 0) {
      ev <- p$fair_joint_prob * (book_parlay_dec - 1) - (1 - p$fair_joint_prob)
      edge_fraction <- ev / (book_parlay_dec - 1)
      kelly_bet <- round(edge_fraction * kelly_mult * bankroll, 2)
    } else {
      ev <- 0
      kelly_bet <- 0
    }

    edges[[length(edges) + 1]] <- data.frame(
      id = p$id,
      home_team = p$home_team,
      away_team = p$away_team,
      commence_time = p$commence_time,
      combo = p$combo,
      total_line = p$total_line,
      book = book_name,
      book_ml = book_ml_american,
      book_total = book_total_american,
      book_parlay_american = book_parlay_american,
      book_parlay_dec = round(book_parlay_dec, 3),
      fair_american = p$fair_american_odds,
      fair_dec = round(fair_dec, 3),
      fair_prob = round(p$fair_joint_prob, 4),
      cf = p$correlation_factor,
      edge_pct = round(edge_pct, 2),
      ev = round(ev, 4),
      kelly_bet = kelly_bet,
      stringsAsFactors = FALSE
    )
  }

  if (length(edges) == 0) return(tibble())
  bind_rows(edges)
}

all_edges <- bind_rows(
  compare_parlays(parlays, wagerzon_odds, "wagerzon"),
  compare_parlays(parlays, hoop88_odds, "hoop88")
)
timer$mark("compare_parlays")

# =============================================================================
# PHASE 7: OUTPUT
# =============================================================================

if (nrow(all_edges) == 0) {
  cat("\nNo parlay comparisons possible (no offshore odds matched).\n")
  cat("=== COLLEGE BASEBALL ANSWER KEY: Complete (no edges) ===\n")
  quit(status = 0)
}

# Filter to +EV only
ev_bets <- all_edges %>%
  filter(edge_pct > 0) %>%
  arrange(desc(edge_pct))

cat(sprintf("\n=== PARLAY EDGE SUMMARY ===\n"))
cat(sprintf("Total comparisons: %d\n", nrow(all_edges)))
cat(sprintf("+EV parlays found: %d\n", nrow(ev_bets)))

if (nrow(ev_bets) > 0) {
  cat(sprintf("Average edge: %.1f%%\n", mean(ev_bets$edge_pct)))
  cat(sprintf("Total Kelly stake: $%.2f\n", sum(ev_bets$kelly_bet)))

  cat("\n=== +EV PARLAYS (sorted by edge) ===\n")
  ev_bets %>%
    transmute(
      game = paste0(away_team, " @ ", home_team),
      combo,
      total_line,
      book,
      book_odds = book_parlay_american,
      fair_odds = fair_american,
      CF = cf,
      edge = paste0("+", edge_pct, "%"),
      kelly = paste0("$", kelly_bet)
    ) %>%
    print(n = 50)

  # Summary by combo type
  cat("\n=== EDGE BY COMBO TYPE ===\n")
  ev_bets %>%
    group_by(combo) %>%
    summarise(
      n = n(),
      avg_edge = mean(edge_pct),
      avg_cf = mean(cf),
      total_kelly = sum(kelly_bet),
      .groups = "drop"
    ) %>%
    arrange(desc(avg_edge)) %>%
    print()

  # Summary by book
  cat("\n=== EDGE BY BOOK ===\n")
  ev_bets %>%
    group_by(book) %>%
    summarise(
      n = n(),
      avg_edge = mean(edge_pct),
      total_kelly = sum(kelly_bet),
      .groups = "drop"
    ) %>%
    print()
}

cat("\n=== COLLEGE BASEBALL ANSWER KEY: Complete ===\n")
timer$report()
