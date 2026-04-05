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
    g.date as game_date,
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
  game_date = as.Date(matched$game_date),
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
    values_from = outcomes_price
  )

# Devig each book's ML
if (!"home" %in% names(ml_odds) || !"away" %in% names(ml_odds)) {
  cat("No moneyline odds found. Exiting.\n")
  quit(status = 0)
}

ml_devigged <- ml_odds %>%
  filter(!is.na(home) & !is.na(away)) %>%
  mutate(
    imp_home = ifelse(home > 0,
                      100 / (home + 100),
                      abs(home) / (abs(home) + 100)),
    imp_away = ifelse(away > 0,
                      100 / (away + 100),
                      abs(away) / (abs(away) + 100)),
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
    consensus_home_ml = median(home),  # American odds for mean_match
    n_books_ml = n(),
    .groups = "drop"
  )

# Build consensus total
totals_raw <- game_odds %>% filter(market_key == "totals")

if (nrow(totals_raw) == 0 || !"outcomes_point" %in% names(game_odds)) {
  cat("No totals found from Odds API. Exiting.\n")
  quit(status = 0)
}

total_odds <- totals_raw %>%
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
  cat("No totals found after pivot. Exiting.\n")
  quit(status = 0)
}

consensus_total <- total_odds %>%
  filter(!is.na(outcomes_point_Over)) %>%
  group_by(id) %>%
  summarise(
    # Use mode (most common line); break ties with median
    total_line = {
      tbl <- table(outcomes_point_Over)
      as.numeric(names(tbl)[which.max(tbl)])
    },
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

# parent_spread must be on the same scale as DT$home_ml_odds (American odds)
# because use_spread_line = FALSE makes mean_match match on home_ml_odds column
targets <- today_odds %>%
  transmute(
    id,
    parent_spread = consensus_home_ml,  # American odds, matches home_ml_odds in DT
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
  use_spread_line = FALSE,
  min_N           = 25  # small pool (~486 games), default 50 too aggressive
)
cat(sprintf("Generated %d samples.\n", length(samples)))
timer$mark("sample_gen")

# =============================================================================
# PHASE 4: COMPUTE FAIR PARLAY ODDS FOR EACH GAME
# =============================================================================

cat("Computing fair parlay odds (4 combos per game)...\n")

# --- Method 1: Answer Key (resampled from historical outcomes) ---
ak_results <- list()

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

    ak_results[[length(ak_results) + 1]] <- data.frame(
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
      method = "answer_key",
      stringsAsFactors = FALSE
    )
  }
}

cat(sprintf("Answer Key: %d parlay prices computed.\n", length(ak_results)))

# --- Method 2: Empirical Correlation Factors ---
# CFs from 486 matched games (all correlation sources captured).
# Medium totals (12-13.5) are the sweet spot: Away+Over CF=1.28, Home+Under CF=1.15.
EMPIRICAL_CFS <- list(
  low    = list(away_over = 0.935, away_under = 1.068, home_over = 1.053, home_under = 0.945),
  medium = list(away_over = 1.280, away_under = 0.799, home_over = 0.791, home_under = 1.150),
  high   = list(away_over = 1.097, away_under = 0.956, home_over = 0.948, home_under = 1.023)
)

get_total_bucket <- function(tl) {
  ifelse(tl < 12, "low", ifelse(tl <= 13.5, "medium", "high"))
}

# Also devig over/under for total odds
total_devigged <- total_odds %>%
  filter(!is.na(outcomes_price_Over) & !is.na(outcomes_price_Under)) %>%
  mutate(
    imp_over = ifelse(outcomes_price_Over > 0,
                      100 / (outcomes_price_Over + 100),
                      abs(outcomes_price_Over) / (abs(outcomes_price_Over) + 100)),
    imp_under = ifelse(outcomes_price_Under > 0,
                       100 / (outcomes_price_Under + 100),
                       abs(outcomes_price_Under) / (abs(outcomes_price_Under) + 100)),
    p_over = imp_over / (imp_over + imp_under)
  ) %>%
  group_by(id) %>%
  summarise(consensus_p_over = median(p_over, na.rm = TRUE), .groups = "drop")

today_odds <- today_odds %>% left_join(total_devigged, by = "id")

cf_results <- list()
for (i in seq_len(nrow(today_odds))) {
  g <- today_odds[i, ]
  tl <- g$total_line
  bucket <- get_total_bucket(tl)
  cfs <- EMPIRICAL_CFS[[bucket]]
  p_home <- g$consensus_p_home
  p_away <- g$consensus_p_away
  p_over <- ifelse(is.na(g$consensus_p_over), 0.5, g$consensus_p_over)
  p_under <- 1 - p_over

  combo_map <- list(
    home_over  = list(p1 = p_home, p2 = p_over,  cf = cfs$home_over),
    home_under = list(p1 = p_home, p2 = p_under, cf = cfs$home_under),
    away_over  = list(p1 = p_away, p2 = p_over,  cf = cfs$away_over),
    away_under = list(p1 = p_away, p2 = p_under, cf = cfs$away_under)
  )

  for (combo_name in names(combo_map)) {
    cm <- combo_map[[combo_name]]
    joint <- cm$p1 * cm$p2 * cm$cf
    joint <- max(0.001, min(0.999, joint))
    fair_dec <- 1 / joint
    fair_am <- ifelse(joint >= 0.5,
                      round(-(joint / (1 - joint)) * 100),
                      round(((1 - joint) / joint) * 100))

    cf_results[[length(cf_results) + 1]] <- data.frame(
      id = g$id,
      home_team = g$home_team,
      away_team = g$away_team,
      commence_time = g$commence_time,
      total_line = tl,
      combo = combo_name,
      fair_joint_prob = joint,
      fair_american_odds = fair_am,
      fair_decimal_odds = fair_dec,
      correlation_factor = cm$cf,
      independent_prob = cm$p1 * cm$p2,
      n_samples = NA_integer_,
      leg1_prob = cm$p1,
      leg2_prob = cm$p2,
      method = "empirical_cf",
      stringsAsFactors = FALSE
    )
  }
}

cat(sprintf("Empirical CF: %d parlay prices computed.\n", length(cf_results)))

# --- Combine both methods ---
parlays <- bind_rows(
  if (length(ak_results) > 0) bind_rows(ak_results) else tibble(),
  if (length(cf_results) > 0) bind_rows(cf_results) else tibble()
)

if (nrow(parlays) == 0) {
  cat("No parlay results generated. Check sample quality.\n")
  quit(status = 0)
}

cat(sprintf("Total: %d parlay fair prices across %d games (2 methods).\n",
            nrow(parlays), n_distinct(parlays$id)))
timer$mark("parlay_fair_odds")

# =============================================================================
# PHASE 5: WAIT FOR SCRAPERS, LOAD BOOK ODDS
# =============================================================================

sentinel <- file.path(getwd(), ".scrapers_done_college_baseball")
waited <- 0
max_wait <- ifelse(file.exists(sentinel), 0,
                   ifelse(any(commandArgs(TRUE) == "--no-wait"), 0, 120))
while (!file.exists(sentinel) && waited < max_wait) {
  Sys.sleep(0.5)
  waited <- waited + 0.5
}
if (file.exists(sentinel)) {
  cat(sprintf("Scrapers done (waited %.1fs for sentinel).\n", waited))
} else if (max_wait > 0) {
  cat("Warning: Scraper sentinel not found after 120s, proceeding anyway.\n")
} else {
  cat("Standalone mode: skipping sentinel wait.\n")
}

# Load scraped odds
wagerzon_odds <- tryCatch(get_wagerzon_odds("college_baseball"), error = function(e) data.frame())
cat(sprintf("Loaded %d Wagerzon records.\n", nrow(wagerzon_odds)))

hoop88_odds <- tryCatch(get_hoop88_odds("college_baseball"), error = function(e) data.frame())
cat(sprintf("Loaded %d Hoop88 records.\n", nrow(hoop88_odds)))

# --- Fuzzy team name resolution for college baseball ---
# Offshore books use short names like "OHIO STATE (BB)" while Odds API uses
# "Ohio State Buckeyes". Strip suffix, normalize, and fuzzy-match.
resolve_college_baseball_teams <- function(offshore_df, api_teams) {
  if (nrow(offshore_df) == 0 || length(api_teams) == 0) return(offshore_df)

  # Normalize: strip (BB) suffix, extra spaces, lowercase, trim
  normalize <- function(x) {
    x <- gsub("\\s*\\(BB\\)\\s*$", "", x, ignore.case = TRUE)
    x <- gsub("\\s+", " ", trimws(tolower(x)))
    x
  }

  # Known aliases (offshore normalized -> api normalized query)
  # Critical: "State" schools MUST be listed explicitly to avoid substring
  # matching "georgia state" -> "georgia" or "kansas state" -> "kansas"
  aliases <- c(
    "miami florida"    = "miami hurricanes",
    "usc upstate"      = "south carolina upstate",
    "connecticut"      = "uconn",
    "cs bakersfield"   = "csu bakersfield",
    "georgia state"    = "georgia st panthers",
    "kansas state"     = "kansas st wildcats",
    "arizona state"    = "arizona st sun devils",
    "central florida"  = "ucf knights",
    "north carolina st"= "nc state wolfpack",
    "n carolina st"    = "nc state wolfpack",
    "nc state"         = "nc state wolfpack",
    "arkansas state"   = "arkansas st red wolves",
    "mississippi state"= "mississippi st bulldogs",
    "ohio state"       = "ohio st buckeyes",
    "oklahoma state"   = "oklahoma st cowboys",
    "oregon state"     = "oregon st beavers",
    "penn state"       = "penn state nittany lions",
    "florida st"       = "florida st seminoles",
    "florida state"    = "florida st seminoles",
    "iowa state"       = "iowa st cyclones",
    "michigan st"      = "michigan st spartans",
    "michigan state"   = "michigan st spartans",
    "texas state"      = "texas st bobcats",
    "ball state"       = "ball st cardinals",
    "boise state"      = "boise st broncos",
    "fresno state"     = "fresno st bulldogs",
    "sam houston state"= "sam houston bearkats",
    "sam houston"      = "sam houston bearkats",
    "wright state"     = "wright st raiders",
    "wichita state"    = "wichita st shockers",
    "kent state"       = "kent st golden flashes",
    "san diego state"  = "san diego st aztecs",
    "san jose state"   = "san jose st spartans",
    "washington state" = "washington st cougars",
    "murray state"     = "murray st racers",
    "louisiana"        = "louisiana ragin' cajuns",
    "louisiana monroe" = "louisiana-monroe warhawks",
    "louisiana lafayette" = "louisiana ragin' cajuns",
    "south carolina"   = "south carolina gamecocks",
    "north carolina"   = "north carolina tar heels",
    "unc wilmington"   = "unc wilmington seahawks",
    "unc greensboro"   = "unc greensboro spartans",
    "ole miss"         = "ole miss rebels",
    "mississippi"      = "ole miss rebels",
    "pitt"             = "pittsburgh panthers",
    "umass"            = "massachusetts minutemen",
    "texas tech"       = "texas tech red raiders",
    "virginia tech"    = "virginia tech hokies",
    "duke"             = "duke blue devils",
    "georgia tech"     = "georgia tech yellow jackets",
    "texas a&m"        = "texas a&m aggies",
    "texas am"         = "texas a&m aggies",
    "texas san antonio"= "utsa roadrunners",
    "utsa"             = "utsa roadrunners",
    "florida atlantic" = "florida atlantic owls",
    "southern miss"    = "southern miss golden eagles"
  )

  # Build index: for each API team, store no-mascot and full lowercase
  api_no_mascot <- setNames(api_teams, sapply(api_teams, function(n) {
    trimws(tolower(sub(" [^ ]+$", "", n)))
  }))
  api_full <- setNames(api_teams, tolower(api_teams))

  resolve_one <- function(raw) {
    norm <- normalize(raw)

    # Check aliases first — these are authoritative and override all fuzzy matching.
    # If alias fires but target not in today's API games, return alias target as-is
    # to prevent false substring matches (e.g., "Georgia State" → "Georgia Bulldogs").
    if (norm %in% names(aliases)) {
      alias_target <- unname(aliases[norm])
      # Find API team matching the alias target (full name)
      hit <- api_full[alias_target]
      if (!is.na(hit)) return(unname(hit))
      # Try no-mascot version of alias target (strip last word = mascot)
      alias_short <- trimws(tolower(sub(" [^ ]+$", "", alias_target)))
      hit <- api_no_mascot[alias_short]
      if (!is.na(hit)) return(unname(hit))
      # Alias is authoritative — return target even if not in today's API games.
      # This prevents false matches; the game simply won't have a comparison.
      return(alias_target)
    }

    # Exact match on no-mascot
    hit <- api_no_mascot[norm]
    if (!is.na(hit)) return(unname(hit))

    # Exact match on full name
    hit <- api_full[norm]
    if (!is.na(hit)) return(unname(hit))

    # Try no-mascot version of input against api_no_mascot keys
    norm_no_mascot <- trimws(sub(" [^ ]+$", "", norm))
    hit <- api_no_mascot[norm_no_mascot]
    if (!is.na(hit)) return(unname(hit))

    # agrep fuzzy match on no-mascot (last resort, strict distance)
    ag <- agrep(norm_no_mascot, names(api_no_mascot), max.distance = 0.15, value = TRUE)
    if (length(ag) == 1) return(unname(api_no_mascot[ag[1]]))  # only if unique match
    return(raw)  # no match
  }

  # Get unique raw names and resolve
  all_raw <- unique(c(offshore_df$home_team, offshore_df$away_team))
  mapping <- setNames(sapply(all_raw, resolve_one), all_raw)

  n_resolved <- sum(mapping != names(mapping))
  cat(sprintf("  Team resolution: %d/%d names resolved\n", n_resolved, length(all_raw)))
  for (raw in names(mapping)) {
    if (mapping[raw] != raw) cat(sprintf("    %s -> %s\n", raw, mapping[raw]))
  }

  offshore_df$home_team <- unname(mapping[offshore_df$home_team])
  offshore_df$away_team <- unname(mapping[offshore_df$away_team])
  offshore_df
}

# Apply to both books using today's Odds API game names
api_team_names <- unique(c(today_odds$home_team, today_odds$away_team))
if (nrow(wagerzon_odds) > 0) {
  wagerzon_odds <- resolve_college_baseball_teams(wagerzon_odds, api_team_names)
}
if (nrow(hoop88_odds) > 0) {
  hoop88_odds <- resolve_college_baseball_teams(hoop88_odds, api_team_names)
}

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
bankroll <- 2000
kelly_mult <- 0.25

# Wagerzon rounds profit to the nearest whole dollar ($63.51 -> $64, $6.49 -> $6).
# Find the bet size near Kelly where raw_win lands just above X.50 to maximize bonus.
optimize_wagerzon_stake <- function(kelly_stake, decimal_odds, search_range = 5) {
  if (kelly_stake <= 0) return(list(stake = 0, win = 0, eff_decimal = decimal_odds, bonus = 0))
  lo <- max(1, floor(kelly_stake) - search_range)
  hi <- ceiling(kelly_stake) + search_range
  candidates <- seq(lo, hi, by = 1)

  best <- list(stake = round(kelly_stake), bonus = -Inf)
  for (s in candidates) {
    raw_win <- s * (decimal_odds - 1)
    rounded_win <- round(raw_win)
    if (rounded_win <= 0) next
    bonus <- rounded_win - raw_win
    eff_dec <- (rounded_win / s) + 1
    if (bonus > best$bonus) {
      best <- list(stake = s, win = rounded_win, eff_decimal = eff_dec, bonus = bonus)
    }
  }
  best
}

compare_parlays <- function(parlays, offshore_odds, book_name) {
  if (nrow(offshore_odds) == 0) return(tibble())

  # Get ML and total odds per game
  # Period can be "Full", "fg", or "Game" depending on book
  full_periods <- c("Full", "fg", "Game")
  ml_odds <- offshore_odds %>%
    filter(market_type == "h2h", period %in% full_periods) %>%
    select(home_team, away_team, odds_home, odds_away)

  total_offshore <- offshore_odds %>%
    filter(market_type == "totals", period %in% full_periods) %>%
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

    # Wagerzon rounding optimization: find stake near Kelly that maximizes
    # the rounding bonus (profit rounded to nearest $5)
    opt_stake <- kelly_bet
    opt_win <- round(kelly_bet * (book_parlay_dec - 1))
    rounding_bonus <- 0
    if (book_name == "wagerzon" && kelly_bet > 0) {
      opt <- optimize_wagerzon_stake(kelly_bet, book_parlay_dec)
      opt_stake <- opt$stake
      opt_win <- opt$win
      rounding_bonus <- round(opt$bonus, 2)
    } else if (kelly_bet > 0) {
      opt_win <- round(kelly_bet * (book_parlay_dec - 1))
    }

    edges[[length(edges) + 1]] <- data.frame(
      id = p$id,
      home_team = p$home_team,
      away_team = p$away_team,
      commence_time = p$commence_time,
      combo = p$combo,
      total_line = p$total_line,
      method = if ("method" %in% names(p)) p$method else "answer_key",
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
      opt_stake = opt_stake,
      opt_win = opt_win,
      rounding_bonus = rounding_bonus,
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

# For each (game, combo, book), take the CONSERVATIVE estimate:
# Use the method that gives the LOWER edge (less likely to overfit).
# Both methods must agree on +EV for us to bet.
conservative_edges <- all_edges %>%
  group_by(id, combo, book) %>%
  summarise(
    home_team = first(home_team),
    away_team = first(away_team),
    commence_time = first(commence_time),
    total_line = first(total_line),
    book_ml = first(book_ml),
    book_total = first(book_total),
    book_parlay_american = first(book_parlay_american),
    book_parlay_dec = first(book_parlay_dec),
    # Conservative: use the higher fair price (lower edge)
    fair_prob = max(fair_prob),
    fair_dec = min(fair_dec),
    fair_american = first(fair_american[which.max(fair_prob)]),
    cf_ak = cf[method == "answer_key"][1],
    cf_empirical = cf[method == "empirical_cf"][1],
    edge_ak = edge_pct[method == "answer_key"][1],
    edge_cf = edge_pct[method == "empirical_cf"][1],
    # Conservative edge: minimum of both methods
    edge_pct = min(edge_pct),
    n_methods = n_distinct(method),
    .groups = "drop"
  ) %>%
  mutate(
    ev = fair_prob * (book_parlay_dec - 1) - (1 - fair_prob),
    kelly_bet = ifelse(edge_pct > 0,
                       round((ev / (book_parlay_dec - 1)) * kelly_mult * bankroll, 2),
                       0)
  ) %>%
  rowwise() %>%
  mutate(
    opt = if (book == "wagerzon" & kelly_bet > 0)
            list(optimize_wagerzon_stake(kelly_bet, book_parlay_dec))
          else list(list(stake = round(kelly_bet),
                         win = round(kelly_bet) * (book_parlay_dec - 1),
                         bonus = 0)),
    opt_stake = opt$stake,
    opt_win = opt$win,
    opt_payout = opt_stake + opt_win,
    rounding_bonus = round(opt$bonus, 2)
  ) %>%
  select(-opt) %>%
  ungroup()

# Filter to +EV only (both methods must agree)
ev_bets_all <- conservative_edges %>%
  filter(edge_pct > 0) %>%
  arrange(desc(edge_pct))

# Dedup: for each (game, combo), keep the book with the higher edge
ev_bets <- ev_bets_all %>%
  group_by(id, combo) %>%
  slice_max(edge_pct, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(edge_pct))

cat(sprintf("\n=== PARLAY EDGE SUMMARY ===\n"))
cat(sprintf("Total comparisons: %d\n", nrow(conservative_edges)))
cat(sprintf("+EV parlays (before dedup): %d\n", nrow(ev_bets_all)))
cat(sprintf("+EV parlays (best book per game): %d\n", nrow(ev_bets)))

if (nrow(ev_bets) > 0) {
  cat(sprintf("Average edge: %.1f%%\n", mean(ev_bets$edge_pct)))
  cat(sprintf("Total Kelly stake: $%.2f\n", sum(ev_bets$opt_stake)))

  cat("\n=== +EV PARLAYS ===\n")
  old_width <- getOption("width")
  options(width = 160)
  ev_bets %>%
    transmute(
      game = paste0(away_team, " @ ", home_team),
      combo,
      TL = total_line,
      book,
      stake = paste0("$", sprintf("%.0f", opt_stake)),
      win = paste0("$", opt_win),
      payout = paste0("$", sprintf("%.0f", opt_payout)),
      book_odds = book_parlay_american,
      fair_odds = fair_american,
      edge = paste0(sprintf("%+.1f", edge_pct), "%")
    ) %>%
    print(n = 50)
  options(width = old_width)

  # Bet card grouped by book
  cat("\n=== BET CARD ===\n")
  for (b in sort(unique(ev_bets$book))) {
    book_bets <- ev_bets %>% filter(book == b) %>% arrange(desc(edge_pct))
    cat(sprintf("\n%s (%d bets, $%.0f total):\n", toupper(b), nrow(book_bets), sum(book_bets$opt_stake)))
    for (i in seq_len(nrow(book_bets))) {
      r <- book_bets[i, ]
      ml_side <- if (grepl("away", r$combo)) r$away_team else r$home_team
      total_side <- if (grepl("over", r$combo)) "Over" else "Under"
      cat(sprintf("  %d. %s ML + %s %.1f  ->  $%.0f  (%+.1f%% edge)\n",
                  i, ml_side, total_side, r$total_line, r$opt_stake, r$edge_pct))
    }
  }

  # Summary by combo type
  cat("\n=== EDGE BY COMBO TYPE ===\n")
  ev_bets %>%
    group_by(combo) %>%
    summarise(
      n = n(),
      avg_edge = sprintf("%.1f%%", mean(edge_pct)),
      avg_cf_empirical = sprintf("%.3f", mean(cf_empirical, na.rm = TRUE)),
      total_stake = sprintf("$%.0f", sum(opt_stake)),
      .groups = "drop"
    ) %>%
    arrange(desc(n)) %>%
    print()

  # Summary by book
  cat("\n=== EDGE BY BOOK ===\n")
  ev_bets %>%
    group_by(book) %>%
    summarise(
      n = n(),
      avg_edge = sprintf("%.1f%%", mean(edge_pct)),
      total_stake = sprintf("$%.0f", sum(opt_stake)),
      .groups = "drop"
    ) %>%
    print()
}

cat("\n=== COLLEGE BASEBALL ANSWER KEY: Complete ===\n")
timer$mark("total")
