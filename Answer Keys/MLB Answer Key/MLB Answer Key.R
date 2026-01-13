library(data.table)
library(oddsapiR)

#Answer Key----
Sys.setenv(ODDS_API_KEY = "4c490c81e9b826798390440845c89d50")
DT <- joined_dt %>% 
  rename(home_ml_odds = "consensus_prob_home",
         away_ml_odds = "consensus_prob_away",
         over_odds = "consensus_over",
         under_odds = "consensus_under",
         actual_cover = "home_winner") %>% 
  mutate(actual_over = ifelse(total_final_score > total_line,1,0))
# Helper: de-vig two American-style odds, return data.frame
devig_american <- function(odd1, odd2) {
  p1_raw <- ifelse(odd1 > 0,
                   100 / (odd1 + 100),
                   -odd1 / (-odd1 + 100))
  p2_raw <- ifelse(odd2 > 0,
                   100 / (odd2 + 100),
                   -odd2 / (-odd2 + 100))
  total_raw <- p1_raw + p2_raw
  data.frame(p1 = p1_raw/total_raw, p2 = p2_raw/total_raw)
}
american_prob <- function(odd1, odd2){
  data.frame(p1 = ifelse(odd1 > 0,
                         100 / (odd1 + 100),
                         -odd1 / (-odd1 + 100)), 
             p2 = ifelse(odd2 > 0,
                         100 / (odd2 + 100),
                         -odd2 / (-odd2 + 100)))
}


qs <- quantile(DT %>% select(home_ml_odds,away_ml_odds) %>% pivot_longer(everything(),names_to = "helper", values_to = "prob") %>% pull(prob)
                    , probs = c(0.05, 0.95))
ss <- qs[[2]] - qs[[1]]
qt <- quantile(DT$total_line,  probs = c(0.05, 0.95))
st <- qt[[2]] - qt[[1]]

distance_index <- function(dt, ps, pt, ss, st) {
  dt[, index := ((home_ml_odds - ps)/ss)^2 + ((total_line - pt)/st)^2]
}

# 1) Mean‐matching loop to select initial N games
library(data.table)

# Function 1: mean‐matching to select initial N games
mean_match <- function(dt, N, parent_spread, parent_total,
                       ss, st, max_iter_mean, tol_mean) {
  dt <- copy(dt)
  adj_spread <- parent_spread
  adj_total  <- parent_total
  for (iter in seq_len(max_iter_mean)) {
    distance_index(dt, adj_spread, adj_total, ss, st)
    setorder(dt, index)
    dt[, included := FALSE]
    dt[1:N, included := TRUE]
    
    # compute current means and errors
    mean_s <- dt[included == TRUE, mean(home_ml_odds)]
    mean_t <- dt[included == TRUE, mean(total_line)]
    err_s  <- mean_s   - parent_spread
    err_t  <- mean_t   - parent_total
    
    # stop if within tolerance
    if (abs(err_s) < tol_mean && abs(err_t) < tol_mean) break
    
    # adjust the parent targets
    adj_spread <- adj_spread - err_s
    adj_total  <- adj_total  - err_t
  }
  # return updated dt and targets
  list(dt = dt,
       parent_spread = parent_spread,
       parent_total  = parent_total)
}

# Function 2: balance sample by greedy remove/add + shrink‐on‐stall
balance_sample <- function(dt, N, target_cover, target_over, tol_error) {
  dt <- copy(dt)
  # initialize errors
  cover_error <- dt[included == TRUE, sum(actual_cover)] -
    round(target_cover * N)
  over_error  <- dt[included == TRUE, sum(actual_over )] -
    round(target_over  * N)
  M <- nrow(dt)
  
  repeat {
    # stop if within tolerance
    if (abs(cover_error) <= tol_error && abs(over_error) <= tol_error) break
    
    # mark failures
    removal_failed  <- TRUE
    addition_failed <- TRUE
    
    # ---- remove one game from 1..N ----
    in1toN    <- dt[1:N,   included]
    cov_help1 <- (cover_error > 0 & dt[1:N, actual_cover] == 1) |
      (cover_error < 0 & dt[1:N, actual_cover] == 0)
    ov_help1  <- (over_error  > 0 & dt[1:N, actual_over ] == 1) |
      (over_error  < 0 & dt[1:N, actual_over ] == 0)
    
    both1   <- which(in1toN & cov_help1 & ov_help1)
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
      cover_error <- dt[included == TRUE, sum(actual_cover)] -
        round(target_cover * N)
      over_error  <- dt[included == TRUE, sum(actual_over )] -
        round(target_over  * N)
    }
    
    # ---- add one game from (N+1)..M ----
    inNplus   <- dt[(N+1):M, included]
    cov_help2 <- (cover_error > 0 & dt[(N+1):M, actual_cover] == 0) |
      (cover_error < 0 & dt[(N+1):M, actual_cover] == 1)
    ov_help2  <- (over_error  > 0 & dt[(N+1):M, actual_over ] == 0) |
      (over_error  < 0 & dt[(N+1):M, actual_over ] == 1)
    
    both2   <- which(!inNplus & cov_help2 & ov_help2)
    either2 <- which(!inNplus & (cov_help2 | ov_help2))
    
    if (length(both2) > 0) {
      j <- max(both2) + N
      dt[j, included := TRUE]
      addition_failed <- FALSE
    } else if (length(either2) > 0) {
      j <- max(either2) + N
      dt[j, included := TRUE]
      addition_failed <- FALSE
    }
    
    if (!addition_failed) {
      cover_error <- dt[included == TRUE, sum(actual_cover)] -
        round(target_cover * N)
      over_error  <- dt[included == TRUE, sum(actual_over )] -
        round(target_over  * N)
    }
    
    # ---- if neither helped, shrink sample size ----
    if (removal_failed && addition_failed) {
      worst_i <- dt[1:N, which.max(index)]
      dt[worst_i, included := FALSE]
      N <- N - 1
      cover_error <- dt[included == TRUE, sum(actual_cover)] -
        round(target_cover * N)
      over_error  <- dt[included == TRUE, sum(actual_over )] -
        round(target_over  * N)
      next
    }
  }
  
  # return final dt and stats
  list(dt          = dt,
       final_N     = N,
       cover_error = cover_error,
       over_error  = over_error)
}


parent_spread <- 0.6
parent_total  <- 8.5
home_sp_odds <- -150; away_sp_odds <- 130
over_odds    <- -110; under_odds    <- +100
N <- round(nrow(DT)*0.05,0)

matched_mean<- mean_match(DT,N,parent_spread,parent_total,ss,st,max_iter_mean = 500, tol_mean = 0.005)
# mean(matched_mean$dt[included == TRUE,home_ml_odds]) # a check to make sure the mean matches

desired_over_prob<- devig_american(over_odds,under_odds)[[1]]
desired_cover_prob<- devig_american(home_sp_odds,away_sp_odds)[[1]]
final_sample <- balance_sample(matched_mean$dt,N,desired_cover_prob,desired_over_prob,tol_error = 1)


#Algoirthm Checks ----
#  mean(matched_mean$dt[included == TRUE,total_line]) checks
dt_final   <- final_sample$dt
final_N    <- final_sample$final_N
cover_err  <- final_sample$cover_error
over_err   <- final_sample$over_error

# 2) Compute final metrics
final_mean_s  <- dt_final[included == TRUE, mean(home_ml_odds)]
final_mean_t  <- dt_final[included == TRUE, mean(total_line)]
final_cover_r <- dt_final[included == TRUE, sum(actual_cover)] / final_N
final_over_r  <- dt_final[included == TRUE, sum(actual_over )] / final_N

# 3) Print everything alongside your targets
cat("===== Sanity Check =====\n")
cat(sprintf("Sample size:      %d  (target was %d)\n", final_N, N))
cat(sprintf("Spread mean:      %.3f  (target %.3f; err %.3f)\n",
            final_mean_s, parent_spread, final_mean_s - parent_spread))
cat(sprintf("Total mean:       %.3f  (target %.3f; err %.3f)\n",
            final_mean_t, parent_total, final_mean_t - parent_total))
cat(sprintf("Cover rate:       %.3f  (target %.3f; err %d games)\n",
            final_cover_r, desired_cover_prob, cover_err))
cat(sprintf("Over rate:        %.3f  (target %.3f; err %d games)\n",
            final_over_r, desired_over_prob, over_err))

#Predictions ----

# get game odds----
game_odds <- toa_sports_odds(
  sport_key = "baseball_mlb",
  regions    = "us,eu",
  markets    = "h2h,totals",
  odds_format = "american",
  date_format = "iso")

book_weights <- tibble::tibble(
  bookmaker_key = c(
    "pinnacle", "circasports", "betonlineag", "lowvig", "matchbook", "bookmaker.eu", "coolbet",
    "everygame", "intertops", "unibet", "unibet_us", "unibet_it", "betmgm", "williamhill_us",
    "caesars", "fanduel", "draftkings", "pointsbetus", "betrivers", "bovada", "sport888", "superbook",
    "wynnbet", "sugarhouse", "betus", "mybookieag", "foxbet", "barstool", "fanatics", "gtbets",
    "tipico_de", "twinspires", "livescorebet_eu", "nordicbet", "betsson", "onexbet"
  ),
  weight = c(
    1.00, 0.95, 0.92, 0.90, 0.88, 0.88, 0.85, 0.80, 0.80, 0.80, 0.75, 0.75, 0.72, 0.72,
    0.70, 0.68, 0.68, 0.65, 0.60, 0.55, 0.52, 0.52, 0.50, 0.50, 0.48, 0.45, 0.45, 0.45, 0.40, 0.40,
    0.38, 0.38, 0.35, 0.35, 0.35, 0.38
  )
)

ml_consensus <- game_odds %>% 
  filter(market_key == "h2h") %>% 
  mutate(odds_type = ifelse(home_team == outcomes_name, "home", ifelse(away_team == outcomes_name, "away", outcomes_name))) %>% 
  pivot_wider(id_cols = c(id,commence_time,home_team, away_team, bookmaker), names_from = odds_type, values_from = outcomes_price, names_prefix = "odds_") %>% 
  mutate(as_tibble(american_prob(odds_home, odds_away))) %>% #consider devigging which can be done by switching devig_american()
  rename(prob_home = p1, prob_away = p2) %>% 
  group_by(id) %>% 
  summarize(
    home_team = first(home_team),
    away_team = first(away_team),
    consensus_prob_home = median(prob_home, na.rm = TRUE),
    consensus_prob_away = median(prob_away, na.rm = TRUE)
  ) %>%
  ungroup()


consensus_total <- game_odds %>% 
  filter(market_key =="totals") %>% 
  mutate(odds_type = ifelse(home_team == outcomes_name, "home", ifelse(away_team == outcomes_name, "away", outcomes_name))) %>% 
  pivot_wider(id_cols = c(id,commence_time,home_team, away_team, bookmaker_key,outcomes_point), names_from = odds_type, values_from = outcomes_price, names_prefix = "odds_") %>%  
  mutate(as_tibble(american_prob(odds_Over, odds_Under))) %>% 
  rename(prob_over = p1, prob_under = p2, total_line = outcomes_point) %>% 
  left_join(book_weights, by = "bookmaker_key") %>% 
  group_by(id,total_line) %>% 
  mutate(total_weight = sum(weight, na.rm = T)) %>% 
  ungroup() %>% 
  group_by(id) %>% 
  filter(total_weight == max(total_weight, na.rm = TRUE)) %>% 
  group_by(id,total_line) %>% 
  summarize(consensus_over = median(prob_over, na.rm = T),
            consensus_under = median(prob_under,na.rm = T))

mlb_odds <- ml_consensus %>% inner_join(consensus_total, by = "id")

#Generate Side predictions----
run_means_for_id <- function(id, parent_spread, parent_total, target_cover, target_over,
                             DT, ss, st, N,
                             max_iter_mean = 500, tol_mean = 0.005, tol_error = 1) {
  mm  <- mean_match(DT, N, parent_spread, parent_total, ss, st, max_iter_mean, tol_mean)
  bal <- balance_sample(mm$dt, N, target_cover, target_over, tol_error)
  
  inc_df <- bal$dt %>% filter(included == TRUE)
  
  # summarise all game_home_ml* columns
  bets_summary <- inc_df %>%
    summarise(across(starts_with("game_home_ml"), ~ mean(.x, na.rm = TRUE)))
  
  bets_summary %>%
    select(everything())
}

targets <- mlb_odds %>%
  transmute(
    id,
    parent_spread = consensus_prob_home,
    parent_total  = total_line,
    target_cover  = consensus_prob_home,
    target_over   = consensus_over
  )

prob_to_american <- function(prob) {
  ifelse(prob >= 0.5,
         # For favorites
         - (prob / (1 - prob)) * 100,
         # For underdogs
         ((1 - prob) / prob) * 100
  ) %>% round(0)  # round to nearest integer
}
final_bets <- targets %>%
  mutate(res = pmap(
    list(id, parent_spread, parent_total, target_cover, target_over),
    ~ run_means_for_id(..1, ..2, ..3, ..4, ..5,
                       DT = DT, ss = ss, st = st, N = N)
  )) %>%
  unnest(res) %>% 
  inner_join(game_odds %>% group_by(id) %>% 
               summarize(home_team = first(home_team),away_team = first(away_team),commence_time = first(commence_time)), by = "id") %>% 
  select(home_team, away_team, commence_time,everything())

predictions <- final_bets %>% 
  mutate(across(starts_with("game_home_ml"),
                ~ prob_to_american(.x),
                .names = "{.col}_american")) %>% 
  relocate(ends_with("_american"), .before = starts_with("game_home_ml"))

#Compare Sides to Actual Odds----

bet_inning <- 5
# market <- "alternate_spreads_1st_5_innings"
market <- "h2h_1st_5_innings" # change this to the market you want to fetch odds for


# 1) Get upcoming MLB events (you can pass regions here, but markets aren’t needed yet)
ev <- GET(
  "https://api.the-odds-api.com/v4/sports/baseball_mlb/events",
  query = list(apiKey = api_key, regions = "us")   # us2 is valid; combine with us if you want more books
)
stop_for_status(ev)
events <- fromJSON(content(ev, "text"), flatten = TRUE) %>% as_tibble() %>% 
  mutate(pt_start_time = with_tz(ymd_hms(commence_time, tz = "UTC"), tzone = "America/Los_Angeles")) %>% 
  filter(pt_start_time > Sys.time() &  date(with_tz(pt_start_time, "America/Los_Angeles")) == date(with_tz(Sys.time(), "America/Los_Angeles"))) #filter to games yet to start

fetch_event_odds <- function(event_id,market) {
  res <- GET(
    paste0("https://api.the-odds-api.com/v4/sports/baseball_mlb/events/", event_id, "/odds"),
    query = list(
      apiKey     = Sys.getenv("ODDS_API_KEY"),
      regions    = "us,us2,us_ex",
      markets    = market,  # change this to whichever market you want
      oddsFormat = "american",
      dateFormat = "iso"
    )
  )
  
  if (http_error(res)) {
    warning(paste("Failed for event:", event_id))
    return(NULL)
  }
  
  fromJSON(content(res, "text"), flatten = TRUE) %>% as_tibble()
}

# Step 3: Loop over all event IDs
all_odds <- map_dfr(events$id, ~fetch_event_odds(.x,market))


flat_betting_odds <- all_odds %>%
  # 1. one row per bookie
  unnest_longer(bookmakers) %>% 
  unnest_wider (bookmakers,   names_sep = "_") %>%
  # 2. one row per market (h2h, totals, …)
  unnest_longer(bookmakers_markets) %>%
  unnest_wider (bookmakers_markets, names_sep = "_") %>%
  # 3. one row per single outcome (e.g. “Detroit Tigers” @ –113)
  unnest_longer(bookmakers_markets_outcomes) %>%
  unnest_wider (bookmakers_markets_outcomes, names_sep = "_") %>%
  # 4. parse your datetimes
  mutate(
    commence_time       = ymd_hms(commence_time,               tz="UTC"),
    market_update       = ymd_hms(bookmakers_markets_last_update, tz="UTC")
  ) %>%
  rename(
    bookmaker_key   = bookmakers_key,
    bookmaker_title = bookmakers_title,
    market_key      = bookmakers_markets_key,
    outcome_name    = bookmakers_markets_outcomes_name,
    closing_odds    = bookmakers_markets_outcomes_price
  ) %>% 
  unnest_wider(outcome_name, names_sep = c("_")) %>% 
  unnest_wider(closing_odds, names_sep = c("_")) %>%
  mutate(
    home_odds = ifelse(home_team == outcome_name_1, closing_odds_1,closing_odds_2),
    away_odds = ifelse(away_team == outcome_name_1, closing_odds_1,closing_odds_2),
  ) %>% 
  group_by(
    id, commence_time, home_team, away_team,
    bookmaker_key, bookmaker_title
  ) %>% 
  summarise(
    book_home_market := first(home_odds),
    book_away_market := first(away_odds)
  )
#Final Side Results----
prediction_set <- flat_betting_odds %>% left_join(predictions %>% select(id,contains(paste("_",bet_inning,sep = ""))), by = "id") %>%
  left_join(mlb_odds %>% select(-home_team,-away_team), by = "id") %>% 
  rename("book_full_game_home_prob" = "consensus_prob_home",
         "book_full_game_away_prob" = "consensus_prob_away",
         "book_full_game_over_prob" = "consensus_over",
         "book_full_game_under_prob" = "consensus_under",
         "book_full_game_total_line" = "total_line") %>% 
  rename(home_predicted_prob = paste0("game_home_ml_inning_",bet_inning,sep = ""),
         home_predicted_american_odds = paste0("game_home_ml_inning_",bet_inning,"_american")) %>% 
  mutate(away_predicted_prob = 1 - home_predicted_prob) %>% 
  mutate(as_tibble(american_prob(book_home_market,book_away_market))) %>% 
  rename(book_market_prob_home = p1, book_market_prob_away = p2) %>% 
  mutate(home_ev = round(home_predicted_prob * ((1/book_market_prob_home)-1) - (1-home_predicted_prob),8),
         away_ev = round(away_predicted_prob * ((1/book_market_prob_away)-1) - (1-away_predicted_prob),3)) %>% 
  mutate(market = market)
bankroll <- 200
kelly_mult <- 0.25

my_bets <- prediction_set %>% filter(bookmaker_key %in% c("fliff","rebet","novig","prophetx")) %>% 
  mutate(pt_start_time = with_tz(ymd_hms(commence_time, tz = "UTC"), tzone = "America/Los_Angeles")) %>% 
  ungroup() %>% 
  select(-id,-commence_time,-bookmaker_title) %>% 
  relocate(c(pt_start_time, bookmaker_key,market,home_ev,away_ev),.after = "away_team") %>% 
  relocate(c(away_predicted_prob,book_market_prob_home,book_market_prob_away), .before = "book_full_game_home_prob") %>% #pivot_wider?
  mutate(home_bet_size = round((home_ev/((1/book_market_prob_home)-1)*kelly_mult)*bankroll,2),
         away_bet_size = round((away_ev/((1/book_market_prob_away)-1)*kelly_mult)*bankroll,2)) %>% 
  relocate(home_bet_size, .after = book_home_market) %>% relocate(away_bet_size, .after = book_away_market)
#Compare Total to Actual Odds----

bet_inning <- 5

total_market <- "alternate_totals_1st_5_innings" # change this to the market you want to fetch odds for
#alternate_totals_1st_5_innings


# Step 3: Loop over all event IDs
total_all_odds <- map_dfr(events$id, ~ fetch_event_odds(.x, total_market)) #uncomment whenver just to hold value


total_flat_betting_odds <- total_all_odds %>%
  # 1. one row per bookie
  unnest_longer(bookmakers) %>% 
  unnest_wider (bookmakers,   names_sep = "_") %>%
  # 2. one row per market (h2h, totals, …)
  unnest_longer(bookmakers_markets) %>%
  unnest_wider (bookmakers_markets, names_sep = "_") %>%
  # 3. one row per single outcome (e.g. “Detroit Tigers” @ –113)
  unnest_longer(bookmakers_markets_outcomes) %>%
  unnest_wider (bookmakers_markets_outcomes, names_sep = "_") %>%
  # 4. parse your datetimes
  mutate(
    commence_time       = ymd_hms(commence_time,               tz="UTC"),
    market_update       = ymd_hms(bookmakers_markets_last_update, tz="UTC")
  ) %>%
  rename(
    bookmaker_key   = bookmakers_key,
    bookmaker_title = bookmakers_title,
    market_key      = bookmakers_markets_key,
    outcome_name    = bookmakers_markets_outcomes_name,
    closing_odds    = bookmakers_markets_outcomes_price
  ) %>% 
  unnest_wider(outcome_name, names_sep = c("_")) %>% 
  unnest_wider(closing_odds, names_sep = c("_")) %>%
  unnest_wider(bookmakers_markets_outcomes_point, names_sep = c("_")) %>% 
  mutate(
  ) %>% 
  group_by(
    id, commence_time, home_team, away_team,
    bookmaker_key, bookmaker_title
  ) %>% 
  summarise(
    book_over_market := first(closing_odds_1),
    book_under_market := first(closing_odds_2),
    book_market_total = first(bookmakers_markets_outcomes_point_1)
  )
#Generate Total predictions----
run_means_for_id_total <- function(id, parent_spread, parent_total, target_cover, target_over, totals, inning,
                             DT, ss, st, N,
                             max_iter_mean = 500, tol_mean = 0.005, tol_error = 1) {
  mm  <- mean_match(DT, N, parent_spread, parent_total, ss, st, max_iter_mean, tol_mean)
  bal <- balance_sample(mm$dt, N, target_cover, target_over, tol_error)
  
  inc_df <- bal$dt %>% filter(included == TRUE)
  
  # summarise all game_home_ml* columns
  vals <- inc_df %>%
    select( "inning_total" = !!sym(paste0("full_game_total_inning_",inning, sep = ""))) %>%
    pull()
  pct_cols <- set_names(
    map(totals, ~ mean(vals > .x, na.rm = TRUE)),
    paste0("pct_over_", gsub("\\.", "_", totals))
  )
  
  tibble(!!!pct_cols)
  
}

targets <- mlb_odds %>%
  transmute(
    id,
    parent_spread = consensus_prob_home,
    parent_total  = total_line,
    target_cover  = consensus_prob_home,
    target_over   = consensus_over
  )

total_final_bets <- targets %>%
  mutate(res = pmap(
    list(id, parent_spread, parent_total, target_cover, target_over),
    ~ run_means_for_id_total(..1, ..2, ..3, ..4, ..5, totals = unique(total_flat_betting_odds$book_market_total), inning = bet_inning,
                       DT = DT, ss = ss, st = st, N = N)
  )) %>%
  unnest(res) %>% 
  inner_join(game_odds %>% group_by(id) %>% 
               summarize(home_team = first(home_team),away_team = first(away_team),commence_time = first(commence_time)), by = "id") %>% 
  select(home_team, away_team, commence_time,everything())

total_predictions <- total_final_bets %>% 
  pivot_longer(starts_with("pct_over"),names_to = "market_total_line", values_to = "over_prediction") %>% 
  mutate(market_total_line = as.numeric(gsub("_",".",sub("^pct_over_", "", market_total_line)))) %>%  #get it ready to join lines
  mutate(across(starts_with("pct_over"),
                ~ prob_to_american(.x),
                .names = "{.col}_american"))

#Final Total Results----
total_prediction_set <- total_flat_betting_odds %>% left_join(total_predictions %>% select(id,market_total_line,over_prediction), by = c("id", "book_market_total" = "market_total_line")) %>%
  left_join(mlb_odds %>% select(-home_team,-away_team), by = "id") %>% 
  rename("book_full_game_home_prob" = "consensus_prob_home",
         "book_full_game_away_prob" = "consensus_prob_away",
         "book_full_game_over_prob" = "consensus_over",
         "book_full_game_under_prob" = "consensus_under",
         "book_full_game_total_line" = "total_line") %>% 
  mutate(under_prediction = 1 - over_prediction) %>% 
  mutate(as_tibble(american_prob(book_over_market,book_under_market))) %>% 
  rename(book_market_prob_over = p1, book_market_prob_under = p2) %>% 
  mutate(over_ev = round(over_prediction * ((1/book_market_prob_over)-1) - (1-over_prediction),3),
         under_ev = round(under_prediction * ((1/book_market_prob_under)-1) - (1-under_prediction),3)) %>% 
  mutate(market = total_market) %>% 
  mutate(over_bet_size = round((over_ev/((1/book_market_prob_over)-1)*kelly_mult)*bankroll,2),
         under_bet_size = round((under_ev/((1/book_market_prob_under)-1)*kelly_mult)*bankroll,2))

total_my_bets <- total_prediction_set %>% filter(bookmaker_key %in% c("fliff","rebet","novig","prophetx")) %>% 
  mutate(pt_start_time = with_tz(ymd_hms(commence_time, tz = "UTC"), tzone = "America/Los_Angeles")) %>% 
  ungroup() %>% 
  select(-id,-commence_time,-bookmaker_title) %>% 
  relocate(c(pt_start_time, bookmaker_key,market,over_ev,under_ev),.after = "away_team") %>% 
  relocate(over_bet_size, .after = book_over_market) %>% relocate(under_bet_size, .after = book_under_market)
  # relocate(c(under_prediction,book_market_over,book_market_under), .before = "book_full_game_home_prob")
