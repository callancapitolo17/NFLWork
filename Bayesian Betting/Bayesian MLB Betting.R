# =======================================================
# 1. SETUP: INSTALL & LOAD REQUIRED PACKAGES
# =======================================================
required_packages <- c("baseballr", "dplyr", "lubridate", "httr", "jsonlite", "rstan", "tidyverse", "pbapply", "tidygeocoder")
for(pkg in required_packages){
    if (!require(pkg, character.only = TRUE)) install.packages(pkg)
}
library(baseballr)
library(dplyr)
library(lubridate)
library(httr)
library(jsonlite)
library(rstan)
library(tidyverse)
library(pbapply)
library(tidygeocoder)

library(dplyr)

# Load Retrosheet play-by-play data
all2023 <- read_csv("all2023.csv")

# Load Retrosheet player roster (downloaded separately from Retrosheet)
roster <- read_csv("allplayers.csv", col_names = c("PLAYER_ID", "LAST_NAME", "FIRST_NAME", "BAT_HAND", "THROW_HAND", "TEAM", "POS"))

# Create full player name
roster <- roster %>%
  mutate(PLAYER_NAME = paste(FIRST_NAME, LAST_NAME)) %>%
  select(PLAYER_ID, PLAYER_NAME) %>% 
  group_by(PLAYER_ID) %>% 
  summarize(PLAYER_NAME = first(PLAYER_NAME))

# Merge player names into Retrosheet dataset
all2023 <- all2023 %>%
  left_join(roster %>% rename("Batter_Name" = PLAYER_NAME), by = c("RESP_BAT_ID" = "PLAYER_ID")) %>% 
  left_join(roster %>% rename("Pitcher_Name" = PLAYER_NAME), by = c("RESP_PIT_ID" = "PLAYER_ID"))

# Team mapping based on Retrosheet abbreviations
team_mapping <- c(
  "TOR" = "Toronto Blue Jays", "WAS" = "Washington Nationals", "KCA" = "Kansas City Royals",
  "OAK" = "Oakland Athletics", "TEX" = "Texas Rangers", "HOU" = "Houston Astros",
  "MIN" = "Minnesota Twins", "BOS" = "Boston Red Sox", "MIA" = "Miami Marlins",
  "CHN" = "Chicago Cubs", "ATL" = "Atlanta Braves", "NYN" = "New York Mets",
  "PHI" = "Philadelphia Phillies", "MIL" = "Milwaukee Brewers", "SLN" = "St. Louis Cardinals",
  "PIT" = "Pittsburgh Pirates", "CIN" = "Cincinnati Reds", "COL" = "Colorado Rockies",
  "SDN" = "San Diego Padres", "LAN" = "Los Angeles Dodgers", "SFN" = "San Francisco Giants",
  "ARI" = "Arizona Diamondbacks", "BAL" = "Baltimore Orioles", "NYA" = "New York Yankees",
  "TBA" = "Tampa Bay Rays", "CHA" = "Chicago White Sox", "DET" = "Detroit Tigers",
  "SEA" = "Seattle Mariners", "ANA" = "Los Angeles Angels", "CLE" = "Cleveland Guardians"
)

# Convert inning format to match MLB_dataRaw.csv
all2023 <- all2023 %>%
  mutate(AWAY_TEAM_NAME = team_mapping[AWAY_TEAM_ID]) %>% 
  mutate(inning = paste0(INN_CT, ifelse(BAT_HOME_ID == 0, "T", "B"))) %>% 
  mutate(HOME_TEAM_ID = substr(GAME_ID, 1, 3),
         HOME_TEAM_NAME = team_mapping[HOME_TEAM_ID],
         batting_team = ifelse(BAT_HOME_ID == 0, AWAY_TEAM_NAME,HOME_TEAM_NAME)) %>% 
  arrange(GAME_ID, INN_CT, BAT_HOME_ID, OUTS_CT) %>% #unsure if this is the right way to group it. Fine for now when not really trying to assign player credit
  group_by(GAME_ID) %>%
  mutate(total_score =HOME_SCORE_CT +AWAY_SCORE_CT,
         run_outcome = lead(total_score, default = last(total_score)) - total_score) %>% 
  ungroup()

# Identify unique teams and pitchers per game
game_info <- all2023 %>%
  group_by(GAME_ID) %>%
  mutate(date = as.Date(substr(GAME_ID,4,11), format ="%Y%m%d")) %>% 
  summarise(
    date = first(date),
    away_team = first(AWAY_TEAM_NAME),
    home_team = first(HOME_TEAM_NAME),
    home_starting_pitcher = first(Pitcher_Name),
    away_starting_pitcher = first(Pitcher_Name[BAT_HOME_ID == 1])
  ) %>% 
  mutate(double_header = ifelse(substr(GAME_ID,12,12) > 0, 1,0))

# Extract final scores per game
score_info <- all2023 %>%
  group_by(GAME_ID) %>%
  summarise(
    home_true_score = max(HOME_SCORE_CT),
    away_true_score = max(AWAY_SCORE_CT)
  )

count_runs <- function(df, team_column, team_name) {
  df %>%
    filter(!!sym(team_column) == team_name) %>%
    group_by(GAME_ID) %>%
    summarise(
      run_1 = sum(run_outcome == 1, na.rm = TRUE),
      run_2 = sum(run_outcome == 2, na.rm = TRUE),
      run_3 = sum(run_outcome == 3, na.rm = TRUE),
      run_4 = sum(run_outcome == 4, na.rm = TRUE)
    )
}

# Initialize lists to store the results
away_runs_list <- list()
home_runs_list <- list()

# Loop through each game to count the runs
for (game in unique(all2023$GAME_ID)) {
  away_team <- all2023$AWAY_TEAM_NAME[all2023$GAME_ID == game]
  home_team <- all2023$HOME_TEAM_NAME[all2023$GAME_ID == game]
  
  away_runs <- count_runs(all2023 %>% filter(GAME_ID == game), "batting_team", away_team)
  home_runs <- count_runs(all2023 %>% filter(GAME_ID == game), "batting_team", home_team)
  
  away_runs_list[[game]] <- away_runs
  home_runs_list[[game]] <- home_runs
}

# Combine all results into single data frames
away_runs <- bind_rows(away_runs_list)
home_runs <- bind_rows(home_runs_list)
final_data <- game_info %>%
  left_join(away_runs, by = "GAME_ID") %>%
  rename(away_run_1 = run_1, away_run_2 = run_2, away_run_3 = run_3, away_run_4 = run_4) %>%
  left_join(home_runs, by = "GAME_ID") %>%
  rename(home_run_1 = run_1, home_run_2 = run_2, home_run_3 = run_3, home_run_4 = run_4) %>%
  left_join(score_info, by = "GAME_ID")

check <- final_data %>% 
  mutate(home_runs_check = home_run_1+home_run_2*2+home_run_3*3+home_run_4*4,
   away_runs_check = away_run_1+away_run_2*2+away_run_3*3+away_run_4*4,
   home_check = home_true_score == home_runs_check,
   away_check = away_runs_check == away_true_score)

help <- all2023 %>% 
  filter(GAME_ID == "ANA202304120") %>% 
  select(AWAY_SCORE_CT,HOME_SCORE_CT,RBI_CT,EVENT_TX,run_outcome)

parse_event_tx <- function(event_tx) {
  case_when(
    str_detect(event_tx, "K") ~ "Strikeout",
    str_detect(event_tx, "W") ~ "Walk",
    str_detect(event_tx, "IW") ~ "Intentional Walk",
    str_detect(event_tx, "HBP") ~ "Hit By Pitch",
    str_detect(event_tx, "E[1-9]") ~ "Error",
    str_detect(event_tx, "FC") ~ "Fielder's Choice",
    str_detect(event_tx, "S[1-9]?") ~ "Single",
    str_detect(event_tx, "D[1-9]?") ~ "Double",
    str_detect(event_tx, "T[1-9]?") ~ "Triple",
    str_detect(event_tx, "HR") ~ "Home Run",
    str_detect(event_tx, "DP") ~ "Double Play",
    str_detect(event_tx, "TP") ~ "Triple Play",
    str_detect(event_tx, "SH") ~ "Sacrifice Hit",
    str_detect(event_tx, "SF") ~ "Sacrifice Fly",
    str_detect(event_tx, "BK") ~ "Balk",
    str_detect(event_tx, "PB") ~ "Passed Ball",
    str_detect(event_tx, "WP") ~ "Wild Pitch",
    str_detect(event_tx, "CS") ~ "Caught Stealing",
    str_detect(event_tx, "PO") ~ "Pickoff",
    str_detect(event_tx, "DI") ~ "Defensive Indifference",
    TRUE ~ "Unknown Play"
  )
}

# =======================================================
# 2. GET MLB DATA FROM ONLINE SOURCES
# =======================================================
# 2.1 Pull MLB Game Schedule & Results (using mlb_schedule)
teams <- c("ANA", "ARI", "ATL", "BAL", "BOS", "CHA", "CHN", "CIN", "CLE", "COL",
           "DET", "HOU", "KCA", "LAN", "MIA", "MIL", "MIN", "NYA", "NYN", "OAK",
           "PHI", "PIT", "SDN", "SEA", "SFN", "SLN", "TBA", "TEX", "TOR", "WAS")

pbp_games <- list()  # Create an empty list to store results

for (team in teams) {
  pbp_games[[team]] <- getRetrosheet("play", 2023, team)  # Store data for each team
}
season <- 2023
game_logs <- mlb_schedule(season = season) 
game_logs <- game_logs %>%
    filter(game_type %in% c("R","W", "L","D","F")) %>%
    mutate(game_date = as.Date(game_date),
           game_key = paste(game_pk, game_date, calendar_event_id, sep = "_")) #fix join issues

# Assume game_logs has a column "game_date" and (if available) "game_datetime"
# If game start time is available, extract the hour (adjust the timezone as needed)
if ("game_datetime" %in% names(game_logs)) {
    game_logs <- game_logs %>%
        mutate(game_datetime = ymd_hms(game_datetime, tz = "America/New_York"),
               game_start_hour = hour(game_datetime))
} else {
    # Otherwise, you can set a default start time (e.g., 13 for 1 PM)
    game_logs <- game_logs %>%
        mutate(game_start_hour = 13)
}

library(pbapply)
library(dplyr)

# Example retry function for one interval:
fetch_interval_data <- function(s_date, e_date, max_attempts = 3) {
    attempt <- 1
    while(attempt <= max_attempts) {
        message("Querying from ", s_date, " to ", e_date, " (Attempt ", attempt, ")")
        result <- tryCatch({
            chunk <- scrape_statcast_savant_pitcher_all(start_date = s_date, end_date = e_date)
            # If no data is found, return NULL so we can retry
            if(nrow(chunk) == 0) stop("No valid data found")
            return(chunk)
        }, error = function(e) {
            message("Attempt ", attempt, " failed: ", e$message)
            return(NULL)
        })
        if(!is.null(result)) {
            return(result)
        }
        Sys.sleep(1)  # Wait a second before retrying
        attempt <- attempt + 1
    }
    warning("No valid data found for interval ", s_date, " to ", e_date)
    return(NULL)
}

# Now apply this function over your date intervals:
days_per_chunk <- 3  # or use a shorter interval if necessary
start_date <- min(game_logs$game_date)
end_date   <- max(game_logs$game_date)

interval_starts <- seq(start_date, end_date, by = paste(days_per_chunk, "days"))
interval_ends   <- interval_starts + days(days_per_chunk - 1)
interval_ends[interval_ends > end_date] <- end_date

pbp_list <- pbapply::pblapply(seq_along(interval_starts), function(i) {
    s_date <- as.character(interval_starts[i])
    e_date <- as.character(interval_ends[i])
    fetch_interval_data(s_date, e_date)
})


# Bind the rows (filtering out NULLs)
pbp_list <- pbp_list[!sapply(pbp_list, is.null)]
pbp_data <- bind_rows(pbp_list)



# NOTE: This example assumes pbp_data includes a column "run_outcome" (with values 1, 2, 3, or 4)
# and relevant team/player columns.

# =======================================================
# 3. GEOCODE STADIUMS & FETCH WEATHER DATA
# =======================================================
# 3.1 Automatically geocode all unique stadiums using tidygeocoder
unique_stadiums <- game_logs %>% distinct(venue_name)
stadium_coords <- unique_stadiums %>%
    geocode(venue_name, method = "osm", lat = latitude, long = longitude)
game_logs <- game_logs %>% left_join(stadium_coords, by = "venue_name")

# 3.2 Define a function to fetch weather during a specified time window.
#     Here we use lubridate to correctly parse times.
get_weather_openmeteo_game_time <- function(lat, lon, date, start_hour, end_hour) {
    url <- paste0(
        "https://archive-api.open-meteo.com/v1/archive?",
        "latitude=", lat,
        "&longitude=", lon,
        "&start_date=", date,
        "&end_date=", date,
        "&hourly=temperature_2m,relativehumidity_2m,windspeed_10m",
        "&timezone=America/New_York"
    )
    resp <- GET(url)
    raw_text <- content(resp, "text", encoding = "UTF-8")
    dat <- fromJSON(raw_text)

    if (is.null(dat$hourly) || length(dat$hourly$time) == 0) {
        warning("No hourly data returned for date: ", date)
        return(tibble(date = as.Date(date), temperature = NA_real_, humidity = NA_real_, wind_speed = NA_real_))
    }

    # Parse times using lubridate to ensure correct conversion from ISO8601
    times <- ymd_hm(dat$hourly$time, tz = "America/New_York")
    hours <- as.numeric(format(times, "%H"))
    message("Available hours for ", date, ": ", paste(hours, collapse = ", "))

    idx <- which(hours >= start_hour & hours <= end_hour)
    if (length(idx) == 0) {
        warning("No hourly data within game time for date: ", date, ". Available hours: ", paste(hours, collapse = ", "))
        return(tibble(date = as.Date(date), temperature = NA_real_, humidity = NA_real_, wind_speed = NA_real_))
    }

    tibble(
        date = as.Date(date),
        temperature = mean(dat$hourly$temperature_2m[idx], na.rm = TRUE),
        humidity = mean(dat$hourly$relativehumidity_2m[idx], na.rm = TRUE),
        wind_speed = mean(dat$hourly$windspeed_10m[idx], na.rm = TRUE)
    )
}

# Compute weather data for each game, storing the game_pk along with venue_name and date.
weather_list <- pbapply::pbsapply(1:nrow(game_logs), function(i) {
    this_date <- as.character(game_logs$game_date[i])
    lat <- game_logs$latitude[i]
    lon <- game_logs$longitude[i]
    venue <- game_logs$venue_name[i]
    game_id <- game_logs$game_pk[i]
    calendar_event_id <- game_logs$calendar_event_id[i]
    start_hr <- ifelse(is.na(game_logs$game_start_hour[i]), 13, game_logs$game_start_hour[i])
    end_hr <- start_hr + 3             # assuming a 3-hour game window

    if (is.na(lat) || is.na(lon)) {
        return(tibble(
            game_pk = game_id,
            date = as.Date(this_date),
            temperature = NA_real_,
            humidity = NA_real_,
            wind_speed = NA_real_,
            venue_name = venue
        ))
    } else {
        tib <- get_weather_openmeteo_game_time(lat, lon, this_date, start_hour = start_hr, end_hour = end_hr)
        tib <- tib %>% mutate(venue_name = venue, game_pk = game_id, calendar_event_id = calendar_event_id)
        return(tib)
    }
}, simplify = FALSE)

weather_data <- bind_rows(weather_list)
clean_weather_data <- weather_data %>%
    mutate(game_date = as.Date(date),
           game_key = paste(game_pk, game_date,calendar_event_id, sep = "_")) %>%
    select(-game_pk,-date,-game_date,-venue_name)



# Now join the weather data to game_logs on the unique game identifier.
weather_game_logs <- game_logs %>%
    left_join(clean_weather_data, by = "game_key")

# =======================================================
# 4. MERGE PLAY-BY-PLAY DATA WITH GAME LOGS
# =======================================================
# Assuming pbp_data has a game_pk column that matches game_logs$game_pk.
final_data <- pbp_data %>%
    distinct(game_pk, at_bat_number,pitch_number, inning, inning_topbot, .keep_all = TRUE) %>%
    mutate(game_join = paste(game_pk, game_date,sep = "_")) %>%
    left_join(weather_game_logs %>%
                  filter(!is.na(teams_away_score)) %>%
                  group_by(game_pk) %>%
                  slice(1) %>%
                  ungroup() %>%
                  mutate(game_join = paste(game_pk, game_date, sep = "_")) %>%
                  select(-game_date),by = "game_pk") %>%
    mutate(
        total_score_pre = home_score + away_score,
        total_score_post = post_home_score + post_away_score,
        run_outcome = total_score_post - total_score_pre
    )
# final_data now includes game metadata (teams, scores, venue, weather)

# =======================================================
# 5. PREPARE DATA FOR THE BAYESIAN MODEL
# =======================================================
# (a) Helper function to count run outcomes per game.
count_runs <- function(df, inning_indicator) {
    df %>%
        filter(inning_topbot == inning_indicator) %>%
        group_by(game_pk) %>%
        summarise(
            run_1 = sum(run_outcome == 1, na.rm = TRUE),
            run_2 = sum(run_outcome == 2, na.rm = TRUE),
            run_3 = sum(run_outcome == 3, na.rm = TRUE),
            run_4 = sum(run_outcome == 4, na.rm = TRUE),
            .groups = "drop"
        )
}
# (b) Extract per-game team info from final_data.
team_info <- final_data %>%
    arrange(game_pk, at_bat_number) %>%
    group_by(game_pk) %>%
    summarise(
        game_date = first(game_date),
        venue_name = first(venue_name),
        home_team = first(home_team),
        away_team = first(away_team),
        # For the away team, take the first pitcher in inning 1 when the inning is "Top"
        away_starting_pitcher = first(pitcher[inning == 1 & inning_topbot == "Top"]),
        # For the home team, take the first pitcher in inning 1 when the inning is "Bot"
        home_starting_pitcher = first(pitcher[inning == 1 & inning_topbot == "Bot"]),
        temperature = first(temperature),
        wind_speed = first(wind_speed),
        humidity = first(humidity),
        home_score = max(teams_home_score),
        away_score = max(teams_away_score),
        .groups = "drop"
    )


# (c) Count run outcomes for each game.
home_runs_list <- list()
away_runs_list <- list()
for(game in unique(final_data$game_pk)){
    home_team <- team_info$home_team[team_info$game_pk == game]
    away_team <- team_info$away_team[team_info$game_pk == game]

    home_runs <- count_runs(final_data %>% filter(game_pk == game), "Bot")
    away_runs <- count_runs(final_data %>% filter(game_pk == game), "Top")

    home_runs_list[[game]] <- home_runs
    away_runs_list[[game]] <- away_runs
}
home_runs <- bind_rows(home_runs_list)
away_runs <- bind_rows(away_runs_list)

# (d) Merge run counts into team_info.
final_data_prepped <- team_info %>%
    left_join(home_runs, by = "game_pk") %>%
    rename(home_run_1 = run_1, home_run_2 = run_2, home_run_3 = run_3, home_run_4 = run_4) %>%
    left_join(away_runs, by = "game_pk") %>%
    rename(away_run_1 = run_1, away_run_2 = run_2, away_run_3 = run_3, away_run_4 = run_4)

check <- final_data_prepped %>% mutate(home_check = home_run_1+home_run_2*2+home_run_3*3+home_run_4*4,
                                                                               away_check = away_run_1+away_run_2*2+away_run_3*3+away_run_4*4,
                                                                               home_tf = home_check != home_score,
                                                                             away_tf = away_check != away_score)

# (e) Create indices for teams, pitchers, and stadiums.
teams <- unique(c(final_data_prepped$home_team, final_data_prepped$away_team))
team_ids <- setNames(seq_along(teams), teams)

pitchers <- unique(c(final_data_prepped$home_starting_pitcher, final_data_prepped$away_starting_pitcher))
pitcher_ids <- setNames(seq_along(pitchers), pitchers)

stadiums <- unique(final_data_prepped$venue_name)
stadium_ids <- setNames(seq_along(stadiums), stadiums)

final_data_prepped <- final_data_prepped %>%
    mutate(
        home_team_id = team_ids[home_team],
        away_team_id = team_ids[away_team],
        home_pitcher_id = pitcher_ids[home_starting_pitcher],
        away_pitcher_id = pitcher_ids[away_starting_pitcher],
        stadium_id = stadium_ids[venue_name]
    )

# (f) Replace any missing run counts with zeros.
final_data_prepped <- final_data_prepped %>%
    mutate_at(vars(home_run_1, home_run_2, home_run_3, home_run_4,
                   away_run_1, away_run_2, away_run_3, away_run_4),
              ~replace_na(., 0))

# (g) Prepare the data list for Stan.
stan_data <- list(
    N = nrow(final_data_prepped),  # number of games
    T = length(teams),             # number of teams
    P = length(pitchers),          # number of pitchers
    S = length(stadiums),          # number of stadiums
    home_team = final_data_prepped$home_team_id,
    away_team = final_data_prepped$away_team_id,
    home_pitcher = final_data_prepped$home_pitcher_id,
    away_pitcher = final_data_prepped$away_pitcher_id,
    stadium = final_data_prepped$stadium_id,
    home_run_1 = final_data_prepped$home_run_1,
    away_run_1 = final_data_prepped$away_run_1,
    home_run_2 = final_data_prepped$home_run_2,
    away_run_2 = final_data_prepped$away_run_2,
    home_run_3 = final_data_prepped$home_run_3,
    away_run_3 = final_data_prepped$away_run_3,
    home_run_4 = final_data_prepped$home_run_4,
    away_run_4 = final_data_prepped$away_run_4,
    temperature = final_data_prepped$temperature,
    wind_speed = final_data_prepped$wind_speed,
    humidity = final_data_prepped$humidity
)

# Optionally, save the prepared Stan data:
saveRDS(stan_data, "stan_input_data.RDS")

# =======================================================
# 6. DEFINE THE VECTORIZED STAN MODEL (ORIGINAL FRAMEWORK + WEATHER & STADIUM)
# =======================================================
stan_model_code <- "
data {
    int<lower=1> N;              // Number of games
    int<lower=1> T;              // Number of teams
    int<lower=1> P;              // Number of pitchers
    int<lower=1> S;              // Number of stadiums
    int home_team[N];            // Home team indices
    int away_team[N];            // Away team indices
    int home_pitcher[N];         // Home starting pitcher indices
    int away_pitcher[N];         // Away starting pitcher indices
    int stadium[N];              // Stadium index for each game
    int home_run_1[N];           // 1-run events for home team
    int away_run_1[N];           // 1-run events for away team
    int home_run_2[N];           // 2-run events for home team
    int away_run_2[N];           // 2-run events for away team
    int home_run_3[N];           // 3-run events for home team
    int away_run_3[N];           // 3-run events for away team
    int home_run_4[N];           // 4-run events for home team
    int away_run_4[N];           // 4-run events for away team
    real temperature[N];         // Game-day temperature (°F or °C depending on API)
    real wind_speed[N];          // Game-day wind speed (mph)
    real humidity[N];            // Game-day relative humidity (%)
}
parameters {
    real<lower=0> theta_run_1;   // Dispersion for 1-run events
    real<lower=0> theta_run_2;   // Dispersion for 2-run events
    real<lower=0> home_advantage; // Global home advantage
    real<lower=0> int_run_1;     // 1-run intercept
    real<lower=0> int_run_2;     // 2-run intercept
    real<lower=0> int_run_3;     // 3-run intercept
    real<lower=0> int_run_4;     // 4-run intercept

    // Raw team ability parameters (offense and defense for each run type)
    vector[T] att_run_1_raw;
    vector[T] def_run_1_raw;
    vector[T] att_run_2_raw;
    vector[T] def_run_2_raw;
    vector[T] att_run_3_raw;
    vector[T] def_run_3_raw;
    vector[T] att_run_4_raw;
    vector[T] def_run_4_raw;
    // Raw starting pitcher ability
    vector[P] pitcher_ability_raw;

    // Weather effect coefficients
    real temp_coeff;
    real wind_coeff;
    real humidity_coeff;

    // Stadium influence effect (applied to both teams)
    vector[S] stadium_effect;
}
transformed parameters {
    vector[T] att_run_1 = att_run_1_raw - mean(att_run_1_raw);
    vector[T] def_run_1 = def_run_1_raw - mean(def_run_1_raw);
    vector[T] att_run_2 = att_run_2_raw - mean(att_run_2_raw);
    vector[T] def_run_2 = def_run_2_raw - mean(def_run_2_raw);
    vector[T] att_run_3 = att_run_3_raw - mean(att_run_3_raw);
    vector[T] def_run_3 = def_run_3_raw - mean(def_run_3_raw);
    vector[T] att_run_4 = att_run_4_raw - mean(att_run_4_raw);
    vector[T] def_run_4 = def_run_4_raw - mean(def_run_4_raw);
    vector[P] pitcher_ability = pitcher_ability_raw - mean(pitcher_ability_raw);
}
model {
    // Global priors
    home_advantage ~ normal(0, 0.1);
    int_run_1 ~ normal(1, 0.2);
    int_run_2 ~ normal(1, 0.2);
    int_run_3 ~ normal(0.5, 0.2);
    int_run_4 ~ normal(0.25, 0.2);
    theta_run_1 ~ gamma(30, 1);
    theta_run_2 ~ gamma(30, 1);

    // Priors for team and pitcher abilities
    att_run_1_raw ~ normal(0, 0.2);
    def_run_1_raw ~ normal(0, 0.2);
    att_run_2_raw ~ normal(0, 0.2);
    def_run_2_raw ~ normal(0, 0.2);
    att_run_3_raw ~ normal(0, 0.2);
    def_run_3_raw ~ normal(0, 0.2);
    att_run_4_raw ~ normal(0, 0.2);
    def_run_4_raw ~ normal(0, 0.2);
    pitcher_ability_raw ~ normal(0, 0.2);

    // Weather effects
    temp_coeff ~ normal(0, 0.1);
    wind_coeff ~ normal(0, 0.1);
    humidity_coeff ~ normal(0, 0.1);

    // Stadium effect
    stadium_effect ~ normal(0, 0.2);

    // Vectorized linear predictors for each run event:
    vector[N] log_mu_home_1 = att_run_1[home_team] + def_run_1[away_team] + home_advantage + int_run_1 +
                               pitcher_ability[away_pitcher] +
                               temp_coeff * temperature + wind_coeff * wind_speed + humidity_coeff * humidity +
                               stadium_effect[stadium];
    vector[N] log_mu_away_1 = att_run_1[away_team] + def_run_1[home_team] + int_run_1 +
                               pitcher_ability[home_pitcher] +
                               temp_coeff * temperature + wind_coeff * wind_speed + humidity_coeff * humidity +
                               stadium_effect[stadium];

    vector[N] log_mu_home_2 = att_run_2[home_team] + def_run_2[away_team] + home_advantage + int_run_2 +
                               pitcher_ability[away_pitcher] +
                               temp_coeff * temperature + wind_coeff * wind_speed + humidity_coeff * humidity +
                               stadium_effect[stadium];
    vector[N] log_mu_away_2 = att_run_2[away_team] + def_run_2[home_team] + int_run_2 +
                               pitcher_ability[home_pitcher] +
                               temp_coeff * temperature + wind_coeff * wind_speed + humidity_coeff * humidity +
                               stadium_effect[stadium];

    vector[N] log_mu_home_3 = att_run_3[home_team] + def_run_3[away_team] + home_advantage + int_run_3 +
                               pitcher_ability[away_pitcher] +
                               temp_coeff * temperature + wind_coeff * wind_speed + humidity_coeff * humidity +
                               stadium_effect[stadium];
    vector[N] log_mu_away_3 = att_run_3[away_team] + def_run_3[home_team] + int_run_3 +
                               pitcher_ability[home_pitcher] +
                               temp_coeff * temperature + wind_coeff * wind_speed + humidity_coeff * humidity +
                               stadium_effect[stadium];

    vector[N] log_mu_home_4 = att_run_4[home_team] + def_run_4[away_team] + home_advantage + int_run_4 +
                               pitcher_ability[away_pitcher] +
                               temp_coeff * temperature + wind_coeff * wind_speed + humidity_coeff * humidity +
                               stadium_effect[stadium];
    vector[N] log_mu_away_4 = att_run_4[away_team] + def_run_4[home_team] + int_run_4 +
                               pitcher_ability[home_pitcher] +
                               temp_coeff * temperature + wind_coeff * wind_speed + humidity_coeff * humidity +
                               stadium_effect[stadium];

    // Likelihoods (vectorized)
    home_run_1 ~ neg_binomial_2_log(log_mu_home_1, theta_run_1);
    away_run_1 ~ neg_binomial_2_log(log_mu_away_1, theta_run_1);
    home_run_2 ~ neg_binomial_2_log(log_mu_home_2, theta_run_2);
    away_run_2 ~ neg_binomial_2_log(log_mu_away_2, theta_run_2);
    home_run_3 ~ poisson_log(log_mu_home_3);
    away_run_3 ~ poisson_log(log_mu_away_3);
    home_run_4 ~ poisson_log(log_mu_home_4);
    away_run_4 ~ poisson_log(log_mu_away_4);
}
"

# =======================================================
# 7. RUN THE STAN MODEL
# =======================================================
# Use the prepared stan_data object (or load it via readRDS if saved)
fit <- stan(
    model_code = stan_model_code,
    data = stan_data,
    iter = 10000,
    warmup = 2000,
    chains = 4,
    cores = 6,
    seed = 1,
    init = "random",
    control = list(max_treedepth = 10)
)

# Print summary of the model fit
print(fit)

# Optionally, produce traceplots for key parameters
traceplot(fit, pars = c("int_run_1", "int_run_2", "int_run_3", "int_run_4"))
