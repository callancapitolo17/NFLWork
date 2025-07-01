library(tidyverse)


# Read all CSV files starting with "epl_"
all_data <- list.files(pattern = "^epl_.*\\.csv$") %>%
  set_names() %>%
  map_dfr(read.csv, .id = "season_file")
write.csv(all_data,"PremBettingHistory.csv")


library(dplyr)
library(stringr)
library(purrr)
library(jsonlite)
library(tidyr)

# assume your raw data is in a data.frame called `df`
df_clean <- all_data %>%
  # 1. Fix quotes so it's valid JSON
  mutate(
    x1x2_json = str_replace_all(X1x2_market, "'", "\"")
  ) %>%
  # 2. Parse each row’s JSON into a list of lists
  mutate(
    odds_list = map(x1x2_json, ~ fromJSON(.x, simplifyDataFrame = FALSE))
  ) %>%
  # 3. Drop the old text columns
  select(-X1x2_market, -x1x2_json) %>%
  # 4. Unnest into long form (one row per bookmaker per match)
  unnest_longer(odds_list) %>%
  unnest_wider(odds_list) %>%
  # 5. Rename and convert to numeric
  mutate(match_id = paste(match_date,home_team,away_team, sep ="-")) %>% 
  transmute(
    # keep whatever identifiers you need, e.g. match_date, teams…
    match_id,
    bookmaker = bookmaker_name,
    home_odds = as.numeric(`1`),
    draw_odds = as.numeric(X),
    away_odds = as.numeric(`2`)
  )
