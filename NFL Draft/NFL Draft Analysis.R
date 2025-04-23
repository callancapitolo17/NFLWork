library(tidyverse)
library(readr)
library(dplyr)
path <- "2025 NFL Mock Draft_Big Board Consensus Tracker - Mock Drafts 2025.csv"




# Read both header rows manually
headers_raw <- read_csv(path, n_max = 2, col_names = FALSE)
header1 <- headers_raw[1, ] %>% unlist() %>% as.character()
header2 <- headers_raw[2, ] %>% unlist() %>% as.character()

# Combine the two header rows
combined_headers <- ifelse(is.na(header2) | header2 == "",
                           header1,
                           paste0(header1, " | ", header2))

# Read full data starting after headers
df <- read_csv(path, skip = 2, col_names = FALSE)
# Fix: assign default names to any empty or NA columns
clean_headers <- ifelse(combined_headers == "" | is.na(combined_headers),
                        paste0("V", seq_along(combined_headers)),
                        combined_headers)
colnames(df) <- clean_headers


# Step 2: Filter only player rows (i.e., rows with player names)
df <- df %>% filter(!is.na(`Consensus ADP current date range: 3/3 - 4/23 | NAME`))

# Step 3: Clean up the columns you want to keep
mock_columns <- df %>% 
  select(matches("SHARP|ESPN|WALTER|NFL Brooks|MOCK|TANK|GRINDING|ATHLETIC|PFF|Kiper|Fantasy Law|Jeremiah|Sikkema|Allbright|Donahue|4FOR4|ITA|Freedman|UD Norris|DEN")) %>%
  mutate_all(~as.numeric(as.character(.)))

df_players <- df %>%
  select(Player = `Consensus ADP current date range: 3/3 - 4/23 | NAME`,
         POS = `Green = one of 13 players invited to draft | POS`,
         COLLEGE = `Italics = need to update column | COLLEGE`) %>%
  bind_cols(mock_columns)

compare_mock_position <- function(df_players, player_a, player_b) {
  require(dplyr)
  require(tidyr)
  
  # Filter for the two players
  mock_df <- df_players %>%
    filter(Player %in% c(player_a, player_b))
  
  # Get only numeric mock columns (exclude Player/POS/etc.)
  mock_cols <- mock_df %>%
    select(where(is.numeric)) %>%
    colnames()
  
  # Reshape long and compare
  result <- mock_df %>%
    select(Player, all_of(mock_cols)) %>%
    pivot_longer(-Player, names_to = "Source", values_to = "Pick") %>%
    filter(!is.na(Pick)) %>%
    pivot_wider(names_from = Player, values_from = Pick)
  
  # Ensure both columns exist
  if (!(player_a %in% names(result)) | !(player_b %in% names(result))) {
    stop("One or both players not found in mock draft columns.")
  }
  
  # Calculate comparison stats
  summary <- result %>%
    mutate(
      A_before_B = !!sym(player_a) < !!sym(player_b),
      Tie = !!sym(player_a) == !!sym(player_b)
    ) %>%
    summarise(
      Player_A = player_a,
      Player_B = player_b,
      Total_Comparisons = sum(!is.na(!!sym(player_a)) & !is.na(!!sym(player_b))),
      A_Before_B = sum(A_before_B, na.rm = TRUE),
      B_Before_A = sum(!A_before_B & !Tie, na.rm = TRUE),
      Ties = sum(Tie, na.rm = TRUE),
      Pct_A_Before_B = A_Before_B / Total_Comparisons
    )
  
  return(summary)
}



# Step 4: Optional — Reshape to long format
df_long <- df_players %>%
  pivot_longer(-c(Player, POS, COLLEGE), names_to = "MockSource", values_to = "MockPick") %>%
  filter(!is.na(MockPick))

# Step 5: Summary statistics for each player (mean, std dev)
df_summary <- df_long %>%
  group_by(Player, POS, COLLEGE) %>%
  summarise(avg_pick = mean(MockPick), sd_pick = sd(MockPick), .groups = "drop")

# Step 6: Join with current betting O/U (if you have another data frame of over/unders)
# Example (you’ll need to add this separately):
# df_summary <- df_summary %>% left_join(df_odds, by = "Player")

# View cleaned data
print(head(df_summary))
