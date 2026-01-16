# Kalshi NFL Coaching Odds Display
# Reads coaching odds from DuckDB and displays heatmap

library(tidyverse)
library(gt)
library(duckdb)
library(scales)

# Path to DuckDB file (created by Python script)
db_path <- "kalshi coaching/kalshi_coaching.duckdb"

# Function to load latest odds from DuckDB
load_latest_odds <- function(db_path) {
  con <- dbConnect(duckdb(), db_path, read_only = TRUE)

  # Get the most recent fetch for each team/candidate
  latest_odds <- dbGetQuery(con, "
    WITH latest_fetch AS (
      SELECT MAX(fetch_time) as max_time
      FROM coaching_odds
    )
    SELECT
      team,
      candidate,
      yes_bid,
      yes_ask,
      last_price,
      volume,
      liquidity,
      fetch_time
    FROM coaching_odds
    WHERE fetch_time = (SELECT max_time FROM latest_fetch)
    ORDER BY team, last_price DESC
  ")

  # Get snapshot count
  snapshot_info <- dbGetQuery(con, "
    SELECT COUNT(DISTINCT fetch_time) as snapshots,
           MIN(fetch_time) as first_fetch,
           MAX(fetch_time) as last_fetch
    FROM coaching_odds
  ")

  dbDisconnect(con, shutdown = TRUE)

  cat("Database has", snapshot_info$snapshots, "snapshots\n")
  cat("Latest data from:", as.character(snapshot_info$last_fetch), "\n\n")

  return(latest_odds)
}

# Function to create a gt table for a single team
create_team_table <- function(df, team_name) {
  df %>%
    filter(team == team_name) %>%
    arrange(desc(last_price)) %>%
    head(10) %>%
    select(candidate, last_price, yes_bid, yes_ask, volume) %>%
    gt() %>%
    tab_header(
      title = team_name,
      subtitle = "NFL Head Coach Candidates"
    ) %>%
    cols_label(
      candidate = "Candidate",
      last_price = "Probability",
      yes_bid = "Bid",
      yes_ask = "Ask",
      volume = "Volume"
    ) %>%
    fmt_percent(
      columns = c(last_price, yes_bid, yes_ask),
      decimals = 0,
      scale_values = FALSE
    ) %>%
    fmt_number(
      columns = volume,
      use_seps = TRUE,
      decimals = 0
    ) %>%
    data_color(
      columns = last_price,
      palette = c("white", "darkgreen"),
      domain = c(0, 100)
    ) %>%
    tab_options(
      table.font.size = 12,
      heading.title.font.size = 16,
      heading.subtitle.font.size = 12
    )
}

# Function to create a summary table of top candidates per team
create_summary_table <- function(df) {
  df %>%
    group_by(team) %>%
    slice_max(last_price, n = 3) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    select(team, rank, candidate, last_price, volume) %>%
    pivot_wider(
      names_from = rank,
      values_from = c(candidate, last_price, volume),
      names_glue = "{.value}_{rank}"
    ) %>%
    gt() %>%
    tab_header(
      title = "NFL Coaching Carousel - Kalshi Odds",
      subtitle = paste("Top 3 candidates per team | Updated:", Sys.time())
    ) %>%
    tab_spanner(
      label = "Favorite",
      columns = ends_with("_1")
    ) %>%
    tab_spanner(
      label = "2nd",
      columns = ends_with("_2")
    ) %>%
    tab_spanner(
      label = "3rd",
      columns = ends_with("_3")
    ) %>%
    cols_label(
      team = "Team",
      candidate_1 = "Name", last_price_1 = "%", volume_1 = "Vol",
      candidate_2 = "Name", last_price_2 = "%", volume_2 = "Vol",
      candidate_3 = "Name", last_price_3 = "%", volume_3 = "Vol"
    ) %>%
    fmt_number(
      columns = starts_with("volume"),
      use_seps = TRUE,
      decimals = 0
    ) %>%
    data_color(
      columns = starts_with("last_price"),
      palette = c("white", "darkgreen"),
      domain = c(0, 100)
    ) %>%
    tab_options(
      table.font.size = 11,
      heading.title.font.size = 16,
      heading.subtitle.font.size = 11,
      column_labels.font.size = 10
    )
}

# Heatmap visualization with ggplot2
# Teams as columns, coaches as rows, color = probability, text size = volume
create_coaching_heatmap <- function(df, min_prob = 5) {
  # Filter to candidates with at least min_prob% chance somewhere
  top_candidates <- df %>%
    group_by(candidate) %>%
    filter(max(last_price) >= min_prob) %>%
    ungroup()

  # Create the full grid (all candidate x team combinations)
  all_combos <- expand_grid(
    candidate = unique(top_candidates$candidate),
    team = unique(top_candidates$team)
  )

  # Join with actual data
  plot_data <- all_combos %>%
    left_join(
      top_candidates %>% select(candidate, team, last_price, volume),
      by = c("candidate", "team")
    )

  # Order candidates by their max probability across all teams
  candidate_order <- plot_data %>%
    group_by(candidate) %>%
    summarize(max_prob = max(last_price, na.rm = TRUE)) %>%
    arrange(desc(max_prob)) %>%
    pull(candidate)

  team_order <- sort(unique(plot_data$team))

  plot_data <- plot_data %>%
    mutate(
      candidate = factor(candidate, levels = rev(candidate_order)),
      team = factor(team, levels = team_order),
      vol_scaled = scales::rescale(sqrt(volume), to = c(2, 5), na.rm = TRUE)
    )

  p <- ggplot(plot_data, aes(x = team, y = candidate)) +
    geom_tile(aes(fill = last_price), color = "white", linewidth = 0.5) +
    geom_text(
      aes(label = ifelse(!is.na(last_price), paste0(last_price, "%"), ""),
          size = vol_scaled),
      color = "white",
      fontface = "bold"
    ) +
    scale_fill_gradient(
      low = "#1a472a",
      high = "#2ecc71",
      na.value = "grey90",
      name = "Probability",
      labels = function(x) paste0(x, "%")
    ) +
    scale_size_identity() +
    scale_x_discrete(position = "top") +
    labs(
      title = "NFL Head Coach Hiring Odds - Kalshi",
      subtitle = paste("Updated:", format(Sys.time(), "%B %d, %Y %I:%M %p")),
      x = NULL,
      y = NULL,
      caption = "Text size reflects trading volume | Source: Kalshi"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0.5, size = 9),
      axis.text.y = element_text(size = 9),
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(color = "grey50", size = 10),
      plot.caption = element_text(color = "grey50", size = 8)
    )

  return(p)
}

# Alternative version with viridis color scale
create_coaching_heatmap_v2 <- function(df, min_prob = 5) {
  top_candidates <- df %>%
    group_by(candidate) %>%
    filter(max(last_price) >= min_prob) %>%
    ungroup()

  candidate_order <- top_candidates %>%
    group_by(candidate) %>%
    summarize(max_prob = max(last_price, na.rm = TRUE)) %>%
    arrange(desc(max_prob)) %>%
    pull(candidate)

  team_order <- sort(unique(top_candidates$team))

  plot_data <- top_candidates %>%
    mutate(
      candidate = factor(candidate, levels = rev(candidate_order)),
      team = factor(team, levels = team_order),
      vol_scaled = scales::rescale(log1p(volume), to = c(2.5, 6), na.rm = TRUE),
      label = paste0(last_price)
    )

  p <- ggplot(plot_data, aes(x = team, y = candidate)) +
    geom_tile(aes(fill = last_price), color = "white", linewidth = 0.8) +
    geom_text(
      aes(label = label, size = vol_scaled),
      color = "white",
      fontface = "bold"
    ) +
    scale_fill_viridis_c(
      option = "D",
      direction = -1,
      na.value = "grey90",
      name = "Prob %"
    ) +
    scale_size_identity() +
    scale_x_discrete(position = "top") +
    labs(
      title = "NFL Coaching Carousel - Kalshi Odds",
      subtitle = paste("Text size = trading volume |", format(Sys.time(), "%b %d, %Y")),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0.5, size = 10, face = "bold"),
      axis.text.y = element_text(size = 9),
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(color = "grey40", size = 11)
    )

  return(p)
}

# =============================================================================
# RUN THIS SECTION - Load latest data from DuckDB and create plots
# =============================================================================

# Load the latest odds from DuckDB
odds <- load_latest_odds(db_path)
cat("Loaded", nrow(odds), "coaching odds records\n")
cat("Teams:", paste(unique(odds$team), collapse = ", "), "\n\n")

# Create and display the heatmap
heatmap_plot <- create_coaching_heatmap_v2(odds, min_prob = 8)
print(heatmap_plot)

# Uncomment below to also show the gt summary table:
# summary_table <- create_summary_table(odds)
# print(summary_table)
