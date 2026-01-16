# Kalshi NFL Coaching Odds Display
# Reads coaching odds from DuckDB and displays heatmap

library(tidyverse)
library(gt)
library(duckdb)
library(scales)

# Path to DuckDB file (created by Python script)
db_path <- "kalshi coaching/kalshi_coaching.duckdb"

# Function to load latest odds from DuckDB with change from previous snapshot
load_latest_odds <- function(db_path) {
  con <- dbConnect(duckdb(), db_path, read_only = TRUE)

  # Get the two most recent fetch times
  fetch_times <- dbGetQuery(con, "
    SELECT DISTINCT fetch_time
    FROM coaching_odds
    ORDER BY fetch_time DESC
    LIMIT 2
  ")

  # Get latest odds
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

  # Get previous odds if we have more than one snapshot
  if (nrow(fetch_times) >= 2) {
    prev_time <- fetch_times$fetch_time[2]
    previous_odds <- dbGetQuery(con, sprintf("
      SELECT
        team,
        candidate,
        last_price as prev_price
      FROM coaching_odds
      WHERE fetch_time = '%s'
    ", prev_time))

    # Join to calculate change
    latest_odds <- latest_odds %>%
      left_join(previous_odds, by = c("team", "candidate")) %>%
      mutate(price_change = last_price - coalesce(prev_price, last_price))
  } else {
    latest_odds$prev_price <- NA
    latest_odds$price_change <- 0
  }

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
    arrange(desc(last_price), desc(volume), .by_group = TRUE) %>%
    slice_head(n = 5) %>%
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
      subtitle = paste("Top 5 candidates per team | Updated:", Sys.time())
    ) %>%
    tab_spanner(
      label = "1st",
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
    tab_spanner(
      label = "4th",
      columns = ends_with("_4")
    ) %>%
    tab_spanner(
      label = "5th",
      columns = ends_with("_5")
    ) %>%
    cols_label(
      team = "Team",
      candidate_1 = "Name", last_price_1 = "%", volume_1 = "Vol",
      candidate_2 = "Name", last_price_2 = "%", volume_2 = "Vol",
      candidate_3 = "Name", last_price_3 = "%", volume_3 = "Vol",
      candidate_4 = "Name", last_price_4 = "%", volume_4 = "Vol",
      candidate_5 = "Name", last_price_5 = "%", volume_5 = "Vol"
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
      axis.text.x = element_text(angle = 75, hjust = 0, vjust = 0.5, size = 9),
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
      axis.text.x = element_text(angle = 0, size = 12, face = "bold"),
      axis.text.y = element_text(size = 12),
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

# Create the heatmap
heatmap_plot <- create_coaching_heatmap_v2(odds, min_prob = 8)

# Create the gt summary table
summary_table <- create_summary_table(odds)

# Check if running from command line (non-interactive) - save to files
if (!interactive()) {
  # Save heatmap as PNG
  ggsave("kalshi coaching/heatmap.png", heatmap_plot, width = 14, height = 10, dpi = 150)
  cat("Saved heatmap to kalshi coaching/heatmap.png\n")

  # Save gt table as HTML
  gtsave(summary_table, "kalshi coaching/table.html")
  cat("Saved table to kalshi coaching/table.html\n")
} else {
  # Interactive mode (RStudio) - display plots
  print(heatmap_plot)
  print(summary_table)
}
