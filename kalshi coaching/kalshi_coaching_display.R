# Kalshi NFL Coaching Odds Display
# Reads coaching odds from DuckDB and displays unified HTML dashboard

library(tidyverse)
library(duckdb)
library(scales)
library(reactable)
library(plotly)
library(htmltools)
library(htmlwidgets)

# Path to DuckDB file (created by Python script)
db_path <- "kalshi coaching/kalshi_coaching.duckdb"

# Function to load latest odds from DuckDB with change from previous snapshot
load_latest_odds <- function(db_path) {
  con <- dbConnect(duckdb(), db_path, read_only = TRUE)

  # Get the two most recent fetch times
  fetch_times <- dbGetQuery(con, "
    SELECT DISTINCT fetch_time
    FROM coaching_odds_v2
    ORDER BY fetch_time DESC
    LIMIT 2
  ")

  # Get latest odds
  latest_odds <- dbGetQuery(con, "
    WITH latest_fetch AS (
      SELECT MAX(fetch_time) as max_time
      FROM coaching_odds_v2
    )
    SELECT
      ticker,
      team,
      candidate,
      yes_bid,
      yes_ask,
      no_bid,
      no_ask,
      last_price,
      volume,
      liquidity,
      fetch_time
    FROM coaching_odds_v2
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
      FROM coaching_odds_v2
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
    FROM coaching_odds_v2
  ")

  dbDisconnect(con, shutdown = TRUE)

  cat("Database has", snapshot_info$snapshots, "snapshots\n")
  cat("Latest data from:", as.character(snapshot_info$last_fetch), "\n\n")

  return(latest_odds)
}

# Function to create an interactive summary table of top candidates per team
create_summary_table <- function(df) {
  # Prepare data: top 10 per team, sorted by last_price then liquidity
  table_data <- df %>%
    group_by(team) %>%
    arrange(desc(last_price), desc(liquidity), .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup() %>%
    mutate(
      yes_spread = paste0(yes_bid, "/", yes_ask)
    ) %>%
    select(team, candidate, last_price, yes_spread, liquidity, volume, price_change)

  # Color scale for probability
  prob_colors <- function(value) {
    if (is.na(value) || value == 0) return(NULL)
    # Scale from white to dark green
    green_intensity <- min(255, round(value * 2.55))
    sprintf("background: rgba(46, 204, 113, %.2f)", value / 100)
  }

  # Create reactable grouped by team
  reactable(
    table_data,
    groupBy = "team",
    searchable = TRUE,
    filterable = TRUE,
    defaultPageSize = 20,
    defaultExpanded = TRUE,
    columns = list(
      team = colDef(name = "Team", minWidth = 100),
      candidate = colDef(
        name = "Candidate",
        minWidth = 150,
        cell = function(value, index) {
          change <- table_data$price_change[index]
          change_text <- if (!is.na(change) && change != 0) {
            sprintf(" (%s%d)", ifelse(change > 0, "+", ""), change)
          } else ""
          color <- if (!is.na(change) && change > 0) "green"
                   else if (!is.na(change) && change < 0) "red"
                   else "inherit"
          div(style = list(color = color, fontWeight = if (change != 0) "bold" else "normal"),
              paste0(value, change_text))
        }
      ),
      last_price = colDef(
        name = "Prob",
        cell = function(value) paste0(value, "%"),
        minWidth = 70
      ),
      yes_spread = colDef(name = "Bid/Ask", minWidth = 80),
      liquidity = colDef(
        name = "Liquidity",
        format = colFormat(currency = "USD", separators = TRUE, digits = 0),
        minWidth = 100
      ),
      volume = colDef(
        name = "Contracts",
        format = colFormat(separators = TRUE),
        minWidth = 100
      ),
      price_change = colDef(show = FALSE)
    ),
    theme = reactableTheme(
      headerStyle = list(fontWeight = "bold", borderBottom = "2px solid #2ecc71"),
      groupHeaderStyle = list(fontWeight = "bold", backgroundColor = "#f0f0f0")
    )
  )
}

# Heatmap visualization with ggplot2
# Teams as columns, coaches as rows, color = probability, text size = liquidity
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
      top_candidates %>% select(candidate, team, last_price, liquidity),
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
      liq_scaled = scales::rescale(sqrt(liquidity), to = c(2, 5), na.rm = TRUE)
    )

  p <- ggplot(plot_data, aes(x = team, y = candidate)) +
    geom_tile(aes(fill = last_price), color = "white", linewidth = 0.5) +
    geom_text(
      aes(label = ifelse(!is.na(last_price), paste0(last_price, "%"), ""),
          size = liq_scaled),
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
      caption = "Text size reflects liquidity | Source: Kalshi"
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

# =============================================================================
# Portfolio Loading Functions
# =============================================================================

# Load positions from DuckDB (latest snapshot) with market descriptions
load_positions <- function(db_path) {
  con <- dbConnect(duckdb(), db_path, read_only = TRUE)
  on.exit(dbDisconnect(con, shutdown = TRUE))

  # Check if positions table exists
  tables <- dbGetQuery(con, "SELECT name FROM sqlite_master WHERE type='table'")
  if (!"positions" %in% tables$name) {
    return(data.frame())
  }

  # Check if market_info table exists
  has_market_info <- "market_info" %in% tables$name

  # Join with coaching_odds and market_info to get names
  # Get both title and subtitle for "title - subtitle" display format
  if (has_market_info) {
    positions <- dbGetQuery(con, "
      WITH latest_pos AS (
        SELECT MAX(fetch_time) as max_time FROM positions
      ),
      latest_odds AS (
        SELECT MAX(fetch_time) as max_time FROM coaching_odds_v2
      )
      SELECT
        p.ticker,
        COALESCE(m.subtitle, o.candidate, p.ticker) as market_name,
        COALESCE(o.team, '') as team,
        m.title as market_title,
        p.position,
        p.market_exposure,
        p.realized_pnl,
        p.total_traded,
        p.fees_paid,
        p.fetch_time
      FROM positions p
      LEFT JOIN (
        SELECT DISTINCT ticker, candidate, team
        FROM coaching_odds_v2
        WHERE fetch_time = (SELECT max_time FROM latest_odds)
      ) o ON p.ticker = o.ticker
      LEFT JOIN market_info m ON p.ticker = m.ticker
      WHERE p.fetch_time = (SELECT max_time FROM latest_pos)
      ORDER BY ABS(p.market_exposure) DESC
    ")
  } else {
    positions <- dbGetQuery(con, "
      WITH latest_pos AS (
        SELECT MAX(fetch_time) as max_time FROM positions
      ),
      latest_odds AS (
        SELECT MAX(fetch_time) as max_time FROM coaching_odds_v2
      )
      SELECT
        p.ticker,
        COALESCE(o.candidate, p.ticker) as market_name,
        COALESCE(o.team, '') as team,
        NULL as market_title,
        p.position,
        p.market_exposure,
        p.realized_pnl,
        p.total_traded,
        p.fees_paid,
        p.fetch_time
      FROM positions p
      LEFT JOIN (
        SELECT DISTINCT ticker, candidate, team
        FROM coaching_odds_v2
        WHERE fetch_time = (SELECT max_time FROM latest_odds)
      ) o ON p.ticker = o.ticker
      WHERE p.fetch_time = (SELECT max_time FROM latest_pos)
      ORDER BY ABS(p.market_exposure) DESC
    ")
  }

  return(positions)
}

# Load resting orders from DuckDB (latest snapshot) with market descriptions
load_resting_orders <- function(db_path) {
  con <- dbConnect(duckdb(), db_path, read_only = TRUE)
  on.exit(dbDisconnect(con, shutdown = TRUE))

  # Check if resting_orders table exists
  tables <- dbGetQuery(con, "SELECT name FROM sqlite_master WHERE type='table'")
  if (!"resting_orders" %in% tables$name) {
    return(data.frame())
  }

  has_market_info <- "market_info" %in% tables$name

  if (has_market_info) {
    orders <- dbGetQuery(con, "
      WITH latest_ord AS (
        SELECT MAX(fetch_time) as max_time FROM resting_orders
      ),
      latest_odds AS (
        SELECT MAX(fetch_time) as max_time FROM coaching_odds_v2
      )
      SELECT
        r.ticker,
        COALESCE(m.subtitle, o.candidate, r.ticker) as market_name,
        COALESCE(o.team, '') as team,
        m.title as market_title,
        r.side,
        r.yes_price,
        r.no_price,
        r.remaining_count,
        r.created_time,
        r.expiration_time,
        r.fetch_time
      FROM resting_orders r
      LEFT JOIN (
        SELECT DISTINCT ticker, candidate, team
        FROM coaching_odds_v2
        WHERE fetch_time = (SELECT max_time FROM latest_odds)
      ) o ON r.ticker = o.ticker
      LEFT JOIN market_info m ON r.ticker = m.ticker
      WHERE r.fetch_time = (SELECT max_time FROM latest_ord)
      ORDER BY r.ticker, r.side
    ")
  } else {
    orders <- dbGetQuery(con, "
      WITH latest_ord AS (
        SELECT MAX(fetch_time) as max_time FROM resting_orders
      ),
      latest_odds AS (
        SELECT MAX(fetch_time) as max_time FROM coaching_odds_v2
      )
      SELECT
        r.ticker,
        COALESCE(o.candidate, r.ticker) as market_name,
        COALESCE(o.team, '') as team,
        NULL as market_title,
        r.side,
        r.yes_price,
        r.no_price,
        r.remaining_count,
        r.created_time,
        r.expiration_time,
        r.fetch_time
      FROM resting_orders r
      LEFT JOIN (
        SELECT DISTINCT ticker, candidate, team
        FROM coaching_odds_v2
        WHERE fetch_time = (SELECT max_time FROM latest_odds)
      ) o ON r.ticker = o.ticker
      WHERE r.fetch_time = (SELECT max_time FROM latest_ord)
      ORDER BY r.ticker, r.side
    ")
  }

  return(orders)
}

# Load position changes (compare current vs previous snapshot) with market descriptions
load_changes <- function(db_path) {
  con <- dbConnect(duckdb(), db_path, read_only = TRUE)
  on.exit(dbDisconnect(con, shutdown = TRUE))

  # Check if positions table exists
  tables <- dbGetQuery(con, "SELECT name FROM sqlite_master WHERE type='table'")
  if (!"positions" %in% tables$name) {
    return(data.frame())
  }

  # Get the two most recent fetch times
  fetch_times <- dbGetQuery(con, "
    SELECT DISTINCT fetch_time
    FROM positions
    ORDER BY fetch_time DESC
    LIMIT 2
  ")

  if (nrow(fetch_times) < 2) {
    return(data.frame())
  }

  current_time <- fetch_times$fetch_time[1]
  prev_time <- fetch_times$fetch_time[2]

  has_market_info <- "market_info" %in% tables$name

  # Compare current vs previous positions with market names
  if (has_market_info) {
    changes <- dbGetQuery(con, sprintf("
      WITH current AS (
        SELECT ticker, position, market_exposure
        FROM positions
        WHERE fetch_time = '%s'
      ),
      previous AS (
        SELECT ticker, position, market_exposure
        FROM positions
        WHERE fetch_time = '%s'
      ),
      latest_odds AS (
        SELECT MAX(fetch_time) as max_time FROM coaching_odds_v2
      ),
      odds_lookup AS (
        SELECT DISTINCT ticker, candidate, team
        FROM coaching_odds_v2
        WHERE fetch_time = (SELECT max_time FROM latest_odds)
      )
      SELECT
        COALESCE(c.ticker, p.ticker) as ticker,
        COALESCE(m.subtitle, o.candidate, COALESCE(c.ticker, p.ticker)) as market_name,
        COALESCE(o.team, '') as team,
        m.title as market_title,
        COALESCE(c.position, 0) as current_pos,
        COALESCE(p.position, 0) as prev_pos,
        COALESCE(c.position, 0) - COALESCE(p.position, 0) as pos_change,
        COALESCE(c.market_exposure, 0) as current_exp,
        COALESCE(p.market_exposure, 0) as prev_exp,
        COALESCE(c.market_exposure, 0) - COALESCE(p.market_exposure, 0) as exp_change,
        CASE
          WHEN c.ticker IS NULL THEN 'CLOSED'
          WHEN p.ticker IS NULL THEN 'NEW'
          ELSE 'CHANGED'
        END as change_type
      FROM current c
      FULL OUTER JOIN previous p ON c.ticker = p.ticker
      LEFT JOIN odds_lookup o ON COALESCE(c.ticker, p.ticker) = o.ticker
      LEFT JOIN market_info m ON COALESCE(c.ticker, p.ticker) = m.ticker
      WHERE c.position != p.position
         OR c.market_exposure != p.market_exposure
         OR c.ticker IS NULL
         OR p.ticker IS NULL
      ORDER BY ABS(COALESCE(c.position, 0) - COALESCE(p.position, 0)) DESC
    ", current_time, prev_time))
  } else {
    changes <- dbGetQuery(con, sprintf("
      WITH current AS (
        SELECT ticker, position, market_exposure
        FROM positions
        WHERE fetch_time = '%s'
      ),
      previous AS (
        SELECT ticker, position, market_exposure
        FROM positions
        WHERE fetch_time = '%s'
      ),
      latest_odds AS (
        SELECT MAX(fetch_time) as max_time FROM coaching_odds_v2
      ),
      odds_lookup AS (
        SELECT DISTINCT ticker, candidate, team
        FROM coaching_odds_v2
        WHERE fetch_time = (SELECT max_time FROM latest_odds)
      )
      SELECT
        COALESCE(c.ticker, p.ticker) as ticker,
        COALESCE(o.candidate, COALESCE(c.ticker, p.ticker)) as market_name,
        COALESCE(o.team, '') as team,
        NULL as market_title,
        COALESCE(c.position, 0) as current_pos,
        COALESCE(p.position, 0) as prev_pos,
        COALESCE(c.position, 0) - COALESCE(p.position, 0) as pos_change,
        COALESCE(c.market_exposure, 0) as current_exp,
        COALESCE(p.market_exposure, 0) as prev_exp,
        COALESCE(c.market_exposure, 0) - COALESCE(p.market_exposure, 0) as exp_change,
        CASE
          WHEN c.ticker IS NULL THEN 'CLOSED'
          WHEN p.ticker IS NULL THEN 'NEW'
          ELSE 'CHANGED'
        END as change_type
      FROM current c
      FULL OUTER JOIN previous p ON c.ticker = p.ticker
      LEFT JOIN odds_lookup o ON COALESCE(c.ticker, p.ticker) = o.ticker
      WHERE c.position != p.position
         OR c.market_exposure != p.market_exposure
         OR c.ticker IS NULL
         OR p.ticker IS NULL
      ORDER BY ABS(COALESCE(c.position, 0) - COALESCE(p.position, 0)) DESC
    ", current_time, prev_time))
  }

  return(changes)
}

# =============================================================================
# Portfolio Table Creation Functions
# =============================================================================

create_positions_table <- function(df) {
  if (nrow(df) == 0) {
    return(NULL)
  }

  table_data <- df %>%
    mutate(
      side = ifelse(position > 0, "YES", "NO"),
      contracts = abs(position),
      # Use title - subtitle format for full context
      display_name = ifelse(
        !is.na(market_title) & nzchar(market_title) & !is.na(market_name) & nzchar(market_name),
        paste0(market_title, " - ", market_name),
        ifelse(
          nzchar(team) & !is.na(team),
          paste0(market_name, " (", team, ")"),
          market_name
        )
      )
    ) %>%
    select(display_name, contracts, side, market_exposure, realized_pnl)

  reactable(
    table_data,
    searchable = TRUE,
    defaultSorted = list(market_exposure = "desc"),
    columns = list(
      display_name = colDef(name = "Market", minWidth = 200),
      contracts = colDef(name = "Contracts", minWidth = 80),
      side = colDef(name = "Side", minWidth = 60),
      market_exposure = colDef(
        name = "Exposure",
        format = colFormat(currency = "USD", digits = 2),
        minWidth = 100
      ),
      realized_pnl = colDef(
        name = "P&L",
        format = colFormat(currency = "USD", digits = 2),
        style = function(value) {
          if (is.na(value)) return(NULL)
          if (value > 0) list(color = "green", fontWeight = "bold")
          else if (value < 0) list(color = "red", fontWeight = "bold")
          else NULL
        },
        minWidth = 100
      )
    ),
    theme = reactableTheme(
      headerStyle = list(fontWeight = "bold", borderBottom = "2px solid #2ecc71")
    )
  )
}

create_orders_table <- function(df) {
  if (nrow(df) == 0) {
    return(NULL)
  }

  table_data <- df %>%
    mutate(
      price = ifelse(side == "yes", yes_price, no_price),
      side = toupper(side),
      display_name = ifelse(
        !is.na(market_title) & nzchar(market_title) & !is.na(market_name) & nzchar(market_name),
        paste0(market_title, " - ", market_name),
        ifelse(
          nzchar(team) & !is.na(team),
          paste0(market_name, " (", team, ")"),
          market_name
        )
      ),
      # Format expiration time
      expires = ifelse(
        !is.na(expiration_time),
        format(as.POSIXct(expiration_time), "%b %d, %I:%M %p"),
        "No expiry"
      )
    ) %>%
    select(display_name, side, price, remaining_count, expires)

  reactable(
    table_data,
    searchable = TRUE,
    columns = list(
      display_name = colDef(name = "Market", minWidth = 200),
      side = colDef(name = "Side", minWidth = 60),
      price = colDef(name = "Price (¢)", minWidth = 80),
      remaining_count = colDef(name = "Contracts", minWidth = 90),
      expires = colDef(name = "Expires", minWidth = 120)
    ),
    theme = reactableTheme(
      headerStyle = list(fontWeight = "bold", borderBottom = "2px solid #2ecc71")
    )
  )
}

create_changes_table <- function(df) {
  if (nrow(df) == 0) {
    return(NULL)
  }

  table_data <- df %>%
    mutate(
      # Determine action based on change type
      action = case_when(
        change_type == "NEW" ~ "Bought",
        change_type == "CLOSED" ~ "Sold",
        change_type == "CHANGED" & abs(current_pos) > abs(prev_pos) ~ "Added",
        change_type == "CHANGED" & abs(current_pos) < abs(prev_pos) ~ "Reduced",
        TRUE ~ "Changed"
      ),
      # Determine side: positive position = YES, negative = NO
      # For changes, look at which direction we moved
      side = case_when(
        change_type == "CLOSED" ~ ifelse(prev_pos > 0, "YES", "NO"),
        pos_change > 0 ~ "YES",
        pos_change < 0 ~ "NO",
        TRUE ~ "—"
      ),
      contracts = abs(pos_change),
      # Calculate average price per contract (exposure change / contracts)
      avg_price = ifelse(
        contracts > 0,
        round(abs(exp_change) / contracts * 100),  # Convert to cents
        NA
      ),
      price_fmt = ifelse(!is.na(avg_price), paste0(avg_price, "¢"), "—"),
      # Format exposure
      exposure_fmt = sprintf("$%.2f", abs(exp_change)),
      # Format before position
      before_fmt = case_when(
        prev_pos == 0 ~ "—",
        prev_pos > 0 ~ paste(abs(prev_pos), "YES"),
        prev_pos < 0 ~ paste(abs(prev_pos), "NO")
      ),
      # Format after position
      after_fmt = case_when(
        current_pos == 0 ~ "—",
        current_pos > 0 ~ paste(abs(current_pos), "YES"),
        current_pos < 0 ~ paste(abs(current_pos), "NO")
      ),
      # Market display name
      display_name = ifelse(
        !is.na(market_title) & nzchar(market_title) & !is.na(market_name) & nzchar(market_name),
        paste0(market_title, " - ", market_name),
        ifelse(
          nzchar(team) & !is.na(team),
          paste0(market_name, " (", team, ")"),
          market_name
        )
      )
    ) %>%
    select(display_name, action, side, contracts, price_fmt, exposure_fmt, before_fmt, after_fmt)

  reactable(
    table_data,
    columns = list(
      display_name = colDef(name = "Market", minWidth = 250),
      action = colDef(
        name = "Action",
        minWidth = 70,
        style = function(value) {
          if (value %in% c("Bought", "Added")) list(color = "green", fontWeight = "bold")
          else if (value %in% c("Sold", "Reduced")) list(color = "red", fontWeight = "bold")
          else NULL
        }
      ),
      side = colDef(name = "Side", minWidth = 50),
      contracts = colDef(name = "Contracts", minWidth = 90),
      price_fmt = colDef(name = "Price", minWidth = 60),
      exposure_fmt = colDef(name = "Spent", minWidth = 70),
      before_fmt = colDef(name = "Before", minWidth = 80),
      after_fmt = colDef(name = "After", minWidth = 80)
    ),
    theme = reactableTheme(
      headerStyle = list(fontWeight = "bold", borderBottom = "2px solid #2ecc71")
    )
  )
}

# =============================================================================
# Unified HTML Report
# =============================================================================

create_unified_report <- function(odds_table, heatmap_widget, positions_table,
                                   orders_table, changes_table, output_path) {
  timestamp <- format(Sys.time(), "%B %d, %Y %I:%M %p")

  # Calculate totals for positions header
  total_exposure <- 0
  total_realized <- 0

  # Build page structure using htmltools
  page <- tagList(
    tags$head(
      tags$meta(charset = "UTF-8"),
      tags$title("NFL Coaching Dashboard"),
      tags$style(HTML("
        * {
          box-sizing: border-box;
        }
        body {
          font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;
          margin: 0;
          padding: 20px;
          padding-left: 200px;
          background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
          min-height: 100vh;
          color: #e0e0e0;
        }

        /* Sidebar */
        .sidebar {
          position: fixed;
          top: 0;
          left: 0;
          width: 180px;
          height: 100vh;
          background: rgba(0, 0, 0, 0.4);
          backdrop-filter: blur(10px);
          -webkit-backdrop-filter: blur(10px);
          border-right: 1px solid rgba(255, 255, 255, 0.1);
          padding: 20px 0;
          z-index: 1000;
          overflow-y: auto;
        }
        .sidebar-title {
          color: #00cec9;
          font-size: 0.75em;
          font-weight: 600;
          text-transform: uppercase;
          letter-spacing: 1px;
          padding: 0 20px;
          margin-bottom: 16px;
        }
        .sidebar a {
          display: flex;
          align-items: center;
          gap: 10px;
          color: #a0a0a0;
          text-decoration: none;
          padding: 12px 20px;
          font-size: 0.9em;
          transition: all 0.2s ease;
          border-left: 3px solid transparent;
        }
        .sidebar a:hover {
          color: #ffffff;
          background: rgba(0, 206, 201, 0.1);
          border-left-color: #00cec9;
        }
        .sidebar a .icon {
          font-size: 1.1em;
        }
        .main-content {
          max-width: 1400px;
        }

        /* Scroll behavior */
        html {
          scroll-behavior: smooth;
        }

        /* Refresh button */
        .refresh-btn {
          background: rgba(255, 255, 255, 0.2);
          border: 1px solid rgba(255, 255, 255, 0.3);
          color: white;
          padding: 10px 20px;
          border-radius: 8px;
          font-size: 0.95em;
          font-weight: 600;
          cursor: pointer;
          transition: all 0.2s ease;
          display: inline-flex;
          align-items: center;
          gap: 8px;
        }
        .refresh-btn:hover {
          background: rgba(255, 255, 255, 0.3);
          transform: translateY(-1px);
        }
        .refresh-btn:active {
          transform: translateY(0);
        }
        .refresh-btn.loading {
          pointer-events: none;
          opacity: 0.7;
        }
        .refresh-btn .spinner {
          display: none;
          width: 16px;
          height: 16px;
          border: 2px solid rgba(255,255,255,0.3);
          border-top-color: white;
          border-radius: 50%;
          animation: spin 0.8s linear infinite;
        }
        .refresh-btn.loading .spinner {
          display: inline-block;
        }
        .refresh-btn.loading .icon {
          display: none;
        }
        @keyframes spin {
          to { transform: rotate(360deg); }
        }
        .header-content {
          display: flex;
          justify-content: space-between;
          align-items: center;
        }
        .header-left h1 {
          margin: 0;
          font-size: 2em;
          font-weight: 700;
          text-shadow: 0 2px 10px rgba(0,0,0,0.2);
        }
        .header-left p {
          margin: 8px 0 0 0;
          opacity: 0.9;
          font-size: 1.1em;
        }
        .toast {
          position: fixed;
          bottom: 30px;
          right: 30px;
          padding: 16px 24px;
          border-radius: 10px;
          color: white;
          font-weight: 500;
          z-index: 2000;
          animation: slideIn 0.3s ease;
          box-shadow: 0 8px 30px rgba(0,0,0,0.3);
        }
        .toast.success {
          background: linear-gradient(135deg, #00b894, #00cec9);
        }
        .toast.error {
          background: linear-gradient(135deg, #e74c3c, #c0392b);
        }
        @keyframes slideIn {
          from { transform: translateX(100%); opacity: 0; }
          to { transform: translateX(0); opacity: 1; }
        }
        .header {
          background: linear-gradient(135deg, #00b894 0%, #00cec9 50%, #0984e3 100%);
          color: white;
          padding: 30px 40px;
          border-radius: 16px;
          margin-bottom: 24px;
          box-shadow: 0 10px 40px rgba(0, 184, 148, 0.3);
          position: relative;
          overflow: hidden;
        }
        .header::before {
          content: '';
          position: absolute;
          top: -50%;
          right: -50%;
          width: 100%;
          height: 200%;
          background: radial-gradient(circle, rgba(255,255,255,0.1) 0%, transparent 70%);
          pointer-events: none;
        }
        .header h1 {
          margin: 0;
          font-size: 2em;
          font-weight: 700;
          text-shadow: 0 2px 10px rgba(0,0,0,0.2);
        }
        .header p {
          margin: 8px 0 0 0;
          opacity: 0.9;
          font-size: 1.1em;
        }
        .section {
          background: rgba(255, 255, 255, 0.05);
          backdrop-filter: blur(10px);
          -webkit-backdrop-filter: blur(10px);
          padding: 24px;
          border-radius: 16px;
          margin-bottom: 20px;
          box-shadow: 0 8px 32px rgba(0, 0, 0, 0.3);
          border: 1px solid rgba(255, 255, 255, 0.1);
          transition: transform 0.2s ease, box-shadow 0.2s ease;
        }
        .section:hover {
          transform: translateY(-2px);
          box-shadow: 0 12px 40px rgba(0, 0, 0, 0.4);
        }
        .section h2 {
          margin-top: 0;
          color: #00cec9;
          border-bottom: 2px solid rgba(0, 206, 201, 0.3);
          padding-bottom: 12px;
          font-size: 1.4em;
          font-weight: 600;
        }
        .section h3, .section p.empty {
          color: #888;
          font-style: italic;
        }
        .plotly { width: 100% !important; }

        /* Reactable table styling - HIGH CONTRAST */
        .reactable {
          background: transparent !important;
        }
        .rt-table {
          background: #1e2a38 !important;
          border-radius: 12px !important;
          overflow: hidden;
          border: 1px solid rgba(255, 255, 255, 0.1);
        }
        .rt-thead {
          background: #0d1b2a !important;
        }
        .rt-th {
          color: #00cec9 !important;
          font-weight: 700 !important;
          font-size: 0.85em !important;
          text-transform: uppercase !important;
          letter-spacing: 0.5px !important;
          border-bottom: 2px solid #00cec9 !important;
          padding: 14px 12px !important;
          white-space: nowrap !important;
        }
        .rt-tbody .rt-tr {
          transition: background 0.15s ease;
          background: #1e2a38 !important;
        }
        .rt-tbody .rt-tr:nth-child(even) {
          background: #243447 !important;
        }
        .rt-tbody .rt-tr:hover {
          background: #2d4a5e !important;
        }
        .rt-td {
          border-bottom: 1px solid rgba(255, 255, 255, 0.08) !important;
          padding: 12px !important;
          color: #ffffff !important;
          font-size: 0.95em !important;
        }
        .rt-search {
          background: #0d1b2a !important;
          border: 1px solid rgba(0, 206, 201, 0.3) !important;
          border-radius: 8px !important;
          color: #ffffff !important;
          padding: 10px 14px !important;
          font-size: 0.95em !important;
        }
        .rt-search:focus {
          border-color: #00cec9 !important;
          outline: none !important;
          box-shadow: 0 0 0 2px rgba(0, 206, 201, 0.2) !important;
        }
        .rt-search::placeholder {
          color: #6b7c8d !important;
        }
        .rt-pagination {
          color: #ffffff !important;
          padding: 12px !important;
          background: #0d1b2a !important;
          border-radius: 0 0 12px 12px !important;
        }
        .rt-pagination button {
          background: #00cec9 !important;
          border: none !important;
          color: #0d1b2a !important;
          font-weight: 600 !important;
          border-radius: 6px !important;
          padding: 8px 14px !important;
          cursor: pointer;
          transition: all 0.15s ease;
        }
        .rt-pagination button:hover {
          background: #00e6df !important;
          transform: translateY(-1px);
        }
        .rt-pagination button:disabled {
          background: #2d4a5e !important;
          color: #6b7c8d !important;
          cursor: not-allowed;
          transform: none;
        }
        .rt-page-info {
          color: #a0aec0 !important;
        }
        .rt-page-size-select {
          background: #1e2a38 !important;
          color: #ffffff !important;
          border: 1px solid rgba(0, 206, 201, 0.3) !important;
          border-radius: 4px !important;
          padding: 4px 8px !important;
        }

        /* Group headers */
        .rt-tr-group-header {
          background: #0d1b2a !important;
          color: #00cec9 !important;
          font-weight: 600 !important;
        }

        /* Filter inputs */
        .rt-filter {
          background: #0d1b2a !important;
          border: 1px solid rgba(0, 206, 201, 0.3) !important;
          border-radius: 4px !important;
          color: #ffffff !important;
          padding: 6px 8px !important;
        }
        .rt-filter:focus {
          border-color: #00cec9 !important;
          outline: none !important;
        }

        /* Make sure inline styles don't break contrast */
        .rt-td div, .rt-td span {
          color: inherit !important;
        }
      "))
    ),
    # Sidebar navigation
    tags$nav(class = "sidebar",
      div(class = "sidebar-title", "Navigation"),
      a(href = "#changes", span(class = "icon", HTML("&#x1F4CA;")), "Changes"),
      a(href = "#heatmap", span(class = "icon", HTML("&#x1F5FA;")), "Heatmap"),
      a(href = "#odds", span(class = "icon", HTML("&#x1F4B0;")), "Team Odds"),
      a(href = "#positions", span(class = "icon", HTML("&#x1F4BC;")), "Positions"),
      a(href = "#orders", span(class = "icon", HTML("&#x23F3;")), "Orders")
    ),
    # Main content
    div(class = "main-content",
      div(class = "header",
        div(class = "header-content",
          div(class = "header-left",
            h1("NFL Coaching Dashboard"),
            p(paste("Updated:", timestamp))
          ),
          tags$button(
            class = "refresh-btn",
            onclick = "refreshData()",
            span(class = "icon", HTML("&#x1F504;")),
            span(class = "spinner"),
            "Refresh Data"
          )
        )
      ),
      # JavaScript for refresh
      tags$script(HTML("
        function refreshData() {
          const btn = document.querySelector('.refresh-btn');
          btn.classList.add('loading');

          fetch('/refresh', { method: 'POST' })
            .then(response => response.json())
            .then(data => {
              btn.classList.remove('loading');
              if (data.success) {
                showToast('Data refreshed! Reloading...', 'success');
                setTimeout(() => location.reload(), 1000);
              } else {
                showToast('Error: ' + data.error, 'error');
              }
            })
            .catch(err => {
              btn.classList.remove('loading');
              showToast('Failed to refresh. Is the server running?', 'error');
            });
        }

        function showToast(message, type) {
          const toast = document.createElement('div');
          toast.className = 'toast ' + type;
          toast.textContent = message;
          document.body.appendChild(toast);
          setTimeout(() => toast.remove(), 4000);
        }
      ")),
      div(class = "section", id = "changes",
        h2("Changes Since Last Run"),
        if (!is.null(changes_table)) changes_table else p(class = "empty", "No changes detected")
      ),
      div(class = "section", id = "heatmap",
        h2("Market Heatmap"),
        p("Hover for details. Drag to zoom, double-click to reset."),
        if (!is.null(heatmap_widget)) heatmap_widget else p("No heatmap data")
      ),
      div(class = "section", id = "odds",
        h2("Detailed Odds by Team"),
        if (!is.null(odds_table)) odds_table else p("No odds data")
      ),
      div(class = "section", id = "positions",
        h2("Open Positions"),
        if (!is.null(positions_table)) positions_table else p(class = "empty", "No open positions")
      ),
      div(class = "section", id = "orders",
        h2("Resting Orders"),
        if (!is.null(orders_table)) orders_table else p(class = "empty", "No resting orders")
      )
    )
  )

  # Save as self-contained HTML
  save_html(page, file = output_path)
  cat("Saved interactive report to", output_path, "\n")
}

# Alternative version with viridis color scale - shows bid/ask spread
# Now includes tooltip text for ggplotly interactivity
create_coaching_heatmap_v2 <- function(df, positions_df = NULL, min_prob = 5) {
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

  # Join with positions if provided
  if (!is.null(positions_df) && nrow(positions_df) > 0) {
    positions_lookup <- positions_df %>%
      filter(position != 0) %>%  # Only include non-zero positions
      select(ticker, position, market_exposure) %>%
      mutate(
        has_position = TRUE,
        pos_side = ifelse(position > 0, "YES", "NO"),
        pos_contracts = abs(position),
        # Calculate average price: exposure / contracts * 100 for cents
        pos_price = round(abs(market_exposure) / abs(position) * 100)
      )

    plot_data <- top_candidates %>%
      left_join(positions_lookup, by = "ticker") %>%
      mutate(has_position = coalesce(has_position, FALSE))
  } else {
    plot_data <- top_candidates %>%
      mutate(has_position = FALSE, pos_side = NA_character_, pos_contracts = NA_integer_, pos_price = NA_integer_)
  }

  plot_data <- plot_data %>%
    mutate(
      candidate = factor(candidate, levels = rev(candidate_order)),
      team = factor(team, levels = team_order),
      liq_scaled = scales::rescale(log1p(liquidity), to = c(2.5, 5), na.rm = TRUE),
      # Show bid/ask spread with star if we have position
      label = ifelse(has_position,
                     paste0("\u2605 ", yes_bid, "/", yes_ask),  # ★ star
                     paste0(yes_bid, "/", yes_ask)),
      # Tooltip text for ggplotly
      tooltip_text = paste0(
        "<b>", candidate, "</b><br>",
        "Team: ", team, "<br>",
        "Probability: ", last_price, "%<br>",
        "Bid/Ask: ", yes_bid, "/", yes_ask, "<br>",
        "Liquidity: $", format(liquidity, big.mark = ",", scientific = FALSE), "<br>",
        "Volume: ", format(volume, big.mark = ",", scientific = FALSE),
        ifelse(has_position,
               paste0("<br><b>Your Position: ", pos_contracts, " ", pos_side, " @ ", pos_price, "¢</b>"),
               "")
      )
    )

  p <- ggplot(plot_data, aes(x = team, y = candidate, text = tooltip_text)) +
    geom_tile(aes(fill = last_price), color = "white", linewidth = 0.8) +
    geom_text(
      aes(label = label, size = liq_scaled),
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
      title = "NFL Coaching Carousel - Kalshi Odds (Yes Bid/Ask)",
      subtitle = paste("Text size = liquidity | \u2605 = your position |", format(Sys.time(), "%b %d, %Y")),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11, face = "bold"),
      axis.text.y = element_text(size = 11),
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(color = "grey40", size = 11),
      plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )

  return(p)
}

# =============================================================================
# RUN THIS SECTION - Load latest data from DuckDB and create unified report
# =============================================================================

# Load the latest odds from DuckDB
odds <- load_latest_odds(db_path)
cat("Loaded", nrow(odds), "coaching odds records\n")
cat("Teams:", paste(unique(odds$team), collapse = ", "), "\n\n")

# Load portfolio data first (needed for heatmap)
positions <- load_positions(db_path)
orders <- load_resting_orders(db_path)
changes <- load_changes(db_path)

cat("Portfolio: ", nrow(positions), " positions, ", nrow(orders), " orders, ",
    nrow(changes), " changes\n\n", sep = "")

# Create the heatmap (with positions overlay)
heatmap_plot <- create_coaching_heatmap_v2(odds, positions_df = positions, min_prob = 8)

# Create the reactable summary table for odds
summary_table <- create_summary_table(odds)

# Create portfolio tables
positions_table <- create_positions_table(positions)
orders_table <- create_orders_table(orders)
changes_table <- create_changes_table(changes)

# Check if running from command line (non-interactive) - save to files
if (!interactive()) {
  # Convert heatmap to interactive plotly widget
  interactive_heatmap <- ggplotly(heatmap_plot, tooltip = "text", height = 600) %>%
    layout(margin = list(t = 120, b = 50, l = 150))
  cat("Generated interactive heatmap\n")

  # Create unified report with interactive widgets
  create_unified_report(
    odds_table = summary_table,
    heatmap_widget = interactive_heatmap,
    positions_table = positions_table,
    orders_table = orders_table,
    changes_table = changes_table,
    output_path = "kalshi coaching/report.html"
  )
} else {
  # Interactive mode (RStudio) - display plots
  print(heatmap_plot)
  print(summary_table)
  if (!is.null(positions_table)) print(positions_table)
  if (!is.null(orders_table)) print(orders_table)
  if (!is.null(changes_table)) print(changes_table)
}
