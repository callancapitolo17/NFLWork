# NFL +EV Betting Dashboard
# Clean reactable-based design with proper styling

library(tidyverse)
library(duckdb)
library(reactable)
library(htmltools)
library(htmlwidgets)
library(jsonlite)
library(digest)
library(lubridate)

# =============================================================================
# CONFIGURATION
# =============================================================================

DASHBOARD_DIR <- normalizePath("~/NFLWork/Answer Keys/NFL Dashboard", mustWork = FALSE)
DB_PATH <- file.path(DASHBOARD_DIR, "nfl_dashboard.duckdb")
OUTPUT_PATH <- file.path(DASHBOARD_DIR, "report.html")

# =============================================================================
# DATA LOADING
# =============================================================================

load_placed_bets <- function(db_path) {
  if (!file.exists(db_path)) {
    return(tibble(bet_hash = character()))
  }
  con <- dbConnect(duckdb(), db_path, read_only = TRUE)
  on.exit(dbDisconnect(con, shutdown = TRUE))
  tryCatch({
    dbGetQuery(con, "SELECT * FROM placed_bets WHERE status = 'pending'")
  }, error = function(e) {
    tibble(bet_hash = character())
  })
}

generate_bet_hash <- function(game_id, market, bet_on, line) {
  line_str <- ifelse(is.na(line), "NA", as.character(line))
  digest(paste(game_id, market, bet_on, line_str, sep = "|"), algo = "sha256")
}

# =============================================================================
# CORRELATION DETECTION
# =============================================================================

load_market_relationships <- function() {
  json_path <- file.path(DASHBOARD_DIR, "market_relationships.json")
  if (file.exists(json_path)) {
    fromJSON(json_path)
  } else {
    list(market_groups = list(), cross_group_correlations = list())
  }
}

get_market_group <- function(market, relationships) {
  for (group_name in names(relationships$market_groups)) {
    if (market %in% relationships$market_groups[[group_name]]$markets) {
      return(group_name)
    }
  }
  return(NA)
}

find_correlated_bets <- function(bet, placed_bets, relationships) {
  if (nrow(placed_bets) == 0) return(list(has_correlation = FALSE, details = NULL))

  # Same game only

  same_game <- placed_bets %>% filter(game_id == bet$id)
  if (nrow(same_game) == 0) return(list(has_correlation = FALSE, details = NULL))

  bet_group <- get_market_group(bet$market, relationships)
  correlations <- list()

  for (i in seq_len(nrow(same_game))) {
    placed <- same_game[i, ]
    placed_group <- get_market_group(placed$market, relationships)

    # Same market group = high correlation
    if (!is.na(bet_group) && !is.na(placed_group) && bet_group == placed_group) {
      correlations <- append(correlations, list(list(
        market = placed$market,
        bet_on = placed$bet_on,
        line = placed$line,
        strength = 0.90,
        level = "high"
      )))
      next
    }

    # Cross-group correlations (JSON loads as data frame)
    cross_df <- relationships$cross_group_correlations
    if (!is.null(cross_df) && nrow(cross_df) > 0) {
      for (j in seq_len(nrow(cross_df))) {
        g1 <- cross_df$group1[j]
        g2 <- cross_df$group2[j]
        strength <- cross_df$strength[j]
        if ((bet_group == g1 && placed_group == g2) ||
            (bet_group == g2 && placed_group == g1)) {
          level <- if (strength >= 0.80) "high" else if (strength >= 0.60) "medium" else "low"
          correlations <- append(correlations, list(list(
            market = placed$market,
            bet_on = placed$bet_on,
            line = placed$line,
            strength = strength,
            level = level
          )))
          break
        }
      }
    }
  }

  if (length(correlations) > 0) {
    max_strength <- max(sapply(correlations, function(x) x$strength))
    max_level <- if (max_strength >= 0.80) "high" else if (max_strength >= 0.60) "medium" else "low"
    return(list(has_correlation = TRUE, level = max_level, details = correlations))
  }

  return(list(has_correlation = FALSE, details = NULL))
}

# =============================================================================
# TABLE CREATION
# =============================================================================

create_placed_bets_table <- function(placed_bets) {
  if (nrow(placed_bets) == 0) return(NULL)

  table_data <- placed_bets %>%
    mutate(
      game = paste(away_team, "@", home_team),
      ev_display = ifelse(model_ev >= 0, sprintf("+%.1f%%", model_ev * 100), sprintf("%.1f%%", model_ev * 100)),
      odds_display = ifelse(odds > 0, paste0("+", odds), as.character(odds)),
      size_display = sprintf("$%.2f", coalesce(actual_size, recommended_size)),
      line_display = case_when(
        is.na(line) ~ "",
        line > 0 ~ paste0("+", line),
        TRUE ~ as.character(line)
      ),
      market_display = market %>%
        str_replace("alternate_", "Alt ") %>%
        str_replace("team_totals", "Team Tot") %>%
        str_replace("totals", "Total") %>%
        str_replace("spreads", "Spread") %>%
        str_replace("h2h_3_way", "ML-3") %>%
        str_replace("h2h", "ML") %>%
        str_replace("_q", " Q") %>%
        str_replace("_h", " H")
    ) %>%
    select(bet_hash, game, market_display, bet_on, line_display, ev_display, odds_display, size_display, bookmaker)

  reactable(
    table_data,
    compact = TRUE,
    columns = list(
      bet_hash = colDef(show = FALSE),
      game = colDef(name = "Game", minWidth = 200),
      market_display = colDef(name = "Market", minWidth = 90),
      bet_on = colDef(name = "Pick", minWidth = 120, cell = function(value, index) {
        line <- table_data$line_display[index]
        if (line != "") {
          div(span(style = "font-weight: 600;", value), span(style = "margin-left: 6px; color: #888;", line))
        } else {
          span(style = "font-weight: 600;", value)
        }
      }),
      line_display = colDef(show = FALSE),
      ev_display = colDef(name = "EV", minWidth = 70, align = "right", style = list(color = "#3fb950", fontWeight = "600")),
      odds_display = colDef(name = "Odds", minWidth = 60, align = "right"),
      size_display = colDef(name = "Size", minWidth = 60, align = "right"),
      bookmaker = colDef(name = "Book", minWidth = 90)
    ),
    theme = reactableTheme(
      color = "#c9d1d9",
      backgroundColor = "#21262d",
      borderColor = "#30363d",
      headerStyle = list(backgroundColor = "#161b22", color = "#8b949e", fontWeight = "600", fontSize = "0.7rem")
    )
  )
}

create_bets_table <- function(all_bets, placed_bets, relationships) {
  placed_hashes <- if (nrow(placed_bets) > 0) placed_bets$bet_hash else character()

  # Check correlations for each bet
  correlation_info <- lapply(seq_len(nrow(all_bets)), function(i) {
    find_correlated_bets(all_bets[i, ], placed_bets, relationships)
  })

  # Prepare table data
  table_data <- all_bets %>%
    mutate(
      bet_hash = pmap_chr(list(id, market, bet_on, line), generate_bet_hash),
      is_placed = bet_hash %in% placed_hashes,
      game = paste(away_team, "@", home_team),
      game_time = format(pt_start_time, "%a %m/%d %I:%M %p"),
      ev_pct = ev * 100,
      ev_display = ifelse(ev >= 0, sprintf("+%.1f%%", ev * 100), sprintf("%.1f%%", ev * 100)),
      line_display = case_when(
        is.na(line) ~ "-",
        line > 0 ~ paste0("+", line),
        TRUE ~ as.character(line)
      ),
      odds_display = ifelse(odds > 0, paste0("+", odds), as.character(odds)),
      size_display = sprintf("$%.2f", bet_size),
      # Correlation warnings
      has_correlation = sapply(correlation_info, function(x) x$has_correlation),
      correlation_level = sapply(correlation_info, function(x) if (x$has_correlation) x$level else "none"),
      correlation_tooltip = sapply(correlation_info, function(x) {
        if (!x$has_correlation) return("")
        details <- x$details
        lines <- sapply(details, function(d) {
          market_name <- d$market %>%
            str_replace("alternate_", "Alt ") %>%
            str_replace("team_totals", "Team Tot") %>%
            str_replace("totals", "Total") %>%
            str_replace("spreads", "Spread") %>%
            str_replace("h2h_3_way", "ML-3") %>%
            str_replace("h2h", "ML") %>%
            str_replace("_q", " Q") %>%
            str_replace("_h", " H")
          line_str <- if (!is.null(d$line) && !is.na(d$line)) {
            if (d$line > 0) paste0(" +", d$line) else paste0(" ", d$line)
          } else ""
          sprintf("%s - %s%s", market_name, d$bet_on, line_str)
        })
        paste("Correlated with:", paste(lines, collapse = ", "))
      }),
      # Simplify market names
      market_display = market %>%
        str_replace("alternate_", "Alt ") %>%
        str_replace("team_totals", "Team Tot") %>%
        str_replace("totals", "Total") %>%
        str_replace("spreads", "Spread") %>%
        str_replace("h2h_3_way", "ML-3") %>%
        str_replace("h2h", "ML") %>%
        str_replace("_q", " Q") %>%
        str_replace("_h", " H")
    ) %>%
    arrange(desc(ev)) %>%
    mutate(warning = "") %>%  # Placeholder for warning column
    select(
      bet_hash, id, warning, game, game_time, market, market_display, bet_on, line, line_display,
      ev_pct, ev_display, odds, odds_display, bet_size, size_display,
      bookmaker_key, is_placed, home_team, away_team, pt_start_time, prob,
      has_correlation, correlation_level, correlation_tooltip
    )

  reactable(
    table_data,
    searchable = TRUE,
    filterable = TRUE,
    defaultPageSize = 50,
    showPageSizeOptions = TRUE,
    pageSizeOptions = c(25, 50, 100),
    striped = TRUE,
    highlight = TRUE,
    compact = TRUE,
    columns = list(
      # Hidden columns for data
      bet_hash = colDef(show = FALSE),
      id = colDef(show = FALSE),
      home_team = colDef(show = FALSE),
      away_team = colDef(show = FALSE),
      pt_start_time = colDef(show = FALSE),
      market = colDef(show = FALSE),
      line = colDef(show = FALSE),
      ev_pct = colDef(show = FALSE),
      odds = colDef(show = FALSE),
      bet_size = colDef(show = FALSE),
      prob = colDef(show = FALSE),
      has_correlation = colDef(show = FALSE),
      correlation_level = colDef(show = FALSE),
      correlation_tooltip = colDef(show = FALSE),

      # Visible columns - Warning indicator first
      warning = colDef(
        name = "",
        minWidth = 40,
        maxWidth = 40,
        align = "center",
        html = TRUE,
        cell = function(value, index) {
          level <- table_data$correlation_level[index]
          tooltip <- table_data$correlation_tooltip[index]
          if (level == "high") {
            sprintf('<span class="warning-icon high" data-tooltip="%s">&#9888;</span>', htmltools::htmlEscape(tooltip))
          } else if (level == "medium") {
            sprintf('<span class="warning-icon medium" data-tooltip="%s">&#9888;</span>', htmltools::htmlEscape(tooltip))
          } else if (level == "low") {
            sprintf('<span class="warning-icon low" data-tooltip="%s">&#9675;</span>', htmltools::htmlEscape(tooltip))
          } else {
            ""
          }
        }
      ),
      game = colDef(
        name = "Game",
        minWidth = 220,
        filterable = TRUE
      ),
      game_time = colDef(
        name = "Time",
        minWidth = 150
      ),
      market_display = colDef(
        name = "Market",
        minWidth = 100,
        filterable = TRUE
      ),
      bet_on = colDef(
        name = "Pick",
        minWidth = 140,
        cell = function(value, index) {
          line <- table_data$line_display[index]
          if (line != "-") {
            div(
              span(style = "font-weight: 600;", value),
              span(style = "margin-left: 8px; color: #888; font-size: 0.9em;", line)
            )
          } else {
            span(style = "font-weight: 600;", value)
          }
        }
      ),
      line_display = colDef(show = FALSE),
      ev_display = colDef(
        name = "EV",
        minWidth = 80,
        align = "right",
        cell = function(value, index) {
          ev <- table_data$ev_pct[index]
          color <- if (ev >= 15) "#3fb950" else if (ev >= 10) "#56d364" else if (ev >= 5) "#7ee787" else "#a5d6a7"
          div(style = list(
            color = color,
            fontWeight = "600"
          ), value)
        }
      ),
      odds_display = colDef(
        name = "Odds",
        minWidth = 70,
        align = "right",
        style = list(fontFamily = "monospace")
      ),
      size_display = colDef(
        name = "Size",
        minWidth = 70,
        align = "right"
      ),
      bookmaker_key = colDef(
        name = "Book",
        minWidth = 100,
        filterable = TRUE
      ),
      is_placed = colDef(
        name = "Action",
        minWidth = 90,
        align = "center",
        filterable = FALSE,
        html = TRUE,
        cell = function(value, index) {
          row <- table_data[index, ]
          data_attrs <- sprintf(
            'data-hash="%s" data-game-id="%s" data-home="%s" data-away="%s" data-time="%s" data-market="%s" data-bet-on="%s" data-line="%s" data-prob="%s" data-ev="%s" data-size="%s" data-odds="%s" data-book="%s"',
            row$bet_hash, row$id, row$home_team, row$away_team,
            as.character(row$pt_start_time), row$market, row$bet_on,
            ifelse(is.na(row$line), "", row$line),
            row$prob, row$ev_pct / 100, row$bet_size, row$odds, row$bookmaker_key
          )

          if (value) {
            sprintf('<button class="btn-placed" onclick="removeBet(this)" %s>Placed</button>', data_attrs)
          } else {
            sprintf('<button class="btn-place" onclick="placeBet(this)" %s>Place</button>', data_attrs)
          }
        }
      )
    ),
    theme = reactableTheme(
      color = "#c9d1d9",
      backgroundColor = "#1c2128",
      borderColor = "#373e47",
      stripedColor = "#22272e",
      highlightColor = "#2d333b",
      headerStyle = list(
        backgroundColor = "#161b22",
        color = "#8b949e",
        fontWeight = "600",
        fontSize = "0.7rem",
        textTransform = "uppercase",
        letterSpacing = "0.5px",
        borderBottom = "1px solid #373e47"
      ),
      searchInputStyle = list(
        backgroundColor = "#161b22",
        color = "#c9d1d9",
        border = "1px solid #373e47",
        borderRadius = "6px",
        padding = "8px 12px"
      ),
      filterInputStyle = list(
        backgroundColor = "#161b22",
        color = "#c9d1d9",
        border = "1px solid #373e47"
      ),
      paginationStyle = list(
        backgroundColor = "#161b22",
        color = "#c9d1d9"
      ),
      pageButtonStyle = list(
        backgroundColor = "#2d333b",
        color = "#c9d1d9"
      ),
      pageButtonActiveStyle = list(
        backgroundColor = "#58a6ff",
        color = "#0d1117"
      )
    )
  )
}

# =============================================================================
# HTML REPORT
# =============================================================================

create_report <- function(bets_table, placed_table, stats, timestamp, games_json, filter_options_json) {
  page <- tagList(
    tags$head(
      tags$meta(charset = "UTF-8"),
      tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
      tags$title("NFL +EV Dashboard"),
      tags$link(
        href = "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap",
        rel = "stylesheet"
      ),
      tags$script(src = "https://cdn.jsdelivr.net/npm/chart.js"),
      # Inject games data for parlay builder
      tags$script(HTML(sprintf("window.GAMES = %s;", games_json))),
      # Inject filter options for bets table
      tags$script(HTML(sprintf("window.FILTER_OPTIONS = %s;", filter_options_json))),
      tags$style(HTML('
        * { box-sizing: border-box; }

        body {
          font-family: "Inter", -apple-system, BlinkMacSystemFont, sans-serif;
          background: #0d1117;
          color: #c9d1d9;
          margin: 0;
          padding: 32px;
          min-height: 100vh;
        }

        .container {
          max-width: 1600px;
          margin: 0 auto;
        }

        .header {
          display: flex;
          justify-content: space-between;
          align-items: center;
          margin-bottom: 28px;
          padding-bottom: 20px;
          border-bottom: 1px solid #21262d;
        }

        .header h1 {
          font-size: 1.4rem;
          font-weight: 600;
          color: #f0f6fc;
          margin: 0;
        }

        .header .subtitle {
          font-size: 0.8rem;
          color: #8b949e;
          margin-top: 4px;
        }

        .refresh-btn {
          background: #238636;
          border: 1px solid #2ea043;
          color: #fff;
          padding: 8px 20px;
          border-radius: 6px;
          font-size: 0.85rem;
          font-weight: 500;
          cursor: pointer;
          transition: all 0.15s;
        }

        .refresh-btn:hover {
          background: #2ea043;
        }

        .stats-row {
          display: flex;
          gap: 12px;
          margin-bottom: 24px;
        }

        .stat-card {
          background: #161b22;
          border: 1px solid #21262d;
          border-radius: 6px;
          padding: 14px 20px;
          min-width: 110px;
        }

        .stat-value {
          font-size: 1.5rem;
          font-weight: 600;
          color: #f0f6fc;
        }

        .stat-label {
          font-size: 0.65rem;
          color: #8b949e;
          text-transform: uppercase;
          letter-spacing: 0.5px;
          margin-top: 2px;
        }

        .table-container {
          background: #161b22;
          border: 1px solid #21262d;
          border-radius: 6px;
          padding: 16px;
        }

        /* Button styles */
        .btn-place {
          background: transparent;
          border: 1px solid #238636;
          color: #3fb950;
          padding: 5px 12px;
          border-radius: 6px;
          font-size: 0.75rem;
          font-weight: 500;
          cursor: pointer;
          transition: all 0.15s;
        }

        .btn-place:hover {
          background: #238636;
          color: #fff;
        }

        .btn-placed {
          background: #21262d;
          border: 1px solid #30363d;
          color: #8b949e;
          padding: 5px 12px;
          border-radius: 6px;
          font-size: 0.75rem;
          font-weight: 500;
          cursor: pointer;
          transition: all 0.15s;
        }

        .btn-placed:hover {
          border-color: #f85149;
          color: #f85149;
        }

        /* Toast */
        .toast {
          position: fixed;
          bottom: 24px;
          right: 24px;
          padding: 12px 16px;
          border-radius: 6px;
          font-weight: 500;
          font-size: 0.85rem;
          z-index: 1000;
          animation: slideIn 0.2s ease;
        }

        .toast.success { background: #238636; color: #fff; }
        .toast.error { background: #da3633; color: #fff; }

        @keyframes slideIn {
          from { transform: translateY(16px); opacity: 0; }
          to { transform: translateY(0); opacity: 1; }
        }

        /* Section headers */
        .section-header {
          font-size: 0.9rem;
          font-weight: 600;
          color: #8b949e;
          text-transform: uppercase;
          letter-spacing: 0.5px;
          margin: 24px 0 12px 0;
        }

        .placed-section {
          margin-bottom: 8px;
        }

        /* Warning icons with tooltips */
        .rt-td { overflow: visible !important; }

        .warning-icon {
          position: relative;
          cursor: help;
          font-size: 1.1em;
          display: inline-block;
        }
        .warning-icon.high { color: #f85149; }
        .warning-icon.medium { color: #d29922; }
        .warning-icon.low { color: #8b949e; font-size: 0.9em; }

        .warning-icon:hover::after {
          content: attr(data-tooltip);
          position: fixed;
          background: #161b22;
          color: #c9d1d9;
          padding: 8px 12px;
          border-radius: 6px;
          font-size: 0.75rem;
          white-space: nowrap;
          z-index: 9999;
          border: 1px solid #30363d;
          box-shadow: 0 4px 12px rgba(0,0,0,0.4);
          margin-top: -40px;
          margin-left: 20px;
        }

        /* Warning colors in rows */
        .correlation-high { background: rgba(248, 81, 73, 0.1) !important; }
        .correlation-medium { background: rgba(210, 153, 34, 0.1) !important; }

        /* Reactable overrides */
        .reactable { font-size: 0.875rem; }
        .rt-search { margin-bottom: 12px !important; }

        /* Filter Bar */
        .filter-bar {
          display: flex;
          gap: 12px;
          margin-bottom: 16px;
          flex-wrap: wrap;
          align-items: flex-end;
        }

        .filter-group {
          position: relative;
        }

        .filter-label {
          display: block;
          font-size: 0.7rem;
          color: #8b949e;
          text-transform: uppercase;
          letter-spacing: 0.5px;
          margin-bottom: 4px;
        }

        .filter-dropdown {
          background: #161b22;
          border: 1px solid #30363d;
          border-radius: 6px;
          padding: 8px 32px 8px 12px;
          color: #c9d1d9;
          font-size: 0.85rem;
          cursor: pointer;
          min-width: 160px;
          position: relative;
        }

        .filter-dropdown::after {
          content: "▼";
          position: absolute;
          right: 10px;
          top: 50%;
          transform: translateY(-50%);
          font-size: 0.6rem;
          color: #8b949e;
          pointer-events: none;
        }

        .filter-menu {
          display: none;
          position: absolute;
          top: 100%;
          left: 0;
          background: #161b22;
          border: 1px solid #30363d;
          border-radius: 6px;
          min-width: 200px;
          max-height: 300px;
          overflow-y: auto;
          z-index: 100;
          box-shadow: 0 4px 12px rgba(0,0,0,0.4);
          margin-top: 4px;
        }

        .filter-menu.open {
          display: block;
        }

        .filter-option {
          display: flex;
          align-items: center;
          padding: 8px 12px;
          cursor: pointer;
          font-size: 0.85rem;
        }

        .filter-option:hover {
          background: #21262d;
        }

        .filter-option input {
          margin-right: 8px;
          accent-color: #58a6ff;
        }

        .filter-option.select-all {
          border-bottom: 1px solid #30363d;
          font-weight: 500;
        }

        .clear-filters-btn {
          background: transparent;
          border: 1px solid #30363d;
          color: #8b949e;
          padding: 8px 16px;
          border-radius: 6px;
          font-size: 0.85rem;
          cursor: pointer;
          height: fit-content;
        }

        .clear-filters-btn:hover {
          border-color: #58a6ff;
          color: #58a6ff;
        }

        .filter-count {
          background: #58a6ff;
          color: #0d1117;
          font-size: 0.7rem;
          padding: 2px 6px;
          border-radius: 10px;
          margin-left: 6px;
        }

        /* Bankroll/Kelly inputs */
        .sizing-controls {
          display: flex;
          gap: 16px;
          align-items: flex-end;
          margin-bottom: 16px;
          padding: 12px 16px;
          background: #161b22;
          border: 1px solid #21262d;
          border-radius: 6px;
        }

        .sizing-group {
          display: flex;
          flex-direction: column;
          gap: 4px;
        }

        .sizing-label {
          font-size: 0.7rem;
          color: #8b949e;
          text-transform: uppercase;
          letter-spacing: 0.5px;
        }

        .sizing-input {
          background: #0d1117;
          border: 1px solid #30363d;
          border-radius: 6px;
          padding: 8px 12px;
          color: #c9d1d9;
          font-size: 0.9rem;
          width: 120px;
        }

        .sizing-input:focus {
          border-color: #58a6ff;
          outline: none;
        }

        .apply-sizing-btn {
          background: #238636;
          border: 1px solid #2ea043;
          color: #fff;
          padding: 8px 16px;
          border-radius: 6px;
          font-size: 0.85rem;
          cursor: pointer;
        }

        .apply-sizing-btn:hover {
          background: #2ea043;
        }

        /* Tab Navigation */
        .tab-nav {
          display: flex;
          gap: 4px;
          margin-bottom: 24px;
          border-bottom: 1px solid #21262d;
          padding-bottom: 0;
        }

        .tab-btn {
          background: transparent;
          border: none;
          color: #8b949e;
          padding: 12px 20px;
          font-size: 0.9rem;
          font-weight: 500;
          cursor: pointer;
          border-bottom: 2px solid transparent;
          margin-bottom: -1px;
          transition: all 0.15s;
        }

        .tab-btn:hover {
          color: #c9d1d9;
        }

        .tab-btn.active {
          color: #58a6ff;
          border-bottom-color: #58a6ff;
        }

        .tab-content {
          display: none;
        }

        .tab-content.active {
          display: block;
        }

        /* Parlay Builder styles */
        .parlay-container {
          display: grid;
          grid-template-columns: 1fr 350px;
          gap: 24px;
        }

        .leg-builder {
          background: #161b22;
          border: 1px solid #21262d;
          border-radius: 6px;
          padding: 20px;
        }

        /* Structured Leg Card */
        .leg-card {
          background: #0d1117;
          border: 1px solid #30363d;
          border-radius: 8px;
          padding: 16px;
          margin-bottom: 12px;
          position: relative;
        }

        .leg-card.has-error {
          border-color: #f85149;
        }

        .leg-header {
          display: flex;
          justify-content: space-between;
          align-items: center;
          margin-bottom: 12px;
        }

        .leg-number {
          font-size: 0.75rem;
          font-weight: 600;
          color: #8b949e;
          text-transform: uppercase;
        }

        .leg-remove-btn {
          background: transparent;
          border: none;
          color: #6e7681;
          cursor: pointer;
          font-size: 1.2rem;
          padding: 4px 8px;
          border-radius: 4px;
        }

        .leg-remove-btn:hover {
          color: #f85149;
          background: rgba(248, 81, 73, 0.1);
        }

        .leg-fields {
          display: grid;
          grid-template-columns: 100px 120px 100px 90px;
          gap: 10px;
          align-items: start;
        }

        .leg-field {
          display: flex;
          flex-direction: column;
        }

        .leg-field-label {
          font-size: 0.65rem;
          color: #6e7681;
          text-transform: uppercase;
          margin-bottom: 4px;
        }

        .leg-select, .leg-line-input {
          background: #161b22;
          border: 1px solid #30363d;
          color: #c9d1d9;
          padding: 8px 10px;
          border-radius: 6px;
          font-size: 0.85rem;
          width: 100%;
        }

        .leg-select:focus, .leg-line-input:focus {
          border-color: #58a6ff;
          outline: none;
        }

        .leg-field.has-error .leg-select,
        .leg-field.has-error .leg-line-input {
          border-color: #f85149;
        }

        .leg-error {
          font-size: 0.7rem;
          color: #f85149;
          margin-top: 4px;
        }

        .leg-field.hidden {
          display: none;
        }

        .add-leg-btn {
          background: transparent;
          border: 1px dashed #30363d;
          color: #8b949e;
          padding: 10px;
          border-radius: 6px;
          cursor: pointer;
          width: 100%;
          margin-bottom: 16px;
        }

        .add-leg-btn:hover {
          border-color: #58a6ff;
          color: #58a6ff;
        }

        .calculate-btn {
          background: #238636;
          border: 1px solid #2ea043;
          color: #fff;
          padding: 12px 24px;
          border-radius: 6px;
          font-size: 0.9rem;
          font-weight: 500;
          cursor: pointer;
          width: 100%;
        }

        .parlay-results {
          background: #161b22;
          border: 1px solid #21262d;
          border-radius: 6px;
          padding: 20px;
        }

        .result-row {
          display: flex;
          justify-content: space-between;
          padding: 12px 0;
          border-bottom: 1px solid #21262d;
        }

        .result-row:last-child {
          border-bottom: none;
        }

        .result-label {
          color: #8b949e;
          font-size: 0.85rem;
        }

        .result-value {
          color: #f0f6fc;
          font-weight: 600;
          font-size: 1rem;
        }

        .result-value.positive { color: #3fb950; }
        .result-value.negative { color: #f85149; }

        /* Props Calculator styles */
        .props-container {
          display: grid;
          grid-template-columns: 400px 1fr;
          gap: 24px;
        }

        .props-form {
          background: #161b22;
          border: 1px solid #21262d;
          border-radius: 6px;
          padding: 20px;
        }

        .form-group {
          margin-bottom: 16px;
        }

        .form-label {
          display: block;
          color: #8b949e;
          font-size: 0.75rem;
          text-transform: uppercase;
          letter-spacing: 0.5px;
          margin-bottom: 6px;
        }

        .form-select, .form-input {
          width: 100%;
          background: #0d1117;
          border: 1px solid #30363d;
          color: #c9d1d9;
          padding: 10px 14px;
          border-radius: 6px;
          font-size: 0.9rem;
        }

        .form-row {
          display: grid;
          grid-template-columns: 1fr 1fr;
          gap: 12px;
        }

        .props-results {
          background: #161b22;
          border: 1px solid #21262d;
          border-radius: 6px;
          padding: 20px;
        }

        .histogram-container {
          height: 300px;
          margin-top: 20px;
        }

        .props-stats {
          display: grid;
          grid-template-columns: repeat(3, 1fr);
          gap: 12px;
          margin-top: 20px;
        }

        .props-stat {
          background: #0d1117;
          border-radius: 6px;
          padding: 12px;
          text-align: center;
        }

        .props-stat-value {
          font-size: 1.2rem;
          font-weight: 600;
          color: #f0f6fc;
        }

        .props-stat-label {
          font-size: 0.7rem;
          color: #8b949e;
          text-transform: uppercase;
          margin-top: 4px;
        }
      '))
    ),

    tags$body(
      tags$div(class = "container",
        # Header
        tags$div(class = "header",
          tags$div(
            tags$h1("NFL +EV Dashboard"),
            tags$div(class = "subtitle", paste("Updated", timestamp))
          ),
          tags$button(class = "refresh-btn", onclick = "refreshData()", "Refresh")
        ),

        # Tab Navigation
        tags$div(class = "tab-nav",
          tags$button(class = "tab-btn active", onclick = "switchTab('bets')", "Bets"),
          tags$button(class = "tab-btn", onclick = "switchTab('parlay')", "Parlay Builder"),
          tags$button(class = "tab-btn", onclick = "switchTab('props')", "Props Calculator")
        ),

        # ============ BETS TAB ============
        tags$div(id = "tab-bets", class = "tab-content active",
          # Stats
          tags$div(class = "stats-row",
            tags$div(class = "stat-card",
              tags$div(class = "stat-value", stats$total_bets),
              tags$div(class = "stat-label", "Total Bets")
            ),
            tags$div(class = "stat-card",
              tags$div(class = "stat-value", stats$placed_count),
              tags$div(class = "stat-label", "Placed")
            ),
            tags$div(class = "stat-card",
              tags$div(class = "stat-value", sprintf("+%.1f%%", stats$avg_ev)),
              tags$div(class = "stat-label", "Avg EV")
            ),
            tags$div(class = "stat-card",
              tags$div(class = "stat-value", sprintf("+%.1f%%", stats$max_ev)),
              tags$div(class = "stat-label", "Best EV")
            )
          ),

          # Placed Bets Section (if any)
          if (!is.null(placed_table)) {
            tagList(
              tags$div(class = "section-header", "Placed Bets"),
              tags$div(class = "table-container placed-section", placed_table)
            )
          },

          # Bankroll/Kelly Controls
          tags$div(class = "sizing-controls",
            tags$div(class = "sizing-group",
              tags$span(class = "sizing-label", "Bankroll ($)"),
              tags$input(id = "bankroll-input", class = "sizing-input", type = "number",
                value = "100", min = "1", step = "10")
            ),
            tags$div(class = "sizing-group",
              tags$span(class = "sizing-label", "Kelly Fraction"),
              tags$input(id = "kelly-input", class = "sizing-input", type = "number",
                value = "0.25", min = "0.01", max = "1", step = "0.05")
            ),
            tags$button(class = "apply-sizing-btn", onclick = "recalculateBetSizes()", "Apply")
          ),

          # Filter Bar
          tags$div(class = "filter-bar",
            tags$div(class = "filter-group",
              tags$span(class = "filter-label", "Game"),
              tags$div(class = "filter-dropdown", id = "filter-game-btn", onclick = "toggleFilter('game')",
                tags$span(id = "filter-game-text", "All Games")
              ),
              tags$div(class = "filter-menu", id = "filter-game-menu")
            ),
            tags$div(class = "filter-group",
              tags$span(class = "filter-label", "Book"),
              tags$div(class = "filter-dropdown", id = "filter-book-btn", onclick = "toggleFilter('book')",
                tags$span(id = "filter-book-text", "All Books")
              ),
              tags$div(class = "filter-menu", id = "filter-book-menu")
            ),
            tags$div(class = "filter-group",
              tags$span(class = "filter-label", "Market"),
              tags$div(class = "filter-dropdown", id = "filter-market-btn", onclick = "toggleFilter('market')",
                tags$span(id = "filter-market-text", "All Markets")
              ),
              tags$div(class = "filter-menu", id = "filter-market-menu")
            ),
            tags$button(class = "clear-filters-btn", onclick = "clearAllFilters()", "Clear Filters")
          ),

          # Available Bets Table
          tags$div(class = "section-header",
            tags$span("Available Bets"),
            tags$span(id = "filtered-count", style = "margin-left: 8px; color: #8b949e; font-size: 0.8rem;")
          ),
          tags$div(class = "table-container", bets_table)
        ),

        # ============ PARLAY TAB ============
        tags$div(id = "tab-parlay", class = "tab-content",
          tags$div(class = "parlay-container",
            # Left: Leg Builder
            tags$div(class = "leg-builder",
              tags$h3(style = "margin: 0 0 16px 0; color: #f0f6fc; font-size: 1rem;", "Build Your Parlay"),

              # Game selector at top
              tags$div(class = "form-group", style = "margin-bottom: 20px;",
                tags$label(class = "form-label", "Game (optional - leave blank for all games)"),
                tags$select(id = "parlay-game", class = "form-select",
                  tags$option(value = "", "All Games")
                )
              ),

              tags$p(style = "color: #8b949e; font-size: 0.85rem; margin-bottom: 16px;",
                "Add legs: period, bet type, and line"),
              tags$div(id = "parlay-legs"),
              tags$button(class = "add-leg-btn", onclick = "addLeg()", "+ Add Leg"),
              tags$button(class = "calculate-btn", onclick = "calculateParlay()", "Calculate Fair Odds")
            ),

            # Right: Results
            tags$div(class = "parlay-results",
              tags$h3(style = "margin: 0 0 16px 0; color: #f0f6fc; font-size: 1rem;", "Results"),
              tags$div(id = "parlay-results-content",
                tags$p(style = "color: #6e7681; text-align: center; padding: 40px 0;", "Add legs and click Calculate")
              )
            )
          )
        ),

        # ============ PROPS TAB ============
        tags$div(id = "tab-props", class = "tab-content",
          tags$div(class = "props-container",
            # Left: Form
            tags$div(class = "props-form",
              tags$h3(style = "margin: 0 0 16px 0; color: #f0f6fc; font-size: 1rem;", "Prop Calculator"),

              # Mode toggle
              tags$div(style = "display: flex; gap: 8px; margin-bottom: 20px;",
                tags$button(id = "props-mode-simple", class = "tab-btn active",
                  style = "flex: 1; padding: 8px;", onclick = "setPropsMode('simple')", "Simple"),
                tags$button(id = "props-mode-custom", class = "tab-btn",
                  style = "flex: 1; padding: 8px;", onclick = "setPropsMode('custom')", "Custom SQL")
              ),

              # Simple mode form
              tags$div(id = "props-simple-form",
                tags$div(class = "form-group",
                  tags$label(class = "form-label", "Stat Column"),
                  tags$select(id = "prop-column", class = "form-select",
                    tags$option(value = "", "Loading columns...")
                  )
                ),

                tags$div(class = "form-group",
                  tags$label(class = "form-label", "Condition"),
                  tags$input(id = "prop-condition", class = "form-input", type = "text",
                    placeholder = "e.g., >= 1, > 0, == TRUE")
                ),

                tags$div(class = "form-row",
                  tags$div(class = "form-group",
                    tags$label(class = "form-label", "Period (optional)"),
                    tags$select(id = "prop-period", class = "form-select",
                      tags$option(value = "", "Full Game"),
                      tags$option(value = "1", "Q1"),
                      tags$option(value = "2", "Q2"),
                      tags$option(value = "3", "Q3"),
                      tags$option(value = "4", "Q4"),
                      tags$option(value = "Half1", "1H"),
                      tags$option(value = "Half2", "2H")
                    )
                  ),
                  tags$div(class = "form-group",
                    tags$label(class = "form-label", "Team (optional)"),
                    tags$select(id = "prop-team", class = "form-select",
                      tags$option(value = "", "Any Team"),
                      tags$option(value = "home", "Home"),
                      tags$option(value = "away", "Away")
                    )
                  )
                ),

                tags$button(class = "calculate-btn", onclick = "calculateProps()", "Calculate Fair Odds")
              ),

              # Custom mode form
              tags$div(id = "props-custom-form", style = "display: none;",
                tags$div(class = "form-group",
                  tags$label(class = "form-label", "SQL WHERE Condition"),
                  tags$textarea(id = "prop-custom-query", class = "form-input",
                    style = "min-height: 100px; font-family: monospace; resize: vertical;",
                    placeholder = "e.g., touchdown >= 1 AND pass_touchdown >= 1")
                ),
                tags$p(style = "color: #6e7681; font-size: 0.75rem; margin: 8px 0 16px 0;",
                  "Query runs against nfl_pbp table. Use column names like: touchdown, pass_touchdown, rush_touchdown, field_goal_attempt, interception, fumble, sack, etc."),
                tags$button(class = "calculate-btn", onclick = "calculatePropsCustom()", "Calculate Fair Odds")
              )
            ),

            # Right: Results + Histogram
            tags$div(class = "props-results",
              tags$h3(style = "margin: 0 0 16px 0; color: #f0f6fc; font-size: 1rem;", "Results"),
              tags$div(id = "props-results-content",
                tags$p(style = "color: #6e7681; text-align: center; padding: 40px 0;", "Select a stat and condition")
              ),
              tags$div(class = "histogram-container",
                tags$canvas(id = "props-histogram")
              )
            )
          )
        )
      ),

      # JavaScript
      tags$script(HTML('
        // ============ BETS FILTERING ============
        const activeFilters = {
          game: new Set(),
          book: new Set(),
          market: new Set()
        };

        document.addEventListener("DOMContentLoaded", function() {
          initBetsFilters();
        });

        function initBetsFilters() {
          const opts = window.FILTER_OPTIONS;
          if (!opts) return;

          // Populate game filter
          populateFilterMenu("game", opts.games);
          populateFilterMenu("book", opts.books);
          populateFilterMenu("market", opts.markets);

          // Close menus when clicking outside
          document.addEventListener("click", function(e) {
            if (!e.target.closest(".filter-group")) {
              document.querySelectorAll(".filter-menu").forEach(m => m.classList.remove("open"));
            }
          });
        }

        function populateFilterMenu(type, options) {
          const menu = document.getElementById("filter-" + type + "-menu");
          if (!menu || !options || !Array.isArray(options)) return;

          let html = "<label class=\\"filter-option select-all\\">" +
            "<input type=\\"checkbox\\" checked onchange=\\"toggleSelectAll(\'" + type + "\', this.checked)\\" /> " +
            "Select All</label>";

          for (var i = 0; i < options.length; i++) {
            var opt = options[i];
            var safeVal = String(opt).replace(/"/g, "&quot;");
            html += "<label class=\\"filter-option\\">" +
              "<input type=\\"checkbox\\" checked data-val=\\"" + safeVal + "\\" onchange=\\"updateFilter(\'" + type + "\')\\"/> " +
              escapeHtml(String(opt)) + "</label>";
          }

          menu.innerHTML = html;
        }

        function toggleFilter(type) {
          const menu = document.getElementById("filter-" + type + "-menu");
          const wasOpen = menu.classList.contains("open");

          // Close all menus first
          document.querySelectorAll(".filter-menu").forEach(m => m.classList.remove("open"));

          // Toggle this one
          if (!wasOpen) {
            menu.classList.add("open");
          }
        }

        function toggleSelectAll(type, isChecked) {
          const menu = document.getElementById("filter-" + type + "-menu");
          menu.querySelectorAll("input[type=checkbox]").forEach(function(cb) {
            cb.checked = isChecked;
          });
          updateFilter(type);
        }

        function updateFilter(type) {
          const menu = document.getElementById("filter-" + type + "-menu");
          const checkboxes = menu.querySelectorAll("input[data-val]");
          const checkedVals = [];

          for (var i = 0; i < checkboxes.length; i++) {
            if (checkboxes[i].checked) {
              checkedVals.push(checkboxes[i].getAttribute("data-val"));
            }
          }

          activeFilters[type] = new Set(checkedVals);

          // Update button text
          const textEl = document.getElementById("filter-" + type + "-text");
          const allLabel = type === "game" ? "All Games" : type === "book" ? "All Books" : "All Markets";

          if (checkedVals.length === checkboxes.length) {
            textEl.textContent = allLabel;
          } else if (checkedVals.length === 0) {
            textEl.textContent = "None selected";
          } else if (checkedVals.length <= 2) {
            textEl.textContent = checkedVals.join(", ");
          } else {
            textEl.textContent = checkedVals.length + " selected";
          }

          // Update select all checkbox
          const selectAllCb = menu.querySelector(".select-all input");
          if (selectAllCb) {
            selectAllCb.checked = checkedVals.length === checkboxes.length;
            selectAllCb.indeterminate = checkedVals.length > 0 && checkedVals.length < checkboxes.length;
          }

          applyFilters();
        }

        function applyFilters() {
          // Only filter the Available Bets table, not Placed Bets
          // The Available Bets table is in the table-container that is NOT in placed-section
          const containers = document.querySelectorAll("#tab-bets .table-container:not(.placed-section)");
          const table = containers.length > 0 ? containers[containers.length - 1].querySelector(".rt-tbody") : null;
          if (!table) return;

          const rows = table.querySelectorAll(".rt-tr-group");
          let visibleCount = 0;
          let totalCount = rows.length;

          rows.forEach(row => {
            const cells = row.querySelectorAll(".rt-td");
            // Find the game, book, and market cells by their content/position
            // Game is typically in the first visible column (index 1 after warning)
            // Market is index 3, Book is index 7

            let gameText = "";
            let bookText = "";
            let marketText = "";

            cells.forEach((cell, idx) => {
              const text = cell.textContent.trim();
              // Game column (contains @)
              if (text.includes("@")) {
                gameText = text;
              }
              // Market column (contains specific patterns)
              if (["ML", "Spread", "Total", "Team Tot", "Alt"].some(m => text.includes(m))) {
                marketText = text;
              }
            });

            // Book is in the last data column before Action
            const bookCell = cells[cells.length - 2];
            if (bookCell) bookText = bookCell.textContent.trim();

            // Determine market type from market text
            let marketType = "Other";
            if (marketText.includes("ML-3") || marketText.includes("ML 3")) {
              marketType = "ML 3-Way";
            } else if (marketText.includes("ML")) {
              marketType = "Moneyline";
            } else if (marketText.includes("Alt")) {
              marketType = "Alternates";
            } else if (marketText.includes("Team Tot")) {
              marketType = "Team Totals";
            } else if (marketText.includes("Spread")) {
              marketType = "Spreads";
            } else if (marketText.includes("Total")) {
              marketType = "Totals";
            }

            // Check filters - if filter has values, row must match one of them
            var gameMatch = activeFilters.game.size === 0 || activeFilters.game.has(gameText);
            var bookMatch = activeFilters.book.size === 0 || activeFilters.book.has(bookText);
            var marketMatch = activeFilters.market.size === 0 || activeFilters.market.has(marketType);

            var visible = gameMatch && bookMatch && marketMatch;
            row.style.display = visible ? "" : "none";
            if (visible) visibleCount++;
          });

          // Update count display
          var countEl = document.getElementById("filtered-count");
          if (countEl) {
            if (visibleCount === totalCount) {
              countEl.textContent = "";
            } else {
              countEl.textContent = "(" + visibleCount + " of " + totalCount + ")";
            }
          }
        }

        function clearAllFilters() {
          ["game", "book", "market"].forEach(function(type) {
            var menu = document.getElementById("filter-" + type + "-menu");
            var checkboxes = menu.querySelectorAll("input[type=checkbox]");
            var vals = [];

            for (var i = 0; i < checkboxes.length; i++) {
              checkboxes[i].checked = true;
              var val = checkboxes[i].getAttribute("data-val");
              if (val) vals.push(val);
            }

            activeFilters[type] = new Set(vals);
            var textEl = document.getElementById("filter-" + type + "-text");
            var allLabel = type === "game" ? "All Games" : type === "book" ? "All Books" : "All Markets";
            textEl.textContent = allLabel;
          });
          applyFilters();
        }

        // ============ BET SIZING ============
        function recalculateBetSizes() {
          var bankroll = parseFloat(document.getElementById("bankroll-input").value) || 100;
          var kellyMult = parseFloat(document.getElementById("kelly-input").value) || 0.25;

          if (bankroll <= 0 || kellyMult <= 0 || kellyMult > 1) {
            showToast("Invalid bankroll or kelly value", "error");
            return;
          }

          // Get all rows in the Available Bets table
          var containers = document.querySelectorAll("#tab-bets .table-container:not(.placed-section)");
          var table = containers.length > 0 ? containers[containers.length - 1] : null;
          if (!table) return;

          var rows = table.querySelectorAll(".rt-tr-group");

          rows.forEach(function(row) {
            var cells = row.querySelectorAll(".rt-td");

            // Find prob and odds from hidden cells or data
            // Based on column order: prob is hidden, odds is visible
            var prob = null;
            var odds = null;
            var sizeCell = null;

            cells.forEach(function(cell, idx) {
              var text = cell.textContent.trim();

              // Odds cell (contains + or - followed by numbers)
              if (/^[+-]\\d+$/.test(text)) {
                odds = parseInt(text);
              }

              // Size cell (contains $)
              if (text.startsWith("$")) {
                sizeCell = cell;
              }
            });

            // Get prob from data attribute on the row or from hidden cell
            // We need to add this - for now, extract from button data
            var btn = row.querySelector("button[data-prob]");
            if (btn) {
              prob = parseFloat(btn.getAttribute("data-prob"));
            }

            // Calculate new bet size if we have the data
            if (prob && odds && sizeCell) {
              var newSize = calculateKellyBet(prob, odds, bankroll, kellyMult);
              sizeCell.textContent = "$" + newSize.toFixed(2);

              // Also update the button data-size attribute
              if (btn) {
                btn.setAttribute("data-size", newSize.toFixed(2));
              }
            }
          });

          showToast("Bet sizes updated", "success");
        }

        function calculateKellyBet(prob, americanOdds, bankroll, kellyMult) {
          // Convert American odds to decimal
          var decimalOdds;
          if (americanOdds > 0) {
            decimalOdds = 1 + (americanOdds / 100);
          } else {
            decimalOdds = 1 + (100 / Math.abs(americanOdds));
          }

          // Calculate edge and kelly fraction
          var edge = (prob * decimalOdds) - 1;

          // If no edge, no bet
          if (edge <= 0) {
            return 0;
          }

          var kellyFraction = edge / (decimalOdds - 1);
          var betSize = bankroll * kellyFraction * kellyMult;

          // Cap at bankroll
          return Math.min(betSize, bankroll);
        }

        function refreshData() {
          const btn = document.querySelector(".refresh-btn");
          btn.textContent = "Refreshing...";
          btn.disabled = true;

          fetch("/refresh", { method: "POST" })
            .then(r => r.json())
            .then(data => {
              if (data.success) {
                showToast("Refreshed!", "success");
                setTimeout(() => location.reload(), 800);
              } else {
                showToast("Error: " + data.error, "error");
                btn.textContent = "Refresh";
                btn.disabled = false;
              }
            })
            .catch(() => {
              showToast("Server error", "error");
              btn.textContent = "Refresh";
              btn.disabled = false;
            });
        }

        function placeBet(btn) {
          const data = {
            bet_hash: btn.dataset.hash,
            game_id: btn.dataset.gameId,
            home_team: btn.dataset.home,
            away_team: btn.dataset.away,
            game_time: btn.dataset.time,
            market: btn.dataset.market,
            bet_on: btn.dataset.betOn,
            line: btn.dataset.line === "" ? null : parseFloat(btn.dataset.line),
            model_prob: parseFloat(btn.dataset.prob),
            model_ev: parseFloat(btn.dataset.ev),
            recommended_size: parseFloat(btn.dataset.size),
            odds: parseInt(btn.dataset.odds),
            bookmaker: btn.dataset.book
          };

          fetch("/api/place-bet", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify(data)
          })
            .then(r => r.json())
            .then(result => {
              if (result.success) {
                btn.className = "btn-placed";
                btn.textContent = "Placed";
                btn.onclick = function() { removeBet(this); };
                showToast("Bet placed!", "success");
              } else {
                showToast(result.error, "error");
              }
            });
        }

        function removeBet(btn) {
          fetch("/api/remove-bet", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ bet_hash: btn.dataset.hash })
          })
            .then(r => r.json())
            .then(result => {
              if (result.success) {
                btn.className = "btn-place";
                btn.textContent = "Place";
                btn.onclick = function() { placeBet(this); };
                showToast("Removed", "success");
              }
            });
        }

        function showToast(msg, type) {
          const t = document.createElement("div");
          t.className = "toast " + type;
          t.textContent = msg;
          document.body.appendChild(t);
          setTimeout(() => t.remove(), 3000);
        }

        // ============ TAB SWITCHING ============
        function switchTab(tabName) {
          // Update buttons
          document.querySelectorAll(".tab-btn").forEach(btn => btn.classList.remove("active"));
          event.target.classList.add("active");

          // Update content
          document.querySelectorAll(".tab-content").forEach(tab => tab.classList.remove("active"));
          document.getElementById("tab-" + tabName).classList.add("active");

          // Load columns for props tab on first visit
          if (tabName === "props" && !window.propsColumnsLoaded) {
            loadPropsColumns();
          }
        }

        // ============ PARLAY BUILDER ============
        let legCounter = 0;

        // Initialize with 2 legs on page load and populate game dropdown
        document.addEventListener("DOMContentLoaded", function() {
          // Populate game dropdown at top
          const gameSelect = document.getElementById("parlay-game");
          if (gameSelect && window.GAMES) {
            window.GAMES.forEach(g => {
              const opt = document.createElement("option");
              opt.value = g.id;
              opt.textContent = g.label;
              gameSelect.appendChild(opt);
            });
          }

          addLeg();
          addLeg();
        });

        function addLeg() {
          legCounter++;
          const container = document.getElementById("parlay-legs");
          const card = document.createElement("div");
          card.className = "leg-card";
          card.dataset.legId = legCounter;

          card.innerHTML = `
            <div class="leg-header">
              <span class="leg-number">Leg ${container.children.length + 1}</span>
              <button class="leg-remove-btn" onclick="removeLeg(this)">&times;</button>
            </div>
            <div class="leg-fields">
              <div class="leg-field" data-field="period">
                <span class="leg-field-label">Period</span>
                <select class="leg-select" name="period">
                  <option value="FG">Full Game</option>
                  <option value="1H">1H</option>
                  <option value="2H">2H</option>
                  <option value="1Q">Q1</option>
                  <option value="2Q">Q2</option>
                  <option value="3Q">Q3</option>
                  <option value="4Q">Q4</option>
                </select>
              </div>
              <div class="leg-field" data-field="betType">
                <span class="leg-field-label">Bet Type</span>
                <select class="leg-select" name="betType" onchange="updateLegFields(this)">
                  <option value="spread">Spread</option>
                  <option value="total">Total</option>
                  <option value="teamTotal">Team Total</option>
                  <option value="ml">Moneyline</option>
                  <option value="ml3">ML 3-Way</option>
                </select>
              </div>
              <div class="leg-field" data-field="side">
                <span class="leg-field-label">Side</span>
                <select class="leg-select" name="side">
                  <option value="home">Home</option>
                  <option value="away">Away</option>
                </select>
                <span class="leg-error"></span>
              </div>
              <div class="leg-field" data-field="line">
                <span class="leg-field-label">Line</span>
                <input class="leg-line-input" type="text" name="line" placeholder="-3">
                <span class="leg-error"></span>
              </div>
            </div>
          `;

          container.appendChild(card);
          updateLegNumbers();
        }

        function removeLeg(btn) {
          const container = document.getElementById("parlay-legs");
          if (container.children.length > 2) {
            btn.closest(".leg-card").remove();
            updateLegNumbers();
          } else {
            showToast("Minimum 2 legs required", "error");
          }
        }

        function updateLegNumbers() {
          const cards = document.querySelectorAll("#parlay-legs .leg-card");
          cards.forEach((card, i) => {
            card.querySelector(".leg-number").textContent = `Leg ${i + 1}`;
          });
        }

        function updateLegFields(selectEl) {
          const card = selectEl.closest(".leg-card");
          const betType = selectEl.value;
          const sideField = card.querySelector("[data-field=\\"side\\"]");
          const lineField = card.querySelector("[data-field=\\"line\\"]");
          const sideSelect = sideField.querySelector("select");

          // Reset visibility
          sideField.classList.remove("hidden");
          lineField.classList.remove("hidden");

          // Update side options based on bet type
          if (betType === "spread" || betType === "ml") {
            sideSelect.innerHTML = `
              <option value="home">Home</option>
              <option value="away">Away</option>
            `;
            if (betType === "ml") lineField.classList.add("hidden");
          } else if (betType === "total") {
            sideSelect.innerHTML = `
              <option value="over">Over</option>
              <option value="under">Under</option>
            `;
          } else if (betType === "teamTotal") {
            sideSelect.innerHTML = `
              <option value="home_over">Home Over</option>
              <option value="home_under">Home Under</option>
              <option value="away_over">Away Over</option>
              <option value="away_under">Away Under</option>
            `;
          } else if (betType === "ml3") {
            sideSelect.innerHTML = `
              <option value="home">Home</option>
              <option value="away">Away</option>
              <option value="tie">Tie</option>
            `;
            lineField.classList.add("hidden");
          }
        }

        function validateLeg(card) {
          let isValid = true;
          const betType = card.querySelector("[name=\\"betType\\"]").value;

          // Clear previous errors
          card.querySelectorAll(".leg-field").forEach(f => {
            f.classList.remove("has-error");
            const errEl = f.querySelector(".leg-error");
            if (errEl) errEl.textContent = "";
          });
          card.classList.remove("has-error");

          // Game is optional - no validation needed

          // Check line for spread, total, teamTotal
          if (["spread", "total", "teamTotal"].includes(betType)) {
            const lineField = card.querySelector("[data-field=\\"line\\"]");
            const lineVal = card.querySelector("[name=\\"line\\"]").value.trim();
            if (!lineVal) {
              lineField.classList.add("has-error");
              lineField.querySelector(".leg-error").textContent = "Required";
              isValid = false;
            } else if (isNaN(parseFloat(lineVal))) {
              lineField.classList.add("has-error");
              lineField.querySelector(".leg-error").textContent = "Invalid";
              isValid = false;
            }
          }

          if (!isValid) card.classList.add("has-error");
          return isValid;
        }

        function buildLegString(card) {
          const period = card.querySelector("[name=\\"period\\"]").value;
          const betType = card.querySelector("[name=\\"betType\\"]").value;
          const side = card.querySelector("[name=\\"side\\"]").value;
          const line = card.querySelector("[name=\\"line\\"]").value.trim();

          // Format: "PERIOD SIDE LINE" matching parlay.R expectations
          // Examples: "1H home -3", "1Q over 10.5", "1H home ML"
          let legStr = period + " ";

          if (betType === "spread") {
            // "1H home -3"
            legStr += side + " " + line;
          } else if (betType === "total") {
            // "1H over 22.5"
            legStr += side + " " + line;
          } else if (betType === "teamTotal") {
            // "1H home over 10.5" - side is like "home_over"
            const parts = side.split("_");
            legStr += parts[0] + " " + parts[1] + " " + line;
          } else if (betType === "ml") {
            // "1H home ML"
            legStr += side + " ML";
          } else if (betType === "ml3") {
            // "1H home ML3" or "1Q tie"
            if (side === "tie") {
              legStr += "tie";
            } else {
              legStr += side + " ML3";
            }
          }

          return legStr.trim();
        }

        function calculateParlay() {
          const cards = document.querySelectorAll("#parlay-legs .leg-card");

          if (cards.length < 2) {
            showToast("Add at least 2 legs", "error");
            return;
          }

          // Validate all legs
          let allValid = true;
          cards.forEach(card => {
            if (!validateLeg(card)) allValid = false;
          });

          if (!allValid) {
            showToast("Please fix errors above", "error");
            return;
          }

          // Build leg strings
          const legs = Array.from(cards).map(buildLegString);

          const resultsDiv = document.getElementById("parlay-results-content");
          resultsDiv.innerHTML = "<p style=\\"color: #8b949e; text-align: center;\\">Calculating...</p>";

          fetch("/api/parlay", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ legs: legs })
          })
            .then(r => r.json())
            .then(data => {
              if (data.success) {
                // Parse the output to extract key values
                const output = data.output;
                const parsed = parseParlayOutput(output);

                if (parsed) {
                  resultsDiv.innerHTML = `
                    <div class="props-stats" style="margin-bottom: 16px;">
                      <div class="props-stat">
                        <div class="props-stat-value">${parsed.probability}</div>
                        <div class="props-stat-label">Probability</div>
                      </div>
                      <div class="props-stat">
                        <div class="props-stat-value" style="color: #3fb950;">${parsed.americanOdds}</div>
                        <div class="props-stat-label">Fair Odds</div>
                      </div>
                      <div class="props-stat">
                        <div class="props-stat-value">${parsed.decimalOdds}</div>
                        <div class="props-stat-label">Decimal</div>
                      </div>
                    </div>
                    <div style="background: #0d1117; padding: 12px; border-radius: 6px; margin-bottom: 12px;">
                      <div style="font-size: 0.75rem; color: #8b949e; margin-bottom: 8px;">LEGS</div>
                      ${parsed.legs.map((leg, i) => `
                        <div style="display: flex; justify-content: space-between; padding: 6px 0; border-bottom: 1px solid #21262d;">
                          <span style="color: #c9d1d9;">${i + 1}. ${leg.description}</span>
                          <span style="color: #8b949e;">${leg.prob}</span>
                        </div>
                      `).join("")}
                    </div>
                    <div style="display: flex; justify-content: space-between; padding: 8px 12px; background: #0d1117; border-radius: 6px;">
                      <span style="color: #8b949e;">Correlation</span>
                      <span style="color: ${parsed.correlationValue > 1 ? "#3fb950" : parsed.correlationValue < 1 ? "#f85149" : "#c9d1d9"};">
                        ${parsed.correlation}
                      </span>
                    </div>
                  `;
                } else {
                  // Fallback to raw output
                  resultsDiv.innerHTML = `
                    <pre style="background: #0d1117; padding: 16px; border-radius: 6px; color: #c9d1d9; font-size: 0.8rem; white-space: pre-wrap; overflow-x: auto;">${escapeHtml(output)}</pre>
                  `;
                }
              } else {
                resultsDiv.innerHTML = `<p style="color: #f85149; text-align: center;">${escapeHtml(data.error)}</p>`;
              }
            })
            .catch(err => {
              resultsDiv.innerHTML = "<p style=\\"color: #f85149; text-align: center;\\">Server error</p>";
            });
        }

        function parseParlayOutput(output) {
          try {
            // Extract legs
            const legMatches = output.match(/\\d+\\.\\s+(.+?)\\s+\\((\\d+\\.?\\d*%?)\\s+individual\\)/g);
            const legs = [];
            if (legMatches) {
              legMatches.forEach(match => {
                const m = match.match(/(\\d+)\\.\\s+(.+?)\\s+\\((\\d+\\.?\\d*%?)\\s+individual\\)/);
                if (m) {
                  legs.push({ description: m[2], prob: m[3] });
                }
              });
            }

            // Extract joint probability
            const probMatch = output.match(/Joint probability:\\s+(\\d+\\.?\\d*%)/);
            const probability = probMatch ? probMatch[1] : "N/A";

            // Extract American odds
            const americanMatch = output.match(/Fair American odds:\\s+([+-]?\\d+)/);
            const americanOdds = americanMatch ? americanMatch[1] : "N/A";

            // Extract decimal odds
            const decimalMatch = output.match(/Fair decimal odds:\\s+(\\d+\\.?\\d*)/);
            const decimalOdds = decimalMatch ? decimalMatch[1] : "N/A";

            // Extract correlation
            const corrMatch = output.match(/Correlation factor:\\s+(\\d+\\.?\\d*)\\s+\\((.+?)\\)/);
            const correlationValue = corrMatch ? parseFloat(corrMatch[1]) : 1;
            const correlation = corrMatch ? `${corrMatch[1]}x (${corrMatch[2]})` : "1.0x";

            if (legs.length === 0) return null;

            return { legs, probability, americanOdds, decimalOdds, correlation, correlationValue };
          } catch (e) {
            return null;
          }
        }

        function escapeHtml(text) {
          const div = document.createElement("div");
          div.textContent = text;
          return div.innerHTML;
        }

        // ============ PROPS CALCULATOR ============
        let propsChart = null;
        window.propsColumnsLoaded = false;
        let propsMode = "simple";

        function setPropsMode(mode) {
          propsMode = mode;
          document.getElementById("props-mode-simple").classList.toggle("active", mode === "simple");
          document.getElementById("props-mode-custom").classList.toggle("active", mode === "custom");
          document.getElementById("props-simple-form").style.display = mode === "simple" ? "block" : "none";
          document.getElementById("props-custom-form").style.display = mode === "custom" ? "block" : "none";

          // Load columns on first visit to simple mode
          if (mode === "simple" && !window.propsColumnsLoaded) {
            loadPropsColumns();
          }
        }

        function loadPropsColumns() {
          fetch("/api/columns")
            .then(r => r.json())
            .then(data => {
              if (data.success) {
                const select = document.getElementById("prop-column");
                select.innerHTML = "<option value=\\\"\\\">Select a column...</option>";

                // Add common columns first
                if (data.common && data.common.length > 0) {
                  const optgroup = document.createElement("optgroup");
                  optgroup.label = "Common Props";
                  data.common.forEach(col => {
                    const opt = document.createElement("option");
                    opt.value = col;
                    opt.textContent = col;
                    optgroup.appendChild(opt);
                  });
                  select.appendChild(optgroup);
                }

                // Add all columns
                const allGroup = document.createElement("optgroup");
                allGroup.label = "All Columns";
                data.columns.forEach(col => {
                  const opt = document.createElement("option");
                  opt.value = col;
                  opt.textContent = col;
                  allGroup.appendChild(opt);
                });
                select.appendChild(allGroup);

                window.propsColumnsLoaded = true;
              }
            });
        }

        function calculateProps() {
          const column = document.getElementById("prop-column").value;
          const condition = document.getElementById("prop-condition").value;
          const period = document.getElementById("prop-period").value;
          const team = document.getElementById("prop-team").value;

          if (!column || !condition) {
            showToast("Select a column and enter a condition", "error");
            return;
          }

          const resultsDiv = document.getElementById("props-results-content");
          resultsDiv.innerHTML = "<p style=\'color: #8b949e; text-align: center;\'>Calculating...</p>";

          fetch("/api/props", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({
              column: column,
              condition: condition,
              period: period || null,
              team: team || null
            })
          })
            .then(r => r.json())
            .then(data => {
              if (data.success && data.result) {
                const r = data.result;
                const probPct = (r.prob * 100).toFixed(1);
                const fairOdds = r.fair_american_odds > 0 ? "+" + r.fair_american_odds : r.fair_american_odds;

                resultsDiv.innerHTML = `
                  <div class="props-stats">
                    <div class="props-stat">
                      <div class="props-stat-value">${probPct}%</div>
                      <div class="props-stat-label">Probability</div>
                    </div>
                    <div class="props-stat">
                      <div class="props-stat-value">${fairOdds}</div>
                      <div class="props-stat-label">Fair Odds</div>
                    </div>
                    <div class="props-stat">
                      <div class="props-stat-value">${r.n_hits}/${r.n_games}</div>
                      <div class="props-stat-label">Hit Rate</div>
                    </div>
                  </div>
                  <div style="margin-top: 16px; padding: 12px; background: #0d1117; border-radius: 6px;">
                    <div style="font-size: 0.75rem; color: #8b949e; text-transform: uppercase; margin-bottom: 8px;">Distribution</div>
                    <div style="display: grid; grid-template-columns: repeat(6, 1fr); gap: 8px; font-size: 0.8rem;">
                      <div><span style="color: #8b949e;">Min:</span> ${r.stat_min?.toFixed(2) ?? "N/A"}</div>
                      <div><span style="color: #8b949e;">Q1:</span> ${r.stat_q1?.toFixed(2) ?? "N/A"}</div>
                      <div><span style="color: #8b949e;">Med:</span> ${r.stat_median?.toFixed(2) ?? "N/A"}</div>
                      <div><span style="color: #8b949e;">Mean:</span> ${r.stat_mean?.toFixed(2) ?? "N/A"}</div>
                      <div><span style="color: #8b949e;">Q3:</span> ${r.stat_q3?.toFixed(2) ?? "N/A"}</div>
                      <div><span style="color: #8b949e;">Max:</span> ${r.stat_max?.toFixed(2) ?? "N/A"}</div>
                    </div>
                  </div>
                `;

                // Update histogram if we have distribution data
                updateHistogram(r);
              } else if (data.output) {
                resultsDiv.innerHTML = `<pre style="background: #0d1117; padding: 16px; border-radius: 6px; color: #c9d1d9; font-size: 0.85rem; white-space: pre-wrap;">${escapeHtml(data.output)}</pre>`;
              } else {
                resultsDiv.innerHTML = `<p style="color: #f85149; text-align: center;">${escapeHtml(data.error || "Unknown error")}</p>`;
              }
            })
            .catch(err => {
              resultsDiv.innerHTML = "<p style=\'color: #f85149; text-align: center;\'>Server error</p>";
            });
        }

        function calculatePropsCustom() {
          const query = document.getElementById("prop-custom-query").value.trim();

          if (!query) {
            showToast("Enter a SQL condition", "error");
            return;
          }

          const resultsDiv = document.getElementById("props-results-content");
          resultsDiv.innerHTML = "<p style=\'color: #8b949e; text-align: center;\'>Calculating...</p>";

          // Hide histogram for custom queries (no distribution data)
          const histContainer = document.querySelector(".histogram-container");
          if (histContainer) histContainer.style.display = "none";

          fetch("/api/props-custom", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ query: query })
          })
            .then(r => r.json())
            .then(data => {
              if (data.success && data.result) {
                const r = data.result;
                const probPct = (r.prob * 100).toFixed(1);
                const fairOdds = r.fair_american_odds ? (r.fair_american_odds > 0 ? "+" + r.fair_american_odds : r.fair_american_odds) : "N/A";

                resultsDiv.innerHTML = `
                  <div class="props-stats">
                    <div class="props-stat">
                      <div class="props-stat-value">${probPct}%</div>
                      <div class="props-stat-label">Probability</div>
                    </div>
                    <div class="props-stat">
                      <div class="props-stat-value">${fairOdds}</div>
                      <div class="props-stat-label">Fair Odds</div>
                    </div>
                    <div class="props-stat">
                      <div class="props-stat-value">${r.n_hits}/${r.n_games}</div>
                      <div class="props-stat-label">Hit Rate</div>
                    </div>
                  </div>
                  <div style="margin-top: 16px; padding: 12px; background: #0d1117; border-radius: 6px;">
                    <div style="font-size: 0.75rem; color: #8b949e; text-transform: uppercase; margin-bottom: 8px;">Query</div>
                    <code style="color: #58a6ff; font-size: 0.8rem;">${escapeHtml(r.query)}</code>
                  </div>
                `;
              } else {
                resultsDiv.innerHTML = `<p style="color: #f85149; text-align: center;">${escapeHtml(data.error || "Unknown error")}</p>`;
              }
            })
            .catch(err => {
              resultsDiv.innerHTML = "<p style=\'color: #f85149; text-align: center;\'>Server error</p>";
            });
        }

        function updateHistogram(result) {
          const ctx = document.getElementById("props-histogram");
          if (!ctx) return;

          if (propsChart) {
            propsChart.destroy();
          }

          // Use actual histogram data if available
          let labels, values, title;
          if (result.histogram && result.histogram.bins && result.histogram.counts) {
            labels = result.histogram.bins.map(b => String(b));
            values = result.histogram.counts;
            title = result.histogram.type === "discrete" ? "Value Frequency" : "Distribution";
          } else {
            // Fallback to summary stats
            labels = ["Min", "Q1", "Median", "Mean", "Q3", "Max"];
            values = [
              result.stat_min || 0,
              result.stat_q1 || 0,
              result.stat_median || 0,
              result.stat_mean || 0,
              result.stat_q3 || 0,
              result.stat_max || 0
            ];
            title = "Distribution Summary";
          }

          propsChart = new Chart(ctx, {
            type: "bar",
            data: {
              labels: labels,
              datasets: [{
                label: "Games",
                data: values,
                backgroundColor: "rgba(88, 166, 255, 0.6)",
                borderColor: "#58a6ff",
                borderWidth: 1
              }]
            },
            options: {
              responsive: true,
              maintainAspectRatio: false,
              plugins: {
                legend: { display: false },
                title: {
                  display: true,
                  text: title,
                  color: "#8b949e"
                }
              },
              scales: {
                x: {
                  ticks: { color: "#8b949e" },
                  grid: { color: "#21262d" }
                },
                y: {
                  ticks: { color: "#8b949e" },
                  grid: { color: "#21262d" }
                }
              }
            }
          });
        }
      '))
    )
  )

  return(page)
}

# =============================================================================
# MAIN
# =============================================================================

setwd("~/NFLWork")

cat("=== NFL +EV Dashboard ===\n\n")

# Load bets from duckdb (saved by NFLCombine.R)
cat("Loading bets from database...\n")
con <- dbConnect(duckdb(), dbdir = "Answer Keys/pbp.duckdb", read_only = TRUE)

# Check if table exists
tables <- dbListTables(con)
if (!"nfl_bets_combined" %in% tables) {
  dbDisconnect(con, shutdown = TRUE)
  stop("No bets found. Run 'python run.py nfl' first to generate predictions.")
}

all_bets <- dbGetQuery(con, "SELECT * FROM nfl_bets_combined")
dbDisconnect(con, shutdown = TRUE)

if (nrow(all_bets) == 0) {
  stop("No bets in database. Run 'python run.py nfl' first.")
}

cat(sprintf("Loaded %d bets\n", nrow(all_bets)))

placed_bets <- load_placed_bets(DB_PATH)
cat(sprintf("Found %d placed bets\n", nrow(placed_bets)))

# Load market relationships for correlation detection
relationships <- load_market_relationships()

# Calculate stats
stats <- list(
  total_bets = nrow(all_bets),
  placed_count = nrow(placed_bets),
  avg_ev = mean(all_bets$ev, na.rm = TRUE) * 100,
  max_ev = max(all_bets$ev, na.rm = TRUE) * 100
)

# Create tables
bets_table <- create_bets_table(all_bets, placed_bets, relationships)
placed_table <- create_placed_bets_table(placed_bets)

# Extract unique games for parlay builder (dedupe by label to avoid duplicates)
games <- all_bets %>%
  distinct(home_team, away_team) %>%
  mutate(
    label = paste(away_team, "@", home_team),
    id = paste0(home_team, "_", away_team)  # Create consistent ID from teams
  ) %>%
  select(id, label, home_team, away_team)
games_json <- toJSON(games, auto_unbox = TRUE)

# Extract filter options for bets table
filter_games <- all_bets %>%
  distinct(home_team, away_team) %>%
  mutate(label = paste(away_team, "@", home_team)) %>%
  pull(label) %>%
  sort()

filter_books <- all_bets %>%
  distinct(bookmaker_key) %>%
  pull(bookmaker_key) %>%
  sort()

filter_markets <- all_bets %>%
  mutate(market_type = case_when(
    grepl("^h2h_3", market) ~ "ML 3-Way",
    grepl("^h2h", market) ~ "Moneyline",
    grepl("^spread", market) ~ "Spreads",
    grepl("^total", market) ~ "Totals",
    grepl("team_total", market) ~ "Team Totals",
    grepl("^alternate", market) ~ "Alternates",
    TRUE ~ "Other"
  )) %>%
  distinct(market_type) %>%
  pull(market_type) %>%
  sort()

# Ensure arrays even for single values (use I() to prevent auto_unbox)
filter_options_json <- toJSON(list(
  games = I(filter_games),
  books = I(filter_books),
  markets = I(filter_markets)
), auto_unbox = FALSE)

# Create report
timestamp <- format(Sys.time(), "%b %d, %Y %I:%M %p")
page <- create_report(bets_table, placed_table, stats, timestamp, games_json, filter_options_json)

# Save
save_html(page, OUTPUT_PATH, libdir = "lib")
cat("\nDashboard saved to:", OUTPUT_PATH, "\n")
cat("Open http://localhost:8081 to view\n")
