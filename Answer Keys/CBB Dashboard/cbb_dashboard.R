# CBB +EV Betting Dashboard
# Bets tab with filtering, Kelly sizing, bet placement, and correlation detection

suppressPackageStartupMessages({
  library(tidyverse)
  library(duckdb)
  library(reactable)
  library(htmltools)
  library(htmlwidgets)
  library(jsonlite)
  library(digest)
  library(lubridate)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

# Resolve DASHBOARD_DIR from script location (works in worktrees)
# R encodes spaces as ~+~ in commandArgs --file=, so decode before using
.args <- commandArgs(trailingOnly = FALSE)
.file_arg <- grep("^--file=", .args, value = TRUE)
DASHBOARD_DIR <- if (length(.file_arg) > 0) {
  .raw <- sub("^--file=", "", .file_arg)
  .raw <- gsub("~\\+~", " ", .raw)  # R encodes spaces as ~+~
  normalizePath(dirname(.raw), mustWork = FALSE)
} else {
  normalizePath("~/NFLWork/Answer Keys/CBB Dashboard", mustWork = FALSE)
}
DB_PATH <- file.path(DASHBOARD_DIR, "cbb_dashboard.duckdb")
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
    dbGetQuery(con, "SELECT * FROM placed_bets WHERE status = 'pending' AND (game_time IS NULL OR game_time > NOW())")
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
        odds = placed$odds,
        size = coalesce(placed$actual_size, placed$recommended_size),
        bookmaker = placed$bookmaker,
        strength = 0.90,
        level = "high"
      )))
      next
    }

    # Cross-group correlations (JSON loads as data frame)
    cross_df <- relationships$cross_group_correlations
    if (!is.null(cross_df) && nrow(cross_df) > 0 && !is.na(bet_group) && !is.na(placed_group)) {
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
            odds = placed$odds,
            size = coalesce(placed$actual_size, placed$recommended_size),
            bookmaker = placed$bookmaker,
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

format_market_name <- function(market) {
  market %>%
    str_replace("alternate_", "Alt ") %>%
    str_replace("team_totals", "Team Tot") %>%
    str_replace("totals", "Total") %>%
    str_replace("spreads", "Spread") %>%
    str_replace("h2h", "ML") %>%
    str_replace("_h", " H")
}

escape_tooltip <- function(text) {
  text %>%
    str_replace_all(fixed('"'), "'") %>%
    str_replace_all(fixed("\\"), "")
}

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
      market_display = format_market_name(market)
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
      game_time = ifelse(is.na(pt_start_time), "",
                         format(pt_start_time, "%Y-%m-%dT%H:%M:%SZ")),
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
          market_name <- format_market_name(d$market)
          line_str <- if (!is.null(d$line) && !is.na(d$line)) {
            if (d$line > 0) paste0(" +", d$line) else paste0(" ", d$line)
          } else ""
          odds_str <- if (!is.null(d$odds) && !is.na(d$odds)) {
            if (d$odds > 0) sprintf(" (%+d)", d$odds) else sprintf(" (%d)", d$odds)
          } else ""
          size_str <- if (!is.null(d$size) && !is.na(d$size)) {
            sprintf(" $%.2f", d$size)
          } else ""
          book_str <- if (!is.null(d$bookmaker) && !is.na(d$bookmaker)) {
            sprintf(" @ %s", d$bookmaker)
          } else ""
          sprintf("%s - %s%s%s%s%s", market_name, d$bet_on, line_str, odds_str, size_str, book_str)
        })
        paste("Correlated with:\n", paste(lines, collapse = "\n"))
      }),
      # Simplify market names
      market_display = format_market_name(market)
    ) %>%
    arrange(desc(ev)) %>%
    mutate(warning = "") %>%
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
          corr_attr <- sprintf('data-corr-level="%s"', level)
          if (level == "high") {
            sprintf('<span class="warning-icon high" %s data-tooltip="%s">&#9888;</span>', corr_attr, escape_tooltip(tooltip))
          } else if (level == "medium") {
            sprintf('<span class="warning-icon medium" %s data-tooltip="%s">&#9888;</span>', corr_attr, escape_tooltip(tooltip))
          } else if (level == "low") {
            sprintf('<span class="warning-icon low" %s data-tooltip="%s">&#9675;</span>', corr_attr, escape_tooltip(tooltip))
          } else {
            sprintf('<span %s></span>', corr_attr)
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
        minWidth = 150,
        cell = JS("function(cellInfo) {
          var val = cellInfo.value;
          if (!val || val === '') return '';
          var d = new Date(val);
          if (isNaN(d)) return val;
          var opts = {weekday:'short', month:'2-digit', day:'2-digit', hour:'numeric', minute:'2-digit'};
          return d.toLocaleString(undefined, opts);
        }")
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

create_report <- function(bets_table, placed_table, stats, timestamp, filter_options_json) {
  page <- tagList(
    tags$head(
      tags$meta(charset = "UTF-8"),
      tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
      tags$title("CBB +EV Dashboard"),
      tags$link(
        href = "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap",
        rel = "stylesheet"
      ),
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

        /* Inline bet confirmation */
        .bet-confirm-group {
          display: flex;
          align-items: center;
          gap: 4px;
        }

        .bet-size-input {
          width: 70px;
          background: #0d1117;
          border: 1px solid #58a6ff;
          border-radius: 4px;
          color: #c9d1d9;
          padding: 4px 6px;
          font-size: 0.75rem;
          text-align: right;
          font-family: monospace;
        }

        .bet-size-input:focus {
          outline: none;
          border-color: #79c0ff;
        }

        .btn-confirm {
          background: #238636;
          border: 1px solid #2ea043;
          color: #fff;
          padding: 4px 8px;
          border-radius: 4px;
          font-size: 0.75rem;
          cursor: pointer;
          line-height: 1;
        }

        .btn-confirm:hover {
          background: #2ea043;
        }

        .btn-cancel {
          background: transparent;
          border: 1px solid #30363d;
          color: #8b949e;
          padding: 4px 8px;
          border-radius: 4px;
          font-size: 0.75rem;
          cursor: pointer;
          line-height: 1;
        }

        .btn-cancel:hover {
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
        .warning-icon {
          cursor: help;
          font-size: 1.1em;
          display: inline-block;
        }
        .warning-icon.high { color: #f85149; }
        .warning-icon.medium { color: #d29922; }
        .warning-icon.low { color: #8b949e; font-size: 0.9em; }

        .corr-tooltip {
          position: fixed;
          background: #161b22;
          color: #c9d1d9;
          padding: 8px 12px;
          border-radius: 6px;
          font-size: 0.75rem;
          white-space: pre-line;
          z-index: 9999;
          border: 1px solid #30363d;
          box-shadow: 0 4px 12px rgba(0,0,0,0.4);
          pointer-events: none;
          max-width: 320px;
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
      '))
    ),

    tags$body(
      tags$div(class = "container",
        # Header
        tags$div(class = "header",
          tags$div(
            tags$h1("CBB +EV Dashboard"),
            tags$div(class = "subtitle", paste("Updated", timestamp))
          ),
          tags$button(class = "refresh-btn", onclick = "refreshData()", "Refresh")
        ),

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
          tags$div(class = "filter-group",
            tags$span(class = "filter-label", "Correlation"),
            tags$div(class = "filter-dropdown", id = "filter-correlation-btn", onclick = "toggleFilter('correlation')",
              tags$span(id = "filter-correlation-text", "All Levels")
            ),
            tags$div(class = "filter-menu", id = "filter-correlation-menu")
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

      # JavaScript
      tags$script(HTML('
        // ============ CORRELATION TOOLTIPS ============
        (function() {
          var tip = document.createElement("div");
          tip.className = "corr-tooltip";
          tip.style.display = "none";
          document.body.appendChild(tip);

          document.addEventListener("mouseover", function(e) {
            var icon = e.target.closest(".warning-icon[data-tooltip]");
            if (!icon) return;
            var text = icon.getAttribute("data-tooltip");
            if (!text) return;
            tip.textContent = text;
            tip.style.display = "block";
            var rect = icon.getBoundingClientRect();
            var tipW = tip.offsetWidth;
            var tipH = tip.offsetHeight;
            // Vertical: prefer above, flip below if clipped
            var top = rect.top - tipH - 6;
            if (top < 4) top = rect.bottom + 6;
            // Horizontal: center on icon, clamp to viewport
            var left = rect.left + rect.width / 2 - tipW / 2;
            left = Math.max(4, Math.min(left, window.innerWidth - tipW - 4));
            tip.style.top = top + "px";
            tip.style.left = left + "px";
          });

          document.addEventListener("mouseout", function(e) {
            var icon = e.target.closest(".warning-icon[data-tooltip]");
            if (!icon) return;
            tip.style.display = "none";
          });
        })();

        // ============ BETS FILTERING ============
        const activeFilters = {
          game: new Set(),
          book: new Set(),
          market: new Set(),
          correlation: new Set()
        };

        document.addEventListener("DOMContentLoaded", function() {
          // Fetch persistent book settings before initializing filters
          // Load sizing settings (bankroll/kelly)
          fetch(\'/api/sizing-settings\')
            .then(function(r) { return r.json(); })
            .then(function(sizing) {
              if (sizing.bankroll) document.getElementById("bankroll-input").value = sizing.bankroll;
              if (sizing.kelly_mult) document.getElementById("kelly-input").value = sizing.kelly_mult;
            })
            .catch(function() {});

          fetch(\'/api/book-settings\')
            .then(function(r) { return r.json(); })
            .then(function(settings) {
              window.BOOK_SETTINGS = settings;
              initBetsFilters();
            })
            .catch(function() {
              window.BOOK_SETTINGS = {};
              initBetsFilters();
            });
        });

        function initBetsFilters() {
          const opts = window.FILTER_OPTIONS;
          if (!opts) return;

          populateFilterMenu("game", opts.games);

          // Merge books from data + book_settings so disabled books still appear in the filter
          var allBooks = new Set(opts.books || []);
          Object.keys(window.BOOK_SETTINGS || {}).forEach(function(b) { allBooks.add(b); });
          var allBooksSorted = Array.from(allBooks).sort();
          populateFilterMenu("book", allBooksSorted);

          populateFilterMenu("market", opts.markets);
          populateFilterMenu("correlation", opts.correlations);

          // Auto-discover: register new books as disabled
          var knownBooks = Object.keys(window.BOOK_SETTINGS || {});
          var newBooks = (opts.books || []).filter(function(b) { return knownBooks.indexOf(b) === -1; });
          if (newBooks.length > 0) {
            fetch(\'/api/book-settings/bulk\', {
              method: \'POST\',
              headers: {\'Content-Type\': \'application/json\'},
              body: JSON.stringify({ books: newBooks })
            });
            newBooks.forEach(function(b) { window.BOOK_SETTINGS[b] = false; });
          }

          // Apply saved book settings to checkboxes
          applyBookSettings();

          document.addEventListener("click", function(e) {
            if (!e.target.closest(".filter-group")) {
              document.querySelectorAll(".filter-menu").forEach(m => m.classList.remove("open"));
            }
          });
        }

        function applyBookSettings() {
          var settings = window.BOOK_SETTINGS || {};
          var menu = document.getElementById("filter-book-menu");
          if (!menu) return;

          var checkboxes = menu.querySelectorAll("input[data-val]");
          checkboxes.forEach(function(cb) {
            var val = cb.getAttribute("data-val");
            cb.checked = settings[val] === true;
          });

          // Update Select All state
          var selectAllCb = menu.querySelector(".select-all input");
          if (selectAllCb) {
            var allEnabled = Object.keys(settings).length > 0 &&
              Object.values(settings).every(function(v) { return v === true; });
            selectAllCb.checked = allEnabled;
          }

          // Sync activeFilters and label
          updateFilter("book");
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

          document.querySelectorAll(".filter-menu").forEach(m => m.classList.remove("open"));

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

          // Persist book settings to server
          if (type === "book") {
            checkboxes.forEach(function(cb) {
              var val = cb.getAttribute("data-val");
              var enabled = cb.checked;
              if (window.BOOK_SETTINGS && window.BOOK_SETTINGS[val] !== enabled) {
                window.BOOK_SETTINGS[val] = enabled;
                fetch(\'/api/book-settings\', {
                  method: \'POST\',
                  headers: {\'Content-Type\': \'application/json\'},
                  body: JSON.stringify({ book: val, enabled: enabled })
                });
              }
            });
          }

          const textEl = document.getElementById("filter-" + type + "-text");
          const allLabel = type === "game" ? "All Games" : type === "book" ? "All Books" :
                           type === "correlation" ? "All Levels" : "All Markets";

          if (checkedVals.length === checkboxes.length) {
            textEl.textContent = allLabel;
          } else if (checkedVals.length === 0) {
            textEl.textContent = "None selected";
          } else if (checkedVals.length <= 2) {
            textEl.textContent = checkedVals.join(", ");
          } else {
            textEl.textContent = checkedVals.length + " selected";
          }

          const selectAllCb = menu.querySelector(".select-all input");
          if (selectAllCb) {
            selectAllCb.checked = checkedVals.length === checkboxes.length;
            selectAllCb.indeterminate = checkedVals.length > 0 && checkedVals.length < checkboxes.length;
          }

          applyFilters();
        }

        function applyFilters() {
          const containers = document.querySelectorAll(".table-container:not(.placed-section)");
          const table = containers.length > 0 ? containers[containers.length - 1].querySelector(".rt-tbody") : null;
          if (!table) return;

          const rows = table.querySelectorAll(".rt-tr-group");
          let visibleCount = 0;
          let totalCount = rows.length;

          rows.forEach(row => {
            const cells = row.querySelectorAll(".rt-td");

            let gameText = "";
            let bookText = "";
            let marketText = "";

            cells.forEach((cell, idx) => {
              const text = cell.textContent.trim();
              if (text.includes("@")) {
                gameText = text;
              }
              if (["ML", "Spread", "Total", "Team Tot", "Alt"].some(m => text.includes(m))) {
                marketText = text;
              }
            });

            const bookCell = cells[cells.length - 2];
            if (bookCell) bookText = bookCell.textContent.trim();

            let marketType = "Other";
            if (marketText.includes("ML")) {
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

            // Correlation level from data attribute
            var corrLevel = "None";
            var corrSpan = row.querySelector("[data-corr-level]");
            if (corrSpan) {
              var raw = corrSpan.getAttribute("data-corr-level");
              if (raw === "high") corrLevel = "High";
              else if (raw === "medium") corrLevel = "Medium";
              else if (raw === "low") corrLevel = "Low";
            }

            var gameMatch = activeFilters.game.size === 0 || activeFilters.game.has(gameText);
            var bookMatch = activeFilters.book.size === 0 || activeFilters.book.has(bookText);
            var marketMatch = activeFilters.market.size === 0 || activeFilters.market.has(marketType);
            var corrMatch = activeFilters.correlation.size === 0 || activeFilters.correlation.has(corrLevel);

            var visible = gameMatch && bookMatch && marketMatch && corrMatch;
            row.style.display = visible ? "" : "none";
            if (visible) visibleCount++;
          });

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
          ["game", "market", "correlation"].forEach(function(type) {
            var menu = document.getElementById("filter-" + type + "-menu");
            if (!menu) return;
            var checkboxes = menu.querySelectorAll("input[type=checkbox]");
            var vals = [];

            for (var i = 0; i < checkboxes.length; i++) {
              checkboxes[i].checked = true;
              var val = checkboxes[i].getAttribute("data-val");
              if (val) vals.push(val);
            }

            activeFilters[type] = new Set(vals);
            var textEl = document.getElementById("filter-" + type + "-text");
            var allLabel = type === "game" ? "All Games" : type === "book" ? "All Books" :
                           type === "correlation" ? "All Levels" : "All Markets";
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

          var containers = document.querySelectorAll(".table-container:not(.placed-section)");
          var table = containers.length > 0 ? containers[containers.length - 1] : null;
          if (!table) return;

          var rows = table.querySelectorAll(".rt-tr-group");

          rows.forEach(function(row) {
            var cells = row.querySelectorAll(".rt-td");

            var prob = null;
            var odds = null;
            var sizeCell = null;

            cells.forEach(function(cell, idx) {
              var text = cell.textContent.trim();

              if (/^[+-]\\d+$/.test(text)) {
                odds = parseInt(text);
              }

              if (text.startsWith("$")) {
                sizeCell = cell;
              }
            });

            var btn = row.querySelector("button[data-prob]");
            if (btn) {
              prob = parseFloat(btn.getAttribute("data-prob"));
            }

            if (prob && odds && sizeCell) {
              var newSize = calculateKellyBet(prob, odds, bankroll, kellyMult);
              sizeCell.textContent = "$" + newSize.toFixed(2);

              if (btn) {
                btn.setAttribute("data-size", newSize.toFixed(2));
              }
            }
          });

          // Persist to server
          fetch(\'/api/sizing-settings\', {
            method: \'POST\',
            headers: {\'Content-Type\': \'application/json\'},
            body: JSON.stringify({ bankroll: bankroll, kelly_mult: kellyMult })
          });

          showToast("Bet sizes updated", "success");
        }

        function calculateKellyBet(prob, americanOdds, bankroll, kellyMult) {
          var decimalOdds;
          if (americanOdds > 0) {
            decimalOdds = 1 + (americanOdds / 100);
          } else {
            decimalOdds = 1 + (100 / Math.abs(americanOdds));
          }

          var edge = (prob * decimalOdds) - 1;

          if (edge <= 0) {
            return 0;
          }

          var kellyFraction = edge / (decimalOdds - 1);
          var betSize = bankroll * kellyFraction * kellyMult;

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
          var cell = btn.closest('.rt-td');
          if (!cell) return;

          // Save original HTML for cancel
          cell.setAttribute('data-original', cell.innerHTML);

          var suggestedSize = parseFloat(btn.dataset.size) || 0;

          cell.innerHTML = '<div class="bet-confirm-group">' +
            '<input type="number" class="bet-size-input" value="' + suggestedSize.toFixed(2) + '" step="0.01" min="0.01">' +
            '<button class="btn-confirm" onclick="confirmBet(this)">&#10003;</button>' +
            '<button class="btn-cancel" onclick="cancelBet(this)">&#10005;</button>' +
            '</div>';

          // Focus and select the input
          var input = cell.querySelector('.bet-size-input');
          input.focus();
          input.select();

          // Allow Enter to confirm, Escape to cancel
          input.addEventListener('keydown', function(e) {
            if (e.key === 'Enter') {
              e.preventDefault();
              confirmBet(cell.querySelector('.btn-confirm'));
            } else if (e.key === 'Escape') {
              e.preventDefault();
              cancelBet(cell.querySelector('.btn-cancel'));
            }
          });
        }

        function confirmBet(confirmBtn) {
          var cell = confirmBtn.closest('.rt-td');
          var input = cell.querySelector('.bet-size-input');
          var actualSize = parseFloat(input.value);

          if (!actualSize || actualSize <= 0) {
            showToast("Enter a valid bet amount", "error");
            input.focus();
            return;
          }

          // Parse the original button to extract data attributes
          var tmp = document.createElement('div');
          tmp.innerHTML = cell.getAttribute('data-original');
          var origBtn = tmp.querySelector('button');

          var data = {
            bet_hash: origBtn.dataset.hash,
            game_id: origBtn.dataset.gameId,
            home_team: origBtn.dataset.home,
            away_team: origBtn.dataset.away,
            game_time: origBtn.dataset.time,
            market: origBtn.dataset.market,
            bet_on: origBtn.dataset.betOn,
            line: origBtn.dataset.line === "" ? null : parseFloat(origBtn.dataset.line),
            model_prob: parseFloat(origBtn.dataset.prob),
            model_ev: parseFloat(origBtn.dataset.ev),
            recommended_size: parseFloat(origBtn.dataset.size),
            actual_size: actualSize,
            odds: parseInt(origBtn.dataset.odds),
            bookmaker: origBtn.dataset.book
          };

          // Disable buttons during request
          confirmBtn.disabled = true;
          cell.querySelector('.btn-cancel').disabled = true;

          fetch("/api/place-bet", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify(data)
          })
            .then(function(r) { return r.json(); })
            .then(function(result) {
              if (result.success) {
                // Restore original button HTML, then restyle as Placed
                cell.innerHTML = cell.getAttribute('data-original');
                var btn = cell.querySelector('button');
                btn.className = 'btn-placed';
                btn.textContent = 'Placed';
                btn.onclick = function() { removeBet(this); };
                showToast("Bet placed: $" + actualSize.toFixed(2), "success");
              } else {
                showToast(result.error, "error");
                confirmBtn.disabled = false;
                cell.querySelector('.btn-cancel').disabled = false;
              }
            })
            .catch(function() {
              showToast("Server error", "error");
              confirmBtn.disabled = false;
              cell.querySelector('.btn-cancel').disabled = false;
            });
        }

        function cancelBet(cancelBtn) {
          var cell = cancelBtn.closest('.rt-td');
          var original = cell.getAttribute('data-original');
          if (original) {
            cell.innerHTML = original;
          }
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

        function escapeHtml(text) {
          const div = document.createElement("div");
          div.textContent = text;
          return div.innerHTML;
        }
      '))
    )
  )

  return(page)
}

# =============================================================================
# MAIN
# =============================================================================

setwd(dirname(dirname(DASHBOARD_DIR)))  # NFLWork root, derived from DASHBOARD_DIR

cat("=== CBB +EV Dashboard ===\n\n")

# Load bets from duckdb (saved by CBBCombine.R)
cat("Loading bets from database...\n")
con <- dbConnect(duckdb(), dbdir = "Answer Keys/cbb.duckdb", read_only = TRUE)

# Check if table exists
tables <- dbListTables(con)
if (!"cbb_bets_combined" %in% tables) {
  all_bets <- tibble()
} else {
  all_bets <- dbGetQuery(con, "SELECT * FROM cbb_bets_combined") %>%
    filter(is.na(pt_start_time) | pt_start_time > Sys.time())
}
dbDisconnect(con, shutdown = TRUE)

cat(sprintf("Loaded %d bets\n", nrow(all_bets)))

placed_bets <- load_placed_bets(DB_PATH)
cat(sprintf("Found %d placed bets\n", nrow(placed_bets)))

# Load market relationships for correlation detection
relationships <- load_market_relationships()

# Calculate stats
stats <- list(
  total_bets = nrow(all_bets),
  placed_count = nrow(placed_bets),
  avg_ev = if (nrow(all_bets) > 0) mean(all_bets$ev, na.rm = TRUE) * 100 else 0,
  max_ev = if (nrow(all_bets) > 0) max(all_bets$ev, na.rm = TRUE) * 100 else 0
)

# Create tables
if (nrow(all_bets) > 0) {
  bets_table <- create_bets_table(all_bets, placed_bets, relationships)
} else {
  bets_table <- tags$div(
    style = "text-align: center; padding: 48px; color: #8b949e;",
    tags$p(style = "font-size: 1.1rem;", "No +EV bets found for upcoming games."),
    tags$p(style = "font-size: 0.85rem;", "Try clicking Refresh when more games are available.")
  )
}
placed_table <- create_placed_bets_table(placed_bets)

# Extract filter options for bets table
filter_games <- if (nrow(all_bets) > 0) {
  all_bets %>%
    distinct(home_team, away_team) %>%
    mutate(label = paste(away_team, "@", home_team)) %>%
    pull(label) %>%
    sort()
} else character()

filter_books <- if (nrow(all_bets) > 0) {
  all_bets %>%
    distinct(bookmaker_key) %>%
    pull(bookmaker_key) %>%
    sort()
} else character()

filter_markets <- if (nrow(all_bets) > 0) {
  all_bets %>%
    mutate(market_type = case_when(
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
} else character()

# Ensure arrays even for single values
filter_options_json <- toJSON(list(
  games = I(filter_games),
  books = I(filter_books),
  markets = I(filter_markets),
  correlations = I(c("None", "Low", "Medium", "High"))
), auto_unbox = FALSE)

# Create report
timestamp <- format(Sys.time(), "%b %d, %Y %I:%M %p")
page <- create_report(bets_table, placed_table, stats, timestamp, filter_options_json)

# Save
save_html(page, OUTPUT_PATH, libdir = "lib")
cat("\nDashboard saved to:", OUTPUT_PATH, "\n")
cat("Open http://localhost:8082 to view\n")
