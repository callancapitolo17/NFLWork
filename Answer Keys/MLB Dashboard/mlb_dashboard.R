# MLB +EV Betting Dashboard
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
  normalizePath("~/NFLWork/Answer Keys/MLB Dashboard", mustWork = FALSE)
}
DB_PATH <- file.path(DASHBOARD_DIR, "mlb_dashboard.duckdb")
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
  digest(paste(game_id, market, bet_on, line_str, sep = "|"), algo = "sha256", serialize = FALSE)
}

# =============================================================================
# SAME-GAME DETECTION
# =============================================================================

find_same_game_bets <- function(row_idx, all_bets, placed_bets) {
  game_id <- all_bets$id[row_idx]
  details <- list()

  # Other pipeline bets on this game (from the table)
  other_idx <- which(all_bets$id == game_id & seq_len(nrow(all_bets)) != row_idx)
  for (j in other_idx) {
    details <- append(details, list(list(
      market = all_bets$market[j],
      bet_on = all_bets$bet_on[j],
      line = all_bets$line[j],
      odds = all_bets$odds[j],
      bet_size = all_bets$bet_size[j],
      bookmaker = all_bets$bookmaker_key[j],
      is_placed = FALSE
    )))
  }

  # Placed bets on this game not already in pipeline output
  if (nrow(placed_bets) > 0) {
    same_game_placed <- placed_bets[placed_bets$game_id == game_id, ]
    for (k in seq_len(nrow(same_game_placed))) {
      p <- same_game_placed[k, ]
      match_idx <- which(all_bets$id == game_id &
                         all_bets$market == p$market &
                         all_bets$bet_on == p$bet_on &
                         all_bets$line == p$line)
      if (length(match_idx) > 0) {
        # Already in table — mark as placed
        for (d in seq_along(details)) {
          if (details[[d]]$market == p$market &&
              details[[d]]$bet_on == p$bet_on &&
              identical(details[[d]]$line, p$line)) {
            details[[d]]$is_placed <- TRUE
            details[[d]]$actual_size <- p$actual_size
            details[[d]]$recommended_size <- p$recommended_size
          }
        }
      } else {
        # Placed bet NOT in pipeline — add as extra entry
        details <- append(details, list(list(
          market = p$market,
          bet_on = p$bet_on,
          line = p$line,
          odds = p$odds,
          bet_size = p$recommended_size,
          bookmaker = p$bookmaker,
          is_placed = TRUE,
          actual_size = p$actual_size,
          recommended_size = p$recommended_size
        )))
      }
    }
  }

  list(has_same_game = length(details) > 0, details = details)
}

# =============================================================================
# TABLE CREATION
# =============================================================================

format_market_name <- function(market) {
  market %>%
    str_replace("_1st_5_innings", " F5") %>%
    str_replace("alternate_", "Alt ") %>%
    str_replace("team_totals", "Team Tot") %>%
    str_replace("totals", "Total") %>%
    str_replace("spreads", "Spread") %>%
    str_replace("h2h", "ML")
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
      actual_display = sprintf("$%.0f", coalesce(actual_size, recommended_size)),
      rec_display = sprintf("$%.0f", recommended_size),
      sizes_differ = !is.na(actual_size) & abs(actual_size - recommended_size) > 0.005,
      line_display = case_when(
        is.na(line) ~ "",
        line > 0 ~ paste0("+", line),
        TRUE ~ as.character(line)
      ),
      market_display = format_market_name(market)
    ) %>%
    select(bet_hash, game, market_display, bet_on, line_display, ev_display, odds_display,
           actual_display, rec_display, sizes_differ, bookmaker)

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
      actual_display = colDef(name = "Actual", minWidth = 70, align = "right",
        cell = function(value, index) {
          if (table_data$sizes_differ[index]) {
            span(class = "size-actual", value)
          } else {
            value
          }
        }
      ),
      rec_display = colDef(name = "Rec", minWidth = 70, align = "right",
        style = list(color = "#8b949e", fontSize = "0.8rem")
      ),
      sizes_differ = colDef(show = FALSE),
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

# =============================================================================
# PARLAY TABLES
# =============================================================================

create_placed_parlays_table <- function(placed_parlays) {
  if (nrow(placed_parlays) == 0) return(NULL)

  table_data <- placed_parlays %>%
    mutate(
      game = paste(away_team, "@", home_team),
      spread_team = ifelse(grepl("Home", combo), home_team, away_team),
      spread_fmt  = ifelse(spread_line > 0, paste0("+", spread_line), as.character(spread_line)),
      ou_prefix   = ifelse(grepl("Over", combo), "O", "U"),
      legs_display = paste0(spread_team, " ", spread_fmt, " \u00b7 ", ou_prefix, total_line),
      odds_display = ifelse(wz_odds > 0, paste0("+", wz_odds), as.character(wz_odds)),
      edge_display = sprintf("+%.1f%%", edge_pct),
      actual_display = sprintf("$%.0f", coalesce(actual_size, kelly_bet)),
      rec_display = sprintf("$%.0f", kelly_bet)
    ) %>%
    select(parlay_hash, game, legs_display, odds_display, edge_display, actual_display, rec_display)

  reactable(
    table_data,
    compact = TRUE,
    columns = list(
      parlay_hash = colDef(show = FALSE),
      game = colDef(name = "Game", minWidth = 180),
      legs_display = colDef(name = "Legs", minWidth = 220),
      odds_display = colDef(name = "WZ Odds", minWidth = 70, align = "right",
        style = list(fontFamily = "monospace")),
      edge_display = colDef(name = "Edge", minWidth = 65, align = "right",
        style = list(color = "#3fb950", fontWeight = "600")),
      actual_display = colDef(name = "Actual", minWidth = 65, align = "right",
        style = list(fontWeight = "600")),
      rec_display = colDef(name = "Rec", minWidth = 65, align = "right",
        style = list(color = "#8b949e", fontSize = "0.8rem"))
    ),
    theme = reactableTheme(
      color = "#c9d1d9",
      backgroundColor = "#21262d",
      borderColor = "#30363d",
      headerStyle = list(backgroundColor = "#161b22", color = "#8b949e", fontWeight = "600", fontSize = "0.7rem")
    )
  )
}

create_parlays_table <- function(parlay_opps, placed_parlays) {
  placed_hashes <- if (nrow(placed_parlays) > 0) placed_parlays$parlay_hash else character()

  # Ensure exact-payout columns exist even on older rows (pre-Stage-2 backfill).
  # parlay_pricer.py --exact-payouts populates these; fall back to Kelly+wz_dec
  # arithmetic when they're NA.
  if (!"exact_wager" %in% names(parlay_opps)) parlay_opps$exact_wager <- NA_integer_
  if (!"exact_to_win" %in% names(parlay_opps)) parlay_opps$exact_to_win <- NA_integer_

  table_data <- parlay_opps %>%
    mutate(
      is_placed = parlay_hash %in% placed_hashes,
      spread_team    = ifelse(grepl("Home", combo), home_team, away_team),
      spread_fmt     = ifelse(spread_line > 0, paste0("+", spread_line), as.character(spread_line)),
      sp_price_fmt   = case_when(
        is.na(spread_price) ~ "?",
        spread_price > 0    ~ paste0("+", spread_price),
        TRUE                ~ as.character(spread_price)
      ),
      ou_prefix      = ifelse(grepl("Over", combo), "O", "U"),
      tot_price_fmt  = case_when(
        is.na(total_price) ~ "?",
        total_price > 0    ~ paste0("+", total_price),
        TRUE               ~ as.character(total_price)
      ),
      legs_display   = paste0(spread_team, " ", spread_fmt, " (", sp_price_fmt, ")",
                              " \u00b7 ", ou_prefix, total_line, " (", tot_price_fmt, ")"),
      fair_display   = ifelse(fair_odds > 0, paste0("+", fair_odds), as.character(fair_odds)),
      wz_display     = ifelse(wz_odds > 0, paste0("+", wz_odds), as.character(wz_odds)),
      edge_display   = sprintf("+%.1f%%", edge_pct),
      # Prefer empirical stake/payout from --exact-payouts. Fall back to Kelly-ideal
      # wager and arithmetic payout for rows the exact-payout step couldn't price
      # (e.g. idgm missing or WZ rejected all candidate stakes).
      size_display   = sprintf("$%.0f", coalesce(as.numeric(exact_wager), kelly_bet)),
      to_win_display = sprintf("$%.0f", coalesce(as.numeric(exact_to_win),
                                                  round(kelly_bet * (wz_dec - 1)))),
      corr_display   = sprintf("%.3f", corr_factor)
    ) %>%
    arrange(desc(edge_pct))

  reactable(
    table_data,
    compact = TRUE,
    striped = TRUE,
    highlight = TRUE,
    defaultPageSize = 25,
    columns = list(
      # Hidden columns for JS data attributes
      parlay_hash = colDef(show = FALSE),
      game_id = colDef(show = FALSE),
      home_team = colDef(show = FALSE),
      away_team = colDef(show = FALSE),
      spread_line = colDef(show = FALSE),
      total_line = colDef(show = FALSE),
      spread_price = colDef(show = FALSE),
      total_price = colDef(show = FALSE),
      fair_odds = colDef(show = FALSE),
      wz_odds = colDef(show = FALSE),
      fair_dec = colDef(show = FALSE),
      wz_dec = colDef(show = FALSE),
      corr_factor = colDef(show = FALSE),
      edge_pct = colDef(show = FALSE),
      kelly_bet = colDef(show = FALSE),
      # Backend-only columns from parlay_pricer.py --exact-payouts and the
      # wagerzon scraper's idgm. They drive size_display / to_win_display
      # (and the Place button's data-size) but shouldn't render as table columns.
      idgm = colDef(show = FALSE),
      exact_wager = colDef(show = FALSE),
      exact_to_win = colDef(show = FALSE),
      model_prob_raw = colDef(show = FALSE),
      blended_prob_raw = colDef(show = FALSE),
      model_prob_pct = colDef(show = FALSE),
      n_samples = colDef(show = FALSE),
      combo = colDef(show = FALSE),
      spread_team = colDef(show = FALSE),
      spread_fmt = colDef(show = FALSE),
      sp_price_fmt = colDef(show = FALSE),
      ou_prefix = colDef(show = FALSE),
      tot_price_fmt = colDef(show = FALSE),

      # Visible columns
      game = colDef(name = "Game", minWidth = 180),
      game_time = colDef(
        name = "Time",
        minWidth = 130,
        cell = JS("function(cellInfo) {
          var val = cellInfo.value;
          if (!val || val === '') return '';
          var d = new Date(val);
          if (isNaN(d)) return val;
          var opts = {weekday:'short', month:'2-digit', day:'2-digit', hour:'numeric', minute:'2-digit'};
          return d.toLocaleString(undefined, opts);
        }")
      ),
      legs_display = colDef(name = "Legs", minWidth = 260),
      fair_display = colDef(name = "Fair", minWidth = 70, align = "right",
        style = list(fontFamily = "monospace", color = "#8b949e")),
      # Per-book devigged fair probabilities (NA when the book didn't price the combo).
      # Muted grey — they're diagnostic detail; the primary "Fair" col is the blend.
      dk_fair_prob = colDef(name = "DK", minWidth = 60, align = "right",
        format = colFormat(percent = TRUE, digits = 1),
        style = list(fontFamily = "monospace", color = "#8b949e")),
      fd_fair_prob = colDef(name = "FD", minWidth = 60, align = "right",
        format = colFormat(percent = TRUE, digits = 1),
        style = list(fontFamily = "monospace", color = "#8b949e")),
      px_fair_prob = colDef(name = "PX", minWidth = 60, align = "right",
        format = colFormat(percent = TRUE, digits = 1),
        style = list(fontFamily = "monospace", color = "#8b949e")),
      nv_fair_prob = colDef(name = "NV", minWidth = 60, align = "right",
        format = colFormat(percent = TRUE, digits = 1),
        style = list(fontFamily = "monospace", color = "#8b949e")),
      wz_display = colDef(name = "WZ", minWidth = 70, align = "right",
        style = list(fontFamily = "monospace")),
      corr_display = colDef(name = "Corr", minWidth = 65, align = "right",
        style = list(color = "#8b949e")),
      edge_display = colDef(name = "Edge %", minWidth = 70, align = "right",
        cell = function(value, index) {
          ep <- table_data$edge_pct[index]
          color <- if (ep >= 15) "#3fb950" else if (ep >= 10) "#56d364" else if (ep >= 5) "#7ee787" else "#a5d6a7"
          span(style = list(color = color, fontWeight = "600"), value)
        }
      ),
      size_display = colDef(name = "Size", minWidth = 65, align = "right"),
      to_win_display = colDef(name = "To Win", minWidth = 65, align = "right",
        style = list(color = "#3fb950")),
      is_placed = colDef(
        name = "Action",
        minWidth = 90,
        align = "center",
        filterable = FALSE,
        html = TRUE,
        cell = function(value, index) {
          row <- table_data[index, ]
          # data-size drives the "Place" action's bet amount. Prefer the nudged
          # exact_wager from Stage 2; fall back to Kelly-ideal if Stage 2
          # couldn't price this row.
          data_size <- if (!is.na(row$exact_wager)) row$exact_wager else row$kelly_bet
          data_attrs <- sprintf(
            'data-hash="%s" data-game-id="%s" data-home="%s" data-away="%s" data-time="%s" data-combo="%s" data-spread="%s" data-total="%s" data-fair-odds="%s" data-wz-odds="%s" data-edge="%s" data-size="%s"',
            row$parlay_hash, row$game_id, row$home_team, row$away_team,
            ifelse(is.na(row$game_time), "", as.character(row$game_time)),
            row$combo, row$spread_line, row$total_line,
            row$fair_odds, row$wz_odds, row$edge_pct, data_size
          )
          if (value) {
            sprintf('<button class="btn-placed" onclick="removeParlay(this)" %s>Placed</button>', data_attrs)
          } else {
            sprintf('<button class="btn-place" onclick="placeParlay(this)" %s>Place</button>', data_attrs)
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
        textTransform = "uppercase"
      )
    )
  )
}

create_bets_table <- function(all_bets, placed_bets) {
  placed_hashes <- if (nrow(placed_bets) > 0) placed_bets$bet_hash else character()
  # Build lookups for placed bet actual_size and recommended_size by hash
  placed_actual <- setNames(
    if (nrow(placed_bets) > 0) placed_bets$actual_size else numeric(),
    if (nrow(placed_bets) > 0) placed_bets$bet_hash else character()
  )
  placed_recommended <- setNames(
    if (nrow(placed_bets) > 0) placed_bets$recommended_size else numeric(),
    if (nrow(placed_bets) > 0) placed_bets$bet_hash else character()
  )

  # Find same-game bets for each bet
  same_game_info <- lapply(seq_len(nrow(all_bets)), function(i) {
    find_same_game_bets(i, all_bets, placed_bets)
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
      odds_display = {
        base <- ifelse(odds > 0, paste0("+", odds), as.character(odds))
        # For Kalshi, show effective cents (fee-adjusted) — use precise value from scraper when available
        eff_cents <- if ("cents" %in% names(.)) cents else NA_real_
        implied_pct <- ifelse(!is.na(eff_cents),
          sprintf("%.2f", eff_cents),
          sprintf("%.2f", ifelse(odds < 0, -odds / (-odds + 100), 100 / (odds + 100)) * 100))
        ifelse(bookmaker_key == "kalshi",
          paste0(base, " (", implied_pct, "\u00A2)"),
          base)
      },
      size_display = sprintf("$%.0f", bet_size),
      # Fill status: compare actual vs current bet_size (not stale placed_rec)
      placed_actual = ifelse(is_placed, placed_actual[bet_hash], NA_real_),
      placed_rec = ifelse(is_placed, placed_recommended[bet_hash], NA_real_),
      fill_status = case_when(
        !is_placed ~ "not_placed",
        is.na(placed_actual) ~ "placed",
        round(placed_actual) >= round(bet_size) ~ "placed",
        TRUE ~ "partial"
      ),
      fill_diff = ifelse(fill_status == "partial", round(bet_size) - round(placed_actual), NA_real_),
      # Same-game indicators
      has_correlation = sapply(same_game_info, function(x) x$has_same_game),
      correlation_level = sapply(same_game_info, function(x) if (x$has_same_game) "same_game" else "none"),
      correlation_tooltip = {
        dot <- intToUtf8(0xB7)     # middle dot separator
        chk <- intToUtf8(0x2713)   # checkmark for placed
        cir <- intToUtf8(0x2013)   # en-dash for unplaced
        sapply(same_game_info, function(x) {
          if (!x$has_same_game) return("")
          details <- x$details
          lines <- sapply(details, function(d) {
            market_name <- format_market_name(d$market)
            line_str <- if (!is.null(d$line) && !is.na(d$line)) {
              if (d$line > 0) paste0(" +", d$line) else paste0(" ", d$line)
            } else ""
            odds_str <- if (!is.null(d$odds) && !is.na(d$odds)) {
              if (d$odds > 0) sprintf(" (%+d)", d$odds) else sprintf(" (%d)", d$odds)
            } else ""
            if (isTRUE(d$is_placed)) {
              act <- d$actual_size
              rec <- d$recommended_size
              size_str <- if (!is.null(act) && !is.na(act)) {
                rec_part <- if (!is.null(rec) && !is.na(rec) && abs(act - rec) > 0.01) {
                  sprintf(" (rec $%.0f)", rec)
                } else ""
                sprintf("$%.0f%s", act, rec_part)
              } else ""
              prefix <- chk
            } else {
              size_str <- if (!is.null(d$bet_size) && !is.na(d$bet_size)) {
                sprintf("$%.0f", d$bet_size)
              } else ""
              prefix <- cir
            }
            book_str <- if (!is.null(d$bookmaker) && !is.na(d$bookmaker)) d$bookmaker else ""
            paste(prefix, paste(Filter(nzchar, c(
              paste0(market_name, " - ", d$bet_on, line_str, odds_str),
              size_str, book_str
            )), collapse = paste0(" ", dot, " ")))
          })
          paste(lines, collapse = "\n")
        })
      },
      # Simplify market names
      market_display = format_market_name(market)
    ) %>%
    arrange(desc(ev)) %>%
    mutate(warning = "") %>%
    select(
      bet_hash, id, warning, game, game_time, market, market_display, bet_on, line, line_display,
      ev_pct, ev_display, odds, odds_display, bet_size, size_display,
      bookmaker_key, is_placed, fill_status, fill_diff, placed_actual, placed_rec,
      home_team, away_team, pt_start_time, prob,
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
      fill_status = colDef(show = FALSE),
      fill_diff = colDef(show = FALSE),
      placed_actual = colDef(show = FALSE),
      placed_rec = colDef(show = FALSE),
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
          if (level == "same_game") {
            sprintf('<span class="warning-icon same-game" %s data-tooltip="%s">&#9679;</span>', corr_attr, escape_tooltip(tooltip))
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
        align = "right",
        html = TRUE,
        cell = function(value, index) {
          sz <- table_data$bet_size[index]
          sprintf('<span data-bet-size="%.2f">%s</span>', sz, value)
        }
      ),
      bookmaker_key = colDef(
        name = "Book",
        minWidth = 100,
        filterable = TRUE
      ),
      is_placed = colDef(
        name = "Action",
        minWidth = 110,
        align = "center",
        filterable = FALSE,
        html = TRUE,
        cell = function(value, index) {
          row <- table_data[index, ]
          data_attrs <- sprintf(
            'data-hash="%s" data-game-id="%s" data-home="%s" data-away="%s" data-time="%s" data-market="%s" data-bet-on="%s" data-line="%s" data-prob="%s" data-ev="%s" data-size="%s" data-odds="%s" data-book="%s" data-actual="%s" data-fill-status="%s"',
            row$bet_hash, row$id, row$home_team, row$away_team,
            as.character(row$pt_start_time), row$market, row$bet_on,
            ifelse(is.na(row$line), "", row$line),
            row$prob, row$ev_pct / 100, row$bet_size, row$odds, row$bookmaker_key,
            ifelse(is.na(row$placed_actual), "", row$placed_actual),
            row$fill_status
          )

          status <- row$fill_status
          if (status == "partial") {
            diff_label <- sprintf("-$%.0f", row$fill_diff)
            sprintf('<button class="btn-partial" onclick="updateBet(this)" %s>Partial %s</button>', data_attrs, diff_label)
          } else if (status == "placed") {
            sprintf('<button class="btn-placed" onclick="removeBet(this)" %s>Placed</button>', data_attrs)
          } else {
            auto_books <- c("wagerzon", "hoop88", "bfa")
            place_btn <- sprintf('<button class="btn-place" onclick="placeBet(this)" %s>Place</button>', data_attrs)
            if (row$bookmaker_key %in% auto_books) {
              auto_btn <- sprintf('<button class="btn-auto" onclick="autoPlaceBet(this)" %s>Auto</button>', data_attrs)
              paste0(place_btn, auto_btn)
            } else {
              place_btn
            }
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

create_report <- function(bets_table, placed_table, stats, timestamp, filter_options_json,
                          parlays_table = NULL, placed_parlays_table = NULL, parlay_opps = tibble(),
                          parlay_filter_options_json = "{}") {
  page <- tagList(
    tags$head(
      tags$meta(charset = "UTF-8"),
      tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
      tags$title("MLB +EV Dashboard"),
      tags$link(
        href = "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap",
        rel = "stylesheet"
      ),
      # Inject filter options for bets table
      tags$script(HTML(sprintf("window.FILTER_OPTIONS = %s;", filter_options_json))),
      # Inject filter options for parlay table
      tags$script(HTML(sprintf("window.PARLAY_FILTER_OPTIONS = %s;", parlay_filter_options_json))),
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

        .tab-bar {
          display: flex;
          gap: 4px;
          margin-bottom: 24px;
          border-bottom: 1px solid #21262d;
          padding-bottom: 0;
        }

        .tab-btn {
          background: transparent;
          border: none;
          border-bottom: 2px solid transparent;
          color: #8b949e;
          padding: 10px 20px;
          font-size: 0.9rem;
          font-weight: 500;
          cursor: pointer;
          transition: all 0.15s;
        }

        .tab-btn:hover {
          color: #c9d1d9;
        }

        .tab-btn.active {
          color: #f0f6fc;
          border-bottom-color: #238636;
        }

        .placed-parlays-live td {
          padding: 8px 12px;
          color: #c9d1d9;
          border-bottom: 1px solid #30363d;
          font-size: 0.85rem;
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

        .btn-partial {
          background: transparent;
          border: 1px solid #d29922;
          color: #d29922;
          padding: 5px 12px;
          border-radius: 6px;
          font-size: 0.75rem;
          font-weight: 500;
          cursor: pointer;
          transition: all 0.15s;
        }

        .btn-partial:hover {
          background: #d29922;
          color: #fff;
        }

        .btn-auto {
          background: transparent;
          border: 1px solid #1f6feb;
          color: #58a6ff;
          padding: 5px 8px;
          border-radius: 6px;
          font-size: 0.7rem;
          font-weight: 500;
          cursor: pointer;
          transition: all 0.15s;
          margin-left: 4px;
        }

        .btn-auto:hover {
          background: #1f6feb;
          color: #fff;
        }

        .btn-navigating {
          background: transparent;
          border: 1px solid #1f6feb;
          color: #58a6ff;
          padding: 5px 8px;
          border-radius: 6px;
          font-size: 0.7rem;
          font-weight: 500;
          cursor: default;
          animation: pulse-border 1.5s ease-in-out infinite;
        }

        .btn-nav-ready {
          background: transparent;
          border: 1px solid #3fb950;
          color: #3fb950;
          padding: 5px 8px;
          border-radius: 6px;
          font-size: 0.7rem;
          font-weight: 500;
          cursor: default;
          animation: pulse-border 1.5s ease-in-out infinite;
        }

        .btn-nav-error {
          background: transparent;
          border: 1px solid #f85149;
          color: #f85149;
          padding: 5px 8px;
          border-radius: 6px;
          font-size: 0.7rem;
          font-weight: 500;
          cursor: pointer;
        }

        @keyframes pulse-border {
          0%, 100% { opacity: 1; }
          50% { opacity: 0.5; }
        }

        /* Place-bet modal */
        .modal-overlay {
          position: fixed;
          top: 0; left: 0; right: 0; bottom: 0;
          background: rgba(0, 0, 0, 0.6);
          z-index: 999;
          display: flex;
          align-items: center;
          justify-content: center;
        }

        .modal-box {
          background: #161b22;
          border: 1px solid #30363d;
          border-radius: 8px;
          padding: 24px;
          min-width: 320px;
          max-width: 400px;
          box-shadow: 0 8px 24px rgba(0, 0, 0, 0.5);
        }

        .modal-title {
          font-size: 1rem;
          font-weight: 600;
          color: #f0f6fc;
          margin: 0 0 16px 0;
        }

        .modal-detail {
          font-size: 0.8rem;
          color: #8b949e;
          margin-bottom: 4px;
        }

        .modal-detail span {
          color: #c9d1d9;
          font-weight: 500;
        }

        .modal-recommended {
          font-size: 0.85rem;
          color: #3fb950;
          font-weight: 600;
          margin: 12px 0 16px 0;
        }

        .modal-input-group {
          margin-bottom: 20px;
        }

        .modal-input-label {
          display: block;
          font-size: 0.75rem;
          color: #8b949e;
          text-transform: uppercase;
          letter-spacing: 0.5px;
          margin-bottom: 6px;
        }

        .modal-input {
          width: 100%;
          background: #0d1117;
          border: 1px solid #30363d;
          border-radius: 6px;
          color: #c9d1d9;
          padding: 10px 12px;
          font-size: 1rem;
          font-family: monospace;
        }

        .modal-input:focus {
          outline: none;
          border-color: #58a6ff;
        }

        .modal-actions {
          display: flex;
          gap: 8px;
          justify-content: flex-end;
        }

        .modal-btn-confirm {
          background: #238636;
          border: 1px solid #2ea043;
          color: #fff;
          padding: 8px 20px;
          border-radius: 6px;
          font-size: 0.85rem;
          font-weight: 500;
          cursor: pointer;
        }

        .modal-btn-confirm:hover { background: #2ea043; }

        .modal-btn-cancel {
          background: transparent;
          border: 1px solid #30363d;
          color: #8b949e;
          padding: 8px 16px;
          border-radius: 6px;
          font-size: 0.85rem;
          cursor: pointer;
        }

        .modal-btn-cancel:hover {
          border-color: #f85149;
          color: #f85149;
        }

        /* Size display in placed bets */
        .size-actual {
          font-weight: 600;
          color: #c9d1d9;
        }

        .size-recommended {
          font-size: 0.75rem;
          color: #8b949e;
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
        .warning-icon.same-game { color: #58a6ff; }

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
            tags$h1("MLB Answer Key Dashboard"),
            tags$div(class = "subtitle", paste("Updated", timestamp))
          ),
          tags$button(class = "refresh-btn", onclick = "refreshData()", "Refresh")
        ),

        # Tab Navigation
        tags$div(class = "tab-bar",
          tags$button(class = "tab-btn active", id = "tab-btn-bets", onclick = "switchTab(\'bets\')", "Bets"),
          tags$button(class = "tab-btn", id = "tab-btn-parlays", onclick = "switchTab(\'parlays\')", "Parlays")
        ),

        # ============ BETS TAB ============
        tags$div(id = "tab-bets", class = "tab-content",

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
            tags$span(class = "filter-label", "Same Game"),
            tags$div(class = "filter-dropdown", id = "filter-correlation-btn", onclick = "toggleFilter('correlation')",
              tags$span(id = "filter-correlation-text", "All")
            ),
            tags$div(class = "filter-menu", id = "filter-correlation-menu")
          ),
          tags$div(class = "filter-group",
            tags$span(class = "filter-label", "Min Size"),
            tags$div(style = "display: flex; align-items: center; background: #161b22; border: 1px solid #30363d; border-radius: 6px; padding: 8px 12px; gap: 4px;",
              tags$span(style = "color: #8b949e; font-size: 0.85rem;", "$"),
              tags$input(id = "filter-size-input", type = "number", min = "0", value = "0",
                style = "width: 50px; padding: 0; background: transparent; border: none; color: #c9d1d9; font-size: 0.85rem; outline: none;",
                onchange = "updateSizeFilter()", oninput = "updateSizeFilter()")
            )
          ),
          tags$div(class = "filter-group",
            tags$span(class = "filter-label", "Status"),
            tags$div(class = "filter-dropdown", id = "filter-status-btn", onclick = "toggleFilter('status')",
              tags$span(id = "filter-status-text", "All Statuses")
            ),
            tags$div(class = "filter-menu", id = "filter-status-menu")
          ),
          tags$button(class = "clear-filters-btn", onclick = "clearAllFilters()", "Clear Filters")
        ),

        # Available Bets Table
        tags$div(class = "section-header",
          tags$span("Available Bets"),
          tags$span(id = "filtered-count", style = "margin-left: 8px; color: #8b949e; font-size: 0.8rem;")
        ),
        tags$div(class = "table-container", id = "bets-table-container", bets_table)
        ), # end tab-bets

        # ============ PARLAYS TAB ============
        tags$div(id = "tab-parlays", class = "tab-content", style = "display: none;",

          # Placed Parlays (always present, JS can append rows)
          tags$div(class = "section-header", id = "placed-parlays-header",
            style = if (is.null(placed_parlays_table)) "display:none;" else "",
            "Placed Parlays"),
          tags$div(class = "table-container placed-section", id = "placed-parlays-section",
            style = if (is.null(placed_parlays_table)) "display:none;" else "",
            if (!is.null(placed_parlays_table)) placed_parlays_table,
            # Plain HTML table for live-injected rows
            tags$table(id = "placed-parlays-live", class = "placed-parlays-live",
              style = "width: 100%; border-collapse: collapse;"
            )
          ),

          # Parlay Sizing Controls
          tags$div(class = "sizing-controls",
            tags$div(class = "sizing-group",
              tags$span(class = "sizing-label", "Parlay Bankroll ($)"),
              tags$input(id = "parlay-bankroll-input", class = "sizing-input", type = "number",
                value = "100", min = "1", step = "10")
            ),
            tags$div(class = "sizing-group",
              tags$span(class = "sizing-label", "Parlay Kelly"),
              tags$input(id = "parlay-kelly-input", class = "sizing-input", type = "number",
                value = "0.25", min = "0.01", max = "1", step = "0.05")
            ),
            tags$div(class = "sizing-group",
              tags$span(class = "sizing-label", "Min Edge (%)"),
              tags$input(id = "parlay-min-edge-input", class = "sizing-input", type = "number",
                value = "0", min = "0", step = "1",
                onchange = "filterParlaysByEdge()", oninput = "filterParlaysByEdge()")
            ),
            tags$button(class = "apply-sizing-btn", onclick = "applyParlaySizing()", "Apply")
          ),

          # Parlay Filter Bar
          tags$div(class = "filter-bar",
            tags$div(class = "filter-group",
              tags$span(class = "filter-label", "Game"),
              tags$div(class = "filter-dropdown", id = "parlay-filter-game-btn",
                onclick = "toggleParlayFilter('game')",
                tags$span(id = "parlay-filter-game-text", "All Games")
              ),
              tags$div(class = "filter-menu", id = "parlay-filter-game-menu")
            ),
            tags$div(class = "filter-group",
              tags$span(class = "filter-label", "Status"),
              tags$div(class = "filter-dropdown", id = "parlay-filter-status-btn",
                onclick = "toggleParlayFilter('status')",
                tags$span(id = "parlay-filter-status-text", "All Statuses")
              ),
              tags$div(class = "filter-menu", id = "parlay-filter-status-menu")
            ),
            tags$div(class = "filter-group",
              tags$span(class = "filter-label", "Min Size"),
              tags$div(
                style = "display: flex; align-items: center; background: #161b22; border: 1px solid #30363d; border-radius: 6px; padding: 8px 12px; gap: 4px;",
                tags$span(style = "color: #8b949e; font-size: 0.85rem;", "$"),
                tags$input(
                  id = "parlay-filter-size-input", type = "number", min = "0", value = "0",
                  style = "width: 50px; padding: 0; background: transparent; border: none; color: #c9d1d9; font-size: 0.85rem; outline: none;",
                  onchange = "updateParlayMinSize()", oninput = "updateParlayMinSize()"
                )
              )
            ),
            tags$button(class = "clear-filters-btn", onclick = "clearParlayFilters()", "Clear Filters")
          ),

          # Parlay Opportunities
          if (!is.null(parlays_table)) {
            tagList(
              tags$div(class = "section-header",
                tags$span("Parlay Opportunities"),
                tags$span(id = "parlay-filtered-count",
                  style = "margin-left: 8px; color: #8b949e; font-size: 0.8rem;",
                  sprintf("(%d)", nrow(parlay_opps)))
              ),
              tags$div(class = "table-container", id = "parlays-table-container", parlays_table)
            )
          },

          # Empty state when no parlays
          if (is.null(parlays_table)) {
            tags$div(style = "text-align: center; padding: 48px; color: #8b949e;",
              tags$p(style = "font-size: 1.1rem;", "No parlay opportunities found."),
              tags$p(style = "font-size: 0.85rem;", "Try clicking Refresh when games are available.")
            )
          }
        ) # end tab-parlays
      ),

      # JavaScript
      tags$script(HTML('
        // ============ TAB SWITCHING ============
        function switchTab(tab) {
          document.getElementById("tab-bets").style.display = tab === "bets" ? "" : "none";
          document.getElementById("tab-parlays").style.display = tab === "parlays" ? "" : "none";
          document.getElementById("tab-btn-bets").className = "tab-btn" + (tab === "bets" ? " active" : "");
          document.getElementById("tab-btn-parlays").className = "tab-btn" + (tab === "parlays" ? " active" : "");
          // Apply parlay edge filter when switching to parlays tab
          if (tab === "parlays") filterParlaysByEdge();
        }

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

        // ============ LIVE CORRELATION RECALCULATION ============
        var _corrRecalcRunning = false;

        // Track bets placed/updated/removed during this session so we can
        // re-apply their visual state after reactable re-renders (React
        // rebuilds the DOM from its initial state, wiping JS modifications).
        var _sessionPlaced = {};

        function reapplyPlacedStates() {
          var table = document.getElementById("bets-table-container");
          if (!table) return;
          table.querySelectorAll("button[data-hash]").forEach(function(btn) {
            var hash = btn.dataset.hash;
            var state = _sessionPlaced[hash];
            if (!state) return;
            btn.className = state.className;
            btn.textContent = state.text;
            btn.setAttribute("data-fill-status", state.fillStatus);
            btn.dataset.actual = state.actual || "";
            if (state.action === "remove") {
              btn.onclick = function() { removeBet(this); };
            } else if (state.action === "update") {
              btn.onclick = function() { updateBet(this); };
            } else if (state.action === "place") {
              btn.onclick = function() { placeBet(this); };
            }
          });
        }
        function formatMarketNameJS(market) {
          return market
            .replace("_1st_5_innings", " F5")
            .replace("alternate_", "Alt ")
            .replace("team_totals", "Team Tot")
            .replace("totals", "Total")
            .replace("spreads", "Spread")
            .replace("h2h", "ML");
        }

        function recalcSameGame(gameId) {
          _corrRecalcRunning = true;
          var table = document.getElementById("bets-table-container");
          if (!table) { _corrRecalcRunning = false; return; }

          var allRows = table.querySelectorAll(".rt-tr-group");
          var gameRows = [];

          allRows.forEach(function(row) {
            var btn = row.querySelector("button[data-game-id]");
            if (!btn || btn.dataset.gameId !== gameId) return;
            gameRows.push({ row: row, btn: btn });
          });

          gameRows.forEach(function(item) {
            var btn = item.btn;
            var row = item.row;
            var myHash = btn.dataset.hash;
            var span = row.querySelector("[data-corr-level]");
            if (!span) return;

            var others = [];
            gameRows.forEach(function(other) {
              if (other.btn.dataset.hash === myHash) return;
              var ob = other.btn;
              var fs = ob.getAttribute("data-fill-status");
              var isPlaced = (fs === "placed" || fs === "partial");
              others.push({
                market: ob.dataset.market, betOn: ob.dataset.betOn,
                line: ob.dataset.line, odds: ob.dataset.odds,
                size: ob.dataset.size, actual: ob.dataset.actual,
                book: ob.dataset.book, isPlaced: isPlaced
              });
            });

            if (others.length === 0) {
              span.setAttribute("data-corr-level", "none");
              span.className = "";
              span.textContent = "";
              span.removeAttribute("data-tooltip");
            } else {
              span.setAttribute("data-corr-level", "same_game");
              span.className = "warning-icon same-game";
              span.innerHTML = "&#9679;";

              var chk = String.fromCharCode(0x2713);
              var cir = String.fromCharCode(0x2013);
              var dot = String.fromCharCode(0xB7);
              var lines = others.map(function(d) {
                var mName = formatMarketNameJS(d.market);
                var lineStr = (d.line && d.line !== "" && d.line !== "NA") ?
                  (parseFloat(d.line) > 0 ? " +" + d.line : " " + d.line) : "";
                var oddsStr = (d.odds && d.odds !== "") ?
                  (parseInt(d.odds) > 0 ? " (+" + d.odds + ")" : " (" + d.odds + ")") : "";
                var sizeStr = "";
                var prefix = cir;
                if (d.isPlaced) {
                  prefix = chk;
                  var act = parseFloat(d.actual);
                  var rec = parseFloat(d.size);
                  if (!isNaN(act) && act > 0) {
                    var recPart = (!isNaN(rec) && Math.abs(act - rec) > 0.01) ?
                      " (rec $" + Math.round(rec) + ")" : "";
                    sizeStr = "$" + Math.round(act) + recPart;
                  }
                } else {
                  var sz = parseFloat(d.size);
                  if (!isNaN(sz)) sizeStr = "$" + Math.round(sz);
                }
                var parts = [mName + " - " + d.betOn + lineStr + oddsStr];
                if (sizeStr) parts.push(sizeStr);
                if (d.book) parts.push(d.book);
                return prefix + " " + parts.join(" " + dot + " ");
              });
              span.setAttribute("data-tooltip", lines.join("\\n"));
            }
          });

          if (typeof applyFilters === "function") applyFilters();
          setTimeout(function() { _corrRecalcRunning = false; }, 150);
        }

        // Re-run same-game indicators after reactable pagination re-renders
        (function() {
          var corrTimer = null;
          function scheduleRecalc() {
            if (_corrRecalcRunning) return;
            if (corrTimer) clearTimeout(corrTimer);
            corrTimer = setTimeout(function() {
              reapplyPlacedStates();
              var table = document.getElementById("bets-table-container");
              if (table) {
                var ids = new Set();
                table.querySelectorAll("button[data-game-id]").forEach(function(btn) { ids.add(btn.dataset.gameId); });
                ids.forEach(function(gid) { recalcSameGame(gid); });
              }
            }, 100);
          }
          document.addEventListener("DOMContentLoaded", function() {
            var table = document.getElementById("bets-table-container");
            if (!table) return;
            var observer = new MutationObserver(scheduleRecalc);
            observer.observe(table, { childList: true, subtree: true });
          });
        })();

        // ============ BETS FILTERING ============
        const activeFilters = {
          game: new Set(),
          book: new Set(),
          market: new Set(),
          correlation: new Set(),
          status: new Set()
        };

        // ============ PARLAY FILTERING ============
        const activeParlayFilters = {
          game: new Set(),
          status: new Set(),
          minEdge: 0,
          minSize: 0
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

          // Load book settings and filter settings in parallel, init filters when both ready
          var bookReady = fetch(\'/api/book-settings\')
            .then(function(r) { return r.json(); })
            .then(function(settings) { window.BOOK_SETTINGS = settings; })
            .catch(function() { window.BOOK_SETTINGS = {}; });

          var filterReady = fetch(\'/api/filter-settings\')
            .then(function(r) { return r.json(); })
            .then(function(fs) { window.FILTER_SETTINGS = fs; })
            .catch(function() { window.FILTER_SETTINGS = {}; });

          Promise.all([bookReady, filterReady]).then(function() {
            initBetsFilters();
            initParlayFilters();
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
          populateFilterMenu("status", opts.statuses);

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

          // Apply saved filter settings (market, correlation, status)
          applyFilterSettings();

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

        function applyFilterSettings() {
          var settings = window.FILTER_SETTINGS || {};
          ["market", "correlation", "status"].forEach(function(type) {
            if (!settings[type]) return;  // no saved state = keep "all checked" default
            var saved = new Set(settings[type]);
            var menu = document.getElementById("filter-" + type + "-menu");
            if (!menu) return;
            var checkboxes = menu.querySelectorAll("input[data-val]");
            checkboxes.forEach(function(cb) {
              cb.checked = saved.has(cb.getAttribute("data-val"));
            });
            var selectAllCb = menu.querySelector(".select-all input");
            if (selectAllCb) {
              selectAllCb.checked = saved.size === checkboxes.length;
              selectAllCb.indeterminate = saved.size > 0 && saved.size < checkboxes.length;
            }
            updateFilter(type);
          });
          // Restore numeric size filter
          if (settings.size && Array.isArray(settings.size) && settings.size.length > 0) {
            var savedSize = parseFloat(settings.size[0]) || 0;
            document.getElementById("filter-size-input").value = savedSize;
          }
          applyFilters();
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

          // Persist filter settings for market/correlation/status
          if (["market", "correlation", "status"].indexOf(type) !== -1) {
            fetch(\'/api/filter-settings\', {
              method: \'POST\',
              headers: {\'Content-Type\': \'application/json\'},
              body: JSON.stringify({ filter_type: type, selected_values: checkedVals })
            });
          }

          const textEl = document.getElementById("filter-" + type + "-text");
          const allLabel = type === "game" ? "All Games" : type === "book" ? "All Books" :
                           type === "correlation" ? "All" :
                           type === "status" ? "All Statuses" : "All Markets";

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

        function updateSizeFilter() {
          var val = parseFloat(document.getElementById("filter-size-input").value) || 0;
          fetch(\'/api/filter-settings\', {
            method: \'POST\',
            headers: {\'Content-Type\': \'application/json\'},
            body: JSON.stringify({ filter_type: "size", selected_values: [val] })
          });
          applyFilters();
        }

        function applyFilters() {
          const betsContainer = document.getElementById("bets-table-container");
          const table = betsContainer ? betsContainer.querySelector(".rt-tbody") : null;
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

            // Same-game level from data attribute
            var corrLevel = "Standalone";
            var corrSpan = row.querySelector("[data-corr-level]");
            if (corrSpan && corrSpan.getAttribute("data-corr-level") === "same_game") {
              corrLevel = "Same Game";
            }

            // Size from data attribute (numeric)
            var betSize = 0;
            var sizeSpan = row.querySelector("[data-bet-size]");
            if (sizeSpan) betSize = parseFloat(sizeSpan.getAttribute("data-bet-size")) || 0;

            // Status from Action button data-fill-status
            var statusLabel = "Not Placed";
            var statusBtn = row.querySelector("button[data-fill-status]");
            if (statusBtn) {
              var fs = statusBtn.getAttribute("data-fill-status");
              if (fs === "placed") statusLabel = "Placed";
              else if (fs === "partial") statusLabel = "Partial Fill";
            }

            var gameMatch = activeFilters.game.size === 0 || activeFilters.game.has(gameText);
            var bookMatch = activeFilters.book.size === 0 || activeFilters.book.has(bookText);
            var marketMatch = activeFilters.market.size === 0 || activeFilters.market.has(marketType);
            var corrMatch = activeFilters.correlation.size === 0 || activeFilters.correlation.has(corrLevel);
            var minSize = parseFloat(document.getElementById("filter-size-input").value) || 0;
            var sizeMatch = betSize >= minSize;
            var statusMatch = activeFilters.status.size === 0 || activeFilters.status.has(statusLabel);

            var visible = gameMatch && bookMatch && marketMatch && corrMatch && sizeMatch && statusMatch;
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
          ["game", "market", "correlation", "status"].forEach(function(type) {
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
                           type === "correlation" ? "All" :
                           type === "status" ? "All Statuses" : "All Markets";
            textEl.textContent = allLabel;

            // Persist "all" state for persistable filters
            if (["market", "correlation", "status"].indexOf(type) !== -1) {
              fetch(\'/api/filter-settings\', {
                method: \'POST\',
                headers: {\'Content-Type\': \'application/json\'},
                body: JSON.stringify({ filter_type: type, selected_values: vals })
              });
            }
          });
          // Reset numeric size filter
          document.getElementById("filter-size-input").value = "0";
          fetch(\'/api/filter-settings\', {
            method: \'POST\',
            headers: {\'Content-Type\': \'application/json\'},
            body: JSON.stringify({ filter_type: "size", selected_values: [0] })
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

          var betsContainer = document.getElementById("bets-table-container");
          if (!betsContainer) return;

          var rows = betsContainer.querySelectorAll(".rt-tr-group");

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
              var sizeSpan = sizeCell.querySelector("[data-bet-size]");
              if (sizeSpan) {
                sizeSpan.textContent = "$" + newSize.toFixed(0);
                sizeSpan.setAttribute("data-bet-size", newSize.toFixed(2));
              } else {
                sizeCell.textContent = "$" + newSize.toFixed(0);
              }

              if (btn) {
                btn.setAttribute("data-size", newSize.toFixed(0));
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
          var recommended = parseFloat(btn.dataset.size) || 0;
          var betOn = btn.dataset.betOn || \'\';
          var odds = btn.dataset.odds || \'\';
          var book = btn.dataset.book || \'\';
          var oddsDisplay = (parseInt(odds) > 0 ? \'+\' : \'\') + odds;

          // Build modal
          var overlay = document.createElement(\'div\');
          overlay.className = \'modal-overlay\';
          overlay.innerHTML =
            \'<div class="modal-box">\' +
              \'<div class="modal-title">Place Bet</div>\' +
              \'<div class="modal-detail">Pick: <span>\' + escapeHtml(betOn) + \' \' + escapeHtml(oddsDisplay) + \'</span></div>\' +
              \'<div class="modal-detail">Book: <span>\' + escapeHtml(book) + \'</span></div>\' +
              \'<div class="modal-recommended">Recommended: $\' + recommended.toFixed(0) + \'</div>\' +
              \'<div class="modal-input-group">\' +
                \'<label class="modal-input-label">Actual Amount ($)</label>\' +
                \'<input type="number" class="modal-input" value="\' + recommended.toFixed(0) + \'" step="1" min="1">\' +
              \'</div>\' +
              \'<div class="modal-actions">\' +
                \'<button class="modal-btn-cancel">Cancel</button>\' +
                \'<button class="modal-btn-confirm">Confirm</button>\' +
              \'</div>\' +
            \'</div>\';

          document.body.appendChild(overlay);

          var input = overlay.querySelector(\'.modal-input\');
          input.focus();
          input.select();

          // Confirm handler
          function doConfirm() {
            var actualSize = parseFloat(input.value);
            if (!actualSize || actualSize <= 0) {
              showToast("Enter a valid bet amount", "error");
              input.focus();
              return;
            }

            var data = {
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
              recommended_size: recommended,
              actual_size: actualSize,
              odds: parseInt(btn.dataset.odds),
              bookmaker: btn.dataset.book
            };

            var confirmBtn = overlay.querySelector(\'.modal-btn-confirm\');
            confirmBtn.disabled = true;
            confirmBtn.textContent = \'Placing...\';

            fetch("/api/place-bet", {
              method: "POST",
              headers: { "Content-Type": "application/json" },
              body: JSON.stringify(data)
            })
              .then(function(r) { return r.json(); })
              .then(function(result) {
                if (result.success) {
                  overlay.remove();
                  btn.dataset.actual = actualSize;
                  if (Math.round(actualSize) >= Math.round(recommended)) {
                    btn.className = \'btn-placed\';
                    btn.textContent = \'Placed\';
                    btn.setAttribute(\'data-fill-status\', \'placed\');
                    btn.onclick = function() { removeBet(this); };
                    _sessionPlaced[btn.dataset.hash] = { className: \'btn-placed\', text: \'Placed\', fillStatus: \'placed\', actual: actualSize, action: \'remove\' };
                  } else {
                    var diff = Math.round(recommended) - Math.round(actualSize);
                    btn.className = \'btn-partial\';
                    btn.textContent = \'Partial -$\' + diff.toFixed(0);
                    btn.setAttribute(\'data-fill-status\', \'partial\');
                    btn.onclick = function() { updateBet(this); };
                    _sessionPlaced[btn.dataset.hash] = { className: \'btn-partial\', text: \'Partial -$\' + diff.toFixed(0), fillStatus: \'partial\', actual: actualSize, action: \'update\' };
                  }
                  showToast("Bet placed: $" + actualSize.toFixed(0), "success");
                  recalcSameGame(btn.dataset.gameId);
                  applyFilters();
                } else {
                  showToast(result.error, "error");
                  confirmBtn.disabled = false;
                  confirmBtn.textContent = \'Confirm\';
                }
              })
              .catch(function() {
                showToast("Server error", "error");
                confirmBtn.disabled = false;
                confirmBtn.textContent = \'Confirm\';
              });
          }

          overlay.querySelector(\'.modal-btn-confirm\').addEventListener(\'click\', doConfirm);
          overlay.querySelector(\'.modal-btn-cancel\').addEventListener(\'click\', function() { overlay.remove(); });
          overlay.addEventListener(\'click\', function(e) { if (e.target === overlay) overlay.remove(); });
          input.addEventListener(\'keydown\', function(e) {
            if (e.key === \'Enter\') { e.preventDefault(); doConfirm(); }
            else if (e.key === \'Escape\') { e.preventDefault(); overlay.remove(); }
          });
        }

        function autoPlaceBet(autoBtn) {
          // Find the Place button (sibling) to get data attributes
          var btn = autoBtn.previousElementSibling;
          if (!btn || !btn.dataset.hash) btn = autoBtn; // fallback if no sibling

          var recommended = parseFloat(btn.dataset.size) || 0;
          var betOn = btn.dataset.betOn || \'\';
          var odds = btn.dataset.odds || \'\';
          var book = btn.dataset.book || \'\';
          var oddsDisplay = (parseInt(odds) > 0 ? \'+\' : \'\') + odds;

          var overlay = document.createElement(\'div\');
          overlay.className = \'modal-overlay\';
          overlay.innerHTML =
            \'<div class="modal-box">\' +
              \'<div class="modal-title">Auto-Place Bet</div>\' +
              \'<div class="modal-detail">Pick: <span>\' + escapeHtml(betOn) + \' \' + escapeHtml(oddsDisplay) + \'</span></div>\' +
              \'<div class="modal-detail">Book: <span>\' + escapeHtml(book) + \'</span></div>\' +
              \'<div class="modal-recommended">Recommended: $\' + recommended.toFixed(0) + \'</div>\' +
              \'<div class="modal-input-group">\' +
                \'<label class="modal-input-label">Actual Amount ($)</label>\' +
                \'<input type="number" class="modal-input" value="\' + recommended.toFixed(0) + \'" step="1" min="1">\' +
              \'</div>\' +
              \'<div class="modal-actions">\' +
                \'<button class="modal-btn-cancel">Cancel</button>\' +
                \'<button class="modal-btn-confirm">Launch Browser</button>\' +
              \'</div>\' +
            \'</div>\';

          document.body.appendChild(overlay);
          var input = overlay.querySelector(\'.modal-input\');
          input.focus();
          input.select();

          function doConfirm() {
            var actualSize = parseFloat(input.value);
            if (!actualSize || actualSize <= 0) {
              showToast("Enter a valid bet amount", "error");
              input.focus();
              return;
            }

            var data = {
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
              recommended_size: recommended,
              actual_size: actualSize,
              odds: parseInt(btn.dataset.odds),
              bookmaker: btn.dataset.book
            };

            var confirmBtn = overlay.querySelector(\'.modal-btn-confirm\');
            confirmBtn.disabled = true;
            confirmBtn.textContent = \'Launching...\';

            // Step 1: Save bet as pending
            fetch("/api/place-bet", {
              method: "POST",
              headers: { "Content-Type": "application/json" },
              body: JSON.stringify(data)
            })
            .then(function(r) { return r.json(); })
            .then(function(result) {
              if (!result.success) {
                showToast(result.error, "error");
                confirmBtn.disabled = false;
                confirmBtn.textContent = \'Launch Browser\';
                return;
              }

              // Step 2: Launch navigator
              return fetch("/api/auto-place", {
                method: "POST",
                headers: { "Content-Type": "application/json" },
                body: JSON.stringify(data)
              })
              .then(function(r2) { return r2.json(); })
              .then(function(navResult) {
                overlay.remove();

                // Update Place button to Placed state
                btn.className = \'btn-placed\';
                btn.textContent = \'Placed\';
                btn.setAttribute(\'data-fill-status\', \'placed\');
                btn.onclick = function() { removeBet(this); };
                btn.dataset.actual = actualSize;
                _sessionPlaced[btn.dataset.hash] = { className: \'btn-placed\', text: \'Placed\', fillStatus: \'placed\', actual: actualSize, action: \'remove\' };

                // Update Auto button to show navigating status
                autoBtn.className = \'btn-navigating\';
                autoBtn.textContent = \'Navigating...\';
                autoBtn.onclick = null;

                showToast("Bet placed, browser launching for " + book, "success");
                recalcSameGame(btn.dataset.gameId);
                applyFilters();

                // Step 3: Poll navigator status
                pollNavStatus(btn.dataset.hash, autoBtn);
              });
            })
            .catch(function() {
              showToast("Server error", "error");
              confirmBtn.disabled = false;
              confirmBtn.textContent = \'Launch Browser\';
            });
          }

          overlay.querySelector(\'.modal-btn-confirm\').addEventListener(\'click\', doConfirm);
          overlay.querySelector(\'.modal-btn-cancel\').addEventListener(\'click\', function() { overlay.remove(); });
          overlay.addEventListener(\'click\', function(e) { if (e.target === overlay) overlay.remove(); });
          input.addEventListener(\'keydown\', function(e) {
            if (e.key === \'Enter\') { e.preventDefault(); doConfirm(); }
            else if (e.key === \'Escape\') { e.preventDefault(); overlay.remove(); }
          });
        }

        function pollNavStatus(betHash, statusBtn) {
          var pollInterval = setInterval(function() {
            fetch("/api/nav-status/" + betHash)
              .then(function(r) { return r.json(); })
              .then(function(data) {
                if (data.status === \'navigating\') {
                  statusBtn.className = \'btn-navigating\';
                  statusBtn.textContent = \'Navigating...\';
                } else if (data.status === \'ready_to_confirm\') {
                  clearInterval(pollInterval);
                  statusBtn.className = \'btn-nav-ready\';
                  statusBtn.textContent = \'Ready\';
                } else if (data.status === \'nav_error\') {
                  clearInterval(pollInterval);
                  statusBtn.className = \'btn-nav-error\';
                  statusBtn.textContent = \'Retry\';
                  statusBtn.onclick = function() { autoPlaceBet(statusBtn); };
                  showToast("Navigator failed — click Retry", "error");
                } else if (data.status === \'nav_timeout\') {
                  clearInterval(pollInterval);
                  statusBtn.className = \'btn-nav-error\';
                  statusBtn.textContent = \'Timeout\';
                  statusBtn.onclick = function() { autoPlaceBet(statusBtn); };
                  showToast("Navigator timed out", "error");
                } else if (data.status === \'pending\') {
                  // Already confirmed or placed normally, stop polling
                  clearInterval(pollInterval);
                  statusBtn.style.display = \'none\';
                }
              })
              .catch(function() {
                // Network error — keep polling
              });
          }, 3000);
        }

        function updateBet(btn) {
          var recommended = parseFloat(btn.dataset.size) || 0;
          var currentActual = parseFloat(btn.dataset.actual) || 0;
          var betOn = btn.dataset.betOn || \'\';
          var odds = btn.dataset.odds || \'\';
          var book = btn.dataset.book || \'\';
          var oddsDisplay = (parseInt(odds) > 0 ? \'+\' : \'\') + odds;

          var overlay = document.createElement(\'div\');
          overlay.className = \'modal-overlay\';
          overlay.innerHTML =
            \'<div class="modal-box">\' +
              \'<div class="modal-title">Update Bet Amount</div>\' +
              \'<div class="modal-detail">Pick: <span>\' + escapeHtml(betOn) + \' \' + escapeHtml(oddsDisplay) + \'</span></div>\' +
              \'<div class="modal-detail">Book: <span>\' + escapeHtml(book) + \'</span></div>\' +
              \'<div class="modal-recommended">Recommended: $\' + recommended.toFixed(0) + \'</div>\' +
              \'<div class="modal-detail">Currently filled: <span>$\' + currentActual.toFixed(0) + \'</span></div>\' +
              \'<div class="modal-input-group">\' +
                \'<label class="modal-input-label">New Amount ($)</label>\' +
                \'<input type="number" class="modal-input" value="\' + currentActual.toFixed(0) + \'" step="1" min="1">\' +
              \'</div>\' +
              \'<div class="modal-actions">\' +
                \'<button class="modal-btn-cancel">Cancel</button>\' +
                \'<button class="modal-btn-confirm">Update</button>\' +
              \'</div>\' +
            \'</div>\';

          document.body.appendChild(overlay);
          var input = overlay.querySelector(\'.modal-input\');
          input.focus();
          input.select();

          function doUpdate() {
            var newAmount = parseFloat(input.value);
            if (!newAmount || newAmount <= 0) {
              showToast("Enter a valid amount", "error");
              input.focus();
              return;
            }

            var confirmBtn = overlay.querySelector(\'.modal-btn-confirm\');
            confirmBtn.disabled = true;
            confirmBtn.textContent = \'Updating...\';

            fetch("/api/update-bet", {
              method: "POST",
              headers: { "Content-Type": "application/json" },
              body: JSON.stringify({ bet_hash: btn.dataset.hash, actual_size: newAmount })
            })
              .then(function(r) { return r.json(); })
              .then(function(result) {
                if (result.success) {
                  overlay.remove();
                  btn.dataset.actual = newAmount;
                  if (Math.round(newAmount) >= Math.round(recommended)) {
                    btn.className = \'btn-placed\';
                    btn.textContent = \'Placed\';
                    btn.setAttribute(\'data-fill-status\', \'placed\');
                    btn.onclick = function() { removeBet(this); };
                    _sessionPlaced[btn.dataset.hash] = { className: \'btn-placed\', text: \'Placed\', fillStatus: \'placed\', actual: newAmount, action: \'remove\' };
                  } else {
                    var diff = Math.round(recommended) - Math.round(newAmount);
                    btn.textContent = \'Partial -$\' + diff.toFixed(0);
                    btn.setAttribute(\'data-fill-status\', \'partial\');
                    _sessionPlaced[btn.dataset.hash] = { className: \'btn-partial\', text: \'Partial -$\' + diff.toFixed(0), fillStatus: \'partial\', actual: newAmount, action: \'update\' };
                  }
                  showToast("Updated to $" + newAmount.toFixed(0), "success");
                  recalcSameGame(btn.dataset.gameId);
                } else {
                  showToast(result.error, "error");
                  confirmBtn.disabled = false;
                  confirmBtn.textContent = \'Update\';
                }
              })
              .catch(function() {
                showToast("Server error", "error");
                confirmBtn.disabled = false;
                confirmBtn.textContent = \'Update\';
              });
          }

          overlay.querySelector(\'.modal-btn-confirm\').addEventListener(\'click\', doUpdate);
          overlay.querySelector(\'.modal-btn-cancel\').addEventListener(\'click\', function() { overlay.remove(); });
          overlay.addEventListener(\'click\', function(e) { if (e.target === overlay) overlay.remove(); });
          input.addEventListener(\'keydown\', function(e) {
            if (e.key === \'Enter\') { e.preventDefault(); doUpdate(); }
            else if (e.key === \'Escape\') { e.preventDefault(); overlay.remove(); }
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
                btn.setAttribute("data-fill-status", "not_placed");
                btn.dataset.actual = "";
                btn.onclick = function() { placeBet(this); };
                _sessionPlaced[btn.dataset.hash] = { className: "btn-place", text: "Place", fillStatus: "not_placed", actual: "", action: "place" };
                showToast("Removed", "success");
                recalcSameGame(btn.dataset.gameId);
              }
            });
        }

        // ============ PARLAY FILTERING ============
        function populateParlayFilterMenu(type, options) {
          const menu = document.getElementById("parlay-filter-" + type + "-menu");
          if (!menu || !options || !Array.isArray(options)) return;

          let html = "<label class=\\"filter-option select-all\\">" +
            "<input type=\\"checkbox\\" checked onchange=\\"toggleParlaySelectAll(\'" + type + "\', this.checked)\\" /> " +
            "Select All</label>";

          for (var i = 0; i < options.length; i++) {
            var opt = options[i];
            var safeVal = String(opt).replace(/"/g, "&quot;");
            html += "<label class=\\"filter-option\\">" +
              "<input type=\\"checkbox\\" checked data-val=\\"" + safeVal + "\\" onchange=\\"updateParlayFilter(\'" + type + "\')\\"/> " +
              escapeHtml(String(opt)) + "</label>";
          }
          menu.innerHTML = html;
        }

        function toggleParlayFilter(type) {
          const menu = document.getElementById("parlay-filter-" + type + "-menu");
          const wasOpen = menu.classList.contains("open");
          document.querySelectorAll(".filter-menu").forEach(m => m.classList.remove("open"));
          if (!wasOpen) menu.classList.add("open");
        }

        function toggleParlaySelectAll(type, isChecked) {
          const menu = document.getElementById("parlay-filter-" + type + "-menu");
          menu.querySelectorAll("input[type=checkbox]").forEach(function(cb) {
            cb.checked = isChecked;
          });
          updateParlayFilter(type);
        }

        function updateParlayFilter(type) {
          const menu = document.getElementById("parlay-filter-" + type + "-menu");
          const checkboxes = menu.querySelectorAll("input[data-val]");
          const checkedVals = [];
          for (var i = 0; i < checkboxes.length; i++) {
            if (checkboxes[i].checked) checkedVals.push(checkboxes[i].getAttribute("data-val"));
          }
          activeParlayFilters[type] = new Set(checkedVals);

          const textEl = document.getElementById("parlay-filter-" + type + "-text");
          const allLabel = type === "game" ? "All Games" : "All Statuses";
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

          // Persist status filter (game not persisted — changes daily)
          if (type === "status") {
            fetch(\'/api/filter-settings\', {
              method: \'POST\',
              headers: {\'Content-Type\': \'application/json\'},
              body: JSON.stringify({ filter_type: "parlay_status", selected_values: checkedVals })
            });
          }

          applyParlayFilters();
        }

        function updateParlayMinSize() {
          var val = parseFloat(document.getElementById("parlay-filter-size-input").value) || 0;
          activeParlayFilters.minSize = val;
          fetch(\'/api/filter-settings\', {
            method: \'POST\',
            headers: {\'Content-Type\': \'application/json\'},
            body: JSON.stringify({ filter_type: "parlay_size", selected_values: [val] })
          });
          applyParlayFilters();
        }

        function initParlayFilters() {
          const opts = window.PARLAY_FILTER_OPTIONS;
          if (!opts) return;
          populateParlayFilterMenu("game", opts.games);
          populateParlayFilterMenu("status", opts.statuses || ["Not Placed", "Placed"]);

          // Restore persisted state
          var settings = window.FILTER_SETTINGS || {};
          if (settings.parlay_status) {
            var saved = new Set(settings.parlay_status);
            var menu = document.getElementById("parlay-filter-status-menu");
            if (menu) {
              menu.querySelectorAll("input[data-val]").forEach(function(cb) {
                cb.checked = saved.has(cb.getAttribute("data-val"));
              });
              var selectAllCb = menu.querySelector(".select-all input");
              if (selectAllCb) {
                selectAllCb.checked = saved.size === menu.querySelectorAll("input[data-val]").length;
                selectAllCb.indeterminate = saved.size > 0 && saved.size < menu.querySelectorAll("input[data-val]").length;
              }
              updateParlayFilter("status");
            }
          }
          if (settings.parlay_size && Array.isArray(settings.parlay_size) && settings.parlay_size.length > 0) {
            var savedSize = parseFloat(settings.parlay_size[0]) || 0;
            document.getElementById("parlay-filter-size-input").value = savedSize;
            activeParlayFilters.minSize = savedSize;
          }
          applyParlayFilters();
        }

        function applyParlayFilters() {
          var container = document.getElementById("parlays-table-container");
          if (!container) return;
          var rows = container.querySelectorAll(".rt-tr-group");
          var minEdge = activeParlayFilters.minEdge;
          var minSize = activeParlayFilters.minSize;
          var visible = 0;

          rows.forEach(function(row) {
            var btn = row.querySelector("button[data-edge]");
            var edge = btn ? parseFloat(btn.dataset.edge) : 0;
            var size = btn ? parseFloat(btn.dataset.size) : 0;
            var gameText = btn ? (btn.dataset.away + " @ " + btn.dataset.home) : "";
            var statusLabel = btn && btn.classList.contains("btn-placed") ? "Placed" : "Not Placed";

            var gameMatch  = activeParlayFilters.game.size === 0   || activeParlayFilters.game.has(gameText);
            var statusMatch = activeParlayFilters.status.size === 0 || activeParlayFilters.status.has(statusLabel);
            var edgeMatch  = edge >= minEdge;
            var sizeMatch  = size >= minSize;

            var show = gameMatch && statusMatch && edgeMatch && sizeMatch;
            row.style.display = show ? "" : "none";
            if (show) visible++;
          });

          var countEl = document.getElementById("parlay-filtered-count");
          if (countEl) {
            if (visible === rows.length) {
              countEl.textContent = "(" + rows.length + ")";
            } else {
              countEl.textContent = "(" + visible + " of " + rows.length + ")";
            }
          }
        }

        function clearParlayFilters() {
          ["game", "status"].forEach(function(type) {
            var menu = document.getElementById("parlay-filter-" + type + "-menu");
            if (!menu) return;
            var checkboxes = menu.querySelectorAll("input[type=checkbox]");
            var vals = [];
            for (var i = 0; i < checkboxes.length; i++) {
              checkboxes[i].checked = true;
              var val = checkboxes[i].getAttribute("data-val");
              if (val) vals.push(val);
            }
            activeParlayFilters[type] = new Set(vals);
            var textEl = document.getElementById("parlay-filter-" + type + "-text");
            textEl.textContent = type === "game" ? "All Games" : "All Statuses";
          });
          document.getElementById("parlay-filter-size-input").value = "0";
          document.getElementById("parlay-min-edge-input").value = "0";
          activeParlayFilters.minSize = 0;
          activeParlayFilters.minEdge = 0;
          // Persist resets
          fetch(\'/api/filter-settings\', {
            method: \'POST\',
            headers: {\'Content-Type\': \'application/json\'},
            body: JSON.stringify({ filter_type: "parlay_status", selected_values: ["Not Placed", "Placed"] })
          });
          fetch(\'/api/filter-settings\', {
            method: \'POST\',
            headers: {\'Content-Type\': \'application/json\'},
            body: JSON.stringify({ filter_type: "parlay_size", selected_values: [0] })
          });
          applyParlayFilters();
        }

        // ============ PARLAY EDGE FILTER ============
        function filterParlaysByEdge() {
          // Sync edge state from DOM before filtering
          var el = document.getElementById("parlay-min-edge-input");
          if (el) activeParlayFilters.minEdge = parseFloat(el.value) || 0;
          applyParlayFilters();
        }

        // ============ PARLAY SIZING ============
        function applyParlaySizing() {
          var bankroll = parseFloat(document.getElementById("parlay-bankroll-input").value);
          var kelly = parseFloat(document.getElementById("parlay-kelly-input").value);
          var minEdge = parseFloat(document.getElementById("parlay-min-edge-input").value) || 0;
          if (!bankroll || bankroll <= 0 || !kelly || kelly <= 0 || kelly > 1) {
            showToast("Invalid sizing values", "error");
            return;
          }

          // Save all parlay settings then refresh to recalculate
          Promise.all([
            fetch("/api/sizing-settings", {
              method: "POST",
              headers: { "Content-Type": "application/json" },
              body: JSON.stringify({ param: "parlay_bankroll", value: bankroll })
            }),
            fetch("/api/sizing-settings", {
              method: "POST",
              headers: { "Content-Type": "application/json" },
              body: JSON.stringify({ param: "parlay_kelly_mult", value: kelly })
            }),
            fetch("/api/sizing-settings", {
              method: "POST",
              headers: { "Content-Type": "application/json" },
              body: JSON.stringify({ param: "parlay_min_edge", value: minEdge })
            })
          ]).then(function() {
            showToast("Parlay settings saved. Refreshing...", "success");
            refreshData();
          });
        }

        // Load parlay settings from server on page load (guard against missing elements)
        (function() {
          var pbi = document.getElementById("parlay-bankroll-input");
          var pki = document.getElementById("parlay-kelly-input");
          var pei = document.getElementById("parlay-min-edge-input");
          if (!pbi || !pki) return;
          fetch("/api/sizing-settings")
            .then(function(r) { return r.json(); })
            .then(function(s) {
              if (s.parlay_bankroll) pbi.value = s.parlay_bankroll;
              if (s.parlay_kelly_mult) pki.value = s.parlay_kelly_mult;
              if (pei && s.parlay_min_edge != null) {
                pei.value = s.parlay_min_edge;
                filterParlaysByEdge();
              }
            })
            .catch(function() {});
        })();

        // ============ PARLAY PLACEMENT ============
        function addPlacedParlayRow(btn, amount) {
          // Show the placed parlays section if hidden
          var header = document.getElementById("placed-parlays-header");
          var section = document.getElementById("placed-parlays-section");
          if (header) header.style.display = "";
          if (section) section.style.display = "";

          // Build leg description: "TeamName -1.5 · O7.0"
          var combo = btn.dataset.combo || "";
          var home = btn.dataset.home || "";
          var away = btn.dataset.away || "";
          var spreadTeam = combo.indexOf("Home") >= 0 ? home : away;
          var spread = btn.dataset.spread || "0";
          var spreadFmt = parseFloat(spread) > 0 ? "+" + spread : spread;
          var ouPrefix = combo.indexOf("Over") >= 0 ? "O" : "U";
          var legs = spreadTeam + " " + spreadFmt + " \\u00b7 " + ouPrefix + btn.dataset.total;

          var odds = btn.dataset.wzOdds || "";
          var oddsFmt = parseInt(odds) > 0 ? "+" + odds : odds;
          var edge = "+" + btn.dataset.edge + "%";
          var game = away + " @ " + home;
          var rec = "$" + Math.round(parseFloat(btn.dataset.size));
          var actual = "$" + Math.round(amount);

          var tr = document.createElement("tr");
          tr.innerHTML =
            \'<td>\' + escapeHtml(game) + \'</td>\' +
            \'<td>\' + escapeHtml(legs) + \'</td>\' +
            \'<td style="text-align:right; font-family:monospace;">\' + escapeHtml(oddsFmt) + \'</td>\' +
            \'<td style="text-align:right; color:#3fb950; font-weight:600;">\' + escapeHtml(edge) + \'</td>\' +
            \'<td style="text-align:right; font-weight:600;">\' + escapeHtml(actual) + \'</td>\' +
            \'<td style="text-align:right; color:#8b949e; font-size:0.8rem;">\' + escapeHtml(rec) + \'</td>\';

          var liveTable = document.getElementById("placed-parlays-live");
          if (liveTable) liveTable.appendChild(tr);
        }

        function placeParlay(btn) {
          var combo = btn.dataset.combo;
          var size = parseFloat(btn.dataset.size);
          var odds = btn.dataset.wzOdds;

          var overlay = document.createElement(\'div\');
          overlay.className = \'modal-overlay\';
          overlay.innerHTML =
            \'<div class="modal-box">\' +
              \'<div class="modal-title">Place Parlay</div>\' +
              \'<div class="modal-detail">Combo: <span>\' + escapeHtml(combo) + \'</span></div>\' +
              \'<div class="modal-detail">Spread: <span>\' + escapeHtml(btn.dataset.spread) + \'</span> | Total: <span>\' + escapeHtml(btn.dataset.total) + \'</span></div>\' +
              \'<div class="modal-detail">WZ Odds: <span>\' + escapeHtml(odds) + \'</span> | Edge: <span>+\' + escapeHtml(btn.dataset.edge) + \'%</span></div>\' +
              \'<div class="modal-recommended">Recommended: $\' + Math.round(size) + \'</div>\' +
              \'<div class="modal-input-group">\' +
                \'<label class="modal-input-label">Actual Amount ($)</label>\' +
                \'<input type="number" class="modal-input" value="\' + Math.round(size) + \'" step="1" min="1">\' +
              \'</div>\' +
              \'<div class="modal-actions">\' +
                \'<button class="modal-btn-cancel">Cancel</button>\' +
                \'<button class="modal-btn-confirm">Confirm</button>\' +
              \'</div>\' +
            \'</div>\';

          document.body.appendChild(overlay);

          overlay.querySelector(\'.modal-btn-cancel\').onclick = function() { overlay.remove(); };
          overlay.querySelector(\'.modal-btn-confirm\').onclick = function() {
            var amount = parseFloat(overlay.querySelector(\'.modal-input\').value);
            if (!amount || amount <= 0) { showToast("Invalid amount", "error"); return; }

            fetch("/api/place-parlay", {
              method: "POST",
              headers: { "Content-Type": "application/json" },
              body: JSON.stringify({
                parlay_hash: btn.dataset.hash,
                game_id: btn.dataset.gameId,
                home_team: btn.dataset.home,
                away_team: btn.dataset.away,
                game_time: btn.dataset.time || null,
                combo: btn.dataset.combo,
                spread_line: parseFloat(btn.dataset.spread) || null,
                total_line: parseFloat(btn.dataset.total) || null,
                fair_odds: parseInt(btn.dataset.fairOdds) || null,
                wz_odds: parseInt(btn.dataset.wzOdds),
                edge_pct: parseFloat(btn.dataset.edge) || null,
                kelly_bet: size,
                actual_size: amount
              })
            })
              .then(function(r) { return r.json(); })
              .then(function(result) {
                if (result.success) {
                  btn.className = "btn-placed";
                  btn.textContent = "Placed";
                  btn.onclick = function() { removeParlay(this); };
                  addPlacedParlayRow(btn, amount);
                  applyParlayFilters();
                  showToast("Parlay placed", "success");
                  applyParlayFilters();
                } else {
                  showToast(result.error || "Failed", "error");
                }
                overlay.remove();
              })
              .catch(function(err) { showToast("Error: " + err, "error"); overlay.remove(); });
          };
        }

        function removeParlay(btn) {
          fetch("/api/remove-parlay", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({ parlay_hash: btn.dataset.hash })
          })
            .then(function(r) { return r.json(); })
            .then(function(result) {
              if (result.success) {
                btn.className = "btn-place";
                btn.textContent = "Place";
                btn.onclick = function() { placeParlay(this); };
                showToast("Parlay removed", "success");
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

# Resolve NFLWork root — if running from a worktree, use main repo for mlb.duckdb
NFLWORK_ROOT <- dirname(dirname(DASHBOARD_DIR))
if (grepl(".claude/worktrees", NFLWORK_ROOT, fixed = TRUE)) {
  NFLWORK_ROOT <- sub("/.claude/worktrees.*", "", NFLWORK_ROOT)
}
setwd(NFLWORK_ROOT)

cat("=== MLB Answer Key Dashboard ===\n\n")

# Load bets from duckdb (saved by MLBCombine.R or equivalent)
cat("Loading bets from database...\n")
con <- dbConnect(duckdb(), dbdir = "Answer Keys/mlb.duckdb", read_only = TRUE)

# Check if table exists
tables <- dbListTables(con)
if (!"mlb_bets_combined" %in% tables) {
  all_bets <- tibble()
} else {
  all_bets <- dbGetQuery(con, "SELECT * FROM mlb_bets_combined") %>%
    filter(is.na(pt_start_time) | pt_start_time > Sys.time())
}
dbDisconnect(con, shutdown = TRUE)

cat(sprintf("Loaded %d bets\n", nrow(all_bets)))

placed_bets <- load_placed_bets(DB_PATH)
cat(sprintf("Found %d placed bets\n", nrow(placed_bets)))

# Calculate stats
stats <- list(
  total_bets = nrow(all_bets),
  placed_count = nrow(placed_bets),
  avg_ev = if (nrow(all_bets) > 0) mean(all_bets$ev, na.rm = TRUE) * 100 else 0,
  max_ev = if (nrow(all_bets) > 0) max(all_bets$ev, na.rm = TRUE) * 100 else 0
)

# Create tables
if (nrow(all_bets) > 0) {
  bets_table <- create_bets_table(all_bets, placed_bets)
} else {
  bets_table <- tags$div(
    style = "text-align: center; padding: 48px; color: #8b949e;",
    tags$p(style = "font-size: 1.1rem;", "No +EV bets found for upcoming games."),
    tags$p(style = "font-size: 0.85rem;", "Try clicking Refresh when more games are available.")
  )
}
placed_table <- create_placed_bets_table(placed_bets)

# Load parlay opportunities
cat("Loading parlay opportunities...\n")
parlay_opps <- tryCatch({
  pcon <- dbConnect(duckdb(), dbdir = "Answer Keys/mlb.duckdb", read_only = TRUE)
  result <- if ("mlb_parlay_opportunities" %in% dbListTables(pcon)) {
    dbGetQuery(pcon, "SELECT * FROM mlb_parlay_opportunities WHERE game_time IS NULL OR CAST(game_time AS TIMESTAMP) > NOW()")
  } else tibble()
  dbDisconnect(pcon, shutdown = TRUE)
  result
}, error = function(e) { cat(sprintf("  Warning: %s\n", e$message)); tibble() })

# Compatibility guard: spread_price / total_price added in feature/mlb-parlay-leg-display.
# If mlb_correlated_parlay.R hasn't been re-run yet, fill NA so dashboard doesn't crash.
if (nrow(parlay_opps) > 0) {
  if (!"spread_price" %in% names(parlay_opps)) parlay_opps$spread_price <- NA_integer_
  if (!"total_price"  %in% names(parlay_opps)) parlay_opps$total_price  <- NA_integer_
}

placed_parlays <- tryCatch({
  ppcon <- dbConnect(duckdb(), dbdir = DB_PATH, read_only = TRUE)
  result <- tryCatch(
    dbGetQuery(ppcon, "SELECT * FROM placed_parlays WHERE status = 'pending' AND (game_time IS NULL OR game_time > NOW())"),
    error = function(e) tibble(parlay_hash = character())
  )
  dbDisconnect(ppcon, shutdown = TRUE)
  result
}, error = function(e) tibble(parlay_hash = character()))

# Load min edge setting and filter parlay opportunities server-side
parlay_min_edge <- 0
tryCatch({
  mecon <- dbConnect(duckdb(), dbdir = DB_PATH, read_only = TRUE)
  me_val <- dbGetQuery(mecon, "SELECT value FROM sizing_settings WHERE param = 'parlay_min_edge'")
  dbDisconnect(mecon, shutdown = TRUE)
  if (nrow(me_val) > 0) parlay_min_edge <- me_val$value[1]
}, error = function(e) NULL)

if (nrow(parlay_opps) > 0 && parlay_min_edge > 0) {
  parlay_opps <- parlay_opps %>% filter(edge_pct >= parlay_min_edge)
}

cat(sprintf("Found %d parlay opportunities (min edge %.0f%%), %d placed parlays\n",
            nrow(parlay_opps), parlay_min_edge, nrow(placed_parlays)))

# Create parlay tables
if (nrow(parlay_opps) > 0) {
  parlays_table <- create_parlays_table(parlay_opps, placed_parlays)
} else {
  parlays_table <- NULL
}
placed_parlays_table <- create_placed_parlays_table(placed_parlays)

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
  correlations = I(c("Standalone", "Same Game")),
  statuses = I(c("Not Placed", "Placed", "Partial Fill"))
), auto_unbox = FALSE)

# Parlay filter options
parlay_games <- if (nrow(parlay_opps) > 0) {
  parlay_opps %>%
    mutate(game_label = paste(away_team, "@", home_team)) %>%
    distinct(game_label) %>%
    arrange(game_label) %>%
    pull(game_label)
} else character()

parlay_filter_options_json <- toJSON(list(
  games    = I(parlay_games),
  statuses = I(c("Not Placed", "Placed"))
), auto_unbox = FALSE)

# Create report
timestamp <- format(Sys.time(), "%b %d, %Y %I:%M %p")
page <- create_report(bets_table, placed_table, stats, timestamp, filter_options_json,
                       parlays_table, placed_parlays_table, parlay_opps,
                       parlay_filter_options_json)

# Save
save_html(page, OUTPUT_PATH, libdir = "lib")
cat("\nDashboard saved to:", OUTPUT_PATH, "\n")
cat("Open http://localhost:8083 to view\n")
