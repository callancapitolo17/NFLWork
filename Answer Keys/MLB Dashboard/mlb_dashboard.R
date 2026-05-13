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

# Pure HTML helpers shared with other dashboards
source(file.path(DASHBOARD_DIR, "..", "books_strip.R"))
# Conditional Kelly residual solver (Task 1) — drives Parlay tab combo residuals
source(file.path(DASHBOARD_DIR, "..", "conditional_kelly.R"))
# Book pill renderer (Task 1) — renders individual book prices
source(file.path(DASHBOARD_DIR, "book_pill.R"))
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

load_trifecta_opps <- function(mlb_db) {
  if (!file.exists(mlb_db)) return(tibble())
  con <- dbConnect(duckdb(), dbdir = mlb_db, read_only = TRUE)
  on.exit(dbDisconnect(con, shutdown = TRUE))
  tryCatch({
    if (!"mlb_trifecta_opportunities" %in% dbListTables(con)) return(tibble())
    dbGetQuery(con, "
      SELECT * FROM mlb_trifecta_opportunities
      WHERE game_time IS NULL OR CAST(game_time AS TIMESTAMP) > NOW()
    ")
  }, error = function(e) tibble())
}

load_placed_trifectas <- function(db_path) {
  if (!file.exists(db_path)) return(tibble())
  con <- dbConnect(duckdb(), dbdir = db_path, read_only = TRUE)
  on.exit(dbDisconnect(con, shutdown = TRUE))
  tryCatch({
    if (!"placed_trifectas" %in% dbListTables(con)) return(tibble())
    dbGetQuery(con, "
      SELECT * FROM placed_trifectas
      WHERE game_time IS NULL OR CAST(game_time AS TIMESTAMP) > NOW()
    ")
  }, error = function(e) tibble())
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

# Compact pill labels for placement statuses. Keep parallel with
# _SHORT_LABEL_BY_STATUS / _SHORT_LABEL_BY_REJECT_KEY in
# Answer Keys/MLB Dashboard/mlb_dashboard_server.py — both must agree
# so a fresh placement (JS render) and an already-stored row (R render)
# show the same pill text.
short_label_for_status <- function(status, error_key = "") {
  if (identical(status, "rejected")) {
    map <- c(
      insufficient_funds = "insufficient $",
      bet_too_large      = "exceeds limit",
      line_unavailable   = "line pulled"
    )
    out <- map[as.character(error_key)]
    return(if (is.na(out)) "rejected" else unname(out))
  }
  map <- c(
    placed         = "placed",
    price_moved    = "price moved",
    auth_error     = "auth fail",
    network_error  = "network err",
    orphaned       = "orphaned"
  )
  out <- map[as.character(status)]
  if (is.na(out)) status else unname(out)
}

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

# Returns parlay_opps with kelly_bet overridden to the conditional residual
# when the row appears in any combo in placed_parlays. Adds a combo_residual_note
# column that the Size column renderer uses for the annotation.
#
# Defensive: if placed_parlays is missing the is_combo column (schema migration
# not yet applied), this is a no-op.
apply_combo_residuals <- function(parlay_opps, placed_parlays, parlay_bankroll, parlay_kelly_mult) {
  parlay_opps$combo_residual_note <- NA_character_

  if (is.null(placed_parlays) || nrow(placed_parlays) == 0) return(parlay_opps)
  if (!"is_combo" %in% names(placed_parlays))                return(parlay_opps)

  combos <- placed_parlays %>% filter(is_combo == TRUE)
  if (nrow(combos) == 0) return(parlay_opps)

  for (i in seq_len(nrow(combos))) {
    leg_ids <- tryCatch(
      jsonlite::fromJSON(combos$combo_leg_ids[i]),
      error = function(e) NULL
    )
    if (is.null(leg_ids) || length(leg_ids) != 2) next

    # combo's wz_odds is American (INTEGER per schema); convert to decimal.
    wz_american <- combos$wz_odds[i]
    if (!is.finite(wz_american)) next
    combo_dec <- if (wz_american >= 0) {
      wz_american / 100 + 1
    } else {
      100 / abs(wz_american) + 1
    }

    s_combo <- combos$actual_size[i]
    if (!is.finite(s_combo) || s_combo <= 0) next

    row_a <- parlay_opps %>% filter(parlay_hash == leg_ids[1])
    row_b <- parlay_opps %>% filter(parlay_hash == leg_ids[2])
    if (nrow(row_a) == 0 || nrow(row_b) == 0) next

    res <- conditional_kelly_residuals(
      p_a     = 1 / row_a$fair_dec,
      d_a     = row_a$wz_dec,
      p_b     = 1 / row_b$fair_dec,
      d_b     = row_b$wz_dec,
      s_combo = s_combo,
      d_combo = combo_dec,
      bankroll   = parlay_bankroll,
      kelly_mult = parlay_kelly_mult
    )

    parlay_opps$kelly_bet[parlay_opps$parlay_hash == leg_ids[1]] <- res$s_a
    parlay_opps$kelly_bet[parlay_opps$parlay_hash == leg_ids[2]] <- res$s_b
    parlay_opps$combo_residual_note[parlay_opps$parlay_hash %in% leg_ids] <-
      sprintf("(residual after combo $%.0f)", s_combo)
  }
  parlay_opps
}

create_parlays_table <- function(parlay_opps, placed_parlays, parlay_bankroll = 100, parlay_kelly_mult = 0.25) {
  parlay_opps <- apply_combo_residuals(parlay_opps, placed_parlays, parlay_bankroll, parlay_kelly_mult)

  # Build a lookup from parlay_hash → placement state columns for the
  # auto-placement Place column. "pending" = manually placed (legacy path);
  # "placed" = auto-placed (show ticket); error statuses = show red pill.
  # Rows with no entry in placed_parlays stay NA.
  placement_cols <- if (nrow(placed_parlays) > 0 && "parlay_hash" %in% names(placed_parlays)) {
    placed_parlays %>%
      select(parlay_hash,
             placement_status     = status,
             ticket_number        = any_of("ticket_number"),
             placement_error      = any_of("error_msg"),
             placement_error_key  = any_of("error_msg_key")) %>%
      # Ensure all four columns exist even on legacy rows that pre-date Task 9 migration
      { if (!"ticket_number"        %in% names(.)) mutate(., ticket_number        = NA_character_) else . } %>%
      { if (!"placement_error"      %in% names(.)) mutate(., placement_error      = NA_character_) else . } %>%
      { if (!"placement_error_key"  %in% names(.)) mutate(., placement_error_key  = NA_character_) else . }
  } else {
    tibble(parlay_hash = character(), placement_status = character(),
           ticket_number = character(), placement_error = character(),
           placement_error_key = character())
  }

  # Filter out the combo rows themselves from the "placed" set — they're not
  # source parlays so they shouldn't be marked as already-placed in the parlay table.
  placed_hashes <- if (nrow(placed_parlays) > 0) {
    if ("is_combo" %in% names(placed_parlays)) {
      is_combo_vec <- placed_parlays$is_combo
      is_combo_vec[is.na(is_combo_vec)] <- FALSE  # treat NA as not-combo
      placed_parlays$parlay_hash[!is_combo_vec]
    } else {
      placed_parlays$parlay_hash
    }
  } else character()

  # Ensure exact-payout columns exist even on older rows (pre-Stage-2 backfill).
  # parlay_pricer.py --exact-payouts populates these; fall back to Kelly+wz_dec
  # arithmetic when they're NA.
  if (!"exact_wager" %in% names(parlay_opps)) parlay_opps$exact_wager <- NA_integer_
  if (!"exact_to_win" %in% names(parlay_opps)) parlay_opps$exact_to_win <- NA_integer_

  table_data <- parlay_opps %>%
    left_join(placement_cols, by = "parlay_hash") %>%
    mutate(
      is_placed = parlay_hash %in% placed_hashes,
      sel = "",
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
      # When this row participates in a placed combo, show the conditional residual
      # (kelly_bet was already overridden by apply_combo_residuals) plus an annotation.
      size_display   = ifelse(
        !is.na(combo_residual_note),
        sprintf("$%.0f<br><span class='combo-note'>%s</span>", kelly_bet, combo_residual_note),
        sprintf("$%.0f", coalesce(as.numeric(exact_wager), kelly_bet))
      ),
      to_win_display = sprintf("$%.0f", coalesce(as.numeric(exact_to_win),
                                                  round(kelly_bet * (wz_dec - 1)))),
      corr_display   = sprintf("%.3f", corr_factor),
      books_strip    = NA  # placeholder; rendered by the colDef's cell function
    ) %>%
    arrange(desc(edge_pct))

  reactable(
    table_data,
    compact = TRUE,
    striped = TRUE,
    highlight = TRUE,
    defaultPageSize = 25,
    columns = list(
      sel = colDef(
        name = "",
        minWidth = 30,
        align = "center",
        filterable = FALSE,
        sortable = FALSE,
        class = "cell-sel",
        html = TRUE,
        cell = function(value, index) {
          row <- table_data[index, ]
          if (isTRUE(row$is_placed) || isTRUE(row$is_combo)) {
            ''  # no checkbox if already placed or this is the combo row
          } else {
            sprintf(
              '<input type="checkbox" class="combo-select" data-hash="%s" data-game-id="%s" onchange="onComboSelectChange(this)">',
              row$parlay_hash, row$game_id
            )
          }
        }
      ),
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
      n_books_blended = colDef(show = FALSE),
      spread_team = colDef(show = FALSE),
      spread_fmt = colDef(show = FALSE),
      sp_price_fmt = colDef(show = FALSE),
      ou_prefix = colDef(show = FALSE),
      tot_price_fmt = colDef(show = FALSE),
      combo_residual_note = colDef(show = FALSE),

      # Visible columns
      game = colDef(
        name = "Game",
        minWidth = 180,
        class = "cell-game",
        html = TRUE,
        cell = JS("function(cellInfo) {
          var matchup = cellInfo.value || '';
          var t = cellInfo.row.game_time;
          if (!t) return matchup;
          var d = new Date(t);
          if (isNaN(d)) return matchup;
          var opts = {weekday:'short', month:'2-digit', day:'2-digit', hour:'numeric', minute:'2-digit'};
          var time = d.toLocaleString(undefined, opts);
          return matchup + '<div style=\"color:#8b949e;font-size:10px\">' + time + '</div>';
        }")
      ),
      game_time = colDef(show = FALSE),
      legs_display = colDef(name = "Legs", minWidth = 260, class = "cell-legs"),
      fair_display = colDef(name = "Fair", minWidth = 70, align = "right",
        class = "cell-fair",
        style = list(fontFamily = "monospace", color = "#8b949e")),
      # Per-book numeric columns folded into the Books pill row below.
      # Hidden (data still in the dataframe) so the cell renderer can read them.
      dk_fair_prob = colDef(show = FALSE),
      fd_fair_prob = colDef(show = FALSE),
      px_fair_prob = colDef(show = FALSE),
      nv_fair_prob = colDef(show = FALSE),
      # Combined Books cell: M / DK / FD / PX / NV / Cons pill row.
      # Reads model_prob_raw (always non-NA in real data — pricer skips no-sample
      # games upstream) + the four per-book devigged probs + blended_prob_raw.
      books_strip = colDef(
        name = "Books (devigged fair %)",
        minWidth = 320,
        html = TRUE,
        sortable = FALSE,
        class = "cell-books",
        cell = function(value, index) {
          row <- table_data[index, ]
          render_books_strip(
            model = row$model_prob_raw,
            dk    = row$dk_fair_prob,
            fd    = row$fd_fair_prob,
            px    = row$px_fair_prob,
            nv    = row$nv_fair_prob,
            cons  = row$blended_prob_raw
          )
        }
      ),
      wz_display = colDef(name = "WZ", minWidth = 70, align = "right",
        class = "cell-wz",
        style = list(fontFamily = "monospace")),
      corr_display = colDef(show = FALSE),
      edge_display = colDef(name = "Edge %", minWidth = 70, align = "right",
        class = "cell-edge",
        cell = function(value, index) {
          ep <- table_data$edge_pct[index]
          color <- if (ep >= 15) "#3fb950" else if (ep >= 10) "#56d364" else if (ep >= 5) "#7ee787" else "#a5d6a7"
          span(style = list(color = color, fontWeight = "600"), value)
        }
      ),
      size_display = colDef(name = "Size", minWidth = 65, align = "right", html = TRUE, class = "cell-size"),
      to_win_display = colDef(name = "To Win", minWidth = 65, align = "right",
        class = "cell-towin",
        style = list(color = "#3fb950")),
      # Auto-placement state columns — hidden data carriers for the cell renderer below
      placement_status     = colDef(show = FALSE),
      ticket_number        = colDef(show = FALSE),
      placement_error      = colDef(show = FALSE),
      placement_error_key  = colDef(show = FALSE),
      is_placed = colDef(
        name = "Action",
        minWidth = 110,
        align = "center",
        filterable = FALSE,
        class = "cell-action",
        html = TRUE,
        cell = function(value, index) {
          row <- table_data[index, ]
          ps  <- row$placement_status
          # Guard: treat NA or "NA" string uniformly
          ps_valid <- !is.na(ps) && nchar(ps) > 0 && ps != "NA"

          # Common data-* attributes mirror what's on the Place button below,
          # so the JS row filter (filterParlaysByEdge) can read edge / size /
          # game from any element with data-edge — button OR placed-span OR
          # error-pill — and classify the row's status uniformly. Without
          # these, placed/failed rows fall through filters as edge=0 + status
          # ="Not Placed" and get triple-excluded.
          common_data_attrs <- sprintf(
            'data-hash="%s" data-home="%s" data-away="%s" data-edge="%s" data-size="%s"',
            htmltools::htmlEscape(row$parlay_hash),
            htmltools::htmlEscape(row$home_team),
            htmltools::htmlEscape(row$away_team),
            row$edge_pct,
            if (!is.na(row$exact_wager)) row$exact_wager else row$kelly_bet
          )

          # Auto-placed: show muted "placed &middot; #<ticket>" label.
          # If ticket is missing (manual log via /api/log-parlay, or auto-place
          # response that didn't carry the identifier), drop the "&middot; #"
          # suffix and just show "placed" — looks intentional rather than
          # broken.
          # Use HTML entity instead of literal U+00B7 to keep the file
          # ASCII-clean — R writes report.html in the platform default
          # encoding (Latin-1 on macOS), which Flask's UTF-8 reader chokes on.
          if (ps_valid && ps == "placed") {
            has_ticket <- !is.na(row$ticket_number) && nchar(as.character(row$ticket_number)) > 0
            label_html <- if (has_ticket) {
              sprintf('placed &middot; #%s', htmltools::htmlEscape(row$ticket_number))
            } else {
              'placed'
            }
            return(sprintf('<span class="placed-parlay-label" %s>%s</span>',
                           common_data_attrs, label_html))
          }

          # Error states: show red pill with a SHORT label for quick scan,
          # full message in the title= tooltip on hover.
          error_statuses <- c("price_moved", "rejected",
                              "auth_error", "network_error", "orphaned")
          if (ps_valid && ps %in% error_statuses) {
            full_msg <- if (!is.na(row$placement_error) && nchar(row$placement_error) > 0)
              row$placement_error else ps
            err_key <- if ("placement_error_key" %in% names(row) &&
                           !is.na(row$placement_error_key)) row$placement_error_key else ""
            short_label <- short_label_for_status(ps, err_key)
            return(sprintf('<span class="pill error" title="%s" %s>%s</span>',
                           htmltools::htmlEscape(full_msg),
                           common_data_attrs,
                           htmltools::htmlEscape(short_label)))
          }

          # Default: Place button (manually placed via "pending" status OR no record)
          # data-size drives the bet amount. Prefer exact_wager from Stage 2;
          # fall back to Kelly-ideal if Stage 2 couldn't price this row.
          # data-risk mirrors data-size and is the attribute the Phase 6
          # insufficient-balance warning sweep reads (decoupled name so a
          # future change to data-size semantics doesn't break the warning).
          data_size <- if (!is.na(row$exact_wager)) row$exact_wager else row$kelly_bet
          data_attrs <- sprintf(
            'data-hash="%s" data-game-id="%s" data-home="%s" data-away="%s" data-time="%s" data-combo="%s" data-spread="%s" data-total="%s" data-fair-odds="%s" data-wz-odds="%s" data-edge="%s" data-size="%s" data-risk="%s"',
            row$parlay_hash, row$game_id, row$home_team, row$away_team,
            ifelse(is.na(row$game_time), "", as.character(row$game_time)),
            row$combo, row$spread_line, row$total_line,
            row$fair_odds, row$wz_odds, row$edge_pct, data_size, data_size
          )
          if (value) {
            sprintf('<button class="btn-placed" onclick="removeParlay(this)" %s>Placed</button>', data_attrs)
          } else {
            # Two buttons side-by-side, plus a transient warning span. The
            # warning is populated by window._wzRecomputeWarnings() whenever
            # the selector changes, balances refresh, or the parlay table
            # re-renders. Empty by default; gets text only when the selected
            # account's available balance is below the row's risk.
            sprintf(paste0(
              '<button class="btn-place" onclick="placeParlay(this)" %s>Place</button>',
              ' ',
              '<button class="btn-log" onclick="logParlay(this)" %s>Log</button>',
              ' ',
              '<span class="wz-insufficient-warning" style="margin-left:6px; color:#f85149; font-size:11px;"></span>'
            ), data_attrs, data_attrs)
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

create_trifectas_table <- function(trifecta_opps, placed_trifectas) {
  if (nrow(trifecta_opps) == 0) {
    return(tags$div(
      style = "text-align: center; padding: 48px; color: #8b949e;",
      tags$p(style = "font-size: 1.1rem;", "No trifectas priced yet."),
      tags$p(style = "font-size: 0.85rem;",
             "Click Refresh after Wagerzon posts today's TRIPLE-PLAY / GRAND-SLAM lines.")
    ))
  }

  # Mark placed rows
  placed_set <- if (nrow(placed_trifectas) > 0) placed_trifectas$trifecta_hash else character(0)

  table_data <- trifecta_opps %>%
    mutate(
      is_placed     = trifecta_hash %in% placed_set,
      pick_display  = sprintf("%s (%s)", target_team, side),
      model_display = ifelse(is.na(model_odds), "—",
                             ifelse(model_odds > 0, paste0("+", model_odds),
                                    as.character(model_odds))),
      dk_display    = ifelse(is.na(dk_odds), "—",
                             ifelse(dk_odds > 0, paste0("+", dk_odds),
                                    as.character(dk_odds))),
      fair_display  = ifelse(is.na(fair_odds), "—",
                             ifelse(fair_odds > 0, paste0("+", fair_odds),
                                    as.character(fair_odds))),
      book_display  = ifelse(book_odds > 0, paste0("+", book_odds),
                             as.character(book_odds)),
      edge_display  = sprintf("%+.1f%%", edge_pct),
      stake_display = ifelse(kelly_bet > 0, sprintf("$%.0f", kelly_bet), "—")
    ) %>%
    arrange(desc(edge_pct))

  # The Action column needs to be a real column on the dataframe so reactable
  # has something to render. Add a placeholder; the cell renderer below builds
  # the actual button HTML from the row's other fields.
  table_data$action <- ""

  reactable(
    table_data,
    searchable = TRUE, filterable = TRUE,
    striped = TRUE, highlight = TRUE, compact = TRUE,
    defaultPageSize = 25,
    columns = list(
      # Hidden helpers (carried for JS data-* attrs and conditional logic)
      trifecta_hash = colDef(show = FALSE),
      game_id       = colDef(show = FALSE),
      game_time     = colDef(show = FALSE),
      target_team   = colDef(show = FALSE),
      side          = colDef(show = FALSE),
      n_samples     = colDef(show = FALSE),
      model_odds    = colDef(show = FALSE),
      dk_odds       = colDef(show = FALSE),
      fair_odds     = colDef(show = FALSE),
      book_odds     = colDef(show = FALSE),
      edge_pct      = colDef(show = FALSE),
      kelly_bet     = colDef(show = FALSE),
      is_placed     = colDef(show = FALSE),

      # Visible columns, in display order
      game = colDef(name = "Game", minWidth = 180, filterable = TRUE),
      pick_display = colDef(
        name = "Pick", minWidth = 140,
        cell = function(value, index) {
          row <- table_data[index, ]
          div(
            span(style = "font-weight: 600;", row$target_team),
            span(style = "margin-left: 8px; color: #888; font-size: 0.9em;",
                 paste0("(", row$side, ")"))
          )
        }
      ),
      prop_type    = colDef(name = "Prop",        minWidth = 110, filterable = TRUE),
      description  = colDef(name = "Description", minWidth = 220, filterable = TRUE),
      model_display = colDef(name = "Model", minWidth = 70, align = "right",
                             style = list(fontFamily = "monospace")),
      dk_display    = colDef(name = "DK",    minWidth = 70, align = "right",
                             style = list(fontFamily = "monospace")),
      fair_display  = colDef(name = "Fair",  minWidth = 70, align = "right",
                             style = list(fontFamily = "monospace", fontWeight = "600")),
      book_display  = colDef(name = "Book",  minWidth = 70, align = "right",
                             style = list(fontFamily = "monospace")),
      edge_display = colDef(
        name = "Edge", minWidth = 80, align = "right",
        cell = function(value, index) {
          ev <- table_data$edge_pct[index]
          color <- if (is.na(ev)) "#8b949e"
                   else if (ev >= 15) "#3fb950"
                   else if (ev >= 10) "#56d364"
                   else if (ev >=  5) "#7ee787"
                   else                "#a5d6a7"
          div(style = list(color = color, fontWeight = "600"), value)
        }
      ),
      stake_display = colDef(name = "Stake", minWidth = 70, align = "right"),
      action = colDef(
        name = "Action", minWidth = 110, align = "center", html = TRUE,
        sortable = FALSE, filterable = FALSE,
        cell = function(value, index) {
          row <- table_data[index, ]
          # data-edge / data-size on every row so filterTrifectasByEdge() can find them
          edge_attrs <- sprintf(
            'data-edge="%.2f" data-size="%.2f"',
            ifelse(is.na(row$edge_pct), 0, row$edge_pct),
            ifelse(is.na(row$kelly_bet), 0, row$kelly_bet)
          )
          if (isTRUE(row$is_placed)) {
            sprintf(
              '<button class="btn-placed" %s data-trifecta-hash="%s" onclick="removeTrifecta(this)">Placed</button>',
              edge_attrs, row$trifecta_hash
            )
          } else if (!is.na(row$kelly_bet) && row$kelly_bet > 0) {
            sprintf(
              '<button class="btn-place" %s data-trifecta-hash="%s" data-actual-wager="%.2f" onclick="placeTrifecta(this)">Place</button>',
              edge_attrs, row$trifecta_hash, row$kelly_bet
            )
          } else {
            # Below threshold: empty marker span so live filter still has data-edge to read
            sprintf('<span class="trifecta-marker" %s></span>', edge_attrs)
          }
        }
      )
    ),
    theme = reactableTheme(
      backgroundColor = "#0d1117", color = "#c9d1d9",
      borderColor = "#30363d", stripedColor = "#161b22"
    )
  )
}

create_bets_table_legacy <- function(all_bets, placed_bets) {
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

# -----------------------------------------------------------------------------
# create_bets_table — card layout (Task 11 of odds-screen rebuild)
# -----------------------------------------------------------------------------
# Replaces the flat per-bet table with a card-per-bet layout that shows the
# pick + opposite side's price across all 8 sportsbooks as pill rows. Each
# pill is rendered server-side via render_book_pill() with three states:
#   - exact line match    -> standard pill
#   - mismatched line     -> amber-tinted pill + line tag
#   - no quote            -> muted dashed em-dash pill
# The pick book's pill gets the green `.pick` override.
#
# Falls back to the legacy flat layout when book_prices_wide is NULL (e.g.,
# when the loader couldn't read mlb_bets_book_prices).
create_bets_table <- function(all_bets, placed_bets, book_prices_wide = NULL) {
  if (is.null(book_prices_wide)) {
    warning("[bets-tab] book_prices_wide is NULL — using legacy flat layout")
    return(create_bets_table_legacy(all_bets, placed_bets))
  }

  placed_hashes <- if (nrow(placed_bets) > 0) placed_bets$bet_hash else character()

  # Lookups for placement state (status + ticket number) used by the action cell.
  placed_status_lookup <- if (nrow(placed_bets) > 0 && "status" %in% names(placed_bets)) {
    setNames(placed_bets$status, placed_bets$bet_hash)
  } else setNames(character(), character())
  placed_ticket_lookup <- if (nrow(placed_bets) > 0 && "ticket_number" %in% names(placed_bets)) {
    setNames(placed_bets$ticket_number, placed_bets$bet_hash)
  } else setNames(character(), character())
  # Keep the legacy placed_actual / placed_rec lookups around so we can still
  # surface "partial fill" data via the action-cell data attrs (Task 12 hook).
  placed_actual_lookup <- setNames(
    if (nrow(placed_bets) > 0 && "actual_size" %in% names(placed_bets)) placed_bets$actual_size else numeric(),
    if (nrow(placed_bets) > 0) placed_bets$bet_hash else character()
  )

  # Same-game correlation lookup (reuse the existing helper unchanged).
  same_game_info <- lapply(seq_len(nrow(all_bets)), function(i) {
    find_same_game_bets(i, all_bets, placed_bets)
  })

  # Re-compute bet_row_id on the bets frame using the same formula MLB.R uses
  # so the join to book_prices_wide matches mlb_bets_book_prices.
  all_bets <- all_bets %>%
    mutate(bet_row_id = vapply(
      paste(id, market, ifelse(is.na(line), "", as.character(line)), bet_on, sep = "|"),
      function(s) digest::digest(s, algo = "md5"), character(1)
    ))

  BOOK_ORDER  <- c("wagerzon", "hoop88", "bfa", "bookmaker", "bet105",
                   "draftkings", "fanduel", "pinnacle")
  BOOK_LABELS <- c(wagerzon = "WZ", hoop88 = "H88", bfa = "BFA",
                   bookmaker = "BKM", bet105 = "B105",
                   draftkings = "DK", fanduel = "FD", pinnacle = "Pinn")

  # Render one full side-row (label + 8 book pills).
  render_side_row <- function(wide_row, side_label_text, is_pick_side,
                              pick_book, side_word) {
    pills <- vapply(BOOK_ORDER, function(b) {
      odds_col  <- paste0(b, "_american_odds")
      lq_col    <- paste0(b, "_line_quoted")
      exact_col <- paste0(b, "_is_exact_line")
      odds  <- if (odds_col  %in% names(wide_row)) wide_row[[odds_col]]  else NA_integer_
      lq    <- if (lq_col    %in% names(wide_row)) wide_row[[lq_col]]    else NA_real_
      exact <- if (exact_col %in% names(wide_row)) wide_row[[exact_col]] else NA
      render_book_pill(
        book          = BOOK_LABELS[[b]],
        american_odds = if (is.na(odds)) NA_integer_ else as.integer(odds),
        line_quoted   = lq,
        is_exact_line = exact,
        is_pick       = is_pick_side && (b == pick_book),
        side          = side_word
      )
    }, character(1))
    paste0(
      '<div class="side-row">',
      sprintf('<span class="side-label">%s</span>',
              htmltools::htmlEscape(side_label_text)),
      paste(pills, collapse = ""),
      '</div>'
    )
  }

  # Display data: one card per bet (no fan-out by book; pills carry that info).
  table_data <- all_bets %>%
    mutate(
      bet_hash = pmap_chr(list(id, market, bet_on, line), generate_bet_hash),
      is_placed = bet_hash %in% placed_hashes,
      game_display = paste(away_team, "@", home_team),
      ev_pct = ev * 100,
      ev_display = ifelse(ev >= 0, sprintf("+%.1f%%", ev * 100),
                                   sprintf("%.1f%%", ev * 100)),
      m_display = sprintf("%.1f%%", prob * 100),
      size_display = sprintf("$%.0f", bet_size),
      towin_display = sprintf("$%.0f",
                              ifelse(odds > 0, bet_size * odds / 100,
                                               bet_size * 100 / abs(odds))),
      pick_display = sprintf("%s %s",
                             toupper(bookmaker_key),
                             ifelse(odds > 0, paste0("+", odds),
                                              as.character(odds))),
      market_display = paste(
        format_market_name(market),
        ifelse(is.na(line), bet_on,
               paste(bet_on,
                     ifelse(line > 0, paste0("+", line), as.character(line))))
      )
    ) %>%
    arrange(desc(ev))

  # Partial-fill detection (mirrors legacy behavior)
  placed_actual <- if (nrow(placed_bets) > 0 && "actual_size" %in% names(placed_bets)) {
    setNames(placed_bets$actual_size, placed_bets$bet_hash)
  } else setNames(numeric(), character())

  table_data <- table_data %>%
    mutate(
      placed_actual = ifelse(is_placed, placed_actual[bet_hash], NA_real_),
      fill_status = case_when(
        !is_placed                                           ~ "not_placed",
        is.na(placed_actual)                                  ~ "placed",
        round(placed_actual) >= round(bet_size)               ~ "placed",
        TRUE                                                  ~ "partial"
      ),
      fill_diff = ifelse(fill_status == "partial",
                         round(bet_size) - round(placed_actual),
                         NA_real_)
    )

  # Pre-render the pick-side HTML for each row. Picks side_word from the
  # bet_on text ("Over X" -> "over"; "Under X" -> "under"; everything else
  # defaults to "over" — only used for the line tag prefix on mismatches).
  table_data$pickside_html <- vapply(seq_len(nrow(table_data)), function(i) {
    bet_id <- table_data$bet_row_id[i]
    wide_pick <- book_prices_wide %>% filter(bet_row_id == bet_id, side == "pick")
    if (nrow(wide_pick) == 0) return("")
    side_word <- if (grepl("^Over",  table_data$bet_on[i], ignore.case = TRUE)) "over"
                 else if (grepl("^Under", table_data$bet_on[i], ignore.case = TRUE)) "under"
                 else "over"
    render_side_row(
      wide_row         = wide_pick[1, ],
      side_label_text  = table_data$bet_on[i],
      is_pick_side     = TRUE,
      pick_book        = table_data$bookmaker_key[i],
      side_word        = side_word
    )
  }, character(1))

  # Pre-render the opposite-side HTML.
  table_data$otherside_html <- vapply(seq_len(nrow(table_data)), function(i) {
    bet_id <- table_data$bet_row_id[i]
    wide_opp <- book_prices_wide %>% filter(bet_row_id == bet_id, side == "opposite")
    if (nrow(wide_opp) == 0) return("")
    opposite_label <- if (grepl("^Over",  table_data$bet_on[i], ignore.case = TRUE)) "Under"
                      else if (grepl("^Under", table_data$bet_on[i], ignore.case = TRUE)) "Over"
                      else "Opp"
    side_word <- if (grepl("^Over", opposite_label, ignore.case = TRUE)) "over" else "under"
    render_side_row(
      wide_row         = wide_opp[1, ],
      side_label_text  = opposite_label,
      is_pick_side     = FALSE,
      pick_book        = table_data$bookmaker_key[i],
      side_word        = side_word
    )
  }, character(1))

  # Same-game correlation corner badge. Reuse the legacy tooltip builder so
  # the hover detail (markets/sizes/books for related legs) is unchanged.
  table_data$corr_html <- vapply(seq_len(nrow(table_data)), function(i) {
    info <- same_game_info[[i]]
    if (!info$has_same_game) return("")
    dot <- intToUtf8(0xB7)
    chk <- intToUtf8(0x2713)
    cir <- intToUtf8(0x2013)
    lines <- vapply(info$details, function(d) {
      market_name <- format_market_name(d$market)
      line_str <- if (!is.null(d$line) && !is.na(d$line)) {
        if (d$line > 0) paste0(" +", d$line) else paste0(" ", d$line)
      } else ""
      odds_str <- if (!is.null(d$odds) && !is.na(d$odds)) {
        if (d$odds > 0) sprintf(" (%+d)", d$odds) else sprintf(" (%d)", d$odds)
      } else ""
      if (isTRUE(d$is_placed)) {
        size_str <- if (!is.null(d$actual_size) && !is.na(d$actual_size)) {
          sprintf("$%.0f", d$actual_size)
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
    }, character(1))
    tooltip <- paste(lines, collapse = "\n")
    sprintf('<span class="corr-badge" title="%s">&#9679;</span>',
            escape_tooltip(tooltip))
  }, character(1))

  # Drop columns that are nominally numeric but contain NAs htmlwidgets
  # serializes as "NA" strings. Mixed types in a single JSON column crash
  # reactable's React widget silently. `cents` is Kalshi-only; `placed_actual`
  # and `fill_diff` are mostly-NA legacy carriers from create_bets_table_legacy.
  # The new card layout doesn't need any of them. `any_of` tolerates absent
  # columns (cents disappears on slates with no Kalshi bet).
  table_data <- table_data %>% select(-any_of(c("cents", "placed_actual", "fill_diff")))

  # The reactable. elementId = "bets-table" so the wrapping container gets
  # id="bets-table-container", matching the CSS scope.
  reactable(
    table_data,
    elementId = "bets-table",
    searchable = TRUE,
    filterable = TRUE,
    defaultPageSize = 25,
    pageSizeOptions = c(25, 50, 100),
    showPageSizeOptions = TRUE,
    columns = list(
      # Hidden data carriers
      bet_hash       = colDef(show = FALSE),
      bet_row_id     = colDef(show = FALSE),
      id             = colDef(show = FALSE),
      home_team      = colDef(show = FALSE),
      away_team      = colDef(show = FALSE),
      market         = colDef(show = FALSE),
      market_type    = colDef(show = FALSE),
      line           = colDef(show = FALSE),
      bet_on         = colDef(show = FALSE),
      bet_size       = colDef(show = FALSE),
      odds           = colDef(show = FALSE),
      prob           = colDef(show = FALSE),
      ev             = colDef(show = FALSE),
      ev_pct         = colDef(show = FALSE),
      bookmaker_key  = colDef(show = FALSE),
      pt_start_time  = colDef(show = FALSE),
      corr_html      = colDef(show = FALSE),
      correlation_adj = colDef(show = FALSE),
      fill_status     = colDef(show = FALSE),

      # Visible cells (ordering via CSS `order:` in #bets-table-container).
      game_display = colDef(
        name = "", html = TRUE, class = "cell-game",
        cell = function(value, index) {
          row <- table_data[index, ]
          time_str <- if (!is.na(row$pt_start_time))
            format(row$pt_start_time, "%a %I:%M %p") else ""
          paste0(
            htmltools::htmlEscape(value),
            sprintf(' <span style="color:#8b949e;font-size:11px;margin-left:8px">%s</span>',
                    htmltools::htmlEscape(time_str)),
            row$corr_html
          )
        }
      ),
      market_display = colDef(
        name = "", class = "cell-market", html = TRUE,
        cell = function(value, index) htmltools::htmlEscape(value)
      ),
      pickside_html = colDef(
        name = "", html = TRUE, class = "cell-pickside",
        cell = function(value, index) value
      ),
      otherside_html = colDef(
        name = "", html = TRUE, class = "cell-otherside",
        cell = function(value, index) value
      ),
      m_display = colDef(
        name = "", class = "cell-m",
        style = list(color = "#7ee787", fontWeight = 600)
      ),
      pick_display = colDef(
        name = "", class = "cell-pick"
      ),
      ev_display = colDef(
        name = "", class = "cell-ev",
        cell = function(value, index) {
          ep <- table_data$ev[index] * 100
          color <- if (ep >= 15) "#3fb950"
                   else if (ep >= 10) "#56d364"
                   else if (ep >= 5) "#7ee787"
                   else "#a5d6a7"
          div(style = list(color = color, fontWeight = "600"), value)
        }
      ),
      size_display = colDef(
        name = "", class = "cell-size",
        html = TRUE,
        cell = function(value, index) {
          sz <- table_data$bet_size[index]
          sprintf('<span data-bet-size="%.2f">%s</span>', sz, value)
        }
      ),
      towin_display = colDef(
        name = "", class = "cell-towin",
        style = list(color = "#3fb950")
      ),
      is_placed = colDef(
        name = "", class = "cell-action", html = TRUE,
        cell = function(value, index) {
          row <- table_data[index, ]
          status <- placed_status_lookup[row$bet_hash]
          ticket <- placed_ticket_lookup[row$bet_hash]
          placed_actual <- placed_actual_lookup[row$bet_hash]
          data_attrs <- sprintf(
            'data-hash="%s" data-game-id="%s" data-home="%s" data-away="%s" data-time="%s" data-market="%s" data-bet-on="%s" data-line="%s" data-prob="%s" data-ev="%s" data-size="%s" data-odds="%s" data-book="%s" data-actual="%s" data-fill-status="%s"',
            row$bet_hash, row$id, row$home_team, row$away_team,
            as.character(row$pt_start_time), row$market, row$bet_on,
            ifelse(is.na(row$line), "", row$line),
            row$prob, row$ev, row$bet_size, row$odds, row$bookmaker_key,
            ifelse(is.na(placed_actual), "", placed_actual),
            row$fill_status
          )
          # Partial fill: show the Partial -$X button so user can update actual_size
          if (!is.na(row$fill_status) && row$fill_status == "partial") {
            diff_label <- sprintf("-$%.0f", row$fill_diff)
            return(sprintf(
              '<button class="btn-partial" onclick="updateBet(this)" %s>Partial %s</button>',
              data_attrs, diff_label))
          }
          if (!is.na(status) && status == "placed") {
            label <- if (!is.na(ticket) && nchar(ticket) > 0)
              sprintf('placed &middot; #%s', htmltools::htmlEscape(ticket))
              else "placed"
            return(sprintf('<span class="placed-bet-label" %s>%s</span>',
                           data_attrs, label))
          }
          if (!is.na(status) && status %in% c("price_moved","rejected","auth_error","network_error","orphaned")) {
            short <- switch(status,
              price_moved   = "drift",
              rejected      = "rejected",
              auth_error    = "auth err",
              network_error = "net err",
              orphaned      = "orphan",
              status)
            return(sprintf(
              '<span class="pill error" %s>%s</span><button class="btn-place" onclick="placeBet(this)" %s>Retry</button><button class="btn-log" onclick="logBet(this)" %s>Log</button>',
              data_attrs, short, data_attrs, data_attrs))
          }
          # Default: not yet placed
          supported_place <- row$bookmaker_key %in% c("wagerzon","hoop88","bfa","betonlineag")
          place_btn <- if (supported_place) {
            sprintf('<button class="btn-place" onclick="placeBet(this)" %s>Place</button>', data_attrs)
          } else {
            sprintf('<button class="btn-place" disabled title="manual log only for this book" %s>Place</button>', data_attrs)
          }
          log_btn <- sprintf('<button class="btn-log" onclick="logBet(this)" %s>Log</button>', data_attrs)
          paste0(place_btn, " ", log_btn)
        }
      )
    ),
    theme = reactableTheme(
      backgroundColor = "transparent",
      borderColor = "transparent",
      stripedColor = "transparent",
      highlightColor = "transparent"
    )
  )
}

# =============================================================================
# HTML REPORT
# =============================================================================

create_report <- function(bets_table, placed_table, stats, timestamp, filter_options_json,
                          parlays_table = NULL, placed_parlays_table = NULL, parlay_opps = tibble(),
                          parlay_filter_options_json = "{}",
                          trifectas_table = NULL, trifecta_opps = tibble()) {
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
          /* Two-row header: .header-row-top (title + Refresh) on top,
             .header-row-accounts (Wagerzon pills) below. Each row owns
             its own internal flex; the parent stacks them vertically. */
          display: flex;
          flex-direction: column;
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

        /* Manual-log button: secondary tone (blue/grey) to avoid competing
           visually with the green Place button. Same size for a clean
           side-by-side row. */
        .btn-log {
          background: transparent;
          border: 1px solid #30363d;
          color: #8b949e;
          padding: 5px 10px;
          border-radius: 6px;
          font-size: 0.75rem;
          font-weight: 500;
          cursor: pointer;
          margin-left: 4px;
          transition: all 0.15s;
        }

        .btn-log:hover {
          border-color: #58a6ff;
          color: #58a6ff;
        }

        .combo-note { color: #8b949e; font-size: 0.7rem; font-style: italic; }

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


        /* === Wagerzon multi-account header pill row (Phase 7 layout merge) ===
           Lives inside .header as a second row. Replaces the old
           full-bleed wz-account-bar that sat above .container. */
        .header-row-top {
          display: flex;
          justify-content: space-between;
          align-items: center;
          padding-bottom: 10px;
        }
        .header-row-accounts {
          display: flex;
          align-items: center;
          gap: 8px;
          padding: 10px 0 4px 0;
        }
        .wz-pills {
          display: flex;
          gap: 6px;
          align-items: center;
        }
        .header-label {
          font-size: 11px;
          color: #8b949e;
          text-transform: uppercase;
          letter-spacing: 0.5px;
          margin-right: 2px;
        }
        .wz-pill {
          padding: 4px 10px;
          border-radius: 14px;
          background: #21262d;
          border: 1px solid #30363d;
          color: #c9d1d9;
          font-size: 13px;
          cursor: pointer;
          user-select: none;
          transition: border-color 0.12s, background 0.12s;
        }
        .wz-pill:hover { border-color: #58a6ff; }
        .wz-pill.selected {
          background: #1f6feb;
          border-color: #1f6feb;
          color: #ffffff;
          font-weight: 600;
        }
        .wz-pill.stale {
          background: #3a1d1d;
          border-color: #4a2a2a;
          color: #ffa198;
        }
        .wz-pill.empty {
          background: transparent;
          border-style: dashed;
          color: #6e7681;
          cursor: default;
        }
        .wz-pill.empty:hover { border-color: #30363d; }
        .wz-icon-btn {
          background: transparent;
          border: 1px solid #30363d;
          color: #8b949e;
          width: 28px;
          height: 28px;
          border-radius: 6px;
          cursor: pointer;
          display: inline-flex;
          align-items: center;
          justify-content: center;
          font-size: 14px;
          padding: 0;
        }
        .wz-icon-btn:hover { color: #c9d1d9; border-color: #58a6ff; }

        /* Parlay tab — books strip (M / DK / FD / PX / NV / Cons pill row) */
        .books-strip {
          display: flex;
          flex-wrap: wrap;
          gap: 6px;
          align-items: center;
        }
        .pill {
          padding: 2px 6px;
          border-radius: 3px;
          font-family: monospace;
          font-size: 10px;
          color: #8b949e;
          background: #21262d;
          white-space: nowrap;
        }
        .pill.model {
          background: #1f3a2c;
          color: #7ee787;
        }
        .pill.book {
          /* default styling; explicit class for selector clarity */
        }
        .pill.cons {
          background: #1c2738;
          border-left: 2px solid #58a6ff;
          color: #79c0ff;
        }
        .pill.dim {
          opacity: 0.4;
        }
        .pill.error {
          background: #3d1c1c;
          color: #f85149;
          border: 1px solid #6e2d2d;
        }
        /* Muted label shown after auto-placement succeeds */
        .placed-parlay-label {
          font-size: 0.72rem;
          color: #8b949e;
          font-family: monospace;
        }

        /* === Parlay tab card layout (scoped — singles tab is unaffected) === */
        /* Flatten the reactable table into a stack of cards. */
        #parlays-table-container .rt-table   { display: block; }
        #parlays-table-container .rt-thead   { display: none; }
        #parlays-table-container .rt-tbody   { display: block; }
        #parlays-table-container .rt-tr-group { display: block; }

        /* Each row becomes a card; cells flex inside so DOM order = visual order. */
        #parlays-table-container .rt-tr {
          display: flex;
          flex-wrap: wrap;
          align-items: baseline;
          gap: 4px;
          background: #161b22;
          border: 1px solid #21262d;
          border-radius: 6px;
          padding: 14px 14px 12px 14px;
          margin-bottom: 10px;
          position: relative;
        }

        #parlays-table-container .rt-td {
          display: block;
          border: 0;
          padding: 0;
          white-space: normal;
          font-size: 14px;
          color: #c9d1d9;
        }

        /* Full-width rows inside the flex card */
        #parlays-table-container .rt-td.cell-game,
        #parlays-table-container .rt-td.cell-legs,
        #parlays-table-container .rt-td.cell-books {
          flex-basis: 100%;
        }

        /* Sel checkbox: top-right corner */
        #parlays-table-container .rt-td.cell-sel {
          position: absolute;
          top: 12px;
          right: 12px;
        }
        #parlays-table-container .combo-select {
          width: 18px;
          height: 18px;
          cursor: pointer;
        }

        /* Game + folded time line */
        #parlays-table-container .rt-td.cell-game {
          font-size: 15px;
          font-weight: 500;
          padding-right: 36px;
        }

        /* Legs */
        #parlays-table-container .rt-td.cell-legs {
          font-size: 14px;
          margin-top: 2px;
          margin-bottom: 10px;
        }

        /* Books pill row */
        #parlays-table-container .rt-td.cell-books {
          margin-bottom: 12px;
        }

        /* Metadata strip — Fair / WZ / Size / To Win flow inline as flex items.
           Parent gap: 4px handles base spacing; margin-right adds breathing room. */
        #parlays-table-container .rt-td.cell-fair,
        #parlays-table-container .rt-td.cell-wz,
        #parlays-table-container .rt-td.cell-size,
        #parlays-table-container .rt-td.cell-towin {
          display: inline-flex;
          align-items: baseline;
          gap: 4px;
          margin-right: 10px;
          font-size: 13px;
          font-family: monospace;
        }

        #parlays-table-container .rt-td.cell-fair::before  { content: "Fair";   color: #8b949e; font-size: 12px; font-family: -apple-system, system-ui, sans-serif; }
        #parlays-table-container .rt-td.cell-wz::before    { content: "WZ";     color: #8b949e; font-size: 12px; font-family: -apple-system, system-ui, sans-serif; }
        #parlays-table-container .rt-td.cell-size::before  { content: "Size";   color: #8b949e; font-size: 12px; font-family: -apple-system, system-ui, sans-serif; }
        #parlays-table-container .rt-td.cell-towin::before { content: "To Win"; color: #8b949e; font-size: 12px; font-family: -apple-system, system-ui, sans-serif; }

        /* Edge — pushed to right via margin-left:auto on first of two right items.
           Action follows in DOM order, so visually:
           [...metadata strip...]                       +12.3%  [ Place ] */
        #parlays-table-container .rt-td.cell-edge {
          margin-left: auto;
          font-size: 15px;
          font-weight: 600;
        }
        #parlays-table-container .rt-td.cell-action {
          margin-left: 12px;
        }

        /* Conditional-Kelly residual note — full-width line below metadata strip */
        #parlays-table-container .combo-note {
          display: block;
          width: 100%;
          margin-top: 4px;
          color: #8b949e;
          font-size: 12px;
          font-style: italic;
        }

        /* Pill upsizing (was 10px in the books-strip round) */
        #parlays-table-container .pill {
          font-size: 13px;
          padding: 3px 8px;
        }

        /* === Responsive fix === */
        /* Reactable injects inline min-width / width / flex-grow on every cell
           derived from the colDef minWidth values (sum ~1340px for the parlay
           table). Inline styles beat external CSS, so the cards stayed 1340px
           wide and overflowed narrow viewports. Override with !important on
           the parlay table so cards fluidly shrink to whatever container
           width is available, and add overflow-x:auto to the placed-parlays
           summary section so it horizontally scrolls inside its own container
           at narrow widths instead of pushing the whole page wider. */
        #parlays-table-container .rt-tbody,
        #parlays-table-container .rt-tr-group,
        #parlays-table-container .rt-tr {
          min-width: 0 !important;
          width: auto !important;
        }
        #parlays-table-container .rt-td {
          min-width: 0 !important;
          width: auto !important;
          flex: 0 1 auto !important;
        }
        #parlays-table-container .rt-td.cell-game,
        #parlays-table-container .rt-td.cell-legs,
        #parlays-table-container .rt-td.cell-books {
          flex-basis: 100% !important;
        }

        /* Visual order in the card — independent of dataframe column order
           reactable preserves. Reads top→bottom: game → legs → books →
           Fair → WZ → Size → To Win → Edge → Action. */
        #parlays-table-container .rt-td.cell-game    { order: 1; }
        #parlays-table-container .rt-td.cell-legs    { order: 2; }
        #parlays-table-container .rt-td.cell-books   { order: 3; }
        #parlays-table-container .rt-td.cell-fair    { order: 4; }
        #parlays-table-container .rt-td.cell-wz      { order: 5; }
        #parlays-table-container .rt-td.cell-size    { order: 6; }
        #parlays-table-container .rt-td.cell-towin   { order: 7; }
        #parlays-table-container .rt-td.cell-edge    { order: 8; }
        #parlays-table-container .rt-td.cell-action  { order: 9; }

        /* Hide any auto-rendered cell that has no explicit cell-*
           class (e.g. blended_prob, which is in the dataframe with no colDef
           entry — would otherwise render as a stray "0.188" cell). */
        #parlays-table-container .rt-td:not(.cell-sel):not(.cell-game):not(.cell-legs):not(.cell-fair):not(.cell-books):not(.cell-wz):not(.cell-edge):not(.cell-size):not(.cell-towin):not(.cell-action) {
          display: none !important;
        }

        /* Placed-parlays summary: keep table layout but contain horizontal
           overflow inside the section instead of pushing page wider. */
        #placed-parlays-section { overflow-x: auto; }
        #placed-parlays-section .rt-thead,
        #placed-parlays-section .rt-tbody,
        #placed-parlays-section .rt-tr-group,
        #placed-parlays-section .rt-tr,
        #placed-parlays-section .rt-tr-header {
          min-width: 0 !important;
        }

        /* Bigger fonts (per design feedback): primary content 15-16px,
           secondary metadata 13-14px. Pill row scales with content. */
        #parlays-table-container .rt-td        { font-size: 15px !important; }
        #parlays-table-container .rt-td.cell-game,
        #parlays-table-container .rt-td.cell-edge { font-size: 16px !important; }
        #parlays-table-container .pill         { font-size: 14px !important; padding: 4px 9px !important; }
        #parlays-table-container .rt-td.cell-fair,
        #parlays-table-container .rt-td.cell-wz,
        #parlays-table-container .rt-td.cell-size,
        #parlays-table-container .rt-td.cell-towin { font-size: 14px !important; }
        #parlays-table-container .rt-td.cell-fair::before,
        #parlays-table-container .rt-td.cell-wz::before,
        #parlays-table-container .rt-td.cell-size::before,
        #parlays-table-container .rt-td.cell-towin::before { font-size: 13px !important; }

        /* === Bets tab odds-screen card layout === */

        /* Each bet renders as a card via reactable row + display:flex on the row container. */
        #bets-table-container .rt-table   { display: block; }
        #bets-table-container .rt-thead   { display: none; }
        #bets-table-container .rt-tbody   { display: block; }
        #bets-table-container .rt-tr-group { display: block; }
        #bets-table-container .rt-tr {
          display: flex;
          flex-wrap: wrap;
          background: #1c2128;
          border: 1px solid #373e47;
          border-radius: 10px;
          padding: 12px 16px;
          margin-bottom: 10px;
          align-items: center;
          position: relative;
        }
        #bets-table-container .rt-td {
          min-width: 0 !important;
          flex: 0 1 auto;
        }

        /* !important is required to override the reactable inline
           style="flex: 100 0 auto" on every .rt-td cell. Without it the
           cells stay at content-width and pills wrap vertically into
           narrow columns. The parlays tab already does this for its
           equivalent cells (see #parlays-table-container rules below). */
        #bets-table-container .rt-td.cell-game,
        #bets-table-container .rt-td.cell-market,
        #bets-table-container .rt-td.cell-pickside,
        #bets-table-container .rt-td.cell-otherside {
          flex-basis: 100% !important;
          width: 100% !important;
          padding: 2px 0;
        }

        #bets-table-container .cell-game {
          order: 1;
          color: #c9d1d9; font-weight: 600;
        }
        #bets-table-container .cell-market {
          order: 2;
          color: #c9d1d9; font-weight: 600;
          border-bottom: 1px solid #2d333b;
          padding-bottom: 6px;
          margin-bottom: 6px;
        }
        #bets-table-container .cell-pickside  { order: 3; }
        #bets-table-container .cell-otherside { order: 4; }

        /* Pill row: each row is a flex container of fixed-width pills */
        #bets-table-container .side-row {
          display: flex; flex-wrap: wrap;
          gap: 6px; margin: 4px 0;
          align-items: center;
        }
        #bets-table-container .side-label {
          flex: 0 0 78px; min-width: 78px;
          color: #c9d1d9; font-weight: 600; font-size: 12px;
        }
        #bets-table-container .pill {
          flex: 0 0 78px; min-width: 78px;
          box-sizing: border-box;
          display: inline-flex; flex-direction: column;
          align-items: center; gap: 1px;
          padding: 5px 6px;
          border-radius: 6px;
          background: #22272e;
          border: 1px solid #373e47;
          color: #c9d1d9;
          font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
          font-size: 12px;
        }
        #bets-table-container .pill .book {
          color: #8b949e; font-size: 9px;
          text-transform: uppercase; letter-spacing: 0.5px;
        }
        #bets-table-container .pill .line-tag {
          color: #f0883e; font-size: 9px; font-weight: 600;
        }
        #bets-table-container .pill.muted {
          color: #6e7681;
          background: transparent;
          border-style: dashed;
        }
        #bets-table-container .pill.mismatched {
          border-color: #5a3d1a;
          background: #2a2317;
        }
        #bets-table-container .pill.pick {
          background: #15321f;
          border-color: #3fb950;
          color: #56d364;
          font-weight: 700;
        }
        #bets-table-container .pill.pick .book { color: #56d364; }

        /* Metadata strip */
        #bets-table-container .cell-m,
        #bets-table-container .cell-pick,
        #bets-table-container .cell-ev,
        #bets-table-container .cell-size,
        #bets-table-container .cell-towin,
        #bets-table-container .cell-action {
          order: 5;
          display: inline-flex; align-items: center;
          gap: 5px;
          padding: 0 18px 0 0;
          font-size: 14px;
        }
        #bets-table-container .cell-m::before    { content: "M";       color: #8b949e; font-size: 11px; text-transform: uppercase; letter-spacing: 0.5px; margin-right: 4px; }
        #bets-table-container .cell-pick::before { content: "Pick";    color: #8b949e; font-size: 11px; text-transform: uppercase; letter-spacing: 0.5px; margin-right: 4px; }
        #bets-table-container .cell-ev::before   { content: "EV";      color: #8b949e; font-size: 11px; text-transform: uppercase; letter-spacing: 0.5px; margin-right: 4px; }
        #bets-table-container .cell-size::before { content: "Size";    color: #8b949e; font-size: 11px; text-transform: uppercase; letter-spacing: 0.5px; margin-right: 4px; }
        #bets-table-container .cell-towin::before{ content: "To Win";  color: #8b949e; font-size: 11px; text-transform: uppercase; letter-spacing: 0.5px; margin-right: 4px; }
        #bets-table-container .cell-action       { margin-left: auto; padding-right: 0; }

        /* Correlation dot — corner badge on the card */
        #bets-table-container .corr-badge {
          position: absolute;
          top: 8px; right: 12px;
          color: #f0883e; font-size: 18px;
          cursor: help;
        }

      '))
    ),

    tags$body(
      tags$div(class = "container",
        # Header — two rows.
        # Row 1: title + subtitle (left) | "Refresh" data button (right).
        # Row 2: "Placing on" caption + clickable Wagerzon account pills +
        #        balance-refresh icon. Replaces the old full-bleed
        #        wz-account-bar that lived above .container.
        tags$div(class = "header",
          tags$div(class = "header-row-top",
            tags$div(
              tags$h1("MLB Answer Key Dashboard"),
              tags$div(class = "subtitle", paste("Updated", timestamp))
            ),
            tags$button(class = "refresh-btn", onclick = "refreshData()", "Refresh")
          ),
          tags$div(class = "header-row-accounts", id = "wz-account-row",
            tags$span(class = "header-label", "Placing on"),
            tags$div(id = "wz-account-pills", class = "wz-pills"),
            tags$button(
              id = "wz-refresh-btn", type = "button", class = "wz-icon-btn",
              title = "Refresh balances",
              HTML("&#x21bb;")
            )
          )
        ),

        # Tab Navigation
        tags$div(class = "tab-bar",
          tags$button(class = "tab-btn active", id = "tab-btn-bets", onclick = "switchTab(\'bets\')", "Bets"),
          tags$button(class = "tab-btn", id = "tab-btn-parlays", onclick = "switchTab(\'parlays\')", "Parlays"),
          tags$button(class = "tab-btn", id = "tab-btn-trifectas", onclick = "switchTab(\'trifectas\')",
                      sprintf("Trifectas (%d)", nrow(trifecta_opps)))
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

          # Combined Parlay banner — hidden until 2 rows checked in the table
          HTML('
<div id="combined-parlay-banner" style="display:none; padding:10px 14px; background:#1f3a5f; border-left:3px solid #4a9eff; margin-bottom:10px; border-radius:4px; color:#c9d1d9; font-family:sans-serif; font-size:13px;">
  <span id="combo-banner-status">Pricing combined ticket…</span>
  <button id="combo-place-btn" onclick="placeCombinedParlay()" style="display:none; float:right; background:#4a9eff; color:#fff; border:none; padding:6px 14px; border-radius:4px; cursor:pointer; font-weight:600;">Place combined →</button>
</div>
'),

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

        ,

        # ============ TRIFECTAS TAB ============
        tags$div(id = "tab-trifectas", class = "tab-content", style = "display: none;",

          # Trifecta Sizing Controls (mirrors parlay sizing)
          tags$div(class = "sizing-controls",
            tags$div(class = "sizing-group",
              tags$span(class = "sizing-label", "Trifecta Bankroll ($)"),
              tags$input(id = "trifecta-bankroll-input", class = "sizing-input", type = "number",
                value = "100", min = "1", step = "10")
            ),
            tags$div(class = "sizing-group",
              tags$span(class = "sizing-label", "Trifecta Kelly"),
              tags$input(id = "trifecta-kelly-input", class = "sizing-input", type = "number",
                value = "0.10", min = "0.01", max = "1", step = "0.05")
            ),
            tags$div(class = "sizing-group",
              tags$span(class = "sizing-label", "Min Edge (%)"),
              tags$input(id = "trifecta-min-edge-input", class = "sizing-input", type = "number",
                value = "5", min = "0", step = "1",
                onchange = "filterTrifectasByEdge()", oninput = "filterTrifectasByEdge()")
            ),
            tags$button(class = "apply-sizing-btn", onclick = "applyTrifectaSizing()", "Apply")
          ),

          # Trifecta Opportunities
          if (!is.null(trifectas_table)) {
            tagList(
              tags$div(class = "section-header",
                tags$span("Trifecta Opportunities"),
                tags$span(id = "trifecta-count",
                  style = "margin-left: 8px; color: #8b949e; font-size: 0.8rem;",
                  sprintf("(%d)", nrow(trifecta_opps)))
              ),
              tags$div(class = "table-container", id = "trifectas-table-container", trifectas_table)
            )
          }

        ) # end tab-trifectas
      ),

      # JavaScript
      tags$script(HTML('
        // ============ TAB SWITCHING ============
        function switchTab(tab) {
          document.getElementById("tab-bets").style.display = tab === "bets" ? "" : "none";
          document.getElementById("tab-parlays").style.display = tab === "parlays" ? "" : "none";
          document.getElementById("tab-trifectas").style.display = tab === "trifectas" ? "" : "none";
          document.getElementById("tab-btn-bets").className = "tab-btn" + (tab === "bets" ? " active" : "");
          document.getElementById("tab-btn-parlays").className = "tab-btn" + (tab === "parlays" ? " active" : "");
          document.getElementById("tab-btn-trifectas").className = "tab-btn" + (tab === "trifectas" ? " active" : "");
          // Apply edge filters when switching to their respective tabs
          if (tab === "parlays") filterParlaysByEdge();
          if (tab === "trifectas") filterTrifectasByEdge();
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

          // Restore Parlays tab if URL has #parlays (used by the
          // hot-swap fallback path in placeCombinedParlay so a full reload
          // doesn\'t dump the user back on the Bets tab).
          if (window.location.hash === "#parlays" &&
              typeof switchTab === "function") {
            switchTab("parlays");
          }

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

            // In the new card layout the second-to-last cell is pickside_html
            // (the full pill row), not a single book name as in the legacy
            // table layout. Read the pick book from the Place/Log button
            // data-book attribute, which carries the bet bookmaker_key
            // (e.g. wagerzon, draftkings).
            const bookBtn = row.querySelector("button[data-book]");
            if (bookBtn) bookText = bookBtn.getAttribute("data-book").trim();

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
          var data = btn.dataset;
          var book = data.book;
          var account = window.WZ_SELECTED_ACCOUNT || null;

          var body = {
            bet_hash:         data.hash,
            bookmaker_key:    book,
            account:          account,
            bet_on:           data.betOn,
            line:             (data.line === \'\' || data.line === undefined) ? null : parseFloat(data.line),
            market:           data.market,
            american_odds:    parseInt(data.odds, 10),
            actual_size:      parseFloat(data.size),
            kelly_bet:        parseFloat(data.size),
            wz_odds_at_place: parseInt(data.odds, 10),
            game_id:          data.gameId,
            home_team:        data.home,
            away_team:        data.away,
            game_time:        data.time,
            model_prob:       data.prob ? parseFloat(data.prob) : 0.0,
            model_ev:         data.ev   ? parseFloat(data.ev)   : 0.0
          };

          if (book === \'wagerzon\' && !account) {
            showToast(\'No Wagerzon account selected — pick one in the header pills\', \'error\');
            return;
          }

          btn.disabled = true;
          var originalLabel = btn.textContent;
          btn.textContent = \'Placing...\';

          fetch(\'/api/place-bet\', {
            method: \'POST\',
            headers: {\'Content-Type\': \'application/json\'},
            body: JSON.stringify(body)
          })
          .then(function(r) {
            return r.json().then(function(j) { return {ok: r.ok, status: r.status, body: j}; });
          })
          .then(function(resp) {
            var ok = resp.ok, status = resp.status, result = resp.body;
            if (status === 409) {
              showToast(\'Bet already in flight: \' + (result.error || result.status), \'warning\');
              btn.disabled = false; btn.textContent = originalLabel;
              return;
            }
            if (result.status === \'placed\') {
              var ticket = result.ticket_number ? \' #\' + result.ticket_number : \'\';
              showToast(\'Placed at \' + book + ticket, \'success\');
              setTimeout(function() { location.reload(); }, 800);
              return;
            }
            if (result.status === \'playwright_launched\') {
              showToast(result.message || \'Browser launching...\', \'info\');
              btn.disabled = false; btn.textContent = originalLabel;
              return;
            }
            if (result.status === \'price_moved\') {
              showToast(\'Price moved — bet not placed\', \'warning\');
              btn.disabled = false; btn.textContent = originalLabel;
              return;
            }
            if (result.error) {
              showToast(result.error, \'error\');
            } else if (result.status) {
              showToast(\'Status: \' + result.status + (result.error_msg ? \' — \' + result.error_msg : \'\'), \'error\');
            } else {
              showToast(\'Unknown response\', \'error\');
            }
            btn.disabled = false; btn.textContent = originalLabel;
          })
          .catch(function(e) {
            showToast(\'Network error: \' + e.message, \'error\');
            btn.disabled = false; btn.textContent = originalLabel;
          });
        }

        function logBet(btn) {
          var data = btn.dataset;
          var account = window.WZ_SELECTED_ACCOUNT || null;
          var body = {
            bet_hash:        data.hash,
            account:         account,
            bookmaker_key:   data.book,
            bet_on:          data.betOn,
            line:            (data.line === \'\' || data.line === undefined) ? null : parseFloat(data.line),
            market:          data.market,
            american_odds:   parseInt(data.odds, 10),
            actual_size:     parseFloat(data.size),
            kelly_bet:       parseFloat(data.size),
            game_id:         data.gameId,
            home_team:       data.home,
            away_team:       data.away,
            game_time:       data.time,
            model_prob:      data.prob ? parseFloat(data.prob) : 0.0,
            model_ev:        data.ev   ? parseFloat(data.ev)   : 0.0
          };

          btn.disabled = true;
          btn.textContent = \'Logging...\';

          fetch(\'/api/log-bet\', {
            method: \'POST\',
            headers: {\'Content-Type\': \'application/json\'},
            body: JSON.stringify(body)
          })
          .then(function(r) { return r.json(); })
          .then(function(result) {
            if (result.success !== false) {
              setTimeout(function() { location.reload(); }, 400);
            } else {
              showToast(\'Log failed: \' + (result.error || \'unknown\'), \'error\');
              btn.disabled = false; btn.textContent = \'Log\';
            }
          })
          .catch(function(e) {
            showToast(\'Log failed: \' + e.message, \'error\');
            btn.disabled = false; btn.textContent = \'Log\';
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
            // Look up the data-bearing element. Three possible shapes per row:
            //   <button data-edge=... class="btn-place">       — un-placed
            //   <button data-edge=... class="btn-placed">      — manually placed (legacy)
            //   <span   data-edge=... class="placed-parlay-label">  — auto-placed
            //   <span   data-edge=... class="pill error">      — auto-place failed
            // Any element with data-edge is the marker.
            var marker = row.querySelector("[data-edge]");
            var edge = marker ? parseFloat(marker.dataset.edge) : 0;
            var size = marker ? parseFloat(marker.dataset.size) : 0;
            var gameText = marker ? (marker.dataset.away + " @ " + marker.dataset.home) : "";
            // Two-bucket classification matching the status dropdown options
            // ("Not Placed" / "Placed"). A failed placement (red pill) is
            // conceptually NOT-placed (the bet didn\'t go through; the user
            // can retry), so it falls into "Not Placed" rather than a third
            // bucket the dropdown can\'t address.
            var statusLabel;
            if (marker && (marker.classList.contains("btn-placed") ||
                           marker.classList.contains("placed-parlay-label"))) {
              statusLabel = "Placed";
            } else {
              statusLabel = "Not Placed";
            }

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

        // ============ TRIFECTA EDGE FILTER ============
        function filterTrifectasByEdge() {
          var container = document.getElementById("trifectas-table-container");
          if (!container) return;
          var rows = container.querySelectorAll(".rt-tr-group");
          var minEdgeEl = document.getElementById("trifecta-min-edge-input");
          var minEdge = minEdgeEl ? (parseFloat(minEdgeEl.value) || 0) : 0;
          var visible = 0;
          rows.forEach(function(row) {
            var marker = row.querySelector("[data-edge]");
            var edge = marker ? parseFloat(marker.dataset.edge) : 0;
            var show = edge >= minEdge;
            row.style.display = show ? "" : "none";
            if (show) visible++;
          });
          var countEl = document.getElementById("trifecta-count");
          if (countEl) {
            if (visible === rows.length) {
              countEl.textContent = "(" + rows.length + ")";
            } else {
              countEl.textContent = "(" + visible + " of " + rows.length + ")";
            }
          }
        }

        // ============ TRIFECTA SIZING ============
        function applyTrifectaSizing() {
          var bankroll = parseFloat(document.getElementById("trifecta-bankroll-input").value);
          var kelly = parseFloat(document.getElementById("trifecta-kelly-input").value);
          var minEdge = parseFloat(document.getElementById("trifecta-min-edge-input").value) || 0;
          if (!bankroll || bankroll <= 0 || !kelly || kelly <= 0 || kelly > 1) {
            showToast("Invalid sizing values", "error");
            return;
          }
          // trifecta_min_edge is stored as a fraction (0-1); UI is percent (0-100).
          Promise.all([
            fetch("/api/sizing-settings", {
              method: "POST",
              headers: { "Content-Type": "application/json" },
              body: JSON.stringify({ param: "trifecta_bankroll", value: bankroll })
            }),
            fetch("/api/sizing-settings", {
              method: "POST",
              headers: { "Content-Type": "application/json" },
              body: JSON.stringify({ param: "trifecta_kelly_mult", value: kelly })
            }),
            fetch("/api/sizing-settings", {
              method: "POST",
              headers: { "Content-Type": "application/json" },
              body: JSON.stringify({ param: "trifecta_min_edge", value: minEdge / 100 })
            })
          ]).then(function() {
            showToast("Trifecta settings saved. Refreshing...", "success");
            refreshData();
          });
        }

        // Load trifecta settings from server on page load (guard against missing elements)
        (function() {
          var tbi = document.getElementById("trifecta-bankroll-input");
          var tki = document.getElementById("trifecta-kelly-input");
          var tei = document.getElementById("trifecta-min-edge-input");
          if (!tbi || !tki) return;
          fetch("/api/sizing-settings")
            .then(function(r) { return r.json(); })
            .then(function(s) {
              if (s.trifecta_bankroll) tbi.value = s.trifecta_bankroll;
              if (s.trifecta_kelly_mult) tki.value = s.trifecta_kelly_mult;
              if (tei && s.trifecta_min_edge != null) {
                // Stored as fraction (0-1); display as integer percent.
                tei.value = (parseFloat(s.trifecta_min_edge) * 100).toFixed(0);
                filterTrifectasByEdge();
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

          // Build leg description: "TeamName -1.5 / O7.0"
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
          var hash = btn.dataset.hash;
          if (!hash) { showToast("Missing parlay_hash", "error"); return; }

          var originalText = btn.textContent;
          btn.disabled = true;
          btn.textContent = "Placing...";

          // Every placement must specify which Wagerzon account it lands on.
          // window.WZ_SELECTED_ACCOUNT (set by the header pill row) is the
          // source of truth; bail loudly if for some reason no account is
          // set (e.g. server returned no accounts at startup).
          if (!window.WZ_SELECTED_ACCOUNT) {
            showToast("No Wagerzon account selected", "error");
            btn.disabled = false;
            btn.textContent = originalText;
            return;
          }

          fetch("/api/place-parlay", {
            method: "POST",
            headers: { "Content-Type": "application/json" },
            body: JSON.stringify({
              parlay_hash: hash,
              account: window.WZ_SELECTED_ACCOUNT
            })
          })
            .then(function(r) { return r.json(); })
            .then(function(result) {
              btn.disabled = false;
              if (result.status === "placed") {
                // Real placement. The Action cell may contain BOTH [Place]
                // and [Log] buttons; replace the entire cell with the
                // placed label so neither is left behind.
                var ticket = result.ticket_number || "";
                var span = document.createElement(\'span\');
                span.className = "placed-parlay-label";
                // No ticket → just "placed" (cleaner than "placed · #")
                span.textContent = ticket ? ("placed \xb7 #" + ticket) : "placed";
                _replaceActionCell(btn, span);
                // Server returns the post-placement balance snapshot for the
                // account it landed on. Push it into the pill immediately
                // (no extra GET) so the user sees the debit reflected.
                if (result.balance_after && window.wzApplyBalanceAfter) {
                  window.wzApplyBalanceAfter(result.balance_after.label, result.balance_after);
                }
                showToast(ticket
                  ? ("Placed on " + window.WZ_SELECTED_ACCOUNT + " (#" + ticket + ")")
                  : ("Placed on " + window.WZ_SELECTED_ACCOUNT),
                  "success");
              } else if (result.transient) {
                // Transient: the bet did not go through (network glitch, line
                // moved, balance / limit refusal, auth blip). The breadcrumb
                // got deleted server-side; restore the [Place] button text so
                // the user can immediately retry. Toast carries the reason.
                btn.textContent = originalText;
                var msg = result.error_msg || result.status || "Failed";
                showToast("Place failed: " + msg + " — try again", "error");
              } else {
                // Persistent failure (orphaned — real money is at WZ but
                // local audit failed). Swap to red pill so the row is
                // visibly different from a clean retryable miss.
                var msg = result.error_msg || result.status || "Failed";
                var label = result.short_label || result.status || "failed";
                var span = document.createElement(\'span\');
                span.className = "pill error";
                span.textContent = label;
                span.title = msg;
                _replaceActionCell(btn, span);
                showToast("Place failed: " + msg, "error");
              }
            })
            .catch(function(err) {
              btn.disabled = false;
              btn.textContent = originalText;
              showToast("Network error: " + err, "error");
            });
        }

        // Replace the parlay row\'s entire Action cell content with `node`.
        // Walk up from the clicked button to the reactable cell wrapper so
        // both buttons ([Place] AND [Log]) get cleared together — preserves
        // the cell\'s data attrs implicitly because the new node carries them.
        function _replaceActionCell(btn, node) {
          var cell = btn.closest(\'.rt-td\') || btn.parentNode;
          while (cell.firstChild) cell.removeChild(cell.firstChild);
          cell.appendChild(node);
        }

        // Manual-log path: user has placed the bet themselves (on Wagerzon\'s
        // UI or anywhere) and just wants the dashboard to track it. Opens a
        // small modal asking for the actual stake (defaults to the row\'s
        // recommended size), then POSTs to /api/log-parlay. NO Wagerzon API
        // call. Result row in placed_parlays carries status=\'placed\' with
        // ticket_number=NULL — the cell renders as just \'placed\' (no #ticket
        // suffix).
        function logParlay(btn) {
          var hash = btn.dataset.hash;
          if (!hash) { showToast("Missing parlay_hash", "error"); return; }

          var combo = btn.dataset.combo || "";
          var spread = btn.dataset.spread || "";
          var total = btn.dataset.total || "";
          var wzOdds = btn.dataset.wzOdds || "";
          var edge = btn.dataset.edge || "";
          var defaultStake = Math.round(parseFloat(btn.dataset.size) || 0);

          var overlay = document.createElement(\'div\');
          overlay.className = \'modal-overlay\';
          overlay.innerHTML =
            \'<div class="modal-box">\' +
              \'<div class="modal-title">Log manual placement</div>\' +
              \'<div class="modal-detail">Combo: <span>\' + escapeHtml(combo) + \'</span></div>\' +
              \'<div class="modal-detail">Spread: <span>\' + escapeHtml(spread) + \'</span> | Total: <span>\' + escapeHtml(total) + \'</span></div>\' +
              \'<div class="modal-detail">WZ Odds: <span>\' + escapeHtml(wzOdds) + \'</span> | Edge: <span>+\' + escapeHtml(edge) + \'%</span></div>\' +
              \'<div class="modal-recommended">Recommended: $\' + defaultStake + \'</div>\' +
              \'<div class="modal-input-group">\' +
                \'<label class="modal-input-label">Stake ($)</label>\' +
                \'<input type="number" class="modal-input" value="\' + defaultStake + \'" step="1" min="1">\' +
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

            fetch("/api/log-parlay", {
              method: "POST",
              headers: { "Content-Type": "application/json" },
              body: JSON.stringify({ parlay_hash: hash, actual_size: amount })
            })
              .then(function(r) { return r.json(); })
              .then(function(result) {
                overlay.remove();
                if (result.status === "placed") {
                  var span = document.createElement(\'span\');
                  span.className = "placed-parlay-label";
                  // No ticket for manual logs — just "placed".
                  span.textContent = "placed";
                  // Carry the same data attrs forward so the row filter still
                  // classifies it as Placed and reads the right edge/size.
                  span.dataset.hash = hash;
                  span.dataset.home = btn.dataset.home || "";
                  span.dataset.away = btn.dataset.away || "";
                  span.dataset.edge = btn.dataset.edge || "";
                  span.dataset.size = String(amount);
                  _replaceActionCell(btn, span);
                  showToast("Logged $" + amount, "success");
                } else {
                  showToast("Log failed: " + (result.error_msg || "unknown"), "error");
                }
              })
              .catch(function(err) {
                overlay.remove();
                showToast("Network error: " + err, "error");
              });
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

        async function placeTrifecta(btn) {
          const hash  = btn.dataset.trifectaHash;
          const wager = parseFloat(btn.dataset.actualWager);
          if (!hash) { console.error(\'placeTrifecta: missing data-trifecta-hash\'); return; }
          try {
            const r = await fetch(\'/api/place-trifecta\', {
              method: \'POST\',
              headers: { \'Content-Type\': \'application/json\' },
              body: JSON.stringify({ trifecta_hash: hash, actual_wager: isNaN(wager) ? null : wager })
            });
            const body = await r.json().catch(() => ({}));
            if (r.ok && body.success !== false) {
              btn.textContent = \'Placed\';
              btn.className = \'btn-placed\';
              btn.removeAttribute(\'data-actual-wager\');
              btn.onclick = () => removeTrifecta(btn);
            } else {
              alert(\'Place failed: \' + (body.error || r.statusText));
            }
          } catch (e) {
            alert(\'Network error: \' + e.message);
          }
        }

        async function removeTrifecta(btn) {
          const hash = btn.dataset.trifectaHash;
          if (!hash) { console.error(\'removeTrifecta: missing data-trifecta-hash\'); return; }
          try {
            const r = await fetch(\'/api/remove-trifecta\', {
              method: \'POST\',
              headers: { \'Content-Type\': \'application/json\' },
              body: JSON.stringify({ trifecta_hash: hash })
            });
            const body = await r.json().catch(() => ({}));
            if (r.ok && body.success !== false) {
              btn.textContent = \'Place\';
              btn.className = \'btn-place\';
              btn.onclick = () => placeTrifecta(btn);
            } else {
              alert(\'Remove failed: \' + (body.error || r.statusText));
            }
          } catch (e) {
            alert(\'Network error: \' + e.message);
          }
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

        // === Combined parlay state ===
        // At most 2 selected at a time, FIFO. Cleared on successful place.
        window.comboSelections = window.comboSelections || [];
        window.comboPricing = null;

        function onComboSelectChange(checkbox) {
          const hash = checkbox.dataset.hash;
          const gameId = checkbox.dataset.gameId;
          if (checkbox.checked) {
            // FIFO: if 2 already selected, uncheck the oldest
            if (window.comboSelections.length >= 2) {
              const oldest = window.comboSelections.shift();
              const oldEl = document.querySelector(`.combo-select[data-hash="${oldest.hash}"]`);
              if (oldEl) oldEl.checked = false;
            }
            window.comboSelections.push({ hash, gameId });
          } else {
            window.comboSelections = window.comboSelections.filter(s => s.hash !== hash);
          }
          refreshComboBanner();
        }

        async function refreshComboBanner() {
          const banner = document.getElementById("combined-parlay-banner");
          const status = document.getElementById("combo-banner-status");
          const placeBtn = document.getElementById("combo-place-btn");
          if (!banner || !status || !placeBtn) return;

          if (window.comboSelections.length !== 2) {
            banner.style.display = "none";
            return;
          }
          // Same-game guard — block at JS layer too so the user gets immediate feedback
          if (window.comboSelections[0].gameId === window.comboSelections[1].gameId) {
            banner.style.display = "block";
            status.textContent = "Same-game combos are already a single row — pick rows from different games.";
            placeBtn.style.display = "none";
            window.comboPricing = null;
            return;
          }

          banner.style.display = "block";
          status.textContent = "Pricing combined ticket…";
          placeBtn.style.display = "none";

          try {
            const res = await fetch("/api/price-combined-parlay", {
              method: "POST",
              headers: { "Content-Type": "application/json" },
              body: JSON.stringify({
                parlay_hash_a: window.comboSelections[0].hash,
                parlay_hash_b: window.comboSelections[1].hash,
              }),
            });
            const body = await res.json();
            if (!body.success) {
              status.textContent = "WZ rejected: " + (body.error || "unknown");
              window.comboPricing = null;
              return;
            }
            window.comboPricing = body;
            const edgePct = (body.joint_edge * 100).toFixed(1);
            const edgeColor = body.joint_edge >= 0 ? "#5fd97a" : "#e57373";
            const edgeSign = body.joint_edge >= 0 ? "+" : "";
            status.innerHTML =
              `<b>2 legs combined</b> · stake <b>$${body.kelly_stake}</b>` +
              ` · WZ pays <b>${body.wz_dec.toFixed(2)}</b>` +
              ` · edge <span style="color:${edgeColor}">${edgeSign}${edgePct}%</span>`;
            placeBtn.style.display = "inline-block";
          } catch (err) {
            status.textContent = "Pricing failed — uncheck and re-check a row to retry.";
            window.comboPricing = null;
          }
        }

        async function placeCombinedParlay() {
          if (!window.comboPricing || window.comboSelections.length !== 2) return;
          const placeBtn = document.getElementById("combo-place-btn");
          placeBtn.disabled = true;
          placeBtn.textContent = "Placing…";

          // Synthetic combo_hash — deterministic from sorted leg hashes
          const sels = window.comboSelections.map(s => s.hash).sort();
          const comboHash = "combo_" + sels.join("_");

          try {
            const res = await fetch("/api/place-combined-parlay", {
              method: "POST",
              headers: { "Content-Type": "application/json" },
              body: JSON.stringify({
                combo_hash: comboHash,
                parlay_hash_a: sels[0],
                parlay_hash_b: sels[1],
                wz_odds: Math.round((window.comboPricing.wz_dec - 1) * 100),
                kelly_bet: window.comboPricing.kelly_stake,
                actual_size: window.comboPricing.kelly_stake,
                combo_label: window.comboPricing.leg_a_combo + " + " + window.comboPricing.leg_b_combo,
              }),
            });
            const body = await res.json();
            if (!body.success) {
              alert("Failed to place combined parlay: " + (body.error || "unknown"));
              placeBtn.disabled = false;
              placeBtn.textContent = "Place combined →";
              return;
            }
            // Hot-swap the source-parlays table in place. /api/parlay-table-
            // fragment runs R with --parlay-fragment, returns the freshly-
            // rendered reactable widget HTML (with apply_combo_residuals
            // already applied to the two source rows). HTMLWidgets.staticRender()
            // re-initializes the React binding after innerHTML replacement.
            // No full page reload → tab state, scroll, and unrelated widgets
            // all persist.
            placeBtn.textContent = "Refreshing…";
            try {
              const fragRes = await fetch("/api/parlay-table-fragment");
              if (!fragRes.ok) throw new Error("fragment HTTP " + fragRes.status);
              const html = await fragRes.text();
              const container = document.getElementById("parlays-table-container");
              if (!container) throw new Error("parlays-table-container not in DOM");
              container.innerHTML = html;
              if (window.HTMLWidgets && typeof HTMLWidgets.staticRender === "function") {
                HTMLWidgets.staticRender();
              }
              // Reset combo selection state since the rows were re-rendered
              if (window.comboSelections) window.comboSelections.length = 0;
              window.comboPricing = null;
              const banner = document.getElementById("combo-banner");
              if (banner) banner.style.display = "none";
              showToast("Combined parlay placed", "success");
            } catch (swapErr) {
              // Fallback: full reload, preserving the parlays tab via URL hash.
              console.error("Hot-swap failed, falling back to reload:", swapErr);
              window.location.hash = "#parlays";
              window.location.reload();
              return;
            }
            placeBtn.disabled = false;
            placeBtn.textContent = "Place combined →";
          } catch (err) {
            alert("Network error placing combined parlay: " + err.message);
            placeBtn.disabled = false;
            placeBtn.textContent = "Place combined →";
          }
        }
      ')),

      tags$script(HTML(r"(
(function() {
  // Wagerzon multi-account header pill row controller (Phase 7).
  // Owns:
  //   * window.WZ_SELECTED_ACCOUNT (string, source-of-truth selection)
  //   * window.WZ_BALANCES         (label -> snapshot from /api/wagerzon/balances)
  //   * window.wzApplyBalanceAfter (called by placeParlay() after a successful
  //                                 placement to refresh the pill without a full GET)
  //   * window._wzRecomputeWarnings (sweeps every Place button with data-risk and
  //                                  populates the sibling .wz-insufficient-warning)
  var PILLS_ID  = 'wz-account-pills';
  var REFRESH_BTN_ID = 'wz-refresh-btn';
  var STALE_SECONDS = 60;   // suffix + .stale class only kick in past this

  window.WZ_SELECTED_ACCOUNT = null;
  window.WZ_BALANCES = {};   // label -> snapshot

  function fmtMoney(n) {
    if (n === null || n === undefined) return '—';  // em dash
    return '$' + Number(n).toLocaleString('en-US', {minimumFractionDigits:2, maximumFractionDigits:2});
  }

  function persistSelection(label) {
    return fetch('/api/wagerzon/last-used', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({label: label})
    });
  }

  function pillTextFor(label, snap) {
    // Healthy:        Label . $1,234.56
    // Errored, fresh: Label . -- (warning sign)
    // Errored, stale: Label . -- (warning sign) followed by (stale Nm ago)
    var text = label + ' · ' + fmtMoney(snap.available);
    if (snap.error) {
      text += ' ⚠';
      // Only label as "stale" when it has been at least a minute since the
      // last successful fetch. Avoids the contradictory "stale 0s ago" text
      // we used to render on every fresh-but-errored fetch.
      if (snap.stale_seconds >= STALE_SECONDS) {
        var sub = Math.floor(snap.stale_seconds / 60) + 'm';
        text += ' (stale ' + sub + ' ago)';
      }
    }
    return text;
  }

  function isStale(snap) {
    return snap.error && snap.stale_seconds >= STALE_SECONDS;
  }

  function renderEmpty(bar) {
    var pill = document.createElement('span');
    pill.className = 'wz-pill empty';
    pill.textContent = 'No Wagerzon accounts configured';
    bar.appendChild(pill);
  }

  function renderPills(orderedLabels) {
    var bar = document.getElementById(PILLS_ID);
    if (!bar) return;
    bar.innerHTML = '';

    if (orderedLabels.length === 0) {
      renderEmpty(bar);
      return;
    }

    orderedLabels.forEach(function(label) {
      var snap = WZ_BALANCES[label] || {available: null, error: 'no_data', stale_seconds: 0};
      var pill = document.createElement('span');
      var classes = ['wz-pill'];
      if (label === window.WZ_SELECTED_ACCOUNT) classes.push('selected');
      if (isStale(snap)) classes.push('stale');
      pill.className = classes.join(' ');
      pill.dataset.label = label;
      pill.textContent = pillTextFor(label, snap);

      pill.addEventListener('click', function() {
        if (label === window.WZ_SELECTED_ACCOUNT) return;
        window.WZ_SELECTED_ACCOUNT = label;
        persistSelection(label).catch(function() {
          // Network blip on persist shouldn't block UI; selection still
          // applies for this session and will retry on next click.
        });
        renderPills(orderedLabels);
        recomputeAllInsufficiencyWarnings();
      });

      bar.appendChild(pill);
    });
  }

  function refreshBalances() {
    return fetch('/api/wagerzon/balances')
      .then(function(r) { return r.json(); })
      .then(function(payload) {
        var incoming = Array.isArray(payload && payload.balances) ? payload.balances : null;

        // Guard against transient empty/malformed payloads wiping good
        // pills mid-session. If we already have rendered balances and the
        // server hands back nothing, keep the old state and warn — a real
        // "no accounts configured" state will still render on first load
        // because WZ_BALANCES is empty then.
        var hadBalances = Object.keys(WZ_BALANCES).length > 0;
        if (!incoming || (incoming.length === 0 && hadBalances)) {
          console.warn('refreshBalances: empty/invalid payload; keeping prior pills', payload);
          return;
        }

        var orderedLabels = [];
        WZ_BALANCES = {};
        incoming.forEach(function(snap) {
          WZ_BALANCES[snap.label] = snap;
          orderedLabels.push(snap.label);
        });

        // Default-to-first behaviour. The old <select> got this for free
        // via browser default-first-option; the div-based pill row needs
        // an explicit default + persist so the next page load is stable.
        if (!window.WZ_SELECTED_ACCOUNT && orderedLabels.length > 0) {
          window.WZ_SELECTED_ACCOUNT = orderedLabels[0];
          persistSelection(orderedLabels[0]).catch(function() {});
        }

        renderPills(orderedLabels);
        recomputeAllInsufficiencyWarnings();
      })
      .catch(function(err) {
        // Network failure or non-JSON response. Don't touch UI state —
        // surface to console so the click doesn't feel like a no-op.
        console.warn('refreshBalances failed:', err);
      });
  }

  function loadLastUsed() {
    return fetch('/api/wagerzon/last-used')
      .then(function(r) { return r.json(); })
      .then(function(payload) {
        window.WZ_SELECTED_ACCOUNT = payload.label || null;
      });
  }

  function recomputeAllInsufficiencyWarnings() {
    if (typeof window._wzRecomputeWarnings === 'function') {
      window._wzRecomputeWarnings();
    }
  }

  // Public hook used by placeParlay() after a successful placement
  // to update the pill without a full re-fetch.
  window.wzApplyBalanceAfter = function(label, snap) {
    if (!snap) return;
    WZ_BALANCES[label] = snap;
    // Caller doesn't have orderedLabels — derive from the current DOM
    // so render order is preserved.
    var bar = document.getElementById(PILLS_ID);
    var labels = [];
    if (bar) {
      bar.querySelectorAll('.wz-pill[data-label]').forEach(function(el) {
        labels.push(el.dataset.label);
      });
    }
    if (labels.length === 0) labels = Object.keys(WZ_BALANCES);
    renderPills(labels);
    recomputeAllInsufficiencyWarnings();
  };

  // Per-parlay insufficient-balance warning. Reads each Place button's
  // data-risk against the currently-selected account's available balance
  // and writes a short message into the sibling .wz-insufficient-warning
  // span. Suppresses output when the account has no successful fetch yet
  // (available === null) to avoid showing misleading text.
  window._wzRecomputeWarnings = function() {
    var label = window.WZ_SELECTED_ACCOUNT;
    var snap  = label ? WZ_BALANCES[label] : null;
    var available = snap ? snap.available : null;

    document.querySelectorAll('button[data-hash][data-risk]').forEach(function(btn) {
      var risk = parseFloat(btn.dataset.risk);
      var warn = btn.parentElement.querySelector('.wz-insufficient-warning');
      if (!warn) return;

      if (available === null || isNaN(risk)) {
        warn.textContent = '';
        return;
      }
      if (risk > available) {
        warn.textContent = '⚠ insufficient on ' + label +
          ' (' + fmtMoney(available) + ' < ' + fmtMoney(risk) + ' risk)';
      } else {
        warn.textContent = '';
      }
    });
  };

  document.addEventListener('DOMContentLoaded', function() {
    var refreshBtn = document.getElementById(REFRESH_BTN_ID);
    if (refreshBtn) {
      refreshBtn.addEventListener('click', refreshBalances);
    }

    loadLastUsed().catch(function() {
      // Network blip on the last-used GET shouldn't block the rest of the
      // initial render. Pills will still populate; default-to-first kicks
      // in inside refreshBalances if WZ_SELECTED_ACCOUNT is still null.
    }).then(refreshBalances);

    // Watch for parlay-table re-renders (reactable pagination, hot-swap
    // after combined placement, etc.) so the warning is recomputed for
    // the new rows. Same pattern as the existing same-game observer on
    // bets-table-container.
    var parlayContainer = document.getElementById('parlays-table-container');
    if (parlayContainer && typeof MutationObserver === 'function') {
      var _wzWarnDebounce = null;
      var obs = new MutationObserver(function() {
        if (_wzWarnDebounce) clearTimeout(_wzWarnDebounce);
        _wzWarnDebounce = setTimeout(function() {
          _wzWarnDebounce = null;
          recomputeAllInsufficiencyWarnings();
        }, 100);
      });
      obs.observe(parlayContainer, {childList: true, subtree: true});
    }
  });
})();
)"))
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

# Fragment-only mode: render JUST the source-parlays reactable (with combo
# residuals applied) and emit its HTML to stdout. Used by /api/parlay-table-
# fragment so the dashboard can hot-swap the parlay table after a combined-
# parlay placement without rebuilding the whole report.
if (any(commandArgs(trailingOnly = TRUE) == "--parlay-fragment")) {
  # mlb_parlay_opportunities was moved to mlb_mm.duckdb in Task 4 to avoid
  # contention with the pipeline's long write lock on mlb.duckdb.
  con_mlb <- dbConnect(duckdb(), dbdir = "Answer Keys/mlb_mm.duckdb", read_only = TRUE)
  on.exit(try(dbDisconnect(con_mlb), silent = TRUE), add = TRUE)
  parlay_opps <- tryCatch(
    dbGetQuery(con_mlb, "SELECT * FROM mlb_parlay_opportunities"),
    error = function(e) tibble()
  )
  dbDisconnect(con_mlb)

  con_dash <- dbConnect(duckdb(), dbdir = DB_PATH, read_only = TRUE)
  on.exit(try(dbDisconnect(con_dash), silent = TRUE), add = TRUE)
  placed_parlays <- tryCatch(
    dbGetQuery(con_dash, "SELECT * FROM placed_parlays"),
    error = function(e) tibble()
  )
  sz <- tryCatch(
    dbGetQuery(con_dash, "SELECT param, value FROM sizing_settings"),
    error = function(e) tibble(param = character(), value = numeric())
  )
  parlay_bankroll <- {
    v <- sz$value[sz$param == "parlay_bankroll"]
    if (length(v) > 0) v[1] else 100
  }
  parlay_kelly_mult <- {
    v <- sz$value[sz$param == "parlay_kelly_mult"]
    if (length(v) > 0) v[1] else 0.25
  }
  dbDisconnect(con_dash)

  parlays_table <- create_parlays_table(parlay_opps, placed_parlays,
                                        parlay_bankroll, parlay_kelly_mult)
  if (!is.null(parlays_table)) {
    # Emit just the widget HTML — no <html>/<head>; dependencies are already
    # loaded on the page. Caller will run HTMLWidgets.staticRender() to
    # re-initialize the widget after DOM-swap.
    cat(as.character(htmltools::as.tags(parlays_table)))
  }
  quit(status = 0)
}

cat("=== MLB Answer Key Dashboard ===\n\n")

# Load bets from duckdb. Source moved from mlb.duckdb → mlb_mm.duckdb
# on 2026-05-05 so the dashboard is no longer blocked by the pipeline's
# long write lock on mlb.duckdb. tryCatch mirrors the parlay/trifecta
# loaders below: a transient lock or missing DB degrades to empty bets
# rather than crashing the whole render.
cat("Loading bets from database...\n")
all_bets <- tryCatch({
  con <- dbConnect(duckdb(), dbdir = "Answer Keys/mlb_mm.duckdb", read_only = TRUE)
  tables <- dbListTables(con)
  result <- if (!"mlb_bets_combined" %in% tables) {
    tibble()
  } else {
    dbGetQuery(con, "SELECT * FROM mlb_bets_combined") %>%
      filter(is.na(pt_start_time) | pt_start_time > Sys.time())
  }
  dbDisconnect(con, shutdown = TRUE)
  result
}, error = function(e) {
  cat(sprintf("  Warning: could not load mlb_bets_combined (%s) — rendering with empty bets\n", e$message))
  tibble()
})

cat(sprintf("Loaded %d bets\n", nrow(all_bets)))

book_prices_long <- tryCatch({
  con <- dbConnect(duckdb(), dbdir = "Answer Keys/mlb_mm.duckdb", read_only = TRUE)
  tables <- dbListTables(con)
  result <- if (!"mlb_bets_book_prices" %in% tables) {
    NULL
  } else {
    dbGetQuery(con, "SELECT * FROM mlb_bets_book_prices")
  }
  dbDisconnect(con, shutdown = TRUE)
  result
}, error = function(e) {
  warning(sprintf(
    "[bets-tab] mlb_bets_book_prices not loaded (%s) — falling back to single-book layout",
    conditionMessage(e)))
  NULL
})

# Pivot long -> wide so each bet has a single row with per-book columns.
# Output columns per side:
#   <book>_line_quoted, <book>_american_odds, <book>_is_exact_line
# for book in (wagerzon, hoop88, bfa, bookmaker, bet105,
#              draftkings, fanduel, pinnacle).
# Two rows per bet (one with side='pick', one with side='opposite').
pivot_book_prices_wide <- function(long_frame) {
  if (is.null(long_frame) || nrow(long_frame) == 0) return(NULL)
  long_frame %>%
    select(bet_row_id, side, bookmaker,
           line_quoted, american_odds, is_exact_line) %>%
    tidyr::pivot_wider(
      names_from  = bookmaker,
      values_from = c(line_quoted, american_odds, is_exact_line),
      names_glue  = "{bookmaker}_{.value}",
      values_fill = list(line_quoted = NA_real_,
                         american_odds = NA_integer_,
                         is_exact_line = NA)
    )
}

book_prices_wide <- pivot_book_prices_wide(book_prices_long)

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
  bets_table <- create_bets_table(all_bets, placed_bets,
                                  book_prices_wide = book_prices_wide)
} else {
  bets_table <- tags$div(
    style = "text-align: center; padding: 48px; color: #8b949e;",
    tags$p(style = "font-size: 1.1rem;", "No +EV bets found for upcoming games."),
    tags$p(style = "font-size: 0.85rem;", "Try clicking Refresh when more games are available.")
  )
}
placed_table <- create_placed_bets_table(placed_bets)

# Load parlay opportunities (from mlb_mm.duckdb — moved in Task 4 to avoid
# contention with the pipeline's long write lock on mlb.duckdb).
cat("Loading parlay opportunities...\n")
parlay_opps <- tryCatch({
  pcon <- dbConnect(duckdb(), dbdir = "Answer Keys/mlb_mm.duckdb", read_only = TRUE)
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
    # Fetch all statuses so the dashboard can render auto-placement state:
    # "pending" = manually placed, "placed" = auto-placed (show ticket),
    # "price_moved"/"rejected"/"auth_error"/"network_error"/"orphaned" = show error.
    dbGetQuery(ppcon, "SELECT * FROM placed_parlays WHERE game_time IS NULL OR game_time > NOW()"),
    error = function(e) tibble(parlay_hash = character())
  )
  dbDisconnect(ppcon, shutdown = TRUE)
  result
}, error = function(e) tibble(parlay_hash = character()))

# Load parlay sizing settings (min edge, bankroll, kelly multiplier).
# Bankroll/kelly_mult feed apply_combo_residuals() so the conditional Kelly
# solver knows the operator's current sizing context. Defaults match the JS
# fallbacks in the dashboard frontend.
parlay_min_edge   <- 0
parlay_bankroll   <- 100
parlay_kelly_mult <- 0.25
tryCatch({
  mecon <- dbConnect(duckdb(), dbdir = DB_PATH, read_only = TRUE)
  ss <- dbGetQuery(mecon, "SELECT param, value FROM sizing_settings WHERE param IN ('parlay_min_edge','parlay_bankroll','parlay_kelly_mult')")
  dbDisconnect(mecon, shutdown = TRUE)
  if (nrow(ss) > 0) {
    me_row <- ss[ss$param == "parlay_min_edge", ]
    if (nrow(me_row) > 0) parlay_min_edge   <- as.numeric(me_row$value[1])
    bk_row <- ss[ss$param == "parlay_bankroll", ]
    if (nrow(bk_row) > 0) parlay_bankroll   <- as.numeric(bk_row$value[1])
    km_row <- ss[ss$param == "parlay_kelly_mult", ]
    if (nrow(km_row) > 0) parlay_kelly_mult <- as.numeric(km_row$value[1])
  }
}, error = function(e) NULL)

if (nrow(parlay_opps) > 0 && parlay_min_edge > 0) {
  parlay_opps <- parlay_opps %>% filter(edge_pct >= parlay_min_edge)
}

cat(sprintf("Found %d parlay opportunities (min edge %.0f%%), %d placed parlays\n",
            nrow(parlay_opps), parlay_min_edge, nrow(placed_parlays)))

# Create parlay tables
if (nrow(parlay_opps) > 0) {
  parlays_table <- create_parlays_table(parlay_opps, placed_parlays,
                                        parlay_bankroll = parlay_bankroll,
                                        parlay_kelly_mult = parlay_kelly_mult)
} else {
  parlays_table <- NULL
}
placed_parlays_table <- create_placed_parlays_table(placed_parlays)

# Load trifecta opportunities + placed trifectas
cat("Loading trifecta opportunities...\n")
trifecta_opps    <- load_trifecta_opps(file.path(NFLWORK_ROOT, "Answer Keys", "mlb_mm.duckdb"))
placed_trifectas <- load_placed_trifectas(DB_PATH)
trifectas_table  <- create_trifectas_table(trifecta_opps, placed_trifectas)
cat(sprintf("Found %d trifecta opportunities, %d placed trifectas\n",
            nrow(trifecta_opps), nrow(placed_trifectas)))

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
                       parlay_filter_options_json,
                       trifectas_table, trifecta_opps)

# Save
save_html(page, OUTPUT_PATH, libdir = "lib")
cat("\nDashboard saved to:", OUTPUT_PATH, "\n")
cat("Open http://localhost:8083 to view\n")
