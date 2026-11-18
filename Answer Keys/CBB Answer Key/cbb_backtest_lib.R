# =============================================================================
# CBB walk-forward backtest: data loading, pricing, and the report metrics.
# =============================================================================
# Sourced by CBB_Backtest.R (the entry script; usage there) and by
# tests/test_cbb_backtest_lib.R. Every function here is pure except:
#   - load_cbb_history() / load_h1_closing_lines(): read cbb.duckdb (read-only)
#   - save_backtest_run(): writes cbb_backtest.duckdb (replaces one label's rows)
#
# What a run measures, per test game priced only off games dated before it:
#   1. Sample centring: does the sample cover the PARENT full-game spread/total
#      at the no-vig target rate (gap, home/over lean), and how far is the
#      sample's median margin/total from the parent line.
#   2. Accuracy: log loss of the sample's 1H spread / total / moneyline
#      probabilities against the book closing fair (median probit devig over
#      books), plus full-game main and alt lines graded on the final score.
#   3. Betting: flat $100 bets at the best closing price wherever the model
#      claims an edge >= the threshold, ROI with a by-game bootstrap interval.
#
# The run needs Tools.R already sourced (run_answer_key_sample,
# compute_dispersion, walk_forward_training_set, has_walk_forward_history).
# =============================================================================

suppressMessages({
  library(data.table)
  library(DBI)
})

# CBB.R sizes N as round(nrow(DT) * 0.02); the backtest sizes it off each test
# game's prior history the same way.
SAMPLE_PCT_DEFAULT <- 0.02
MIN_BOOKS_PER_LINE <- 2          # a 1H line needs 2+ books for a closing fair
MIN_NONPUSH_SAMPLE_GAMES <- 10   # fewer non-push sample games -> no probability
EDGE_THRESHOLD_DEFAULT <- 0.03
FLAT_STAKE <- 100
BOOTSTRAP_DRAWS <- 2000
LOGLOSS_PROB_CLIP <- 1e-4
H1_MARKETS <- c("spreads_h1", "totals_h1", "h2h_h1")
# Full-game lines graded on the final score around the parent line
FG_SPREAD_OFFSETS <- c(-6, -3, 0, 3, 6)
FG_TOTAL_OFFSETS  <- c(-10, -5, 0, 5, 10)

PREDICTIONS_TABLE <- "backtest_predictions"
RUNS_TABLE <- "backtest_runs"


# -----------------------------------------------------------------------------
# Data loading
# -----------------------------------------------------------------------------

#' Historical CBB games, built the way CBB.R builds its DT.
#'
#' The parent lines and no-vig targets come from cbb_betting_pbp's own
#' home_spread / total_line / consensus_devig_* columns (written by
#' build_betting_pbp.R). The old backtests re-joined cbb_closing_odds on top,
#' which duplicated home_spread / total_line into .x/.y and crashed.
#'
#' @param con DBI connection to cbb.duckdb
#' @return data.table, one row per game with game_date, odds_id, parent lines,
#'   targets, full-game and 1H scores, actual_cover / actual_over (push = 0.5)
load_cbb_history <- function(con) {
  games <- as.data.table(dbGetQuery(con, "
    SELECT odds_id, game_id, game_date, home_team, away_team,
           home_spread, total_line,
           consensus_devig_home_odds, consensus_devig_over_odds,
           game_home_margin_fg, game_total_fg,
           game_home_margin_h1, game_total_h1
    FROM cbb_betting_pbp
    WHERE home_spread IS NOT NULL AND total_line IS NOT NULL
      AND consensus_devig_home_odds IS NOT NULL AND consensus_devig_over_odds IS NOT NULL
      AND game_home_margin_fg IS NOT NULL AND game_total_fg IS NOT NULL
      AND game_home_margin_h1 IS NOT NULL AND game_total_h1 IS NOT NULL
  "))
  if (nrow(games) == 0) stop("load_cbb_history: cbb_betting_pbp returned no complete games")
  games[, game_date := as.Date(game_date)]
  # Same column names and push grading as CBB.R
  games[, `:=`(home_margin = game_home_margin_fg, total_final_score = game_total_fg,
               home_spread_odds = consensus_devig_home_odds, over_odds = consensus_devig_over_odds)]
  games[, actual_cover := ifelse(home_margin > -home_spread, 1, ifelse(home_margin == -home_spread, 0.5, 0))]
  games[, actual_over  := ifelse(total_final_score > total_line, 1, ifelse(total_final_score == total_line, 0.5, 0))]
  setorder(games, game_date, odds_id)
  games[]
}

#' Two-way probit devig, vectorised; matches Tools.R::devig_american's p1.
#' @return fair probability of side 1 (NA where either price is NA or 0)
devig_probit_side1 <- function(side1_american_odds, side2_american_odds, eps = 1e-9) {
  implied <- function(odds) ifelse(odds > 0, 100 / (odds + 100), -odds / (-odds + 100))
  bad <- is.na(side1_american_odds) | is.na(side2_american_odds) |
    side1_american_odds == 0 | side2_american_odds == 0
  z1 <- qnorm(pmin(pmax(implied(side1_american_odds), eps), 1 - eps))
  z2 <- qnorm(pmin(pmax(implied(side2_american_odds), eps), 1 - eps))
  fair <- pnorm(z1 - (z1 + z2) / 2)
  fair[bad] <- NA_real_
  fair
}

#' 1H closing lines with a book fair and the best price on each side.
#'
#' Side 1 is home (spreads_h1 / h2h_h1) or over (totals_h1); `line` is the
#' home spread, the total, or 0 for the moneyline.
#'
#' @param con DBI connection to cbb.duckdb (reads cbb_derivative_closing_odds)
#' @param books NULL for every book, else bookmaker_keys to keep (e.g.
#'   "betonline"); with one book a line needs only that book
#' @return data.table: event_id, market, line, book_fair_prob (median over
#'   books), n_books, best_side1_american_odds, best_side2_american_odds
load_h1_closing_lines <- function(con, books = NULL, min_books = MIN_BOOKS_PER_LINE) {
  # One price per book, market, outcome and point (a handful of duplicate rows)
  book_pairs <- as.data.table(dbGetQuery(con, "
    WITH quotes AS (
      SELECT event_id, bookmaker_key, market, outcome_name, outcome_price, outcome_point,
             home_team, away_team
      FROM cbb_derivative_closing_odds
      WHERE market IN ('spreads_h1', 'totals_h1', 'h2h_h1') AND outcome_price IS NOT NULL
      QUALIFY row_number() OVER (
        PARTITION BY event_id, bookmaker_key, market, outcome_name, outcome_point
        ORDER BY outcome_price) = 1
    )
    SELECT s1.event_id, s1.market, s1.bookmaker_key,
           CASE WHEN s1.market = 'h2h_h1' THEN 0 ELSE s1.outcome_point END AS line,
           s1.outcome_price AS side1_american_odds,
           s2.outcome_price AS side2_american_odds
    FROM quotes s1
    JOIN quotes s2
      ON s1.event_id = s2.event_id AND s1.bookmaker_key = s2.bookmaker_key
     AND s1.market = s2.market
    WHERE (s1.market = 'spreads_h1' AND s1.outcome_name = s1.home_team
           AND s2.outcome_name = s2.away_team AND s2.outcome_point = -s1.outcome_point)
       OR (s1.market = 'h2h_h1' AND s1.outcome_name = s1.home_team
           AND s2.outcome_name = s2.away_team)
       OR (s1.market = 'totals_h1' AND s1.outcome_name = 'Over'
           AND s2.outcome_name = 'Under' AND s2.outcome_point = s1.outcome_point)
  "))
  if (!is.null(books)) {
    unknown <- setdiff(books, book_pairs$bookmaker_key)
    if (length(unknown) > 0) {
      stop("load_h1_closing_lines: no 1H prices for book(s) ", paste(unknown, collapse = ", "),
           "; books present: ", paste(sort(unique(book_pairs$bookmaker_key)), collapse = ", "))
    }
    book_pairs <- book_pairs[bookmaker_key %in% books]
    min_books <- min(min_books, length(books))
  }
  book_pairs[, side1_fair_prob := devig_probit_side1(side1_american_odds, side2_american_odds)]
  lines <- book_pairs[is.finite(side1_fair_prob), .(
    book_fair_prob = median(side1_fair_prob),
    n_books = .N,
    best_side1_american_odds = max(side1_american_odds),  # higher American = better payout
    best_side2_american_odds = max(side2_american_odds)
  ), by = .(event_id, market, line)]
  lines <- lines[n_books >= min_books]
  setkey(lines, event_id)
  lines[]
}


# -----------------------------------------------------------------------------
# Pricing
# -----------------------------------------------------------------------------

#' Share of non-push sample games above `threshold`; NA when too few.
sample_prob_over <- function(values, threshold, min_nonpush = MIN_NONPUSH_SAMPLE_GAMES) {
  non_push <- values[!is.na(values) & values != threshold]
  if (length(non_push) < min_nonpush) return(NA_real_)
  mean(non_push > threshold)
}

#' 1 above the threshold, 0 below, NA on a push.
outcome_vs_threshold <- function(value, threshold) {
  ifelse(is.na(value) | value == threshold, NA_real_, as.numeric(value > threshold))
}

#' Test games: enough walk-forward history behind them and 1H closing lines.
select_test_games <- function(history, h1_lines, max_games = 0, seed = 42) {
  eligible <- has_walk_forward_history(history, history$game_date) &
    history$odds_id %in% h1_lines$event_id
  test_idx <- which(eligible)
  if (max_games > 0 && max_games < length(test_idx)) {
    set.seed(seed)
    test_idx <- sort(sample(test_idx, max_games))
  }
  test_idx
}

#' Price one test game off the games dated before it.
#'
#' @param test_game one row of load_cbb_history()
#' @param history all of load_cbb_history()
#' @param game_h1_lines this game's rows of load_h1_closing_lines()
#' @param dispersion compute_dispersion() output (list with ss, st)
#' @return data.table of predictions (one row per line priced) with the game's
#'   sample diagnostics repeated on each row, or NULL when the game cannot be
#'   priced (not enough history, or the sampler failed)
price_test_game <- function(test_game, history, game_h1_lines, dispersion,
                            sample_pct = SAMPLE_PCT_DEFAULT) {
  train <- walk_forward_training_set(history, test_game$game_date)
  if (is.null(train)) return(NULL)
  target_N <- round(nrow(train) * sample_pct)
  result <- tryCatch(
    suppressWarnings(run_answer_key_sample(
      test_game$odds_id, test_game$home_spread, test_game$total_line,
      test_game$consensus_devig_home_odds, test_game$consensus_devig_over_odds,
      train, dispersion$ss, dispersion$st, target_N, use_spread_line = TRUE)),
    error = function(e) NULL)
  if (is.null(result) || is.null(result$sample) || nrow(result$sample) == 0) return(NULL)
  sample <- result$sample

  rows <- list()
  add_row <- function(market, line, model_prob, book_fair_prob, best1, best2, n_books, outcome) {
    rows[[length(rows) + 1]] <<- data.table(
      market = market, line = line, model_prob = model_prob, book_fair_prob = book_fair_prob,
      best_side1_american_odds = best1, best_side2_american_odds = best2,
      n_books = n_books, outcome = outcome)
  }

  for (j in seq_len(nrow(game_h1_lines))) {
    quote <- game_h1_lines[j]
    if (quote$market == "spreads_h1") {
      values <- sample$game_home_margin_h1; actual <- test_game$game_home_margin_h1; threshold <- -quote$line
    } else if (quote$market == "totals_h1") {
      values <- sample$game_total_h1; actual <- test_game$game_total_h1; threshold <- quote$line
    } else {  # h2h_h1: home leads at the half; a tied half is a push
      values <- sample$game_home_margin_h1; actual <- test_game$game_home_margin_h1; threshold <- 0
    }
    add_row(quote$market, quote$line, sample_prob_over(values, threshold), quote$book_fair_prob,
            quote$best_side1_american_odds, quote$best_side2_american_odds, quote$n_books,
            outcome_vs_threshold(actual, threshold))
  }

  # Full-game lines: book fair only at the parent line (the no-vig target)
  for (offset in FG_SPREAD_OFFSETS) {
    line <- test_game$home_spread + offset
    add_row(sprintf("fg_spread_%+d", offset), line, sample_prob_over(sample$home_margin, -line),
            if (offset == 0) test_game$consensus_devig_home_odds else NA_real_, NA, NA, NA,
            outcome_vs_threshold(test_game$home_margin, -line))
  }
  for (offset in FG_TOTAL_OFFSETS) {
    line <- test_game$total_line + offset
    add_row(sprintf("fg_total_%+d", offset), line, sample_prob_over(sample$total_final_score, line),
            if (offset == 0) test_game$consensus_devig_over_odds else NA_real_, NA, NA, NA,
            outcome_vs_threshold(test_game$total_final_score, line))
  }

  predictions <- rbindlist(rows)
  predictions[, `:=`(
    odds_id = test_game$odds_id, game_date = test_game$game_date,
    parent_spread = test_game$home_spread, parent_total = test_game$total_line,
    target_N = target_N, final_N = result$final_N, converged = isTRUE(result$converged),
    sample_median_margin_miss = median(sample$home_margin) + test_game$home_spread,
    sample_median_total_miss = median(sample$total_final_score) - test_game$total_line,
    sample_mean_spread_drift = mean(sample$home_spread) - test_game$home_spread)]
  setcolorder(predictions, c("odds_id", "game_date", "market", "line"))
  predictions[]
}

#' Price every test game, in parallel (forked; Tools.R must be sourced).
#' @return list(predictions, unpriced_games = games the sampler could not price)
price_test_games <- function(history, h1_lines, test_idx, sample_pct = SAMPLE_PCT_DEFAULT,
                             cores = 4) {
  dispersion <- compute_dispersion(history, moneyline = FALSE,
                                   spread_col = "home_spread", total_col = "total_line")
  price_one <- function(i) {
    game <- history[i]
    price_test_game(game, history, h1_lines[.(game$odds_id), nomatch = NULL], dispersion, sample_pct)
  }
  per_game <- parallel::mclapply(test_idx, price_one, mc.cores = cores, mc.preschedule = TRUE)
  failed <- vapply(per_game, function(x) inherits(x, "try-error"), logical(1))
  if (any(failed)) {
    stop("price_test_games: ", sum(failed), " workers errored; first: ",
         as.character(per_game[[which(failed)[1]]]))
  }
  unpriced <- sum(vapply(per_game, is.null, logical(1)))
  list(predictions = rbindlist(per_game), unpriced_games = unpriced)
}


# -----------------------------------------------------------------------------
# Storage (cbb_backtest.duckdb)
# -----------------------------------------------------------------------------

#' Save one run's predictions under `label`, replacing that label's old rows.
#'
#' Side effect: writes PREDICTIONS_TABLE and RUNS_TABLE in the DB at `con`
#' (DELETE that label's rows, then append), inside one transaction.
save_backtest_run <- function(con, label, predictions, run_info) {
  predictions <- copy(predictions)
  predictions[, label := label]
  setcolorder(predictions, "label")
  run_row <- as.data.table(c(list(label = label), run_info))
  dbWithTransaction(con, {
    for (table in c(PREDICTIONS_TABLE, RUNS_TABLE)) {
      if (dbExistsTable(con, table)) {
        dbExecute(con, sprintf("DELETE FROM %s WHERE label = ?", table), params = list(label))
      }
    }
    dbWriteTable(con, PREDICTIONS_TABLE, predictions, append = TRUE)
    dbWriteTable(con, RUNS_TABLE, run_row, append = TRUE)
  })
  invisible(nrow(predictions))
}

#' Read saved predictions for the given labels; stops naming any missing label.
load_backtest_runs <- function(con, labels) {
  if (!dbExistsTable(con, PREDICTIONS_TABLE)) stop("no saved backtest runs yet; run `run` first")
  saved <- dbGetQuery(con, sprintf("SELECT DISTINCT label FROM %s", PREDICTIONS_TABLE))$label
  missing <- setdiff(labels, saved)
  if (length(missing) > 0) {
    stop("no saved run for label(s): ", paste(missing, collapse = ", "),
         "; saved labels: ", paste(saved, collapse = ", "))
  }
  placeholders <- paste(rep("?", length(labels)), collapse = ", ")
  predictions <- as.data.table(dbGetQuery(
    con, sprintf("SELECT * FROM %s WHERE label IN (%s)", PREDICTIONS_TABLE, placeholders),
    params = as.list(labels)))
  predictions[, game_date := as.Date(game_date)]
  runs <- as.data.table(dbGetQuery(
    con, sprintf("SELECT * FROM %s WHERE label IN (%s)", RUNS_TABLE, placeholders),
    params = as.list(labels)))
  list(predictions = predictions, runs = runs)
}


# -----------------------------------------------------------------------------
# Report metrics (pure; input is predictions with a `label` column)
# -----------------------------------------------------------------------------

#' Keep lines every label priced with a graded (non-push) outcome, so labels
#' are compared on the same games and lines.
common_graded_lines <- function(predictions) {
  graded <- predictions[!is.na(model_prob) & !is.na(outcome)]
  graded[, line_key := paste(odds_id, market, line)]
  n_labels <- uniqueN(graded$label)
  keep <- graded[, .(n = uniqueN(label)), by = line_key][n == n_labels, line_key]
  graded <- graded[line_key %in% keep]
  graded[, market_group := fifelse(market %in% H1_MARKETS, market,
                           fifelse(market == "fg_spread_+0", "fg_spread_main",
                           fifelse(market == "fg_total_+0", "fg_total_main",
                           fifelse(startsWith(market, "fg_spread"), "fg_spread_alts", "fg_total_alts"))))]
  graded[]
}

log_loss <- function(prob, outcome) {
  prob <- pmin(pmax(prob, LOGLOSS_PROB_CLIP), 1 - LOGLOSS_PROB_CLIP)
  -(outcome * log(prob) + (1 - outcome) * log(1 - prob))
}

american_to_decimal <- function(american_odds) {
  ifelse(american_odds > 0, 1 + american_odds / 100, 1 + 100 / abs(american_odds))
}

#' By-game bootstrap of sum(numerator) / sum(denominator); rows of one game
#' are resampled together because they share one sample.
bootstrap_ratio_by_game <- function(game, numerator, denominator, draws = BOOTSTRAP_DRAWS, seed = 1) {
  per_game <- data.table(game, numerator, denominator)[, .(num = sum(numerator), den = sum(denominator)), by = game]
  set.seed(seed)
  resampled <- replicate(draws, {
    pick <- sample.int(nrow(per_game), replace = TRUE)
    sum(per_game$num[pick]) / sum(per_game$den[pick])
  })
  c(lo = unname(quantile(resampled, 0.025)), hi = unname(quantile(resampled, 0.975)))
}

#' Sample centring per label: does the sample land on the parent line?
#' cover_gap = sample's non-push cover rate of the parent spread minus the
#' no-vig target; mean_cover_gap > 0 is a home lean, mean_over_gap > 0 an over lean.
summarise_sample_centring <- function(graded) {
  per_game <- unique(graded[, .(label, odds_id, target_N, final_N, converged,
                                sample_median_margin_miss, sample_median_total_miss,
                                sample_mean_spread_drift)])
  diagnostics <- per_game[, .(
    games = .N,
    mean_abs_median_margin_miss = mean(abs(sample_median_margin_miss)),
    mean_abs_median_total_miss = mean(abs(sample_median_total_miss)),
    mean_abs_spread_drift = mean(abs(sample_mean_spread_drift)),
    pct_converged = mean(converged),
    mean_N_ratio = mean(final_N / target_N)), by = label]
  main <- graded[market_group %in% c("fg_spread_main", "fg_total_main")]
  gaps <- main[, .(lines = .N,
                   mean_abs_gap = mean(abs(model_prob - book_fair_prob)),
                   mean_gap = mean(model_prob - book_fair_prob),
                   whole_line_mean_gap = mean((model_prob - book_fair_prob)[line == round(line)])),
               by = .(label, market_group)]
  list(diagnostics = diagnostics, main_line_gaps = gaps)
}

#' Log loss of the model and the book close by label and market group.
summarise_accuracy <- function(graded) {
  graded[, .(
    lines = .N, games = uniqueN(odds_id),
    model_logloss = mean(log_loss(model_prob, outcome)),
    book_logloss = if (all(is.na(book_fair_prob))) NA_real_ else mean(log_loss(book_fair_prob, outcome), na.rm = TRUE),
    mean_abs_vs_book = mean(abs(model_prob - book_fair_prob)),
    mean_model_minus_book = mean(model_prob - book_fair_prob)
  ), by = .(label, market_group)][order(market_group, label)]
}

#' Flat bets at the best closing price where the model's edge >= threshold.
#' One 1H line can give a bet on either side, never both.
flat_bets <- function(graded, edge_threshold = EDGE_THRESHOLD_DEFAULT) {
  h1 <- graded[market_group %in% H1_MARKETS]
  h1[, `:=`(side1_ev = model_prob * american_to_decimal(best_side1_american_odds) - 1,
            side2_ev = (1 - model_prob) * american_to_decimal(best_side2_american_odds) - 1)]
  side1 <- h1[side1_ev >= edge_threshold, .(
    label, market_group, odds_id, claimed_edge = side1_ev,
    pnl = FLAT_STAKE * ifelse(outcome == 1, american_to_decimal(best_side1_american_odds) - 1, -1))]
  side2 <- h1[side2_ev >= edge_threshold, .(
    label, market_group, odds_id, claimed_edge = side2_ev,
    pnl = FLAT_STAKE * ifelse(outcome == 0, american_to_decimal(best_side2_american_odds) - 1, -1))]
  rbind(side1, side2)
}

summarise_flat_bets <- function(bets, by = c("label", "market_group")) {
  if (nrow(bets) == 0) return(data.table())
  bets[, {
    interval <- bootstrap_ratio_by_game(odds_id, pnl, rep(FLAT_STAKE, .N))
    .(bets = .N, claimed_edge = mean(claimed_edge), roi = sum(pnl) / (FLAT_STAKE * .N),
      roi_lo = interval[["lo"]], roi_hi = interval[["hi"]])
  }, by = by]
}

#' Paired change in log loss from label `from` to label `to` by market group,
#' with a by-game bootstrap 95% interval (negative = `to` is better).
compare_logloss <- function(graded, from, to) {
  pair <- graded[label %in% c(from, to)]
  pair[, model_logloss := log_loss(model_prob, outcome)]
  wide <- dcast(pair, odds_id + line_key + market_group ~ label, value.var = "model_logloss")
  wide[, change := get(to) - get(from)]
  wide[, {
    interval <- bootstrap_ratio_by_game(odds_id, change, rep(1, .N))
    .(change = paste(from, "->", to), lines = .N, logloss_change = mean(change),
      lo = interval[["lo"]], hi = interval[["hi"]])
  }, by = market_group][order(market_group)]
}

#' Print the whole report for the labels in `predictions` (first label = baseline).
print_backtest_report <- function(predictions, runs, edge_threshold = EDGE_THRESHOLD_DEFAULT) {
  labels <- unique(runs$label)
  graded <- common_graded_lines(predictions)
  if (nrow(graded) == 0) stop("print_backtest_report: no graded lines common to all labels")
  order_labels <- function(dt) if ("label" %in% names(dt)) dt[order(match(label, labels))] else dt

  cat("\n== Runs ==\n");                  print(runs)
  cat(sprintf("\nCompared on %d games, %s to %s (lines every label priced, pushes dropped)\n",
              uniqueN(graded$odds_id), format(min(graded$game_date)), format(max(graded$game_date))))
  centring <- summarise_sample_centring(graded)
  cat("\n== Sample centring (per game) ==\n"); print(order_labels(centring$diagnostics), digits = 3)
  cat("\n== Main line: sample cover/over rate minus no-vig target (mean_gap > 0 = home / over lean) ==\n")
  print(centring$main_line_gaps[order(market_group, match(label, labels))], digits = 3)
  cat("\n== Log loss by market (lower is better) ==\n"); print(summarise_accuracy(graded)[order(market_group, match(label, labels))], digits = 4)
  if (length(labels) > 1) {
    cat("\n== Paired log loss change vs", labels[1], "(by-game bootstrap 95%; negative = better) ==\n")
    print(rbindlist(lapply(labels[-1], function(l) compare_logloss(graded, labels[1], l))), digits = 3)
  }
  bets <- flat_bets(graded, edge_threshold)
  cat(sprintf("\n== Flat $%d bets at best closing price, claimed edge >= %.0f%% (1H) ==\n",
              FLAT_STAKE, 100 * edge_threshold))
  print(order_labels(summarise_flat_bets(bets))[order(market_group, match(label, labels))], digits = 3)
  print(order_labels(summarise_flat_bets(bets, by = "label")), digits = 3)
  invisible(graded)
}
