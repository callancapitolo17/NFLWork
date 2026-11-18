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
  library(tidyverse)  # Tools.R's predict_*_from_sample use purrr / tibble, as CBB.R loads them
})

# CBB.R sizes N as round(nrow(DT) * 0.02); the backtest sizes it off each test
# game's prior history the same way.
SAMPLE_PCT_DEFAULT <- 0.02
MIN_BOOKS_PER_LINE <- 2          # a 1H line needs 2+ books for a closing fair
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
    SELECT odds_id, game_date, home_spread, total_line,
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
  games[, `:=`(home_margin = game_home_margin_fg, total_final_score = game_total_fg)]
  games[, actual_cover := ifelse(home_margin > -home_spread, 1, ifelse(home_margin == -home_spread, 0.5, 0))]
  games[, actual_over  := ifelse(total_final_score > total_line, 1, ifelse(total_final_score == total_line, 0.5, 0))]
  setorder(games, game_date, odds_id)
  games[]
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
  # One price per book, market, outcome and point. About 250 of ~800k rows are
  # duplicates; the table has no fetch time, so the lowest price is kept,
  # which can pair the two sides from different snapshots of one book.
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
  book_pairs[, side1_fair_prob := devig_american(side1_american_odds, side2_american_odds)$p1]
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

#' Sample probability that the home side covers each home spread in `lines`,
#' through Tools.R's predict_spreads_from_sample (as live prices), so a change
#' to it in --code-dir shows up. Pushes are left out; a line no non-push sample
#' game decides is NA. `period` "h1" / "fg" picks game_home_margin_<period>.
sample_home_cover_probs <- function(sample, lines, period) {
  # The predictor names a column per line, so it is asked once per distinct line
  distinct_lines <- unique(lines)
  probs <- unlist(predict_spreads_from_sample(sample, distinct_lines, margin_col = "game_home_margin",
                                              period = period), use.names = FALSE)
  probs[!is.finite(probs)] <- NA_real_
  probs[match(lines, distinct_lines)]
}

#' Sample probability of over each total in `lines` (predict_totals_from_sample).
sample_over_probs <- function(sample, lines, period) {
  distinct_lines <- unique(lines)
  probs <- unlist(predict_totals_from_sample(sample, distinct_lines, total_col = "game_total",
                                             period = period), use.names = FALSE)
  probs[!is.finite(probs)] <- NA_real_
  probs[match(lines, distinct_lines)]
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

# One row per priced line; save_backtest_run checks a saved table against this
# before pricing starts.
PREDICTION_COLUMNS <- c(
  "odds_id", "game_date", "market", "line", "model_prob", "book_fair_prob",
  "best_side1_american_odds", "best_side2_american_odds", "outcome",
  "target_N", "final_N", "converged",
  "sample_median_margin_miss", "sample_median_total_miss", "sample_mean_spread_drift")

#' Prediction rows for one game's sample: its 1H closing lines, and full-game
#' lines around the parent (book fair only at the parent line itself, the
#' no-vig target).
prediction_rows <- function(test_game, sample, game_h1_lines) {
  h1 <- copy(game_h1_lines)
  h1[, model_prob := NA_real_]
  h1[, outcome := NA_real_]
  # h2h_h1 is the home side at spread 0: a tied half is a push
  is_spread <- h1$market %in% c("spreads_h1", "h2h_h1")
  h1[is_spread, model_prob := sample_home_cover_probs(sample, line, "h1")]
  h1[is_spread, outcome := outcome_vs_threshold(test_game$game_home_margin_h1, -line)]
  h1[market == "totals_h1", model_prob := sample_over_probs(sample, line, "h1")]
  h1[market == "totals_h1", outcome := outcome_vs_threshold(test_game$game_total_h1, line)]

  spread_lines <- test_game$home_spread + FG_SPREAD_OFFSETS
  total_lines <- test_game$total_line + FG_TOTAL_OFFSETS
  fg <- rbind(
    data.table(market = sprintf("fg_spread_%+d", FG_SPREAD_OFFSETS), line = spread_lines,
               model_prob = sample_home_cover_probs(sample, spread_lines, "fg"),
               book_fair_prob = ifelse(FG_SPREAD_OFFSETS == 0, test_game$consensus_devig_home_odds, NA_real_),
               outcome = outcome_vs_threshold(test_game$home_margin, -spread_lines)),
    data.table(market = sprintf("fg_total_%+d", FG_TOTAL_OFFSETS), line = total_lines,
               model_prob = sample_over_probs(sample, total_lines, "fg"),
               book_fair_prob = ifelse(FG_TOTAL_OFFSETS == 0, test_game$consensus_devig_over_odds, NA_real_),
               outcome = outcome_vs_threshold(test_game$total_final_score, total_lines)))
  rbind(h1[, .(market, line, model_prob, book_fair_prob, best_side1_american_odds,
               best_side2_american_odds, outcome)], fg, fill = TRUE)
}

#' Price one test game off the games dated before it.
#'
#' The distance weights (compute_dispersion) come from the training games
#' too, so nothing after the test date reaches the sample.
#'
#' @param test_game one row of load_cbb_history()
#' @param history all of load_cbb_history()
#' @param game_h1_lines this game's rows of load_h1_closing_lines()
#' @return list(status, predictions, message). status is "priced",
#'   "no_history" (under half the history behind it), "not_converged" (the
#'   sampler returned nothing, as live skips such a game) or "error" (the
#'   sampler threw; message holds why)
price_test_game <- function(test_game, history, game_h1_lines, sample_pct = SAMPLE_PCT_DEFAULT) {
  train <- walk_forward_training_set(history, test_game$game_date)
  if (is.null(train)) return(list(status = "no_history"))
  dispersion <- compute_dispersion(train, moneyline = FALSE,
                                   spread_col = "home_spread", total_col = "total_line")
  target_N <- round(nrow(train) * sample_pct)
  result <- tryCatch(
    suppressWarnings(run_answer_key_sample(
      test_game$odds_id, test_game$home_spread, test_game$total_line,
      test_game$consensus_devig_home_odds, test_game$consensus_devig_over_odds,
      train, dispersion$ss, dispersion$st, target_N, use_spread_line = TRUE)),
    error = function(e) e)
  if (inherits(result, "error")) return(list(status = "error", message = conditionMessage(result)))
  if (is.null(result) || is.null(result$sample) || nrow(result$sample) == 0) {
    return(list(status = "not_converged"))
  }
  sample <- as.data.table(result$sample)

  predictions <- prediction_rows(test_game, sample, game_h1_lines)
  predictions[, `:=`(
    odds_id = test_game$odds_id, game_date = test_game$game_date,
    target_N = target_N, final_N = result$final_N, converged = isTRUE(result$converged),
    sample_median_margin_miss = median(sample$home_margin) + test_game$home_spread,
    sample_median_total_miss = median(sample$total_final_score) - test_game$total_line,
    sample_mean_spread_drift = mean(sample$home_spread) - test_game$home_spread)]
  setcolorder(predictions, PREDICTION_COLUMNS)
  list(status = "priced", predictions = predictions)
}

#' Price every test game, in parallel (forked; Tools.R must be sourced).
#'
#' Stops when any game's sampler threw or a worker died, naming the first
#' error, so a run never saves a survivor subset as if it were complete.
#' @return list(predictions, not_converged_games)
price_test_games <- function(history, h1_lines, test_idx, sample_pct = SAMPLE_PCT_DEFAULT,
                             cores = 4) {
  price_one <- function(i) {
    game <- history[i]
    price_test_game(game, history, h1_lines[.(game$odds_id), nomatch = NULL], sample_pct)
  }
  per_game <- parallel::mclapply(test_idx, price_one, mc.cores = cores, mc.preschedule = TRUE)
  status <- vapply(per_game, function(x) if (is.list(x) && !is.null(x$status)) x$status else "worker_died",
                   character(1))
  if (any(status == "worker_died")) {
    stop("price_test_games: ", sum(status == "worker_died"), " games came back from a worker that ",
         "errored or died (out of memory?); first: ", as.character(per_game[[which(status == "worker_died")[1]]]))
  }
  if (any(status == "error")) {
    first <- per_game[[which(status == "error")[1]]]
    stop("price_test_games: the sampler threw on ", sum(status == "error"), " games; first (",
         history$odds_id[test_idx[which(status == "error")[1]]], "): ", first$message)
  }
  if (any(status == "no_history")) {
    stop("price_test_games: ", sum(status == "no_history"), " test games have under half the history ",
         "behind them; select them with select_test_games()")
  }
  priced <- per_game[status == "priced"]
  list(predictions = rbindlist(lapply(priced, `[[`, "predictions")),
       not_converged_games = sum(status == "not_converged"))
}


# -----------------------------------------------------------------------------
# Storage (cbb_backtest.duckdb)
# -----------------------------------------------------------------------------

# One row per run; run_info passed to save_backtest_run has these names, in order.
RUN_COLUMNS <- c("label", "run_at", "code_dir", "code_commit", "sample_pct", "max_games",
                 "books", "games_priced", "games_not_converged", "minutes")

# Appends are positional, so a table saved by an older version of this script
# must be rebuilt rather than appended to. Call before pricing, so a mismatch
# costs nothing.
stop_if_columns_changed <- function(con, table, columns) {
  if (!dbExistsTable(con, table)) return(invisible(TRUE))
  saved_columns <- dbListFields(con, table)
  if (!identical(saved_columns, columns)) {
    stop(table, " was saved with columns (", paste(saved_columns, collapse = ", "),
         ") but this run has (", paste(columns, collapse = ", "),
         "); delete cbb_backtest.duckdb and rerun the labels you need")
  }
  invisible(TRUE)
}

#' Save one run's predictions under `label`, replacing that label's old rows.
#'
#' Side effect: writes PREDICTIONS_TABLE and RUNS_TABLE in the DB at `con`
#' (DELETE that label's rows, then append), inside one transaction.
check_saved_schema <- function(con) {
  stop_if_columns_changed(con, PREDICTIONS_TABLE, c("label", PREDICTION_COLUMNS))
  stop_if_columns_changed(con, RUNS_TABLE, RUN_COLUMNS)
}

save_backtest_run <- function(con, label, predictions, run_info) {
  predictions <- copy(predictions)
  predictions[, label := label]
  setcolorder(predictions, c("label", PREDICTION_COLUMNS))
  run_row <- as.data.table(c(list(label = label), run_info))
  if (!identical(names(run_row), RUN_COLUMNS)) {
    stop("save_backtest_run: run_info must have ", paste(RUN_COLUMNS[-1], collapse = ", "),
         " in that order; got ", paste(names(run_info), collapse = ", "))
  }
  check_saved_schema(con)
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
    con, sprintf("SELECT %s FROM %s WHERE label IN (%s)",
                 paste(c("label", PREDICTION_COLUMNS), collapse = ", "), PREDICTIONS_TABLE, placeholders),
    params = as.list(labels)))
  predictions[, game_date := as.Date(game_date)]
  runs <- as.data.table(dbGetQuery(
    con, sprintf("SELECT %s FROM %s WHERE label IN (%s)",
                 paste(RUN_COLUMNS, collapse = ", "), RUNS_TABLE, placeholders),
    params = as.list(labels)))
  if (!setequal(runs$label, labels)) {
    stop("backtest_runs has no row for label(s): ", paste(setdiff(labels, runs$label), collapse = ", "))
  }
  list(predictions = predictions, runs = runs[order(match(label, labels))])
}


# -----------------------------------------------------------------------------
# Report metrics (pure; input is predictions with a `label` column)
# -----------------------------------------------------------------------------

#' Keep lines every label priced with a graded (non-push) outcome, so labels
#' are compared on the same games and lines.
common_graded_lines <- function(predictions) {
  graded <- predictions[!is.na(model_prob) & !is.na(outcome)]
  empty <- setdiff(unique(predictions$label), graded$label)
  if (length(empty) > 0) stop("no graded lines for label(s): ", paste(empty, collapse = ", "))
  graded[, line_key := paste(odds_id, market, line)]
  n_labels <- uniqueN(predictions$label)
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
#' gap = sample's non-push cover (over) rate of the parent spread (total)
#' minus the no-vig target; mean_gap > 0 is a home (over) lean.
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
  ), by = .(label, market_group)]
}

#' Flat bets at the best closing price where the model's edge >= threshold.
#' At most one side per line: when best prices across books let both sides
#' clear it, only the bigger edge is bet. EV is per graded bet, so a line that
#' pushed is neither staked nor counted (pushes are dropped upstream).
flat_bets <- function(graded, edge_threshold = EDGE_THRESHOLD_DEFAULT) {
  h1 <- graded[market_group %in% H1_MARKETS]
  h1[, `:=`(side1_ev = model_prob * american_to_decimal(best_side1_american_odds) - 1,
            side2_ev = (1 - model_prob) * american_to_decimal(best_side2_american_odds) - 1)]
  side1 <- h1[side1_ev >= edge_threshold & side1_ev >= side2_ev, .(
    label, market_group, odds_id, claimed_edge = side1_ev,
    pnl = FLAT_STAKE * ifelse(outcome == 1, american_to_decimal(best_side1_american_odds) - 1, -1))]
  side2 <- h1[side2_ev >= edge_threshold & side2_ev > side1_ev, .(
    label, market_group, odds_id, claimed_edge = side2_ev,
    pnl = FLAT_STAKE * ifelse(outcome == 0, american_to_decimal(best_side2_american_odds) - 1, -1))]
  rbind(side1, side2)
}

summarise_flat_bets <- function(bets, by = c("label", "market_group")) {
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

#' Print the whole report; runs' label order sets the order (first = baseline).
print_backtest_report <- function(predictions, runs, edge_threshold = EDGE_THRESHOLD_DEFAULT) {
  labels <- runs$label
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
  cat(sprintf("\n== Flat $%d bets at best closing price, claimed edge >= %.1f%% (1H) ==\n",
              FLAT_STAKE, 100 * edge_threshold))
  if (nrow(bets) == 0) {
    cat("no line cleared the edge threshold\n")
  } else {
    print(summarise_flat_bets(bets)[order(market_group, match(label, labels))], digits = 3)
    print(order_labels(summarise_flat_bets(bets, by = "label")), digits = 3)
  }
  invisible(graded)
}
