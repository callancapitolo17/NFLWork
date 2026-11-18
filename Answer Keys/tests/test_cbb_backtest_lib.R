# Answer Keys/tests/test_cbb_backtest_lib.R
# Guards the CBB walk-forward backtest (CBB Answer Key/cbb_backtest_lib.R):
# data loading without the old .x/.y join clash, 1H line pairing and devig,
# per-game pricing rows, and the report metrics.
#
# Run from "Answer Keys/" directory:
#   Rscript -e 'testthat::test_file("tests/test_cbb_backtest_lib.R")'

library(testthat)
suppressMessages({ library(data.table); library(DBI); library(duckdb) })
source("../Tools.R")
source("../CBB Answer Key/cbb_backtest_lib.R")

.tiny_db <- function() {
  path <- tempfile(fileext = ".duckdb")
  con <- dbConnect(duckdb(), dbdir = path)
  dbWriteTable(con, "cbb_betting_pbp", data.frame(
    odds_id = c("e1", "e2", "e3"), game_id = c("g1", "g2", "g3"),
    game_date = as.Date(c("2024-01-01", "2024-01-02", "2024-01-03")),
    home_team = c("A", "C", "E"), away_team = c("B", "D", "F"),
    home_spread = c(-3, -4.5, 2), total_line = c(140, 150.5, 133),
    consensus_devig_home_odds = c(0.5, 0.52, 0.48), consensus_devig_away_odds = c(0.5, 0.48, 0.52),
    consensus_devig_over_odds = c(0.5, 0.49, 0.51), consensus_devig_under_odds = c(0.5, 0.51, 0.49),
    game_home_margin_fg = c(3, 10, NA), game_total_fg = c(140, 150, 130),
    game_home_margin_h1 = c(0, 4, 1), game_total_h1 = c(70, 72, 60)))
  quote <- function(event, book, market, outcome, price, point, home, away) data.frame(
    event_id = event, home_team = home, away_team = away, commence_time = "2024-01-01T00:00:00Z",
    bookmaker_key = book, market = market, outcome_name = outcome,
    outcome_price = as.integer(price), outcome_point = point)
  dbWriteTable(con, "cbb_derivative_closing_odds", rbind(
    quote("e1", "bk1", "spreads_h1", "A", -110, -1.5, "A", "B"),
    quote("e1", "bk1", "spreads_h1", "B", -110,  1.5, "A", "B"),
    quote("e1", "bk2", "spreads_h1", "A", -120, -1.5, "A", "B"),
    quote("e1", "bk2", "spreads_h1", "B",  100,  1.5, "A", "B"),
    quote("e1", "bk1", "totals_h1", "Over", -105, 68.5, "A", "B"),
    quote("e1", "bk1", "totals_h1", "Under", -115, 68.5, "A", "B"),  # one book only: dropped
    quote("e1", "bk1", "h2h_h1", "A", -150, NA, "A", "B"),
    quote("e1", "bk1", "h2h_h1", "B",  130, NA, "A", "B"),
    quote("e1", "bk2", "h2h_h1", "A", -140, NA, "A", "B"),
    quote("e1", "bk2", "h2h_h1", "B",  120, NA, "A", "B")))
  list(con = con, path = path)
}

test_that("history loads off cbb_betting_pbp's own lines, one row per complete game", {
  db <- .tiny_db(); on.exit(dbDisconnect(db$con, shutdown = TRUE))
  history <- load_cbb_history(db$con)
  expect_equal(history$odds_id, c("e1", "e2"))  # e3 has no final margin
  expect_true(all(c("home_spread", "total_line") %in% names(history)))
  expect_false(any(grepl("\\.(x|y)$", names(history))))
  # e1: margin 3 on -3 is a push; total 140 on 140 is a push
  expect_equal(history$actual_cover, c(0.5, 1))
  expect_equal(history$actual_over, c(0.5, 0))
  expect_s3_class(history$game_date, "Date")
})

test_that("1H lines pair each book's two sides, devig by probit, need 2+ books", {
  db <- .tiny_db(); on.exit(dbDisconnect(db$con, shutdown = TRUE))
  lines <- load_h1_closing_lines(db$con)
  expect_setequal(lines$market, c("spreads_h1", "h2h_h1"))
  spread <- lines[market == "spreads_h1"]
  expect_equal(spread$line, -1.5)
  expect_equal(spread$n_books, 2L)
  fair_bk1 <- devig_american(-110, -110)$p1
  fair_bk2 <- devig_american(-120, 100)$p1
  expect_equal(spread$book_fair_prob, median(c(fair_bk1, fair_bk2)), tolerance = 1e-12)
  expect_equal(spread$best_side1_american_odds, -110)
  expect_equal(spread$best_side2_american_odds, 100)
  expect_equal(lines[market == "h2h_h1"]$line, 0)
})

test_that("a single named book keeps its one-book lines; an unknown book fails loudly", {
  db <- .tiny_db(); on.exit(dbDisconnect(db$con, shutdown = TRUE))
  lines <- load_h1_closing_lines(db$con, books = "bk1")
  expect_setequal(lines$market, c("spreads_h1", "totals_h1", "h2h_h1"))
  expect_error(load_h1_closing_lines(db$con, books = "nope"), "no 1H prices for book")
})

test_that("sample probabilities come from Tools.R's predictors and drop pushes", {
  sample <- data.table(game_home_margin_h1 = c(-3, 0, 2, 3, 5), game_total_h1 = c(60, 65, 65, 70, 71))
  # home -2.5 covers when margin > 2.5: {3, 5} of 5; home 0 (moneyline): {2, 3, 5} of 4 non-ties
  expect_equal(sample_home_cover_probs(sample, c(-2.5, 0, -2.5), "h1"), c(2 / 5, 3 / 4, 2 / 5))
  # over 65 drops the two 65s: {70, 71} of 3
  expect_equal(sample_over_probs(sample, c(65, 80), "h1"), c(2 / 3, 0))
  expect_true(is.na(sample_home_cover_probs(data.table(game_home_margin_h1 = c(1, 1)), -1, "h1")))
  expect_equal(outcome_vs_threshold(c(4, 3, 2, NA), 3), c(1, NA, 0, NA))
})

.mock_sampler <- function(result_fn) {
  real_sampler <- run_answer_key_sample
  assign("run_answer_key_sample", result_fn, envir = globalenv())
  function() assign("run_answer_key_sample", real_sampler, envir = globalenv())
}

.four_game_history <- function() {
  data.table(
    odds_id = sprintf("e%d", 1:4), game_date = as.Date("2024-01-01") + 0:3,
    home_spread = c(-2, -3, -1, -2), total_line = c(140, 150, 135, 140),
    consensus_devig_home_odds = 0.5, consensus_devig_over_odds = 0.5,
    home_margin = c(5, -5, 1, 3), total_final_score = c(150, 130, 141, 139),
    game_home_margin_h1 = c(2, -1, 0, 0), game_total_h1 = c(75, 64, 70, 69))
}

test_that("a test game is priced only off earlier games and grades each line", {
  history <- .four_game_history()
  history[, `:=`(game_home_margin_fg = home_margin, game_total_fg = total_final_score)]
  seen_train_dates <- NULL
  restore <- .mock_sampler(function(id, parent_spread, parent_total, target_cover, target_over,
                                    DT, ss, st, N, ...) {
    seen_train_dates <<- DT$game_date
    list(sample = DT, final_N = nrow(DT), converged = TRUE)
  })
  on.exit(restore())
  h1 <- data.table(event_id = "e4", market = c("spreads_h1", "h2h_h1", "totals_h1"),
                   line = c(-0.5, 0, 69.5), book_fair_prob = c(0.5, 0.55, 0.5), n_books = 2L,
                   best_side1_american_odds = c(-105, -120, -110),
                   best_side2_american_odds = c(-105, 110, -110))
  priced <- price_test_game(history[4], history, h1, sample_pct = 1)
  expect_equal(priced$status, "priced")
  rows <- priced$predictions
  expect_identical(names(rows), PREDICTION_COLUMNS)
  expect_true(all(seen_train_dates < as.Date("2024-01-04")))
  spread <- rows[market == "spreads_h1"]
  expect_equal(spread$model_prob, 1 / 3)  # 1H margins 2, -1, 0: only 2 covers -0.5
  expect_equal(spread$outcome, 0)         # test game's 1H margin 0 does not cover -0.5
  ml <- rows[market == "h2h_h1"]
  expect_equal(ml$model_prob, 1 / 2)      # 2 and -1; the tied half drops out
  expect_equal(ml$outcome, NA_real_)      # tied half = push
  total <- rows[market == "totals_h1"]
  expect_equal(total$model_prob, 2 / 3)   # 75, 70 over 69.5; 64 under
  expect_equal(total$outcome, 0)          # 69 < 69.5
  main <- rows[market == "fg_spread_+0"]
  expect_equal(main$book_fair_prob, 0.5)
  expect_equal(main$outcome, 1)           # margin 3 beats -2
  expect_true(all(is.na(rows[market == "fg_spread_+3"]$book_fair_prob)))
  expect_equal(nrow(rows), 3 + length(FG_SPREAD_OFFSETS) + length(FG_TOTAL_OFFSETS))
  expect_equal(rows$target_N[1], 3)
  expect_equal(rows$sample_median_margin_miss[1], median(c(5, -5, 1)) - 2)
})

test_that("sampler failures are told apart: no convergence is skipped, an error is kept", {
  history <- .four_game_history()
  history[, `:=`(game_home_margin_fg = home_margin, game_total_fg = total_final_score)]
  restore <- .mock_sampler(function(...) NULL)
  expect_equal(price_test_game(history[4], history, data.table(), sample_pct = 1)$status, "not_converged")
  restore()
  restore <- .mock_sampler(function(...) stop("boom"))
  on.exit(restore())
  failed <- price_test_game(history[4], history, data.table(), sample_pct = 1)
  expect_equal(failed[c("status", "message")], list(status = "error", message = "boom"))
  expect_equal(price_test_game(history[1], history, data.table())$status, "no_history")
})

test_that("a sampler error anywhere stops the whole run, naming the game", {
  history <- .four_game_history()
  history[, `:=`(game_home_margin_fg = home_margin, game_total_fg = total_final_score)]
  restore <- .mock_sampler(function(id, ...) if (id == "e4") stop("bad game") else NULL)
  on.exit(restore())
  h1 <- data.table(event_id = c("e3", "e4"), market = "h2h_h1", line = 0, book_fair_prob = 0.5,
                   n_books = 2L, best_side1_american_odds = -110, best_side2_american_odds = -110)
  setkey(h1, event_id)
  expect_error(price_test_games(history, h1, 3:4, sample_pct = 1, cores = 1), "e4\\): bad game")
})

test_that("test games need half the history behind them and 1H lines; max_games subsets", {
  history <- data.table(odds_id = sprintf("e%d", 1:10), game_date = as.Date("2024-01-01") + 0:9)
  h1 <- data.table(event_id = c("e2", "e7", "e8", "e9"))
  expect_equal(select_test_games(history, h1), c(7L, 8L, 9L))
  expect_length(select_test_games(history, h1, max_games = 2), 2)
  expect_equal(select_test_games(history, h1, max_games = 50), c(7L, 8L, 9L))
})

.graded_fixture <- function() {
  base <- data.table(odds_id = c("e1", "e1", "e2"), game_date = as.Date("2025-01-01"),
                     market = c("spreads_h1", "fg_spread_+0", "totals_h1"), line = c(-1.5, -3, 70),
                     model_prob = c(0.60, 0.51, 0.40), book_fair_prob = c(0.50, 0.50, 0.50),
                     best_side1_american_odds = c(100, NA, -110), best_side2_american_odds = c(-120, NA, 110),
                     outcome = c(1, 1, 0), target_N = 100, final_N = 100, converged = TRUE,
                     sample_median_margin_miss = 0.5, sample_median_total_miss = -0.5,
                     sample_mean_spread_drift = 0)
  setcolorder(base, PREDICTION_COLUMNS)
  rbind(copy(base)[, label := "a"], copy(base)[, label := "b"][market != "totals_h1"])
}

.run_info <- function() list(run_at = Sys.time(), code_dir = "x", code_commit = "abc", sample_pct = 0.02,
                             max_games = 0L, books = "all", games_priced = 2L,
                             games_not_converged = 0L, minutes = 1)

test_that("labels are compared only on lines every label priced", {
  graded <- common_graded_lines(.graded_fixture())
  expect_setequal(unique(graded$market), c("spreads_h1", "fg_spread_+0"))
  expect_equal(graded[market == "fg_spread_+0"]$market_group[1], "fg_spread_main")
})

test_that("a label with no graded lines stops the report by name", {
  preds <- .graded_fixture()
  preds[label == "b", model_prob := NA_real_]
  expect_error(common_graded_lines(preds), "no graded lines for label\\(s\\): b")
})

test_that("flat bets take the side with the edge at the best price", {
  graded <- common_graded_lines(.graded_fixture()[label == "a"])
  bets <- flat_bets(graded, edge_threshold = 0.03)
  # spreads_h1: 0.60 * 2.0 - 1 = 0.20 on side 1, wins +100
  # totals_h1:  0.60 * 2.1 - 1 = 0.26 on side 2 (under), outcome 0 wins +110
  expect_equal(nrow(bets), 2)
  expect_equal(bets[market_group == "spreads_h1"]$pnl, 100)
  expect_equal(bets[market_group == "totals_h1"]$pnl, 110)
  summary <- summarise_flat_bets(bets, by = "label")
  expect_equal(summary$roi, 1.05)
  expect_equal(summary$bets, 2L)
})

test_that("a line where best prices let both sides clear the edge bets only the bigger one", {
  graded <- common_graded_lines(.graded_fixture()[label == "a"])
  graded[market == "spreads_h1", `:=`(model_prob = 0.52, best_side1_american_odds = 110,
                                      best_side2_american_odds = 110)]
  bets <- flat_bets(graded, edge_threshold = 0.03)[market_group == "spreads_h1"]
  expect_equal(nrow(bets), 1)
  expect_equal(bets$claimed_edge, 0.52 * 2.1 - 1)
})

test_that("log loss change is to minus from, so negative means the second label is better", {
  preds <- .graded_fixture()[market != "totals_h1"]
  preds[label == "b" & market == "spreads_h1", model_prob := 0.9]  # outcome 1: b is better
  graded <- common_graded_lines(preds)
  change <- compare_logloss(graded, "a", "b")
  expect_equal(change[market_group == "spreads_h1"]$logloss_change, log(0.6) - log(0.9))
  expect_equal(change[market_group == "fg_spread_main"]$logloss_change, 0)
})

test_that("centring and accuracy summaries report the gap to the no-vig target", {
  graded <- common_graded_lines(.graded_fixture()[label == "a"])
  centring <- summarise_sample_centring(graded)
  expect_equal(centring$main_line_gaps$mean_gap, 0.01)
  expect_equal(centring$diagnostics$mean_abs_median_margin_miss, 0.5)
  accuracy <- summarise_accuracy(graded)
  expect_equal(accuracy[market_group == "spreads_h1"]$model_logloss, -log(0.6))
  expect_equal(accuracy[market_group == "spreads_h1"]$book_logloss, -log(0.5))
})

test_that("the report prints for two labels and with no bets over the edge", {
  preds <- .graded_fixture()
  runs <- data.table(label = c("a", "b"))
  expect_output(print_backtest_report(preds, runs, edge_threshold = 0.03), "Paired log loss change vs a")
  expect_output(print_backtest_report(preds, runs, edge_threshold = 0.9), "no line cleared the edge threshold")
})

test_that("saving a label replaces only that label's rows", {
  path <- tempfile(fileext = ".duckdb")
  con <- dbConnect(duckdb(), dbdir = path); on.exit(dbDisconnect(con, shutdown = TRUE))
  preds <- .graded_fixture()[label == "a"][, label := NULL]
  save_backtest_run(con, "a", preds, .run_info())
  save_backtest_run(con, "b", preds, .run_info())
  save_backtest_run(con, "a", preds[1], .run_info())
  saved <- load_backtest_runs(con, c("b", "a"))
  expect_equal(nrow(saved$predictions[label == "a"]), 1)
  expect_equal(nrow(saved$predictions[label == "b"]), 3)
  expect_equal(saved$runs$label, c("b", "a"))
  expect_error(load_backtest_runs(con, "zzz"), "no saved run for label")
})

test_that("a saved table with other columns stops the run before pricing", {
  path <- tempfile(fileext = ".duckdb")
  con <- dbConnect(duckdb(), dbdir = path); on.exit(dbDisconnect(con, shutdown = TRUE))
  dbWriteTable(con, RUNS_TABLE, data.frame(label = "old", run_at = Sys.time()))
  expect_error(check_saved_schema(con), "backtest_runs was saved with columns")
  info <- .run_info(); info$books <- NULL
  expect_error(save_backtest_run(dbConnect(duckdb(), dbdir = tempfile(fileext = ".duckdb")), "a",
                                 .graded_fixture()[label == "a"][, label := NULL], info),
               "run_info must have")
})
