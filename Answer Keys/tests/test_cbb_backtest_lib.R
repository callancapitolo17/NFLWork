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

test_that("devig_probit_side1 matches Tools.R::devig_american", {
  odds1 <- c(-110, 150, -250, 105); odds2 <- c(-110, -170, 200, -125)
  expect_equal(devig_probit_side1(odds1, odds2), devig_american(odds1, odds2)$p1, tolerance = 1e-12)
  expect_true(is.na(devig_probit_side1(0, -110)))
})

test_that("sample probabilities drop pushes and need enough games", {
  expect_equal(sample_prob_over(c(1, 2, 3, 3, 4), 3, min_nonpush = 1), 1 / 3)  # 4 of {1, 2, 4}
  expect_true(is.na(sample_prob_over(c(1, 5), 3, min_nonpush = 3)))
  expect_equal(outcome_vs_threshold(c(4, 3, 2, NA), 3), c(1, NA, 0, NA))
})

test_that("a test game is priced only off earlier games and grades each line", {
  history <- data.table(
    odds_id = sprintf("e%d", 1:4), game_date = as.Date("2024-01-01") + 0:3,
    home_spread = -2, total_line = 140, consensus_devig_home_odds = 0.5, consensus_devig_over_odds = 0.5,
    home_margin = c(5, -5, 1, 3), total_final_score = c(150, 130, 141, 139),
    game_home_margin_h1 = c(2, -1, 0, 0), game_total_h1 = c(75, 64, 70, 69))
  seen_train_dates <- NULL
  local_mocked <- function(id, parent_spread, parent_total, target_cover, target_over, DT, ss, st, N, ...) {
    seen_train_dates <<- DT$game_date
    list(sample = DT, final_N = nrow(DT), converged = TRUE)
  }
  real_sampler <- run_answer_key_sample
  assign("run_answer_key_sample", local_mocked, envir = globalenv())
  on.exit(assign("run_answer_key_sample", real_sampler, envir = globalenv()))
  h1 <- data.table(event_id = "e4", market = c("spreads_h1", "h2h_h1"), line = c(-0.5, 0),
                   book_fair_prob = c(0.5, 0.55), n_books = 2L,
                   best_side1_american_odds = c(-105, -120), best_side2_american_odds = c(-105, 110))
  rows <- price_test_game(history[4], history, h1, list(ss = 1, st = 1), sample_pct = 1)
  expect_true(all(seen_train_dates < as.Date("2024-01-04")))
  spread <- rows[market == "spreads_h1"]
  expect_equal(spread$model_prob, NA_real_)  # 3 sample games < MIN_NONPUSH_SAMPLE_GAMES
  expect_equal(spread$outcome, 0)            # home -0.5 in the half: margin 0 does not cover
  expect_equal(rows[market == "h2h_h1"]$outcome, NA_real_)  # tied half = push
  main <- rows[market == "fg_spread_+0"]
  expect_equal(main$book_fair_prob, 0.5)
  expect_equal(main$outcome, 1)  # margin 3 beats -2
  expect_equal(rows$target_N[1], 3)
  expect_equal(rows$sample_median_margin_miss[1], median(c(5, -5, 1)) - 2)
})

.graded_fixture <- function() {
  base <- data.table(odds_id = c("e1", "e1", "e2"), game_date = as.Date("2025-01-01"),
                     market = c("spreads_h1", "fg_spread_+0", "totals_h1"), line = c(-1.5, -3, 70),
                     model_prob = c(0.60, 0.51, 0.40), book_fair_prob = c(0.50, 0.50, 0.50),
                     best_side1_american_odds = c(100, NA, -110), best_side2_american_odds = c(-120, NA, 110),
                     n_books = 2L, outcome = c(1, 1, 0), parent_spread = -3, parent_total = 140,
                     target_N = 100, final_N = 100, converged = TRUE, sample_median_margin_miss = 0.5,
                     sample_median_total_miss = -0.5, sample_mean_spread_drift = 0)
  rbind(copy(base)[, label := "a"], copy(base)[, label := "b"][market != "totals_h1"])
}

test_that("labels are compared only on lines every label priced", {
  graded <- common_graded_lines(.graded_fixture())
  expect_setequal(unique(graded$market), c("spreads_h1", "fg_spread_+0"))
  expect_equal(graded[market == "fg_spread_+0"]$market_group[1], "fg_spread_main")
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
})

test_that("saving a label replaces only that label's rows", {
  path <- tempfile(fileext = ".duckdb")
  con <- dbConnect(duckdb(), dbdir = path); on.exit(dbDisconnect(con, shutdown = TRUE))
  preds <- .graded_fixture()[label == "a"][, label := NULL]
  info <- list(run_at = Sys.time(), code_dir = "x", code_commit = "abc", sample_pct = 0.02,
               max_games = 0, books = "all", games_priced = 2L, games_unpriced = 0L, minutes = 1)
  save_backtest_run(con, "a", preds, info)
  save_backtest_run(con, "b", preds, info)
  save_backtest_run(con, "a", preds[1], info)
  saved <- load_backtest_runs(con, c("a", "b"))
  expect_equal(nrow(saved$predictions[label == "a"]), 1)
  expect_equal(nrow(saved$predictions[label == "b"]), 3)
  expect_equal(nrow(saved$runs), 2)
  expect_error(load_backtest_runs(con, "zzz"), "no saved run for label")
})
