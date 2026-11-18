# =============================================================================
# CBB walk-forward backtest of the answer key
# =============================================================================
# Prices each test game off only the games dated before it (as live does),
# saves every prediction under a label, and prints the report: sample centring
# on the parent line, 1H spread / total / moneyline log loss vs the book close,
# full-game main and alt lines, and flat-bet ROI at the best closing price.
# Metrics and loading live in cbb_backtest_lib.R.
#
# Usage (from anywhere):
#   Rscript "Answer Keys/CBB Answer Key/CBB_Backtest.R" run --label main
#   Rscript "Answer Keys/CBB Answer Key/CBB_Backtest.R" run --label before --code-dir /path/to/worktree/"Answer Keys"
#   Rscript "Answer Keys/CBB Answer Key/CBB_Backtest.R" report --labels before,main
#
# run options:
#   --label NAME        required; a rerun replaces that label's saved rows
#   --code-dir DIR      "Answer Keys" folder whose Tools.R is the code under
#                       test (default: this repo's). Point it at a worktree of
#                       another commit to compare code versions on the same data.
#   --max-games N       random subset of the eligible games (default 0 = all,
#                       ~9,600 games, ~3 min on 4 cores)
#   --sample-pct P      N as a share of each game's prior history (default 0.02, as CBB.R)
#   --cores N           parallel workers (default 4)
#   --books a,b         1H fair and best price from these bookmaker_keys only
#                       (e.g. betonline); default every book, a line needing 2+
#   --edge E            claimed-edge threshold for the flat-bet table (default 0.03)
# report options: --labels a,b,... (first is the baseline), --edge E
#
# Reads (read-only): this repo's "Answer Keys/cbb.duckdb" — cbb_betting_pbp
#   (history, parent lines, no-vig targets) and cbb_derivative_closing_odds
#   (1H closing prices; a line needs 2+ books).
# Writes: "Answer Keys/cbb_backtest.duckdb" (gitignored) — backtest_predictions
#   and backtest_runs; `run` deletes then appends the rows of its own label.
# =============================================================================

suppressMessages({
  library(data.table)
  library(DBI)
  library(duckdb)
})

script_path <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0) stop("run CBB_Backtest.R with Rscript so it can find its own folder")
  normalizePath(gsub("~\\+~", " ", sub("^--file=", "", file_arg[1])))
})
answer_keys_dir <- dirname(dirname(script_path))
DATA_DB_PATH <- file.path(answer_keys_dir, "cbb.duckdb")
BACKTEST_DB_PATH <- file.path(answer_keys_dir, "cbb_backtest.duckdb")
source(file.path(dirname(script_path), "cbb_backtest_lib.R"))

parse_cli <- function(args) {
  if (length(args) == 0 || !args[1] %in% c("run", "report")) {
    stop("usage: CBB_Backtest.R run --label NAME [options] | report --labels a,b [--edge E]")
  }
  opts <- list(command = args[1], code_dir = answer_keys_dir, max_games = 0,
               sample_pct = SAMPLE_PCT_DEFAULT, cores = 4, edge = EDGE_THRESHOLD_DEFAULT)
  rest <- args[-1]
  if (length(rest) %% 2 != 0) stop("every option needs a value, got: ", paste(rest, collapse = " "))
  for (i in seq(1, length(rest), by = 2)) {
    key <- sub("^--", "", rest[i]); value <- rest[i + 1]
    switch(key,
      "label"      = opts$label <- value,
      "labels"     = opts$labels <- strsplit(value, ",", fixed = TRUE)[[1]],
      "code-dir"   = opts$code_dir <- normalizePath(value, mustWork = TRUE),
      "max-games"  = opts$max_games <- as.integer(value),
      "sample-pct" = opts$sample_pct <- as.numeric(value),
      "cores"      = opts$cores <- as.integer(value),
      "books"      = opts$books <- strsplit(value, ",", fixed = TRUE)[[1]],
      "edge"       = opts$edge <- as.numeric(value),
      stop("unknown option --", key))
  }
  if (opts$command == "run" && is.null(opts$label)) stop("run needs --label")
  if (opts$command == "report" && is.null(opts$labels)) stop("report needs --labels a,b,...")
  opts
}

# Tools.R compiles src/sampling.cpp relative to its own folder.
source_code_under_test <- function(code_dir) {
  tools_path <- file.path(code_dir, "Tools.R")
  if (!file.exists(tools_path)) stop("no Tools.R in --code-dir ", code_dir)
  old_wd <- setwd(code_dir)
  on.exit(setwd(old_wd))
  source(tools_path)
  needed <- c("run_answer_key_sample", "compute_dispersion",
              "walk_forward_training_set", "has_walk_forward_history")
  missing <- needed[!vapply(needed, exists, logical(1))]
  if (length(missing) > 0) {
    stop(tools_path, " lacks ", paste(missing, collapse = ", "),
         "; walk-forward needs a version at or after PR #332")
  }
}

git_commit_of <- function(dir) {
  commit <- tryCatch(system2("git", c("-C", shQuote(dir), "rev-parse", "--short", "HEAD"),
                             stdout = TRUE, stderr = FALSE), error = function(e) character(0))
  if (length(commit) == 0) NA_character_ else commit[1]
}

run_command <- function(opts) {
  source_code_under_test(opts$code_dir)
  data_con <- dbConnect(duckdb(), dbdir = DATA_DB_PATH, read_only = TRUE)
  history <- load_cbb_history(data_con)
  h1_lines <- load_h1_closing_lines(data_con, books = opts$books)
  dbDisconnect(data_con, shutdown = TRUE)

  test_idx <- select_test_games(history, h1_lines, opts$max_games)
  cat(sprintf("%s: %d games of history, pricing %d walk-forward games with 1H closing lines\n",
              opts$label, nrow(history), length(test_idx)))
  started <- Sys.time()
  priced <- price_test_games(history, h1_lines, test_idx, opts$sample_pct, opts$cores)
  minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  cat(sprintf("priced %d games (%d the sampler could not price) in %.1f min\n",
              uniqueN(priced$predictions$odds_id), priced$unpriced_games, minutes))
  if (nrow(priced$predictions) == 0) stop("no game could be priced; nothing saved")

  run_info <- list(run_at = Sys.time(), code_dir = opts$code_dir,
                   code_commit = git_commit_of(opts$code_dir), sample_pct = opts$sample_pct,
                   max_games = opts$max_games,
                   books = if (is.null(opts$books)) "all" else paste(opts$books, collapse = ","),
                   games_priced = uniqueN(priced$predictions$odds_id),
                   games_unpriced = priced$unpriced_games, minutes = round(minutes, 1))
  out_con <- dbConnect(duckdb(), dbdir = BACKTEST_DB_PATH)
  on.exit(dbDisconnect(out_con, shutdown = TRUE))
  save_backtest_run(out_con, opts$label, priced$predictions, run_info)
  cat(sprintf("saved label '%s' to %s\n", opts$label, BACKTEST_DB_PATH))
  saved <- load_backtest_runs(out_con, opts$label)
  print_backtest_report(saved$predictions, saved$runs, opts$edge)
}

report_command <- function(opts) {
  con <- dbConnect(duckdb(), dbdir = BACKTEST_DB_PATH, read_only = TRUE)
  on.exit(dbDisconnect(con, shutdown = TRUE))
  saved <- load_backtest_runs(con, opts$labels)
  saved$runs <- saved$runs[order(match(label, opts$labels))]
  print_backtest_report(saved$predictions, saved$runs, opts$edge)
}

opts <- parse_cli(commandArgs(trailingOnly = TRUE))
if (opts$command == "run") run_command(opts) else report_command(opts)
