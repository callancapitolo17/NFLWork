# Answer Keys/tests/test_nflwork_root.R
library(testthat)
source("../Tools.R")

test_that("nflwork_root falls back to ~/NFLWork when unset", {
  set_nflwork_root(NULL)
  expect_equal(nflwork_root(), normalizePath(path.expand("~/NFLWork"), mustWork = FALSE))
})

test_that("set_nflwork_root overrides the root", {
  tmp <- tempfile("root_"); dir.create(tmp)
  on.exit(set_nflwork_root(NULL))
  set_nflwork_root(tmp)
  expect_equal(nflwork_root(), normalizePath(tmp, mustWork = FALSE))
})

test_that("derive_repo_root walks up to the dir containing 'Answer Keys'", {
  root <- tempfile("repo_"); dir.create(root)
  deep <- file.path(root, "Answer Keys", "MLB Answer Key")
  dir.create(deep, recursive = TRUE)
  script <- file.path(deep, "MLB.R")
  expect_equal(
    derive_repo_root(script_path = script),
    normalizePath(root, mustWork = FALSE)
  )
})

test_that("derive_repo_root falls back to ~/NFLWork when no marker found", {
  orphan <- file.path(tempfile("orphan_"))
  dir.create(orphan, recursive = TRUE)
  expect_equal(
    derive_repo_root(script_path = file.path(orphan, "x.R")),
    normalizePath(path.expand("~/NFLWork"), mustWork = FALSE)
  )
})

test_that("derive_repo_root falls back when no --file= in commandArgs", {
  # Called with no script_path under testthat, there is no --file= arg, so it
  # must hit the length(fa) == 0 fallback. (The ~+~ space-decode path is verified
  # manually via a real Rscript launch; it cannot be unit-tested without a subprocess.)
  expect_equal(
    derive_repo_root(),
    normalizePath(path.expand("~/NFLWork"), mustWork = FALSE)
  )
})

test_that("get_*_odds defaults route through nflwork_root()", {
  # Inspect each function's db_path DEFAULT source — deterministic, no dependency
  # on DB existence or warning text. (nflwork_root() correctness is proven above.)
  book_dirs <- list(
    get_wagerzon_odds  = "wagerzon_odds",
    get_hoop88_odds    = "hoop88_odds",
    get_bfa_odds       = "bfa_odds",
    get_bookmaker_odds = "bookmaker_odds",
    get_dk_odds        = "dk_odds",
    get_fd_odds        = "fd_odds",
    get_bet105_odds    = "bet105_odds",
    get_kalshi_odds    = "kalshi_odds"
  )
  for (fn in names(book_dirs)) {
    default_src <- paste(deparse(formals(get(fn))$db_path), collapse = " ")
    expect_match(default_src, "nflwork_root()", fixed = TRUE, info = fn)
    expect_match(default_src, book_dirs[[fn]], fixed = TRUE, info = fn)
  }
})
