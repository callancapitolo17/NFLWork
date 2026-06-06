# Answer Keys/tests/test_nflwork_root.R
library(testthat)
source("../Tools.R")

test_that("nflwork_root falls back to ~/NFLWork when unset", {
  set_nflwork_root(NULL)
  expect_equal(nflwork_root(), normalizePath(path.expand("~/NFLWork"), mustWork = FALSE))
})

test_that("set_nflwork_root overrides the root", {
  tmp <- tempfile("root_"); dir.create(tmp)
  set_nflwork_root(tmp)
  expect_equal(nflwork_root(), normalizePath(tmp, mustWork = FALSE))
  set_nflwork_root(NULL)  # reset for other tests
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
