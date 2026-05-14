# Answer Keys/tests/test_book_cell.R
# Tests for book_cell.R: .devig_american_pair (probit math) and
# render_book_cell (HTML output with raw + fair spans).
#
# Run from "Answer Keys/" directory:
#   Rscript -e 'testthat::test_file("tests/test_book_cell.R")'
#
# NOTE: Expected American-odds values here match Tools.R::devig_american()
# rounded to the nearest integer. At p = 0.5 exactly, the convention is
# -100 (not +100); at p = 0.5645 (the asymmetric case below), the unrounded
# American odds are -129.60 / +129.60 which round to -130 / +130. The
# parity test in test_devig_pair_matches_tools.R independently enforces
# equivalence with Tools.R across a wider set of inputs.

library(testthat)
source("../MLB Dashboard/book_cell.R")

# Helper: convert American odds to implied probability for the assertion math.
.amer_to_prob <- function(o) ifelse(o > 0, 100 / (o + 100), -o / (-o + 100))

test_that(".devig_american_pair returns devigged American odds (no-vig coin flip)", {
  # -110 / -110 -> both implied 52.38%, devigged -> exactly 50% / 50%.
  # Convention: p >= 0.5 emits negative-sign American odds, so p=0.5 -> -100.
  out <- .devig_american_pair(-110, -110)
  expect_equal(out$fair1, -100, tolerance = 0.5)
  expect_equal(out$fair2, -100, tolerance = 0.5)
})

test_that(".devig_american_pair returns devigged American odds (asymmetric)", {
  # -140 / +120 -> implied 58.3% / 45.5%, sum 103.8% (vig)
  # Probit-devigged -> 56.4468% / 43.5532% -> unrounded American
  # -129.604 / +129.604 -> round() -> -130 / +130
  out <- .devig_american_pair(-140, +120)
  expect_lt(abs(out$fair1 - (-130)), 1)
  expect_lt(abs(out$fair2 - (+130)), 1)
})

test_that(".devig_american_pair returns NA on bad inputs", {
  expect_true(is.na(.devig_american_pair(NA, +120)$fair1))
  expect_true(is.na(.devig_american_pair(0,  +120)$fair1))
  expect_true(is.na(.devig_american_pair(-110, NA)$fair2))
})
