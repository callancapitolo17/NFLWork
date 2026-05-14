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

# ----------------------------------------------------------------------------
# render_book_cell HTML output tests (PR B Task 2)
# ----------------------------------------------------------------------------
# These tests verify that render_book_cell emits BOTH a <span class="raw">
# (the book's actual American odds) and a <span class="fair"> (the
# probit-devigged fair American odds) when the opposite side's odds are
# supplied. Later PR-B tasks add CSS + JS to hide one or the other based
# on a class on the parent .price-grid container (the toggle).
#
# Backwards-compat: when opposite_american_odds is NA (legacy callers
# like render_book_pill), no fair span should be emitted — only raw.
#
# NOTE on expected fair values: .devig_american_pair(+120, -140) returns
# fair1 = +130 (NOT +129 as the original PR-B plan documented). The plan's
# literal values were off by 1 due to rounding (+129.604 rounds to +130).
# The asymmetric devig math test above already documents this.

test_that("render_book_cell with both odds emits raw + fair spans", {
  html <- render_book_cell(
    american_odds          = +120,
    opposite_american_odds = -140,
    line_quoted            = -0.5,
    is_exact_line          = TRUE,
    is_pick                = FALSE,
    side_word              = "over",
    is_totals              = FALSE
  )
  expect_match(html, '<span class="raw">\\+120</span>',  fixed = FALSE)
  expect_match(html, '<span class="fair">\\+130</span>', fixed = FALSE)
  # No alt-line tag because is_exact_line = TRUE
  expect_false(grepl('alt-line', html))
})

test_that("render_book_cell without opposite_american_odds emits raw + fair em-dash (no opposite available)", {
  html <- render_book_cell(
    american_odds = +120,
    line_quoted   = -0.5,
    is_exact_line = TRUE,
    is_pick       = FALSE,
    side_word     = "over",
    is_totals     = FALSE
  )
  # Raw span present; fair span renders the em-dash fallback (no devig
  # computable without the opposite side's odds). Guarantees the Task 5
  # CSS toggle (.show-fair hides .raw) never produces a blank cell.
  expect_match(html, '<span class="raw">\\+120</span>', fixed = FALSE)
  expect_match(html, '<span class="fair">&mdash;</span>', fixed = TRUE)
})

test_that("render_book_cell pick cell with both odds emits green pick class and both spans, no alt tag", {
  html <- render_book_cell(
    american_odds          = +120,
    opposite_american_odds = -140,
    line_quoted            = -0.5,
    is_exact_line          = TRUE,
    is_pick                = TRUE,
    side_word              = "over",
    is_totals              = FALSE
  )
  expect_match(html, 'class="cell pick"',                  fixed = TRUE)
  expect_match(html, '<span class="raw">\\+120</span>',    fixed = FALSE)
  expect_match(html, '<span class="fair">\\+130</span>',   fixed = FALSE)
  expect_false(grepl('alt-line', html))
})

test_that("render_book_cell formats negative fair odds without a + prefix", {
  # Inverse asymmetric input: -140 raw should devig fair to -130 (negative side
  # of the same pair as the +130 test above). Confirms the fair_str else-branch
  # used for negative odds emits "-130" rather than "+-130" or "+130".
  html <- render_book_cell(
    american_odds          = -140,
    opposite_american_odds = +120,
    line_quoted            = -0.5,
    is_exact_line          = TRUE,
    is_pick                = FALSE,
    side_word              = "under",
    is_totals              = FALSE
  )
  expect_match(html, '<span class="raw">-140</span>',  fixed = FALSE)
  expect_match(html, '<span class="fair">-130</span>', fixed = FALSE)
})

test_that("render_book_cell on alt line keeps amber tag in both views", {
  html <- render_book_cell(
    american_odds          = -117,
    opposite_american_odds = +102,
    line_quoted            = 0,
    is_exact_line          = FALSE,
    is_pick                = FALSE,
    side_word              = "over",
    is_totals              = FALSE
  )
  expect_match(html, 'alt-line', fixed = FALSE)
  expect_match(html, '<span class="raw">-117</span>',  fixed = FALSE)
  expect_match(html, '<span class="fair">', fixed = FALSE)  # devig present
})
