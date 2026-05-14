library(testthat)
source("../book_cell.R")  # adjust path: from tests/ dir up to book_cell.R

test_that("render_book_cell empty state (no quote)", {
  out <- render_book_cell(american_odds = NA_integer_, line_quoted = NA_real_,
                          is_exact_line = NA, is_pick = FALSE,
                          side_word = "over", is_totals = TRUE)
  expect_match(out, 'class="cell empty"')
  expect_match(out, '&mdash;')
})

test_that("render_book_cell exact state (book matches model line, positive odds)", {
  out <- render_book_cell(american_odds = 160L, line_quoted = -2.5,
                          is_exact_line = TRUE, is_pick = FALSE,
                          side_word = "over", is_totals = FALSE)
  expect_match(out, 'class="cell exact"')
  expect_match(out, '\\+160')
  # No alt-line tag rendered when is_exact_line == TRUE
  expect_no_match(out, 'class="alt-line"')
})

test_that("render_book_cell alt state for spread (signed line tag)", {
  out <- render_book_cell(american_odds = -174L, line_quoted = 1.5,
                          is_exact_line = FALSE, is_pick = FALSE,
                          side_word = "under", is_totals = FALSE)
  expect_match(out, 'class="cell alt"')
  expect_match(out, '<span class="alt-line">\\+1.5</span>')
  expect_match(out, '-174')
})

test_that("render_book_cell alt state for totals (O/U prefix on line tag)", {
  out <- render_book_cell(american_odds = -115L, line_quoted = 7.5,
                          is_exact_line = FALSE, is_pick = FALSE,
                          side_word = "over", is_totals = TRUE)
  expect_match(out, '<span class="alt-line">O7.5</span>')
})

test_that("render_book_cell pick state overrides others", {
  out <- render_book_cell(american_odds = 160L, line_quoted = -2.5,
                          is_exact_line = TRUE, is_pick = TRUE,
                          side_word = "over", is_totals = FALSE)
  expect_match(out, 'class="cell pick"')
  expect_match(out, '\\+160')
})
