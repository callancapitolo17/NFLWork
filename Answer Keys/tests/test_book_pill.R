# Answer Keys/tests/test_book_pill.R
# Tests for render_book_pill — bets-tab pill renderer with 3 visual states + pick override.
# Run from "Answer Keys/" directory:
#   Rscript -e 'testthat::test_file("tests/test_book_pill.R")'

library(testthat)
source("../MLB Dashboard/book_pill.R")

test_that("exact-line pill renders book label + price only", {
  out <- render_book_pill(book = "WZ", american_odds = 120,
                         line_quoted = 5.5, is_exact_line = TRUE, is_pick = FALSE)
  expect_match(out, 'class="pill"', fixed = TRUE)
  expect_match(out, '<span class="book">WZ</span>', fixed = TRUE)
  expect_match(out, '>+120<', fixed = TRUE)
  expect_false(grepl("line-tag", out))
})

test_that("negative odds render with minus sign, not plus", {
  out <- render_book_pill("FD", -115, 5.5, is_exact_line = TRUE, FALSE)
  expect_match(out, '>-115<', fixed = TRUE)
  expect_false(grepl(">+-115<", out, fixed = TRUE))
})

test_that("mismatched-line pill renders amber line tag and amber pill class", {
  out <- render_book_pill("BFA", -115, line_quoted = 5.0, is_exact_line = FALSE,
                         is_pick = FALSE)
  expect_match(out, 'class="pill mismatched"', fixed = TRUE)
  expect_match(out, '<span class="line-tag">O5</span>', fixed = TRUE)
  expect_match(out, '>-115<', fixed = TRUE)
})

test_that("mismatched-line under-side renders U<line> tag", {
  out <- render_book_pill("BFA", -105, line_quoted = 5.0, is_exact_line = FALSE,
                         is_pick = FALSE, side = "under")
  expect_match(out, '<span class="line-tag">U5</span>', fixed = TRUE)
})

test_that("NA odds render as muted dashed pill with em-dash", {
  out <- render_book_pill("DK", NA_integer_, NA_real_, is_exact_line = TRUE, FALSE)
  expect_match(out, 'class="pill muted"', fixed = TRUE)
  expect_match(out, '&mdash;', fixed = TRUE)
})

test_that("pick override applies green class on exact line", {
  out <- render_book_pill("H88", 125, 5.5, is_exact_line = TRUE, is_pick = TRUE)
  expect_match(out, 'class="pill pick"', fixed = TRUE)
  expect_match(out, '<span class="book">H88</span>', fixed = TRUE)
})

test_that("pick override applies green class even on mismatched line", {
  out <- render_book_pill("H88", 125, line_quoted = 5.0, is_exact_line = FALSE,
                         is_pick = TRUE)
  # green pick takes precedence visually; class string explicitly notes both
  expect_match(out, 'class="pill pick mismatched"', fixed = TRUE)
  expect_match(out, '<span class="line-tag">O5</span>', fixed = TRUE)
})

test_that("integer line renders without trailing .0 in the tag", {
  out <- render_book_pill("BFA", -115, 5.0, is_exact_line = FALSE, FALSE)
  expect_match(out, '<span class="line-tag">O5</span>', fixed = TRUE)
  expect_false(grepl("O5.0", out))
})

test_that("half-run line renders with .5", {
  out <- render_book_pill("BFA", -115, 5.5, is_exact_line = FALSE, FALSE)
  expect_match(out, '<span class="line-tag">O5.5</span>', fixed = TRUE)
})

test_that("spread mismatch shows signed line value with no O/U prefix", {
  out <- render_book_pill("DK", -110, line_quoted = -1.5, is_exact_line = FALSE,
                          is_pick = FALSE, side = "over", is_totals = FALSE)
  expect_match(out, '<span class="line-tag">-1.5</span>', fixed = TRUE)
  expect_false(grepl("O-1.5", out, fixed = TRUE))
})

test_that("spread mismatch on dog side shows positive signed line", {
  out <- render_book_pill("DK", +120, line_quoted = 1.5, is_exact_line = FALSE,
                          is_pick = FALSE, side = "under", is_totals = FALSE)
  expect_match(out, '<span class="line-tag">+1.5</span>', fixed = TRUE)
  expect_false(grepl("U1.5", out, fixed = TRUE))
})

test_that("spread exact-line does not render any tag (is_totals=FALSE)", {
  out <- render_book_pill("DK", -110, line_quoted = -1.5, is_exact_line = TRUE,
                          is_pick = FALSE, is_totals = FALSE)
  expect_false(grepl("line-tag", out))
})
