# Answer Keys/tests/test_compute_bet_row_id.R
# compute_bet_row_id() must hash the CANONICAL bet identity so the model and
# market paths produce identical ids for the same wager.
# Run from "Answer Keys/":
#   Rscript -e 'testthat::test_file("tests/test_compute_bet_row_id.R")'

library(testthat)
source("../MLB Answer Key/odds_screen.R")

test_that("alternate_ and main collapse to the same id", {
  main <- compute_bet_row_id("g1", "totals", "F5", 8.5, "Over")
  alt  <- compute_bet_row_id("g1", "alternate_totals", "F5", 8.5, "Over")
  expect_equal(main, alt)
})

test_that("line is normalized (8.5 == 8.50) and NA -> empty", {
  expect_equal(compute_bet_row_id("g1", "totals", "F5", 8.5,  "Over"),
               compute_bet_row_id("g1", "totals", "F5", 8.50, "Over"))
  # moneyline (NA line) is stable and distinct from a 0 line
  id_na <- compute_bet_row_id("g1", "h2h", "F5", NA, "NYY")
  expect_type(id_na, "character")
  expect_false(is.na(id_na))
  expect_false(id_na == compute_bet_row_id("g1", "h2h", "F5", 0, "NYY"))
})

test_that("different identity fields produce different ids", {
  base <- compute_bet_row_id("g1", "totals", "F5", 8.5, "Over")
  expect_false(base == compute_bet_row_id("g2", "totals", "F5", 8.5, "Over"))
  expect_false(base == compute_bet_row_id("g1", "spreads", "F5", 8.5, "Over"))
  expect_false(base == compute_bet_row_id("g1", "totals", "FG", 8.5, "Over"))
  expect_false(base == compute_bet_row_id("g1", "totals", "F5", 9.5, "Over"))
  expect_false(base == compute_bet_row_id("g1", "totals", "F5", 8.5, "Under"))
})

test_that("vectorized over inputs", {
  ids <- compute_bet_row_id(c("g1","g2"), c("totals","spreads"),
                            c("F5","F5"), c(8.5, -1.5), c("Over","NYY"))
  expect_length(ids, 2)
  expect_false(ids[1] == ids[2])
})
