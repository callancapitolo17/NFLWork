# Answer Keys/tests/test_devig_pair_matches_tools.R
# Parity guard: book_cell.R::.devig_american_pair must match
# Tools.R::devig_american() for the 2-way case across representative inputs.
#
# Run from "Answer Keys/" directory:
#   Rscript -e 'testthat::test_file("tests/test_devig_pair_matches_tools.R")'

library(testthat)
source("../Tools.R")
source("../MLB Dashboard/book_cell.R")

# Convert a probability to American odds using the same convention as
# .devig_american_pair internally.
.prob_to_amer <- function(p) {
  if (is.na(p)) return(NA_real_)
  if (p >= 0.5) -100 * p / (1 - p)
  else           100 * (1 - p) / p
}

cases <- tibble::tribble(
  ~odd1, ~odd2,
  -110,  -110,
  -140,  +120,
  -200,  +170,
  +150,  -180,
  +105,  -135,
  -300,  +260
)

test_that(".devig_american_pair matches Tools.R::devig_american (2-way) within rounding", {
  for (i in seq_len(nrow(cases))) {
    o1 <- cases$odd1[i]; o2 <- cases$odd2[i]
    pair    <- .devig_american_pair(o1, o2)
    tools_p <- devig_american(o1, o2)            # data.frame(p1, p2)
    tools_amer1 <- .prob_to_amer(tools_p$p1)
    tools_amer2 <- .prob_to_amer(tools_p$p2)
    # .devig_american_pair rounds to the nearest integer; Tools.R returns
    # raw probabilities, so the round-trip can be off by up to 0.5 of an
    # American-odds unit. Use < 1 as a generous-but-tight bound.
    expect_true(abs(pair$fair1 - tools_amer1) < 1,
                label = sprintf("input (%d, %d) side 1: pair=%s tools=%.4f",
                                o1, o2, pair$fair1, tools_amer1))
    expect_true(abs(pair$fair2 - tools_amer2) < 1,
                label = sprintf("input (%d, %d) side 2: pair=%s tools=%.4f",
                                o1, o2, pair$fair2, tools_amer2))
  }
})
