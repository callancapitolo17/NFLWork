# Answer Keys/tests/test_parse_legs.R
library(testthat)
source("../parse_legs.R")

test_that("parse_fraction handles integer input", {
  expect_equal(parse_fraction("3"), 3)
})

test_that("parse_fraction handles decimal input", {
  expect_equal(parse_fraction("2.5"), 2.5)
})

test_that("parse_fraction handles ½ unicode half", {
  expect_equal(parse_fraction("2½"), 2.5)
})

test_that("parse_fraction handles ¼ ¾ unicode quarters", {
  expect_equal(parse_fraction("1¼"), 1.25)
  expect_equal(parse_fraction("3¾"), 3.75)
})

test_that("parse_fraction returns NA on unparseable input", {
  expect_true(is.na(parse_fraction("abc")))
  expect_true(is.na(parse_fraction("")))
})

test_that("parse_legs handles triple-play (3 legs)", {
  expect_equal(
    parse_legs("GIANTS TRIPLE-PLAY (SCR 1ST, 1H & GM)"),
    list(
      list(type = "scores_first"),
      list(type = "wins_period", period = "F5"),
      list(type = "wins_period", period = "FG")
    )
  )
})

test_that("parse_legs handles grand slam with unicode fraction", {
  expect_equal(
    parse_legs("GIANTS GRAND-SLAM (SCR 1ST, 1H, GM & SCR U2½)"),
    list(
      list(type = "scores_first"),
      list(type = "wins_period", period = "F5"),
      list(type = "wins_period", period = "FG"),
      list(type = "team_total_under", line = 2.5)
    )
  )
})

test_that("parse_legs handles grand slam with decimal fraction", {
  expect_equal(
    parse_legs("GIANTS GRAND-SLAM (SCR 1ST, 1H, GM & SCR U2.5)"),
    list(
      list(type = "scores_first"),
      list(type = "wins_period", period = "F5"),
      list(type = "wins_period", period = "FG"),
      list(type = "team_total_under", line = 2.5)
    )
  )
})

test_that("parse_legs handles F3 and F7 period tokens", {
  expect_equal(
    parse_legs("GIANTS MEGA (SCR 1ST, F3, 1H, F7 & GM)"),
    list(
      list(type = "scores_first"),
      list(type = "wins_period", period = "F3"),
      list(type = "wins_period", period = "F5"),
      list(type = "wins_period", period = "F7"),
      list(type = "wins_period", period = "FG")
    )
  )
})

test_that("parse_legs handles SCR O<N> team-total over", {
  expect_equal(
    parse_legs("GIANTS BIG (SCR 1ST, GM & SCR O4.5)"),
    list(
      list(type = "scores_first"),
      list(type = "wins_period", period = "FG"),
      list(type = "team_total_over", line = 4.5)
    )
  )
})

test_that("parse_legs returns NULL and warns on unknown token", {
  expect_warning(
    result <- parse_legs("GIANTS MYSTERY (SCR 1ST, XYZ & GM)"),
    "Unknown leg token"
  )
  expect_null(result)
})

test_that("parse_legs returns NULL on description with no parenthetical", {
  expect_null(parse_legs("GIANTS TRIPLE-PLAY"))
})

test_that("eval_leg scores_first home reads home_scored_in_1st", {
  samples <- data.frame(
    home_scored_in_1st = c(1L, 0L, 1L, NA),
    away_scored_in_1st = c(1L, 1L, 0L, 0L)   # noise — should be ignored
  )
  result <- eval_leg(list(type = "scores_first"),
                    samples, side = "home",
                    team_runs = NULL, opp_runs = NULL)
  expect_equal(result, c(TRUE, FALSE, TRUE, NA))
})

test_that("eval_leg scores_first away reads away_scored_in_1st", {
  samples <- data.frame(
    home_scored_in_1st = c(1L, 0L, 1L, NA),   # noise — should be ignored
    away_scored_in_1st = c(0L, 1L, 1L, 1L)
  )
  result <- eval_leg(list(type = "scores_first"),
                    samples, side = "away",
                    team_runs = NULL, opp_runs = NULL)
  expect_equal(result, c(FALSE, TRUE, TRUE, TRUE))
})

test_that("eval_leg wins_period F5 home checks home_margin_f5 > 0", {
  samples <- data.frame(home_margin_f5 = c(2, 0, -1, 3))
  result <- eval_leg(list(type = "wins_period", period = "F5"),
                     samples, side = "home",
                     team_runs = c(NA, NA, NA, NA), opp_runs = c(NA, NA, NA, NA))
  expect_equal(result, c(TRUE, FALSE, FALSE, TRUE))
})

test_that("eval_leg wins_period F5 away checks home_margin_f5 < 0 (strict)", {
  samples <- data.frame(home_margin_f5 = c(2, 0, -1, -3))
  result <- eval_leg(list(type = "wins_period", period = "F5"),
                     samples, side = "away",
                     team_runs = c(NA, NA, NA, NA), opp_runs = c(NA, NA, NA, NA))
  expect_equal(result, c(FALSE, FALSE, TRUE, TRUE))
})

test_that("eval_leg wins_period FG uses home_margin column", {
  samples <- data.frame(home_margin = c(1, -1, 0))
  result <- eval_leg(list(type = "wins_period", period = "FG"),
                     samples, side = "home",
                     team_runs = c(NA, NA, NA), opp_runs = c(NA, NA, NA))
  expect_equal(result, c(TRUE, FALSE, FALSE))
})

test_that("eval_leg wins_period F3 uses home_margin_f3 column", {
  samples <- data.frame(home_margin_f3 = c(2, -1, 0))
  result <- eval_leg(list(type = "wins_period", period = "F3"),
                     samples, side = "home",
                     team_runs = c(NA, NA, NA), opp_runs = c(NA, NA, NA))
  expect_equal(result, c(TRUE, FALSE, FALSE))
})

test_that("eval_leg wins_period F7 uses home_margin_f7 column", {
  samples <- data.frame(home_margin_f7 = c(2, -1, 0))
  result <- eval_leg(list(type = "wins_period", period = "F7"),
                     samples, side = "home",
                     team_runs = c(NA, NA, NA), opp_runs = c(NA, NA, NA))
  expect_equal(result, c(TRUE, FALSE, FALSE))
})

test_that("eval_leg team_total_under uses team_runs (side-dependent)", {
  samples <- data.frame(x = c(1, 2, 3))  # unused
  result <- eval_leg(list(type = "team_total_under", line = 2.5),
                     samples, side = "home",
                     team_runs = c(1, 2, 3), opp_runs = c(5, 5, 5))
  expect_equal(result, c(TRUE, TRUE, FALSE))
})

test_that("eval_leg team_total_over uses team_runs", {
  samples <- data.frame(x = c(1, 2, 3))
  result <- eval_leg(list(type = "team_total_over", line = 2.5),
                     samples, side = "home",
                     team_runs = c(1, 2, 3), opp_runs = c(5, 5, 5))
  expect_equal(result, c(FALSE, FALSE, TRUE))
})

test_that("eval_leg unknown type raises error", {
  samples <- data.frame(x = 1)
  expect_error(
    eval_leg(list(type = "nonexistent"),
             samples, side = "home",
             team_runs = c(1), opp_runs = c(1)),
    "Unknown leg type"
  )
})

test_that("compute_prop_fair triple-play home returns joint P(all 3 legs)", {
  # 10 synthetic rows. home_triple = scored_first==1 AND margin_f5>0 AND home_margin>0
  # Expected hits (same fixture as test_triple_play's home-side test): 5/10 = 0.5
  samples <- data.frame(
    home_margin        = c( 2,  3, -1,  4,  1,  5, -2,  2,  6, -3),
    home_margin_f5     = c( 1,  2,  0,  3,  0,  2, -1,  1,  4, -2),
    home_scored_in_1st = c( 1L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 1L, 0L),
    away_scored_in_1st = c( 0L, 0L, 1L, 0L, 0L, 0L, 1L, 1L, 0L, 1L),
    total_final_score  = rep(10, 10)  # irrelevant for triple-play
  )
  legs <- list(
    list(type = "scores_first"),
    list(type = "wins_period", period = "F5"),
    list(type = "wins_period", period = "FG")
  )
  expect_equal(compute_prop_fair(samples, "home", legs), 0.5)
})

test_that("compute_prop_fair grand slam (team_total_under) home", {
  # 4-leg: scores_first AND wins_f5 AND wins_game AND team_under_2.5
  # home_runs = (total + margin) / 2
  samples <- data.frame(
    home_margin        = c( 2,  3,  1,  2,  2),
    home_margin_f5     = c( 1,  2,  1,  1,  1),
    home_scored_in_1st = c( 1L, 1L, 1L, 1L, 1L),
    away_scored_in_1st = c( 0L, 0L, 0L, 0L, 0L),
    total_final_score  = c( 4,  6,  2,  4,  4)
  )
  # home_runs: (4+2)/2=3, (6+3)/2=4.5, (2+1)/2=1.5, (4+2)/2=3, (4+2)/2=3
  # team_under_2.5: FALSE, FALSE, TRUE, FALSE, FALSE → 1/5 hits
  # All other legs pass on all rows
  legs <- list(
    list(type = "scores_first"),
    list(type = "wins_period", period = "F5"),
    list(type = "wins_period", period = "FG"),
    list(type = "team_total_under", line = 2.5)
  )
  expect_equal(compute_prop_fair(samples, "home", legs), 0.2)
})

test_that("compute_prop_fair returns NA_real_ on NULL legs", {
  samples <- data.frame(
    home_margin = 1, home_margin_f5 = 1,
    home_scored_in_1st = 1L, away_scored_in_1st = 0L,
    total_final_score = 2
  )
  expect_true(is.na(compute_prop_fair(samples, "home", NULL)))
})

test_that("compute_prop_fair returns NA_real_ on empty legs list", {
  samples <- data.frame(
    home_margin = 1, home_margin_f5 = 1,
    home_scored_in_1st = 1L, away_scored_in_1st = 0L,
    total_final_score = 2
  )
  expect_true(is.na(compute_prop_fair(samples, "home", list())))
})

test_that("compute_prop_fair returns NA_real_ on empty samples", {
  samples <- data.frame(
    home_margin = integer(0), home_margin_f5 = integer(0),
    home_scored_in_1st = integer(0), away_scored_in_1st = integer(0),
    total_final_score = integer(0)
  )
  legs <- list(list(type = "scores_first"))
  expect_true(is.na(compute_prop_fair(samples, "home", legs)))
})

test_that("compute_prop_fair excludes rows where the side-specific scored-in-1st column is NA", {
  # NA filter is now per-side, applied after AND-reduce.
  # The home side and away side independently drop different rows based on
  # which inning-1 column is NA — they do NOT share a global NA filter.
  #
  # Fixture: 4 rows with margins set so home wins F5 + FG everywhere
  #   (home_margin > 0 and home_margin_f5 > 0). For the AWAY-side test we
  #   use a separate fixture with mirrored margins so away wins F5 + FG.
  # That way the only leg that gates hits is scores_first.
  #
  # HOME-side fixture: row 4 has home_scored_in_1st = NA → drop row 4.
  #   Remaining: c(1, 0, 0) → hits c(T, F, F) → 1/3.
  home_samples <- data.frame(
    home_margin        = c(1, 1, 1, 1),
    home_margin_f5     = c(1, 1, 1, 1),
    home_scored_in_1st = c(1L, 0L, 0L, NA),
    away_scored_in_1st = c(0L, 1L, 1L, 1L),  # noise for HOME side
    total_final_score  = c(2, 2, 2, 2)
  )
  legs <- list(
    list(type = "scores_first"),
    list(type = "wins_period", period = "F5"),
    list(type = "wins_period", period = "FG")
  )
  expect_equal(compute_prop_fair(home_samples, "home", legs), 1/3)

  # AWAY-side fixture: away wins F5 + FG everywhere; row 1 has
  # away_scored_in_1st = NA → drop row 1.
  #   Remaining away_scored_in_1st: c(1, 1, 0) → hits c(T, T, F) → 2/3.
  away_samples <- data.frame(
    home_margin        = c(-1, -1, -1, -1),
    home_margin_f5     = c(-1, -1, -1, -1),
    home_scored_in_1st = c( 1L,  0L,  0L,  0L),  # noise for AWAY side
    away_scored_in_1st = c(NA,  1L,  1L,  0L),
    total_final_score  = c( 2,   2,   2,   2)
  )
  expect_equal(compute_prop_fair(away_samples, "away", legs), 2/3)
})

test_that("compute_prop_fair away side with team_total_under uses away_runs", {
  # away_runs = (total - margin) / 2
  # For side=away to win: home_margin_f5 < 0, home_margin < 0, away_scored_in_1st == 1
  samples <- data.frame(
    home_margin        = c(-2, -3, -1),
    home_margin_f5     = c(-1, -2, -1),
    home_scored_in_1st = c( 0L, 0L, 0L),
    away_scored_in_1st = c( 1L, 1L, 1L),
    total_final_score  = c( 4,  6,  2)
  )
  # away_runs: (4-(-2))/2=3, (6-(-3))/2=4.5, (2-(-1))/2=1.5
  # under 2.5: FALSE, FALSE, TRUE → 1/3
  legs <- list(
    list(type = "scores_first"),
    list(type = "wins_period", period = "F5"),
    list(type = "wins_period", period = "FG"),
    list(type = "team_total_under", line = 2.5)
  )
  expect_equal(compute_prop_fair(samples, "away", legs), 1/3)
})
