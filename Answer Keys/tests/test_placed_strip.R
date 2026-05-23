# Answer Keys/tests/test_placed_strip.R
library(testthat)
source("../MLB Dashboard/placed_strip.R")

test_that("renders one chip per placement with account, risk and ticket", {
  chips <- list(
    list(account = "Wagerzon",  risk = 200, ticket = "W1234"),
    list(account = "WagerzonJ", risk = 150, ticket = "W5678")
  )
  html <- render_placed_strip(
    chips,
    all_wz_accounts = c("Wagerzon", "WagerzonJ", "WagerzonC"),
    bet_hash        = "abc",
    book            = "wagerzon"
  )
  # Two chips with expected text
  expect_match(html, '"placement-chip"', fixed = TRUE)
  expect_match(html, '>WZ<', fixed = TRUE)
  expect_match(html, '>WZJ<', fixed = TRUE)
  expect_match(html, '$200', fixed = TRUE)
  expect_match(html, '$150', fixed = TRUE)
  expect_match(html, '#W1234', fixed = TRUE)
  expect_match(html, '#W5678', fixed = TRUE)
})

test_that("+ another appears (enabled) when an untouched WZ account exists", {
  chips <- list(list(account = "Wagerzon", risk = 200, ticket = "W1234"))
  html <- render_placed_strip(
    chips,
    all_wz_accounts = c("Wagerzon", "WagerzonJ"),
    bet_hash        = "abc",
    book            = "wagerzon"
  )
  expect_match(html, "add-another", fixed = TRUE)
  expect_false(grepl("disabled", html, fixed = TRUE))
})

test_that("+ another is disabled once every WZ account has a chip", {
  chips <- list(
    list(account = "Wagerzon",  risk = 200, ticket = "W1"),
    list(account = "WagerzonJ", risk = 150, ticket = "W2")
  )
  html <- render_placed_strip(
    chips,
    all_wz_accounts = c("Wagerzon", "WagerzonJ"),
    bet_hash        = "abc",
    book            = "wagerzon"
  )
  expect_match(html, "add-another", fixed = TRUE)
  expect_match(html, "disabled", fixed = TRUE)
})

test_that("+ another is omitted entirely for non-wagerzon books", {
  chips <- list(list(account = "Hoop88", risk = 100, ticket = "H1"))
  html <- render_placed_strip(
    chips,
    all_wz_accounts = c("Wagerzon", "WagerzonJ"),
    bet_hash        = "abc",
    book            = "hoop88"
  )
  expect_false(grepl("add-another", html, fixed = TRUE))
})

test_that("chips carry data-account, data-risk, data-ticket for JS handlers", {
  chips <- list(list(account = "Wagerzon", risk = 200, ticket = "W1234"))
  html <- render_placed_strip(
    chips, all_wz_accounts = c("Wagerzon"),
    bet_hash = "abc", book = "wagerzon"
  )
  expect_match(html, 'data-account="Wagerzon"', fixed = TRUE)
  expect_match(html, 'data-risk="200"',          fixed = TRUE)
  expect_match(html, 'data-ticket="W1234"',      fixed = TRUE)
  expect_match(html, 'data-bet-hash="abc"',      fixed = TRUE)
})
