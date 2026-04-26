# Answer Keys/tests/test_books_strip.R
library(testthat)
source("../books_strip.R")

test_that("render_books_strip emits all six pills with values", {
  out <- render_books_strip(0.274, 0.269, 0.271, 0.278, 0.267, 0.272)
  expect_match(out, '<span class="books-strip">', fixed = TRUE)
  expect_match(out, '<span class="pill model">M 27.4</span>', fixed = TRUE)
  expect_match(out, '<span class="pill book">DK 26.9</span>', fixed = TRUE)
  expect_match(out, '<span class="pill book">FD 27.1</span>', fixed = TRUE)
  expect_match(out, '<span class="pill book">PX 27.8</span>', fixed = TRUE)
  expect_match(out, '<span class="pill book">NV 26.7</span>', fixed = TRUE)
  expect_match(out, '<span class="pill cons">Cons 27.2</span>', fixed = TRUE)
})

test_that("render_books_strip dims and dashes NA values", {
  out <- render_books_strip(0.274, NA_real_, 0.271, NA_real_, 0.267, 0.270)
  expect_match(out, '<span class="pill book dim">DK &mdash;</span>', fixed = TRUE)
  expect_match(out, '<span class="pill book dim">PX &mdash;</span>', fixed = TRUE)
  expect_match(out, '<span class="pill model">M 27.4</span>', fixed = TRUE)
  expect_match(out, '<span class="pill book">FD 27.1</span>', fixed = TRUE)
})

test_that("render_books_strip handles all NA except model", {
  out <- render_books_strip(0.30, NA_real_, NA_real_, NA_real_, NA_real_, 0.30)
  expect_match(out, '<span class="pill model">M 30.0</span>', fixed = TRUE)
  expect_match(out, '<span class="pill cons">Cons 30.0</span>', fixed = TRUE)
  expect_match(out, '<span class="pill book dim">DK &mdash;</span>', fixed = TRUE)
  expect_match(out, '<span class="pill book dim">FD &mdash;</span>', fixed = TRUE)
  expect_match(out, '<span class="pill book dim">PX &mdash;</span>', fixed = TRUE)
  expect_match(out, '<span class="pill book dim">NV &mdash;</span>', fixed = TRUE)
})

test_that("render_books_strip rounds to one decimal", {
  out <- render_books_strip(0.27449, 0.26951, 0.27101, 0.27801, 0.26701, 0.27201)
  expect_match(out, '>M 27.4</span>', fixed = TRUE)
  expect_match(out, '>DK 27.0</span>', fixed = TRUE)  # 26.951 → 27.0 (sprintf rounds away from zero at .5)
  expect_match(out, '>PX 27.8</span>', fixed = TRUE)
})

test_that("render_books_strip output is wrapped in container", {
  out <- render_books_strip(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
  expect_true(startsWith(out, '<span class="books-strip">'))
  expect_true(endsWith(out, '</span>'))
})
