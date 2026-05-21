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

# ---------------------------------------------------------------------------
# compute_k_within() — counts voices within +/- band of the median.
# Used by the parlay tab to surface "how many books agree" as a single number.
# ---------------------------------------------------------------------------

test_that("compute_k_within: tight cluster + 1 outlier returns k=4 n=5", {
  out <- compute_k_within(c(0.40, 0.402, 0.403, 0.403, 0.55))
  expect_equal(out$k, 4L)
  expect_equal(out$n, 5L)
})

test_that("compute_k_within: scattered values return k=1 n=5", {
  out <- compute_k_within(c(0.40, 0.50, 0.55, 0.60, 0.70))
  expect_equal(out$k, 1L)
  expect_equal(out$n, 5L)
})

test_that("compute_k_within: full agreement returns k=n=5", {
  out <- compute_k_within(c(0.400, 0.402, 0.403, 0.403, 0.405))
  expect_equal(out$k, 5L)
  expect_equal(out$n, 5L)
})

test_that("compute_k_within: NA entries are excluded from both k and n", {
  out <- compute_k_within(c(0.40, 0.402, 0.403, NA_real_, NA_real_))
  expect_equal(out$k, 3L)
  expect_equal(out$n, 3L)
})

test_that("compute_k_within: fewer than 2 non-NA voices returns NA k", {
  out1 <- compute_k_within(c(0.40, NA_real_, NA_real_, NA_real_, NA_real_))
  expect_true(is.na(out1$k))
  expect_equal(out1$n, 1L)

  out2 <- compute_k_within(c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_))
  expect_true(is.na(out2$k))
  expect_equal(out2$n, 0L)
})

test_that("compute_k_within: custom band widens the consensus window", {
  # median = 0.45; default band 0.01 would give k=2 (the two 0.45s).
  # With band=0.05, the 0.40 and 0.42 also fall in; 0.55 stays out.
  out <- compute_k_within(c(0.40, 0.42, 0.45, 0.45, 0.55), band = 0.05)
  expect_equal(out$k, 4L)
  expect_equal(out$n, 5L)
})
