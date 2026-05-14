# MLB Bets Tab — PR B (display features) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add two display features to the MLB Dashboard bets tab: a manual Kelly Calculator widget below the bankroll/Kelly settings, and a per-cell devigged-odds toggle (RAW / FAIR) on each bet card's price grid that defaults to FAIR and renders devigged values in American odds via probit (z-shift) math.

**Architecture:** All changes are server-render + JS. No new Python, no new server endpoints, no schema changes. The devig math runs in R during card render using a small inline helper (`.devig_american_pair`) that mirrors `Tools.R::devig_american` for the 2-way closed-form case — chosen over importing Tools.R wholesale (heavy, side-effect-laden) and over a JS reimplementation (risk of drift). Each cell emits *both* a raw and a fair `<span>`; CSS hides one or the other based on a class on the grid container, so the toggle is a one-line JS class swap. Kelly calculator reuses the existing `calculateKellyBet` JS function at `mlb_dashboard.R:3746` and reads the dashboard's existing Bankroll + Kelly Fraction inputs live.

**Tech Stack:** R 4.x, testthat (existing); plain DOM JavaScript (existing pattern).

**Worktree:** Before Task 1, create a fresh worktree off `main` named `mlb-bets-tab-pr-b-display` (e.g. via `EnterWorktree({name: "mlb-bets-tab-pr-b-display"})` if you have the Superpowers worktree skill, or `git worktree add .claude/worktrees/mlb-bets-tab-pr-b-display -b worktree-mlb-bets-tab-pr-b-display main` otherwise). All paths in this plan assume that worktree as the working directory. **PR A must be merged to `main` first** so the `_fg` suffix canonicalization and the spread-opposite line-tag fix are present — otherwise the FAIR-view screenshots in Step 9 won't match expectations.

---

## File Structure

```
Answer Keys/
├── MLB Dashboard/
│   ├── mlb_dashboard.R                       (modify: 4 sections)
│   └── book_cell.R                           (modify: render_book_cell signature)
├── tests/
│   ├── test_book_cell.R                      (CREATE: probit-pair + cell HTML tests)
│   └── test_devig_pair_matches_tools.R       (CREATE: parity test vs Tools.R)
├── CLAUDE.md                                 (modify: note Kelly + devig toggle)
└── MLB Dashboard/README.md                   (modify: feature list)
```

Each file's responsibility:

- **`book_cell.R`** — adds a private `.devig_american_pair(odd1, odd2)` helper (2-way probit closed form, ~15 lines). `render_book_cell()` gains an `opposite_american_odds` parameter. When both odds are non-NA, the function computes the devigged pair and emits both `<span class="raw">…</span>` and `<span class="fair">…</span>` inside the cell. CSS in `mlb_dashboard.R` toggles which span is visible.
- **`mlb_dashboard.R`** — four sections change:
  1. `render_price_grid_row()` gains `other_side_wide_row` parameter; for each book it looks up the *other* side's `<book>_american_odds` and passes it to `render_book_cell()`.
  2. `render_bet_card()` calls `render_price_grid_row()` with both wide rows so each row has access to its sibling.
  3. The `grid_html` template in `render_bet_card()` gets a new toggle-buttons block (`<div class="grid-toggle">RAW | FAIR</div>`) and the `<div class="price-grid">` gets a default `show-fair` class.
  4. The bets-tab page template gains: (a) the Kelly Calculator widget HTML block right after the existing `.sizing-controls` strip, and (b) a small `<script>` block with the toggle click handler and Kelly calc input handlers.
- **`tests/test_book_cell.R`** — new test file. Three sets of tests: `.devig_american_pair` math correctness; `render_book_cell` HTML output (both spans present, correct values); the legacy `render_book_pill` shim still works (no opposite_american_odds, falls back to old behavior).
- **`tests/test_devig_pair_matches_tools.R`** — new test file. Asserts that `book_cell.R`'s `.devig_american_pair()` output matches `Tools.R::devig_american()` output across a table of representative inputs. Catches drift if Tools.R changes or if our reimplementation has a math error.
- **`Answer Keys/CLAUDE.md`** — add bullets describing the Kelly Calculator widget placement and the per-cell devig toggle (default FAIR, probit math, mirrors Tools.R).
- **`MLB Dashboard/README.md`** — add Kelly Calculator + per-cell devig toggle to the feature list.

---

## Task 1 — `.devig_american_pair` helper in book_cell.R

**Why first.** The math primitive everything else depends on. TDD-friendly with assert-against-Tools.R parity.

**Files:**
- Modify: `Answer Keys/MLB Dashboard/book_cell.R` (add helper near top)
- Create: `Answer Keys/tests/test_book_cell.R`
- Create: `Answer Keys/tests/test_devig_pair_matches_tools.R`

- [ ] **Step 1: Write the failing tests for the helper math**

Create `Answer Keys/tests/test_book_cell.R`:

```r
# Answer Keys/tests/test_book_cell.R
# Tests for book_cell.R: .devig_american_pair (probit math) and
# render_book_cell (HTML output with raw + fair spans).
#
# Run from "Answer Keys/" directory:
#   Rscript -e 'testthat::test_file("tests/test_book_cell.R")'

library(testthat)
source("MLB Dashboard/book_cell.R")

# Helper: convert American odds to implied probability for the assertion math.
.amer_to_prob <- function(o) ifelse(o > 0, 100 / (o + 100), -o / (-o + 100))

test_that(".devig_american_pair returns devigged American odds (no-vig coin flip)", {
  # -110 / -110 → both implied 52.38%, devigged → 50% / 50% → +100 / +100
  out <- .devig_american_pair(-110, -110)
  expect_equal(out$fair1, +100, tolerance = 0.5)
  expect_equal(out$fair2, +100, tolerance = 0.5)
})

test_that(".devig_american_pair returns devigged American odds (asymmetric)", {
  # -140 / +120 → implied 58.3% / 45.5%, sum 103.8% (vig)
  # Probit-devigged → ~56.4% / ~43.6% → about -129 / +129
  out <- .devig_american_pair(-140, +120)
  expect_lt(abs(out$fair1 - (-129)), 1)
  expect_lt(abs(out$fair2 - (+129)), 1)
})

test_that(".devig_american_pair returns NA on bad inputs", {
  expect_true(is.na(.devig_american_pair(NA, +120)$fair1))
  expect_true(is.na(.devig_american_pair(0,  +120)$fair1))
  expect_true(is.na(.devig_american_pair(-110, NA)$fair2))
})
```

Create `Answer Keys/tests/test_devig_pair_matches_tools.R`:

```r
# Answer Keys/tests/test_devig_pair_matches_tools.R
# Parity guard: book_cell.R::.devig_american_pair must match
# Tools.R::devig_american() for the 2-way case across representative inputs.
#
# Run from "Answer Keys/" directory:
#   Rscript -e 'testthat::test_file("tests/test_devig_pair_matches_tools.R")'

library(testthat)
source("Tools.R")
source("MLB Dashboard/book_cell.R")

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
    expect_lt(abs(pair$fair1 - tools_amer1), 1,
              info = sprintf("input (%d, %d) side 1 mismatch", o1, o2))
    expect_lt(abs(pair$fair2 - tools_amer2), 1,
              info = sprintf("input (%d, %d) side 2 mismatch", o1, o2))
  }
})
```

- [ ] **Step 2: Run the tests; confirm they fail (helper not defined yet)**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-pr-b-display/Answer Keys" && \
  Rscript -e 'testthat::test_file("tests/test_book_cell.R")' && \
  Rscript -e 'testthat::test_file("tests/test_devig_pair_matches_tools.R")'
```

Expected: errors like `could not find function ".devig_american_pair"` for both files.

- [ ] **Step 3: Add the helper to book_cell.R**

Open `Answer Keys/MLB Dashboard/book_cell.R`. Add this helper near the top of the file, immediately after `.format_line_value` (around line 23):

```r
#' 2-way probit devig: given two American odds, return the no-vig
#' (fair) American odds for both sides.
#'
#' Mirrors `Tools.R::devig_american()` for the 2-way case using the
#' closed-form solution `c* = -(z1 + z2) / 2`. Inlined here to avoid
#' sourcing all of Tools.R into the dashboard. A parity test in
#' `Answer Keys/tests/test_devig_pair_matches_tools.R` guards against
#' drift across a table of representative inputs.
#'
#' @param odd1 American odds for side 1 (integer or numeric).
#' @param odd2 American odds for side 2.
#' @return list(fair1 = numeric, fair2 = numeric). Returns NA pair when
#'   either input is NA or zero.
.devig_american_pair <- function(odd1, odd2) {
  if (is.na(odd1) || is.na(odd2) || odd1 == 0 || odd2 == 0) {
    return(list(fair1 = NA_real_, fair2 = NA_real_))
  }
  # Implied probabilities from American odds
  p1 <- if (odd1 > 0) 100 / (odd1 + 100) else -odd1 / (-odd1 + 100)
  p2 <- if (odd2 > 0) 100 / (odd2 + 100) else -odd2 / (-odd2 + 100)
  # Probit z-shift: z' = z + c, with c chosen so p1' + p2' = 1
  eps <- 1e-9
  z1 <- qnorm(min(max(p1, eps), 1 - eps))
  z2 <- qnorm(min(max(p2, eps), 1 - eps))
  c_star <- -(z1 + z2) / 2
  q1 <- pnorm(z1 + c_star)
  q2 <- pnorm(z2 + c_star)
  # Convert devigged probabilities back to American odds
  to_amer <- function(p) {
    if (p >= 0.5) round(-100 * p / (1 - p))
    else          round( 100 * (1 - p) / p)
  }
  list(fair1 = to_amer(q1), fair2 = to_amer(q2))
}
```

- [ ] **Step 4: Run the tests; confirm green**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-pr-b-display/Answer Keys" && \
  Rscript -e 'testthat::test_file("tests/test_book_cell.R")' && \
  Rscript -e 'testthat::test_file("tests/test_devig_pair_matches_tools.R")'
```

Expected: every test passes (the third test in test_book_cell.R about NA inputs, the asymmetric +/-129 test, the parity assertions).

- [ ] **Step 5: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-pr-b-display" && \
  git add "Answer Keys/MLB Dashboard/book_cell.R" \
          "Answer Keys/tests/test_book_cell.R" \
          "Answer Keys/tests/test_devig_pair_matches_tools.R" && \
  git commit -m "$(cat <<'EOF'
feat(bets-tab): add .devig_american_pair helper to book_cell.R

Inlines a 2-way probit devig (closed-form c* = -(z1+z2)/2) used by
the upcoming per-cell devig toggle. Mirrors Tools.R::devig_american
exactly for the 2-way case; full Tools.R is too heavy to source into
the dashboard. Parity test enforces equivalence across representative
inputs so any future drift fails CI.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2 — Update `render_book_cell` to emit raw + fair spans

**Files:**
- Modify: `Answer Keys/MLB Dashboard/book_cell.R:38-68` (`render_book_cell` body)
- Modify: `Answer Keys/tests/test_book_cell.R` (add HTML output tests)

- [ ] **Step 1: Append HTML output tests to `test_book_cell.R`**

Add to the bottom of `Answer Keys/tests/test_book_cell.R`:

```r
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
  expect_match(html, '<span class="fair">\\+129</span>', fixed = FALSE)
  # No alt-line tag because is_exact_line = TRUE
  expect_false(grepl('alt-line', html))
})

test_that("render_book_cell without opposite_american_odds emits raw only (legacy)", {
  html <- render_book_cell(
    american_odds = +120,
    line_quoted   = -0.5,
    is_exact_line = TRUE,
    is_pick       = FALSE,
    side_word     = "over",
    is_totals     = FALSE
  )
  # Raw span present; fair span omitted (or empty) when opposite missing
  expect_match(html, '<span class="raw">\\+120</span>', fixed = FALSE)
  expect_false(grepl('<span class="fair">\\+', html))
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
```

- [ ] **Step 2: Run the new tests; confirm they fail**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-pr-b-display/Answer Keys" && \
  Rscript -e 'testthat::test_file("tests/test_book_cell.R")'
```

Expected: the three new tests fail because `render_book_cell` doesn't accept `opposite_american_odds` and doesn't emit `<span class="raw">` / `<span class="fair">`.

- [ ] **Step 3: Update `render_book_cell` signature and body**

Open `Answer Keys/MLB Dashboard/book_cell.R`. Replace the `render_book_cell` function (lines 38-68) with:

```r
#' Render one bets-tab grid cell.
#'
#' @param american_odds Integer odds (e.g., 125, -110). NA -> empty state.
#' @param opposite_american_odds Integer odds for the OTHER side at the same
#'   book (e.g., the Under price when this cell is the Over). Used to compute
#'   probit-devigged fair odds for the FAIR view of the toggle. NA -> no
#'   fair span emitted (cell shows raw only, behaves like legacy).
#' @param line_quoted Numeric line the book is showing on this side.
#' @param is_exact_line Boolean: TRUE when book's line matches the model line
#'   exactly; FALSE -> alt state.
#' @param is_pick TRUE if this is the picked book on the pick side; overrides
#'   exact/alt with the green pick state.
#' @param side_word "over" or "under" — only used for the O/U prefix on a
#'   mismatched totals line tag.
#' @param is_totals TRUE for totals markets (line tag gets O/U prefix);
#'   FALSE for spreads (signed line value, e.g. "-1.5").
#' @return HTML string for the cell. Contains both <span class="raw"> and
#'   (when devig is computable) <span class="fair">; CSS on the parent
#'   .price-grid container determines which is visible (toggle).
render_book_cell <- function(american_odds, line_quoted, is_exact_line,
                              is_pick = FALSE, side_word = "over",
                              is_totals = TRUE,
                              opposite_american_odds = NA_integer_) {
  # State 1: empty (no quote)
  if (is.na(american_odds)) {
    return('<div class="cell empty"><span class="raw">&mdash;</span><span class="fair">&mdash;</span></div>')
  }

  raw_str <- if (american_odds > 0) paste0("+", american_odds) else as.character(american_odds)

  # Compute devigged American for the FAIR span (if we have both sides).
  fair_html <- ""
  if (!is.na(opposite_american_odds)) {
    pair <- .devig_american_pair(american_odds, opposite_american_odds)
    if (!is.na(pair$fair1)) {
      fair_str <- if (pair$fair1 > 0) paste0("+", as.integer(pair$fair1))
                  else as.character(as.integer(pair$fair1))
      fair_html <- sprintf('<span class="fair">%s</span>', fair_str)
    }
  }

  is_mismatched <- !isTRUE(is_exact_line)

  cell_class <- if (is_pick) "cell pick"
                else if (is_mismatched) "cell alt"
                else "cell exact"

  tag_html <- ""
  if (is_mismatched && !is_pick) {
    if (is_totals) {
      prefix <- if (side_word == "under") "U" else "O"
      tag_html <- sprintf('<span class="alt-line">%s%s</span>',
                          prefix, .format_line_value(line_quoted))
    } else {
      tag_html <- sprintf('<span class="alt-line">%s</span>',
                          .format_line_value(line_quoted, signed = TRUE))
    }
  }

  sprintf('<div class="%s">%s<span class="raw">%s</span>%s</div>',
          cell_class, tag_html, raw_str, fair_html)
}
```

- [ ] **Step 4: Run the test suite; confirm green**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-pr-b-display/Answer Keys" && \
  Rscript -e 'testthat::test_file("tests/test_book_cell.R")' && \
  Rscript -e 'testthat::test_file("tests/test_devig_pair_matches_tools.R")' && \
  Rscript -e 'testthat::test_file("tests/test_odds_screen.R")'
```

Expected: all green. The legacy `render_book_pill` shim still calls `render_book_cell` without `opposite_american_odds`, so it falls into the "no fair span" branch and behaves exactly as before.

- [ ] **Step 5: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-pr-b-display" && \
  git add "Answer Keys/MLB Dashboard/book_cell.R" "Answer Keys/tests/test_book_cell.R" && \
  git commit -m "$(cat <<'EOF'
feat(bets-tab): render_book_cell emits raw + fair spans for devig toggle

Adds opposite_american_odds parameter (default NA, backward compatible
with legacy render_book_pill shim). When both odds are present, computes
the probit-devigged fair American odds via .devig_american_pair and
emits both <span class="raw"> and <span class="fair">. CSS on the
.price-grid container will hide one or the other based on the
toggle state. Legacy single-odds callers fall through to raw-only
output (no behavior change).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3 — `render_price_grid_row` passes opposite-side odds

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R:1161-1185` (`render_price_grid_row`)

- [ ] **Step 1: Update `render_price_grid_row` to accept and pass `other_side_wide_row`**

Open `Answer Keys/MLB Dashboard/mlb_dashboard.R`. Replace `render_price_grid_row` (lines 1161-1185) with:

```r
#' Render one row of the price grid (1 row label + 8 book cells).
#'
#' @param wide_row Single-row tibble from book_prices_wide (THIS side's data),
#'   or NULL if no row exists for this side.
#' @param other_side_wide_row Single-row tibble from book_prices_wide (the
#'   OTHER side's data — pass oppside when rendering the pick row, and vice
#'   versa). NULL is allowed; the cells will fall back to raw-only display.
#' @param side_label Text for the leftmost cell, e.g. "BOS -2.5" or "Over 5.5"
#' @param is_pick_side TRUE for the pick row; FALSE for the opposite row
#' @param pick_book Bookmaker key (e.g. "wagerzon") of the picked book
#' @param side_word "over" or "under" — drives O/U prefix on alt-line tag
#' @param is_totals Boolean — drives spread vs totals line-tag formatting
#' @return HTML string: row-label cell + 8 grid cells.
render_price_grid_row <- function(wide_row, side_label, is_pick_side,
                                   pick_book, side_word, is_totals,
                                   other_side_wide_row = NULL) {
  cells <- vapply(BOOK_ORDER_V8, function(b) {
    odds_col  <- paste0(b, "_american_odds")
    lq_col    <- paste0(b, "_line_quoted")
    exact_col <- paste0(b, "_is_exact_line")
    odds  <- if (!is.null(wide_row) && odds_col  %in% names(wide_row)) wide_row[[odds_col]]  else NA_integer_
    lq    <- if (!is.null(wide_row) && lq_col    %in% names(wide_row)) wide_row[[lq_col]]    else NA_real_
    exact <- if (!is.null(wide_row) && exact_col %in% names(wide_row)) wide_row[[exact_col]] else NA
    # Other-side odds for THIS book — used for devig FAIR view.
    opp_odds <- if (!is.null(other_side_wide_row) &&
                    odds_col %in% names(other_side_wide_row)) {
      other_side_wide_row[[odds_col]]
    } else NA_integer_
    render_book_cell(
      american_odds          = if (is.na(odds))     NA_integer_ else as.integer(odds),
      opposite_american_odds = if (is.na(opp_odds)) NA_integer_ else as.integer(opp_odds),
      line_quoted            = lq,
      is_exact_line          = exact,
      is_pick                = is_pick_side && (b == pick_book),
      side_word              = side_word,
      is_totals              = is_totals
    )
  }, character(1))

  paste0(
    sprintf('<div class="row-hdr">%s</div>',
            htmltools::htmlEscape(side_label)),
    paste(cells, collapse = "")
  )
}
```

- [ ] **Step 2: Spot-check the function still loads cleanly**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-pr-b-display/Answer Keys" && \
  Rscript -e 'source("MLB Dashboard/book_cell.R"); source("MLB Dashboard/mlb_dashboard.R", echo = FALSE)' 2>&1 | head -20
```

Expected: any sourcing errors that appear here would also appear when loading the dashboard. (mlb_dashboard.R has top-level side effects so this might print other status — that's fine; you're looking for syntax/parse errors specific to the new function.)

- [ ] **Step 3: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-pr-b-display" && \
  git add "Answer Keys/MLB Dashboard/mlb_dashboard.R" && \
  git commit -m "$(cat <<'EOF'
feat(bets-tab): render_price_grid_row threads opposite-side odds to cells

Adds other_side_wide_row parameter so each cell can compute devig
against the same book's other-side price. Backward compatible (default
NULL falls back to raw-only display). Caller wiring in render_bet_card
follows in the next commit.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4 — `render_bet_card` passes both wide rows

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R:1287-1300` (the two `render_price_grid_row` calls inside `render_bet_card`)

- [ ] **Step 1: Update both calls to pass the other side's row**

Open `Answer Keys/MLB Dashboard/mlb_dashboard.R`. Locate the `render_bet_card` body around lines 1287-1300:

```r
  pick_row <- render_price_grid_row(pickside_wide_row, pickside_label,
                                     is_pick_side = TRUE,
                                     pick_book = pick_book,
                                     side_word = side_word,
                                     is_totals = is_totals)

  opp_row <- if (!is.null(oppside_wide_row) && !is.na(oppside_label) &&
                 nchar(oppside_label) > 0) {
    render_price_grid_row(oppside_wide_row, oppside_label,
                           is_pick_side = FALSE,
                           pick_book = pick_book,
                           side_word = if (side_word == "over") "under" else "over",
                           is_totals = is_totals)
  } else ""
```

Replace with:

```r
  pick_row <- render_price_grid_row(pickside_wide_row, pickside_label,
                                     is_pick_side = TRUE,
                                     pick_book = pick_book,
                                     side_word = side_word,
                                     is_totals = is_totals,
                                     other_side_wide_row = oppside_wide_row)

  opp_row <- if (!is.null(oppside_wide_row) && !is.na(oppside_label) &&
                 nchar(oppside_label) > 0) {
    render_price_grid_row(oppside_wide_row, oppside_label,
                           is_pick_side = FALSE,
                           pick_book = pick_book,
                           side_word = if (side_word == "over") "under" else "over",
                           is_totals = is_totals,
                           other_side_wide_row = pickside_wide_row)
  } else ""
```

- [ ] **Step 2: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-pr-b-display" && \
  git add "Answer Keys/MLB Dashboard/mlb_dashboard.R" && \
  git commit -m "$(cat <<'EOF'
feat(bets-tab): pass both wide rows to render_price_grid_row for devig

Wires render_bet_card so each row's grid cells have access to the
other side's per-book quote. Pick row receives oppside_wide_row;
opposite row receives pickside_wide_row. Each cell can now emit its
devigged FAIR span alongside the RAW span.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5 — Toggle markup + CSS in the V8 card grid

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` — `render_bet_card`'s `grid_html` template (around line 1302), plus the CSS block in the bets-tab page template.

- [ ] **Step 1: Add toggle buttons + default `show-fair` class to the grid div**

Open `Answer Keys/MLB Dashboard/mlb_dashboard.R`. Locate line 1302:

```r
  grid_html <- paste0('<div class="price-grid">', col_hdrs, pick_row, opp_row, '</div>')
```

Replace with:

```r
  # Toggle buttons sit just above the grid; clicking "RAW" or "FAIR"
  # flips the `show-raw` / `show-fair` class on the .price-grid below it
  # (handled by the page-level <script> initialized in setupDevigToggles).
  # Default is FAIR — matches user's manual no-vig workflow. Class is
  # `show-fair` on every render; per-card state, no global preference.
  toggle_html <- paste0(
    '<div class="grid-toggle">',
      '<button type="button" data-view="raw"  onclick="toggleDevigView(this, \'raw\')">RAW</button>',
      '<button type="button" data-view="fair" class="on" onclick="toggleDevigView(this, \'fair\')">FAIR</button>',
    '</div>'
  )
  grid_html <- paste0(toggle_html,
                       '<div class="price-grid show-fair">',
                       col_hdrs, pick_row, opp_row, '</div>')
```

- [ ] **Step 2: Add the CSS rules for the toggle and the raw/fair span visibility**

In the same file (`mlb_dashboard.R`), find the CSS block for the V8 card layout. Search for the existing class definitions like `.price-grid` — they live inside an inline `<style>` block in the bets-tab page template. Add these rules anywhere in that block (group with existing `.price-grid` rules for readability):

```css
/* Per-card toggle for raw vs devigged American odds. Default: show fair. */
.grid-toggle {
  display: inline-flex;
  background: #161b22;
  border: 1px solid #2a3442;
  border-radius: 7px;
  overflow: hidden;
  font-size: 11px;
  letter-spacing: 0.04em;
  font-family: 'SF Mono', monospace;
  margin: 8px 0 6px auto;
}
.grid-toggle button {
  padding: 5px 12px;
  background: transparent;
  border: none;
  color: #7d8590;
  cursor: pointer;
  font-family: inherit;
}
.grid-toggle button.on {
  background: rgba(63, 185, 80, 0.15);
  color: #3fb950;
}
.price-grid .cell .raw,
.price-grid .cell .fair { display: none; }
.price-grid.show-raw  .cell .raw  { display: inline; }
.price-grid.show-fair .cell .fair { display: inline; }
/* Empty cells render a dash for both spans; show whichever the parent class wants. */
.price-grid .cell.empty .raw,
.price-grid .cell.empty .fair { color: #7d8590; }
```

> If you can't find the existing `<style>` block, grep:
> ```bash
> grep -n '.price-grid\|<style>\|</style>' "Answer Keys/MLB Dashboard/mlb_dashboard.R" | head
> ```
> Add the new rules inside the same `<style>...</style>` block where `.price-grid` is already defined.

- [ ] **Step 3: Commit (markup + CSS, JS handler comes in Task 6)**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-pr-b-display" && \
  git add "Answer Keys/MLB Dashboard/mlb_dashboard.R" && \
  git commit -m "$(cat <<'EOF'
feat(bets-tab): per-card RAW/FAIR toggle markup + CSS

Adds the toggle buttons and the .price-grid CSS rules that hide raw or
fair spans based on the .show-raw / .show-fair class on the grid.
Default rendered class is .show-fair so devigged American odds are
visible on first paint. JS click handler comes in the next commit.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 6 — JS click handler for the toggle

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` — add `<script>` block on the bets-tab page

- [ ] **Step 1: Add the click-handler JS to the bets-tab page template**

Find the existing inline `<script>` blocks at the bottom of the bets-tab page template (search for `function placeBet` or `DOMContentLoaded` in `mlb_dashboard.R`). Append this near the other DOMContentLoaded handlers:

```html
<script>
  // Per-card RAW/FAIR toggle for the devig view.
  // Each .grid-toggle has two buttons; clicking one flips the parent
  // .price-grid's class between .show-raw and .show-fair. Cells emit
  // both <span class="raw"> and <span class="fair"> server-side; CSS
  // hides one based on the grid's class.
  function toggleDevigView(btn, view) {
    var toggle = btn.parentElement;            // .grid-toggle
    var grid   = toggle.nextElementSibling;    // .price-grid
    if (!grid || !grid.classList.contains('price-grid')) return;
    toggle.querySelectorAll('button').forEach(function(b) {
      b.classList.toggle('on', b.dataset.view === view);
    });
    grid.classList.remove('show-raw', 'show-fair');
    grid.classList.add(view === 'raw' ? 'show-raw' : 'show-fair');
  }
  // Expose for inline onclick handlers.
  window.toggleDevigView = toggleDevigView;
</script>
```

> Make sure this `<script>` is rendered into the page that contains the bets cards (not the parlay or trifecta tabs). The cleanest spot is alongside the existing `placeBet` JS — search for `function placeBet(btn)` to find that block, and add this script in the same `<script>...</script>` tag or directly above it.

- [ ] **Step 2: Manually verify the toggle works**

Start the dashboard server (or rely on a running instance) and open the bets tab. Pick any card. By default every cell should show devigged American odds (e.g. WZ shows `+129` instead of `+120`). Click `RAW` on the toggle — the same cells should switch to raw American odds (`+120`). Click `FAIR` — back to devigged.

If the toggle doesn't work, open the browser console and check for `toggleDevigView is not defined` or selector errors. The most likely cause is the `<script>` block being added to the wrong tab's template.

- [ ] **Step 3: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-pr-b-display" && \
  git add "Answer Keys/MLB Dashboard/mlb_dashboard.R" && \
  git commit -m "$(cat <<'EOF'
feat(bets-tab): wire RAW/FAIR toggle click handler

JS handler swaps the grid's .show-raw / .show-fair class on click,
and updates the active button's .on state. Default-FAIR rendered
state is preserved across refreshes; per-card state, no persistence.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 7 — Kelly Calculator widget

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` — add HTML block + CSS + JS in the bets-tab page template

- [ ] **Step 1: Add the Kelly Calculator HTML block right after `.sizing-controls`**

Open `Answer Keys/MLB Dashboard/mlb_dashboard.R`. Locate the existing `Bankroll/Kelly Controls` block around line 2891-2904 (the `tags$div(class = "sizing-controls", ...)` block). Immediately AFTER the closing of `.sizing-controls` (after line 2904) and BEFORE the `# Filter Bar` block (line 2907), insert:

```r
          # Kelly Calculator widget — manual no-vig Kelly sizing tool.
          # Reads dashboard's Bankroll + Kelly Fraction live; computes
          # via the existing calculateKellyBet JS function. Two text
          # inputs (Odds, Fair); live output (Risk, To Win, edge, Kelly%,
          # breakeven). Read-only — does not mutate any bet card.
          tags$div(class = "kelly-calc",
            tags$div(class = "kelly-label",
              tags$div(class = "h", "KELLY CALCULATOR"),
              tags$div(class = "s", "manual sizing tool")
            ),
            tags$div(class = "kelly-fields",
              tags$div(class = "field",
                tags$label("Offered Odds"),
                tags$input(id = "kc-odds", type = "text",
                           class = "kc-input", placeholder = "+120")
              ),
              tags$span(class = "arrow", intToUtf8(0xB7)),  # middle dot
              tags$div(class = "field",
                tags$label("Your Fair"),
                tags$input(id = "kc-fair", type = "text",
                           class = "kc-input", placeholder = "52.5% or -110")
              )
            ),
            tags$div(class = "kelly-out",
              tags$div(class = "risk-row",
                tags$span(class = "risk-label", "Risk"),
                tags$span(class = "risk-value", id = "kc-risk", "$0")
              ),
              tags$div(class = "towin",
                "to win ", tags$span(id = "kc-towin", "$0")
              ),
              tags$div(class = "detail",
                tags$span("edge"),
                tags$span(id = "kc-edge", class = "chip", "-"),
                tags$span(class = "sep", intToUtf8(0xB7)),
                tags$span("Kelly"),
                tags$span(id = "kc-kelly", "-"),
                tags$span(class = "sep", intToUtf8(0xB7)),
                tags$span("BE"),
                tags$span(id = "kc-be", "-")
              )
            )
          ),
```

(Note: the trailing comma is intentional; the next sibling in this `tagList` is the filter bar.)

- [ ] **Step 2: Add the CSS for the widget**

In the same file, inside the same `<style>...</style>` block where `.sizing-controls` is defined, add:

```css
.kelly-calc {
  display: grid;
  grid-template-columns: auto 1fr auto;
  gap: 24px;
  align-items: center;
  padding: 14px 18px;
  background: linear-gradient(135deg,
    rgba(63, 185, 80, 0.06) 0%,
    rgba(63, 185, 80, 0.02) 50%,
    #161b22 100%);
  border: 1px solid #2a3442;
  border-left: 3px solid #3fb950;
  border-radius: 10px;
  margin: 0 0 14px 0;
}
.kelly-calc .kelly-label .h {
  font-family: 'SF Mono', monospace;
  font-size: 11px;
  letter-spacing: 0.14em;
  color: #e6edf3;
  font-weight: 600;
}
.kelly-calc .kelly-label .s {
  color: #7d8590;
  font-size: 11px;
  margin-top: 2px;
}
.kelly-calc .kelly-fields {
  display: flex; align-items: center; gap: 14px; flex-wrap: wrap;
}
.kelly-calc .field {
  display: flex; flex-direction: column; gap: 4px;
}
.kelly-calc .field label {
  font-size: 10px; letter-spacing: 0.12em;
  text-transform: uppercase; color: #7d8590;
}
.kelly-calc .kc-input {
  background: #0d1117;
  border: 1px solid #2a3442;
  color: #e6edf3;
  font-family: 'SF Mono', monospace;
  font-size: 16px;
  font-variant-numeric: tabular-nums;
  font-weight: 500;
  padding: 8px 12px;
  border-radius: 7px;
  width: 110px;
  text-align: center;
  transition: all 0.15s ease;
}
.kelly-calc .kc-input:focus {
  outline: none;
  border-color: #3fb950;
  background: #1c222b;
  box-shadow: 0 0 0 3px rgba(63, 185, 80, 0.15);
}
.kelly-calc .kc-input.invalid {
  border-color: #f85149;
  box-shadow: 0 0 0 3px rgba(248, 81, 73, 0.12);
}
.kelly-calc .arrow {
  color: #7d8590;
  font-family: 'SF Mono', monospace;
  font-size: 18px;
  margin-top: 14px;
}
.kelly-calc .kelly-out {
  text-align: right;
  display: flex; flex-direction: column; gap: 4px;
}
.kelly-calc .risk-row {
  display: flex; align-items: baseline; justify-content: flex-end; gap: 8px;
}
.kelly-calc .risk-label {
  font-size: 10px; letter-spacing: 0.14em;
  text-transform: uppercase; color: #7d8590;
}
.kelly-calc .risk-value {
  font-family: 'SF Mono', monospace;
  font-size: 28px; font-weight: 600;
  color: #3fb950; font-variant-numeric: tabular-nums;
}
.kelly-calc .risk-value.neg { color: #7d8590; }
.kelly-calc .towin {
  font-family: 'SF Mono', monospace;
  font-size: 11px; color: #7d8590;
}
.kelly-calc .detail {
  font-family: 'SF Mono', monospace;
  font-size: 11px; color: #7d8590;
  font-variant-numeric: tabular-nums;
}
.kelly-calc .detail .chip {
  padding: 2px 7px; border-radius: 4px; margin-left: 6px;
}
.kelly-calc .detail .chip.pos { color: #3fb950; background: rgba(63, 185, 80, 0.10); }
.kelly-calc .detail .chip.neg { color: #f85149; background: rgba(248, 81, 73, 0.10); }
.kelly-calc .detail .sep { color: #3a4658; margin: 0 6px; }
```

- [ ] **Step 3: Add the JS that wires the widget to the dashboard's bankroll/kelly inputs**

Find the existing `<script>` block in the bets-tab page template that defines `placeBet`, `recalculateBetSizes`, etc. Add this Kelly-calc setup (wrap in DOMContentLoaded if there isn't already a wrapper):

```html
<script>
  (function setupKellyCalc() {
    var $ = function(id) { return document.getElementById(id); };
    function parseAmerican(s) {
      if (s == null) return null;
      s = String(s).trim().replace(/\s/g, '');
      if (!s) return null;
      var n = Number(s);
      if (!isFinite(n) || n === 0) return null;
      return n;
    }
    function parseFair(s) {
      if (s == null) return null;
      s = String(s).trim();
      if (!s) return null;
      if (s.endsWith('%')) {
        var p = Number(s.slice(0, -1));
        return isFinite(p) ? p / 100 : null;
      }
      var n = Number(s);
      if (!isFinite(n)) return null;
      if (n > 0 && n < 1) return n;          // 0.525 → 52.5%
      if (n > 1 && n < 100) return n / 100;  // 52.5  → 52.5%
      // Otherwise: American fair odds (e.g. -110, +120)
      return n > 0 ? 100 / (n + 100) : -n / (-n + 100);
    }
    function fmtMoney(x) {
      if (x == null || !isFinite(x) || x <= 0) return '$0';
      return '$' + Math.round(x).toLocaleString();
    }
    function fmtPct(x, signed) {
      if (x == null || !isFinite(x)) return '-';
      var v = (x * 100).toFixed(1);
      return (signed && x > 0 ? '+' : '') + v + '%';
    }

    function recompute() {
      var oddsEl = $('kc-odds');
      var fairEl = $('kc-fair');
      if (!oddsEl || !fairEl) return;
      var american = parseAmerican(oddsEl.value);
      var fairProb = parseFair(fairEl.value);
      oddsEl.classList.toggle('invalid', oddsEl.value !== '' && american == null);
      fairEl.classList.toggle('invalid', fairEl.value !== '' && fairProb == null);

      var bankroll  = Number(($('bankroll-input') || {}).value) || 0;
      var kellyMult = Number(($('kelly-input')     || {}).value) || 0;

      var risk = 0, towin = 0, edge = NaN, kellyPct = NaN, be = NaN;
      if (american != null && fairProb != null && bankroll > 0 && kellyMult > 0) {
        var dec = american > 0 ? american / 100 + 1 : 100 / Math.abs(american) + 1;
        be = american > 0 ? 100 / (american + 100) : Math.abs(american) / (Math.abs(american) + 100);
        var b = dec - 1, p = fairProb, q = 1 - p;
        kellyPct = Math.max(0, (b * p - q) / b);
        risk = Math.min(bankroll, bankroll * kellyMult * kellyPct);
        towin = american > 0 ? risk * american / 100 : risk * 100 / Math.abs(american);
        edge = p - be;
      }

      $('kc-risk').textContent  = fmtMoney(risk);
      $('kc-risk').classList.toggle('neg', !(risk > 0));
      $('kc-towin').textContent = fmtMoney(towin);
      var edgeEl = $('kc-edge');
      edgeEl.textContent = isNaN(edge) ? '-' : fmtPct(edge, true);
      edgeEl.className   = 'chip' + (isFinite(edge) ? (edge > 0 ? ' pos' : ' neg') : '');
      $('kc-kelly').textContent = isNaN(kellyPct) ? '-' : fmtPct(kellyPct, false);
      $('kc-be').textContent    = isNaN(be) ? '-' : fmtPct(be, false);
    }

    ['kc-odds', 'kc-fair', 'bankroll-input', 'kelly-input'].forEach(function(id) {
      var el = document.getElementById(id);
      if (el) el.addEventListener('input', recompute);
    });
    recompute();
  })();
</script>
```

- [ ] **Step 4: Manually verify the calculator works**

Reload the bets tab. The Kelly Calculator strip should appear directly under the Bankroll/Kelly Fraction strip. Test:
- Enter `+120` in Odds, `52.5%` in Fair → Risk should compute (with bankroll $5000, kelly 0.5 → ~$76).
- Enter `-110` in Fair → should also work (American fair odds path).
- Enter `40%` in Fair against `+120` Odds → Risk should grey to `$0` (fair < breakeven).
- Tweak the dashboard's Bankroll input → Risk should recompute live without re-typing odds/fair.

- [ ] **Step 5: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-pr-b-display" && \
  git add "Answer Keys/MLB Dashboard/mlb_dashboard.R" && \
  git commit -m "$(cat <<'EOF'
feat(bets-tab): Kelly Calculator widget

Manual no-vig Kelly sizing tool below the bankroll/kelly settings strip.
Two text inputs (Odds, Fair); live output (Risk, To Win, edge, Kelly%,
breakeven). Reads dashboard bankroll + kelly_mult live, so a single
source of truth for sizing math. Fair input accepts %, decimal, or
American fair odds. Negative-EV greys Risk to $0; invalid input red-tints
the offending field. Read-only — does not mutate any bet card.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 8 — Documentation updates

**Files:**
- Modify: `Answer Keys/CLAUDE.md`
- Modify: `Answer Keys/MLB Dashboard/README.md`

- [ ] **Step 1: Update `Answer Keys/CLAUDE.md`**

Find the "MLB Dashboard — Odds screen + WZ single-bet placer" section. Append a subsection:

```markdown
### Per-cell devig toggle + Kelly Calculator widget (PR B, 2026-05)

- Bets-tab cards have a per-card `RAW` / `FAIR` toggle in the top
  corner of the price grid. Default is `FAIR`: each cell shows the
  book's own probit-devigged American odds (math via
  `book_cell.R::.devig_american_pair`, parity-tested against
  `Tools.R::devig_american` in `tests/test_devig_pair_matches_tools.R`).
  RAW shows the original book quotes. Per-card state — no persistence,
  no global preference.
- A "Kelly Calculator" strip sits below the Bankroll/Kelly Fraction
  settings on the bets tab. Manual sizing tool: type Odds + Fair, get
  recommended Risk and To Win. Reuses the existing
  `calculateKellyBet` JS and reads the dashboard's bankroll/kelly_mult
  inputs live.
```

- [ ] **Step 2: Update `Answer Keys/MLB Dashboard/README.md`**

Find the "Features" or equivalent section. Add bullets:

```markdown
- **Kelly Calculator** — manual sizing tool below the settings strip on the
  bets tab. Type American odds + a fair % (or American fair odds), get a
  recommended Risk based on your bankroll × Kelly fraction. Useful when
  you've devigged a market manually and want to size at your own fair.

- **Per-cell devig toggle** — every bet card's price grid has a
  `RAW / FAIR` toggle. FAIR is the default view: each book cell shows
  its own probit-devigged American odds (computed against the book's
  own two-sided quote). Click `RAW` to see the original book prices.
```

- [ ] **Step 3: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-pr-b-display" && \
  git add "Answer Keys/CLAUDE.md" "Answer Keys/MLB Dashboard/README.md" && \
  git commit -m "$(cat <<'EOF'
docs: Kelly Calculator + per-cell devig toggle

Notes the new bets-tab features in CLAUDE.md (architecture pointers,
math source) and README.md (user-facing description).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 9 — Final verification + handoff

**Files:** none (verification only).

- [ ] **Step 1: Run the full test suite**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-pr-b-display/Answer Keys" && \
  for t in tests/test_book_cell.R tests/test_devig_pair_matches_tools.R tests/test_odds_screen.R; do
    echo "=== $t ==="; Rscript -e "testthat::test_file(\"$t\")" || break
  done
```

Expected: all green.

- [ ] **Step 2: End-to-end dashboard sanity check**

Reload the dashboard. On the bets tab:

1. **Kelly Calculator** is visible below the Bankroll/Kelly settings; basic typed inputs produce sensible output.
2. **Every card defaults to FAIR view** — cells show devigged American odds (e.g. WZ on a -140/+120 spread shows `+129` instead of `+120`).
3. **Toggle works** — click RAW on any card; that card's cells switch to raw odds; other cards remain on FAIR (per-card state).
4. **Pick book** still has the green border in both views.
5. **Off-line cells** (e.g. B105 on `0` for a `-0.5` bet) still show the amber background + tiny line tag in both views.
6. **Empty cells** show `—` in both views.
7. **No JS errors** in the browser console.

- [ ] **Step 3: Summarize the branch state**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-pr-b-display" && \
  git log --oneline main..HEAD && \
  git diff --stat main..HEAD
```

Expected: ~8 commits ahead of main. Diff stat shows `book_cell.R`, `mlb_dashboard.R`, the two new test files, `Answer Keys/CLAUDE.md`, `MLB Dashboard/README.md`.

Hand off to the user with: "PR B is ready. Kelly Calculator widget + per-cell RAW/FAIR devig toggle (default FAIR, probit math). Tests green, dashboard verified. Ready to merge to main when you give the word."

---

## Out-of-band notes

- **Do NOT merge to main without explicit user approval.** Project policy.
- **PR A must already be on main** before this PR runs end-to-end verification — Task 9 step 2 sub-bullet 5 ("off-line cells show amber tag in both views") relies on PR A's bug-3 fix (otherwise the *opposite-row* cells will all show spurious amber tags from the unrelated bug).
- **The legacy `render_book_pill` shim is preserved** — the new `opposite_american_odds` parameter has a default of NA so backward callers (the legacy `create_bets_table_legacy` fallback path) keep working unchanged.
- **No persistence of toggle state.** Refresh resets every card to FAIR. If you decide later you want to remember the toggle per-card across refreshes, that's a separate change (would need `dashboard_settings` table entries).
