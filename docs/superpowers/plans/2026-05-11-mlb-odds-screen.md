# MLB Odds-Screen + WZ Singles Placer — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the MLB Dashboard's flat bets table with a card-layout odds screen that shows every bet at every tracked sportsbook (5 bettable + 3 reference), with amber-tagged mismatched-line pills, plus a Wagerzon direct-API auto-placer for singles mirroring what parlays already have.

**Architecture:** Pipeline (`MLB.R` + new `odds_screen.R` helper) stops discarding the per-book bet comparison frames it already computes and writes them to a new `mlb_bets_book_prices` DuckDB table. Dashboard (`mlb_dashboard.R`) replaces `create_bets_table()` with a card-shaped reactable using a new `render_book_pill()` helper. Server (`mlb_dashboard_server.py`) makes `/api/place-bet` a dispatcher: Wagerzon → new `wagerzon_odds/single_placer.py` (REST); Hoop88/BFA/BetOnline → existing Playwright path; others → 400.

**Tech Stack:** R (testthat, reactable, htmltools, DuckDB), Python (Flask, DuckDB, requests), Wagerzon REST API.

**Spec:** `Answer Keys/MLB Dashboard/PLAN_odds_screen.md` (commits `bdeefba` + `39f2cd7` on this branch).

---

## File structure

### Files created

| Path | Responsibility |
|---|---|
| `Answer Keys/MLB Answer Key/odds_screen.R` | Pure helper `expand_bets_to_book_prices(bets, book_odds_by_book)` returning the long-format frame written to `mlb_bets_book_prices`. ±1 unit matching logic lives here. |
| `Answer Keys/MLB Dashboard/book_pill.R` | Pure helper `render_book_pill(book, american_odds, line_quoted, model_line, is_pick)`. Three pill states + pick highlight. Mirrors the precedent set by `Answer Keys/MLB Dashboard/books_strip.R`. |
| `Answer Keys/tests/test_odds_screen.R` | testthat coverage for the matching logic — exact / mismatched-within-1 / mismatched-outside-1 / tiebreak. |
| `Answer Keys/tests/test_book_pill.R` | testthat coverage for `render_book_pill` HTML output across the 5 visual states. |
| `Answer Keys/MLB Dashboard/migrations/002_single_placer_columns.py` | Idempotent ADD COLUMN migration for `placed_bets` (account, status, ticket_number, error_msg, error_msg_key, wz_odds_at_place). |
| `wagerzon_odds/single_placer.py` | `place_single(account, bet)` — load session, build 1-leg payload, preflight + drift check, submit, extract ticket, return structured result. |
| `wagerzon_odds/tests/test_single_placer.py` | Unit tests for sel encoding, preflight drift detection, ticket extraction, error taxonomy. Mocks WZ HTTP. |
| `Answer Keys/MLB Dashboard/tests/test_place_bet_dispatch.py` | Flask test client coverage for `/api/place-bet` dispatch by `bookmaker_key`. |

### Files modified

| Path | What changes |
|---|---|
| `Answer Keys/MLB Answer Key/MLB.R` | Expand Odds API `markets=` (~line 155). Add `bet_row_id` hash on `all_bets_combined`. Call `expand_bets_to_book_prices` between dedup (line 787) and Kelly correlation (~line 820). DROP+WRITE `mlb_bets_book_prices` alongside `mlb_bets_combined` (~line 851). |
| `Answer Keys/MLB Dashboard/mlb_dashboard.R` | New loader for `mlb_bets_book_prices` with try/catch fallback. New pivot long→wide. `create_bets_table()` body replaced with card layout. New CSS classes in the inline style block. New JS dispatch `placeBet()`. Correlation dot relocated to corner badge. |
| `Answer Keys/MLB Dashboard/mlb_dashboard_server.py` | `/api/place-bet` becomes a dispatcher (WZ / Playwright / 400 manual-log). Imports `single_placer` from `wagerzon_odds`. |
| `Answer Keys/MLB Dashboard/README.md` | Bets-tab section rewritten: card layout sketch, pill states, Place/Log dispatch behavior, schema notes. |
| `Answer Keys/CLAUDE.md` | New `mlb_bets_book_prices` table + new `placed_bets` columns listed in the DuckDB databases section. |

### Files not changed (intentional)

- `bet_placer/placer.py` — existing Playwright path stays as-is for Hoop88/BFA/BetOnline.
- `wagerzon_odds/parlay_pricer.py` — `single_placer.py` imports its shared session/preflight/error-classification helpers but doesn't modify it.
- `Answer Keys/MLB Answer Key/Tools.R` — `expand_bets_to_book_prices` lives in its own file (`odds_screen.R`) to keep Tools.R focused.

### TDD note

Pure functions (`expand_bets_to_book_prices`, `render_book_pill`, `single_placer.place_single`, the `/api/place-bet` dispatcher) are TDD-driven. The dashboard's `create_bets_table()` body and the JS dispatch are verified by browser smoke test (no headless harness exists for the reactable; manual verification + the existing parlay tab's CSS pattern is the precedent).

---

## Task 1: `render_book_pill` helper + tests

**Files:**
- Create: `Answer Keys/MLB Dashboard/book_pill.R`
- Create: `Answer Keys/tests/test_book_pill.R`

- [ ] **Step 1: Write the failing tests**

Path: `Answer Keys/tests/test_book_pill.R`

```r
# Answer Keys/tests/test_book_pill.R
# Tests for render_book_pill — bets-tab pill renderer with 5 visual states.
# Run from "Answer Keys/" directory:
#   Rscript -e 'testthat::test_file("tests/test_book_pill.R")'

library(testthat)
source("../MLB Dashboard/book_pill.R")

test_that("exact-line pill renders book label + price only", {
  out <- render_book_pill(book = "WZ", american_odds = 120,
                         line_quoted = 5.5, model_line = 5.5, is_pick = FALSE)
  expect_match(out, 'class="pill"', fixed = TRUE)
  expect_match(out, '<span class="book">WZ</span>', fixed = TRUE)
  expect_match(out, '>+120<', fixed = TRUE)
  expect_false(grepl("line-tag", out))
})

test_that("negative odds render with minus sign, not plus", {
  out <- render_book_pill("FD", -115, 5.5, 5.5, FALSE)
  expect_match(out, '>-115<', fixed = TRUE)
  expect_false(grepl(">+-115<", out, fixed = TRUE))
})

test_that("mismatched-line pill renders amber line tag and amber pill class", {
  out <- render_book_pill("BFA", -115, line_quoted = 5.0, model_line = 5.5,
                         is_pick = FALSE)
  expect_match(out, 'class="pill mismatched"', fixed = TRUE)
  expect_match(out, '<span class="line-tag">O5</span>', fixed = TRUE)
  expect_match(out, '>-115<', fixed = TRUE)
})

test_that("mismatched-line under-side renders U<line> tag", {
  out <- render_book_pill("BFA", -105, line_quoted = 5.0, model_line = 5.5,
                         is_pick = FALSE, side = "under")
  expect_match(out, '<span class="line-tag">U5</span>', fixed = TRUE)
})

test_that("NA odds render as muted dashed pill with em-dash", {
  out <- render_book_pill("DK", NA_integer_, NA_real_, 5.5, FALSE)
  expect_match(out, 'class="pill muted"', fixed = TRUE)
  expect_match(out, '&mdash;', fixed = TRUE)
})

test_that("pick override applies green class on exact line", {
  out <- render_book_pill("H88", 125, 5.5, 5.5, is_pick = TRUE)
  expect_match(out, 'class="pill pick"', fixed = TRUE)
  expect_match(out, '<span class="book">H88</span>', fixed = TRUE)
})

test_that("pick override applies green class even on mismatched line", {
  out <- render_book_pill("H88", 125, line_quoted = 5.0, model_line = 5.5,
                         is_pick = TRUE)
  # green pick takes precedence visually; class string explicitly notes both
  expect_match(out, 'class="pill pick mismatched"', fixed = TRUE)
  expect_match(out, '<span class="line-tag">O5</span>', fixed = TRUE)
})

test_that("integer line renders without trailing .0 in the tag", {
  out <- render_book_pill("BFA", -115, 5.0, 5.5, FALSE)
  expect_match(out, '<span class="line-tag">O5</span>', fixed = TRUE)
  expect_false(grepl("O5.0", out))
})

test_that("half-run line renders with .5", {
  out <- render_book_pill("BFA", -115, 5.5, 6.0, FALSE)
  expect_match(out, '<span class="line-tag">O5.5</span>', fixed = TRUE)
})
```

- [ ] **Step 2: Run tests to verify they fail**

Run:
```bash
cd "Answer Keys"
Rscript -e 'testthat::test_file("tests/test_book_pill.R")'
```

Expected: all 9 tests FAIL with `could not find function "render_book_pill"`.

- [ ] **Step 3: Write minimal implementation**

Path: `Answer Keys/MLB Dashboard/book_pill.R`

```r
# Answer Keys/MLB Dashboard/book_pill.R
# Render one book pill for the bets-tab odds-screen card.
#
# Three visual states (driven by line_quoted vs model_line):
#   - exact-line match    -> standard pill, just book + price
#   - mismatched (within +/-1 unit) -> amber-tinted pill + amber line tag
#   - no quote (NA odds)  -> muted dashed pill with em-dash
# Pick override (is_pick = TRUE) applies the green pick class on top of
# whichever state the row is in.
#
# Companion file: Answer Keys/MLB Dashboard/books_strip.R (parlays-tab pill row)

#' Format a line value for the line tag (e.g., 5.0 -> "5", 5.5 -> "5.5").
#' Returns the bare line; the side prefix (O/U) is prepended by the caller.
.format_line_value <- function(x) {
  if (is.na(x)) return("")
  if (x == round(x)) format(as.integer(x))
  else sub("\\.?0+$", "", sprintf("%.2f", x))
}

#' Render one book pill.
#'
#' @param book Display label for the book (e.g., "WZ", "H88").
#' @param american_odds Integer odds (e.g., 125, -110). NA -> no-quote pill.
#' @param line_quoted The line this book is actually showing on this side.
#' @param model_line The model's line (comparison anchor).
#' @param is_pick TRUE if this is the pick book on the pick side; applies green.
#' @param side Either "over" (default for totals over / favorite spread) or
#'   "under" (totals under / dog spread). Drives O/U prefix on mismatched
#'   line tag. For spreads, pass "over" for the favored side and "under" for
#'   the dog side; the prefix becomes the sign on the line ("+1" / "-1").
#' @return HTML string for the pill.
render_book_pill <- function(book, american_odds, line_quoted, model_line,
                             is_pick = FALSE, side = "over") {
  # State 1: no quote
  if (is.na(american_odds)) {
    return(sprintf('<span class="pill muted"><span class="book">%s</span>&mdash;</span>',
                   htmltools::htmlEscape(book)))
  }

  price_str <- if (american_odds > 0) paste0("+", american_odds) else as.character(american_odds)

  is_mismatched <- !is.na(line_quoted) && !is.na(model_line) && line_quoted != model_line

  classes <- "pill"
  if (is_pick) classes <- paste(classes, "pick")
  if (is_mismatched) classes <- paste(classes, "mismatched")

  tag_html <- ""
  if (is_mismatched) {
    prefix <- if (side == "under") "U" else "O"
    tag_html <- sprintf('<span class="line-tag">%s%s</span>',
                        prefix, .format_line_value(line_quoted))
  }

  sprintf('<span class="%s"><span class="book">%s%s</span>%s</span>',
          classes,
          htmltools::htmlEscape(book),
          ifelse(nzchar(tag_html), paste0(" ", tag_html), ""),
          price_str)
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run:
```bash
cd "Answer Keys"
Rscript -e 'testthat::test_file("tests/test_book_pill.R")'
```

Expected: 9 passed, 0 failed.

- [ ] **Step 5: Commit**

```bash
git add "Answer Keys/MLB Dashboard/book_pill.R" "Answer Keys/tests/test_book_pill.R"
git commit -m "feat(mlb-dashboard): render_book_pill helper with three pill states"
```

---

## Task 2: `expand_bets_to_book_prices` helper + tests

**Files:**
- Create: `Answer Keys/MLB Answer Key/odds_screen.R`
- Create: `Answer Keys/tests/test_odds_screen.R`

- [ ] **Step 1: Write the failing tests**

Path: `Answer Keys/tests/test_odds_screen.R`

```r
# Answer Keys/tests/test_odds_screen.R
# Tests for expand_bets_to_book_prices — long-format expansion of bets
# into per-book pill rows with +/-1 unit nearest-line matching.
#
# Run from "Answer Keys/" directory:
#   Rscript -e 'testthat::test_file("tests/test_odds_screen.R")'

library(testthat)
library(dplyr)
library(tibble)
source("../MLB Answer Key/odds_screen.R")

# A canonical 1-bet input frame matching the columns mlb_bets_combined
# carries: id, market, line, bet_on, market_type, plus the helper's
# computed bet_row_id (the caller hashes it in).
make_bet_row <- function(game_id = "g1",
                         market = "totals_1st_5_innings",
                         line = 5.5,
                         bet_on = "Over",
                         market_type = "totals") {
  tibble(
    bet_row_id  = "fakehash",
    game_id     = game_id,
    market      = market,
    market_type = market_type,
    period      = "F5",
    line        = line,
    bet_on      = bet_on,
    pick_side   = "pick"  # the bet itself is always the pick side
  )
}

# Per-book odds frame matching the shape scrapers + game_odds expose.
# market and side mirror Odds API conventions; line is the book's quote.
make_book_odds <- function(rows) {
  bind_rows(rows)
}

book_row <- function(game_id, market, period, side, line, american_odds,
                     fetch_time = as.POSIXct("2026-05-11 12:00:00", tz = "UTC")) {
  tibble(
    game_id = game_id, market = market, period = period, side = side,
    line = line, american_odds = american_odds, fetch_time = fetch_time
  )
}

test_that("exact-line quote at one book emits two rows (pick + opposite)", {
  bets <- make_bet_row()
  book_odds <- list(wagerzon = bind_rows(
    book_row("g1", "totals", "F5", "Over",  5.5, +120),
    book_row("g1", "totals", "F5", "Under", 5.5, -140)
  ))
  out <- expand_bets_to_book_prices(bets, book_odds)
  expect_equal(nrow(out), 2)
  expect_setequal(out$side, c("pick", "opposite"))
  expect_true(all(out$is_exact_line))
  expect_equal(out$line_quoted[out$side == "pick"], 5.5)
  expect_equal(out$american_odds[out$side == "pick"], 120)
})

test_that("missing book emits zero rows for that book", {
  bets <- make_bet_row()
  book_odds <- list(wagerzon = bind_rows(
    book_row("g1", "totals", "F5", "Over", 5.5, +120)
  ))
  # hoop88 absent from book_odds entirely
  out <- expand_bets_to_book_prices(bets, book_odds)
  expect_false("hoop88" %in% out$bookmaker)
})

test_that("nearest-line within +/-1 unit is taken with is_exact_line=FALSE", {
  bets <- make_bet_row(line = 5.5)
  book_odds <- list(bfa = bind_rows(
    book_row("g1", "totals", "F5", "Over",  5.0, -115),
    book_row("g1", "totals", "F5", "Under", 5.0, -105)
  ))
  out <- expand_bets_to_book_prices(bets, book_odds)
  expect_equal(nrow(out), 2)
  expect_false(any(out$is_exact_line))
  expect_equal(unique(out$line_quoted), 5.0)
})

test_that("lines outside +/-1 unit emit no row", {
  bets <- make_bet_row(line = 5.5)
  book_odds <- list(bfa = bind_rows(
    book_row("g1", "totals", "F5", "Over", 3.5, -115)  # 2 units off
  ))
  out <- expand_bets_to_book_prices(bets, book_odds)
  expect_equal(nrow(out), 0)
})

test_that("when both 0.5-unit and 1.0-unit lines exist, the closer one wins", {
  bets <- make_bet_row(line = 5.5)
  book_odds <- list(bfa = bind_rows(
    book_row("g1", "totals", "F5", "Over", 5.0, -115),   # 0.5 off
    book_row("g1", "totals", "F5", "Over", 4.5, -150),   # 1.0 off
    book_row("g1", "totals", "F5", "Under", 5.0, -105),
    book_row("g1", "totals", "F5", "Under", 4.5, +130)
  ))
  out <- expand_bets_to_book_prices(bets, book_odds)
  expect_equal(unique(out$line_quoted), 5.0)
})

test_that("equidistant tiebreak prefers the line worse for the bettor (over)", {
  # For an Over bet at 5.5, "worse for the bettor" on the Over side is
  # the HIGHER line (over 6 is harder to hit than over 5). When both 5
  # and 6 are equidistant (0.5 each), pick 6.
  bets <- make_bet_row(line = 5.5, bet_on = "Over")
  book_odds <- list(bfa = bind_rows(
    book_row("g1", "totals", "F5", "Over", 5.0, +110),
    book_row("g1", "totals", "F5", "Over", 6.0, -130),
    book_row("g1", "totals", "F5", "Under", 5.0, -130),
    book_row("g1", "totals", "F5", "Under", 6.0, +110)
  ))
  out <- expand_bets_to_book_prices(bets, book_odds) %>%
    filter(side == "pick")
  expect_equal(out$line_quoted, 6.0)
})

test_that("bet_row_id is preserved on every emitted row", {
  bets <- make_bet_row()
  book_odds <- list(wagerzon = bind_rows(
    book_row("g1", "totals", "F5", "Over", 5.5, +120),
    book_row("g1", "totals", "F5", "Under", 5.5, -140)
  ))
  out <- expand_bets_to_book_prices(bets, book_odds)
  expect_true(all(out$bet_row_id == "fakehash"))
})

test_that("multiple books on the same bet expand independently", {
  bets <- make_bet_row()
  book_odds <- list(
    wagerzon = bind_rows(
      book_row("g1", "totals", "F5", "Over",  5.5, +120),
      book_row("g1", "totals", "F5", "Under", 5.5, -140)
    ),
    bfa = bind_rows(
      book_row("g1", "totals", "F5", "Over",  5.0, -115),
      book_row("g1", "totals", "F5", "Under", 5.0, -105)
    )
  )
  out <- expand_bets_to_book_prices(bets, book_odds)
  expect_setequal(out$bookmaker, c("wagerzon", "bfa"))
  expect_equal(nrow(out), 4)
})

test_that("normalize_book_odds_frame splits market into (market, period)", {
  raw <- tibble(
    game_id = "g1",
    market_name = c("totals_1st_5_innings", "spreads_1st_3_innings",
                    "h2h", "alternate_totals"),
    bet_on = c("Over", "Home", "Home", "Over"),
    line = c(5.5, -1.5, NA, 8.5),
    american_odds = c(-110L, +120L, -150L, +105L),
    fetch_time = rep(as.POSIXct("2026-05-11 12:00", tz="UTC"), 4)
  )
  out <- normalize_book_odds_frame(raw)
  expect_setequal(out$market, c("totals", "spreads", "h2h", "alternate_totals"))
  expect_setequal(out$period, c("F5", "F3", "FG", "FG"))
  expect_equal(out$side, c("Over", "Home", "Home", "Over"))
})

test_that("normalize_book_odds_frame handles F7 and alts", {
  raw <- tibble(
    game_id = "g1",
    market_name = c("totals_1st_7_innings", "alternate_spreads_1st_5_innings"),
    bet_on = c("Under", "Away"),
    line = c(6.5, 1.5),
    american_odds = c(-115L, -120L),
    fetch_time = rep(as.POSIXct("2026-05-11 12:00", tz="UTC"), 2)
  )
  out <- normalize_book_odds_frame(raw)
  expect_equal(out$market, c("totals", "alternate_spreads"))
  expect_equal(out$period, c("F7", "F5"))
})
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd "Answer Keys"
Rscript -e 'testthat::test_file("tests/test_odds_screen.R")'
```

Expected: 8 tests FAIL with `could not find function "expand_bets_to_book_prices"`.

- [ ] **Step 3: Write minimal implementation**

Path: `Answer Keys/MLB Answer Key/odds_screen.R`

```r
# Answer Keys/MLB Answer Key/odds_screen.R
# Pure helpers for the bets-tab odds-screen pipeline path:
#   - normalize_book_odds_frame()    -- raw scraper rows -> canonical shape
#   - expand_bets_to_book_prices()   -- +/-1 unit nearest-line matching
#
# Sourced by MLB.R between the dedup step and the table write.
# Tested by tests/test_odds_screen.R.

library(dplyr)
library(tibble)
library(stringr)

LINE_MATCH_TOLERANCE <- 1.0  # max abs(line_quoted - model_line) we'll emit

#' Normalize a raw scraper / Odds API frame to the canonical shape consumed
#' by expand_bets_to_book_prices.
#'
#' Maps the full Odds API market name (e.g., "totals_1st_5_innings") to a
#' (market, period) tuple (e.g., market="totals", period="F5"). Full-game
#' markets ("h2h", "totals", "spreads", "alternate_totals", "alternate_spreads")
#' map to period="FG".
#'
#' @param raw A tibble with columns: game_id, market_name, bet_on, line,
#'   american_odds, fetch_time. (Some scrapers use `market` instead of
#'   `market_name` — the caller renames first.)
#' @return Tibble with columns game_id, market, period, side, line,
#'   american_odds, fetch_time.
normalize_book_odds_frame <- function(raw) {
  if (nrow(raw) == 0) {
    return(tibble(game_id = character(), market = character(),
                  period = character(), side = character(),
                  line = numeric(), american_odds = integer(),
                  fetch_time = as.POSIXct(character())))
  }
  raw %>%
    mutate(
      # Period extraction: look for _1st_<N>_innings suffix; default FG.
      period = case_when(
        str_detect(market_name, "_1st_3_innings$") ~ "F3",
        str_detect(market_name, "_1st_5_innings$") ~ "F5",
        str_detect(market_name, "_1st_7_innings$") ~ "F7",
        TRUE                                       ~ "FG"
      ),
      # Strip period suffix to get the base market type.
      market = str_replace(market_name, "_1st_[357]_innings$", "")
    ) %>%
    transmute(
      game_id, market, period,
      side = bet_on,
      line, american_odds = as.integer(american_odds),
      fetch_time
    )
}

#' Pick the closest matching line from a candidate set of book quotes.
#'
#' @param candidates A subset of one book's odds for one (game, market,
#'   period, side). Must have columns line and american_odds.
#' @param model_line The model's line (comparison anchor).
#' @param bet_on Direction string from the bet ("Over"/"Under" for totals,
#'   "Home"/"Away" or canonical team name for spreads). Used only for the
#'   equidistant tiebreaker.
#' @return A 1-row tibble of the chosen quote (with is_exact_line,
#'   line_quoted, american_odds) or an empty tibble if no quote is within
#'   LINE_MATCH_TOLERANCE.
.pick_closest_line <- function(candidates, model_line, bet_on) {
  if (nrow(candidates) == 0) {
    return(tibble(line_quoted = numeric(), american_odds = integer(),
                  is_exact_line = logical()))
  }
  c2 <- candidates %>%
    mutate(.dist = abs(line - model_line)) %>%
    filter(.dist <= LINE_MATCH_TOLERANCE)
  if (nrow(c2) == 0) {
    return(tibble(line_quoted = numeric(), american_odds = integer(),
                  is_exact_line = logical()))
  }
  min_dist <- min(c2$.dist)
  c3 <- c2 %>% filter(.dist == min_dist)
  if (nrow(c3) > 1) {
    # Equidistant tiebreaker: prefer the line worse for the bettor.
    # Over side: higher line is harder to hit -> pick the higher line.
    # Under side: lower line is harder to hit -> pick the lower line.
    # Spread favorite: more negative line (e.g. -2 vs -1) is harder -> pick lower.
    # Spread dog: more positive line is easier; "worse for bettor" is lower.
    # Implementation: for over/favorite-style picks we pick max(line),
    # for under/dog-style we pick min(line). We approximate this from
    # the bet_on string: starts with "Over" or is the favored side ->
    # pick higher line. Otherwise pick lower line. Caller hands us
    # whatever convention exists in their data.
    pick_high <- grepl("^(Over|over)", bet_on) || (model_line < 0)  # spread favorite has negative line
    if (pick_high) c3 <- c3 %>% slice_max(line, n = 1)
    else c3 <- c3 %>% slice_min(line, n = 1)
  }
  row <- c3[1, , drop = FALSE]
  tibble(
    line_quoted   = row$line,
    american_odds = row$american_odds,
    is_exact_line = abs(row$line - model_line) < 1e-9
  )
}

#' Map a bet's bet_on string to the "pick" and "opposite" side labels used
#' in the per-book odds frames. For totals, the side names match exactly
#' (Over/Under). For moneyline and spreads, the convention here lines up
#' with what the per-book scrapers write to mlb_odds.
.sides_for_bet <- function(bet_on, market_type) {
  if (market_type == "totals") {
    list(pick = bet_on,
         opposite = if (grepl("^Over", bet_on)) "Under" else "Over")
  } else if (market_type %in% c("spreads", "moneyline")) {
    # For spreads/h2h, the per-book frames use the team name in `side`.
    # The expansion writes opposite = the OTHER team name on the same
    # row; that mapping is computed upstream and passed in as the bet's
    # opposite_side column. If absent, we fall back to "opposite".
    list(pick = bet_on, opposite = NA_character_)  # caller wires up
  } else {
    list(pick = bet_on, opposite = NA_character_)
  }
}

#' Main entry point.
#'
#' @param bets A frame with one row per recommended bet. Required columns:
#'   bet_row_id, game_id, market, market_type, period, line, bet_on.
#'   For spreads/moneyline, also opposite_side (the other team's canonical
#'   label).
#' @param book_odds_by_book Named list. Keys are bookmaker_key values
#'   (e.g., "wagerzon", "hoop88", "bfa", "bookmaker", "bet105",
#'   "draftkings", "fanduel", "pinnacle"). Each value is a tibble with
#'   columns game_id, market, period, side, line, american_odds,
#'   fetch_time. The `market` value should be the normalized market type
#'   (e.g., "totals", "spreads", "moneyline") matching what the caller
#'   passes in `bets$market_type`.
#' @return Long-format tibble with columns matching the
#'   mlb_bets_book_prices schema.
expand_bets_to_book_prices <- function(bets, book_odds_by_book) {
  if (nrow(bets) == 0) {
    return(tibble(bet_row_id = character(), game_id = character(),
                  market = character(), period = character(),
                  side = character(), bookmaker = character(),
                  line = numeric(), line_quoted = numeric(),
                  is_exact_line = logical(), american_odds = integer(),
                  fetch_time = as.POSIXct(character())))
  }

  out_rows <- vector("list", length(book_odds_by_book) * nrow(bets) * 2)
  k <- 0L

  for (i in seq_len(nrow(bets))) {
    bet <- bets[i, , drop = FALSE]
    sides <- .sides_for_bet(bet$bet_on, bet$market_type)

    side_labels <- c(pick = sides$pick, opposite = sides$opposite)

    for (book_name in names(book_odds_by_book)) {
      book_frame <- book_odds_by_book[[book_name]]
      if (is.null(book_frame) || nrow(book_frame) == 0) next

      for (slot in names(side_labels)) {
        side_value <- side_labels[[slot]]
        if (is.na(side_value)) next

        candidates <- book_frame %>%
          filter(game_id == bet$game_id,
                 market  == bet$market_type,
                 period  == bet$period,
                 side    == side_value)

        chosen <- .pick_closest_line(candidates, bet$line, bet$bet_on)
        if (nrow(chosen) == 0) next

        # fetch_time from the chosen row (re-find since pick_closest_line
        # returns a thin tibble without it).
        ft <- candidates %>%
          filter(abs(line - chosen$line_quoted) < 1e-9) %>%
          slice_head(n = 1) %>%
          pull(fetch_time)

        k <- k + 1L
        out_rows[[k]] <- tibble(
          bet_row_id    = bet$bet_row_id,
          game_id       = bet$game_id,
          market        = bet$market_type,
          period        = bet$period,
          side          = slot,
          bookmaker     = book_name,
          line          = bet$line,
          line_quoted   = chosen$line_quoted,
          is_exact_line = chosen$is_exact_line,
          american_odds = as.integer(chosen$american_odds),
          fetch_time    = ft
        )
      }
    }
  }

  if (k == 0L) {
    return(tibble(bet_row_id = character(), game_id = character(),
                  market = character(), period = character(),
                  side = character(), bookmaker = character(),
                  line = numeric(), line_quoted = numeric(),
                  is_exact_line = logical(), american_odds = integer(),
                  fetch_time = as.POSIXct(character())))
  }
  bind_rows(out_rows[1:k])
}
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
cd "Answer Keys"
Rscript -e 'testthat::test_file("tests/test_odds_screen.R")'
```

Expected: 8 passed, 0 failed.

- [ ] **Step 5: Commit**

```bash
git add "Answer Keys/MLB Answer Key/odds_screen.R" "Answer Keys/tests/test_odds_screen.R"
git commit -m "feat(mlb): expand_bets_to_book_prices helper with ±1 unit line matching"
```

---

## Task 3: Expand Odds API `markets=` in MLB.R

**Files:**
- Modify: `Answer Keys/MLB Answer Key/MLB.R` (the Odds API pull, ~line 150-160)

- [ ] **Step 1: Locate the current `markets=` string**

Run:
```bash
grep -n "markets=" "Answer Keys/MLB Answer Key/MLB.R" | head -5
```

Find the line that builds the Odds API URL. The spec calls out line 155 as the area; the exact line is whatever assignment includes the `markets=` query param.

- [ ] **Step 2: Read 40 lines around the markets= line to understand the URL structure**

Run:
```bash
sed -n '130,180p' "Answer Keys/MLB Answer Key/MLB.R"
```

You're looking for either a single concatenated `markets=h2h_1st_5_innings,totals_1st_5_innings,...` string, or a vector that gets `paste(collapse=",")`-ed in.

- [ ] **Step 3: Modify the markets list**

Replace the existing markets enumeration (currently F5 mains + F5 alt totals) with the expanded set. The new list:

```
h2h, totals, spreads,
alternate_totals, alternate_spreads,
h2h_1st_3_innings, totals_1st_3_innings, spreads_1st_3_innings,
h2h_1st_5_innings, totals_1st_5_innings, spreads_1st_5_innings,
alternate_totals_1st_5_innings, alternate_spreads_1st_5_innings,
h2h_1st_7_innings, totals_1st_7_innings, spreads_1st_7_innings
```

(That's 4 FG markets + 3 F3 + 5 F5 + 3 F7 = 15 markets total. The old call had 4; new call has 15.)

The exact edit will look like this (the actual surrounding code may differ — adapt the markets string in place):

Before:
```r
markets <- "h2h_1st_5_innings,totals_1st_5_innings,spreads_1st_5_innings,alternate_totals_1st_5_innings"
```

After:
```r
markets <- paste(c(
  "h2h", "totals", "spreads",
  "alternate_totals", "alternate_spreads",
  "h2h_1st_3_innings", "totals_1st_3_innings", "spreads_1st_3_innings",
  "h2h_1st_5_innings", "totals_1st_5_innings", "spreads_1st_5_innings",
  "alternate_totals_1st_5_innings", "alternate_spreads_1st_5_innings",
  "h2h_1st_7_innings", "totals_1st_7_innings", "spreads_1st_7_innings"
), collapse = ",")
```

- [ ] **Step 4: Run the pipeline locally and verify the new markets populated**

```bash
cd "Answer Keys/MLB Answer Key"
Rscript MLB.R 2>&1 | tail -40
```

Wait 3-8 minutes. After completion:

```bash
duckdb "/Users/callancapitolo/NFLWork/Answer Keys/mlb.duckdb" \
  "SELECT DISTINCT market FROM mlb_odds_temp ORDER BY market" 2>/dev/null
```

Expected: see the new market keys (e.g., `h2h_1st_3_innings`, `alternate_spreads`, `spreads`) in the output. If only the original 4 appear, the markets= edit didn't take effect — re-check the URL build.

- [ ] **Step 5: Check Odds API credit usage**

The Odds API exposes the remaining-credits header. Check today's pipeline log:
```bash
grep -i "remaining\|credits\|x-requests" /tmp/mlb_pipeline.log 2>/dev/null | tail -5
```
or `Answer Keys/MLB Answer Key/*.log` — wherever the run wrote its output.

Expected: a credit count consistent with ~135 additional credits used (vs. the prior baseline). If credit usage exceeds the daily budget, drop F3/F7 markets from the list and re-run.

- [ ] **Step 6: Commit**

```bash
git add "Answer Keys/MLB Answer Key/MLB.R"
git commit -m "feat(mlb): expand Odds API markets= for reference-book coverage on derivatives"
```

---

## Task 4: Wire `expand_bets_to_book_prices` into MLB.R + write the new table

**Files:**
- Modify: `Answer Keys/MLB Answer Key/MLB.R` (insertion point: between dedup at line 787 and Kelly correlation at ~line 820; also write block at ~line 851)

- [ ] **Step 1: Read the surrounding code**

```bash
sed -n '780,860p' "Answer Keys/MLB Answer Key/MLB.R"
```

Confirm:
- The dedup `group_by(...) %>% filter(ev == max(ev)) %>% slice_head(n=1)` produces `all_bets_combined`.
- The pre-dedup frame (the union of `wagerzon_bets`, `hoop88_bets`, etc.) is available — you may need to save it as `all_bets_combined_pre_dedup <- bind_rows(...)` before the dedup runs.
- The `mlb_bets_combined` DROP+WRITE block sits at ~line 851.

- [ ] **Step 2: Source the new helper near the top of MLB.R**

Find the existing `source(...)` block (near the top, where Tools.R and other helpers are loaded) and add:

```r
source(file.path(this_dir, "odds_screen.R"))
```

(Match the existing source-path convention. If existing sources use bare filenames, use `source("odds_screen.R")` instead.)

- [ ] **Step 3: Add `bet_row_id` to `all_bets_combined` before the dedup**

Just before the `group_by(...) %>% filter(ev == max(ev)) %>% slice_head(n=1)` line, add:

```r
all_bets_combined <- all_bets_combined %>%
  mutate(bet_row_id = vapply(
    paste(id, market, ifelse(is.na(line), "", line), bet_on, sep = "|"),
    digest::digest,
    character(1),
    algo = "md5"
  ))
```

(Use `library(digest)` if not already loaded; otherwise call `tools::md5sum` or any other deterministic hash you have on hand. The key is determinism across pipeline runs.)

Then after dedup, the `bet_row_id` carries through naturally.

- [ ] **Step 4: Build the per-book odds list for the helper**

The helper expects each book's odds in canonical shape:
`(game_id, market, period, side, line, american_odds, fetch_time)`.
Today MLB.R has raw scraper frames (`wagerzon_odds`, `hoop88_odds`,
`bfa_odds`, `bookmaker_odds`, `bet105_odds`) in scope, plus the Odds
API frame `game_odds`. We normalize each through
`normalize_book_odds_frame()` (defined in `odds_screen.R` from Task 2).

First, locate the raw frames in MLB.R:

```bash
grep -nE "wagerzon_odds|hoop88_odds|bfa_odds|bookmaker_odds|bet105_odds|game_odds" \
  "Answer Keys/MLB Answer Key/MLB.R" | head -20
```

Confirm the column names on each frame. Most scrapers' `mlb_odds` table
has columns like `(game_id, market, bet_on, line, american_odds, fetch_time)`.
If a scraper's column is `market_name` rather than `market`, rename in
the `transmute` below.

Add right before the `expand_bets_to_book_prices()` call:

```r
# Normalize each book's raw frame to the canonical shape.
# normalize_book_odds_frame expects (game_id, market_name, bet_on, line,
# american_odds, fetch_time). Adapt the rename() calls if a scraper's
# column happens to be named differently.
normalize <- function(raw) {
  raw %>% rename(market_name = market) %>% normalize_book_odds_frame()
}

book_odds_by_book <- list(
  wagerzon   = normalize(wagerzon_odds),
  hoop88     = normalize(hoop88_odds),
  bfa        = normalize(bfa_odds),
  bookmaker  = normalize(bookmaker_odds),
  bet105     = normalize(bet105_odds),
  draftkings = game_odds %>% filter(bookmaker_key == "draftkings") %>%
                 rename(market_name = market) %>%
                 normalize_book_odds_frame(),
  fanduel    = game_odds %>% filter(bookmaker_key == "fanduel") %>%
                 rename(market_name = market) %>%
                 normalize_book_odds_frame(),
  pinnacle   = game_odds %>% filter(bookmaker_key == "pinnacle") %>%
                 rename(market_name = market) %>%
                 normalize_book_odds_frame()
)
```

If any of those frame variable names don't exist in MLB.R's local
scope at this point (the scraper code may save to different names),
adjust to whatever's in scope. The normalization treats them
identically once they reach `normalize_book_odds_frame()`.

- [ ] **Step 5: Call the helper and write the table**

Right before the `mlb_bets_combined` DROP+WRITE block, add:

```r
book_prices_long <- expand_bets_to_book_prices(all_bets_combined,
                                              book_odds_by_book)

# Data-bug guard: pick book should always be on the exact line.
pick_mismatches <- book_prices_long %>%
  inner_join(all_bets_combined %>% select(bet_row_id, pick_book = bookmaker_key),
             by = "bet_row_id") %>%
  filter(side == "pick" & bookmaker == pick_book & !is_exact_line)
if (nrow(pick_mismatches) > 0) {
  warning(sprintf("[odds_screen] %d pick-book rows on mismatched line — investigate",
                  nrow(pick_mismatches)))
}
```

Then in the writer block (around line 851), add a sibling DROP+WRITE:

```r
con <- dbConnect(duckdb::duckdb(), MLB_MM_DB_PATH, read_only = FALSE)
dbExecute(con, "DROP TABLE IF EXISTS mlb_bets_book_prices")
dbWriteTable(con, "mlb_bets_book_prices", as.data.frame(book_prices_long))
dbDisconnect(con, shutdown = TRUE)
```

(If the existing writer block already wraps a single `dbConnect/dbDisconnect`, fold the new DROP+WRITE into it rather than opening a second connection.)

- [ ] **Step 6: Run the pipeline and verify the table**

```bash
cd "Answer Keys/MLB Answer Key"
Rscript MLB.R 2>&1 | tail -20
```

After completion (~3-8 min):

```bash
duckdb "/Users/callancapitolo/NFLWork/Answer Keys/mlb_mm.duckdb" "DESCRIBE mlb_bets_book_prices"
duckdb "/Users/callancapitolo/NFLWork/Answer Keys/mlb_mm.duckdb" \
  "SELECT COUNT(*) AS rows, COUNT(DISTINCT bet_row_id) AS bets, COUNT(DISTINCT bookmaker) AS books FROM mlb_bets_book_prices"
```

Expected output:
- DESCRIBE shows 11 columns matching the schema (bet_row_id, game_id, market, period, side, bookmaker, line, line_quoted, is_exact_line, american_odds, fetch_time).
- COUNT shows rows > 0, bets > 0, books between 1 and 8 (depending on slate coverage).

- [ ] **Step 7: Run the threshold-violation guard**

```bash
duckdb "/Users/callancapitolo/NFLWork/Answer Keys/mlb_mm.duckdb" \
  "SELECT COUNT(*) FROM mlb_bets_book_prices WHERE is_exact_line = FALSE AND ABS(line_quoted - line) > 1.0"
```

Expected: `0`. Any row > 0 indicates a bug in `.pick_closest_line` and should block the commit.

- [ ] **Step 8: Commit**

```bash
git add "Answer Keys/MLB Answer Key/MLB.R"
git commit -m "feat(mlb): write mlb_bets_book_prices with per-book pill rows"
```

This completes the pipeline-side work (Phases 1 + Odds API expansion from the spec).

---

## Task 5: Schema migration `002_single_placer_columns.py`

**Files:**
- Create: `Answer Keys/MLB Dashboard/migrations/002_single_placer_columns.py`

- [ ] **Step 1: Write the migration**

Path: `Answer Keys/MLB Dashboard/migrations/002_single_placer_columns.py`

```python
"""One-shot migration: add WZ single-placer columns to placed_bets.

Mirrors 001_combined_parlay_columns.py. Idempotent — safe to re-run.

Run once against the live mlb_dashboard.duckdb:

    python "Answer Keys/MLB Dashboard/migrations/002_single_placer_columns.py" \
        "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb"
"""
from __future__ import annotations
import sys
import duckdb


COLUMNS = [
    ("account",          "VARCHAR"),
    ("status",           "VARCHAR DEFAULT 'placed'"),
    ("ticket_number",    "VARCHAR"),
    ("error_msg",        "VARCHAR"),
    ("error_msg_key",    "VARCHAR"),
    ("wz_odds_at_place", "INTEGER"),
]


def run(db_path: str) -> None:
    con = duckdb.connect(db_path)
    try:
        existing = {row[0] for row in con.execute("DESCRIBE placed_bets").fetchall()}
        for name, ddl in COLUMNS:
            if name not in existing:
                con.execute(f"ALTER TABLE placed_bets ADD COLUMN {name} {ddl}")
    finally:
        con.close()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <path-to-mlb_dashboard.duckdb>", file=sys.stderr)
        sys.exit(1)
    run(sys.argv[1])
```

- [ ] **Step 2: Run against a copy first, then the live DB**

Make a backup and dry-run on the copy:

```bash
cp "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb" /tmp/mlb_dashboard_backup.duckdb
python3 "Answer Keys/MLB Dashboard/migrations/002_single_placer_columns.py" /tmp/mlb_dashboard_backup.duckdb
duckdb /tmp/mlb_dashboard_backup.duckdb "DESCRIBE placed_bets" | grep -E "account|status|ticket_number|error_msg|wz_odds_at_place"
```

Expected: 6 new columns listed.

Then run idempotently to verify safety:
```bash
python3 "Answer Keys/MLB Dashboard/migrations/002_single_placer_columns.py" /tmp/mlb_dashboard_backup.duckdb
```
Expected: completes silently (no error from re-adding an existing column).

Once verified, run against the live DB:
```bash
python3 "Answer Keys/MLB Dashboard/migrations/002_single_placer_columns.py" \
  "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb"
```

- [ ] **Step 3: Commit**

```bash
git add "Answer Keys/MLB Dashboard/migrations/002_single_placer_columns.py"
git commit -m "feat(mlb-dashboard): migration 002 — placed_bets columns for WZ single placer"
```

---

## Task 6: `wagerzon_odds/single_placer.py` + tests

**Files:**
- Create: `wagerzon_odds/single_placer.py`
- Create: `wagerzon_odds/tests/test_single_placer.py`

Before writing this task's code, read the parlay placer to understand its public helpers:

```bash
grep -n "def \|class \|^from\|^import" wagerzon_odds/parlay_pricer.py | head -40
```

You're looking for:
- The session-loader function (loads WZ auth cookies from disk).
- The preflight/ConfirmWagerHelper call signature.
- The MakeWagerHelper submit signature.
- The error-classification helpers (rejected/auth_error/network_error/orphaned).
- The sel encoding helper.

`wagerzon_odds/tests/test_scraper_specials.py` shows the import-from-tests convention: `sys.path.insert(0, str(Path(__file__).parent.parent))`.

- [ ] **Step 1: Write failing tests**

Path: `wagerzon_odds/tests/test_single_placer.py`

```python
"""Tests for wagerzon_odds.single_placer

Mocks the WZ HTTP session and verifies:
- sel string for a 1-leg payload matches parlay placer's encoding for one leg
- preflight drift detection: returns price_moved when WZ shows different odds
- successful path: returns {status: 'placed', ticket_number: ...}
- error classification: auth_error / rejected / network_error / orphaned
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import pytest
from unittest.mock import MagicMock, patch

import single_placer


def make_bet(**overrides):
    """Default 1-leg bet payload for tests."""
    bet = dict(
        bet_hash         = "abc123",
        bookmaker_key    = "wagerzon",
        idgm             = 100001,    # WZ game id
        play             = 1,         # 1 = home spread
        line             = -1.5,
        american_odds    = 110,
        kelly_bet        = 42,
        actual_size      = 42,
        wz_odds_at_place = 110,       # what client thought WZ was showing
    )
    bet.update(overrides)
    return bet


def test_sel_encoding_single_leg():
    bet = make_bet()
    sel = single_placer.build_sel_for_single(bet)
    # Same shape as the parlay placer's per-leg encoding
    assert sel == "1_100001_-1.5_110"


def test_preflight_drift_returns_price_moved():
    """If WZ's ConfirmWagerHelper returns odds different from wz_odds_at_place,
    placer must NOT call MakeWagerHelper and must return price_moved."""
    bet = make_bet(wz_odds_at_place=110)
    fake_session = MagicMock()
    fake_session.post.return_value.json.return_value = {
        "result": {"details": [{"Win": 4400, "Risk": 4000, "Odds": 100}]}
        # Odds 100 differs from client's 110 -> drift
    }
    result = single_placer.place_single(account="primary", bet=bet,
                                        session=fake_session)
    assert result["status"] == "price_moved"
    assert result.get("ticket_number") is None
    # MakeWagerHelper should NOT have been called
    make_calls = [c for c in fake_session.post.call_args_list
                  if "MakeWagerHelper" in str(c)]
    assert len(make_calls) == 0


def test_successful_placement_returns_ticket():
    bet = make_bet(wz_odds_at_place=110)
    fake_session = MagicMock()
    # ConfirmWagerHelper: no drift
    confirm_response = MagicMock()
    confirm_response.json.return_value = {
        "result": {"details": [{"Win": 4620, "Risk": 4200, "Odds": 110}]}
    }
    # MakeWagerHelper: success with ticket
    make_response = MagicMock()
    make_response.json.return_value = {
        "result": {"WagerNumber": "T-9876", "AvailBalance": 153.50}
    }
    fake_session.post.side_effect = [confirm_response, make_response]

    result = single_placer.place_single("primary", bet, session=fake_session)
    assert result["status"] == "placed"
    assert result["ticket_number"] == "T-9876"
    assert result["balance_after"] == pytest.approx(153.50)


def test_wz_rejection_returns_rejected_with_reason():
    bet = make_bet()
    fake_session = MagicMock()
    confirm_response = MagicMock()
    confirm_response.json.return_value = {
        "result": {"details": [{"Win": 4620, "Risk": 4200, "Odds": 110}]}
    }
    make_response = MagicMock()
    make_response.json.return_value = {
        "result": {"ErrorCode": "INSUFFICIENT_BALANCE",
                   "ErrorMessage": "Not enough funds"}
    }
    fake_session.post.side_effect = [confirm_response, make_response]

    result = single_placer.place_single("primary", bet, session=fake_session)
    assert result["status"] == "rejected"
    assert "INSUFFICIENT_BALANCE" in (result.get("error_msg_key") or "") + \
           (result.get("error_msg") or "")


def test_network_error_returns_network_error():
    import requests
    bet = make_bet()
    fake_session = MagicMock()
    fake_session.post.side_effect = requests.ConnectionError("connection reset")

    result = single_placer.place_single("primary", bet, session=fake_session)
    assert result["status"] == "network_error"


def test_auth_error_returns_auth_error():
    bet = make_bet()
    fake_session = MagicMock()
    auth_response = MagicMock()
    auth_response.status_code = 401
    auth_response.json.return_value = {"error": "session expired"}
    fake_session.post.return_value = auth_response

    result = single_placer.place_single("primary", bet, session=fake_session)
    assert result["status"] == "auth_error"


def test_empty_details_treated_as_rejected():
    """ConfirmWagerHelper sometimes returns details=[] (line pulled) — mirror
    the parlay placer's behavior of converting this into 'rejected'."""
    bet = make_bet()
    fake_session = MagicMock()
    confirm_response = MagicMock()
    confirm_response.json.return_value = {"result": {"details": []}}
    fake_session.post.return_value = confirm_response

    result = single_placer.place_single("primary", bet, session=fake_session)
    assert result["status"] == "rejected"
    assert "details" in (result.get("error_msg") or "").lower() or \
           result.get("error_msg_key") == "empty_details"
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd wagerzon_odds
python3 -m pytest tests/test_single_placer.py -v 2>&1 | tail -20
```

Expected: 7 tests FAIL with `ModuleNotFoundError: single_placer`.

- [ ] **Step 3: Implement `single_placer.py`**

Path: `wagerzon_odds/single_placer.py`

```python
"""Wagerzon direct-API placer for SINGLE bets.

Mirrors parlay_pricer.py's placement path but for 1-leg payloads. Reuses
the shared session loader, preflight pattern, and error taxonomy.

Public API:
    place_single(account: str, bet: dict, session=None) -> dict

Returns a dict with keys:
    status:         'placed' | 'price_moved' | 'rejected' | 'auth_error'
                    | 'network_error' | 'orphaned'
    ticket_number:  WZ ticket string when status == 'placed', else None
    balance_after:  Available balance snapshot from WZ when status == 'placed'
    error_msg:      Human-readable error text (None on success)
    error_msg_key:  Short key for pill rendering (None on success)
"""
from __future__ import annotations
from typing import Optional
import requests

# Reuse the parlay placer's session loader. Before implementing this
# import line, run:
#   grep -nE "def .*session|class.*Session" wagerzon_odds/parlay_pricer.py
# to find the actual function or class name. Common candidates in the
# existing parlay code path:
#   - parlay_pricer.load_session_for_account(account)
#   - parlay_pricer.get_session(account_label)
#   - wagerzon_accounts.load_session(label)  (in wagerzon_accounts.py)
# Replace the import below with whichever exists.
from parlay_pricer import load_session_for_account  # type: ignore

DRIFT_TOLERANCE_AMERICAN = 1   # +/- 1 cent on american odds triggers price_moved


def build_sel_for_single(bet: dict) -> str:
    """Encode a single-leg sel string matching parlay_pricer's per-leg format:
    "<play>_<idgm>_<points>_<odds>" where points may include sign."""
    play   = bet["play"]
    idgm   = bet["idgm"]
    points = bet["line"]
    odds   = bet["american_odds"]
    pts_str = points if isinstance(points, str) else _format_points(points)
    return f"{play}_{idgm}_{pts_str}_{odds}"


def _format_points(p: float) -> str:
    """Match parlay_pricer.py's points formatting: signed string, no
    trailing .0 for integers."""
    if p == int(p):
        return f"{int(p):+d}" if p < 0 else str(int(p))
    s = f"{p:+g}"
    return s if p < 0 else s.lstrip("+")


def _classify_response(resp_json: dict) -> Optional[dict]:
    """If the response is a known error shape, return a structured error.
    Otherwise return None (caller proceeds to the success path)."""
    result = resp_json.get("result") or {}
    if "ErrorCode" in result:
        return {
            "status": "rejected",
            "ticket_number": None,
            "balance_after": None,
            "error_msg": result.get("ErrorMessage") or result["ErrorCode"],
            "error_msg_key": result["ErrorCode"],
        }
    return None


def place_single(account: str, bet: dict,
                 session=None) -> dict:
    """Submit a single bet at Wagerzon.

    The `session` kwarg is for testing only; production code should let
    the function load the session from the account label.
    """
    sess = session if session is not None else load_session_for_account(account)

    # 1) ConfirmWagerHelper (preflight + drift check)
    confirm_payload = {
        "sel": build_sel_for_single(bet),
        "Risk": int(bet["actual_size"]),
    }
    try:
        confirm = sess.post("https://wagerzon.com/ConfirmWagerHelper.aspx",
                            data=confirm_payload)
    except requests.ConnectionError as e:
        return _network_error(str(e))
    except requests.RequestException as e:
        return _network_error(str(e))

    if getattr(confirm, "status_code", 200) == 401:
        return {"status": "auth_error", "ticket_number": None,
                "balance_after": None, "error_msg": "session expired",
                "error_msg_key": "session_expired"}

    confirm_json = confirm.json() if hasattr(confirm, "json") else confirm
    err = _classify_response(confirm_json)
    if err is not None:
        return err

    details = (confirm_json.get("result") or {}).get("details") or []
    if not details:
        return {"status": "rejected", "ticket_number": None,
                "balance_after": None,
                "error_msg": "Wagerzon returned empty details (line pulled?)",
                "error_msg_key": "empty_details"}

    wz_odds_now = details[0].get("Odds")
    if wz_odds_now is None:
        return {"status": "rejected", "ticket_number": None,
                "balance_after": None,
                "error_msg": "Wagerzon details missing Odds field",
                "error_msg_key": "missing_odds"}

    expected = bet["wz_odds_at_place"]
    if abs(int(wz_odds_now) - int(expected)) > DRIFT_TOLERANCE_AMERICAN:
        return {"status": "price_moved", "ticket_number": None,
                "balance_after": None,
                "error_msg": f"Odds drifted from {expected} to {wz_odds_now}",
                "error_msg_key": "drift"}

    # 2) MakeWagerHelper (real submission)
    make_payload = dict(confirm_payload)  # same payload shape; WZ accepts
    try:
        make = sess.post("https://wagerzon.com/MakeWagerHelper.aspx",
                         data=make_payload)
    except requests.RequestException as e:
        # Network failure AFTER preflight is an orphan candidate — log
        # forensics + return network_error so the dashboard surfaces it.
        return _network_error(str(e))

    make_json = make.json() if hasattr(make, "json") else make
    err = _classify_response(make_json)
    if err is not None:
        return err

    result = make_json.get("result") or {}
    ticket = result.get("WagerNumber")
    if not ticket:
        # WZ confirmed but didn't return a ticket; treat as orphan candidate.
        return {"status": "orphaned", "ticket_number": None,
                "balance_after": None,
                "error_msg": "submitted to WZ but no WagerNumber returned",
                "error_msg_key": "no_ticket"}

    return {"status": "placed", "ticket_number": str(ticket),
            "balance_after": result.get("AvailBalance"),
            "error_msg": None, "error_msg_key": None}


def _network_error(msg: str) -> dict:
    return {"status": "network_error", "ticket_number": None,
            "balance_after": None,
            "error_msg": msg, "error_msg_key": "network_error"}
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
cd wagerzon_odds
python3 -m pytest tests/test_single_placer.py -v
```

Expected: 7 passed.

- [ ] **Step 5: Commit**

```bash
git add wagerzon_odds/single_placer.py wagerzon_odds/tests/test_single_placer.py
git commit -m "feat(wagerzon): single_placer module with preflight, drift check, error taxonomy"
```

---

## Task 7: `/api/place-bet` dispatcher + tests

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard_server.py` (`/api/place-bet` route)
- Create: `Answer Keys/MLB Dashboard/tests/test_place_bet_dispatch.py`

- [ ] **Step 1: Read the current `/api/place-bet` implementation**

```bash
grep -n "place-bet\|place_bet" "Answer Keys/MLB Dashboard/mlb_dashboard_server.py" | head -10
```

Read 40 lines around the current handler. Today it inserts a row into `placed_bets` for the manual-log path. We're adding a dispatch branch on `bookmaker_key`.

- [ ] **Step 2: Write the failing test**

Path: `Answer Keys/MLB Dashboard/tests/test_place_bet_dispatch.py`

```python
"""Tests for /api/place-bet dispatcher.

Three branches:
  - wagerzon -> call single_placer.place_single, return its response
  - hoop88 / bfa / betonlineag -> spawn Playwright via existing auto-place
  - everything else -> 400 with manual-log message
"""
import json
import sys
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest

# The dashboard module imports its own modules from project root
HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))            # MLB Dashboard/
sys.path.insert(0, str(HERE.parents[3]))            # NFLWork/

import mlb_dashboard_server as srv


@pytest.fixture
def client():
    srv.app.config["TESTING"] = True
    with srv.app.test_client() as c:
        yield c


def _wz_bet(**overrides):
    payload = dict(
        bet_hash="abc123",
        bookmaker_key="wagerzon",
        account="primary",
        kelly_bet=42, actual_size=42,
        idgm=100001, play=1, line=-1.5,
        american_odds=110, wz_odds_at_place=110,
        game_id="g1", market="spreads_1st_5_innings",
        bet_on="Home", home_team="NYY", away_team="BOS",
    )
    payload.update(overrides)
    return payload


def test_dispatch_wagerzon_calls_single_placer(client):
    fake_result = {"status": "placed", "ticket_number": "T-9876",
                   "balance_after": 153.50,
                   "error_msg": None, "error_msg_key": None}
    with patch("mlb_dashboard_server.single_placer.place_single",
               return_value=fake_result) as mock_place:
        resp = client.post("/api/place-bet", json=_wz_bet())
    assert resp.status_code == 200
    body = resp.get_json()
    assert body["status"] == "placed"
    assert body["ticket_number"] == "T-9876"
    mock_place.assert_called_once()


def test_dispatch_hoop88_routes_to_playwright(client):
    with patch("mlb_dashboard_server._spawn_playwright_placer") as mock_pw:
        mock_pw.return_value = {"success": True, "message": "Browser launching..."}
        resp = client.post("/api/place-bet",
                           json=_wz_bet(bookmaker_key="hoop88"))
    assert resp.status_code == 200
    mock_pw.assert_called_once()


def test_dispatch_unsupported_book_returns_400(client):
    resp = client.post("/api/place-bet",
                       json=_wz_bet(bookmaker_key="draftkings"))
    assert resp.status_code == 400
    assert "manual log" in resp.get_json().get("error", "").lower()


def test_wagerzon_dispatch_writes_breadcrumb_before_call(client):
    """The breadcrumb (status='placing') should be inserted before the placer
    call so that even if the process crashes mid-call we have a record."""
    fake_result = {"status": "placed", "ticket_number": "T-1",
                   "balance_after": 100.0, "error_msg": None, "error_msg_key": None}
    with patch("mlb_dashboard_server.single_placer.place_single",
               return_value=fake_result), \
         patch("mlb_dashboard_server._insert_placement_breadcrumb") as mock_bc:
        client.post("/api/place-bet", json=_wz_bet())
    mock_bc.assert_called_once()
    # breadcrumb must be called BEFORE place_single returns
    assert mock_bc.call_args.kwargs.get("status") == "placing" or \
           "placing" in str(mock_bc.call_args)
```

- [ ] **Step 3: Run tests to verify they fail**

```bash
cd "Answer Keys/MLB Dashboard"
python3 -m pytest tests/test_place_bet_dispatch.py -v 2>&1 | tail -20
```

Expected: 4 tests fail with import / attribute errors for `single_placer`, `_spawn_playwright_placer`, `_insert_placement_breadcrumb`.

- [ ] **Step 4: Add the import + dispatcher**

In `mlb_dashboard_server.py`, near the top of the file (with other imports):

```python
# Project-root sibling for single_placer
sys.path.insert(0, str(_REPO_ROOT / "wagerzon_odds"))
import single_placer  # type: ignore
```

(Use whatever import idiom matches the existing code's pattern for cross-module imports from `wagerzon_odds`.)

Add a constant near the existing `SUPPORTED_AUTO_BOOKS`:

```python
WZ_DIRECT_API_BOOKS = {"wagerzon"}
PLAYWRIGHT_BOOKS = {"hoop88", "bfa", "betonlineag"}
```

Add the helper for the Playwright path (extract the existing logic from `/api/auto-place` so it can be called internally):

```python
def _spawn_playwright_placer(data: dict) -> dict:
    """Spawn the existing Playwright bet placer subprocess.
    Returns the same shape as the old /api/auto-place success/error."""
    bookmaker = data.get("bookmaker_key")
    if bookmaker not in PLAYWRIGHT_BOOKS:
        return {"success": False,
                "error": f"Playwright path not configured for {bookmaker}"}
    try:
        local_root = PROJECT_ROOT
        repo_root = _REPO_ROOT
        placer_script = str(local_root / "bet_placer" / "placer.py")
        bet_json = json.dumps(data)
        placer_python = str(repo_root / "wagerzon_odds" / "venv" / "bin" / "python3")
        if not Path(placer_python).exists():
            placer_python = sys.executable
        log_path = str(local_root / "bet_placer" / "placer.log")
        with open(log_path, "a") as log_file:
            subprocess.Popen(
                [placer_python, placer_script, bet_json],
                cwd=str(local_root / "bet_placer"),
                stdout=log_file, stderr=log_file,
                start_new_session=True,
            )
        return {"success": True, "status": "playwright_launched",
                "message": f"Browser launching for {bookmaker}..."}
    except Exception as e:
        return {"success": False, "error": str(e)}


def _insert_placement_breadcrumb(bet_hash: str, account: str,
                                  bet_meta: dict, status: str = "placing"):
    """Insert (or upsert) a placed_bets row with the given status. Used to
    create a breadcrumb before calling the WZ placer so we never lose a
    ticket to a mid-call crash."""
    con = duckdb.connect(str(DB_PATH))
    try:
        existing = con.execute(
            "SELECT bet_hash FROM placed_bets WHERE bet_hash = ?",
            [bet_hash]).fetchone()
        if existing is None:
            con.execute("""
                INSERT INTO placed_bets
                  (bet_hash, game_id, home_team, away_team, market, bet_on,
                   line, odds, actual_size, recommended_size, bookmaker,
                   account, status, wz_odds_at_place)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """, [
                bet_hash, bet_meta.get("game_id"),
                bet_meta.get("home_team"), bet_meta.get("away_team"),
                bet_meta.get("market"), bet_meta.get("bet_on"),
                bet_meta.get("line"), bet_meta.get("american_odds"),
                bet_meta.get("actual_size"), bet_meta.get("kelly_bet"),
                bet_meta.get("bookmaker_key"),
                account, status,
                bet_meta.get("wz_odds_at_place"),
            ])
        else:
            con.execute("""
                UPDATE placed_bets
                SET status = ?, account = ?, wz_odds_at_place = ?
                WHERE bet_hash = ?
            """, [status, account,
                  bet_meta.get("wz_odds_at_place"), bet_hash])
    finally:
        con.close()


def _finalize_placement(bet_hash: str, result: dict):
    """Upsert the final status from the placer's result onto the breadcrumb row."""
    con = duckdb.connect(str(DB_PATH))
    try:
        con.execute("""
            UPDATE placed_bets
            SET status = ?, ticket_number = ?,
                error_msg = ?, error_msg_key = ?
            WHERE bet_hash = ?
        """, [
            result.get("status"), result.get("ticket_number"),
            result.get("error_msg"), result.get("error_msg_key"),
            bet_hash,
        ])
    finally:
        con.close()
```

Then replace the body of `/api/place-bet` (after reading the existing handler — preserve the manual-log path for unsupported books):

```python
@app.route("/api/place-bet", methods=["POST"])
def place_bet():
    data = request.json or {}
    bookmaker_key = data.get("bookmaker_key")
    bet_hash      = data.get("bet_hash")
    account       = data.get("account")

    if not bet_hash:
        return jsonify({"success": False, "error": "bet_hash required"}), 400

    # WZ direct-API path
    if bookmaker_key in WZ_DIRECT_API_BOOKS:
        if not account:
            return jsonify({"success": False,
                            "error": "account required for wagerzon placement"}), 400
        _insert_placement_breadcrumb(bet_hash, account, data, status="placing")
        result = single_placer.place_single(account=account, bet=data)
        _finalize_placement(bet_hash, result)
        return jsonify(result)

    # Playwright path (Hoop88 / BFA / BetOnline)
    if bookmaker_key in PLAYWRIGHT_BOOKS:
        return jsonify(_spawn_playwright_placer(data))

    # Reference / sharp books — no placement integration
    return jsonify({"success": False,
                    "error": f"manual log only for '{bookmaker_key}'; "
                             f"use the [Log] button"}), 400
```

If there's a legacy "manual log" semantic in the old `/api/place-bet` (writing a row to placed_bets without placement), move that to a new `/api/log-bet` endpoint — the spec calls for a separate `[Log]` button that always logs without any placement attempt:

```python
@app.route("/api/log-bet", methods=["POST"])
def log_bet():
    """Manual log: record a placement without contacting any book."""
    data = request.json or {}
    bet_hash = data.get("bet_hash")
    if not bet_hash:
        return jsonify({"success": False, "error": "bet_hash required"}), 400
    _insert_placement_breadcrumb(bet_hash, data.get("account"),
                                  data, status="placed")
    return jsonify({"success": True, "status": "placed",
                    "ticket_number": None})
```

- [ ] **Step 5: Run tests to verify they pass**

```bash
cd "Answer Keys/MLB Dashboard"
python3 -m pytest tests/test_place_bet_dispatch.py -v
```

Expected: 4 passed.

- [ ] **Step 6: Commit**

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard_server.py" \
        "Answer Keys/MLB Dashboard/tests/test_place_bet_dispatch.py"
git commit -m "feat(mlb-dashboard): /api/place-bet dispatcher + /api/log-bet"
```

---

## Task 8: Dashboard loader for `mlb_bets_book_prices` (defensive)

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` (around line 4476 where `mlb_bets_combined` is loaded)

- [ ] **Step 1: Read the existing loader**

```bash
sed -n '4470,4500p' "Answer Keys/MLB Dashboard/mlb_dashboard.R"
```

You're looking for the `dbGetQuery(con, "SELECT * FROM mlb_bets_combined")` call.

- [ ] **Step 2: Add the new table load alongside it**

Right after the existing `mlb_bets_combined` load, insert:

```r
book_prices_long <- tryCatch({
  dbGetQuery(con, "SELECT * FROM mlb_bets_book_prices")
}, error = function(e) {
  warning(sprintf(
    "[bets-tab] mlb_bets_book_prices not loaded (%s) — falling back to single-book layout",
    conditionMessage(e)))
  NULL
})
```

(`book_prices_long` will be either a populated frame or `NULL`. Downstream renderers check `is.null(book_prices_long)` and degrade.)

- [ ] **Step 3: Source the new pill helper**

Near the top of `mlb_dashboard.R` where the file sources `books_strip.R`, also source the new helper:

```r
source("book_pill.R")
```

- [ ] **Step 4: Manual verify by running the dashboard once**

```bash
cd "Answer Keys/MLB Dashboard"
python3 mlb_dashboard_server.py &
# Wait 3 seconds, then check
sleep 3
curl -sf http://localhost:8083/ > /dev/null && echo "loaded OK" || echo "FAILED"
# Stop the server
lsof -ti:8083 | xargs kill 2>/dev/null
```

Expected: `loaded OK`. (No regression — the dashboard renders today's table since the new card layout isn't wired yet.)

- [ ] **Step 5: Commit**

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "feat(mlb-dashboard): defensive loader for mlb_bets_book_prices"
```

---

## Task 9: Pivot book_prices long→wide for the renderer

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` (near the new loader from Task 8)

- [ ] **Step 1: Write the pivot helper**

Add to `mlb_dashboard.R`, right after the `book_prices_long <- tryCatch(...)` block:

```r
# Pivot long -> wide so each bet has a single row with per-book columns.
# Output columns per side:
#   <book>_line_quoted, <book>_american_odds, <book>_is_exact_line
# for book in (wagerzon, hoop88, bfa, bookmaker, bet105,
#              draftkings, fanduel, pinnacle).
# Two rows per bet (one with side='pick', one with side='opposite').
pivot_book_prices_wide <- function(long_frame) {
  if (is.null(long_frame) || nrow(long_frame) == 0) return(NULL)
  long_frame %>%
    select(bet_row_id, side, bookmaker, line_quoted, american_odds, is_exact_line) %>%
    tidyr::pivot_wider(
      names_from  = bookmaker,
      values_from = c(line_quoted, american_odds, is_exact_line),
      names_glue  = "{bookmaker}_{.value}",
      values_fill = list(line_quoted = NA_real_,
                         american_odds = NA_integer_,
                         is_exact_line = NA)
    )
}

book_prices_wide <- pivot_book_prices_wide(book_prices_long)
```

- [ ] **Step 2: Spot-check the pivot output**

Add a temporary `cat()` after the pivot to inspect:

```r
if (!is.null(book_prices_wide)) {
  cat("[bets-tab] pivoted:", nrow(book_prices_wide), "rows,",
      ncol(book_prices_wide), "cols\n")
  cat("[bets-tab] columns:", paste(colnames(book_prices_wide), collapse = ", "), "\n")
}
```

Run the dashboard once and check the stdout. Expected: rows = 2 × bet count, columns include `wagerzon_american_odds`, `hoop88_line_quoted`, `bookmaker_is_exact_line`, etc.

Remove the temporary `cat()` lines after verification.

- [ ] **Step 3: Commit**

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "feat(mlb-dashboard): pivot mlb_bets_book_prices long->wide for card renderer"
```

---

## Task 10: CSS for cards + pills

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` (the inline `<style>` block inside `create_report`)

- [ ] **Step 1: Locate the inline style block**

```bash
grep -n "parlays-table-container" "Answer Keys/MLB Dashboard/mlb_dashboard.R" | head -5
```

The CSS lives in `create_report` around lines 1879-2070 (parlays-tab card styles).

- [ ] **Step 2: Add the bets-tab card classes**

Add this block after the existing parlays CSS, inside the same `<style>` element:

```css
/* === Bets tab odds-screen card layout === */

/* Each bet renders as a card via reactable row + display:flex on the row container. */
#bets-table-container .rt-table   { display: block; }
#bets-table-container .rt-thead   { display: none; }
#bets-table-container .rt-tbody   { display: block; }
#bets-table-container .rt-tr-group { display: block; }
#bets-table-container .rt-tr {
  display: flex;
  flex-wrap: wrap;
  background: #1c2128;
  border: 1px solid #373e47;
  border-radius: 10px;
  padding: 12px 16px;
  margin-bottom: 10px;
  align-items: center;
}
#bets-table-container .rt-td {
  min-width: 0 !important;
  flex: 0 1 auto;
}

#bets-table-container .cell-game,
#bets-table-container .cell-market,
#bets-table-container .cell-pickside,
#bets-table-container .cell-otherside {
  flex-basis: 100%;
  width: 100%;
  padding: 2px 0;
}

#bets-table-container .cell-game {
  order: 1;
  color: #c9d1d9; font-weight: 600;
}
#bets-table-container .cell-market {
  order: 2;
  color: #c9d1d9; font-weight: 600;
  border-bottom: 1px solid #2d333b;
  padding-bottom: 6px;
  margin-bottom: 6px;
}
#bets-table-container .cell-pickside  { order: 3; }
#bets-table-container .cell-otherside { order: 4; }

/* Pill row: each row is a flex container of fixed-width pills */
#bets-table-container .side-row {
  display: flex; flex-wrap: wrap;
  gap: 6px; margin: 4px 0;
  align-items: center;
}
#bets-table-container .side-label {
  flex: 0 0 78px; min-width: 78px;
  color: #c9d1d9; font-weight: 600; font-size: 12px;
}
#bets-table-container .pill {
  flex: 0 0 78px; min-width: 78px;
  box-sizing: border-box;
  display: inline-flex; flex-direction: column;
  align-items: center; gap: 1px;
  padding: 5px 6px;
  border-radius: 6px;
  background: #22272e;
  border: 1px solid #373e47;
  color: #c9d1d9;
  font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
  font-size: 12px;
}
#bets-table-container .pill .book {
  color: #8b949e; font-size: 9px;
  text-transform: uppercase; letter-spacing: 0.5px;
}
#bets-table-container .pill .line-tag {
  color: #f0883e; font-size: 9px; font-weight: 600;
}
#bets-table-container .pill.muted {
  color: #6e7681;
  background: transparent;
  border-style: dashed;
}
#bets-table-container .pill.mismatched {
  border-color: #5a3d1a;
  background: #2a2317;
}
#bets-table-container .pill.pick {
  background: #15321f;
  border-color: #3fb950;
  color: #56d364;
  font-weight: 700;
}
#bets-table-container .pill.pick .book { color: #56d364; }

/* Metadata strip */
#bets-table-container .cell-m,
#bets-table-container .cell-pick,
#bets-table-container .cell-ev,
#bets-table-container .cell-size,
#bets-table-container .cell-towin,
#bets-table-container .cell-action {
  order: 5;
  display: inline-flex; align-items: center;
  gap: 5px;
  padding: 0 18px 0 0;
  font-size: 14px;
}
#bets-table-container .cell-m::before    { content: "M";       color: #8b949e; font-size: 11px; text-transform: uppercase; letter-spacing: 0.5px; margin-right: 4px; }
#bets-table-container .cell-pick::before { content: "Pick";    color: #8b949e; font-size: 11px; text-transform: uppercase; letter-spacing: 0.5px; margin-right: 4px; }
#bets-table-container .cell-ev::before   { content: "EV";      color: #8b949e; font-size: 11px; text-transform: uppercase; letter-spacing: 0.5px; margin-right: 4px; }
#bets-table-container .cell-size::before { content: "Size";    color: #8b949e; font-size: 11px; text-transform: uppercase; letter-spacing: 0.5px; margin-right: 4px; }
#bets-table-container .cell-towin::before{ content: "To Win";  color: #8b949e; font-size: 11px; text-transform: uppercase; letter-spacing: 0.5px; margin-right: 4px; }
#bets-table-container .cell-action       { margin-left: auto; padding-right: 0; }
#bets-table-container .meta-strip-row {
  width: 100%;
  display: flex; align-items: center; flex-wrap: wrap; gap: 14px;
  margin-top: 10px;
  padding-top: 10px;
  border-top: 1px solid #2d333b;
}

/* Correlation dot — corner badge on the card */
#bets-table-container .corr-badge {
  position: absolute;
  top: 8px; right: 12px;
  color: #f0883e; font-size: 18px;
  cursor: help;
}
#bets-table-container .rt-tr { position: relative; }
```

- [ ] **Step 3: Verify the dashboard still loads**

```bash
cd "Answer Keys/MLB Dashboard"
python3 mlb_dashboard_server.py &
sleep 3
curl -sf http://localhost:8083/ > /dev/null && echo "loaded OK" || echo "FAILED"
lsof -ti:8083 | xargs kill 2>/dev/null
```

Expected: `loaded OK`. (Visual change isn't visible yet — the bets table still uses the old reactable; CSS is dormant until Task 11 wires the cell classes.)

- [ ] **Step 4: Commit**

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "feat(mlb-dashboard): CSS for bets-tab card layout (pills, meta strip, corr dot)"
```

---

## Task 11: Replace `create_bets_table()` body with card layout

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` (function at line 817)

This is the largest single change. The function builds a reactable; we're changing the column definitions to render the card layout via the CSS from Task 10 and the helper from Task 1.

- [ ] **Step 1: Read the current function (200 lines)**

```bash
sed -n '817,1120p' "Answer Keys/MLB Dashboard/mlb_dashboard.R"
```

Note the existing data attribute scheme on the Place button (`data-hash`, `data-game-id`, etc.) — we'll keep all of those.

- [ ] **Step 2: Replace the function**

Replace the entire `create_bets_table()` function with the version below. Read each segment carefully and adapt variable names to whatever `all_bets` carries in your environment.

```r
create_bets_table <- function(all_bets, placed_bets, book_prices_wide = NULL) {
  placed_hashes <- if (nrow(placed_bets) > 0) placed_bets$bet_hash else character()

  placed_status_lookup <- if (nrow(placed_bets) > 0 && "status" %in% names(placed_bets)) {
    setNames(placed_bets$status, placed_bets$bet_hash)
  } else setNames(character(), character())
  placed_ticket_lookup <- if (nrow(placed_bets) > 0 && "ticket_number" %in% names(placed_bets)) {
    setNames(placed_bets$ticket_number, placed_bets$bet_hash)
  } else setNames(character(), character())

  same_game_info <- lapply(seq_len(nrow(all_bets)), function(i) {
    find_same_game_bets(i, all_bets, placed_bets)
  })

  # Build the long-format display frame: two rows per bet (pick + opposite).
  if (is.null(book_prices_wide)) {
    # Fall back to today's flat table if the new table isn't available.
    warning("[bets-tab] book_prices_wide is NULL — falling back to legacy renderer")
    return(create_bets_table_legacy(all_bets, placed_bets))
  }

  # Compute bet_row_id on the bets frame the same way MLB.R does, so the join
  # matches mlb_bets_book_prices.
  all_bets <- all_bets %>%
    mutate(bet_row_id = vapply(
      paste(id, market, ifelse(is.na(line), "", line), bet_on, sep = "|"),
      digest::digest, character(1), algo = "md5"
    ))

  # Wide frame has 2 rows per bet_row_id (side = pick / opposite).
  card_rows <- all_bets %>%
    inner_join(book_prices_wide, by = "bet_row_id") %>%
    arrange(desc(ev), bet_row_id, desc(side))  # pick row above opposite

  # Helper: render a side row's pills given the row data
  BOOK_ORDER <- c("wagerzon", "hoop88", "bfa", "bookmaker", "bet105",
                  "draftkings", "fanduel", "pinnacle")
  BOOK_LABELS <- c(wagerzon = "WZ", hoop88 = "H88", bfa = "BFA",
                   bookmaker = "BKM", bet105 = "B105",
                   draftkings = "DK", fanduel = "FD", pinnacle = "Pinn")

  render_side_row <- function(row, model_line, side_label_text, is_pick_side, pick_book) {
    pills <- vapply(BOOK_ORDER, function(b) {
      odds <- row[[paste0(b, "_american_odds")]]
      lq   <- row[[paste0(b, "_line_quoted")]]
      render_book_pill(
        book          = BOOK_LABELS[[b]],
        american_odds = if (is.na(odds)) NA_integer_ else as.integer(odds),
        line_quoted   = lq,
        model_line    = model_line,
        is_pick       = is_pick_side && (b == pick_book),
        side          = side_label_text
      )
    }, character(1))
    paste0(
      '<div class="side-row">',
      sprintf('<span class="side-label">%s</span>', htmltools::htmlEscape(side_label_text)),
      paste(pills, collapse = ""),
      '</div>'
    )
  }

  # Build a single string of HTML per bet representing the whole card body.
  # reactable renders each row as a flex container; we package the side rows
  # under cell-pickside / cell-otherside, and use cell-meta-strip-row for
  # the bottom metadata.
  card_data <- all_bets %>%
    mutate(
      ev_display = ifelse(ev >= 0, sprintf("+%.1f%%", ev * 100), sprintf("%.1f%%", ev * 100)),
      m_display  = sprintf("%.1f%%", prob * 100),
      size_display = sprintf("$%.0f", bet_size),
      towin_display = sprintf("$%.0f",
                              ifelse(odds > 0, bet_size * odds / 100,
                                              bet_size * 100 / abs(odds))),
      bet_hash = pmap_chr(list(id, market, bet_on, line), generate_bet_hash),
      is_placed = bet_hash %in% placed_hashes,
      game_display = paste(away_team, "@", home_team),
      pick_display = sprintf("%s %s",
                             toupper(bookmaker_key),
                             ifelse(odds > 0, paste0("+", odds), as.character(odds))),
      market_display = sprintf("%s %s", format_market_name(market),
                               ifelse(is.na(line), bet_on,
                                      paste(bet_on, ifelse(line > 0, paste0("+", line), line))))
    )

  # The reactable still expects one row per BET (not per side). We pre-render
  # both side rows as HTML strings stored in cell-pickside / cell-otherside.
  card_data$pickside_html <- vapply(seq_len(nrow(card_data)), function(i) {
    bet_id <- card_data$bet_row_id[i]
    wide_pick <- book_prices_wide %>% filter(bet_row_id == bet_id, side == "pick")
    if (nrow(wide_pick) == 0) return("")
    render_side_row(wide_pick[1, ], card_data$line[i],
                    side_label_text = card_data$bet_on[i],
                    is_pick_side = TRUE,
                    pick_book = card_data$bookmaker_key[i])
  }, character(1))

  card_data$otherside_html <- vapply(seq_len(nrow(card_data)), function(i) {
    bet_id <- card_data$bet_row_id[i]
    wide_opp <- book_prices_wide %>% filter(bet_row_id == bet_id, side == "opposite")
    if (nrow(wide_opp) == 0) return("")
    opposite_label <- if (grepl("^Over", card_data$bet_on[i])) "Under"
                      else if (grepl("^Under", card_data$bet_on[i])) "Over"
                      else "Opp"
    render_side_row(wide_opp[1, ], card_data$line[i],
                    side_label_text = opposite_label,
                    is_pick_side = FALSE,
                    pick_book = card_data$bookmaker_key[i])
  }, character(1))

  reactable(
    card_data,
    elementId = "bets-table",
    searchable = TRUE,
    filterable = TRUE,
    defaultPageSize = 25,
    columns = list(
      bet_hash      = colDef(show = FALSE),
      bet_row_id    = colDef(show = FALSE),
      id            = colDef(show = FALSE),
      home_team     = colDef(show = FALSE),
      away_team     = colDef(show = FALSE),
      market        = colDef(show = FALSE),
      line          = colDef(show = FALSE),
      bet_on        = colDef(show = FALSE),
      bet_size      = colDef(show = FALSE),
      odds          = colDef(show = FALSE),
      prob          = colDef(show = FALSE),
      ev            = colDef(show = FALSE),
      bookmaker_key = colDef(show = FALSE),
      pt_start_time = colDef(show = FALSE),
      market_type   = colDef(show = FALSE),
      correlation_adj = colDef(show = FALSE),

      game_display = colDef(
        name = "", html = TRUE, class = "cell-game",
        cell = function(value, index) {
          row <- card_data[index, ]
          time_str <- if (!is.na(row$pt_start_time))
            format(row$pt_start_time, "%a %I:%M %p") else ""
          corr_badge <- ""
          if (same_game_info[[index]]$has_same_game) {
            corr_badge <- sprintf(
              '<span class="corr-badge" title="%s">&#9679;</span>',
              htmltools::htmlEscape(paste("Other bets on this game",
                                         collapse = "; ")))
          }
          sprintf('%s<span style="color:#8b949e;font-size:11px;margin-left:8px">%s</span>%s',
                  htmltools::htmlEscape(value), time_str, corr_badge)
        }
      ),
      market_display = colDef(
        name = "", html = TRUE, class = "cell-market",
        cell = function(value, index) htmltools::htmlEscape(value)
      ),
      pickside_html = colDef(
        name = "", html = TRUE, class = "cell-pickside",
        cell = function(value, index) value
      ),
      otherside_html = colDef(
        name = "", html = TRUE, class = "cell-otherside",
        cell = function(value, index) value
      ),
      m_display = colDef(
        name = "", class = "cell-m",
        style = list(color = "#7ee787", fontWeight = 600)
      ),
      pick_display = colDef(
        name = "", class = "cell-pick"
      ),
      ev_display = colDef(
        name = "", class = "cell-ev",
        cell = function(value, index) {
          ep <- card_data$ev[index] * 100
          color <- if (ep >= 15) "#3fb950" else if (ep >= 10) "#56d364"
                   else if (ep >= 5) "#7ee787" else "#a5d6a7"
          htmltools::div(style = list(color = color, fontWeight = 600), value)
        }
      ),
      size_display = colDef(
        name = "", class = "cell-size"
      ),
      towin_display = colDef(
        name = "", class = "cell-towin",
        style = list(color = "#3fb950")
      ),
      is_placed = colDef(
        name = "", class = "cell-action", html = TRUE,
        cell = function(value, index) {
          row <- card_data[index, ]
          status <- placed_status_lookup[row$bet_hash]
          ticket <- placed_ticket_lookup[row$bet_hash]
          data_attrs <- sprintf(
            'data-hash="%s" data-game-id="%s" data-home="%s" data-away="%s" data-market="%s" data-bet-on="%s" data-line="%s" data-prob="%s" data-ev="%s" data-size="%s" data-odds="%s" data-book="%s" data-account-required="%s"',
            row$bet_hash, row$id, row$home_team, row$away_team,
            row$market, row$bet_on,
            ifelse(is.na(row$line), "", row$line),
            row$prob, row$ev, row$bet_size, row$odds,
            row$bookmaker_key,
            ifelse(row$bookmaker_key == "wagerzon", "true", "false")
          )
          if (!is.na(status) && status == "placed") {
            label <- if (!is.na(ticket) && nchar(ticket) > 0)
                       sprintf('placed &middot; #%s', htmltools::htmlEscape(ticket))
                     else "placed"
            return(sprintf('<span class="placed-bet-label" %s>%s</span>',
                          data_attrs, label))
          }
          if (!is.na(status) && status %in% c("price_moved", "rejected",
                                              "auth_error", "network_error",
                                              "orphaned")) {
            short <- switch(status,
              price_moved = "drift", rejected = "rejected",
              auth_error = "auth err", network_error = "net err",
              orphaned = "orphan", status)
            return(sprintf('<span class="pill error" %s>%s</span><button class="btn-place" onclick="placeBet(this)" %s>Retry</button><button class="btn-log" onclick="logBet(this)" %s>Log</button>',
                          data_attrs, short, data_attrs, data_attrs))
          }
          # Not placed yet
          place_disabled <- !(row$bookmaker_key %in% c("wagerzon", "hoop88",
                                                       "bfa", "betonlineag"))
          place_btn <- if (place_disabled) {
            sprintf('<button class="btn-place" disabled title="manual log only for this book" %s>Place</button>', data_attrs)
          } else {
            sprintf('<button class="btn-place" onclick="placeBet(this)" %s>Place</button>', data_attrs)
          }
          log_btn <- sprintf('<button class="btn-log" onclick="logBet(this)" %s>Log</button>', data_attrs)
          paste0(place_btn, " ", log_btn)
        }
      )
    ),
    theme = reactableTheme(
      backgroundColor = "transparent",
      borderColor = "transparent",
      stripedColor = "transparent",
      highlightColor = "transparent"
    )
  )
}

# Legacy fallback: today's flat table renderer, kept verbatim. Called when
# mlb_bets_book_prices is unavailable (defensive). Move the OLD function body
# (lines 817-1126 in the prior version of this file) into this function.
create_bets_table_legacy <- function(all_bets, placed_bets) {
  # [paste the existing create_bets_table function body here, unchanged]
  stop("create_bets_table_legacy not yet wired — copy the prior implementation")
}
```

(Read the existing `create_bets_table` body once and paste it into `create_bets_table_legacy()` — the fallback only fires if `mlb_bets_book_prices` is missing.)

- [ ] **Step 3: Update the caller to pass `book_prices_wide`**

Find the call site (where `create_report` invokes `create_bets_table`):

```bash
grep -n "create_bets_table(" "Answer Keys/MLB Dashboard/mlb_dashboard.R" | head -3
```

Change the call to pass the wide frame:

```r
bets_table <- create_bets_table(all_bets, placed_bets, book_prices_wide = book_prices_wide)
```

- [ ] **Step 4: Manual browser verify**

```bash
cd "Answer Keys/MLB Dashboard"
./run.sh &
sleep 5
open http://localhost:8083
```

Verify in browser:
- Each card shows game header + market header + two side-pill rows + metadata strip.
- Pick book pill is green-tinted + bold.
- Some pills are amber (mismatched line) — verify the line tag.
- Some pills are dashed `—` (no quote).
- `[Place]` button shows for WZ/H88/BFA picks; disabled for DK/FD/Pinn/BKM/B105 picks; tooltip explains why.
- `[Log]` button shows on every card.
- Correlation dot appears on cards where you have placed bets on the same game.

Stop the server:
```bash
lsof -ti:8083 | xargs kill 2>/dev/null
```

- [ ] **Step 5: Commit**

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "feat(mlb-dashboard): card-layout bets tab with pill rows + meta strip"
```

---

## Task 12: JS wiring — `placeBet()` dispatcher + `logBet()`

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` (inline `<script>` block inside `create_report`)

- [ ] **Step 1: Locate the existing JS block**

```bash
grep -n "function placeBet\|function autoPlaceBet\|function placeParlay" "Answer Keys/MLB Dashboard/mlb_dashboard.R" | head -10
```

The inline `<script>` block lives in `create_report`. Find where `placeBet`, `autoPlaceBet`, etc. are defined.

- [ ] **Step 2: Replace `placeBet()` with the dispatcher**

In the JS block, replace the existing `placeBet` and `autoPlaceBet` functions with:

```javascript
function placeBet(btn) {
  const data = btn.dataset;
  const book = data.book;
  const account = window._wzActiveAccount || null;

  // Build the request body. Keys match what the server's
  // _insert_placement_breadcrumb expects.
  const body = {
    bet_hash:       data.hash,
    bookmaker_key:  book,
    account:        account,
    bet_on:         data.betOn,
    line:           parseFloat(data.line || 'NaN') || null,
    market:         data.market,
    american_odds:  parseInt(data.odds),
    actual_size:    parseFloat(data.size),
    kelly_bet:      parseFloat(data.size),
    wz_odds_at_place: parseInt(data.odds),
    game_id:        data.gameId,
    home_team:      data.home,
    away_team:      data.away,
    // For WZ direct API, idgm + play must be passed through. They aren't
    // in the bet_hash; the server can look them up by bet_hash from the
    // pipeline-side bet frame (mlb_bets_combined or a sibling).
    // For now they're absent here; the server fills them in.
  };

  if (book === 'wagerzon' && !account) {
    showToast('No Wagerzon account selected — pick one in the header pills', 'error');
    return;
  }

  btn.disabled = true;
  btn.textContent = 'Placing...';

  fetch('/api/place-bet', {
    method: 'POST',
    headers: {'Content-Type': 'application/json'},
    body: JSON.stringify(body)
  })
  .then(r => r.json())
  .then(result => {
    if (result.status === 'placed') {
      showToast(`Placed at ${book} (#${result.ticket_number || 'logged'})`, 'success');
      setTimeout(() => location.reload(), 800);
    } else if (result.status === 'playwright_launched') {
      showToast(result.message, 'info');
    } else if (result.status === 'price_moved') {
      showToast('Price moved — bet not placed', 'warning');
      btn.disabled = false; btn.textContent = 'Place';
    } else if (result.error) {
      showToast(result.error, 'error');
      btn.disabled = false; btn.textContent = 'Place';
    } else {
      showToast(`Status: ${result.status}`, 'error');
      btn.disabled = false; btn.textContent = 'Place';
    }
  })
  .catch(e => {
    showToast('Network error: ' + e.message, 'error');
    btn.disabled = false; btn.textContent = 'Place';
  });
}

function logBet(btn) {
  const data = btn.dataset;
  fetch('/api/log-bet', {
    method: 'POST',
    headers: {'Content-Type': 'application/json'},
    body: JSON.stringify({
      bet_hash: data.hash,
      account: window._wzActiveAccount,
      bookmaker_key: data.book,
      bet_on: data.betOn, line: parseFloat(data.line || 'NaN') || null,
      market: data.market,
      american_odds: parseInt(data.odds),
      actual_size: parseFloat(data.size),
      kelly_bet: parseFloat(data.size),
      game_id: data.gameId, home_team: data.home, away_team: data.away,
    })
  })
  .then(r => r.json())
  .then(() => location.reload())
  .catch(e => showToast('Log failed: ' + e.message, 'error'));
}
```

(If `showToast` doesn't exist in the current JS, swap it for `alert(...)` calls; the parlay tab has a real toast helper — re-use whatever pattern is already there.)

- [ ] **Step 3: Manual verify with a fake WZ bet**

Start the dashboard:
```bash
cd "Answer Keys/MLB Dashboard"
./run.sh &
sleep 5
```

Open browser → bets tab → find a WZ-pick card. Click `[Place]`. Verify:

1. **Dry-run-style first verification (no real money):** use a tiny bet by editing the slider to size=$1, OR set a developer flag in single_placer that returns a fake-success without touching WZ.
2. Watch the network tab — POST to `/api/place-bet`, response contains a status.
3. If status is `placed`, card flips to `placed · #<ticket>`.
4. If status is `price_moved` or `rejected`, toast appears with the message.

Then click `[Log]` on a different card. Verify it instantly flips to `[Placed]` without any network call to WZ.

Stop the server:
```bash
lsof -ti:8083 | xargs kill 2>/dev/null
```

- [ ] **Step 4: Commit**

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "feat(mlb-dashboard): placeBet() dispatcher + logBet() JS hooks"
```

---

## Task 13: End-to-end verification + documentation

**Files:**
- Modify: `Answer Keys/MLB Dashboard/README.md`
- Modify: `Answer Keys/CLAUDE.md`

- [ ] **Step 1: Run the full pipeline end-to-end**

```bash
cd "Answer Keys/MLB Dashboard"
./run.sh
```

Wait ~10 minutes. Watch for errors in the output. Expected: pipeline completes, dashboard starts on :8083.

- [ ] **Step 2: Run the verification queries from the spec**

```bash
DB="/Users/callancapitolo/NFLWork/Answer Keys/mlb_mm.duckdb"

duckdb "$DB" "DESCRIBE mlb_bets_book_prices"
duckdb "$DB" "SELECT COUNT(*) AS rows, COUNT(DISTINCT bet_row_id) AS bets, COUNT(DISTINCT bookmaker) AS books FROM mlb_bets_book_prices"
duckdb "$DB" "SELECT market, period, is_exact_line, COUNT(*) AS rows FROM mlb_bets_book_prices GROUP BY 1,2,3 ORDER BY 1,2,3"
duckdb "$DB" "SELECT COUNT(*) FROM mlb_bets_book_prices WHERE is_exact_line = FALSE AND ABS(line_quoted - line) > 1.0"
duckdb "$DB" "SELECT market, period, SUM(CASE WHEN bookmaker IN ('draftkings','fanduel','pinnacle') THEN 1 ELSE 0 END) AS ref_rows, COUNT(*) AS total FROM mlb_bets_book_prices GROUP BY 1,2 ORDER BY 1,2"
```

Expected:
- DESCRIBE: 11 columns
- COUNT: positive rows / bets / books
- is_exact_line distribution: most F5 mains TRUE, F3/F7 should show some FALSE
- Threshold guard: zero rows
- Reference rows: non-zero for FG/F3/F5/F7 mains + alts; zero for team_totals + odd_even

- [ ] **Step 3: Browser smoke test**

Open http://localhost:8083 → bets tab. Verify:
- Cards render with the locked layout (game header, market header, two pill rows aligned, meta strip).
- Pick book pill is green-tinted.
- Mismatched-line pills are amber with O5/U6 line tags.
- `—` pills appear on at least some cards (no quote within ±1).
- Click `[Place]` on a WZ card with size=$1 — expect placed pill with ticket #.
- Click `[Place]` on an H88 card — expect Playwright window to open.
- Click `[Log]` on a DK card — expect `[Placed]` flip.
- Switch WZ accounts via header pills — verify the next Place uses the new account.

- [ ] **Step 4: Update README**

Path: `Answer Keys/MLB Dashboard/README.md`

Locate the bets-tab section and replace with:

```markdown
- **Bets tab — odds-screen card layout** — Each opportunity renders as a
  card mirroring the parlays tab. Card body shows two pill rows (pick
  side + opposite side), each row laid out as one fixed-width pill per
  book in a stable left-to-right order: WZ, H88, BFA, BKM, B105, DK, FD,
  Pinn. Three pill states:
    - exact-line price (default)
    - mismatched within ±1 unit (amber background + amber line tag like
      `O5` over the price)
    - no quote (muted dashed `—`)
  Pick book pill is green-tinted + bold. Metadata strip at the bottom:
  M (model fair %) / Pick / EV / Size / To Win / [Place] [Log].
- **`[Place]` dispatches by book:**
    - Wagerzon → direct REST API (no browser, returns ticket #), reusing
      the parlay-placer's session/preflight/drift/orphan code paths
    - Hoop88 / BFA / BetOnlineAG → existing Playwright browser flow
    - DraftKings / FanDuel / Pinnacle / Bookmaker / Bet105 → button
      disabled; use `[Log]` to record placement made externally
- **`[Log]` always records a manual placement** without contacting any
  book. Mirrors the parlays tab's Log button.
```

Update the Data Storage section to mention the new table:

```markdown
- Pipeline bets read from `Answer Keys/mlb_mm.duckdb`:
  `mlb_bets_combined`, `mlb_bets_book_prices` (new — per-book pill
  rows for the bets-tab odds screen), `mlb_parlay_opportunities`,
  `mlb_trifecta_opportunities`.
```

Update the API Endpoints table to mention `/api/log-bet`:

```markdown
| `/api/log-bet` | POST | Manual placement log (no book contact) |
```

- [ ] **Step 5: Update `Answer Keys/CLAUDE.md`**

In the "DuckDB Databases" section, update the `mlb_mm.duckdb` line:

```markdown
- `mlb_mm.duckdb` — MLB MM consumer tables (`mlb_bets_combined`,
  `mlb_bets_book_prices` [new 2026-05-11 — per-book pill rows for the
  bets-tab odds screen], `mlb_game_samples`, `mlb_samples_meta`,
  `mlb_sgp_odds`, `mlb_parlay_lines`, `mlb_parlay_opportunities`,
  `mlb_trifecta_opportunities`). ...
```

Note the schema change to `placed_bets` (after the existing Wagerzon
multi-account schema notes):

```markdown
### `placed_bets` schema additions (2026-05-11)

- `placed_bets.account` — WZ account label; NULL for non-WZ placements
- `placed_bets.status` — placing / placed / price_moved / rejected /
  auth_error / network_error / orphaned
- `placed_bets.ticket_number` — WZ ticket for direct-API placements
- `placed_bets.error_msg`, `error_msg_key` — error rendering
- `placed_bets.wz_odds_at_place` — drift forensics

Migration: `Answer Keys/MLB Dashboard/migrations/002_single_placer_columns.py`
```

- [ ] **Step 6: Commit docs**

```bash
git add "Answer Keys/MLB Dashboard/README.md" "Answer Keys/CLAUDE.md"
git commit -m "docs(mlb-dashboard): document odds-screen bets tab + WZ singles placer"
```

---

## Pre-merge

Run the pre-merge checklist from the spec
(`Answer Keys/MLB Dashboard/PLAN_odds_screen.md` § Pre-merge executive
review checklist). All items in three categories (Pipeline / Dashboard /
Auto-placer) plus General must pass before merging.

```bash
git diff main..HEAD --stat
```

If everything passes and the user explicitly approves:

```bash
# from main
cd /Users/callancapitolo/NFLWork
git checkout main
git merge --no-ff claude/mlb-sportsbook-comparison-98fM8
git push origin main

# Cleanup
git worktree remove .worktrees/mlb-odds-screen-rev
git branch -d claude/mlb-sportsbook-comparison-98fM8
git push origin --delete claude/mlb-sportsbook-comparison-98fM8
```

If any task's commit failed pre-merge review, revisit the relevant task,
fix in place, append a fixup commit, and re-run the checklist.

---

## Notes for the implementer

- **Sandbox / paths:** This plan assumes work happens in the worktree at
  `/Users/callancapitolo/NFLWork/.worktrees/mlb-odds-screen-rev/`. All
  paths above are repo-relative.
- **Pipeline runtime:** Each `Rscript MLB.R` run is 3-8 minutes. Don't
  iterate on the pipeline task by re-running the full pipeline — write
  TDD tests against fixtures and only run the full pipeline at the end
  of each pipeline-modifying task.
- **DuckDB locks:** `mlb.duckdb` is held under a long write lock during
  pipeline runs. The new table lives in `mlb_mm.duckdb` to avoid this.
  Don't move it.
- **Wagerzon live testing:** Test the WZ single placer first with a $1
  bet on a real game. Verify the ticket appears in `placed_bets` AND in
  Wagerzon's web ticket history before scaling up.
- **Fallback path matters:** the `create_bets_table_legacy()` function
  must contain a complete copy of the prior `create_bets_table` body.
  Without it, a deploy with a missing `mlb_bets_book_prices` table will
  fail to render the bets tab.
