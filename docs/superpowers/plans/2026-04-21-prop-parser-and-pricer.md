# Generic Prop Parser + Pricer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the per-prop-type fair-probability functions in `Answer Keys/mlb_triple_play.R` with a generic description-driven parser and evaluator, so adding a new prop type (grand slam, future 4-leggers) requires only a token-registry entry — not a new function. Also ship `wagerzon_odds/recon_specials.py`, ready to run tomorrow when new MLB specials post, to capture the Wagerzon page structure that a follow-up plan will turn into a production scraper.

**Architecture:** Add `Answer Keys/parse_legs.R` containing `parse_legs()`, `eval_leg()`, `compute_prop_fair()`, and a token registry that maps Wagerzon description tokens (`SCR 1ST`, `1H`, `F3`, `F7`, `GM`, `SCR U<N>`, `SCR O<N>`) to structured leg specs. Modify `mlb_triple_play.R` to source the parser, extend the samples query to include full-game + F3/F7 columns, and pipe each tribble row through `parse_legs(description) -> compute_prop_fair(samples, side, legs)`. The tribble keeps its current team/odds fields and adds a single `description` column — so the user still types the Wagerzon label verbatim once per day, but the parser now drives all leg logic. `compute_triple_play_fair()` is removed (dead after rewire). Recon script is a one-off Playwright artifact with no runtime dependencies.

**Tech Stack:** R (testthat, dplyr, data.table, DBI, duckdb), Python 3 (Playwright) for recon. No new dependencies.

---

## File Structure

**Created:**
- `Answer Keys/parse_legs.R` — pure functions: `parse_fraction()`, `parse_legs()`, `eval_leg()`, `compute_prop_fair()`, and a top-level `TOKEN_REGISTRY` list. One file so a pricer can `source()` one thing.
- `Answer Keys/tests/test_parse_legs.R` — testthat file exercising the parser + evaluator + generic pricer on synthetic inputs.
- `wagerzon_odds/recon_specials.py` — Playwright recon script (sibling of `recon_wagerzon.py`), captures network traffic from the Wagerzon MLB-specials page and dumps to `recon_specials.json`. Manually run, not part of any pipeline.

**Modified:**
- `Answer Keys/mlb_triple_play.R` — remove `compute_triple_play_fair()`; `source("parse_legs.R")`; extend the sample query to include `total_final_score, home_margin_f3, home_margin_f7`; add a `description` column to `todays_lines`; rewire the main-block `mutate` to call `compute_prop_fair(game_samples, side, parse_legs(description))`; add `prop_type` to the display table (parsed from description).
- `Answer Keys/tests/test_triple_play.R` — remove the now-dead tests for `compute_triple_play_fair`; keep helper tests (`determine_home_scored_first*`) untouched; keep `prob_to_american` / `american_to_prob` tests (they're still used).
- `Answer Keys/CLAUDE.md` — extend the "Triple-Play Data Flow" subsection with the new parser layer.

---

## Worktree & Version Control Plan

- Branch: `feature/prop-parser-pricer`
- Worktree: `.worktrees/prop-parser-pricer` (created before any code changes)
  ```bash
  cd /Users/callancapitolo/NFLWork
  git worktree add .worktrees/prop-parser-pricer -b feature/prop-parser-pricer main
  cd .worktrees/prop-parser-pricer
  git branch  # confirm feature/prop-parser-pricer
  ```
- DuckDB files: **not needed** for this plan — all tests use synthetic in-memory frames. The integration smoke test at the end of Task 9 runs against `main`'s `mlb.duckdb` in read-only mode AFTER merge.
- Commits (one per task completion):
  1. Task 1 → `feat(mlb): add parse_fraction helper for Wagerzon ½ notation`
  2. Task 2 → `feat(mlb): add parse_legs token registry for triple-play + grand slam`
  3. Task 3 → `feat(mlb): add eval_leg evaluator for scores_first + wins_period + team totals`
  4. Task 4 → `feat(mlb): add compute_prop_fair generic joint-probability pricer`
  5. Task 5 → `refactor(mlb): extend mlb_game_samples query with total + F3/F7`
  6. Task 6 → `refactor(mlb): rewire mlb_triple_play.R to use generic parser`
  7. Task 7 → `refactor(mlb): remove compute_triple_play_fair (superseded by parse_legs + compute_prop_fair)`
  8. Task 8 → `feat(mlb): add prop_type column to triple-play output table`
  9. Task 9 → `feat(wagerzon): add recon_specials.py for MLB specials page`
  10. Task 10 → `docs(mlb): document prop parser in Triple-Play Data Flow`
- Cleanup after merge: `git worktree remove .worktrees/prop-parser-pricer && git branch -d feature/prop-parser-pricer`
- Never merge to main without explicit user approval

---

## Documentation Plan

Task 10 extends `Answer Keys/CLAUDE.md`'s existing "Triple-Play Data Flow" subsection (added 2026-04-21) with the new parser layer. No other README updates required for this plan — the scraper README update belongs to the follow-up plan.

---

## Pre-Merge Review Checklist

Run `git diff main..HEAD` in the worktree and verify:

- **Behavior preservation**: for the 8 triple-plays we priced today, the NEW generic path must produce the same fair probabilities as the OLD `compute_triple_play_fair` — verified by Task 6's regression test which hard-codes expected values from today's run.
- **Grand slam coverage**: Task 4's tests exercise a 4-leg grand slam fixture end-to-end.
- **NA safety**: `parse_legs()` returns NULL on unknown tokens; `compute_prop_fair()` returns `NA_real_` on NULL legs — Task 2 + Task 4 cover both.
- **Unicode ½**: Task 1 has a unit test asserting `parse_fraction("U2½")` returns 2.5. R regex must use `perl = TRUE`.
- **Dead code**: `compute_triple_play_fair` is fully removed in Task 7; no callers remain after Task 6's rewire.
- **Recon script**: Task 9 creates the script but does NOT run it (that's user-manual). Verify the script has `if __name__ == "__main__"` guard and clean Playwright teardown.
- **No DB writes**: this plan makes zero writes to `wagerzon.duckdb`, `mlb.duckdb`, or any pipeline DB. All new DuckDB usage is read-only.

Fix any issues. Document findings as ISSUES-TO-FIX vs ACCEPTABLE-RISKS. Get explicit user approval before merging.

---

## Task 1: Pure helper — parse_fraction

**Files:**
- Create: `Answer Keys/parse_legs.R`
- Create: `Answer Keys/tests/test_parse_legs.R`

### - [ ] Step 1: Write failing tests

```r
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
```

### - [ ] Step 2: Run tests to verify they fail

```bash
cd "Answer Keys/tests" && Rscript -e 'library(testthat); test_file("test_parse_legs.R")'
```
Expected: failure with `could not find function "parse_fraction"` or `cannot open file '../parse_legs.R'`.

### - [ ] Step 3: Implement parse_fraction

Create `Answer Keys/parse_legs.R` with only this content (rest of the file is added in later tasks):

```r
# Answer Keys/parse_legs.R
# Wagerzon prop-description parser and generic leg-based joint probability pricer.

#' Convert a Wagerzon fraction string to a numeric.
#'
#' Handles three formats Wagerzon uses interchangeably:
#'   "3"   → 3
#'   "2.5" → 2.5
#'   "2½"  → 2.5 (using the U+00BD unicode fraction)
#' Also handles ¼ (U+00BC → 0.25) and ¾ (U+00BE → 0.75) for future-proofing.
#'
#' Returns NA_real_ on unparseable input.
parse_fraction <- function(s) {
  if (is.null(s) || is.na(s) || !nzchar(s)) return(NA_real_)
  frac_map <- c("¼" = 0.25, "½" = 0.5, "¾" = 0.75)
  # Try pure numeric first (handles "3", "2.5")
  n <- suppressWarnings(as.numeric(s))
  if (!is.na(n)) return(n)
  # Check for unicode fraction suffix
  last_char <- substr(s, nchar(s), nchar(s))
  if (last_char %in% names(frac_map)) {
    int_part <- suppressWarnings(as.numeric(substr(s, 1L, nchar(s) - 1L)))
    if (is.na(int_part)) int_part <- 0
    return(int_part + frac_map[[last_char]])
  }
  NA_real_
}
```

### - [ ] Step 4: Run tests to verify all pass

```bash
cd "Answer Keys/tests" && Rscript -e 'library(testthat); test_file("test_parse_legs.R")'
```
Expected: 6 tests pass.

### - [ ] Step 5: Commit

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/prop-parser-pricer
git add "Answer Keys/parse_legs.R" "Answer Keys/tests/test_parse_legs.R"
git commit -m "feat(mlb): add parse_fraction helper for Wagerzon ½ notation"
```

---

## Task 2: parse_legs + TOKEN_REGISTRY

**Files:**
- Modify: `Answer Keys/parse_legs.R` — append the registry + parser
- Modify: `Answer Keys/tests/test_parse_legs.R` — append parser tests

### - [ ] Step 1: Append failing tests

Append to `Answer Keys/tests/test_parse_legs.R`:

```r

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
```

### - [ ] Step 2: Run tests to verify they fail

```bash
cd "Answer Keys/tests" && Rscript -e 'library(testthat); test_file("test_parse_legs.R")'
```
Expected: 7 new failures with `could not find function "parse_legs"` (the 6 parse_fraction tests still pass).

### - [ ] Step 3: Append parse_legs + registry to parse_legs.R

Append to `Answer Keys/parse_legs.R`:

```r

#' Token registry: regex → function(match) returning a leg spec.
#'
#' Each entry has:
#'   pattern: POSIX-extended regex (anchored with ^ and $) matched against the
#'            stripped, uppercased token string.
#'   spec:    function(character vector of capture groups) returning a named list.
#'
#' Add new tokens here to teach the parser new prop types.
TOKEN_REGISTRY <- list(
  list(pattern = "^F3$",
       spec    = function(m) list(type = "wins_period", period = "F3")),
  list(pattern = "^1H$",
       spec    = function(m) list(type = "wins_period", period = "F5")),
  list(pattern = "^F7$",
       spec    = function(m) list(type = "wins_period", period = "F7")),
  list(pattern = "^GM$",
       spec    = function(m) list(type = "wins_period", period = "FG")),
  list(pattern = "^SCR 1ST$",
       spec    = function(m) list(type = "scores_first")),
  # For patterns with a capture group, m[[1]] is the full match and
  # m[[2]] is the first capture. parse_fraction needs the capture only.
  list(pattern = "^SCR U(.+)$",
       spec    = function(m) list(type = "team_total_under",
                                  line = parse_fraction(m[[2]]))),
  list(pattern = "^SCR O(.+)$",
       spec    = function(m) list(type = "team_total_over",
                                  line = parse_fraction(m[[2]])))
)

#' Match a single token against TOKEN_REGISTRY. Returns leg spec or NULL.
.match_token <- function(token) {
  for (entry in TOKEN_REGISTRY) {
    m <- regmatches(token, regexec(entry$pattern, token, perl = TRUE))[[1]]
    if (length(m) > 0) return(entry$spec(m))
  }
  NULL
}

#' Parse a Wagerzon description string into a list of leg specs.
#'
#' "GIANTS TRIPLE-PLAY (SCR 1ST, 1H & GM)" → 3 leg specs.
#' Unknown tokens trigger a warning and return NULL (pricer treats NULL as
#' "cannot price" → fair_prob = NA).
#' Descriptions with no parenthetical also return NULL.
parse_legs <- function(description) {
  if (is.null(description) || is.na(description) || !nzchar(description)) {
    return(NULL)
  }
  # Extract the parenthetical content (e.g. "SCR 1ST, 1H & GM")
  inside <- regmatches(description,
                       regexec("\\((.+)\\)", description, perl = TRUE))[[1]]
  if (length(inside) < 2) return(NULL)
  raw <- inside[[2]]
  # Tokenize: split on "," and "&", trim whitespace
  tokens <- trimws(unlist(strsplit(raw, "[,&]")))
  tokens <- tokens[nzchar(tokens)]
  # Map each token via registry
  legs <- vector("list", length(tokens))
  for (i in seq_along(tokens)) {
    spec <- .match_token(tokens[[i]])
    if (is.null(spec)) {
      warning(sprintf("Unknown leg token: '%s' in description '%s'",
                      tokens[[i]], description))
      return(NULL)
    }
    legs[[i]] <- spec
  }
  legs
}
```

### - [ ] Step 4: Run tests to verify all pass

```bash
cd "Answer Keys/tests" && Rscript -e 'library(testthat); test_file("test_parse_legs.R")'
```
Expected: 13 tests pass (6 from Task 1 + 7 new).

### - [ ] Step 5: Commit

```bash
git add "Answer Keys/parse_legs.R" "Answer Keys/tests/test_parse_legs.R"
git commit -m "feat(mlb): add parse_legs token registry for triple-play + grand slam"
```

---

## Task 3: eval_leg evaluator

**Files:**
- Modify: `Answer Keys/parse_legs.R` — append evaluator
- Modify: `Answer Keys/tests/test_parse_legs.R` — append evaluator tests

### - [ ] Step 1: Append failing tests

```r

test_that("eval_leg scores_first home matches home_scored_first == 1L", {
  samples <- data.frame(
    home_scored_first = c(1L, 0L, 1L, NA)
  )
  team_runs <- c(NA, NA, NA, NA)
  opp_runs  <- c(NA, NA, NA, NA)
  result <- eval_leg(list(type = "scores_first"),
                     samples, side = "home",
                     team_runs = team_runs, opp_runs = opp_runs)
  expect_equal(result, c(TRUE, FALSE, TRUE, NA))
})

test_that("eval_leg scores_first away matches home_scored_first == 0L", {
  samples <- data.frame(home_scored_first = c(1L, 0L, 1L))
  result <- eval_leg(list(type = "scores_first"),
                     samples, side = "away",
                     team_runs = c(NA, NA, NA), opp_runs = c(NA, NA, NA))
  expect_equal(result, c(FALSE, TRUE, FALSE))
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
```

### - [ ] Step 2: Run tests to verify they fail

```bash
cd "Answer Keys/tests" && Rscript -e 'library(testthat); test_file("test_parse_legs.R")'
```
Expected: 10 new failures with `could not find function "eval_leg"`.

### - [ ] Step 3: Append eval_leg to parse_legs.R

Append to `Answer Keys/parse_legs.R`:

```r

#' Evaluate a single leg spec against a samples frame. Returns a logical vector
#' of length nrow(samples).
#'
#' @param leg Named list returned by parse_legs() (e.g. list(type="scores_first")).
#' @param samples data.frame with columns matching what the leg type needs.
#'   scores_first → home_scored_first
#'   wins_period  → home_margin, home_margin_f3, home_margin_f5, or home_margin_f7
#'   team_total_* → nothing from samples directly (uses team_runs arg)
#'   opp_total_*  → uses opp_runs arg
#' @param side "home" or "away" — direction of inequalities and indicator value.
#' @param team_runs numeric vector of length nrow(samples), runs scored by the
#'   prop's subject team. Computed once per game by caller from home_margin +
#'   total_final_score.
#' @param opp_runs  numeric vector, runs scored by the opposing team.
eval_leg <- function(leg, samples, side, team_runs, opp_runs) {
  switch(leg$type,
    scores_first = {
      target <- if (side == "home") 1L else 0L
      samples$home_scored_first == target
    },
    wins_period = {
      col <- switch(leg$period,
                    F3 = samples$home_margin_f3,
                    F5 = samples$home_margin_f5,
                    F7 = samples$home_margin_f7,
                    FG = samples$home_margin,
                    stop("Unknown period: ", leg$period))
      if (side == "home") col > 0 else col < 0
    },
    team_total_under = team_runs <  leg$line,
    team_total_over  = team_runs >  leg$line,
    opp_total_under  = opp_runs  <  leg$line,
    opp_total_over   = opp_runs  >  leg$line,
    stop("Unknown leg type: ", leg$type)
  )
}
```

### - [ ] Step 4: Run tests to verify all pass

```bash
cd "Answer Keys/tests" && Rscript -e 'library(testthat); test_file("test_parse_legs.R")'
```
Expected: 23 tests pass (13 prior + 10 new).

### - [ ] Step 5: Commit

```bash
git add "Answer Keys/parse_legs.R" "Answer Keys/tests/test_parse_legs.R"
git commit -m "feat(mlb): add eval_leg evaluator for scores_first + wins_period + team totals"
```

---

## Task 4: compute_prop_fair generic pricer

**Files:**
- Modify: `Answer Keys/parse_legs.R` — append pricer
- Modify: `Answer Keys/tests/test_parse_legs.R` — append pricer tests

### - [ ] Step 1: Append failing tests

```r

test_that("compute_prop_fair triple-play home returns joint P(all 3 legs)", {
  # 10 synthetic rows. home_triple = scored_first==1 AND margin_f5>0 AND home_margin>0
  # Expected hits (same fixture as test_triple_play's home-side test): 5/10 = 0.5
  samples <- data.frame(
    home_margin       = c( 2,  3, -1,  4,  1,  5, -2,  2,  6, -3),
    home_margin_f5    = c( 1,  2,  0,  3,  0,  2, -1,  1,  4, -2),
    home_scored_first = c( 1L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 1L, 0L),
    total_final_score = rep(10, 10)  # irrelevant for triple-play
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
    home_margin       = c( 2,  3,  1,  2,  2),
    home_margin_f5    = c( 1,  2,  1,  1,  1),
    home_scored_first = c( 1L, 1L, 1L, 1L, 1L),
    total_final_score = c( 4,  6,  2,  4,  4)
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
    home_scored_first = 1L, total_final_score = 2
  )
  expect_true(is.na(compute_prop_fair(samples, "home", NULL)))
})

test_that("compute_prop_fair returns NA_real_ on empty legs list", {
  samples <- data.frame(
    home_margin = 1, home_margin_f5 = 1,
    home_scored_first = 1L, total_final_score = 2
  )
  expect_true(is.na(compute_prop_fair(samples, "home", list())))
})

test_that("compute_prop_fair returns NA_real_ on empty samples", {
  samples <- data.frame(
    home_margin = integer(0), home_margin_f5 = integer(0),
    home_scored_first = integer(0), total_final_score = integer(0)
  )
  legs <- list(list(type = "scores_first"))
  expect_true(is.na(compute_prop_fair(samples, "home", legs)))
})

test_that("compute_prop_fair excludes NA home_scored_first rows", {
  # 4 rows: 3 valid with 1 hit, 1 NA. Expected: 1/3, not 1/4.
  samples <- data.frame(
    home_margin       = c(1, 1, 1, 1),
    home_margin_f5    = c(1, 1, 1, 1),
    home_scored_first = c(1L, 0L, 0L, NA),
    total_final_score = c(2, 2, 2, 2)
  )
  legs <- list(
    list(type = "scores_first"),
    list(type = "wins_period", period = "F5"),
    list(type = "wins_period", period = "FG")
  )
  expect_equal(compute_prop_fair(samples, "home", legs), 1/3)
})

test_that("compute_prop_fair away side with team_total_under uses away_runs", {
  # away_runs = (total - margin) / 2
  # For side=away to win: home_margin_f5 < 0, home_margin < 0, home_scored_first == 0
  samples <- data.frame(
    home_margin       = c(-2, -3, -1),
    home_margin_f5    = c(-1, -2, -1),
    home_scored_first = c( 0L, 0L, 0L),
    total_final_score = c( 4,  6,  2)
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
```

### - [ ] Step 2: Run tests to verify they fail

```bash
cd "Answer Keys/tests" && Rscript -e 'library(testthat); test_file("test_parse_legs.R")'
```
Expected: 7 new failures with `could not find function "compute_prop_fair"`.

### - [ ] Step 3: Append compute_prop_fair to parse_legs.R

Append to `Answer Keys/parse_legs.R`:

```r

#' Generic joint-probability pricer for any multi-leg prop.
#'
#' Drops rows with NA home_scored_first, derives team/opp runs from margin +
#' total, then AND-reduces the leg evaluators across samples and returns the
#' mean of that logical vector. Returns NA_real_ if no valid rows or no legs.
#'
#' @param samples data.frame with at minimum: home_margin, home_margin_f5,
#'   home_scored_first, total_final_score. wins_period legs for F3/F7 also
#'   need home_margin_f3 / home_margin_f7 respectively.
#' @param side "home" or "away" — which team is the prop's subject.
#' @param legs list of leg specs as returned by parse_legs().
compute_prop_fair <- function(samples, side, legs) {
  if (is.null(legs) || length(legs) == 0) return(NA_real_)
  samples <- samples[!is.na(samples$home_scored_first), ]
  if (nrow(samples) == 0) return(NA_real_)

  # Derive team-specific run totals once from margin + total.
  home_runs <- (samples$total_final_score + samples$home_margin) / 2
  away_runs <- (samples$total_final_score - samples$home_margin) / 2
  team_runs <- if (side == "home") home_runs else away_runs
  opp_runs  <- if (side == "home") away_runs else home_runs

  hits <- rep(TRUE, nrow(samples))
  for (leg in legs) {
    hits <- hits & eval_leg(leg, samples, side, team_runs, opp_runs)
  }
  mean(hits)
}
```

### - [ ] Step 4: Run tests to verify all pass

```bash
cd "Answer Keys/tests" && Rscript -e 'library(testthat); test_file("test_parse_legs.R")'
```
Expected: 30 tests pass (23 prior + 7 new).

### - [ ] Step 5: Commit

```bash
git add "Answer Keys/parse_legs.R" "Answer Keys/tests/test_parse_legs.R"
git commit -m "feat(mlb): add compute_prop_fair generic joint-probability pricer"
```

---

## Task 5: Extend samples query in mlb_triple_play.R

**Files:**
- Modify: `Answer Keys/mlb_triple_play.R` — update the SQL string

### - [ ] Step 1: Read current state

The file currently has:
```r
samples_df <- dbGetQuery(con,
  "SELECT game_id, home_margin, home_margin_f5, home_scored_first
   FROM mlb_game_samples")
```

### - [ ] Step 2: Extend the SELECT to include total + F3/F7

Edit `Answer Keys/mlb_triple_play.R`. Find the `samples_df <- dbGetQuery(...)` call. Replace with:

```r
  samples_df <- dbGetQuery(con,
    "SELECT game_id, home_margin, total_final_score,
            home_margin_f3, home_margin_f5, home_margin_f7,
            home_scored_first
     FROM mlb_game_samples")
```

### - [ ] Step 3: Verify syntax

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/prop-parser-pricer
Rscript -e 'parse(file = "Answer Keys/mlb_triple_play.R"); cat("SYNTAX OK\n")'
```
Expected: `SYNTAX OK`.

### - [ ] Step 4: Run the existing triple-play test suite

```bash
cd "Answer Keys/tests" && Rscript -e 'library(testthat); test_file("test_triple_play.R")'
```
Expected: all passing tests from Task 1-4 of the prior plan still pass (none of those tests touch the SQL query).

### - [ ] Step 5: Commit

```bash
git add "Answer Keys/mlb_triple_play.R"
git commit -m "refactor(mlb): extend mlb_game_samples query with total + F3/F7"
```

---

## Task 6: Rewire pricer to use generic parser

**Files:**
- Modify: `Answer Keys/mlb_triple_play.R` — source parse_legs.R, extend tribble, replace per-leg logic
- Modify: `Answer Keys/tests/test_triple_play.R` — add regression test asserting grand slam can be priced

### - [ ] Step 1: Source parse_legs.R from the main block

Edit `Answer Keys/mlb_triple_play.R`. Inside the `if (!interactive() && sys.nframe() == 0L)` guard, find this existing line:

```r
  setwd("~/NFLWork/Answer Keys")
```

Immediately AFTER that line, add one new line:

```r
  source("parse_legs.R")
```

That's the only source-loading edit in this task. Do NOT add a `source()` call at the top of the file — `setwd()` establishes the working directory that the relative `source("parse_legs.R")` path depends on, and that only happens inside the main block.

### - [ ] Step 2: Extend todays_lines tribble with description column

Edit the `todays_lines <- tribble(...)` call. Replace entirely with:

```r
  # Today's book lines — edit this tribble whenever new props post.
  # home_team / away_team must match Odds API canonical names in mlb_consensus_temp.
  # description is the Wagerzon label verbatim; parse_legs() derives leg logic from it.
  todays_lines <- tribble(
    ~home_team,              ~away_team,             ~target_team, ~side,   ~book_odds, ~description,
    "Colorado Rockies",      "San Diego Padres",     "Rockies",    "home",  +530,       "ROCKIES TRIPLE-PLAY (SCR 1ST, 1H & GM)",
    "Colorado Rockies",      "San Diego Padres",     "Padres",     "away",  +190,       "PADRES TRIPLE-PLAY (SCR 1ST, 1H & GM)",
    "San Francisco Giants",  "Los Angeles Dodgers",  "Giants",     "home",  +750,       "GIANTS TRIPLE-PLAY (SCR 1ST, 1H & GM)",
    "San Francisco Giants",  "Los Angeles Dodgers",  "Dodgers",    "away",  +155,       "DODGERS TRIPLE-PLAY (SCR 1ST, 1H & GM)",
    "Seattle Mariners",      "Athletics",            "Mariners",   "home",  +215,       "MARINERS TRIPLE-PLAY (SCR 1ST, 1H & GM)",
    "Seattle Mariners",      "Athletics",            "Athletics",  "away",  +455,       "ATHLETICS TRIPLE-PLAY (SCR 1ST, 1H & GM)",
    "Arizona Diamondbacks",  "Chicago White Sox",    "DBacks",     "home",  +240,       "DBACKS TRIPLE-PLAY (SCR 1ST, 1H & GM)",
    "Arizona Diamondbacks",  "Chicago White Sox",    "White Sox",  "away",  +415,       "WHITE SOX TRIPLE-PLAY (SCR 1ST, 1H & GM)"
  )
```

### - [ ] Step 3: Replace the rowwise mutate block to use parse_legs + compute_prop_fair

Find the `priced <- matched %>% rowwise() %>% mutate(...)` block. Replace the mutate call entirely with:

```r
  priced <- matched %>%
    rowwise() %>%
    mutate(
      game_samples = list(samples_df[samples_df$game_id == id, ]),
      n_samples    = nrow(game_samples),
      legs         = list(parse_legs(description)),
      fair_prob    = compute_prop_fair(game_samples, side, legs),
      fair_odds    = prob_to_american(fair_prob),
      book_prob    = american_to_prob(book_odds),
      edge_pct     = (fair_prob / book_prob - 1) * 100
    ) %>%
    ungroup() %>%
    select(target_team, side, n_samples, fair_prob, fair_odds, book_odds, edge_pct) %>%
    arrange(desc(edge_pct))
```

### - [ ] Step 4: Add a regression test for the triple-play behavior

In `Answer Keys/tests/test_triple_play.R`, append a new test at the end:

```r

test_that("compute_prop_fair on triple-play legs matches prior compute_triple_play_fair output", {
  # Source both the old function (still present until Task 7) and the new generic pricer
  source("../parse_legs.R")

  # Same 10-row fixture used in the original compute_triple_play_fair test
  samples <- data.frame(
    home_margin       = c( 2,  3, -1,  4,  1,  5, -2,  2,  6, -3),
    home_margin_f5    = c( 1,  2,  0,  3,  0,  2, -1,  1,  4, -2),
    home_scored_first = c( 1L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 1L, 0L),
    total_final_score = rep(10, 10)
  )
  triple_legs <- list(
    list(type = "scores_first"),
    list(type = "wins_period", period = "F5"),
    list(type = "wins_period", period = "FG")
  )
  new_result <- compute_prop_fair(samples, "home", triple_legs)
  old_result <- compute_triple_play_fair(samples, "home")
  expect_equal(new_result, old_result)
})
```

### - [ ] Step 5: Run the test file

```bash
cd "Answer Keys/tests" && Rscript -e 'library(testthat); test_file("test_triple_play.R")'
```
Expected: all tests pass including the new regression test.

### - [ ] Step 6: Copy main's mlb.duckdb into the worktree for regression verification

Per project rule, DuckDB files must be COPIED, not symlinked:

```bash
cp "/Users/callancapitolo/NFLWork/Answer Keys/mlb.duckdb" \
   "/Users/callancapitolo/NFLWork/.worktrees/prop-parser-pricer/Answer Keys/mlb.duckdb"
```

Note: use double-quoted paths, never backslash-escaped spaces (project rule). The copied file is gitignored by the existing `*.duckdb` rule, so no risk of accidental commit. Remove after verification.

### - [ ] Step 7: Run the rewired pricer against the copied mlb.duckdb

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/prop-parser-pricer/Answer Keys"
Rscript mlb_triple_play.R 2>&1 | tail -20
```

Expected: table prints with 8 rows (matched 8/8 posted lines), `fair_odds` populated for all 8 triple-plays. Compare the `fair_odds` column against today's committed run:

| Team | Today's fair | Acceptable range |
|---|---:|---:|
| Giants (home) | +369 | +360 to +380 |
| Athletics (away) | +347 | +340 to +355 |
| White Sox (away) | +320 | +315 to +325 |
| Rockies (home) | +447 | +440 to +455 |
| Dodgers (away) | +137 | +133 to +141 |
| Mariners (home) | +196 | +192 to +200 |
| DBacks (home) | +227 | +222 to +232 |
| Padres (away) | +186 | +182 to +190 |

Identical values are ideal (same samples, same math). Deviation beyond the acceptable range indicates a regression — stop and investigate before proceeding to Step 8.

### - [ ] Step 8: Remove the copied DuckDB from the worktree

```bash
rm "/Users/callancapitolo/NFLWork/.worktrees/prop-parser-pricer/Answer Keys/mlb.duckdb"
```

This keeps the worktree clean and prevents stale data from lingering. `git status` should still show a clean working tree after this delete because the file was gitignored.

### - [ ] Step 9: Commit

```bash
git add "Answer Keys/mlb_triple_play.R" "Answer Keys/tests/test_triple_play.R"
git commit -m "refactor(mlb): rewire mlb_triple_play.R to use generic parser"
```

---

## Task 7: Remove compute_triple_play_fair (dead after rewire)

**Files:**
- Modify: `Answer Keys/mlb_triple_play.R` — delete the old function
- Modify: `Answer Keys/tests/test_triple_play.R` — delete tests for the removed function + the regression test from Task 6

### - [ ] Step 1: Remove compute_triple_play_fair from mlb_triple_play.R

Delete the entire function definition from `Answer Keys/mlb_triple_play.R`:

```r
compute_triple_play_fair <- function(samples, side = c("home", "away")) {
  side <- match.arg(side)
  samples <- samples[!is.na(samples$home_scored_first), ]
  if (nrow(samples) == 0) return(NA_real_)
  if (side == "home") {
    mean(samples$home_scored_first == 1L &
         samples$home_margin_f5   > 0 &
         samples$home_margin      > 0)
  } else {
    mean(samples$home_scored_first == 0L &
         samples$home_margin_f5   < 0 &
         samples$home_margin      < 0)
  }
}
```

Also delete its roxygen docstring (the `#' Fair probability for a team's triple-play...` block above it).

### - [ ] Step 2: Remove compute_triple_play_fair tests from test_triple_play.R

Delete ALL of these test_that blocks from `Answer Keys/tests/test_triple_play.R`:

1. `"compute_triple_play_fair (home) returns joint P of 3 legs"`
2. `"compute_triple_play_fair (away) returns joint P for away side"`
3. `"compute_triple_play_fair excludes NA home_scored_first rows"`
4. `"compute_triple_play_fair returns NA_real_ on empty samples"`
5. `"compute_triple_play_fair returns 1.0 when all rows pass all 3 legs"`
6. `"compute_triple_play_fair returns 0.0 when no rows pass"`
7. `"compute_prop_fair on triple-play legs matches prior compute_triple_play_fair output"` (the Task 6 regression test — its purpose was to bridge the rewire; now it references a deleted function)

Keep these tests unchanged:
- All `determine_home_scored_first*` tests (6 total from Task 1 of prior plan)
- `"prob_to_american handles favorites and dogs"`
- `"american_to_prob is the inverse of prob_to_american"`
- `"american_to_prob(0) returns NA_real_"` (this tests the surviving function)
- `"mlb_game_samples has home_scored_first column populated"` smoke test

Coverage of the boundary cases (empty samples, all-pass, no-pass) is preserved because `test_parse_legs.R` has equivalent tests against `compute_prop_fair` (added in Task 4 of this plan).

### - [ ] Step 3: Run tests

```bash
cd "Answer Keys/tests" && Rscript -e 'library(testthat); test_file("test_triple_play.R")'
```
Expected: remaining tests pass. No references to `compute_triple_play_fair` anywhere.

### - [ ] Step 4: Syntax-check the pricer

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/prop-parser-pricer
Rscript -e 'parse(file = "Answer Keys/mlb_triple_play.R"); cat("SYNTAX OK\n")'
```
Expected: `SYNTAX OK`.

### - [ ] Step 5: Grep for lingering references

```bash
grep -rn "compute_triple_play_fair" "Answer Keys/" || echo "NO REFERENCES FOUND"
```
Expected: `NO REFERENCES FOUND`. If any match, delete or replace them.

### - [ ] Step 6: Commit

```bash
git add "Answer Keys/mlb_triple_play.R" "Answer Keys/tests/test_triple_play.R"
git commit -m "refactor(mlb): remove compute_triple_play_fair (superseded by parse_legs + compute_prop_fair)"
```

---

## Task 8: Add prop_type column to output display

**Files:**
- Modify: `Answer Keys/mlb_triple_play.R` — parse prop_type from description, add to select + display

### - [ ] Step 1: Add prop_type derivation

In `Answer Keys/mlb_triple_play.R`, find the `priced <- matched %>% rowwise() %>% mutate(...)` block from Task 6. Add a `prop_type` field inside the `mutate`:

```r
  priced <- matched %>%
    rowwise() %>%
    mutate(
      game_samples = list(samples_df[samples_df$game_id == id, ]),
      n_samples    = nrow(game_samples),
      legs         = list(parse_legs(description)),
      fair_prob    = compute_prop_fair(game_samples, side, legs),
      fair_odds    = prob_to_american(fair_prob),
      book_prob    = american_to_prob(book_odds),
      edge_pct     = (fair_prob / book_prob - 1) * 100,
      prop_type    = {
        # Anchor on known prop-type tokens so multi-word team names work
        # ("WHITE SOX TRIPLE-PLAY ..." and "GIANTS GRAND-SLAM ..." both parse).
        # Extend this pattern as new prop types are added to TOKEN_REGISTRY.
        known_props <- "(TRIPLE-PLAY|GRAND-SLAM)"
        m <- regmatches(description,
                        regexec(paste0("\\b", known_props, "\\b"), description))[[1]]
        if (length(m) >= 2) m[[2]] else NA_character_
      }
    ) %>%
    ungroup() %>%
    select(target_team, prop_type, side, n_samples,
           fair_prob, fair_odds, book_odds, edge_pct) %>%
    arrange(desc(edge_pct))
```

### - [ ] Step 2: Update the display formatting block

Find the `display <- priced %>% mutate(...)` block. Update to include prop_type (no format change — it's already a string):

```r
  display <- priced %>%
    mutate(
      fair_prob = sprintf("%.3f", fair_prob),
      fair_odds = ifelse(fair_odds > 0,
                         paste0("+", fair_odds), as.character(fair_odds)),
      book_odds = ifelse(book_odds > 0,
                         paste0("+", book_odds), as.character(book_odds)),
      edge_pct  = sprintf("%+.1f%%", edge_pct)
    )
  print(as.data.frame(display), row.names = FALSE)
```

### - [ ] Step 3: Syntax check

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/prop-parser-pricer
Rscript -e 'parse(file = "Answer Keys/mlb_triple_play.R"); cat("SYNTAX OK\n")'
```
Expected: `SYNTAX OK`.

### - [ ] Step 4: Run tests

```bash
cd "Answer Keys/tests" && Rscript -e 'library(testthat); test_file("test_parse_legs.R") ; test_file("test_triple_play.R")'
```
Expected: all tests pass.

### - [ ] Step 5: Commit

```bash
git add "Answer Keys/mlb_triple_play.R"
git commit -m "feat(mlb): add prop_type column to triple-play output table"
```

---

## Task 9: recon_specials.py — Wagerzon reconnaissance script

**Files:**
- Create: `wagerzon_odds/recon_specials.py`

### - [ ] Step 1: Create the recon script

Create `wagerzon_odds/recon_specials.py`:

```python
"""
Wagerzon MLB Specials Recon Script

Sibling of recon_wagerzon.py. Launches Chromium with the existing
.wagerzon_profile, navigates to the MLB specials page (league 4899),
captures every network request + response, and dumps to recon_specials.json.

Usage:
    python3 wagerzon_odds/recon_specials.py

Prerequisites:
    - Playwright installed: pip install playwright && playwright install chromium
    - .wagerzon_profile directory populated (log in via recon_wagerzon.py once
      if starting from scratch).

Workflow:
    1. Script launches persistent-context Chromium.
    2. Script navigates to the candidate specials URL (see CANDIDATE_URLS).
    3. You manually verify the specials page rendered correctly — if the URL
       is wrong, close the browser, fix CANDIDATE_URLS, re-run.
    4. Press ENTER in the terminal when the page is fully loaded with all
       triple-play / grand-slam rows visible.
    5. Captured traffic is saved to wagerzon_odds/recon_specials.json.
    6. Inspect the JSON to determine: JSON API endpoint? HTML-only? URL
       structure? Section-header format?

The follow-up plan uses this recon output to design the production scraper.
"""

from playwright.sync_api import sync_playwright
import json
import os

WAGERZON_URL = "https://backend.wagerzon.com"

# Candidate URLs for MLB specials. If the first fails, adjust based on what
# you see in the browser's devtools Network tab when you manually navigate
# to the specials page.
CANDIDATE_URLS = [
    f"{WAGERZON_URL}/wager/NewSchedule.aspx?WT=0&lg=4899",
]

PROFILE_DIR = os.path.join(os.path.dirname(__file__), ".wagerzon_profile")
OUT_FILE = os.path.join(os.path.dirname(__file__), "recon_specials.json")


def run_recon():
    captured = []

    def handle_request(request):
        if "wagerzon.com" not in request.url:
            return
        captured.append({
            "phase": "request",
            "url": request.url,
            "method": request.method,
            "headers": dict(request.headers),
            "post_data": request.post_data,
            "resource_type": request.resource_type,
        })
        rtype = request.resource_type
        tag = {"document": "[PAGE]", "xhr": "[API ]", "fetch": "[API ]"}.get(rtype, "[____]")
        print(f"  {tag} {request.method} {request.url[:130]}")

    def handle_response(response):
        if "wagerzon.com" not in response.url:
            return
        try:
            content_type = response.headers.get("content-type", "")
            body = None
            if "json" in content_type:
                body = response.text()[:20000]
            elif "html" in content_type and response.request.resource_type == "document":
                body = response.text()[:40000]
            captured.append({
                "phase": "response",
                "url": response.url,
                "status": response.status,
                "content_type": content_type,
                "body_preview": body,
                "body_length": len(response.text()) if body is not None else None,
            })
        except Exception as e:
            captured.append({
                "phase": "response",
                "url": response.url,
                "status": response.status,
                "error": str(e),
            })

    with sync_playwright() as p:
        context = p.chromium.launch_persistent_context(
            PROFILE_DIR,
            headless=False,
            viewport={"width": 1440, "height": 900},
        )
        page = context.new_page()
        page.on("request", handle_request)
        page.on("response", handle_response)

        for url in CANDIDATE_URLS:
            print(f"\n=== Navigating to candidate URL: {url} ===")
            try:
                page.goto(url, wait_until="networkidle", timeout=60000)
                print("Page loaded. Verify in browser that specials are visible.")
            except Exception as e:
                print(f"Navigation failed: {e}")

        input("\nPress ENTER when page is fully loaded with all MLB specials visible...")

        out = {
            "candidate_urls": CANDIDATE_URLS,
            "captured_events": captured,
        }
        with open(OUT_FILE, "w") as f:
            json.dump(out, f, indent=2, default=str)
        print(f"\nSaved {len(captured)} captured events to {OUT_FILE}")

        input("Press ENTER to close the browser...")
        context.close()


if __name__ == "__main__":
    run_recon()
```

### - [ ] Step 2: Syntax check

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/prop-parser-pricer
python3 -c "import ast; ast.parse(open('wagerzon_odds/recon_specials.py').read()); print('SYNTAX OK')"
```
Expected: `SYNTAX OK`.

### - [ ] Step 3: Import check (does not execute the browser)

```bash
python3 -c "
import sys
sys.path.insert(0, 'wagerzon_odds')
import recon_specials
assert hasattr(recon_specials, 'run_recon'), 'run_recon function missing'
assert recon_specials.CANDIDATE_URLS, 'CANDIDATE_URLS empty'
print('IMPORT OK')
"
```
Expected: `IMPORT OK`.

### - [ ] Step 4: Commit

```bash
git add "wagerzon_odds/recon_specials.py"
git commit -m "feat(wagerzon): add recon_specials.py for MLB specials page"
```

---

## Task 10: Documentation — parser layer in CLAUDE.md

**Files:**
- Modify: `Answer Keys/CLAUDE.md` — extend the existing "Triple-Play Data Flow" section

### - [ ] Step 1: Update the data-flow diagram

In `Answer Keys/CLAUDE.md`, find the "Triple-Play Data Flow" section (added 2026-04-21). Replace the `mlb_triple_play.R (standalone pricer)` block of the diagram and the bullets beneath with:

Find this existing content:
```
mlb_triple_play.R (standalone pricer)
  ├── Reads mlb_game_samples + hardcoded today's lines
  ├── compute_triple_play_fair(samples, side) → joint P across the 3 legs
  └── Prints fair odds vs book + edge per posted line
```
- `home_scored_first` is NA for games with no scoring in innings 1–5 (~5% of historical games); excluded from the mean before the ratio is computed.
- F5 3-way: leg passes only with strict lead (`margin_f5 > 0` for home, `< 0` for away). F5 ties kill the parlay.
- Helper lives in `triple_play_helpers.R` so both MLB.R and the pricer can source it without duplication.

Replace with:
```
parse_legs.R (generic prop parser)
  ├── TOKEN_REGISTRY maps Wagerzon description tokens to structured leg specs:
  │     SCR 1ST → scores_first;  1H → wins_period(F5);  GM → wins_period(FG);
  │     F3/F7 → wins_period(F3/F7);  SCR U<N>/O<N> → team_total_under/over
  ├── parse_legs(description) → list of leg specs (or NULL + warning on unknown)
  ├── eval_leg(leg, samples, side, team_runs, opp_runs) → logical vector
  └── compute_prop_fair(samples, side, legs) → AND-reduce across legs,
        NA rows on home_scored_first excluded before the mean

mlb_triple_play.R (standalone pricer)
  ├── Reads mlb_game_samples (total_final_score + margin at F3/F5/F7/FG + scored_first)
  ├── For each tribble row: parse_legs(description) → compute_prop_fair(...)
  └── Prints fair odds vs book + edge per posted line, grouped by prop_type
```
- `home_scored_first` is NA for games with no scoring in innings 1–5 (~5% of historical games); excluded from the mean before the ratio is computed.
- F5 3-way: `wins_period` with period=F5 passes only with strict lead (`margin_f5 > 0` for home, `< 0` for away). F5 ties kill the parlay.
- `SCR U<N>` means the listed team's team-total under N (team_total_under leg). Numeric parser handles `"2"`, `"2.5"`, and unicode `"2½"`.
- Helpers: `triple_play_helpers.R` (determine_home_scored_first*) is sourced by MLB.R Phase 1; `parse_legs.R` is sourced by mlb_triple_play.R's main block. Both files are pure — no DB or network side effects.
- Adding a new prop type = add a token to `TOKEN_REGISTRY`. No new function needed.

### - [ ] Step 2: Verify the file renders sensibly

```bash
grep -A 40 "Triple-Play Data Flow" "Answer Keys/CLAUDE.md" | head -60
```
Expected: the updated section with the parser layer visible.

### - [ ] Step 3: Commit

```bash
git add "Answer Keys/CLAUDE.md"
git commit -m "docs(mlb): document prop parser in Triple-Play Data Flow"
```

---

## Self-Review

**Spec coverage (from `docs/superpowers/specs/2026-04-21-wagerzon-specials-scraper-design.md`):**

| Spec section | Covered by |
|---|---|
| Description parser with token registry | Task 2 (`parse_legs` + `TOKEN_REGISTRY`) |
| Generic `compute_prop_fair` + `eval_leg` | Tasks 3, 4 |
| Unicode fraction handling (`2½`) | Task 1 |
| Pricer wired to parser, sample query extended | Tasks 5, 6 |
| `prop_type` in output display | Task 8 |
| Recon script for specials page | Task 9 |
| Documentation | Task 10 |
| Scraper + DB schema + `run.py` integration | **Deferred to follow-up plan** (post-recon) |

The scraper + DB integration parts of the spec are deferred — noted in the plan's opening paragraph. This plan ships a parser that operates on a tribble today and is ready to operate on scraped DuckDB rows tomorrow with no parser-side changes.

**Placeholder scan:** No TBDs, no "implement later" without concrete code, no vague error-handling hand-waves. Each step has either exact code or an exact command.

**Type consistency:**
- `parse_legs()` returns `list` (of leg specs) or `NULL`. Consumer `compute_prop_fair` handles both.
- `eval_leg()` returns a logical vector of length `nrow(samples)`. Consumer ANDs it elementwise — correct shape.
- `compute_prop_fair()` returns scalar numeric or `NA_real_`. Consumer `prob_to_american()` handles NA.
- Column names `home_margin`, `home_margin_f3/f5/f7`, `home_scored_first`, `total_final_score` consistent across query (Task 5), samples exposure (Task 6), evaluator (Task 3), and pricer (Task 4).

No issues found.

---

## Execution Handoff

**Plan complete and saved to `docs/superpowers/plans/2026-04-21-prop-parser-and-pricer.md`. Two execution options:**

**1. Subagent-Driven (recommended)** — I dispatch a fresh subagent per task, review between tasks, fast iteration.

**2. Inline Execution** — execute tasks in this session with checkpoints.

**Which approach?**
