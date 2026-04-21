# MLB Triple-Play Pricer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Price MLB "triple-play" props — team scores first in the game AND wins the 1st half (F5, strict lead) AND wins the game — using the dispersion-matched samples already produced by the MLB answer-key pipeline.

**Architecture:** Compute a `home_scored_first` indicator per historical game during MLB.R Phase 1, plumb it through `generate_all_samples()` (which preserves all DT columns through resampling) and the `mlb_game_samples` export. Then build a standalone pricer `mlb_triple_play.R` that mirrors `mlb_correlated_parlay.R`'s pattern: read `mlb_game_samples`, compute the 3-leg joint probability empirically on the resampled rows, compare to a hardcoded table of today's 8 book lines.

**Tech Stack:** R (data.table, dplyr, duckdb, DBI), testthat for unit tests. Reuses the existing Tools.R sampling engine — no changes to `run_answer_key_sample()` or `generate_all_samples()`.

---

## File Structure

**Created:**
- `Answer Keys/triple_play_helpers.R` — pure functions: `determine_home_scored_first()`, `determine_home_scored_first_vec()`. Isolated so Phase 1 of MLB.R and the pricer can both `source()` it without duplicating logic.
- `Answer Keys/mlb_triple_play.R` — standalone pricer, mirrors `mlb_correlated_parlay.R` structure. Loads samples + hardcoded lines, prints fair-vs-book comparison.
- `Answer Keys/tests/test_triple_play.R` — unit tests for the helper, the pricer core function, and a smoke test against `mlb_game_samples`.

**Modified:**
- `Answer Keys/MLB.R` — Phase 1: compute `home_scored_first` on DT. Phase 4 samples export: add the column to the `sample_rows` tibble.
- `Answer Keys/CLAUDE.md` — new "Triple-Play Data Flow" subsection documenting how the indicator flows through the pipeline.

---

## Worktree & Version Control Plan

- **Before any code change**: create worktree and switch to feature branch.
  ```bash
  git worktree add .worktrees/mlb-triple-play -b feature/mlb-triple-play-pricer main
  cd .worktrees/mlb-triple-play
  git branch  # confirm feature/mlb-triple-play-pricer
  ```
- **DuckDB**: do NOT symlink `mlb.duckdb` into the worktree (per project CLAUDE.md — WAL corruption risk). Either copy it once for the Task 2 pipeline run, or run MLB.R from the worktree directory so it writes its own fresh `mlb.duckdb`. We'll copy, since a full MLB.R run takes ~10 minutes.
  ```bash
  cp "../../Answer Keys/mlb.duckdb" "Answer Keys/mlb.duckdb"
  cp "../../Answer Keys/pbp.duckdb" "Answer Keys/pbp.duckdb"
  ```
- **Commits** (one per task as it's completed):
  1. Task 1 → `feat(mlb): add first-scorer helper for triple-play pricing`
  2. Task 2 → `feat(mlb): plumb home_scored_first through samples pipeline`
  3. Task 3 → `feat(mlb): add triple-play fair probability computation`
  4. Task 4 → `feat(mlb): print triple-play fair odds vs today's book lines`
  5. Task 5 → `docs(mlb): document triple-play data flow`
- **Pre-merge review**: run through the checklist below, present findings to user, wait for explicit approval before merging.
- **Merge + cleanup** (only after approval):
  ```bash
  cd ../..
  git checkout main
  git merge --no-ff feature/mlb-triple-play-pricer
  git worktree remove .worktrees/mlb-triple-play
  git branch -d feature/mlb-triple-play-pricer
  ```

---

## Documentation Plan

Task 5 adds a "Triple-Play Data Flow" subsection to `Answer Keys/CLAUDE.md` modeled on the existing "Race-to-X Data Flow" subsection. It documents: where `home_scored_first` is computed, what NA means (~5% of games — no scoring in innings 1–5), how `mlb_game_samples` carries it, and how `mlb_triple_play.R` consumes it.

No top-level `README.md` or `CLAUDE.md` update needed — this is scoped entirely to Answer Keys.

---

## Pre-Merge Review Checklist

Run `git diff main..HEAD` in the worktree and verify:

- **Data integrity**: `home_scored_first` values are only 0, 1, or NA in `mlb_game_samples`. Spot-check 5 random sample rows by joining back to the original pbp row and recomputing by hand.
- **NA handling**: `compute_triple_play_fair()` excludes `home_scored_first == NA` rows before computing the mean. Without the exclusion, `NA == 1` would poison the mean.
- **F5 3-way ties**: fair prob uses strict inequality (`margin_f5 > 0` / `< 0`), not `>=` — ties kill the parlay.
- **Resource safety**: `mlb_triple_play.R` uses `on.exit(dbDisconnect(con))`. No new lock files.
- **Log hygiene**: output is console only — no log files created.
- **Dead code**: no unused functions or imports introduced.
- **Team-name matching**: `inner_join(consensus, by = c("home_team", "away_team"))` — if any of the 8 lines don't match (typo in `todays_lines`), they silently drop. Print `nrow(priced)` vs `nrow(todays_lines)` and warn if they differ.
- **Pipeline impact**: the MLB.R change adds one column to `DT` and one column to `mlb_game_samples`. Verify `mlb_correlated_parlay.R` and the dashboard still load samples without error (extra column is harmless, but confirm).

Document findings as ISSUES-TO-FIX vs ACCEPTABLE-RISKS. Fix issues. Then ask user explicitly: "Ready to merge to main?"

---

## Task 1: Pure helper — determine_home_scored_first

**Files:**
- Create: `Answer Keys/triple_play_helpers.R`
- Create: `Answer Keys/tests/test_triple_play.R`

### - [ ] Step 1: Write failing tests

```r
# Answer Keys/tests/test_triple_play.R
library(testthat)
source("../triple_play_helpers.R")

test_that("home scored first when away has 0 runs in inning 1", {
  # Inning 1: home +2, total 2 → home scored 2, away 0
  expect_equal(
    determine_home_scored_first(m1=2, t1=2,
                                m2=NA, t2=NA, m3=NA, t3=NA,
                                m4=NA, t4=NA, m5=NA, t5=NA),
    1L
  )
})

test_that("away scored first when away scored any runs in first scoring inning", {
  # Inning 1: total 3, margin -1 → away 2, home 1. Away bats first → away scored first.
  expect_equal(
    determine_home_scored_first(m1=-1, t1=3,
                                m2=NA, t2=NA, m3=NA, t3=NA,
                                m4=NA, t4=NA, m5=NA, t5=NA),
    0L
  )
})

test_that("falls through to inning 2 when inning 1 is scoreless", {
  # Inning 1: 0-0. Inning 2: total 1, margin +1 → home scored.
  expect_equal(
    determine_home_scored_first(m1=0, t1=0,
                                m2=1, t2=1, m3=NA, t3=NA,
                                m4=NA, t4=NA, m5=NA, t5=NA),
    1L
  )
})

test_that("returns NA when no scoring in innings 1 through 5", {
  expect_true(is.na(
    determine_home_scored_first(m1=0, t1=0, m2=0, t2=0,
                                m3=0, t3=0, m4=0, t4=0,
                                m5=0, t5=0)
  ))
})

test_that("returns NA when data is NA before any scoring", {
  expect_true(is.na(
    determine_home_scored_first(m1=NA, t1=NA,
                                m2=NA, t2=NA, m3=NA, t3=NA,
                                m4=NA, t4=NA, m5=NA, t5=NA)
  ))
})

test_that("vectorized version matches scalar over multiple rows", {
  result <- determine_home_scored_first_vec(
    m1 = c(2, -1, 0), t1 = c(2, 3, 0),
    m2 = c(NA, NA, 1), t2 = c(NA, NA, 1),
    m3 = c(NA, NA, NA), t3 = c(NA, NA, NA),
    m4 = c(NA, NA, NA), t4 = c(NA, NA, NA),
    m5 = c(NA, NA, NA), t5 = c(NA, NA, NA)
  )
  expect_equal(result, c(1L, 0L, 1L))
})
```

### - [ ] Step 2: Run tests to confirm they fail

Run from the worktree root:
```bash
cd "Answer Keys/tests" && Rscript -e 'library(testthat); test_file("test_triple_play.R")'
```
Expected: failure with `could not find function "determine_home_scored_first"`

### - [ ] Step 3: Implement the helper

```r
# Answer Keys/triple_play_helpers.R

#' Determine if the home team scored first in the game.
#'
#' Baseball half-inning order: away bats top (first), home bats bottom (second).
#' So if the away team scored any runs in the first inning where scoring
#' occurred, away scored first — even if home also scored in the bottom half.
#'
#' Inputs are cumulative inning margins (home-minus-away) and totals after
#' innings 1 through 5.
#'
#' @return 1L if home scored first, 0L if away scored first, NA_integer_ if
#'   no scoring through inning 5 (rare — ~5% of games).
determine_home_scored_first <- function(m1, t1, m2, t2, m3, t3, m4, t4, m5, t5) {
  m <- c(m1, m2, m3, m4, m5)
  t <- c(t1, t2, t3, t4, t5)
  prev_m <- 0
  prev_t <- 0
  for (i in seq_len(5)) {
    if (is.na(t[i]) || is.na(m[i])) return(NA_integer_)
    if (t[i] > prev_t) {
      dt <- t[i] - prev_t
      dm <- m[i] - prev_m
      away_runs <- (dt - dm) / 2
      home_runs <- (dt + dm) / 2
      if (away_runs > 0) return(0L)
      if (home_runs > 0) return(1L)
      return(NA_integer_)  # shouldn't hit — dt > 0 means someone scored
    }
    prev_m <- m[i]
    prev_t <- t[i]
  }
  NA_integer_
}

#' Vectorized wrapper for use inside data.table/dplyr mutates.
determine_home_scored_first_vec <- function(m1, t1, m2, t2, m3, t3, m4, t4, m5, t5) {
  mapply(determine_home_scored_first,
         m1, t1, m2, t2, m3, t3, m4, t4, m5, t5,
         USE.NAMES = FALSE)
}
```

### - [ ] Step 4: Run tests to confirm all pass

Run:
```bash
cd "Answer Keys/tests" && Rscript -e 'library(testthat); test_file("test_triple_play.R")'
```
Expected: all 6 tests pass.

### - [ ] Step 5: Commit

```bash
cd ../..  # back to worktree root
git add "Answer Keys/triple_play_helpers.R" "Answer Keys/tests/test_triple_play.R"
git commit -m "feat(mlb): add first-scorer helper for triple-play pricing"
```

---

## Task 2: Plumb home_scored_first through MLB.R

**Files:**
- Modify: `Answer Keys/MLB.R` — two edits: Phase 1 (after DT load) and Phase 4 (samples export)
- Modify: `Answer Keys/tests/test_triple_play.R` — add smoke test

### - [ ] Step 1: Append a smoke test to test_triple_play.R

```r
# Append to Answer Keys/tests/test_triple_play.R

test_that("mlb_game_samples has home_scored_first column populated", {
  skip_if_not(file.exists("../mlb.duckdb"),
              "mlb.duckdb not present — run MLB.R first")
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = "../mlb.duckdb", read_only = TRUE)
  on.exit(DBI::dbDisconnect(con))
  cols <- DBI::dbGetQuery(con, "PRAGMA table_info('mlb_game_samples')")$name
  expect_true("home_scored_first" %in% cols)
  # Values should be 0, 1, or NA only
  vals <- DBI::dbGetQuery(con,
    "SELECT DISTINCT home_scored_first FROM mlb_game_samples")$home_scored_first
  expect_true(all(vals %in% c(0L, 1L, NA_integer_)))
  # Coverage should be ~95% non-NA per historical data
  non_na_pct <- DBI::dbGetQuery(con,
    "SELECT AVG(CASE WHEN home_scored_first IS NOT NULL THEN 1.0 ELSE 0.0 END) AS p
     FROM mlb_game_samples")$p
  expect_gt(non_na_pct, 0.90)
})
```

### - [ ] Step 2: Run tests — smoke test should fail

Run:
```bash
cd "Answer Keys/tests" && Rscript -e 'library(testthat); test_file("test_triple_play.R")'
```
Expected: smoke test FAILS with `home_scored_first not in cols` (or skipped if no mlb.duckdb yet).

### - [ ] Step 3: Edit MLB.R Phase 1 — compute home_scored_first on DT

In `Answer Keys/MLB.R`, locate line 60 (the closing `as.data.table()` of the DT pipeline). Immediately after that line, insert:

```r
# Compute home_scored_first indicator for triple-play prop pricing.
# Carried through sampling by preserving as an extra column on DT.
source("triple_play_helpers.R")
DT[, home_scored_first := determine_home_scored_first_vec(
  game_home_margin_inning_inning_1, game_total_inning_inning_1,
  game_home_margin_inning_inning_2, game_total_inning_inning_2,
  game_home_margin_inning_inning_3, game_total_inning_inning_3,
  game_home_margin_inning_inning_4, game_total_inning_inning_4,
  game_home_margin_inning_inning_5, game_total_inning_inning_5
)]
cat(sprintf("home_scored_first coverage: %.1f%% of games determinable in innings 1-5.\n",
            100 * mean(!is.na(DT$home_scored_first))))
```

### - [ ] Step 4: Edit MLB.R Phase 4 — add column to samples export

In `Answer Keys/MLB.R`, locate the `sample_rows <- imap_dfr(...)` block at lines 309–323. Replace it with:

```r
sample_rows <- imap_dfr(samples, function(s, game_id) {
  samp <- s$sample
  tibble(
    game_id            = game_id,
    sim_idx            = seq_len(nrow(samp)),
    home_margin        = samp$game_home_margin_period_FG,
    total_final_score  = samp$game_total_period_FG,
    home_margin_f3     = samp$game_home_margin_period_F3,
    total_f3           = samp$game_total_period_F3,
    home_margin_f5     = samp$game_home_margin_period_F5,
    total_f5           = samp$game_total_period_F5,
    home_margin_f7     = samp$game_home_margin_period_F7,
    total_f7           = samp$game_total_period_F7,
    home_scored_first  = samp$home_scored_first
  )
})
```

### - [ ] Step 5: Run MLB.R end-to-end from the worktree

```bash
cd "Answer Keys" && Rscript MLB.R
```
Expected: pipeline completes without error. Console should include a line like `home_scored_first coverage: 95.2% of games determinable in innings 1-5.`

### - [ ] Step 6: Run smoke test to confirm column present

```bash
cd "Answer Keys/tests" && Rscript -e 'library(testthat); test_file("test_triple_play.R")'
```
Expected: all tests including the smoke test pass.

### - [ ] Step 7: Commit

```bash
cd ../..
git add "Answer Keys/MLB.R" "Answer Keys/tests/test_triple_play.R"
git commit -m "feat(mlb): plumb home_scored_first through samples pipeline"
```

---

## Task 3: Triple-play fair-prob core function

**Files:**
- Create: `Answer Keys/mlb_triple_play.R`
- Modify: `Answer Keys/tests/test_triple_play.R` — add pricer tests

### - [ ] Step 1: Append pricer tests to test_triple_play.R

```r
# Append to Answer Keys/tests/test_triple_play.R
source("../mlb_triple_play.R", local = TRUE)  # safe: main block is guarded

test_that("compute_triple_play_fair (home) returns joint P of 3 legs", {
  # 10 synthetic rows. home_triple = scored_first==1 AND margin_f5>0 AND home_margin>0.
  # Row-by-row: 1:(1,1,2)✓  2:(1,2,3)✓  3:(0,0,-1)✗  4:(1,3,4)✓
  #             5:(1,0,1)✗f5tie  6:(1,2,5)✓  7:(0,-1,-2)✗
  #             8:(0,1,2)✗scored_first=0  9:(1,4,6)✓  10:(0,-2,-3)✗
  # home_triple hits: rows 1,2,4,6,9 → 5/10 = 0.5
  samples <- data.frame(
    home_margin       = c( 2,  3, -1,  4,  1,  5, -2,  2,  6, -3),
    home_margin_f5    = c( 1,  2,  0,  3,  0,  2, -1,  1,  4, -2),
    home_scored_first = c( 1L, 1L, 0L, 1L, 1L, 1L, 0L, 0L, 1L, 0L)
  )
  expect_equal(compute_triple_play_fair(samples, side = "home"), 0.5)
})

test_that("compute_triple_play_fair (away) returns joint P for away side", {
  # away_triple = scored_first==0 AND margin_f5<0 AND home_margin<0
  samples <- data.frame(
    home_margin       = c(-2, -3,  1, -4,  1, -5,  2, -2, -6,  3),
    home_margin_f5    = c(-1, -2,  0, -3,  0, -2,  1, -1, -4,  2),
    home_scored_first = c( 0L, 0L, 1L, 0L, 0L, 0L, 1L, 1L, 0L, 1L)
  )
  # hits: 1,2,4,6,9 → 5/10
  expect_equal(compute_triple_play_fair(samples, side = "away"), 0.5)
})

test_that("compute_triple_play_fair excludes NA home_scored_first rows", {
  # 6 rows — 3 valid with 1 hit, 3 NA. Fair should be 1/3, not 1/6.
  samples <- data.frame(
    home_margin       = c( 2,  3,  1,  4,  5,  6),
    home_margin_f5    = c( 1,  2,  0,  3,  4,  5),
    home_scored_first = c( 1L, 0L, 1L, NA, NA, NA)
  )
  expect_equal(compute_triple_play_fair(samples, side = "home"), 1/3)
})

test_that("prob_to_american handles favorites and dogs", {
  expect_equal(prob_to_american(0.5), -100L)
  expect_equal(prob_to_american(0.6), -150L)   # -100 * 0.6/0.4
  expect_equal(prob_to_american(0.25), 300L)   # 100 * 0.75/0.25
  expect_true(is.na(prob_to_american(0)))
  expect_true(is.na(prob_to_american(1)))
  expect_true(is.na(prob_to_american(NA_real_)))
})

test_that("american_to_prob is the inverse of prob_to_american", {
  expect_equal(american_to_prob(-150), 0.6)
  expect_equal(american_to_prob(+300), 0.25)
  expect_equal(american_to_prob(+100), 0.5)
})
```

### - [ ] Step 2: Run tests — should fail with function-not-found

```bash
cd "Answer Keys/tests" && Rscript -e 'library(testthat); test_file("test_triple_play.R")'
```
Expected: failures with `could not find function "compute_triple_play_fair"`, `prob_to_american`, etc.

### - [ ] Step 3: Create mlb_triple_play.R with core functions (no main block yet)

```r
#!/usr/bin/env Rscript
# MLB Triple-Play Pricer
#
# Prices "SCR 1ST, 1H & GM" props: team scores first in the game AND wins the
# 1st half (F5, strict lead) AND wins the game. Uses the dispersion-matched
# samples from mlb_game_samples — same framework as mlb_correlated_parlay.R.
#
# Usage: Rscript mlb_triple_play.R

suppressPackageStartupMessages({
  library(duckdb)
  library(dplyr)
  library(tibble)
  library(DBI)
})

# =============================================================================
# CORE PRICER FUNCTIONS (pure — unit-tested)
# =============================================================================

#' Fair probability for a team's triple-play, computed empirically on sample rows.
#'
#' @param samples data.frame with columns home_margin, home_margin_f5,
#'   home_scored_first (0/1/NA)
#' @param side "home" or "away"
#' @return scalar probability in [0, 1], or NA if no valid rows
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

#' Probability → American odds (integer). Returns NA for p <= 0 or p >= 1.
prob_to_american <- function(p) {
  if (is.na(p) || p <= 0 || p >= 1) return(NA_integer_)
  if (p >= 0.5) return(as.integer(round(-100 * p / (1 - p))))
  as.integer(round(100 * (1 - p) / p))
}

#' American odds → implied probability (includes vig).
american_to_prob <- function(o) {
  if (is.na(o)) return(NA_real_)
  if (o > 0) return(100 / (o + 100))
  (-o) / ((-o) + 100)
}

# Main block goes here in Task 4 — guarded so `source()` from tests is safe.
if (!interactive() && sys.nframe() == 0L) {
  # Task 4 will fill this in
  invisible(NULL)
}
```

### - [ ] Step 4: Run tests — all should pass

```bash
cd "Answer Keys/tests" && Rscript -e 'library(testthat); test_file("test_triple_play.R")'
```
Expected: all tests pass (helper, smoke, pricer, conversion functions).

### - [ ] Step 5: Commit

```bash
cd ../..
git add "Answer Keys/mlb_triple_play.R" "Answer Keys/tests/test_triple_play.R"
git commit -m "feat(mlb): add triple-play fair probability computation"
```

---

## Task 4: Wire today's lines and print the comparison table

**Files:**
- Modify: `Answer Keys/mlb_triple_play.R` — replace the main block

### - [ ] Step 1: Replace the guarded main block with the real script

Replace the `if (!interactive() && ...)` block at the bottom of `mlb_triple_play.R` with:

```r
# =============================================================================
# MAIN
# =============================================================================

setwd("~/NFLWork/Answer Keys")
MLB_DB <- "mlb.duckdb"

# Today's book lines — edit this tribble whenever new triple-plays post.
# home_team / away_team must match Odds API canonical names in mlb_consensus_temp.
todays_lines <- tribble(
  ~home_team,              ~away_team,             ~target_team, ~side,   ~book_odds,
  "San Diego Padres",      "Colorado Rockies",     "Padres",     "home",  +190,
  "San Diego Padres",      "Colorado Rockies",     "Rockies",    "away",  +530,
  "San Francisco Giants",  "Los Angeles Dodgers",  "Giants",     "home",  +750,
  "San Francisco Giants",  "Los Angeles Dodgers",  "Dodgers",    "away",  +155,
  "Seattle Mariners",      "Athletics",            "Mariners",   "home",  +215,
  "Seattle Mariners",      "Athletics",            "Athletics",  "away",  +455,
  "Arizona Diamondbacks",  "Chicago White Sox",    "DBacks",     "home",  +240,
  "Arizona Diamondbacks",  "Chicago White Sox",    "White Sox",  "away",  +415
)

con <- dbConnect(duckdb(), dbdir = MLB_DB, read_only = TRUE)
on.exit(tryCatch(dbDisconnect(con), error = function(e) NULL), add = TRUE)

samples_df <- dbGetQuery(con,
  "SELECT game_id, home_margin, home_margin_f5, home_scored_first
   FROM mlb_game_samples")
consensus  <- dbGetQuery(con,
  "SELECT id, home_team, away_team FROM mlb_consensus_temp")

if (!"home_scored_first" %in% names(samples_df)) {
  stop("mlb_game_samples is missing home_scored_first. Re-run MLB.R to regenerate.")
}

# Join today's lines to their game_id via team names
matched <- todays_lines %>%
  inner_join(consensus, by = c("home_team", "away_team"))

n_matched <- nrow(matched)
n_posted  <- nrow(todays_lines)
if (n_matched < n_posted) {
  dropped <- anti_join(todays_lines, consensus, by = c("home_team", "away_team"))
  warning(sprintf("Dropped %d/%d lines with no matching consensus game:\n%s",
                  n_posted - n_matched, n_posted,
                  paste(capture.output(print(dropped)), collapse = "\n")))
}

# Price each line
priced <- matched %>%
  rowwise() %>%
  mutate(
    game_samples = list(samples_df[samples_df$game_id == id, ]),
    n_samples    = nrow(game_samples),
    fair_prob    = compute_triple_play_fair(game_samples, side),
    fair_odds    = prob_to_american(fair_prob),
    book_prob    = american_to_prob(book_odds),
    edge_pct     = (fair_prob / book_prob - 1) * 100
  ) %>%
  ungroup() %>%
  select(target_team, side, n_samples, fair_prob, fair_odds, book_odds, edge_pct) %>%
  arrange(desc(edge_pct))

cat("\n=== MLB Triple-Play Fair Prices (SCR 1ST + F5 + GM) ===\n")
cat(sprintf("Matched %d / %d posted lines.\n\n", n_matched, n_posted))

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

dbDisconnect(con)
```

### - [ ] Step 2: Run the pricer end-to-end

```bash
cd "Answer Keys" && Rscript mlb_triple_play.R
```
Expected: 8-row table printed, sorted by edge descending, with fair_prob / fair_odds / book_odds / edge_pct populated. No warnings about dropped lines (all 4 games match).

### - [ ] Step 3: Sanity check against the earlier Python ML-bucket prototype

Ballpark expectations (dispersion-matched samples should be close but not identical to the ±3% ML-bucket prototype):

| Team | Prototype fair | Expected from pricer |
|---|---|---|
| Dodgers | +151 | +140 to +165 |
| Mariners | +204 | +190 to +220 |
| DBacks | +211 | +200 to +230 |
| Padres | +217 | +200 to +240 |
| Athletics | +338 | +310 to +370 |
| White Sox | +321 | +300 to +350 |
| Rockies | +312 | +290 to +340 |
| Giants | +394 | +370 to +430 |

If any fair number is wildly outside this range, investigate before committing.

### - [ ] Step 4: Commit

```bash
git add "Answer Keys/mlb_triple_play.R"
git commit -m "feat(mlb): print triple-play fair odds vs today's book lines"
```

---

## Task 5: Documentation

**Files:**
- Modify: `Answer Keys/CLAUDE.md` — add "Triple-Play Data Flow" section

### - [ ] Step 1: Insert the new section after "Race-to-X Data Flow"

In `Answer Keys/CLAUDE.md`, locate the end of the `### Race-to-X Data Flow` subsection (the line containing `Backward compat: \`first_to_10_h1\` still extracted alongside...`). After that line's closing, add:

````markdown

### Triple-Play Data Flow
Triple-play props = team scores first in the game AND wins F5 (3-way, strict lead) AND wins the game.
```
MLB.R Phase 1 (DT load)
  └── determine_home_scored_first_vec() → home_scored_first column on DT
      (baseball half-inning order: away bats first → away scored first if
       away_runs_in_first_scoring_inning > 0)
MLB.R Phase 4 (samples export)
  └── home_scored_first carried into mlb_game_samples as 0/1/NA column
      (preserved automatically because run_answer_key_sample returns a row
       subset of DT with all columns intact)
mlb_triple_play.R (standalone pricer)
  ├── Reads mlb_game_samples + hardcoded today's lines
  ├── compute_triple_play_fair(samples, side) → joint P across the 3 legs
  └── Prints fair odds vs book + edge per posted line
```
- `home_scored_first` is NA for games with no scoring in innings 1–5 (~5% of historical games); excluded from the mean before the ratio is computed.
- F5 3-way: leg passes only with strict lead (`margin_f5 > 0` for home, `< 0` for away). F5 ties kill the parlay.
- Helper lives in `triple_play_helpers.R` so both MLB.R and the pricer can source it without duplication.

````

### - [ ] Step 2: Commit

```bash
git add "Answer Keys/CLAUDE.md"
git commit -m "docs(mlb): document triple-play data flow"
```

---

## Self-Review

**Spec coverage (user said: "plan similar to correlated parlays, basic to start, just use today's lines"):**
- Mirrors `mlb_correlated_parlay.R` structure — standalone script, reads `mlb_game_samples`, computes joint prob on resampled rows ✓ (Tasks 3–4)
- Basic scope: no DuckDB persistence, no dashboard integration, no Kelly sizing ✓
- Today's 8 lines hardcoded in a tribble ✓ (Task 4)
- Uses the real answer-key sampling engine ✓ (Task 2 plumbs through samples)

**Placeholder scan:** No TBDs, no "implement later". All code complete.

**Type consistency:**
- `determine_home_scored_first` returns `integer` (0L/1L/NA_integer_); `home_scored_first := ...` stores integer; comparisons use `== 1L` / `== 0L` ✓
- Sample column names `home_margin`, `home_margin_f5`, `home_scored_first` match across export (Task 2), tests (Task 3), and main block (Task 4) ✓
- `side` values `"home"` / `"away"` consistent everywhere ✓

**Worktree / version control / docs / pre-merge review:** all present per CLAUDE.md requirements.

---

## Execution Handoff

**Plan saved to `docs/superpowers/plans/2026-04-21-mlb-triple-play-pricer.md`. Two execution options:**

**1. Subagent-Driven (recommended)** — fresh subagent per task, review between tasks, fast iteration.

**2. Inline Execution** — execute tasks in this session with checkpoints for review.

**Which approach?**
