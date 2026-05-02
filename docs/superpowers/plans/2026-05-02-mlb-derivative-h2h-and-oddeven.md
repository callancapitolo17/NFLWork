# MLB Derivative h2h + Odd/Even Matching Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add straight-bet matching for two market families that are scraped today but silently dropped on the MLB dashboard: F3/F7 moneyline (`h2h_1st_3_innings`, `h2h_1st_7_innings`) and full-game odd/even total runs (`odd_even_runs`). This is "Bucket A" of the broader MLB matching expansion (gap #3 + odd/even from the trace).

**Architecture:** Two new `else if` branches inside `compare_alts_to_samples()` in `Answer Keys/Tools.R`. The h2h branch uses the existing `game_home_margin_period_F3` / `game_home_margin_period_F7` sample columns (already populated by `MLB.R` Phase 1 DT prep). The odd/even branch uses parity of `game_total_period_FG`. No new sample-data work, no MLB.R changes, no PBP joins, no team-name dictionary additions.

**Tech Stack:** R, testthat, DuckDB

**Branch:** `feature/mlb-h2h-derivative-and-oddeven`
**Worktree:** `/Users/callancapitolo/NFLWork/.worktrees/mlb-h2h-derivative-and-oddeven`

---

## Background — what's broken today

`compare_alts_to_samples()` in `Tools.R:4504` already filters offshore-odds rows to derivatives:

```r
alt_odds <- offshore_odds %>%
  filter(
    grepl("^alternate_", market) |
    grepl("^team_totals_", market) |
    grepl("1st_3_innings", market) |
    grepl("1st_7_innings", market)
  )
```

So `h2h_1st_3_innings` and `h2h_1st_7_innings` rows pass the filter. But the per-row `if/else if` chain only handles `spread`, non-team `total`, and `team_totals_*`. h2h falls through silently and no bet is emitted. `odd_even_runs` doesn't even pass the filter.

The Wagerzon scraper at `wagerzon_odds/scraper_v2.py:365-379, 446-453` already produces both market types — they just have no consumer.

---

## File Structure

| File | Role | Touch |
|------|------|-------|
| `Answer Keys/Tools.R` | Hosts `compare_alts_to_samples()` (and ~150 other helpers) — single source of truth for all matching logic | **Modify** — extend filter regex + add 2 new branches in the per-row chain (~80 lines added inside one function) |
| `Answer Keys/tests/test_compare_alts_to_samples.R` | New unit-test file using `testthat`, mirrors the convention from `tests/test_books_strip.R`. Builds tiny synthetic `samples` list + `offshore_odds` tibble + `consensus_odds` tibble, then asserts `compare_alts_to_samples()` emits the expected bet rows. | **Create** |
| `Answer Keys/CLAUDE.md` | Pipeline overview — currently says "F5 (first 5 innings) markets only" | **Modify** — extend the bullet to note F3/F7 h2h + odd/even are now matched too |
| `Answer Keys/MLB Dashboard/README.md` | Markets section currently lists only F5 | **Modify** — add `h2h_1st_3_innings`, `h2h_1st_7_innings`, `odd_even_runs` |

No changes to `MLB.R`, no changes to scrapers, no changes to dashboard JS or HTML. The dashboard's market filter is built dynamically from `mlb_bets_combined.market` — new markets show up automatically.

---

## Worktree Section

- **Worktree:** `.worktrees/mlb-h2h-derivative-and-oddeven` (already created)
- **Branch:** `feature/mlb-h2h-derivative-and-oddeven` (already created off `main`)
- **DuckDB rule:** Test data must NOT use symlinks. The end-to-end smoke test (Task 5) copies `mlb.duckdb`, `mlb_mm.duckdb`, `wagerzon.duckdb` into the worktree before running `python run.py mlb` — never symlinks them. Per project rule (`/Users/callancapitolo/NFLWork/CLAUDE.md` housekeeping #5), symlinks lose WAL data on worktree removal.
- **Cleanup after merge:** `git worktree remove .worktrees/mlb-h2h-derivative-and-oddeven && git branch -d feature/mlb-h2h-derivative-and-oddeven`

---

## Version Control Section

- All commits land on `feature/mlb-h2h-derivative-and-oddeven`.
- Commit cadence: one commit per task (5 commits total).
- Final merge: `git checkout main && git merge --no-ff feature/mlb-h2h-derivative-and-oddeven` — only after explicit user approval.
- Never use `--no-verify`. No force-push.

Commit messages:
1. `test(mlb): add fixture-driven tests for compare_alts_to_samples`
2. `feat(mlb): match F3/F7 derivative moneylines in compare_alts_to_samples`
3. `feat(mlb): match odd_even_runs in compare_alts_to_samples`
4. `test(mlb): end-to-end smoke for new derivative + odd/even bets`
5. `docs(mlb): note F3/F7 h2h and odd/even matching in CLAUDE.md + dashboard README`

---

## Documentation Section

After code is finalized and reviewed:

- **`Answer Keys/CLAUDE.md`** — under "Pipeline Flow (MLB)", change the bullet `F5 (first 5 innings) markets only: h2h, totals, spreads` to also list `+ F3/F7 h2h via compare_alts_to_samples; FG odd/even runs via compare_alts_to_samples`.
- **`Answer Keys/MLB Dashboard/README.md`** — under the "Markets" section, append the three new market keys.
- No `MLB Answer Key/README.md` change needed (it's user-facing setup, not algorithmic detail).

---

## Tasks

### Task 1: Test scaffolding for compare_alts_to_samples (synthetic fixture)

**Files:**
- Create: `Answer Keys/tests/test_compare_alts_to_samples.R`

**Why this task first:** Establishes a synthetic samples + odds fixture that Tasks 2-3 reuse. Without this, the only way to exercise the new branches is via end-to-end pipeline runs, which are slow and dependent on whatever Wagerzon happens to be posting today. A fixture decouples the unit logic from live data.

- [ ] **Step 1: Write the test file with a fixture builder**

```r
# Answer Keys/tests/test_compare_alts_to_samples.R
# Fixture-driven tests for compare_alts_to_samples derivative-h2h + odd/even branches.
# Run from "Answer Keys/" directory: Rscript -e 'testthat::test_file("tests/test_compare_alts_to_samples.R")'

library(testthat)
library(dplyr)
library(tibble)
library(data.table)
source("Tools.R")

# Helper: build a synthetic samples list with one game.
# margins_f3 / margins_f7 / totals_fg are vectors of 1000 sim outcomes.
make_synthetic_samples <- function(
    game_id = "test-game-1",
    margins_f3 = NULL,
    margins_f7 = NULL,
    totals_fg = NULL
) {
  n <- 1000L
  if (is.null(margins_f3)) margins_f3 <- rep(0L, n)
  if (is.null(margins_f7)) margins_f7 <- rep(0L, n)
  if (is.null(totals_fg)) totals_fg <- rep(8L, n)
  stopifnot(length(margins_f3) == n, length(margins_f7) == n, length(totals_fg) == n)
  samp <- data.frame(
    game_home_margin_period_F3 = margins_f3,
    game_total_period_F3       = abs(margins_f3) + 6L,  # plausible total
    game_home_margin_period_F5 = margins_f3,            # placeholder
    game_total_period_F5       = abs(margins_f3) + 7L,
    game_home_margin_period_F7 = margins_f7,
    game_total_period_F7       = abs(margins_f7) + 8L,
    game_home_margin_period_FG = margins_f7,
    game_total_period_FG       = totals_fg
  )
  list(setNames(list(list(sample = samp)), game_id))[[1]]
}

make_consensus <- function(game_id = "test-game-1") {
  tibble(
    id = game_id,
    home_team = "Test Home",
    away_team = "Test Away",
    commence_time = as.POSIXct("2026-05-02 19:00:00", tz = "UTC")
  )
}
```

- [ ] **Step 2: Add the first failing test (h2h F3 fair-prob)**

Append to the same file:

```r
test_that("compare_alts_to_samples emits h2h_1st_3_innings bets when home wins F3 60% of the time", {
  # 600 home wins (margin > 0), 300 away wins (margin < 0), 100 ties
  margins <- c(rep(2L, 600), rep(-2L, 300), rep(0L, 100))
  samples <- make_synthetic_samples(margins_f3 = margins)
  consensus <- make_consensus()

  # Wagerzon-shaped row: away_ml -110 / home_ml -110 (so devigged fair = 50/50)
  # Model says home wins 60% of non-push → big edge on home
  offshore <- tibble(
    bookmaker_key = "wagerzon",
    home_team = "Test Home",
    away_team = "Test Away",
    market = "h2h_1st_3_innings",
    line = NA_real_,
    odds_away = -110L,
    odds_home = -110L,
    odds_over = NA_integer_,
    odds_under = NA_integer_,
    home_spread = NA_real_,
    away_spread = NA_real_,
    commence_time = as.POSIXct("2026-05-02 19:00:00", tz = "UTC")
  )

  bets <- compare_alts_to_samples(
    samples = samples,
    offshore_odds = offshore,
    consensus_odds = consensus,
    bankroll = 100,
    kelly_mult = 0.25,
    ev_threshold = 0.02
  )

  expect_true(nrow(bets) >= 1)
  home_bet <- bets[bets$bet_on == "Test Home", ]
  expect_equal(nrow(home_bet), 1)
  # 600/(600+300) = 0.6667 — exclude pushes from denominator
  expect_equal(home_bet$prob, 600 / 900, tolerance = 1e-9)
  expect_equal(home_bet$market, "h2h_1st_3_innings")
})
```

- [ ] **Step 3: Run the test to verify it fails for the right reason**

Run from the worktree's `Answer Keys/` directory:

```bash
cd "Answer Keys" && Rscript -e 'testthat::test_file("tests/test_compare_alts_to_samples.R")'
```

Expected: 1 FAIL on the first test. Failure message will be either "nrow(bets) is 0" (today's behavior — h2h falls through silently) or an error about an unhandled regex match. Both are acceptable failure modes — both confirm the branch doesn't exist yet.

- [ ] **Step 4: Commit**

```bash
git add "Answer Keys/tests/test_compare_alts_to_samples.R"
git commit -m "test(mlb): add fixture-driven tests for compare_alts_to_samples"
```

---

### Task 2: Implement h2h derivative branch

**Files:**
- Modify: `Answer Keys/Tools.R:4691` — add a new `else if` inside the per-row chain inside `compare_alts_to_samples()`. The existing chain is `if (grepl("spread", row$market) ...) {} else if (grepl("total", row$market) ...) {} else if (grepl("team_totals_", row$market) ...) {}`. The h2h branch slots in as a fourth `else if`.

**Why this insertion point:** The existing `spread` branch (`Tools.R:4565`) already loads `margins <- sample_df[[col_name]]` from the period-mapped column. The h2h branch can mirror that exactly — same column lookup, just sign instead of comparison-to-spread.

- [ ] **Step 1: Read the current end of the per-row chain to confirm exact insertion point**

Run:
```bash
sed -n '4685,4695p' "Answer Keys/Tools.R"
```

Expected output: the closing `}` of the `team_totals_` branch followed by `}` for the `for` loop. The new `else if` goes BEFORE the closing `}` of the `for` loop (i.e., immediately after the `team_totals_` branch's closing `}`).

- [ ] **Step 2: Add the h2h branch**

Insert after the closing `}` of the `team_totals_` branch in `Tools.R` (currently line 4691, but verify with the grep above — line numbers shift):

```r
    } else if (grepl("^h2h_", row$market) && !is.na(row$odds_home) && !is.na(row$odds_away)) {
      # Derivative moneyline (e.g. h2h_1st_3_innings, h2h_1st_7_innings).
      # 2-way pricing: home wins iff margin > 0, away wins iff margin < 0,
      # margin == 0 is a push (refund) and excluded from the denominator —
      # mirrors how Wagerzon settles F-period MLs in practice.
      col_name <- paste0(margin_col, "_", period)
      if (!col_name %in% names(sample_df)) next
      margins <- sample_df[[col_name]]
      margins <- margins[!is.na(margins)]
      if (length(margins) == 0) next

      non_push <- margins[margins != 0]
      if (length(non_push) == 0) next
      p_home <- sum(non_push > 0) / length(non_push)
      p_away <- 1 - p_home

      probs <- american_prob(row$odds_away, row$odds_home)
      if (any(is.na(probs)) || any(probs == 0)) next

      home_ev <- compute_ev(p_home, probs$p2)
      away_ev <- compute_ev(p_away, probs$p1)
      home_size <- kelly_stake(home_ev, probs$p2, bankroll, kelly_mult)
      away_size <- kelly_stake(away_ev, probs$p1, bankroll, kelly_mult)

      if (home_ev >= ev_threshold) {
        all_bets[[length(all_bets) + 1]] <- tibble(
          id = game_id, home_team = row$home_team, away_team = row$away_team,
          pt_start_time = pt_start_time, bookmaker_key = book_key,
          market = row$market, bet_on = row$home_team,
          line = NA_real_, bet_size = home_size, ev = home_ev,
          odds = row$odds_home, prob = p_home
        )
      }
      if (away_ev >= ev_threshold) {
        all_bets[[length(all_bets) + 1]] <- tibble(
          id = game_id, home_team = row$home_team, away_team = row$away_team,
          pt_start_time = pt_start_time, bookmaker_key = book_key,
          market = row$market, bet_on = row$away_team,
          line = NA_real_, bet_size = away_size, ev = away_ev,
          odds = row$odds_away, prob = p_away
        )
      }
```

(The `american_prob`, `compute_ev`, `kelly_stake` helpers and the `all_bets`, `book_key`, `pt_start_time`, `period` variables are already in scope inside the loop — same pattern as the spread branch above.)

- [ ] **Step 3: Run the test from Task 1 to verify it passes**

```bash
cd "Answer Keys" && Rscript -e 'testthat::test_file("tests/test_compare_alts_to_samples.R")'
```

Expected: PASS on the first test. The fair prob should be exactly `600/900 = 0.6667`.

- [ ] **Step 4: Add a second test — F7 case**

Append to `tests/test_compare_alts_to_samples.R`:

```r
test_that("compare_alts_to_samples emits h2h_1st_7_innings bets using F7 margins", {
  # 700 home wins, 250 away wins, 50 ties — at F7
  margins <- c(rep(3L, 700), rep(-3L, 250), rep(0L, 50))
  samples <- make_synthetic_samples(margins_f7 = margins)
  consensus <- make_consensus()

  offshore <- tibble(
    bookmaker_key = "wagerzon",
    home_team = "Test Home",
    away_team = "Test Away",
    market = "h2h_1st_7_innings",
    line = NA_real_,
    odds_away = +200L,   # implied 33.3% — way too long if home is 73.7%
    odds_home = -150L,   # implied 60.0% — undervaluing home
    odds_over = NA_integer_,
    odds_under = NA_integer_,
    home_spread = NA_real_,
    away_spread = NA_real_,
    commence_time = as.POSIXct("2026-05-02 19:00:00", tz = "UTC")
  )

  bets <- compare_alts_to_samples(
    samples = samples,
    offshore_odds = offshore,
    consensus_odds = consensus,
    bankroll = 100,
    kelly_mult = 0.25,
    ev_threshold = 0.02
  )

  home_bet <- bets[bets$bet_on == "Test Home" & bets$market == "h2h_1st_7_innings", ]
  expect_equal(nrow(home_bet), 1)
  expect_equal(home_bet$prob, 700 / 950, tolerance = 1e-9)
})
```

- [ ] **Step 5: Add a third test — degenerate-input safety**

Append:

```r
test_that("compare_alts_to_samples returns no bets when h2h row has NA odds", {
  margins <- c(rep(2L, 600), rep(-2L, 400))
  samples <- make_synthetic_samples(margins_f3 = margins)
  consensus <- make_consensus()

  offshore <- tibble(
    bookmaker_key = "wagerzon",
    home_team = "Test Home",
    away_team = "Test Away",
    market = "h2h_1st_3_innings",
    line = NA_real_,
    odds_away = NA_integer_,   # missing odds — must skip cleanly
    odds_home = NA_integer_,
    odds_over = NA_integer_,
    odds_under = NA_integer_,
    home_spread = NA_real_,
    away_spread = NA_real_,
    commence_time = as.POSIXct("2026-05-02 19:00:00", tz = "UTC")
  )

  bets <- compare_alts_to_samples(
    samples = samples, offshore_odds = offshore, consensus_odds = consensus,
    bankroll = 100, kelly_mult = 0.25, ev_threshold = 0.02
  )

  expect_equal(nrow(bets), 0)
})
```

- [ ] **Step 6: Run all three tests**

```bash
cd "Answer Keys" && Rscript -e 'testthat::test_file("tests/test_compare_alts_to_samples.R")'
```

Expected: 3 PASS, 0 FAIL.

- [ ] **Step 7: Commit**

```bash
git add "Answer Keys/Tools.R" "Answer Keys/tests/test_compare_alts_to_samples.R"
git commit -m "feat(mlb): match F3/F7 derivative moneylines in compare_alts_to_samples"
```

---

### Task 3: Implement odd/even runs branch

**Files:**
- Modify: `Answer Keys/Tools.R:4521-4525` — extend the filter regex
- Modify: `Answer Keys/Tools.R` — add a new `else if` branch immediately after the h2h branch from Task 2

**Why both edits are needed:** The current filter regex doesn't match `odd_even_runs`, so rows are dropped before the per-row chain. We must add it to the filter AND handle it in the chain. Bundling these two edits in one task keeps the filter-vs-handler invariant intact.

- [ ] **Step 1: Add the failing odd/even test**

Append to `Answer Keys/tests/test_compare_alts_to_samples.R`:

```r
test_that("compare_alts_to_samples emits odd_even_runs bet when totals are odd 60% of the time", {
  totals <- c(rep(7L, 600), rep(8L, 400))   # 60% odd, 40% even
  samples <- make_synthetic_samples(totals_fg = totals)
  consensus <- make_consensus()

  # Wagerzon convention: away_ml = ODD price, home_ml = EVEN price
  # (scraper_v2.py:446-453 — vtm = "TOTAL RUNS ODD", htm = "TOTAL RUNS EVEN")
  offshore <- tibble(
    bookmaker_key = "wagerzon",
    home_team = "Test Home",
    away_team = "Test Away",
    market = "odd_even_runs",
    line = NA_real_,
    odds_away = -110L,   # ODD — fair if 52.4%
    odds_home = -110L,   # EVEN — fair if 52.4%
    odds_over = NA_integer_,
    odds_under = NA_integer_,
    home_spread = NA_real_,
    away_spread = NA_real_,
    commence_time = as.POSIXct("2026-05-02 19:00:00", tz = "UTC")
  )

  bets <- compare_alts_to_samples(
    samples = samples, offshore_odds = offshore, consensus_odds = consensus,
    bankroll = 100, kelly_mult = 0.25, ev_threshold = 0.02
  )

  odd_bet <- bets[bets$bet_on == "Odd" & bets$market == "odd_even_runs", ]
  expect_equal(nrow(odd_bet), 1)
  expect_equal(odd_bet$prob, 0.60, tolerance = 1e-9)
})
```

- [ ] **Step 2: Run it to verify it fails**

```bash
cd "Answer Keys" && Rscript -e 'testthat::test_file("tests/test_compare_alts_to_samples.R")'
```

Expected: 1 FAIL ("nrow(odd_bet) is 0") because the filter regex drops `odd_even_runs` rows.

- [ ] **Step 3: Extend the filter regex**

In `Tools.R:4519-4525`, change:

```r
  alt_odds <- offshore_odds %>%
    filter(
      grepl("^alternate_", market) |
      grepl("^team_totals_", market) |
      grepl("1st_3_innings", market) |
      grepl("1st_7_innings", market)
    )
```

to:

```r
  alt_odds <- offshore_odds %>%
    filter(
      grepl("^alternate_", market) |
      grepl("^team_totals_", market) |
      grepl("1st_3_innings", market) |
      grepl("1st_7_innings", market) |
      market == "odd_even_runs"
    )
```

- [ ] **Step 4: Add the odd/even handler branch**

Insert the new `else if` immediately after the h2h branch from Task 2 (and before the closing `}` of the `for` loop):

```r
    } else if (row$market == "odd_even_runs" && !is.na(row$odds_home) && !is.na(row$odds_away)) {
      # Wagerzon convention (scraper_v2.py:446-453):
      #   away side = ODD total runs, home side = EVEN total runs
      # Pricing: parity of game_total_period_FG. Baseball totals are integers,
      # so no push case to handle.
      col_name <- paste0(total_col, "_FG")
      if (!col_name %in% names(sample_df)) next
      totals <- sample_df[[col_name]]
      totals <- totals[!is.na(totals)]
      if (length(totals) == 0) next

      p_odd <- sum(totals %% 2 == 1) / length(totals)
      p_even <- 1 - p_odd

      probs <- american_prob(row$odds_away, row$odds_home)
      if (any(is.na(probs)) || any(probs == 0)) next

      odd_ev  <- compute_ev(p_odd,  probs$p1)
      even_ev <- compute_ev(p_even, probs$p2)
      odd_size  <- kelly_stake(odd_ev,  probs$p1, bankroll, kelly_mult)
      even_size <- kelly_stake(even_ev, probs$p2, bankroll, kelly_mult)

      if (odd_ev >= ev_threshold) {
        all_bets[[length(all_bets) + 1]] <- tibble(
          id = game_id, home_team = row$home_team, away_team = row$away_team,
          pt_start_time = pt_start_time, bookmaker_key = book_key,
          market = row$market, bet_on = "Odd",
          line = NA_real_, bet_size = odd_size, ev = odd_ev,
          odds = row$odds_away, prob = p_odd
        )
      }
      if (even_ev >= ev_threshold) {
        all_bets[[length(all_bets) + 1]] <- tibble(
          id = game_id, home_team = row$home_team, away_team = row$away_team,
          pt_start_time = pt_start_time, bookmaker_key = book_key,
          market = row$market, bet_on = "Even",
          line = NA_real_, bet_size = even_size, ev = even_ev,
          odds = row$odds_home, prob = p_even
        )
      }
```

- [ ] **Step 5: Run all tests**

```bash
cd "Answer Keys" && Rscript -e 'testthat::test_file("tests/test_compare_alts_to_samples.R")'
```

Expected: 4 PASS, 0 FAIL.

- [ ] **Step 6: Commit**

```bash
git add "Answer Keys/Tools.R" "Answer Keys/tests/test_compare_alts_to_samples.R"
git commit -m "feat(mlb): match odd_even_runs in compare_alts_to_samples"
```

---

### Task 4: End-to-end smoke test against copied production DBs

**Why this task exists:** Unit tests prove the new branches behave correctly on synthetic input. This task proves they also fire correctly when called inside the real MLB pipeline against real Wagerzon data — catches integration bugs the unit tests can't (team-name resolution, `nearest_game_match` joining, `map_scraper_markets_mlb` interaction, `mlb_bets_combined` insert).

**Files (test artifacts only, never committed):**
- Copy: `mlb.duckdb`, `mlb_mm.duckdb`, the live scraper DBs into the worktree

- [ ] **Step 1: Copy live DuckDB files into the worktree (NEVER symlink)**

```bash
cp /Users/callancapitolo/NFLWork/Answer\ Keys/mlb.duckdb \
   "/Users/callancapitolo/NFLWork/.worktrees/mlb-h2h-derivative-and-oddeven/Answer Keys/mlb.duckdb"
cp /Users/callancapitolo/NFLWork/Answer\ Keys/mlb_mm.duckdb \
   "/Users/callancapitolo/NFLWork/.worktrees/mlb-h2h-derivative-and-oddeven/Answer Keys/mlb_mm.duckdb"
cp /Users/callancapitolo/NFLWork/Answer\ Keys/pbp.duckdb \
   "/Users/callancapitolo/NFLWork/.worktrees/mlb-h2h-derivative-and-oddeven/Answer Keys/pbp.duckdb"
cp /Users/callancapitolo/NFLWork/wagerzon_odds/wagerzon.duckdb \
   "/Users/callancapitolo/NFLWork/.worktrees/mlb-h2h-derivative-and-oddeven/wagerzon_odds/wagerzon.duckdb"
cp /Users/callancapitolo/NFLWork/bookmaker_odds/bookmaker.duckdb \
   "/Users/callancapitolo/NFLWork/.worktrees/mlb-h2h-derivative-and-oddeven/bookmaker_odds/bookmaker.duckdb"
cp /Users/callancapitolo/NFLWork/bet105_odds/bet105.duckdb \
   "/Users/callancapitolo/NFLWork/.worktrees/mlb-h2h-derivative-and-oddeven/bet105_odds/bet105.duckdb"
cp /Users/callancapitolo/NFLWork/bfa_odds/bfa.duckdb \
   "/Users/callancapitolo/NFLWork/.worktrees/mlb-h2h-derivative-and-oddeven/bfa_odds/bfa.duckdb"
```

(These DBs are gitignored — no commit needed. They're throwaway test fixtures.)

- [ ] **Step 2: Run the MLB pipeline end-to-end on the worktree's copies**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/mlb-h2h-derivative-and-oddeven/Answer Keys"
python3 run.py mlb 2>&1 | tee /tmp/mlb_pipeline_smoke.log
```

Expected: exit 0. Look for these new log lines from `compare_alts_to_samples`:
- `Generated N wagerzon alt/team-total bets from samples` where N includes new h2h_1st_3/7_innings rows when edge exists.

If Wagerzon happens to not be posting `odd_even_runs` at scrape time (often the case before evening — see session note), the odd/even portion of the smoke test won't have rows to exercise. The unit tests in Task 3 cover that path; this is acceptable.

- [ ] **Step 3: Inspect mlb_bets_combined for new market types**

```bash
python3 -c "
import duckdb
con = duckdb.connect('/Users/callancapitolo/NFLWork/.worktrees/mlb-h2h-derivative-and-oddeven/Answer Keys/mlb.duckdb', read_only=True)
print(con.execute('SELECT market, bookmaker_key, COUNT(*) n FROM mlb_bets_combined GROUP BY 1,2 ORDER BY 1,2').fetchdf().to_string())
"
```

Expected: presence of one or more of `h2h_1st_3_innings`, `h2h_1st_7_innings`, `odd_even_runs` from `wagerzon` or `bookmaker` IF edges exist for tonight's slate. Absence of those rows is a SOFT failure — could just mean no edge today. To distinguish hard from soft failure:

```bash
# Are there any spreads_f3 / spreads_f7 wagerzon rows in the raw scraper data?
python3 -c "
import duckdb
con = duckdb.connect('/Users/callancapitolo/NFLWork/.worktrees/mlb-h2h-derivative-and-oddeven/wagerzon_odds/wagerzon.duckdb', read_only=True)
print(con.execute(\"SELECT market, COUNT(*) n FROM mlb_odds WHERE market IN ('spreads_f3','spreads_f7','odd_even_runs') GROUP BY 1\").fetchdf().to_string())
"
```

If raw scraper rows exist but no matched bets emerge AND there's no `nearest_game_match` warning in the pipeline log → hard failure, debug. If raw rows exist and matched bets emerge → integration confirmed.

- [ ] **Step 4: Spot-check one matched bet for sanity**

If the smoke produced any `h2h_1st_3_innings` row, pick one and hand-verify. Example query:

```bash
python3 -c "
import duckdb
con = duckdb.connect('/Users/callancapitolo/NFLWork/.worktrees/mlb-h2h-derivative-and-oddeven/Answer Keys/mlb.duckdb', read_only=True)
df = con.execute(\"\"\"
  SELECT id, home_team, away_team, market, bet_on, odds, prob, ev, bet_size
  FROM mlb_bets_combined
  WHERE market LIKE 'h2h_1st_%' OR market = 'odd_even_runs'
  ORDER BY ev DESC
  LIMIT 5
\"\"\").fetchdf()
print(df.to_string())
"
```

Cross-check one row by hand:
- `prob` should be in [0, 1]
- For h2h: prob ≈ (sample fraction with margin > 0 at that period, excluding ties)
- `odds` should match the raw Wagerzon row for that game and side
- `ev = (prob × profit) − ((1 − prob) × stake)` per the project's `compute_ev` definition

If any sanity check fails, debug before the docs commit.

- [ ] **Step 5: Commit nothing (DBs are gitignored — this task produces no committed artifacts)**

Confirm cleanliness:

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/mlb-h2h-derivative-and-oddeven && git status
```

Expected output: `nothing to commit, working tree clean` (the copied .duckdb files don't show because they're gitignored).

---

### Task 5: Documentation updates

**Files:**
- Modify: `Answer Keys/CLAUDE.md` — Pipeline Flow (MLB) section
- Modify: `Answer Keys/MLB Dashboard/README.md` — Markets section

- [ ] **Step 1: Update `Answer Keys/CLAUDE.md`**

Find the MLB pipeline bullets (around the line `- F5 (first 5 innings) markets only: h2h, totals, spreads`) and change that bullet to:

```markdown
- F5 (first 5 innings) markets: h2h, totals, spreads
- Derivative markets matched via `compare_alts_to_samples` in `Tools.R`:
  - F3 / F7 spreads + totals + h2h (Wagerzon, Bookmaker)
  - FG alt spreads + alt totals (Wagerzon, Bet105)
  - FG odd/even total runs (Wagerzon)
- Historical data: `pbp.duckdb/mlb_betting_pbp` (12,719 games)
```

- [ ] **Step 2: Update `Answer Keys/MLB Dashboard/README.md`**

Find the "Markets" section (around `Currently F5 (first 5 innings) only:`) and replace with:

```markdown
## Markets

F5 mains (priced via `build_*_from_samples` against the Odds API):
- `h2h_1st_5_innings` — F5 moneyline
- `totals_1st_5_innings` — F5 total + alternate totals
- `spreads_1st_5_innings` — F5 run line

Derivatives (priced via `compare_alts_to_samples` against scraped Wagerzon/Bookmaker/BFA/Bet105 data):
- `h2h_1st_3_innings`, `h2h_1st_7_innings` — F3/F7 moneyline (Wagerzon, Bookmaker for F3 only)
- `spreads_1st_3_innings`, `spreads_1st_7_innings` — F3/F7 run line
- `totals_1st_3_innings`, `totals_1st_7_innings` — F3/F7 total
- `alternate_spreads_fg`, `alternate_totals_fg` — full-game alt lines
- `alternate_totals_f5`, `alternate_spreads_f5` — F5 alt lines (h1-suffix scrapers remap to f5)
- `odd_even_runs` — full-game total runs odd vs even (Wagerzon only)
```

- [ ] **Step 3: Commit docs**

```bash
git add "Answer Keys/CLAUDE.md" "Answer Keys/MLB Dashboard/README.md"
git commit -m "docs(mlb): note F3/F7 h2h and odd/even matching in CLAUDE.md + dashboard README"
```

---

### Task 6: Pre-merge review (REQUIRED before merge)

Per `/Users/callancapitolo/NFLWork/CLAUDE.md` "Pre-merge review" section — must complete before requesting merge approval.

- [ ] **Step 1: Review the full diff**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/mlb-h2h-derivative-and-oddeven && git diff main..HEAD --stat
git diff main..HEAD
```

- [ ] **Step 2: Walk the executive checklist**

| Check | Status |
|-------|--------|
| Data integrity — no duplicate writes; pushes excluded from h2h denominator (margin == 0); odd/even has no push case (integer totals) | verify |
| Resource safety — no new DB connections opened (we only modify pure functions) | verify |
| Edge cases — h2h with all-tie sample, h2h with NA odds, odd/even with all-same-parity totals, missing sample columns | covered by tests |
| Dead code — no leftover comments, no commented-out blocks, no unused parameters introduced | verify |
| Log/disk hygiene — `compare_alts_to_samples` already prints a summary line; no new logging that could grow unbounded | verify |
| Security — no new env vars, no new credentials, no new external calls | verify |

- [ ] **Step 3: Document the review summary in the PR/merge message**

Compose the merge commit body to include the issues-found-and-fixed list (or "no issues found").

- [ ] **Step 4: Surface to user with the diff and ask for explicit approval to merge to `main`**

Per project rule and user preference — never merge without explicit approval. Even after all tests pass.

---

### Task 7: Merge + cleanup (only after user approves)

- [ ] **Step 1: Merge to main**

```bash
cd /Users/callancapitolo/NFLWork
git checkout main
git merge --no-ff feature/mlb-h2h-derivative-and-oddeven -m "Merge feature/mlb-h2h-derivative-and-oddeven: F3/F7 h2h + odd/even matching"
```

- [ ] **Step 2: Re-run unit tests on main to confirm clean merge**

```bash
cd "/Users/callancapitolo/NFLWork/Answer Keys" && Rscript -e 'testthat::test_file("tests/test_compare_alts_to_samples.R")'
```

Expected: 4 PASS, 0 FAIL.

- [ ] **Step 3: Clean up worktree + branch**

```bash
git worktree remove .worktrees/mlb-h2h-derivative-and-oddeven
git branch -d feature/mlb-h2h-derivative-and-oddeven
```

- [ ] **Step 4: Confirm worktree list is clean**

```bash
git worktree list
```

Expected: only `main` worktree remains.

---

## Self-Review (post-write checklist)

- [x] **Spec coverage** — All three goal markets (h2h_1st_3_innings, h2h_1st_7_innings, odd_even_runs) have a task that implements them. Tasks 2 + 3.
- [x] **Placeholder scan** — No "TBD", "TODO", "implement later", "similar to Task N", "add appropriate error handling". Every code block is complete.
- [x] **Type consistency** — `prob`, `bet_size`, `ev`, `odds`, `bet_on`, `line` columns match the existing tibble schema returned by the spread/total/team_totals branches at `Tools.R:4589, 4630, 4674`. New rows slot into the same `bind_rows()` collection.
- [x] **No new dependencies** — Reuses `american_prob`, `compute_ev`, `kelly_stake`, `nearest_game_match` already in `Tools.R`. No new packages.
- [x] **DRY** — h2h branch duplicates ~5 lines of column-lookup boilerplate from the spread branch. Acceptable given the variable bindings differ (margin vs. spread comparison) and extracting a helper would obscure the per-branch math. The user-facing review pattern is "read all four branches side-by-side"; mechanical similarity helps that.
- [x] **YAGNI** — No new sample columns, no new helper functions, no new MLB.R changes, no new dashboard JS, no new config knobs.
- [x] **TDD** — Tests come before implementation in every task.
- [x] **Frequent commits** — 5 commits, one per logical unit.
