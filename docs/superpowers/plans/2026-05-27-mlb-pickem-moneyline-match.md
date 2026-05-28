# MLB Dashboard — Pick'em ↔ Moneyline Matching Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** On the MLB bets-tab odds screen, a bet that is a spread at line 0 (a pick'em / draw-no-bet) shows each book's true pick'em price — derived from the book's 2-way period winner, or from a 3-way winner (devig, drop tie, renormalize) — instead of the wrong-bet ±0.5 run line fallback. Display-only; the model is untouched.

**Architecture:** Approach A (per spec). One pure helper `derive_pickem_american()` in `odds_screen.R` does the DNB math (2-way direct, 3-way devig-drop-tie). The matcher (`expand_bets_to_book_prices`) gains a `line == 0` branch that looks up each book's period winner rows and calls the helper. A new `derived_fair_odds` column on `mlb_bets_book_prices` carries the precomputed FAIR; `book_cell.R` displays stored RAW/FAIR for derived rows without re-devigging. v1 wires DraftKings and FanDuel.

**Tech Stack:** R (dplyr, tibble, testthat-style asserts inline), Python (curl_cffi-based DK/FD scrapers, pytest), DuckDB (atomic table writes), probit devig (`Tools.R::devig_american` / `devig_american_3way`).

**Branch / Worktree:** This plan executes in worktree `.claude/worktrees/feat+pickem-moneyline-match` on branch `feature/pickem-moneyline-match` (already created from local `main`). On completion: pre-merge review → user approval → merge to `main` → `git worktree remove` + `git branch -d`.

---

## File Structure

| File | Action | Responsibility |
|---|---|---|
| `Answer Keys/MLB Answer Key/odds_screen.R` | Modify | Add `derive_pickem_american()`; extend `.related_market_types` + `scraper_to_canonical` for 3-way; add `line == 0` branch to `expand_bets_to_book_prices` |
| `Answer Keys/tests/test_pickem_match.R` | Create | Unit tests for helper + matcher pick'em branch (2-way and 3-way) |
| `Answer Keys/MLB Answer Key/MLB.R` | Modify (Phase 8) | Add `derived_fair_odds DOUBLE` column to `mlb_bets_book_prices` `CREATE TABLE` |
| `Answer Keys/Tools.R` | Modify | `get_dk_odds()` + `get_fd_odds()`: surface F-period 2-way moneylines as canonical h2h rows with period-tagged market name |
| `mlb_sgp/scraper_draftkings_singles.py` | Modify | `classify_market`: recognize bare "1st N Innings" name as a 2-way period winner (`h2h_period`) |
| `mlb_sgp/tests/test_dk_singles_parser.py` | Modify | New test: bare "1st 3 Innings" classifies as 2-way winner and produces a moneyline row |
| `mlb_sgp/scraper_fanduel_singles.py` | Modify | FD analogue: classify FD's F-period winner market (name confirmed in Task 1) |
| `mlb_sgp/tests/test_fd_singles_parser.py` (or per-file test) | Modify or Create | Mirror test for FD |
| `Answer Keys/MLB Dashboard/book_cell.R` | Modify | `render_book_cell`: when `derived_fair_odds` is non-NULL, display stored values per RAW/FAIR toggle without pair-devig |
| `Answer Keys/MLB Dashboard/mlb_dashboard.R` | Modify | Load `derived_fair_odds` from `mlb_bets_book_prices`; pass through to renderer |
| `Answer Keys/CLAUDE.md` | Modify | Document the pick'em→winner matching rule + `derived_fair_odds` |
| `mlb_sgp/README.md` | Modify | Note DK/FD period-winner capture |

---

## Task 0: Setup verification

**Files:** none (read-only).

- [ ] **Step 0.1: Confirm worktree state**

Run from worktree root:
```bash
git -C /Users/callancapitolo/NFLWork/.claude/worktrees/feat+pickem-moneyline-match status
git -C /Users/callancapitolo/NFLWork/.claude/worktrees/feat+pickem-moneyline-match log --oneline -3
```
Expected: clean working tree; HEAD = the spec commit on `feature/pickem-moneyline-match`, branched from local `main` (parent: `20dfcb3`).

- [ ] **Step 0.2: Verify R unit-test pattern in repo**

Run:
```bash
ls "/Users/callancapitolo/NFLWork/Answer Keys/tests/"
```
Note the existing test style (e.g. `test_devig_pair_matches_tools.R`, `test_books_strip.R`). New tests follow the same `stopifnot()`-based or `testthat`-style pattern already present.

---

## Task 1: Probe FD's F-period winner markets (discovery)

**Files:** `$CLAUDE_JOB_DIR/probe_fd_winners.py` (scratch, not committed).

This is a discovery task — its output decides FD's classify_market change in Task 8. We don't want to guess FD's market names.

- [ ] **Step 1.1: Write the probe script**

Create `$CLAUDE_JOB_DIR/probe_fd_winners.py`:
```python
import sys, json
sys.path.insert(0, "/Users/callancapitolo/NFLWork")
sys.path.insert(0, "/Users/callancapitolo/NFLWork/mlb_sgp")
from mlb_sgp.scraper_fanduel_singles import fetch_event_markets, list_mlb_events
events = list_mlb_events()
ev = events[0]  # any in-progress MLB event
print(f"probing event {ev}")
mkts = fetch_event_markets(ev["event_id"])  # adapt to FD client's actual function
keep = []
for m in mkts:
    nm = (m.get("name") or m.get("market_name") or "").strip()
    low = nm.lower()
    if "1st 3" in low or "1st 5" in low or "1st 7" in low or "first 3" in low \
       or "first 5" in low or "first 7" in low:
        keep.append(nm)
for nm in sorted(set(keep)):
    print(" ", nm)
```

Note: the actual FD client function names differ — adjust to match `mlb_sgp/scraper_fanduel_singles.py`'s entrypoints. Read the file first if needed.

- [ ] **Step 1.2: Run the probe**

```bash
cd /Users/callancapitolo/NFLWork && python3 "$CLAUDE_JOB_DIR/probe_fd_winners.py"
```

- [ ] **Step 1.3: Record findings**

Add a comment block at the top of `mlb_sgp/scraper_fanduel_singles.py` (in the docstring) listing the FD F-period winner market names found. If FD posts no F-period 2-way winner at all, halt Task 8 and report to the user — v1 may need to redefine FD's role (e.g. FD shows "—" on pick'em cards in v1, fast-follow the run-line→DNB conversion later).

- [ ] **Step 1.4: No commit** (discovery only; findings travel as comments in Task 8's commit).

---

## Task 2: `derive_pickem_american()` helper — 2-way path

**Files:**
- Modify: `Answer Keys/MLB Answer Key/odds_screen.R` (append the new helper near the existing devig helpers)
- Create: `Answer Keys/tests/test_pickem_match.R`

- [ ] **Step 2.1: Write the failing test**

Create `Answer Keys/tests/test_pickem_match.R`:
```r
# Tests for derive_pickem_american() — pick'em DNB math.
# Run from repo root: Rscript -e 'source("Answer Keys/tests/test_pickem_match.R")'
suppressPackageStartupMessages({
  source("Answer Keys/Tools.R")          # devig_american, devig_american_3way
  source("Answer Keys/MLB Answer Key/odds_screen.R")
})

# 2-way: a winner market already excludes ties.
# DK 1st-3-innings COL +140 / LA -180. Probit 2-way devig should produce
# matched (p_away + p_home == 1) probs > 0.
r <- derive_pickem_american(home_raw = -180, away_raw = 140, tie_raw = NA)
stopifnot(is.list(r))
stopifnot(all(c("home_raw_dnb","away_raw_dnb",
                "home_fair_dnb","away_fair_dnb") %in% names(r)))
# Raw equals input for 2-way (no transformation).
stopifnot(r$home_raw_dnb == -180)
stopifnot(r$away_raw_dnb == 140)
# Fair: devigged probs sum to 1, so the two American-odds round-trip back to
# probs that sum to 1.
p_home <- if (r$home_fair_dnb < 0) -r$home_fair_dnb / (-r$home_fair_dnb + 100) else 100 / (r$home_fair_dnb + 100)
p_away <- if (r$away_fair_dnb < 0) -r$away_fair_dnb / (-r$away_fair_dnb + 100) else 100 / (r$away_fair_dnb + 100)
stopifnot(abs((p_home + p_away) - 1) < 1e-6)
cat("OK 2-way path\n")
```

- [ ] **Step 2.2: Run test to verify it fails**

```bash
cd /Users/callancapitolo/NFLWork && Rscript "Answer Keys/tests/test_pickem_match.R"
```
Expected: error `could not find function "derive_pickem_american"`.

- [ ] **Step 2.3: Implement 2-way path**

Append to `Answer Keys/MLB Answer Key/odds_screen.R`. **Dependency note:** this
helper calls `devig_american`, `devig_american_3way`, and `prob_to_american`
— all defined in `Tools.R` (already loaded alongside `odds_screen.R` in both
the pipeline and the test). Do NOT redefine `prob_to_american`; it exists at
`Tools.R:937`. Verified return shapes (controller pre-checked):
`devig_american(odd1, odd2)` → 1-row df with `$p1` (↔odd1), `$p2` (↔odd2);
`prob_to_american(p)` → rounded American (uses `%>%`, fine once dplyr is
attached, which `odds_screen.R` does).
```r
#' Derive the draw-no-bet (pick'em) American odds for a book's period winner.
#'
#' Two source shapes:
#'  - 2-way winner (no tie outcome): the market already excludes ties, so the
#'    devigged probabilities ARE the DNB probabilities. Raw = input unchanged.
#'  - 3-way winner: probit-devig home/away/tie, drop the tie, renormalize.
#'
#' Depends on Tools.R (devig_american, devig_american_3way, prob_to_american).
#'
#' @param home_raw American odds, home side (numeric, e.g. -180)
#' @param away_raw American odds, away side
#' @param tie_raw  American odds for the tie outcome, or NA for 2-way
#' @return list with home_raw_dnb, away_raw_dnb, home_fair_dnb, away_fair_dnb
derive_pickem_american <- function(home_raw, away_raw, tie_raw = NA) {
  na_result <- list(home_raw_dnb = NA_real_, away_raw_dnb = NA_real_,
                    home_fair_dnb = NA_real_, away_fair_dnb = NA_real_)
  if (is.na(home_raw) || is.na(away_raw)) return(na_result)

  # Guarded prob -> American: reuse Tools.R::prob_to_american for valid probs.
  to_amer <- function(p) if (is.na(p) || p <= 0 || p >= 1) NA_real_ else prob_to_american(p)

  if (is.na(tie_raw)) {
    # 2-way: raw_dnb = input; fair_dnb = probit 2-way devig.
    # devig_american(odd1, odd2) -> df with $p1 (↔odd1), $p2 (↔odd2).
    devig <- devig_american(away_raw, home_raw)   # p1 = away, p2 = home
    if (is.null(devig) || any(is.na(c(devig$p1, devig$p2)))) {
      return(list(home_raw_dnb = home_raw, away_raw_dnb = away_raw,
                  home_fair_dnb = NA_real_, away_fair_dnb = NA_real_))
    }
    return(list(
      home_raw_dnb  = home_raw,
      away_raw_dnb  = away_raw,
      home_fair_dnb = to_amer(devig$p2),
      away_fair_dnb = to_amer(devig$p1)
    ))
  }
  # 3-way path implemented in Task 3 (replaces this stop()).
  stop("derive_pickem_american 3-way path not implemented yet")
}
```

- [ ] **Step 2.4: Verify `devig_american` return shape**

Run interactively to confirm the field names (`p1` = away, `p2` = home per `compare_alts_to_samples` calling convention at Tools.R:5097):
```bash
cd /Users/callancapitolo/NFLWork && Rscript -e 'source("Answer Keys/Tools.R"); print(devig_american(140, -180))'
```
Expected: `$p1` (away) ≈ 0.40, `$p2` (home) ≈ 0.60, sum ≈ 1.

If field names differ, fix Step 2.3 accordingly.

- [ ] **Step 2.5: Run test to verify it passes**

```bash
cd /Users/callancapitolo/NFLWork && Rscript "Answer Keys/tests/test_pickem_match.R"
```
Expected: `OK 2-way path`.

- [ ] **Step 2.6: Commit**

```bash
git add "Answer Keys/MLB Answer Key/odds_screen.R" "Answer Keys/tests/test_pickem_match.R"
git commit -m "feat(odds-screen): derive_pickem_american() 2-way path

Pure helper for converting a book's 2-way period-winner odds into the
draw-no-bet American price. A 2-way market already excludes ties, so
raw_dnb = input and fair_dnb = probit 2-way devig. 3-way path follows.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: `derive_pickem_american()` — 3-way path

**Files:**
- Modify: `Answer Keys/MLB Answer Key/odds_screen.R`
- Modify: `Answer Keys/tests/test_pickem_match.R`

- [ ] **Step 3.1: Add failing 3-way test**

Append to `Answer Keys/tests/test_pickem_match.R`:
```r
# 3-way: synthetic Home -150 / Tie +400 / Away +180.
# Raw implied: q_h = 150/250 = 0.60, q_t = 100/500 = 0.20, q_a = 100/280 ≈ 0.357
# Sum ≈ 1.157 (vig ~15.7%). Raw DNB (drop tie, renorm on raw):
#   p_h_raw = 0.60 / (0.60 + 0.357) ≈ 0.627
#   p_a_raw = 0.357 / 0.957 ≈ 0.373
# Fair DNB uses probit 3-way devig first; assert just that sums to 1.
r3 <- derive_pickem_american(home_raw = -150, away_raw = 180, tie_raw = 400)
stopifnot(!is.na(r3$home_fair_dnb), !is.na(r3$away_fair_dnb))
p_h_fair <- if (r3$home_fair_dnb < 0) -r3$home_fair_dnb / (-r3$home_fair_dnb + 100) else 100 / (r3$home_fair_dnb + 100)
p_a_fair <- if (r3$away_fair_dnb < 0) -r3$away_fair_dnb / (-r3$away_fair_dnb + 100) else 100 / (r3$away_fair_dnb + 100)
stopifnot(abs((p_h_fair + p_a_fair) - 1) < 1e-6)
# Raw DNB checks: sum to 1 and home favored.
p_h_raw <- if (r3$home_raw_dnb < 0) -r3$home_raw_dnb / (-r3$home_raw_dnb + 100) else 100 / (r3$home_raw_dnb + 100)
p_a_raw <- if (r3$away_raw_dnb < 0) -r3$away_raw_dnb / (-r3$away_raw_dnb + 100) else 100 / (r3$away_raw_dnb + 100)
stopifnot(abs((p_h_raw + p_a_raw) - 1) < 1e-6)
stopifnot(p_h_raw > p_a_raw)   # home was -150, should be favored after drop-tie
cat("OK 3-way path\n")

# Degenerate input: NA tie returns 2-way result; NA home returns all NA.
stopifnot(is.na(derive_pickem_american(NA, 140, NA)$home_fair_dnb))
cat("OK degenerate inputs\n")
```

- [ ] **Step 3.2: Run test, confirm 3-way assertion fails**

```bash
cd /Users/callancapitolo/NFLWork && Rscript "Answer Keys/tests/test_pickem_match.R"
```
Expected: error from the `stop("...3-way path not implemented yet")` in the helper.

- [ ] **Step 3.3: Implement 3-way path**

Replace the `stop(...)` line in `derive_pickem_american()` with (note:
`devig_american_3way(odd_home, odd_away, odd_tie)` returns a 1-row df with
`$p_home`, `$p_away`, `$p_tie` — controller-verified signature):
```r
  # 3-way: devig home/away/tie -> drop the tie -> renormalize home/away.
  devig3 <- devig_american_3way(home_raw, away_raw, tie_raw)
  if (is.null(devig3) || any(is.na(c(devig3$p_home, devig3$p_away)))) return(na_result)
  denom_f <- devig3$p_home + devig3$p_away
  if (!is.finite(denom_f) || denom_f <= 0) return(na_result)
  # Raw DNB: same drop-tie-renormalize on RAW implied probs (no devig).
  implied <- function(american) {
    if (is.na(american)) return(NA_real_)
    if (american < 0) -american / (-american + 100) else 100 / (american + 100)
  }
  q_h <- implied(home_raw); q_a <- implied(away_raw); denom_r <- q_h + q_a
  list(
    home_raw_dnb  = to_amer(q_h / denom_r),
    away_raw_dnb  = to_amer(q_a / denom_r),
    home_fair_dnb = to_amer(devig3$p_home / denom_f),
    away_fair_dnb = to_amer(devig3$p_away / denom_f)
  )
```
(`to_amer` and `na_result` are in scope — both defined at the top of
`derive_pickem_american` in Task 2.3.)

- [ ] **Step 3.4: Verify `devig_american_3way` return shape**

```bash
cd /Users/callancapitolo/NFLWork && Rscript -e 'suppressWarnings(suppressMessages(source("Answer Keys/MLB Answer Key/odds_screen.R"))); source("Answer Keys/Tools.R"); print(devig_american_3way(-150, 180, 400))'
```
Expected: 1-row data frame with columns `p_home`, `p_away`, `p_tie` summing
to 1 (controller already confirmed this shape; sourcing odds_screen.R first
attaches dplyr so `prob_to_american`'s `%>%` resolves).

- [ ] **Step 3.5: Run test to verify all assertions pass**

```bash
cd /Users/callancapitolo/NFLWork && Rscript "Answer Keys/tests/test_pickem_match.R"
```
Expected: `OK 2-way path`, `OK 3-way path`, `OK degenerate inputs`.

- [ ] **Step 3.6: Commit**

```bash
git add "Answer Keys/MLB Answer Key/odds_screen.R" "Answer Keys/tests/test_pickem_match.R"
git commit -m "feat(odds-screen): derive_pickem_american() 3-way path

Probit 3-way devig -> drop the tie -> renormalize home/away to draw-no-bet.
Raw DNB is the same drop-tie-renormalize on raw implied probs (no devig).
Handles NA inputs and degenerate denominators by returning NA pairs.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: Canonical 3-way market shape

**Files:**
- Modify: `Answer Keys/MLB Answer Key/odds_screen.R` (extend `.derive_period`/`.derive_market_type`, `.related_market_types`, and `scraper_to_canonical`'s h2h branch to emit a Tie row when present)
- Modify: `Answer Keys/tests/test_pickem_match.R`

- [ ] **Step 4.1: Add failing test for 3-way canonicalization**

Append to `Answer Keys/tests/test_pickem_match.R`:
```r
# A wide-scraper row with odds_home/odds_away/odds_tie and line == NA should
# canonicalize to three h2h_3way rows (home, away, Tie) with the right period.
library(tibble)
raw <- tibble(
  market = "h2h_3way_1st_5_innings",
  home_team = "LA Dodgers", away_team = "COL Rockies",
  line = NA_real_,
  odds_home = -150, odds_away = 180, odds_tie = 400,
  fetch_time = Sys.time()
)
lookup <- tibble(id = "g1", home_team = "LA Dodgers", away_team = "COL Rockies")
canon <- scraper_to_canonical(raw, lookup, book_name = "test")
stopifnot(nrow(canon) == 3)
stopifnot(all(canon$market == "h2h_3way"))
stopifnot(all(canon$period == "F5"))
stopifnot(setequal(canon$side, c("LA Dodgers", "COL Rockies", "Tie")))
cat("OK 3-way canonicalization\n")
```

- [ ] **Step 4.2: Run test, confirm it fails**

```bash
cd /Users/callancapitolo/NFLWork && Rscript "Answer Keys/tests/test_pickem_match.R"
```
Expected: `nrow(canon) == 3` assertion fails (current canonical only emits 2 h2h rows, ignores tie).

- [ ] **Step 4.3: Extend `.derive_market_type` to recognize `h2h_3way_*`**

In `odds_screen.R::.derive_market_type`, update the regex strip so `h2h_3way_1st_5_innings` → `h2h_3way` (and similar for F3/F7/FG suffixes):
```r
.derive_market_type <- function(market_name) {
  market_name %>%
    str_replace("_1st_[357]_innings$", "") %>%
    str_replace("_(fg|f3|f5|f7|h2)$", "")
}
```

Verify: `.derive_market_type("h2h_3way_1st_5_innings")` returns `"h2h_3way"`. Same function already handles the `_1st_[357]_innings$` strip — no change needed if the input naming follows that pattern. Confirm with a quick `cat(.derive_market_type("h2h_3way_1st_5_innings"))`.

- [ ] **Step 4.4: Extend `scraper_to_canonical()` h2h branch to emit a Tie row**

In the per-row loop, after the existing two `is.na(row$line)` branches that emit home_team / away_team h2h rows, add:
```r
    # 3-way moneyline (h2h_3way_*): emit a "Tie" row alongside home/away.
    if (!is.na(row$odds_home) && is.na(row$line) &&
        "odds_tie" %in% names(row) && !is.na(row$odds_tie)) {
      rows[[length(rows) + 1]] <- tibble(
        game_id = gid, market_name = mkt,
        bet_on = "Tie", line = NA_real_,
        american_odds = as.integer(row$odds_tie), fetch_time = ft
      )
    }
```

Place this immediately after the existing away-team h2h emission. The two existing branches (home and away) already fire on `is.na(row$line)`; the new branch ALSO fires only when an `odds_tie` column is present on this scraper row — older 2-way scraper frames are unaffected (column absent).

- [ ] **Step 4.5: Extend `.related_market_types` so pick'em can find both 2-way and 3-way winners**

```r
.related_market_types <- function(mt) {
  if (mt == "spreads"           || mt == "alternate_spreads")
    return(c("spreads", "alternate_spreads"))
  if (mt == "totals"            || mt == "alternate_totals")
    return(c("totals", "alternate_totals"))
  if (mt == "h2h" || mt == "h2h_3way")
    return(c("h2h", "h2h_3way"))
  mt
}
```

(The pick'em branch in Task 6 will *additionally* expand a `spreads`-bet at line 0 to include `h2h` + `h2h_3way`. Adding h2h↔h2h_3way here means a future h2h bet would also see both — fine, since the matcher's side filter (team name vs "Tie") already excludes the tie row from matching to a real bet side.)

- [ ] **Step 4.6: Run test to verify it passes**

```bash
cd /Users/callancapitolo/NFLWork && Rscript "Answer Keys/tests/test_pickem_match.R"
```
Expected: all prior tests still pass plus `OK 3-way canonicalization`.

- [ ] **Step 4.7: Commit**

```bash
git add "Answer Keys/MLB Answer Key/odds_screen.R" "Answer Keys/tests/test_pickem_match.R"
git commit -m "feat(odds-screen): canonical 3-way market shape (h2h_3way)

scraper_to_canonical emits a third 'Tie' row when a scraper wide-frame
provides odds_tie alongside odds_home/odds_away (and line == NA). The new
market_type 'h2h_3way' joins 'h2h' in .related_market_types so the matcher
can reach either. No effect on existing 2-way scrapers (odds_tie absent).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: Schema — add `derived_fair_odds` to `mlb_bets_book_prices`

**Files:**
- Modify: `Answer Keys/MLB Answer Key/MLB.R` (Phase 8, `CREATE TABLE` for `mlb_bets_book_prices`)

- [ ] **Step 5.1: Locate the CREATE TABLE block**

```bash
grep -n "CREATE TABLE mlb_bets_book_prices" "/Users/callancapitolo/NFLWork/Answer Keys/MLB Answer Key/MLB.R"
```
Expected match around line 977 (existing CREATE TABLE with explicit schema, immediately after `DROP TABLE IF EXISTS mlb_bets_book_prices`).

- [ ] **Step 5.2: Add the new column to the schema**

Edit the CREATE TABLE to add `derived_fair_odds DOUBLE` after `american_odds`:
```r
dbExecute(con_bets, "
  CREATE TABLE mlb_bets_book_prices (
    bet_row_id      VARCHAR,
    game_id         VARCHAR,
    market          VARCHAR,
    period          VARCHAR,
    side            VARCHAR,
    bookmaker       VARCHAR,
    line            DOUBLE,
    line_quoted     DOUBLE,
    is_exact_line   BOOLEAN,
    american_odds   INTEGER,
    derived_fair_odds DOUBLE,
    fetch_time      TIMESTAMPTZ,
    game_start_time TIMESTAMPTZ
  )
")
```

The column is added before `fetch_time` so the order matches the order columns appear in the matcher's output tibble (Task 6).

- [ ] **Step 5.3: Confirm `dbAppendTable` still aligns**

`dbAppendTable` matches by column name, not position, so no further change here. Just confirm the matcher output (after Task 6) emits a `derived_fair_odds` column — Task 6's step 6.4 covers that.

- [ ] **Step 5.4: Commit**

```bash
git add "Answer Keys/MLB Answer Key/MLB.R"
git commit -m "schema(mlb_bets_book_prices): add derived_fair_odds DOUBLE column

NULL for all rows except derived pick'em cells, which precompute the
draw-no-bet FAIR price in the matcher so book_cell.R can render it
directly instead of re-devigging an artificial pair.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 6: Matcher `line == 0` branch in `expand_bets_to_book_prices`

**Files:**
- Modify: `Answer Keys/MLB Answer Key/odds_screen.R`
- Modify: `Answer Keys/tests/test_pickem_match.R`

- [ ] **Step 6.1: Add failing matcher test for a line-0 spread bet vs a 2-way winner**

Append to `Answer Keys/tests/test_pickem_match.R`:
```r
# Bet: PHI 0 spread, period F3. Book: DK with a 2-way h2h_1st_3_innings winner.
# Matcher should emit a pick row and an opposite row using derive_pickem_american.
bets <- tibble(
  bet_row_id = "b1", id = "g1", home_team = "SD Padres", away_team = "PHI Phillies",
  market = "spreads_1st_3_innings", market_type = "spreads", period = "F3",
  line = 0, bet_on = "PHI Phillies",
  opposite_side = "SD Padres",
  pt_start_time = as.POSIXct(NA, tz = "UTC")
)
dk_canonical <- normalize_book_odds_frame(tibble(
  game_id = "g1", market_name = "h2h_1st_3_innings",
  bet_on = c("SD Padres", "PHI Phillies"),
  line = c(NA_real_, NA_real_),
  american_odds = c(-180L, 140L),
  fetch_time = Sys.time()
))
out <- expand_bets_to_book_prices(bets, list(dk = dk_canonical))
# Two rows: pick + opposite, both at line 0, both is_exact_line=TRUE, both
# carry derived_fair_odds.
stopifnot(nrow(out) == 2)
stopifnot(all(out$line == 0))
stopifnot(all(out$is_exact_line))
stopifnot(all(!is.na(out$derived_fair_odds)))
# Pick (PHI/away) raw == 140; opposite (SD/home) raw == -180.
pick <- out[out$side == "pick", ]; opp <- out[out$side == "opposite", ]
stopifnot(pick$american_odds == 140)
stopifnot(opp$american_odds == -180)
cat("OK matcher 2-way pick'em\n")
```

- [ ] **Step 6.2: Run test, confirm failure**

```bash
cd /Users/callancapitolo/NFLWork && Rscript "Answer Keys/tests/test_pickem_match.R"
```
Expected: matcher returns 0 rows (no line-0 branch) OR a fallback row without `derived_fair_odds`. Either way, asserts fail.

- [ ] **Step 6.3: Add the pick'em branch to `expand_bets_to_book_prices`**

Inside the `for (book_name in names(book_odds_by_book))` loop, BEFORE the existing `for (slot in names(side_labels))` block, branch on the pick'em case:
```r
      # Pick'em (line-0 spread): derive draw-no-bet from the book's period
      # winner instead of falling back to its nearest run line. h2h (2-way)
      # is used directly; h2h_3way is devigged + tie dropped + renormalized.
      is_pickem <- bet$market_type %in% c("spreads", "alternate_spreads") &&
                   !is.na(bet$line) && bet$line == 0
      if (is_pickem) {
        winners <- book_frame %>%
          filter(game_id == bet$game_id,
                 period  == bet$period,
                 market  %in% c("h2h", "h2h_3way", "spreads", "alternate_spreads"))
        # For spread rows in the pick'em candidate set, only true 0-handicaps
        # count — anything else is a real run line and was already excluded
        # by the bet's line == 0 trigger.
        winners <- winners %>% filter(
          market %in% c("h2h", "h2h_3way") |
          (market %in% c("spreads", "alternate_spreads") &
             !is.na(line) & abs(line) < 1e-9)
        )
        if (nrow(winners) == 0) next   # cell is "—"

        # Decide source shape. Priority: 2-way h2h > 3-way > true 0-handicap.
        h2h_rows  <- winners %>% filter(market == "h2h")
        h3w_rows  <- winners %>% filter(market == "h2h_3way")
        sp0_rows  <- winners %>% filter(market %in% c("spreads", "alternate_spreads"))
        use_3way  <- nrow(h2h_rows) < 2 && nrow(h3w_rows) >= 3
        use_sp0   <- nrow(h2h_rows) < 2 && nrow(h3w_rows) < 3 && nrow(sp0_rows) >= 2

        if (use_3way) {
          home_row <- h3w_rows %>% filter(side == home_team_value) %>% slice_head(n = 1)
          away_row <- h3w_rows %>% filter(side == away_team_value) %>% slice_head(n = 1)
          tie_row  <- h3w_rows %>% filter(side == "Tie")            %>% slice_head(n = 1)
          if (nrow(home_row) == 0 || nrow(away_row) == 0 || nrow(tie_row) == 0) next
          dnb <- derive_pickem_american(home_row$american_odds,
                                        away_row$american_odds,
                                        tie_row$american_odds)
          ft  <- home_row$fetch_time
        } else if (use_sp0) {
          # True 0-handicap: treat the spread rows as a 2-way pair.
          home_row <- sp0_rows %>% filter(side == home_team_value) %>% slice_head(n = 1)
          away_row <- sp0_rows %>% filter(side == away_team_value) %>% slice_head(n = 1)
          if (nrow(home_row) == 0 || nrow(away_row) == 0) next
          dnb <- derive_pickem_american(home_row$american_odds,
                                        away_row$american_odds,
                                        NA)
          ft  <- home_row$fetch_time
        } else {
          home_row <- h2h_rows %>% filter(side == home_team_value) %>% slice_head(n = 1)
          away_row <- h2h_rows %>% filter(side == away_team_value) %>% slice_head(n = 1)
          if (nrow(home_row) == 0 || nrow(away_row) == 0) next
          dnb <- derive_pickem_american(home_row$american_odds,
                                        away_row$american_odds,
                                        NA)
          ft  <- home_row$fetch_time
        }
        if (is.na(dnb$home_raw_dnb) || is.na(dnb$away_raw_dnb)) next

        # Pick is the team the bet is on; opposite is the other team.
        bet_is_home <- identical(bet$bet_on, home_team_value)
        pick_raw  <- if (bet_is_home) dnb$home_raw_dnb  else dnb$away_raw_dnb
        opp_raw   <- if (bet_is_home) dnb$away_raw_dnb  else dnb$home_raw_dnb
        pick_fair <- if (bet_is_home) dnb$home_fair_dnb else dnb$away_fair_dnb
        opp_fair  <- if (bet_is_home) dnb$away_fair_dnb else dnb$home_fair_dnb

        for (slot in c("pick", "opposite")) {
          american <- if (slot == "pick") pick_raw  else opp_raw
          fair_odd <- if (slot == "pick") pick_fair else opp_fair
          if (is.na(american)) next
          k <- k + 1L
          out_rows[[k]] <- tibble(
            bet_row_id        = bet$bet_row_id,
            game_id           = bet$game_id,
            market            = bet$market_type,
            period            = bet$period,
            side              = slot,
            bookmaker         = book_name,
            line              = bet$line,
            line_quoted       = 0,
            is_exact_line     = TRUE,
            american_odds     = as.integer(round(american)),
            derived_fair_odds = as.numeric(fair_odd),
            fetch_time        = ft,
            game_start_time   = bet_game_start
          )
        }
        next   # do not fall through to the run-line matching path
      }
```

Also: in the **default branch** (the existing `for (slot in names(side_labels))` block), every emitted `out_rows[[k]] <- tibble(...)` must add `derived_fair_odds = NA_real_,` so the column is uniformly present. Place it adjacent to `american_odds`.

- [ ] **Step 6.4: Update the zero-rows fallback in `expand_bets_to_book_prices`**

Both empty-return `tibble(...)` constructors at the top and bottom of the function gain `derived_fair_odds = numeric()`. Match the column ordering used in the CREATE TABLE.

- [ ] **Step 6.5: Run test to verify the 2-way matcher case passes**

```bash
cd /Users/callancapitolo/NFLWork && Rscript "Answer Keys/tests/test_pickem_match.R"
```
Expected: all prior tests pass plus `OK matcher 2-way pick'em`.

- [ ] **Step 6.6: Add 3-way matcher test**

Append to `Answer Keys/tests/test_pickem_match.R`:
```r
# Same bet, but the book only has a 3-way h2h_3way winner.
wz_canonical <- normalize_book_odds_frame(tibble(
  game_id = "g1", market_name = "h2h_3way_1st_3_innings",
  bet_on = c("SD Padres", "PHI Phillies", "Tie"),
  line = c(NA_real_, NA_real_, NA_real_),
  american_odds = c(-150L, 180L, 400L),
  fetch_time = Sys.time()
))
out2 <- expand_bets_to_book_prices(bets, list(wz = wz_canonical))
stopifnot(nrow(out2) == 2)
stopifnot(all(out2$is_exact_line))
stopifnot(all(!is.na(out2$derived_fair_odds)))
cat("OK matcher 3-way pick'em\n")
```

```bash
cd /Users/callancapitolo/NFLWork && Rscript "Answer Keys/tests/test_pickem_match.R"
```
Expected: `OK matcher 3-way pick'em`.

- [ ] **Step 6.7: Add regression test that a NON-zero spread bet is unaffected**

Append to `Answer Keys/tests/test_pickem_match.R`:
```r
# A non-zero spread bet must NOT enter the pick'em branch.
bets_alt <- tibble(
  bet_row_id = "b2", id = "g1", home_team = "SD Padres", away_team = "PHI Phillies",
  market = "alternate_spreads_fg", market_type = "alternate_spreads", period = "FG",
  line = -2.5, bet_on = "PHI Phillies",
  opposite_side = "SD Padres",
  pt_start_time = as.POSIXct(NA, tz = "UTC")
)
fake_book <- normalize_book_odds_frame(tibble(
  game_id = "g1", market_name = "alternate_spreads_fg",
  bet_on = c("PHI Phillies"),
  line = -2.5, american_odds = 240L, fetch_time = Sys.time()
))
out3 <- expand_bets_to_book_prices(bets_alt, list(test = fake_book))
stopifnot(nrow(out3) >= 1)
stopifnot(all(is.na(out3$derived_fair_odds)))
cat("OK non-zero spread unaffected\n")
```

```bash
cd /Users/callancapitolo/NFLWork && Rscript "Answer Keys/tests/test_pickem_match.R"
```
Expected: `OK non-zero spread unaffected`.

- [ ] **Step 6.8: Commit**

```bash
git add "Answer Keys/MLB Answer Key/odds_screen.R" "Answer Keys/tests/test_pickem_match.R"
git commit -m "feat(odds-screen): line-0 spread matches book period winner (DNB)

expand_bets_to_book_prices(): when the bet is a spread/alt-spread at line
exactly 0, route per-book to its h2h or h2h_3way winner rows and emit cells
with derive_pickem_american (raw + fair DNB). Non-zero spreads are
untouched. Books with no period winner produce no row (cell -> '—').

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 7: DK scraper — classify "1st N Innings" as 2-way period winner

**Files:**
- Modify: `mlb_sgp/scraper_draftkings_singles.py::classify_market`
- Modify: `mlb_sgp/tests/test_dk_singles_parser.py`

- [ ] **Step 7.1: Write failing parser test**

Append to `mlb_sgp/tests/test_dk_singles_parser.py`:
```python
def test_bare_period_name_classifies_as_two_way_winner():
    """DK names the period winner exactly '1st 3 Innings' (and 1st 5/7).
    classify_market must return (period, 'main') with ML semantics so the
    selections (no line, team names) become away_ml/home_ml."""
    assert classify_market("1st 3 Innings") == ("F3", "main")
    assert classify_market("1st 5 Innings") == ("F5", "main")
    assert classify_market("1st 7 Innings") == ("F7", "main")


def test_period_winner_selections_yield_moneyline_row():
    event = Event("e1", "Yankees", "Red Sox", "2026-05-12T22:00:00Z")
    selections = [
        Selection("s1", "m_pw", "Yankees", None, -180),
        Selection("s2", "m_pw", "Red Sox", None, 140),
    ]
    market_meta = {"m_pw": ("F3", "main")}
    rows = parse_selections_to_wide_rows(event, selections, market_meta,
                                          fetch_time=datetime(2026, 5, 12, 14, 0))
    assert len(rows) == 1
    r = rows[0]
    assert r["period"] == "F3"
    assert r["market"] == "main"
    assert r["home_ml"] == -180
    assert r["away_ml"] == 140
```

Run:
```bash
cd /Users/callancapitolo/NFLWork && PYTHONPATH="$PWD:$PWD/mlb_sgp" python3 -m pytest mlb_sgp/tests/test_dk_singles_parser.py -q
```
Expected: `test_bare_period_name_classifies_as_two_way_winner` fails — `classify_market("1st 3 Innings")` returns `None`.

- [ ] **Step 7.2: Locate the bare-period-name branch in `classify_market`**

`classify_market` already does period detection (`"1st 3 innings" in n → period = "F3"`) and market-type detection. The fix: a bare period name (no "run line"/"moneyline"/"total"/"alternate") IS a winner moneyline.

Add a new clause just before the final `return None`:
```python
    # Bare period name like "1st 3 Innings" / "1st 5 Innings" / "1st 7 Innings"
    # is DK's 2-way period winner (no tie outcome on first-N-innings). The
    # selections have no line — they'll flow through the moneyline detection
    # branch in parse_selections_to_wide_rows as home_ml/away_ml.
    if n in ("1st 3 innings", "1st 5 innings", "1st 7 innings",
             "first 3 innings", "first 5 innings", "first 7 innings"):
        return (period, "main")
```

- [ ] **Step 7.3: Run the parser test suite**

```bash
cd /Users/callancapitolo/NFLWork && PYTHONPATH="$PWD:$PWD/mlb_sgp" python3 -m pytest mlb_sgp/tests/test_dk_singles_parser.py -q
```
Expected: all tests pass.

- [ ] **Step 7.4: Live smoke-check: a real DK event yields F3 ML rows**

```bash
cd /Users/callancapitolo/NFLWork && python3 -c "
import sys; sys.path.insert(0, '.'); sys.path.insert(0, 'mlb_sgp')
from datetime import datetime, timezone
from mlb_sgp.dk_client import DraftKingsClient
from mlb_sgp.scraper_draftkings_sgp import DK_SGP_PARLAYS_URL
from mlb_sgp.scraper_draftkings_singles import classify_market, parse_selections_to_wide_rows
c = DraftKingsClient()
ev = c.list_events()[0]
mkts = c.fetch_event_markets(ev.event_id)
sels = c.fetch_event_selections(ev.event_id)
meta = {m.market_id: classify_market(m.name) for m in mkts if classify_market(m.name)}
r = c.session.get(f'{DK_SGP_PARLAYS_URL}/{ev.event_id}', timeout=60)
for m in ((r.json() or {}).get('data') or {}).get('markets') or []:
    mid = str(m.get('id',''))
    if mid and mid not in meta:
        cl = classify_market(m.get('name','') or '')
        if cl: meta[mid] = cl
rows = parse_selections_to_wide_rows(ev, sels, meta, datetime.now(timezone.utc))
f3_ml = [x for x in rows if x['period']=='F3' and x['market']=='main' and (x['home_ml'] is not None or x['away_ml'] is not None)]
print('F3 ML rows for', ev.away_team, '@', ev.home_team, ':', len(f3_ml))
for r in f3_ml:
    print(' ', r['away_team'], '@', r['away_ml'], '/', r['home_team'], '@', r['home_ml'])
"
```
Expected: 1 F3 ML row with both team prices populated (matching the DK board).

- [ ] **Step 7.5: Commit**

```bash
git add mlb_sgp/scraper_draftkings_singles.py mlb_sgp/tests/test_dk_singles_parser.py
git commit -m "feat(dk-scraper): capture bare-period-name 2-way period winners

DK names the first-N-innings winner exactly '1st 3 Innings' / '1st 5
Innings' / '1st 7 Innings' — no 'run line' / 'moneyline' keyword, so
classify_market returned None and the selections were dropped. Recognize
the bare period name as (period, 'main'); the existing moneyline branch in
parse_selections_to_wide_rows then captures the team prices as
away_ml/home_ml. This is the DNB-equivalent the dashboard needs to render
pick'em (line-0 spread) cards.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 8: `get_dk_odds()` — flow F-period ML to canonical `h2h_<period>`

**Files:**
- Modify: `Answer Keys/Tools.R` (`get_dk_odds`)

The DK scraper now writes F-period ML rows (`market = "main"`, `period = "F3"/"F5"/"F7"`, `away_ml`/`home_ml` populated, no spread/total). `get_dk_odds` reads those rows into the wide canonical frame; `scraper_to_canonical` then produces h2h rows. The frame's `market` column must carry a name that `.derive_period` will resolve to the right period (e.g. `h2h_1st_3_innings`).

- [ ] **Step 8.1: Read the current `get_dk_odds()` to see how it labels markets**

```bash
grep -n "get_dk_odds <- function\|h2h\|market =" "/Users/callancapitolo/NFLWork/Answer Keys/Tools.R" | head -30
```

- [ ] **Step 8.2: Verify the current behavior for an F-period ML row**

Write a one-off probe:
```r
# Rscript -e ...
source("Answer Keys/Tools.R")
dk <- get_dk_odds()
dk_f3_ml <- subset(dk, period == "F3" & is.na(line) &
                   (!is.na(odds_home) | !is.na(odds_away)))
print(head(dk_f3_ml[, c("home_team","away_team","market","period",
                        "line","odds_home","odds_away")]))
```
- If `market` is already a value like `h2h_1st_3_innings`, no change is needed in Step 8.3.
- If `market` is `main` or `h2h` without a period suffix, you must remap it before `scraper_to_canonical` consumes it.

- [ ] **Step 8.3: If a remap is needed, add it inside `get_dk_odds()`**

```r
  # F-period moneylines (no line, no spread) need a period-tagged market label
  # so scraper_to_canonical's canonicalization (.derive_period) puts them in
  # the right period bucket.
  result <- result %>%
    mutate(market = case_when(
      period == "F3" & market %in% c("h2h", "main") & is.na(line) ~ "h2h_1st_3_innings",
      period == "F5" & market %in% c("h2h", "main") & is.na(line) ~ "h2h_1st_5_innings",
      period == "F7" & market %in% c("h2h", "main") & is.na(line) ~ "h2h_1st_7_innings",
      TRUE ~ market
    ))
```

Place this just before the function's `return(result)`. Apply the same remap in `get_fd_odds()` in Task 10 — the exact location varies per helper; align with their existing post-processing patterns.

- [ ] **Step 8.4: Confirm `.derive_period("h2h_1st_3_innings")` returns "F3"**

```bash
Rscript -e 'source("Answer Keys/MLB Answer Key/odds_screen.R"); cat(.derive_period("h2h_1st_3_innings"), "\n")'
```
Expected: `F3`.

- [ ] **Step 8.5: Run the existing R test suite to confirm nothing else broke**

```bash
cd /Users/callancapitolo/NFLWork && Rscript "Answer Keys/tests/test_pickem_match.R"
Rscript "Answer Keys/tests/test_devig_pair_matches_tools.R" 2>/dev/null || true
Rscript "Answer Keys/tests/test_books_strip.R" 2>/dev/null || true
```
Expected: all existing tests still pass.

- [ ] **Step 8.6: Commit**

```bash
git add "Answer Keys/Tools.R"
git commit -m "feat(tools): get_dk_odds flows F-period MLs as h2h_1st_N_innings

Period-tag the market label so scraper_to_canonical's .derive_period puts
the F3/F5/F7 winner rows in the right period bucket. Only fires on rows
with no line and a populated ML side; existing labels untouched.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 9: FD scraper — classify FD's period winner market

**Files:**
- Modify: `mlb_sgp/scraper_fanduel_singles.py::classify_market`
- Modify: FD parser tests (path determined by reading the file)

This task depends on Task 1's discovery output. Steps below assume FD names its period winner with a recognizable pattern; substitute the actual name(s) from Task 1.

- [ ] **Step 9.1: Open `scraper_fanduel_singles.py` and find `classify_market` (or equivalent)**

```bash
grep -n "def classify_market\|def classify\|def _classify\|market.*1st 3\|first 3 innings" "/Users/callancapitolo/NFLWork/mlb_sgp/scraper_fanduel_singles.py"
```

- [ ] **Step 9.2: Add the FD-specific period-winner recognition**

Pick the branch matching Task 1's findings.

**Branch B1 — FD posts a bare period name** ("1st 3 Innings", "1st 5 Innings", "1st 7 Innings"). Add to FD's classifier near where it does period detection:
```python
if n in ("1st 3 innings", "1st 5 innings", "1st 7 innings",
         "first 3 innings", "first 5 innings", "first 7 innings"):
    return (period, "main")   # 2-way winner (no tie outcome)
```

**Branch B2 — FD names it "Moneyline - 1st N Innings"** (FD's main FG ML is `Moneyline`, so this is plausible). Add a separate recognizer:
```python
if "moneyline" in n and ("1st 3 innings" in n or "first 3 innings" in n):
    return ("F3", "main")
if "moneyline" in n and ("1st 5 innings" in n or "first 5 innings" in n):
    return ("F5", "main")
if "moneyline" in n and ("1st 7 innings" in n or "first 7 innings" in n):
    return ("F7", "main")
```
(Place these BEFORE the existing FG-only `moneyline` clause so the period-suffixed match wins.)

**Branch B3 — FD posts no F-period 2-way winner.** Add this comment near FD's classifier and stop Task 9:
```python
# FD does not currently post a 2-way period winner for 1st 3/5/7 innings;
# pick'em cards will show "—" for FD in v1. Revisit if FD adds the market
# or if a 3-way derivation path is added.
```
Then skip Step 9.4 and proceed to Task 10's smoke check noting FD as absent.

Whichever branch you pick, **document the chosen branch (and the FD market names you saw in Task 1) in the FD scraper module docstring** so the next reader knows why the classifier looks the way it does.

- [ ] **Step 9.3: Write FD parser test**

Mirror the DK pattern: assert FD's winner market name classifies to `(period, ...)` with ML semantics and that selections (no line, team names) flow as moneyline. The test file is whichever FD parser tests already exist; if none exist, create `mlb_sgp/tests/test_fd_singles_parser.py` modeled after `test_dk_singles_parser.py`.

- [ ] **Step 9.4: Run FD parser tests**

```bash
cd /Users/callancapitolo/NFLWork && PYTHONPATH="$PWD:$PWD/mlb_sgp" python3 -m pytest mlb_sgp/tests/test_fd_singles_parser.py -q
```
Expected: tests pass.

- [ ] **Step 9.5: Live smoke-check**

Adapt the Task 7.4 script for FD: list events, classify, parse, count F3 ML rows. Confirm at least one MLB event yields an F3 ML row with both team prices.

- [ ] **Step 9.6: Commit**

```bash
git add mlb_sgp/scraper_fanduel_singles.py mlb_sgp/tests/test_fd_singles_parser.py
git commit -m "feat(fd-scraper): capture FD F-period 2-way winners as moneylines

Recognize FD's first 3/5/7 innings winner market in classify_market so
the selections flow as moneylines into get_fd_odds() and produce
canonical h2h_1st_N_innings rows for the dashboard pick'em path. Exact
branch chosen documented in the FD scraper module docstring.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 10: `get_fd_odds()` — flow F-period ML to canonical `h2h_<period>`

**Files:**
- Modify: `Answer Keys/Tools.R` (`get_fd_odds`)

Mirror Task 8 for FD.

- [ ] **Step 10.1: Apply the same period-tagged label remap inside `get_fd_odds()`**

Use the exact `mutate(market = case_when(...))` block from Task 8.3, placed at the end of `get_fd_odds()`.

- [ ] **Step 10.2: Probe to confirm F-period ML rows flow as h2h_1st_N_innings**

```r
Rscript -e 'source("Answer Keys/Tools.R"); fd <- get_fd_odds(); fd_pw <- subset(fd, grepl("^h2h_1st_[357]_innings$", market)); print(head(fd_pw[, c("home_team","away_team","market","odds_home","odds_away")]))'
```
Expected: at least one row.

- [ ] **Step 10.3: Run all R tests**

```bash
cd /Users/callancapitolo/NFLWork && Rscript "Answer Keys/tests/test_pickem_match.R"
```

- [ ] **Step 10.4: Commit**

```bash
git add "Answer Keys/Tools.R"
git commit -m "feat(tools): get_fd_odds flows F-period MLs as h2h_1st_N_innings

Mirror get_dk_odds. Period-tag the market label so scraper_to_canonical
routes FD's F3/F5/F7 winner rows to the right period bucket.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 11: `book_cell.R` — render derived rows without re-devigging

**Files:**
- Modify: `Answer Keys/MLB Dashboard/book_cell.R`

- [ ] **Step 11.1: Read `render_book_cell` to find where RAW/FAIR is computed**

```bash
grep -n "render_book_cell\|devig_american_pair\|fair\|raw" "/Users/callancapitolo/NFLWork/Answer Keys/MLB Dashboard/book_cell.R" | head
```

- [ ] **Step 11.2: Add the `derived_fair_odds` arg and the bypass branch**

`render_book_cell()` already accepts the row's American odds. Extend its signature with `derived_fair_odds = NA_real_` and short-circuit the devig path when non-NULL:
```r
render_book_cell <- function(american_odds, opposite_american_odds, ...,
                              derived_fair_odds = NA_real_,
                              ...) {
  # ...existing setup...

  if (!is.na(derived_fair_odds)) {
    # Derived pick'em cell. RAW = stored american_odds (raw DNB), FAIR = the
    # precomputed devigged DNB. No pair-devig — the matcher did the math.
    raw_disp  <- format_american(american_odds)
    fair_disp <- format_american(derived_fair_odds)
    return(htmltools::tagList(
      span(class = "cell-raw",  raw_disp),
      span(class = "cell-fair", fair_disp)
    ))
  }

  # ...existing on-the-fly pair-devig path unchanged...
}
```

Adapt the exact HTML/tag structure to the existing patterns in `book_cell.R`. Key invariant: when `derived_fair_odds` is provided, **do not call `.devig_american_pair`** for this cell.

- [ ] **Step 11.3: Update all `render_book_cell()` call sites in `mlb_dashboard.R`**

```bash
grep -n "render_book_cell\|render_book_pill" "/Users/callancapitolo/NFLWork/Answer Keys/MLB Dashboard/mlb_dashboard.R"
```
Each call site must pass `derived_fair_odds = <row>$derived_fair_odds` (or named field if pivoted). Use `NA_real_` for any call that doesn't have a derived value.

- [ ] **Step 11.4: Manual render check (uses copied live DBs per the verify-UI policy)**

```bash
# Copy the live DBs into a scratch dir.
mkdir -p "$CLAUDE_JOB_DIR/render_test"
cp "/Users/callancapitolo/NFLWork/Answer Keys/mlb_mm.duckdb"        "$CLAUDE_JOB_DIR/render_test/"
cp "/Users/callancapitolo/NFLWork/Answer Keys/MLB Dashboard/mlb_dashboard.duckdb" "$CLAUDE_JOB_DIR/render_test/" 2>/dev/null || true
# Render the report against the copy. (Adapt to the project's render entrypoint
# — see Answer Keys/MLB Dashboard/README.md or the report rmarkdown.)
```

For the first-3-innings PHI/SD pick'em card in the rendered report:
- DK should show its 2-way DNB price (no line tag).
- FD should show its DNB price (no line tag) if FD's winner capture landed.
- Other books should show "—" (no run-line fallback).
- Non-pick'em cards (e.g. the ARI -2.5 alt card from the FG fix branch) must look identical to before.

- [ ] **Step 11.5: Commit**

```bash
git add "Answer Keys/MLB Dashboard/book_cell.R" "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "feat(dashboard): render derived pick'em cells without pair-devig

book_cell.R::render_book_cell short-circuits the devig path when the row
carries a non-NULL derived_fair_odds: RAW shows stored american_odds (raw
DNB), FAIR shows the precomputed devigged DNB. All other cells render
exactly as before. mlb_dashboard.R passes derived_fair_odds through.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 12: End-to-end pipeline verification

**Files:** none (verification only).

- [ ] **Step 12.1: Trigger a full MLB pipeline run** (or the scraper-only path, whichever the dashboard's `/refresh` uses today)

```bash
# Adapt to the project's invocation — see run.py mlb or the dashboard refresh
cd /Users/callancapitolo/NFLWork && python3 run.py mlb
```

- [ ] **Step 12.2: Query the resulting `mlb_bets_book_prices` for a pick'em bet**

```bash
python3 -c "
import duckdb
con = duckdb.connect('/Users/callancapitolo/NFLWork/Answer Keys/mlb_mm.duckdb', read_only=True)
# Find a pick'em bet (line == 0 spread) and inspect its book-price rows.
df = con.execute('''
  SELECT b.bet_row_id, b.bookmaker, b.side, b.market, b.period,
         b.line, b.line_quoted, b.is_exact_line, b.american_odds,
         b.derived_fair_odds
  FROM mlb_bets_book_prices b
  JOIN mlb_bets_combined c USING (bet_row_id)
  WHERE c.line = 0 AND c.market LIKE \"%spread%\"
  ORDER BY b.bet_row_id, b.side, b.bookmaker
  LIMIT 30
''').fetchdf()
print(df.to_string())
con.close()
"
```
Expected: For each pick'em bet, DK + FD rows (if they post winners) carry `is_exact_line = TRUE`, `derived_fair_odds` populated, `line_quoted = 0`. Books without a winner produce no row.

- [ ] **Step 12.3: Visual confirmation on the rendered dashboard**

Open the rendered report against the live DBs. Verify:
- A first-3-innings PK card now shows DK + FD prices on the FAIR toggle (and raw DNB on RAW), no `-0.5` line tag.
- A non-PK card (FG run line, FG alt spread, totals) is unchanged.

---

## Task 13: Documentation

**Files:**
- Modify: `Answer Keys/CLAUDE.md`
- Modify: `mlb_sgp/README.md`

- [ ] **Step 13.1: Update `Answer Keys/CLAUDE.md`**

Add a subsection under "MLB Dashboard — Odds screen" describing the pick'em→winner matching rule, source priority (2-way → 3-way devig-drop-tie → 0-handicap → "—"), and the `derived_fair_odds` column. Reference this spec by path.

- [ ] **Step 13.2: Update `mlb_sgp/README.md`**

Add a line under DK/FD coverage noting that bare "1st N Innings" period winners are now captured and flowed as moneylines so the dashboard pick'em path can render them.

- [ ] **Step 13.3: Commit docs**

```bash
git add "Answer Keys/CLAUDE.md" mlb_sgp/README.md
git commit -m "docs: pick'em (line-0 spread) -> moneyline matching

Document the new matching rule, source priority, derived_fair_odds column,
and the DK/FD bare-period-name winner capture.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 14: Pre-merge review

**Files:** none (review only).

- [ ] **Step 14.1: Diff against `main`**

```bash
cd /Users/callancapitolo/NFLWork/.claude/worktrees/feat+pickem-moneyline-match && git diff main..HEAD --stat
git diff main..HEAD
```

- [ ] **Step 14.2: Run the executive-engineer review checklist (per CLAUDE.md)**

For each category, write a short verdict (PASS or ISSUE):
- Data integrity (no duplicate writes; `derived_fair_odds` NULL on all non-PK rows).
- Resource safety (no new long-held DB connections).
- Edge cases (no winner markets at all; degenerate 3-way odds; non-MLB sports unaffected; first run with the new column).
- Dead code (any unused arg/branch in `render_book_cell`?).
- Log/disk hygiene (no new unbounded logs).
- Security (no new secrets surfaced in logs).

- [ ] **Step 14.3: Stop and report to user**

Do not merge. Output the diff summary and the checklist verdicts. Ask the user to approve merge.

---

## Version control & worktree lifecycle

- **Branch:** `feature/pickem-moneyline-match` (cut from local `main` `20dfcb3`).
- **Worktree:** `.claude/worktrees/feat+pickem-moneyline-match`.
- **Commit cadence:** one commit per task (sizes above). All commit messages end with `Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>`.
- **Merge:** only after Task 14 user approval. After merge, clean up: `git worktree remove .claude/worktrees/feat+pickem-moneyline-match` + `git branch -d feature/pickem-moneyline-match`.
- **Parallel branch:** `worktree-fix+dk-alt-spread-bucket-collapse` (FG abs-collapse fix) remains parked. Before merging *that* branch, re-apply its 2-file edit onto local `main` since its current base (`origin/main`) is stale.

## Documentation updates

Covered in Task 13. Required edits:
- `Answer Keys/CLAUDE.md`: pick'em matching rule + `derived_fair_odds`.
- `mlb_sgp/README.md`: DK + FD period-winner capture.

No new files; no schema migration script needed (existing `DROP TABLE IF EXISTS` + `CREATE TABLE` in MLB.R Phase 8 rebuilds the table each run, picking up the new column automatically).
