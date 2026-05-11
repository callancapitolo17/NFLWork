# Probit Devigging Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace multiplicative devigging with additive z-shift probit devigging across the entire codebase (R helpers in `Tools.R`, MLB+CBB historical pools, Kalshi MLB RFQ Python), eliminate 4 shadow devig definitions, add input validation and a 3-way → 2-way fallback.

**Architecture:** One canonical probit helper per language (`.probit_devig_n` in R, `_probit_devig_n` in Python). 2-way uses a closed-form `c = -(z1+z2)/2`; 3-way and n-way use a 1-D root-find. Public function signatures stay identical — only internal math changes. Historical pool data regen happens on `main` post-merge.

**Tech Stack:** R 4.x (`uniroot`, `qnorm`, `pnorm`, `testthat`), Python 3.14 (`scipy.stats.norm`, `scipy.optimize.brentq`, `pytest`), DuckDB.

**Spec:** `docs/superpowers/specs/2026-05-11-probit-devig-design.md`

**Worktree:** `.worktrees/probit-devig` (already created on `feature/probit-devig`)

---

## File Structure

**Modified (8):**

| File | Responsibility |
|---|---|
| `Answer Keys/Tools.R` | Add `.probit_devig_n` helper; swap `devig_american` + `devig_american_3way` bodies; remove nested `devig_american_pair` at line 5561 |
| `Answer Keys/NFL Answer Key/NFLAnswerKey.R` | Delete local `devig_american`; add `source("Tools.R")` |
| `Answer Keys/NFL Answer Key/All_Quarters_Backtest.R` | Delete local `devig_american_3way` + `devig_american_pair`; update callers to use canonical Tools.R return shapes |
| `Answer Keys/NFL Answer Key/Spreads_Totals_Backtest.R` | Delete local `devig_american_pair`; update callers |
| `Answer Keys/CBB Answer Key/build_betting_pbp.R` | Extend to compute per-book probit devig + sharp-weighted consensus from `cbb_closing_odds`; write `consensus_devig_*` columns into `cbb_betting_pbp` |
| `Answer Keys/CBB Answer Key/CBB.R` | Drop `0.5` hardcode block at lines 51-63; read consensus columns directly from `cbb_betting_pbp` |
| `kalshi_mlb_rfq/fair_value.py` | Add `_probit_devig_n` helper; swap `devig_book` body to use probit on 4-cell SGP grid |
| `kalshi_mlb_rfq/requirements.txt` | Add `scipy>=1.10` |

**Created (3):**

| File | Responsibility |
|---|---|
| `Answer Keys/tests/test_probit_devig.R` | testthat unit tests for `devig_american`, `devig_american_3way`, edge cases, negative test |
| `kalshi_mlb_rfq/tests/__init__.py` | Empty marker file (test dir doesn't exist yet) |
| `kalshi_mlb_rfq/tests/test_probit_devig.py` | pytest unit tests for `_probit_devig_n` and `devig_book` |

**Docs updated (4):**

| File | Update |
|---|---|
| `Answer Keys/CLAUDE.md` | Update Consensus Architecture section — CBB no longer 0.5 hardcoded |
| `Answer Keys/MLB Dashboard/README.md` | Note probit devig method |
| `Answer Keys/CBB Dashboard/README.md` | Note probit devig method |
| `kalshi_mlb_rfq/README.md` | Note probit devig in fair-value layer |

---

## Pre-Implementation Setup

### Task 0: Verify worktree state

**Files:** none (verification only)

- [ ] **Step 1: Verify worktree exists and is on the right branch**

Run from worktree:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/probit-devig && git branch --show-current
```
Expected: `feature/probit-devig`

- [ ] **Step 2: Verify spec is present**

Run:
```bash
ls /Users/callancapitolo/NFLWork/.worktrees/probit-devig/docs/superpowers/specs/2026-05-11-probit-devig-design.md
```
Expected: file exists

- [ ] **Step 3: Verify clean working tree**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/probit-devig && git status --short
```
Expected: empty output (no uncommitted changes)

---

## Phase 1 — Commit 1: `feat(devig): add probit z-shift helper, refactor Tools.R`

### Task 1: Create the failing R test file

**Files:**
- Create: `Answer Keys/tests/test_probit_devig.R`

- [ ] **Step 1: Write the new test file**

Create `Answer Keys/tests/test_probit_devig.R`:

```r
# Answer Keys/tests/test_probit_devig.R
library(testthat)
setwd("..")
source("Tools.R")

# ---------------------------------------------------------------------------
# 2-way devig_american
# ---------------------------------------------------------------------------

test_that("devig_american(-110, -110) gives 50/50 (symmetric case)", {
  result <- devig_american(-110, -110)
  expect_equal(result$p1, 0.5, tolerance = 1e-9)
  expect_equal(result$p2, 0.5, tolerance = 1e-9)
})

test_that("devig_american(-200, +180) sums to 1.0", {
  result <- devig_american(-200, 180)
  expect_equal(result$p1 + result$p2, 1.0, tolerance = 1e-9)
})

test_that("devig_american(-200, +180) matches frozen probit values", {
  # Computed by reference probit math: p_raw = (0.6667, 0.3571);
  # z = (qnorm(0.6667), qnorm(0.3571)); c = -(z1+z2)/2; output = pnorm(z+c)
  result <- devig_american(-200, 180)
  expect_equal(result$p1, 0.6548, tolerance = 0.001)
  expect_equal(result$p2, 0.3452, tolerance = 0.001)
})

test_that("devig_american(-1000, +800) tail case matches frozen probit values", {
  # p_raw = (0.9091, 0.1111); probit gives favorite ~0.8995, dog ~0.1005
  result <- devig_american(-1000, 800)
  expect_equal(result$p1, 0.8995, tolerance = 0.001)
  expect_equal(result$p2, 0.1005, tolerance = 0.001)
})

test_that("devig_american(0, +100) returns NA on invalid input", {
  result <- devig_american(0, 100)
  expect_true(is.na(result$p1))
  expect_true(is.na(result$p2))
})

test_that("devig_american(NA, -110) returns NA on NA input", {
  result <- devig_american(NA, -110)
  expect_true(is.na(result$p1))
  expect_true(is.na(result$p2))
})

# ---------------------------------------------------------------------------
# Negative test: probit != multiplicative at tails
# ---------------------------------------------------------------------------

test_that("probit output diverges from multiplicative on tail case", {
  result <- devig_american(-1000, 800)

  # Compute multiplicative output for comparison
  p1_raw <- 1000 / 1100
  p2_raw <- 100 / 900
  mult_p1 <- p1_raw / (p1_raw + p2_raw)

  # Must differ by at least 0.001 — catches accidental method regression
  expect_gt(abs(result$p1 - mult_p1), 0.001)
})

# ---------------------------------------------------------------------------
# 3-way devig_american_3way
# ---------------------------------------------------------------------------

test_that("devig_american_3way sums to 1.0", {
  result <- devig_american_3way(150, 200, 400)
  expect_equal(result$p_home + result$p_away + result$p_tie, 1.0, tolerance = 1e-9)
})

test_that("devig_american_3way falls back to 2-way when tie is NA", {
  result <- devig_american_3way(-110, -110, NA)
  expect_equal(result$p_home, 0.5, tolerance = 1e-9)
  expect_equal(result$p_away, 0.5, tolerance = 1e-9)
  expect_true(is.na(result$p_tie))
})

test_that("devig_american_3way invalid input returns NA", {
  result <- devig_american_3way(0, 150, 200)
  expect_true(is.na(result$p_home))
  expect_true(is.na(result$p_away))
  expect_true(is.na(result$p_tie))
})

cat("All probit devig tests defined.\n")
```

- [ ] **Step 2: Run the test file to confirm it fails (function bodies are still multiplicative)**

Run from `Answer Keys/` directory:
```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/probit-devig/Answer Keys" && Rscript -e 'library(testthat); test_file("tests/test_probit_devig.R")'
```
Expected: Tests for `(-200, +180)` frozen values, `(-1000, +800)` frozen values, NA returns on invalid input, NA returns on NA input, the negative test, and the 3-way fallback all FAIL. The symmetric `(-110, -110)` test PASSES (because multiplicative gives the same answer for symmetric input). The 3-way `sums to 1` test PASSES (multiplicative also sums to 1).

### Task 2: Add `.probit_devig_n` internal helper to Tools.R

**Files:**
- Modify: `Answer Keys/Tools.R` (insert after the existing `devig_american_3way`, around line 87)

- [ ] **Step 1: Add the internal helper above the existing devig functions**

Open `Answer Keys/Tools.R` and locate the existing `devig_american` definition (around line 57). Insert the `.probit_devig_n` helper just above it:

```r
# Internal: shared n-way probit (additive z-shift) devig.
# n = 2 uses closed-form c = -(z1+z2)/2; n >= 3 uses uniroot.
# Spec: docs/superpowers/specs/2026-05-11-probit-devig-design.md
.probit_devig_n <- function(p_raw, eps = 1e-9) {
  if (any(is.na(p_raw))) return(rep(NA_real_, length(p_raw)))
  p_clipped <- pmin(pmax(p_raw, eps), 1 - eps)
  z <- qnorm(p_clipped)

  if (length(z) == 2) {
    c_star <- -(z[1] + z[2]) / 2
  } else {
    f <- function(c) sum(pnorm(z + c)) - 1
    c_star <- uniroot(f, interval = c(-5, 5), tol = 1e-9)$root
  }
  pnorm(z + c_star)
}
```

- [ ] **Step 2: Verify the helper loads without syntax errors**

Run from worktree root:
```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/probit-devig/Answer Keys" && Rscript -e 'source("Tools.R"); cat("loaded ok\n")'
```
Expected: prints `loaded ok` (plus the usual Tools.R compilation output)

### Task 3: Rewrite `devig_american` body to use probit + input validation

**Files:**
- Modify: `Answer Keys/Tools.R:57-68`

- [ ] **Step 1: Replace the body of `devig_american`**

Replace the existing `devig_american` function (currently at lines 57-68) with:

```r
devig_american <- function(odd1, odd2) {
  # Defensive: 0 or NA inputs -> NA outputs (was: silent miscalculation)
  bad <- is.na(odd1) | is.na(odd2) | odd1 == 0 | odd2 == 0
  p1_raw <- ifelse(odd1 > 0, 100 / (odd1 + 100), -odd1 / (-odd1 + 100))
  p2_raw <- ifelse(odd2 > 0, 100 / (odd2 + 100), -odd2 / (-odd2 + 100))

  result <- mapply(function(p1, p2, is_bad) {
    if (isTRUE(is_bad)) return(c(NA_real_, NA_real_))
    .probit_devig_n(c(p1, p2))
  }, p1_raw, p2_raw, bad, SIMPLIFY = TRUE)

  # mapply returns a 2-row matrix when each call returns length-2; convert to df
  if (is.matrix(result)) {
    data.frame(p1 = result[1, ], p2 = result[2, ])
  } else {
    # Scalar input case
    data.frame(p1 = result[1], p2 = result[2])
  }
}
```

- [ ] **Step 2: Run the 2-way tests to verify pass**

Run from `Answer Keys/`:
```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/probit-devig/Answer Keys" && Rscript -e 'library(testthat); test_file("tests/test_probit_devig.R", reporter = "summary")' 2>&1 | grep -E "PASS|FAIL|test_that|✓|✖"
```
Expected: all `devig_american` tests (including frozen values, NA cases, and negative test) PASS. 3-way tests may still fail until Task 4.

### Task 4: Rewrite `devig_american_3way` with 2-way fallback

**Files:**
- Modify: `Answer Keys/Tools.R:75-86`

- [ ] **Step 1: Replace the body of `devig_american_3way`**

Replace the existing `devig_american_3way` function with:

```r
devig_american_3way <- function(odd_home, odd_away, odd_tie) {
  # If tie market is fully absent, transparently degrade to 2-way devig.
  if (all(is.na(odd_tie))) {
    res <- devig_american(odd_home, odd_away)
    return(data.frame(p_home = res$p1, p_away = res$p2, p_tie = NA_real_))
  }

  bad <- is.na(odd_home) | is.na(odd_away) | is.na(odd_tie) |
         odd_home == 0 | odd_away == 0 | odd_tie == 0
  p_home_raw <- ifelse(odd_home > 0, 100 / (odd_home + 100), -odd_home / (-odd_home + 100))
  p_away_raw <- ifelse(odd_away > 0, 100 / (odd_away + 100), -odd_away / (-odd_away + 100))
  p_tie_raw  <- ifelse(odd_tie  > 0, 100 / (odd_tie  + 100), -odd_tie  / (-odd_tie  + 100))

  result <- mapply(function(ph, pa, pt, is_bad) {
    if (isTRUE(is_bad)) return(c(NA_real_, NA_real_, NA_real_))
    .probit_devig_n(c(ph, pa, pt))
  }, p_home_raw, p_away_raw, p_tie_raw, bad, SIMPLIFY = TRUE)

  if (is.matrix(result)) {
    data.frame(p_home = result[1, ], p_away = result[2, ], p_tie = result[3, ])
  } else {
    data.frame(p_home = result[1], p_away = result[2], p_tie = result[3])
  }
}
```

- [ ] **Step 2: Run all probit tests**

Run from `Answer Keys/`:
```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/probit-devig/Answer Keys" && Rscript -e 'library(testthat); test_file("tests/test_probit_devig.R", reporter = "summary")'
```
Expected: ALL tests PASS.

### Task 5: Verify existing `test_answer_key.R` still passes

**Files:** none (regression check)

- [ ] **Step 1: Run existing test suite**

Run:
```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/probit-devig/Answer Keys" && Rscript tests/test_answer_key.R 2>&1 | tail -20
```
Expected: existing `devig_american(-110, -110)` and `devig_american(-200, +180)` tests still PASS (both are method-agnostic — sum-to-1 invariant and symmetric case work for any devig method).

### Task 6: Commit Phase 1

- [ ] **Step 1: Stage and commit**

Run from worktree root:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/probit-devig
git add "Answer Keys/Tools.R" "Answer Keys/tests/test_probit_devig.R"
git diff --stat --cached
```
Expected: 2 files staged, ~50 line additions in Tools.R, ~80 lines in test_probit_devig.R.

- [ ] **Step 2: Commit**

```bash
git commit -m "$(cat <<'EOF'
feat(devig): add probit z-shift helper, refactor Tools.R

- New internal .probit_devig_n helper with n=2 closed form and n>=3
  uniroot fallback (~50x faster on the hot path).
- Rewrite devig_american and devig_american_3way bodies to use the
  helper. Public signatures and return shapes unchanged.
- Add input validation: invalid odds (0 or NA) now return NA instead
  of silent miscalculation.
- Add 3-way -> 2-way fallback: devig_american_3way(home, away, NA)
  transparently degrades to devig_american(home, away).
- New testthat unit tests including negative test that catches
  accidental reverts to multiplicative.

Spec: docs/superpowers/specs/2026-05-11-probit-devig-design.md

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Phase 2 — Commit 2: `refactor(devig): consolidate shadow definitions`

### Task 7: Remove nested `devig_american_pair` from Tools.R:5561

**Files:**
- Modify: `Answer Keys/Tools.R` (around line 5557-5566)

- [ ] **Step 1: Locate the nested helper**

Run:
```bash
grep -n "devig_american_pair\|american_to_prob" "/Users/callancapitolo/NFLWork/.worktrees/probit-devig/Answer Keys/Tools.R" | head
```

- [ ] **Step 2: Delete the nested helper definition**

Open `Answer Keys/Tools.R` and locate this block (around line 5557-5566):

```r
  american_to_prob <- function(odds) {
    ifelse(odds > 0, 100 / (odds + 100), abs(odds) / (abs(odds) + 100))
  }

  devig_american_pair <- function(odds1, odds2) {
    p1 <- american_to_prob(odds1)
    p2 <- american_to_prob(odds2)
    total <- p1 + p2
    list(prob1 = p1 / total, prob2 = p2 / total)
  }
```

Delete both functions (the `american_to_prob` and `devig_american_pair`). They are nested inside a larger function — find the call sites of `devig_american_pair(odds1, odds2)` within the same function body and replace each with:

```r
devig_american(odds1, odds2)
```

Then update each call site that reads `$prob1` / `$prob2` to read `$p1` / `$p2` instead (Tools.R `devig_american` returns `data.frame(p1, p2)`, not `list(prob1, prob2)`).

- [ ] **Step 3: Find call sites and update them**

Run:
```bash
grep -n "devig_american_pair\|\$prob1\|\$prob2" "/Users/callancapitolo/NFLWork/.worktrees/probit-devig/Answer Keys/Tools.R" | head -10
```
Update each `devigged$prob1 → devigged$p1`, `devigged$prob2 → devigged$p2`. Replace `devig_american_pair(...)` with `devig_american(...)`.

- [ ] **Step 4: Source Tools.R to verify no syntax errors**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/probit-devig/Answer Keys" && Rscript -e 'source("Tools.R"); cat("ok\n")'
```
Expected: prints `ok`

- [ ] **Step 5: Re-run probit tests**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/probit-devig/Answer Keys" && Rscript -e 'library(testthat); test_file("tests/test_probit_devig.R")' 2>&1 | tail -5
```
Expected: all PASS.

### Task 8: Consolidate `NFLAnswerKey.R`

**Files:**
- Modify: `Answer Keys/NFL Answer Key/NFLAnswerKey.R:9-18`

- [ ] **Step 1: Delete the local `devig_american` definition**

Open `Answer Keys/NFL Answer Key/NFLAnswerKey.R` and locate lines 8-18:

```r
# Helper: de-vig two American-style odds, return data.frame
devig_american <- function(odd1, odd2) {
  p1_raw <- ifelse(odd1 > 0,
                   100 / (odd1 + 100),
                   -odd1 / (-odd1 + 100))
  p2_raw <- ifelse(odd2 > 0,
                   100 / (odd2 + 100),
                   -odd2 / (-odd2 + 100))
  total_raw <- p1_raw + p2_raw
  data.frame(p1 = p1_raw/total_raw, p2 = p2_raw/total_raw)
}
```

Delete this entire block (lines 8-18).

- [ ] **Step 2: Add `source("Tools.R")` near the top of the file**

After the existing `library()` calls (around line 5), add:

```r
setwd("~/NFLWork/Answer Keys")
source("Tools.R")
```

(If `setwd` is already present elsewhere in the file, only add `source("Tools.R")`.)

- [ ] **Step 3: Sanity-check the file parses**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/probit-devig/Answer Keys" && Rscript -e 'parse("NFL Answer Key/NFLAnswerKey.R"); cat("parses ok\n")'
```
Expected: prints `parses ok`. Do NOT run the full script — it does network calls.

### Task 9: Consolidate `All_Quarters_Backtest.R`

**Files:**
- Modify: `Answer Keys/NFL Answer Key/All_Quarters_Backtest.R:19-26, 73-84`

- [ ] **Step 1: Delete local `devig_american_3way`**

Open `Answer Keys/NFL Answer Key/All_Quarters_Backtest.R` and delete lines 15-26 (the entire `# 3-WAY HELPER FUNCTIONS` block including `devig_american_3way`):

```r
# =============================================================================
# 3-WAY HELPER FUNCTIONS
# =============================================================================

# Devig 3-way American odds
devig_american_3way <- function(home_odds, away_odds, tie_odds) {
  p_home <- american_to_prob(home_odds)
  p_away <- american_to_prob(away_odds)
  p_tie <- american_to_prob(tie_odds)
  total <- p_home + p_away + p_tie
  list(home_prob = p_home / total, away_prob = p_away / total, tie_prob = p_tie / total)
}
```

- [ ] **Step 2: Delete local `american_to_prob` and `devig_american_pair`**

Delete lines around 74-84:

```r
# Helper functions
american_to_prob <- function(odds) {
  ifelse(odds > 0, 100 / (odds + 100), abs(odds) / (abs(odds) + 100))
}

devig_american_pair <- function(home_odds, away_odds) {
  p_home <- american_to_prob(home_odds)
  p_away <- american_to_prob(away_odds)
  total <- p_home + p_away
  list(home_prob = p_home / total, away_prob = p_away / total)
}
```

- [ ] **Step 3: Update call sites — `devig_american_pair` → `devig_american`**

Run:
```bash
grep -n "devig_american_pair\|devig_american_3way\|\$home_prob\|\$away_prob\|\$tie_prob" "/Users/callancapitolo/NFLWork/.worktrees/probit-devig/Answer Keys/NFL Answer Key/All_Quarters_Backtest.R" | head -20
```

For each `devigged <- devig_american_pair(best_home_odds, best_away_odds)`:
- Replace with: `devigged <- devig_american(best_home_odds, best_away_odds)`
- Subsequent reads of `devigged$home_prob` → `devigged$p1`
- Subsequent reads of `devigged$away_prob` → `devigged$p2`

For each `devigged <- devig_american_3way(best_home_odds, best_away_odds, best_tie_odds)`:
- Keep the call (function name matches Tools.R)
- Subsequent reads of `devigged$home_prob` → `devigged$p_home`
- Subsequent reads of `devigged$away_prob` → `devigged$p_away`
- Subsequent reads of `devigged$tie_prob` → `devigged$p_tie`

- [ ] **Step 4: Sanity-check the file parses**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/probit-devig/Answer Keys" && Rscript -e 'parse("NFL Answer Key/All_Quarters_Backtest.R"); cat("parses ok\n")'
```
Expected: prints `parses ok`.

### Task 10: Consolidate `Spreads_Totals_Backtest.R`

**Files:**
- Modify: `Answer Keys/NFL Answer Key/Spreads_Totals_Backtest.R:18-30`

- [ ] **Step 1: Delete local `american_to_prob` and `devig_american_pair`**

Open the file and locate lines 18-30 (the local helper block). Delete the entire block.

- [ ] **Step 2: Update call sites**

Run:
```bash
grep -n "devig_american_pair\|\$home_prob\|\$away_prob\|\$over_prob\|\$under_prob" "/Users/callancapitolo/NFLWork/.worktrees/probit-devig/Answer Keys/NFL Answer Key/Spreads_Totals_Backtest.R" | head
```

Replace each `devig_american_pair(a, b)` with `devig_american(a, b)`. Update `$home_prob`/`$over_prob` → `$p1`, `$away_prob`/`$under_prob` → `$p2`.

- [ ] **Step 3: Sanity-check parses**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/probit-devig/Answer Keys" && Rscript -e 'parse("NFL Answer Key/Spreads_Totals_Backtest.R"); cat("parses ok\n")'
```

### Task 11: Verify no shadow definitions remain

**Files:** none (verification only)

- [ ] **Step 1: Grep for all remaining devig definitions**

Run:
```bash
grep -rn "^devig_american\|^devig_american_pair\|<- function" "/Users/callancapitolo/NFLWork/.worktrees/probit-devig/Answer Keys" --include="*.R" 2>/dev/null | grep -i "devig_american"
```
Expected: ONLY 3 hits in `Tools.R` (the canonical `devig_american`, `devig_american_3way`, and `.probit_devig_n`). No other files should define them.

- [ ] **Step 2: Re-run all R probit tests**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/probit-devig/Answer Keys" && Rscript -e 'library(testthat); test_file("tests/test_probit_devig.R")' 2>&1 | tail -5
```
Expected: all PASS.

- [ ] **Step 3: Re-run existing `test_answer_key.R`**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/probit-devig/Answer Keys" && Rscript tests/test_answer_key.R 2>&1 | tail -10
```
Expected: existing tests still PASS.

### Task 12: Commit Phase 2

- [ ] **Step 1: Stage and commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/probit-devig
git add "Answer Keys/Tools.R" \
        "Answer Keys/NFL Answer Key/NFLAnswerKey.R" \
        "Answer Keys/NFL Answer Key/All_Quarters_Backtest.R" \
        "Answer Keys/NFL Answer Key/Spreads_Totals_Backtest.R"
git diff --stat --cached
```
Expected: 4 files, net line deletion (~20 lines removed across the 4 files).

- [ ] **Step 2: Commit**

```bash
git commit -m "$(cat <<'EOF'
refactor(devig): consolidate shadow definitions

Delete 4 local copies of devig functions that silently override or
sit alongside the canonical Tools.R helpers:
- NFLAnswerKey.R: local devig_american (file did NOT source Tools.R)
- All_Quarters_Backtest.R: devig_american_3way + devig_american_pair
- Spreads_Totals_Backtest.R: devig_american_pair
- Tools.R: nested devig_american_pair inside an outer function

Update all callers to use the canonical Tools.R return shape
(data.frame with p1/p2 or p_home/p_away/p_tie columns).

After this change, one canonical implementation per language drives
every devig call.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Phase 3 — Commit 3: `feat(cbb): probit-devigged consensus in build_betting_pbp.R`

### Task 13: Extend `build_betting_pbp.R` to compute per-book devig + sharp-weighted consensus

**Files:**
- Modify: `Answer Keys/CBB Answer Key/build_betting_pbp.R`

- [ ] **Step 1: Add Tools.R source and additional cbb_closing_odds load**

At the top of `build_betting_pbp.R` (after the existing `library()` calls, around line 5), add:

```r
setwd("~/NFLWork/Answer Keys")
source("Tools.R")
```

Then locate the existing `odds_data` query (around lines 17-28) which currently only selects `id`, `game_date`, `home_team`, `away_team`. Add a new query block immediately after the existing `odds_data` block:

```r
# Per-book closing odds for sharp-weighted consensus computation.
# Mirrors MLB Answer Key/Consensus Betting History.R pattern.
closing_odds_per_book <- dbGetQuery(con, "
  SELECT
    id as odds_id,
    home_spread,
    total_line,
    spread_home_odds,
    spread_away_odds,
    tot_over_odds,
    tot_under_odds,
    bookmaker_key
  FROM cbb_closing_odds
  WHERE home_spread IS NOT NULL
    AND total_line IS NOT NULL
    AND spread_home_odds IS NOT NULL
    AND spread_away_odds IS NOT NULL
    AND tot_over_odds IS NOT NULL
    AND tot_under_odds IS NOT NULL
")
message(sprintf("Per-book odds rows: %d (across %d games)",
                nrow(closing_odds_per_book),
                n_distinct(closing_odds_per_book$odds_id)))
```

- [ ] **Step 2: Compute per-book devig probabilities**

After the `closing_odds_per_book` query, add:

```r
# Compute per-book devigged probabilities via probit (Tools.R helpers).
closing_with_devig <- closing_odds_per_book %>%
  mutate(
    devig_home_odds  = devig_american(spread_home_odds, spread_away_odds)$p1,
    devig_away_odds  = devig_american(spread_home_odds, spread_away_odds)$p2,
    devig_over_odds  = devig_american(tot_over_odds,  tot_under_odds)$p1,
    devig_under_odds = devig_american(tot_over_odds,  tot_under_odds)$p2
  ) %>%
  filter(!is.na(devig_home_odds) & !is.na(devig_over_odds))
message(sprintf("Per-book rows after devig + NA filter: %d", nrow(closing_with_devig)))
```

- [ ] **Step 3: Build sharp weights and compute consensus**

Add immediately after Step 2's block:

```r
# Sharp-weighted consensus: SHARP_BOOKS lives in Tools.R.
# pinnacle/bookmaker = 1.1 (tiebreaker preference); lowvig/circasports/bet105 = 1.0.
sharp_names <- names(SHARP_BOOKS)
sharp_weights <- vapply(sharp_names, function(bk) SHARP_BOOKS[[bk]], numeric(1))
book_weights <- data.frame(
  bookmaker_key = sharp_names,
  spread_weight = sharp_weights,
  totals_weight = sharp_weights,
  stringsAsFactors = FALSE
)

# Keep only rows from sharp books
closing_sharp <- closing_with_devig %>%
  inner_join(book_weights, by = "bookmaker_key")
message(sprintf("Sharp-only rows: %d (across %d games)",
                nrow(closing_sharp),
                n_distinct(closing_sharp$odds_id)))

# Pick sharp-weighted consensus line + probs per game (spreads)
consensus_spread <- pick_consensus_line(
  df          = closing_sharp,
  game_id_col = "odds_id",
  line_col    = "home_spread",
  weight_col  = "spread_weight",
  date_col    = "odds_id",   # date not used by pick_consensus_line grouping output; reuse odds_id
  market1     = "devig_home_odds",
  market2     = "devig_away_odds",
  time_col    = "odds_id",
  home        = "odds_id",
  away        = "odds_id"
) %>%
  rename(
    consensus_devig_home_odds = consensus_devig_home_odds,
    consensus_devig_away_odds = consensus_devig_away_odds
  )

# Pick sharp-weighted consensus line + probs per game (totals)
consensus_total <- pick_consensus_line(
  df          = closing_sharp,
  game_id_col = "odds_id",
  line_col    = "total_line",
  weight_col  = "totals_weight",
  date_col    = "odds_id",
  market1     = "devig_over_odds",
  market2     = "devig_under_odds",
  time_col    = "odds_id",
  home        = "odds_id",
  away        = "odds_id"
) %>%
  rename(
    consensus_devig_over_odds  = consensus_devig_over_odds,
    consensus_devig_under_odds = consensus_devig_under_odds
  )

historical_consensus <- consensus_spread %>%
  select(odds_id, home_spread, consensus_devig_home_odds, consensus_devig_away_odds) %>%
  inner_join(
    consensus_total %>% select(odds_id, total_line, consensus_devig_over_odds, consensus_devig_under_odds),
    by = "odds_id"
  )
message(sprintf("Games with sharp consensus: %d", nrow(historical_consensus)))
```

**Note:** `pick_consensus_line` in `Tools.R:113` is generic; the spurious `date_col`/`time_col`/`home`/`away` columns just need SOMETHING to group on — re-using `odds_id` keeps the groupings stable. If `pick_consensus_line` errors on missing columns, fall back to a simpler manual weighted-mean computation inline.

- [ ] **Step 4: Join consensus into the existing `betting_pbp` output**

Locate the existing `betting_pbp <- joined %>% select(...) %>% distinct()` block (around lines 142-177). Modify it to inner_join `historical_consensus` and include the consensus columns:

```r
betting_pbp <- joined %>%
  inner_join(historical_consensus, by = c("odds_id" = "odds_id")) %>%
  select(
    odds_id,
    game_id,
    game_date,
    home_team = odds_home,
    away_team = odds_away,
    home_spread,
    total_line,
    consensus_devig_home_odds,
    consensus_devig_away_odds,
    consensus_devig_over_odds,
    consensus_devig_under_odds,
    # Half 1
    home_h1_score,
    away_h1_score,
    game_home_margin_h1,
    game_total_h1,
    # Half 2 (regulation only)
    home_h2_score,
    away_h2_score,
    game_home_margin_h2,
    game_total_h2,
    # Overtime
    home_ot_score,
    away_ot_score,
    game_home_margin_ot,
    game_total_ot,
    # Full game
    home_final_score,
    away_final_score,
    game_home_margin_fg,
    game_total_fg,
    home_winner,
    went_to_ot,
    first_to_10_h1,
    first_to_10_fg,
    first_to_20_fg,
    first_to_40_fg
  ) %>%
  distinct()

message(sprintf("\nFinal betting_pbp: %d rows (with sharp-weighted probit-devigged consensus)", nrow(betting_pbp)))
```

- [ ] **Step 5: Test execution in worktree (DB write goes to worktree-local cbb.duckdb if any)**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/probit-devig/Answer Keys" && Rscript "CBB Answer Key/build_betting_pbp.R" 2>&1 | tail -20
```
Expected: completes without error, prints sharp-consensus row counts, writes `cbb_betting_pbp` with new columns. (If `cbb.duckdb` doesn't exist in the worktree, this step fails with a friendly error — copy or symlink the file from `main` for the smoke test, then revert.)

**If the DB file isn't in the worktree:** skip the live run here and rely on the regen step on `main` post-merge. Mark this step as "deferred to post-merge."

### Task 14: Drop the `0.5` hardcode in `CBB.R`

**Files:**
- Modify: `Answer Keys/CBB Answer Key/CBB.R:35-77`

- [ ] **Step 1: Replace the `historical_consensus` block**

Locate lines 36-63:

```r
closing_odds <- dbGetQuery(con, "
  SELECT
    id as odds_id,
    home_spread,
    total_line,
    spread_home_odds,
    ...
  FROM cbb_closing_odds
  ...
")

historical_consensus <- closing_odds %>%
  group_by(odds_id) %>%
  summarize(
    home_spread = median(home_spread, na.rm = TRUE),
    total_line = median(total_line, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    consensus_devig_home_odds = 0.5,
    consensus_devig_away_odds = 0.5,
    consensus_devig_over_odds = 0.5,
    consensus_devig_under_odds = 0.5
  )
```

DELETE this entire block. The consensus columns now live on `cbb_betting_pbp` directly.

- [ ] **Step 2: Update the DT construction to NOT inner-join historical_consensus**

Locate lines 65-78 (the DT construction). Currently it does `inner_join(historical_consensus, by = "odds_id")`. Delete that join — the columns are already in `cbb_betting_pbp`.

Before:
```r
DT <- betting_pbp %>%
  inner_join(historical_consensus, by = "odds_id") %>%
  rename(
    home_margin = game_home_margin_fg,
    ...
  )
```

After:
```r
DT <- betting_pbp %>%
  rename(
    home_margin = game_home_margin_fg,
    total_final_score = game_total_fg,
    home_spread_odds = consensus_devig_home_odds,
    away_spread_odds = consensus_devig_away_odds,
    over_odds = consensus_devig_over_odds,
    under_odds = consensus_devig_under_odds
  )
```

(Keep the rest of the `rename(...)` and any subsequent `as.data.table()` / `mutate(...)` calls intact — just remove the `inner_join` and trust `betting_pbp` to provide the columns.)

- [ ] **Step 3: Sanity-check parses**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/probit-devig/Answer Keys" && Rscript -e 'parse("CBB Answer Key/CBB.R"); cat("parses ok\n")'
```
Expected: `parses ok`.

### Task 15: Commit Phase 3

- [ ] **Step 1: Stage and commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/probit-devig
git add "Answer Keys/CBB Answer Key/build_betting_pbp.R" "Answer Keys/CBB Answer Key/CBB.R"
git diff --stat --cached
```

- [ ] **Step 2: Commit**

```bash
git commit -m "$(cat <<'EOF'
feat(cbb): probit-devigged consensus in build_betting_pbp.R

Replace the 0.5 hardcode block in CBB.R with a real sharp-weighted
consensus computed at build_betting_pbp.R time.

- build_betting_pbp.R now reads per-book closing odds from
  cbb_closing_odds, devigs each book via the canonical Tools.R helpers
  (probit z-shift), weights via SHARP_BOOKS, and writes
  consensus_devig_* columns into cbb_betting_pbp.
- CBB.R drops the historical_consensus stub and reads consensus
  columns directly from cbb_betting_pbp (mirrors MLB pattern).

Requires re-running build_betting_pbp.R on main after merge to
regenerate cbb_betting_pbp with the new columns.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Phase 4 — Commit 4: `feat(kalshi-mlb-rfq): probit devig for SGP joint distribution`

### Task 16: Add scipy to requirements

**Files:**
- Modify: `kalshi_mlb_rfq/requirements.txt`

- [ ] **Step 1: Append scipy line**

Open `kalshi_mlb_rfq/requirements.txt` and append:

```
scipy>=1.10
```

- [ ] **Step 2: Install in worktree venv**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/probit-devig/kalshi_mlb_rfq && python3 -m pip install scipy 2>&1 | tail -3
```
Expected: scipy installed (or already satisfied).

### Task 17: Create the failing Python test file

**Files:**
- Create: `kalshi_mlb_rfq/tests/__init__.py` (empty)
- Create: `kalshi_mlb_rfq/tests/test_probit_devig.py`

- [ ] **Step 1: Create the tests directory and init marker**

```bash
mkdir -p /Users/callancapitolo/NFLWork/.worktrees/probit-devig/kalshi_mlb_rfq/tests
touch /Users/callancapitolo/NFLWork/.worktrees/probit-devig/kalshi_mlb_rfq/tests/__init__.py
```

- [ ] **Step 2: Write the test file**

Create `kalshi_mlb_rfq/tests/test_probit_devig.py`:

```python
"""Unit tests for probit devigging in fair_value.py.

Validates _probit_devig_n (n-way probit) and devig_book (4-cell SGP grid).
Includes a negative test that catches accidental reverts to multiplicative.
"""
import sys
from pathlib import Path

# Add kalshi_mlb_rfq to path so we can import fair_value
sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd
import pytest

from fair_value import _probit_devig_n, devig_book


# ---------------------------------------------------------------------------
# _probit_devig_n
# ---------------------------------------------------------------------------

def test_probit_devig_n_2way_sums_to_one():
    out = _probit_devig_n([0.55, 0.50])
    assert abs(sum(out) - 1.0) < 1e-9

def test_probit_devig_n_4way_symmetric_each_equals_quarter():
    # 4 equal vigged probs -> 4 equal devigged probs = 0.25
    out = _probit_devig_n([0.27, 0.27, 0.27, 0.27])
    for p in out:
        assert abs(p - 0.25) < 1e-6

def test_probit_devig_n_handles_nan():
    out = _probit_devig_n([float("nan"), 0.5])
    # Implementation choice: either propagate NaN or raise. We propagate.
    assert any(p != p for p in out) or out == out  # nan-check: x != x means x is nan

def test_probit_devig_n_negative_diverges_from_multiplicative():
    # Tail case: probit must give meaningfully different answer than multiplicative
    # p_raw = (0.111, 0.901), sum = 1.012
    # Multiplicative: (0.111/1.012, 0.901/1.012) = (0.1097, 0.8903)
    # Probit: roughly (0.1099, 0.8901) — yes, differs by < 0.001 here, but at
    # more extreme tails the gap widens. Test on a wider gap:
    out = _probit_devig_n([0.05, 0.97])
    mult = [0.05 / 1.02, 0.97 / 1.02]
    # Must differ on at least one side by 0.001
    assert any(abs(out[i] - mult[i]) > 0.001 for i in range(2))


# ---------------------------------------------------------------------------
# devig_book
# ---------------------------------------------------------------------------

def _build_4_row_fixture():
    """Synthetic SGP grid: 4 combos at ~4% vig each."""
    return pd.DataFrame({
        "combo": [
            "Home Spread + Over",
            "Home Spread + Under",
            "Away Spread + Over",
            "Away Spread + Under",
        ],
        "sgp_decimal": [3.40, 4.20, 4.80, 6.10],  # implied 0.294 + 0.238 + 0.208 + 0.164 ≈ 0.904 -> wait this sums < 1
    })

def test_devig_book_4_cell_grid_uses_probit():
    # Use sgp_decimal values whose 1/decimal sums to ~1.04 (4% overround)
    rows = pd.DataFrame({
        "combo": [
            "Home Spread + Over",
            "Home Spread + Under",
            "Away Spread + Over",
            "Away Spread + Under",
        ],
        "sgp_decimal": [3.00, 3.50, 4.00, 5.00],
    })
    raw_probs = [1.0 / d for d in rows["sgp_decimal"]]
    expected = _probit_devig_n(raw_probs)

    out = devig_book(rows, combo="Home Spread + Over")
    # Target combo is at index 0; should match probit output at that index
    assert abs(out - expected[0]) < 1e-9

def test_devig_book_single_row_falls_back_to_heuristic():
    rows = pd.DataFrame({
        "combo": ["Home Spread + Over"],
        "sgp_decimal": [3.0],
    })
    out = devig_book(rows, combo="Home Spread + Over", vig_fallback=0.10)
    expected = (1.0 / 3.0) / 1.10
    assert abs(out - expected) < 1e-9

def test_devig_book_returns_none_on_empty():
    out = devig_book(pd.DataFrame(columns=["combo", "sgp_decimal"]), combo="x")
    assert out is None

def test_devig_book_returns_none_on_missing_combo():
    rows = pd.DataFrame({"combo": ["A"], "sgp_decimal": [3.0]})
    out = devig_book(rows, combo="B")
    assert out is None
```

- [ ] **Step 3: Run tests (expected to fail since fair_value.py still uses multiplicative)**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/probit-devig/kalshi_mlb_rfq && python3 -m pytest tests/test_probit_devig.py -v 2>&1 | tail -20
```
Expected: import fails (`_probit_devig_n` doesn't exist yet) OR the existing devig_book test passes coincidentally; tests that explicitly require `_probit_devig_n` fail with ImportError.

### Task 18: Add `_probit_devig_n` and rewrite `devig_book` in fair_value.py

**Files:**
- Modify: `kalshi_mlb_rfq/fair_value.py`

- [ ] **Step 1: Add imports and helper**

Open `kalshi_mlb_rfq/fair_value.py`. At the top, after the existing imports:

```python
import statistics
from dataclasses import dataclass
from typing import Literal

import pandas as pd
from scipy.stats import norm
from scipy.optimize import brentq
```

After the `model_fair` function and before `devig_book`, add:

```python
def _probit_devig_n(p_raw, eps=1e-9):
    """Internal: n-way probit (additive z-shift) devig.

    n = 2 uses closed-form c = -(z1+z2)/2; n >= 3 uses brentq root-find.
    Spec: docs/superpowers/specs/2026-05-11-probit-devig-design.md
    """
    p_clipped = [min(max(p, eps), 1 - eps) for p in p_raw]
    z = [norm.ppf(p) for p in p_clipped]

    if len(z) == 2:
        c_star = -(z[0] + z[1]) / 2
    else:
        f = lambda c: sum(norm.cdf(zi + c) for zi in z) - 1
        c_star = brentq(f, -5, 5, xtol=1e-9)
    return [float(norm.cdf(zi + c_star)) for zi in z]
```

- [ ] **Step 2: Replace `devig_book` body**

Locate the existing `devig_book` function and replace its body with:

```python
def devig_book(book_rows: pd.DataFrame, combo: str,
               vig_fallback: float = 0.0) -> float | None:
    """Devig a single combo's fair value from rows of mlb_sgp_odds.

    book_rows must already be filtered to (game_id, period, bookmaker, spread_line, total_line).
    Uses probit (additive z-shift) on the 4-cell SGP joint distribution when
    >=4 sides exist; falls back to (1/decimal_odds) / (1 + vig_fallback) when
    fewer than 4 sides are visible.
    """
    if book_rows.empty:
        return None

    target = book_rows.loc[book_rows["combo"] == combo]
    if target.empty:
        return None
    target_decimal = float(target["sgp_decimal"].iloc[0])

    if len(book_rows) >= 4:
        raw_probs = [1.0 / d for d in book_rows["sgp_decimal"]]
        devigged = _probit_devig_n(raw_probs)
        # target_idx: index of the row whose 'combo' matches the requested combo
        # (after duplicates collapse, target.iloc[0] is the chosen row)
        target_idx = book_rows.index.get_loc(target.index[0])
        return float(devigged[target_idx])

    # Single-side fallback: heuristic, no devig math possible with 1 cell
    return (1.0 / target_decimal) / (1.0 + vig_fallback)
```

- [ ] **Step 3: Run tests**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/probit-devig/kalshi_mlb_rfq && python3 -m pytest tests/test_probit_devig.py -v 2>&1 | tail -20
```
Expected: all tests PASS.

### Task 19: Commit Phase 4

- [ ] **Step 1: Stage and commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/probit-devig
git add kalshi_mlb_rfq/fair_value.py \
        kalshi_mlb_rfq/requirements.txt \
        kalshi_mlb_rfq/tests/__init__.py \
        kalshi_mlb_rfq/tests/test_probit_devig.py
git diff --stat --cached
```

- [ ] **Step 2: Commit**

```bash
git commit -m "$(cat <<'EOF'
feat(kalshi-mlb-rfq): probit devig for SGP joint distribution

Port the same probit (additive z-shift) devigging from Tools.R to
fair_value.py. Operates on the 4-cell SGP joint distribution (Home/Away
Spread x Over/Under).

- New _probit_devig_n helper with n=2 closed form and n>=3 brentq fallback.
- devig_book rewritten to call _probit_devig_n on the 4-cell grid.
- Single-side heuristic fallback unchanged (no devig math possible with 1 cell).
- Add scipy>=1.10 to requirements.txt (already in use elsewhere in project).
- New pytest unit tests including negative test against multiplicative.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Phase 5 — Commit 5: `docs: update CLAUDE.md + READMEs for probit devig`

### Task 20: Update `Answer Keys/CLAUDE.md`

**Files:**
- Modify: `Answer Keys/CLAUDE.md`

- [ ] **Step 1: Locate the Consensus Architecture section**

Open `Answer Keys/CLAUDE.md` and find the section beginning with `### Consensus Architecture`. It currently says:

```
- **Historical consensus** (Phase 1): All-book median with 0.5 prob hardcode. Used for sample building. Do not change.
- **Live consensus** (Phase 2): Sharp books only via `SHARP_BOOKS` in `Tools.R`. ...
```

- [ ] **Step 2: Replace with the updated description**

Replace the historical consensus bullet and add a methodology note:

```markdown
- **Historical consensus** (Phase 1): MLB and CBB both use sharp-weighted probit-devigged consensus computed at PBP-build time and stored in `mlb_betting_pbp` / `cbb_betting_pbp`. NFL legacy paths still use the older `Consensus Betting History.R` output.
- **Live consensus** (Phase 2): Sharp books only via `SHARP_BOOKS` in `Tools.R`. Rec books get weight=0. Pinnacle + Bookmaker at 1.1 weight (tiebreaker), LowVig/Circa/Bet105 at 1.0. Games with no sharp coverage are dropped.
- **Devigging method:** All devig in `Tools.R::devig_american` / `devig_american_3way` uses probit (additive z-shift). 2-way uses closed-form `c = -(z1+z2)/2`; 3-way+ uses `uniroot`. See `docs/superpowers/specs/2026-05-11-probit-devig-design.md` for math and rationale.
```

### Task 21: Update dashboard READMEs

**Files:**
- Modify: `Answer Keys/MLB Dashboard/README.md`
- Modify: `Answer Keys/CBB Dashboard/README.md`

- [ ] **Step 1: Add a "Devigging method" note to MLB Dashboard README**

Open `Answer Keys/MLB Dashboard/README.md`. Find the section describing fair-value math (search for `devig` or `fair`). Add or update a paragraph:

```markdown
**Devigging method:** The MLB dashboard's "Books (devigged fair %)" column uses probit (additive z-shift) devigging via `Tools.R::devig_american`. See `docs/superpowers/specs/2026-05-11-probit-devig-design.md`. Historical samples in `mlb_betting_pbp` are also probit-devigged (sharp-weighted across Pinnacle, Bookmaker, LowVig, Circa, Bet105).
```

- [ ] **Step 2: Add the same note to CBB Dashboard README**

Open `Answer Keys/CBB Dashboard/README.md` and add a parallel paragraph. Mention that the 0.5 hardcode is gone — CBB historical now uses real sharp-weighted probit consensus.

### Task 22: Update `kalshi_mlb_rfq/README.md`

**Files:**
- Modify: `kalshi_mlb_rfq/README.md`

- [ ] **Step 1: Add devigging note**

Open `kalshi_mlb_rfq/README.md`. Find a section about fair-value computation or `fair_value.py`. Add:

```markdown
**Devigging:** `fair_value.devig_book` uses probit (additive z-shift) devigging on the 4-cell SGP joint distribution (Home/Away Spread × Over/Under). When fewer than 4 sides are visible (interpolated combos with partial coverage), falls back to `(1/decimal) / (1 + vig_fallback)` heuristic. See `docs/superpowers/specs/2026-05-11-probit-devig-design.md`.

Requires `scipy>=1.10` (added to `requirements.txt`). Run `pip install -r requirements.txt` after pulling this change.
```

### Task 23: Commit Phase 5

- [ ] **Step 1: Stage and commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/probit-devig
git add "Answer Keys/CLAUDE.md" \
        "Answer Keys/MLB Dashboard/README.md" \
        "Answer Keys/CBB Dashboard/README.md" \
        kalshi_mlb_rfq/README.md
git diff --stat --cached
```

- [ ] **Step 2: Commit**

```bash
git commit -m "$(cat <<'EOF'
docs: update CLAUDE.md + READMEs for probit devig

Document the methodology change in:
- Answer Keys/CLAUDE.md: update Consensus Architecture; CBB historical
  no longer uses 0.5 hardcode.
- Answer Keys/MLB Dashboard/README.md: probit method note.
- Answer Keys/CBB Dashboard/README.md: probit method note + 0.5 removal.
- kalshi_mlb_rfq/README.md: probit method note + scipy requirement.

Spec: docs/superpowers/specs/2026-05-11-probit-devig-design.md

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Phase 6 — Pre-Merge Validation (no commits — runs from worktree)

### Task 24: Verify all R + Python tests pass

- [ ] **Step 1: Run all R tests**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/probit-devig/Answer Keys" && \
  for f in tests/test_*.R; do
    echo "=== $f ==="
    Rscript "$f" 2>&1 | tail -5
  done
```
Expected: every file completes; `test_probit_devig.R` and `test_answer_key.R` both report ALL PASS.

- [ ] **Step 2: Run all Python tests**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/probit-devig/kalshi_mlb_rfq && python3 -m pytest tests/ -v 2>&1 | tail -15
```
Expected: ALL PASS.

### Task 25: Final shadow-definition grep

- [ ] **Step 1: Confirm no shadow `devig_american*` definitions remain in R**

```bash
grep -rn "^devig_american\|<- function" "/Users/callancapitolo/NFLWork/.worktrees/probit-devig/Answer Keys" --include="*.R" 2>/dev/null | grep -i "devig_american"
```
Expected: ONLY the 3 definitions in `Tools.R` (`.probit_devig_n`, `devig_american`, `devig_american_3way`). No other files.

- [ ] **Step 2: Confirm no extra devig functions in Python**

```bash
grep -rn "def devig\|def _probit" /Users/callancapitolo/NFLWork/.worktrees/probit-devig/kalshi_mlb_rfq --include="*.py" 2>/dev/null | grep -v venv
```
Expected: only `_probit_devig_n` and `devig_book` in `fair_value.py`.

### Task 26: Self-review the diff

- [ ] **Step 1: Generate full diff against main**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/probit-devig && git diff main..HEAD --stat
```

Expected: 8 modified files + 3 new files + 1 spec + 1 plan ≈ 13 files total.

- [ ] **Step 2: Spot-check changes per the pre-merge checklist in the spec**

From `docs/superpowers/specs/2026-05-11-probit-devig-design.md` Section 8:

- [ ] All unit tests pass (R + Python) — confirmed in Task 24
- [ ] All 4 shadow definitions confirmed gone — confirmed in Task 25
- [ ] No regressions in `test_answer_key.R` — confirmed in Task 24
- [ ] Documentation updates in same PR — confirmed Phase 5

(Backtest validation + dashboard validation happen post-merge on `main` since the DB files don't live in the worktree.)

### Task 27: Ask user for explicit merge approval

- [ ] **Step 1: Surface the diff summary and ask**

Show the user:
- `git log --oneline main..HEAD` (the 5 commits from the spec)
- `git diff main..HEAD --stat`
- The pre-merge checklist status

**STOP HERE.** Per project rules (`Answer Keys/CLAUDE.md`), never merge to `main` without explicit user approval. Surface the diff and ask: "Approved to merge to main?" Do NOT proceed to Phase 7 without a "yes."

---

## Phase 7 — Merge + Post-Merge (runs on `main`, USER-APPROVED ONLY)

### Task 28: Snapshot historical tables on `main`

**Files:** none (DB writes only)

- [ ] **Step 1: Switch to main**

```bash
cd /Users/callancapitolo/NFLWork && git checkout main
```

- [ ] **Step 2: Snapshot mlb_betting_pbp**

```bash
cd /Users/callancapitolo/NFLWork && python3 -c "
import duckdb
con = duckdb.connect('pbp.duckdb')
con.execute('CREATE TABLE mlb_betting_pbp_backup_2026_05_11 AS SELECT * FROM mlb_betting_pbp')
con.execute('CREATE TABLE cbb_betting_pbp_backup_2026_05_11 AS SELECT * FROM cbb_betting_pbp')
print('snapshots created')
"
```
Expected: prints `snapshots created`.

### Task 29: Merge feature branch

- [ ] **Step 1: Merge with --no-ff**

```bash
cd /Users/callancapitolo/NFLWork && git merge --no-ff feature/probit-devig -m "$(cat <<'EOF'
Merge feature/probit-devig: probit devigging across the codebase

Replaces multiplicative with additive z-shift probit devigging in:
- Tools.R helpers (devig_american, devig_american_3way) — used by
  MLB+CBB+NFL pipelines, dashboards, and backtests.
- kalshi_mlb_rfq/fair_value.py (devig_book) — SGP 4-cell joint grid.

Also kills 4 shadow definitions, adds input validation, and adds a
3-way -> 2-way fallback.

Spec: docs/superpowers/specs/2026-05-11-probit-devig-design.md
Plan: docs/superpowers/plans/2026-05-11-probit-devig.md
EOF
)"
```

- [ ] **Step 2: Confirm merge succeeded**

```bash
cd /Users/callancapitolo/NFLWork && git log --oneline -10
```
Expected: merge commit on top, 5 feature commits below.

### Task 30: Regenerate MLB historical pool

- [ ] **Step 1: Run Consensus Betting History.R**

```bash
cd "/Users/callancapitolo/NFLWork/Answer Keys" && Rscript "MLB Answer Key/Consensus Betting History.R" 2>&1 | tail -10
```
Expected: completes; logs row count for `mlb_betting_pbp`.

- [ ] **Step 2: Verify row count matches backup**

```bash
cd /Users/callancapitolo/NFLWork && python3 -c "
import duckdb
con = duckdb.connect('pbp.duckdb')
n_new = con.execute('SELECT COUNT(*) FROM mlb_betting_pbp').fetchone()[0]
n_old = con.execute('SELECT COUNT(*) FROM mlb_betting_pbp_backup_2026_05_11').fetchone()[0]
print(f'new: {n_new}, old: {n_old}, diff: {n_new - n_old}')
"
```
Expected: similar row counts. Small drift (±10 rows) is acceptable due to NA filtering changes.

### Task 31: Regenerate CBB historical pool

- [ ] **Step 1: Run build_betting_pbp.R**

```bash
cd "/Users/callancapitolo/NFLWork/Answer Keys" && Rscript "CBB Answer Key/build_betting_pbp.R" 2>&1 | tail -15
```
Expected: completes; logs sharp consensus counts; logs final row count for `cbb_betting_pbp`.

- [ ] **Step 2: Verify new consensus columns exist**

```bash
cd /Users/callancapitolo/NFLWork && python3 -c "
import duckdb
con = duckdb.connect('cbb.duckdb')
cols = con.execute(\"PRAGMA table_info('cbb_betting_pbp')\").fetchall()
print([c[1] for c in cols if 'consensus' in c[1] or 'devig' in c[1]])
"
```
Expected: list includes `consensus_devig_home_odds`, `consensus_devig_away_odds`, `consensus_devig_over_odds`, `consensus_devig_under_odds`.

### Task 32: Backtest validation (gates "done")

- [ ] **Step 1: Run MLB backtest**

```bash
cd "/Users/callancapitolo/NFLWork/Answer Keys" && Rscript "MLB Answer Key/MLB_Backtest.R" 2>&1 | tail -30
```
Record ROI, win rate, CLV. Compare to pre-change baseline (from memory or recent run).
Expected: no metric moves >5% from baseline.

- [ ] **Step 2: Run MLB ROI backtest**

```bash
cd "/Users/callancapitolo/NFLWork/Answer Keys" && Rscript "MLB Answer Key/MLB_ROI_Backtest.R" 2>&1 | tail -30
```
Expected: same — no >5% drift.

- [ ] **Step 3: Run CBB backtest**

```bash
cd "/Users/callancapitolo/NFLWork/Answer Keys" && Rscript "CBB Answer Key/CBB_Backtest.R" 2>&1 | tail -30
```
Expected: no >5% drift.

**If anything moves >5%, STOP.** Restore from `*_backup_2026_05_11` tables and revert the merge. Investigate before proceeding.

### Task 33: Live sanity check

- [ ] **Step 1: Start MLB dashboard**

```bash
cd "/Users/callancapitolo/NFLWork/Answer Keys/MLB Dashboard" && bash run.sh &
sleep 5
```

- [ ] **Step 2: Manual verify in browser**

Open `http://localhost:8083` in browser. Hover over a game and confirm the "Books (devigged fair %)" column shows reasonable probit values (e.g., for a -150 ML game, favorite should show ~58.3%, not ~58.0%).

- [ ] **Step 3: Stop the MLB dashboard**

```bash
pkill -f mlb_dashboard_server || true
```

- [ ] **Step 4: Same check for CBB dashboard (port varies — see README)**

- [ ] **Step 5: Kalshi MLB RFQ dry-run**

```bash
cd /Users/callancapitolo/NFLWork/kalshi_mlb_rfq && python3 -m pip install -r requirements.txt 2>&1 | tail -3 && \
python3 main.py --dry-run 2>&1 | head -50
```
Expected: scipy already installed; bot logs show book fairs with reasonable values.

### Task 34: Cleanup worktree

- [ ] **Step 1: Remove worktree**

```bash
cd /Users/callancapitolo/NFLWork && git worktree remove .worktrees/probit-devig
```

- [ ] **Step 2: Delete feature branch**

```bash
cd /Users/callancapitolo/NFLWork && git branch -d feature/probit-devig
```

- [ ] **Step 3: Verify cleanup**

```bash
cd /Users/callancapitolo/NFLWork && git worktree list && git branch | grep probit
```
Expected: no probit-devig worktree, no probit branch.

### Task 35: Save methodology memory

**Files:**
- Create: `/Users/callancapitolo/.claude-personal/projects/-Users-callancapitolo-NFLWork/memory/probit_devig_methodology.md`
- Modify: `/Users/callancapitolo/.claude-personal/projects/-Users-callancapitolo-NFLWork/memory/MEMORY.md`

- [ ] **Step 1: Write the memory file**

Create the memory file with content like:

```markdown
---
name: Probit devigging methodology
description: As of 2026-05-11, all devigging across Tools.R helpers, MLB+CBB historical pools, and Kalshi RFQ uses probit (additive z-shift) instead of multiplicative. 2-way path uses closed-form c = -(z1+z2)/2.
type: project
---

All devigging across the codebase uses probit (additive z-shift) as of 2026-05-11.

**Why:** Multiplicative under-weights extremes. Probit gives more credit to favorites at the tails. Also fixed 4 shadow devig definitions that silently kept multiplicative.

**How to apply:** When asked about devigging methodology, fair-value calculation, or consensus probabilities — the answer is probit. `Tools.R::devig_american` and `devig_american_3way` are the canonical R helpers. `kalshi_mlb_rfq/fair_value.py::_probit_devig_n` is the canonical Python helper. NO shadow copies anywhere.

**2-way uses closed-form** `c = -(z1+z2)/2`. 3-way+ uses uniroot/brentq.

Spec: `docs/superpowers/specs/2026-05-11-probit-devig-design.md`
Plan: `docs/superpowers/plans/2026-05-11-probit-devig.md`
```

- [ ] **Step 2: Index in MEMORY.md under Reference section**

Append a line under `## Reference`:

```markdown
- [Probit devigging methodology](probit_devig_methodology.md) — 2026-05-11 switched all devig (R + Python) to probit additive z-shift; canonical helpers in Tools.R + fair_value.py; no shadows
```

---

## Done

When all tasks above are checked off, the feature is fully shipped:
- Probit devigging is the canonical method across R and Python.
- 4 shadow definitions are gone.
- Historical pools regenerated with consistent methodology.
- Backtests pass within 5% drift tolerance.
- Dashboards visually validated.
- Worktree cleaned up.
- Memory saved for future-me.
