# MLB Parlay Tab — Books-Agree Count Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an "Agree k/n" summary pill to the end of the existing Books strip on the MLB dashboard's Parlay tab. The pill shows how many of the five reference fair-prob voices (DK, FD, ProphetX, Novig, model) lie within ±1pp of the median of the five — a caution signal that the user reads alongside Edge %.

**Architecture:** All work happens in two existing files. `Answer Keys/books_strip.R` (the pure HTML helper) gets a new `compute_k_within()` function and an extended `render_books_strip()` signature. `Answer Keys/MLB Dashboard/mlb_dashboard.R` gets a tiny edit in the Books cell renderer to compute k/n and pass them to the renderer, plus one new CSS rule for the `.pill.agree` class. No DB schema changes, no Python changes, no new files.

**Tech Stack:** R, testthat, reactable (HTML cell renderers).

**Spec:** [`docs/superpowers/specs/2026-05-21-parlay-books-agree-count-design.md`](../specs/2026-05-21-parlay-books-agree-count-design.md)

---

## File Structure

**Modified files (3):**
- `Answer Keys/books_strip.R` — adds `compute_k_within()` and extends `render_books_strip()` to accept optional `k_agree` / `n_agree` arguments. Both functions stay pure (no DB, no Shiny refs).
- `Answer Keys/MLB Dashboard/mlb_dashboard.R` — two edits: (1) add `.pill.agree` CSS rule alongside the existing `.pill.cons` rule; (2) update the Books column cell function inside `create_parlays_table()` to compute k/n and pass them into `render_books_strip()`.
- `Answer Keys/tests/test_books_strip.R` — adds testthat cases for `compute_k_within()` and for the new pill output in `render_books_strip()`.
- `Answer Keys/CLAUDE.md` — adds a documentation paragraph describing the new pill.

**No new files. No DB schema changes. No Python changes.**

**Why edit `books_strip.R` rather than `mlb_dashboard.R` for the helper (deviation from spec):** The spec body says "add `compute_k_within()` to `mlb_dashboard.R` near `apply_combo_residuals()`". On inspection, `render_books_strip()` already lives in a separate, pure, testable file (`Answer Keys/books_strip.R`) with an existing test file. Putting the new helper in the same place keeps related code together AND lets us write fast unit tests via the existing test harness without sourcing the entire dashboard (which has top-level side effects). The user has accepted similar minor structural deviations on past PRs to enable testing; flagging here for transparency.

---

## Task 1: Add `compute_k_within()` helper (TDD)

**Files:**
- Test: `Answer Keys/tests/test_books_strip.R`
- Implementation: `Answer Keys/books_strip.R`

### Step 1: Write the failing tests

Append the following block to `Answer Keys/tests/test_books_strip.R` (after the existing `render_books_strip` tests at the bottom of the file):

```r
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
```

### Step 2: Run the tests to verify they fail

Run from the worktree root:

```bash
cd "Answer Keys/tests" && Rscript -e "testthat::test_file('test_books_strip.R')"
```

Expected: the 6 new `compute_k_within` tests fail with an error like `could not find function "compute_k_within"`. The 5 existing `render_books_strip` tests still pass.

### Step 3: Implement `compute_k_within()` in `books_strip.R`

Append the following function to `Answer Keys/books_strip.R` (after the existing `render_books_strip()` definition):

```r

# compute_k_within(): counts how many entries in `probs` fall within +/- `band`
# of the median of the non-NA entries. NAs are excluded from both the median
# and the count. Returns list(k, n) where n is the count of non-NA entries.
#
# When fewer than 2 non-NA entries exist, k is NA (no meaningful consensus
# check is possible with a single voice). The caller decides whether to
# render anything in that case.
#
# `probs` is on the [0, 1] scale (probabilities, not percentages). `band`
# is on the same scale, so 0.01 == 1 percentage point.
compute_k_within <- function(probs, band = 0.01) {
  probs <- probs[!is.na(probs)]
  n <- length(probs)
  if (n < 2) return(list(k = NA_integer_, n = n))
  med <- median(probs)
  list(k = sum(abs(probs - med) <= band), n = n)
}
```

### Step 4: Run the tests to verify they pass

```bash
cd "Answer Keys/tests" && Rscript -e "testthat::test_file('test_books_strip.R')"
```

Expected: all 11 tests pass (5 existing + 6 new). Output ends with `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 11 ]` or similar.

### Step 5: Commit

```bash
git add "Answer Keys/books_strip.R" "Answer Keys/tests/test_books_strip.R"
git commit -m "feat(parlay-tab): add compute_k_within() consensus helper

Counts how many of a vector of fair-probability voices fall within
+/- band of the median. NA-tolerant: drops NAs from both numerator
and denominator so a quiet book doesn't penalize the row.

Used by the next commit to surface 'Agree k/n' on the parlay tab."
```

---

## Task 2: Extend `render_books_strip()` to render the Agree pill (TDD)

**Files:**
- Test: `Answer Keys/tests/test_books_strip.R`
- Implementation: `Answer Keys/books_strip.R`

### Step 1: Write the failing tests

Append the following block to `Answer Keys/tests/test_books_strip.R` (after the `compute_k_within` tests):

```r
# ---------------------------------------------------------------------------
# render_books_strip(): optional Agree k/n suffix pill.
# ---------------------------------------------------------------------------

test_that("render_books_strip omits Agree pill when k_agree is not provided", {
  out <- render_books_strip(0.40, 0.41, 0.42, 0.40, 0.55, 0.41)
  expect_false(grepl("Agree", out, fixed = TRUE))
})

test_that("render_books_strip renders Agree pill when k/n provided", {
  out <- render_books_strip(0.40, 0.41, 0.42, 0.40, 0.55, 0.41,
                            k_agree = 4L, n_agree = 5L)
  expect_match(out, '<span class="pill agree">Agree 4/5</span>', fixed = TRUE)
})

test_that("render_books_strip omits Agree pill when n_agree < 2", {
  out <- render_books_strip(0.40, NA_real_, NA_real_, NA_real_, NA_real_, 0.40,
                            k_agree = NA_integer_, n_agree = 1L)
  expect_false(grepl("Agree", out, fixed = TRUE))
})

test_that("render_books_strip renders partial-denominator Agree pill (k/n=3/4)", {
  # Four books quoted, one missing — denominator drops to 4.
  out <- render_books_strip(0.40, 0.41, 0.42, NA_real_, 0.55, 0.41,
                            k_agree = 3L, n_agree = 4L)
  expect_match(out, '<span class="pill agree">Agree 3/4</span>', fixed = TRUE)
})

test_that("render_books_strip Agree pill appears AFTER the Cons pill", {
  out <- render_books_strip(0.40, 0.41, 0.42, 0.40, 0.55, 0.41,
                            k_agree = 4L, n_agree = 5L)
  cons_pos <- regexpr('Cons 41.0', out, fixed = TRUE)
  agree_pos <- regexpr('Agree 4/5', out, fixed = TRUE)
  expect_true(cons_pos > 0 && agree_pos > 0)
  expect_true(agree_pos > cons_pos)
})
```

### Step 2: Run the tests to verify they fail

```bash
cd "Answer Keys/tests" && Rscript -e "testthat::test_file('test_books_strip.R')"
```

Expected: the 5 new tests fail because `render_books_strip()` doesn't accept `k_agree` / `n_agree` arguments and doesn't emit any `Agree` pill. The other tests (11 from Task 1 + existing) still pass.

### Step 3: Extend `render_books_strip()` in `books_strip.R`

Replace the existing `render_books_strip()` definition in `Answer Keys/books_strip.R` with this version (the function body grows to accept `k_agree` / `n_agree` and append a final pill when both are usable):

```r
render_books_strip <- function(model, dk, fd, px, nv, cons,
                                k_agree = NULL, n_agree = NULL) {
  fmt_pill <- function(label, prob, css_class) {
    if (is.na(prob)) {
      sprintf('<span class="pill %s dim">%s &mdash;</span>', css_class, label)
    } else {
      sprintf('<span class="pill %s">%s %.1f</span>', css_class, label, prob * 100)
    }
  }

  # Optional Agree k/n suffix pill — see compute_k_within() above.
  # Omit entirely when n_agree < 2 (no meaningful consensus possible) or
  # when k_agree is NA. The caller is responsible for computing k/n via
  # compute_k_within() over whichever voices are in scope.
  agree_pill <- if (!is.null(k_agree) && !is.null(n_agree) &&
                    !is.na(k_agree) && n_agree >= 2L) {
    sprintf('<span class="pill agree">Agree %d/%d</span>',
            as.integer(k_agree), as.integer(n_agree))
  } else {
    ""
  }

  paste0(
    '<span class="books-strip">',
    fmt_pill("M",    model, "model"),
    fmt_pill("DK",   dk,    "book"),
    fmt_pill("FD",   fd,    "book"),
    fmt_pill("PX",   px,    "book"),
    fmt_pill("NV",   nv,    "book"),
    fmt_pill("Cons", cons,  "cons"),
    agree_pill,
    '</span>'
  )
}
```

### Step 4: Run the tests to verify they pass

```bash
cd "Answer Keys/tests" && Rscript -e "testthat::test_file('test_books_strip.R')"
```

Expected: all tests pass (5 existing + 6 compute_k_within + 5 new render tests = 16 total). Output: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 16 ]`.

### Step 5: Commit

```bash
git add "Answer Keys/books_strip.R" "Answer Keys/tests/test_books_strip.R"
git commit -m "feat(parlay-tab): render_books_strip emits optional Agree k/n pill

Adds two optional args (k_agree, n_agree). When both are provided and
n_agree>=2, appends a 'pill agree' span at the end of the books strip
showing 'Agree k/n'. Existing callers that don't pass the new args
behave identically — no behavior change for any current consumer."
```

---

## Task 3: Add `.pill.agree` CSS to `mlb_dashboard.R`

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` (insert new CSS rule after the existing `.pill.dim` block)

### Step 1: Locate the insertion point

The dashboard inlines all parlay-tab CSS as a long string. The pill rules currently live around lines 2435–2470. Find the `.pill.dim` rule:

```css
.pill.dim {
  opacity: 0.4;
}
```

The new rule goes immediately after it, before the `.pill.error` rule.

### Step 2: Insert the new CSS rule

Open `Answer Keys/MLB Dashboard/mlb_dashboard.R`. After the closing `}` of `.pill.dim`, insert:

```css
        .pill.agree {
          background: #1f6feb22;
          color: #8ab4f8;
          font-weight: 700;
          margin-left: 2px;
        }
```

Match the existing 8-space indentation used by the surrounding rules. The CSS lives inside an `HTML(...)` string in the R file, so the indentation is part of the string literal — preserve it exactly.

### Step 3: Smoke-check the syntax

Make sure no R syntax was broken by the edit:

```bash
Rscript -e 'parse("Answer Keys/MLB Dashboard/mlb_dashboard.R")' >/dev/null && echo OK
```

Expected: prints `OK`. If it errors, the CSS block was inserted outside the surrounding R string — fix and re-run.

### Step 4: Commit

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "style(parlay-tab): add .pill.agree CSS for consensus-count pill

Blue-tinted variant of the existing book pills so the eye can find it
while still reading as part of the books strip. Matches the B3 mockup
from the design spec."
```

---

## Task 4: Wire `compute_k_within()` into the Books cell renderer

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` (the Books column cell function inside `create_parlays_table()`, around line 544–561)

### Step 1: Locate the current cell function

The existing Books column definition reads:

```r
      books_strip = colDef(
        name = "Books (devigged fair %)",
        minWidth = 320,
        html = TRUE,
        sortable = FALSE,
        class = "cell-books",
        cell = function(value, index) {
          row <- table_data[index, ]
          render_books_strip(
            model = row$model_prob_raw,
            dk    = row$dk_fair_prob,
            fd    = row$fd_fair_prob,
            px    = row$px_fair_prob,
            nv    = row$nv_fair_prob,
            cons  = row$blended_prob_raw
          )
        }
      ),
```

### Step 2: Edit the cell function to compute and pass k/n

Replace the cell function body with this version. Note that `Cons` is intentionally excluded from the count — it's a blended derivative of model + books, not an independent voice; including it would double-count the cluster:

```r
      books_strip = colDef(
        name = "Books (devigged fair %)",
        minWidth = 320,
        html = TRUE,
        sortable = FALSE,
        class = "cell-books",
        cell = function(value, index) {
          row <- table_data[index, ]
          # Agree k/n is computed over the 5 independent voices (DK / FD / PX /
          # NV / model). Cons is excluded because it's a blended derivative of
          # model+books. NA voices reduce the denominator rather than counting
          # as "disagree" — see books_strip.R::compute_k_within.
          agree <- compute_k_within(c(
            row$dk_fair_prob, row$fd_fair_prob,
            row$px_fair_prob, row$nv_fair_prob,
            row$model_prob_raw
          ))
          render_books_strip(
            model = row$model_prob_raw,
            dk    = row$dk_fair_prob,
            fd    = row$fd_fair_prob,
            px    = row$px_fair_prob,
            nv    = row$nv_fair_prob,
            cons  = row$blended_prob_raw,
            k_agree = agree$k,
            n_agree = agree$n
          )
        }
      ),
```

### Step 3: Smoke-check the syntax

```bash
Rscript -e 'parse("Answer Keys/MLB Dashboard/mlb_dashboard.R")' >/dev/null && echo OK
```

Expected: prints `OK`.

### Step 4: Run the full books-strip test suite one more time to confirm nothing regressed

```bash
cd "Answer Keys/tests" && Rscript -e "testthat::test_file('test_books_strip.R')"
```

Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 16 ]`.

### Step 5: Commit

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "feat(parlay-tab): wire compute_k_within into Books cell renderer

Each parlay row now passes its 5 independent voices (DK / FD / PX /
NV / model) through compute_k_within() and feeds k/n into
render_books_strip() so an 'Agree k/5' pill renders at the end of
the existing books pill row. Cons is excluded from the count to
avoid double-counting the cluster."
```

---

## Task 5: Live dashboard smoke verification

**Goal:** Confirm the new pill renders correctly on real parlay rows and that the existing UI hasn't regressed.

### Step 1: Render an HTML sample directly via Rscript

This catches any obvious string/escape bug without launching the full dashboard. Run from the worktree root (single line — bash backslash-newlines inside single quotes don't work in R):

```bash
Rscript -e 'source("Answer Keys/books_strip.R"); a <- compute_k_within(c(0.40,0.41,0.42,0.40,0.55)); cat(render_books_strip(0.40, 0.41, 0.42, 0.40, 0.55, 0.41, k_agree = a$k, n_agree = a$n), "\n")'
```

Expected output (one line, no errors):

```html
<span class="books-strip"><span class="pill model">M 40.0</span><span class="pill book">DK 41.0</span><span class="pill book">FD 42.0</span><span class="pill book">PX 40.0</span><span class="pill book">NV 55.0</span><span class="pill cons">Cons 41.0</span><span class="pill agree">Agree 4/5</span></span>
```

Verify the `Agree 4/5` pill is the last span before `</span>` closing `books-strip`.

### Step 2: Try a "no-pill" case (n=1)

```bash
Rscript -e 'source("Answer Keys/books_strip.R"); a <- compute_k_within(c(0.40, NA, NA, NA, NA)); cat(render_books_strip(0.40, NA_real_, NA_real_, NA_real_, NA_real_, 0.40, k_agree = a$k, n_agree = a$n), "\n")'
```

Expected: no `Agree` pill in the output (just the dim DK/FD/PX/NV pills and the Cons pill).

### Step 3: Start the dashboard locally and view the parlay tab

From the worktree root, run whatever the local convention is for starting the MLB dashboard (e.g., `cd "Answer Keys/MLB Dashboard" && Rscript mlb_dashboard.R`, or the project's normal launch command). Then visit `http://localhost:8083` and click the Parlay tab.

**Verification checklist:**
- [ ] At least one parlay row shows a blue `Agree k/5` pill at the end of the Books cell.
- [ ] The pill is visually distinct from the existing M / DK / FD / PX / NV / Cons pills (blue tint vs. their colors).
- [ ] When a book is missing (dim pill `—`), the Agree pill shows a reduced denominator like `Agree 3/4`.
- [ ] No JS console errors, no R errors in the launching shell.
- [ ] Existing columns (Game / Legs / Fair / WZ / Edge / Size / Action) still render correctly.

If anything is off, the most likely culprits are (1) CSS not applying (refresh browser, hard reload), (2) the column's `minWidth = 320` is now too tight — bump to `360` if pills wrap onto two lines.

### Step 4: Commit if any tweaks were needed

If you bumped `minWidth` or made any other tiny adjustments in this step, commit them:

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "chore(parlay-tab): tweak books column width for Agree pill"
```

Otherwise skip this step.

---

## Task 6: Document the new pill in `Answer Keys/CLAUDE.md`

**Files:**
- Modify: `Answer Keys/CLAUDE.md`

### Step 1: Find the insertion point

Open `Answer Keys/CLAUDE.md`. Scroll to the bottom of the existing "MLB Dashboard — Odds screen + WZ single-bet placer" section (it ends around the "Editable Risk + WZ verified-quote (PR C, 2026-05)" subsection).

### Step 2: Append a new subsection

After the last paragraph of that section, append the following Markdown:

````markdown

## MLB Dashboard — Parlay tab "Agree" pill

The parlay-tab Books cell ends with an `Agree k/n` pill (e.g. `Agree 4/5`)
that summarizes how many of the five reference voices — DK, FD, ProphetX,
Novig, and the model — sit within ±1pp of the median of the five. It is a
display-only caution signal: high `k` means the market is in price-discovery
consensus; low `k` means the books disagree and the edge against any single
book is less trustworthy.

- **Implementation:** `Answer Keys/books_strip.R::compute_k_within()` does the
  count; `render_books_strip()` accepts optional `k_agree` / `n_agree` args
  and appends the suffix pill. The caller (in `mlb_dashboard.R`'s
  `create_parlays_table()`, around line 552) builds the 5-element voice
  vector and passes the result through.
- **Cons is intentionally excluded** from the count because it's a blended
  derivative of model + books, not an independent voice.
- **NA voices reduce the denominator** rather than counting as "disagree";
  a row with only 4 quoting voices shows `Agree k/4`.
- **No filter / no color / no auto-bet gate** in v1. Layer those on once
  live values show what threshold actually feels right.
- **Known blind spot:** when the model is the lone dissenter from a tight
  book cluster, the metric reads identically to "a book is the lone dissenter"
  — see scenario H in the design spec. Easy mitigation later: a separate
  `model_in_cluster` boolean.
- **Tests:** `Answer Keys/tests/test_books_strip.R` (both `compute_k_within`
  and the `render_books_strip` extension).
````

### Step 3: Commit

```bash
git add "Answer Keys/CLAUDE.md"
git commit -m "docs(parlay-tab): document the Books Agree k/n pill"
```

---

## Task 7: Hand off for merge approval

### Step 1: Review the full diff against `main`

```bash
git diff main..HEAD --stat
git diff main..HEAD
```

Sanity checks before reporting ready-to-merge:

- [ ] Only 4 files changed: `Answer Keys/books_strip.R`, `Answer Keys/tests/test_books_strip.R`, `Answer Keys/MLB Dashboard/mlb_dashboard.R`, `Answer Keys/CLAUDE.md`.
- [ ] No DB files, no temp files, no scratch outputs.
- [ ] All commits have meaningful messages.
- [ ] No `TODO`/`FIXME` strings left in the changes.

### Step 2: Run all tests one final time

```bash
cd "Answer Keys/tests" && Rscript -e "testthat::test_file('test_books_strip.R')"
```

Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 16 ]`.

### Step 3: Report status and wait for explicit merge approval

Report back to the user:

> "Implementation complete on worktree `parlay-agree-count`. 4 files changed, 16 tests passing, dashboard smoke-verified locally. Awaiting your approval to merge into `main` and clean up the worktree."

Per project policy (`/Users/callancapitolo/NFLWork/CLAUDE.md`), **do NOT merge to `main` without explicit user approval**, even though tests pass.

### Step 4 (only after user approves merge)

```bash
git checkout main
git merge --no-ff worktree-parlay-agree-count
git worktree remove "/Users/callancapitolo/NFLWork/.claude/worktrees/parlay-agree-count"
git branch -d worktree-parlay-agree-count
```

Report final status.
