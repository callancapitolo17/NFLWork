# MLB Parlay Dashboard — Books Strip Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Surface our model's joint probability and replace the four per-book probability columns with a single inline pill row (`M · DK · FD · PX · NV · Cons`) so the parlay table reads cleanly on phone and laptop without hiding any data.

**Architecture:** Extract the pill-row HTML rendering into a pure R helper (`render_books_strip()`) with `testthat` unit coverage. Wire it into the parlay reactable as a single `colDef` `cell` function. Hide the four obsolete per-book `colDef`s and the standalone Time column. Fold the game time into the Game column via a JS cell renderer reading `cellInfo.row.game_time`. Add inline CSS to the dashboard's existing `<style>` block.

**Tech Stack:** R, `tidyverse`, `reactable`, `htmltools`, `testthat`. Single dashboard file (`mlb_dashboard.R`); no schema, server, R-pricer, or Python changes.

**Spec:** [`docs/superpowers/specs/2026-04-25-mlb-parlay-dashboard-redesign-design.md`](../specs/2026-04-25-mlb-parlay-dashboard-redesign-design.md)

---

## File Structure

| Action | Path | Responsibility |
|--------|------|----------------|
| Create | `Answer Keys/books_strip.R` | Pure helper `render_books_strip()` returning the pill-row HTML string |
| Create | `Answer Keys/tests/test_books_strip.R` | `testthat` tests for the helper |
| Modify | `Answer Keys/MLB Dashboard/mlb_dashboard.R` | Source helper, swap per-book `colDef`s for one Books `colDef`, hide model/time, add CSS to inline `<style>` block, fold time into Game cell |
| Modify | `Answer Keys/MLB Dashboard/README.md` | One-paragraph note under Features describing the strip + Cons |

The helper lives in `Answer Keys/` (sibling to `parse_legs.R`) so the existing test pattern (`tests/test_*.R` sources sibling files via `source("../helper.R")`) continues to work.

---

## Version Control

- **Branch:** `feature/parlay-books-strip` (created from `main`)
- **Worktree:** `.worktrees/parlay-books-strip` — created in Task 0, removed in Task 6
- **Commits:**
  1. `feat(mlb-dash): pure render_books_strip helper + tests` (Task 1)
  2. `feat(mlb-dash): swap per-book columns for Books pill row` (Task 2)
  3. `style(mlb-dash): pill styling + phone responsive Corr hide` (Task 3)
  4. `docs(mlb-dash): README note on books strip + Cons` (Task 5)
- **Merge to main:** Task 6, after pre-merge review (`git diff main..HEAD`) and explicit user approval.

---

## Task 0: Worktree setup

**Files:** none modified — branch creation only.

- [ ] **Step 1: Confirm clean main**

```bash
cd /Users/callancapitolo/NFLWork
git status --short
git rev-parse --abbrev-ref HEAD
```

Expected: empty (or only untracked) status, `main` branch.

- [ ] **Step 2: Create worktree on a new branch**

```bash
cd /Users/callancapitolo/NFLWork
git worktree add -b feature/parlay-books-strip .worktrees/parlay-books-strip main
```

Expected: `Preparing worktree (new branch 'feature/parlay-books-strip')`

- [ ] **Step 3: Verify worktree branch**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-books-strip"
git branch --show-current
```

Expected: `feature/parlay-books-strip`

**All subsequent commands run from the worktree:** `/Users/callancapitolo/NFLWork/.worktrees/parlay-books-strip`. Never edit files in the main checkout while the worktree exists.

---

## Task 1: Pure helper + unit tests (TDD)

**Files:**
- Create: `Answer Keys/books_strip.R`
- Create: `Answer Keys/tests/test_books_strip.R`

- [ ] **Step 1: Write the failing tests**

Create `Answer Keys/tests/test_books_strip.R` with:

```r
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
  # Non-NA pills still render their values
  expect_match(out, '<span class="pill model">M 27.4</span>', fixed = TRUE)
  expect_match(out, '<span class="pill book">FD 27.1</span>', fixed = TRUE)
})

test_that("render_books_strip handles all NA except model", {
  out <- render_books_strip(0.30, NA_real_, NA_real_, NA_real_, NA_real_, 0.30)
  expect_match(out, '<span class="pill model">M 30.0</span>', fixed = TRUE)
  expect_match(out, '<span class="pill book dim">DK &mdash;</span>', fixed = TRUE)
  expect_match(out, '<span class="pill book dim">FD &mdash;</span>', fixed = TRUE)
  expect_match(out, '<span class="pill book dim">PX &mdash;</span>', fixed = TRUE)
  expect_match(out, '<span class="pill book dim">NV &mdash;</span>', fixed = TRUE)
})

test_that("render_books_strip rounds to one decimal", {
  out <- render_books_strip(0.27449, 0.26951, 0.27101, 0.27801, 0.26701, 0.27201)
  # 27.449 → 27.4, 27.801 → 27.8 (HALF_EVEN or HALF_UP both fine here — boundary cases excluded)
  expect_match(out, '>M 27.4</span>', fixed = TRUE)
  expect_match(out, '>DK 27.0</span>', fixed = TRUE)  # 26.951 rounds up to 27.0
  expect_match(out, '>PX 27.8</span>', fixed = TRUE)
})

test_that("render_books_strip output is wrapped in container", {
  out <- render_books_strip(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
  expect_true(startsWith(out, '<span class="books-strip">'))
  expect_true(endsWith(out, '</span>'))
})
```

- [ ] **Step 2: Run the tests and confirm they fail**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-books-strip/Answer Keys/tests"
Rscript -e 'library(testthat); test_file("test_books_strip.R")'
```

Expected: error sourcing `../books_strip.R` (`cannot open file`) — the helper doesn't exist yet.

- [ ] **Step 3: Write the minimal helper**

Create `Answer Keys/books_strip.R`:

```r
# Answer Keys/books_strip.R
#
# Pure HTML renderer for the parlay-tab "books strip" — a single inline pill row
# showing the model + four book devigged fair probabilities + blended consensus.
# Used by `MLB Dashboard/mlb_dashboard.R` as a reactable cell renderer.
#
# All inputs are probabilities on [0, 1] or NA. Output is a character string.

render_books_strip <- function(model, dk, fd, px, nv, cons) {
  fmt_pill <- function(label, prob, css_class) {
    if (is.na(prob)) {
      sprintf('<span class="pill %s dim">%s &mdash;</span>', css_class, label)
    } else {
      sprintf('<span class="pill %s">%s %.1f</span>', css_class, label, prob * 100)
    }
  }

  paste0(
    '<span class="books-strip">',
    fmt_pill("M",    model, "model"),
    fmt_pill("DK",   dk,    "book"),
    fmt_pill("FD",   fd,    "book"),
    fmt_pill("PX",   px,    "book"),
    fmt_pill("NV",   nv,    "book"),
    fmt_pill("Cons", cons,  "cons"),
    '</span>'
  )
}
```

- [ ] **Step 4: Run the tests and confirm they pass**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-books-strip/Answer Keys/tests"
Rscript -e 'library(testthat); test_file("test_books_strip.R")'
```

Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 5 ]` (or whatever the testthat success summary looks like — five `test_that` blocks, all green).

- [ ] **Step 5: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-books-strip"
git add "Answer Keys/books_strip.R" "Answer Keys/tests/test_books_strip.R"
git commit -m "$(cat <<'EOF'
feat(mlb-dash): pure render_books_strip helper + tests

Renders the parlay tab's pill row (M / DK / FD / PX / NV / Cons) as an HTML
string. NA probabilities render as a dimmed em-dash; numeric values format
to one decimal percent. testthat covers all-values, NA dimming, all-but-model
NA, rounding, and container wrapping.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Wire helper into the dashboard parlay table

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` — source the helper, replace four per-book `colDef`s with one Books `colDef`, hide `model_prob_raw` (already hidden), hide `game_time` separately, fold time into Game cell.

- [ ] **Step 1: Source the helper at the top of the dashboard script**

Find the `library(...)` block at lines 4–13 of `Answer Keys/MLB Dashboard/mlb_dashboard.R`. Immediately after the closing `})` (line 13), add:

```r

# Pure HTML helpers shared with other dashboards
source(file.path(DASHBOARD_DIR, "..", "books_strip.R"))
```

Wait — `DASHBOARD_DIR` is defined at lines 19–29, *after* the library block. Move the `source(...)` call to immediately after `DASHBOARD_DIR <- ...` (line 29) instead. The exact insertion point: after line 30 (`DB_PATH <- ...`), add a blank line then:

```r
# Pure HTML helpers shared with other dashboards
source(file.path(DASHBOARD_DIR, "..", "books_strip.R"))
```

- [ ] **Step 2: Hide `game_time` and update Game cell to include time**

In `create_parlays_table()` (starts at line 246), find the Visible columns block (line 326):

```r
      # Visible columns
      game = colDef(name = "Game", minWidth = 180),
      game_time = colDef(
        name = "Time",
        minWidth = 130,
        cell = JS("function(cellInfo) {
          var val = cellInfo.value;
          if (!val || val === '') return '';
          var d = new Date(val);
          if (isNaN(d)) return val;
          var opts = {weekday:'short', month:'2-digit', day:'2-digit', hour:'numeric', minute:'2-digit'};
          return d.toLocaleString(undefined, opts);
        }")
      ),
```

Replace with:

```r
      # Visible columns
      game = colDef(
        name = "Game",
        minWidth = 180,
        html = TRUE,
        cell = JS("function(cellInfo) {
          var matchup = cellInfo.value || '';
          var t = cellInfo.row.game_time;
          if (!t) return matchup;
          var d = new Date(t);
          if (isNaN(d)) return matchup;
          var opts = {weekday:'short', month:'2-digit', day:'2-digit', hour:'numeric', minute:'2-digit'};
          var time = d.toLocaleString(undefined, opts);
          return matchup + '<div style=\"color:#8b949e;font-size:10px\">' + time + '</div>';
        }")
      ),
      game_time = colDef(show = FALSE),
```

- [ ] **Step 3: Replace the four per-book `colDef`s with a single Books `colDef`**

In the same `columns = list(...)` block, find lines 343–356:

```r
      # Per-book devigged fair probabilities (NA when the book didn't price the combo).
      # Muted grey — they're diagnostic detail; the primary "Fair" col is the blend.
      dk_fair_prob = colDef(name = "DK", minWidth = 60, align = "right",
        format = colFormat(percent = TRUE, digits = 1),
        style = list(fontFamily = "monospace", color = "#8b949e")),
      fd_fair_prob = colDef(name = "FD", minWidth = 60, align = "right",
        format = colFormat(percent = TRUE, digits = 1),
        style = list(fontFamily = "monospace", color = "#8b949e")),
      px_fair_prob = colDef(name = "PX", minWidth = 60, align = "right",
        format = colFormat(percent = TRUE, digits = 1),
        style = list(fontFamily = "monospace", color = "#8b949e")),
      nv_fair_prob = colDef(name = "NV", minWidth = 60, align = "right",
        format = colFormat(percent = TRUE, digits = 1),
        style = list(fontFamily = "monospace", color = "#8b949e")),
```

Replace with:

```r
      # Hide the per-book numeric columns — they're now rendered inline in the Books pill row.
      dk_fair_prob = colDef(show = FALSE),
      fd_fair_prob = colDef(show = FALSE),
      px_fair_prob = colDef(show = FALSE),
      nv_fair_prob = colDef(show = FALSE),
      # Combined Books cell: M / DK / FD / PX / NV / Cons pill row.
      # Reads model_prob_raw (always non-NA in real data — pricer skips no-sample
      # games upstream) + the four per-book devigged probs + blended_prob_raw.
      books_strip = colDef(
        name = "Books (devigged fair %)",
        minWidth = 320,
        html = TRUE,
        sortable = FALSE,
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

The `books_strip` colDef is a **synthetic column** — it doesn't correspond to a column in `table_data`. Reactable allows this if the `cell` function returns content; the column header still renders. To make this work, we need to add a placeholder column to `table_data`.

In `create_parlays_table()`, find the `table_data <- parlay_opps %>% mutate(...)` block (starts at line 255). At the end of the `mutate()` call (just before `arrange(desc(edge_pct))` at line 284), add a new column:

```r
      books_strip   = NA  # placeholder; rendered by the colDef's cell function
```

That is, change:

```r
      corr_display   = sprintf("%.3f", corr_factor)
    ) %>%
    arrange(desc(edge_pct))
```

to:

```r
      corr_display   = sprintf("%.3f", corr_factor),
      books_strip    = NA  # placeholder; rendered by the colDef's cell function
    ) %>%
    arrange(desc(edge_pct))
```

- [ ] **Step 4: Smoke-test the dashboard generates without error**

The dashboard reads from a live DuckDB. As long as `mlb_parlay_opportunities` has at least one row, this will exercise the pill renderer:

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-books-strip/Answer Keys/MLB Dashboard"
Rscript mlb_dashboard.R 2>&1 | tail -20
```

Expected: no R errors. The script normally writes `report.html` next to itself; final lines should mention rows written and `report.html` saved.

If `mlb_parlay_opportunities` is empty (off-season / no MLB games today), this check is **inconclusive** — proceed to visual verification in Task 4 against any saved fixture, and skip this step. Do not fake data into the table to force the smoke test.

- [ ] **Step 5: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-books-strip"
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "$(cat <<'EOF'
feat(mlb-dash): swap per-book columns for Books pill row

Replace the four narrow DK/FD/PX/NV probability columns and the standalone
Time column with one combined Books cell containing an inline pill row
(M / DK / FD / PX / NV / Cons). Time folds into the Game cell as a second
line. Net column count: 14 → 10.

Sources Answer Keys/books_strip.R for the pure HTML renderer.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: CSS for pill styling + phone-responsive Corr hide

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` — add rules to the inline `<style>` block at line 746.

- [ ] **Step 1: Locate the inline style block**

Open `Answer Keys/MLB Dashboard/mlb_dashboard.R`, find the line:

```r
      tags$style(HTML('
        * { box-sizing: border-box; }
```

(currently at line 746–747). The HTML string continues for many lines before its closing `'))`. Find the closing — search for the line `        '))` after the opening of the style block. We're going to insert new CSS just before that closing.

Run from the worktree to find the exact closing line number:

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-books-strip"
awk 'NR>=746 && /^        \047\)\),?$/ {print NR; exit}' "Answer Keys/MLB Dashboard/mlb_dashboard.R"
```

Expected: a line number — the line containing the closing `'))` of the style block. Call this `$STYLE_END`.

- [ ] **Step 2: Insert pill + responsive CSS before the closing `'))`**

Immediately before the line at `$STYLE_END`, insert these CSS rules (still inside the HTML(' ... ') string):

```css
        /* Parlay tab — books strip (M / DK / FD / PX / NV / Cons pill row) */
        .books-strip {
          display: flex;
          flex-wrap: wrap;
          gap: 6px;
          align-items: center;
        }
        .pill {
          padding: 2px 6px;
          border-radius: 3px;
          font-family: monospace;
          font-size: 10px;
          color: #8b949e;
          background: #21262d;
          white-space: nowrap;
        }
        .pill.model {
          background: #1f3a2c;
          color: #7ee787;
        }
        .pill.book {
          /* default styling; explicit class for selector clarity */
        }
        .pill.cons {
          background: #1c2738;
          border-left: 2px solid #58a6ff;
          color: #79c0ff;
        }
        .pill.dim {
          opacity: 0.4;
        }

        /* Hide Corr column on phones (low-information; Edge stays visible) */
        @media (max-width: 700px) {
          /* Reactable column data attribute is the dataframe column name */
          [data-col="corr_display"] { display: none !important; }
        }
```

The mobile selector relies on reactable's auto-applied `data-col` attribute. If reactable doesn't emit it in this version, the safer fallback is to use the column index — but reactable >= 0.4 does emit `data-col`. The dashboard's existing `library(reactable)` is sufficient.

- [ ] **Step 3: Smoke-run the dashboard again**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-books-strip/Answer Keys/MLB Dashboard"
Rscript mlb_dashboard.R 2>&1 | tail -10
```

Expected: still no errors; `report.html` regenerates.

- [ ] **Step 4: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-books-strip"
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "$(cat <<'EOF'
style(mlb-dash): pill styling + phone responsive Corr hide

Adds .books-strip flex-wrap container and .pill / .pill.model / .pill.book /
.pill.cons / .pill.dim rules to the dashboard's inline style block. Adds a
media query at <=700px that hides the Corr column so the row stays scrollable
on phone width without horizontal overflow.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Visual verification across viewports

**Files:** none modified — testing only.

- [ ] **Step 1: Make sure `report.html` exists and is fresh**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-books-strip/Answer Keys/MLB Dashboard"
ls -la report.html
```

If it's stale or missing, regenerate:

```bash
Rscript mlb_dashboard.R
```

If `mlb_parlay_opportunities` is empty (off-season), copy the most recent saved `report.html` from `main` and run the dashboard with the same fixture. If neither path is available, document the inconclusive state in the merge PR and defer visual check until live data exists.

- [ ] **Step 2: Open `report.html` in a browser and switch to the Parlays tab**

```bash
open "/Users/callancapitolo/NFLWork/.worktrees/parlay-books-strip/Answer Keys/MLB Dashboard/report.html"
```

Click the **Parlays** tab. Confirm:

- The Parlay Opportunities table renders.
- Each row shows a pill row with six pills: `M` (green tint), `DK`/`FD`/`PX`/`NV` (grey), `Cons` (blue tint with left border).
- The Game column shows the matchup on line 1 and the formatted time on line 2 in muted grey.
- No standalone Time, DK, FD, PX, NV columns are visible.

- [ ] **Step 3: Resize the browser and confirm responsive behavior**

Resize to each of the following widths (use DevTools device toolbar or just drag the window). Confirm at each:

| Width | Expected behavior |
|-------|-------------------|
| 1400px | All six pills on one line per row. All 10 columns visible. No horizontal scroll. |
| 1024px (laptop) | Pills may begin wrapping to a 2nd line in the Books cell. All columns still visible. |
| 768px (tablet) | Pills wrap onto 2 lines. Layout still cleanly readable. |
| 414px (phone) | Pills wrap; **Corr column hidden** (verify by checking column headers). Container scrolls horizontally only as a last resort, not for the pill row itself. |

- [ ] **Step 4: Sanity-check the consensus pill matches the existing Fair column**

Pick any row. Convert the Fair column's American odds to a probability, e.g. `+265 → 100/(265+100) = 27.4%`. Compare to the Cons pill on that row — they should agree to within 0.1pp (rounding).

If they disagree by more, the pill cell is reading the wrong column — double-check `cons = row$blended_prob_raw` in Task 2 Step 3.

- [ ] **Step 5: Place a test parlay**

In the dashboard, click **Place** on any parlay. In the modal, enter a small size (e.g., $1) and submit. Confirm:

- A toast appears: "Parlay placed".
- The button label flips to "Placed".
- The row appears in the **Placed Parlays** section above the table.

Then click **Placed** to remove it (so this is non-destructive). Confirm the toast "Parlay removed" appears and the button flips back. This validates the row's data attributes and Place/Remove flow are unchanged.

- [ ] **Step 6: NA-pill check (optional, only if a real NA exists)**

Look for any row in the table where one of the four book pills shows a dimmed em-dash (`—`). If found, take a screenshot for the merge note. If no NA pills are present in the current data, skip — the unit test in Task 1 covers the rendering path.

No commit for this task — verification only. If any check fails, fix and re-commit under Task 2 or 3 as appropriate, then re-run verification.

---

## Task 5: README update + commit

**Files:**
- Modify: `Answer Keys/MLB Dashboard/README.md`

- [ ] **Step 1: Read the current Features section**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-books-strip"
sed -n '23,32p' "Answer Keys/MLB Dashboard/README.md"
```

Expected output (lines 23–32):

```
## Features

- Filterable bet table (book, market, EV threshold, correlation status)
- **Parlay tab** — MLB-specific: correlated 2-leg parlays (spread + total) priced via `mlb_correlated_parlay.R` with conditional Kelly sizing
- Kelly sizing with configurable bankroll + multiplier
- Same-game correlation detection with visual tooltips
- Bet placement tracking (placed vs recommended)
- CLV (Closing Line Value) computation post-game via `MLB Answer Key/clv_compute.py`
- Auto-place integration via Playwright for supported books (wagerzon, hoop88, bfa, betonlineag)
```

- [ ] **Step 2: Replace the Parlay tab bullet with an expanded version**

Find the line:

```
- **Parlay tab** — MLB-specific: correlated 2-leg parlays (spread + total) priced via `mlb_correlated_parlay.R` with conditional Kelly sizing
```

Replace with:

```
- **Parlay tab** — MLB-specific: correlated 2-leg parlays (spread + total) priced via `mlb_correlated_parlay.R` with conditional Kelly sizing. Each row shows a single "Books" pill cell with our model's joint probability (M), the four per-book devigged fair probabilities (DK / FD / PX / NV), and the blended consensus (Cons) — making model-vs-market disagreement visible at a glance and keeping the table readable on phone and laptop widths.
```

- [ ] **Step 3: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-books-strip"
git add "Answer Keys/MLB Dashboard/README.md"
git commit -m "$(cat <<'EOF'
docs(mlb-dash): README note on books strip + Cons

Document the new Parlay tab "Books" cell — model + four book devigged probs +
blended consensus — and the readability win on phone/laptop viewports.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Pre-merge review, merge, cleanup

**Files:** none — git workflow only.

- [ ] **Step 1: Re-run unit tests from the worktree**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-books-strip/Answer Keys/tests"
Rscript -e 'library(testthat); test_file("test_books_strip.R")'
```

Expected: all 5 `test_that` blocks pass.

- [ ] **Step 2: Re-run the dashboard end-to-end**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-books-strip/Answer Keys/MLB Dashboard"
Rscript mlb_dashboard.R 2>&1 | tail -20
```

Expected: no errors. Confirm `report.html` modification time is fresh.

- [ ] **Step 3: Executive engineer review of the full diff**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-books-strip"
git diff main..HEAD --stat
git diff main..HEAD
```

Verify against the CLAUDE.md pre-merge checklist:

- **Data integrity:** No new writes to DuckDB. Read paths unchanged. ✓ (only presentation layer)
- **Resource safety:** No new DB connections. ✓
- **Edge cases:** NA pills covered; off-season run produces empty table without error (existing behavior, untouched). ✓
- **Dead code:** Confirm `dk_fair_prob`, `fd_fair_prob`, `px_fair_prob`, `nv_fair_prob` are still referenced (just hidden) — do **not** drop them from `table_data` because the `books_strip` cell renderer reads them by row.
- **Log/disk hygiene:** No new files written. ✓
- **Security:** No new env/secret/log content. ✓

If any check fails, fix on the feature branch before proceeding.

- [ ] **Step 4: Hand off to user for explicit merge approval**

Per `CLAUDE.md` Branch hygiene rules, do **not** merge without explicit user approval. Report:

> "Feature branch `feature/parlay-books-strip` is ready. Tests pass, dashboard runs cleanly, viewport checks done at 1400/1024/768/414px. Diff summary: [paste the `git diff main..HEAD --stat` output]. May I merge to `main`?"

Wait for explicit "yes" / "merge" / equivalent before continuing.

- [ ] **Step 5: Merge to main (only after approval)**

```bash
cd /Users/callancapitolo/NFLWork
git checkout main
git merge --no-ff feature/parlay-books-strip -m "Merge feat/parlay-books-strip: surface model prob + collapse per-book columns into pill row"
```

- [ ] **Step 6: Cleanup worktree and feature branch**

```bash
cd /Users/callancapitolo/NFLWork
git worktree remove .worktrees/parlay-books-strip
git branch -d feature/parlay-books-strip
git worktree list
```

Expected: only the main checkout listed; no `.worktrees/parlay-books-strip` entry.

- [ ] **Step 7: Confirm clean state**

```bash
git status
git log --oneline -6
```

Expected: clean working tree on `main`; the merge commit appears at the top of the log alongside the four feature commits.

---

## Self-Review

**Spec coverage check:**

- ✅ Surface model probability → Task 1 (helper) + Task 2 (wiring)
- ✅ Pill row with M / DK / FD / PX / NV / Cons → Task 1 (helper) + Task 3 (CSS)
- ✅ Fold Time into Game cell → Task 2 Step 2
- ✅ Hide per-book columns → Task 2 Step 3
- ✅ Responsive wrap (no horizontal scroll on phone) → Task 3 (`flex-wrap: wrap`) + Task 4 verification
- ✅ Hide Corr on phones → Task 3 Step 2 + Task 4 Step 3
- ✅ NA pills dimmed with em-dash → Task 1 helper, covered in tests
- ✅ Sort default unchanged (`desc(edge_pct)`) → Task 2 leaves arrange clause intact, books_strip colDef sets `sortable = FALSE`
- ✅ Place button regression → Task 4 Step 5
- ✅ README updated → Task 5
- ✅ Worktree lifecycle → Task 0 (create) + Task 6 (cleanup)
- ✅ Pre-merge review + explicit user approval → Task 6 Steps 3–5
- ✅ Cons matches Fair column logic → Task 4 Step 4 sanity check

**Type/name consistency:** `render_books_strip(model, dk, fd, px, nv, cons)` signature is identical in Task 1 (definition + tests), Task 2 Step 3 (call site), and the README mention in Task 5. CSS classes (`books-strip`, `pill`, `pill.model`, `pill.book`, `pill.cons`, `pill.dim`) match between Task 1 (HTML output), Task 3 (CSS rules), and the unit tests.

**Placeholder scan:** No "TODO", "TBD", "fill in details". Every step contains either an exact command or a complete code block. The one explicit unknown (`$STYLE_END` line number) has a deterministic awk command that produces it.
