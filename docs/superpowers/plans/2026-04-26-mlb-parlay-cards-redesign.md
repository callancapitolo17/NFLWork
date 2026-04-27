# MLB Parlay Tab — Cards Redesign Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Convert the Parlay Opportunities reactable from a horizontal table into a vertical list of cards via scoped CSS, drop the Corr and n_books_blended columns, and adopt a 2-size font system (14px / 12px) so the parlay tab reads cleanly at any viewport from a wide laptop down through split-screen to phone width without horizontal scroll.

**Architecture:** CSS-only flip. Reactable still emits a real `<table>`; a CSS block scoped under `#parlays-table-container` overrides `display: table` → `display: flex`-card on `<tr>` and stacks `<td>` cells inside. Per-column targeting comes from `class = "cell-<name>"` parameters added to each visible `colDef`. No R logic changes, no JS rewriting, no data flow changes. Singles tab is untouched (different container ID).

**Tech Stack:** R, reactable, htmltools. Single dashboard file (`Answer Keys/MLB Dashboard/mlb_dashboard.R`) plus its README.

**Spec:** [`docs/superpowers/specs/2026-04-26-mlb-parlay-cards-redesign-design.md`](../specs/2026-04-26-mlb-parlay-cards-redesign-design.md)

---

## File Structure

| Action | Path | Responsibility |
|--------|------|----------------|
| Modify | `Answer Keys/MLB Dashboard/mlb_dashboard.R` | (a) hide three columns from `create_parlays_table()`, (b) add `class = "cell-<name>"` to each visible parlay `colDef`, (c) drop the obsolete `.corr-col` rule and add the cards CSS block to the inline `<style>` block |
| Modify | `Answer Keys/MLB Dashboard/README.md` | Update the Features bullet describing the parlay tab to mention card layout and universal-viewability |

No new files.

---

## Version Control

- **Branch:** `feature/parlay-cards-redesign` (already created in `.worktrees/parlay-cards-redesign`, rebased onto current main `b31dbdb`)
- **Worktree:** `.worktrees/parlay-cards-redesign` — created during brainstorm; **do NOT recreate it**. Cleanup happens in Task 5.
- **Commits (4 atomic):**
  1. `feat(mlb-dash): hide Corr / n_books_blended / is_combo from parlay table` (Task 1)
  2. `feat(mlb-dash): add cell-<name> classes for card-layout targeting` (Task 2)
  3. `style(mlb-dash): card layout for parlay opportunities` (Task 3)
  4. `docs(mlb-dash): README note on parlay card layout` (Task 4)
- **Merge to main:** Task 5, after the executive engineer review and explicit user approval.

---

## Task 1: Hide unwanted columns + drop obsolete `.corr-col` rule

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` — remove `class`/`headerClass` from `corr_display`, add three `colDef(show = FALSE)` entries, remove the `.corr-col` media query.

- [ ] **Step 1: Confirm the worktree branch is correct**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-cards-redesign"
git branch --show-current
```

Expected: `feature/parlay-cards-redesign`. If anything else, stop and ask.

- [ ] **Step 2: Hide `corr_display` (was visible, now hidden — Corr column gone from view)**

In `Answer Keys/MLB Dashboard/mlb_dashboard.R`, locate the `corr_display` colDef in `create_parlays_table()` (currently around line 525). It looks like:

```r
      corr_display = colDef(
        name = "Corr",
        minWidth = 65,
        align = "right",
        class = "corr-col",
        headerClass = "corr-col",
        style = list(color = "#8b949e")
      ),
```

Replace the entire colDef with a single line:

```r
      corr_display = colDef(show = FALSE),
```

- [ ] **Step 3: Hide `n_books_blended` and `is_combo` (auto-rendering today)**

Both columns are in the dataframe but have no explicit colDef, so reactable auto-renders them. Add explicit hides. Find the block of hidden columns near the top of the `columns = list(...)` in `create_parlays_table()` (currently around lines 393–425, the comment header is `# Hidden columns for JS data attributes`).

Add these two lines anywhere inside that hidden-columns block (place them after `parlay_hash = colDef(show = FALSE),` for grouping; the order within the list does not affect rendering):

```r
      n_books_blended = colDef(show = FALSE),
      is_combo = colDef(show = FALSE),
```

- [ ] **Step 4: Remove the obsolete `.corr-col` media query from the inline style block**

The phone-hide rule was added in the books-strip round to hide the Corr column on phones. Now that the column is hidden at all widths, the rule is dead code. Locate it in the inline `<style>` block (currently around lines 1608–1614):

```css
        /* Hide Corr column on phones (low-information; Edge stays visible).
           reactable applies class+headerClass from colDef to both data cells
           and the header cell, so .corr-col matches both. !important overrides
           reactable inline display:table-cell on cells. */
        @media (max-width: 700px) {
          .corr-col { display: none !important; }
        }
```

Delete those 7 lines (the comment block and the media query). Leave the surrounding lines untouched.

- [ ] **Step 5: Syntax check**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-cards-redesign"
Rscript -e 'parse("Answer Keys/MLB Dashboard/mlb_dashboard.R"); cat("SYNTAX OK\n")'
```

Expected: `SYNTAX OK`.

- [ ] **Step 6: Smoke run (against a copied DB)**

The worktree won't have `mlb.duckdb` since CLAUDE.md prohibits symlinking DuckDB files. Copy fresh copies in for testing:

```bash
cd "/Users/callancapitolo/NFLWork"
cp "Answer Keys/mlb.duckdb" ".worktrees/parlay-cards-redesign/Answer Keys/mlb.duckdb"
cp "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb" ".worktrees/parlay-cards-redesign/Answer Keys/MLB Dashboard/mlb_dashboard.duckdb"
```

Then run the dashboard end-to-end:

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-cards-redesign/Answer Keys/MLB Dashboard"
Rscript mlb_dashboard.R 2>&1 | tail -10
```

Expected: no R errors. Final lines mention rows loaded and `report.html` saved.

If `mlb_parlay_opportunities` is empty (off-season / no MLB games today), the run is **inconclusive but not a failure** — proceed and note it. Visual check is in Task 4.

- [ ] **Step 7: Confirm Corr / n_books / is_combo are absent from `report.html`**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-cards-redesign/Answer Keys/MLB Dashboard"
grep -oE '"name":"Corr"|"name":"n_books_blended"|"name":"is_combo"' report.html | sort -u
```

Expected: empty output (none of those headers should appear in the rendered HTML). If any appear, the corresponding `colDef(show = FALSE)` didn't take effect — re-check Step 2 or Step 3.

- [ ] **Step 8: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-cards-redesign"
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "$(cat <<'EOF'
feat(mlb-dash): hide Corr / n_books_blended / is_combo from parlay table

Drops the Corr column from view (low signal), and explicitly hides
n_books_blended + is_combo which were auto-rendering because they had
no colDef entry. Also removes the .corr-col media query — dead now that
the column is hidden at all widths.

Prep for the parlay-card layout change in the next commit.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Add `class = "cell-<name>"` to each visible parlay colDef

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` — add `class` parameter to each of the 10 visible parlay colDefs.

These class names give the CSS in Task 3 stable selectors to position each cell inside the card. Reactable applies the `class` parameter to the `<td>` element so the selector `.rt-td.cell-game` matches each Game cell.

- [ ] **Step 1: Add `class = "cell-sel"` to the `sel` colDef**

In `create_parlays_table()`, find the `sel = colDef(` block (currently around line 424). It currently looks like:

```r
      sel = colDef(
        name = "",
        minWidth = 30,
        align = "center",
        filterable = FALSE,
        sortable = FALSE,
        html = TRUE,
        cell = function(value, index) { ... }
      ),
```

Add `class = "cell-sel",` immediately after `sortable = FALSE,`. Result:

```r
      sel = colDef(
        name = "",
        minWidth = 30,
        align = "center",
        filterable = FALSE,
        sortable = FALSE,
        class = "cell-sel",
        html = TRUE,
        cell = function(value, index) { ... }
      ),
```

(Don't touch the `cell = function(...) {...}` body. Leave it as-is; only insert the new line.)

- [ ] **Step 2: Add `class = "cell-game"` to the `game` colDef**

Find the `game = colDef(` in `create_parlays_table()` (the one with `html = TRUE` and the JS cell renderer that folds the time line in — currently around line 478, NOT the simpler `game = colDef(name = "Game", minWidth = 200)` at line 168 which is a different function). Add `class = "cell-game",` after `minWidth = 180,`. Result:

```r
      game = colDef(
        name = "Game",
        minWidth = 180,
        class = "cell-game",
        html = TRUE,
        cell = JS("function(cellInfo) { ... }")
      ),
```

- [ ] **Step 3: Add `class = "cell-legs"` to `legs_display`**

Currently at line 494:

```r
      legs_display = colDef(name = "Legs", minWidth = 260),
```

Add `class = "cell-legs"`:

```r
      legs_display = colDef(name = "Legs", minWidth = 260, class = "cell-legs"),
```

- [ ] **Step 4: Add `class = "cell-fair"` to `fair_display`**

Currently at line 495:

```r
      fair_display = colDef(name = "Fair", minWidth = 70, align = "right",
        style = list(fontFamily = "monospace", color = "#8b949e")),
```

Add `class = "cell-fair"`:

```r
      fair_display = colDef(name = "Fair", minWidth = 70, align = "right",
        class = "cell-fair",
        style = list(fontFamily = "monospace", color = "#8b949e")),
```

- [ ] **Step 5: Add `class = "cell-books"` to `books_strip`**

Currently around line 506:

```r
      books_strip = colDef(
        name = "Books (devigged fair %)",
        minWidth = 320,
        html = TRUE,
        sortable = FALSE,
        cell = function(value, index) { ... }
      ),
```

Add `class = "cell-books",` after `sortable = FALSE,`:

```r
      books_strip = colDef(
        name = "Books (devigged fair %)",
        minWidth = 320,
        html = TRUE,
        sortable = FALSE,
        class = "cell-books",
        cell = function(value, index) { ... }
      ),
```

- [ ] **Step 6: Add `class = "cell-wz"` to `wz_display`**

Currently at line 523:

```r
      wz_display = colDef(name = "WZ", minWidth = 70, align = "right",
        style = list(fontFamily = "monospace")),
```

Add the class:

```r
      wz_display = colDef(name = "WZ", minWidth = 70, align = "right",
        class = "cell-wz",
        style = list(fontFamily = "monospace")),
```

- [ ] **Step 7: Add `class = "cell-edge"` to `edge_display`**

Currently around line 533:

```r
      edge_display = colDef(name = "Edge %", minWidth = 70, align = "right",
        cell = function(value, index) {
          ep <- table_data$edge_pct[index]
          color <- if (ep >= 15) "#3fb950" else if (ep >= 10) "#56d364" else if (ep >= 5) "#7ee787" else "#a5d6a7"
          span(style = list(color = color, fontWeight = "600"), value)
        }
      ),
```

Add `class = "cell-edge",` after `align = "right",`:

```r
      edge_display = colDef(name = "Edge %", minWidth = 70, align = "right",
        class = "cell-edge",
        cell = function(value, index) {
          ep <- table_data$edge_pct[index]
          color <- if (ep >= 15) "#3fb950" else if (ep >= 10) "#56d364" else if (ep >= 5) "#7ee787" else "#a5d6a7"
          span(style = list(color = color, fontWeight = "600"), value)
        }
      ),
```

- [ ] **Step 8: Add `class = "cell-size"` to `size_display`**

Currently at line 540:

```r
      size_display = colDef(name = "Size", minWidth = 65, align = "right", html = TRUE),
```

Add the class:

```r
      size_display = colDef(name = "Size", minWidth = 65, align = "right", html = TRUE, class = "cell-size"),
```

- [ ] **Step 9: Add `class = "cell-towin"` to `to_win_display`**

Currently at lines 541–542:

```r
      to_win_display = colDef(name = "To Win", minWidth = 65, align = "right",
        style = list(color = "#3fb950")),
```

Add the class:

```r
      to_win_display = colDef(name = "To Win", minWidth = 65, align = "right",
        class = "cell-towin",
        style = list(color = "#3fb950")),
```

- [ ] **Step 10: Add `class = "cell-action"` to `is_placed`**

The Action cell — find the `is_placed = colDef(` block (currently around line 548). Looks like:

```r
      is_placed = colDef(
        name = "Action",
        minWidth = 110,
        align = "center",
        filterable = FALSE,
        html = TRUE,
        cell = function(value, index) { ... long body ... }
      )
```

Add `class = "cell-action",` after `filterable = FALSE,`:

```r
      is_placed = colDef(
        name = "Action",
        minWidth = 110,
        align = "center",
        filterable = FALSE,
        class = "cell-action",
        html = TRUE,
        cell = function(value, index) { ... long body ... }
      )
```

(Leave the `cell = function(...) {...}` body untouched.)

- [ ] **Step 11: Syntax + smoke run**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-cards-redesign"
Rscript -e 'parse("Answer Keys/MLB Dashboard/mlb_dashboard.R"); cat("SYNTAX OK\n")'
cd "Answer Keys/MLB Dashboard"
Rscript mlb_dashboard.R 2>&1 | tail -10
```

Expected: SYNTAX OK; dashboard renders without errors; `report.html` is fresh.

- [ ] **Step 12: Confirm the cell classes appear in the rendered HTML**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-cards-redesign/Answer Keys/MLB Dashboard"
grep -oE 'cell-(sel|game|legs|fair|books|wz|edge|size|towin|action)' report.html | sort -u
```

Expected: 10 distinct class names, one per visible colDef.

- [ ] **Step 13: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-cards-redesign"
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "$(cat <<'EOF'
feat(mlb-dash): add cell-<name> classes for card-layout targeting

Adds class = "cell-<name>" to each of the 10 visible parlay colDefs
(sel, game, legs, fair, books, wz, edge, size, towin, action).
reactable propagates the class onto each <td>, giving the upcoming
card-layout CSS stable selectors to position individual cells inside
the per-row card container.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Card-layout CSS

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` — insert the cards CSS block into the inline `<style>` block, just before its closing `'))` line.

- [ ] **Step 1: Locate the closing `'))` of the inline style block**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-cards-redesign"
grep -nE "^[[:space:]]+'\)\)" "Answer Keys/MLB Dashboard/mlb_dashboard.R" | head -3
```

Expected: at least two hits — the first one is the close of the dashboard's main inline style block (this is the one we want); the second is the close of a different `tags$style(HTML(...))` further down (don't touch). Note the **first** line number — call it `$STYLE_END`.

(After Task 1's deletions, this line will have shifted up by ~7 lines from the previous round. After Task 2's additions it doesn't shift further since Task 2 didn't change line counts inside the style block.)

- [ ] **Step 2: Insert the cards CSS just before `$STYLE_END`**

Use the Edit tool to insert these rules. The insertion goes inside the HTML(' ... ') string, immediately before the closing `'))` line. Match the existing 8-space indentation of CSS rules in the block.

Insert this block:

```css
        /* === Parlay tab card layout (scoped — singles tab is unaffected) === */
        /* Flatten the reactable table into a stack of cards. */
        #parlays-table-container .rt-table   { display: block; }
        #parlays-table-container .rt-thead   { display: none; }
        #parlays-table-container .rt-tbody   { display: block; }
        #parlays-table-container .rt-tr-group { display: block; }

        /* Each row becomes a card; cells flex inside so DOM order = visual order. */
        #parlays-table-container .rt-tr {
          display: flex;
          flex-wrap: wrap;
          align-items: baseline;
          gap: 4px;
          background: #161b22;
          border: 1px solid #21262d;
          border-radius: 6px;
          padding: 14px 14px 12px 14px;
          margin-bottom: 10px;
          position: relative;
        }

        #parlays-table-container .rt-td {
          display: block;
          border: 0;
          padding: 0;
          white-space: normal;
          font-size: 14px;
          color: #c9d1d9;
        }

        /* Full-width rows inside the flex card */
        #parlays-table-container .rt-td.cell-game,
        #parlays-table-container .rt-td.cell-legs,
        #parlays-table-container .rt-td.cell-books {
          flex-basis: 100%;
        }

        /* Sel checkbox: top-right corner */
        #parlays-table-container .rt-td.cell-sel {
          position: absolute;
          top: 12px;
          right: 12px;
        }
        #parlays-table-container .combo-select {
          width: 18px;
          height: 18px;
          cursor: pointer;
        }

        /* Game + folded time line */
        #parlays-table-container .rt-td.cell-game {
          font-size: 15px;
          font-weight: 500;
          padding-right: 36px;
        }

        /* Legs */
        #parlays-table-container .rt-td.cell-legs {
          font-size: 14px;
          margin-top: 2px;
          margin-bottom: 10px;
        }

        /* Books pill row */
        #parlays-table-container .rt-td.cell-books {
          margin-bottom: 12px;
        }

        /* Metadata strip — Fair / WZ / Size / To Win flow inline as flex items.
           Parent gap: 4px handles base spacing; margin-right adds breathing room. */
        #parlays-table-container .rt-td.cell-fair,
        #parlays-table-container .rt-td.cell-wz,
        #parlays-table-container .rt-td.cell-size,
        #parlays-table-container .rt-td.cell-towin {
          display: inline-flex;
          align-items: baseline;
          gap: 4px;
          margin-right: 10px;
          font-size: 13px;
          font-family: monospace;
        }

        #parlays-table-container .rt-td.cell-fair::before  { content: "Fair";   color: #8b949e; font-size: 12px; font-family: -apple-system, system-ui, sans-serif; }
        #parlays-table-container .rt-td.cell-wz::before    { content: "WZ";     color: #8b949e; font-size: 12px; font-family: -apple-system, system-ui, sans-serif; }
        #parlays-table-container .rt-td.cell-size::before  { content: "Size";   color: #8b949e; font-size: 12px; font-family: -apple-system, system-ui, sans-serif; }
        #parlays-table-container .rt-td.cell-towin::before { content: "To Win"; color: #8b949e; font-size: 12px; font-family: -apple-system, system-ui, sans-serif; }

        /* Edge — pushed to right via margin-left:auto on first of two right items.
           Action follows in DOM order, so visually:
           [...metadata strip...]                       +12.3%  [ Place ] */
        #parlays-table-container .rt-td.cell-edge {
          margin-left: auto;
          font-size: 15px;
          font-weight: 600;
        }
        #parlays-table-container .rt-td.cell-action {
          margin-left: 12px;
        }

        /* Conditional-Kelly residual note — full-width line below metadata strip */
        #parlays-table-container .combo-note {
          display: block;
          width: 100%;
          margin-top: 4px;
          color: #8b949e;
          font-size: 12px;
          font-style: italic;
        }

        /* Pill upsizing (was 10px in the books-strip round) */
        #parlays-table-container .pill {
          font-size: 13px;
          padding: 3px 8px;
        }
```

- [ ] **Step 3: Syntax + smoke run**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-cards-redesign"
Rscript -e 'parse("Answer Keys/MLB Dashboard/mlb_dashboard.R"); cat("SYNTAX OK\n")'
cd "Answer Keys/MLB Dashboard"
Rscript mlb_dashboard.R 2>&1 | tail -10
```

Expected: SYNTAX OK; `report.html` regenerates without errors.

- [ ] **Step 4: Confirm the new selectors appear in the HTML output**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-cards-redesign/Answer Keys/MLB Dashboard"
grep -oE '#parlays-table-container \.rt-tr|cell-fair::before|cell-action' report.html | sort -u
```

Expected: at least the three substrings present in the inlined CSS.

- [ ] **Step 5: Visual verification (manual)**

Open `report.html` in a browser:

```bash
open "/Users/callancapitolo/NFLWork/.worktrees/parlay-cards-redesign/Answer Keys/MLB Dashboard/report.html"
```

Switch to the **Parlays** tab and verify:

1. Each parlay opportunity renders as a card (no table rows visible). Singles tab still looks like a table.
2. Card top: matchup line + time line + legs line. Top-right corner of each card has the Sel checkbox (only on cards that aren't already placed or aren't a combo).
3. Middle: pills row M / DK / FD / PX / NV / Cons (with Cons in blue).
4. Bottom: a metadata strip with `Fair`, `WZ`, `Size`, `To Win` labels and values. The `(residual after combo $X)` annotation, when present, sits below the metadata strip on its own line.
5. Bottom-right: edge percentage (color-coded green) followed by the Action cell (Place button OR "placed · #ticket" muted label OR red error pill).
6. Resize the browser to ~1400 / ~860 / ~400 px. Cards should scale cleanly. No horizontal scroll on phone width. Pills should wrap to a 2nd line at narrow widths.
7. Tick two Sel checkboxes on different cards. The combo-banner above the cards should appear with combined pricing. Click "Place Combined Parlay" and confirm it lands; click Remove and confirm it reverts.
8. Click Place on any card at $1. Toast appears, button flips to "Placed", placed-parlays section above shows the row. Click Remove. Reverts.

If any check fails, debug; do NOT commit.

- [ ] **Step 6: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-cards-redesign"
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "$(cat <<'EOF'
style(mlb-dash): card layout for parlay opportunities

Adds ~80 lines of CSS scoped under #parlays-table-container that flips
the reactable's table-display rules so each row renders as a vertical
card. Each card stacks game + legs + books pill row + metadata strip
+ edge + action, with the Sel checkbox absolutely positioned top-right.
The CSS is cell-layout only — reactable's R logic, JS handlers, and
data flow are all untouched. Singles tab uses a different container ID
and is unaffected.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: README update

**Files:**
- Modify: `Answer Keys/MLB Dashboard/README.md`

- [ ] **Step 1: Read the current Features section**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-cards-redesign"
sed -n '20,35p' "Answer Keys/MLB Dashboard/README.md"
```

Find the bullet that mentions the Parlay tab — it currently mentions the books strip and per-book devigged probabilities (added in commit `1cc850f`).

- [ ] **Step 2: Replace the Parlay tab bullet**

Find this bullet (the wording was set in the previous round; if it's been edited, match the prefix and rewrite the rest):

```
- **Parlay tab** — MLB-specific: correlated 2-leg parlays (spread + total) priced via `mlb_correlated_parlay.R` with conditional Kelly sizing. Each row shows a single "Books" pill cell with our model's joint probability (M), the four per-book devigged fair probabilities (DK / FD / PX / NV), and the blended consensus (Cons) — making model-vs-market disagreement visible at a glance and keeping the table readable on phone and laptop widths.
```

Replace with:

```
- **Parlay tab** — MLB-specific: correlated 2-leg parlays (spread + total) priced via `mlb_correlated_parlay.R` with conditional Kelly sizing. Each opportunity renders as a card containing the matchup, legs, a Books pill row (model M plus per-book devigged fair probabilities for DK / FD / PX / NV plus blended consensus Cons), and a metadata strip (Fair / WZ / Size / To Win) with edge percentage and the Place / placed-label / error-pill action. The card layout reads identically across laptop, split-screen, and phone — no column hiding, no horizontal scroll. Combined-parlay selection (the Sel checkbox in the top-right corner of each card) and auto-placement still work unchanged.
```

- [ ] **Step 3: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-cards-redesign"
git add "Answer Keys/MLB Dashboard/README.md"
git commit -m "$(cat <<'EOF'
docs(mlb-dash): README note on parlay card layout

Replaces the Parlay tab bullet to describe the card layout (game, legs,
books strip, metadata strip, edge + action), the universal-viewability
across laptop / split-screen / phone, and that combined-parlay selection
+ auto-placement remain unchanged.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Pre-merge review, merge, cleanup

**Files:** none — git workflow only.

- [ ] **Step 1: Re-run the dashboard end-to-end one final time**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-cards-redesign/Answer Keys/MLB Dashboard"
Rscript mlb_dashboard.R 2>&1 | tail -10
```

Expected: no errors. Confirm `report.html` modification time is fresh.

- [ ] **Step 2: Executive engineer review**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/parlay-cards-redesign"
git log --oneline main..HEAD
git diff main..HEAD --stat
git diff main..HEAD
```

Verify against the CLAUDE.md pre-merge checklist:

- **Data integrity:** No new writes. Read paths unchanged. ✓ (presentation-layer only)
- **Resource safety:** No new DB connections. ✓
- **Edge cases:** NA pills handled by existing `render_books_strip` helper from the previous round. Singles tab untouched. Empty table = no `.rt-tr` elements rendered (existing empty-state above the table is unchanged). ✓
- **Dead code:** Confirm `corr_display`, `n_books_blended`, `is_combo` are still in the dataframe (NOT dropped from `table_data`) — they're hidden via `colDef(show = FALSE)` only. The `.corr-col` media query was deleted (now dead, removed). The `class="corr-col"` and `headerClass="corr-col"` parameters on `corr_display` are gone (column hidden, those params irrelevant).
- **Log/disk hygiene:** No new files written. ✓
- **Security:** No new env / secret / log content. ✓

If any check fails, fix on the feature branch before proceeding.

- [ ] **Step 3: Hand off to user for explicit merge approval**

Per `CLAUDE.md` Branch hygiene rules, do **NOT** merge without explicit user approval. Report:

> "Feature branch `feature/parlay-cards-redesign` is ready. Smoke run clean, dashboard renders cards as designed, viewport / Place button / combined-parlay flow all verified. Diff summary: [paste `git diff main..HEAD --stat` output]. May I merge to `main`?"

Wait for explicit "yes" / "merge" / equivalent before continuing.

- [ ] **Step 4: Merge to main (only after approval)**

```bash
cd /Users/callancapitolo/NFLWork
git checkout main
git merge --no-ff feature/parlay-cards-redesign -m "Merge feature/parlay-cards-redesign: card layout for MLB parlay opportunities"
```

- [ ] **Step 5: Push (only if explicitly requested by the user)**

If the user approved a push along with the merge in Step 3, run:

```bash
cd /Users/callancapitolo/NFLWork
git push origin main 2>&1 | tail -5
```

If the user did NOT request a push, stop here and let them decide later.

- [ ] **Step 6: Remove worktree and feature branch (only after merge succeeds)**

```bash
cd /Users/callancapitolo/NFLWork
git worktree remove .worktrees/parlay-cards-redesign
git branch -d feature/parlay-cards-redesign
git worktree list
```

Expected: only the main checkout listed (plus other unrelated stale worktrees if any exist); no `.worktrees/parlay-cards-redesign` entry.

- [ ] **Step 7: Final sanity check**

```bash
cd /Users/callancapitolo/NFLWork
git status
git log --oneline -8
```

Expected: clean working tree on `main`; the merge commit at the top of the log alongside the four feature commits.

---

## Self-Review

**Spec coverage:**

- ✅ Universal viewability across viewports → Task 3 (CSS flex layout that adapts naturally) + Task 3 Step 5 (manual viewport sweep at 1400 / 860 / 400 px)
- ✅ Drop Corr column → Task 1 Step 2
- ✅ Drop n_books_blended (and is_combo) → Task 1 Step 3
- ✅ 2-size font system (14px primary, 12px secondary) → Task 3 CSS (`font-size: 14px` on `.rt-td`, `font-size: 12px` on `::before` labels and the time-line div emitted by `cell-game`)
- ✅ Preserve Sel checkbox + combined-parlay flow → Task 2 (`cell-sel` class added) + Task 3 (CSS positions it absolute top-right) + Task 3 Step 5 verification
- ✅ Preserve Place button + auto-placement Action variants → Task 2 (`cell-action` class added) + existing `.pill.error` and `.placed-parlay-label` CSS untouched
- ✅ Preserve combo_residual_note → Task 3 (`.combo-note` CSS rule sets `display: block; width: 100%`)
- ✅ Preserve edge color thresholds → cell renderer left untouched in Task 2 Step 7
- ✅ Preserve to-win green → cell style left untouched in Task 2 Step 9
- ✅ Singles tab untouched → CSS is scoped to `#parlays-table-container`
- ✅ README update → Task 4
- ✅ Worktree lifecycle → already created during brainstorm; Task 5 Step 6 cleans it up
- ✅ Pre-merge review + explicit user approval → Task 5 Steps 2–4

**Type / name consistency:** Cell class names match between Task 2 (R-side) and Task 3 (CSS): `cell-sel`, `cell-game`, `cell-legs`, `cell-books`, `cell-fair`, `cell-wz`, `cell-edge`, `cell-size`, `cell-towin`, `cell-action`. Container ID `#parlays-table-container` matches what's used in `mlb_dashboard.R` (line 1833 in current main) and the existing JS handlers.

**Placeholder scan:** No "TODO", "TBD", "fill in details", or "similar to". Every step contains either an exact command or a complete code block.

**Risk note for the implementer:** If reactable is upgraded to a future major version that changes its internal class names (`.rt-tr` / `.rt-td` / `.rt-thead` etc.), the CSS will silently stop matching and the layout will revert to a normal table. This is acceptable for now; flag it in a future tech-debt ticket if the team upgrades reactable.
