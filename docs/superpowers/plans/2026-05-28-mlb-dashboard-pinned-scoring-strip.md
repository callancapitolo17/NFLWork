# MLB Dashboard — Pinned Scoring Strip Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Freeze the `Placing on` account pills (all tabs) and the Kelly Calculator (Bets tab) to the top of the MLB Dashboard so they stay visible while scrolling the bet-card list.

**Architecture:** Pure CSS `position: sticky`. The account-pills row is relocated out of the short `.header` to a page-spanning sibling under `.container` so its sticky freeze holds through the whole list; the Kelly Calculator is wrapped in a full-width opaque sticky container that pins just below the account bar, with a JS-measured `top` offset so it stays flush when pills wrap. No server, DB, schema, or Kelly-math changes.

**Tech Stack:** R (`htmltools`/`tags` HTML generation in `mlb_dashboard.R`), inline CSS, vanilla JS. Verification by rendering `report.html` and inspecting in a browser (the repo's established UI-verification pattern — there is no pure function to unit-test here).

**Spec:** `docs/superpowers/specs/2026-05-27-mlb-dashboard-pinned-scoring-strip-design.md`

---

## A note on testing approach

This change is entirely presentation-layer: a DOM relocation, CSS additions, and a
3-line offset helper. There is no pure function to assert on, so classic unit-TDD
does not apply. Instead each task ends with **render `report.html` + visual
inspection in a browser**, which is how dashboard UI is verified in this repo
(see the `verify_ui_features_by_rendering` convention). Do not claim the feature
works from a static diff — render and look.

## File structure

Only one source file changes:

- **Modify:** `Answer Keys/MLB Dashboard/mlb_dashboard.R`
  - Body markup (~`:3122-3139`): relocate the account row; (~`:3452`): wrap the Kelly Calculator.
  - Inline `<style>` block (~`:2372` and ~`:3271`): add `.pinned-account-bar` and `.kelly-calc-pin` rules; add `position: sticky` to the Kelly wrapper.
  - Inline `<script>` (end of `renderPills`, ~`:5374`, plus a new helper): add `syncPinOffset()` and its listeners.
- **Modify (docs, Task 4):** `Answer Keys/MLB Dashboard/README.md`, `Answer Keys/CLAUDE.md`.

Reference facts (verified in the current code, post-rebase onto local `main`):
- `.container` (`:3116`) is the page-spanning wrapper; `.header` (`:3122`) is its short first child; `.tab-bar` (`:3142`) and `#tab-bets` (`:3150`) follow.
- Account row markup is `#wz-account-row` / `.header-row-accounts` (`:3130-3138`); pills render into `#wz-account-pills`.
- `.header-row-accounts` CSS at `:2372`; `.kelly-calc` CSS at `:3271` (`display: inline-flex; max-width: fit-content; margin: 0 0 14px 0`).
- `renderPills()` is defined at `:5341` and ends at `:5374` (then `refreshBalances()` at `:5376`); it is the function that populates the pills (so the bar's height is only correct *after* it runs).
- Existing z-index layers: filter dropdown menu `100` (`:2251`), modal overlay `999` (`:2016`), toast `1000`/`9999` (`:2140`,`:2182`). Our sticky elements use `40` and `30` — above cards, below all of those.

---

## Task 1: Relocate the account row and wrap the Kelly Calculator (markup only)

This task changes DOM structure but adds **no** sticky behavior yet. After it,
the dashboard should look almost identical (the account bar moves a few pixels
down, out of the header border) and all pills/calc behavior is unchanged.

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R:3122-3139` (header / account row)
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R:3452-3483` (Kelly Calculator wrapper)

- [ ] **Step 1: Relocate the account row out of `.header`**

Replace the current header block (`:3122-3139`):

```r
        tags$div(class = "header",
          tags$div(class = "header-row-top",
            tags$div(
              tags$h1("MLB Answer Key Dashboard"),
              tags$div(class = "subtitle", paste("Updated", timestamp))
            ),
            tags$button(class = "refresh-btn", onclick = "refreshData()", "Refresh")
          ),
          tags$div(class = "header-row-accounts", id = "wz-account-row",
            tags$span(class = "header-label", "Placing on"),
            tags$div(id = "wz-account-pills", class = "wz-pills"),
            tags$button(
              id = "wz-refresh-btn", type = "button", class = "wz-icon-btn",
              title = "Refresh balances",
              HTML("&#x21bb;")
            )
          )
        ),
```

with this (account row pulled out to be a sibling of `.header`, before `.tab-bar`,
and given the `pinned-account-bar` class; ids and inner structure preserved):

```r
        tags$div(class = "header",
          tags$div(class = "header-row-top",
            tags$div(
              tags$h1("MLB Answer Key Dashboard"),
              tags$div(class = "subtitle", paste("Updated", timestamp))
            ),
            tags$button(class = "refresh-btn", onclick = "refreshData()", "Refresh")
          )
        ),

        # Pinned account bar — relocated out of .header so its sticky parent is
        # the page-spanning .container; stays frozen through the whole card list
        # on every tab. Ids/inner structure preserved so the pill JS is untouched.
        tags$div(class = "header-row-accounts pinned-account-bar", id = "wz-account-row",
          tags$span(class = "header-label", "Placing on"),
          tags$div(id = "wz-account-pills", class = "wz-pills"),
          tags$button(
            id = "wz-refresh-btn", type = "button", class = "wz-icon-btn",
            title = "Refresh balances",
            HTML("&#x21bb;")
          )
        ),
```

- [ ] **Step 2: Wrap the Kelly Calculator in a full-width sticky container**

The Kelly Calculator div begins at `:3452` with `tags$div(class = "kelly-calc",`.
Because `.kelly-calc` is `display: inline-flex; max-width: fit-content`, pinning
it directly would let cards show through beside the narrow strip. Wrap it in a
full-width container that will carry the sticky + opaque background.

Find the opening line:

```r
          tags$div(class = "kelly-calc",
            tags$span(class = "title", "Kelly Calculator"),
```

Change it to open a wrapper first:

```r
          tags$div(class = "kelly-calc-pin",
          tags$div(class = "kelly-calc",
            tags$span(class = "title", "Kelly Calculator"),
```

Then find the **closing** of the `.kelly-calc` div. It is the `)` that closes the
`tags$div(class = "kelly-calc", ...)` call — currently at `:3483`, immediately
before the `# Filter Bar` comment. The exact current text (note the comment is
indented 8 spaces) is:

```r
            )
          ),

        # Filter Bar
```

Add one more closing paren for the new wrapper (so `.kelly-calc` closes with a
bare `)` and `.kelly-calc-pin` closes with `),`):

```r
            )
          )
          ),

        # Filter Bar
```

(If line numbers have drifted, anchor on the `# Filter Bar` comment: the wrapper's
extra `)` goes on the line just before it, closing `.kelly-calc-pin`.)

- [ ] **Step 3: Render and inspect (structure unchanged, behavior intact)**

The R script resolves its DB path next to itself, so copy the live DBs into the
worktree first (copy, never symlink — DuckDB WAL safety), then render:

```bash
WT="/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-dashboard-pin-kelly-book"
SRC="/Users/callancapitolo/NFLWork/Answer Keys"
cp "$SRC/MLB Dashboard/mlb_dashboard.duckdb" "$WT/Answer Keys/MLB Dashboard/" 2>/dev/null || true
cp "$SRC/mlb_mm.duckdb"  "$WT/Answer Keys/" 2>/dev/null || true
cp "$SRC/mlb.duckdb"     "$WT/Answer Keys/" 2>/dev/null || true
cp "$SRC/pbp.duckdb"     "$WT/Answer Keys/" 2>/dev/null || true
cd "$WT" && Rscript "Answer Keys/MLB Dashboard/mlb_dashboard.R"
open "Answer Keys/MLB Dashboard/report.html"
```

Expected: the script prints its normal generation output and exits 0; `report.html`
is regenerated. In the browser (no scrolling yet), confirm: the `Placing on` pills
render with balances, clicking a pill changes the selected (blue) account, and the
Kelly Calculator still computes Risk/To-Win/EV/Kelly when you type Odds + Fair.
Nothing should be pinned yet.

- [ ] **Step 4: Commit**

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-dashboard-pin-kelly-book"
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "refactor(mlb-dashboard): relocate account row + wrap Kelly Calc for pinning"
```

---

## Task 2: Make the strip sticky (CSS + offset JS)

This task adds the actual freeze behavior. After it, the feature is complete.

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` (CSS block ~`:2376`, `:3270`; JS ~`:5498`)

- [ ] **Step 1: Add `.pinned-account-bar` CSS**

Locate the end of the `.header-row-accounts` rule (the block starting at `:2372`).
Immediately after its closing `}`, add:

```css
        /* Pinned scoring strip — account bar frozen at the top on every tab.
           Parent is the page-spanning .container, so the freeze holds through
           the whole card list. */
        .pinned-account-bar {
          position: sticky;
          top: 0;
          z-index: 40;              /* above cards; below filter-menu(100)/modal(999)/toast */
          background: #0d1117;       /* opaque so cards scroll cleanly underneath */
          border-bottom: 1px solid #21262d;
        }
```

- [ ] **Step 2: Add `.kelly-calc-pin` CSS**

Locate the end of the `.kelly-calc` rule (the block starting at `:3271`, closing
`}` at ~`:3291`). Immediately after it, add:

```css
        /* Full-width opaque sticky surface that carries the Kelly Calculator.
           Pins just below the account bar; --pin-top is set by syncPinOffset()
           to the account bar's measured height so it sits flush even when the
           pills wrap on narrow widths. Fallback 48px covers first paint. */
        .kelly-calc-pin {
          position: sticky;
          top: var(--pin-top, 48px);
          z-index: 30;              /* below the account bar, above cards */
          background: #0d1117;
          padding: 8px 0 2px 0;
        }
```

- [ ] **Step 3: Add the offset-measurement JS helper**

Find the end of `renderPills()` (its closing `}` at `:5374`, right before
`function refreshBalances()` at `:5376`). Insert the helper and a call so the
offset reflects the real rendered pill height. Replace:

```js
      bar.appendChild(pill);
    });
  }

  function refreshBalances() {
```

with:

```js
      bar.appendChild(pill);
    });
    syncPinOffset();   // pills just changed height; refresh the sticky offset
  }

  // Keep the Kelly Calc pinned flush below the account bar by measuring the
  // bar's real height into --pin-top. Recomputed on load, on resize, and after
  // pills (re)render (they are JS-populated, so height is only known then).
  function syncPinOffset() {
    var barEl = document.getElementById('wz-account-row');
    if (!barEl) return;
    document.documentElement.style.setProperty('--pin-top', barEl.offsetHeight + 'px');
  }
  window.addEventListener('load', syncPinOffset);
  window.addEventListener('resize', syncPinOffset);

  function refreshBalances() {
```

(If `renderPills`/`refreshBalances` have moved, anchor on the
`bar.appendChild(pill); });` that closes `renderPills`'s `forEach`, and on the
`function refreshBalances()` declaration.)

- [ ] **Step 4: Render and run the full verification**

Re-render (DBs already copied from Task 1):

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-dashboard-pin-kelly-book"
Rscript "Answer Keys/MLB Dashboard/mlb_dashboard.R"
open "Answer Keys/MLB Dashboard/report.html"
```

In the browser, confirm the spec's 7 checks:
1. On **Bets**, scroll into the cards: the pills bar is frozen at the very top and
   the Kelly Calculator is frozen directly below it — no gap, no overlap.
2. Type Odds/Fair in the Kelly Calc **while scrolled** → it still computes Risk /
   To Win / EV / Kelly (it reads the scrolled-away Bankroll + Kelly Fraction from
   the DOM, so the math is unaffected).
3. Click a pill **while scrolled** → account selection + balances still update,
   and the calc stays pinned.
4. Switch to **Parlays** / **Trifectas** → only the pills bar is frozen; the Kelly
   Calc is absent there; parlay placement still reads the selected account.
5. Resize the window narrow enough that the pills wrap onto two lines → the Kelly
   Calc re-pins flush below the now-taller bar (offset recomputed).
6. Cards scroll cleanly **under** both pinned elements (nothing peeks through
   beside the Kelly Calc strip).
7. Open a modal (e.g. a confirm dialog) / trigger a toast → both render **above**
   the pinned strip.

Expected: all 7 pass. If #1 shows a gap/overlap, `--pin-top` isn't matching the
bar height — confirm `syncPinOffset()` is being called after `renderPills` and on
load.

- [ ] **Step 5: Commit**

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "feat(mlb-dashboard): pin account pills + Kelly Calc while scrolling (sticky strip)"
```

---

## Task 3: Documentation

**Files:**
- Modify: `Answer Keys/MLB Dashboard/README.md`
- Modify: `Answer Keys/CLAUDE.md`

- [ ] **Step 1: Add a README note**

In `Answer Keys/MLB Dashboard/README.md`, under the bets-tab / odds-screen
section, add:

```markdown
### Pinned scoring strip

While scrolling the bet-card list, two controls stay frozen at the top:
- the **`Placing on` account pills** — frozen on every tab (Bets/Parlays/Trifectas).
- the **Kelly Calculator** — frozen directly below the pills, on the Bets tab only.

Always-on (no toggle). Implemented with CSS `position: sticky`: the account row
is a direct child of `.container` (a `.pinned-account-bar`), and the Kelly
Calculator is wrapped in `.kelly-calc-pin`. A small `syncPinOffset()` helper
measures the account bar's height into the `--pin-top` CSS variable so the calc
pins flush even when the pills wrap. The page title/Refresh row and the tab bar
are intentionally not pinned (scroll up to reach them).
```

- [ ] **Step 2: Add a CLAUDE.md maintenance note**

In `Answer Keys/CLAUDE.md`, in the MLB Dashboard area, add:

```markdown
### Pinned scoring strip (sticky account bar + Kelly Calc)

- The `Placing on` account row (`#wz-account-row`) was **relocated out of
  `.header`** to be a direct child of `.container` (class `pinned-account-bar`)
  so its `position: sticky; top: 0` freeze holds through the whole card list on
  every tab. Do not move it back inside `.header` — sticky only holds while its
  parent is on screen, and `.header` is too short.
- The Kelly Calculator is wrapped in `.kelly-calc-pin` (full-width opaque sticky
  surface) because `.kelly-calc` itself is `inline-flex; fit-content` and would
  let cards peek through beside it.
- `.kelly-calc-pin`'s `top` is `var(--pin-top)`, set by `syncPinOffset()` (called
  on load, on resize, and at the end of `renderPills`) to the account bar's
  measured height. If you change the account bar's height/padding, the offset
  self-corrects; no constant to update.
- z-index: account bar 40, Kelly wrapper 30 — above cards, below filter menu
  (100) / modal (999) / toast (1000/9999).
```

- [ ] **Step 3: Commit**

```bash
git add "Answer Keys/MLB Dashboard/README.md" "Answer Keys/CLAUDE.md"
git commit -m "docs(mlb-dashboard): document pinned scoring strip"
```

---

## Pre-merge cleanup & review

- [ ] **Remove the copied DBs from the worktree** (they are gitignored, so this is
  just disk hygiene — confirm `git status` is clean of any `.duckdb`):

```bash
cd "/Users/callancapitolo/NFLWork/.claude/worktrees/mlb-dashboard-pin-kelly-book"
rm -f "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb" \
      "Answer Keys/mlb_mm.duckdb" "Answer Keys/mlb.duckdb" "Answer Keys/pbp.duckdb"
rm -f "Answer Keys/MLB Dashboard/report.html"
git status --short   # expect: only the 3 commits' source changes, no data files
```

- [ ] **Executive engineer review** of `git diff main..HEAD` (per repo convention):
  confirm no DB/schema/server changes leaked in, ids/structure for the pill JS are
  intact, no dead CSS/JS, z-index values are below the modal/toast layer.

- [ ] **Get explicit user approval, then merge** to `main`. After merge: remove the
  worktree and delete the branch.

---

## Self-review (author check against the spec)

- **Spec coverage:** sticky strip (Task 2) ✓; minimal contents pills+calc, no
  bankroll/filters (no task adds them) ✓; pills all tabs / calc Bets-only (account
  bar is global sibling, `.kelly-calc-pin` lives in `#tab-bets`) ✓; relocate pills
  out of `.header` (Task 1 Step 1) ✓; always-on, no toggle/persistence (nothing
  added) ✓; offset JS (Task 2 Step 3) ✓; full-width opaque surfaces so cards don't
  peek (`.pinned-account-bar` bg, `.kelly-calc-pin` wrapper) ✓; test plan = spec's
  7 checks (Task 2 Step 4) ✓; docs (Task 3) ✓; version-control/worktree cleanup
  (pre-merge section) ✓.
- **Placeholder scan:** every code step shows the exact before/after text; no TBD.
- **Name consistency:** `pinned-account-bar`, `kelly-calc-pin`, `--pin-top`,
  `syncPinOffset`, `#wz-account-row`, `#wz-account-pills` used identically across
  all tasks and docs.
