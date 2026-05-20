# Kelly Calculator Redesign Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Visual redesign of the Kelly Calculator widget on the MLB Dashboard bets tab — tight inline strip sized to content (not full-width), uniform 15px SF Mono throughout, hierarchy via color + weight (not size), title "Kelly Calculator", same green left-bar accent as the bet-card hero strips. JS untouched.

**Architecture:** Server-render only — replace one CSS block and one HTML markup block in `Answer Keys/MLB Dashboard/mlb_dashboard.R`. The existing JS (`setupKellyCalc()` IIFE) keys off seven element IDs (`kc-odds`, `kc-fair`, `kc-risk`, `kc-towin`, `kc-edge`, `kc-kelly`, `kc-be`) plus three classes (`.invalid` on inputs, `.neg` on `#kc-risk`, `.pos`/`.neg` modifiers on `#kc-edge`'s `.chip` class). All IDs + class operations preserved verbatim in the redesigned HTML so no JS edit is needed.

**Tech Stack:** R 4.x, htmltools (existing). No new dependencies.

**Worktree:** This plan executes in `/Users/callancapitolo/NFLWork/.claude/worktrees/kelly-calc-redesign` on branch `worktree-kelly-calc-redesign`, already rebased onto current `main` (which contains the merged PR B Kelly Calculator). All paths below are relative to that worktree.

**Spec:** `docs/superpowers/specs/2026-05-14-kelly-calc-redesign-design.md`

---

## File Structure

```
Answer Keys/
└── MLB Dashboard/
    └── mlb_dashboard.R         (modify: 1 CSS block + 1 HTML block)
```

One file, two contiguous regions replaced:

- **CSS block** at `mlb_dashboard.R:2996-3095` — the `/* ====== KELLY CALCULATOR WIDGET ====== */` comment + all `.kelly-calc *` rules.
- **HTML block** at `mlb_dashboard.R:3180-3222` — the `tags$div(class = "kelly-calc", ...)` tagList (with its leading R comment).

The JS region (`setupKellyCalc()` IIFE at ~`mlb_dashboard.R:5308-5390`) is **not** touched.

---

## Task 1 — Replace the `.kelly-calc` CSS rules

**Why first.** CSS without the new HTML still works — the *old* HTML's classes (`.kelly-label`, `.kelly-fields`, `.kelly-out`, `.kc-input`, etc.) will simply lose their styling. The widget will render unstyled but functional, which is fine for a single uncommitted step. Task 2 brings the markup into alignment and commits both changes together.

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R:2996-3095` (CSS block — 100 lines total)

- [ ] **Step 1: Locate the CSS block to replace**

Run:
```bash
grep -n '====== KELLY CALCULATOR WIDGET ======\|.kelly-calc .detail .sep' "Answer Keys/MLB Dashboard/mlb_dashboard.R"
```

Expected output (line numbers may have drifted slightly):
```
2996:        /* ====== KELLY CALCULATOR WIDGET ====== */
3095:        .kelly-calc .detail .sep { color: #3a4658; margin: 0 6px; }
```

The block starts at the comment line and ends at the `.sep` rule (inclusive). 100 lines.

- [ ] **Step 2: Replace the CSS block**

In `Answer Keys/MLB Dashboard/mlb_dashboard.R`, replace the entire block from the `/* ====== KELLY CALCULATOR WIDGET ====== */` comment through (and including) the `.kelly-calc .detail .sep` rule with this new block. The indentation (8 spaces) must match the surrounding lines so the embedded CSS stays inside the `tags$style(HTML('...'))` block.

```r
        /* ====== KELLY CALCULATOR WIDGET ====== */
        /* Tight inline strip sized to content. Uniform 15px SF Mono;     */
        /* hierarchy from color + weight only. Same green-bar accent as   */
        /* the V8 bet-card hero strips. JS keys off seven element IDs     */
        /* (kc-odds, kc-fair, kc-risk, kc-towin, kc-edge, kc-kelly, kc-be)*/
        /* and the .invalid / .risk.neg / .chip.pos / .chip.neg classes.  */
        .kelly-calc {
          display: inline-flex;
          align-items: center;
          gap: 20px;
          padding: 14px 22px 14px 20px;
          background: linear-gradient(135deg,
            rgba(63, 185, 80, 0.10) 0%,
            rgba(63, 185, 80, 0.03) 50%,
            #161b22 100%);
          border: 1px solid #2a3442;
          border-left: 4px solid #3fb950;
          border-radius: 10px;
          box-shadow: inset 1px 0 0 rgba(63, 185, 80, 0.18);
          max-width: fit-content;
          color: #7d8590;
          font-family: "SF Mono", SFMono-Regular, Menlo, Consolas, monospace;
          font-size: 15px;
          font-variant-numeric: tabular-nums;
          line-height: 1.2;
          margin: 0 0 14px 0;
        }
        .kelly-calc * {
          font-family: inherit;
          font-size: inherit;
          font-variant-numeric: tabular-nums;
          line-height: 1.2;
        }
        .kelly-calc .title {
          color: #e6edf3;
          font-weight: 600;
          padding-right: 16px;
          border-right: 1px solid #232b36;
        }
        .kelly-calc .pair {
          display: flex; align-items: center; gap: 10px;
        }
        .kelly-calc .pair label {
          color: #7d8590; font-weight: 500;
        }
        .kelly-calc .pair input {
          background: #0d1117;
          border: 1px solid #2a3442;
          color: #e6edf3;
          font-weight: 500;
          padding: 7px 13px;
          border-radius: 6px;
          width: 102px;
          text-align: center;
          transition: border-color 0.12s ease, box-shadow 0.12s ease;
        }
        .kelly-calc .pair input::placeholder { color: #484f58; }
        .kelly-calc .pair input:focus {
          outline: none;
          border-color: #3fb950;
          box-shadow: 0 0 0 3px rgba(63, 185, 80, 0.18);
        }
        .kelly-calc .pair input.invalid {
          border-color: #f85149;
          box-shadow: 0 0 0 3px rgba(248, 81, 73, 0.12);
        }
        .kelly-calc .arrow { color: #3a4658; }
        .kelly-calc .risk-grp {
          display: flex; align-items: baseline; gap: 9px;
          padding-left: 18px; border-left: 1px solid #232b36;
        }
        .kelly-calc .risk-grp .lbl   { color: #7d8590; font-weight: 500; }
        .kelly-calc .risk-grp .risk  { color: #3fb950; font-weight: 700; }
        .kelly-calc .risk-grp .risk.neg { color: #7d8590; }
        .kelly-calc .risk-grp .towin { color: #7d8590; }
        .kelly-calc .metrics {
          display: flex; gap: 16px;
          padding-left: 18px; border-left: 1px solid #232b36;
        }
        .kelly-calc .metrics .grp {
          display: flex; align-items: baseline; gap: 7px;
        }
        .kelly-calc .metrics .k   { color: #7d8590; font-weight: 500; }
        .kelly-calc .metrics .val { color: #c9d1d9; font-weight: 600; }
        .kelly-calc .metrics .chip {
          color: #5d6470;
          font-weight: 600;
        }
        .kelly-calc .metrics .chip.pos {
          color: #3fb950;
          background: rgba(63, 185, 80, 0.10);
          padding: 2px 9px;
          border-radius: 4px;
        }
        .kelly-calc .metrics .chip.neg {
          color: #f85149;
          background: rgba(248, 81, 73, 0.10);
          padding: 2px 9px;
          border-radius: 4px;
        }
```

> **Indentation note:** every line of CSS inside the `tags$style(HTML('...'))` block in `mlb_dashboard.R` is indented with **8 spaces** (matching the surrounding R `tags$style(...)` wrapper). When you paste the new block, preserve that 8-space indent. If your editor auto-indents differently, fix it before saving.

> **Quote choice:** the CSS uses **double quotes** for `font-family: "SF Mono", ...` because the outer R wrapper is `HTML('...')` (single-quoted). A single-quoted `'SF Mono'` inside would terminate the R string early. This matches the established pattern in this file (the existing `.kelly-calc` block used the same trick).

- [ ] **Step 3: Run a parse check**

Run from the worktree root:
```bash
Rscript -e 'parse("Answer Keys/MLB Dashboard/mlb_dashboard.R"); cat("PARSE OK\n")'
```

Expected: `PARSE OK`. If anything else (e.g. `Error: ... unexpected ...`), you broke a string literal — most likely an unescaped `'` or a mismatched indent inside the `HTML('...')` block. Fix and re-run before continuing.

> **Do NOT commit yet.** The new CSS still expects new class names (`.title`, `.pair`, `.risk-grp`, `.metrics`, `.grp`, `.k`, `.val`, `.chip`) that the existing HTML doesn't emit. Task 2 fixes the HTML and commits both changes together.

---

## Task 2 — Replace the `.kelly-calc` HTML block + verify + commit

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R:3180-3222` (HTML tagList block — ~43 lines)

- [ ] **Step 1: Locate the HTML block to replace**

Run:
```bash
grep -n 'Kelly Calculator widget — manual no-vig\|tags\$div(class = "kelly-calc"\|class = "detail"' "Answer Keys/MLB Dashboard/mlb_dashboard.R" | head -10
```

Expected output (line numbers may have drifted by a few):
```
3180:          # Kelly Calculator widget — manual no-vig Kelly sizing tool.
3185:          tags$div(class = "kelly-calc",
3211:              tags$div(class = "detail",
```

The block starts at the leading R comment (`# Kelly Calculator widget — manual ...`) on line 3180 and ends at the closing `),` of the `tags$div(class = "kelly-calc", ...)` block (currently line 3222). The next sibling in the tagList is `# Filter Bar` (currently line 3224 or 3225).

- [ ] **Step 2: Replace the HTML block**

In `Answer Keys/MLB Dashboard/mlb_dashboard.R`, replace the block from the `# Kelly Calculator widget — manual no-vig Kelly sizing tool.` comment line through (and including) the closing `),` of the `tags$div(class = "kelly-calc", ...)` block with this new block. The leading indent must be **10 spaces** (matching the surrounding `tagList` indent — the existing code is at 10 spaces too).

```r
          # Kelly Calculator widget — manual no-vig Kelly sizing tool.
          # Tight inline strip sized to content. Reads dashboard Bankroll
          # + Kelly Fraction live; computes risk / to-win / edge / Kelly%
          # / breakeven from two text inputs. Element IDs and class
          # operations (.invalid / .neg / .pos) preserved so the
          # setupKellyCalc() JS handler below runs unchanged.
          tags$div(class = "kelly-calc",
            tags$span(class = "title", "Kelly Calculator"),
            tags$div(class = "pair",
              tags$label("Odds"),
              tags$input(id = "kc-odds", type = "text", placeholder = "+120")
            ),
            tags$span(class = "arrow", HTML("&rarr;")),
            tags$div(class = "pair",
              tags$label("Fair"),
              tags$input(id = "kc-fair", type = "text", placeholder = "52.5%")
            ),
            tags$div(class = "risk-grp",
              tags$span(class = "lbl", "Risk"),
              tags$span(id = "kc-risk", class = "risk", "$0"),
              tags$span(class = "towin", "/ to win "),
              tags$span(id = "kc-towin", class = "towin", "$0")
            ),
            tags$div(class = "metrics",
              tags$span(class = "grp",
                tags$span(class = "k", "Edge"),
                tags$span(id = "kc-edge", class = "chip", HTML("&mdash;"))
              ),
              tags$span(class = "grp",
                tags$span(class = "k", "Kelly"),
                tags$span(id = "kc-kelly", class = "val", HTML("&mdash;"))
              ),
              tags$span(class = "grp",
                tags$span(class = "k", "BE"),
                tags$span(id = "kc-be", class = "val", HTML("&mdash;"))
              )
            )
          ),
```

> **Trailing comma:** the block ends with `),` (close-paren plus comma). The comma is intentional — the next sibling in the `tagList` is the filter bar. Without the comma, the file will parse but the filter bar will be passed as a positional argument to `tags$div(class = "kelly-calc", ...)` instead of being its own sibling.

> **`HTML("&rarr;")` and `HTML("&mdash;")`:** the `HTML()` wrapper tells `htmltools` not to escape the entity. Without it, the literal text `&rarr;` would render in the browser as `&amp;rarr;`. Same pattern is used elsewhere in this file for `&#x21bb;` (the refresh icon).

- [ ] **Step 3: Parse check**

Run:
```bash
Rscript -e 'parse("Answer Keys/MLB Dashboard/mlb_dashboard.R"); cat("PARSE OK\n")'
```

Expected: `PARSE OK`. Most likely failure: missed the trailing `),` or a typo in `tags$span(...)`. Fix and re-run.

- [ ] **Step 4: Smoke render check**

Run from the worktree root:
```bash
Rscript "Answer Keys/MLB Dashboard/mlb_dashboard.R" 2>&1 | tail -5
```

Expected last 1–2 lines:
```
Dashboard saved to: /Users/callancapitolo/NFLWork/.claude/worktrees/kelly-calc-redesign/Answer Keys/MLB Dashboard/report.html
Open http://localhost:8083 to view
```

A clean exit. Any earlier `Error` or `dyn.load` failure means the R-side rendering broke — investigate before continuing.

- [ ] **Step 5: Verify the rendered HTML contains the new markup**

The rendered file is HTML on one very long line, so `grep -c` counts lines (always 0 or 1). Use `grep -ao ... | wc -l` for actual occurrence counts.

Run:
```bash
REPORT="Answer Keys/MLB Dashboard/report.html"
for kw in 'kelly-calc' 'class="title">Kelly Calculator' 'id="kc-odds"' 'id="kc-fair"' 'id="kc-risk"' 'id="kc-towin"' 'id="kc-edge"' 'id="kc-kelly"' 'id="kc-be"' 'class="risk-grp"' 'class="metrics"' 'class="arrow"'; do
  count=$(grep -ao "$kw" "$REPORT" | wc -l | tr -d ' ')
  printf "%-40s %s\n" "$kw" "$count"
done
```

Expected: every count ≥ 1. If any is 0, the HTML block didn't replace cleanly — diff vs Step 2 and try again. Also confirm that none of the **old** class names still appear:

```bash
for kw in 'class="kelly-label"' 'class="kelly-fields"' 'class="kelly-out"' 'class="kc-input"' 'class="field"' 'class="risk-row"' 'class="risk-label"' 'class="risk-value"' 'class="detail"' 'class="sep"'; do
  count=$(grep -ao "$kw" "$REPORT" | wc -l | tr -d ' ')
  printf "%-30s %s\n" "$kw" "$count"
done
```

Expected: every count = 0. If any is ≥ 1, an old reference survived somewhere — most likely the comment block or an accidental partial replacement. Find and remove it.

- [ ] **Step 6: Manual browser sanity check (visual contract)**

> **Note for subagent implementers:** the live-interaction checks below (typing into inputs, clicking around) require either a human or computer-use tooling to perform. If you don't have either, do the static checks you *can* do (Steps 1–5 above) and explicitly hand the live-interaction checks off to the user in your final report. Do NOT claim "verified" if you haven't observed the browser.

Open the rendered report in a browser:
```bash
open "Answer Keys/MLB Dashboard/report.html"
```

On the bets tab, locate the Kelly Calculator strip (directly under the Bankroll / Kelly Fraction row). Verify all five state contracts from the spec:

1. **Empty state** (initial render): grey `$0` for Risk, em-dash for Edge (dim grey), em-dash for Kelly (near-white), em-dash for BE (near-white). Strip sized to content, not stretched full-width. Green 4px left border, subtle green inner glow.
2. **Live positive-EV state**: type `+120` in Odds and `52.5%` in Fair. Risk goes green (e.g. `$76`); "to win" shows a real dollar value; Edge chip turns green-filled with `+X.X%`; Kelly and BE show real percentages.
3. **Live negative-EV state**: change Fair to `40%`. Risk greys back to `$0`; Edge chip turns red-filled with a negative percentage.
4. **Invalid input state**: type garbage (e.g. `abc`) in Odds. The Odds input border turns red (the `.invalid` class).
5. **Live coupling**: change the Bankroll input (above the widget) — Risk should recompute without re-typing Odds/Fair.

Also confirm:
- The font is SF Mono everywhere (labels, values, title, "to win").
- All text appears at the same size; no large Risk number jumping above the rest.
- No JavaScript console errors (open browser dev tools to check).
- The rest of the bets tab (bet cards, per-cell RAW/FAIR toggle, parlay tab) still renders normally.

If anything fails, fix it and re-run from Step 3.

**Static fallback** (if you can't observe the browser): Steps 1–5 above already confirmed the new markup is in the rendered HTML and the old markup is gone. The JS contract preservation is also static-verifiable — `grep -ao 'id="kc-' "Answer Keys/MLB Dashboard/report.html" | sort -u` should list exactly the seven IDs (`kc-be`, `kc-edge`, `kc-fair`, `kc-kelly`, `kc-odds`, `kc-risk`, `kc-towin`). If so, the JS will function correctly because every reference it makes resolves. The remaining unknowns (CSS visual rendering, click-and-type behavior) require human eyes.

- [ ] **Step 7: Commit**

Run from the worktree root:
```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R" && \
  git commit -m "$(cat <<'EOF'
feat(bets-tab): redesign Kelly Calculator widget — tight inline strip

Replaces the full-width grid layout that left a large dead zone in the
middle with a tight inline-flex strip sized to content. Uniform 15px
SF Mono throughout — hierarchy now comes from color and weight, not
size jumps. Title reads "Kelly Calculator" (title case). Same 4px green
left-bar accent and gradient as the V8 bet-card hero strips, so the
widget reads as part of the same visual family.

Element IDs and class operations (.invalid on inputs, .neg on #kc-risk,
.pos/.neg on #kc-edge.chip) preserved verbatim, so the setupKellyCalc()
JS handler runs unchanged.

Spec: docs/superpowers/specs/2026-05-14-kelly-calc-redesign-design.md

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

Expected: a single new commit on `worktree-kelly-calc-redesign` touching only `Answer Keys/MLB Dashboard/mlb_dashboard.R`.

---

## Task 3 — Final verification + handoff

**Files:** none (verification only).

- [ ] **Step 1: Branch summary**

Run:
```bash
git log --oneline main..HEAD && \
  echo "---DIFF STAT---" && \
  git diff --stat main..HEAD
```

Expected output:
```
<sha> feat(bets-tab): redesign Kelly Calculator widget — tight inline strip
<sha> docs(spec): Kelly Calculator visual redesign
---DIFF STAT---
 Answer Keys/MLB Dashboard/mlb_dashboard.R          | <±N> +/-/-
 docs/superpowers/specs/2026-05-14-kelly-calc-redesign-design.md | +342 ++
 2 files changed, ...
```

Two commits ahead of `main`: the spec (already on this branch from earlier) and the new feature commit. One source file changed.

- [ ] **Step 2: Re-run the existing R test suite to confirm zero regression**

The redesign is HTML + CSS only, but the dashboard's R-side tests should still pass cleanly. Run:

```bash
cd "Answer Keys" && \
  for t in tests/test_book_cell.R tests/test_devig_pair_matches_tools.R tests/test_odds_screen.R; do
    echo "=== $t ==="
    Rscript -e "testthat::test_file('$t')" 2>&1 | grep -E "^\[" | tail -1
  done
```

Expected: every test file reports `FAIL 0`. The pass counts should match what `main` produces (currently 21 / 12 / 75). Any new failures mean the R rewrite accidentally broke something — investigate before handoff.

- [ ] **Step 3: Hand off to the user**

Report back with:
- Branch name (`worktree-kelly-calc-redesign`) and head SHA.
- Diff stat against `main`.
- Test pass counts.
- Confirmation that all four state contracts (empty / positive-EV / negative-EV / invalid input) and live bankroll coupling rendered correctly in the browser.
- Open question: **does the user want to merge to `main` now**, or first review the live dashboard themselves? Per project policy, **do NOT merge to `main` without explicit user approval**.

---

## Out-of-band notes

- **No tests** were added or modified by this plan. The redesign is visual; the existing JS computation is unchanged; there is no per-test contract to assert that's different from "the dashboard renders cleanly". Smoke render + manual browser check are sufficient.
- **No README / CLAUDE.md edits** are strictly required. The behavior is unchanged, so the existing `Answer Keys/CLAUDE.md` "Per-cell devig toggle + Kelly Calculator widget (PR B, 2026-05)" subsection still accurately describes the widget. An optional one-line "visual refresh 2026-05-15" note could go in that subsection if the user wants a future grep landmark — flag this in handoff but don't add it unsolicited.
- **Worktree cleanup** happens after merge approval, not as part of this plan. The implementer's job ends at "ready to merge".
