# Kelly Calculator — Visual Redesign Spec

> Scope: visual-only redesign of the Kelly Calculator widget on the MLB Dashboard bets tab (`Answer Keys/MLB Dashboard/mlb_dashboard.R`). HTML markup + CSS only. JS logic and computation untouched.

---

## Review Pack

**What we're building** — A visual redesign of the Kelly Calculator widget that sits below the Bankroll/Kelly settings strip on the MLB Dashboard bets tab. Same inputs (Odds + Your Fair), same outputs (Risk + To Win + Edge + Kelly% + Breakeven), same JS computation. New shape: a tight inline strip sized to its content (not stretched full-width), uniform 15px SF Mono throughout, hierarchy via color and weight (not size), one bigger green left-bar accent matching the bet-card hero treatment.

**Key decisions:**

1. **Inline-flex strip, sized to content** — not full-width. The current widget uses a 3-zone grid `auto 1fr auto` that stretches across the page and leaves a giant dead-zone between inputs and outputs. The redesign uses `display: inline-flex; max-width: fit-content`, so the strip is only as wide as needed. Alternative rejected: re-balance the full-width grid. That still feels like a status strip, not a tool.
2. **Uniform 15px font size across every text element** — hierarchy comes from color (greys / whites / greens / red) and weight (500 / 600 / 700). Alternative rejected: size hierarchy with a 24px Risk number. User specifically prefers single-size for visual consistency.
3. **Title is "Kelly Calculator" (title case), no subtitle** — drops the "manual sizing tool" subtitle from the V8 design. At 15px uniform sizing it added caption noise without informing.
4. **HTML restructure required, not just CSS** — current DOM is `.kelly-calc > {.kelly-label, .kelly-fields, .kelly-out}` with metrics nested *inside* `.kelly-out`. The new flat inline-flex needs a flatter DOM: `title · pair · arrow · pair · risk-grp · metrics` as direct children. So this is markup + CSS rewrite, not CSS-only.
5. **JS untouched** — the `setupKellyCalc()` IIFE references seven element IDs (`kc-odds`, `kc-fair`, `kc-risk`, `kc-towin`, `kc-edge`, `kc-kelly`, `kc-be`) plus the `.invalid` class on inputs and the `.neg` / `.pos` / `.chip` classes on `#kc-edge`. All preserved exactly in the new HTML; the JS file region runs unchanged.

**Risks / push back here:**

- **Removing the "manual sizing tool" subtitle.** Choice not constraint — if you want a hover tooltip or a small `?` info icon explaining what the widget is for, that's a different design. Tell me if you want it back.
- **Inline-flex `max-width: fit-content` left-aligns and stops at content width.** If you'd prefer it centered horizontally in the page or stretching to a max-width but not full-width, that's a different CSS choice we can make now.
- **Section divider style** (1px solid `#232b36`). If the dividers read too dark in the live dashboard, we can switch to dashed, lighter grey, or remove them and rely on padding alone.
- **Removing the existing chip-style for the empty-state "—"**. The current widget shows the dash inside the green chip (`<span id="kc-edge" class="chip">—</span>`). In the redesign, the chip styling (background-fill) only applies when `.pos` or `.neg` is set; the bare `.chip` with no modifier renders as plain text. This is a small JS-CSS contract change — the JS already sets the modifier class on real values, so empty state degrades cleanly.

**Worth understanding** (1 concept):

- **CSS specificity and "single font size everywhere"** — the trick is `.kelly-calc { font-size: 15px }` on the parent and `.kelly-calc * { font-size: inherit }` on the universal child selector. The `*` selector with `inherit` overrides default browser styling (notably `<input>`, which ships with its own smaller default font-size). This is similar to R's `inherit.aes = TRUE` in ggplot2: children take the parent's setting unless they override it explicitly.

---

## Goal

Rebuild the Kelly Calculator widget on the MLB Dashboard bets tab so it:

1. Reads as a deliberate, scannable tool rather than a stretched-out status strip with empty space.
2. Uses one font, one size, one rhythm — hierarchy from color and weight alone.
3. Matches the visual family of the V8 bet-card hero strips (green left-bar accent + green-tinted gradient).
4. Preserves all existing functionality: same inputs, same outputs, same JS, same dashboard inputs (`bankroll-input`, `kelly-input`) live-coupling, same validation feedback (red border on invalid input), same negative-EV greying.

## Non-goals

- **No JS changes.** `setupKellyCalc()` and `recompute()` run unchanged.
- **No new dependencies** (no Tailwind, no CSS framework, no icon font).
- **No persistence.** Inputs still reset on page reload.
- **No re-positioning** of the widget relative to the page. Stays directly below `.sizing-controls`, directly above the filter bar.
- **No layout change to the dashboard's `bankroll-input` / `kelly-input` controls.** They stay as they are.

## Scope

**Files affected:**

- `Answer Keys/MLB Dashboard/mlb_dashboard.R` — replace the `.kelly-calc` HTML block (currently at lines 3180-3222) and the `.kelly-calc *` CSS rules (currently at lines 2997-3094).
- `Answer Keys/CLAUDE.md` — minor description tweak in the "Per-cell devig toggle + Kelly Calculator widget (PR B, 2026-05)" subsection if its description is now stale.

**Files NOT affected:**

- The `setupKellyCalc()` IIFE in `mlb_dashboard.R` (~lines 5308-5390) — not touched.
- All other dashboard files.
- All tests — there are no Kelly-Calc-specific tests today, and the redesign is visual, so smoke render + manual browser verification is sufficient.

## Visual Specification

### Container

```css
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
  border-left: 4px solid #3fb950;            /* signature green bar */
  border-radius: 10px;

  box-shadow: inset 1px 0 0 rgba(63, 185, 80, 0.18);  /* subtle inner glow */

  max-width: fit-content;                    /* don't stretch full-width */

  color: #7d8590;                            /* default descendant text */

  font-family: 'SF Mono', SFMono-Regular, Menlo, Consolas, monospace;
  font-size: 15px;
  font-variant-numeric: tabular-nums;
  line-height: 1.2;

  margin: 0 0 14px 0;                        /* preserve current vertical spacing */
}

.kelly-calc * {
  font-family: inherit;
  font-size: inherit;                        /* force inputs to also be 15px */
  font-variant-numeric: tabular-nums;
  line-height: 1.2;
}
```

### Color tokens (reused from V8 palette)

| Token | Hex | Used for |
|---|---|---|
| Bright white | `#e6edf3` | Title |
| Near-white | `#c9d1d9` | Input values, Kelly/BE metric values |
| Medium grey | `#7d8590` | Labels, "to win", neutral-state Risk |
| Dim grey | `#5d6470` | Empty-state em-dashes |
| Green accent | `#3fb950` | Risk (live), positive Edge chip, left border |
| Red accent | `#f85149` | Negative Edge chip |
| Divider | `#232b36` | 1px vertical section dividers |
| Borders | `#2a3442`, `#3a4658` | Input borders, arrow |
| Surfaces | `#161b22`, `#0d1117` | Gradient stop, input background |

### Element rules

**Title:**
```css
.kelly-calc .title {
  color: #e6edf3;
  font-weight: 600;
  padding-right: 16px;
  border-right: 1px solid #232b36;
}
```

**Input pair (label + input):**
```css
.kelly-calc .pair {
  display: flex;
  align-items: center;
  gap: 10px;
}

.kelly-calc .pair label {
  color: #7d8590;
  font-weight: 500;
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
```

**Arrow (the visual separator between input pairs):**
```css
.kelly-calc .arrow { color: #3a4658; }
```

**Risk group:**
```css
.kelly-calc .risk-grp {
  display: flex;
  align-items: baseline;
  gap: 9px;
  padding-left: 18px;
  border-left: 1px solid #232b36;
}

.kelly-calc .risk-grp .lbl { color: #7d8590; font-weight: 500; }
.kelly-calc .risk-grp .risk { color: #3fb950; font-weight: 700; }
.kelly-calc .risk-grp .risk.neg { color: #7d8590; }   /* JS toggles .neg on negative-EV */
.kelly-calc .risk-grp .towin { color: #7d8590; }
```

**Metrics group (Edge / Kelly / BE):**
```css
.kelly-calc .metrics {
  display: flex;
  gap: 16px;
  padding-left: 18px;
  border-left: 1px solid #232b36;
}

.kelly-calc .metrics .grp {
  display: flex;
  align-items: baseline;
  gap: 7px;
}

.kelly-calc .metrics .k { color: #7d8590; font-weight: 500; }
.kelly-calc .metrics .val { color: #c9d1d9; font-weight: 600; }

.kelly-calc .metrics .chip {
  color: #5d6470;          /* dim grey for empty-state em-dash */
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

### HTML structure (R `tags$div` form, replaces lines 3185-3222)

```r
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
)
```

### State variations (visual contract)

| State | Risk color | Risk class toggle | Edge chip class | Kelly / BE metrics |
|---|---|---|---|---|
| Empty (initial render) | grey `#7d8590` | `.risk` (no toggle) | `.chip` plain → em-dash in dim grey `#5d6470` | em-dash in `.val` near-white `#c9d1d9` |
| Live, positive EV | green `#3fb950` | `.risk` (no toggle) | `.chip.pos` (green fill, green text) | live values in near-white |
| Live, negative EV | grey `#7d8590` | `.risk.neg` | `.chip.neg` (red fill, red text) | live values in near-white |
| Invalid input | (unchanged) | (unchanged) | (unchanged) | input gets `.invalid` red border |

> **Minor cosmetic note (accepted trade-off):** in the empty state, the Edge em-dash renders in dim grey (`#5d6470`) while the Kelly/BE em-dashes render in near-white (`#c9d1d9`). This is because the JS toggles a modifier class on `#kc-edge` (`.chip` → `.chip.pos`/`.chip.neg`) but only updates the `textContent` of `#kc-kelly` and `#kc-be` without touching their class. Making all three em-dashes match would require a small JS edit; we accept the inconsistency to keep JS untouched. Users rarely see the empty state — most renders carry live values.

### JS contract (preserved verbatim)

The existing `setupKellyCalc()` IIFE at `mlb_dashboard.R:5308-5390` already:

- Reads `#kc-odds.value` and `#kc-fair.value`.
- Writes `#kc-risk.textContent`, `#kc-towin.textContent`, `#kc-edge.textContent`, `#kc-kelly.textContent`, `#kc-be.textContent`.
- Toggles `.invalid` on `#kc-odds` and `#kc-fair` based on parse success.
- Toggles `.neg` on `#kc-risk` when EV ≤ 0.
- Sets `#kc-edge.className` to `'chip'` or `'chip pos'` or `'chip neg'`.
- Reads `#bankroll-input.value` and `#kelly-input.value` and re-runs on their `input` event.

Every ID and every class operation in the redesigned HTML matches one-to-one. The new CSS has rules for `.invalid`, `.risk.neg`, `.chip.pos`, `.chip.neg`. No JS edit required.

## Out of scope (explicitly)

- **Tooltips / help text** — if we want the "manual sizing tool" subtitle back in some form, that's a separate design decision.
- **Keyboard shortcuts** (e.g., tab order, arrow keys to step) — existing tab order is fine.
- **Validation messages** — `.invalid` red border is the only feedback; no inline error text.
- **Mobile / narrow viewport** — the dashboard is desktop-only; we don't redesign for narrow screens.
- **Theming / dark-light mode** — dashboard is dark-only.

## Testing

No new logic; no new tests. Verification is:

1. **Parse** — `Rscript -e 'parse("Answer Keys/MLB Dashboard/mlb_dashboard.R"); cat("PARSE OK\n")'`.
2. **Smoke render** — `Rscript "Answer Keys/MLB Dashboard/mlb_dashboard.R"`. Grep the resulting `report.html` for: `kelly-calc`, `kc-odds`, `kc-fair`, `kc-risk`, `kc-edge` — all present.
3. **Manual browser check** — open the rendered dashboard in a browser:
   - Confirm empty state renders (grey `$0`, em-dashes for Edge/Kelly/BE).
   - Type `+120` in Odds and `52.5%` in Fair — confirm Risk goes green with a real dollar value, Edge chip turns green with a `+X.X%` value.
   - Type `40%` in Fair — confirm Risk greys to `$0` and Edge chip turns red.
   - Type garbage in Odds — confirm input border turns red.
   - Tweak the dashboard's Bankroll input — confirm Risk recomputes live.
4. **Regression** — confirm the rest of the bets tab still renders normally and the per-cell devig toggle still works (PR B feature unaffected).

## Rollback

Single-commit revert. The redesign is HTML + CSS in one file, no migrations, no state changes. If something looks wrong post-merge, `git revert <commit-sha>` restores the prior widget.

## Implementation order

The implementation plan (next step) will sequence:

1. CSS — replace the `.kelly-calc *` rule block.
2. HTML — replace the `tags$div(class = "kelly-calc", ...)` block.
3. Parse + smoke render verification.
4. Update `Answer Keys/CLAUDE.md` if the description is stale.
5. Single commit.

## Version control

- **Worktree**: `worktree-kelly-calc-redesign` at `/Users/callancapitolo/NFLWork/.claude/worktrees/kelly-calc-redesign`, already rebased onto current `main`. ✓
- **Branch**: `worktree-kelly-calc-redesign`.
- **Commits expected**: 1 (small visual change, one file).
- **Merge plan**: after user approval, merge `--no-ff` to local `main` (matching project pattern). No push to remote without explicit user approval.
- **Worktree cleanup**: remove the worktree + delete the branch after merge.

## Documentation

- `Answer Keys/CLAUDE.md` — review the "Per-cell devig toggle + Kelly Calculator widget (PR B, 2026-05)" subsection. The description ("A 'Kelly Calculator' strip sits below the Bankroll/Kelly Fraction settings on the bets tab. Manual sizing tool: type Odds + Fair, get recommended Risk and To Win.") is still accurate after the redesign — no edit strictly required. If we want to record the redesign date for future archaeology, add a one-line note.
- `Answer Keys/MLB Dashboard/README.md` — the "Per-cell devig toggle" bullet describing the Kelly Calculator stays accurate.

No README/CLAUDE.md changes are strictly required by the redesign; both files described the *behavior*, which is unchanged. Optional one-line "visual refresh 2026-05-14" note in CLAUDE.md if useful for future grep.
