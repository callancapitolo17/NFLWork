# MLB Parlay Dashboard — Books Strip Redesign

**Date:** 2026-04-25
**Scope:** Parlay tab of the MLB dashboard only
**Files touched:** `Answer Keys/MLB Dashboard/mlb_dashboard.R`, `Answer Keys/MLB Dashboard/README.md`
**Out of scope:** Singles tab, R pricer, server, scrapers, schema, CBB dashboard

## Problem

The Parlay Opportunities table on the MLB dashboard has 14 visible columns. On phones and laptops the table is hard to read — too many narrow columns competing for horizontal space, requiring sideways scrolling. Separately, the dashboard does not surface our own model's joint probability for each parlay, even though it is already computed and stored — only DK / FD / PX / NV per-book devigged probabilities are visible.

## Goals

1. **Surface our model's joint probability** alongside the four book probabilities so we can see model-vs-market disagreement.
2. **Reduce column count** so the table is readable on a phone or laptop without horizontal scrolling.
3. **Preserve all existing information** — no per-book probability is hidden behind a click or a toggle.
4. **Zero impact** on the R pricer, the data schema, the dashboard server, the singles tab, or the bet placement / CLV flow.

## Solution

Replace the five separate per-book probability columns (Model + DK + FD + PX + NV) with **one combined cell per row** containing an inline pill row, plus a **consensus pill** (model + books blended) at the end:

```
[ M 27.4 ] [ DK 26.9 ] [ FD 27.1 ] [ PX 27.8 ] [ NV 26.7 ] [ Cons 27.2 ]
```

- **Model pill (M)** — green tint background, brighter text. Visually flags "this is ours."
- **Book pills (DK / FD / PX / NV)** — neutral grey, monospace.
- **Consensus pill (Cons)** — blue tint with a thin left border. Equals the mean of the available inputs (model + books); identical to the probability that already drives the existing "Fair" American-odds column.

The strip uses CSS `flex-wrap: wrap` so the six pills sit on one line on wide screens and break to two lines on narrow ones (no horizontal scroll, no information hidden).

## Final Column Layout

The table goes from 14 visible columns → 10:

| # | Column      | Min width | Notes                                                         |
|---|-------------|-----------|---------------------------------------------------------------|
| 1 | Game        | 180px     | Matchup on first line; game time on second line (Time merged) |
| 2 | Legs        | 260px     | Unchanged                                                     |
| 3 | **Books**   | 320px     | New combined pill row: M / DK / FD / PX / NV / Cons           |
| 4 | Fair        | 70px      | American odds, unchanged                                      |
| 5 | WZ          | 70px      | Unchanged                                                     |
| 6 | Corr        | 65px      | Unchanged (hidden via media query at <700px width)            |
| 7 | Edge %      | 70px      | Unchanged                                                     |
| 8 | Size        | 65px      | Unchanged                                                     |
| 9 | To Win      | 65px      | Unchanged                                                     |
| 10 | Action     | 90px      | Unchanged                                                     |

Columns removed from view: separate `dk_fair_prob`, `fd_fair_prob`, `px_fair_prob`, `nv_fair_prob` (folded into Books); standalone `Time` (folded into Game). The underlying data columns remain in the dataframe — only their `colDef` rendering changes.

## Data Flow

No R-pricer changes. Every value the new cell needs is already written to `mlb_parlay_opportunities` by `mlb_correlated_parlay.R`:

| Pill   | Source column        | Currently in dashboard?                |
|--------|----------------------|----------------------------------------|
| M      | `model_prob_raw`     | Stored, hidden (`colDef(show = FALSE)`) |
| DK     | `dk_fair_prob`       | Visible (own column)                    |
| FD     | `fd_fair_prob`       | Visible (own column)                    |
| PX     | `px_fair_prob`       | Visible (own column)                    |
| NV     | `nv_fair_prob`       | Visible (own column)                    |
| Cons   | `blended_prob_raw`   | Stored, surfaced today only as American odds in the Fair column |

All values are already devigged probabilities on a `[0, 1]` scale; they are formatted as one-decimal percentages in the pill renderer.

## Pill Cell Renderer

The Books column uses an HTML cell renderer (`html = TRUE` in reactable) emitting:

```html
<span class="books-strip">
  <span class="pill model">M  27.4</span>
  <span class="pill book">DK 26.9</span>
  <span class="pill book">FD 27.1</span>
  <span class="pill book">PX 27.8</span>
  <span class="pill book">NV 26.7</span>
  <span class="pill cons">Cons 27.2</span>
</span>
```

Supporting CSS goes into the dashboard's existing inline `<style>` block:

- `.books-strip` — `display: flex; flex-wrap: wrap; gap: 6px;`
- `.pill` — `padding: 2px 6px; border-radius: 3px; font-family: monospace; font-size: 10px; color: #8b949e; background: #21262d;`
- `.pill.model` — green tint (`background: #1f3a2c; color: #7ee787`)
- `.pill.book` — default grey (no override)
- `.pill.cons` — blue tint with left accent (`background: #1c2738; border-left: 2px solid #58a6ff; color: #79c0ff`)
- `.pill.dim` — dimmed for NA values (`opacity: 0.4`)

The `%` glyph is omitted inside pills to save width; the column header reads `Books (devigged fair %)`.

## Edge Cases

| Case                                            | Behavior                                                                                       |
|-------------------------------------------------|------------------------------------------------------------------------------------------------|
| Model probability is NA                         | Cannot occur in practice — the R pricer skips games with no samples (`mlb_correlated_parlay.R:601-605`) so every row in `mlb_parlay_opportunities` has a real model prob. Defensive guard only: if `model_prob_raw` is NA the pill renders `M —` dimmed, same as a missing book pill. |
| All books missing prices for a combo            | Existing R pricer skips the row. No change.                                                    |
| Single book missing (e.g., DK NA, others priced) | That book's pill shows `DK —`, dimmed. Cons = mean of available inputs, exactly as today.      |
| FG vs F5 combos                                 | Period stays embedded in legs text (e.g. "F5 NYY -1.5 …"). No badge.                            |
| Existing placed parlays                         | Pill row is presentation-only. All `data-*` attributes the Place button reads are unchanged.    |

## Responsive Behavior

The phone/laptop fix is in three layers:

1. **Pill row wraps within its cell** — `flex-wrap: wrap`, so the six pills break onto a second line as the cell narrows. Row height grows; nothing overflows.
2. **Corr column hides on phones** — CSS media query at `max-width: 700px` sets `.parlay-table .corr-col { display: none; }`. Lowest-information column; Edge stays.
3. **Container `overflow-x: auto`** — safety net for very narrow viewports; scrolls instead of breaking. Already present in the reactable wrapper.

Manual viewport checks at 1400px / 1024px / 768px / 414px during testing (see Verification).

## Sort Behavior

Default sort stays `desc(edge_pct)` (unchanged). The Books cell is not sortable — it contains six values per row, and per-book sort was not useful pre-change either. The hidden `model_prob_pct` column stays in the dataframe so sort-by-model can be added later with a one-line `colDef` change if requested.

## What Is Not Changing

- **R pricer (`mlb_correlated_parlay.R`)** — no changes. All probabilities already produced.
- **Schema (`mlb_parlay_opportunities`)** — no changes.
- **Server (`mlb_dashboard_server.py`)** — no changes. The Place button payload uses the same `data-*` attributes.
- **Singles tab** — left as-is unless explicitly asked.
- **Filter bar, sizing controls, placed-parlays box** — all unchanged.
- **CBB dashboard** — not touched.
- **CLV computation, auto-place flow, scheduled captures** — all unchanged.

## Verification

Don't merge until all four pass on a fixture that includes (a) at least one game where one book is NA, (b) phone-width viewport rendered cleanly. (Item 4 below covers the synthetic NA-model-prob check separately, since real data never has model NA.)

1. **Visual regression.** Generate `report.html` from `main` and from the feature branch using the same `mlb_parlay_opportunities` snapshot. Open both in a browser, scan side by side. Check: pills render, NAs show dimmed `—`, Cons matches the existing Fair column's implied probability within rounding.
2. **Responsive check.** Resize browser to 1400 / 1024 / 768 / 414 px. Confirm pills wrap cleanly at 1024px and below, Corr hides at <700px, no horizontal scroll on phone.
3. **Action button regression.** Place a test parlay from the new layout. Confirm `placed_bets` row has correct `parlay_hash`, `game_id`, `fair_odds`, `wz_odds`, `edge_pct`, size — exactly as today.
4. **Empty-data run.** Use a fixture with one row where every book is NA except the model. Confirm pills render with five dimmed `—` and one model value, no JS errors in the console.

## Risks & Mitigations

- **Risk:** The HTML pill renderer breaks reactable filter/sort behavior in subtle ways (free-text column filter on the old per-book numeric columns).
  **Mitigation:** Per-book columns weren't filterable in the existing UI. Confirmed during exploration. Game / status / minSize / minEdge filters all read from other columns.
- **Risk:** Very long pill rows force the column wider than 320px, pushing other columns off-screen on laptops.
  **Mitigation:** `flex-wrap` is the primary defense. Pills max out at ~7 chars (`Cons 99.9`) so worst-case six pills × ~50px each ≈ 300px, fits the 320px min.
- **Risk:** Action button JS (which reads `data-prob` from the old DK pill column for some bookkeeping) breaks.
  **Mitigation:** `data-prob` is read from the Place button's own attributes (set in the `is_placed` cell renderer), not from the prob columns. Verified at `mlb_dashboard.R` lines 2266 / 2378. No coupling to remove.

## Documentation Updates

Same commit as the code change:

- `Answer Keys/MLB Dashboard/README.md` — under "Features," add a sentence describing the books strip + Cons. The existing "Parlay tab" bullet expands to mention model-vs-market visibility.

No CLAUDE.md updates are needed — the architectural facts (which columns exist in `mlb_parlay_opportunities`, where the dashboard lives, how reactable cells are rendered) are unchanged.
