# MLB Parlay Tab — Cards Redesign

**Date:** 2026-04-26
**Branch:** `feature/parlay-cards-redesign`
**Scope:** Parlay Opportunities table on the MLB dashboard's Parlay tab only
**Files touched:** `Answer Keys/MLB Dashboard/mlb_dashboard.R`, `Answer Keys/MLB Dashboard/README.md`
**Out of scope:** Singles (Bets) tab, placed-parlays section, combo banner placement, R pricer, server, schema, CBB dashboard

## Problem

The shipped books-strip redesign (commit `14c3721`, 2026-04-25) collapsed five per-book probability columns into one inline pill row but left the rest of the table cramped. The font sizes are too small (10–13px varying), the `Corr` column carries little signal, the `n_books_blended` column auto-renders without a `colDef` and adds noise, and the responsive strategy (one media query hiding `Corr` at <700px, written against a non-existent reactable attribute and later patched to `class="corr-col"`) does not produce a phone-usable table at any width below ~720px. The combined-parlay feature shipped between then and 2026-04-26 (commits `8675faa..d3803ba`) added a `Sel` checkbox column and a combined-parlay banner, increasing the column count back to 11 and tightening the layout further.

## Goals

1. **Universally viewable** — the Parlay Opportunities list reads cleanly at any viewport width from a wide laptop (1400px+) down through split-screen (~720px) to phone (~390px) without horizontal scroll, mode-switch hand-off, or unreadable cell sizes.
2. **Drop low-signal columns** — `Corr` and `n_books_blended` removed from view.
3. **Establish a 2-size font system** — 14px primary content, 12px secondary metadata. Replace the current ~5 size variants with two.
4. **Preserve everything else** — Sel checkbox / combined-parlay flow, Place button, conditional Kelly residual annotation, edge color thresholds, To-Win green, sort-by-edge default — all unchanged in behavior.
5. **Single-file change** — only `mlb_dashboard.R` and its README. No schema, server, R pricer, or scraper changes.

## Solution

Convert the Parlay Opportunities reactable from a horizontal table into a vertical list of cards, **using CSS only**. Each parlay opportunity becomes a self-contained card (game / time / legs / pills / metadata strip / edge + place button), rendered identically at every viewport width. Cards just get narrower as the viewport narrows; pills wrap to a second line; nothing hides, nothing switches modes.

The reactable component still emits a real `<table>` under the hood. A scoped CSS block targeting the parlay container (`#parlays-table-container`) overrides the default table-display behavior so each `<tr>` becomes a card and each `<td>` becomes a content block within the card. The cells are positioned within the card via per-column `class` attributes that we add to each `colDef`.

This is not a refactor — no R logic, no JS, no data flow changes. The reactable component, its `cell` functions, its filter/sort JS, and the Place button machinery all keep working exactly as today.

## Card Anatomy

```
┌───────────────────────────────────────────────────────────────┐
│  CHC @ SD                                                ☐    │  Game (15px) + Sel checkbox (top-right)
│  Sat 7:05p                                                    │  time (12px muted)
│  SD -1.5 (+165) · O 8.5 (-110)                                │  Legs (14px)
│                                                               │
│  M 27.4  DK 26.9  FD 27.1  PX 27.8  NV 26.7  Cons 27.2        │  Books pill row (13px pills, wraps)
│                                                               │
│  Fair +265   WZ +330   Size $24   To Win $79                  │  Metadata strip (13px, labels via ::before)
│  (residual after combo $50)                                   │  combo_residual_note (12px italic) — when present
│                                          +12.3%  [ Place ]    │  Edge (15px, color-coded) + Place button
└───────────────────────────────────────────────────────────────┘
```

## Final Visible Columns

The reactable still has columns; they just lay out as content blocks inside the card. Visible columns shrink from 11 → 10 (Corr removed; n_books_blended and is_combo were auto-rendering and now get explicit `colDef(show = FALSE)`):

| Cell class    | Source colDef    | Card position                    |
|---------------|------------------|----------------------------------|
| `cell-sel`    | `sel`            | Absolute top-right corner        |
| `cell-game`   | `game`           | Top of card (matchup + time line)|
| `cell-legs`   | `legs_display`   | Below game                        |
| `cell-books`  | `books_strip`    | Middle row                        |
| `cell-fair`   | `fair_display`   | Metadata strip — labeled inline   |
| `cell-wz`     | `wz_display`     | Metadata strip — labeled inline   |
| `cell-size`   | `size_display`   | Metadata strip — labeled inline   |
| `cell-towin`  | `to_win_display` | Metadata strip — labeled inline   |
| `cell-edge`   | `edge_display`   | Bottom-right of card              |
| `cell-action` | `is_placed`      | Bottom-right of card (right of edge) |

**Removed from view:** `corr_display` (no longer wanted) and `n_books_blended` + `is_combo` (auto-rendering today; getting explicit `colDef(show = FALSE)`).

## Data Flow

No R-pricer changes. No schema changes. No new dataframe columns.

The only R-side change in the parlay reactable is:
- Three `colDef(show = FALSE)` additions: `corr_display`, `n_books_blended`, `is_combo`.
- A `class = "cell-<name>"` parameter added to each visible `colDef`.
- The existing `class = "corr-col"` and `headerClass = "corr-col"` on `corr_display` are removed (column is gone).

Everything else in `create_parlays_table()` stays.

## CSS Specification

All rules go into the dashboard's existing inline `<style>` block (`tags$style(HTML('...'))` near line 750), scoped under `#parlays-table-container`. Approximate length: ~80 lines.

```css
/* === Parlay tab card layout (scoped — singles tab is unaffected) === */

/* Flatten the reactable table into a stack of cards */
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

/* These cells take full row width inside the flex card */
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

/* Game + time (Game cell already renders matchup + <div> with time inside) */
#parlays-table-container .rt-td.cell-game {
  font-size: 15px;
  font-weight: 500;
  padding-right: 36px;  /* room for the absolutely-positioned checkbox */
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

/* Metadata strip — Fair / WZ / Size / To Win sit inline as flex items.
   Parent `gap: 4px` handles base spacing; `margin-right` adds breathing room. */
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

/* Edge — pushed to the right via margin-left:auto on the first of the
   two right-aligned items. Action follows it in DOM order, so visually:
   [ ...metadata strip... ]                  +12.3%  [ Place ] */
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

/* Pill upsizing (was 10px in the shipped design) */
#parlays-table-container .pill {
  font-size: 13px;
  padding: 3px 8px;
}
```

## Edge Cases

| Case | Behavior |
|------|----------|
| **Empty parlay opportunities** | Reactable renders nothing inside `.rt-tbody`. The existing empty-state above the table is unchanged. |
| **Filtering hides rows** | `applyParlayFilters()` toggles `display: none` on `.rt-tr-group`. Cards hide the same way. |
| **Pagination** | `defaultPageSize = 25` still applies. Cards paginate as rows did. |
| **Sort UI** | `.rt-thead` is hidden, so column sort buttons aren't user-clickable. No UI was relying on them; default sort `desc(edge_pct)` still applies. |
| **Sel checkbox on placed / combo rows** | Existing R logic in the `sel` cell function (`isTRUE(row$is_placed) || isTRUE(row$is_combo) → ''`) keeps the cell empty on those rows. The absolute-positioned slot stays empty; nothing renders. |
| **`combo_residual_note` annotation** | Already emitted today as `'<br><span class="combo-note">...</span>'` inside the Size cell. New `.combo-note` rule with `display: block; width: 100%` lifts it onto its own line below the metadata strip. |
| **Edge color thresholds** | Computed inside `edge_display`'s existing cell function (≥15 / ≥10 / ≥5 / default). CSS doesn't touch color logic. |
| **To-Win green** | `style = list(color = "#3fb950")` on `to_win_display` colDef preserved. |
| **Place button data attributes** | Unchanged. Button still emitted with all 12 `data-*` attributes; JS still finds them via `[data-hash]` and `[data-game-id]` regardless of card vs. table layout. |
| **Combined-parlay banner above table** | Lives outside `#parlays-table-container`. Untouched. |
| **Singles tab** | Not under `#parlays-table-container`. Renders as a normal table. |
| **JS that walks `.rt-tr-group`** (e.g., `applyParlayFilters`) | Selectors still match — `display: block` doesn't change DOM structure, only layout. |
| **Pill row wrapping at narrow widths** | `.books-strip` already uses `flex-wrap: wrap` from the previous round. Continues to work. |

## What Is Not Changing

- **R pricer (`mlb_correlated_parlay.R`)** — no changes.
- **Schema (`mlb_parlay_opportunities`)** — no changes.
- **Dashboard server (`mlb_dashboard_server.py`)** — no changes.
- **`books_strip.R`** — no changes (already produces the right HTML).
- **Combined-parlay flow** — Sel checkbox, banner, `onComboSelectChange`, `/api/price-combined-parlay`, `/api/place-combined-parlay` all unchanged.
- **CLV computation, auto-place, scheduled captures** — all unchanged.
- **Singles (Bets) tab** — untouched. Still renders as a standard reactable table.
- **CBB dashboard** — not touched.

## Verification

Don't merge until all five pass.

1. **Smoke run.** Copy live `mlb.duckdb` + `mlb_dashboard.duckdb` into the worktree, run `Rscript mlb_dashboard.R`, confirm `report.html` regenerates without R errors and includes at least 1 `.rt-tr` element on the parlay tab.
2. **Card layout sanity.** Open `report.html`. Parlay tab. Confirm each parlay opportunity renders as a card with the layout from the anatomy diagram above. Singles tab still looks like a table.
3. **Resize sweep.** Drag the browser window to ~1400px / ~860px (split-screen) / ~400px (phone). Cards proportionally scale. Pills wrap to 2 lines below ~700px. No horizontal scroll on phone.
4. **Combined-parlay flow regression.** Tick two Sel checkboxes on different cards. Combined-parlay banner appears with combined pricing. Click "Place Combined Parlay." Combo lands in placed-parlays section. Click Remove. Combo disappears. Validates the recent-commit machinery still works after layout flip.
5. **Place + Remove single parlay.** Click Place on any card at $1. Toast appears, button flips to Placed, row appears in placed-parlays. Click Remove. Reverts.

## Risks & Mitigations

- **Risk:** A future reactable version changes its internal class names (`.rt-tr`, `.rt-td`, etc.). The CSS would silently stop matching and the layout would revert to a table.
  **Mitigation:** Pin the reactable version in `renv.lock` or document the expected class names. (Out of scope for this design — flag for future work if reactable bumps.)
- **Risk:** Float-based positioning of Edge + Place can cause clearance bugs at extreme narrow widths.
  **Mitigation:** Verify in the resize sweep at 400px. If it breaks, swap floats for a flex container on the card itself (handled in implementation if needed).
- **Risk:** `padding-right: 36px` on `.cell-game` may not be enough room for very long matchup names + a checkbox.
  **Mitigation:** The matchup format is "Away Team @ Home Team" using full team names; width tested at 400px during verification. Bump to 44px if any team name overflows.
- **Risk:** Existing JS (`applyParlayFilters`, `addPlacedParlayRow`) relies on querying `.rt-tr-group` and reading button data attributes. The CSS flip doesn't change DOM structure, only `display`, but combined-parlay JS that walks specific cells could break.
  **Mitigation:** Verification step 4 explicitly exercises the combined-parlay flow end-to-end. Failure caught there.

## Documentation Updates

Same commit (or a follow-up commit) as the implementation:

- `Answer Keys/MLB Dashboard/README.md` — replace the existing "books strip" sentence under Features with one describing the card layout and universal-viewability across viewports.

No CLAUDE.md updates needed — the architectural facts are unchanged.
