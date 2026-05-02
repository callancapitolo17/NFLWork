# MLB Dashboard — Wagerzon Account Bar Merge into Header

**Date:** 2026-05-02
**Branch:** `feature/mlb-dashboard-header-merge`
**Worktree:** `.worktrees/mlb-dash-header-merge`
**Scope:** Layout/UI only. No backend, schema, or API changes.

## Problem

The Wagerzon multi-account header bar (shipped in `b31f7fb`, see
`fee03cf`/`eda227e`) sits as a full-width strip above `.container`,
while the rest of the dashboard is centered inside a 1600px container.
This creates two stacked "shelves" before any data is visible:

1. WZ bar — full-bleed, edge-to-edge.
2. `.header` — `MLB Answer Key Dashboard` h1 + subtitle + green Refresh.

Concrete issues observed on the live dashboard:

- **Bar misaligns with the container.** Pills sit far-left of the page
  title, breaking visual flow.
- **Title feels demoted.** The bar above has more visual weight than
  the page H1 underneath.
- **"Placing on:" + selector are far from the pills.** The selector
  is on the right, the thing it selects is on the left.
- **The selected-pill indicator is a 2px border** against a 1px default
  border — easy to miss.
- **Two distinct refresh controls** (↻ for balances, green "Refresh"
  for data) read as one ambiguous "refresh".
- **Bar uses inline styles**, inconsistent with the rest of the file's
  CSS-class system.
- **Error rendering is alarmist:** when the balance fetch fails, every
  pill shows `WagerzonX: — ⚠ (stale 0s ago)`. The "stale 0s ago" text
  is contradictory (fresh data labelled stale).

The underlying `wz_error` on all three accounts is a separate
auth/scraper bug. It is **out of scope** for this spec — the UI just
has to degrade gracefully.

## Solution — option B3: dedicated pill row inside `.header`

Replace the standalone `wz-account-bar` with a pill row that lives
inside the existing dashboard `.header`, as its own second row. One
unified header. Two rows:

- **Row 1:** Title + subtitle on the left. Green "Refresh" button on
  the right. (Same as today.)
- **Row 2:** `Placing on` caption + clickable account pills + balance
  refresh icon. (New.)

Pills become click-to-select — the `<select>` dropdown is removed.
The active pill renders with a filled blue background, so it pops
against the inactive pills.

### Layout (after)

```
.container (max 1600px, centered)
└── .header
    ├── row 1 (flex row, justify-between)
    │   ├── h1  "MLB Answer Key Dashboard"
    │   │   .subtitle "Updated YYYY-MM-DD HH:MM"
    │   └── .refresh-btn  "Refresh"
    └── row 2 (flex row, align-center)
        ├── .header-label "Placing on"
        ├── .wz-pill[.selected] WagerzonJ · $1,250.42
        ├── .wz-pill            Wagerzon  · $874.10
        ├── .wz-pill            WagerzonC · $402.55
        └── .wz-icon-btn  ↻
.tab-bar
.tab-content
...
```

### What gets removed

- The entire `tags$div(id = "wz-account-bar", ...)` block above
  `.container` (current `mlb_dashboard.R` lines ~1992–2016).
- `tags$select(id = "wz-account-select", ...)` — the dropdown.
- The `<select>` change-listener inside the JS controller block.

### What gets added

- Inside the existing `.header`, a sibling `tags$div(class = "header")`
  reorganises into:
  - `tags$div(class = "header-row-top", ...)` — current title/refresh.
  - `tags$div(class = "header-row-accounts", id = "wz-account-row", ...)`
    — caption + `#wz-account-pills` + `#wz-refresh-btn`.
- Pill rendering becomes click-driven. The renderer attaches a click
  handler that:
  1. Sets `window.WZ_SELECTED_ACCOUNT = label`.
  2. POSTs to `/api/wagerzon/last-used` with `{label}`.
  3. Re-renders pills (selected one toggles the `.selected` class).
  4. Calls `window._wzRecomputeWarnings()` so insufficient-balance
     warnings on parlay rows recompute against the new account.

### CSS classes (new)

Sit in the existing `<style>` block alongside `.header`, `.tab-bar`,
`.stat-card`. Replace the current inline `pill.style.cssText` lines.

```css
.header-row-top {
  display: flex; justify-content: space-between; align-items: center;
  padding-bottom: 10px;
}
.header-row-accounts {
  display: flex; align-items: center; gap: 8px;
  padding: 10px 0 4px 0;
}
.header-label {
  font-size: 11px; color: #8b949e;
  text-transform: uppercase; letter-spacing: 0.5px;
}
.wz-pill {
  padding: 4px 10px; border-radius: 14px;
  background: #21262d; border: 1px solid #30363d;
  color: #c9d1d9; font-size: 13px;
  cursor: pointer; user-select: none;
  transition: border-color 0.12s, background 0.12s;
}
.wz-pill:hover { border-color: #58a6ff; }
.wz-pill.selected {
  background: #1f6feb; border-color: #1f6feb;
  color: #ffffff; font-weight: 600;
}
.wz-pill.stale {
  background: #3a1d1d; border-color: #4a2a2a; color: #ffa198;
}
.wz-pill.empty {
  background: transparent; border-style: dashed;
  color: #6e7681; cursor: default;
}
.wz-icon-btn {
  background: transparent; border: 1px solid #30363d;
  color: #8b949e; width: 28px; height: 28px;
  border-radius: 6px; cursor: pointer;
  display: inline-flex; align-items: center; justify-content: center;
  font-size: 14px;
}
.wz-icon-btn:hover { color: #c9d1d9; border-color: #58a6ff; }
```

### Component behaviour

| State | Pill rendering | Click | Notes |
|---|---|---|---|
| Healthy + selected | `.wz-pill.selected` text `Label · $X,XXX.XX` | no-op | one per account |
| Healthy + inactive | `.wz-pill` text `Label · $X,XXX.XX` | sets active | hover shows blue border |
| Error, fresh (`stale_seconds < 60`) | `.wz-pill` text `Label · — ⚠` | sets active | no "stale 0s ago" |
| Error, stale (`stale_seconds >= 60`) | `.wz-pill.stale` text `Label · — ⚠ (stale Nm ago)` | sets active | red-tinted |
| 0 accounts configured | Single `.wz-pill.empty` text "No Wagerzon accounts configured" | no-op | row still rendered for layout consistency |
| 1 account | One pill, always `.selected` | no-op | refresh icon still works |

The "stale 0s ago" bug is fixed: the stale suffix only appears when
`stale_seconds >= 60`. Fresh-but-errored pills show just `— ⚠` —
neutral, not alarmist.

### What stays the same

- `GET /api/wagerzon/balances` — unchanged.
- `GET/POST /api/wagerzon/last-used` — unchanged.
- `POST /api/place-parlay` — still receives `account` from
  `WZ_SELECTED_ACCOUNT`, no contract change.
- `placed_parlays.account` schema — unchanged.
- `dashboard_settings.wagerzon_last_used` persistence — unchanged.
- `window._wzRecomputeWarnings()` for insufficient-balance pings on
  parlay rows — same logic, just driven by class changes instead of
  `<select>` change events.
- `MutationObserver` on `parlays-table-container` — unchanged.
- The dashboard-wide green "Refresh" button — unchanged in placement,
  styling, and behaviour. Only the WZ ↻ moves with the pills.

## Files Affected

| File | Change |
|---|---|
| `Answer Keys/MLB Dashboard/mlb_dashboard.R` | Layout (R `tags$*`), CSS block, JS controller. ~80 lines edited, ~30 net removed. |

No other files touched. No backend, no schema, no tests added (the
existing combined-parlay tests in `tests/test_combined_parlay.py`
don't exercise the layout).

## Edge Cases & Resilience

- **Empty state** (0 WZ accounts in `.env`): `/api/wagerzon/balances`
  returns `{"balances": []}`. Pill row shows the single empty pill.
  Insufficient-balance warnings stay silent — `WZ_SELECTED_ACCOUNT`
  is `null`, the warning helper short-circuits on `available === null`.
- **Single account** (most common configuration): pill is always
  selected, click is a no-op, refresh still works.
- **All accounts errored** (today's live state): pills render with the
  truncated error text. Selected pill keeps its `.selected` class even
  while errored — the user can still place if they have stale-but-valid
  knowledge of their balance.
- **Reactable re-renders**: existing `MutationObserver` on
  `parlays-table-container` runs `_wzRecomputeWarnings()` after table
  hot-swaps. No change.
- **Page-load race**: `loadLastUsed()` runs before `refreshBalances()`.
  If `last-used` returns `{"label": null}` on a fresh DB,
  `WZ_SELECTED_ACCOUNT` is `null` until pills render. The new pill
  renderer must explicitly default to the first label in render order
  and persist it via `POST /api/wagerzon/last-used` so the next page
  load is stable. (Today's `<select>` got this for free via browser
  default-first-option behaviour; a `<div>`-based pill row does not.)

## Out of Scope

- The underlying `wz_error` on all 3 accounts (separate auth/scraper
  bug). The UI degrades gracefully but does not fix the source.
- Sticky header on scroll (not currently sticky; out of scope).
- Insufficient-balance warning copy or placement on parlay rows.
- The `Refresh` button label clarification — keeping it as "Refresh"
  is fine once the WZ ↻ no longer competes for the same word at the
  top of the page (it's the only "Refresh" left in the global header).

## Version Control

- **Branch:** `feature/mlb-dashboard-header-merge`
- **Worktree:** `.worktrees/mlb-dash-header-merge`
- **Commit structure:** one commit, per `CLAUDE.md` rule that docs
  ship with the feature:
  - `feat(mlb-dashboard): merge WZ account bar into header pill row` —
    `mlb_dashboard.R` layout/CSS/JS + README update in the same commit.
- **Files created:** none (this spec lives at
  `docs/superpowers/specs/2026-05-02-mlb-dashboard-header-merge-design.md`).
- **Files modified:** `Answer Keys/MLB Dashboard/mlb_dashboard.R`,
  `Answer Keys/MLB Dashboard/README.md`.
- **Files deleted:** none.

## Worktree Lifecycle

1. Worktree created at `.worktrees/mlb-dash-header-merge` on branch
   `feature/mlb-dashboard-header-merge` (already done).
2. Implementation + manual visual test in worktree (open
   `http://localhost:8083` after running the dashboard from the
   worktree path).
3. Pre-merge review: `git diff main..HEAD` reviewed against the
   pre-merge checklist in `CLAUDE.md`.
4. User approves merge.
5. Merge to `main` with `--no-ff` to preserve branch history (matches
   recent merges like `b31f7fb`, `864542b`).
6. Worktree removed: `git worktree remove .worktrees/mlb-dash-header-merge`.
7. Branch deleted: `git branch -d feature/mlb-dashboard-header-merge`.

## Documentation

- **`Answer Keys/MLB Dashboard/README.md`** — needs a one-liner under
  the multi-account section noting that pills are now click-to-select
  and live in the dashboard header (no separate bar). The current
  README references the bar as a separate component.
- **`Answer Keys/CLAUDE.md`** — the "MLB Dashboard — Wagerzon
  multi-account" section talks about endpoints/schema, not layout. No
  change needed.
- **`wagerzon_odds/CLAUDE.md`** — unchanged (account discovery is
  unaffected).

## Manual Test Plan

After implementation, verify in browser at `localhost:8083`:

1. **Multiple accounts, healthy:** All pills render with balances. The
   one matching `wagerzon_last_used` is filled blue. Click another;
   it becomes filled, prior one returns to grey, and refreshing the
   page persists the selection.
2. **All accounts errored** (current live state): Pills show
   `Label · — ⚠` without "stale 0s ago" suffix. Selected pill is
   still visually distinct.
3. **Stale > 60s:** Force a stale state, confirm the suffix
   `(stale Nm ago)` appears and the pill goes red-tinted.
4. **Single account:** Pill is always selected, hover does nothing
   special, click is a no-op.
5. **0 accounts** (rename `.env` entries away temporarily): "No
   Wagerzon accounts configured" pill renders, no JS errors.
6. **Place parlay:** Picks the active account label, persists to
   `placed_parlays.account`, balance pill updates after success.
7. **Insufficient balance:** With a parlay risk > active account
   available, the row's `.wz-insufficient-warning` populates. Switch
   to an account with sufficient balance, warning clears.
8. **Reactable re-render:** Trigger a parlay-table hot-swap (place a
   combined parlay), confirm warnings recompute on the new rows.
9. **Container alignment:** At default and resized widths, the header
   pill row aligns with the title and tabs (no full-bleed misalignment).
