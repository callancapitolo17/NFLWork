# MLB Dashboard — Pinned Scoring Strip

_Design spec · 2026-05-27 · branch `worktree-mlb-dashboard-pin-kelly-book`_

## Review Pack

**What we're building**
A "frozen" strip at the top of the MLB Dashboard so the **account pills**
(`Placing on …`) and the **Kelly Calculator** stay visible while you scroll
down the long list of bet cards ("scoring"). The account pills freeze on every
tab; the Kelly Calculator freezes (just below the pills) on the Bets tab only.

**Key decisions**
1. **Sticky strip, not a floating panel** — uses CSS `position: sticky` so the
   bar lives in normal document flow and freezes at the top on scroll. Rejected
   a fixed/floating widget (approach B) because it overlaps cards and needs
   drag/collapse JS. Rejected a collapsible click-to-expand bar (approach C)
   because you want the calc *visible*, not one click away.
2. **Minimal contents: pills + Kelly Calc only** — Bankroll/Kelly-Fraction
   inputs and the filter bar are *not* pinned. Rejected the larger variants to
   keep the frozen footprint small (your call: "Minimal").
3. **Pills pinned on all tabs; Kelly Calc only on Bets** — the pills already
   live in the global header and drive parlay/trifecta placement too, so
   freezing them globally is both useful and avoids duplicating account state.
   The Kelly Calc is a Bets-tab tool, so it only freezes there.
4. **Relocate the pills row out of `.header`** — `position: sticky` only holds
   while the element's *parent* is on screen. The pills row sits inside the
   short `.header`, so on its own it would unfreeze immediately. We move it up
   one level to be a direct child of `.container` (which spans the whole page)
   so it stays frozen through the entire card list. The big title/Refresh row is
   deliberately left behind to scroll away (keeps the footprint minimal).
5. **Always-on, no toggle, no persistence** — rejected a "📌 Pin" button; you
   want it always frozen, so there's no state to store.

**Risks / push back here**
- **The tab bar scrolls away.** Once you scroll into the cards, the pinned
  region is `[pills]` (+ Kelly Calc on Bets) — the tab-switcher buttons scroll
  off behind the pinned bar. To change tabs you scroll back up. If you switch
  tabs mid-scroll often, say so and we'll pin the tab bar too (costs ~36px of
  height).
- **The pinned strip eats vertical screen height** on the Bets tab (pills row +
  Kelly Calc row ≈ 90–110px frozen). On a laptop that's a few fewer cards
  visible at once. This is the inherent cost of approach A; flagged so it's a
  conscious trade.
- **Title row leaves the screen on scroll.** The "MLB Answer Key Dashboard" /
  "Updated …" / Refresh row is not pinned. If you want the Refresh button always
  reachable, that's a change.

**Worth understanding** _(opt-in)_
- **`position: sticky` and its parent-boundary rule.** A sticky element behaves
  like normal flow until its edge hits the offset you give it (`top: 0`), then
  it "pins." But it only stays pinned while its *containing block* (its parent)
  is still scrolling past — once you scroll beyond the parent, it lets go. This
  is the single most surprising thing about sticky and it's the whole reason we
  relocate the pills row to a tall parent. Closest R analogy: think of a
  `kableExtra` table with a frozen header row — the freeze only "works" within
  the table's own region; scroll past the table and the header is gone. Same
  idea: the freeze is scoped to the parent.

---

## Goal

When working the Bets tab, the user scrolls through many bet cards and wants two
controls always in view:
- the **`Placing on` account pills** (which Wagerzon account a placement hits), and
- the **Kelly Calculator** (manual no-vig sizing tool: type Odds + Fair → Risk / To Win / EV / Kelly%).

Today both scroll off the top: the pills sit in the global header, the Kelly
Calculator sits below the Bankroll/Kelly controls near the top of the Bets tab.

## Non-goals

- No floating/draggable panel.
- No pinning of the Bankroll/Kelly-Fraction inputs or the filter bar.
- No pin on/off toggle and no saved preference.
- No changes to placement logic, the Kelly math, server endpoints, or the DB schema.

## Current layout (as built)

DOM order inside `.container` (from `mlb_dashboard.R`):

```
.container
  .header
    .header-row-top        (title + subtitle + Refresh)
    .header-row-accounts    #wz-account-row   ← the Placing-on pills (GLOBAL)
  .tab-bar                  (Bets / Parlays / Trifectas)
  #tab-bets .tab-content
    [Placed Bets section]   (conditional)
    .sizing-controls        (Bankroll + Kelly Fraction + Apply)
    .kelly-calc             ← the Kelly Calculator widget
    .filter-bar             (Game / Book / Market / … dropdowns)
    #bets-table-container   (the scrolling bet cards)
  #tab-parlays  .tab-content (hidden unless active)
  #tab-trifectas .tab-content (hidden unless active)
```

Key facts that drive the design:
- `.header` is a **direct child of `.container`**, but `.header` itself is short,
  so a sticky element *inside* it cannot stay pinned past the header.
- `.container` spans the full page height (it wraps every tab's content).
- `#tab-bets` spans the full height of the Bets tab's content (including all cards).
- The pills are rendered by JS into `#wz-account-pills` (inside `#wz-account-row`).
  Selection + balances JS keys off these ids — see `renderPills()` etc. around
  `mlb_dashboard.R:5416+`.

## Target layout

Relocate the account row up one level so its parent is the page-spanning
`.container`, and add sticky positioning:

```
.container
  .header
    .header-row-top         (title + subtitle + Refresh)   ← scrolls away
  #wz-account-row   .pinned-account-bar   position: sticky; top: 0   ← pins on ALL tabs
  .tab-bar                                  ← scrolls away (under the pinned bar)
  #tab-bets .tab-content
    [Placed Bets section]
    .sizing-controls                        ← scrolls away
    .kelly-calc   position: sticky; top: var(--pin-top)   ← pins on the Bets tab
    .filter-bar                             ← scrolls away (under the pinned calc)
    #bets-table-container                   ← cards scroll under the strip
  #tab-parlays / #tab-trifectas             ← only the account bar pins here
```

When fully scrolled on the Bets tab the frozen region reads, top to bottom:
`[ account pills ]` then `[ Kelly Calculator ]`, with cards sliding underneath.
On the Parlays/Trifectas tabs only the account pills freeze.

## Implementation

All changes are in **`Answer Keys/MLB Dashboard/mlb_dashboard.R`** — the inline
`<style>` block plus the body markup. No server, DB, schema, or JS-logic changes;
the only JS added is a tiny height-measurement helper for the offset.

### 1. Move the account row out of the header

In the body builder (around `mlb_dashboard.R:3121`), take the existing
`#wz-account-row` div (currently the second child of `.header`) and emit it as a
**sibling of `.header`**, immediately after it and before `.tab-bar`. Keep the
`id="wz-account-row"`, the inner `#wz-account-pills` container, the
`.header-label` ("Placing on"), and the refresh button **exactly as-is** so all
existing pill JS keeps working untouched. Add a class (e.g. `pinned-account-bar`)
for the new sticky styling.

### 2. Sticky CSS for the account bar

```css
.pinned-account-bar {
  position: sticky;
  top: 0;
  z-index: 40;                 /* above cards, below modal/toast layer */
  background: #0d1117;          /* opaque so cards scroll cleanly under it */
  border-bottom: 1px solid #21262d;
  /* keep existing flex layout from .header-row-accounts */
}
```
(Fold in or keep the existing `.header-row-accounts` flex rules so the pills,
label, and refresh button still lay out the same.)

### 3. Sticky CSS for the Kelly Calculator (Bets tab)

```css
.kelly-calc {
  position: sticky;
  top: var(--pin-top, 48px);   /* sits flush below the account bar */
  z-index: 30;                  /* below the account bar, above cards */
  /* existing .kelly-calc styles unchanged */
}
```
`--pin-top` must equal the account bar's rendered height so the calc pins flush
below it with no gap or overlap. The account bar's height varies slightly when
the pills wrap on narrow widths, so we measure it rather than hardcode.

### 4. Offset-measurement JS (the only new logic)

A small handler that sets `--pin-top` to the account bar's height on load and on
resize:

```js
function syncPinOffset() {
  var bar = document.getElementById('wz-account-row');
  if (!bar) return;
  document.documentElement.style.setProperty('--pin-top', bar.offsetHeight + 'px');
}
window.addEventListener('load', syncPinOffset);
window.addEventListener('resize', syncPinOffset);
// Also call once after the pills first render (they're populated by JS),
// e.g. at the end of the existing renderPills() so the height reflects real pills.
```

This is the minimal robust approach; a `ResizeObserver` on the bar would also
work but `load`/`resize` + a post-render call covers the cases that matter.

### 5. z-index sanity

Pick `z-index` values for the two sticky elements that sit **above** the bet
cards but **below** the modal overlay and toast layer (both already use high
`z-index` in the existing CSS). The values above (40 / 30) are placeholders to
verify against the existing modal/toast `z-index` during implementation.

## Edge cases & details

- **Body padding.** `body` has `padding: 32px`. A `sticky; top: 0` element pins
  flush to the viewport top (y=0), i.e. into that 32px band. Verify the pinned
  bar looks intentional there (its own background + bottom border should read as
  a toolbar). If it looks cramped against the very top edge, give the bar a small
  internal `padding` rather than changing `top`.
- **Tab switching while scrolled.** The tab bar is not pinned; switching tabs
  deep in the list requires scrolling up. Accepted (see Review Pack risks).
- **Parlays/Trifectas tabs.** `.kelly-calc` is not present there, so only the
  account bar pins — no offset needed for those tabs.
- **Placed Bets section** (conditional, top of Bets tab) scrolls normally; the
  Kelly Calc pins only once you scroll past it and the sizing controls.
- **Narrow widths / pill wrap.** Handled by `syncPinOffset()` recomputing
  `--pin-top` on resize and after render.
- **Existing pill JS.** Unaffected — ids and inner structure are preserved; only
  the row's DOM parent and CSS change.

## Testing / verification

Per project convention (`verify_ui_features_by_rendering`), do **not** rely on a
static diff. Render and look:
1. Regenerate `report.html` (run the dashboard's R render against copied live DBs)
   and open it.
2. On the **Bets** tab, scroll into the card list and confirm: pills bar frozen
   at top, Kelly Calculator frozen directly below it (no gap/overlap), cards
   scrolling cleanly underneath.
3. Type into the Kelly Calc's Odds/Fair fields **while scrolled** and confirm it
   still computes Risk / To Win / EV / Kelly (reads Bankroll + Kelly Fraction
   from the now-scrolled-away inputs — those values persist in the DOM, so the
   math is unaffected).
4. Click a pill **while scrolled** and confirm account selection + balances still
   update.
5. Switch to **Parlays** / **Trifectas** and confirm only the pills bar freezes
   and parlay placement still reads the selected account.
6. Resize the window narrow enough that pills wrap; confirm the Kelly Calc still
   pins flush below the (now taller) pills bar (offset recomputed).
7. Open a modal / trigger a toast and confirm they still render **above** the
   pinned strip.

## Version control

- **Branch / worktree:** `worktree-mlb-dashboard-pin-kelly-book` (already created;
  this spec is being authored in it).
- **Files modified:** `Answer Keys/MLB Dashboard/mlb_dashboard.R` only.
- **Commits:** (1) this design doc; (2) the implementation (CSS + markup move +
  offset JS) as one focused commit.
- **Cleanup:** after merge to `main`, remove the worktree and delete the branch.

## Documentation

- Update `Answer Keys/MLB Dashboard/README.md` with a short "Pinned scoring strip"
  note (what freezes, on which tabs, always-on).
- Update `Answer Keys/CLAUDE.md`'s "MLB Dashboard" section to record that the
  account row was relocated out of `.header` to a sticky `.pinned-account-bar`
  sibling, and that `.kelly-calc` is sticky on the Bets tab with a JS-measured
  `--pin-top` offset (so a future editor doesn't "fix" the relocation).
- Docs land in the **same merge** as the code.
