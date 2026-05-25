# MLB Dashboard — Multi-Account Re-Place (v1)

**Status:** Draft for review
**Date:** 2026-05-23
**Author:** brainstorming session w/ claude
**Scope:** MLB Dashboard bets tab, Wagerzon only

---

## Review Pack

**What we're building** — When a Wagerzon account hits its per-wager max
on a bet you want to place big on, today you have to manually re-pick
the account in the header and click Place again — but the original
card is locked into a "placed" state, so you can't. This spec adds a
"+ another" affordance on placed Wagerzon cards that re-opens the
hero strip for a second placement on a different WZ account, with each
placement tracked as its own row in `placed_bets`.

**Key decisions**

1. **Re-open & re-place, not a split panel.** Considered an inline
   "split panel" that lets you type stakes for N accounts and fires
   them in parallel from one click. Rejected for v1: re-open is a
   strict subset of the same schema work, ships in 1-2 days vs 1-2
   weeks, and we don't yet know how often you actually split. If
   re-open friction proves real after 2-3 weeks of use, the panel
   builds on top with zero schema rework.
2. **Inline chips for the placed display, not stacked rows or
   summary.** Inline chips (`placed · WZ $200 #W1234 · WZJ $150 #W5678`)
   keep the placed strip the same single-line shape as today's
   single-ticket label and scale horizontally as more chips appear.
   Stacked rows grow the card vertically (worse on a dense bets list);
   summary-only ("Placed $350 across 2 accts") hides ticket numbers.
3. **v1 spreads across accounts; does NOT stack on the same account.**
   The composite PK `(bet_hash, account)` enforces this — you can't
   place the same bet twice on Wagerzon. The stated use case is hitting
   per-wager max, not splitting one account into N tickets. If we ever
   want stacking, the PK needs to include a sequence number.
4. **Explicit "+ another" button + header pill switch, not
   pill-switch-auto-unlocks.** Considered: switching the header pill
   automatically un-locks any placed WZ card that doesn't have a chip
   for the new account. Rejected because it would silently "wake up"
   unrelated placed cards on the page whenever you switch the pill —
   jarring and easy to misread. Explicit button = explicit intent.

**Risks / push back here**

- **Frequency assumption:** v1 assumes multi-account splits happen <5×/week. If you
  start doing 4+ account splits every day, the manual click cadence
  (per account: switch pill → edit Risk → click Place) becomes a real
  CLV drag and the split panel becomes the right next step. Worth
  reviewing after 2-3 weeks of use.
- **No stacking on same account:** if WZ caps you at $50/wager and you
  want to put $500 down on Wagerzon alone, this feature doesn't help
  — you'd need stacking, which is a future schema change.
- **Header pill switch clears `data-expected-win` on every card on the
  page.** Minor side-effect: if you were mid-review on a totally
  unrelated card with a custom Risk amount, that card's verified quote
  also clears. Invisible in practice (the next Risk edit re-verifies
  silently before the Place button fires), but worth flagging.
- **WZ-only scope:** Hoop88 / BFA / DK / FD do not have a
  multi-account registry, so "+ another" never renders for those
  picks. If a future book grows a registry, the rendering rule needs
  to generalize (probably a `books_with_multi_account` set).

**Worth understanding** *(opt-in, anchored to R)*

- **Composite primary keys.** In R you might dedupe a `data.frame` by
  combining columns: `df[!duplicated(paste0(bet_hash, account)), ]`.
  A SQL `PRIMARY KEY (bet_hash, account)` is the formal version — the
  database itself enforces uniqueness across the *pair*, not either
  column alone. So two rows with the same `bet_hash` are perfectly
  legal as long as their `account` differs. This is the single
  enabling change for "multiple placements on the same bet" — every
  other feature change in this spec depends on it.
- **In-place DOM swap vs full re-render.** Today's `_replaceActionCell`
  helper mutates the existing card element to show "placed" without
  reloading the page — cheap and instant. When we add a second chip,
  we just append another `<span class="placement-chip">` to the
  existing hero-placed strip instead of regenerating the card. It's
  analogous to mutating one column of a `data.frame` in R vs
  reassigning the whole frame: smaller, faster, and avoids losing
  unrelated UI state (e.g. a Risk edit you just typed).

---

## 1. Background

The MLB Dashboard bets tab (port 8083) currently supports single-account
Wagerzon placement: the "Placing on" header pills are single-select, and
clicking Place on a bet card fires one bet against the currently-selected
account via `/api/place-bet` → `single_placer.place_single()`. Each bet
is keyed by `bet_hash` in `placed_bets` (PRIMARY KEY).

When a Wagerzon account hits its per-wager limit, WZ rejects with
`bet_too_large` ("exceeds limit"). Today's recovery: manually switch the
header pill to another WZ account and click Place again — but the
original card is locked in a "placed" or "error" state, and the
`bet_hash`-as-PK constraint means the new placement would overwrite the
breadcrumb row from the first attempt. The user must instead manually
log the second placement as a separate bet entirely, losing the link
between the two halves of the same logical bet.

### Existing plumbing this spec builds on

- **WZ account registry**: `wagerzon_odds/wagerzon_accounts.py` discovers
  accounts from env vars (`WAGERZON_USERNAME`, `WAGERZONJ_USERNAME`,
  `WAGERZONC_USERNAME`, etc.) and exposes `list_accounts()` /
  `get_account(label)`.
- **Header pill row**: `mlb_dashboard.R` line ~3153 renders the "Placing
  on" pills. Selection is persisted via
  `/api/wagerzon/last-used` and exposed to JS as
  `window.WZ_SELECTED_ACCOUNT`.
- **Place flow**: `placeBet(btn)` (line ~4177) reads
  `window.WZ_SELECTED_ACCOUNT`, POSTs to `/api/place-bet` with body
  containing `account`. Server dispatches to
  `single_placer.place_single(account, bet)` and returns a status
  envelope.
- **Click-to-edit Risk + WZ-verified quote**: the hero strip's Risk
  cell is editable in-place (`risk-value` span, hero strip lines 1267-
  1273). On edit, the JS coordinator fires `/api/wz-quote-single` with
  `{bet_hash, amount, account}` and writes the returned win value to
  `data-expected-win`, which the placer uses for its Win-on-Win drift
  check.
- **In-place DOM swap**: on successful placement, `_replaceActionCell`
  swaps the Place button for a `<span class="placed-bet-label">` chip
  showing "placed · #ticket". The R-side regen on next dashboard
  refresh produces an identical chip, so the JS swap is visually
  consistent with the server-rendered state.

### What's missing

There is no path to:
1. Add a second placement row for the same `bet_hash` with a different
   `account`.
2. Visually accumulate multiple ticket chips on one card.
3. Re-open the hero strip for editing after a successful placement.

## 2. Goals & non-goals

### Goals (v1)

- Allow placing a Wagerzon bet on up to N different WZ accounts where
  N = number of WZ accounts registered.
- Each placement is independently tracked in `placed_bets` (separate
  row, separate ticket number, separate risk amount).
- Per-placement Risk is user-editable (default = previous Risk value or
  Kelly suggestion) and verified against the live WZ quote via the
  existing `/api/wz-quote-single` endpoint before firing.
- No regression in single-account behavior: if you only ever place
  once per bet, the card flow looks exactly like today plus a small
  dashed "+ another" button on the placed strip.

### Non-goals (v1)

- Stacking multiple placements on the **same** account (e.g. $200 on
  Wagerzon now, $300 on Wagerzon again later). Composite PK
  `(bet_hash, account)` blocks this; deliberate v1 constraint.
- Batch / parallel "place across all accounts in one click." That's
  the split panel; revisit if click cadence proves painful.
- Multi-account support for non-WZ books (Hoop88, BFA, DK, FD). None
  of those have a multi-account registry today.
- Automatic fall-back ("WZ rejected with exceeds-limit, try WZJ
  automatically"). Manual control preferred for v1.
- Visual cross-card linking ("these two ticket chips belong to the
  same logical bet on different cards"). Each chip already lives on
  the same card, so the link is implicit.

## 3. UX flow

The full 4-step walkthrough is mirrored in
`.superpowers/brainstorm/97806-1779557989/content/walkthrough.html`.
Summary:

| Step | State | User action |
|---|---|---|
| 1 | Card in "place mode", header pill on Wagerzon | Click Place |
| 2 | Card flips to placed state with 1 chip + dashed `+ another` button | Click `+ another` |
| 3 | Hero strip re-opens; chip remains visible above; Place is **disabled** because header pill still matches the chip | Click WagerzonJ pill |
| 4 | Place re-enables; edit Risk → quote re-verifies under WagerzonJ → click Place; second chip lands | (optionally repeat for WZC) |

### State machine for the placed strip

```
[empty]
   |
   | first Place succeeds
   v
[placed, 1 chip] -- click "+another" --> [re-opened, 1 chip + editable hero]
                                                |
                                                | header pill switched
                                                | to untouched WZ account,
                                                | Risk edited, Place clicked,
                                                | WZ accepts
                                                v
                                          [placed, 2 chips]  (loop)
                                                |
                                                | all WZ accounts now in chip list
                                                v
                                          [placed, N chips, "+another" disabled]
```

### Rendering rules

- **"+ another" visibility**: rendered on the placed strip iff
  `book == "wagerzon"` AND the bet has at least one WZ account not
  already in its chip list.
- **Place button enable state in re-opened mode**: disabled when
  `window.WZ_SELECTED_ACCOUNT ∈ chip_accounts(bet_hash)`. Disabled
  state shows tooltip "Switch the header pill to a different account."
- **Chip ordering**: chronological (oldest first, left to right).
  Matches the placement order, which usually matches the registry
  order anyway.
- **Risk default on re-open**: the previous successful Risk amount
  (e.g. $200), not Kelly — assumption is the user is splitting the
  same intended stake across accounts.
- **Toast on `+ another` click**: a small informational toast: "Switch
  the header pill to add another account." Auto-dismisses on header
  pill click.

## 4. Data model changes

### `placed_bets` table

Current PK: `bet_hash TEXT PRIMARY KEY`
New PK: `PRIMARY KEY (bet_hash, account)`

`account` is already a column (added by migration `002_single_placer_columns.py`,
2026-05). Today's pre-multi-account rows have `account = NULL` or `'Wagerzon'`;
both must be backfilled before the PK migration runs.

Backfill rule:
- If `account IS NULL`: set to `'Wagerzon'` (the legacy single-account
  default).
- If `account = ''`: same.
- Otherwise: leave as-is.

### Migration script

New file: `Answer Keys/MLB Dashboard/migrations/003_placed_bets_composite_pk.py`

Idempotent steps:
1. Verify all rows have non-null, non-empty `account` (backfill if not).
2. `DROP TABLE IF EXISTS placed_bets_new`.
3. `CREATE TABLE placed_bets_new (... PRIMARY KEY (bet_hash, account))`
   with identical column definitions.
4. `INSERT INTO placed_bets_new SELECT * FROM placed_bets`.
5. `DROP TABLE placed_bets`.
6. `ALTER TABLE placed_bets_new RENAME TO placed_bets`.
7. Print row counts before/after for sanity.

Run as part of the worktree merge, before restarting the dashboard
server. Manual command documented in the README update.
**Restart the dashboard server after migration** — the running R
process holds a connection to the old schema and will not see the
PK change without a restart.

### `/api/place-bet` 409 in-flight check

Current behavior (`_insert_placement_breadcrumb`): if a row exists with
status `placing`, return 409 "bet already in flight."

New behavior: scope the in-flight check by `(bet_hash, account)`.
Placement on WZ for a bet that has an in-flight `WagerzonJ` row is
allowed (different account). Placement on WZ for a bet with an
in-flight `Wagerzon` row is still blocked (same account = duplicate
attempt).

Implementation: one-line SQL change in `_insert_placement_breadcrumb`
to add `AND account = ?` to the existence check.

## 5. Server changes

Minimal — most of the work was already done when `account` became a
required field of `/api/place-bet` (migration 002 timeframe).

Changes:
- `_insert_placement_breadcrumb` 409 check scoped by `(bet_hash, account)` as above.
- `_finalize_placement` (and any UPSERT that touches `placed_bets`)
  scoped by composite key — DuckDB will reject the upsert if it
  doesn't match the new PK. Audit all writers and update WHERE
  clauses.

No new endpoints. `/api/place-bet` and `/api/wz-quote-single` work
unchanged at the wire level — both already accept `account` in the
body.

## 6. Frontend changes (`mlb_dashboard.R`)

### New rendering helper

`render_placed_strip(chips)` — produces the hero-placed strip HTML
with one `<span class="placement-chip">` per chip, plus the `+ another`
dashed button if the rendering rule above is satisfied.

Chip shape (mirrors today's single placed-bet-label, multiplied):
```html
<span class="placement-chip"
      data-account="Wagerzon"
      data-risk="200"
      data-ticket="W1234">
  <span class="acct">WZ</span>$200 <span class="ticket">#W1234</span>
</span>
```

### CSS additions

New rules under `.bet-card-v8` (in the existing inline `<style>` block
at line ~2891):
- `.hero-placed` — green-tinted flex container, mirrors `.hero` but
  shorter (no Fair/EV/Risk/To Win stats, just chips).
- `.placement-chip` — green pill with account label + amount + ticket.
- `.add-another` — dashed-border ghost button.

### New JS handlers

- `addAnother(btn)` — handler for the `+ another` button click:
  - Swap the placed strip out, swap an editable hero strip in (mirror of the
    pre-placement hero, with Risk pre-filled from the previous chip).
  - Show toast "Switch the header pill to add another account."
  - Disable the Place button if `window.WZ_SELECTED_ACCOUNT` matches
    any chip's `data-account`.

- `placeAnother(btn)` — invoked by Place click in re-opened mode:
  - Same body shape as today's `placeBet` (uses
    `window.WZ_SELECTED_ACCOUNT`, current Risk, `data-expected-win`).
  - On success, append a new `placement-chip` to the existing strip
    (don't regenerate); recompute `+ another` enable state; collapse
    back to placed-strip view.

- **Header pill click handler addition**: on every pill click, walk all
  cards on the page and clear `data-expected-win` (set to empty
  string). Also re-evaluate Place button enable state on any re-opened
  card.

### Card render path

`create_bets_table()` in `mlb_dashboard.R` currently passes either the
editable hero (un-placed) or a `placed-bet-label` chip (placed). On
page load, it needs to:

1. Query `placed_bets` for all rows with `bet_hash IN (...)` (already
   done; just add `account` to the SELECT).
2. Group by `bet_hash` → list of chips.
3. If `length(chips) > 0`: render `render_placed_strip(chips)`.
4. If `length(chips) == 0`: render the editable hero as today.

Visual continuity: the JS chip-append in `placeAnother` produces the
same HTML the R render would on a page refresh. Tested by reloading
the page mid-flow.

## 7. Edge cases

| Case | Behavior |
|---|---|
| First WZ account placement gets rejected (line pulled, exceeds limit) | No chip added. Card stays in "place mode" with an error toast. Same as today. |
| User clicks `+ another` but never switches the header pill, then clicks elsewhere | Re-opened state persists on the card (does not auto-collapse). Next page refresh server-renders only the chips that exist in `placed_bets`, so the re-opened editable hero goes away on refresh. Acceptable. |
| User clicks `+ another`, switches header pill, then switches it back to the already-placed account | Place button disables again. Toast doesn't repeat. |
| All N WZ accounts have been placed on for a bet | `+ another` is rendered but disabled with tooltip "All WZ accounts placed." |
| Two browser tabs both on the dashboard, both place on different accounts simultaneously | Composite PK + 409 in-flight check by `(bet_hash, account)` handle this naturally. No race. |
| WZ pill switch happens while a Risk edit is in-flight (`/api/wz-quote-single` pending) | Existing handler debounces; new account context will re-verify on the next edit. Pending response from the old account is discarded by the existing handler. |
| Non-WZ pick (Hoop88, BFA, etc.) gets placed | No "+ another" button rendered. Card looks exactly like today. |
| Migration runs against a `placed_bets` table with NULL `account` rows | Backfill step sets NULL → 'Wagerzon' before PK creation. Idempotent: re-runs are no-ops. |

## 8. Version control plan

- **Branch**: `worktree-feature+dashboard-multi-account-re-place-spec`
  (already created via worktree).
- **Worktree**: `/Users/callancapitolo/NFLWork/.claude/worktrees/feature+dashboard-multi-account-re-place-spec`
  (already entered; this spec is being written inside it).
- **Files modified**:
  - `Answer Keys/MLB Dashboard/mlb_dashboard.R` — CSS, JS handlers,
    `create_bets_table` chip-list rendering, header pill click hook.
  - `Answer Keys/MLB Dashboard/mlb_dashboard_server.py` — 409 check
    scoping, finalize-placement UPSERT scoping.
  - `Answer Keys/MLB Dashboard/migrations/003_placed_bets_composite_pk.py`
    (new).
  - `docs/superpowers/specs/2026-05-23-mlb-dashboard-multi-account-re-place-design.md`
    (this file).
- **Commit structure** (suggested):
  1. spec doc (this file)
  2. migration script + manual run instructions
  3. server changes (409 + UPSERT scoping)
  4. frontend rendering helper + CSS
  5. frontend JS handlers (`addAnother`, `placeAnother`, header pill
     hook)
  6. README updates
- **Pre-merge review**: full diff review per CLAUDE.md pre-merge
  checklist before requesting merge approval.
- **Worktree cleanup**: after merge, `git worktree remove` +
  `git branch -d` the feature branch.

## 9. Documentation updates

- **`Answer Keys/CLAUDE.md`** — the "MLB Dashboard — Wagerzon
  multi-account" section currently describes single-account placement
  semantics. Add a paragraph noting:
  - `placed_bets` PK is now composite `(bet_hash, account)`.
  - Multiple placements per bet across accounts are supported via the
    "+ another" affordance.
  - v1 does not support stacking on the same account.
- **`Answer Keys/MLB Dashboard/README.md`** (if it exists; create if
  not) — usage note for `+ another` button.
- **Migration README** — add `003_placed_bets_composite_pk.py` to the
  "manual migrations" list with the run command.

## 10. Testing

- **Migration test**: create a throwaway DuckDB with the old PK,
  insert a few rows with NULL account, run the migration, verify
  backfill + new PK + same row count. New script:
  `Answer Keys/MLB Dashboard/migrations/test_003_placed_bets_composite_pk.py`.
- **Manual end-to-end test plan** (no automated UI test for the
  dashboard yet):
  1. Place a WZ bet on Wagerzon at $50. Verify chip + `+ another`
     appears.
  2. Click `+ another`. Verify Place is disabled.
  3. Switch to WagerzonJ. Verify Place enables.
  4. Edit Risk to $40. Verify To Win updates from WZ quote.
  5. Place. Verify second chip lands.
  6. Repeat for WagerzonC. Verify `+ another` becomes disabled after
     all 3 are placed.
  7. Refresh page. Verify all 3 chips persist (server-rendered from
     `placed_bets`).
  8. Place a Hoop88 bet. Verify no `+ another` button appears (non-WZ
     book).
- **Regression spot-checks**:
  - Single-account WZ placement (today's flow) still works end-to-end.
  - 409 in-flight check still triggers on rapid double-click of Place
    on the same account.
  - WZ-verified quote (Win-on-Win drift check) still triggers correctly
    on the second placement.
  - Removing a placement (existing `/api/remove-bet` flow, if it's
    used) still works under composite PK — needs SQL update to
    `WHERE bet_hash = ? AND account = ?`.

## 11. Open questions

- **Remove a single chip vs remove all chips?** If you want to delete
  one ticket (e.g. miscategorized), today's "Placed" button toggles to
  un-place. Under multi-chip, should each chip get its own click-to-
  remove, or does the whole placed strip remove all at once? Suggest
  per-chip remove (each chip is a row, removal is row-scoped).
  Out of scope for v1 — flag for v1.1 if needed.
- **CLV tracking under multi-placement.** CLV is computed per `bet_hash`
  today. Under multi-placement, do we want per-chip CLV (each
  placement's own closing price diff) or one rolled-up CLV for the
  logical bet? Probably per-chip is the right answer (sometimes you
  place at different prices), but the existing CLV pipeline assumes
  one row per bet. Out of scope for v1 — document the asymmetry.
- **Header pill switch UX during re-open.** Should the dashboard
  highlight WZ accounts that don't have chips yet (green border on
  WagerzonJ, dim WagerzonC if it's been placed)? Nice-to-have, not
  required. Out of scope for v1.

---
