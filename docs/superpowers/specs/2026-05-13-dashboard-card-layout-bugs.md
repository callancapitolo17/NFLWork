# Dashboard Card-Layout Bugs — Root-Cause Analysis & Follow-up Plan

**Date:** 2026-05-13
**Context:** During the singles-scrapers smoke session, we surfaced multiple
pre-existing bugs in the odds-screen card layout (merged in via `6056fd5`).
Some were fixed in-line; others remain. This doc captures the trace and a
plan for the cleanup.

---

## Bugs FIXED in this session

### 1. ✅ Mixed-type NA columns crash reactable

**Symptom:** Empty bets-tab card area; "0 OF N" counter shows but no cards render.

**Root cause:** `cents` (Kalshi-only) and `placed_actual` (NA for unplaced bets)
are nominally `numeric` columns but `htmlwidgets/jsonlite` serializes `NA_real_`
as the string `"NA"`. When a column has BOTH numeric and string values
(`cents: [3.2, "NA"]`), the React `reactable` widget can't parse the column's
declared type and silently fails to render any rows.

**Fix:** `select(-any_of(c("cents", "placed_actual", "fill_diff")))` before
passing `table_data` to `reactable()`. The new card layout doesn't need any
of those columns — they were legacy carriers from `create_bets_table_legacy()`.

**Commits:** `8889b74`, `70537ad`

### 2. ✅ Bets-tab book filter hides all rows

**Symptom:** Even with "All Books" selected, BOOK filter shows "4 selected"
and all 22 rows hidden by `applyFilters()`.

**Root cause:** The JS reads book from `cells[cells.length - 2]` —
second-to-last cell. In the *legacy* table layout that was a single-book
cell. In the *new card layout* that position is `pickside_html` — the
entire pill row text ("Over WZ -120 H88 — BFA — ..."). The book filter
then compares that long string against `activeFilters.book` (which contains
canonical book keys like `"wagerzon"`); never matches; every row hidden.

**Fix:** Read pick book from the Place/Log button's `data-book` attribute
instead. (Already present on the button per the existing cell-action HTML.)

**Commit:** `e09ddec`

### 3. ✅ DK/FD canonical rows had wrong market-name convention

**Symptom:** Zero DK and FD pills in `mlb_bets_book_prices` for any bet
despite 500+ raw DK records and 400+ FD records loaded.

**Root cause:** `get_dk_odds()` / `get_fd_odds()` emitted `market="main"` /
`"alternate_spreads"` / `"alternate_totals"` with `period` as a separate
column. But `normalize_book_odds_frame()` derives `period` from the market
string via `.derive_period(market_name)` (regex `_1st_[357]_innings$`).
"main" doesn't match the regex → period defaults to "FG" → `expand_bets_to_book_prices`'s
filter (`period == bet$period`) never matches an F5/F7 bet.

**Fix:** Add `.singles_market_name(market_type, period)` helper in `Tools.R`
that emits `"totals_1st_5_innings"`-style strings. Call from both `get_*_odds()`.

**Commit:** `5566784`

### 4. ✅ CSS pill row stacks vertically

**Symptom:** Each book pill renders on its own line in a narrow column
instead of side-by-side rows across the card width.

**Root cause:** `#bets-table-container .cell-pickside { flex-basis: 100%; }`
without `!important`. Reactable applies an inline `style="flex: 100 0 auto"`
on every `.rt-td` that wins over external CSS without `!important`. The cell
stays at content-width; the inner `.side-row` flex container has no horizontal
space; pills wrap vertically. The parlays-tab CSS already uses `!important`
on its equivalent rule — bets-tab missed it.

**Fix:** `flex-basis: 100% !important; width: 100% !important;` plus tighter
selector `.rt-td.cell-*`.

**Commit:** `f3d6875`

### 5. ✅ Line tolerance too strict for F-period totals display

**Symptom:** DK posts F7 totals at 7.5 but the model bets Over 5.5; the pill
showed `—` because |5.5 - 7.5| = 2.0 exceeded `LINE_MATCH_TOLERANCE = 1.0`.

**Not actually a bug — UX choice.** The pill renderer (`book_pill.R`) already
paints amber + adds a `O7.5` line tag when `is_exact_line = FALSE`, so a wider
display tolerance just exposes more book quotes with clear visual indication
of the line gap.

**Fix:** `LINE_MATCH_TOLERANCE = 3.0` in `odds_screen.R`. Pick-book line is
still required to be exact by the data-bug guard in `MLB.R` (line 919) — this
only affects comparison pill display.

**Commit:** `f3d6875` (bundled with #4)

---

## Bugs IDENTIFIED but NOT yet fixed

### 6. ⚠️ DK F7 alt-line totals not making it to per-book DB

**Symptom (user verified in DK's UI on 2026-05-13):** DK's "Total Runs - 1st 7 Innings"
section shows Over/Under at lines 5.5 / 6.5 / 7.5. Our `dk_odds/dk.duckdb` has
only ONE F7 total per game (the "main" line, e.g. 7.5 for Texas/Arizona).

**Hypothesis:** Multiple alt lines come back as separate `Total Alternate - 1st 7 Innings`
markets in DK's parlays API (one market per alt line, each with one Over/Under pair).
Either:

  - (a) `classify_market` isn't recognizing today's market names (verified our captured
        fixture only had ONE alt market for F7; live API today likely has multiple), OR
  - (b) the scraper IS capturing them but the parser/DB-write step coalesces them, OR
  - (c) DK actually returns them under one umbrella market with multiple selections
        and our parser drops all but one.

**Investigation needed:** capture today's live DK parlays payload for a game with
upcoming start time, dump all F7 markets + their selections, verify how DK is
structuring the data right now (the fixture is from 2026-05-12 and the live
shape may differ).

**Impact:** Most F-period totals bets show only ONE DK line in the dashboard
pill row. With the line-tolerance fix (#5) at 3.0, that one line shows even
when far from the model's line — but multiple alt lines are missing entirely.

### 7. ⚠️ Spread line tag shows "O-1.5" or "0-1.5" instead of just "-1.5"

**Symptom (user screenshot of Kansas City Royals card):** WZ pill on a spread
bet shows `0-1.5` instead of the expected `-1.5`.

**Hypothesis:** Memory note from 2026-05-11 (`mlb_odds_screen_pending_merge.md`)
says `commit b2c6b3b` was supposed to fix the spread tag (no O/U prefix). Either:

  - (a) That fix didn't propagate correctly into `book_pill.R` (which always
        prepends `O` or `U`), OR
  - (b) The `is_totals` parameter isn't being passed through `render_book_pill()`
        when called from `create_bets_table()`.

**Investigation needed:** read `book_pill.R` and the call site in `mlb_dashboard.R`;
verify `is_totals_market` is computed and passed; check whether the tag suppression
logic exists.

**Impact:** Cosmetic — odds are correct, just the line-tag label is wrong on
spread pills.

---

## Plan for follow-up session

### Step 1: Investigate bug #6 (F7 alt totals missing)

Capture a fresh DK parlays payload for an upcoming game (one that hasn't tipped
off and DK hasn't pulled markets yet). Dump all `1st 7 Innings` markets +
selections. Compare to today's `dk_odds/dk.duckdb` content for the same game.

Decision tree:
- If DK returns multiple "Total Alternate - 1st 7 Innings" markets → debug
  `classify_market` / `parse_selections_to_wide_rows` for why they're dropped.
- If DK returns one market with multiple Over/Under pairs at different lines →
  fix the parser to bucket totals by line even when `market_type == "main"`.

### Step 2: Investigate bug #7 (spread line tag format)

Read `book_pill.R` and `mlb_dashboard.R` call site. Confirm whether `is_totals`
flag is plumbed and respected. Add it if missing.

### Step 3: Re-render dashboard + visual smoke

After both fixes, regenerate `report.html`, hard-refresh browser, verify:
- F7 totals cards show DK at multiple lines (e.g., 5.5 / 6.5 / 7.5)
- Spread cards show line tag as `-1.5` (no `O` or `0` prefix)

### Step 4: Pre-merge gate (revised)

Update the pre-merge review in `2026-05-12-mlb-dk-fd-singles-scrapers-design.md`
to cover all 7 bugs discovered today (5 fixed, 2 deferred to that follow-up
session, OR fixed during this follow-up).

### Step 5: Merge with user approval

Per CLAUDE.md feedback_always_ask_merge.

---

## Architectural takeaway

The odds-screen card layout was merged with at least 7 latent bugs, all surfacing
only when the new singles-scrapers data flowed in. The card layout's pre-existing
"pending in-browser smoke" note understated how much polish remained. Future
multi-feature branches with cross-cutting display changes should get a more
thorough visual smoke per-feature before merge.
