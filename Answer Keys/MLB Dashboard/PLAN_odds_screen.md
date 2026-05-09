# MLB Bets Tab → Odds Screen

**Branch:** `claude/mlb-sportsbook-comparison-98fM8`
**Status:** Plan, not yet implemented
**Authored:** 2026-05-09

---

## Quick start (terminal)

```bash
# 1. Get the branch locally
cd /Users/callancapitolo/NFLWork
git fetch origin claude/mlb-sportsbook-comparison-98fM8
git checkout claude/mlb-sportsbook-comparison-98fM8

# 2. View / edit the plan
less "Answer Keys/MLB Dashboard/PLAN_odds_screen.md"
# or open in your editor:
code "Answer Keys/MLB Dashboard/PLAN_odds_screen.md"
open -a "Visual Studio Code" "Answer Keys/MLB Dashboard/PLAN_odds_screen.md"

# 3. Key files this plan touches (open all three to start)
code "Answer Keys/MLB Answer Key/MLB.R"
code "Answer Keys/MLB Dashboard/mlb_dashboard.R"
code "Answer Keys/MLB Answer Key/Tools.R"

# 4. Run the pipeline (writes mlb_bets_combined + new mlb_bets_book_prices)
cd "Answer Keys/MLB Answer Key"
Rscript MLB.R                  # ~3-8 min depending on slate size

# 5. Start the dashboard (port 8083)
cd "/Users/callancapitolo/NFLWork/Answer Keys/MLB Dashboard"
./run.sh                        # serves http://localhost:8083
# or directly:
python3 mlb_dashboard_server.py

# 6. Inspect the new table after a pipeline run
duckdb "/Users/callancapitolo/NFLWork/Answer Keys/mlb_mm.duckdb" \
  "SELECT * FROM mlb_bets_book_prices LIMIT 20"

# 7. View the GitHub branch in browser (if PR helpful)
# https://github.com/callancapitolo17/NFLWork/tree/claude/mlb-sportsbook-comparison-98fM8
```

**Path note:** This plan was authored from a sandboxed environment at
`/home/user/NFLWork`. Your local clone is at `/Users/callancapitolo/NFLWork`.
All file paths in this plan are repo-relative — prepend either prefix
depending on where you're running.

---

## Goal

Replace the current bets tab in the MLB dashboard with an odds-screen layout
that shows every bet alongside the same line+price at every other tracked
sportsbook (bettable + reference), so the user can verify the model's edge in
one glance instead of manually opening DraftKings to cross-check.

The model's recommendation (Pick book + side + price + EV + Kelly size) and
the placement flow stay identical to today. Only the *display* changes.

## Why this is a low-risk change

The cross-book comparison data is **already computed every pipeline run**
and immediately discarded. `MLB.R:778-792` builds per-book bets frames
(`wagerzon_bets`, `hoop88_bets`, `bfa_bets`, `bookmaker_bets`,
`bet105_bets`), unions them, then `group_by(id, base_market, bet_on) %>%
filter(ev == max(ev)) %>% slice_head(n=1)` collapses them to one row per
market keeping only the best book. The dashboard never sees the
discarded rows.

This plan stops discarding them — that's the whole data-side change. The
display side is a card-layout rebuild of `create_bets_table()` mirroring
the existing parlays-tab pattern.

Implication: nothing in the placement flow, EV math, Kelly sizing, or
pipeline cadence changes. The risk surface is contained to (a) one new
DuckDB table and (b) one rebuilt reactable in the dashboard.

## Non-goals (v1)

- No per-cell place buttons. Single Place button per bet, places at Pick.
- No alt ladder. Each bet block anchors on the model's exact line; books
  that don't quote that line show their nearest available line.
- No EV recompute at non-pick prices.
- No new tab. This *replaces* the existing bets tab.
- No changes to the parlays tab, trifectas tab, RFQ bot, scrapers, or
  pipeline cadence.
- No model fair column.

## Final design (locked)

**Card layout, mirroring the parlays tab.** Each bet renders as a card
(reactable row with `display:flex` + per-cell `order:` overrides), not as
a flat table row. Inside the card the visual stack is:

1. **Game header** (full width) — matchup + date/time
2. **Market header** (full width) — e.g., `Spread NYY -1.5`
3. **Pick-side pills row** (full width) — side label + one pill per book
4. **Other-side pills row** (full width) — side label + one pill per book
5. **Metadata strip** (inline-flex) — `Pick / EV / Size / To Win / [Place]`

```
┌────────────────────────────────────────────────────────────────────────────┐
│ NYY @ BOS  Fri 7:05 PM                                                     │
│ Spread NYY -1.5                                                            │
│ NYY -1.5:  [WZ +110] [▓H88 +115 ★▓] [BFA +108] [BKM ⟨-2⟩ -120] [B105 +112] │
│            [DK -110] [FD -108] [Pinn -105]                                 │
│ BOS +1.5:  [WZ -130] [H88 -135] [BFA -128] [BKM ⟨+2⟩ +100] [B105 -132]    │
│            [DK -130] [FD -132] [Pinn -135]                                 │
│ Pick H88 +115   EV 4.2%   Size $42   To Win $48                   [Place]  │
└────────────────────────────────────────────────────────────────────────────┘
```

**CSS pattern** mirrors parlays directly (`mlb_dashboard.R:1907-2070`):
- Row container `display:flex; flex-wrap:wrap`.
- Full-width cells (`cell-game`, `cell-market`, `cell-pickside`,
  `cell-otherside`) carry `flex-basis:100%`.
- Metadata cells (`cell-pick`, `cell-ev`, `cell-size`, `cell-towin`,
  `cell-action`) are `display:inline-flex` and flow horizontally.
- Each metadata cell has its label baked in via `::before` pseudo-element
  (e.g., `.cell-pick::before { content: "Pick" }`) — same pattern as
  parlays' Fair / WZ / Size / To Win labels.
- Place button pushed right via `margin-left:auto` on the preceding flex
  item, identical to parlays' Edge → Action pattern.
- `order:` overrides give the visual stack independent of the dataframe
  column order reactable preserves.
- Same responsive override block (`min-width:0 !important; flex:0 1 auto`)
  so cards reflow at narrow widths.

**Pill rendering** — new helper `render_book_odds_strip(side_label, model_line, books)`
in `mlb_dashboard.R`, analogous to the existing `render_books_strip`:
- Always renders all 8 books in fixed order (WZ, H88, BFA, BKM, B105, DK,
  FD, Pinn) so the eye learns positions across cards.
- Exact-line pill: `BookCode  +Odds` (e.g., `H88  +115`).
- Mismatched-line pill: `BookCode  ⟨line⟩  Odds` with the bracketed line
  in amber (e.g., `BKM  ⟨-2⟩  -120`).
- Missing-quote pill: `BookCode  —` in muted grey.
- Pick pill on the pick side: green/yellow background tint + ★ glyph,
  reusing the `cell-pick-highlight` class.

**Place button:** unchanged behavior. Single `[Place]` per card in the
metadata strip, places at the Pick book at the Pick price. Same JS hooks
and `data-*` attributes as today.

## Code trace (file paths + line numbers)

These are the load-bearing locations a fresh session would need to
re-discover. Captured here so the implementer doesn't need to spelunk.

### Dashboard (the consumer)

| What | File | Line(s) | Notes |
|---|---|---|---|
| Bets table loader | `Answer Keys/MLB Dashboard/mlb_dashboard.R` | 4476-4485 | `dbGetQuery(con, "SELECT * FROM mlb_bets_combined")` from `mlb_mm.duckdb` |
| `create_bets_table()` | `Answer Keys/MLB Dashboard/mlb_dashboard.R` | search for `create_bets_table` | Builds the reactable for the bets tab — this gets replaced |
| Parlays card layout (CSS template to mirror) | `Answer Keys/MLB Dashboard/mlb_dashboard.R` | 1907-2070 | `display:flex` + `flex-basis:100%` + `order:` + `::before` labels + responsive overrides |
| Parlay `render_books_strip()` (helper to copy/adapt) | `Answer Keys/MLB Dashboard/mlb_dashboard.R` | called at 543-552; defined elsewhere — grep `render_books_strip` | Renders the per-book pill row for parlays |
| Parlays table builder (column order + cell-classes) | `Answer Keys/MLB Dashboard/mlb_dashboard.R` | ~440-570 | `create_parlays_table()` — reference for the cell class assignments |
| Place button JS hooks | `Answer Keys/MLB Dashboard/mlb_dashboard.R` | search for `data-game-id`, `data-book`, `placeBet(` | Keep these data attributes on the pick-side row |
| Server (Flask) | `Answer Keys/MLB Dashboard/mlb_dashboard_server.py` | 2334 lines | `/api/place-bet`, `/api/place-parlay`, `/api/wagerzon/*` — **no changes for this plan** |
| `run_closing_capture()` (proves cross-book read pattern) | `Answer Keys/MLB Dashboard/mlb_dashboard_server.py` | 548 | Already shops per-book DBs for closing snapshots — same pattern to reuse |

### Pipeline (the producer)

| What | File | Line(s) | Notes |
|---|---|---|---|
| Per-book bets construction | `Answer Keys/MLB Answer Key/MLB.R` | 585-770 | `wagerzon_bets`, `wz_alt_bets`, `hoop88_bets`, `bfa_bets`, `bfa_alt_bets`, `bookmaker_bets`, `bkm_alt_bets`, `bet105_bets`, `bet105_alt_bets`, `kalshi_bets` |
| Comparison helpers (already produce per-book rows) | `Answer Keys/MLB Answer Key/Tools.R` | search `compare_spreads_to_wagerzon`, `compare_totals_to_wagerzon`, `compare_moneylines_to_wagerzon`, `compare_alts_to_samples` | These return per-book bet rows; we just need to keep them |
| **Dedup that throws away comparison data** | `Answer Keys/MLB Answer Key/MLB.R` | 783-787 | `group_by(id, base_market, bet_on) %>% filter(ev == max(ev)) %>% slice_head(n=1)` — **insertion point for the new expansion is right before this** so we still have all rows |
| `mlb_bets_combined` write | `Answer Keys/MLB Answer Key/MLB.R` | 846-857 | `DROP TABLE IF EXISTS` + `dbWriteTable` to `mlb_mm.duckdb`. Same pattern for the new table. |
| Odds API frame (DK, FD, Pinn, etc. live here) | `Answer Keys/MLB Answer Key/MLB.R` | 155 | `all_books <- unique(game_odds$bookmaker_key)` |
| Reference book key list (full universe) | `Answer Keys/MLB Answer Key/Back Testing Odds.R` | 141-144 | `pinnacle, circasports, betonlineag, lowvig, matchbook, bookmaker.eu, coolbet, everygame, intertops, unibet, ..., betmgm, ..., caesars, fanduel, draftkings, ...` |
| Sharp consensus selection | `Answer Keys/MLB Answer Key/MLB.R` | 127-183 | `SHARP_BOOKS` filter for the consensus model — not used for the new table but useful context |

### Team name resolution (already canonical, no work needed)

| What | File |
|---|---|
| Python → Odds API canonical | `Answer Keys/canonical_match.py` |
| R → Odds API canonical | `resolve_offshore_teams()` in `Answer Keys/MLB Answer Key/Tools.R` |

All scrapers already write canonical names. The new join key
`(game_id, market_type, period, side, line)` is unambiguous across books.

### Existing helpers to reuse (don't reinvent)

| Helper | Why useful here |
|---|---|
| `render_books_strip()` (parlays) | Direct template for `render_book_odds_strip()` — copy structure, swap content (line+price instead of devigged %) |
| `compare_*_to_wagerzon()` family | Already produces per-book bet rows with `prob`, `ev`, `bet_size` — no new pricing math needed |
| `compare_alts_to_samples()` | Same, for alt lines |
| `apply_combo_residuals()` (parlays) | Not reusable — parlay-specific |
| `canonical_match.py` / `resolve_offshore_teams()` | Names already resolved upstream — just join |

### Conversation history (design decisions and what was rejected)

The card layout was the third design considered. Trail of decisions:

1. **First proposal — "show other books for each bet":** Three options:
   - **A**: persist pre-dedup frame, expandable row in reactable. Cheap.
   - **B**: live on-demand `/api/compare?...` endpoint hitting per-book DBs. Live but adds lock contention.
   - **C**: precomputed best/median/devigged-fair summary columns inline.
   - User reframed: *trust* in model is the issue, not shopping for prices.

2. **Second proposal — three layouts** mirroring the trust framing:
   - **A flat**: books as columns, one row per recommended bet.
   - **B flat**: both sides shown per market (true odds-screen). User picked B.
   - **C flat**: today's table + click-to-expand mini odds-screen.

3. **Highlight options for the picked book:**
   - **(a)** cell-only highlight, **(b)** cell + "PICK" badge, **(c)** cell + dedicated Pick column, **(d)** side-row stripe + cell.
   - User picked **(c)** — cell + Pick column.

4. **Place button:** initially considered click-any-cell-to-place. User
   chose: **single Place button at the recommended Pick** (today's
   behavior). Per-cell place is rejected for v1.

5. **Reference vs bettable books:** user clarified the actual workflow —
   model says +250 on WZ on Giants -2.5; user manually opens DraftKings
   to verify. The verification book is a sharp/mainstream book (DK / FD /
   Pinnacle), not just "another offshore book." Settled on: 5 bettable
   (WZ, H88, BFA, BKM, B105) + 3 reference (DK, FD, Pinn) = 8 books total.

6. **Final layout switch — flat → card:** user noted the parlays tab uses
   a richer "dynamic format" — not a flat reactable row. Found the
   parlays-tab card pattern at `mlb_dashboard.R:1907-2070`. Adopted it
   wholesale: same `display:flex` + `order:` + `::before` labels +
   responsive overrides. Now the bets tab visually matches the parlays
   tab convention.

7. **Dropped from v1:** model fair column, bettable/reference visual
   divider (Pick column already conveys "where to bet"), per-cell place,
   live `/api/compare` endpoint, alt ladder.

## Data sources (already exist)

| Source | Path | Books |
|---|---|---|
| `wagerzon.duckdb` → `mlb_odds` | per-book scraper DBs | WZ |
| `hoop88.duckdb` → `mlb_odds` | | H88 |
| `bfa.duckdb` → `mlb_odds` | | BFA |
| `bookmaker.duckdb` → `mlb_odds` | | BKM |
| `bet105.duckdb` → `mlb_odds` | | B105 |
| `game_odds` frame in `MLB.R` | from Odds API, already loaded at `MLB.R:155` | DK, FD, Pinn (and others available if needed) |

All of these already use canonical Odds API team names (resolved by
`canonical_match.py`) and a shared `(game_id, market, bet_on, line, period)`
key. No new normalization needed.

## Schema — new table `mlb_bets_book_prices`

Long format, written to `mlb_mm.duckdb` alongside `mlb_bets_combined`. One
row per `(bet × book × side)`:

| Column | Type | Notes |
|---|---|---|
| `bet_row_id` | VARCHAR | Stable hash of `(game_id, market, line, bet_on)` from the parent bet — joins back to `mlb_bets_combined` |
| `game_id` | VARCHAR | Odds API id |
| `market` | VARCHAR | Normalized: `spreads`, `totals`, `moneyline`, plus alt variants |
| `period` | VARCHAR | `FG`, `F5`, `F3` |
| `side` | VARCHAR | `pick` or `opposite` (relative to the bet the model recommended) |
| `bookmaker` | VARCHAR | `wagerzon`, `hoop88`, `bfa`, `bookmaker`, `bet105`, `draftkings`, `fanduel`, `pinnacle` |
| `line_quoted` | DOUBLE | The line *this book* offers on this side (may differ from model's line) |
| `american_odds` | INTEGER | Price at `line_quoted` |
| `is_exact_line` | BOOLEAN | True if `line_quoted == model_line` |
| `fetch_time` | TIMESTAMP | Source row's `fetch_time` |

Grain: `(bet_row_id, bookmaker, side)`. For a slate of ~250 bets × 8 books ×
2 sides ≈ 4000 rows. Trivial for DuckDB.

## Implementation

### Phase 1 — pipeline (`MLB.R`)

**File:** `Answer Keys/MLB Answer Key/MLB.R`
**Insertion point:** between line 792 (end of `slice_head(n=1)` dedup) and
line 820 (`adjust_kelly_for_correlation`).

1. Add `bet_row_id` column to `all_bets_combined` via deterministic hash of
   `(game_id, market, line, bet_on)`. Same hash logic gets used by the
   expansion below, so bets and per-book rows share the join key.
2. Build `book_prices_long`:
   - For each row in `all_bets_combined`, expand into one row per
     `(bookmaker, side ∈ {pick, opposite})`.
   - For offshore books, source from per-book `mlb_odds` tables already
     loaded in memory (`wagerzon_odds`, `hoop88_odds`, `bfa_odds`,
     `bookmaker_odds`, `bet105_odds`).
   - For DK / FD / Pinn, source from `game_odds` filtered to those
     `bookmaker_key`s.
   - Join key: `(game_id, market_type, period, side)` first, then within
     that subset pick the row with `line == model_line` if present;
     otherwise the row minimizing `abs(line - model_line)` — tiebreaker:
     prefer the line worse for the bettor.
3. Write `mlb_bets_book_prices` to `mlb_mm.duckdb` in the same writer block
   as `mlb_bets_combined` (`MLB.R:851-857`). Drop-and-replace, same as the
   bets table.

**New helper** in `Answer Keys/MLB Answer Key/Tools.R` (or a new file
`Answer Keys/MLB Answer Key/odds_screen.R` if Tools.R is getting heavy):

```r
expand_bets_to_book_prices <- function(bets, book_odds_by_book) { ... }
```

Pure function — takes the bets frame + a named list of per-book odds
frames, returns the long-format table. Easy to unit-test.

### Phase 2 — dashboard (`mlb_dashboard.R`)

**File:** `Answer Keys/MLB Dashboard/mlb_dashboard.R`

1. **Loader (~line 4476):** load `mlb_bets_book_prices` from `mlb_mm.duckdb`
   alongside `mlb_bets_combined`. If missing, fall back to today's display
   (defensive — don't break the dashboard during deploy).
2. **Pivot** book_prices long → wide on `(bet_row_id, side)`. Resulting
   columns per side: `wz_line`, `wz_odds`, `wz_is_exact`, …, `pinn_odds`,
   `pinn_is_exact`. Two rows per bet (one per side) joined back to the bet's
   metadata (game, market, EV, Size, Pick).
3. **Replace `create_bets_table()`** (currently single-row-per-bet) with
   new column definitions matching the locked layout.
4. **Cell renderer for each book column:** small reactable JS function that
   stacks line on top, odds below; greys out the whole cell if line+odds
   are NA; amber-colors the line text if `*_is_exact == FALSE`; green-tints
   the cell background when `(book == pick_book) && (side == pick_side)`.
5. **Place button column:** only render `[Bet]` on the pick-side row of each
   bet block. Empty cell on the opposite-side row. Reuse existing JS data
   attributes (`data-book`, `data-odds`, `data-line`, etc.) so the
   placement endpoint is unchanged.
6. **CSS additions** in the inline style block: `.cell-pick-highlight`,
   `.cell-line-mismatch`, `.cell-no-quote`. Match existing dashboard
   palette (greys/blues, no neon).

### Phase 3 — server (`mlb_dashboard_server.py`)

No changes. The placement payload is identical to today.

## Implementation order (suggested)

Smallest-blast-radius commits first, so the tree stays runnable
throughout:

1. **Helper + test (no behavior change yet):**
   - Add `expand_bets_to_book_prices()` to `Tools.R` (or new
     `odds_screen.R`).
   - Hand-call it from a scratch R session to verify output schema on a
     real `mlb_bets_combined` snapshot before wiring it into `MLB.R`.
   - Commit: `feat(mlb): helper to expand bets to per-book prices`.

2. **Pipeline write (data flowing, dashboard still old):**
   - Wire `expand_bets_to_book_prices()` into `MLB.R` between
     line 792 and line 820.
   - Add `bet_row_id` hash to `all_bets_combined` *before* the dedup so
     both tables share the join key.
   - Add the DROP+WRITE block alongside the existing `mlb_bets_combined`
     write at lines 851-857.
   - Run the pipeline once. Confirm `mlb_bets_book_prices` populated.
   - Commit: `feat(mlb): write mlb_bets_book_prices in pipeline`.

3. **Dashboard rebuild (behavior change visible):**
   - Add `render_book_odds_strip()` helper in `mlb_dashboard.R`.
   - Add `cell-pick`, `cell-ev`, `cell-size`, `cell-towin`,
     `cell-action`, `cell-pickside`, `cell-otherside`, `cell-market`,
     `cell-pick-highlight`, `cell-line-mismatch`, `cell-no-quote` CSS
     classes inside the existing inline style block. Mirror the
     parlays selectors at `mlb_dashboard.R:1907-2070`.
   - Replace `create_bets_table()` body — load
     `mlb_bets_book_prices`, pivot long→wide, expand to two rows per
     bet, render via the new column defs.
   - Defensive try/catch in the loader so the dashboard falls back to
     today's table if `mlb_bets_book_prices` is missing.
   - Test in a browser. Verify Place button still works on standard +
     alt bets.
   - Commit: `feat(mlb-dashboard): card-layout odds-screen for bets tab`.

4. **Documentation:**
   - Update `Answer Keys/MLB Dashboard/README.md` with the new card
     layout (paste the ASCII sketch from this plan), explanation of
     pill colors (highlight / amber / muted-grey), and note the new
     `mlb_bets_book_prices` table.
   - Update `Answer Keys/CLAUDE.md` Database section to list the new
     table.
   - Commit: `docs(mlb-dashboard): document odds-screen layout`.

5. **Pre-merge review** per checklist below; user explicit approval;
   merge to `main`.

## How to verify (post-implementation)

After implementing and running `Rscript MLB.R`:

```bash
# 1. New table exists with expected schema
duckdb "Answer Keys/mlb_mm.duckdb" \
  "DESCRIBE mlb_bets_book_prices"

# 2. Row count looks reasonable (~bets × 8 books × 2 sides)
duckdb "Answer Keys/mlb_mm.duckdb" \
  "SELECT COUNT(*) AS rows,
          COUNT(DISTINCT bet_row_id) AS bets,
          COUNT(DISTINCT bookmaker) AS books
   FROM mlb_bets_book_prices"

# 3. Spot-check one bet across all 8 books, both sides
duckdb "Answer Keys/mlb_mm.duckdb" \
  "SELECT bookmaker, side, line_quoted, american_odds, is_exact_line
   FROM mlb_bets_book_prices
   WHERE bet_row_id = (SELECT bet_row_id FROM mlb_bets_book_prices LIMIT 1)
   ORDER BY side, bookmaker"

# 4. Sanity: every bet in mlb_bets_combined should have at least one
#    row in mlb_bets_book_prices (the pick-book pick-side row at minimum)
duckdb "Answer Keys/mlb_mm.duckdb" \
  "SELECT c.bet_row_id, COUNT(p.bookmaker) AS book_rows
   FROM mlb_bets_combined c
   LEFT JOIN mlb_bets_book_prices p USING (bet_row_id)
   GROUP BY 1
   HAVING book_rows = 0"
# (should return zero rows)
```

In the dashboard:
- Open http://localhost:8083, switch to bets tab.
- Each card shows: game header, market header, two pill rows, metadata
  strip with Pick / EV / Size / To Win / Place.
- Pick book pill is highlighted with ★.
- At least one card should show a mismatched line (BKM frequently hangs
  alt totals at half-run intervals different from main books) — verify
  amber rendering.
- At least one card should show `—` for some books (Pinnacle doesn't
  always quote alts, BFA misses team totals, etc.) — verify muted grey
  rendering.
- Click Place on a card → bet appears in placed_bets. Same flow as today.

## Rollback

If anything breaks after merging:

```bash
# Revert the dashboard rebuild commit (data writes can stay — they're harmless)
git revert <dashboard-commit-sha>
git push origin main
# Restart the dashboard
cd "Answer Keys/MLB Dashboard"
./run.sh
```

If the pipeline write itself is causing issues (DB lock, schema mismatch):

```bash
git revert <pipeline-commit-sha>
# Drop the orphaned table
duckdb "Answer Keys/mlb_mm.duckdb" "DROP TABLE IF EXISTS mlb_bets_book_prices"
```

The defensive try/catch in the dashboard loader means it falls back to
today's table automatically if `mlb_bets_book_prices` is missing — so a
partial rollback (revert dashboard only, leave pipeline) is also safe.

## Version control

- **Branch:** `claude/mlb-sportsbook-comparison-98fM8` (already current).
- **Worktree:** session is already running directly on the branch in the
  primary working directory; no separate worktree needed.
- **Commits, in order:**
  1. `feat(mlb): expand per-book odds for each bet in MLB.R` — adds the
     helper, the expansion logic, the new table write.
  2. `feat(mlb-dashboard): odds-screen layout for bets tab` — replaces
     `create_bets_table`, adds CSS, adds loader for `mlb_bets_book_prices`.
  3. `docs(mlb-dashboard): document odds-screen bets tab` — README updates
     (see Documentation section below).
- **No merge to `main`** until: pipeline runs end-to-end with the new
  table, dashboard renders without errors on a real slate, executive review
  per `CLAUDE.md` complete, user explicitly approves.

## Documentation

Files to update in the same PR as the code:

- `Answer Keys/MLB Dashboard/README.md` — update bets-tab section with
  screenshot/sketch of new layout, list of book columns, explanation of
  highlight + amber tint conventions, explanation that nearest-line
  fallback is used when exact line isn't quoted.
- `Answer Keys/MLB Answer Key/CLAUDE.md` (if it exists; check) — note that
  `MLB.R` now writes `mlb_bets_book_prices` to `mlb_mm.duckdb` in addition
  to `mlb_bets_combined`.
- Root `CLAUDE.md` — no changes expected; the change is contained inside
  the MLB dashboard subsystem.

## Edge cases / decisions

- **Doubleheaders:** `game_id` is unique per game, so the join key handles
  this naturally. No special logic.
- **Period:** F5 / F3 / FG are filtered as part of the join key — a bet on
  F5 spread won't accidentally pull FG odds from another book.
- **Markets the Odds API doesn't track** (e.g., team totals): reference
  columns show `—`. Bettable columns still populated from per-book scraper
  DBs.
- **Book outage / missed scrape:** column shows `—` for affected games.
- **Alts in the new view:** alt bets get their own block in the bets table
  (already the case today). Within each block, books that don't offer that
  exact alt show their nearest available line on the same side.
- **First run after deploy with no `mlb_bets_book_prices` table yet:**
  dashboard loader has try/catch that falls back to today's column set,
  prints a warning. Pipeline run repopulates it within one cycle.

## Pre-merge executive review checklist

Per `CLAUDE.md` pre-merge requirements, before merging to `main`:

- [ ] `git diff main..HEAD --stat` reviewed
- [ ] No duplicate writes to `mlb_bets_book_prices` (DROP + WRITE pattern,
      same as `mlb_bets_combined`)
- [ ] DB connections all wrapped in `on.exit(dbDisconnect(...))`
- [ ] Off-season behavior: pipeline produces empty bets → empty
      `mlb_bets_book_prices` → dashboard renders empty bets tab without
      error
- [ ] Doubleheader behavior: two games same teams → distinct `bet_row_id`s
      → no row collision
- [ ] No new dependencies added to `requirements.txt` / `DESCRIPTION`
- [ ] Dashboard renders end-to-end with real pipeline output
- [ ] Place button still works for both standard and alt bets
- [ ] No secrets / API keys logged
- [ ] README updates committed in the same merge
- [ ] User explicit approval to merge

## Open questions (none blocking)

None. All design decisions resolved in conversation.
