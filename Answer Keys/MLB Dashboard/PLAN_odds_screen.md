# MLB Bets Tab → Odds Screen

**Branch:** `claude/mlb-sportsbook-comparison-98fM8`
**Status:** Plan, not yet implemented
**Authored:** 2026-05-09

## Goal

Replace the current bets tab in the MLB dashboard with an odds-screen layout
that shows every bet alongside the same line+price at every other tracked
sportsbook (bettable + reference), so the user can verify the model's edge in
one glance instead of manually opening DraftKings to cross-check.

The model's recommendation (Pick book + side + price + EV + Kelly size) and
the placement flow stay identical to today. Only the *display* changes.

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

Layout — two rows per bet (both sides of the market), one block per
`(game × market × line the model picked)`:

```
┌────────────────────┬─────────┬──────┬─────┬─────┬──────────────┬────────────────────────────────────────┬───────┐
│ Game               │ Market  │ Side │ EV% │ Size│ Pick         │ WZ   H88  BFA  BKM  B105  DK   FD  Pin │ Place │
├────────────────────┼─────────┼──────┼─────┼─────┼──────────────┼────────────────────────────────────────┼───────┤
│ NYY @ BOS  7:05 PM │ Spread  │ NYY  │     │     │              │ -1.5 -1.5 -1.5 -2   -1.5 -1.5 -1.5 -1.5│       │
│                    │         │      │ 4.2 │ $42 │ ★ H88  +115  │ +110 ▓+115▓+108 -120 +112 -110 -108-105│ [Bet] │
│                    │         │ BOS  │     │     │              │ +1.5 +1.5 +1.5 +2   +1.5 +1.5 +1.5 +1.5│       │
│                    │         │      │     │     │              │ -130 -135 -128 +100 -132 -130 -132 -135│       │
└────────────────────┴─────────┴──────┴─────┴─────┴──────────────┴────────────────────────────────────────┴───────┘
```

- **Highlight:** picked book's cell on the picked side gets a green/yellow
  background (style `(c)` from design discussion).
- **Pick column:** explicit `★ {Book}  {line} / {odds}` so a long slate
  scans easily.
- **Line mismatch:** when a book doesn't quote the model's exact line, show
  that book's nearest line + its price for the same side, with the line
  rendered in amber so the mismatch is visually obvious.
- **Missing quote:** book has no line on that side at all → `—`.
- **Column order:** by-account books first (WZ, H88, BFA, BKM, B105),
  reference books after (DK, FD, Pinnacle). No visual divider.
- **Place button:** unchanged behavior. Single button per bet, places at the
  Pick book at the Pick price.

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
