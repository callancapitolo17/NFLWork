# MLB Bets Tab → Odds Screen

**Branch:** `claude/mlb-sportsbook-comparison-98fM8`
**Status:** Plan, not yet implemented
**Authored:** 2026-05-09
**Revised:** 2026-05-11 (added auto-placer, dropped nearest-line fallback, locked
fixed-width pill columns, added M pill in metadata strip, expanded Odds API
coverage, added worktree lifecycle)

---

## Quick start (terminal)

```bash
# 1. Get the branch locally
cd /Users/callancapitolo/NFLWork
git fetch origin claude/mlb-sportsbook-comparison-98fM8
git worktree add .worktrees/mlb-odds-screen claude/mlb-sportsbook-comparison-98fM8
cd .worktrees/mlb-odds-screen

# 2. View / edit the plan
less "Answer Keys/MLB Dashboard/PLAN_odds_screen.md"

# 3. Key files this plan touches
code "Answer Keys/MLB Answer Key/MLB.R"
code "Answer Keys/MLB Dashboard/mlb_dashboard.R"
code "Answer Keys/MLB Dashboard/mlb_dashboard_server.py"
code "Answer Keys/MLB Answer Key/Tools.R"

# 4. Run the pipeline (writes mlb_bets_combined + new mlb_bets_book_prices)
cd "Answer Keys/MLB Answer Key"
Rscript MLB.R                  # ~3-8 min depending on slate size

# 5. Start the dashboard (port 8083)
cd "/Users/callancapitolo/NFLWork/.worktrees/mlb-odds-screen/Answer Keys/MLB Dashboard"
./run.sh

# 6. Inspect the new table after a pipeline run
duckdb "/Users/callancapitolo/NFLWork/Answer Keys/mlb_mm.duckdb" \
  "SELECT * FROM mlb_bets_book_prices LIMIT 20"
```

---

## Goal

Replace the current bets tab in the MLB dashboard with an odds-screen layout
that shows every bet alongside the same line+price at every other tracked
sportsbook (bettable + reference), so the user can verify the model's edge in
one glance instead of manually opening DraftKings to cross-check. Bring the
Wagerzon direct-API auto-placer to singles (mirroring what parlays already
have), so a single click places the bet at Wagerzon without a browser.

The model's recommendation (Pick book + side + price + EV + Kelly size) is
unchanged. The layout, the per-card book comparison, and the WZ direct-API
placement path are new.

## Why this is a low-risk change

The cross-book comparison data is **already computed every pipeline run**
and immediately discarded. `MLB.R:778-792` builds per-book bets frames
(`wagerzon_bets`, `hoop88_bets`, `bfa_bets`, `bookmaker_bets`,
`bet105_bets`), unions them, then `group_by(id, base_market, bet_on) %>%
filter(ev == max(ev)) %>% slice_head(n=1)` collapses them to one row per
market keeping only the best book. The dashboard never sees the
discarded rows.

This plan stops discarding them — that's the data-side change. The
display side is a card-layout rebuild of `create_bets_table()` mirroring
the parlays-tab pattern. The auto-placer side reuses the existing
Wagerzon parlay placer's session/preflight/drift code paths (already
production-hardened) — singles are a narrower call than parlays so the
risk surface is strictly smaller.

Implication: nothing in the EV math, Kelly sizing, scraper cadence, or
parlay/trifecta flows changes. The risk surface is contained to (a) one
new DuckDB table, (b) one rebuilt reactable in the dashboard, and (c) a
dispatch branch in `/api/place-bet` that calls the existing WZ placer
when the pick book is Wagerzon.

## Non-goals (v1)

- No alt ladder. Each bet block anchors on the model's exact line; books
  that don't quote that line render `—` (no nearest-line fallback).
- No EV recompute at non-pick prices.
- No new tab. This *replaces* the existing bets tab.
- No changes to the parlays tab, trifectas tab, RFQ bot, or scrapers.
- No dry-run toggle. `[Place]` is always real placement (parlays have
  dry-run; the simpler bets path doesn't need it).
- No H88 / BFA / BetOnline direct-API placement. Those books keep the
  existing Playwright browser flow. WZ is the only direct-API book v1.
- No DK / FD / Pinnacle / BKM / B105 placement. These are reference and
  sharp books — read-only on the bets tab.
- No replacement of the placed-bets summary table at the top of the
  bets tab (stays flat).

## Final design (locked)

**Card layout, mirroring the parlays tab with fixed-width pill columns.**
Each bet renders as a card (reactable row with `display:flex` + per-cell
`order:` overrides), not as a flat table row. Inside the card the visual
stack is:

1. **Game header** (full width) — matchup + date/time
2. **Market header** (full width) — e.g., `Spread F5 · NYY -1.5`
3. **Pick-side pills row** (full width) — side label + one pill per book
4. **Opposite-side pills row** (full width) — side label + one pill per book
5. **Metadata strip** (inline-flex) — `M / Pick / EV / Size / To Win / [Place] [Log]`

```
┌─────────────────────────────────────────────────────────────────────────┐
│ NYY @ BOS  Fri 7:05 PM PT                                               │
│ Spread F5 · NYY -1.5                                                    │
│                                                                         │
│ NYY -1.5  [WZ +110] [H88 +115] [BFA +108] [BKM —]                       │
│           [B105 +112] [DK -110] [FD -108] [Pinn -105]                   │
│ BOS +1.5  [WZ -130] [H88 -135] [BFA -128] [BKM —]                       │
│           [B105 -132] [DK -130] [FD -132] [Pinn -135]                   │
│                                                                         │
│ M 52.4%  Pick H88 +115  EV +4.2%  Size $42  To Win $48  [Place] [Log]   │
└─────────────────────────────────────────────────────────────────────────┘
```

Pick-book pill (H88 in the example) renders with green background + green
border + bold weight. Same shape as other pills — no star, no extra glyph
that would break column alignment. Books that don't quote the model's line
render as a muted dashed `—` pill ("no quote at this line" on hover).

**CSS pattern** mirrors parlays directly (`mlb_dashboard.R:1907-2070`):

- Row container `display:flex; flex-wrap:wrap`.
- Full-width cells (`cell-game`, `cell-market`, `cell-pickside`,
  `cell-otherside`) carry `flex-basis:100%`.
- Side label and every pill in the pill row have
  `flex: 0 0 78px; min-width: 78px; box-sizing: border-box;`. Identical
  widths and equal gaps mean pills sit in the same columns on both side
  rows. When the card narrows, pill rows wrap at the same column
  boundary because every pill has the same shape.
- Pill content is stacked: book label (small, muted) on top, price below.
  Inline-flex column inside each pill.
- Metadata cells (`cell-m`, `cell-pick`, `cell-ev`, `cell-size`,
  `cell-towin`, `cell-action`) are `display:inline-flex` and flow
  horizontally. Each cell's label is baked in via `::before`
  pseudo-element (e.g., `.cell-pick::before { content: "Pick" }`) — same
  pattern as parlays' Fair / WZ / Size / To Win labels.
- Metadata text bumped relative to parlays: values 14px, labels 11px,
  button text 13px. Bets are the highest-frequency surface; the
  metadata strip is the actionable strip; bumped sizing keeps it
  scannable.
- Place button pushed right via `margin-left:auto` on the preceding flex
  item.

**Pill rendering** — new helper `render_book_pill(book, american_odds)` in
`mlb_dashboard.R`:

- 8 books in fixed left-to-right order: `WZ, H88, BFA, BKM, B105, DK, FD,
  Pinn`. Same order on every card so the eye learns positions.
- Quoted: `BookCode` (small, top) over `+Odds` (mono, below).
- No quote at the model's line: muted dashed pill with `—`.
- Pick highlight: green tint, green border, bold weight (no star).

**No M pill in the pill row.** M (model fair %) lives in the metadata
strip, next to Pick, so the pill row stays book-only and column-aligned.
This is the only deviation from the parlays-card pattern, and it's
deliberate — parlay cards have 6 pills (M + 4 books + Cons), bets cards
have 8 books, so M-in-pill-row would make alignment-with-the-other-side
harder.

**Action column behavior:**

- Not yet placed → `[Place]` + `[Log]` buttons
- Successful auto-place → green `placed · #<ticket>` pill (clickable to
  remove the placement record locally; does NOT void the ticket)
- After manual log → `[Placed]` button (click to remove the log record)
- After error → red pill with short status label + `[Place]` and `[Log]`
  still available for retry

`[Place]` dispatches by pick book:

- **WZ pick** → Wagerzon direct REST API (no browser, returns ticket #)
- **H88 / BFA / BetOnline pick** → existing Playwright browser flow
  (today's `[Auto]` behavior, now unified under `[Place]`)
- **DK / FD / Pinn / BKM / B105 pick** → button disabled with tooltip
  "no placement integration; use [Log]"

`[Log]` always records a manual placement, identical to today's behavior
for the bets tab and the parlays tab.

## Auto-placer for WZ singles

Mirrors the parlay direct-API placer. The Wagerzon parlay placer is
production-hardened: it owns the auth session refresh, the preflight
that re-fetches WZ's current price and rejects if it has drifted by more
than $0.01, the ticket extraction, and the orphan-detection path
(WZ-confirmed-but-local-write-failed). Singles re-use all of that.

### Endpoint

- `POST /api/place-bet` becomes a dispatcher with body
  `{bet_hash, account, kelly_bet, actual_size, ...row metadata}`.
  Server inspects `bookmaker_key`:
  - `bookmaker_key == "wagerzon"` → invoke
    `wagerzon_odds.single_placer.place_single(...)`, return same
    response shape as `/api/place-parlay` (status, ticket_number,
    balance_after, error_msg, error_msg_key).
  - `bookmaker_key in {"hoop88", "bfa", "betonlineag"}` → spawn
    Playwright via the existing `bet_placer/placer.py` machinery
    (today's `/api/auto-place` behavior; can be inlined or kept as a
    separate route called internally).
  - any other book → 400, "manual log only for this book."

A new `wagerzon_odds/single_placer.py` Python module:

1. Loads the active WZ account session via the existing helper
2. Calls `ConfirmWagerHelper` for the single leg (same endpoint the
   parlay placer uses, with a 1-leg payload)
3. Compares returned price to `wz_odds_at_place` from the request body;
   if drift > $0.01 → return `price_moved` status
4. Calls `MakeWagerHelper` to submit; extracts ticket on success
5. Returns ticket # + balance snapshot on success, or status + error_msg
   on failure (rejected / auth_error / network_error / orphaned)

The dashboard reads the returned status and renders the appropriate pill
on the card.

### Schema migration

Add to `placed_bets` (mirror `placed_parlays`):

| Column | Type | Notes |
|---|---|---|
| `account` | VARCHAR | WZ account label; NULL for non-WZ placements |
| `status` | VARCHAR | `placing` / `placed` / `price_moved` / `rejected` / `auth_error` / `network_error` / `orphaned`; existing rows backfill to `placed` |
| `ticket_number` | VARCHAR | WZ ticket; NULL for non-WZ |
| `error_msg` | VARCHAR | Full error text; surfaces in pill tooltip |
| `error_msg_key` | VARCHAR | Short key (e.g., `balance_low`, `line_pulled`) for short label rendering |
| `wz_odds_at_place` | INTEGER | American odds WZ was showing at click time, for drift forensics |

Migration script: `Answer Keys/MLB Dashboard/migrations/002_single_placer_columns.py`
— idempotent (`ADD COLUMN IF NOT EXISTS` pattern like
`001_combined_parlay_columns.py`).

`placement_orphans` table from the parlay placer is re-used as-is
(orphaned single = same forensic record format).

### Account selector

Re-use the existing Wagerzon account pills at the top of the dashboard
(per `wagerzon_multi_account.md`). The Place button reads the active
account from `dashboard_settings.wagerzon_last_used` exactly as the
parlay Place button does. Balance gating: if active account's available
balance < bet size, the Place button shows the same red warning pill
("Insufficient balance on <account>") that parlay rows already show.

### Out of scope

- H88 / BFA direct-API placement (their auth flows are partially
  documented in `navigator_learnings.md` but never wired for direct
  REST; Playwright stays the path for v1).
- BetOnlineAG direct-API placement.
- Partial-fill detection beyond today's `fill_status` column (Wagerzon
  fills singles atomically — no observed partial fills on singles in
  production logs).
- Auto-place driven by closing-snapshot triggers.

## Code trace (file paths + line numbers)

### Dashboard (the consumer)

| What | File | Line(s) | Notes |
|---|---|---|---|
| Bets table loader | `Answer Keys/MLB Dashboard/mlb_dashboard.R` | 4476-4485 | `dbGetQuery(con, "SELECT * FROM mlb_bets_combined")` from `mlb_mm.duckdb` |
| `create_bets_table()` | `Answer Keys/MLB Dashboard/mlb_dashboard.R` | 817 | Builds the reactable for the bets tab — this gets replaced |
| Parlays card layout (CSS template to mirror) | `Answer Keys/MLB Dashboard/mlb_dashboard.R` | 1907-2070 | `display:flex` + `flex-basis:100%` + `order:` + `::before` labels + responsive overrides |
| Parlay `render_books_strip()` (helper to copy/adapt) | `Answer Keys/MLB Dashboard/mlb_dashboard.R` | called at 543-552 | Renders the per-book pill row for parlays — copy the shape, drop the M pill, swap content (line+price instead of devigged %) |
| Parlays table builder (column order + cell-classes) | `Answer Keys/MLB Dashboard/mlb_dashboard.R` | 360-680 | `create_parlays_table()` — reference for cell class assignments |
| Place button JS hooks | `Answer Keys/MLB Dashboard/mlb_dashboard.R` | search `data-game-id`, `placeBet(` | Keep these data attributes on the pick-side row |
| WZ account selector + balance gating JS | `Answer Keys/MLB Dashboard/mlb_dashboard.R` | search `_wzRecomputeWarnings`, `wagerzon_last_used` | Re-use for singles' Place button |
| Server (Flask) | `Answer Keys/MLB Dashboard/mlb_dashboard_server.py` | 2334 lines | `/api/place-bet` dispatch branch is the entry point for the WZ direct-API path |
| `auto_place()` (existing Playwright path) | `Answer Keys/MLB Dashboard/mlb_dashboard_server.py` | 2014-2057 | Keep for H88/BFA/BetOnline; WZ no longer routes here |
| `place_parlay()` (reference for the WZ singles dispatch) | `Answer Keys/MLB Dashboard/mlb_dashboard_server.py` | search `/api/place-parlay` | Same wiring template: read account, call placer, write breadcrumb, finalize on response |

### Pipeline (the producer)

| What | File | Line(s) | Notes |
|---|---|---|---|
| Per-book bets construction | `Answer Keys/MLB Answer Key/MLB.R` | 585-770 | `wagerzon_bets`, `wz_alt_bets`, `hoop88_bets`, `bfa_bets`, `bfa_alt_bets`, `bookmaker_bets`, `bkm_alt_bets`, `bet105_bets`, `bet105_alt_bets`, `kalshi_bets` |
| Comparison helpers (already produce per-book rows) | `Answer Keys/MLB Answer Key/Tools.R` | search `compare_spreads_to_wagerzon`, `compare_totals_to_wagerzon`, `compare_moneylines_to_wagerzon`, `compare_alts_to_samples` | Return per-book bet rows; we just keep them |
| **Dedup that throws away comparison data** | `Answer Keys/MLB Answer Key/MLB.R` | 783-787 | `group_by(id, base_market, bet_on) %>% filter(ev == max(ev)) %>% slice_head(n=1)` — **insertion point for the new expansion is right before this** |
| `mlb_bets_combined` write | `Answer Keys/MLB Answer Key/MLB.R` | 846-857 | `DROP TABLE IF EXISTS` + `dbWriteTable` to `mlb_mm.duckdb`. Same pattern for the new table. |
| Odds API frame (DK, FD, Pinn live here) | `Answer Keys/MLB Answer Key/MLB.R` | 155 | `all_books <- unique(game_odds$bookmaker_key)` |
| Odds API `markets=` configuration | `Answer Keys/MLB Answer Key/MLB.R` | search `markets <-` or `odds_api_url` | **Expansion point for reference-book coverage on F3/F7/FG markets** |
| Sharp consensus selection | `Answer Keys/MLB Answer Key/MLB.R` | 127-183 | `SHARP_BOOKS` filter for consensus model — unchanged |

### Team name resolution (already canonical, no work needed)

| What | File |
|---|---|
| Python → Odds API canonical | `Answer Keys/canonical_match.py` |
| R → Odds API canonical | `resolve_offshore_teams()` in `Answer Keys/MLB Answer Key/Tools.R` |

All scrapers already write canonical names. The join key
`(game_id, market, period, side, line)` is unambiguous across books.

### Wagerzon placer (reuse path)

| What | File |
|---|---|
| Shared WZ session helper | `wagerzon_odds/` — used by parlay placer; will be imported by the new single placer |
| Parlay placer (reference implementation) | `wagerzon_odds/parlay_pricer.py` and related — the single placer's preflight/drift/ticket logic mirrors this exactly |
| Account discovery | `wagerzon_odds/wagerzon_accounts.py` |

## Conversation history (design decisions and what was rejected)

### Original plan (2026-05-09)

1. **First proposal — "show other books for each bet":** Three options:
   A persist pre-dedup frame + expandable row; B live `/api/compare`
   endpoint; C precomputed summary columns inline. User reframed: trust
   in model is the issue, not shopping for prices.
2. **Second proposal — three layouts**: flat books-as-columns; flat
   both-sides-per-market (true odds-screen, user picked); click-to-expand
   mini odds-screen.
3. **Highlight options:** cell-only / cell+badge / cell+Pick column /
   side-stripe+cell. User picked cell + Pick column.
4. **Place button:** rejected per-cell-place; single Place button at
   recommended Pick.
5. **Reference books:** model says +250 at WZ; user opens DraftKings to
   verify. Verification book = sharp/mainstream. Settled on 5 bettable
   (WZ, H88, BFA, BKM, B105) + 3 reference (DK, FD, Pinn).
6. **Final layout switch — flat → card:** parlays tab uses card layout;
   adopted wholesale from `mlb_dashboard.R:1907-2070`.

### Revision pass (2026-05-11)

7. **Critique surfaced four gaps:** plan ignored the user's "set up the
   auto placer" ask entirely; dropped the model fair (M) pill that
   parlay cards have; reference books would be `—` on most derivative
   markets because Odds API only pulled F5 mains; smaller hygiene
   issues (worktree, hash collision risk, nearest-line fallback
   contradicting "no alt ladder" non-goal).
8. **M pill placement:** considered (A) leftmost on pick-side pill row,
   matching parlay cards exactly, vs (B) in the metadata strip next to
   EV. Picked B — keeps pill rows book-only and column-aligned.
9. **Column alignment:** considered a CSS grid with column headers vs
   fixed-width flex pills. Picked fixed-width pills (`flex: 0 0 78px`)
   because they wrap responsively without the grid's overflow problem.
   Removed the pick-pill star — green tint + bold conveys the same
   signal without breaking pill width.
10. **Metadata strip text size:** bumped values 12 → 14px, labels
    10 → 11px, button 12 → 13px so the actionable strip is scannable.
11. **Nearest-line fallback dropped:** when a book doesn't quote the
    exact line, render `—` (no row in `mlb_bets_book_prices`) rather
    than amber-colored "kind-of comparison." Removes the
    `is_exact_line` column from the schema.
12. **Odds API expansion approved:** add FG mains + F3/F7 mains + FG
    alts + F5 alt spreads to `markets=`. Future scrapers reserved for
    team-totals / odd-even-runs if needed.
13. **Auto-placer scope locked:** WZ direct-API for singles only;
    H88/BFA/BetOnline keep Playwright; reference books are manual-log
    only. No dry-run.
14. **Manual `[Log]` button stays** for every card alongside `[Place]`,
    mirroring parlays.

## Data sources (already exist)

| Source | Path | Books |
|---|---|---|
| `wagerzon.duckdb` → `mlb_odds` | `wagerzon_odds/` | WZ |
| `hoop88.duckdb` → `mlb_odds` | `hoop88_odds/` | H88 |
| `bfa.duckdb` → `mlb_odds` | (per-book scraper DB) | BFA |
| `bookmaker.duckdb` → `mlb_odds` | `bookmaker_odds/` | BKM |
| `bet105.duckdb` → `mlb_odds` | `bet105_odds/` | B105 |
| `game_odds` frame in `MLB.R` | from Odds API, line 155 | DK, FD, Pinn (and others available — see Odds API expansion below) |

All canonical-team-named (resolved by `canonical_match.py`) and on the
shared key `(game_id, market, bet_on, line, period)`. No new
normalization needed.

### Odds API `markets=` expansion

Today's `MLB.R` Odds API pull covers (F5 mains):

- `h2h_1st_5_innings, totals_1st_5_innings, spreads_1st_5_innings`
- `alternate_totals_1st_5_innings`

Adding for reference-book coverage on more market types:

- `h2h, totals, spreads` (FG mains)
- `alternate_totals, alternate_spreads` (FG alts)
- `alternate_spreads_1st_5_innings` (F5 alt spreads — F5 alt totals
  already present)
- `h2h_1st_3_innings, totals_1st_3_innings, spreads_1st_3_innings`
- `h2h_1st_7_innings, totals_1st_7_innings, spreads_1st_7_innings`

**Cost estimate:** ~9 additional markets × ~15 events per pipeline run
≈ ~135 extra credit-hits on top of today's ~60. Each market is one HTTP
call regardless of book count — DK, FD, Pinn (and any others in the
default region) all come free once a market is requested. If credits
become a constraint, F3/F7 reference coverage can be dropped first (the
offshore scrapers still populate those bettable books).

**Markets the Odds API doesn't carry** (and where reference pills
render `—`): team totals, odd/even runs. v2 future work: targeted DK/FD
single-leg scrapers if verification on these markets becomes important.

## Schema — new table `mlb_bets_book_prices`

Long format, written to `mlb_mm.duckdb` alongside `mlb_bets_combined`.
One row per `(bet × book × side)`. Only exact-line quotes are written —
no nearest-line fallback.

| Column | Type | Notes |
|---|---|---|
| `bet_row_id` | VARCHAR | Stable hash of `(game_id, market, line, bet_on)` — joins back to `mlb_bets_combined`. **Does NOT include `bookmaker_key`** so all books for the same bet share an id |
| `game_id` | VARCHAR | Odds API id |
| `market` | VARCHAR | Normalized: `spreads`, `totals`, `moneyline`, plus alt variants |
| `period` | VARCHAR | `FG`, `F5`, `F3`, `F7` |
| `side` | VARCHAR | `pick` or `opposite` (relative to the bet the model recommended) |
| `bookmaker` | VARCHAR | One of: `wagerzon`, `hoop88`, `bfa`, `bookmaker`, `bet105`, `draftkings`, `fanduel`, `pinnacle` |
| `line` | DOUBLE | Always equals the model's line (since we only write exact-line rows) |
| `american_odds` | INTEGER | Price at `line` |
| `fetch_time` | TIMESTAMP | Source row's `fetch_time` (per-book freshness signal) |

Grain: `(bet_row_id, bookmaker, side)`. ~250 bets × 8 books × 2 sides ≈
4000 rows. Trivial for DuckDB.

**Hash separation (intentional):** the existing per-book `bet_hash` in
`mlb_dashboard.R:55` (game + market + line + bet_on + **book**) keeps
owning `placed_bets` keys. `bet_row_id` (no book) is purely a join key
into `mlb_bets_book_prices`. Placing at WZ writes a row keyed on the
per-book `bet_hash` for `placed_bets`; it does NOT mark BFA's pill on the
same card as "placed."

## Implementation

### Phase 1 — pipeline (`MLB.R` + `Tools.R`)

**File:** `Answer Keys/MLB Answer Key/MLB.R`
**Insertion point:** between the dedup `slice_head(n=1)` (line 787) and
the `adjust_kelly_for_correlation` step (~line 820).

1. **Expand Odds API `markets=`** at the existing odds-pull call (line
   155 area). Add the 9 markets listed above to the `markets=` string.
   Verify the Odds API budget after one production run.
2. Add `bet_row_id` column to `all_bets_combined` via deterministic hash
   of `(game_id, market, line, bet_on)`.
3. Build `book_prices_long`:
   - For each row in `all_bets_combined`, expand into one row per
     `(bookmaker, side ∈ {pick, opposite})` **only when that book quotes
     the exact line+market+side**.
   - For offshore books: source from per-book `mlb_odds` tables already
     loaded (`wagerzon_odds`, `hoop88_odds`, `bfa_odds`, `bookmaker_odds`,
     `bet105_odds`).
   - For DK / FD / Pinn: source from `game_odds` filtered to those
     `bookmaker_key`s.
   - Join key: `(game_id, market, period, side, line)` — strict equality
     on `line`. No book contributes a row if it doesn't quote that
     exact line.
4. Write `mlb_bets_book_prices` to `mlb_mm.duckdb` in the same writer
   block as `mlb_bets_combined` (`MLB.R:851-857`). Drop-and-replace,
   same as the bets table.

**New helper** in `Answer Keys/MLB Answer Key/Tools.R`:

```r
expand_bets_to_book_prices <- function(bets, book_odds_by_book) { ... }
```

Pure function — takes the bets frame + a named list of per-book odds
frames, returns the long-format table. Unit-testable in isolation.

### Phase 2 — dashboard (`mlb_dashboard.R`)

**File:** `Answer Keys/MLB Dashboard/mlb_dashboard.R`

1. **Loader (~line 4476):** load `mlb_bets_book_prices` from
   `mlb_mm.duckdb` alongside `mlb_bets_combined`. Defensive try/catch:
   if missing, fall back to today's table renderer.
2. **Pivot** book_prices long → wide on `(bet_row_id, side, bookmaker)`.
   Resulting columns per side: `wz_line`, `wz_odds`, …, `pinn_odds`.
   Two rows per bet (one per side) joined back to the bet metadata.
3. **Replace `create_bets_table()`** with a card-layout reactable
   matching the locked layout. Use the same column-classes pattern as
   the parlays table.
4. **New helper `render_book_pill(book, american_odds)`:** renders one
   pill. NA odds → muted dashed pill with `—`. Pick book (when
   `book == pick_book && side == pick_side`) → green-tint pill with
   bold weight.
5. **Metadata strip cells** (`cell-m`, `cell-pick`, `cell-ev`,
   `cell-size`, `cell-towin`, `cell-action`) with bumped font sizing
   per the design.
6. **Place button + Log button** in the action cell. Place button JS
   reads pick book from `data-book`, dispatches:
   - `wagerzon` → POST `/api/place-bet` (server branches to the WZ
     single placer).
   - `hoop88` / `bfa` / `betonlineag` → POST `/api/auto-place`
     (existing Playwright path).
   - others → button disabled with tooltip "no placement integration."
7. **Correlation dot preserved** as a small corner badge on the card
   (top-right of the game header). Same tooltip listing other bets on
   the same game.
8. **Filters preserved.** Book filter / market filter / EV threshold /
   correlation-status filter all apply to which cards render. No
   filter UI changes.
9. **Placed-bets summary** at the top of the bets tab stays as a flat
   reactable (separate from the cards below).
10. **CSS additions** in the inline style block — match existing
    dashboard palette. Document new classes (`cell-pickside`,
    `cell-otherside`, `cell-m`, `cell-pick-highlight`, `cell-no-quote`).

### Phase 3 — server + placer (`mlb_dashboard_server.py` + new `wagerzon_odds/single_placer.py`)

**File:** `Answer Keys/MLB Dashboard/mlb_dashboard_server.py`

1. `/api/place-bet` becomes a dispatcher:
   - Read `bookmaker_key` from the bet metadata (or accept it in the
     body for explicitness).
   - WZ → call `wagerzon_odds.single_placer.place_single(...)`. Mirror
     `/api/place-parlay`'s response shape: status, ticket_number,
     balance_after, error_msg, error_msg_key.
   - H88 / BFA / BetOnline → spawn Playwright via existing
     `/api/auto-place` machinery (now invoked internally rather than
     via a separate JS button).
   - Other books → 400 with "manual log only for this book."
2. Reuse the parlay placer's account-selector pattern: read active
   account from `dashboard_settings.wagerzon_last_used`; require
   `account` in the request body so server-side state is explicit.
3. Write a breadcrumb row to `placed_bets` with `status='placing'`
   before invoking the placer; upsert to the final status after the
   placer returns. Same anti-orphan pattern as parlays.

**File:** `wagerzon_odds/single_placer.py` (new)

1. `place_single(account, bet)` — load session, build single-leg
   payload (one selection encoded the same way the parlay placer does),
   call `ConfirmWagerHelper` for preflight + drift check, call
   `MakeWagerHelper` for submission, extract ticket, return structured
   result.
2. Error taxonomy matches parlays: `placed`, `price_moved`, `rejected`,
   `auth_error`, `network_error`, `orphaned`.
3. Import shared session/preflight helpers from the existing parlay
   placer module — no duplicate code.

**Migration script:** `Answer Keys/MLB Dashboard/migrations/002_single_placer_columns.py`
— idempotent ADD COLUMN for `account`, `status`, `ticket_number`,
`error_msg`, `error_msg_key`, `wz_odds_at_place` on `placed_bets`.

## Implementation order (suggested)

Two commits on the worktree branch, smallest blast radius first, docs
ship with the feature commit per `CLAUDE.md`:

1. **Pipeline (data flowing, dashboard still old):**
   - Add `bet_row_id` hash + `expand_bets_to_book_prices()` helper in
     `Tools.R`.
   - Expand Odds API `markets=` in `MLB.R`.
   - Wire the expansion into `MLB.R` between line 792 and line 820.
   - Add DROP+WRITE for `mlb_bets_book_prices` alongside
     `mlb_bets_combined`.
   - Run pipeline once; confirm table populated; spot-check a slate
     across multiple games.
   - Commit: `feat(mlb): write mlb_bets_book_prices in pipeline`.

2. **Dashboard rebuild + WZ single placer (behavior change visible) +
   docs (same commit):**
   - Run migration `002_single_placer_columns.py` on
     `mlb_dashboard.duckdb`.
   - Add `wagerzon_odds/single_placer.py`.
   - Update `/api/place-bet` dispatch in `mlb_dashboard_server.py`.
   - Replace `create_bets_table()` body — new card layout, helpers,
     CSS, Place + Log button wiring.
   - Defensive loader fallback if `mlb_bets_book_prices` is missing.
   - Test in browser: standard + alt bets, WZ Place flow end-to-end,
     non-WZ Place flow (Playwright still works), manual Log flow,
     filters, correlation dot.
   - Update `Answer Keys/MLB Dashboard/README.md` with new layout
     sketch, button behaviors, schema notes.
   - Update `Answer Keys/CLAUDE.md` to list `mlb_bets_book_prices` and
     the new `placed_bets` columns.
   - Commit: `feat(mlb-dashboard): card-layout odds-screen bets tab + WZ single-bet auto-placer`.

3. **Pre-merge review** per checklist below; explicit user approval;
   merge to `main`; remove worktree + delete local branch.

## How to verify (post-implementation)

After implementing and running `Rscript MLB.R`:

```bash
# 1. New table exists with expected schema
duckdb "Answer Keys/mlb_mm.duckdb" \
  "DESCRIBE mlb_bets_book_prices"

# 2. Row count looks reasonable (~bets × 8 books × 2 sides, but lower
#    because not every book quotes every line)
duckdb "Answer Keys/mlb_mm.duckdb" \
  "SELECT COUNT(*) AS rows,
          COUNT(DISTINCT bet_row_id) AS bets,
          COUNT(DISTINCT bookmaker) AS books
   FROM mlb_bets_book_prices"

# 3. Spot-check one bet across books, both sides
duckdb "Answer Keys/mlb_mm.duckdb" \
  "SELECT bookmaker, side, line, american_odds
   FROM mlb_bets_book_prices
   WHERE bet_row_id = (SELECT bet_row_id FROM mlb_bets_book_prices LIMIT 1)
   ORDER BY side, bookmaker"

# 4. Reference-book coverage spot-check by market+period
duckdb "Answer Keys/mlb_mm.duckdb" \
  "SELECT market, period,
          SUM(CASE WHEN bookmaker IN ('draftkings','fanduel','pinnacle') THEN 1 ELSE 0 END) AS ref_rows,
          COUNT(*) AS total
   FROM mlb_bets_book_prices
   GROUP BY market, period
   ORDER BY market, period"
# Reference rows should be non-zero for FG/F3/F5/F7 mains + alts.
# Should be zero for team_totals_* and odd_even_runs (Odds API doesn't carry).
```

In the dashboard:

- Open `http://localhost:8083`, switch to bets tab.
- Each card: game header, market header, two pill rows with aligned
  columns, metadata strip with M / Pick / EV / Size / To Win / Place +
  Log.
- Pick book pill is green-tinted + bold (no star).
- At least one card should show `—` pills for some books — verify
  muted dashed rendering.
- Click `[Place]` on a WZ-pick card → toast `placed · #<ticket>`; card
  flips to placed pill; ticket appears in `placed_bets` with
  `status='placed'` and `ticket_number` populated.
- Click `[Place]` on an H88-pick card → Playwright opens with bet
  pre-filled (existing flow); click confirm on the book; manually click
  `[Log]` to record locally.
- Click `[Log]` on any card → instant log, no placement attempt.
- Switch WZ accounts via the top-bar selector → next Place uses the
  new account; account label appears in `placed_bets.account`.
- Filters (book / market / EV / correlation status) all reduce the
  card set.

## Rollback

If the dashboard rebuild breaks:

```bash
# Revert just the dashboard+placer commit
git revert <dashboard-commit-sha>
git push origin main
cd "Answer Keys/MLB Dashboard"
./run.sh
```

The pipeline write commit can stay — it's harmless dataset growth. The
defensive loader fallback means the dashboard reverts to today's table
automatically if `mlb_bets_book_prices` is missing.

If the pipeline write itself is causing issues:

```bash
git revert <pipeline-commit-sha>
duckdb "Answer Keys/mlb_mm.duckdb" "DROP TABLE IF EXISTS mlb_bets_book_prices"
```

If the single placer is causing WZ orphan rows: disable the WZ branch
in the `/api/place-bet` dispatcher (force fallthrough to manual log)
until the root cause is fixed. `placement_orphans` table will contain
forensics for any orphans created.

## Version control

- **Branch:** `claude/mlb-sportsbook-comparison-98fM8` (already exists
  with the 3 plan commits).
- **Worktree:** `.worktrees/mlb-odds-screen` based on
  `claude/mlb-sportsbook-comparison-98fM8`. Created with
  `git worktree add .worktrees/mlb-odds-screen claude/mlb-sportsbook-comparison-98fM8`.
  Removed after merge with
  `git worktree remove .worktrees/mlb-odds-screen` and
  `git branch -d claude/mlb-sportsbook-comparison-98fM8`.
- **Commits, in order:**
  1. `feat(mlb): expand per-book odds + write mlb_bets_book_prices in MLB.R`
  2. `feat(mlb-dashboard): card-layout odds-screen bets tab + WZ single-bet auto-placer` (includes docs)
- **No merge to `main`** until: pipeline runs end-to-end with the new
  table, dashboard renders without errors on a real slate, WZ single
  Place flow places + receives ticket against the live Wagerzon API in
  one verified bet, executive review per `CLAUDE.md` complete, user
  explicitly approves.

## Documentation

Files to update in the same PR as the code:

- `Answer Keys/MLB Dashboard/README.md` — update bets-tab section with
  the new card layout sketch, list of book columns, explanation of the
  Place / Log button dispatch by book, schema note for the new
  `mlb_bets_book_prices` and `placed_bets` columns.
- `Answer Keys/CLAUDE.md` — note that `MLB.R` writes
  `mlb_bets_book_prices` to `mlb_mm.duckdb`, and that `placed_bets`
  has new columns from the WZ single placer.
- Root `CLAUDE.md` — no changes expected; the change is contained
  inside the MLB dashboard subsystem.

## Edge cases / decisions

- **Doubleheaders:** `game_id` is unique per game, so the join key
  handles this naturally. No special logic.
- **Period:** F5 / F3 / F7 / FG are filtered as part of the join key.
- **Markets the Odds API doesn't track** (team totals, odd/even runs):
  reference columns show `—`. Bettable columns still populated from
  per-book scraper DBs.
- **Book outage / missed scrape:** column shows `—` for affected games.
- **Alts:** alt bets get their own card; within each card, books that
  don't offer that exact alt show `—`.
- **First run after deploy with no `mlb_bets_book_prices` table yet:**
  dashboard loader try/catch falls back to today's column set, prints
  a warning. Pipeline run repopulates it within one cycle.
- **WZ session expired during Place click:** placer attempts one
  re-auth (existing parlay placer behavior). If still failing → return
  `auth_error` status, render red pill on card.
- **WZ Place succeeds but local write fails:** orphan written to
  `placement_orphans`; card shows `orphaned` red pill. User checks WZ
  ticket history to confirm placement before retrying.
- **Account switched mid-session:** Place button re-reads
  `wagerzon_last_used` on each click; no stale account state.
- **Non-WZ Place with no Playwright integration:** button disabled with
  tooltip; user clicks Log instead.
- **Correlation dot:** displayed only when the user has placed bets on
  the same game (existing behavior). Tooltip lists other bets with
  market / side / size / book.

## Pre-merge executive review checklist

Per `CLAUDE.md` pre-merge requirements, before merging to `main`:

**Pipeline (Phase 1):**

- [ ] `git diff main..HEAD --stat` reviewed
- [ ] No duplicate writes to `mlb_bets_book_prices` (DROP + WRITE
      pattern, same as `mlb_bets_combined`)
- [ ] DB connections wrapped in `on.exit(dbDisconnect(...))`
- [ ] Off-season behavior: pipeline produces empty bets → empty
      `mlb_bets_book_prices` → dashboard renders empty bets tab
      without error
- [ ] Doubleheader behavior: two games same teams → distinct
      `bet_row_id`s → no row collision
- [ ] Odds API expansion didn't push the daily credit budget over a
      red line; check usage after one full pipeline run

**Dashboard (Phase 2):**

- [ ] No new dependencies added to `requirements.txt` / `DESCRIPTION`
- [ ] Dashboard renders end-to-end with real pipeline output
- [ ] Place + Log buttons work for both standard and alt bets
- [ ] Filters (book / market / EV / correlation status) all still
      filter the card set
- [ ] Correlation dot still renders + tooltip lists other bets
- [ ] Placed-bets summary table at top of tab still renders + accepts
      remove
- [ ] No secrets / API keys logged

**Auto-placer (Phase 3):**

- [ ] Migration `002_single_placer_columns.py` is idempotent against
      an existing `mlb_dashboard.duckdb`
- [ ] One verified live WZ single-bet placement returns a ticket and
      `placed · #<ticket>` pill renders correctly
- [ ] One forced `price_moved` test (mock WZ response with drifted
      price) returns `price_moved` status without placement
- [ ] WZ session re-auth path is exercised at least once (let the
      session expire, click Place, confirm re-auth succeeds)
- [ ] Account selector switches accounts mid-session without losing
      state; `placed_bets.account` populated correctly
- [ ] `placement_orphans` inspectable; one forced orphan test (mock
      WZ confirms but local write raises) lands cleanly

**General:**

- [ ] README updates committed in the same commit as the feature code
- [ ] User explicit approval to merge

## Open questions

None blocking. All design decisions resolved through the 2026-05-09
brainstorm + 2026-05-11 revision pass.
