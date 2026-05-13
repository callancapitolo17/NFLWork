# MLB Dashboard Bets Tab — Layout & Data Improvements

**Status:** Design · 2026-05-13
**Branch / worktree:** `worktree-mlb-bets-tab-improvements`
**Owner:** Callan
**Mockups:** `.superpowers/brainstorm/37414-1778712316/content/bet-card-layouts-v8.html` (V7 Variant B + centered stats)

---

## Review Pack

**What we're building**

A redesign of the "Available Bets" cards on the MLB Dashboard (port 8083). Each card switches from today's loose pill row to a strict 2×8 grid (sides × books) with a small "what to do" hero strip — pick book, edge, risk, to-win, action buttons — sitting in one band at the top. Two data bugs that motivated the redesign are fixed at the same time: FanDuel's First-7-Innings markets are missing entirely from the scraper whitelist, and the FanDuel alternate-run-line scraper isn't capturing all alt lines on at least some games.

**Key decisions**

- **Decision:** Render each bet as a 2×N price grid (rows = sides, columns = books) with cell-level alignment. **Rejected:** keep current pill rows; add a sidebar summary; promote stats to a header chip strip. **Why:** the underlying data *is* a sides × books matrix, and column alignment makes cross-book price comparison the primary affordance — that is the whole point of the bets tab.
- **Decision:** Always render both sides on spreads (e.g., BOS -2.5 row + PHI +2.5 row), mirroring how totals already render Over + Under. **Rejected:** only render the picked side. **Why:** user explicitly asked to see opponent + other side odds; symmetry with totals removes a special case.
- **Decision:** Replace "Market %" with "Fair" (fair de-vigged American odds, e.g. `+199`) as a standalone stat in the hero strip alongside EV / Risk / To Win. **Rejected:** keep "Market %" as a probability; tuck Fair under the pick odds in muted text; drop Fair entirely. **Why:** fair price in American odds is the most legible side-by-side comparison to the actual quoted odds — "I'm getting +160, fair is +199" reads at a glance, while "Market 33.4%" did not. Giving Fair its own slot rather than tucking it under the pick odds also makes it scannable across cards.
- **Decision:** Keep the label "EV". **Rejected:** rename to "Edge". **Why:** user prefers the existing label; "Edge" was speculative on my part.
- **Decision:** Recon FanDuel's actual market availability for F7 + alt-run-lines *before* coding the fixes. **Rejected:** ship the whitelist extension blind. **Why:** if FD doesn't post F7 totals or paginates alt lines server-side, the fix path is different.

**Risks / push back here**

- The grid degrades below ~400px wide (8 book columns get cramped). MLB Dashboard is desktop-only today, but if you ever want it on a phone, we'd need a stacked layout breakpoint.
- The FD alt-spread fix has unknown scope until recon — it could be a one-line whitelist tweak or it could require pagination logic in the scraper.
- Some bets won't have a clean "opposite" — e.g., when the model is on an alt-spread line that's not symmetric across books. The opposite row will show alt-line fallbacks (amber) rather than exact matches, same as the pick side already does. Flagging in case you want a stricter "hide opposite when no exact match" behavior.
- Fair price is computed by de-vigging the consensus / Pinnacle line — make sure the pipeline that already produces "Market %" continues to expose enough info to derive a fair American-odds number cleanly. If the existing pricing path drops the fair-prob → fair-odds conversion downstream of the model, we'll need to surface it in `mlb_bets_combined`.

**Worth understanding**

- **CSS Grid `1fr` columns.** The grid uses `grid-template-columns: minmax(85px, 125px) repeat(8, minmax(0, 1fr))`. This means: one fixed-ish row label column, then 8 book columns that share whatever width is left, equally. The `minmax(0, 1fr)` is important — without it, cells refuse to shrink below their content size. Closest R analogue: `par(mfrow = c(2, 8))` divides plotting space proportionally, but CSS Grid lets you express that *with constraints* on each track.
- **Long → wide pivot for the cards.** `mlb_bets_book_prices` is stored as long-form (one row per bet × book × side). The dashboard pivots it to wide format (`{book}_american_odds`, `{book}_line_quoted`, etc.) before rendering. This is exactly `tidyr::pivot_wider()`. If you ever see a card with a phantom column or a misaligned cell, the join key in `expand_bets_to_book_prices()` is the place to look.

---

## Background

Today's bets tab renders each `mlb_bets_combined` row as a card with:

- A single header line: `{matchup} · {tipoff}` + a market title line below it.
- A row of pills (one per book) showing line + odds, with the picked book highlighted in green and alt-line mismatches in amber.
- A footer strip with EV / Market % / Size / To Win / Pick / Place + Log buttons in left-to-right flex order.

Four issues motivated this redesign:

1. **FanDuel F7 totals are silently dropped.** The `_FD_MARKET_WHITELIST` in `mlb_sgp/scraper_fanduel_singles.py` covers FG and F5 only. F7 markets ("First 7 Innings Total Runs", etc.) are not listed, so they never reach `mlb_bets_book_prices`. The screenshot user shared shows an F7 totals bet where FD's column is empty even though FD does post that market.
2. **FanDuel alt-spread fallback to closest line, not exact match.** On the Phillies @ Red Sox alt-spread BOS -2.5 bet, FanDuel's pill shows `-1.5 +150` (closest line within ±3.0 tolerance), not `-2.5` — even though FD posts -2.5. Cause unknown until recon: likely either (a) the scraper doesn't paginate through all alt lines, or (b) `get_fd_odds()` in `Tools.R` collapses multi-row alt frames.
3. **Spread cards don't render the opposite side.** Totals already render both Over and Under rows; spreads only render the picked side. User wants symmetric behavior so they can see what the *other* team's odds look like across books.
4. **Visual feel.** The footer flex strip with EV / M / SIZE / TO WIN / PICK feels "spread out" and the bet description gets buried below the matchup line. The pill row isn't column-aligned across books, so cross-book price comparison takes work.

## Final design (locked in V6)

Visual reference: `bet-card-layouts-v6.html` in the brainstorm directory. Mockup served at http://localhost:60595 during brainstorming.

### Card structure (top-to-bottom)

1. **Bet title row** — uniform 18px. Primary text (bright): the bet itself, e.g. `BOS -2.5`. Secondary text (muted): the market type, e.g. `· Alt Spread · Full Game`. Same font size on the whole row; color contrast does the hierarchy.
2. **Matchup row** — uniform 15px. Primary (bright): team names, e.g. `Philadelphia Phillies @ Boston Red Sox`. Secondary (muted): tipoff time, e.g. `· Wed 10:46 PM`. Same size-on-row rule applies.
3. **Hero strip** — bordered green-tinted band containing six elements in a flex row. The four stats (Fair / EV / Risk / To Win) are **center-aligned within their slots** (label and value both centered, each slot has a `min-width: 60px` so they don't squeeze together); the pick block stays left-anchored and the action buttons stay grouped immediately after To Win.
   - **Pick block:** book name (small all-caps green) + odds (big white, 22px tabular-nums)
   - **Vertical divider** (1px, low opacity)
   - **Fair stat:** label "FAIR" (small) + fair de-vigged American odds (18px white), e.g. `+199`
   - **EV stat:** label "EV" (small) + percentage (18px, **always green** — every bet on this tab is pre-filtered to +EV, so the value never goes negative; no conditional coloring)
   - **Risk stat:** label "RISK" (small) + dollar value (18px white)
   - **To Win stat:** label "TO WIN" (small) + dollar value (18px green)
   - **Action buttons:** Place Bet (filled green-tint, bold) + Log (outlined). Sit immediately after To Win with an 8px gap — no auto-spacer.
4. **Price grid** — strict 2×8 CSS grid:
   - **Columns:** row-label column + 8 book columns (WZ, H88, BFA, BKM, B105, DK, FD, PINN). Columns share remaining width equally via `minmax(0, 1fr)`. Book labels are 9px uppercase headers separated from cells by a 1px bottom border.
   - **Rows:** one per side. For spreads/ML: pick-side row + opposite-side row. For totals: Over row + Under row.
   - **Cell states:**
     - `pick`: green border + green-tint fill, price in green. Exactly one cell per card.
     - `exact`: book has this exact line; subtle dark-gray fill, price in neutral white.
     - `alt`: book's closest line differs from the bet line. Amber-tint fill, with the book's actual line (e.g., `-1.5`) shown small above the price in amber.
     - `empty`: book doesn't post this market. Transparent fill, dim "—" placeholder.

### Responsive behavior

Grid columns are `1fr` so cells resize with viewport width. Hero strip uses `flex-wrap: wrap` so buttons drop to a second line below the stats on narrow widths. **Known limitation:** below ~400px wide, book column headers begin to overlap and alt-line cells get cramped. Acceptable for a desktop dashboard; would need a stacked breakpoint for mobile (out of scope).

### Removed / replaced from current UI

- **"Market %"** — replaced by **"Fair"** in the hero strip. Same underlying number (de-vigged probability), expressed as American odds instead of a percentage so it sits next to the actual quoted odds visually.
- **Empty pills as visible placeholders** — replaced by dim "—" cells in the grid. Still visible (so you can tell at a glance which books aren't posting), but they don't dominate visually because they're transparent rather than bordered.

## Scope of changes (code)

Five workstreams. Each is independently shippable.

### 1. Backend: write `opposite` rows for spread bets

**File:** `Answer Keys/MLB Answer Key/odds_screen.R` (function `expand_bets_to_book_prices`)

Today `expand_bets_to_book_prices` emits `side = "pick"` for every bet, plus `side = "opposite"` only for totals (where the "Over 5.5" pick implies an "Under 5.5" opposite at the same line). Extend this so spreads / moneylines also emit an opposite row:

- For spread bets (`market` starts with `spreads` or `alt_spreads`): opposite has flipped line (`+2.5` ↔ `-2.5`) and references the other team. Use the existing `away_spread` / `home_spread` columns in each book's odds frame to pull the actual quoted price.
- For ML bets (`market = h2h`): opposite is the other team's moneyline on the same book.
- The opposite row uses the same line-matching logic (closest within `LINE_MATCH_TOLERANCE = 3.0`) — so opposite cells also get the `is_exact_line` flag and can be flagged as `alt` in the UI when the book's nearest line differs from the bet line.

Schema change to `mlb_bets_book_prices`: none (the `side` column already exists, we're just populating more rows). Row count roughly doubles for spread + ML bets.

### 2. Dashboard layout: rewrite `create_bets_table()`

**File:** `Answer Keys/MLB Dashboard/mlb_dashboard.R` (currently lines ~1143–1506)
**File:** `Answer Keys/MLB Dashboard/book_pill.R` → likely rename to `book_cell.R`

Replace the pill-row rendering with the grid layout from V6. Concrete substeps:

- Replace `BOOK_ORDER` flex loop with a CSS Grid template. Move per-row rendering into a `render_price_row(side_label, prices_for_side, pick_book, is_totals_market)` helper.
- Replace `render_book_pill()` with `render_book_cell()` — same inputs, but emits a grid cell (one of `pick` / `exact` / `alt` / `empty`) rather than a fixed-width pill.
- Build the hero strip as a separate component: takes `(pick_book, pick_odds, fair_odds, ev_pct, risk_dollars, towin_dollars, place_button_html, log_button_html)` and renders the green-tinted band with the four stats centered in their slots.
- Header: switch from two distinct font sizes per line to uniform-per-row sizing with color-only hierarchy.
- Replace the Market % cell with the new Fair stat. Compute fair American odds from the existing fair probability (`if p >= 0.5: -100 * p / (1 - p)` else `100 * (1 - p) / p`). If `mlb_bets_combined` already exposes a fair-odds column, use it directly; otherwise add the conversion in the pricing path upstream.
- Inline CSS for the new card class (single style block at the top of `mlb_dashboard.R`, near the existing card CSS around line 2544).

### 3. Data fix: FanDuel First-7-Innings market coverage

**File:** `mlb_sgp/scraper_fanduel_singles.py`

Extend `_FD_MARKET_WHITELIST` to include F7 markets:
```python
# First 7 Innings
"First 7 Innings Run Line",
"First 7 Innings Total Runs",
"First 7 Innings Money Line",
"First 7 Innings Alternate Run Lines",
"First 7 Innings Alternate Total Runs",
```

**Recon step first:** before adding these, run the scraper in debug/dry-run mode against a live MLB slate and confirm FanDuel actually returns these market names. If FD uses different strings ("1st 7 Innings", "F7 Total Runs", etc.), use the actual strings observed.

`classify_market()` in `scraper_draftkings_singles.py` already accepts F7 via its period detection (it just doesn't see F3) — so DK F7 markets should already work. Confirm during recon that the DK column does in fact populate for F7 cards in production.

### 4. Data fix: FanDuel alternate-run-line exact-line capture

**File:** `mlb_sgp/scraper_fanduel_singles.py` and possibly `Answer Keys/Tools.R::get_fd_odds`

**Recon required.** Pull a live FD `Alternate Run Lines` market response on a game we have a -2.5 bet on, and check:

- Does the API response include the -2.5 selection? If no → FD doesn't post it on that game and there's no bug.
- If yes, does the scraper persist it to `fd_odds/fd.duckdb::mlb_odds`? If no → scraper is filtering or paginating wrong.
- If persisted, does `get_fd_odds()` in `Tools.R` pivot it through? If no → loader collapses alt rows.

Fix scope depends on what recon finds. The lightweight outcome is "the FD API returns all alt lines in one response and the scraper persists them, but `get_fd_odds()` only emits the main line per game." The heavier outcome is "FD paginates alt lines and the scraper only fetches the first page."

### 5. Label & copy polish

- Keep label "EV" as-is (no rename). EV value always renders green since bets-tab is pre-filtered to +EV; no `if ev > 0` conditional in the rendering code.
- Time/date format on the matchup row: use the existing `Wed 10:46 PM` format, just shifted to smaller secondary text.
- Place button copy: "Place Bet" (current is "Place"). Worth checking that the button handler still wires up — should be a pure text change.

## File-by-file change list

| File | Change | Workstream |
|------|--------|------------|
| `Answer Keys/MLB Answer Key/odds_screen.R` | `expand_bets_to_book_prices`: emit opposite rows for spreads + ML | 1 |
| `Answer Keys/MLB Dashboard/mlb_dashboard.R` | Rewrite `create_bets_table()`; new card CSS; replace Market % with Fair | 2 |
| `Answer Keys/MLB Dashboard/book_pill.R` → `book_cell.R` | Convert pill renderer → grid cell renderer | 2 |
| `mlb_sgp/scraper_fanduel_singles.py` | Extend `_FD_MARKET_WHITELIST` w/ F7 markets after recon | 3 |
| `mlb_sgp/scraper_fanduel_singles.py` | Alt-line fix (scope TBD by recon) | 4 |
| `Answer Keys/Tools.R` (possibly) | Alt-line fix in `get_fd_odds` (scope TBD by recon) | 4 |
| `Answer Keys/MLB Dashboard/README.md` | Update screenshot + bets-tab section | docs |
| `Answer Keys/CLAUDE.md` (or top-level) | Update if the `mlb_bets_book_prices` row-count semantics changed enough to flag | docs |

## Worktree & branch lifecycle

- Already in worktree `worktree-mlb-bets-tab-improvements` at `.claude/worktrees/mlb-bets-tab-improvements/`. Branch of same name.
- Implementation plan will run inside this worktree.
- Before merge to main: full pipeline run in worktree (`Rscript MLB.R` → restart dashboard) and verify all four issues resolved on the live dashboard.
- After merge: `git worktree remove .claude/worktrees/mlb-bets-tab-improvements && git branch -d worktree-mlb-bets-tab-improvements`.
- Do **not** merge to main until I've explicitly approved the implementation result.

## Documentation updates

In the same merge commit as the code changes:

- `Answer Keys/MLB Dashboard/README.md` — update bets-tab section to describe the new grid layout, refresh any screenshot.
- `Answer Keys/CLAUDE.md` — confirm the existing "Odds screen" section is still accurate; update the "DraftKings and FanDuel pill data" wording since we no longer render pills.
- `mlb_sgp/README.md` — if F7 support is being added to the FD scraper, add a one-line note under the FanDuel singles section noting F7 coverage.

## Risks

- **Layout-only changes have UI regression surface area.** The bets tab is the primary action surface during a slate. Mitigate by previewing in the worktree dashboard (a separate port) before merge.
- **`expand_bets_to_book_prices` change doubles row count for spread bets.** Should be fine — `mlb_bets_book_prices` is small (< 10k rows typical) — but flag if downstream consumers (combined parlay, bet logger) read this table and don't filter by `side = 'pick'`.
- **Recon-then-fix workstreams (3, 4) carry timing risk.** If FD's market names don't match what we guessed, the spec's whitelist additions need updating before code lands.
- **CSS in R Shiny / shinydashboard isn't hot-reload-friendly.** Each iteration cycle is "edit R → reload dashboard → eyeball." Plan for ~5–10 visual iterations after the structural change lands.

## Out of scope

- Mobile / phone layout (under-400px breakpoint).
- Bets-tab filtering, sort order, or search-bar behavior — untouched.
- Other tabs (Pricing, Slate, Performance) — untouched.
- DraftKings F7 verification (DK already covers F7 in `classify_market`; confirm during recon but don't expect a fix here).
- Pinnacle data pipeline — out of scope; PINN cells will continue to populate from the Odds API path as today.
