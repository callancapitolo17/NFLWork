# MLB Dashboard

Interactive web dashboard displaying +EV MLB F5 betting opportunities and correlated parlays.

## Architecture

```
run.py mlb (pipeline)  →  mlb_mm.duckdb/mlb_bets_combined          (dashboard refresh)
                       →  mlb_mm.duckdb/mlb_parlay_opportunities    (place-bet loader)
                       →  mlb_mm.duckdb/mlb_trifecta_opportunities
                                ↓
                      mlb_dashboard.R (HTML generation)  →  report.html
                                                                  ↓
                                                   mlb_dashboard_server.py (Flask, :8083)
                                                                  ↓
                                                            Browser (:8083)
```

- **mlb_dashboard.R**: Reads pipeline bets + parlays, generates reactable HTML
- **mlb_dashboard_server.py**: Flask server + REST API for placement/CLV state
- **run.sh**: Full pipeline launcher (scrapers → R model → parlay pricer → correlated parlay finder → dashboard → Flask)

## Features

- **Bets tab odds screen** — card-based layout (mirrors the parlays-tab style) showing every recommended single bet at every tracked sportsbook. The card body renders a strict **2×8 price grid** (2 sides × 8 books) — books in fixed order: WZ, H88, BFA, BKM, B105, DK, FD, Pinn. Spreads, alt-spreads, and moneyline bets now render **both sides** (e.g. `BOS -2.5` and `PHI +2.5`), mirroring how totals already showed Over + Under. Each grid cell is in one of three states:
  - **Exact-line price** — book offers the bet at the model's exact line; shows the American odds.
  - **Mismatched-line** — book offers the bet at a nearby line (within ±3.0 units); the cell shows the line tag in amber.
  - **No-quote** — book has no line for this bet; cell renders muted with a dashed border.
  Above the grid, a green-tinted **hero strip** surfaces the pick book / Fair (de-vigged American odds) / EV / Risk / To Win, plus the `[Place]` and `[Log]` buttons.
  - `[Place]` dispatches by book: Wagerzon → direct REST API (no browser, returns ticket #); Hoop88 / BFA / BetOnlineAG → existing Playwright browser flow; DraftKings / FanDuel / Pinnacle / Bookmaker / Bet105 → button disabled (use `[Log]` instead).
  - `[Log]` records a manual placement without contacting any book — works for any sportsbook.
- **Parlay tab** — MLB-specific: correlated 2-leg parlays (spread + total) priced via `mlb_correlated_parlay.R` with conditional Kelly sizing. Each opportunity renders as a card containing the matchup, legs, a Books pill row (model M plus per-book devigged fair probabilities for DK / FD / PX / NV plus blended consensus Cons), and a metadata strip (Fair / WZ / Size / To Win) with edge percentage and the Place / placed-label / error-pill action. The card layout reads identically across laptop, split-screen, and phone — no column hiding, no horizontal scroll. Combined-parlay selection (the Sel checkbox in the top-right corner of each card) and auto-placement still work unchanged.
- **Kelly Calculator** — manual sizing tool below the settings strip on the
  bets tab. Type American odds + a fair % (or American fair odds), get a
  recommended Risk based on your bankroll × Kelly fraction. Useful when
  you've devigged a market manually and want to size at your own fair.
- **Per-cell devig toggle** — every bet card's price grid has a
  `RAW / FAIR` toggle. FAIR is the default view: each book cell shows
  its own probit-devigged American odds (computed against the book's
  own two-sided quote). Click `RAW` to see the original book prices.
- Kelly sizing with configurable bankroll + multiplier
- Same-game correlation detection with visual tooltips
- Bet placement tracking (placed vs recommended)
- CLV (Closing Line Value) computation post-game via `MLB Answer Key/clv_compute.py`
- Auto-place integration via Playwright for supported books (wagerzon, hoop88, bfa, betonlineag)

## Devigging method

The MLB dashboard's "Books (devigged fair %)" column uses probit (additive z-shift) devigging via `Tools.R::devig_american`. See `docs/superpowers/specs/2026-05-11-probit-devig-design.md`. Historical samples in `mlb_betting_pbp` are also probit-devigged (sharp-weighted across Pinnacle, Bookmaker, LowVig, Circa, Bet105).

## API Endpoints

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/` | GET | Serve dashboard HTML |
| `/api/place-bet` | POST | Dispatch by bookmaker: WZ direct API, Playwright for offshore, 400 for reference books |
| `/api/log-bet` | POST | Manual placement log (no book contact) |
| `/api/remove-bet` | POST | Remove a placed bet |
| `/api/update-bet` | POST | Update actual_size (partial fills) |
| `/api/placed-bets` | GET | Retrieve all placed bets |
| `/api/place-parlay` | POST | Record a placed parlay (MLB-specific) |
| `/api/remove-parlay` | POST | Remove a placed parlay |
| `/api/exposure` | GET | Sum exposure by game |
| `/api/book-settings` | GET/POST | Book enable/disable state |
| `/api/book-settings/bulk` | POST | Enable/disable all books at once |
| `/api/sizing-settings` | GET/POST | Bankroll + Kelly multiplier (separate for singles vs parlays) |
| `/api/filter-settings` | GET/POST | Market / correlation / status filters |
| `/api/clv-compute` | POST | Trigger CLV computation |
| `/api/clv-summary` | GET | CLV aggregated by market/book |
| `/api/clv-details` | GET | Per-bet CLV detail |
| `/api/scheduled-captures` | GET | Upcoming closing-odds snapshot schedule |
| `/api/auto-place` | POST | Legacy Playwright spawn (kept for backwards-compat with the fallback table) |
| `/api/nav-status/<bet_hash>` | GET | Poll navigator status |
| `/refresh` | POST | Re-run the MLB pipeline end-to-end |

## Data Storage

- **mlb_dashboard.duckdb** — dashboard state:
  - `placed_bets` — user-placed singles with status + fill tracking
  - `placed_parlays` — user-placed parlays
  - `book_settings` — which books are enabled
  - `sizing_settings` — bankroll, kelly_mult, parlay_bankroll, parlay_kelly_mult, parlay_min_edge
  - `filter_settings` — UI filter state
  - `closing_snapshots` — odds captured 15 min before first pitch (for CLV)
  - `bet_clv` — post-game CLV computations
- Pipeline bets read from `Answer Keys/mlb_mm.duckdb`: `mlb_bets_combined`, `mlb_bets_book_prices` (new — per-book pill rows for the bets-tab odds screen), `mlb_parlay_opportunities`, `mlb_trifecta_opportunities` — all consumer tables in one DB to avoid contention with the pipeline's long write lock on `mlb.duckdb`
- Historical PBP + odds in `Answer Keys/pbp.duckdb` (`mlb_betting_pbp`, `mlb_betting_history`, `mlb_pbp_all`)

## Running

```bash
# Full pipeline + dashboard (from this directory)
bash run.sh

# Dashboard server only (if pipeline already ran)
cd "Answer Keys/MLB Dashboard"
python3 mlb_dashboard_server.py
# → http://localhost:8083
```

`run.sh` runs the following sequence:
1. Kills any existing process on port 8083
2. Runs `python3 run.py mlb` (sharp scrapers → parallel rec scrapers + MLB.R)
3. Runs `wagerzon_odds/parlay_pricer.py mlb` (fetches exact Wagerzon parlay prices)
4. Runs `mlb_correlated_parlay.R` (correlated parlay edge finder)
5. Runs `mlb_dashboard.R` (generates report.html)
6. Starts `mlb_dashboard_server.py` on port 8083

## CLV Tracking

15 minutes before each game's first pitch, the server's scheduler runs offshore scrapers to snapshot closing odds. After the game completes, `MLB Answer Key/clv_compute.py` fetches Pinnacle's closing snapshot (sharp reference) and computes two CLVs per placed bet:

- **Market CLV**: vs sharp market (Pinnacle at T-15min)
- **Book CLV**: vs the offshore book's own closing odds

Both use normal CDF re-pricing when the closing line differs from the placement line (σ=2.5 for F5 totals, 2.0 for spreads).

## Markets

The Odds API call in `MLB Answer Key/MLB.R` requests 16 markets via `all_deriv_markets` (defined at MLB.R:305):

Full game (5):
- `h2h`, `totals`, `spreads` — FG mains
- `alternate_totals`, `alternate_spreads` — FG alt lines (Wagerzon, Bet105)

First 3 innings (3):
- `h2h_1st_3_innings`, `totals_1st_3_innings`, `spreads_1st_3_innings`

First 5 innings (5):
- `h2h_1st_5_innings` — F5 moneyline
- `totals_1st_5_innings`, `spreads_1st_5_innings` — F5 total + run line
- `alternate_totals_1st_5_innings`, `alternate_spreads_1st_5_innings` — F5 alt lines

First 7 innings (3):
- `h2h_1st_7_innings`, `totals_1st_7_innings`, `spreads_1st_7_innings`

Additional derivative markets (scraped Wagerzon/Bookmaker/BFA/Bet105 data, not via Odds API):
- `odd_even_runs` — full-game total runs odd vs even (Wagerzon only)
- `h2h_3way_1st_5_innings` — F5 3-way moneyline (Wagerzon only) with home/away/tie outcomes; bet_on label is the team name for home/away or `"Tie"` for the draw

Notes:
- F3/F7 moneyline fires only when a book posts F-period MLs; Wagerzon currently posts spread+total only at F3/F7
- F5 alt lines (h1-suffix scrapers like BFA/Bet105 are remapped to f5)

## Configuration

| Setting | Default | Where |
|---------|---------|-------|
| `bankroll` | $100 | `sizing_settings` table |
| `kelly_mult` | 0.25 | `sizing_settings` table |
| `parlay_bankroll` | $100 | `sizing_settings` table |
| `parlay_kelly_mult` | 0.25 | `sizing_settings` table |
| `parlay_min_edge` | 3% | `sizing_settings` table |
| Server port | 8083 | `mlb_dashboard_server.py` |

## Auto-Placement of Correlated Parlays

One-click parlay placement from the Parlays tab. No manual confirmation, no browser opens.

**Workflow:**
1. Parlays tab shows recommended parlay rows with a "Place" button
2. Toggle the amber "Dry run" bar at the top (checked = dry run, unchecked = real placement)
3. Click "Place" → backend calls Wagerzon's REST API → status appears in the Place column

**Dry Run Mode:**
- **Checked (default):** Runs the preflight + drift check but does NOT place. Toast shows `"Dry run OK — would win $X.XX"`
- **Unchecked:** Real placement. Toast shows `"Placed at Wagerzon (#TICKET)"`

**Status Values (visible in Place column):**
- `placed · #<ticket>` — Bet accepted; ticket number from Wagerzon
- `price_moved` — Wagerzon's current price differed by >$0.01; aborted (no money at risk)
- `rejected: <reason>` — Wagerzon refused (balance, size limit, line pulled, etc.)
- `auth_error` — Session expired; re-login retry also failed
- `network_error` — Request incomplete; **VERIFY** Wagerzon's ticket history before retrying (may have been placed)
- `orphaned` — Wagerzon confirmed but local DB write failed; forensics in `placement_orphans` table
- `would_place` — Dry run passed (no money placed)

**Data Storage:**
- `placed_parlays` table in `mlb_dashboard.duckdb` — all placements + status
- `placement_orphans` table — orphaned bets (confirmed at Wagerzon, local write failed)

**Sheets Integration:**
Auto-placed bets are picked up by the existing `bet_logger/scraper_wagerzon.py` on its next run (HistoryHelper feed). No special handling needed — they log like manually-placed bets.

## Wagerzon account selector

The dashboard header includes a row of pills — one per configured Wagerzon
account (discovered by `wagerzon_odds/wagerzon_accounts.py`). Each pill shows
the account label and current available balance. Click a pill to switch the
active placement account; the selection is persisted to
`dashboard_settings.wagerzon_last_used` via `POST /api/wagerzon/last-used`
and used by `POST /api/place-parlay`.

- The selected pill is filled blue. Click an inactive pill (dark grey) to
  switch.
- A pill rendered as `Label · — ⚠` indicates the latest balance fetch
  failed but is still considered fresh (under one minute old). After one
  minute, it gains a `(stale Nm ago)` suffix and turns red-tinted.
- The small `↻` icon next to the pills refetches balances on demand. The
  green **Refresh** button on the right of the title row re-runs the
  dashboard pipeline.
- With zero accounts configured, the row renders a single dashed
  "No Wagerzon accounts configured" pill. Placement is disabled.

## Troubleshooting

- **"Port 8083 already in use"** — `run.sh` kills existing processes, but if you started Flask manually, `lsof -ti:8083 | xargs kill`
- **Parlays tab empty** — `mlb_parlay_opportunities` hasn't been populated; ensure `mlb_correlated_parlay.R` ran successfully (check `run.sh` output)
- **No bets shown** — check `mlb_bets_combined` has rows in `mlb_mm.duckdb`; if empty, check MLB.R output in pipeline logs
- **CLV not computed** — ensure `clv_compute.py` ran post-game; the closing snapshots must exist in `closing_snapshots` table
- **Auto-place fails** — single-leg auto-queue supports `wagerzon`, `hoop88`, `bfa`, `betonlineag`; parlay auto-placement supports `wagerzon` only. Other books must be placed manually.
- **Place button doesn't show** — Browser cache; hard refresh (Cmd+Shift+R)

## Combined Parlay (Wagerzon cash-efficiency)

The Parlay tab supports combining **two recommended parlay rows from different games** into a single 4-leg Wagerzon ticket. Useful when WZ balance is the binding constraint — the combined ticket needs far less cash than placing both source parlays separately.

### How to use

1. On the Parlay tab, check the leftmost checkbox on **two rows from different games**.
2. The "Combined Parlay" banner appears above the table with:
   - Joint fair odds (= product of the two `fair_dec`s, since cross-game legs are independent)
   - Wagerzon's exact 4-leg payout (live `ConfirmWagerHelper` call)
   - Joint edge and recommended Kelly stake (computed from your `parlay_bankroll` and `parlay_kelly_mult` settings)
3. Click **Place Combined →** to record the placement. The page reloads.
4. After reload, the two source rows' Kelly column shows the **conditional residual** — the optimal additional stake on each single given the combo is already placed. Place those residuals as singles too if you want more exposure.

### When to use it

- Wagerzon balance is below the cost of placing both source parlays separately
- You're OK trading EV per dollar bet for cash efficiency

### When NOT to use it

- Wagerzon balance is plentiful — placing both as singles captures more EV
- Same-game combos: blocked, since each game's spread+total combo is already a single Parlay row

### Math

- Combined Kelly stake: `kelly_stake(joint_edge, wz_decimal_payout, parlay_bankroll, parlay_kelly_mult)`
- Conditional residual: numerical max-log-growth optimization (L-BFGS-B) over the 4-outcome joint distribution `(p_A * p_B, p_A * (1-p_B), (1-p_A) * p_B, (1-p_A) * (1-p_B))` with the combo stake fixed
- Implementation:
  - R: `Answer Keys/conditional_kelly.R::conditional_kelly_residuals()`
  - R: `Answer Keys/Tools.R::compute_combined_parlay_pricing()`
  - Python: `Answer Keys/MLB Dashboard/combined_parlay.py::joint_pricing()`
  - Python: `wagerzon_odds/parlay_pricer.py::get_combined_parlay_price()`

### Server endpoints

- `POST /api/price-combined-parlay` — body `{parlay_hash_a, parlay_hash_b}` returns joint pricing + WZ exact payout. 60s in-memory TTL cache keyed on sorted hash tuple.
- `POST /api/place-combined-parlay` — body `{combo_hash, parlay_hash_a, parlay_hash_b, wz_odds, kelly_bet, actual_size, combo_label}` inserts a row in `placed_parlays` with `is_combo = TRUE` and `combo_leg_ids = [hash_a, hash_b]` (JSON).

### One-time setup

After pulling this feature for the first time, run the schema migration once against the live `mlb_dashboard.duckdb`:

```bash
python "Answer Keys/MLB Dashboard/migrations/001_combined_parlay_columns.py" \
    "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb"
```

This adds three columns to `placed_parlays`: `is_combo BOOLEAN`, `combo_leg_ids VARCHAR`, `parent_combo_id INTEGER`. The migration is idempotent — running it twice is safe.

### Out-of-scope (v1)

- N>2 leg combinations
- Combining a parlay with a single from the Bets tab
- Cash-budget-aware portfolio Kelly (set "available WZ balance" → dashboard solves for optimal split)
- Auto-place via WZ API — for now, the dashboard records intent in `placed_parlays`; you place the actual bet at WZ manually. (Single parlays use the auto-placement flow above; combined parlays remain manual for v1.)

## Trifectas tab

Shows priced TRIPLE-PLAY and GRAND-SLAM specials (from Wagerzon) blended with DraftKings SGP fair odds. Mirrors the bets-tab UI (filterable, color-coded EV, in-place Place button) on the parlay-tab plumbing (pricer-writes-table on each refresh, separate placed_X table for dedup).

### Data flow

1. `wagerzon_odds/scraper_specials.py` populates `wagerzon_specials` (sport='mlb', prop_type IN ('TRIPLE-PLAY','GRAND-SLAM'))
2. `Answer Keys/mlb_triple_play.R` runs in parallel with `mlb_correlated_parlay.R` during refresh:
   - reads posted lines from `wagerzon_specials`
   - invokes `mlb_sgp/scraper_draftkings_trifecta.py` for live DK SGP odds
   - blends model fair × DK fair (vig 1.25)
   - writes `mlb_trifecta_opportunities` to `Answer Keys/mlb_mm.duckdb`
3. Dashboard reads `mlb_trifecta_opportunities` and renders the Trifectas tab.

### Placement (manual log)

- Click **Place** to log a trifecta: `POST /api/place-trifecta` with the row's hash; the server fetches the full opportunity row by hash and inserts into `placed_trifectas` (in `mlb_dashboard.duckdb`). Idempotent: re-clicking is a no-op.
- Click **Placed** to undo: `POST /api/remove-trifecta` deletes the row.
- Auto-placement (direct submission to Wagerzon) is **not** wired in this version. Place the bet on Wagerzon yourself, then click Place here to log it.

### Sizing settings

Three rows in `sizing_settings` control trifecta sizing (separate from parlay sizing because trifectas are 3-4 legs at higher payouts with single-book DK blending):

- `trifecta_bankroll` — default 100
- `trifecta_kelly_mult` — default 0.10 (10% Kelly, half the parlay default; trifecta vig 1.25 is conservative until validated)
- `trifecta_min_edge` — default 0.05 (5%; below this the row stays visible but no Place button is shown)

To change these defaults, update the `sizing_settings` table in `mlb_dashboard.duckdb` directly (a UI panel for trifecta sizing has not been built yet — the parlay sliders only persist parlay sizing).
