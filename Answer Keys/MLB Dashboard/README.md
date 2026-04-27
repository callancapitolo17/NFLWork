# MLB Dashboard

Interactive web dashboard displaying +EV MLB F5 betting opportunities and correlated parlays.

## Architecture

```
run.py mlb (pipeline)  →  mlb.duckdb/mlb_bets_combined
                      →  mlb.duckdb/mlb_parlay_opportunities
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

- Filterable bet table (book, market, EV threshold, correlation status)
- **Parlay tab** — MLB-specific: correlated 2-leg parlays (spread + total) priced via `mlb_correlated_parlay.R` with conditional Kelly sizing. Each row shows a single "Books" pill cell with our model's joint probability (M), the four per-book devigged fair probabilities (DK / FD / PX / NV), and the blended consensus (Cons) — making model-vs-market disagreement visible at a glance and keeping the table readable on phone and laptop widths.
- Kelly sizing with configurable bankroll + multiplier
- Same-game correlation detection with visual tooltips
- Bet placement tracking (placed vs recommended)
- CLV (Closing Line Value) computation post-game via `MLB Answer Key/clv_compute.py`
- Auto-place integration via Playwright for supported books (wagerzon, hoop88, bfa, betonlineag)

## API Endpoints

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/` | GET | Serve dashboard HTML |
| `/api/place-bet` | POST | Record a placed bet |
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
| `/api/auto-place` | POST | Launch Playwright bet placer |
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
- Pipeline bets read from `Answer Keys/mlb.duckdb` (`mlb_bets_combined`, `mlb_parlay_opportunities`)
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

Currently F5 (first 5 innings) only:
- `h2h_1st_5_innings` — F5 moneyline
- `totals_1st_5_innings` — F5 total + alternate totals
- `spreads_1st_5_innings` — F5 run line

## Configuration

| Setting | Default | Where |
|---------|---------|-------|
| `bankroll` | $100 | `sizing_settings` table |
| `kelly_mult` | 0.25 | `sizing_settings` table |
| `parlay_bankroll` | $4,000 | `sizing_settings` table |
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

## Troubleshooting

- **"Port 8083 already in use"** — `run.sh` kills existing processes, but if you started Flask manually, `lsof -ti:8083 | xargs kill`
- **Parlays tab empty** — `mlb_parlay_opportunities` hasn't been populated; ensure `mlb_correlated_parlay.R` ran successfully (check `run.sh` output)
- **No bets shown** — check `mlb_bets_combined` has rows; if empty, check MLB.R output in pipeline logs
- **CLV not computed** — ensure `clv_compute.py` ran post-game; the closing snapshots must exist in `closing_snapshots` table
- **Auto-place fails** — only `wagerzon` is currently supported for parlay placement; other books must be placed manually
- **Place button doesn't show** — Browser cache; hard refresh (Cmd+Shift+R)
