# CBB Dashboard

Interactive web dashboard displaying +EV college basketball betting opportunities.

## Architecture

```
CBB.R (pipeline output)  →  cbb_dashboard.R (HTML generation)  →  cbb_dashboard_server.py (Flask)
                                                                         ↓
                                                                    Browser (:8082)
```

- **cbb_dashboard.R**: Reads pipeline bets, generates reactive HTML with reactable tables
- **cbb_dashboard_server.py**: Flask server serving static HTML + REST API for state management
- **run.sh**: Full pipeline launcher (scrapers → R model → dashboard → Flask)

## Features

- Filterable bet table (by book, market, EV threshold, correlation status)
- Kelly sizing with configurable bankroll + multiplier
- Same-game correlation detection with visual tooltips
- Bet placement tracking (placed vs recommended)
- CLV (Closing Line Value) computation post-game
- Auto-place integration via Playwright (spawns bet_placer)

## API Endpoints

| Endpoint | Method | Purpose |
|----------|--------|---------|
| `/api/place-bet` | POST | Record a placed bet |
| `/api/remove-bet` | POST | Remove a placed bet |
| `/api/update-bet` | POST | Update bet details |
| `/api/placed-bets` | GET | Retrieve all placed bets |
| `/api/exposure` | GET | Sum exposure by game |
| `/api/book-settings` | GET/POST | Book enable/disable state |
| `/api/sizing-settings` | GET/POST | Bankroll + Kelly multiplier |
| `/api/filter-settings` | GET/POST | Market/correlation/status filters |
| `/api/clv-compute` | POST | Trigger CLV computation |
| `/api/clv-summary` | GET | CLV by market/book |
| `/api/clv-details` | GET | Per-bet CLV detail |
| `/api/auto-place` | POST | Launch Playwright bet placer |

## Data Storage

- **cbb_dashboard.duckdb**: placed_bets, book_settings, sizing_settings, filter_settings, closing_snapshots
- Pipeline bets read from `cbb.duckdb` via R

## Running

```bash
# Full pipeline + dashboard
bash run.sh

# With Kalshi tournament prop edge finder
bash run.sh --kalshi-edges

# Dashboard server only (if pipeline already ran)
python cbb_dashboard_server.py
# → http://localhost:8082
```

## Kalshi Edges Tab

The `--kalshi-edges` flag enables tournament prop analysis:

- **kalshi_edges.R**: Runs 10,000 Monte Carlo bracket simulations, fetches all Kalshi March Madness prop markets (round advancement, seed sum, highest seed, upset count, seed win), maps each to a sim probability, and computes EV after Kalshi's 7% taker fee
- EV calculation matches `kalshi_mm/taker.py` exactly: `fee = 0.07 * P * (1-P) * 100`, `ev% = (fair - price - fee) / price`
- Markets covered: team round advancement (R32 through Championship), seed sum (F4/Title), seed count per round, highest seed, upset count, seed wins

## CLV Tracking

The dashboard schedules offshore scraper runs 15 minutes before tipoff to capture closing odds. Post-game, `clv_compute.py` compares placed bet odds against closing consensus to measure CLV — the ultimate validation metric for edge.
