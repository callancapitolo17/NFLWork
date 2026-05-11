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

## Devigging method

The CBB dashboard's "Books (devigged fair %)" column uses probit (additive z-shift) devigging via `Tools.R::devig_american`. See `docs/superpowers/specs/2026-05-11-probit-devig-design.md`. Historical samples in `cbb_betting_pbp` are also probit-devigged (sharp-weighted across Pinnacle, Bookmaker, LowVig, Circa, Bet105) — the prior 0.5 prob hardcode is gone, so historical samples now match the live consensus methodology.

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

# Dashboard server only (if pipeline already ran)
python cbb_dashboard_server.py
# → http://localhost:8082
```

## CLV Tracking

The dashboard schedules offshore scraper runs 15 minutes before tipoff to capture closing odds. Post-game, `clv_compute.py` compares placed bet odds against closing consensus to measure CLV — the ultimate validation metric for edge.
