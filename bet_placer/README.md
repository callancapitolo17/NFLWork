# Bet Placer

Automated bet placement system using Playwright browser automation. Places bets across 4 sportsbooks from the dashboard with a single click.

## Architecture

```
Dashboard "Place Bet" button
        ↓
cbb_dashboard_server.py (/api/auto-place)
        ↓
placer.py (CLI entry point)
        ↓
BaseNavigator.lookup_game()  →  DuckDB (find game ID)
        ↓
BookNavigator.place_bets()   →  Playwright browser (visible, not headless)
        ↓
User confirms in browser     →  Status updated in dashboard
```

## Supported Books

| Book | Navigator | Auth Method | Browser Profile |
|------|-----------|-------------|-----------------|
| BFA Gaming | `navigator_bfa.py` | Keycloak cookies | `.bfa_nav_profile` |
| Hoop88 | `navigator_hoop88.py` | REST JWT + cookies | `.hoop88_nav_profile` |
| Wagerzon | `navigator_wagerzon.py` | ASP.NET form POST | New session each run |
| BetOnline | `navigator_betonline.py` | Keycloak OAuth2 | `.betonline_nav_profile` |

## Usage

```bash
# Single bet
python placer.py '{"bookmaker": "bfa", "home_team": "Duke", "away_team": "UNC", "market": "spreads_fg", "bet_on": "Duke", "line": -4, "odds": -110, "recommended_size": 50}'

# Batch (same book)
python placer.py '[{...}, {...}]'

# Custom timeout (default 300s)
python placer.py '<json>' --timeout 120
```

## Bet Status Flow

```
dashboard creates  →  "navigating"  →  "ready_to_confirm"  →  "pending" (user confirmed)
                                    →  "nav_error" (couldn't find game/odds)
                                    →  "nav_timeout" (user didn't confirm in time)
```

## Key Design Decisions

- **Visible browser**: User must manually confirm the bet (human-in-the-loop)
- **Scoring system**: Each navigator scores potential odds buttons (line match, position, text match) and only clicks if confidence >= 3 points
- **Batch support**: Groups bets by period to minimize page navigation within a single browser session
- **Persistent profiles**: Browser profiles survive across runs (no re-login each time)

## Dependencies

- `playwright` (Chromium automation)
- `requests` (API pre-fetch for game lookup)
- `duckdb` (read game data from odds databases)
- `python-dotenv` (credentials from `bet_logger/.env`)
- Real Chrome required for BetOnline (Cloudflare bypass)

## Files

| File | Purpose |
|------|---------|
| `placer.py` | CLI entry point, dispatches to correct navigator |
| `base_navigator.py` | Shared utilities: game lookup, status updates, market parsing |
| `navigator_bfa.py` | BFA Gaming placement logic |
| `navigator_hoop88.py` | Hoop88 placement logic |
| `navigator_wagerzon.py` | Wagerzon placement logic |
| `navigator_betonline.py` | BetOnline placement logic |
