# Bet105 Odds Scraper

Scrapes live odds from Bet105.ag (LinePros white-label platform) via WebSocket.

## Method

Connects to `wss://pandora.ganchrow.com/socket.io/` using Socket.IO v4 protocol. Receives gzip-compressed binary payloads containing event data and coefficient updates.

## Markets Captured

- Main spreads, totals, moneylines (full game + 1H)
- Alt spreads (±15 pts), alt totals (±25 pts)
- Team totals (home + away, main + alts)

## Sports

- CBB (NCAA Div I + Extra)
- NBA

## Usage

```bash
python scraper.py cbb
python scraper.py nba
```

## Auth

Requires in `.env`:
- `BET105_PARTNER_ID`
- `BET105_PREMATCH_KEY`
- `BET105_USER_ID`
- `BET105_GROUP_ID`

Prematch key rotates periodically. Run `recon_bet105.py` to capture fresh params from browser.

## Storage

DuckDB: `bet105.duckdb` → tables: `cbb_odds`, `nba_odds` (18-column standard schema)
