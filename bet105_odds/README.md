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

Prematch key rotates periodically. Run `recon_bet105.py` to capture fresh params:
it opens a plain Chrome on the persistent profile (`.bet105_profile/`), you log
in yourself, and it attaches over CDP only afterwards (a Playwright launch
crawled: every odds-socket frame was relayed to Python). It also records the
bets API calls My Plays makes to `.bet105_recon_api.json` — the shapes
`unabated_ticket/bets_service/sources/bet105.py` is pinned to. It saves no
cookies; the bets service reads Bet105 through the Unabated Ticket extension.

## Storage

DuckDB: `bet105.duckdb` → tables: `cbb_odds`, `nba_odds` (17-column standard schema — game timing is carried by `game_start_time TIMESTAMPTZ` (UTC); the legacy `game_date VARCHAR` + `game_time VARCHAR` pair was retired. Bet105's source timestamps are already UTC, so the scraper stores them verbatim — no TZ conversion needed.)
