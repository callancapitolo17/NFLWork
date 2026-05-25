# Hoop88 Odds Scraper

Scrapes live odds from Hoop88 (private offshore book) via REST API with JWT authentication.

## Method

REST API:
- Auth: POST `/cloud/api/System/authenticateCustomer` → JWT token
- Odds: POST `/cloud/api/Lines/Get_LeagueLines2` with Bearer token

## Markets Captured

- Main spreads, totals, moneylines
- Team totals (home + away, main only — no alts)
- Multiple periods: FG, 1H, 1Q, 2Q, 3Q, 4Q (varies by sport)

## Sports

- NFL, NCAAF, CBB, NBA

## Usage

```bash
python scraper.py nfl
python scraper.py ncaaf
python scraper.py cbb
python scraper.py nba
```

## Auth

Requires in `.env`:
- `HOOP88_USERNAME`
- `HOOP88_PASSWORD`
- `HOOP88_URL` (default: `https://hoop88.com`)

## Storage

DuckDB: `hoop88.duckdb` → tables: `nfl_odds`, `ncaaf_odds`, `cbb_odds`, `nba_odds` (17-column standard schema — game timing is carried by `game_start_time TIMESTAMPTZ` (UTC); the legacy `game_date VARCHAR` + `game_time VARCHAR` pair was retired).

## Timezone handling

Hoop88's response timestamps are Pacific Time (a 2026-05 audit found PT,
NOT ET as the original code claimed — see
`tools/TZ_AUDIT_FINDINGS.md`). The scraper now converts PT → UTC at
write-time via `ZoneInfo("America/Los_Angeles")` before storing into
`game_start_time TIMESTAMPTZ`, so all downstream consumers can compare
game times directly without worrying about the source TZ.
