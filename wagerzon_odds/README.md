# Wagerzon Odds Scraper

Scrapes live odds from Wagerzon (private offshore book) via REST API with ASP.NET form authentication.

## Method

REST API:
- Auth: ASP.NET form POST (username + password + hidden ViewState fields)
- Odds: `https://backend.wagerzon.com/wager/NewScheduleHelper.aspx`
- Session maintained via `ASP.NET_SessionId` cookie

## Markets Captured

- Main spreads, totals, moneylines (full game + 1H)
- Alt spreads, alt totals (paired, stored separately)
- Team totals (full game + 1H)

## Sports

- NFL, CBB, NBA

## Usage

```bash
python scraper_v2.py nfl
python scraper_v2.py cbb
python scraper_v2.py nba
```

## Auth

Requires in `.env`:
- `WAGERZON_USERNAME`
- `WAGERZON_PASSWORD`

## Files

| File | Purpose |
|------|---------|
| `scraper_v2.py` | Main scraper (current version) |
| `config.py` | Sport configurations (league IDs, API params, table names) |
| `team_mapping.py` | Team name mappings (e.g., "SEA SEAHAWKS" → "Seattle Seahawks") |
| `transform.py` | Transform raw Wagerzon records to standard format |

## Storage

DuckDB: `wagerzon.duckdb` → tables: `nfl_odds`, `cbb_odds`, `nba_odds` (18-column standard schema)

Also: `wagerzon_cbb.duckdb` (CBB-specific historical data)
