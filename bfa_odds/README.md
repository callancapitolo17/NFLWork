# BFA Gaming Odds Scraper

Scrapes live odds from BFA Gaming via public REST API. No authentication required.

## Method

REST API endpoints:
- `/oddsservice/events/popular/{sport}` (primary)
- `/oddsservice/events/leagues` (fallback)

## Markets Captured

- Main spreads, totals, moneylines
- Alt spreads (±15 pts), alt totals (±25 pts)
- Team totals (home + away)
- Multiple periods: FG, 1H, 1Q, 2Q, 3Q, 4Q (varies by sport)

## Sports

- NFL, NBA, CBB

## Usage

```bash
python scraper.py nfl
python scraper.py nba
python scraper.py cbb
```

## Auth

None required (public API).

## Storage

DuckDB: `bfa.duckdb` → tables: `nfl_odds`, `nba_odds`, `cbb_odds` (18-column standard schema)
