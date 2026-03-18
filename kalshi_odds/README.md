# Kalshi Odds Scraper

Scrapes live odds from Kalshi prediction market exchange via public REST API.

## Method

REST API: `https://api.elections.kalshi.com/trade-api/v2/markets` — no authentication required.

## Markets Captured

- 1H spreads (from individual "Team wins by >X" contracts)
- 1H totals (from "Total >X" contracts)
- 1H moneyline / winner (3-way: team A, team B, tie)

Kalshi prices are probabilities (0-100 cents). Converted to American odds with 7% taker fee baked in.

## Sports

- CBB only (currently)

## Usage

```bash
python scraper.py cbb
```

## Auth

None required (public API).

## Storage

DuckDB: `kalshi.duckdb` → table: `cbb_odds` (26-column extended schema — includes `*_cents` fields for effective probability after fee)

## Series Tickers

- `KXNCAAMB1HSPREAD` — 1H spreads
- `KXNCAAMB1HTOTAL` — 1H totals
- `KXNCAAMB1HWINNER` — 1H winner

## Notes

- Liquidity filter: bid-ask spread must be ≤20 cents
- Each spread contract is a separate market; R pipeline picks best line per event
