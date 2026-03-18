# Bookmaker.eu Odds Scraper

Scrapes live odds from Bookmaker.eu using curl_cffi for Cloudflare bypass.

## Method

HTTP POST to `https://be.bookmaker.eu/gateway/BetslipProxy.aspx` with Chrome TLS fingerprint impersonation via curl_cffi.

## Markets Captured

- Main spreads, totals, moneylines only (no alt lines, no team totals)
- Two periods: full game + 1st half

## Sports

- CBB, NBA

## Usage

```bash
python scraper.py cbb
python scraper.py nba
```

## Auth

Requires in `.env`:
- `BOOKMAKER_USERNAME`
- `BOOKMAKER_PASSWORD`

Cookies cached in `.bookmaker_cookies.json`. If blocked/expired, runs `recon_bookmaker.py` to refresh.

## Storage

DuckDB: `bookmaker.duckdb` → tables: `cbb_odds`, `nba_odds` (18-column standard schema)
