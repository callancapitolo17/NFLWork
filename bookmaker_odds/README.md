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

Cookies cached in `.bookmaker_cookies.json`. If Cloudflare blocks the request AND stdin is a TTY AND the last recon attempt was > 1h ago, the scraper launches `recon_bookmaker.py` in a real Chrome window to refresh `cf_clearance`.

Otherwise (piped subprocess from `run.py`, recent recon attempt, or healthy session with no games posted) the scraper clears stale data, logs the reason, and exits 0 without opening a browser. To refresh cookies manually when the pipeline has skipped it:

```bash
cd bookmaker_odds
./venv/bin/python recon_bookmaker.py
```

## Storage

DuckDB: `bookmaker.duckdb` → tables: `cbb_odds`, `nba_odds` (18-column standard schema)
