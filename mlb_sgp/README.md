# MLB SGP Odds Scrapers

Fetch Same Game Parlay (SGP) odds from multiple sources for MLB correlated parlay edge finding.

## Sources

| Source | Method | Speed | Status |
|--------|--------|-------|--------|
| **DraftKings** | CDP click-and-capture | ~15s/game | Working (verified) |
| **Pikkit Pro** | API intercept | ~10s/game | Ready to test |

## DraftKings SGP Scraper

Uses Chrome DevTools Protocol (CDP) to connect to your running Chrome browser. Clicks the Run Line + Total cells in DK's SGP builder and intercepts the `calculateBets` API response.

### Setup

```bash
# 1. Start Chrome with remote debugging
/Applications/Google\ Chrome.app/Contents/MacOS/Google\ Chrome --remote-debugging-port=9222

# 2. Navigate to DK MLB page in that Chrome window
# https://sportsbook.draftkings.com/leagues/baseball/mlb

# 3. Run the scraper
cd mlb_sgp
source venv/bin/activate
python scraper_draftkings_sgp.py --verbose
```

### Why CDP?

DraftKings uses Akamai Bot Manager. The SGP pricing endpoint returns 404 for:
- Direct HTTP requests (even with curl_cffi Chrome impersonation)
- `page.evaluate(fetch())` from inside Playwright
- Cookie transfer from browser to requests session

Only DK's own JavaScript can call the SGP endpoint — it includes Akamai sensor tokens. The click-and-capture approach lets DK's JS make the call, we just intercept the response.

### Verified Result

Cubs -1.5 + Over 7.5 → `trueOdds=3.64` (+264)
- Individual legs: -1.5 (+129), O 7.5 (-115)
- Independent multiply would give ~+340
- DK's correlation adjustment: ~20% tighter

## Pikkit Pro Scraper

Pikkit is an odds aggregator that shows SGP odds from multiple books (FanDuel, DraftKings, Novig, ProphetX) in a single API response.

### Setup

```bash
# First time: SMS login required
cd mlb_sgp
source venv/bin/activate
python scraper_pikkit_mlb.py --visible
# Login manually, session saves automatically

# Subsequent runs
python scraper_pikkit_mlb.py
```

## DK Recon Tool

Network traffic capture tool used to discover DK's SGP API endpoints.

```bash
python recon_draftkings_sgp.py
```

## Output

Both scrapers write to `mlb_sgp_odds` table in `Answer Keys/mlb.duckdb`:

| Column | Type | Description |
|--------|------|-------------|
| game_id | VARCHAR | Event ID |
| combo | VARCHAR | e.g., "Home Spread + Over" |
| period | VARCHAR | "FG" |
| bookmaker | VARCHAR | "draftkings", "fanduel", etc. |
| sgp_decimal | DOUBLE | Decimal odds |
| sgp_american | INTEGER | American odds |
| fetch_time | TIMESTAMP | When scraped |
| source | VARCHAR | "draftkings_direct" or "pikkit" |

## Files

| File | Purpose |
|------|---------|
| `scraper_draftkings_sgp.py` | DK SGP scraper via CDP click-and-capture |
| `scraper_pikkit_mlb.py` | Pikkit MLB SGP scraper |
| `pikkit_common.py` | Reusable Pikkit functions (session, API intercept) |
| `recon_draftkings_sgp.py` | DK network recon tool |
| `db.py` | DuckDB helpers for mlb_sgp_odds |
