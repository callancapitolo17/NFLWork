# MLB SGP Odds Scrapers

Fetch Same Game Parlay (SGP) odds from multiple sources for MLB correlated parlay edge finding.

## The Idea

SGP engines (DraftKings, FanDuel, etc.) price spread+total parlays with built-in correlation adjustments. By comparing their prices to our sample-based fair odds, we get a second independent estimate of true parlay value.

## Sources

| Source | Method | Books Covered | Login Required |
|--------|--------|---------------|----------------|
| **Pikkit Pro** | Browser automation + API intercept | FanDuel, DraftKings, Novig, ProphetX | Yes (SMS) |
| **DraftKings** | TBD (needs recon first) | DraftKings only | No |

## Setup

```bash
cd mlb_sgp
pip install playwright python-dotenv duckdb
playwright install chromium
```

### Pikkit Login (first time)

Pikkit requires phone number + SMS verification:

```bash
python scraper_pikkit_mlb.py --visible
```

1. Browser opens to Pikkit
2. Login with your phone number
3. Complete SMS verification
4. Session saves automatically after 120 seconds

Re-run if session expires.

## Usage

### Pikkit MLB Scraper

```bash
# Scrape all today's games (headless)
python scraper_pikkit_mlb.py

# Show browser for debugging
python scraper_pikkit_mlb.py --visible

# Single game only
python scraper_pikkit_mlb.py --game-id abc123

# Custom trusted books
python scraper_pikkit_mlb.py --books fanduel draftkings
```

### DraftKings Recon (one-time discovery)

Before building the DK scraper, we need to find their SGP pricing API endpoint:

```bash
python recon_draftkings_sgp.py
```

1. Browser opens to DK MLB page
2. Follow the prompts: navigate to a game, add SGP legs
3. Tool captures all network traffic at each step
4. Results saved to `recon_dk_sgp.json`
5. Look for POST endpoints with SGP/parlay keywords in the JSON responses

## Output

Both scrapers write to the `mlb_sgp_odds` table in `Answer Keys/mlb.duckdb`:

| Column | Type | Description |
|--------|------|-------------|
| game_id | VARCHAR | Odds API event ID (joins to mlb_parlay_opportunities) |
| combo | VARCHAR | e.g., "Home Spread + Over" (matches R pipeline combo names) |
| period | VARCHAR | "FG" or "F5" |
| bookmaker | VARCHAR | e.g., "draftkings", "fanduel" |
| sgp_decimal | DOUBLE | Decimal odds |
| sgp_american | INTEGER | American odds |
| fetch_time | TIMESTAMP | When scraped |
| source | VARCHAR | "pikkit" or "draftkings_direct" |

## Files

| File | Purpose |
|------|---------|
| `pikkit_common.py` | Reusable Pikkit functions (session, API intercept, odds parsing) |
| `scraper_pikkit_mlb.py` | MLB Pikkit scraper |
| `recon_draftkings_sgp.py` | DK network recon tool |
| `scraper_draftkings_sgp.py` | DK scraper (TBD after recon) |
| `db.py` | DuckDB read/write helpers for mlb_sgp_odds |

## How Pikkit Works

Pikkit is an odds aggregator. When you build a parlay on their site, their backend calls each sportsbook's SGP engine and returns all the prices in a single API response.

The scraper exploits this by:
1. Clicking spread + total legs on the Pikkit game page (DOM automation)
2. Intercepting the `prod-website.pikkit.app/betslip` API response (Playwright response listener)
3. Parsing the structured JSON for each book's SGP odds
4. Filtering out books that altered the requested line (adjusted leg detection)

The API intercept is the reliable part. The DOM clicking is the fragile part.
