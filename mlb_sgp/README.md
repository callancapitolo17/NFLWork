# MLB SGP Odds Scrapers

Fetch Same Game Parlay (SGP) odds from DraftKings and FanDuel for MLB correlated parlay edge finding. Both write to the same `mlb_sgp_odds` table so the scanner can shop combos across books.

## Quick Start

```bash
cd mlb_sgp
source venv/bin/activate
python scraper_draftkings_sgp.py           # all games, ~30s
python scraper_fanduel_sgp.py              # all games, ~1s
python scraper_draftkings_sgp.py --verbose
python scraper_fanduel_sgp.py --verbose
```

Requirements: `pip install curl_cffi duckdb` and the MLB pipeline must have run (`mlb_parlay_opportunities` table populated).

## How It Works

Pure REST API — no browser, no Chrome, no clicking. Uses `curl_cffi` with Chrome TLS impersonation to bypass DraftKings' Akamai bot protection.

```
DK League API          → event list (team names, event IDs)
DK Event Markets API   → main market number (Run Line ID)
DK SGP Parlays API     → ALL selection IDs (2MB response, main + alt lines)
Match Wagerzon total   → exact Over/Under selection ID
DK calculateBets       → correlation-adjusted SGP trueOdds
                       → mlb_sgp_odds table in DuckDB
```

## DraftKings API Endpoints

| Endpoint | Auth | Purpose |
|----------|------|---------|
| `sportsbook-nash.../league/leagueSubcategory/v1/markets` | None | List MLB events |
| `sportsbook-nash.../event/eventSubcategory/v1/markets` | None | Main market IDs per game |
| `sportsbook-nash.../parlays/v1/sgp/events/{id}` | curl_cffi | **All selection IDs** (2MB response) |
| `gaming-us-nj.../api/wager/v1/calculateBets` | curl_cffi | **SGP pricing** (POST, returns trueOdds) |
| `sportsbook-nash.../sgp/dkusnj/sportsdata/v2/sgp` | Full Akamai | SGP pricing (DK frontend only — **inaccessible** via REST) |

### Why curl_cffi?

DraftKings uses Akamai Bot Manager. Plain `requests` gets blocked by TLS fingerprinting. `curl_cffi` impersonates Chrome's TLS signature, which is enough to bypass Akamai on `calculateBets` and `parlays/v1/sgp/events`. The `sportsdata/v2/sgp` endpoint has stricter protection and is inaccessible — that's why we use `calculateBets` instead.

**Things that DON'T work:** direct HTTP requests, page.evaluate(fetch()), cookie transfer from browser to requests, Playwright stealth plugins. All tested extensively.

## Selection ID Format

```
Spread: 0HC{market_num}{N|P}{line*100}_{suffix}
Total:  0OU{market_num}{O|U}{line*100}_{suffix}

Examples:
  0HC84191361N150_1   → Home team -1.5, market 84191361, suffix _1
  0HC84191361P150_3   → Away team +1.5, market 84191361, suffix _3
  0OU84203528O750_1   → Over 7.5, alt market 84203528, suffix _1
  0OU84203528U750_3   → Under 7.5, alt market 84203528, suffix _3
```

- **Suffix** (`_1` vs `_3`) varies per game — can't predict, must read from SGP parlays data
- **N** = negative spread (favorite), **P** = positive spread (underdog)
- **Market number** comes from market ID (e.g., `2_84191361` → `84191361`)
- Main market and alt market have **different** market numbers

## DK Market Structure

- **Main market** (subcategory 4519): Run Line (±1.5) + one Total (DK's main line) + Moneyline
- **Alt market**: Alt spreads (±1.0, ±2.5, etc.) + alt totals (every 0.5 from 5.0 to 13.0+)
- Both appear in the SGP parlays response
- **Game line markets** have BOTH spread AND total selections — inning/prop markets have only one. The scraper uses this to filter out non-game markets.

## Known DK Restrictions

### Cross-Market Blocking Near Main Line
DK won't combine the main run line with alt totals within ±0.5 of their main total. Example: if DK's main total is O/U 8.5, you can SGP with O7.5 or O9.5, but NOT O8.0 or O9.0. This is confirmed on DK's website too ("Sorry, your picks cannot be parlayed"). The scraper returns no price for these games rather than using a different total.

### Transient Rejections
Some games temporarily return `SelectionsCannotBeCombined` then work minutes later. The scraper retries once after a 2-second delay.

### SelectionClosed
Games near first pitch may close SGP pricing entirely.

## calculateBets Request/Response

### Request
```json
{
  "selections": [],
  "selectionsForYourBet": [
    {"id": "0HC84191361N150_1", "yourBetGroup": 0},
    {"id": "0OU84203528O750_1", "yourBetGroup": 0}
  ],
  "selectionsForCombinator": [],
  "selectionsForProgressiveParlay": [],
  "oddsStyle": "american"
}
```

### Response (success)
```json
{
  "selectionsForYourBet": [
    {"id": "0HC84191361N150_1", "trueOdds": 2.59, "displayOdds": "+159", "points": -1.5},
    {"id": "0OU84203528O750_1", "trueOdds": 1.87, "displayOdds": "−115", "points": 7.5}
  ],
  "bets": [
    {
      "type": "YourBet",
      "selectionsMapped": [{"id": "0HC84191361N150_1"}, {"id": "0OU84203528O750_1"}],
      "trueOdds": 4.0,
      "displayOdds": "+300"
    }
  ]
}
```

### Error responses
- **422 `SelectionsCannotBeCombined`** — DK won't combine these selections (cross-market restriction)
- **422 `SelectionClosed`** — Game near first pitch, SGP unavailable
- **200 with `combinabilityRestrictions`** — Selections recognized but can't be parlayed

## Output

Writes to `mlb_sgp_odds` table in `Answer Keys/mlb.duckdb`:

| Column | Type | Description |
|--------|------|-------------|
| game_id | VARCHAR | Odds API event ID (joins to mlb_parlay_opportunities) |
| combo | VARCHAR | e.g., "Home Spread + Over" |
| period | VARCHAR | "FG" |
| bookmaker | VARCHAR | "draftkings" |
| sgp_decimal | DOUBLE | Decimal odds |
| sgp_american | INTEGER | American odds |
| fetch_time | TIMESTAMP | When scraped |
| source | VARCHAR | "draftkings_direct" |

## FanDuel Scraper

`scraper_fanduel_sgp.py` mirrors the DK scraper but is much simpler since FD's selection IDs are plain integers tied to marketIds (no DK-style `0HC...` decoding) and FD doesn't lock its pricing endpoint behind Akamai.

### How it works

```
FD scan API           → MLB events list (competitionId 11196870)
canonical_match       → resolve FD team names to internal game_ids
FD event-page API     → Run Line + Total Runs runners (marketId, selectionId)
FD implyBets API      → POST 4 combos per game, parse the betCombinations
                        entry where isSGM=true → that's the correlated SGP price
                      → mlb_sgp_odds (bookmaker='fanduel')
```

### Required headers

FD's API silently strips the SGP combination from `implyBets` responses (returning only single-leg prices) if any of these are missing:

| Header | Value | Notes |
|---|---|---|
| `x-application` | `FhMFpcPWXMeyZxOx` | FD's API key, also passed as `?_ak=` query param |
| `x-sportsbook-region` | `NJ` | Geo header. The `nj.` hostnames are FD's backend routing — works from any state, no VPN required. |
| `x-px-context` | `_pxvid=...;pxcts=...;` | PerimeterX visitor token. Hardcoded in the scraper; long-lived but rotates eventually. |

### Combo logic

Per game, 4 combos × 2 periods (FG + F5) = up to 8 prices. Matches Wagerzon's exact lines:

- Home Spread + Over / Under
- Away Spread + Over / Under
- F5 Home Spread + Over / Under
- F5 Away Spread + Over / Under

**Exact line matching only.** The scraper fetches main + alt spreads and totals from FD's SGP tab, then looks up the exact Wagerzon spread and total. If FD doesn't have the precise line, that game is skipped (no approximate matching). Alt spreads (e.g., ±2.5, ±3.5) are resolved by matching runner team names against the event's home/away teams.

### When the PerimeterX token expires

If FD starts returning 400s with empty bodies, or `implyBets` stops returning the `isSGM=true` entry, the `x-px-context` token has rotated. To refresh:

1. Open `sportsbook.fanduel.com/navigation/mlb` in Chrome with DevTools open
2. Find any `event-page` request in the Network tab
3. Copy the `x-px-context` request header
4. Paste into `FD_PX_CONTEXT` constant in `scraper_fanduel_sgp.py`

A future v2 may bootstrap this automatically via headless Chrome.

## Files

| File | Purpose |
|------|---------|
| `scraper_draftkings_sgp.py` | DK SGP scraper (pure REST, curl_cffi) |
| `scraper_fanduel_sgp.py` | FD SGP scraper (pure REST, curl_cffi) |
| `scraper_pikkit_mlb.py` | Pikkit MLB SGP scraper (fallback) |
| `pikkit_common.py` | Reusable Pikkit functions |
| `recon_draftkings_sgp.py` | DK network recon tool |
| `quick_recon.py` | Lightweight CDP recon (attaches to running Chrome) |
| `db.py` | DuckDB helpers for mlb_sgp_odds |

## Troubleshooting

**"No mlb_parlay_opportunities table"** — Run the MLB pipeline first (`cd "Answer Keys" && python run.py --sport mlb`).

**All games return "no price"** — curl_cffi session may have expired. The scraper auto-reinits after 3 consecutive failures, but if all games fail, try running again.

**"Total X not found in DK selection IDs"** — Wagerzon total doesn't exist on DK for this game. Rare — DK offers totals from 5.0 to 13.0+.

**"Spread ±1.5 not found"** — Game may have been removed from DK or SGP not yet available (too far from game time).
