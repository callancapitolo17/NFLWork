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

## SGP Scraping Playbook (for adding new books)

Lessons learned from building DK and FD scrapers. Use this when expanding to new sites (BetMGM, Caesars, ESPN BET, etc.).

### The 3 endpoints you need to find

Every book has these, just named differently:

| What | DK | FD | What to look for |
|---|---|---|---|
| **Event listing** | `league/.../v1/markets` | `scan/.../facet/.../search` | Returns today's games with event IDs + team names |
| **Market catalog** | `parlays/v1/sgp/events/{id}` | `sbapi/event-page?eventId=X&tab=same-game-parlay-` | Returns ALL selection IDs (main + alt) for one game |
| **SGP pricing** | `wager/v1/calculateBets` | `fixedodds/transactional/v1/implyBets` | POST with 2 selection IDs, returns correlated price |

### Recon process

1. Open the book in Chrome with DevTools → Network tab
2. Build a 2-leg SGP manually (spread + total)
3. Capture every network request during "add 2nd leg" — the SGP pricing call fires here
4. **Capture request HEADERS, not just URLs.** Missing headers was a multi-hour debugging session on FD. The API returned 200 with degraded data (singles only) instead of an error.
5. Try replaying the call with `curl_cffi` → if it works, you're done. If 400/403, check for missing headers or bot protection tokens.

### Common gotchas

**Silent degradation > explicit errors.** FD returns 200 with single-leg prices when headers are wrong, instead of 403. DK returns 422 with a clear error code. Assume the worst: always verify the response contains the SGP combined entry, not just a 200 status.

**Alt markets may hide behind a different tab/category.** FD's default event-page returns 3-5 markets. Adding `&tab=same-game-parlay-` returns 156+. Always check if there's a tab/category parameter that unlocks more markets.

**Selection ID formats vary wildly.** DK encodes market number + sign + line + suffix into strings like `0HC84191361N150_1`. FD uses plain integer `selectionId` tied to a `marketId`. Don't assume one format.

**Team disambiguation in alt spread markets.** Alt run line markets list BOTH teams at each line (e.g., "Reds -1.5" AND "Marlins -1.5"). You need to match runner team names against known home/away to avoid pricing the wrong team's spread. Key by `("home", line)` / `("away", line)`, not by handicap sign.

**Live games return adjusted handicaps.** A game in progress shows Run Line at +2.5/-2.5 instead of ±1.5. Filter events with `openDate < now` before scraping. Two events can exist for the same matchup (live + tomorrow's pre-game).

**PerimeterX / bot protection tokens.** Some books require a `x-px-context` or similar token that's set by JavaScript on page load. `curl_cffi` impersonation alone isn't enough — you need the token too. These tokens are semi-persistent (days/weeks) but eventually rotate. Hardcode for v1, add auto-refresh later if needed.

### Vig measurement

For any book, compute vig from the 4 mutually-exclusive combos per game:

```
vig = sum(1/D for all 4 combos: HomeSpread+Over, HomeSpread+Under, AwaySpread+Over, AwaySpread+Under)
```

This sum exceeds 1.0 by the book's vig charge. Measured values:
- **DK:** stable ~1.125 (12.5% vig), consistent across time-to-game
- **FD:** bimodal — ~1.13 for games >21h out, ~1.22 for games <16h out

**Don't hardcode a single vig constant** when the book's vig varies. Use per-game devigging: divide each combo's implied prob by the per-game sum. Falls back to a constant when <4 combos are available.

**FG and F5 are separate partitions.** Never sum across both periods — that doubles the measured vig. Group by `(game_id, period)`.

### Exact line matching

The scraper must match Wagerzon's exact spread and total lines. If the book doesn't have the precise line (e.g., Wagerzon has total 8.5 but the book only offers 8.0), skip the game entirely. No rounding, no closest-line fallback. The parlay fair odds are calibrated to Wagerzon's specific lines — a different line is a different bet.

### Integration into the R scanner

New books write to `mlb_sgp_odds` with their own `bookmaker` and `source` values. The R scanner (`mlb_correlated_parlay.R`) reads all books, computes per-game vig for each, and blends their devigged fair probs equally with the model. Adding a third book requires:

1. Add `'newbook_direct'` to the `WHERE source IN (...)` query
2. Add a `NEWBOOK_SGP_VIG_DEFAULT` fallback constant
3. Add a `nb_vig_lookup` table (same pattern as DK/FD)
4. Add the `nb_row` / `nb_fair_prob` block in the per-combo loop (copy DK or FD block)
5. Optionally add `nb_fair_prob` to the output tibble

The blend automatically scales: `mean(model, dk, fd, nb)` when all present, falls back gracefully when books are missing.

---

## FanDuel Scraper

`scraper_fanduel_sgp.py` mirrors the DK scraper but is much simpler since FD's selection IDs are plain integers tied to marketIds (no DK-style `0HC...` decoding) and FD doesn't lock its pricing endpoint behind Akamai.

### How it works

```
FD scan API           → MLB events list (competitionId 11196870)
canonical_match       → resolve FD team names to internal game_ids
FD event-page API     → all SGP-eligible runners (main + alt, FG + F5)
                        via &tab=same-game-parlay- (156+ markets)
Match Wagerzon lines  → exact spread + total lookup
FD implyBets API      → POST each combo, parse isSGM=true entry
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

### FD-specific lessons learned

- **Event name format:** `"Away Team (P Pitcher) @ Home Team (P Pitcher)"` — strip the pitcher parens before team matching.
- **`nj.` hostnames are NOT geo-restricted.** `sib.nj.sportsbook.fanduel.com` works from California with no VPN. The `nj` prefix is FD's backend routing, not a geo gate.
- **`implyBets` returns 3 entries for a 2-leg combo:** 2 `SINGLE` (one per leg, `isSGM=false`) and 1 `DOUBLE` (`isSGM=true` — this is the SGP price). Parse `winAvgOdds.trueOdds.decimalOdds.decimalOdds` for the decimal value.
- **FD's SGP vig is bimodal by time-to-game:** ~13% for games >21h out, ~21% for games <16h out. The step change happens around 16-21h before first pitch — possibly tied to lineup posting windows.
- **Alt total runner names use parens:** `"Over (8.5)"`, `"Under (7.5)"`. Main total runners are just `"Over"` / `"Under"` with the line in the `handicap` field.
- **Alt spread runner names embed the team:** `"Cincinnati Reds +3.5"`. The `handicap` field is 0 for alts. Parse team name + signed line from the string, match team to home/away.
- **F5 totals at integer values (e.g., 5.0) may not exist.** FD's F5 alt totals jump in 1.0 increments (2.5, 3.5, 4.5, 5.5...). If Wagerzon has F5 total 5.0, FD won't have it and the game is correctly skipped.

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
