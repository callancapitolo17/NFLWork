# New Sportsbook Scraper

Build a new odds scraper for sportsbook: $ARGUMENTS

## Process

### 1. Reconnaissance
- Create `<book>_odds/recon_<book>.py` to explore the sportsbook's API/HTML structure
- Determine approach: **REST API interception** (preferred) or Playwright browser automation
- Identify auth mechanism (token, cookie, API key)
- Map available sports and market types
- Run recon and share findings before proceeding

### 1b. Akamai / Bot Protection
Most US sportsbooks (DraftKings, FanDuel, etc.) use Akamai Bot Manager. Escalation path:
1. **`curl_cffi` with `impersonate="chrome"`** — try this first. Spoofs Chrome TLS fingerprint, bypasses most Akamai. Works on DK's `calculateBets` and `parlays/v1/sgp/events`.
2. **`curl_cffi` + load homepage first** — some endpoints need cookies established by visiting the site: `session.get(homepage)` before API calls.
3. **CDP to real Chrome** — `playwright.chromium.connect_over_cdp("http://localhost:9222")`. Inherits all browser state. Use for endpoints that reject curl_cffi.
4. **Click-and-capture** — last resort. Automate clicking in the browser, intercept API responses via `page.on("response")`.

**Don't waste time on:** cookie transfer from browser to requests, `page.evaluate(fetch())`, Playwright stealth plugins, undetected-chromedriver. These don't work against serious Akamai deployments.

**Key insight:** Different endpoints on the same site can have **different protection levels**. Always test each endpoint individually. DK's `calculateBets` works with curl_cffi but `sportsdata/v2/sgp` doesn't.

### 2. Architecture Decision
Based on recon, choose approach:
- **REST API** (preferred): Use `requests` or `curl_cffi` with intercepted auth tokens
- **Playwright** (fallback): Only if API is not accessible or requires browser-level auth
- Reference existing scrapers for patterns:
  - REST API style: `bfa_odds/scraper.py`, `bookmaker_odds/scraper.py`, `bet105_odds/scraper.py`
  - curl_cffi style: `mlb_sgp/scraper_draftkings_sgp.py` (Akamai bypass)
  - Playwright style: `hoop88_odds/scraper.py`, `wagerzon_odds/scraper_v2.py`

### 2b. SGP / Parlay Scrapers (if applicable)
SGP pricing endpoints are typically the **most protected** — they may need different auth than regular odds endpoints.
- **Find the selection listing endpoint** — every book has an API that returns all available selection IDs. Don't guess ID formats. For DK it's `parlays/v1/sgp/events/{id}` (2MB response with all IDs).
- **Selection IDs vary per game** — suffixes, market numbers, encoding can all change. Build discovery logic, not hardcoded patterns.
- **Test cross-market combos explicitly** — books may restrict which markets can be combined in SGP. A combo that works for one game may fail for another.
- **Restrictions are game-specific and time-dependent** — a combo rejected at 2pm may work at 5pm. Add retry logic.
- **Return no price rather than wrong line** — if the exact Wagerzon total isn't available, don't fall back to a different total. Honest data beats approximate data.
- **Validate with implied probability sum** — across all 4 combos for a game, the sum should be consistent (~1.10 for DK). If it varies wildly, IDs are wrong.

### 3. Build Scraper
Create `<book>_odds/scraper.py` with:
- Sport argument support (`--sport nfl|cbb|nba`)
- DuckDB storage in `<book>_odds/<book>.duckdb`
- Team name normalization using canonical names from The Odds API
- Deduplication on (game, market, book, timestamp)
- Proper error handling and logging
- `.env` for credentials, `.env.example` committed

### 4. Team Name Normalization
- The canonical team names come from The Odds API (used in `run.py` via `fetch_canonical_games()`)
- Map the new book's team names to these canonical names
- For CBB: use fuzzy matching or build explicit mapping table
- Store mapping in `<book>_odds/team_mapping.py` or inline in scraper

### 5. Integrate into Pipeline
- Add entry to `SCRAPER_CONFIGS` in `Answer Keys/run.py`:
  ```python
  "<book>": {
      "script": "../<book>_odds/scraper.py",
      "sports": ["cbb"]  # adjust as appropriate
  }
  ```
- Test with: `cd "Answer Keys" && python run.py --sport cbb`

### 6. Add to Dashboard
- Register the new book in `cbb_dashboard_server.py` book_settings
- Add book column to R dashboard if needed

### 7. Validate
- Run the full pipeline end-to-end
- Verify odds appear in the dashboard
- Check team names match across all books
- Verify no duplicate entries in DuckDB

## Rules
- Use a worktree for this work
- Never commit `.duckdb` files or `.env` files
- Create `.env.example` with placeholder values
- Normalize team names on ingest, not at display time
- Store raw + normalized names in DuckDB for debugging
- Always use `canonical_match.py` for team name resolution — it handles abbreviations, punctuation, state suffixes across all books
- Use `fetchall()` not `fetchdf()` to avoid numpy/pandas dependencies
- Batch DuckDB operations (single DELETE + single INSERT) instead of N individual operations
- Dynamic repo root resolution for worktree compatibility (detect `.worktrees` in path)
- Add session re-init logic after consecutive failures (catches auth/cookie expiry)
- Add retry logic for transient API rejections (sportsbooks temporarily block then allow)
