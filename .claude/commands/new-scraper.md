# New Sportsbook Scraper

Build a new odds scraper for sportsbook: $ARGUMENTS

## Process

### 1. Reconnaissance
- Create `<book>_odds/recon_<book>.py` to explore the sportsbook's API/HTML structure
- Determine approach: **REST API interception** (preferred) or Playwright browser automation
- Identify auth mechanism (token, cookie, API key)
- Map available sports and market types
- Run recon and share findings before proceeding

### 2. Architecture Decision
Based on recon, choose approach:
- **REST API** (preferred): Use `requests` or `curl_cffi` with intercepted auth tokens
- **Playwright** (fallback): Only if API is not accessible or requires browser-level auth
- Reference existing scrapers for patterns:
  - REST API style: `bfa_odds/scraper.py`, `bookmaker_odds/scraper.py`, `bet105_odds/scraper.py`
  - Playwright style: `hoop88_odds/scraper.py`, `wagerzon_odds/scraper_v2.py`

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
