# Bet Logger

Multi-platform bet history scraper. Scrapes settled bets from Wagerzon, Hoop88, BFA Gaming, and BetOnline, then uploads to Google Sheets with automatic duplicate detection.

## Setup

### 1. Python Environment

```bash
cd "/Users/callancapitolo/NFLWork/bet_logger"
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
playwright install chromium
```

### 2. Google Sheets Credentials

1. Go to [Google Cloud Console](https://console.cloud.google.com/)
2. Enable the **Google Sheets API**
3. Create a **Service Account** and download the JSON key as `credentials.json` in this folder
4. Share your Google Sheet with the service account email (Editor access)

### 3. Environment Variables

```bash
cp .env.example .env
```

Edit `.env` with credentials for each platform (Wagerzon, Hoop88, BFA) and your Google Sheet ID.

## Usage

### Run all scrapers

```bash
./run_all_scrapers.sh
```

Runs Wagerzon, Hoop88, BFA, and BetOnline sequentially. Continues to the next scraper if one fails.

### Run individual scrapers

```bash
./venv/bin/python3 scraper_wagerzon.py
./venv/bin/python3 scraper_hoop88.py    [--weeks 1] [--visible] [--dry-run]
./venv/bin/python3 scraper_bfa.py       [--weeks 1] [--visible] [--dry-run]
./venv/bin/python3 scraper_betonline.py [--days 7] [--since-last] [--dry-run]
```

**Common flags:**
- `--weeks N` — How many weeks back to fetch (0=this week, 1=last week, default: 1)
- `--visible` — Show the browser window instead of running headless
- `--dry-run` — Scrape but don't upload to Google Sheets
- `--test` — Parse from a saved HTML/JSON file (offline testing)

**BetOnline-specific flags:**
- `--days N` — Days of history to fetch (default: 7)
- `--since-last` — Scrape from the date of the last BetOnline bet in Google Sheets
- `--refresh-only` — Just refresh the auth token and exit (keeps token alive)

### Automated Monday runs (launchd)

The scrapers are configured to run automatically every Monday at 5 AM via macOS launchd. Logs go to `bet_logger/logs/`. To check status:

```bash
launchctl list | grep betlogger
```

To manually load/unload:

```bash
launchctl load ~/Library/LaunchAgents/com.callancapitolo.betlogger.plist
launchctl unload ~/Library/LaunchAgents/com.callancapitolo.betlogger.plist
```

## Output Columns (Google Sheets)

| Column | Field | Example |
|--------|-------|---------|
| A | Date | 1/6/26 |
| B | Platform | Wagerzon / Hoop88 / BFA / BetOnline |
| C | Sport | NFL, NBA, NHL, NCAAF, NCAAM |
| D | Description | Houston Texans -0.5 +130 - 1st Quarter |
| E | Bet Type | Parlay, Straight, Prop, Teaser, Contest |
| F | Line | +3.5, Over 45.5 |
| G | Odds | +450, -110 (American format) |
| H | Bet Amount | $100.00 |
| I | Decimal Odds | 5.50 |
| J | Result | win / loss / push |

Duplicate detection uses (date + description + bet amount) as the key.

## Platform Notes

**Wagerzon** — Scrapes the HTML history table at `backend.wagerzon.com/wager/History.aspx`. Auto-login with credentials, falls back to 60-second manual login window if needed.

**Hoop88** — Navigates the bet history UI by clicking the Balance box, selecting a week from the dropdown, and expanding bet details. Handles parlays by combining legs with `|`. Converts fractional lines (1/2, 1/4, 3/4).

**BFA Gaming** — Uses Keycloak SSO authentication (redirects to `auth.bfagaming.com`). Blazor SPA with date range filters. Parses teasers with adjusted lines.

**BetOnline** — Browser-free API scraper. Uses a one-time recon script (`recon_betonline.py`) to capture auth tokens via Chrome, then hits BetOnline's REST API directly with `requests`. Refresh token rotates on each use and expires after 3 days of inactivity. A daily LaunchAgent (`com.callancapitolo.betonline-token-refresh`) runs `--refresh-only` at noon to keep the token alive. If the token expires, re-run `recon_betonline.py` to re-authenticate.

## Troubleshooting

- **Login fails** — Check credentials in `.env`. If auto-login fails, run with `--visible` to log in manually.
- **No bets found** — Check the time period filter. Verify bets exist on the site for that week.
- **Duplicates still appearing** — Detection matches on exact (date, description, amount). Any small difference bypasses it.
- **Permission denied on Sheets** — Share the spreadsheet with the service account email from `credentials.json`.

## Architecture

```
run_all_scrapers.sh          # Entry point (runs all 4)
  scraper_wagerzon.py        # Wagerzon scraper (headless Chromium)
  scraper_hoop88.py          # Hoop88 scraper (headless Chromium)
  scraper_bfa.py             # BFA Gaming scraper (headless Chromium)
  scraper_betonline.py       # BetOnline scraper (REST API, no browser)
recon_betonline.py           # One-time token capture via Chrome
utils.py                     # Shared: odds calc, sport detection, parsing
sheets.py                    # Google Sheets upload + duplicate detection
```
