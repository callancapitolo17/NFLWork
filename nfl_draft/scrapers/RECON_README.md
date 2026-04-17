# NFL Draft Reconnaissance Scripts

Four scripts that capture the raw API responses powering each book's NFL Draft
futures page. The captured JSON is saved to
`nfl_draft/tests/fixtures/<book>/draft_markets.json` so parser subagents can
work **offline** from fixtures (no live book calls during parser development).

## Python

All four scripts run under:
```bash
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python nfl_draft/scrapers/recon_<book>.py
```

Dependencies installed in that venv: `curl_cffi`, `playwright`, `requests`,
`python-dotenv`.

## Credentials

The BM and WZ scripts read from `bet_logger/.env`:
```
BOOKMAKER_USERNAME=...
BOOKMAKER_PASSWORD=...
WAGERZON_USERNAME=...
WAGERZON_PASSWORD=...
```
DK and FD require no credentials (public APIs + TLS fingerprint bypass).

## Per-Book Workflow

### DraftKings — `recon_dk.py`

- **Auth:** none (public DK content API, curl_cffi `impersonate='chrome'`).
- **Endpoint:** `sportsbook-nash.draftkings.com/sites/US-SB/api/sportscontent/dkusnj/v1/leagues/88808`
  and each `2026 Draft` subcategory (categoryId=1803).
- **Run:** `python nfl_draft/scrapers/recon_dk.py`
- **Expected output:** ~600KB JSON, 14 subcategories, 59 markets. Fixture at
  `nfl_draft/tests/fixtures/draftkings/draft_markets.json`.
- **Fallback:** `--browser` opens headed Chrome and captures the biggest
  sportscontent response during page load (useful if Akamai throttles REST).

### FanDuel — `recon_fd.py`

- **Auth:** none, but requires three headers (`x-application`,
  `x-sportsbook-region`, `x-px-context`). The PX token in the script is valid
  as of 2026-04-17; rotate if requests start returning HTTP 400 + empty body.
- **Endpoint:** `api.sportsbook.fanduel.com/sbapi/content-managed-page?page=CUSTOM&customPageId=nfl`
- **Run:** `python nfl_draft/scrapers/recon_fd.py`
- **Expected output:** ~750KB JSON. The script finds `layout.tabs['391']`
  (title `NFL Draft`) and reports 16 draft markets of 69 total on the NFL
  page. Fixture at `nfl_draft/tests/fixtures/fanduel/draft_markets.json`.
- **Fallback:** `--browser` opens Chrome; user can manually reach the draft
  page and the script grabs the largest event-page/content-managed-page response.
- **Refreshing x-px-context:** open `sportsbook.fanduel.com` in Chrome
  DevTools, filter requests by `event-page`, copy the `x-px-context` header
  from any one of them into `FD_PX_CONTEXT` at the top of `recon_fd.py`.

### Wagerzon — `recon_wz.py`

- **Auth:** ASP.NET `__VIEWSTATE` POST (same as `wagerzon_odds/scraper_v2.py`).
- **Endpoint:** `backend.wagerzon.com/wager/NewScheduleHelper.aspx?WT=0&lg=<comma-list>`
- **Run:** `python nfl_draft/scrapers/recon_wz.py`
- **Expected output:** ~870KB JSON with 32 active "NFL DRAFT 2026" league
  buckets. The script reports each `IdLeague` and its game count so you can
  see which draft markets the book is actually carrying today. Fixture at
  `nfl_draft/tests/fixtures/wagerzon/draft_markets.json`.
- **Confirmed active lg IDs (2026-04-17):**
  Overall-pick props: `1270,2867,1370,1867,2611,3699,3700,3701,3702,3703`
  (#1-10); `2614` (#11-15); `2618` (Mr. Irrelevant).
  Bundles: `4242` (Top 5), `4243` (Top 10), `4245` (1st Round), `4246`
  (Draft Specials). Live: `1269`. Plus position/team props: `2579, 478, 590,
  977, 1231, 3237, 2516, 2560-2565, 2580-2581, 4537`.
- **Options:**
  - `--probe N` — single-ID debug: fetches just lg=N and dumps league descriptions.
  - `--browser` — headed Chrome fallback; user navigates and the script captures
    the NewScheduleHelper response fired by the page.

### Bookmaker — `recon_bm.py`

- **Auth:** cookie file `bookmaker_odds/.bookmaker_cookies.json` (managed by
  `bookmaker_odds/recon_bookmaker.py`) + ASP.NET Login POST.
- **Endpoint:** `be.bookmaker.eu/gateway/BetslipProxy.aspx/GetSchedule`
- **Status:** **BM's REST probe does not expose NFL Draft to the dummy
  account** — a parallel sweep of IDs 100-1500 (all IDs, no gaps) surfaced
  only 16 non-empty leagues, of which just `104/105` (NFL 1H/2H) and
  `117/208/211` (TNT "ODDS TO WIN", all empty) were NFL-adjacent. No leagues
  had "NFL DRAFT" in the Description. BM may gate futures behind a funded
  account, or NFL Draft markets may be offline until closer to draft day.
- **Recommended workflow:**
  1. Refresh cookies: `python bookmaker_odds/recon_bookmaker.py`
  2. Run `python nfl_draft/scrapers/recon_bm.py` — this will try the short
     REST probe, then fall back to a headed Chrome window.
  3. **In the browser: log in if prompted, then click Football -> NFL ->
     Odds to Win / NFL Futures / NFL Draft.** The script listens for
     `GetSchedule` responses with >2KB bodies and real league data. When
     the draft markets are visible, press ENTER in the terminal.
  4. The script prints the discovered `league_id` extracted from the POST
     body's `LeaguesIdList` — record this ID so the real scraper can hardcode
     it (similar to how `bookmaker_odds/scraper.py` hardcodes cbb=4, nba=3, etc.).
- **Options:**
  - `--probe N` — probe a single league ID (e.g. `--probe 117` for the
    `ODDS TO WIN` TNT league observed during recon).
  - `--browser` — skip REST, go straight to the browser.
- **Fixture at** `nfl_draft/tests/fixtures/bookmaker/draft_markets.json`.

## Troubleshooting

| Symptom | Book | Fix |
|---|---|---|
| HTTP 403 from DK REST | DK | Akamai rate-limit; wait 5 minutes, or run `--browser` |
| HTTP 400 empty body from FD | FD | `x-px-context` token expired — refresh from Chrome DevTools into `FD_PX_CONTEXT` |
| "WZ login failed" | WZ | Check `WAGERZON_USERNAME`/`WAGERZON_PASSWORD` in `bet_logger/.env` |
| BM REST probe returns no hits | BM | Expected — the dummy account doesn't expose draft leagues via REST. Use browser fallback. |
| BM browser captures `GetConfig` or empty schedule | BM | The filter now rejects responses <2KB and responses with empty `League[]` arrays. If you still get nothing, you didn't navigate to a page that actually shows draft odds. |
| Playwright not installed | any | `pip install playwright && playwright install chromium` into the kalshi_draft venv |

## Re-running

Every script is idempotent: re-running overwrites the fixture. The fixture
envelope includes `captured_at` and `meta` fields so you can see when and how
each fixture was produced:

```json
{
  "book": "draftkings",
  "captured_at": "2026-04-17T20:06:00Z",
  "meta": { "method": "rest_controldata", "league_id": 88808, ... },
  "data": { ... raw book response ... }
}
```

## Guardrails

- Fixtures **do not** contain auth cookies, tokens, or user-identifying data —
  only market JSON.
- Credentials are read from `bet_logger/.env` (gitignored), never hardcoded.
- Playwright profiles (`nfl_draft/.cookies/…`, `bookmaker_odds/.bookmaker_profile`,
  `wagerzon_odds/.wagerzon_profile`) are already gitignored.
