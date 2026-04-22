# NFL Draft EV Portal

Trader's-cockpit web portal for surfacing +EV NFL Draft bets across
Kalshi and 6 sportsbooks (DraftKings, FanDuel, Bookmaker, Wagerzon, Hoop88, BetOnline).

**Current venue status** (2026-04-22):
- Kalshi, DraftKings, Wagerzon, Hoop88, Bookmaker, BetOnline — posting draft markets, all scraping live.
- FanDuel — draft page temporarily offline as of 2026-04-20 (they pulled the "NFL Draft" tab from `customPageId=nfl`'s layout). Scraper is intact and will pick up automatically when FD reposts; the dashboard's staleness filter hides FD rows while they're gone.

## Architecture

Single DuckDB at `nfl_draft/nfl_draft.duckdb` (legacy `kalshi_draft.duckdb` migrated in and retired).
Python orchestrator `nfl_draft/run.py` invoked by cron.

- `--mode scrape --book all` — pulls odds from all 7 venues (kalshi, draftkings, fanduel, bookmaker, wagerzon, hoop88, betonline), devigs, writes to `draft_odds`, triggers legacy `kalshi_draft/edge_detector.py` + `consensus.py`
- `--mode trades` — polls Kalshi trade tape with cursor + dedup

Dashboard extends `kalshi_draft/app.py` with 4 new tabs under "Portal" section:
Cross-Book Grid, +EV Candidates, Trade Tape, Bet Log.

### Cross-Book Grid — Kalshi vs sportsbook semantics

Each cell is a venue's fair-value estimate:

- Sportsbooks (DK, FD, Bookmaker, Wagerzon, Hoop88): devigged implied probability.
- Kalshi: mid of buy/sell (falls back to last trade when one-sided).

The outlier flag (⚑) fires when a venue's **take price** sits at least
`threshold_pp` away from the cross-venue median:

- Sportsbooks use their devigged probability for the comparison.
- Kalshi uses the actual **buy price** (what you'd pay to take Yes). Because
  Kalshi has no vig, the buy price is close to fair, so offsets from median are
  real +EV signal rather than vig noise.

Hover a Kalshi cell to see the raw buy, sell, and last-trade prices plus the
source Kalshi ticker.

See `docs/superpowers/specs/2026-04-17-nfl-draft-portal-design.md` for the full design.

## Setup

1. Dependencies (should already be installed from the existing kalshi_draft venv):
   ```bash
   /Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/pip install -r nfl_draft/requirements.txt
   /Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/playwright install chromium
   ```

2. Credentials. Each scraper has a different auth model — there is NO single
   repo-root `.env`. If you already have the production `bookmaker_odds` /
   `wagerzon_odds` scrapers working, you are already configured.

   - **DraftKings** (`scrapers/draftkings.py`): public JSON API, no credentials
     required.
   - **FanDuel** (`scrapers/fanduel.py`): public JSON API, no credentials
     required.
   - **Bookmaker** (`scrapers/bookmaker.py` -> `scrapers/recon_bm.py`):
     credentials live in `bet_logger/.env` (shared with the production
     `bookmaker_odds` scraper). Path is resolved by
     `_recon_util.load_env()`. Env vars (see `scrapers/recon_bm.py`):
     ```
     BOOKMAKER_USERNAME=...
     BOOKMAKER_PASSWORD=...
     ```
     On expired cookies, the scraper auto-re-logs-in via the same HTTP
     flow the production scraper uses; only falls back to the headful
     Playwright path (`recon_bm.py --browser`) if HTTP login also fails
     (e.g., Cloudflare JS challenge). League IDs `12273` (PROPOSITIONS)
     and `13425` (ODDS TO WIN) were discovered via XHR interception.
   - **Wagerzon** (`scrapers/wagerzon.py` -> `scrapers/recon_wz.py`):
     credentials live in `bet_logger/.env` (same shared file). Env vars
     (see `scrapers/recon_wz.py`):
     ```
     WAGERZON_USERNAME=...
     WAGERZON_PASSWORD=...
     ```
   - **Kalshi** (`scrapers/kalshi.py`): API key stays in `kalshi_draft/.env`
     (unchanged from the existing `kalshi_draft` setup). Env vars:
     ```
     KALSHI_API_KEY_ID=...
     KALSHI_PRIVATE_KEY_PATH=/path/to/kalshi-private-key.pem
     ```
   - **Hoop88** (`scrapers/hoop88.py`): JWT login via
     `/cloud/api/System/authenticateCustomer`, then
     `/cloud/api/Lines/Get_LeagueLines2` with `sportType=FOOTBALL`,
     `sportSubType='NFL Draft 2026'`, `period='Prop'`, and one call per
     `propDescription` (e.g. `'2nd Overall Pick'`, `'1st Wide Receiver'`).
     The scraper enumerates propDescriptions live via
     `/cloud/api/League/Get_SportsLeagues`, so new markets appear
     automatically. Credentials live in `bet_logger/.env`:
     ```
     HOOP88_URL=https://hoop88.com
     HOOP88_USERNAME=...
     HOOP88_PASSWORD=...
     ```
   - **BetOnline** (`scrapers/betonline.py` -> `scrapers/recon_betonline.py`):
     NFL Draft markets are anonymous-readable, so no sportsbook login is
     required at runtime. The scraper does need a valid Cloudflare
     `cf_clearance` cookie — we reuse the cookie jar at
     `bet_logger/recon_betonline_cookies.json` maintained by the
     `bet_logger/recon_betonline.py` pipeline. Two-step auth at runtime:
     `GET https://www.betonline.ag/get-token` mints an anonymous JWT; we
     use it as `Authorization: Bearer` + seven BetOnline-specific headers
     (`gsetting`, `contests`, `gmt-offset`, `utc-offset`, `actual-time`,
     `iso-time`, `utc-time`) when POSTing per NFL Draft bucket slug to
     `api-offering.betonline.ag/api/offering/Sports/get-contests-by-contest-type2`.
     On HTTP 403 (Cloudflare) or 401 (BetOnline app rejects), refresh:
     ```
     python bet_logger/recon_betonline.py
     ```
     then re-run the scraper. Refresh tokens rotate automatically on each
     authenticated call, so cookies stay warm as long as the bet_logger
     pipeline runs at least daily.

3. Run the one-time migration (preserves historical Kalshi draft odds):
   ```bash
   /Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m nfl_draft.lib.migrate_from_kalshi_draft
   ```

4. Install pre-draft cron:
   ```bash
   crontab nfl_draft/crontab.pre
   ```

5. Keep the laptop awake for the draft window:
   ```bash
   caffeinate -i &
   ```

## Usage

- **Manual scrape**: `python -m nfl_draft.run --mode scrape --book all`
- **Single book (debug)**: `python -m nfl_draft.run --mode scrape --book draftkings`
- **Trade tape**: `python -m nfl_draft.run --mode trades`
- **Dashboard**: `python kalshi_draft/app.py` -> http://127.0.0.1:8090/
  - Override the port with `NFL_DRAFT_DASHBOARD_PORT=9001 python kalshi_draft/app.py`.
  - On launch, the dashboard spawns a detached background scrape so data is populated within ~2 minutes without any manual action. A daemon thread then re-scrapes every 15 minutes while the dashboard is running. Logs to `/tmp/nfl_draft_startup_scrape.log` and `/tmp/nfl_draft_periodic_scrape.log`.
  - A second daemon thread polls the Kalshi trade tape **in-process every 15 seconds** while the dashboard is running, so the Trade Tape tab stays live without any external scheduler. This makes `crontab.pre` / `crontab.draft`'s `--mode trades` line **redundant whenever the dashboard is up** — if you install cron _and_ run the dashboard, both writers will contend for the DuckDB file lock (per-ticker failures are caught, but Kalshi API calls are wasted). Pick one.
  - The "Refresh Data" button in the header kicks an on-demand scrape (non-blocking; returns instantly with a toast).
  - The Cross-Book Grid hides rows older than 2 hours, so venues that go silent (cookies expired, site changes) drop out automatically instead of showing stale data.

## Draft-day mode

Morning of April 23, swap to faster cadence:
```bash
crontab nfl_draft/crontab.draft
```
And toggle the dashboard mode header from "Pre-draft" to "Draft-day" (changes auto-refresh cadence).

## Maintenance

When a scraper produces unmapped players or markets (visible in the dashboard footer):
1. Edit `nfl_draft/config/players.py` or `nfl_draft/config/markets.py`
2. Save - the next cron tick reseeds automatically (or run `python -m nfl_draft.lib.seed` to apply immediately)

## Reconnaissance

If a sportsbook changes its API shape, re-capture the fixture:
```bash
python nfl_draft/scrapers/recon_dk.py
python nfl_draft/scrapers/recon_fd.py
python nfl_draft/scrapers/recon_bm.py --browser   # requires manual login
python nfl_draft/scrapers/recon_wz.py
```
See `nfl_draft/scrapers/RECON_README.md` for troubleshooting.

## Troubleshooting

- **Auth failures** surface in the dashboard footer per book.
- **DuckDB lock errors** should be zero; if observed, the spec describes a JSONL-tail fallback (Phase 2).
- **Unmapped market count growing** - add market_map entries in `config/markets.py`, re-run seed.

## Testing

```bash
/Users/callancapitolo/NFLWork/kalshi_draft/venv/bin/python -m pytest nfl_draft/tests/ -v
```

Three tiers: `unit/` (fast, no I/O), `integration/` (temp DuckDB, fixtures, no network), `live/` (hits real APIs, run manually before merge).
