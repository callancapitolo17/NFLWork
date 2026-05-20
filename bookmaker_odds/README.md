# Bookmaker.eu Odds Scraper

Scrapes live odds from Bookmaker.eu using curl_cffi for Cloudflare bypass.

## Method

HTTP POST to `https://be.bookmaker.eu/gateway/BetslipProxy.aspx` with Chrome
TLS fingerprint impersonation via curl_cffi.

## Sports + markets

| Sport | Periods captured | BKM leagueDescEn labels |
|---|---|---|
| MLB | full game, 1st 5 innings, 2nd halves | `GAME LINES`, `1ST 5 INNINGS`, `2ND HALVES` |
| NBA | full game, 1st halves | `GAME LINES`, `1ST HALVES` |
| CBB | full game, 1st halves | `GAME LINES`, `1ST HALVES` (off-season at last verify) |

Main spreads, totals, and moneylines per market — no alt lines, no team totals.

## How league IDs are resolved (no more hardcoded magic numbers)

League IDs are **not hardcoded** in `SPORT_CONFIGS`. Instead, each market is
declared by its human-readable BKM name (`leagueDescEn`) plus the
`(sportId, region)` qualifier. At scrape time, the scraper:

1. Calls `GetRoutingInfo` to fetch BKM's live league catalog.
2. Calls `resolve_leagues()` to look each wanted market up by
   `(sportId, region, leagueDescEn)`.
3. Iterates the resolved `leagueId`s and fetches their schedules.

**Why:** an earlier version hardcoded `leagueId=503` as "1st 3 Innings", but
BKM actually labels that league as `2ND HALVES`. The scraper happily fetched
2nd-half data and tagged it as F3 — silent corruption that polluted the MLB
dashboard for days. Dynamic discovery turns that failure mode into a
**visible break**: if BKM renames a league or the config has a typo, the
scraper logs `WARNING: no BKM league matching leagueDescEn='X'` and lists
the patterns that ARE available at that scope — actionable diagnostic, not
silent fetching of the wrong endpoint.

**Disambiguation matters:** BKM groups some leagues under shared `sportId`s
(e.g., both NBA and WNBA sit under `sportId=NBA`, distinguished only by
`region`). The lookup key is always the triple
`(sportId, region, leagueDescEn)`. Pin all three in `SPORT_CONFIGS`.

### Adding a new market

1. Run a fresh `GetRoutingInfo` (or inspect `recon_bookmaker_api.json`) to
   find the league.
2. Add an entry to `SPORT_CONFIGS[<sport>]["markets"]`:
   ```python
   {"league_pattern": "<BKM leagueDescEn>", "market": "<internal_market>",
    "period": "<internal_period>"},
   ```
3. The scraper will resolve the live `leagueId` automatically on next run.

## Usage

```bash
python scraper.py mlb
python scraper.py nba
python scraper.py cbb
```

## Auth

Requires in `.env`:
- `BOOKMAKER_USERNAME`
- `BOOKMAKER_PASSWORD`

Cookies cached in `.bookmaker_cookies.json`. If Cloudflare blocks the request
AND stdin is a TTY AND the last recon attempt was > 1h ago, the scraper
launches `recon_bookmaker.py` in a real Chrome window to refresh
`cf_clearance`.

Otherwise (piped subprocess from `run.py`, recent recon attempt, or session
healthy but BKM returned no leagues for this sport) the scraper clears stale
data, logs the reason, and exits 0 without opening a browser. To refresh
cookies manually when the pipeline has skipped it:

```bash
cd bookmaker_odds
./venv/bin/python recon_bookmaker.py
```

## Storage

DuckDB: `bookmaker.duckdb` → tables: `mlb_odds`, `nba_odds`, `cbb_odds`
(18-column standard schema).

## Tests

Pure-Python helper tests in `tests/` cover the league-resolution contract
(NBA/WNBA disambiguation, missing-market warnings, empty-routing defenses)
and the existing recon/auth helpers. No network, no DuckDB:

```bash
cd bookmaker_odds && ./venv/bin/python -m pytest tests/ -v
```
