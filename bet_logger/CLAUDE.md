# Bet Logger — AI Context

## What This Is
Scrapes bet history from 4 sportsbooks and logs them to Google Sheets for P&L tracking.

## Architecture
- `scraper_*.py` — One per sportsbook (Wagerzon, Hoop88, BFA, BetOnline)
- `sheets.py` — Google Sheets API integration (append rows, read existing)
- `utils.py` — Shared utilities (date parsing, deduplication)
- `create_summary.py` — Generate summary reports
- `run_all_scrapers.sh` — Run all scrapers sequentially

## Multi-account scrapers

`scraper_bfa.py` supports a second account via `--account j`. `scraper_wagerzon.py`
still defines `--account j` / `--account c` entries, but as of 2026-06-26 only the
primary Wagerzon login is active — see "Wagerzon account consolidation" below.

- **BFAJ** (BFA second account) — flat `bet_adjustment: -15` (subtract $15 per bet). Uploads adjusted to Sheet1 and raw to `Shared` tab.
- **WagerzonJ** (Wagerzon second account) — multiplicative `bet_multiplier: 0.875` (user holds 87.5% of risk). Uploads adjusted to Sheet1 and raw to `Shared` tab. **Inactive: login lost 2026-06-26 (`WAGERZONJ_*` env vars unset).**
- **WagerzonC** (Wagerzon third account) — `bet_multiplier: 1.0` (full risk to user, no adjustment). Uploads adjusted to Sheet1 only; no Shared tab. **This is now the only live Wagerzon account; its credentials moved to the primary `WAGERZON_*` slot (see below).**

The two scrapers each define their own `ACCOUNTS` dict; they're not yet sharing a generic helper (intentional — see spec `docs/superpowers/specs/2026-04-29-wagerzon-second-account-design.md`, Approaches 1 vs 2).

### Wagerzon account consolidation (2026-06-26)

The original primary Wagerzon login and WagerzonJ were lost. The surviving
account (formerly **WagerzonC**) was promoted into the primary
`WAGERZON_USERNAME` / `WAGERZON_PASSWORD` slot in `.env`. To keep its P&L under
the existing **WagerzonC** spreadsheet column, `ACCOUNTS['default']` now carries
`platform: 'WagerzonC'` while still reading the primary env slot, and
`run_all_scrapers.sh` runs a single Wagerzon scrape (no `--account` flag). The
`j` and `c` `ACCOUNTS` entries are retained for reference but are inert — their
`WAGERZONJ_*` / `WAGERZONC_*` env vars are unset, so the cron no longer invokes
them. Note the odds/placement side (`wagerzon_odds/`, dashboard registry) still
labels this single login **"Wagerzon"** internally; only the bet-history sheet
calls it "WagerzonC".

## Critical Details

### Credentials
All credentials live in `.env` in this directory. This is the SHARED `.env` — odds scrapers and bet_placer also read from here.

### BetOnline EndDate Bug
The BetOnline API returns 0 results if EndDate is "tomorrow" in BetOnline's timezone. Use `datetime.now()` (local time), NOT `datetime.now(timezone.utc)`.

### Google Sheets
- `credentials.json` — Service account key (gitignored)
- Sheets are append-only; deduplication happens in the scraper before writing
- All Sheets API calls go through `_execute_with_retry()` in `sheets.py`, which
  retries transient DNS/connection errors and 429/5xx responses (4 attempts,
  linear backoff). This guards the 5 AM cron run against network blips — a DNS
  failure used to crash uncaught mid-upload and silently drop a scraper's bets.
  Genuine outages still fail loudly (non-zero exit → FAILED notification).

### Recon Scripts
`recon_*.py` and `recon_*.json` files capture API endpoints and auth tokens from browser sessions. Run these when a scraper breaks due to expired tokens.

## When Making Changes
- Test with `--dry-run` flag if available
- Check deduplication logic before running — duplicate rows in Sheets are hard to clean up
- Never commit `.env` or `credentials.json`
