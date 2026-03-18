# Bet Logger — AI Context

## What This Is
Scrapes bet history from 4 sportsbooks and logs them to Google Sheets for P&L tracking.

## Architecture
- `scraper_*.py` — One per sportsbook (Wagerzon, Hoop88, BFA, BetOnline)
- `sheets.py` — Google Sheets API integration (append rows, read existing)
- `utils.py` — Shared utilities (date parsing, deduplication)
- `create_summary.py` — Generate summary reports
- `run_all_scrapers.sh` — Run all scrapers sequentially

## Critical Details

### Credentials
All credentials live in `.env` in this directory. This is the SHARED `.env` — odds scrapers and bet_placer also read from here.

### BetOnline EndDate Bug
The BetOnline API returns 0 results if EndDate is "tomorrow" in BetOnline's timezone. Use `datetime.now()` (local time), NOT `datetime.now(timezone.utc)`.

### Google Sheets
- `credentials.json` — Service account key (gitignored)
- Sheets are append-only; deduplication happens in the scraper before writing

### Recon Scripts
`recon_*.py` and `recon_*.json` files capture API endpoints and auth tokens from browser sessions. Run these when a scraper breaks due to expired tokens.

## When Making Changes
- Test with `--dry-run` flag if available
- Check deduplication logic before running — duplicate rows in Sheets are hard to clean up
- Never commit `.env` or `credentials.json`
