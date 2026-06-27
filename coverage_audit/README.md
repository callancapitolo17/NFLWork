# MLB Scraper Coverage Audit

Daily check that every MLB book you scrape is still producing the markets it
used to, is fresh, and (for the 5 pill-rendered books) reaches the odds screen.
Detection is a deterministic, read-only Python script; fixing is done by a daily
Claude Desktop scheduled task following `AGENT_PLAYBOOK.md`.

## What it checks (deterministic, network-free)
- **Regression** — a `(market, period)` the book posted in ≥80% of the last 7
  days of runs but not in today's run (parse break / auth expiry). Skipped for
  tables with no `market`/`period` columns (e.g. `wagerzon_specials`).
- **Freshness** — latest run older than 26h, or DB unreadable/empty. Naive
  timestamps (e.g. `wagerzon_specials.scraped_at`) are treated as UTC.
- **Row-count** — today's run below 50% of the trailing median (partial break).
- **Screen presence** — an `expected_on_screen` book with no rows in the latest
  ~26h of `mlb_bets_book_prices`. (Window-based, because each book's rows carry
  their own scraper `fetch_time` — there is no single shared run timestamp.)

Fuzzy scraped-vs-screen label matching and live "never-seen market" discovery
are intentionally the agent's job at runtime — see `AGENT_PLAYBOOK.md`.

## The 9 books
`wagerzon, wagerzon_specials, hoop88, bfa, bookmaker, bet105, kalshi,
draftkings_singles, fanduel_singles`. Only 5 render on the odds screen
(`wagerzon, bookmaker, bet105, draftkings, fanduel`); the rest feed pricing /
consensus and are not expected on the screen (`expected_on_screen=False`).

## Run it manually
```
python -m coverage_audit.audit            # writes a row, notifies on NEW gaps, prints JSON
python -m coverage_audit.audit --no-notify --json-only
```
Run from the repo root so the package resolves and the per-book DuckDB paths
(derived from the repo root) point at the real data. Output table:
`coverage_audit/coverage.duckdb::coverage_gaps` (one row per gap per run).

## Register the daily Desktop scheduled task
This runs locally on your Mac (NOT a `/schedule` cloud routine — that runs on
Anthropic infra with only a git clone, so it can't reach your DBs / auth / IP).
In the **Claude Code Desktop app**:
1. Sidebar → **Routines** → **New routine** → choose **Local** (not Remote).
2. Working folder: the repo root.
3. Schedule: daily ~9:00 AM Pacific.
4. Instructions: "Follow `coverage_audit/AGENT_PLAYBOOK.md`."

Caveat: the Desktop app must be running for the task to fire (it does not need
an open chat session). If you ever want detection to run even with the app
closed, schedule `python -m coverage_audit.audit` alone via launchd and keep
only the wiring as a Desktop task.

## Add a book or change thresholds
- New book → add a `Book(...)` row to `coverage_audit/registry.py`. If its table
  uses a different timestamp column, set `ts_col=`. If it has no `market`/
  `period` columns, detection automatically falls back to freshness + rowcount.
- Thresholds are constants/arguments: `BASELINE_DAYS` (audit.py), and
  `presence_threshold` / `max_age_hours` / `min_ratio` / `within_hours` in
  `detectors.py`.

## Tests
```
python -m pytest coverage_audit/ -v
```
Unit tests use synthetic temp DuckDBs (fast, isolated). A real read-only run
against the live book DBs is the integration smoke test (`python -m
coverage_audit.audit --no-notify`).
