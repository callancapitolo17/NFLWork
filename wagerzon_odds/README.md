# Wagerzon Odds Scraper

Scrapes live odds from Wagerzon (private offshore book) via REST API with ASP.NET form authentication.

## Method

REST API:
- Auth: ASP.NET form POST (username + password + hidden ViewState fields)
- Odds: `https://backend.wagerzon.com/wager/NewScheduleHelper.aspx`
- Session maintained via `ASP.NET_SessionId` cookie

## Markets Captured

- Main spreads, totals, moneylines (full game + 1H)
- Alt spreads, alt totals (paired, stored separately)
- Team totals (full game + 1H)

## Sports

- NFL, CBB, NBA

## Usage

```bash
python scraper_v2.py nfl
python scraper_v2.py cbb
python scraper_v2.py nba
```

## Auth

Requires in `.env`:
- `WAGERZON_USERNAME`
- `WAGERZON_PASSWORD`

## Files

| File | Purpose |
|------|---------|
| `scraper_v2.py` | Main scraper (current version) |
| `config.py` | Sport configurations (league IDs, API params, table names) |
| `team_mapping.py` | Team name mappings (e.g., "SEA SEAHAWKS" → "Seattle Seahawks") |
| `transform.py` | Transform raw Wagerzon records to standard format |
| `parlay_pricer.py` | Exact parlay pricing via `ConfirmWagerHelper` (see below) |
| `parlay_placer.py` | One-click parlay placement via REST API (see below) |
| `recon_place_parlay.py` | Captures placement HTTP requests from browser (see below) |
| `migrate_placed_parlays.py` | Schema migration for placement tables (see below) |
| `recon_parlay_slip.py` | Playwright recon for the bet-slip flow — captures network traffic to reverse-engineer new endpoints |

## Parlay Pricer

`parlay_pricer.py` calls Wagerzon's `ConfirmWagerHelper.aspx` endpoint to get exact
parlay payouts. It runs in two modes:

### Stage 1 — populate `mlb_parlay_prices`

```bash
python parlay_pricer.py mlb
```

For each MLB game, queries all 4 spread+total parlay combos (FG and F5) and stores
`{wz_decimal, wz_american, wz_win}` in `wagerzon.duckdb/mlb_parlay_prices`. Queries
at `$10000` by default (max decimal precision); falls back to `$100` when a
parlay's max-risk ceiling is below `$10000`.

### Stage 2 — empirical nudge + exact payout per Kelly-sized parlay

```bash
python parlay_pricer.py mlb --exact-payouts
```

Reads `Answer Keys/mlb.duckdb/mlb_parlay_opportunities` (output of
`mlb_correlated_parlay.R`), sweeps stakes `kelly_bet ± NUDGE_RANGE` around each
sized parlay, picks the stake with the best `win/stake` ratio (empirical round-up
nudge — no math assumptions about WZ rounding), and writes `(exact_wager,
exact_to_win)` back onto each row. The dashboard displays those columns directly
so "To Win" matches the WZ slip to the dollar.

### API notes

- `RiskWin: "2"` in the POST body marks the call as a preview quote; WZ skips the
  balance check. Never pass `RiskWin: 0` — it turns every call into a real-bet
  intent and the API returns `BALANCEEXCEED` whenever balance < amount.
- Odds in the `sel` string are **unsigned** (e.g. `_136`, not `_-136`).
- `MINWAGERONLINE` still returns valid `Risk/Win` (just a warning that the amount
  is below WZ's minimum for placement); `MAXPARLAYRISKEXCEED` returns `Risk=0,
  Win=0` (no usable data).

## Storage

DuckDB: `wagerzon.duckdb` → tables: `nfl_odds`, `cbb_odds`, `nba_odds`,
`mlb_odds` (18-column standard schema), `mlb_parlay_prices` (Stage 1 output).

Also: `wagerzon_cbb.duckdb` (CBB-specific historical data).

## Parlay Placer

`parlay_placer.py` is a pure REST client for one-click parlay placement on Wagerzon.

**What it does:**
- Takes a list of `ParlaySpec` objects (legs, amount, expected payout)
- Validates pricing against Wagerzon's `ConfirmWagerHelper` endpoint (preflight drift check)
- Places via `PostWagerMultipleHelper` if prices match
- Returns placement status + ticket number or failure reason

**Public API:**
```python
from parlay_placer import place_parlays, ParlaySpec, PlacementResult

specs = [
    ParlaySpec(legs=[(game_id, "spread", -1.5), ...], wager_amount=100.00, expected_payout=250.00),
    ...
]
results: list[PlacementResult] = place_parlays(specs, dry_run=False)
```

**Status values in `PlacementResult`:**
- `placed` — Bet accepted; ticket number in result
- `price_moved` — Wagerzon's current price ≥ $0.01 off dashboard; aborted (no money at risk)
- `rejected` — Wagerzon refused (insufficient balance, bet too large, line pulled, etc.)
- `auth_error` — Session expired; re-login retry also failed
- `network_error` — Request failed to complete; **AMBIGUOUS** — check Wagerzon's ticket history to confirm
- `orphaned` — Wagerzon confirmed but local DB insert failed; full record in `placement_orphans` table
- `would_place` — Dry run mode; preflight passed (no actual placement)

**Usage:**
```bash
# From repo root; requires WAGERZON_USERNAME + WAGERZON_PASSWORD in .env
python3 wagerzon_odds/parlay_placer.py --help
```

**Tests:** `pytest wagerzon_odds/test_parlay_placer.py -v` (22 test cases)

**Used by:** MLB dashboard's `/api/place-parlay` endpoint

## Recon: Place Parlay

`recon_place_parlay.py` captures HTTP traffic when you click "Place Bet" on a Wagerzon bet slip in your browser.

**What it's for:**
- Discovers endpoint shapes and parameters for new market types or new sportsbooks
- Used during development to reverse-engineer API contracts

**Security note:**
The captured `recon_place_parlay.json` contains your plaintext password in the `confirmPassword` field. The committed recon file has been **scrubbed to `REDACTED`**. If you re-record, **always scrub the file before sharing**.

**How to run:**
1. Start the capture proxy: `python3 wagerzon_odds/recon_place_parlay.py` from repo root
2. Open your browser to `http://localhost:8888`; agree to install the CA certificate
3. Log in to Wagerzon and place a parlay manually
4. Proxy saves the request to `recon_place_parlay.json`

**Storage:** `.gitignore` excludes `recon_*.json` files (not committed).

## Schema Migration

`migrate_placed_parlays.py` ensures `placed_parlays` table has all required columns and creates `placement_orphans` table if missing.

**What it does:**
- Idempotent: safe to re-run without side effects
- Creates/alters `placed_parlays` in `Answer Keys/MLB Dashboard/mlb_dashboard.duckdb`
- Creates `placement_orphans` table for forensics on failed DB writes

**Run once per environment:**
```bash
python3 wagerzon_odds/migrate_placed_parlays.py
```
