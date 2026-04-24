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
