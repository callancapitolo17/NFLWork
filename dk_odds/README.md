# dk_odds/

DraftKings single-leg MLB odds DuckDB. Written by
`mlb_sgp/scraper_draftkings_singles.py`.

## Tables

- `mlb_odds` — single-leg odds in the offshore 18-column wide schema. One
  row per `(game_id, period, market_tier, line)` group, where:
  - `period` is `FG`, `F5`, or `F7`
  - `market` (the row's market tier) is `main`, `alternate_spreads`, or
    `alternate_totals`
  - Within a row, `away_spread`/`home_spread`/`total` + their prices share
    one quote pair, and `away_ml`/`home_ml` carry moneylines

Each cycle the table is fully rewritten via `BEGIN TRANSACTION` →
`DELETE FROM mlb_odds` → `INSERT` → `COMMIT` (atomic). No history is kept
inside this DB; MLB.R snapshots whatever it needs into `mlb.duckdb`.

## Read pattern (R)

```r
con <- dbConnect(duckdb::duckdb(), "dk_odds/dk.duckdb", read_only = TRUE)
dk <- dbGetQuery(con, "SELECT * FROM mlb_odds")
dbDisconnect(con, shutdown = TRUE)
```

Then feed `dk` into `scraper_to_canonical()` in
`Answer Keys/MLB Answer Key/odds_screen.R`.

## Team-name canonicalization

DK uses abbreviated city prefixes (`"CLE Guardians"`, `"LA Angels"`,
`"STL Cardinals"`, …). The scraper translates these to canonical
Odds-API names (`"Cleveland Guardians"`, etc.) via the `DK_TEAM_MAP` dict
in `scraper_draftkings_singles.py` before writing rows. If DK ever
introduces a new team string, the scraper logs a WARNING listing the
unmapped names so the dict can be updated.
