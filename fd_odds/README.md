# fd_odds/

FanDuel single-leg MLB odds DuckDB. Written by
`mlb_sgp/scraper_fanduel_singles.py`.

## Tables

- `mlb_odds` — single-leg odds in the offshore 18-column wide schema (same
  as `dk_odds/mlb_odds`). One row per `(game_id, period, market_tier, line)`
  group, where:
  - `period` is `FG` or `F5` (FD does not post F3 or F7 spread/total markets)
  - `market` (the row's market tier) is `main`, `alternate_spreads`, or
    `alternate_totals`
  - Within a row, `away_spread`/`home_spread`/`total` + their prices share
    one quote pair, and `away_ml`/`home_ml` carry moneylines

Each cycle the table is fully rewritten via `BEGIN TRANSACTION` →
`DELETE FROM mlb_odds` → `INSERT` → `COMMIT` (atomic). No history is kept
inside this DB; MLB.R snapshots whatever it needs into `mlb.duckdb`.

## Read pattern (R)

```r
con <- dbConnect(duckdb::duckdb(), "fd_odds/fd.duckdb", read_only = TRUE)
fd <- dbGetQuery(con, "SELECT * FROM mlb_odds")
dbDisconnect(con, shutdown = TRUE)
```

Then feed `fd` into `scraper_to_canonical()` in
`Answer Keys/MLB Answer Key/odds_screen.R`.

## Team names

FD team names already match canonical Odds-API forms (`"Athletics"`,
`"St. Louis Cardinals"`, `"New York Yankees"`, …), so no translation map
is needed — the scraper writes FD names directly. This was verified during
Phase 1 of the singles-scraper plan.

## Market whitelist

The scraper only ingests FD's canonical 10 main/alt markets:

| FD market name                             | period | tier              |
|--------------------------------------------|--------|-------------------|
| `Run Line`                                 | FG     | main              |
| `Total Runs`                               | FG     | main              |
| `Moneyline`                                | FG     | main              |
| `Alternate Run Lines`                      | FG     | alternate_spreads |
| `Alternate Total Runs`                     | FG     | alternate_totals  |
| `First 5 Innings Run Line`                 | F5     | main              |
| `First 5 Innings Total Runs`               | F5     | main              |
| `First 5 Innings Money Line`               | F5     | main              |
| `First 5 Innings Alternate Run Lines`      | F5     | alternate_spreads |
| `First 5 Innings Alternate Total Runs`     | F5     | alternate_totals  |

A whitelist is necessary because FD posts ~150 markets per event and many
of them — `Run Line / Total Runs Parlay`, `Line / Total Parlay 1..12`,
`Moneyline Home Listed`, `Total Runs (Bands)`, `<Team> Total Runs`, etc. —
look like main/alt markets to a naive keyword classifier but are not in
our single-bet schema.

## FD line-encoding quirks

FD's alt-market runners encode the line in the runner NAME instead of on
the `handicap` field. The parser regex-parses those names:

- Alt spread: `"St. Louis Cardinals +1.5"` (handicap=0 → line parsed from name)
- Alt total: `"Over (8.5)"` (handicap=0 → line parsed from parens in name)
- Main spread: `"St. Louis Cardinals"` (handicap=1.5 → line from handicap)
- Main total: `"Over"` (handicap=10.0 → line from handicap)
- Moneyline: team name only, no line at all
