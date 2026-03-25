# Answer Keys — AI Context

## What This Is
The pricing engine for the entire betting operation. Generates fair odds for any market by matching current game lines to historically similar games.

## Critical Architecture

### Pipeline Flow (CBB — primary active sport)
```
run.py (orchestrator)
  ├── [first]    Sharp scrapers (bookmaker, bet105) → their DuckDB files
  ├── [parallel] Other scrapers (wagerzon, hoop88, bfa, kalshi) + CBB.R
  ├── Sentinel: .scrapers_done_cbb signals all scrapers complete
  └── CBB.R reads scraper DBs, generates fair prices, writes pipeline output
```

### Consensus Architecture
- **Historical consensus** (Phase 1): All-book median with 0.5 prob hardcode. Used for sample building. Do not change.
- **Live consensus** (Phase 2): Sharp books only via `SHARP_BOOKS` in `Tools.R`. Rec books get weight=0. Pinnacle + Bookmaker at 1.1 weight (tiebreaker), LowVig/Circa/Bet105 at 1.0. Games with no sharp coverage are dropped.
- **Scraper integration**: `scraper_to_odds_api_format()` converts offshore scraper output to Odds API long format so sharp scrapers (bookmaker, bet105) can participate in consensus.
- Sharp scrapers run before R starts so their DuckDB data is fresh when CBB.R Phase 2 reads it.

### Team Name Resolution (THE #1 source of bugs)
- Python side: `canonical_match.py` — two-layer matching (dict → game-level fallback)
- R side: `resolve_offshore_teams()` in `Tools.R`
- Both must agree on canonical names (Odds API format)
- When adding a new scraper or sport, team name mismatches WILL happen. Test exhaustively.

### DuckDB Databases
- `cbb.duckdb` — Pipeline data (historical odds, betting PBP, dashboard bets)
- `cbb_mm.duckdb` — MM export (predictions, game samples). Separate from `cbb.duckdb` to avoid lock contention with the Kalshi MM bot's 30s polling.
- `pbp.duckdb` — Historical play-by-play data
- `cbb_dashboard.duckdb` — Dashboard state (placed_bets, settings, CLV)
- Never symlink DuckDB files. Always copy if needed in worktrees.

### Race-to-10 Data Flow
```
Acquire CBB Data.R --daily-pbp
  └── extract_first_to_10_h1(pbp) → first_to_10_h1 column in cbb_pbp_v2
build_betting_pbp.R
  └── Propagates first_to_10_h1 into cbb_betting_pbp
CBB.R
  ├── first_to_10_h1 flows into DT via betting_pbp join (no rename needed)
  ├── Exported as 3rd column in cbb_game_samples
  └── race_to_10_h1 predictions added to cbb_raw_predictions
```
- Backfill existing games: `Rscript "Acquire CBB Data.R" --backfill-first10`
- Fair value = `mean(first_to_10_h1)` across resampled historical games

## Common Pitfalls

1. **Odds API commence_time is a string, not POSIXct** — Must parse with `ymd_hms(commence_time, tz = "UTC")` before any time comparison
2. **Scraper auth tokens expire** — Recon scripts refresh them, but check before debugging "0 games found"
3. **Tools.R Rcpp compilation** — Optional but 10-20x faster. Falls back gracefully if Rcpp unavailable.
4. **Sample staleness** — parlay.R and props.R auto-regenerate if >5 min old. Don't cache samples.
5. **18-column schema** — All scrapers MUST write this exact schema. Kalshi is the exception (26 columns with probability fields).

## When Making Changes

- **Adding a new scraper**: Add to `run.py` scraper config, create DuckDB table with 18-column schema, add team name mappings
- **Adding a new sport**: Update `run.py` targets, create sport-specific R script, add to canonical_match.py
- **Adding a new market type**: Update R pricing logic in the sport's answer key, add to dashboard filters
- **Modifying Tools.R**: Changes affect ALL sports. Test CBB and NFL at minimum.
