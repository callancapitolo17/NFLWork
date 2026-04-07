# Answer Keys — AI Context

## What This Is
The pricing engine for the entire betting operation. Generates fair odds for any market by matching current game lines to historically similar games.

## Critical Architecture

### Pipeline Flow (CBB — primary active sport)
```
run.py cbb (orchestrator)
  ├── [first]    Sharp scrapers (bookmaker, bet105) → their DuckDB files
  ├── [parallel] Other scrapers (wagerzon, hoop88, bfa, kalshi) + CBB.R
  ├── Sentinel: .scrapers_done_cbb signals all scrapers complete
  └── CBB.R reads scraper DBs, generates fair prices, writes pipeline output
```

### Pipeline Flow (MLB)
```
run.py mlb (orchestrator)
  ├── [first]    Sharp scrapers (bookmaker, bet105) → their DuckDB files
  ├── [parallel] Other scrapers (wagerzon, hoop88, bfa) + MLB.R
  ├── Sentinel: .scrapers_done_mlb signals all scrapers complete
  └── MLB.R reads scraper DBs, generates F5 fair prices, writes pipeline output
```
- MLB uses moneyline-based matching (`use_spread_line = FALSE`)
- F5 (first 5 innings) markets only: h2h, totals, spreads
- Historical data: `pbp.duckdb/mlb_betting_pbp` (12,719 games)
- Dashboard: port 8083

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
- `cbb.duckdb` — CBB pipeline data (historical odds, betting PBP, dashboard bets)
- `cbb_mm.duckdb` — CBB MM export (predictions, game samples). Separate to avoid lock contention.
- `mlb.duckdb` — MLB pipeline output (bets_combined, team_dict)
- `mlb_dashboard.duckdb` — MLB dashboard state (placed_bets, settings, CLV)
- `pbp.duckdb` — Historical play-by-play data (shared: MLB + others)
- `cbb_dashboard.duckdb` — CBB dashboard state (placed_bets, settings, CLV)
- Never symlink DuckDB files. Always copy if needed in worktrees.

### Race-to-X Data Flow
All race-to-X props are **full game** — "which team reaches X points first" at any point.
```
Acquire CBB Data.R --daily-pbp
  └── extract_race_to_fg(pbp, threshold) → first_to_{10,20,40}_fg in cbb_pbp_v2
build_betting_pbp.R
  └── Propagates first_to_{10,20,40}_fg into cbb_betting_pbp
CBB.R
  ├── first_to_X_fg flows into DT via betting_pbp join
  ├── Exported as columns 3-5 in cbb_game_samples
  └── race_to_{10,20,40} predictions added to cbb_raw_predictions
```
- Backfill: `Rscript "Acquire CBB Data.R" --backfill-race-fg <threshold>` (10, 20, or 40)
- Fair value = `mean(first_to_X_fg)` across resampled historical games
- Backward compat: `first_to_10_h1` still extracted alongside FG columns during transition

## Common Pitfalls

1. **Odds API commence_time is a string, not POSIXct** — Must parse with `ymd_hms(commence_time, tz = "UTC")` before any time comparison
2. **Scraper auth tokens expire** — Recon scripts refresh them, but check before debugging "0 games found"
3. **Tools.R Rcpp compilation** — Optional but 10-20x faster. Falls back gracefully if Rcpp unavailable.
4. **Sample staleness** — parlay.R and props.R auto-regenerate if >5 min old. Don't cache samples.
5. **18-column schema** — All scrapers MUST write this exact schema. Kalshi is the exception (26 columns with probability fields).
6. **`tol_error` auto-scales to 0.2% of N** — Per Feustel's spec ("+/-1 means you are within 0.2% of your target"), `run_answer_key_sample()` and `generate_all_samples()` default `tol_error = NULL` which computes `max(1L, round(0.002 * current_N))` inside the shrink loop. At N=500 this is 1 (Feustel's canonical value); at N=1272 it's 3. **Never hardcode `tol_error = 1` at the caller** — that silently tightens the rule by 2.5× at MLB's N and causes pathological shrinkage of sparse-region games.

## When Making Changes

- **Adding a new scraper**: Add to `run.py` scraper config, create DuckDB table with 18-column schema, add team name mappings
- **Adding a new sport**: Update `run.py` targets, create sport-specific R script, add to canonical_match.py
- **Adding a new market type**: Update R pricing logic in the sport's answer key, add to dashboard filters
- **Modifying Tools.R**: Changes affect ALL sports. Test CBB, NFL, and MLB at minimum.
- **Adding tables/sections to dashboards**: NEVER use generic selectors like `querySelectorAll(".table-container")[last]` to find a specific table. Always use `getElementById("specific-id")`. Adding a new table shifts positional indices and breaks every function that relied on "grab the last one." This bug has caused silent failures three separate times. Every table container must have a unique ID and all JS must reference it by ID.
