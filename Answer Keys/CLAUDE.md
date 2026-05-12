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
- F5 (first 5 innings) markets matched via Odds API: h2h, totals, spreads, alternate_totals
- Derivative markets matched via `compare_alts_to_samples` against scraped offshore odds:
  - F3 / F7 spreads + totals + h2h (Wagerzon, Bookmaker for F3 only). Note: Wagerzon currently posts spread+total at F3/F7 but rarely posts moneyline — h2h_1st_*_innings branch is built but typically idle until a book posts F-period MLs.
  - FG alt spreads + alt totals (Wagerzon, Bet105)
  - FG odd/even total runs (Wagerzon — single-game prop, away_ml side = ODD, home_ml side = EVEN per scraper convention)
  - F5 3-way moneyline `h2h_3way_1st_5_innings` (Wagerzon only — `idgmtyp=29` league `lg=1280`). Tie is a real outcome, not push refund. Coexists with the existing 2-way F5 ML matched via `compare_moneylines_to_wagerzon`.
- Team totals (`team_totals_*_fg`, `team_totals_*_h1`) are scraped but NOT yet matched — MLB samples lack per-team score columns. Tracked as gap #2 in the broader matching plan.
- Historical data: `pbp.duckdb/mlb_betting_pbp` (12,719 games)
- Dashboard: port 8083

### Consensus Architecture
- **Historical consensus** (Phase 1): MLB and CBB both use sharp-weighted probit-devigged consensus computed at PBP-build time and stored in `mlb_betting_pbp` / `cbb_betting_pbp`. CBB no longer uses the 0.5 prob hardcode — historical samples now match the live consensus methodology. NFL legacy paths still use the older `Consensus Betting History.R` output.
- **Live consensus** (Phase 2): Sharp books only via `SHARP_BOOKS` in `Tools.R`. Rec books get weight=0. Pinnacle + Bookmaker at 1.1 weight (tiebreaker), LowVig/Circa/Bet105 at 1.0. Games with no sharp coverage are dropped.
- **Devigging method:** All devig in `Tools.R::devig_american` / `devig_american_3way` uses probit (additive z-shift). 2-way uses closed-form `c = -(z1+z2)/2`; 3-way+ uses `uniroot`. See `docs/superpowers/specs/2026-05-11-probit-devig-design.md` for math and rationale.
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
- `mlb.duckdb` — MLB pipeline working tables (`mlb_team_dict`, `mlb_consensus_temp`, `mlb_odds_temp`, `mlb_trifecta_sgp_odds`, scraper raw odds, historical PBP joins). Held under a long write lock during a pipeline run; nothing the dashboard or RFQ bot needs lives here anymore.
- `mlb_mm.duckdb` — MLB MM consumer tables (`mlb_bets_combined`, `mlb_game_samples`, `mlb_samples_meta`, `mlb_sgp_odds`, `mlb_parlay_lines`, `mlb_parlay_opportunities`, `mlb_trifecta_opportunities`). Separate from `mlb.duckdb` so the pipeline's write lock on `mlb.duckdb` never blocks the RFQ bot or the dashboard. `mlb_bets_combined` moved here on 2026-05-05; the others moved 2026-04-30. Mirrors the CBB pattern.
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

### Triple-Play Data Flow
Triple-play props = team scores ≥1 run in the 1st inning AND wins F5 (3-way,
strict lead) AND wins the game.
```
MLB.R Phase 1 (DT load)
  └── determine_inning_1_scoring_vec() → home_scored_in_1st + away_scored_in_1st
      columns on DT (each is 1 if the team scored ≥1 run in inning 1,
      derived from m1 = home_margin_inning_1 and t1 = total_inning_1).
MLB.R Phase 4 (samples export)
  └── home_scored_in_1st + away_scored_in_1st carried into mlb_game_samples
      as 0/1/NA columns (preserved automatically because run_answer_key_sample
      returns a row subset of DT with all columns intact)
wagerzon_odds/scraper_specials.py (NEW)
  ├── Authenticated GET to NewScheduleHelper.aspx?WT=0&lg=4899
  ├── Filters to single-team TRIPLE-PLAY + GRAND-SLAM (cross-game props skipped)
  └── Appends one row per prop to wagerzon_specials in wagerzon.duckdb,
       keyed by scraped_at for snapshot history

parse_legs.R (generic prop parser)
  ├── TOKEN_REGISTRY maps Wagerzon description tokens to structured leg specs:
  │     SCR 1ST → scores_first;  1H → wins_period(F5);  GM → wins_period(FG);
  │     F3/F7 → wins_period(F3/F7);  SCR U<N>/O<N> → team_total_under/over
  ├── parse_legs(description) → list of leg specs (or NULL + warning on unknown)
  ├── eval_leg(leg, samples, side, team_runs, opp_runs) → logical vector
  └── compute_prop_fair(samples, side, legs) → AND-reduce across legs,
        NA rows on home_scored_in_1st / away_scored_in_1st excluded before the mean

mlb_triple_play.R (standalone pricer)
  ├── Reads wagerzon_specials (latest scraped_at) for posted lines
  ├── Reads mlb_game_samples (total_final_score + margin at F3/F5/F7/FG + home/away_scored_in_1st)
  ├── Reads mlb_trifecta_sgp_odds (latest fetch_time) for DK SGP fair odds
  ├── Maps Wagerzon team names → Odds API canonical via WZ_TO_CANONICAL dict
  ├── Joins to consensus_temp for game_id + side (home/away)
  ├── For each row:
  │     parse_legs(description) → compute_prop_fair(...) → model_fair_prob
  │     dk_sgp lookup → devig with DK_SGP_VIG_DEFAULT=1.25 → dk_fair_prob
  │     ↑ populated by mlb_sgp/scraper_draftkings_trifecta.py via system2():
  │       - Map Odds API id → DK event id (last-word team match)
  │       - Resolve each leg via dk_leg_resolvers.LEG_RESOLVERS
  │         (<Team> Run Scored - 1st Inning?, 1st 5 Innings, Moneyline, <Team>: Team Total Runs)
  │       - POST to DK calculateBets with selectionsForYourBet payload
  │       - Write per-prop rows to mlb_trifecta_sgp_odds (atomic snapshot)
  │     blend_dk_with_model(model, dk, vig) → fair_prob (the published fair)
  ├── Prints model_odds, dk_odds, fair_odds (blended), book_odds, edge_pct
  └── Writes mlb_trifecta_opportunities (game_id, hash, kelly_bet, ...) for the dashboard
```
- `home_scored_in_1st` / `away_scored_in_1st` are NA only when inning-1 PBP
  data is missing — extremely rare. Each is independent of the other (both
  teams can score in inning 1). The pricer drops NA rows from the leg's
  AND-reduce inside `compute_prop_fair`.
- F5 3-way: `wins_period` with period=F5 passes only with strict lead (`margin_f5 > 0` for home, `< 0` for away). F5 ties kill the parlay.
- `SCR U<N>` means the listed team's team-total under N (team_total_under leg). Numeric parser handles `"2"`, `"2.5"`, and unicode `"2½"`.
- Helpers: `triple_play_helpers.R` (determine_inning_1_scoring*) is sourced by MLB.R Phase 1; `parse_legs.R` is sourced by mlb_triple_play.R's main block. Both files are pure — no DB or network side effects.
- Adding a new prop type = add a token to `TOKEN_REGISTRY`. No new function needed.
- Adding new MLB teams (none today, but if expansion happens) requires updating `WZ_TO_CANONICAL` in the pricer.
- DK trifecta SGP blend mirrors mlb_correlated_parlay.R but uses `DK_SGP_VIG_DEFAULT = 1.25` (vs. 1.10 for parlays) because trifectas are 3-4 legs and DK shaves harder on additional legs. Only 2 obs/game (home + away) so per-game vig fitting isn't possible; revisit the constant once Plan #2 collects real DK data.
- Blend = mean of available probs (model + DK). When DK is unavailable (table empty, missing leg, scraper failure), blend reduces to model-only.
- DK trifecta SGP odds populated by `mlb_sgp/scraper_draftkings_trifecta.py`,
  invoked from the pricer via `system2()`. Resolvers in
  `mlb_sgp/dk_leg_resolvers.py` map each `parse_legs` leg type to a DK
  selection ID using primitive markets: `<TEAM> Run Scored - 1st Inning?`
  (scores_first → 'Yes' selection), `Moneyline` (FG ML), `1st 5 Innings`
  (F5 ML), `<TEAM>: Team Total Runs` (team totals). The inning-1 Y/N market
  may not be SGP-eligible at all events (it was absent from the 2026-04 SGP
  event_state dump); when missing, the resolver returns None and the blend
  degrades to model-only — same graceful-degrade path as F3/F7 legs.
- F3 / F7 wins_period legs and opp_total_* legs return NULL from the resolver — DK doesn't reliably post these as 2-way primitives. Props containing them get NULL `sgp_decimal` and the R blend correctly degrades to model-only.
- 4-leg GRAND-SLAM SGPs (scores_first + F5 ML + FG ML + team_total_under) are accepted by DK as resolvers but combinability tests have shown DK rejects the 4-leg combination at calculateBets time. Those rows currently get NULL `sgp_decimal` and degrade to model-only fair. TRIPLE-PLAY 3-leg SGPs are accepted cleanly.
- Dashboard tab: `mlb_dashboard.R` reads `mlb_trifecta_opportunities` + `placed_trifectas` and renders a Trifectas tab next to Parlays. Manual-log only — `/api/place-trifecta` and `/api/remove-trifecta` toggle a row in `placed_trifectas` (in `mlb_dashboard.duckdb`). Sizing settings: `trifecta_bankroll`, `trifecta_kelly_mult`, `trifecta_min_edge` rows in `sizing_settings`.

## Common Pitfalls

1. **Odds API commence_time is a string, not POSIXct** — Must parse with `ymd_hms(commence_time, tz = "UTC")` before any time comparison
2. **Scraper auth tokens expire** — Recon scripts refresh them, but check before debugging "0 games found"
3. **Tools.R Rcpp compilation** — Optional but 10-20x faster. Falls back gracefully if Rcpp unavailable.
4. **Sample staleness** — parlay.R and props.R auto-regenerate if >5 min old. Don't cache samples.
5. **18-column schema** — All scrapers MUST write this exact schema. Kalshi is the exception (26 columns with probability fields).
6. **`tol_error` auto-scales to 0.2% of N** — Per Feustel's spec ("+/-1 means you are within 0.2% of your target"), `run_answer_key_sample()` and `generate_all_samples()` default `tol_error = NULL` which computes `max(1L, round(0.002 * current_N))` inside the shrink loop. At N=500 this is 1 (Feustel's canonical value); at N=1272 it's 3. **Never hardcode `tol_error = 1` at the caller** — that silently tightens the rule by 2.5× at MLB's N and causes pathological shrinkage of sparse-region games.

## Known model biases

**MLB correlated parlay — FG `-1.5` fav + over at total ≤ 7 is structurally -EV.**
- Evidence: n=33 placed bets, ROI -73%, bootstrap 95% CI [-97%, -38%] (fully below zero). Confirmed 2026-05-11.
- Mechanism: `mlb_game_samples` overstates the marginal probability of over at low totals. The correlation factor itself (~1.04) is mild — the issue is the underlying over-leg simulator, not the joint structure.
- Mitigation: filter Home Spread + Over and Away Spread + Over combos in `mlb_correlated_parlay.R` when `row$total_line <= 7.0` for FG slice.
- Details: `docs/superpowers/analysis/2026-05-11-fg-fav-leak/findings.md`. Memory: `memory/mlb_parlay_edge_overestimation.md`.

## When Making Changes

- **Adding a new scraper**: Add to `run.py` scraper config, create DuckDB table with 18-column schema, add team name mappings
- **Adding a new sport**: Update `run.py` targets, create sport-specific R script, add to canonical_match.py
- **Adding a new market type**: Update R pricing logic in the sport's answer key, add to dashboard filters
- **Modifying Tools.R**: Changes affect ALL sports. Test CBB, NFL, and MLB at minimum.
- **Adding tables/sections to dashboards**: NEVER use generic selectors like `querySelectorAll(".table-container")[last]` to find a specific table. Always use `getElementById("specific-id")`. Adding a new table shifts positional indices and breaks every function that relied on "grab the last one." This bug has caused silent failures three separate times. Every table container must have a unique ID and all JS must reference it by ID.

## MLB Dashboard — Wagerzon multi-account

The MLB correlated parlay dashboard (port 8083) supports multiple
Wagerzon accounts. Configuration lives in `bet_logger/.env` and is
discovered by `wagerzon_odds/wagerzon_accounts.py`. See
`wagerzon_odds/CLAUDE.md` for the discovery pattern and how to add
an account.

### New endpoints (`mlb_dashboard_server.py`)

- `GET /api/wagerzon/balances` — list all account snapshots.
- `GET /api/wagerzon/last-used` — persisted selector value, or
  `{"label": null}` when no accounts configured.
- `POST /api/wagerzon/last-used` — body `{"label": "<account>"}`,
  returns `{"ok": true}` or `{"success": false, "error": "..."}`.
- `POST /api/place-parlay` — now requires `{"account": "<label>"}`
  in the body in addition to `parlay_hash`. Returns `balance_after`
  snapshot on success.

### Schema

- `placed_parlays.account` — TEXT, label of the WZ account the parlay
  was placed on. NULL for pre-feature rows (treat as primary). Set on
  the breadcrumb INSERT (status='placing') and preserved through the
  finalising upsert.
- `dashboard_settings(key, value, updated_at)` — generic key/value
  preferences. Currently only `wagerzon_last_used`.
