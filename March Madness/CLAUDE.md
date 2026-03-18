# March Madness — AI Context

## What This Is
NCAA tournament simulator using power ratings and Monte Carlo methods. Outputs championship probabilities and Survivor Value rankings.

## Key Files
- `shared.R` — Everything shared: team name fixes, rating fetchers, `simulate_game()`, bracket utilities
- `Basketball Model.R` — Pre-tournament: simulates all 63 games × 10,000 runs
- `Dynamic Simulator.R` — Mid-tournament: only simulates remaining games, respects actual results
- `espn_bracket.R` — Fetches bracket structure from ESPN API (not HTML scraping)

## Critical Details

### Team Name Resolution
`TEAM_NAME_FIXES` in `shared.R` is the single source of truth. All 4 rating sources (BPI, KenPom, Torvik, EvanMiya) map through this. When a name doesn't resolve, add it here — NOT in individual fetcher functions.

### Bracket Structure
ESPN API returns seeding order. `get_bracket_matchups()` enforces correct bracket structure (1v16, 8v9, etc.) and region assignment.

### Simulation Method
Residual sampling: `margin = (rating_A - rating_B) + random_residual`. Uses actual historical variance. This is stochastic — run enough simulations (10k minimum).

## When Making Changes
- Adding a rating source: add fetcher in `shared.R`, merge into composite in `fetch_power_ratings()`
- Team name issues: always fix in `TEAM_NAME_FIXES`, never in individual fetchers
- Mid-tournament bugs: check that eliminated teams are properly filtered in Dynamic Simulator
