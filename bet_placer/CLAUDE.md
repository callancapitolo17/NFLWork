# Bet Placer — AI Context

## What This Is
Automated bet placement via Playwright browser automation. Launched from the dashboard, opens a visible browser, navigates to the sportsbook, and pre-fills the betslip. User must manually confirm.

## Architecture
- `placer.py` — CLI entry point, dispatches to the correct navigator
- `base_navigator.py` — Shared: game lookup (DuckDB), market parsing, status updates, profile cleanup
- `navigator_*.py` — One per sportsbook with book-specific DOM selectors and login flows

## Critical Patterns

### Status Flow
```
"navigating" → "ready_to_confirm" → "pending" (user confirmed in dashboard)
           → "nav_error" (game/odds not found)
           → "nav_timeout" (browser closed without confirmation)
```
Status lives in `cbb_dashboard.duckdb` table `placed_bets`. Updates use `only_if` guards to prevent race conditions.

### Game Lookup
Each navigator queries its book's DuckDB to find the game by canonical team names. If the book uses non-canonical names, the navigator does fuzzy matching (first-word, substring, abbreviation expansion).

### DOM Scoring
Navigators score potential odds buttons before clicking:
- Exact text match = +5, partial = +1-2
- Position match = +2 bonus
- Minimum score 3 to click (prevents wrong-bet placement)

## Common Pitfalls

1. **Singleton Lock**: Orphaned Chrome processes leave `SingletonLock` files. `cleanup_singleton_lock()` handles this.
2. **SPA re-renders**: Hoop88's SPA clears betslip values when adding new bets. Solution: add all odds first, then fill amounts.
3. **Cloudflare**: BetOnline requires real Chrome (not Playwright Chromium). Uses persistent profile.
4. **Worktree paths**: `base_navigator.py` detects worktree context and resolves DuckDB paths to main repo.

## When Making Changes
- Test with a single small bet first
- Never run headless (user must see and confirm)
- DOM selectors change frequently — recon scripts capture current structure
