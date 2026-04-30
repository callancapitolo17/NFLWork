# WagerzonJ — Second Wagerzon Account in Bet Logger

**Date:** 2026-04-29
**Branch:** `feature/wagerzon-second-account`
**Worktree:** `~/NFLWork/.worktrees/wagerzon-second-account`

## Goal

Add a second Wagerzon account ("WagerzonJ") to the bet_logger pipeline. The user has an 87.5% share of the risk on this account (a partner covers the remaining 12.5%). Sheet1 should record the user's adjusted share so existing summary tabs compute correct P&L automatically; raw (unadjusted) bets go to the existing `Shared` tab for week-over-week balance reconciliation.

## Approach

Mirror the existing BFA two-account pattern (`scraper_bfa.py:38` `ACCOUNTS` dict + `--account` flag) inside `scraper_wagerzon.py`. Add a parallel invocation in `run_all_scrapers.sh`. Reuse the existing `Shared` Google Sheets tab (already used by BFAJ) for raw-bet verification.

The only novel element vs. BFA's pattern is the adjustment type: BFA subtracts a flat `$15` per bet (additive). Wagerzon needs `× 0.875` (multiplicative). The new account config field is `bet_multiplier` rather than `bet_adjustment`.

## Code changes

### `bet_logger/scraper_wagerzon.py`

1. Add a module-level `ACCOUNTS` dict:
   ```python
   ACCOUNTS = {
       'default': {
           'username_env': 'WAGERZON_USERNAME',
           'password_env': 'WAGERZON_PASSWORD',
           'platform': 'Wagerzon',
           'bet_multiplier': 1.0,
           'shared_sheet': None,
       },
       'j': {
           'username_env': 'WAGERZONJ_USERNAME',
           'password_env': 'WAGERZONJ_PASSWORD',
           'platform': 'WagerzonJ',
           'bet_multiplier': 0.875,
           'shared_sheet': 'Shared',
       },
   }
   ```

2. Refactor `login(session)` → `login(session, username, password)` so credentials are passed in instead of read from module-level globals.

3. Refactor `scrape_wagerzon(weeks_back, all_weeks)` → `scrape_wagerzon(weeks_back, all_weeks, account_name='default')`. Inside, look up the account config and pass per-account credentials and `platform` / `bet_multiplier` down to parsing.

4. In `parse_api_bets(history, platform='Wagerzon', bet_multiplier=1.0)`:
   - Compute `adjusted_risk = round(risk * bet_multiplier, 2)`
   - Set `bet['platform'] = platform`
   - Set `bet['bet_amount'] = adjusted_risk`
   - Add `bet['_raw_risk'] = risk` (matches BFA's field name, so the `Shared` upload code uses identical structure)
   - **Do not** scale `american_odds`, `decimal_odds`, or `result` — these are ratios/labels independent of stake size. P&L automatically scales because `0.875 × risk × (decimal_odds − 1) = your share of winnings`.

5. CLI: add `--account` argument with `choices=list(ACCOUNTS.keys())`, default `'default'`.

6. `__main__` block: mirror `scraper_bfa.py:498-523` structure:
   - Always upload adjusted bets to `Sheet1`
   - If `acct['shared_sheet']` is set, also upload raw bets (using `_raw_risk` as `bet_amount`) to that tab.

### `bet_logger/run_all_scrapers.sh`

After the existing Wagerzon block (current lines 19–28), insert:

```bash
# Run Wagerzon scraper — WagerzonJ account
echo "[$(date '+%H:%M:%S')] Running Wagerzon scraper (WagerzonJ)..."
echo "----------------------------------------"
if ./venv/bin/python3 scraper_wagerzon.py --account j; then
    echo "[$(date '+%H:%M:%S')] WagerzonJ: done"
else
    echo "[$(date '+%H:%M:%S')] WagerzonJ: FAILED (exit $?)"
    FAILED=$((FAILED + 1))
    FAILED_NAMES="${FAILED_NAMES}WagerzonJ, "
fi
echo ""
```

Update the success message string from `"All 5 scrapers completed successfully."` to `"All 6 scrapers completed successfully."`.

### `bet_logger/.env.example`

Add the new variables (commented placeholder, no real values):

```
# Wagerzon second account (WagerzonJ — 87.5% of risk attributable to user)
WAGERZONJ_USERNAME=
WAGERZONJ_PASSWORD=
```

The user populates the real `.env` themselves; `.env` is gitignored.

## Documentation changes

### `bet_logger/README.md`

- Top-line description: `"Wagerzon (2 accounts), Hoop88, BFA Gaming (2 accounts), and BetOnline"`.
- `.env` setup: add `WAGERZONJ_USERNAME` / `WAGERZONJ_PASSWORD` to the variables list.
- "Run individual scrapers" usage: add `[--account j]` flag to the Wagerzon line.
- "Run all scrapers": update count from "Runs Wagerzon, Hoop88, BFA (×2), BetOnline" to "Runs Wagerzon (×2), Hoop88, BFA (×2), BetOnline".
- "Output Columns" → Platform example list: append `WagerzonJ`.
- Add a "Wagerzon — two accounts" subsection in Platform Notes explaining:
  - WagerzonJ is the second account, 87.5% of risk attributed to user
  - Adjusted bets go to Sheet1, raw bets go to `Shared` tab
  - Identical adjustment pattern to BFAJ but multiplicative (×0.875) instead of additive (−$15)

### `bet_logger/CLAUDE.md`

Add a brief note in the "Architecture" section:
> Both BFA and Wagerzon support two accounts (`--account j`). BFA's BFAJ adjusts by flat $-15/bet; Wagerzon's WagerzonJ adjusts by ×0.875. Adjusted bets → Sheet1; raw bets → `Shared` tab.

## Version control plan

- **Worktree:** `~/NFLWork/.worktrees/wagerzon-second-account` on branch `feature/wagerzon-second-account`. Created before any spec/code work.
- **Files modified:**
  - `bet_logger/scraper_wagerzon.py`
  - `bet_logger/run_all_scrapers.sh`
  - `bet_logger/.env.example`
  - `bet_logger/README.md`
  - `bet_logger/CLAUDE.md`
  - `docs/superpowers/specs/2026-04-29-wagerzon-second-account-design.md` (this file)
- **Files NOT touched:**
  - `bet_logger/.env` — user adds real credentials locally; never committed
  - `bet_logger/sheets.py`, `bet_logger/utils.py` — no changes needed; existing `append_bets_to_sheet(bets, sheet_name=...)` already supports tab targeting
  - `bet_logger/scraper_bfa.py` — kept untouched per Approach 1 (minimal blast radius)
- **Commits:** single commit covering code + docs together (per CLAUDE.md "documentation in same commit as feature").
- **Pre-merge testing (from inside worktree):**
  1. `./venv/bin/python3 scraper_wagerzon.py --account j --dry-run` — verify login, parsing, multiplier (×0.875), `_raw_risk` preservation
  2. `./venv/bin/python3 scraper_wagerzon.py --account default --dry-run` — regression check that primary account still works
  3. Confirm dry-run output shows correct adjusted vs. raw amounts
  4. (Optional, with user OK) live run with `--account j` and check Sheet1 + Shared tab end-to-end
- **Pre-merge review:** full `git diff main..HEAD` against the executive engineer review checklist (data integrity, edge cases, dead code, no secrets in logs, etc.). I will surface findings and get explicit approval before merging.
- **Merge:** only after explicit user approval.
- **Post-merge cleanup:** `git worktree remove .worktrees/wagerzon-second-account` + `git branch -d feature/wagerzon-second-account`.

## Edge cases & risks

1. **Floating-point rounding** — `round(risk * 0.875, 2)` ensures sheet-friendly amounts (no `$87.49999999`).
2. **Dedup safety** — existing key is `(date, platform, description, amount)`. Distinct platform label (`WagerzonJ` vs `Wagerzon`) prevents collision with primary.
3. **Historical data not retro-relabeled** — bets previously placed on the second account but logged under `Wagerzon` will not be auto-corrected. Out of scope for this change. User can manually relabel if needed.
4. **Summary tabs** — MLB Summary and CBB Summary filter by sport, not platform. WagerzonJ bets in those sports will appear in summaries at the user's adjusted (87.5%) share, which is the correct tracking of user's actual P&L exposure.
5. **First scrape** — `--account j` reads the same `--weeks 1` default as primary. No `--since-last` complications because the existing Wagerzon scraper doesn't use `--since-last`.
6. **Auth failure on either account does not block the other** — `run_all_scrapers.sh` already runs each block with `if … then … else FAILED++`, continuing through failures.
7. **No new credentials files** — Wagerzon uses pure HTTP login per request, no token caching, so no `recon_*` files are needed (unlike BFA's Keycloak refresh tokens).

## Out of scope

- Refactoring BFA and Wagerzon to share a generic multi-account helper module (Approach 2 considered and rejected as YAGNI).
- Backfilling historical WagerzonJ bets that were logged as `Wagerzon`.
- Adjusting MLB Summary / CBB Summary tab filters to exclude or separately report WagerzonJ.
- Any changes to BFA's two-account handling.
