# WagerzonC — Third Wagerzon Account in Bet Logger

**Date:** 2026-04-30
**Branch:** `feature/wagerzon-third-account`
**Worktree:** `~/NFLWork/.worktrees/wagerzon-third-account`
**Builds on:** `docs/superpowers/specs/2026-04-29-wagerzon-second-account-design.md`

## Goal

Add a third Wagerzon account (`WagerzonC`) to the bet_logger pipeline. Full risk attributed to user (no multiplier, no partner share). Behaves exactly like the primary account except for separate credentials and a distinct platform label so bets can be tracked independently in `Sheet1`.

## Approach

Same pattern as WagerzonJ, simpler config. Add one new entry to the `ACCOUNTS` dict with `bet_multiplier: 1.0` and `shared_sheet: None`. The multi-account scaffolding shipped with WagerzonJ already supports any number of `ACCOUNTS` entries — `--account` argparse choices, login/scrape parameterization, and the dual-tab upload guard are all gated on per-account config. No code changes outside the dict.

## Code changes

### `bet_logger/scraper_wagerzon.py`

Add a third entry to the `ACCOUNTS` dict (after the existing `'j'` entry):

```python
'c': {
    'username_env': 'WAGERZONC_USERNAME',
    'password_env': 'WAGERZONC_PASSWORD',
    'platform': 'WagerzonC',
    'bet_multiplier': 1.0,
    'shared_sheet': None,
},
```

No other code edits. The existing argparse line `parser.add_argument('--account', choices=list(ACCOUNTS.keys()), default='default', ...)` automatically picks up `'c'`. The `__main__` block's `if acct['shared_sheet']:` guard skips the Shared-tab upload (since `None` is falsy). The `parse_api_bets` multiplier branch is a no-op when `bet_multiplier == 1.0`.

### `bet_logger/test_scraper_wagerzon.py`

Two test changes:

1. Update `test_accounts_has_default_and_j` to also include `'c'`:
   ```python
   def test_accounts_keys():
       assert set(ACCOUNTS.keys()) == {'default', 'j', 'c'}
   ```
   (Rename and update the existing test rather than adding a parallel one.)

2. Add `test_accounts_c_fields` mirroring the existing `test_accounts_default_fields`:
   ```python
   def test_accounts_c_fields():
       acct = ACCOUNTS['c']
       assert acct['username_env'] == 'WAGERZONC_USERNAME'
       assert acct['password_env'] == 'WAGERZONC_PASSWORD'
       assert acct['platform'] == 'WagerzonC'
       assert acct['bet_multiplier'] == 1.0
       assert acct['shared_sheet'] is None
   ```

The existing 5 multiplier tests (default + j cases) continue to cover the math; no new multiplier tests needed since WagerzonC has the same `bet_multiplier: 1.0` as the primary.

## Env (`bet_logger/.env.example`)

Append two new placeholder lines after the existing `WAGERZONJ_PASSWORD=`:

```
# Wagerzon third account (WagerzonC — full risk to user, no adjustment)
WAGERZONC_USERNAME=
WAGERZONC_PASSWORD=
```

The user populates real values in `.env` separately. `.env` is gitignored.

## Scheduler (`bet_logger/run_all_scrapers.sh`)

Insert a third Wagerzon block immediately after the existing WagerzonJ block (right before the Hoop88 block):

```bash
# Run Wagerzon scraper — WagerzonC account
echo "[$(date '+%H:%M:%S')] Running Wagerzon scraper (WagerzonC)..."
echo "----------------------------------------"
if ./venv/bin/python3 scraper_wagerzon.py --account c; then
    echo "[$(date '+%H:%M:%S')] WagerzonC: done"
else
    echo "[$(date '+%H:%M:%S')] WagerzonC: FAILED (exit $?)"
    FAILED=$((FAILED + 1))
    FAILED_NAMES="${FAILED_NAMES}WagerzonC, "
fi
echo ""
```

Update the success message: `"All 6 scrapers completed successfully."` → `"All 7 scrapers completed successfully."`.

## Documentation

### `bet_logger/README.md`

- Top-line description: `Wagerzon (2 accounts)` → `Wagerzon (3 accounts)`.
- `.env` setup: add `WagerzonC` to the credentials list.
- Run-all description: `(primary + WagerzonJ)` → `(primary + WagerzonJ + WagerzonC)`.
- Wagerzon flags subsection: extend `--account` help to mention `c`.
- Output Columns Platform example: append `WagerzonC`.
- Wagerzon Platform Notes: add a third bullet under the "Two accounts supported" block for WagerzonC describing it as full-risk-to-user, Sheet1 only.
- Architecture diagram: `runs all 6` → `runs all 7`; update Wagerzon×2 → Wagerzon×3.

### `bet_logger/CLAUDE.md`

Add one bullet to the Multi-account scrapers section:

> - **WagerzonC** (Wagerzon third account) — `bet_multiplier: 1.0` (full risk to user, no adjustment). Uploads adjusted to Sheet1 only; no Shared tab.

## Version control plan

- **Worktree:** `~/NFLWork/.worktrees/wagerzon-third-account` on branch `feature/wagerzon-third-account`. Already created.
- **Files modified:**
  - `bet_logger/scraper_wagerzon.py`
  - `bet_logger/test_scraper_wagerzon.py`
  - `bet_logger/.env.example`
  - `bet_logger/run_all_scrapers.sh`
  - `bet_logger/README.md`
  - `bet_logger/CLAUDE.md`
  - `docs/superpowers/specs/2026-04-30-wagerzon-third-account-design.md` (this file)
- **Files NOT touched:** `sheets.py`, `utils.py`, `scraper_bfa.py`, BFA-related files, `.env`.
- **Commits:** one commit per task as in the WagerzonJ plan.
- **Pre-merge testing (in worktree):**
  1. `pytest test_scraper_wagerzon.py -v` — 8 tests pass (7 existing + 1 new + 1 renamed)
  2. `bash -n run_all_scrapers.sh` — shell syntax OK
  3. `--help` shows `--account {default,j,c}`
  4. `--account c --dry-run` — login + parse + correct platform label, no Shared-tab upload
  5. `--account default --dry-run` — regression: 120 primary bets unchanged
- **Pre-merge review:** full `git diff main..HEAD` against the executive engineer checklist.
- **Merge:** only after explicit user approval.
- **Post-merge cleanup:** `git worktree remove`, `git branch -d`.
- **One-off live run after merge:** mirror the WagerzonJ flow — run `--account c` on `main` once, no dry-run, to backfill last week's bets.

## Edge cases & risks

1. **Dedup safety** — `WagerzonC` is a new platform label, so existing dedup key `(date, platform, description, amount)` won't collide with `Wagerzon` or `WagerzonJ`.
2. **Shared tab not touched** — `shared_sheet: None` skips the upload entirely. WagerzonC bets land in Sheet1 only.
3. **No multiplier math** — `bet_multiplier: 1.0` makes `round(risk * 1.0, 2) == risk`, and the print-format polish suppresses the `(raw $...)` suffix. Output identical to primary.
4. **Summary tabs** — MLB Summary and CBB Summary filter by sport, not platform. WagerzonC bets in those sports get included automatically at full amount, which is the correct full-exposure tracking.
5. **Run order** — WagerzonC runs after WagerzonJ in `run_all_scrapers.sh`. If primary or WagerzonJ auth fails, WagerzonC still runs (existing failure-tolerant `if/else` pattern).

## Out of scope

- Refactoring ACCOUNTS into a generic helper module shared with BFA — still YAGNI even with three Wagerzon accounts; the dict is small enough that duplication is cheaper than abstraction.
- Backfilling historical WagerzonC bets that were logged under `Wagerzon` — out of scope; user can manually relabel if needed.
- Changes to BFA's two-account handling.
