# WagerzonC Third Account Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a third Wagerzon account (`WagerzonC`) to bet_logger that scrapes alongside the primary and WagerzonJ. Full risk to user — no multiplier, no Shared tab. Behaves identically to the primary except for distinct credentials and platform label.

**Architecture:** Add one new entry to the `ACCOUNTS` dict in `scraper_wagerzon.py` with `bet_multiplier: 1.0` and `shared_sheet: None`. Multi-account scaffolding shipped with WagerzonJ already supports any number of entries — argparse choices, login parameterization, and dual-tab guard are all driven from the dict. Surrounding work is wiring (env vars, shell scheduler, docs) and a one-off live backfill.

**Tech Stack:** Python 3 (`requests`, `python-dotenv`, pytest 9 in `bet_logger/venv`), Bash, Google Sheets API (already wired).

**Worktree:** `~/NFLWork/.worktrees/wagerzon-third-account`
**Branch:** `feature/wagerzon-third-account`
**Spec:** `docs/superpowers/specs/2026-04-30-wagerzon-third-account-design.md`

---

## File Structure

**Modified:**
- `bet_logger/scraper_wagerzon.py` — add `'c'` entry to `ACCOUNTS` dict
- `bet_logger/test_scraper_wagerzon.py` — rename keys-test, add c-fields test
- `bet_logger/.env.example` — `WAGERZONC_USERNAME` / `WAGERZONC_PASSWORD` placeholders
- `bet_logger/run_all_scrapers.sh` — third Wagerzon block + count bump 6 → 7
- `bet_logger/README.md` — count, usage, platform notes (7 small edits)
- `bet_logger/CLAUDE.md` — bullet under Multi-account scrapers

**Not touched:** `sheets.py`, `utils.py`, `scraper_bfa.py`, `.env` (user adds creds locally; gitignored).

---

## Task 1: Add `'c'` entry to ACCOUNTS dict (TDD cycle)

**Files:**
- Modify: `bet_logger/scraper_wagerzon.py` — `ACCOUNTS` dict (around lines 36-51)
- Modify: `bet_logger/test_scraper_wagerzon.py` — keys test + new c-fields test

- [ ] **Step 1.1: Update test file**

In `bet_logger/test_scraper_wagerzon.py`:

a) Rename the existing `test_accounts_has_default_and_j` to `test_accounts_keys` and update the asserted set:

Find:
```python
def test_accounts_has_default_and_j():
    assert set(ACCOUNTS.keys()) == {'default', 'j'}
```

Replace with:
```python
def test_accounts_keys():
    assert set(ACCOUNTS.keys()) == {'default', 'j', 'c'}
```

b) Add a new test function after `test_accounts_j_fields`:

```python
def test_accounts_c_fields():
    acct = ACCOUNTS['c']
    assert acct['username_env'] == 'WAGERZONC_USERNAME'
    assert acct['password_env'] == 'WAGERZONC_PASSWORD'
    assert acct['platform'] == 'WagerzonC'
    assert acct['bet_multiplier'] == 1.0
    assert acct['shared_sheet'] is None
```

- [ ] **Step 1.2: Run tests to verify they fail**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-third-account/bet_logger
/Users/callancapitolo/NFLWork/bet_logger/venv/bin/python3 -m pytest test_scraper_wagerzon.py -v
```

Expected: 2 failures
- `test_accounts_keys` — `AssertionError: {'default', 'j'} != {'default', 'j', 'c'}` (or similar)
- `test_accounts_c_fields` — `KeyError: 'c'`

If anything else fails or no tests fail at all, stop and investigate before proceeding.

- [ ] **Step 1.3: Add `'c'` entry to `ACCOUNTS`**

Open `bet_logger/scraper_wagerzon.py`. Find the existing `ACCOUNTS` dict (currently the `'j'` entry ends with `'shared_sheet': 'Shared',` followed by `},` and then `}`). Insert the new `'c'` entry between the closing `},` of `'j'` and the dict-closing `}`:

Find:
```python
    'j': {
        'username_env': 'WAGERZONJ_USERNAME',
        'password_env': 'WAGERZONJ_PASSWORD',
        'platform': 'WagerzonJ',
        'bet_multiplier': 0.875,
        'shared_sheet': 'Shared',
    },
}
```

Replace with:
```python
    'j': {
        'username_env': 'WAGERZONJ_USERNAME',
        'password_env': 'WAGERZONJ_PASSWORD',
        'platform': 'WagerzonJ',
        'bet_multiplier': 0.875,
        'shared_sheet': 'Shared',
    },
    'c': {
        'username_env': 'WAGERZONC_USERNAME',
        'password_env': 'WAGERZONC_PASSWORD',
        'platform': 'WagerzonC',
        'bet_multiplier': 1.0,
        'shared_sheet': None,
    },
}
```

- [ ] **Step 1.4: Run tests to verify they pass**

```bash
/Users/callancapitolo/NFLWork/bet_logger/venv/bin/python3 -m pytest test_scraper_wagerzon.py -v
```

Expected: 8 PASS (the renamed `test_accounts_keys`, the new `test_accounts_c_fields`, plus the 6 unchanged tests from prior work).

- [ ] **Step 1.5: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-third-account
git add bet_logger/scraper_wagerzon.py bet_logger/test_scraper_wagerzon.py
git commit -m "feat(bet-logger): add WagerzonC third account to ACCOUNTS"
```

---

## Task 2: Add `WAGERZONC_*` env vars to `.env.example`

**Files:**
- Modify: `bet_logger/.env.example`

- [ ] **Step 2.1: Read the file to find the WagerzonJ section**

```bash
grep -n -i "wagerzon" bet_logger/.env.example
```

Locate the line beginning `WAGERZONJ_PASSWORD=`.

- [ ] **Step 2.2: Append the new section**

Right after the existing `WAGERZONJ_PASSWORD=` line, append:

```
# Wagerzon third account (WagerzonC — full risk to user, no adjustment)
WAGERZONC_USERNAME=
WAGERZONC_PASSWORD=
```

If the existing file uses different formatting between sections (blank-line separators, etc.), match the style of the surrounding WagerzonJ block. Report what convention was matched.

- [ ] **Step 2.3: Commit**

```bash
git add bet_logger/.env.example
git commit -m "feat(bet-logger): document WAGERZONC_* env vars in .env.example"
```

---

## Task 3: Wire WagerzonC into `run_all_scrapers.sh`

**Files:**
- Modify: `bet_logger/run_all_scrapers.sh`

- [ ] **Step 3.1: Insert WagerzonC block after the WagerzonJ block**

Find the existing WagerzonJ block (begins `# Run Wagerzon scraper — WagerzonJ account`, ends with `echo ""`). Right after that block's `echo ""` line and BEFORE the `# Run Hoop88 scraper` comment, insert:

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

- [ ] **Step 3.2: Update success message count**

Find:
```bash
    MSG="All 6 scrapers completed successfully."
```

Replace with:
```bash
    MSG="All 7 scrapers completed successfully."
```

- [ ] **Step 3.3: Lint the shell script**

```bash
bash -n bet_logger/run_all_scrapers.sh && echo "OK"
```
Expected: `OK`.

- [ ] **Step 3.4: Commit**

```bash
git add bet_logger/run_all_scrapers.sh
git commit -m "feat(bet-logger): run WagerzonC in run_all_scrapers.sh"
```

---

## Task 4: Update `bet_logger/README.md`

**Files:**
- Modify: `bet_logger/README.md` (multiple sections)

Read the file first to confirm the current text matches each "from" version. If it doesn't, stop and report the discrepancy.

- [ ] **Step 4.1: Top-line description**

- FROM: `Multi-platform bet history scraper. Scrapes settled bets from Wagerzon (2 accounts), Hoop88, BFA Gaming (2 accounts), and BetOnline, then uploads to Google Sheets with automatic duplicate detection.`
- TO: `Multi-platform bet history scraper. Scrapes settled bets from Wagerzon (3 accounts), Hoop88, BFA Gaming (2 accounts), and BetOnline, then uploads to Google Sheets with automatic duplicate detection.`

- [ ] **Step 4.2: `.env` setup mention**

- FROM: `Edit \`.env\` with credentials for each platform (Wagerzon, WagerzonJ, Hoop88, BFA primary, BFAJ) and your Google Sheet ID.`
- TO: `Edit \`.env\` with credentials for each platform (Wagerzon, WagerzonJ, WagerzonC, Hoop88, BFA primary, BFAJ) and your Google Sheet ID.`

- [ ] **Step 4.3: "Run all scrapers" description**

- FROM: `Runs Wagerzon (primary + WagerzonJ), Hoop88, BFA (primary + BFAJ), and BetOnline sequentially. Continues to the next scraper if one fails.`
- TO: `Runs Wagerzon (primary + WagerzonJ + WagerzonC), Hoop88, BFA (primary + BFAJ), and BetOnline sequentially. Continues to the next scraper if one fails.`

- [ ] **Step 4.4: Wagerzon usage line in "Run individual scrapers"**

- FROM: `./venv/bin/python3 scraper_wagerzon.py     [--weeks 1] [--all-weeks] [--account j] [--dry-run]`
- TO: `./venv/bin/python3 scraper_wagerzon.py     [--weeks 1] [--all-weeks] [--account j|c] [--dry-run]`

- [ ] **Step 4.5: Wagerzon flags subsection**

Find the existing block:
```
**Wagerzon flags:**
- `--weeks N` — Which week to fetch (0=This Week, 1=Last Week, default: 1)
- `--all-weeks` — Fetch every week back-to-back until an empty one (catches up after missed runs)
- `--account j` — Scrape the WagerzonJ second account instead of primary
```

Replace the last bullet line so the block reads:
```
**Wagerzon flags:**
- `--weeks N` — Which week to fetch (0=This Week, 1=Last Week, default: 1)
- `--all-weeks` — Fetch every week back-to-back until an empty one (catches up after missed runs)
- `--account j|c` — Scrape the WagerzonJ (87.5% of risk) or WagerzonC (full risk) account instead of primary
```

- [ ] **Step 4.6: Output Columns Platform example**

- FROM: `| B | Platform | Wagerzon / WagerzonJ / Hoop88 / Betfastaction / BFAJ / BetOnline |`
- TO: `| B | Platform | Wagerzon / WagerzonJ / WagerzonC / Hoop88 / Betfastaction / BFAJ / BetOnline |`

- [ ] **Step 4.7: Wagerzon Platform Notes — add WagerzonC bullet**

Find the existing two-bullet block:

```markdown
- **Primary (`Wagerzon`)** — full risk attributed to user. Bets land in `Sheet1` only.
- **WagerzonJ (`--account j`)** — partner-shared account. User holds 87.5% of risk. Adjusted bets (`risk × 0.875`, rounded to cents) land in `Sheet1`; raw (full) bets land in the `Shared` tab for week-over-week balance reconciliation. Same pattern as BFAJ but with a multiplicative adjustment instead of BFAJ's flat $-15.
```

Append a third bullet AFTER the WagerzonJ bullet so the block reads:

```markdown
- **Primary (`Wagerzon`)** — full risk attributed to user. Bets land in `Sheet1` only.
- **WagerzonJ (`--account j`)** — partner-shared account. User holds 87.5% of risk. Adjusted bets (`risk × 0.875`, rounded to cents) land in `Sheet1`; raw (full) bets land in the `Shared` tab for week-over-week balance reconciliation. Same pattern as BFAJ but with a multiplicative adjustment instead of BFAJ's flat $-15.
- **WagerzonC (`--account c`)** — third account, full risk attributed to user (no multiplier). Bets land in `Sheet1` only — no `Shared` tab upload. Behaves identically to the primary except for separate credentials and the `WagerzonC` platform label.
```

- [ ] **Step 4.8: Architecture diagram count**

- FROM: `run_all_scrapers.sh          # Entry point (runs all 6: Wagerzon×2, Hoop88, BFA×2, BetOnline)`
- TO: `run_all_scrapers.sh          # Entry point (runs all 7: Wagerzon×3, Hoop88, BFA×2, BetOnline)`

- [ ] **Step 4.9: Confirm changes**

```bash
git diff bet_logger/README.md
```
Verify exactly the eight targeted edits, no whitespace noise.

- [ ] **Step 4.10: Commit**

```bash
git add bet_logger/README.md
git commit -m "docs(bet-logger): document WagerzonC third account in README"
```

---

## Task 5: Add WagerzonC bullet to `bet_logger/CLAUDE.md`

**Files:**
- Modify: `bet_logger/CLAUDE.md`

- [ ] **Step 5.1: Locate the Multi-account scrapers section**

```bash
grep -n "Multi-account\|WagerzonJ\|BFAJ" bet_logger/CLAUDE.md
```

- [ ] **Step 5.2: Add the WagerzonC bullet**

Find the existing bullet block in the "Multi-account scrapers" section:

```markdown
- **BFAJ** (BFA second account) — flat `bet_adjustment: -15` (subtract $15 per bet). Uploads adjusted to Sheet1 and raw to `Shared` tab.
- **WagerzonJ** (Wagerzon second account) — multiplicative `bet_multiplier: 0.875` (user holds 87.5% of risk). Uploads adjusted to Sheet1 and raw to `Shared` tab.
```

Append a third bullet (right after the WagerzonJ line, before the next paragraph) so the block reads:

```markdown
- **BFAJ** (BFA second account) — flat `bet_adjustment: -15` (subtract $15 per bet). Uploads adjusted to Sheet1 and raw to `Shared` tab.
- **WagerzonJ** (Wagerzon second account) — multiplicative `bet_multiplier: 0.875` (user holds 87.5% of risk). Uploads adjusted to Sheet1 and raw to `Shared` tab.
- **WagerzonC** (Wagerzon third account) — `bet_multiplier: 1.0` (full risk to user, no adjustment). Uploads adjusted to Sheet1 only; no Shared tab.
```

- [ ] **Step 5.3: Commit**

```bash
git add bet_logger/CLAUDE.md
git commit -m "docs(bet-logger): note WagerzonC third account in CLAUDE.md"
```

---

## Task 6: Dry-run verification (manual gate)

**Files:** none (verification only)

**Prerequisites:** The user must have populated `WAGERZONC_USERNAME` / `WAGERZONC_PASSWORD` in `bet_logger/.env`. If not, **stop here and ask the user to add the credentials before continuing.**

The `.env` is in the main worktree (`/Users/callancapitolo/NFLWork/bet_logger/.env`), not the feature worktree. Either symlink it into the feature worktree (the WagerzonJ work used this approach: `ln -s /Users/callancapitolo/NFLWork/bet_logger/.env /Users/callancapitolo/NFLWork/.worktrees/wagerzon-third-account/bet_logger/.env`) or run the verification from the main worktree's bet_logger directory pointing at the feature-worktree's code via PYTHONPATH (clunkier — symlink wins).

- [ ] **Step 6.1: Confirm `.env` has the new vars (do not print values)**

```bash
grep -c "^WAGERZONC_" /Users/callancapitolo/NFLWork/bet_logger/.env
```
Expected: `2`. If `0` or `1`, ask the user to fill in `.env` before proceeding.

- [ ] **Step 6.2: Symlink `.env` into the feature worktree**

```bash
ln -s /Users/callancapitolo/NFLWork/bet_logger/.env /Users/callancapitolo/NFLWork/.worktrees/wagerzon-third-account/bet_logger/.env
```

(If a symlink already exists from prior work, this command will error; ignore.)

- [ ] **Step 6.3: Dry-run primary account (regression check)**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-third-account/bet_logger
/Users/callancapitolo/NFLWork/bet_logger/venv/bin/python3 scraper_wagerzon.py --account default --dry-run
```

Expected:
- Login succeeds as primary
- Each bet line shows `$X.XX - {result}` with **no** `(raw $X.XX)` suffix (multiplier is 1.0; suppressed by the polish from prior work)
- Final summary `Successfully scraped N bets from Wagerzon`
- "Dry run — skipping upload" message
- No Shared-tab upload attempted

If the bet count seems off vs. typical weekly volume, stop and report.

- [ ] **Step 6.4: Dry-run WagerzonC**

```bash
/Users/callancapitolo/NFLWork/bet_logger/venv/bin/python3 scraper_wagerzon.py --account c --dry-run
```

Expected:
- Login succeeds as WagerzonC (`Logging in to Wagerzon (WagerzonC)...`)
- Each bet line shows `$X.XX - {result}` with **no** `(raw $X.XX)` suffix
- Final summary `Successfully scraped N bets from WagerzonC`
- "Dry run — skipping upload" message
- No multiplier summary block (only prints when `bet_multiplier != 1.0`)
- No Shared-tab upload attempted

- [ ] **Step 6.5: Verify `--help` shows the new choice**

```bash
/Users/callancapitolo/NFLWork/bet_logger/venv/bin/python3 scraper_wagerzon.py --help
```

Expected: `--account {default,j,c}` line.

- [ ] **Step 6.6: Run all unit tests one final time**

```bash
/Users/callancapitolo/NFLWork/bet_logger/venv/bin/python3 -m pytest test_scraper_wagerzon.py test_mlb_summary.py -v
```

Expected: 23 PASS (8 wagerzon + 15 mlb_summary).

- [ ] **Step 6.7: Document verification results**

Report to the user:
- Primary regression: bet count, plausibility check
- WagerzonC: bet count, login confirmation, no multiplier suffix in output
- `--help` shows the new account choice
- All tests pass

No commit — verification only.

---

## Task 7: Pre-merge review and merge approval

**Files:** none (review only)

- [ ] **Step 7.1: Show the full feature diff**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-third-account
MERGE_BASE=$(git merge-base main HEAD)
git log --oneline $MERGE_BASE..HEAD
git diff --stat $MERGE_BASE..HEAD
```

(Use merge-base because main has progressed since the worktree was created.)

- [ ] **Step 7.2: Run executive engineer review checklist**

Confirm in writing back to the user:

- **Data integrity:**
  - `WagerzonC` is a new platform label → no dedup collision with `Wagerzon` or `WagerzonJ`
  - Sheet1 written once per bet; no Shared upload (because `shared_sheet: None`)
  - `parse_api_bets` multiplier branch is no-op (1.0)
- **Resource safety:** N/A (pure HTTP, no DB connections)
- **Edge cases:**
  - First-run with no WagerzonC bets → "No bets found, exit 1" (existing path)
  - Auth failure on WagerzonC doesn't block primary or WagerzonJ (failure-tolerant `if/else` in run_all_scrapers.sh)
  - Empty-week stop logic preserved
  - argparse `choices=list(ACCOUNTS.keys())` auto-picks `'c'`
- **Dead code:** No new unused params, imports, branches
- **Log/disk hygiene:** stdout only
- **Security:** Credentials read from env, never logged. Missing-cred error shows env-var names, not values.

- [ ] **Step 7.3: Surface findings**

Post the diffstat, the checklist results, and any open issues. List ISSUES TO FIX vs ACCEPTABLE RISKS. Do NOT merge yet.

- [ ] **Step 7.4: If issues found, fix and re-test**

Fix on the same feature branch, re-run Task 6 dry-runs, return to Step 7.3.

- [ ] **Step 7.5: Get explicit user approval to merge**

Wait for user to say "yes, merge" or equivalent. **Never merge without explicit approval.**

- [ ] **Step 7.6: Merge to main and clean up**

Once approved:

```bash
cd /Users/callancapitolo/NFLWork
git checkout main
git status            # confirm clean (only unrelated runtime artifacts)
git merge --no-ff feature/wagerzon-third-account -m "Merge feat/wagerzon-third-account: WagerzonC full-risk account"
git worktree remove .worktrees/wagerzon-third-account
git branch -d feature/wagerzon-third-account
```

- [ ] **Step 7.7: Confirm post-merge state**

```bash
git log --oneline -5
git worktree list
git branch
```

Expected: feature branch gone, no stale worktree, latest commit is the merge.

---

## Task 8: One-off live backfill for last week

**Files:** none (live operation)

User explicitly requested a one-off live run after merge to populate last week's WagerzonC bets, mirroring the WagerzonJ flow.

- [ ] **Step 8.1: Confirm we're on `main` after merge**

```bash
cd /Users/callancapitolo/NFLWork
git branch --show-current
```
Expected: `main`.

- [ ] **Step 8.2: Live run for last week**

```bash
cd /Users/callancapitolo/NFLWork/bet_logger
./venv/bin/python3 scraper_wagerzon.py --account c
```

(No `--dry-run`. No `--weeks` override — default is 1, which is last week.)

Expected:
- Login as WagerzonC
- Bets scraped and printed
- "Checking for duplicates..." then upload to Sheet1
- "Successfully updated N cells" / "Added N new bets to sheet"
- No Shared-tab upload (`shared_sheet: None`)
- Exit 0

- [ ] **Step 8.3: Document outcome**

Report bet count, sheet rows added, and any duplicates skipped. If 0 bets are found, that's fine — the account may simply have had no activity last week.

---

## Self-review

**Spec coverage:** Every section of `2026-04-30-wagerzon-third-account-design.md` maps to a task:
- ACCOUNTS dict + tests → Task 1
- `.env.example` → Task 2
- `run_all_scrapers.sh` → Task 3
- `README.md` (8 edits) → Task 4
- `CLAUDE.md` → Task 5
- Pre-merge testing → Task 6
- Pre-merge review + merge → Task 7
- One-off live run → Task 8 (added because the user explicitly requested it for WagerzonJ and asked for the same here in conversation)

**Placeholder scan:** No "TBD" / "TODO" / "implement later" / "add appropriate handling". Every step has actual code or commands.

**Type consistency:** Field names (`username_env`, `password_env`, `platform`, `bet_multiplier`, `shared_sheet`) are identical across Tasks 1–7. The platform label `WagerzonC` is consistent in Tasks 1, 4, 5. The CLI value `c` (paired with `--account c`) is consistent in Tasks 1, 3, 4, 6, 8.
