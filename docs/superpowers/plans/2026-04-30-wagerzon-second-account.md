# WagerzonJ Second Account Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a second Wagerzon account (`WagerzonJ`) to bet_logger that scrapes alongside the primary, applies a 0.875 risk multiplier to user-share P&L in `Sheet1`, and uploads raw amounts to the existing `Shared` tab for week-over-week reconciliation.

**Architecture:** Mirror the existing BFA two-account pattern — module-level `ACCOUNTS` dict, `--account` CLI flag, per-account credentials/platform/multiplier/shared_sheet config. WagerzonJ is the only entry that needs a multiplier; the primary uses 1.0. Adjusted bets land in `Sheet1`, raw bets in `Shared`. Single second invocation in `run_all_scrapers.sh`.

**Tech Stack:** Python 3 (`requests`, `python-dotenv`, `pytest` 9 in `bet_logger/venv`), Bash, Google Sheets API (already wired in `sheets.py`).

**Worktree:** `~/NFLWork/.worktrees/wagerzon-second-account`
**Branch:** `feature/wagerzon-second-account`
**Spec:** `docs/superpowers/specs/2026-04-29-wagerzon-second-account-design.md`

---

## File Structure

**Created:**
- `bet_logger/test_scraper_wagerzon.py` — pytest unit tests for `ACCOUNTS` and `parse_api_bets` multiplier semantics

**Modified:**
- `bet_logger/scraper_wagerzon.py` — `ACCOUNTS` dict, refactored `login()` / `scrape_wagerzon()` / `parse_api_bets()`, `--account` CLI flag, dual-tab upload in `__main__`
- `bet_logger/run_all_scrapers.sh` — second WZ block, success-message count bump
- `bet_logger/.env.example` — `WAGERZONJ_USERNAME` / `WAGERZONJ_PASSWORD` placeholders
- `bet_logger/README.md` — usage, .env, platform notes, output columns
- `bet_logger/CLAUDE.md` — architecture note on two-account variants

**Not touched:** `bet_logger/sheets.py`, `bet_logger/utils.py`, `bet_logger/scraper_bfa.py`, `bet_logger/.env` (user adds credentials locally; gitignored).

---

## Task 1: Add `ACCOUNTS` dict and lock its shape with tests

**Files:**
- Modify: `bet_logger/scraper_wagerzon.py:30-31`
- Create: `bet_logger/test_scraper_wagerzon.py`

- [ ] **Step 1.1: Write failing tests for `ACCOUNTS` shape**

Create `bet_logger/test_scraper_wagerzon.py`:

```python
"""Unit tests for scraper_wagerzon.py — account config and bet parsing.

These tests run without hitting Wagerzon's live API. Network-dependent
behaviors (login, fetch_history_json) are validated end-to-end via
--dry-run runs in Task 9.
"""
import copy
import pytest
from scraper_wagerzon import ACCOUNTS, parse_api_bets


def test_accounts_has_default_and_j():
    assert set(ACCOUNTS.keys()) == {'default', 'j'}


def test_accounts_default_fields():
    acct = ACCOUNTS['default']
    assert acct['username_env'] == 'WAGERZON_USERNAME'
    assert acct['password_env'] == 'WAGERZON_PASSWORD'
    assert acct['platform'] == 'Wagerzon'
    assert acct['bet_multiplier'] == 1.0
    assert acct['shared_sheet'] is None


def test_accounts_j_fields():
    acct = ACCOUNTS['j']
    assert acct['username_env'] == 'WAGERZONJ_USERNAME'
    assert acct['password_env'] == 'WAGERZONJ_PASSWORD'
    assert acct['platform'] == 'WagerzonJ'
    assert acct['bet_multiplier'] == 0.875
    assert acct['shared_sheet'] == 'Shared'
```

- [ ] **Step 1.2: Run tests to verify they fail**

Run from the worktree:
```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-second-account/bet_logger
./venv/bin/python3 -m pytest test_scraper_wagerzon.py -v
```
Expected: ImportError on `ACCOUNTS` because it doesn't exist yet.

- [ ] **Step 1.3: Add `ACCOUNTS` dict to `scraper_wagerzon.py`**

Replace lines 30–31 in `bet_logger/scraper_wagerzon.py` (which currently read):

```python
WAGERZON_USERNAME = os.getenv("WAGERZON_USERNAME")
WAGERZON_PASSWORD = os.getenv("WAGERZON_PASSWORD")
```

with:

```python
# Account configurations.
# Wagerzon supports two accounts using the same credentials flow but
# different sheet labels and stake adjustments. The primary uses a
# multiplier of 1.0 (full risk attributed to user). WagerzonJ is a
# partner account where the user holds 87.5% of the risk; raw bets
# also land on the Shared tab for week-over-week reconciliation.
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

- [ ] **Step 1.4: Run tests to verify they pass**

```bash
./venv/bin/python3 -m pytest test_scraper_wagerzon.py -v
```
Expected: 3 PASS.

- [ ] **Step 1.5: Commit**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-second-account
git add bet_logger/test_scraper_wagerzon.py bet_logger/scraper_wagerzon.py
git commit -m "feat(bet-logger): ACCOUNTS dict + tests for Wagerzon two-account scaffold"
```

---

## Task 2: Refactor `parse_api_bets` to apply `bet_multiplier`

**Files:**
- Modify: `bet_logger/scraper_wagerzon.py:236-324` (`parse_api_bets` function body)
- Modify: `bet_logger/test_scraper_wagerzon.py` (add multiplier tests)

- [ ] **Step 2.1: Write failing tests for multiplier behavior**

Append to `bet_logger/test_scraper_wagerzon.py`:

```python
# Minimal fixture matching the HistoryHelper.aspx JSON shape.
# - Single straight bet, MLB, $100 risk, lost, with American odds -110.
# - Sufficient to exercise platform/multiplier/raw-risk paths.
SAMPLE_HISTORY = {
    'details': [
        {
            'wager': [
                {
                    'WagerOrTrans': 'WAGER',
                    'PlacedDate': '04/29/2026',
                    'RiskAmount': '100.00',
                    'WinAmount': '90.91',
                    'WinLoss': '-100.00',
                    'Result': 'LOSE',
                    'HeaderDesc': 'STRAIGHT BET',
                    'details': [
                        {
                            'IdSport': 'MLB',
                            'DetailDesc': '[967] NY YANKEES -1.5-110 (NY YANKEES vrs SF GIANTS)',
                        }
                    ],
                }
            ]
        }
    ]
}


def test_parse_api_bets_default_platform_and_no_multiplier():
    bets = parse_api_bets(SAMPLE_HISTORY, platform='Wagerzon', bet_multiplier=1.0)
    assert len(bets) == 1
    bet = bets[0]
    assert bet['platform'] == 'Wagerzon'
    assert bet['bet_amount'] == 100.0
    assert bet['_raw_risk'] == 100.0


def test_parse_api_bets_j_applies_multiplier_to_bet_amount_only():
    bets = parse_api_bets(SAMPLE_HISTORY, platform='WagerzonJ', bet_multiplier=0.875)
    assert len(bets) == 1
    bet = bets[0]
    assert bet['platform'] == 'WagerzonJ'
    assert bet['bet_amount'] == 87.5
    assert bet['_raw_risk'] == 100.0


def test_parse_api_bets_multiplier_does_not_change_odds_or_result():
    """Odds, decimal odds, and result must be invariant to multiplier."""
    default = parse_api_bets(SAMPLE_HISTORY, platform='Wagerzon', bet_multiplier=1.0)[0]
    j = parse_api_bets(SAMPLE_HISTORY, platform='WagerzonJ', bet_multiplier=0.875)[0]
    assert default['odds'] == j['odds']
    assert default['dec'] == j['dec']
    assert default['result'] == j['result']


def test_parse_api_bets_rounds_to_two_decimals():
    """0.875 × 33.33 = 29.16375 should round to 29.16, not write a long float."""
    history = copy.deepcopy(SAMPLE_HISTORY)
    history['details'][0]['wager'][0]['RiskAmount'] = '33.33'
    bets = parse_api_bets(history, platform='WagerzonJ', bet_multiplier=0.875)
    assert bets[0]['bet_amount'] == 29.16
    assert bets[0]['_raw_risk'] == 33.33
```

- [ ] **Step 2.2: Run tests to verify they fail**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-second-account/bet_logger
./venv/bin/python3 -m pytest test_scraper_wagerzon.py -v
```
Expected: 4 new tests FAIL (the function doesn't accept `platform`/`bet_multiplier` kwargs yet, and `_raw_risk` isn't set).

- [ ] **Step 2.3: Refactor `parse_api_bets` to accept `platform` and `bet_multiplier`**

In `bet_logger/scraper_wagerzon.py`, replace the `parse_api_bets` function (currently at line 236) with:

```python
def parse_api_bets(history: dict, platform: str = 'Wagerzon',
                   bet_multiplier: float = 1.0) -> list:
    """Convert HistoryHelper JSON to the standard bet dict format.

    The JSON has a 'details' list where each item is a day with a
    list of wagers. Each wager has legs in its own 'details' list.

    Args:
        history: The 'result' dict from HistoryHelper.aspx.
        platform: Sheet platform label (e.g. 'Wagerzon', 'WagerzonJ').
        bet_multiplier: Fraction of the wagered amount attributable to
                        the user. 1.0 for the primary account; 0.875
                        for the partner-shared WagerzonJ account.
                        The original risk is preserved in '_raw_risk'
                        so it can be uploaded to a verification tab.
    """
    bets = []
    days = history.get('details', [])

    for day in days:
        for wager in day.get('wager', []):
            try:
                # Skip non-wager transactions (cash in/out, transfers)
                if wager.get('WagerOrTrans') != 'WAGER':
                    continue

                risk = float(wager.get('RiskAmount', 0))
                win_potential = float(wager.get('WinAmount', 0))
                win_loss = float(wager.get('WinLoss', 0))
                result = map_result(wager.get('Result', ''))

                # Skip pending/open bets
                if not result:
                    continue

                # Use the actual settlement (WinLoss) to calculate odds.
                # WinLoss reflects the real payout — e.g. parlays with pushed
                # legs pay reduced odds, and WinLoss captures that.
                if result == 'win' and risk > 0 and win_loss > 0:
                    american_odds = calculate_american_odds(risk, win_loss)
                elif risk > 0 and win_potential > 0:
                    # For losses and pushes, use potential win to show the
                    # odds when the bet was placed.
                    american_odds = calculate_american_odds(risk, win_potential)
                else:
                    american_odds = 0
                decimal_odds = calculate_decimal_odds_from_american(american_odds) if american_odds else 0

                # Build description from legs
                legs = wager.get('details', [])
                leg_descs = []
                line = ''
                sport = ''

                for leg in legs:
                    raw = leg.get('DetailDesc', '')
                    cleaned = clean_desc(raw)
                    leg_descs.append(cleaned)

                    # Use first leg for line and sport
                    if not line:
                        line = parse_line_from_desc(cleaned)
                    if not sport:
                        # IdSport from the API is the primary source
                        id_sport = leg.get('IdSport', '')
                        sport = map_sport(id_sport) if id_sport else ''

                # Fallback: parse sport from description text if API didn't provide it
                if not sport:
                    sport = parse_sport(' '.join(leg_descs))

                # Prepend HeaderDesc to match old format:
                # "PARLAY (2 TEAMS) | leg1 | leg2" or "STRAIGHT BET | desc"
                header = wager.get('HeaderDesc', '').strip()
                leg_text = ' | '.join(leg_descs)
                description = f"{header} | {leg_text}" if header else leg_text
                bet_type = parse_bet_type(header)

                # Apply per-account stake adjustment.
                # bet_amount = user's share for Sheet1; _raw_risk = original
                # for the Shared verification tab.
                adjusted_risk = round(risk * bet_multiplier, 2)

                bet = {
                    'date': parse_date(wager.get('PlacedDate', '')),
                    'platform': platform,
                    'sport': sport,
                    'description': description,
                    'bet_type': bet_type,
                    'line': line,
                    'odds': american_odds,
                    'bet_amount': adjusted_risk,
                    'dec': decimal_odds,
                    'result': result,
                    '_raw_risk': risk,
                }

                bets.append(bet)
                print(f"  {bet['date']} - {bet_type} - {description[:60]} - ${adjusted_risk:.2f} (raw ${risk:.2f}) - {result}")

            except Exception as e:
                print(f"  Error parsing wager: {e}")
                continue

    return bets
```

- [ ] **Step 2.4: Run tests to verify they pass**

```bash
./venv/bin/python3 -m pytest test_scraper_wagerzon.py -v
```
Expected: 7 tests PASS (3 from Task 1 + 4 new).

- [ ] **Step 2.5: Commit**

```bash
git add bet_logger/scraper_wagerzon.py bet_logger/test_scraper_wagerzon.py
git commit -m "feat(bet-logger): apply bet_multiplier in Wagerzon parse, preserve _raw_risk"
```

---

## Task 3: Refactor `login()` and `scrape_wagerzon()` to be account-aware

**Files:**
- Modify: `bet_logger/scraper_wagerzon.py:37-74` (`login` function)
- Modify: `bet_logger/scraper_wagerzon.py:330-383` (`scrape_wagerzon` function)

This is a pure refactor: behavior of the `default` account must remain unchanged. No new unit tests — exercised end-to-end in Task 9 dry-runs.

- [ ] **Step 3.1: Replace `login()` to take credentials as arguments**

Replace the existing `login` function (currently `def login(session: requests.Session):`) with:

```python
def login(session: requests.Session, username: str, password: str):
    """Login to Wagerzon via ASP.NET form POST.

    Same auth approach as wagerzon_odds/scraper_v2.py:
    1. GET the login page to capture __VIEWSTATE and other hidden fields
    2. POST credentials with those hidden fields
    3. Session cookie (ASP.NET_SessionId) maintains auth for subsequent requests
    """
    resp = session.get(WAGERZON_BASE_URL, timeout=15)
    resp.raise_for_status()

    # If already redirected to a logged-in page, session is still valid
    if "History" in resp.url or "NewSchedule" in resp.url or "Welcome" in resp.url:
        print("Already authenticated")
        return

    html = resp.text

    # Extract ASP.NET hidden fields from the login form.
    # ASP.NET uses these to validate form submissions — they're generated
    # server-side and must be posted back exactly as received.
    fields = {}
    for name in ["__VIEWSTATE", "__VIEWSTATEGENERATOR", "__EVENTVALIDATION",
                 "__EVENTTARGET", "__EVENTARGUMENT"]:
        match = re.search(rf'(?:name|id)="{name}"[^>]*value="([^"]*)"', html)
        if match:
            fields[name] = match.group(1)

    if "__VIEWSTATE" not in fields:
        raise RuntimeError("Could not find __VIEWSTATE on login page — page structure may have changed")

    fields["Account"] = username
    fields["Password"] = password
    fields["BtnSubmit"] = ""

    resp = session.post(WAGERZON_BASE_URL, data=fields, timeout=15)
    resp.raise_for_status()
    print("Logged in successfully")
```

- [ ] **Step 3.2: Refactor `scrape_wagerzon()` to look up account config**

Replace the entire `scrape_wagerzon` function (currently `def scrape_wagerzon(weeks_back: int = 1, all_weeks: bool = False) -> list:`) with:

```python
def scrape_wagerzon(weeks_back: int = 1, all_weeks: bool = False,
                    account_name: str = 'default') -> list:
    """
    Log into Wagerzon via HTTP and fetch bet history from the JSON API.

    Uses HistoryHelper.aspx — the same endpoint the React frontend calls.
    Returns structured data with sport IDs, so no HTML parsing or team
    name guessing needed.

    Args:
        weeks_back: Which week to fetch (0 = current, 1 = last week, etc.)
        all_weeks: If True, fetch weeks 0 through N until an empty week is
                   found. Useful for catching up after missed weeks.
        account_name: Account key from ACCOUNTS ('default' or 'j').

    Returns:
        List of parsed bet dictionaries with platform/multiplier applied.
    """
    acct = ACCOUNTS[account_name]
    username = os.getenv(acct['username_env'])
    password = os.getenv(acct['password_env'])
    if not username or not password:
        raise ValueError(
            f"{acct['username_env']} and {acct['password_env']} must be set in .env file"
        )

    session = requests.Session()
    session.headers.update({
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36",
    })

    # Step 1: Authenticate
    print(f"Logging in to Wagerzon ({acct['platform']})...")
    login(session, username, password)

    # Step 2: Determine which weeks to fetch
    if all_weeks:
        weeks_to_fetch = range(0, 10)  # cap at 10 weeks back
    else:
        weeks_to_fetch = [weeks_back]

    # Step 3: Fetch and parse each week
    all_bets = []
    for week in weeks_to_fetch:
        week_label = 'This Week' if week == 0 else f'Week {week}'
        print(f"Fetching bet history ({week_label})...")
        history = fetch_history_json(session, week=week)

        days = history.get('details', [])
        total_wagers = sum(len(d.get('wager', [])) for d in days)
        print(f"  {total_wagers} wagers across {len(days)} days")

        if total_wagers == 0 and all_weeks and week > 0:
            print(f"  Empty week — stopping.\n")
            break

        bets = parse_api_bets(
            history,
            platform=acct['platform'],
            bet_multiplier=acct['bet_multiplier'],
        )
        all_bets.extend(bets)
        print()

    return all_bets
```

- [ ] **Step 3.3: Run unit tests to confirm refactor didn't break anything**

```bash
./venv/bin/python3 -m pytest test_scraper_wagerzon.py -v
```
Expected: 7 PASS. (Tests don't exercise login/scrape_wagerzon, but they import the module — confirms no syntax errors.)

- [ ] **Step 3.4: Commit**

```bash
git add bet_logger/scraper_wagerzon.py
git commit -m "refactor(bet-logger): make Wagerzon login/scrape account-parameterized"
```

---

## Task 4: Add `--account` CLI flag and dual-tab upload in `__main__`

**Files:**
- Modify: `bet_logger/scraper_wagerzon.py:386-427` (the `__main__` block)

- [ ] **Step 4.1: Replace the `__main__` block**

Replace the entire `if __name__ == "__main__":` block (currently at line 386) with:

```python
if __name__ == "__main__":
    import argparse
    import copy

    parser = argparse.ArgumentParser(description='Scrape bet history from Wagerzon')
    parser.add_argument('--weeks', type=int, default=1,
                        help='Weeks back to fetch (0=This Week, 1=Last Week, default: 1)')
    parser.add_argument('--all-weeks', action='store_true',
                        help='Fetch all weeks until an empty one (catches up after missed weeks)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Scrape but do not upload to Google Sheets')
    parser.add_argument('--account', choices=list(ACCOUNTS.keys()), default='default',
                        help='Which Wagerzon account to scrape (default or j)')
    args = parser.parse_args()

    acct = ACCOUNTS[args.account]

    print("=" * 60)
    print(f"WAGERZON BET HISTORY SCRAPER — {acct['platform']}")
    print("=" * 60)

    try:
        bets = scrape_wagerzon(weeks_back=args.weeks,
                               all_weeks=args.all_weeks,
                               account_name=args.account)
    except Exception as e:
        print(f"\n❌ Error: {e}")
        sys.exit(1)

    print(f"\n{'=' * 60}")
    print(f"Successfully scraped {len(bets)} bets from {acct['platform']}")
    print(f"{'=' * 60}\n")

    if bets and not args.dry_run:
        from sheets import append_bets_to_sheet

        # Adjusted bets (user's share) → main Sheet1.
        result = append_bets_to_sheet(bets)

        if result['status'] == 'success':
            print(f"\n✅ SUCCESS! Added {result['rows_added']} new bets to sheet")
            print(f"   Rows {result['start_row']} to {result['end_row']}")
        elif result['status'] == 'skipped':
            print(f"\n⚠️  {result['message']}")
        else:
            print(f"\n❌ Error uploading to sheets: {result.get('message', 'Unknown error')}")

        # Raw bets (original full risk) → Shared verification tab, if configured.
        if acct['shared_sheet']:
            print(f"\nUploading raw bets to '{acct['shared_sheet']}' tab...")
            raw_bets = []
            for bet in bets:
                raw = copy.copy(bet)
                raw['bet_amount'] = raw.pop('_raw_risk', raw['bet_amount'])
                raw_bets.append(raw)
            raw_result = append_bets_to_sheet(raw_bets, sheet_name=acct['shared_sheet'])
            if raw_result['status'] == 'success':
                print(f"Added {raw_result['rows_added']} raw bets to {acct['shared_sheet']}")
            elif raw_result['status'] == 'skipped':
                print(f"{raw_result['message']}")
            else:
                print(f"Error uploading raw bets: {raw_result.get('message', 'Unknown error')}")

    elif args.dry_run:
        print("Dry run — skipping upload to Google Sheets")
        if acct['bet_multiplier'] != 1.0:
            print(f"\nBet multiplier: ×{acct['bet_multiplier']}")
            raw_total = sum(b.get('_raw_risk', 0) for b in bets)
            adj_total = sum(b['bet_amount'] for b in bets)
            print(f"Raw total wagered:      ${raw_total:,.2f}")
            print(f"Adjusted total wagered: ${adj_total:,.2f}")
    elif not bets:
        print("No bets found to upload")
        sys.exit(1)
```

- [ ] **Step 4.2: Confirm `--help` shows the new flag**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-second-account/bet_logger
./venv/bin/python3 scraper_wagerzon.py --help
```
Expected: output includes `--account {default,j}` line.

- [ ] **Step 4.3: Confirm tests still pass**

```bash
./venv/bin/python3 -m pytest test_scraper_wagerzon.py -v
```
Expected: 7 PASS.

- [ ] **Step 4.4: Commit**

```bash
git add bet_logger/scraper_wagerzon.py
git commit -m "feat(bet-logger): --account flag + dual-tab upload for Wagerzon"
```

---

## Task 5: Add `WAGERZONJ_*` env vars to `.env.example`

**Files:**
- Modify: `bet_logger/.env.example`

- [ ] **Step 5.1: Read current `.env.example` to find the Wagerzon section**

```bash
grep -n -i "wager" bet_logger/.env.example
```

- [ ] **Step 5.2: Append the new variables**

After the existing `WAGERZON_PASSWORD=` line, add:

```
# Wagerzon second account (WagerzonJ — user holds 87.5% of risk;
# remaining 12.5% covered by partner). Adjusted bets go to Sheet1,
# raw bets to the Shared tab.
WAGERZONJ_USERNAME=
WAGERZONJ_PASSWORD=
```

If the file uses a different formatting style for sections (review the file first), match the existing style.

- [ ] **Step 5.3: Commit**

```bash
git add bet_logger/.env.example
git commit -m "feat(bet-logger): document WAGERZONJ_* env vars in .env.example"
```

> **Note for the engineer:** The real credentials go in `bet_logger/.env`, which is gitignored. The user will populate that file separately — do not edit `.env` from the worktree.

---

## Task 6: Wire WagerzonJ into `run_all_scrapers.sh`

**Files:**
- Modify: `bet_logger/run_all_scrapers.sh:19-29` (after the existing Wagerzon block)
- Modify: `bet_logger/run_all_scrapers.sh` (success-message count)

- [ ] **Step 6.1: Insert the WagerzonJ block right after the existing Wagerzon block**

Find the existing Wagerzon block (lines 19–28, ending with `echo ""`). Right after that `echo ""` line, insert:

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

- [ ] **Step 6.2: Update the success message count**

Find the line:
```bash
    MSG="All 5 scrapers completed successfully."
```
Replace with:
```bash
    MSG="All 6 scrapers completed successfully."
```

- [ ] **Step 6.3: Lint the shell script**

```bash
bash -n bet_logger/run_all_scrapers.sh && echo "OK"
```
Expected: `OK`.

- [ ] **Step 6.4: Commit**

```bash
git add bet_logger/run_all_scrapers.sh
git commit -m "feat(bet-logger): run WagerzonJ in run_all_scrapers.sh"
```

---

## Task 7: Update `bet_logger/README.md`

**Files:**
- Modify: `bet_logger/README.md` (multiple sections — see steps)

- [ ] **Step 7.1: Update top-line description (line 3)**

Change:
```
Multi-platform bet history scraper. Scrapes settled bets from Wagerzon, Hoop88, BFA Gaming (2 accounts), and BetOnline, then uploads to Google Sheets with automatic duplicate detection.
```
to:
```
Multi-platform bet history scraper. Scrapes settled bets from Wagerzon (2 accounts), Hoop88, BFA Gaming (2 accounts), and BetOnline, then uploads to Google Sheets with automatic duplicate detection.
```

- [ ] **Step 7.2: Update `.env` setup mention (line 30)**

Change:
```
Edit `.env` with credentials for each platform (Wagerzon, Hoop88, BFA primary, BFAJ) and your Google Sheet ID.
```
to:
```
Edit `.env` with credentials for each platform (Wagerzon, WagerzonJ, Hoop88, BFA primary, BFAJ) and your Google Sheet ID.
```

- [ ] **Step 7.3: Update "Run all scrapers" description (line 40)**

Change:
```
Runs Wagerzon, Hoop88, BFA (primary + BFAJ), and BetOnline sequentially. Continues to the next scraper if one fails.
```
to:
```
Runs Wagerzon (primary + WagerzonJ), Hoop88, BFA (primary + BFAJ), and BetOnline sequentially. Continues to the next scraper if one fails.
```

- [ ] **Step 7.4: Add `--account j` to the Wagerzon usage example (line 45)**

Change:
```
./venv/bin/python3 scraper_wagerzon.py
```
to:
```
./venv/bin/python3 scraper_wagerzon.py     [--weeks 1] [--all-weeks] [--account j] [--dry-run]
```

- [ ] **Step 7.5: Add a Wagerzon-specific flags subsection**

Right above the existing `**Hoop88 flags:**` block, insert:

```
**Wagerzon flags:**
- `--weeks N` — Which week to fetch (0=This Week, 1=Last Week, default: 1)
- `--all-weeks` — Fetch every week back-to-back until an empty one (catches up after missed runs)
- `--account j` — Scrape the WagerzonJ second account instead of primary
```

- [ ] **Step 7.6: Update the Output Columns Platform example**

Change:
```
| B | Platform | Wagerzon / Hoop88 / Betfastaction / BFAJ / BetOnline |
```
to:
```
| B | Platform | Wagerzon / WagerzonJ / Hoop88 / Betfastaction / BFAJ / BetOnline |
```

- [ ] **Step 7.7: Replace the Wagerzon Platform Notes paragraph**

Find the line beginning with `**Wagerzon** —` and replace that paragraph (one line) with:

```
**Wagerzon** — Scrapes the JSON history endpoint at `backend.wagerzon.com/wager/HistoryHelper.aspx`. Auto-login via ASP.NET form POST (same flow as `wagerzon_odds/scraper_v2.py`). Pure HTTP, no browser. Two accounts supported via `--account`:

- **Primary (`Wagerzon`)** — full risk attributed to user. Bets land in `Sheet1` only.
- **WagerzonJ (`--account j`)** — partner-shared account. User holds 87.5% of risk. Adjusted bets (`risk × 0.875`, rounded to cents) land in `Sheet1`; raw (full) bets land in the `Shared` tab for week-over-week balance reconciliation. Same pattern as BFAJ but with a multiplicative adjustment instead of BFAJ's flat $-15.
```

- [ ] **Step 7.8: Update Architecture diagram**

In the architecture block at the bottom, change:
```
  scraper_wagerzon.py        # Wagerzon scraper (headless Chromium)
```
to:
```
  scraper_wagerzon.py        # Wagerzon scraper (HTTP, supports --account j)
```

(Note: existing line says "headless Chromium" but the scraper is actually HTTP-only; this is a pre-existing doc accuracy fix.)

- [ ] **Step 7.9: Confirm changes look right**

```bash
git diff bet_logger/README.md
```
Expected: only the targeted edits, no whitespace noise.

- [ ] **Step 7.10: Commit**

```bash
git add bet_logger/README.md
git commit -m "docs(bet-logger): document WagerzonJ second account in README"
```

---

## Task 8: Add architecture note to `bet_logger/CLAUDE.md`

**Files:**
- Modify: `bet_logger/CLAUDE.md`

- [ ] **Step 8.1: Read the current CLAUDE.md to locate the Architecture section**

```bash
grep -n "Architecture\|## " bet_logger/CLAUDE.md
```

- [ ] **Step 8.2: Add a multi-account note**

Append after the existing Architecture bullet list (or after `run_all_scrapers.sh — Run all scrapers sequentially`):

```
## Multi-account scrapers

Both `scraper_bfa.py` and `scraper_wagerzon.py` support a second account via `--account j`:

- **BFAJ** (BFA second account) — flat `bet_adjustment: -15` (subtract $15 per bet). Uploads adjusted to Sheet1 and raw to `Shared` tab.
- **WagerzonJ** (Wagerzon second account) — multiplicative `bet_multiplier: 0.875` (user holds 87.5% of risk). Uploads adjusted to Sheet1 and raw to `Shared` tab.

The two scrapers each define their own `ACCOUNTS` dict; they're not yet sharing a generic helper (intentional — see spec 2026-04-29-wagerzon-second-account-design.md, Approaches 1 vs 2).
```

- [ ] **Step 8.3: Commit**

```bash
git add bet_logger/CLAUDE.md
git commit -m "docs(bet-logger): note Wagerzon two-account variant in CLAUDE.md"
```

---

## Task 9: End-to-end dry-run verification

**Files:** none (verification only)

**Prerequisites:** The user must have populated `WAGERZONJ_USERNAME` / `WAGERZONJ_PASSWORD` in `bet_logger/.env`. If not, **stop here and ask the user to add the credentials before continuing.**

- [ ] **Step 9.1: Confirm `.env` has the new vars (do not print values)**

```bash
grep -c "^WAGERZONJ_" bet_logger/.env
```
Expected: `2`. If `0` or `1`, ask the user to fill in `.env` before proceeding.

- [ ] **Step 9.2: Dry-run the primary account (regression check)**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-second-account/bet_logger
./venv/bin/python3 scraper_wagerzon.py --account default --dry-run
```
Expected:
- Logs in successfully as the primary
- Prints `Successfully scraped N bets from Wagerzon`
- Each bet line shows `$X.XX (raw $X.XX)` where the two values are equal
- No upload attempted

If the dry-run fails or the bet count seems implausible vs. the user's normal weekly volume, stop and report.

- [ ] **Step 9.3: Dry-run WagerzonJ**

```bash
./venv/bin/python3 scraper_wagerzon.py --account j --dry-run
```
Expected:
- Logs in successfully as WagerzonJ
- Prints `Successfully scraped N bets from WagerzonJ`
- Each bet line shows `$adj (raw $raw)` where `adj ≈ raw × 0.875`
- Final summary block prints `Bet multiplier: ×0.875`, `Raw total wagered`, `Adjusted total wagered`
- Adjusted total ≈ raw total × 0.875 (within rounding)
- No upload attempted

- [ ] **Step 9.4: Spot-check a single bet's math by hand**

Pick any bet from the WagerzonJ dry-run output. Verify:
- `adjusted_amount == round(raw_amount * 0.875, 2)`
- `odds` and `dec` match the raw description's odds (multiplier should not affect them)

If the math is off on any bet, stop and investigate.

- [ ] **Step 9.5: Run all unit tests one final time**

```bash
./venv/bin/python3 -m pytest test_scraper_wagerzon.py test_mlb_summary.py -v
```
Expected: all PASS. (Including the existing MLB summary tests as a regression check that nothing else was disturbed.)

- [ ] **Step 9.6: Document verification results in the conversation**

Report to the user:
- Primary dry-run: bet count, plausibility check
- WagerzonJ dry-run: bet count, raw total, adjusted total, multiplier confirmation
- Hand-checked bet math on at least one row

No commit — this task is verification only.

---

## Task 10: Pre-merge review and merge approval

**Files:** none (review only)

- [ ] **Step 10.1: Show the full feature diff**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/wagerzon-second-account
git log --oneline main..HEAD
git diff --stat main..HEAD
```

- [ ] **Step 10.2: Run executive engineer review checklist (per `~/NFLWork/CLAUDE.md`)**

For each item, confirm in writing (post a checklist back to the user):

- **Data integrity:**
  - `bet_multiplier` applied exactly once per bet (not twice)
  - `_raw_risk` always preserved before any adjustment
  - Dedup unaffected (distinct platform labels)
  - No risk of double-uploading WagerzonJ bets to Sheet1
- **Resource safety:** N/A (pure HTTP; no DB connections introduced)
- **Edge cases:**
  - First-run with no bets in either account → handled (`if bets and not args.dry_run`)
  - Empty week stop condition still works (`all_weeks` path unchanged)
  - Floating-point: `round(..., 2)` covers fractional cents
- **Dead code:** No unused params, imports, or branches introduced
- **Log/disk hygiene:** Only stdout prints; no new log files
- **Security:** Credentials read from env, never logged. No secrets printed.

- [ ] **Step 10.3: Surface findings to the user**

Post the diffstat, the checklist results, and any open questions. List any ISSUES TO FIX vs ACCEPTABLE RISKS. Do NOT merge yet.

- [ ] **Step 10.4: If issues found, fix and re-test**

Fix on the same feature branch, re-run Task 9 dry-runs after each fix, return to Step 10.3.

- [ ] **Step 10.5: Get explicit user approval to merge**

Wait for the user to say "yes, merge" or equivalent.

- [ ] **Step 10.6: Merge to main and clean up**

Once approved:

```bash
cd /Users/callancapitolo/NFLWork
# Verify main is clean of unrelated work first
git checkout main
git status
# Merge feature branch
git merge --no-ff feature/wagerzon-second-account -m "Merge feat/wagerzon-second-account: WagerzonJ partner-shared account"
# Clean up worktree and branch
git worktree remove .worktrees/wagerzon-second-account
git branch -d feature/wagerzon-second-account
```

- [ ] **Step 10.7: Confirm post-merge state**

```bash
git log --oneline -5
git worktree list
git branch
```
Expected: feature branch gone, no stale worktrees, latest commit is the merge.

- [ ] **Step 10.8: Live verification next Monday morning** (no action this run)

The first scheduled launchd run with the new account will be the next Monday at 5 AM. After that run, check:
- `bet_logger/logs/` for any errors in the WagerzonJ block
- Sheet1 for new `WagerzonJ` rows with adjusted amounts
- `Shared` tab for matching raw rows

If any issue surfaces, open a follow-up.

---

## Self-review notes

**Spec coverage:** All sections of `2026-04-29-wagerzon-second-account-design.md` map to tasks:
- ACCOUNTS dict → Task 1
- parse_api_bets multiplier → Task 2
- login/scrape_wagerzon refactor → Task 3
- --account flag + dual-tab upload → Task 4
- .env.example → Task 5
- run_all_scrapers.sh → Task 6
- README.md → Task 7
- CLAUDE.md → Task 8
- Pre-merge dry-run testing → Task 9
- Pre-merge review + merge → Task 10

**Placeholder scan:** No "TBD"/"TODO"/"add appropriate handling" — every step shows the exact code or command to run. Steps that depend on user action (Task 9 prerequisite, Task 10.5 approval) are explicit gates with stop conditions.

**Type consistency:** `bet_multiplier` (float), `platform` (str), `shared_sheet` (str|None), `_raw_risk` (float, leading underscore matches BFA's convention) — used consistently across Tasks 1, 2, 3, 4. CLI flag `--account` with `choices=list(ACCOUNTS.keys())` references the dict defined in Task 1.
