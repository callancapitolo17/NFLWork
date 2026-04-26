# MLB Summary Tab Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an `MLB Summary` tab to the bet-logging Google Sheet that auto-updates from `Sheet1` and shows aggregate performance of MLB spread+total correlated parlays placed on Wagerzon.

**Architecture:** A one-shot Python setup script (`bet_logger/create_mlb_summary.py`) creates the tab, writes pure spreadsheet formulas (no hidden helper sheet, no scheduled job), and adds two embedded charts via the Sheets API. After the initial run, the tab is self-maintaining — every formula references `Sheet1` directly with `SUMPRODUCT + REGEXMATCH`, so totals update live as new bets are scraped in. Charts re-render automatically when their backing data ranges change.

**Tech Stack:** Python 3, `google-api-python-client` (Sheets API v4), service-account auth via `google.oauth2.service_account`, `pytest` for unit tests. All existing and already installed in `bet_logger/venv`.

**Spec:** `docs/superpowers/specs/2026-04-24-mlb-summary-tab-design.md`

---

## File Structure

**Create:**
- `bet_logger/create_mlb_summary.py` — the setup script. Single file, mirrors `bet_logger/create_summary.py`'s structure but does not use a helper sheet and uses `SORT(UNIQUE(FILTER(…)))` for dynamic dates.
- `bet_logger/test_mlb_summary.py` — pytest file for the pure-Python filter helpers (regex classification, FG/F5 detection). Top-level in the `bet_logger/` directory, matching the existing `hoop88_correlation/test_*.py` convention elsewhere in the repo.

**Modify:**
- `bet_logger/README.md` — add a short section describing the MLB Summary tab and how to run the setup script.

**Touch:** nothing else. No changes to `sheets.py`, `scraper_wagerzon.py`, `create_summary.py`, `run_all_scrapers.sh`, or any DuckDB layer.

---

## Version Control

- **Branch:** `feature/mlb-summary-tab` (off `main`)
- **Worktree:** `../NFLWork-mlb-summary` (to avoid conflicts with in-flight work in the main tree)
- **Commits:** one per task where something shippable is produced (usually after a passing test + working code)
- **Merge:** after the script has been run successfully against the live sheet and the user has eyeballed the resulting tab; user approval required before merge to `main`; delete branch and worktree after merge

---

## Task 1: Create worktree and feature branch

**Files:** none yet — this is just environment setup.

- [ ] **Step 1: Verify on main, clean tree**

Run: `git -C /Users/callancapitolo/NFLWork status && git -C /Users/callancapitolo/NFLWork branch --show-current`
Expected: clean working tree, current branch = `main`

- [ ] **Step 2: Create worktree**

Run: `git -C /Users/callancapitolo/NFLWork worktree add -b feature/mlb-summary-tab ../NFLWork-mlb-summary main`
Expected: `Preparing worktree (new branch 'feature/mlb-summary-tab')` then `HEAD is now at ...`

- [ ] **Step 3: Confirm worktree is live**

Run: `ls /Users/callancapitolo/NFLWork-mlb-summary && git -C /Users/callancapitolo/NFLWork-mlb-summary branch --show-current`
Expected: directory contents identical to main repo; branch = `feature/mlb-summary-tab`

All subsequent tasks happen inside `/Users/callancapitolo/NFLWork-mlb-summary`.

---

## Task 2: Pure-Python filter helpers + tests (TDD)

We test the regex logic because it is the one thing most likely to break silently. Everything else in the script is formula strings that can only be validated by running against the live sheet.

**Files:**
- Create: `bet_logger/test_mlb_summary.py`
- Create: `bet_logger/create_mlb_summary.py` (initial skeleton — constants + helper functions only)

- [ ] **Step 1: Write the failing test file**

Create `bet_logger/test_mlb_summary.py`:

```python
"""Unit tests for create_mlb_summary.py filter helpers."""
import pytest
from create_mlb_summary import (
    is_mlb_correlated_parlay,
    is_f5_parlay,
    SPREAD_REGEX,
    TOTAL_REGEX,
    F5_REGEX,
)


# Real Wagerzon description samples (format: "PARLAY (2 TEAMS) | leg1 | leg2")
SPREAD_PLUS_TOTAL_FG = (
    "PARLAY (2 TEAMS) | Detroit Tigers -1.5 +140 | Over 8.5 -110"
)
SPREAD_PLUS_TOTAL_F5 = (
    "PARLAY (2 TEAMS) | Detroit Tigers 1st 5 Innings -0.5 +115 | "
    "1st 5 Innings Over 4.5 -105"
)
ML_PLUS_ML = (
    "PARLAY (2 TEAMS) | Detroit Tigers ML +100 | Chicago Cubs ML -120"
)
TOTAL_PLUS_TOTAL = (
    "PARLAY (2 TEAMS) | Over 8.5 -110 | Over 4.5 -105 1st 5 Innings"
)
THREE_LEGS = (
    "PARLAY (3 TEAMS) | Detroit Tigers -1.5 +140 | Over 8.5 -110 | "
    "Yankees ML -200"
)
STRAIGHT_BET = "STRAIGHT BET | Detroit Tigers -1.5 +140"


@pytest.mark.parametrize(
    "desc, sport, bet_type, expected",
    [
        (SPREAD_PLUS_TOTAL_FG, "MLB", "Parlay", True),
        (SPREAD_PLUS_TOTAL_F5, "MLB", "Parlay", True),
        (ML_PLUS_ML, "MLB", "Parlay", False),
        (TOTAL_PLUS_TOTAL, "MLB", "Parlay", False),
        (THREE_LEGS, "MLB", "Parlay", False),
        (STRAIGHT_BET, "MLB", "Straight", False),
        (SPREAD_PLUS_TOTAL_FG, "NFL", "Parlay", False),
    ],
)
def test_is_mlb_correlated_parlay(desc, sport, bet_type, expected):
    assert is_mlb_correlated_parlay(desc, sport, bet_type) is expected


@pytest.mark.parametrize(
    "desc, expected",
    [
        (SPREAD_PLUS_TOTAL_FG, False),
        (SPREAD_PLUS_TOTAL_F5, True),
        ("PARLAY (2 TEAMS) | Yankees F5 -0.5 +110 | Over 4.5 -105", True),
        ("PARLAY (2 TEAMS) | Yankees First 5 -0.5 +110 | Over 4.5 -105", True),
        ("PARLAY (2 TEAMS) | Yankees -1.5 +140 | Over 8.5 -110", False),
    ],
)
def test_is_f5_parlay(desc, expected):
    assert is_f5_parlay(desc) is expected
```

- [ ] **Step 2: Run the tests and verify they fail**

Run: `cd /Users/callancapitolo/NFLWork-mlb-summary/bet_logger && ./venv/bin/python3 -m pytest test_mlb_summary.py -v`
Expected: `ImportError` (no `create_mlb_summary` module yet).

- [ ] **Step 3: Create the skeleton with the helpers**

Create `bet_logger/create_mlb_summary.py`:

```python
#!/usr/bin/env python3
"""
Create the MLB Summary tab in the bet-logging Google Sheet.

One-shot setup script. Writes pure spreadsheet formulas that auto-update
from Sheet1 as new Wagerzon bets are scraped in. Also creates two embedded
charts (equity curve + weekly P&L bars). No helper sheet. No cron.

Filter: sport=MLB + bet_type=Parlay + description starts with
"PARLAY (2 TEAMS)" + description has a spread token + has a total token.
"""

import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from sheets import (  # noqa: E402
    get_sheets_service,
    SPREADSHEET_ID,
    SHEET_NAME,
    ensure_sheet_exists,
)
from googleapiclient.errors import HttpError  # noqa: E402


# ── Constants ─────────────────────────────────────────────────────────────
TAB = "MLB Summary"
SRC = SHEET_NAME  # "Sheet1"
N = 10000         # Max row the formulas scan in Sheet1


# ── Filter regex ──────────────────────────────────────────────────────────
# A spread leg contains "<team-name> <signed decimal>". The signed decimal
# must be the spread, not the American odds, so we require a decimal point
# OR a bare "+/- 0.5 / 1.5" — the shapes Wagerzon emits for MLB spreads.
SPREAD_REGEX = re.compile(r"(?<!\d)[+\-](?:0\.5|1\.5|2\.5)(?!\d)")

# A total leg contains "Over N.N" or "Under N.N". We accept the word forms
# only (the bare "O 8.5" / "U 8.5" shorthand does not appear in Wagerzon's
# MLB output — confirmed against the scraper's clean_desc output).
TOTAL_REGEX = re.compile(r"\b(Over|Under)\s+\d+(?:\.\d+)?\b", re.IGNORECASE)

# F5 markers in leg descriptions.
F5_REGEX = re.compile(r"\b(?:1st\s*5|F5|First\s*5)\b", re.IGNORECASE)


def is_mlb_correlated_parlay(description: str, sport: str, bet_type: str) -> bool:
    """True iff this Sheet1 row is a tracked MLB correlated parlay."""
    if sport != "MLB":
        return False
    if bet_type != "Parlay":
        return False
    if not description.startswith("PARLAY (2 TEAMS)"):
        return False
    if not SPREAD_REGEX.search(description):
        return False
    if not TOTAL_REGEX.search(description):
        return False
    return True


def is_f5_parlay(description: str) -> bool:
    """True iff this parlay is a First-5-Innings parlay."""
    return bool(F5_REGEX.search(description))


def main():
    # Filled in by later tasks
    raise NotImplementedError("main() not implemented yet — see later tasks")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests and verify they pass**

Run: `cd /Users/callancapitolo/NFLWork-mlb-summary/bet_logger && ./venv/bin/python3 -m pytest test_mlb_summary.py -v`
Expected: all 12 tests pass.

- [ ] **Step 5: Commit**

```bash
cd /Users/callancapitolo/NFLWork-mlb-summary
git add bet_logger/create_mlb_summary.py bet_logger/test_mlb_summary.py
git commit -m "feat(bet_logger): add MLB correlated-parlay filter helpers with tests"
```

---

## Task 3: Formula builders (filter mask + metric formulas)

Now write the Google Sheets-formula helpers that mirror the Python filter. These are pure string builders. No network, no tests — verified by visually inspecting the final tab in Task 9.

**Files:**
- Modify: `bet_logger/create_mlb_summary.py`

- [ ] **Step 1: Add the formula builders**

Append to `bet_logger/create_mlb_summary.py`, after the regex constants but before `main()`:

```python
# ── Formula builders ──────────────────────────────────────────────────────
# Column map in Sheet1:
#   A=Date  B=Platform  C=Sport  D=Description  E=Bet Type
#   F=Line  G=Odds      H=Stake  I=Dec Odds     J=Result

S = f"'{SRC}'"
COL_DATE = f"{S}!A$2:A${N}"
COL_SPORT = f"{S}!C$2:C${N}"
COL_DESC = f"{S}!D$2:D${N}"
COL_TYPE = f"{S}!E$2:E${N}"
COL_STAKE = f"{S}!H$2:H${N}"
COL_DEC = f"{S}!I$2:I${N}"
COL_RESULT = f"{S}!J$2:J${N}"

# Filter mask as a Sheets expression (no leading `=`). Five conditions,
# multiplied together (boolean AND). Leaves a 1/0 array suitable for
# SUMPRODUCT. Mirrors is_mlb_correlated_parlay() exactly.
FILTER_MASK = (
    f'({COL_SPORT}="MLB")'
    f'*({COL_TYPE}="Parlay")'
    f'*(LEFT({COL_DESC},16)="PARLAY (2 TEAMS)")'
    f'*IFERROR(REGEXMATCH({COL_DESC},'
    r'"(?<!\d)[+\-](?:0\.5|1\.5|2\.5)(?!\d)"'
    f'),FALSE)'
    f'*IFERROR(REGEXMATCH({COL_DESC},'
    r'"\b(?i)(Over|Under)\s+\d+(?:\.\d+)?\b"'
    f'),FALSE)'
)

# The FG/F5 modifiers — multiplied INTO the filter mask to split the set.
F5_COND = (
    f'IFERROR(REGEXMATCH({COL_DESC},'
    r'"\b(?i)(?:1st\s*5|F5|First\s*5)\b"'
    f'),FALSE)'
)
FG_ADD = f"*(1-{F5_COND})"   # NOT F5
F5_ADD = f"*{F5_COND}"       # IS F5

SETTLED = f'({COL_RESULT}<>"")'
WIN = f'({COL_RESULT}="win")'
LOSS = f'({COL_RESULT}="loss")'
PUSH = f'({COL_RESULT}="push")'


def placed_f(mask=FILTER_MASK, extra=""):
    """Count of qualifying bets (settled + pending)."""
    return f"=SUMPRODUCT({mask}{extra})"


def settled_f(mask=FILTER_MASK, extra=""):
    """Count of qualifying settled bets."""
    return f"=SUMPRODUCT({mask}{extra}*{SETTLED})"


def wagered_f(mask=FILTER_MASK, extra=""):
    """Sum of stake across settled qualifying bets."""
    return f"=SUMPRODUCT({mask}{extra}*{SETTLED}*{COL_STAKE})"


def _profit_expr(mask=FILTER_MASK, extra=""):
    """Profit expression WITHOUT leading `=` (for embedding)."""
    wins_payout = (
        f"SUMPRODUCT({mask}{extra}*{WIN}*({COL_DEC}<>\"\")"
        f"*{COL_STAKE}*({COL_DEC}-1))"
    )
    losses = f"SUMPRODUCT({mask}{extra}*{LOSS}*{COL_STAKE})"
    return f"{wins_payout}-{losses}"


def profit_f(mask=FILTER_MASK, extra=""):
    return f"={_profit_expr(mask, extra)}"


def wins_f(mask=FILTER_MASK, extra=""):
    return f"=SUMPRODUCT({mask}{extra}*{WIN})"


def losses_f(mask=FILTER_MASK, extra=""):
    return f"=SUMPRODUCT({mask}{extra}*{LOSS})"


def pushes_f(mask=FILTER_MASK, extra=""):
    return f"=SUMPRODUCT({mask}{extra}*{PUSH})"


def record_f(mask=FILTER_MASK, extra=""):
    w = f"SUMPRODUCT({mask}{extra}*{WIN})"
    l = f"SUMPRODUCT({mask}{extra}*{LOSS})"
    p = f"SUMPRODUCT({mask}{extra}*{PUSH})"
    return f'={w}&"-"&{l}&"-"&{p}'


def avg_odds_f(mask=FILTER_MASK, extra=""):
    num = f"SUMPRODUCT({mask}{extra}*{SETTLED}*{COL_DEC})"
    den = f"SUMPRODUCT({mask}{extra}*{SETTLED})"
    return f"=IF({den}=0,0,{num}/{den})"
```

- [ ] **Step 2: Confirm the file still imports cleanly**

Run: `cd /Users/callancapitolo/NFLWork-mlb-summary/bet_logger && ./venv/bin/python3 -c "import create_mlb_summary; print('OK')"`
Expected: `OK`

- [ ] **Step 3: Commit**

```bash
cd /Users/callancapitolo/NFLWork-mlb-summary
git add bet_logger/create_mlb_summary.py
git commit -m "feat(bet_logger): add MLB Summary formula builders (filter mask, metrics)"
```

---

## Task 4: `build_rows()` — Block 1 (Overall) + Block 2 (FG vs F5 split)

**Files:**
- Modify: `bet_logger/create_mlb_summary.py`

- [ ] **Step 1: Add `build_rows()` covering Blocks 1 & 2**

Insert before `main()`:

```python
# ── Row builder ───────────────────────────────────────────────────────────

def build_rows():
    """Return the 2D list of cell values for the MLB Summary tab.

    Layout (row indices 1-based):
      1  Title
      2  (blank)
      3  "OVERALL STATS" header
      4  column labels
      5-13 overall metrics (anchor row numbers used by ROI/win-rate formulas)
      14 (blank)
      15 "BY GAME WINDOW" header
      16 column labels
      17 FG row
      18 F5 row
      19 (blank)
      20 "EQUITY CURVE (Daily P&L)" header
      21 column labels ("Date","Bets","Wagered","P&L","Cumulative P&L")
      22+ daily data rows (dynamic, auto-extend)
      ...
      (a later section for weekly — filled by Task 5)
    """
    rows = []

    # Row 1: title
    rows.append(["MLB CORRELATED PARLAYS — WAGERZON"])
    # Row 2: blank
    rows.append([])

    # Row 3: section header
    rows.append(["OVERALL STATS"])
    # Row 4: column labels
    rows.append(["", "Value"])
    # Row 5: Total bets placed
    rows.append(["Total bets placed", placed_f()])
    # Row 6: Settled bets
    rows.append(["Settled bets (W+L+P)", settled_f()])
    # Row 7: Total wagered
    rows.append(["Total wagered", wagered_f()])
    # Row 8: Net P&L
    rows.append(["Net P&L", profit_f()])
    # Row 9: ROI = P&L / wagered  (references anchors B7, B8)
    rows.append(["ROI", "=IF(B7=0,0,B8/B7)"])
    # Row 10: Win rate = wins / (wins + losses), NOT wins / settled
    w_expr = f"SUMPRODUCT({FILTER_MASK}*{WIN})"
    l_expr = f"SUMPRODUCT({FILTER_MASK}*{LOSS})"
    rows.append([
        "Win rate",
        f"=IF(({w_expr}+{l_expr})=0,0,{w_expr}/({w_expr}+{l_expr}))",
    ])
    # Row 11: Record
    rows.append(["Record (W-L-P)", record_f()])
    # Row 12: Avg decimal odds
    rows.append(["Avg decimal odds", avg_odds_f()])
    # Row 13: blank spacer
    rows.append([])

    # Row 14: "BY GAME WINDOW" header
    rows.append(["BY GAME WINDOW"])
    # Row 15: column labels
    rows.append([
        "Window", "Bets", "Settled", "Wagered",
        "P&L", "ROI", "Win Rate", "Record", "Avg Odds",
    ])
    # Rows 16-17: FG and F5 rows
    for label, extra in (("Full Game (FG)", FG_ADD), ("First 5 (F5)", F5_ADD)):
        r = len(rows) + 1
        w_here = f"SUMPRODUCT({FILTER_MASK}{extra}*{WIN})"
        l_here = f"SUMPRODUCT({FILTER_MASK}{extra}*{LOSS})"
        rows.append([
            label,
            placed_f(extra=extra),
            settled_f(extra=extra),
            wagered_f(extra=extra),
            profit_f(extra=extra),
            f"=IF(D{r}=0,0,E{r}/D{r})",     # ROI
            f"=IF(({w_here}+{l_here})=0,0,{w_here}/({w_here}+{l_here}))",
            record_f(extra=extra),
            avg_odds_f(extra=extra),
        ])

    # Row 18: blank spacer (Task 5 continues from here)
    rows.append([])

    return rows
```

- [ ] **Step 2: Sanity-check that `build_rows()` produces the expected shape**

Run:
```
cd /Users/callancapitolo/NFLWork-mlb-summary/bet_logger && \
./venv/bin/python3 -c "import create_mlb_summary as m; r = m.build_rows(); print(len(r), 'rows'); print(r[4]); print(r[15]); print(r[16])"
```
Expected: `18 rows`, then the "Total bets placed" row, then the column-labels row, then the FG row — each visibly showing formulas prefixed with `=`.

- [ ] **Step 3: Commit**

```bash
cd /Users/callancapitolo/NFLWork-mlb-summary
git add bet_logger/create_mlb_summary.py
git commit -m "feat(bet_logger): add Overall + FG/F5 blocks for MLB Summary tab"
```

---

## Task 5: Block 3 — dynamic daily + weekly tables

Unlike the CBB tab, which pre-populates every date from season start to April 15, we let the sheet generate the date list dynamically via `SORT(UNIQUE(FILTER(...)))`. The table auto-extends whenever a bet on a new date lands in `Sheet1`.

**Files:**
- Modify: `bet_logger/create_mlb_summary.py`

- [ ] **Step 1: Extend `build_rows()` with Blocks 3a (daily) and 3b (weekly)**

Append to `build_rows()` — inside the function, before the final `return rows`:

```python
    # ── Block 3a: Daily P&L (dynamic, feeds the equity curve chart) ──
    rows.append(["DAILY P&L"])
    rows.append(["Date", "Bets", "Wagered", "P&L", "Cumulative P&L"])
    daily_header_row = len(rows)          # 1-based row of column labels
    daily_start_row = daily_header_row + 1

    # Column A: a single dynamic formula. The spilled dates fill downward.
    # We put the formula in the top cell; rows below it will be filled by
    # the ARRAYFORMULA spill when the sheet is opened.
    rows.append([
        (
            f"=IFERROR(SORT(UNIQUE(FILTER({COL_DATE},{FILTER_MASK}=1))),\"\")"
        ),
        # The other four columns use ARRAYFORMULA over the spilled date column
        # in A — they must reference A{start}:A (whole column from start down).
        f"=ARRAYFORMULA(IF(A{daily_start_row}:A=\"\",\"\","
        f"SUMPRODUCT(({COL_DATE}=A{daily_start_row}:A)*{FILTER_MASK}*{SETTLED})))",
        f"=ARRAYFORMULA(IF(A{daily_start_row}:A=\"\",\"\","
        f"SUMPRODUCT(({COL_DATE}=A{daily_start_row}:A)*{FILTER_MASK}*{SETTLED}*{COL_STAKE})))",
        # P&L per day — wins payout minus losses stake
        (
            f"=ARRAYFORMULA(IF(A{daily_start_row}:A=\"\",\"\","
            f"SUMPRODUCT(({COL_DATE}=A{daily_start_row}:A)*{FILTER_MASK}*{WIN}"
            f"*({COL_DEC}<>\"\")*{COL_STAKE}*({COL_DEC}-1))"
            f"-SUMPRODUCT(({COL_DATE}=A{daily_start_row}:A)*{FILTER_MASK}*{LOSS}*{COL_STAKE})))"
        ),
        # Cumulative P&L — running sum of column D starting from daily_start_row.
        # ARRAYFORMULA of MMULT gives a running sum without per-row formulas.
        (
            f"=ARRAYFORMULA(IF(A{daily_start_row}:A=\"\",\"\","
            f"MMULT(--(ROW(D{daily_start_row}:D)>=TRANSPOSE(ROW(D{daily_start_row}:D))),"
            f"IFERROR(D{daily_start_row}:D*1,0))))"
        ),
    ])
    rows.append([])  # spacer after daily block

    # ── Block 3b: Weekly P&L (dynamic, feeds the bar chart) ──
    rows.append(["WEEKLY P&L"])
    rows.append(["Week of (Mon)", "Bets", "Wagered", "P&L", "Cumulative P&L"])
    weekly_header_row = len(rows)
    weekly_start_row = weekly_header_row + 1

    # Week-start bucket = A - WEEKDAY(A, 2) + 1  (Monday of that date's week).
    # Wrap the FILTER in ARRAYFORMULA so the transform is applied element-wise.
    week_of_date = (
        f"({COL_DATE}-WEEKDAY({COL_DATE},2)+1)"
    )
    rows.append([
        (
            f"=IFERROR(SORT(UNIQUE("
            f"FILTER(ARRAYFORMULA({week_of_date}),{FILTER_MASK}=1))),\"\")"
        ),
        f"=ARRAYFORMULA(IF(A{weekly_start_row}:A=\"\",\"\","
        f"SUMPRODUCT(({week_of_date}=A{weekly_start_row}:A)*{FILTER_MASK}*{SETTLED})))",
        f"=ARRAYFORMULA(IF(A{weekly_start_row}:A=\"\",\"\","
        f"SUMPRODUCT(({week_of_date}=A{weekly_start_row}:A)*{FILTER_MASK}*{SETTLED}*{COL_STAKE})))",
        (
            f"=ARRAYFORMULA(IF(A{weekly_start_row}:A=\"\",\"\","
            f"SUMPRODUCT(({week_of_date}=A{weekly_start_row}:A)*{FILTER_MASK}*{WIN}"
            f"*({COL_DEC}<>\"\")*{COL_STAKE}*({COL_DEC}-1))"
            f"-SUMPRODUCT(({week_of_date}=A{weekly_start_row}:A)*{FILTER_MASK}*{LOSS}*{COL_STAKE})))"
        ),
        (
            f"=ARRAYFORMULA(IF(A{weekly_start_row}:A=\"\",\"\","
            f"MMULT(--(ROW(D{weekly_start_row}:D)>=TRANSPOSE(ROW(D{weekly_start_row}:D))),"
            f"IFERROR(D{weekly_start_row}:D*1,0))))"
        ),
    ])

    # Also stash the anchor row numbers so Task 6 (charts) can reference
    # them without re-computing. Attach as an attribute on the returned list.
    rows.__anchor_daily_start__ = daily_start_row          # type: ignore[attr-defined]
    rows.__anchor_weekly_start__ = weekly_start_row        # type: ignore[attr-defined]

    return rows
```

(Note: the trailing `return rows` already in Task 4's version should be *replaced* by the version above — keep the earlier block content, just extend before the return. Re-open the file and confirm only one `return rows` exists.)

- [ ] **Step 2: Sanity-check the output**

Run:
```
cd /Users/callancapitolo/NFLWork-mlb-summary/bet_logger && \
./venv/bin/python3 -c "import create_mlb_summary as m; r = m.build_rows(); print(len(r), 'rows'); print('daily_start:', r.__anchor_daily_start__); print('weekly_start:', r.__anchor_weekly_start__); print('Daily header row content:', r[r.__anchor_daily_start__ - 2])"
```
Expected: rows count ≥ 25; `daily_start` and `weekly_start` printed as integers; the daily header row is `['Date', 'Bets', 'Wagered', 'P&L', 'Cumulative P&L']`.

- [ ] **Step 3: Commit**

```bash
cd /Users/callancapitolo/NFLWork-mlb-summary
git add bet_logger/create_mlb_summary.py
git commit -m "feat(bet_logger): add dynamic daily and weekly P&L blocks"
```

---

## Task 6: Chart creation (equity curve + weekly bars)

The Google Sheets API's `addChart` request embeds a chart tied to a cell range. When the range's values change, the chart re-renders automatically.

**Files:**
- Modify: `bet_logger/create_mlb_summary.py`

- [ ] **Step 1: Add `get_sheet_id()` helper and `build_chart_requests()`**

Append to `bet_logger/create_mlb_summary.py`, before `main()`:

```python
# ── Chart builders ────────────────────────────────────────────────────────

def get_sheet_id(service, tab_name):
    meta = service.spreadsheets().get(spreadsheetId=SPREADSHEET_ID).execute()
    for s in meta.get("sheets", []):
        if s["properties"]["title"] == tab_name:
            return s["properties"]["sheetId"]
    return None


def build_chart_requests(sheet_id, daily_start_row, weekly_start_row, last_row):
    """Return Sheets API batchUpdate requests for the two charts.

    Both charts read a generous range (through `last_row`) so they pick up
    new dates automatically as the dynamic tables extend downward.
    """
    # Convert 1-based row numbers to API's 0-based indices.
    # The daily data block occupies columns A-E (indices 0-5 exclusive).
    daily_header_idx = daily_start_row - 2       # 0-based index of the header row
    daily_data_start = daily_start_row - 1       # 0-based index of first data row

    weekly_header_idx = weekly_start_row - 2
    weekly_data_start = weekly_start_row - 1

    # EQUITY CURVE (line chart) — X: Date (col A), Y: Cumulative P&L (col E)
    equity_chart = {
        "addChart": {
            "chart": {
                "spec": {
                    "title": "Equity Curve — Cumulative P&L",
                    "basicChart": {
                        "chartType": "LINE",
                        "legendPosition": "NO_LEGEND",
                        "axis": [
                            {"position": "BOTTOM_AXIS", "title": "Date"},
                            {"position": "LEFT_AXIS", "title": "Cumulative P&L"},
                        ],
                        "domains": [{
                            "domain": {"sourceRange": {"sources": [{
                                "sheetId": sheet_id,
                                "startRowIndex": daily_data_start,
                                "endRowIndex": last_row,
                                "startColumnIndex": 0, "endColumnIndex": 1,
                            }]}}
                        }],
                        "series": [{
                            "series": {"sourceRange": {"sources": [{
                                "sheetId": sheet_id,
                                "startRowIndex": daily_data_start,
                                "endRowIndex": last_row,
                                "startColumnIndex": 4, "endColumnIndex": 5,
                            }]}},
                            "targetAxis": "LEFT_AXIS",
                        }],
                        "headerCount": 0,
                    },
                },
                "position": {
                    "overlayPosition": {
                        "anchorCell": {
                            "sheetId": sheet_id,
                            "rowIndex": daily_header_idx - 12,  # above DAILY P&L block
                            "columnIndex": 6,                   # column G
                        },
                        "widthPixels": 600,
                        "heightPixels": 300,
                    }
                },
            }
        }
    }

    # WEEKLY BARS — X: Week of (col A), Y: P&L (col D)
    weekly_chart = {
        "addChart": {
            "chart": {
                "spec": {
                    "title": "Weekly P&L",
                    "basicChart": {
                        "chartType": "COLUMN",
                        "legendPosition": "NO_LEGEND",
                        "axis": [
                            {"position": "BOTTOM_AXIS", "title": "Week of"},
                            {"position": "LEFT_AXIS", "title": "P&L"},
                        ],
                        "domains": [{
                            "domain": {"sourceRange": {"sources": [{
                                "sheetId": sheet_id,
                                "startRowIndex": weekly_data_start,
                                "endRowIndex": last_row,
                                "startColumnIndex": 0, "endColumnIndex": 1,
                            }]}}
                        }],
                        "series": [{
                            "series": {"sourceRange": {"sources": [{
                                "sheetId": sheet_id,
                                "startRowIndex": weekly_data_start,
                                "endRowIndex": last_row,
                                "startColumnIndex": 3, "endColumnIndex": 4,
                            }]}},
                            "targetAxis": "LEFT_AXIS",
                        }],
                        "headerCount": 0,
                    },
                },
                "position": {
                    "overlayPosition": {
                        "anchorCell": {
                            "sheetId": sheet_id,
                            "rowIndex": weekly_header_idx - 6,
                            "columnIndex": 6,                   # column G, below equity curve
                        },
                        "widthPixels": 600,
                        "heightPixels": 300,
                    }
                },
            }
        }
    }

    return [equity_chart, weekly_chart]
```

- [ ] **Step 2: Import sanity check**

Run: `cd /Users/callancapitolo/NFLWork-mlb-summary/bet_logger && ./venv/bin/python3 -c "import create_mlb_summary; print(create_mlb_summary.build_chart_requests(0, 22, 30, 1000))" | head -c 400`
Expected: a truncated printed list containing two dicts, each with an `addChart` key.

- [ ] **Step 3: Commit**

```bash
cd /Users/callancapitolo/NFLWork-mlb-summary
git add bet_logger/create_mlb_summary.py
git commit -m "feat(bet_logger): add chart requests for equity curve + weekly bars"
```

---

## Task 7: Formatting (currency, percent, conditional colors)

**Files:**
- Modify: `bet_logger/create_mlb_summary.py`

- [ ] **Step 1: Add formatting request builders**

Append to `bet_logger/create_mlb_summary.py`, before `main()`:

```python
# ── Formatting ────────────────────────────────────────────────────────────

def _range(sid, r1, r2, c1, c2):
    return {"sheetId": sid, "startRowIndex": r1, "endRowIndex": r2,
            "startColumnIndex": c1, "endColumnIndex": c2}


def _text_fmt(sid, r1, r2, c1, c2, bold=False, size=None):
    tf = {"bold": bold}
    if size:
        tf["fontSize"] = size
    return {"repeatCell": {
        "range": _range(sid, r1, r2, c1, c2),
        "cell": {"userEnteredFormat": {"textFormat": tf}},
        "fields": "userEnteredFormat.textFormat",
    }}


def _bg_fmt(sid, r1, r2, c1, c2, rgb):
    return {"repeatCell": {
        "range": _range(sid, r1, r2, c1, c2),
        "cell": {"userEnteredFormat": {"backgroundColor":
            {"red": rgb[0], "green": rgb[1], "blue": rgb[2]}}},
        "fields": "userEnteredFormat.backgroundColor",
    }}


def _num_fmt(sid, r1, r2, c1, c2, ntype, pattern):
    return {"repeatCell": {
        "range": _range(sid, r1, r2, c1, c2),
        "cell": {"userEnteredFormat": {"numberFormat":
            {"type": ntype, "pattern": pattern}}},
        "fields": "userEnteredFormat.numberFormat",
    }}


def _cond_fmt(sid, r1, r2, c1, c2, cond_type, value, rgb):
    return {"addConditionalFormatRule": {"rule": {
        "ranges": [_range(sid, r1, r2, c1, c2)],
        "booleanRule": {
            "condition": {"type": cond_type,
                          "values": [{"userEnteredValue": value}]},
            "format": {"textFormat": {"foregroundColorStyle":
                {"rgbColor": {"red": rgb[0], "green": rgb[1], "blue": rgb[2]}}}},
        },
    }, "index": 0}}


def build_format_requests(sid, daily_start_row, weekly_start_row):
    """Return the batchUpdate requests to format the MLB Summary tab."""
    reqs = []
    GREEN = (0.13, 0.55, 0.13)
    RED = (0.80, 0.13, 0.13)
    HEADER_BG = (0.85, 0.92, 1.00)
    LABEL_BG = (0.93, 0.93, 0.93)

    # Title (row 1)
    reqs.append(_text_fmt(sid, 0, 1, 0, 10, bold=True, size=16))

    # Section headers: rows 3 (OVERALL), 14 (BY GAME WINDOW),
    # (daily_header_row - 1), and (weekly_header_row - 1).
    for section_row_1based in (3, 14, daily_start_row - 1, weekly_start_row - 1):
        r = section_row_1based - 1
        reqs.append(_bg_fmt(sid, r, r + 1, 0, 10, HEADER_BG))
        reqs.append(_text_fmt(sid, r, r + 1, 0, 10, bold=True, size=12))

    # Column-label rows (the one right below each section header)
    for label_row_1based in (4, 15, daily_start_row, weekly_start_row):
        r = label_row_1based - 1
        reqs.append(_bg_fmt(sid, r, r + 1, 0, 10, LABEL_BG))
        reqs.append(_text_fmt(sid, r, r + 1, 0, 10, bold=True))

    # OVERALL STATS — rows 5-12 (0-based 4-12), column B (index 1)
    #   Row indices of: 5=placed 6=settled 7=wagered 8=pnl 9=roi 10=winrate 11=record 12=avgodds
    reqs.append(_num_fmt(sid, 6, 8, 1, 2, "CURRENCY", "$#,##0.00"))   # wagered, P&L
    reqs.append(_num_fmt(sid, 7, 8, 1, 2, "CURRENCY", "$#,##0.00"))   # P&L (overlaps intentionally)
    reqs.append(_num_fmt(sid, 8, 9, 1, 2, "PERCENT", "0.00%"))        # ROI
    reqs.append(_num_fmt(sid, 9, 10, 1, 2, "PERCENT", "0.0%"))        # Win Rate
    reqs.append(_num_fmt(sid, 11, 12, 1, 2, "NUMBER", "0.00"))        # Avg odds

    # FG / F5 rows — rows 16-17 (0-based 15-16)
    reqs.append(_num_fmt(sid, 15, 17, 3, 5, "CURRENCY", "$#,##0.00"))  # wagered, P&L
    reqs.append(_num_fmt(sid, 15, 17, 5, 6, "PERCENT", "0.00%"))       # ROI
    reqs.append(_num_fmt(sid, 15, 17, 6, 7, "PERCENT", "0.0%"))        # Win Rate
    reqs.append(_num_fmt(sid, 15, 17, 8, 9, "NUMBER", "0.00"))         # Avg odds
    reqs.append(_cond_fmt(sid, 15, 17, 4, 5, "NUMBER_GREATER", "0", GREEN))
    reqs.append(_cond_fmt(sid, 15, 17, 4, 5, "NUMBER_LESS", "0", RED))

    # Daily and weekly blocks — generous range (1000 rows) for dynamic data
    for start_1based in (daily_start_row, weekly_start_row):
        s = start_1based - 1
        e = s + 1000
        reqs.append(_num_fmt(sid, s, e, 2, 5, "CURRENCY", "$#,##0.00"))  # wagered, P&L, cum
        reqs.append(_cond_fmt(sid, s, e, 3, 4, "NUMBER_GREATER", "0", GREEN))
        reqs.append(_cond_fmt(sid, s, e, 3, 4, "NUMBER_LESS", "0", RED))
        reqs.append(_cond_fmt(sid, s, e, 4, 5, "NUMBER_GREATER", "0", GREEN))
        reqs.append(_cond_fmt(sid, s, e, 4, 5, "NUMBER_LESS", "0", RED))

    # Column widths (A-I)
    for i, w in enumerate([155, 80, 85, 110, 115, 75, 85, 130, 90]):
        reqs.append({"updateDimensionProperties": {
            "range": {"sheetId": sid, "dimension": "COLUMNS",
                      "startIndex": i, "endIndex": i + 1},
            "properties": {"pixelSize": w}, "fields": "pixelSize",
        }})

    return reqs
```

- [ ] **Step 2: Commit**

```bash
cd /Users/callancapitolo/NFLWork-mlb-summary
git add bet_logger/create_mlb_summary.py
git commit -m "feat(bet_logger): add formatting requests for MLB Summary tab"
```

---

## Task 8: Wire everything together in `main()`

**Files:**
- Modify: `bet_logger/create_mlb_summary.py`

- [ ] **Step 1: Replace the placeholder `main()`**

Find the stub `def main(): raise NotImplementedError...` and replace with:

```python
# ── Main ──────────────────────────────────────────────────────────────────

def main():
    print("Setting up MLB Summary tab...")
    service = get_sheets_service()
    ensure_sheet_exists(service, TAB)

    # Clear any prior contents
    try:
        service.spreadsheets().values().clear(
            spreadsheetId=SPREADSHEET_ID,
            range=f"'{TAB}'!A:Z",
        ).execute()
    except HttpError as e:
        print(f"  (clear failed — continuing: {e})")

    # Build rows
    rows = build_rows()
    daily_start = rows.__anchor_daily_start__          # type: ignore[attr-defined]
    weekly_start = rows.__anchor_weekly_start__        # type: ignore[attr-defined]

    # Pad rows to 9 columns so the update range is rectangular
    max_cols = max((len(r) for r in rows), default=1)
    padded = [list(r) + [""] * (max_cols - len(r)) for r in rows]

    print(f"Writing {len(padded)} rows ({max_cols} cols)...")
    service.spreadsheets().values().update(
        spreadsheetId=SPREADSHEET_ID,
        range=f"'{TAB}'!A1",
        valueInputOption="USER_ENTERED",
        body={"values": padded},
    ).execute()

    # Apply formatting and add charts in one batchUpdate
    sid = get_sheet_id(service, TAB)
    if sid is None:
        print("  (could not resolve sheet id — skipping formatting/charts)")
        return

    requests = []
    # Remove any existing charts on this tab first (in case of re-run)
    meta = service.spreadsheets().get(spreadsheetId=SPREADSHEET_ID).execute()
    for s in meta.get("sheets", []):
        if s["properties"]["sheetId"] == sid:
            for ch in s.get("charts", []):
                requests.append({"deleteEmbeddedObject": {"objectId": ch["chartId"]}})
            break

    requests.extend(build_format_requests(sid, daily_start, weekly_start))
    requests.extend(build_chart_requests(sid, daily_start, weekly_start, last_row=5000))

    service.spreadsheets().batchUpdate(
        spreadsheetId=SPREADSHEET_ID,
        body={"requests": requests},
    ).execute()

    print(f"Done. Open the sheet and look at the '{TAB}' tab.")
```

- [ ] **Step 2: Import sanity check (no network yet)**

Run: `cd /Users/callancapitolo/NFLWork-mlb-summary/bet_logger && ./venv/bin/python3 -c "import create_mlb_summary; print('main is', create_mlb_summary.main)"`
Expected: `main is <function main at 0x...>` (no ImportError).

- [ ] **Step 3: Run the tests one more time to make sure nothing regressed**

Run: `cd /Users/callancapitolo/NFLWork-mlb-summary/bet_logger && ./venv/bin/python3 -m pytest test_mlb_summary.py -v`
Expected: 12 passing.

- [ ] **Step 4: Commit**

```bash
cd /Users/callancapitolo/NFLWork-mlb-summary
git add bet_logger/create_mlb_summary.py
git commit -m "feat(bet_logger): wire up create_mlb_summary main() with charts + formatting"
```

---

## Task 9: Run against the live sheet and verify

This is where we find out whether the formulas, regex, and charts actually work. Visual verification — open the sheet and look.

**Files:** none changed.

- [ ] **Step 1: Run the script**

Run: `cd /Users/callancapitolo/NFLWork-mlb-summary/bet_logger && ./venv/bin/python3 create_mlb_summary.py`
Expected output includes:
```
Setting up MLB Summary tab...
Writing N rows (M cols)...
Done. Open the sheet and look at the 'MLB Summary' tab.
```
No tracebacks. No HttpError.

- [ ] **Step 2: Open the Google Sheet and visually verify**

Open: `https://docs.google.com/spreadsheets/d/1t9_7HmsrQAu34HI_gMIDPErC2kThz2MHBwwzZIYU6H4/edit#gid=0`
Switch to the `MLB Summary` tab.

Check:
- Block 1 shows reasonable numbers for "Total bets placed" and "Settled bets" (these should match roughly what you remember placing).
- Net P&L is a number, not `#ERROR!` or `#N/A`.
- FG + F5 "Total bets placed" roughly sums to the overall count (off by any unsettled bets).
- Daily P&L table has at least one row of dates with stake/P&L/cumulative filled in.
- Weekly P&L table has at least one row (week-of Monday date, stake, P&L).
- Equity curve chart is visible, with the cumulative line climbing/dipping over time.
- Weekly bar chart is visible, bars colored by sign.

- [ ] **Step 3: If verification fails, iterate**

Common failure modes and fixes:

| Symptom | Likely cause | Fix |
|---|---|---|
| "Total bets placed" = 0 | Filter regex didn't match real descriptions | Open the sheet, inspect a real MLB parlay row in Sheet1 column D, adjust `SPREAD_REGEX` / `TOTAL_REGEX` in create_mlb_summary.py, update test fixtures, re-run |
| "Total bets placed" > expected | False positives (maybe ML+total slipped through) | Tighten `SPREAD_REGEX` — likely matching odds rather than spread; add a stricter test case first |
| Daily table shows only one row but you have >1 date | `SORT(UNIQUE(FILTER(...)))` fell through to IFERROR `""` | Inspect `COL_DATE` range — Sheet1 column A should be real dates not strings; may need `DATEVALUE` wrap |
| `#REF!` in Cumulative P&L | MMULT range geometry mismatch | Check that the MMULT inner dimensions line up; ensure no blank rows inside the data |
| Charts are empty | Data range indices off | Double-check `daily_data_start` and `weekly_data_start` in `build_chart_requests` against the actual 0-based row of the first data row |

Iterate: edit the script, re-run (`./venv/bin/python3 create_mlb_summary.py` — it clears and rewrites the tab). Commit once it looks right.

- [ ] **Step 4: Commit any fixes**

```bash
cd /Users/callancapitolo/NFLWork-mlb-summary
git add bet_logger/create_mlb_summary.py bet_logger/test_mlb_summary.py
git commit -m "fix(bet_logger): tune MLB Summary regex/ranges after live-sheet verification"
```

(If no fixes were needed, skip this commit.)

---

## Task 10: Update documentation

**Files:**
- Modify: `bet_logger/README.md`

- [ ] **Step 1: Add an "MLB Summary tab" section**

Find the "Output Columns (Google Sheets)" section in `bet_logger/README.md`. Just above it, insert:

```markdown
## Summary Tabs

The sheet also includes auto-updating summary tabs generated by one-shot setup scripts. Run each once; after that, the tab's formulas refresh live whenever new bets are scraped into `Sheet1`.

### MLB Summary (`MLB Summary`)

Tracks performance of MLB spread+total correlated parlays placed on Wagerzon (both full-game and F5). Auto-populates from `Sheet1` via spreadsheet formulas — no helper sheet, no cron.

**What it tracks:** rows where `Sport = MLB`, `Bet Type = Parlay`, description starts with `PARLAY (2 TEAMS)`, and description contains both a spread token and a total token.

**What's on the tab:**
- Overall stats (bets placed, settled, wagered, P&L, ROI, win rate, record, avg odds)
- FG vs F5 split (same metrics, two rows)
- Daily P&L table (auto-extends as new dates arrive)
- Weekly P&L table (auto-extends)
- Equity curve chart (cumulative P&L over time)
- Weekly P&L bar chart (green/red bars)

**Setup (one time):**

```bash
./venv/bin/python3 create_mlb_summary.py
```

Re-run only if you need to regenerate the tab (it clears and rewrites).

### CBB Summary (`CBB Summary`)

See `create_summary.py`. Tracks NCAAM half-game answer-key bets across all platforms.
```

- [ ] **Step 2: Commit**

```bash
cd /Users/callancapitolo/NFLWork-mlb-summary
git add bet_logger/README.md
git commit -m "docs(bet_logger): document MLB Summary tab and setup script"
```

---

## Task 11: Merge, clean up, ask for approval

**Files:** none changed.

- [ ] **Step 1: Show the user the full diff for approval**

Run: `git -C /Users/callancapitolo/NFLWork-mlb-summary log main..HEAD --stat` and `git -C /Users/callancapitolo/NFLWork-mlb-summary diff main..HEAD`
Share the output with the user.

- [ ] **Step 2: Wait for explicit user approval before merging to main**

The user must say "yes, merge" or similar. Do not proceed without it.

- [ ] **Step 3: Merge to main**

```bash
cd /Users/callancapitolo/NFLWork
git checkout main
git pull --ff-only  # only if a remote is configured; skip otherwise
git merge --no-ff feature/mlb-summary-tab -m "feat(bet_logger): add MLB Summary tab"
```

- [ ] **Step 4: Clean up worktree and branch**

```bash
git -C /Users/callancapitolo/NFLWork worktree remove /Users/callancapitolo/NFLWork-mlb-summary
git -C /Users/callancapitolo/NFLWork branch -d feature/mlb-summary-tab
```

- [ ] **Step 5: Verify clean state**

Run: `git -C /Users/callancapitolo/NFLWork worktree list && git -C /Users/callancapitolo/NFLWork branch`
Expected: only the main worktree; `feature/mlb-summary-tab` gone.

---

## Self-Review Check

- **Spec coverage:**
  - Filter (Sheet1 sport/bet_type/header/spread/total) → Task 2 (`is_mlb_correlated_parlay`) + Task 3 (`FILTER_MASK`).
  - FG vs F5 classification → Task 2 (`is_f5_parlay`) + Task 3 (`F5_COND`, `FG_ADD`, `F5_ADD`).
  - Settlement treatment (win/loss/push, void→push upstream, pending excluded) → Task 3 (`WIN`/`LOSS`/`PUSH`/`SETTLED`).
  - Block 1 (Overall) → Task 4 rows 5-12.
  - Block 2 (FG/F5 split) → Task 4 rows 16-17.
  - Block 3 (daily + weekly, dynamic) → Task 5.
  - Equity curve + weekly bar charts → Task 6.
  - No helper sheet, no script on an ongoing basis → honored; only the one-shot setup script.
  - Version control (worktree, branch, single merge) → Tasks 1 and 11.
  - README update → Task 10.
- **Placeholder scan:** no TBDs, TODOs, or hand-wavy "handle edge cases" steps.
- **Type/name consistency:** `FILTER_MASK`, `FG_ADD`, `F5_ADD`, `SETTLED`, `WIN`, `LOSS`, `PUSH` defined once in Task 3, referenced identically in Tasks 4 and 5. `daily_start_row` / `weekly_start_row` stashed on the rows list in Task 5, consumed with identical names in Task 8. `get_sheet_id` and `build_chart_requests` / `build_format_requests` defined before `main()` uses them.
- **Known unknown (intentionally deferred to live-sheet verification in Task 9):** the exact `SPREAD_REGEX` and `TOTAL_REGEX` patterns — test fixtures cover expected Wagerzon shapes, but real Wagerzon output may include edge cases (fractional lines, team names with digits in them) not represented in the fixtures. Task 9 step 3 explicitly covers the iteration loop.
