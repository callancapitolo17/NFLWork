# MLB Bets Tab Improvements — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the MLB Dashboard bets-tab card with a 2×8 price-grid + hero-strip layout, render opponent/other-side rows on spreads, and fix two FanDuel scraper gaps (missing F7 markets, missing alt-spread exact lines).

**Architecture:**
- Backend: extend `expand_bets_to_book_prices()` in `Answer Keys/MLB Answer Key/odds_screen.R` so it emits an `opposite` row for spread / moneyline bets (it already does so for totals).
- Frontend: rewrite `create_bets_table()` in `Answer Keys/MLB Dashboard/mlb_dashboard.R` to render each bet as a card with (1) bet title, (2) matchup line, (3) hero strip with `WAGERZON +160 | FAIR +199 | EV +13% | RISK $41 | TO WIN $98 | [Place] [Log]`, (4) 2×8 price grid (sides × books). Rename `book_pill.R` → `book_cell.R` and replace `render_book_pill()` with `render_book_cell()`.
- Data fixes: extend `_FD_MARKET_WHITELIST` in `mlb_sgp/scraper_fanduel_singles.py` to include First-7-Innings markets (after recon), and resolve the FanDuel alt-run-line "closest line vs exact line" gap (scope determined by recon).

**Tech Stack:** R (Shiny + reactable + htmltools), Python 3 (FanDuel scraper), DuckDB, HTML/CSS for the card layout.

**Spec:** `docs/superpowers/specs/2026-05-13-mlb-bets-tab-improvements-design.md`
**Visual reference (locked design):** `.superpowers/brainstorm/37414-1778712316/content/bet-card-layouts-v8.html` (open via `http://localhost:60595` if companion server is still running)
**Worktree:** `.claude/worktrees/mlb-bets-tab-improvements/` (branch `worktree-mlb-bets-tab-improvements`) — rebased onto local `main` at start of plan-writing.

---

## File Structure

Files this plan creates / modifies:

| File | Action | Responsibility |
|------|--------|----------------|
| `Answer Keys/MLB Answer Key/odds_screen.R` | Modify (in `expand_bets_to_book_prices` + `.sides_for_bet`) | Emit opposite-side rows for spreads / moneyline. |
| `Answer Keys/tests/test_odds_screen.R` | Modify (add tests) | Cover new opposite-row behavior. |
| `Answer Keys/MLB Dashboard/book_pill.R` | Rename → `book_cell.R` and rewrite | Grid-cell renderer (4 states: `pick` / `exact` / `alt` / `empty`). |
| `Answer Keys/MLB Dashboard/mlb_dashboard.R` | Modify (lines ~1131–1506 + CSS block) | New `create_bets_table()` body that emits the V8 card. Inline CSS for `.bet-card-v8`. |
| `Answer Keys/MLB Dashboard/tests/test_book_cell.R` | Create | Unit tests for `render_book_cell()` HTML output. |
| `mlb_sgp/scraper_fanduel_singles.py` | Modify (`_FD_MARKET_WHITELIST`) | Add F7 market names after recon confirms strings. |
| `mlb_sgp/tests/test_fanduel_singles.py` | Modify (add tests) | Cover F7 classify_market + alt-run-line capture. |
| `Answer Keys/MLB Dashboard/README.md` | Modify | Update "bets tab" section to describe the new grid layout. |
| `Answer Keys/CLAUDE.md` | Modify | Update the "Odds screen" section (no more pill row; book_pill.R → book_cell.R). |
| `mlb_sgp/README.md` | Modify (if F7 added) | Note FD F7 coverage. |
| `docs/superpowers/research/2026-05-13-fd-recon-findings.md` | Create | Recon notes from Task 2 (one-pager; informs Tasks 3 + 4). |

---

## Task 1: Verify worktree environment + run baseline dashboard

**Files:**
- Read only — no edits.

**Why:** Before changing anything we need to confirm we can run the existing dashboard inside the worktree on a non-conflicting port, so we have a "before" baseline to compare against and a hot-loop for visual verification later.

- [ ] **Step 1: Confirm worktree branch state**

Run:
```bash
git status
git log --oneline -3
git branch --show-current
```

Expected: clean working tree; branch `worktree-mlb-bets-tab-improvements`; HEAD shows `docs(mlb-dashboard): spec updates from V8 design lock-in` then `docs(mlb-dashboard): bets tab improvements design spec` then `Merge feature/mlb-dk-fd-singles-scrapers into main`.

- [ ] **Step 2: Confirm key files exist with expected contents**

Run:
```bash
test -f "Answer Keys/MLB Answer Key/odds_screen.R" && echo "odds_screen OK"
test -f "Answer Keys/MLB Dashboard/book_pill.R" && echo "book_pill OK"
test -f "Answer Keys/MLB Dashboard/mlb_dashboard.R" && echo "mlb_dashboard OK"
test -f "mlb_sgp/scraper_fanduel_singles.py" && echo "fd scraper OK"
grep -n "_FD_MARKET_WHITELIST" "mlb_sgp/scraper_fanduel_singles.py" | head -3
grep -n "^create_bets_table " "Answer Keys/MLB Dashboard/mlb_dashboard.R"
```

Expected: 4 "OK" lines; `_FD_MARKET_WHITELIST` defined around line 61; `create_bets_table` defined around line 1143.

- [ ] **Step 3: Copy live `mlb_mm.duckdb` into the worktree (do NOT symlink)**

Per `CLAUDE.md`: "Never symlink DuckDB files. Always copy."

Run:
```bash
# Copy the live DB into a worktree-local path so we don't disturb production
mkdir -p "Answer Keys"
cp -p "/Users/callancapitolo/NFLWork/Answer Keys/mlb_mm.duckdb" "Answer Keys/mlb_mm.duckdb"
cp -p "/Users/callancapitolo/NFLWork/Answer Keys/mlb_dashboard.duckdb" "Answer Keys/MLB Dashboard/mlb_dashboard.duckdb" 2>/dev/null || true
cp -p "/Users/callancapitolo/NFLWork/Answer Keys/mlb.duckdb" "Answer Keys/mlb.duckdb" 2>/dev/null || true
ls -l "Answer Keys/mlb_mm.duckdb"
```

Expected: file copied successfully.

- [ ] **Step 4: Snapshot pre-change schema of `mlb_bets_book_prices`**

Run:
```bash
duckdb "Answer Keys/mlb_mm.duckdb" "DESCRIBE mlb_bets_book_prices" 2>&1 | head -20
duckdb "Answer Keys/mlb_mm.duckdb" "SELECT side, COUNT(*) FROM mlb_bets_book_prices GROUP BY side" 2>&1
```

Expected output (write what you see down for later comparison):
- Columns: `bet_row_id, game_id, market, period, side, bookmaker, line, line_quoted, is_exact_line, american_odds, fetch_time`
- Side breakdown: `pick` rows and (only for totals) `opposite` rows. After Task 5 this will include `opposite` for spreads too — that's how we'll know the change took effect.

- [ ] **Step 5: Confirm `fair_odds` column exists on `mlb_bets_combined`**

Run:
```bash
duckdb "Answer Keys/mlb_mm.duckdb" "DESCRIBE mlb_bets_combined" 2>&1 | grep -E "fair_odds|fair_prob|prob"
```

Expected: a `fair_odds` column (INTEGER) is present. If only `prob` is present and no `fair_odds`, **STOP** — the plan assumes `fair_odds` is upstream of the dashboard. In that case, before continuing, add a small derivation in MLB.R: `fair_odds = ifelse(prob >= 0.5, round(-100 * prob / (1 - prob)), round(100 * (1 - prob) / prob))`. (Per code reading at `mlb_dashboard.R:424`, `fair_odds` is already present and used in another tab; this step is a sanity check.)

- [ ] **Step 6: Start dashboard on alternate port to confirm baseline renders**

Run (in a terminal you can leave open):
```bash
cd "Answer Keys/MLB Dashboard"
PORT=8086 Rscript -e 'options(shiny.port=8086); source("mlb_dashboard.R"); shiny::runApp(shiny::shinyApp(ui=ui, server=server), port=8086, host="127.0.0.1", launch.browser=FALSE)'
```

If `mlb_dashboard.R` doesn't expose `ui` / `server` directly, use the existing runner. Check `mlb_dashboard.R` for the bottom-of-file `shinyApp(...)` invocation and figure out the right command — typical pattern is just `Rscript mlb_dashboard.R` and the port lives inside the script.

Expected: dashboard starts on `http://127.0.0.1:8086`. Open the bets tab. Confirm the existing pill rows render. Take note of one card (e.g., the Phillies @ Red Sox alt-spread) so you can compare visuals after.

- [ ] **Step 7: Stop the dashboard, commit nothing**

Kill the Rscript process. No changes to commit — Task 1 is read-only verification.

---

## Task 2: Recon FanDuel F7 + alternate-run-line markets

**Files:**
- Create: `docs/superpowers/research/2026-05-13-fd-recon-findings.md`
- Read only: `mlb_sgp/scraper_fanduel_singles.py`, `mlb_sgp/fd_client.py`

**Why:** Workstreams 3 and 4 are gated on what FD actually returns. The spec explicitly says: "if FD doesn't post F7 totals or paginates alt lines server-side, the fix path is different." We capture the findings in a tiny one-pager that Tasks 3 + 4 reference.

- [ ] **Step 1: Find a live FD MLB event we can probe**

Run:
```bash
duckdb fd_odds/fd.duckdb "SELECT DISTINCT game_id, away_team, home_team, game_time
  FROM mlb_odds WHERE game_time >= NOW() - INTERVAL '2 hours' LIMIT 5" 2>&1
```

Expected: a few rows with `game_id`. Note one as `$FD_EVENT_ID` for the next steps. If no rows, the slate is finished — pick yesterday's slate and live with the fact that markets won't be priced (we only need the structure of the API response).

- [ ] **Step 2: Pull the raw FD market list for that event**

Run a small Python recon snippet:
```bash
cd mlb_sgp
. venv/bin/activate 2>/dev/null || true   # if the venv exists
python -c "
from fd_client import FanDuelClient
c = FanDuelClient(verbose=True)
events = c.list_events()
target = events[0]   # first live event
print('event:', target.event_id, target.away_team, '@', target.home_team)
markets = c.fetch_event_markets(target.event_id)
print('TOTAL MARKETS:', len(markets))
print()
print('=== ALL F7 / 7-INNING MARKET NAMES ===')
for m in markets:
    if '7 Inning' in m.name or 'F7' in m.name or '1st 7' in m.name.lower() or 'first 7' in m.name.lower():
        print(' ', repr(m.name))
print()
print('=== ALL ALTERNATE-RUN-LINE MARKET NAMES ===')
for m in markets:
    if 'alternate' in m.name.lower() and 'run' in m.name.lower():
        print(' ', repr(m.name))
"
```

Expected output: a list of exact FD market-name strings for F7. **Write down each exact string verbatim** — they go in Task 3's whitelist entries. Common cases:
- "First 7 Innings Run Line"
- "First 7 Innings Total Runs"
- "First 7 Innings Money Line"
- "First 7 Innings Alternate Run Lines"
- "First 7 Innings Alternate Total Runs"

If FD uses different strings (e.g. "1st 7 Innings ..."), **use what you observed** — do not assume.

- [ ] **Step 3: Inspect the FD `Alternate Run Lines` payload for the alt-spread bug**

For a game where the model has a -2.5 bet, dump every runner returned by the alt-spread market:
```bash
python -c "
from fd_client import FanDuelClient
c = FanDuelClient(verbose=True)
events = c.list_events()
# Find an event where you know the model bets -2.5; if unsure, just dump the first event
target = events[0]
markets = c.fetch_event_markets(target.event_id)
alt_rl = [m for m in markets if m.name == 'Alternate Run Lines']
print('alt_rl markets:', len(alt_rl))
runners = c.fetch_event_runners(target.event_id)
for r in runners:
    for m in alt_rl:
        if r.market_id == m.market_id:
            print('  runner:', r.name, 'line=', r.line, 'odds=', r.american_odds)
"
```

Expected: a list of runners like:
```
runner: Yankees -1.5 line=-1.5 odds=...
runner: Red Sox +1.5 line=+1.5 odds=...
runner: Yankees -2.5 line=None odds=...     # <-- look at this!
runner: Red Sox +2.5 line=None odds=...
```

The critical thing to verify: does the payload include the -2.5 (or whatever line the model is betting)? If YES → the bug is in our scraper/loader, not FD. If NO → FD just doesn't post it on this game (no bug to fix, document as such).

- [ ] **Step 4: If -2.5 is in the payload, trace why it isn't in `fd_odds/fd.duckdb`**

Run:
```bash
duckdb fd_odds/fd.duckdb "SELECT period, market, home_spread, away_spread, home_spread_price, away_spread_price, fetch_time
  FROM mlb_odds
  WHERE game_id = '$FD_EVENT_ID' AND market = 'alternate_spreads'
  ORDER BY ABS(home_spread)" 2>&1
```

Expected: rows for every alt-line the scraper captured. Compare against Step 3's runner list. The gap is the bug.

Likely findings:
- (a) **All alt lines present in DB**: bug is downstream — see Step 5.
- (b) **Only ±1.5 in DB, ±2.5 missing**: the FD parser is dropping some lines. Look at `parse_runners_to_wide_rows` in `scraper_fanduel_singles.py` lines 89–196.
- (c) **No alt rows at all**: classify_market may be rejecting "Alternate Run Lines" (unlikely since it's whitelisted), or `client.fetch_event_markets` is paginating and missing pages.

- [ ] **Step 5: If alt lines are in `fd_odds/fd.duckdb` but not making it to `mlb_bets_book_prices`, trace `get_fd_odds()` in Tools.R**

Run R:
```r
source("Answer Keys/Tools.R")
fd <- get_fd_odds()
fd[fd$game_id == "$FD_EVENT_ID" & fd$market == "alternate_spreads", c("home_team","away_team","line","home_spread","away_spread","odds_home","odds_away")]
```

Expected: one row per alt-line. If only one line shows, `get_fd_odds()` is collapsing rows.

- [ ] **Step 6: Write recon findings to a one-pager**

Create `docs/superpowers/research/2026-05-13-fd-recon-findings.md`:
```markdown
# FanDuel recon — 2026-05-13

## F7 markets
FanDuel posts these market names for First-7-Innings:
- `<exact string>`
- `<exact string>`
- ...

(Or: "FanDuel does NOT post F7 markets for MLB" — if that turns out to be true.)

## Alt-run-line capture
Verified on event `<event_id>` (`<away_team>` @ `<home_team>`):
- FD API returns alt-spread runners at: ±0.5, ±1.5, ±2.5, ±3.5, ±4.5 (etc.)
- `fd_odds/fd.duckdb::mlb_odds` contains: <list what's actually persisted>
- `get_fd_odds()` in Tools.R returns: <list what's exposed>

## Bug location
<one of: scraper drops X / Tools.R collapses Y / FD itself doesn't post Z>

## Fix path
<one paragraph: what we'll change and where>
```

- [ ] **Step 7: Commit**

```bash
git add docs/superpowers/research/2026-05-13-fd-recon-findings.md
git commit -m "$(cat <<'EOF'
docs(mlb-bets-tab): FanDuel F7 + alt-run-line recon findings

Captures exact FD market-name strings for F7 markets and the trace of
the alt-spread -2.5 capture gap. Informs the FD whitelist extension
(Task 3) and the alt-line fix scope (Task 4).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Extend FanDuel whitelist with F7 markets

**Files:**
- Modify: `mlb_sgp/scraper_fanduel_singles.py:61` (the `_FD_MARKET_WHITELIST` dict)
- Modify: `mlb_sgp/tests/test_fanduel_singles.py` (or create if it doesn't exist)

**Why:** Workstream 3 from the spec. The whitelist gates which FD markets enter our pipeline; today it has FG + F5 but no F7, which is why FanDuel pills are blank on F7 totals cards.

- [ ] **Step 1: Check whether the test file exists**

Run:
```bash
ls mlb_sgp/tests/test_fanduel_singles.py 2>&1
```

If it exists, read it to understand the style. If not, you'll create one in Step 2.

- [ ] **Step 2: Write the failing test for F7 classification**

Add to `mlb_sgp/tests/test_fanduel_singles.py` (creating the file if needed):

```python
"""Tests for FanDuel singles scraper market classification."""
from mlb_sgp.scraper_fanduel_singles import classify_market


def test_classify_market_fg_main():
    assert classify_market("Run Line") == ("FG", "main")
    assert classify_market("Total Runs") == ("FG", "main")
    assert classify_market("Moneyline") == ("FG", "main")


def test_classify_market_fg_alt():
    assert classify_market("Alternate Run Lines") == ("FG", "alternate_spreads")
    assert classify_market("Alternate Total Runs") == ("FG", "alternate_totals")


def test_classify_market_f5_main():
    assert classify_market("First 5 Innings Run Line") == ("F5", "main")
    assert classify_market("First 5 Innings Total Runs") == ("F5", "main")
    assert classify_market("First 5 Innings Money Line") == ("F5", "main")


def test_classify_market_f5_alt():
    assert classify_market("First 5 Innings Alternate Run Lines") == ("F5", "alternate_spreads")
    assert classify_market("First 5 Innings Alternate Total Runs") == ("F5", "alternate_totals")


def test_classify_market_f7_main():
    # NEW — these should pass once the whitelist is extended
    assert classify_market("First 7 Innings Run Line") == ("F7", "main")
    assert classify_market("First 7 Innings Total Runs") == ("F7", "main")
    assert classify_market("First 7 Innings Money Line") == ("F7", "main")


def test_classify_market_f7_alt():
    assert classify_market("First 7 Innings Alternate Run Lines") == ("F7", "alternate_spreads")
    assert classify_market("First 7 Innings Alternate Total Runs") == ("F7", "alternate_totals")


def test_classify_market_rejects_unknown():
    assert classify_market("Random Garbage Market") is None
    assert classify_market("Yankees Total Runs") is None  # team total, out of scope
    assert classify_market("Total Runs (Bands)") is None
```

**Important:** if the recon in Task 2 found different exact strings (e.g. "1st 7 Innings Run Line"), use those instead of "First 7 Innings ..." in the F7 tests above.

- [ ] **Step 3: Run the tests to verify F7 tests fail**

Run:
```bash
cd mlb_sgp
python -m pytest tests/test_fanduel_singles.py -v
```

Expected: all FG / F5 tests PASS; all `test_classify_market_f7_*` tests FAIL with assertion errors (e.g., `assert None == ("F7", "main")`).

- [ ] **Step 4: Extend `_FD_MARKET_WHITELIST`**

In `mlb_sgp/scraper_fanduel_singles.py`, replace lines 61–76 (the existing whitelist dict) with:

```python
_FD_MARKET_WHITELIST: dict[str, tuple[str, str]] = {
    # FG main
    "Run Line":                              ("FG", "main"),
    "Total Runs":                            ("FG", "main"),
    "Moneyline":                             ("FG", "main"),
    # FG alt
    "Alternate Run Lines":                   ("FG", "alternate_spreads"),
    "Alternate Total Runs":                  ("FG", "alternate_totals"),
    # F5 main
    "First 5 Innings Run Line":              ("F5", "main"),
    "First 5 Innings Total Runs":            ("F5", "main"),
    "First 5 Innings Money Line":            ("F5", "main"),
    # F5 alt
    "First 5 Innings Alternate Run Lines":   ("F5", "alternate_spreads"),
    "First 5 Innings Alternate Total Runs":  ("F5", "alternate_totals"),
    # F7 main  (added 2026-05-13)
    "First 7 Innings Run Line":              ("F7", "main"),
    "First 7 Innings Total Runs":            ("F7", "main"),
    "First 7 Innings Money Line":            ("F7", "main"),
    # F7 alt
    "First 7 Innings Alternate Run Lines":   ("F7", "alternate_spreads"),
    "First 7 Innings Alternate Total Runs":  ("F7", "alternate_totals"),
}
```

**If recon found different strings**, use those exact strings instead of `First 7 Innings ...`.

- [ ] **Step 5: Re-run the tests to verify they pass**

Run:
```bash
cd mlb_sgp && python -m pytest tests/test_fanduel_singles.py -v
```

Expected: all tests PASS, including the four `test_classify_market_f7_*`.

- [ ] **Step 6: Run a live FD scrape against today's slate**

Run:
```bash
cd mlb_sgp
python scraper_fanduel_singles.py --verbose 2>&1 | tail -40
```

Expected: a non-zero row count printed; no errors. The verbose output should show F7 markets being classified now.

- [ ] **Step 7: Verify F7 rows landed in `fd_odds/fd.duckdb`**

Run:
```bash
duckdb fd_odds/fd.duckdb "SELECT period, market, COUNT(*)
  FROM mlb_odds
  WHERE fetch_time >= NOW() - INTERVAL '5 minutes'
  GROUP BY period, market
  ORDER BY period, market" 2>&1
```

Expected: rows for `F7` × `main`, `F7` × `alternate_spreads`, `F7` × `alternate_totals` in addition to existing FG / F5 rows.

- [ ] **Step 8: Commit**

```bash
git add mlb_sgp/scraper_fanduel_singles.py mlb_sgp/tests/test_fanduel_singles.py
git commit -m "$(cat <<'EOF'
feat(mlb-sgp/fd): add First-7-Innings markets to scraper whitelist

FD posts F7 main + alt run-line + alt total markets but the singles
scraper's strict whitelist had only FG + F5 entries, dropping F7
silently. Bets-tab F7 totals cards had blank FD pills as a result.

Tests cover classify_market() for the 5 new F7 strings plus the
existing FG/F5 entries.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Fix FanDuel alt-spread exact-line capture

**Files:**
- **Modify location depends on Task 2 recon findings**. The three possible cases below cover the realistic outcomes.
- Update tests in `mlb_sgp/tests/test_fanduel_singles.py` regardless.

**Why:** Workstream 4 from the spec. User flagged that FD shows `-1.5 +150` (closest alt line) on a bet where the model is on `-2.5` and FD actually posts `-2.5`. Fix path is recon-driven.

- [ ] **Step 1: Read the recon findings**

Read `docs/superpowers/research/2026-05-13-fd-recon-findings.md` from Task 2. Identify which case applies:
- **Case A:** FD parser is dropping some alt lines → fix in `parse_runners_to_wide_rows`.
- **Case B:** Tools.R `get_fd_odds()` is collapsing alt rows → fix in `Answer Keys/Tools.R`.
- **Case C:** FD doesn't post the line on that game → no code change; document in recon notes and skip Steps 2–7.

- [ ] **Step 2 (Case A): Write a failing test for the parser**

Add to `mlb_sgp/tests/test_fanduel_singles.py`:

```python
from datetime import datetime
from mlb_sgp.scraper_fanduel_singles import parse_runners_to_wide_rows
from mlb_sgp.fd_client import Event, Runner


def test_parse_runners_alt_spread_emits_all_lines():
    """FD alt-spread runners across +/-0.5, +/-1.5, +/-2.5 must all produce
    separate rows. Regression test for the bug where -2.5 was dropped."""
    event = Event(
        event_id="fd-test-1",
        away_team="Boston Red Sox",
        home_team="New York Yankees",
        start_time=datetime(2026, 5, 13, 19, 0),
    )
    runners = [
        # ±0.5
        Runner(market_id="mkt-alt", name="New York Yankees -0.5",
               line=None, american_odds=-180),
        Runner(market_id="mkt-alt", name="Boston Red Sox +0.5",
               line=None, american_odds=+150),
        # ±1.5
        Runner(market_id="mkt-alt", name="New York Yankees -1.5",
               line=None, american_odds=+130),
        Runner(market_id="mkt-alt", name="Boston Red Sox +1.5",
               line=None, american_odds=-160),
        # ±2.5
        Runner(market_id="mkt-alt", name="New York Yankees -2.5",
               line=None, american_odds=+220),
        Runner(market_id="mkt-alt", name="Boston Red Sox +2.5",
               line=None, american_odds=-280),
    ]
    market_meta = {"mkt-alt": ("FG", "alternate_spreads")}
    rows = parse_runners_to_wide_rows(event, runners, market_meta,
                                       fetch_time=datetime(2026, 5, 13))
    # Three buckets expected: |0.5|, |1.5|, |2.5|
    assert len(rows) == 3
    by_line = {abs(r["home_spread"]): r for r in rows}
    assert sorted(by_line) == [0.5, 1.5, 2.5]
    # Spot-check the -2.5 row
    r25 = by_line[2.5]
    assert r25["home_spread"] == -2.5
    assert r25["away_spread"] == +2.5
    assert r25["home_spread_price"] == +220
    assert r25["away_spread_price"] == -280
```

- [ ] **Step 3 (Case A): Run the test to verify it fails**

Run:
```bash
cd mlb_sgp && python -m pytest tests/test_fanduel_singles.py::test_parse_runners_alt_spread_emits_all_lines -v
```

Expected: FAIL. The failure mode tells you whether the bucket key is collapsing (`len(rows) < 3`) or one of the lines is being skipped (`KeyError: 2.5`).

- [ ] **Step 4 (Case A): Fix the parser**

The most likely bug is in the bucket-key construction at `scraper_fanduel_singles.py:151–157`:
```python
if market_type == "main":
    bucket_line: float | None = None
elif market_type == "alternate_spreads" and effective_line is not None:
    bucket_line = abs(effective_line)
else:
    bucket_line = effective_line
key = (period, market_type, bucket_line)
```

If `effective_line` is None for the -2.5 runner (because the regex `_FD_ALT_SPREAD_RE` failed to match), the bucket collapses. Check `r.name` against the regex:
```python
import re
_FD_ALT_SPREAD_RE = re.compile(r"^(?P<team>.+?)\s+(?P<line>[+-]\d+(?:\.\d+)?)\s*$")
print(_FD_ALT_SPREAD_RE.match("New York Yankees -2.5"))  # should match
print(_FD_ALT_SPREAD_RE.match("Yankees -2.5"))           # should match
```

If the regex doesn't match the real runner names, broaden it (e.g., allow unicode whitespace). Apply the minimal fix that makes the test pass.

- [ ] **Step 5 (Case A): Run the test to verify it passes**

```bash
cd mlb_sgp && python -m pytest tests/test_fanduel_singles.py -v
```

Expected: all tests PASS.

- [ ] **Step 2 (Case B): Write a failing test for `get_fd_odds()` in Tools.R**

(Skip if Case A applied.) Tools.R doesn't have unit tests in the same style, so write a smoke check instead:
```r
source("Answer Keys/Tools.R")
fd <- get_fd_odds()
n_alt <- fd %>% filter(market == "alternate_spreads", game_id == "<event_id>") %>% nrow()
stopifnot(n_alt >= 4)   # expect at least 4 alt-line rows for a typical game
```

- [ ] **Step 3 (Case B): Trace `get_fd_odds()` — find where the collapse happens**

Read `Answer Keys/Tools.R::get_fd_odds` (search for `get_fd_odds <-`). The function reads `mlb_odds` from `fd_odds/fd.duckdb` and explodes each wide row into three records (spread / totals / ML — see lines 3656–3697 per spec). The likely collapse: a `distinct()` or `slice_head()` somewhere is dropping rows that share `(game_id, market, period)` but differ in `line`.

Fix by adjusting the de-dup to use `line` as part of the key, or removing the collapse entirely. **Apply the minimal change.**

- [ ] **Step 4 (Case B): Verify the fix**

Re-run the smoke check from Step 2 — should now return ≥ 4 rows.

- [ ] **Step 6: Commit (works for either case)**

```bash
git add mlb_sgp/scraper_fanduel_singles.py mlb_sgp/tests/test_fanduel_singles.py "Answer Keys/Tools.R" 2>/dev/null
# Only add files you actually modified
git commit -m "$(cat <<'EOF'
fix(mlb-sgp/fd): capture all alternate-run-line lines, not just the closest

FD posts ±0.5 / ±1.5 / ±2.5 / ... alt-spread runners but our pipeline
was dropping the further-from-zero lines, causing the bets-tab FD pill
to show e.g. -1.5 +150 (closest match) on a bet for -2.5 that FD does
actually post.

See docs/superpowers/research/2026-05-13-fd-recon-findings.md for the
trace and exact bug location.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Emit opposite-side rows for spreads + moneyline in `expand_bets_to_book_prices`

**Files:**
- Modify: `Answer Keys/MLB Answer Key/odds_screen.R:106–114` (`.sides_for_bet`)
- Modify: `Answer Keys/MLB Answer Key/odds_screen.R:161–204` (the loop in `expand_bets_to_book_prices`)
- Modify: `Answer Keys/tests/test_odds_screen.R` (add tests for new behavior)

**Why:** Workstream 1 from the spec. Today the function emits `side = "opposite"` only for totals (Over/Under auto-derived). For spreads and ML, it requires the caller to pass `opposite_side` as a column on the bets frame, which MLB.R doesn't do. We're going to make the function derive `opposite_side` itself from `home_team`/`away_team` columns that are already on `mlb_bets_combined`.

- [ ] **Step 1: Read existing tests to learn the style**

```bash
sed -n '1,80p' "Answer Keys/tests/test_odds_screen.R"
```

Note the test framework (likely `testthat`), the file structure, and how existing tests build a `bets` frame + a `book_odds_by_book` list.

- [ ] **Step 2: Write a failing test for the spread opposite-row**

Append to `Answer Keys/tests/test_odds_screen.R`:

```r
test_that("expand_bets_to_book_prices emits opposite row for spread bets", {
  bets <- tibble(
    bet_row_id   = "row1",
    id           = "g1",
    home_team    = "New York Yankees",
    away_team    = "Boston Red Sox",
    market       = "spreads",
    market_type  = "spreads",
    period       = "FG",
    line         = -1.5,
    bet_on       = "New York Yankees"
  )

  # WZ has -1.5 on home and +1.5 on away.
  wz <- tibble(
    game_id       = "g1",
    market        = "spreads",
    period        = "FG",
    side          = c("New York Yankees", "Boston Red Sox"),
    line          = c(-1.5, 1.5),
    american_odds = c(-120L, 100L),
    fetch_time    = as.POSIXct("2026-05-13 12:00", tz = "UTC")
  )

  out <- expand_bets_to_book_prices(bets, list(wz = wz))

  expect_equal(nrow(out), 2)  # one pick row + one opposite row
  expect_setequal(out$side, c("pick", "opposite"))

  pick <- out[out$side == "pick", ]
  opp  <- out[out$side == "opposite", ]
  expect_equal(pick$american_odds, -120L)
  expect_equal(pick$line_quoted, -1.5)
  expect_equal(opp$american_odds,  100L)
  expect_equal(opp$line_quoted,    1.5)
})

test_that("expand_bets_to_book_prices emits opposite row for moneyline bets", {
  bets <- tibble(
    bet_row_id   = "row2",
    id           = "g2",
    home_team    = "New York Yankees",
    away_team    = "Boston Red Sox",
    market       = "h2h",
    market_type  = "h2h",
    period       = "FG",
    line         = NA_real_,
    bet_on       = "Boston Red Sox"
  )

  wz <- tibble(
    game_id       = "g2",
    market        = "h2h",
    period        = "FG",
    side          = c("Boston Red Sox", "New York Yankees"),
    line          = c(NA_real_, NA_real_),
    american_odds = c(150L, -170L),
    fetch_time    = as.POSIXct("2026-05-13 12:00", tz = "UTC")
  )

  out <- expand_bets_to_book_prices(bets, list(wz = wz))
  expect_equal(nrow(out), 2)
  expect_setequal(out$side, c("pick", "opposite"))
  expect_equal(out[out$side == "pick", ]$american_odds,  150L)
  expect_equal(out[out$side == "opposite", ]$american_odds, -170L)
})
```

- [ ] **Step 3: Run the tests to confirm they fail**

Run:
```bash
cd "Answer Keys" && Rscript -e 'testthat::test_file("tests/test_odds_screen.R")'
```

Expected: the two new tests FAIL with "nrow(out) == 1" (only pick, no opposite). All pre-existing tests still PASS.

- [ ] **Step 4: Update `.sides_for_bet` and `expand_bets_to_book_prices`**

In `Answer Keys/MLB Answer Key/odds_screen.R`:

Replace `.sides_for_bet` (lines 106–114) with:
```r
#' Map a bet's bet_on string to the "pick" and "opposite" side labels used
#' in the per-book odds frames.
#'
#' For totals, opposite is derived automatically (Over <-> Under).
#' For spreads, moneyline, and h2h, the opposite side is the OTHER team —
#' resolved from home_team / away_team if present, else from the explicit
#' opposite_side column, else NA (no opposite-side row emitted).
.sides_for_bet <- function(bet_on, market_type,
                            opposite_side = NA_character_,
                            home_team = NA_character_,
                            away_team = NA_character_) {
  if (market_type == "totals" || grepl("^alternate_totals", market_type) ||
      grepl("^totals_1st_[357]_innings$", market_type)) {
    list(pick = bet_on,
         opposite = if (grepl("^Over", bet_on, ignore.case = TRUE)) "Under" else "Over")
  } else {
    # Spreads / alt_spreads / h2h / moneyline: opposite is the other team.
    derived_opp <- if (!is.na(home_team) && !is.na(away_team)) {
      if (identical(bet_on, home_team)) away_team
      else if (identical(bet_on, away_team)) home_team
      else NA_character_
    } else NA_character_
    final_opp <- if (!is.na(derived_opp)) derived_opp else opposite_side
    list(pick = bet_on, opposite = final_opp)
  }
}
```

Then update the call site inside `expand_bets_to_book_prices` (currently around lines 162–165):

Replace:
```r
    opposite_col_value <- if ("opposite_side" %in% names(bet)) bet$opposite_side else NA_character_
    sides <- .sides_for_bet(bet$bet_on, bet$market_type, opposite_col_value)
    side_labels <- c(pick = sides$pick, opposite = sides$opposite)
```

With:
```r
    opposite_col_value <- if ("opposite_side" %in% names(bet)) bet$opposite_side else NA_character_
    home_team_value    <- if ("home_team"     %in% names(bet)) bet$home_team     else NA_character_
    away_team_value    <- if ("away_team"     %in% names(bet)) bet$away_team     else NA_character_
    sides <- .sides_for_bet(bet$bet_on, bet$market_type, opposite_col_value,
                            home_team = home_team_value, away_team = away_team_value)
    side_labels <- c(pick = sides$pick, opposite = sides$opposite)
```

- [ ] **Step 5: Re-run the tests**

```bash
cd "Answer Keys" && Rscript -e 'testthat::test_file("tests/test_odds_screen.R")'
```

Expected: ALL tests PASS, including the two new ones.

- [ ] **Step 6: Run a full pipeline pass and verify row counts**

```bash
cd "Answer Keys"
Rscript "MLB Answer Key/MLB.R" 2>&1 | tail -30
```

Expected: no errors. Then check row counts:
```bash
duckdb "Answer Keys/mlb_mm.duckdb" "SELECT side, COUNT(*) FROM mlb_bets_book_prices GROUP BY side" 2>&1
```

Expected: a meaningful number of `opposite` rows (more than before — for every spread bet × every book that quotes the opposite team). Compare against the baseline you snapshotted in Task 1 Step 4.

- [ ] **Step 7: Sanity-check downstream consumers don't break on the doubled rows**

Search for code that reads `mlb_bets_book_prices` without filtering on `side`:
```bash
grep -rn "mlb_bets_book_prices" --include='*.R' --include='*.py' 2>&1
```

Each match should either (a) filter on `side = 'pick'` or (b) be the dashboard renderer (which already handles both sides). If anything else reads this table un-filtered (combined parlay, bet logger, etc.), it may now double-count. Patch by adding a `WHERE side = 'pick'` filter.

- [ ] **Step 8: Commit**

```bash
git add "Answer Keys/MLB Answer Key/odds_screen.R" "Answer Keys/tests/test_odds_screen.R"
git commit -m "$(cat <<'EOF'
feat(mlb-bets-tab): emit opposite-side rows for spread + ML bets

Today expand_bets_to_book_prices only emits side='opposite' for totals
(Over <-> Under auto-derived). Spreads and moneyline bets get pick-only
because the function required the caller to pass opposite_side on the
bet row, which MLB.R doesn't do.

Internalize the derivation: when market_type is spreads / alt_spreads /
h2h, look at home_team / away_team on the bet row and resolve the
opposite team there. mlb_bets_combined already carries those columns,
so no upstream change needed.

Bets-tab spread cards will now render both BOS -2.5 and PHI +2.5
rows, mirroring how totals already render Over + Under.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Rename `book_pill.R` → `book_cell.R` and rewrite to a grid-cell renderer

**Files:**
- Move + rewrite: `Answer Keys/MLB Dashboard/book_pill.R` → `Answer Keys/MLB Dashboard/book_cell.R`
- Create: `Answer Keys/MLB Dashboard/tests/test_book_cell.R`

**Why:** Workstream 2 from the spec — the pill renderer becomes a grid-cell renderer with 4 states (`pick` / `exact` / `alt` / `empty`) and no inline book-label (the label moves to a column header above the grid).

- [ ] **Step 1: `git mv` the file to its new name**

Run:
```bash
git mv "Answer Keys/MLB Dashboard/book_pill.R" "Answer Keys/MLB Dashboard/book_cell.R"
```

- [ ] **Step 2: Write the failing test for `render_book_cell`**

Create `Answer Keys/MLB Dashboard/tests/test_book_cell.R`:

```r
library(testthat)
source("../book_cell.R")  # adjust path: from tests/ dir up to book_cell.R

test_that("render_book_cell empty state (no quote)", {
  out <- render_book_cell(american_odds = NA_integer_, line_quoted = NA_real_,
                          is_exact_line = NA, is_pick = FALSE,
                          side_word = "over", is_totals = TRUE)
  expect_match(out, 'class="cell empty"')
  expect_match(out, '&mdash;')
})

test_that("render_book_cell exact state (book matches model line, positive odds)", {
  out <- render_book_cell(american_odds = 160L, line_quoted = -2.5,
                          is_exact_line = TRUE, is_pick = FALSE,
                          side_word = "over", is_totals = FALSE)
  expect_match(out, 'class="cell exact"')
  expect_match(out, '\\+160')
  # No alt-line tag rendered when is_exact_line == TRUE
  expect_no_match(out, 'class="alt-line"')
})

test_that("render_book_cell alt state for spread (signed line tag)", {
  out <- render_book_cell(american_odds = -174L, line_quoted = 1.5,
                          is_exact_line = FALSE, is_pick = FALSE,
                          side_word = "under", is_totals = FALSE)
  expect_match(out, 'class="cell alt"')
  expect_match(out, '<span class="alt-line">\\+1.5</span>')
  expect_match(out, '-174')
})

test_that("render_book_cell alt state for totals (O/U prefix on line tag)", {
  out <- render_book_cell(american_odds = -115L, line_quoted = 7.5,
                          is_exact_line = FALSE, is_pick = FALSE,
                          side_word = "over", is_totals = TRUE)
  expect_match(out, '<span class="alt-line">O7.5</span>')
})

test_that("render_book_cell pick state overrides others", {
  out <- render_book_cell(american_odds = 160L, line_quoted = -2.5,
                          is_exact_line = TRUE, is_pick = TRUE,
                          side_word = "over", is_totals = FALSE)
  expect_match(out, 'class="cell pick"')
  expect_match(out, '\\+160')
})
```

- [ ] **Step 3: Run the test to verify it fails**

Run:
```bash
cd "Answer Keys/MLB Dashboard"
Rscript -e 'testthat::test_file("tests/test_book_cell.R")'
```

Expected: all 5 tests FAIL because `render_book_cell` doesn't exist yet.

- [ ] **Step 4: Replace the contents of `book_cell.R`**

Write the entire new file `Answer Keys/MLB Dashboard/book_cell.R`:

```r
# Answer Keys/MLB Dashboard/book_cell.R
# Render one grid cell for the bets-tab V8 card.
#
# Four visual states:
#   - empty:   book has no quote on this side -> dim em-dash
#   - exact:   book's line matches the bet's model line -> plain price
#   - alt:     book's closest line differs from model line -> amber background +
#              line tag above the price ("-1.5" for spreads, "O7.5" / "U5.5"
#              for totals)
#   - pick:    overrides exact/alt -> green border + green price (applied when
#              is_pick = TRUE; pick is exact by construction in expand_bets_to_book_prices)
#
# The book label (WZ, H88, etc.) is NOT rendered here — it sits in a column
# header above the grid. This keeps the cell compact.

#' Format a line value for the alt-line tag.
#' (Carried over from the old book_pill.R unchanged.)
.format_line_value <- function(x, signed = FALSE) {
  if (is.na(x)) return("")
  base <- if (x == round(x)) format(as.integer(x))
          else sub("\\.?0+$", "", sprintf("%.2f", x))
  if (signed && x > 0) paste0("+", base) else base
}

#' Render one bets-tab grid cell.
#'
#' @param american_odds Integer odds (e.g., 125, -110). NA -> empty state.
#' @param line_quoted Numeric line the book is showing on this side.
#' @param is_exact_line Boolean: TRUE when book's line matches the model line
#'   exactly; FALSE -> alt state.
#' @param is_pick TRUE if this is the picked book on the pick side; overrides
#'   exact/alt with the green pick state.
#' @param side_word "over" or "under" — only used for the O/U prefix on a
#'   mismatched totals line tag.
#' @param is_totals TRUE for totals markets (line tag gets O/U prefix);
#'   FALSE for spreads (signed line value, e.g. "-1.5").
#' @return HTML string for the cell (a single <div class="cell ..."> ... </div>).
render_book_cell <- function(american_odds, line_quoted, is_exact_line,
                              is_pick = FALSE, side_word = "over",
                              is_totals = TRUE) {
  # State 1: empty (no quote)
  if (is.na(american_odds)) {
    return('<div class="cell empty"><span class="price">&mdash;</span></div>')
  }

  price_str <- if (american_odds > 0) paste0("+", american_odds) else as.character(american_odds)

  is_mismatched <- !isTRUE(is_exact_line)

  cell_class <- if (is_pick) "cell pick"
                else if (is_mismatched) "cell alt"
                else "cell exact"

  tag_html <- ""
  if (is_mismatched && !is_pick) {
    if (is_totals) {
      prefix <- if (side_word == "under") "U" else "O"
      tag_html <- sprintf('<span class="alt-line">%s%s</span>',
                          prefix, .format_line_value(line_quoted))
    } else {
      tag_html <- sprintf('<span class="alt-line">%s</span>',
                          .format_line_value(line_quoted, signed = TRUE))
    }
  }

  sprintf('<div class="%s">%s<span class="price">%s</span></div>',
          cell_class, tag_html, price_str)
}

# Backwards-compat shim: old code may still source book_pill.R via legacy paths.
# Provide render_book_pill as a thin wrapper so a stale source() doesn't crash
# during the transition. Remove after Task 7 lands.
render_book_pill <- function(book, american_odds, line_quoted, is_exact_line,
                              is_pick = FALSE, side = "over", is_totals = TRUE) {
  warning("render_book_pill is deprecated; use render_book_cell. The book ",
          "label is no longer part of the cell — it lives in the column header.",
          call. = FALSE)
  render_book_cell(american_odds, line_quoted, is_exact_line, is_pick,
                   side_word = side, is_totals = is_totals)
}
```

- [ ] **Step 5: Re-run the tests to verify they pass**

```bash
cd "Answer Keys/MLB Dashboard"
Rscript -e 'testthat::test_file("tests/test_book_cell.R")'
```

Expected: all 5 tests PASS.

- [ ] **Step 6: Update `mlb_dashboard.R` to source the renamed file**

```bash
grep -n 'book_pill\.R' "Answer Keys/MLB Dashboard/mlb_dashboard.R"
```

Replace any `source(".../book_pill.R")` with `source(".../book_cell.R")`. (The shim provides `render_book_pill` so the legacy `create_bets_table_legacy` path still works without further changes.)

- [ ] **Step 7: Commit**

```bash
git add "Answer Keys/MLB Dashboard/book_cell.R" \
        "Answer Keys/MLB Dashboard/tests/test_book_cell.R" \
        "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "$(cat <<'EOF'
refactor(mlb-bets-tab): rename book_pill.R -> book_cell.R; new cell renderer

Drop-in renamed file for the new V8 grid layout. render_book_cell()
emits a <div class="cell ..."> with 4 states (empty / exact / alt /
pick) and no inline book label — the label lives in the column header
above the grid instead.

Legacy render_book_pill is preserved as a deprecation shim so the
create_bets_table_legacy fallback still works during transition.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: Add card render helpers to `mlb_dashboard.R`

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` (insert new helpers just above `create_bets_table`, around line 1130)

**Why:** Three new helpers that the new `create_bets_table` will use: `render_price_grid_row`, `render_hero_strip`, `render_bet_card`. Putting them above the function keeps the call graph readable.

- [ ] **Step 1: Find the insertion point**

Run:
```bash
grep -n "^# create_bets_table " "Answer Keys/MLB Dashboard/mlb_dashboard.R"
```

Expected: a line like `1131:# create_bets_table — card layout (Task 11 of odds-screen rebuild)`. Insert the new helpers just above this line.

- [ ] **Step 2: Insert the helpers**

Open `Answer Keys/MLB Dashboard/mlb_dashboard.R` and insert the following block immediately before line 1131 (`# create_bets_table — card layout ...`):

```r
# =============================================================================
# Bets-tab V8 card helpers
#
# The V8 card has three components stacked top-to-bottom:
#   1. Bet title row (uniform 18px) + matchup row (uniform 15px)
#   2. Hero strip — pick book/odds, Fair, EV, Risk, To Win, [Place] [Log]
#   3. Price grid — 2 rows (pick side / opposite side) × 8 book columns
#
# Render helpers below build one card given precomputed inputs.
# CSS class names match the `.bet-card-v8` block added to mlb_dashboard.R
# inline-CSS (see Task 8).
# =============================================================================

BOOK_ORDER_V8  <- c("wagerzon", "hoop88", "bfa", "bookmaker", "bet105",
                    "draftkings", "fanduel", "pinnacle")
BOOK_LABELS_V8 <- c(wagerzon = "WZ", hoop88 = "H88", bfa = "BFA",
                    bookmaker = "BKM", bet105 = "B105",
                    draftkings = "DK", fanduel = "FD", pinnacle = "PINN")

#' Render one row of the price grid (1 row label + 8 book cells).
#'
#' @param wide_row Single-row tibble from book_prices_wide (the pivoted side's data)
#' @param side_label Text for the leftmost cell, e.g. "BOS -2.5" or "Over 5.5"
#' @param is_pick_side TRUE for the pick row; FALSE for the opposite row
#' @param pick_book Bookmaker key (e.g. "wagerzon") of the picked book
#' @param side_word "over" or "under" — drives O/U prefix on alt-line tag
#' @param is_totals Boolean — drives spread vs totals line-tag formatting
#' @return HTML string: row-label cell + 8 grid cells.
render_price_grid_row <- function(wide_row, side_label, is_pick_side,
                                   pick_book, side_word, is_totals) {
  cells <- vapply(BOOK_ORDER_V8, function(b) {
    odds_col  <- paste0(b, "_american_odds")
    lq_col    <- paste0(b, "_line_quoted")
    exact_col <- paste0(b, "_is_exact_line")
    odds  <- if (!is.null(wide_row) && odds_col  %in% names(wide_row)) wide_row[[odds_col]]  else NA_integer_
    lq    <- if (!is.null(wide_row) && lq_col    %in% names(wide_row)) wide_row[[lq_col]]    else NA_real_
    exact <- if (!is.null(wide_row) && exact_col %in% names(wide_row)) wide_row[[exact_col]] else NA
    render_book_cell(
      american_odds = if (is.na(odds)) NA_integer_ else as.integer(odds),
      line_quoted   = lq,
      is_exact_line = exact,
      is_pick       = is_pick_side && (b == pick_book),
      side_word     = side_word,
      is_totals     = is_totals
    )
  }, character(1))

  paste0(
    sprintf('<div class="row-hdr">%s</div>',
            htmltools::htmlEscape(side_label)),
    paste(cells, collapse = "")
  )
}

#' Format an American-odds integer with explicit sign ("-120" or "+160").
.format_odds_signed <- function(odds) {
  if (is.na(odds)) return("--")
  if (odds > 0) sprintf("+%d", as.integer(odds))
  else sprintf("%d", as.integer(odds))
}

#' Render the green-tinted hero strip — pick block, divider, four centered
#' stats, action buttons.
#'
#' @param pick_book Display label (e.g. "WAGERZON")
#' @param pick_odds Integer American odds
#' @param fair_odds Integer American odds (de-vigged fair)
#' @param ev_pct Percentage value (e.g. 13.0 for +13.0%)
#' @param risk_dollars Numeric (e.g. 41)
#' @param towin_dollars Numeric (e.g. 98)
#' @param action_html Pre-rendered HTML for the right-aligned action buttons.
#' @return HTML string for the entire <div class="hero">.
render_hero_strip <- function(pick_book, pick_odds, fair_odds,
                               ev_pct, risk_dollars, towin_dollars,
                               action_html) {
  ev_str   <- sprintf("+%.1f%%", ev_pct)
  risk_str <- sprintf("$%.0f", risk_dollars)
  win_str  <- sprintf("$%.0f", towin_dollars)
  fair_str <- .format_odds_signed(fair_odds)
  odds_str <- .format_odds_signed(pick_odds)

  sprintf(
    '<div class="hero">
       <div class="pick">
         <span class="book">%s</span>
         <span class="odds">%s</span>
       </div>
       <div class="divider"></div>
       <div class="stat"><span class="lbl">Fair</span><span class="val fair">%s</span></div>
       <div class="stat"><span class="lbl">EV</span><span class="val ev">%s</span></div>
       <div class="stat"><span class="lbl">Risk</span><span class="val risk">%s</span></div>
       <div class="stat"><span class="lbl">To Win</span><span class="val win">%s</span></div>
       <div class="actions">%s</div>
     </div>',
    htmltools::htmlEscape(pick_book),
    odds_str, fair_str, ev_str, risk_str, win_str, action_html
  )
}

#' Render the bet title row + matchup row + hero strip + price grid for one bet.
#'
#' All inputs are pre-computed by the caller. This function does no data
#' lookups — keeping it pure makes testing trivial.
#'
#' @param bet_title Text for H1, e.g. "BOS -2.5"
#' @param market_label Secondary text on same line, e.g. "Alt Spread · Full Game"
#' @param matchup Text for H2, e.g. "Philadelphia Phillies @ Boston Red Sox"
#' @param tipoff Smaller secondary text on same line, e.g. "Wed 10:46 PM"
#' @param hero_html Output of render_hero_strip
#' @param pickside_label e.g. "BOS -2.5" (left-most cell of pick row)
#' @param pickside_wide_row Single-row tibble from book_prices_wide
#' @param oppside_label e.g. "PHI +2.5" (NULL/NA to skip the opposite row,
#'   which happens when no opposite data exists)
#' @param oppside_wide_row Single-row tibble or NULL
#' @param pick_book Bookmaker key for green-highlighting in the grid
#' @param side_word "over"/"under"
#' @param is_totals Boolean
#' @param corr_badge_html Optional same-game-corr badge HTML (or "")
#' @return HTML string — the entire card.
render_bet_card <- function(bet_title, market_label, matchup, tipoff, hero_html,
                             pickside_label, pickside_wide_row,
                             oppside_label, oppside_wide_row,
                             pick_book, side_word, is_totals,
                             corr_badge_html = "") {
  bet_line <- sprintf(
    '<div class="bet-line">
       <span class="primary">%s</span>
       <span class="sep">&middot;</span>
       <span class="secondary">%s</span>
       %s
     </div>',
    htmltools::htmlEscape(bet_title),
    htmltools::htmlEscape(market_label),
    corr_badge_html
  )

  matchup_line <- sprintf(
    '<div class="matchup-line">
       <span class="primary">%s</span>
       <span class="sep">&middot;</span>
       <span class="secondary">%s</span>
     </div>',
    htmltools::htmlEscape(matchup),
    htmltools::htmlEscape(tipoff)
  )

  # Column header row (book labels)
  col_hdrs <- paste(
    "<div></div>",  # blank corner cell over the row-label column
    paste(sprintf('<div class="col-hdr">%s</div>', BOOK_LABELS_V8[BOOK_ORDER_V8]),
          collapse = ""),
    sep = ""
  )

  pick_row <- render_price_grid_row(pickside_wide_row, pickside_label,
                                     is_pick_side = TRUE,
                                     pick_book = pick_book,
                                     side_word = side_word,
                                     is_totals = is_totals)

  opp_row <- if (!is.null(oppside_wide_row) && !is.na(oppside_label) &&
                 nchar(oppside_label) > 0) {
    render_price_grid_row(oppside_wide_row, oppside_label,
                           is_pick_side = FALSE,
                           pick_book = pick_book,
                           side_word = if (side_word == "over") "under" else "over",
                           is_totals = is_totals)
  } else ""

  grid_html <- paste0('<div class="price-grid">', col_hdrs, pick_row, opp_row, '</div>')

  paste0('<div class="bet-card-v8">', bet_line, matchup_line, hero_html, grid_html, '</div>')
}
```

- [ ] **Step 3: Smoke-test the helpers manually**

Run a quick R smoke check:
```bash
cd "Answer Keys/MLB Dashboard"
Rscript -e '
source("book_cell.R")
source("mlb_dashboard.R", local = new.env())  # WARNING: this might run all of the dashboard;
# if it does, abort. Instead, source just the helpers section by copying it into a tmp file
# or extracting the relevant lines.
'
```

If sourcing the whole file is too heavy (it runs the Shiny app), copy the helpers block from Step 2 into a temp file `/tmp/v8_helpers.R`, then:
```r
source("book_cell.R")
source("/tmp/v8_helpers.R")
hero <- render_hero_strip(pick_book = "WAGERZON", pick_odds = 160L,
                          fair_odds = 199L, ev_pct = 13.0,
                          risk_dollars = 41, towin_dollars = 98,
                          action_html = '<button>Place</button>')
cat(hero, "\n")
```

Expected: clean HTML output with the four stats visible.

- [ ] **Step 4: Commit**

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "$(cat <<'EOF'
feat(mlb-bets-tab): add V8 card render helpers

Pure helper functions that build the new card HTML given precomputed
inputs:
- render_price_grid_row: one row of the 8-book price grid
- render_hero_strip: green-tinted band with pick + Fair + EV + Risk +
  To Win + action buttons (stats centered)
- render_bet_card: the whole card assembled

These plug into create_bets_table in a follow-up commit.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: Add new card CSS to `mlb_dashboard.R`

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R` (CSS block — find the existing bets-table CSS around line 2544 and add the new block adjacent to it)

**Why:** The HTML emitted by the new render helpers needs CSS classes defined. We add a single self-contained block so all the V8 styles live together.

- [ ] **Step 1: Find the existing CSS block**

```bash
grep -n "bets-table-container\|cell-game\|cell-ev" "Answer Keys/MLB Dashboard/mlb_dashboard.R" | head -10
```

Expected: a block around line 2544 containing CSS-in-R-string for the existing card. Add the new block immediately after that one (so the order in the rendered `<style>` tag is: old block, new block).

- [ ] **Step 2: Insert the V8 CSS block**

Append this CSS string (R `paste0` or `HTML` style — match whatever the existing block uses, likely a long string passed to `tags$style` or written into an `HTML(...)` block). The V8 styles:

```css
/* === V8 BET CARD =========================================================== */
.bet-card-v8 {
  background: #0d1117;
  border: 1px solid #30363d;
  border-radius: 12px;
  padding: 18px 20px;
  font-family: -apple-system, BlinkMacSystemFont, "SF Pro Text", Inter, sans-serif;
  color: #c9d1d9;
  margin-bottom: 14px;
}

.bet-card-v8 .bet-line {
  font-size: 18px;
  font-weight: 600;
  letter-spacing: -0.01em;
  line-height: 1.25;
}
.bet-card-v8 .bet-line .primary    { color: #f0f6fc; font-weight: 700; }
.bet-card-v8 .bet-line .sep        { color: #484f58; margin: 0 8px; }
.bet-card-v8 .bet-line .secondary  { color: #8b949e; font-weight: 500; }

.bet-card-v8 .matchup-line {
  font-size: 15px;
  font-weight: 500;
  line-height: 1.3;
  margin-top: 6px;
}
.bet-card-v8 .matchup-line .primary    { color: #c9d1d9; }
.bet-card-v8 .matchup-line .sep        { color: #484f58; margin: 0 8px; }
.bet-card-v8 .matchup-line .secondary  { color: #8b949e; }

.bet-card-v8 .hero {
  display: flex;
  align-items: center;
  gap: 24px;
  flex-wrap: wrap;
  margin: 14px 0 16px;
  padding: 12px 14px;
  background: linear-gradient(90deg, rgba(63,185,80,0.10) 0%, rgba(63,185,80,0.04) 60%, transparent 100%);
  border: 1px solid rgba(63,185,80,0.30);
  border-radius: 8px;
}
.bet-card-v8 .hero .pick { display: flex; flex-direction: column; gap: 2px; }
.bet-card-v8 .hero .pick .book {
  color: #3fb950; font-weight: 700; font-size: 13px; letter-spacing: 0.04em;
}
.bet-card-v8 .hero .pick .odds {
  color: #f0f6fc; font-weight: 700; font-size: 22px;
  font-variant-numeric: tabular-nums; line-height: 1;
}
.bet-card-v8 .hero .divider {
  width: 1px; height: 32px; background: rgba(255,255,255,0.08);
}
.bet-card-v8 .hero .stat {
  display: flex; flex-direction: column; gap: 3px;
  align-items: center; text-align: center; min-width: 60px;
}
.bet-card-v8 .hero .stat .lbl {
  color: #8b949e; font-size: 10px; letter-spacing: 0.08em; text-transform: uppercase;
}
.bet-card-v8 .hero .stat .val {
  color: #f0f6fc; font-weight: 700; font-size: 18px;
  font-variant-numeric: tabular-nums; line-height: 1;
}
.bet-card-v8 .hero .stat .val.ev    { color: #3fb950; }
.bet-card-v8 .hero .stat .val.win   { color: #3fb950; }
.bet-card-v8 .hero .stat .val.risk  { color: #c9d1d9; }
.bet-card-v8 .hero .stat .val.fair  { color: #c9d1d9; }
.bet-card-v8 .hero .actions {
  display: flex; gap: 8px; margin-left: 8px;
}
.bet-card-v8 .hero .actions button {
  border: 1px solid #30363d; border-radius: 6px; padding: 9px 18px;
  color: #c9d1d9; background: #161b22;
  font-size: 12px; font-weight: 600; cursor: pointer;
}
.bet-card-v8 .hero .actions button.btn-place {
  border-color: #3fb950; color: #ffffff;
  background: rgba(63,185,80,0.28); font-weight: 700;
}

.bet-card-v8 .price-grid {
  display: grid;
  grid-template-columns: minmax(85px, 125px) repeat(8, minmax(0, 1fr));
  gap: 3px;
  font-variant-numeric: tabular-nums;
}
.bet-card-v8 .price-grid .col-hdr {
  text-align: center; color: #8b949e; font-size: 9px;
  letter-spacing: 0.08em; text-transform: uppercase;
  padding: 0 2px 6px 2px; border-bottom: 1px solid #21262d;
}
.bet-card-v8 .price-grid .row-hdr {
  color: #c9d1d9; font-weight: 600; font-size: 12px;
  display: flex; align-items: center; padding: 8px 8px 8px 0;
}
.bet-card-v8 .price-grid .cell {
  display: flex; flex-direction: column;
  align-items: center; justify-content: center;
  padding: 6px 2px;
  background: #0d1117; border: 1px solid transparent; border-radius: 5px;
  min-height: 36px; font-size: 12px;
}
.bet-card-v8 .price-grid .cell .price       { color: #c9d1d9; font-weight: 600; }
.bet-card-v8 .price-grid .cell .alt-line    { color: #d29922; font-size: 10px; margin-bottom: 1px; }
.bet-card-v8 .price-grid .cell.empty .price { color: #484f58; }
.bet-card-v8 .price-grid .cell.exact        { background: #161b22; }
.bet-card-v8 .price-grid .cell.alt          { background: rgba(210,153,34,0.06); }
.bet-card-v8 .price-grid .cell.pick {
  background: rgba(63,185,80,0.18); border-color: #3fb950;
}
.bet-card-v8 .price-grid .cell.pick .price { color: #3fb950; }
```

Wrap into a raw R string. If the existing CSS lives in a `tags$style(HTML("..."))` call, append to that string. Otherwise, add a new `tags$style(HTML("..."))` block to the UI.

- [ ] **Step 3: Restart the dashboard and confirm CSS loads**

```bash
cd "Answer Keys/MLB Dashboard"
Rscript mlb_dashboard.R 2>&1 | head -20
# (use whatever the existing entry-point command is)
```

Open the bets tab in the browser. The old layout will still render (we haven't wired the new helpers in yet), but inspect the page source — search for `.bet-card-v8` in the `<style>` tag. Expected: the CSS block is in the page, just unused.

- [ ] **Step 4: Commit**

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "$(cat <<'EOF'
style(mlb-bets-tab): add V8 card CSS block

CSS for .bet-card-v8 with all sub-classes (bet-line, matchup-line, hero,
price-grid, cell states). Adjacent to the existing bets-table CSS in
mlb_dashboard.R. Not wired to any HTML yet — that lands in Task 9.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 9: Wire `create_bets_table` to use the V8 helpers

**Files:**
- Modify: `Answer Keys/MLB Dashboard/mlb_dashboard.R:1131–1506` (replace body of `create_bets_table`)

**Why:** This is the integration step. Replace the existing per-cell `colDef` rendering with a single card-html column driven by the new helpers from Task 7.

- [ ] **Step 1: Add a `format_market_for_card` helper**

Look at the existing `format_market_name` function used in the current implementation. Find it:
```bash
grep -n "^format_market_name " "Answer Keys/MLB Dashboard/mlb_dashboard.R"
```

Read its definition. We'll reuse it for the "Alt Spread · Full Game" market label string. If it produces a too-long string, write a thin wrapper that strips redundancy.

Insert just above `create_bets_table` (so above the helper block from Task 7):

```r
#' Produce the human-readable market label for the V8 card title row.
#'
#' E.g.:  market="alternate_spreads",   bet_on="Boston Red Sox", line=-2.5
#'        -> bet_title="BOS -2.5",  market_label="Alt Spread · Full Game"
#'
#' Returns list(bet_title, market_label).
format_market_for_card <- function(market, bet_on, line, home_team, away_team) {
  is_total <- grepl("^Over|^Under", bet_on, ignore.case = TRUE)

  # Bet title — what is the bet?
  bet_title <- if (is_total) {
    line_str <- if (!is.na(line)) {
      if (line == round(line)) format(as.integer(line))
      else sub("\\.?0+$", "", sprintf("%.2f", line))
    } else ""
    paste(bet_on, line_str)
  } else if (is.na(line)) {
    # Moneyline — just the team
    .team_abbr(bet_on)
  } else {
    # Spread — team + signed line
    line_str <- if (line > 0) paste0("+", line) else as.character(line)
    paste(.team_abbr(bet_on), line_str)
  }

  # Market label — what kind of bet?
  period <- if (grepl("_1st_3_innings$", market)) "1st 3 Innings"
            else if (grepl("_1st_5_innings$", market)) "1st 5 Innings"
            else if (grepl("_1st_7_innings$", market)) "1st 7 Innings"
            else "Full Game"

  market_kind <- if (grepl("^alternate_spreads", market))  "Alt Spread"
                 else if (grepl("^alternate_totals", market)) "Alt Total"
                 else if (grepl("^spreads", market))           "Spread"
                 else if (grepl("^totals", market))            "Total"
                 else if (grepl("^h2h", market))               "Moneyline"
                 else gsub("_", " ", market, fixed = TRUE)

  list(bet_title    = bet_title,
       market_label = paste(market_kind, period, sep = " · "))
}

#' 3-letter team abbreviation used in bet-title and side-label.
.team_abbr <- function(full_name) {
  # Simple map — extend if needed. Falls back to first three letters.
  abbr_map <- c(
    "Arizona Diamondbacks" = "ARI", "Atlanta Braves" = "ATL",
    "Baltimore Orioles" = "BAL",   "Boston Red Sox" = "BOS",
    "Chicago Cubs" = "CHC",        "Chicago White Sox" = "CHW",
    "Cincinnati Reds" = "CIN",     "Cleveland Guardians" = "CLE",
    "Colorado Rockies" = "COL",    "Detroit Tigers" = "DET",
    "Houston Astros" = "HOU",      "Kansas City Royals" = "KCR",
    "Los Angeles Angels" = "LAA",  "Los Angeles Dodgers" = "LAD",
    "Miami Marlins" = "MIA",       "Milwaukee Brewers" = "MIL",
    "Minnesota Twins" = "MIN",     "New York Mets" = "NYM",
    "New York Yankees" = "NYY",    "Athletics" = "ATH",
    "Oakland Athletics" = "OAK",   "Philadelphia Phillies" = "PHI",
    "Pittsburgh Pirates" = "PIT",  "San Diego Padres" = "SDP",
    "San Francisco Giants" = "SFG","Seattle Mariners" = "SEA",
    "St. Louis Cardinals" = "STL", "Tampa Bay Rays" = "TBR",
    "Texas Rangers" = "TEX",       "Toronto Blue Jays" = "TOR",
    "Washington Nationals" = "WSN"
  )
  abbr_map[[full_name]] %||% toupper(substr(full_name, 1, 3))
}

# %||% null-coalesce — define if not already defined
if (!exists("%||%")) `%||%` <- function(x, y) if (is.null(x)) y else x
```

- [ ] **Step 2: Replace the body of `create_bets_table`**

Find `create_bets_table <- function(all_bets, placed_bets, book_prices_wide = NULL) {` (around line 1143) and replace from there through the closing `}` of the function (around line 1506) with:

```r
create_bets_table <- function(all_bets, placed_bets, book_prices_wide = NULL) {
  if (is.null(book_prices_wide)) {
    warning("[bets-tab] book_prices_wide is NULL — using legacy flat layout")
    return(create_bets_table_legacy(all_bets, placed_bets))
  }

  placed_hashes <- if (nrow(placed_bets) > 0) placed_bets$bet_hash else character()

  placed_status_lookup <- if (nrow(placed_bets) > 0 && "status" %in% names(placed_bets)) {
    setNames(placed_bets$status, placed_bets$bet_hash)
  } else setNames(character(), character())
  placed_ticket_lookup <- if (nrow(placed_bets) > 0 && "ticket_number" %in% names(placed_bets)) {
    setNames(placed_bets$ticket_number, placed_bets$bet_hash)
  } else setNames(character(), character())
  placed_actual_lookup <- setNames(
    if (nrow(placed_bets) > 0 && "actual_size" %in% names(placed_bets)) placed_bets$actual_size else numeric(),
    if (nrow(placed_bets) > 0) placed_bets$bet_hash else character()
  )

  same_game_info <- lapply(seq_len(nrow(all_bets)), function(i) {
    find_same_game_bets(i, all_bets, placed_bets)
  })

  # Re-compute bet_row_id to match the join key in mlb_bets_book_prices
  all_bets <- all_bets %>%
    mutate(bet_row_id = vapply(
      paste(id, market, ifelse(is.na(line), "", as.character(line)), bet_on, sep = "|"),
      function(s) digest::digest(s, algo = "md5"), character(1)
    ))

  # Compute display data for each bet
  table_data <- all_bets %>%
    mutate(
      bet_hash       = pmap_chr(list(id, market, bet_on, line), generate_bet_hash),
      is_placed      = bet_hash %in% placed_hashes,
      matchup        = paste(away_team, "@", home_team),
      tipoff         = ifelse(is.na(pt_start_time), "",
                              format(pt_start_time, "%a %I:%M %p")),
      ev_pct         = ev * 100,
      risk_amt       = bet_size,
      towin_amt      = ifelse(odds > 0, bet_size * odds / 100,
                                        bet_size * 100 / abs(odds))
    ) %>%
    arrange(desc(ev))

  # Partial-fill detection
  placed_actual <- if (nrow(placed_bets) > 0 && "actual_size" %in% names(placed_bets)) {
    setNames(placed_bets$actual_size, placed_bets$bet_hash)
  } else setNames(numeric(), character())

  table_data <- table_data %>%
    mutate(
      placed_actual = ifelse(is_placed, placed_actual[bet_hash], NA_real_),
      fill_status = case_when(
        !is_placed                                   ~ "not_placed",
        is.na(placed_actual)                          ~ "placed",
        round(placed_actual) >= round(bet_size)       ~ "placed",
        TRUE                                          ~ "partial"
      ),
      fill_diff = ifelse(fill_status == "partial",
                         round(bet_size) - round(placed_actual), NA_real_)
    )

  # Render one card per row
  table_data$card_html <- vapply(seq_len(nrow(table_data)), function(i) {
    row <- table_data[i, ]
    bet_id <- row$bet_row_id

    # Look up book prices for both sides
    wide_pick <- book_prices_wide %>% filter(bet_row_id == bet_id, side == "pick")
    wide_opp  <- book_prices_wide %>% filter(bet_row_id == bet_id, side == "opposite")
    pickside_wide_row <- if (nrow(wide_pick) > 0) wide_pick[1, ] else NULL
    oppside_wide_row  <- if (nrow(wide_opp)  > 0) wide_opp[1,  ] else NULL

    # Bet title + market label (e.g., "BOS -2.5" + "Alt Spread · Full Game")
    title_parts <- format_market_for_card(row$market, row$bet_on, row$line,
                                           row$home_team, row$away_team)

    # Pick row label (same as title for spreads/ML; "Over 5.5" for totals)
    pick_label <- if (grepl("^Over|^Under", row$bet_on, ignore.case = TRUE)) {
      line_disp <- if (!is.na(row$line)) {
        if (row$line == round(row$line)) format(as.integer(row$line))
        else sub("\\.?0+$", "", sprintf("%.2f", row$line))
      } else ""
      paste(row$bet_on, line_disp)
    } else if (is.na(row$line)) {
      .team_abbr(row$bet_on)
    } else {
      line_disp <- if (row$line > 0) paste0("+", row$line) else as.character(row$line)
      paste(.team_abbr(row$bet_on), line_disp)
    }

    # Opposite row label
    is_totals_market <- grepl("^(Over|Under)", row$bet_on, ignore.case = TRUE)
    opp_label <- if (is_totals_market) {
      other_side <- if (grepl("^Over", row$bet_on, ignore.case = TRUE)) "Under" else "Over"
      line_disp <- if (!is.na(row$line)) {
        if (row$line == round(row$line)) format(as.integer(row$line))
        else sub("\\.?0+$", "", sprintf("%.2f", row$line))
      } else ""
      paste(other_side, line_disp)
    } else if (is.na(row$line)) {
      # ML — opposite team
      other_team <- if (identical(row$bet_on, row$home_team)) row$away_team else row$home_team
      .team_abbr(other_team)
    } else {
      # Spread — flip team + flip sign
      other_team <- if (identical(row$bet_on, row$home_team)) row$away_team else row$home_team
      flipped <- -row$line
      line_disp <- if (flipped > 0) paste0("+", flipped) else as.character(flipped)
      paste(.team_abbr(other_team), line_disp)
    }

    side_word <- if (grepl("^Over",  row$bet_on, ignore.case = TRUE)) "over"
                 else if (grepl("^Under", row$bet_on, ignore.case = TRUE)) "under"
                 else "over"

    # Action-button HTML
    status        <- placed_status_lookup[row$bet_hash]
    ticket        <- placed_ticket_lookup[row$bet_hash]
    placed_actual <- placed_actual_lookup[row$bet_hash]
    data_attrs <- sprintf(
      'data-hash="%s" data-game-id="%s" data-home="%s" data-away="%s" data-time="%s" data-market="%s" data-bet-on="%s" data-line="%s" data-prob="%s" data-ev="%s" data-size="%s" data-odds="%s" data-book="%s" data-actual="%s" data-fill-status="%s"',
      row$bet_hash, row$id, row$home_team, row$away_team,
      as.character(row$pt_start_time), row$market, row$bet_on,
      ifelse(is.na(row$line), "", row$line),
      row$prob, row$ev, row$bet_size, row$odds, row$bookmaker_key,
      ifelse(is.na(placed_actual), "", placed_actual),
      row$fill_status
    )

    action_html <- if (!is.na(row$fill_status) && row$fill_status == "partial") {
      sprintf('<button class="btn-partial" onclick="updateBet(this)" %s>Partial -$%.0f</button>',
              data_attrs, row$fill_diff)
    } else if (!is.na(status) && status == "placed") {
      label <- if (!is.na(ticket) && nchar(ticket) > 0)
        sprintf("placed &middot; #%s", htmltools::htmlEscape(ticket)) else "placed"
      sprintf('<span class="placed-bet-label" %s>%s</span>', data_attrs, label)
    } else if (!is.na(status) && status %in%
               c("price_moved","rejected","auth_error","network_error","orphaned")) {
      short <- switch(status, price_moved="drift", rejected="rejected",
                              auth_error="auth err", network_error="net err",
                              orphaned="orphan", status)
      sprintf(
        '<span class="pill error" %s>%s</span><button class="btn-place" onclick="placeBet(this)" %s>Retry</button><button class="btn-log" onclick="logBet(this)" %s>Log</button>',
        data_attrs, short, data_attrs, data_attrs)
    } else {
      supported_place <- row$bookmaker_key %in% c("wagerzon","hoop88","bfa","betonlineag")
      place_btn <- if (supported_place) {
        sprintf('<button class="btn-place" onclick="placeBet(this)" %s>Place Bet</button>', data_attrs)
      } else {
        sprintf('<button class="btn-place" disabled title="manual log only for this book" %s>Place Bet</button>', data_attrs)
      }
      log_btn <- sprintf('<button class="btn-log" onclick="logBet(this)" %s>Log</button>', data_attrs)
      paste0(place_btn, " ", log_btn)
    }

    # Same-game-corr badge (re-uses the existing tooltip builder)
    info <- same_game_info[[i]]
    corr_html <- if (info$has_same_game) {
      dot <- intToUtf8(0xB7); chk <- intToUtf8(0x2713); cir <- intToUtf8(0x2013)
      lines_text <- vapply(info$details, function(d) {
        # (Identical body to the legacy one — copied verbatim to avoid drift.)
        market_name <- format_market_name(d$market)
        line_str <- if (!is.null(d$line) && !is.na(d$line)) {
          if (d$line > 0) paste0(" +", d$line) else paste0(" ", d$line)
        } else ""
        odds_str <- if (!is.null(d$odds) && !is.na(d$odds)) {
          if (d$odds > 0) sprintf(" (%+d)", d$odds) else sprintf(" (%d)", d$odds)
        } else ""
        size_str <- if (isTRUE(d$is_placed)) {
          if (!is.null(d$actual_size) && !is.na(d$actual_size))
            sprintf("$%.0f", d$actual_size) else ""
        } else {
          if (!is.null(d$bet_size) && !is.na(d$bet_size))
            sprintf("$%.0f", d$bet_size) else ""
        }
        prefix <- if (isTRUE(d$is_placed)) chk else cir
        book_str <- if (!is.null(d$bookmaker) && !is.na(d$bookmaker)) d$bookmaker else ""
        paste(prefix, paste(Filter(nzchar, c(
          paste0(market_name, " - ", d$bet_on, line_str, odds_str), size_str, book_str
        )), collapse = paste0(" ", dot, " ")))
      }, character(1))
      tooltip <- paste(lines_text, collapse = "\n")
      sprintf('<span class="corr-badge" title="%s">&#9679;</span>',
              escape_tooltip(tooltip))
    } else ""

    hero_html <- render_hero_strip(
      pick_book      = toupper(row$bookmaker_key),
      pick_odds      = row$odds,
      fair_odds      = row$fair_odds,
      ev_pct         = row$ev_pct,
      risk_dollars   = row$risk_amt,
      towin_dollars  = row$towin_amt,
      action_html    = action_html
    )

    render_bet_card(
      bet_title          = title_parts$bet_title,
      market_label       = title_parts$market_label,
      matchup            = row$matchup,
      tipoff             = row$tipoff,
      hero_html          = hero_html,
      pickside_label     = pick_label,
      pickside_wide_row  = pickside_wide_row,
      oppside_label      = opp_label,
      oppside_wide_row   = oppside_wide_row,
      pick_book          = row$bookmaker_key,
      side_word          = side_word,
      is_totals          = is_totals_market,
      corr_badge_html    = corr_html
    )
  }, character(1))

  # Drop columns that aren't needed by reactable to avoid serialization issues
  table_data <- table_data %>%
    select(card_html, bet_hash, ev, market, bet_on, away_team, home_team,
           any_of(c("matchup", "tipoff")))

  reactable(
    table_data,
    elementId = "bets-table",
    searchable = TRUE,
    filterable = TRUE,
    defaultPageSize = 25,
    pageSizeOptions = c(25, 50, 100),
    showPageSizeOptions = TRUE,
    columns = list(
      bet_hash   = colDef(show = FALSE),
      ev         = colDef(show = FALSE),
      market     = colDef(show = FALSE),
      bet_on     = colDef(show = FALSE),
      away_team  = colDef(show = FALSE),
      home_team  = colDef(show = FALSE),
      matchup    = colDef(show = FALSE),
      tipoff     = colDef(show = FALSE),
      card_html  = colDef(
        name = "", html = TRUE,
        cell = function(value, index) value
      )
    ),
    theme = reactableTheme(
      backgroundColor = "transparent",
      borderColor     = "transparent",
      stripedColor    = "transparent",
      highlightColor  = "transparent"
    )
  )
}
```

- [ ] **Step 3: Reload the dashboard**

Stop any running dashboard, then:
```bash
cd "Answer Keys/MLB Dashboard"
Rscript mlb_dashboard.R 2>&1 | head -30
```

Open the bets tab in the browser. Expected: V8 cards visible. Each card has the new structure (bet title H1 + matchup line + green hero strip + 2×8 price grid).

- [ ] **Step 4: Eyeball-compare against the V8 mockup**

The brainstorm server may have stopped (30-min auto-exit). If you want to re-start it for side-by-side comparison:
```bash
/Users/callancapitolo/.claude-personal/plugins/cache/claude-plugins-official/superpowers/5.1.0/skills/brainstorming/scripts/start-server.sh --project-dir /Users/callancapitolo/NFLWork/.claude/worktrees/mlb-bets-tab-improvements
```

The mockup HTML lives at `.superpowers/brainstorm/37414-1778712316/content/bet-card-layouts-v8.html`. Open both side-by-side. Things to check:
- [ ] Bet title and matchup line uniform-size-per-row
- [ ] Pick book + odds visible in green
- [ ] Fair / EV / Risk / To Win centered in their slots
- [ ] Place + Log buttons sit right after To Win
- [ ] Price grid: row labels left, 8 book columns share remaining width
- [ ] Pick cell has green background and green text
- [ ] Alt-line cells (where book's line differs from bet line) show the line tag in amber

- [ ] **Step 5: Commit**

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "$(cat <<'EOF'
feat(mlb-bets-tab): wire V8 card layout into create_bets_table

Replace per-cell reactable colDefs (game / market / pickside /
otherside / m / pick / ev / size / towin / action) with a single
card_html column rendered via render_bet_card() helpers.

The card now matches the V8 mockup:
- Bet title H1 (uniform 18px, color-only hierarchy)
- Matchup line (uniform 15px)
- Green hero strip: WAGERZON +160 | FAIR | EV | RISK | TO WIN | Place / Log
- 2x8 price grid: sides x books, with centered alignment

Same-game-corr badge and partial-fill / status states preserved
unchanged in the action-button area.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 10: End-to-end visual verification + small polish iterations

**Files:**
- Modify (likely): small CSS tweaks in `Answer Keys/MLB Dashboard/mlb_dashboard.R`
- Modify (possibly): copy adjustments in `format_market_for_card` if labels read awkward on real data

**Why:** The spec calls out that "CSS in R Shiny isn't hot-reload-friendly. Plan for ~5–10 visual iterations after the structural change lands."

- [ ] **Step 1: Run the full pipeline against a live slate**

```bash
cd "Answer Keys"
Rscript "MLB Answer Key/MLB.R" 2>&1 | tail -20
```

Expected: no errors. Output mentions `mlb_bets_book_prices` row count.

- [ ] **Step 2: Restart dashboard and walk through the bets tab**

```bash
cd "Answer Keys/MLB Dashboard"
Rscript mlb_dashboard.R
```

Verify all four user-flagged issues are resolved:

- [ ] **(Issue 1) F7 totals card shows a FanDuel price** (not "—") on a card where the model is on an F7 Over/Under. Find an F7 card in the bets tab; FD column should have a price.

- [ ] **(Issue 2) Alt-spread BOS -2.5 card shows FanDuel exact match** (no amber tag, no "-1.5" line tag on the FD cell — just the price). Find an alt-spread card; FD column should show the exact line.

- [ ] **(Issue 3) Spread card renders both sides.** The Phillies @ Red Sox alt-spread card has both a `BOS -2.5` row and a `PHI +2.5` row in the grid.

- [ ] **(Issue 4) Cards feel right.** The hero strip with Fair / EV / Risk / To Win is centered; Place + Log buttons sit immediately after To Win; bet title row is uniform 18px; matchup row is uniform 15px.

- [ ] **Step 3: Polish iteration — gather any mismatches**

For each visual mismatch vs the mockup, note the change needed (e.g., "the matchup row is too close to the bet line", "the EV value isn't perfectly centered", "the place button looks dim"). Apply CSS tweaks in `mlb_dashboard.R`'s CSS block, restart dashboard, re-verify. Repeat until visuals match the mockup.

- [ ] **Step 4: Test placing a bet from the new card**

Click `Place Bet` on a card backed by Wagerzon (WZ — the only book with auto-placement support per row 1488 logic). Verify:
- The button click handler fires (look for `placeBet` JS call).
- The status updates to `placed` and the card shows `placed · #<ticket>`.
- No JS errors in the browser console.

If anything breaks, the data attrs on the button (`data-hash`, `data-game-id`, etc.) are the most likely culprit — they must match what `mlb_dashboard_server.py::/api/place-bet` expects.

- [ ] **Step 5: Test "Log" button**

Click `Log` on any card. Verify the manual placement log endpoint is hit and the status flips to logged.

- [ ] **Step 6: Test search + filter still work**

Use the reactable search box at the top of the bets table. Search for a team name (e.g., "Red Sox"). Confirm cards filter correctly.

- [ ] **Step 7: Commit polish tweaks**

```bash
git add "Answer Keys/MLB Dashboard/mlb_dashboard.R"
git commit -m "style(mlb-bets-tab): polish V8 card visuals after end-to-end check"
```

(Skip this commit if Steps 3+ required no changes.)

---

## Task 11: Update documentation

**Files:**
- Modify: `Answer Keys/MLB Dashboard/README.md`
- Modify: `Answer Keys/CLAUDE.md`
- Modify: `mlb_sgp/README.md`

**Why:** Spec requires docs to be updated in the same merge.

- [ ] **Step 1: Update `Answer Keys/MLB Dashboard/README.md`**

Find the bets-tab section. Update its description so:
- The "pill row" wording is replaced with "price grid" or "2×8 grid".
- The screenshot (if there is one) is refreshed.
- Mention that spreads and ML bets now render both sides.

If there's no current screenshot, no need to add one — the README is reference, not a tutorial.

- [ ] **Step 2: Update `Answer Keys/CLAUDE.md`**

Find the section "MLB Dashboard — Odds screen + WZ single-bet placer". The "Data flow" subsection mentions pills (point 1 in the data flow):
> 1. MLB.R writes `mlb_bets_book_prices` to `mlb_mm.duckdb` alongside `mlb_bets_combined`. Each row is one (bet × book × side) at the model's exact line OR the closest line within ±1 unit.

Update the "±1 unit" to "±3.0 units" (per `LINE_MATCH_TOLERANCE = 3.0` in `odds_screen.R:18`).

In the "Helpers" subsection, update:
> - `Answer Keys/MLB Dashboard/book_pill.R` — `render_book_pill()` HTML helper

to:
> - `Answer Keys/MLB Dashboard/book_cell.R` — `render_book_cell()` grid-cell HTML helper. Replaces the legacy `book_pill.R` as part of the V8 card redesign (2026-05-13).

Add a short note about the V8 layout — one paragraph saying the card is now a strict 2×8 grid with a hero strip surfacing Fair, EV, Risk, To Win.

- [ ] **Step 3: Update `mlb_sgp/README.md`**

Find the FanDuel singles section. Add a line under "Supported markets" (or wherever the market list lives):
> First-7-Innings markets (run line, total runs, money line, alt run lines, alt totals) added to the whitelist 2026-05-13.

- [ ] **Step 4: Commit**

```bash
git add "Answer Keys/MLB Dashboard/README.md" "Answer Keys/CLAUDE.md" mlb_sgp/README.md
git commit -m "$(cat <<'EOF'
docs(mlb-bets-tab): update READMEs + CLAUDE.md for V8 layout + FD F7

- MLB Dashboard README: replace "pill row" with "2x8 price grid"; note
  that spreads + ML now render both sides
- Answer Keys CLAUDE.md: update Odds-screen section to reference
  book_cell.R / render_book_cell() and the V8 grid layout; correct
  the line-match tolerance to 3.0
- mlb_sgp/README: note F7 market coverage added to FD whitelist

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 12: Pre-merge executive review

**Files:**
- Read only — no edits, just a review pass.

**Why:** Project CLAUDE.md requires an executive-engineer review of the full diff before any merge to main. This is the gate.

- [ ] **Step 1: Generate the cumulative diff**

```bash
git diff main..HEAD --stat
git log --oneline main..HEAD
```

Expected: a clean list of commits, each focused. Each commit should compile / not break tests when checked out individually.

- [ ] **Step 2: Run the full review checklist**

Per the project CLAUDE.md pre-merge review checklist:

- [ ] **Data integrity:** no duplicate writes; `mlb_bets_book_prices` row-count doubling for spread bets is intentional (Task 5).
- [ ] **Resource safety:** any new DB connections use `on.exit(dbDisconnect(...))`. (Check by grepping for `dbConnect` in the modified files.)
- [ ] **Edge cases:** first-run with no bets; bets with NA `fair_odds`; bets with no `opposite` data; off-season behavior.
- [ ] **Dead code:** the `render_book_pill` shim in `book_cell.R` is the only one — leave for one release, then remove.
- [ ] **Log / disk hygiene:** no new log files created.
- [ ] **Security:** no secrets in logs; data-attrs on action buttons are HTML-escaped.

- [ ] **Step 3: Run the full test suite**

```bash
cd "Answer Keys" && Rscript -e 'testthat::test_dir("tests")' 2>&1 | tail -15
cd "Answer Keys/MLB Dashboard" && Rscript -e 'testthat::test_dir("tests")' 2>&1 | tail -15
cd mlb_sgp && python -m pytest tests/ -v 2>&1 | tail -15
```

Expected: all PASS, no FAIL.

- [ ] **Step 4: Write a one-page summary of findings**

In chat (not in a file), summarize:
- Issues found: (list any).
- Acceptable risks: (the spec's known risks — row-count doubling, alt-line scope uncertainty).
- Test results: pass/fail counts.

Hand summary back to the user. **Do not merge yet** — wait for explicit user approval per `[[feedback_always_ask_merge]]`.

- [ ] **Step 5: Wait for user approval, then merge**

Once user approves:

```bash
git checkout main
git merge --no-ff worktree-mlb-bets-tab-improvements -m "Merge bets-tab V8 layout + FD F7 + alt-spread fixes"
```

Run the full pipeline one more time on main to confirm nothing broke during the merge:
```bash
cd "Answer Keys"
Rscript "MLB Answer Key/MLB.R" 2>&1 | tail -10
```

Restart the dashboard on its production port (8083) and visually confirm.

- [ ] **Step 6: Clean up the worktree**

```bash
cd /Users/callancapitolo/NFLWork
git worktree remove .claude/worktrees/mlb-bets-tab-improvements
git branch -d worktree-mlb-bets-tab-improvements
git worktree list  # confirm gone
```

---

## Self-Review

**1. Spec coverage:**
- ✓ Workstream 1 (opposite rows for spreads) — Task 5
- ✓ Workstream 2 (dashboard layout) — Tasks 6, 7, 8, 9
- ✓ Workstream 3 (FD F7 whitelist) — Tasks 2, 3
- ✓ Workstream 4 (FD alt-line fix) — Tasks 2, 4
- ✓ Workstream 5 (label & copy polish) — Task 9 (Place Bet copy; "EV" kept; Fair stat added)
- ✓ Documentation — Task 11
- ✓ Worktree lifecycle — Tasks 1 + 12

**2. Placeholder scan:** none. Every step has full code or full commands.

**3. Type consistency:**
- `render_book_cell(american_odds, line_quoted, is_exact_line, is_pick, side_word, is_totals)` — same signature in Task 6 and Task 7. ✓
- `render_hero_strip(pick_book, pick_odds, fair_odds, ev_pct, risk_dollars, towin_dollars, action_html)` — same in Task 7 and Task 9. ✓
- `render_bet_card(bet_title, market_label, matchup, tipoff, hero_html, pickside_label, pickside_wide_row, oppside_label, oppside_wide_row, pick_book, side_word, is_totals, corr_badge_html)` — same in Task 7 and Task 9. ✓
- `BOOK_ORDER_V8` and `BOOK_LABELS_V8` — used in Task 7 (definition) and Task 7 (consumer). ✓
- `.team_abbr` — used in Task 9 (definition + 3 call sites). ✓

**4. Known risks captured in plan:** alt-line fix scope unknown until Task 2 recon; FD F7 may not be posted at all (Task 2 covers that path); CSS visual polish is iterative (Task 10 explicit).
