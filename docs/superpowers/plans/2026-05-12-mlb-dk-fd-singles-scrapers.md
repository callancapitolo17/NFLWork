# MLB DK/FD Single-Leg Scrapers Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace Odds API as the source for DraftKings and FanDuel single-leg pill data on the MLB Dashboard bets tab by adding two new scrapers (`scraper_draftkings_singles.py`, `scraper_fanduel_singles.py`) that share API client classes with the existing SGP scrapers.

**Architecture:** Extract DK/FD API code into per-book client classes (`dk_client.py`, `fd_client.py`). Refactor existing SGP scrapers to use them (behavior-preserving). Build singles scrapers on the same clients. Write to per-book DuckDB files in the offshore 18-column `mlb_odds` schema; MLB.R reads via `scraper_to_canonical()`.

**Tech Stack:** Python 3 + `curl_cffi` (Chrome TLS, Akamai bypass) + DuckDB; R for MLB.R integration; pytest for Python tests.

**Design spec:** `docs/superpowers/specs/2026-05-12-mlb-dk-fd-singles-scrapers-design.md`

---

## Version control & worktree lifecycle

- **Worktree:** `/Users/callancapitolo/NFLWork/.claude/worktrees/feature+mlb-dk-fd-singles-scrapers`
- **Branch:** `worktree-feature+mlb-dk-fd-singles-scrapers`
- **Base:** `main`
- **Commit cadence:** one commit per task (each task ends in a commit step).
- **Pre-merge gate:** all tasks complete + full pipeline smoke + executive review (Task 17).
- **Merge:** explicit user approval before `git checkout main && git merge` (per CLAUDE.md feedback_always_ask_merge).
- **Cleanup after merge:** `git worktree remove <path>` + `git branch -d worktree-feature+mlb-dk-fd-singles-scrapers`.

## Documentation updates (in Task 16, same merge as code)

- `mlb_sgp/README.md` — new "Singles scrapers" section.
- `Answer Keys/CLAUDE.md` — Pipeline Flow diagram updated to show DK/FD in parallel-scrapers tier.
- Project-root `CLAUDE.md` "MLB Dashboard — Odds screen" section — note DK/FD source change.
- `dk_odds/README.md`, `fd_odds/README.md` — new (one-line description + duckdb table reference).
- `mlb_sgp/CLAUDE.md` (if exists) — note client-class refactor.

---

## File structure (locked from design)

**New files:**
```
mlb_sgp/dk_client.py                          # DraftKingsClient + dataclasses
mlb_sgp/fd_client.py                          # FanDuelClient + dataclasses
mlb_sgp/scraper_draftkings_singles.py         # Singles scrape → dk_odds/dk.duckdb
mlb_sgp/scraper_fanduel_singles.py            # Singles scrape → fd_odds/fd.duckdb
mlb_sgp/tests/test_dk_client.py               # Fixture-based parser tests
mlb_sgp/tests/test_fd_client.py
mlb_sgp/tests/test_dk_singles_parser.py       # parse_selections_to_wide_rows tests
mlb_sgp/tests/test_fd_singles_parser.py
mlb_sgp/tests/fixtures/dk_event_markets.json  # Captured API responses
mlb_sgp/tests/fixtures/dk_selections.json
mlb_sgp/tests/fixtures/fd_runners.json
dk_odds/README.md
fd_odds/README.md
```

**Modified files:**
```
mlb_sgp/scraper_draftkings_sgp.py             # Import from dk_client (no behavior change)
mlb_sgp/scraper_fanduel_sgp.py                # Import from fd_client (no behavior change)
mlb_sgp/README.md                             # +Singles scrapers section
Answer Keys/run.py                            # +2 SCRAPER_CONFIGS entries
Answer Keys/MLB Answer Key/MLB.R              # book_odds_by_book DK/FD swap
Answer Keys/CLAUDE.md                         # Pipeline diagram
CLAUDE.md (project root)                      # MLB Dashboard — Odds screen section
```

**Untouched:** All offshore scrapers, dashboard code, Kalshi/Novig/ProphetX SGP scrapers, Pikkit.

---

## Phase 1 — Verify open uncertainties (BLOCKING — must complete before client work)

### Task 1: Verify API surface assumptions

**Goal:** confirm three uncertainties from the design spec before committing to the client class surface. Investigation only — no production code changes.

**Files:**
- Create (transient): `mlb_sgp/tests/fixtures/dk_selections.json`, `mlb_sgp/tests/fixtures/dk_event_markets.json`, `mlb_sgp/tests/fixtures/fd_runners.json`

- [ ] **Step 1: Capture a live DK SGP payload for one event**

```bash
cd mlb_sgp
source venv/bin/activate
python -c "
from scraper_draftkings_sgp import init_session, fetch_dk_events, _fetch_subcat_markets, fetch_main_market_nums, fetch_selection_ids
import json
s = init_session()
events = fetch_dk_events(s)
eid = events[0]['event_id']
main_nums = fetch_main_market_nums(s, eid)
sel_ids = fetch_selection_ids(s, eid, main_nums, verbose=True)
# Dump the raw structure so we can see if prices live in selection rows
with open('tests/fixtures/dk_selections.json', 'w') as f:
    json.dump({'event_id': eid, 'main_nums': main_nums, 'sel_ids': sel_ids}, f, indent=2, default=str)
print('Dumped', len(sel_ids), 'selection groups')
"
```

Expected: a JSON file with selection groups. Inspect to confirm each selection has an `odds`/`displayOdds`/`american_odds` field (or whatever DK calls the price). If prices ARE in the payload → singles scraper is one-call-per-event. If NOT → singles scraper needs to additionally hit `event/eventSubcategory/v1/markets` per event (which `fetch_main_market_nums` already calls).

- [ ] **Step 2: Inspect for F3/F7 markets on DK**

```bash
python -c "
import json
data = json.load(open('tests/fixtures/dk_selections.json'))
names = set()
for grp in data['sel_ids'].values() if isinstance(data['sel_ids'], dict) else data['sel_ids']:
    if isinstance(grp, dict):
        names.add(grp.get('name', ''))
print('Distinct selection names with inning markers:')
for n in sorted(names):
    if any(k in n.lower() for k in ['1st 3', '1st 7', 'f3', 'f7', 'innings']):
        print(' -', n)
"
```

Expected: list of inning-period selection names. If F3 and F7 appear → singles scraper writes those periods. If only F5 → record limitation in `mlb_sgp/README.md` (FD-style "no F3/F7 yet").

- [ ] **Step 3: Capture an FD runners payload + check F3/F7**

```bash
python -c "
from scraper_fanduel_sgp import init_session, fetch_fd_events, fetch_event_runners
import json
s = init_session()
events = fetch_fd_events(s)
eid = events[0]['fd_event_id']
home = events[0].get('home_team', 'home')
away = events[0].get('away_team', 'away')
runners = fetch_event_runners(s, eid, home, away)
with open('tests/fixtures/fd_runners.json', 'w') as f:
    json.dump({'event_id': eid, 'runners': runners}, f, indent=2, default=str)
names = {r.get('runnerName', r.get('name', '')) for r in runners if isinstance(r, dict)} if isinstance(runners, list) else set()
print('Distinct runner names with inning markers:')
for n in sorted(names):
    if any(k in n.lower() for k in ['1st 3', '1st 7', 'f3', 'f7', 'innings']):
        print(' -', n)
"
```

Expected: list of F-period runner names from FD. If empty → FD has no F3/F7; record in design spec's open-uncertainty resolution.

- [ ] **Step 4: Verify DK/FD team-name canonicalization**

```bash
python -c "
from scraper_draftkings_sgp import fetch_dk_events, init_session
from scraper_fanduel_sgp import fetch_fd_events
import duckdb

con = duckdb.connect('../Answer Keys/mlb.duckdb', read_only=True)
canonical = {row[0] for row in con.execute('SELECT DISTINCT home_team FROM mlb_consensus_temp UNION SELECT DISTINCT away_team FROM mlb_consensus_temp').fetchall()}
con.close()

s = init_session()
dk_teams = set()
for e in fetch_dk_events(s):
    dk_teams.add(e.get('home_team', ''))
    dk_teams.add(e.get('away_team', ''))
print('DK teams NOT in canonical:', dk_teams - canonical)

fd_teams = set()
for e in fetch_fd_events(s):
    fd_teams.add(e.get('home_team', ''))
    fd_teams.add(e.get('away_team', ''))
print('FD teams NOT in canonical:', fd_teams - canonical)
"
```

Expected: empty sets, or a short list of mismatches. Any mismatches go into the team-name dictionary in Task 5 / Task 7.

- [ ] **Step 5: Document findings**

Append a short "Phase 1 findings" block to the design spec:

```bash
cat >> docs/superpowers/specs/2026-05-12-mlb-dk-fd-singles-scrapers-design.md <<'EOF'

---

## Phase 1 findings (2026-05-12)

- DK SGP payload includes single-leg prices: <YES/NO + price field name>
- DK F3/F7 markets: <PRESENT/ABSENT — list specific selection names if present>
- FD F3/F7 markets: <PRESENT/ABSENT>
- DK team-name drift: <NONE / list mismatches>
- FD team-name drift: <NONE / list mismatches>
EOF
```

Replace `<...>` placeholders with actual findings from steps 1-4.

- [ ] **Step 6: Commit findings**

```bash
cd ../..
git add docs/superpowers/specs/2026-05-12-mlb-dk-fd-singles-scrapers-design.md mlb_sgp/tests/fixtures/
git commit -m "investigation: verify DK/FD API surface for singles scrapers

Findings appended to design spec. Fixtures captured for parser unit tests.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Phase 2 — Build dk_client.py with tests

### Task 2: Create dk_client.py skeleton with dataclasses

**Files:**
- Create: `mlb_sgp/dk_client.py`
- Create: `mlb_sgp/tests/test_dk_client.py`

- [ ] **Step 1: Write the failing test**

`mlb_sgp/tests/test_dk_client.py`:

```python
"""Unit tests for DraftKingsClient — fixture-based, no network."""
import json
from pathlib import Path
import pytest
from mlb_sgp.dk_client import DraftKingsClient, Event, Market, Selection

FIXTURES = Path(__file__).parent / "fixtures"


def test_event_dataclass_has_required_fields():
    e = Event(event_id="123", home_team="Yankees", away_team="Red Sox",
              start_time="2026-05-12T22:00:00Z")
    assert e.event_id == "123"
    assert e.home_team == "Yankees"
    assert e.away_team == "Red Sox"
    assert e.start_time == "2026-05-12T22:00:00Z"


def test_selection_dataclass_has_required_fields():
    s = Selection(selection_id="0HC1234N150_1", market_id="1234",
                  name="Yankees -1.5", line=-1.5, american_odds=120)
    assert s.line == -1.5
    assert s.american_odds == 120


def test_market_dataclass_has_required_fields():
    m = Market(market_id="1234", name="Run Line", subcategory="game_lines")
    assert m.subcategory == "game_lines"
```

- [ ] **Step 2: Run test to verify it fails**

```bash
cd mlb_sgp
source venv/bin/activate
pip install pytest  # if not already
python -m pytest tests/test_dk_client.py -v
```

Expected: FAIL with `ImportError: cannot import name 'DraftKingsClient' from 'mlb_sgp.dk_client'`

- [ ] **Step 3: Write minimal dk_client.py skeleton**

`mlb_sgp/dk_client.py`:

```python
"""DraftKings API client — extracted from scraper_draftkings_sgp.py.

Owns the curl_cffi Chrome-TLS session and exposes the three operations
both the SGP scraper and the singles scraper need:
  - list_events()           — all MLB events today
  - fetch_event_markets()   — market metadata (FG vs F5 vs alt)
  - fetch_event_selections()— all selections with prices, in one call
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Any


@dataclass
class Event:
    event_id: str
    home_team: str
    away_team: str
    start_time: str  # ISO UTC string


@dataclass
class Market:
    market_id: str
    name: str
    subcategory: str  # "game_lines" | "alt_lines" | "innings" etc.


@dataclass
class Selection:
    selection_id: str
    market_id: str
    name: str            # e.g. "Yankees -1.5"
    line: float | None   # spread/total value; None for moneylines
    american_odds: int


class DraftKingsClient:
    def __init__(self, verbose: bool = False) -> None:
        from scraper_draftkings_sgp import init_session
        self.session = init_session()
        self.verbose = verbose

    def list_events(self) -> list[Event]:
        raise NotImplementedError("Task 3")

    def fetch_event_markets(self, event_id: str) -> list[Market]:
        raise NotImplementedError("Task 4")

    def fetch_event_selections(self, event_id: str) -> list[Selection]:
        raise NotImplementedError("Task 4")
```

- [ ] **Step 4: Run test to verify it passes**

```bash
python -m pytest tests/test_dk_client.py -v
```

Expected: 3 passed.

- [ ] **Step 5: Commit**

```bash
cd ../..
git add mlb_sgp/dk_client.py mlb_sgp/tests/test_dk_client.py
git commit -m "feat(mlb_sgp): scaffold DraftKingsClient + dataclasses

Empty methods raise NotImplementedError; tests cover dataclass shape only.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

### Task 3: Implement DraftKingsClient.list_events()

**Files:**
- Modify: `mlb_sgp/dk_client.py`
- Modify: `mlb_sgp/tests/test_dk_client.py`

- [ ] **Step 1: Write the failing test**

Append to `mlb_sgp/tests/test_dk_client.py`:

```python
def test_list_events_parses_fixture(monkeypatch):
    """Verifies list_events parses a captured DK league API response."""
    # Fixture captured in Task 1. Load it and monkeypatch the HTTP layer.
    import json
    from mlb_sgp import dk_client

    fixture_path = FIXTURES / "dk_league_response.json"
    with open(fixture_path) as f:
        fake_response_body = json.load(f)

    class FakeResponse:
        status_code = 200
        def json(self):
            return fake_response_body

    class FakeSession:
        def get(self, url, **kwargs):
            return FakeResponse()

    client = DraftKingsClient.__new__(DraftKingsClient)  # bypass __init__
    client.session = FakeSession()
    client.verbose = False

    events = client.list_events()
    assert len(events) > 0
    assert all(isinstance(e, Event) for e in events)
    assert all(e.event_id and e.home_team and e.away_team for e in events)
```

Also create the fixture file by extracting the league response from the existing scraper. In Phase 1, dump it explicitly:

```bash
cd mlb_sgp
python -c "
from scraper_draftkings_sgp import init_session
import json, curl_cffi.requests as cr
s = init_session()
url = 'https://sportsbook-nash.draftkings.com/api/sportscontent/dkusnj/v1/controldata/league/leagueSubcategory/v1/markets?leagueIds=84240&subCategoryIds=4519'
r = s.get(url)
with open('tests/fixtures/dk_league_response.json', 'w') as f:
    json.dump(r.json(), f, indent=2)
"
```

- [ ] **Step 2: Run test to verify it fails**

```bash
python -m pytest tests/test_dk_client.py::test_list_events_parses_fixture -v
```

Expected: FAIL with `NotImplementedError("Task 3")`.

- [ ] **Step 3: Implement list_events**

Replace `list_events` body in `dk_client.py`. Use the existing `fetch_dk_events` logic from `scraper_draftkings_sgp.py` lines 93–129 (read its content with `sed -n '93,129p' mlb_sgp/scraper_draftkings_sgp.py` to see the current implementation), then wrap each dict result in an `Event` dataclass:

```python
    def list_events(self) -> list[Event]:
        """Returns one Event per MLB game today."""
        from scraper_draftkings_sgp import fetch_dk_events
        raw = fetch_dk_events(self.session)
        return [Event(
            event_id=str(e["event_id"]),
            home_team=e["home_team"],
            away_team=e["away_team"],
            start_time=e.get("start_time", "")
        ) for e in raw]
```

Note: this is a thin wrapper. Task 6 (SGP refactor) will eventually move the *logic* of `fetch_dk_events` into `dk_client.py` itself and `scraper_draftkings_sgp.py` will import from `dk_client`. For now we wrap to keep behavior identical and unit-test against a fixture.

- [ ] **Step 4: Run test to verify it passes**

```bash
python -m pytest tests/test_dk_client.py::test_list_events_parses_fixture -v
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
cd ../..
git add mlb_sgp/dk_client.py mlb_sgp/tests/test_dk_client.py mlb_sgp/tests/fixtures/dk_league_response.json
git commit -m "feat(mlb_sgp): DraftKingsClient.list_events with fixture test

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

### Task 4: Implement fetch_event_markets + fetch_event_selections

**Files:**
- Modify: `mlb_sgp/dk_client.py`
- Modify: `mlb_sgp/tests/test_dk_client.py`

- [ ] **Step 1: Write the failing tests**

Append to `mlb_sgp/tests/test_dk_client.py`:

```python
def test_fetch_event_markets_parses_fixture():
    """Verifies fetch_event_markets parses a captured DK markets response."""
    import json
    from mlb_sgp.dk_client import DraftKingsClient, Market

    with open(FIXTURES / "dk_event_markets.json") as f:
        fixture = json.load(f)

    class FakeResponse:
        status_code = 200
        def __init__(self, body):
            self._body = body
        def json(self):
            return self._body

    class FakeSession:
        def get(self, url, **kwargs):
            return FakeResponse(fixture)

    client = DraftKingsClient.__new__(DraftKingsClient)
    client.session = FakeSession()
    client.verbose = False

    markets = client.fetch_event_markets(event_id="dummy")
    assert len(markets) > 0
    assert all(isinstance(m, Market) for m in markets)
    assert all(m.market_id and m.name for m in markets)


def test_fetch_event_selections_returns_priced_rows():
    """Verifies selections include american_odds parsed correctly."""
    import json
    from mlb_sgp.dk_client import DraftKingsClient, Selection

    with open(FIXTURES / "dk_selections.json") as f:
        fixture = json.load(f)

    class FakeResponse:
        status_code = 200
        def json(self):
            return fixture

    class FakeSession:
        def get(self, url, **kwargs):
            return FakeResponse()

    client = DraftKingsClient.__new__(DraftKingsClient)
    client.session = FakeSession()
    client.verbose = False

    selections = client.fetch_event_selections(event_id="dummy")
    assert len(selections) > 0
    assert all(isinstance(s, Selection) for s in selections)
    # Every selection must have an integer american_odds
    assert all(isinstance(s.american_odds, int) for s in selections)
```

- [ ] **Step 2: Capture fixture if not already present**

```bash
cd mlb_sgp
ls tests/fixtures/dk_event_markets.json tests/fixtures/dk_selections.json
# Should exist from Task 1; if not, re-run the Task 1 Step 1 / Step 2 capture scripts.
```

- [ ] **Step 3: Run tests to verify failures**

```bash
python -m pytest tests/test_dk_client.py::test_fetch_event_markets_parses_fixture tests/test_dk_client.py::test_fetch_event_selections_returns_priced_rows -v
```

Expected: 2 failures with `NotImplementedError("Task 4")`.

- [ ] **Step 4: Implement both methods**

Use Phase 1's findings to know exactly where prices live in the payload. Replace bodies in `mlb_sgp/dk_client.py`:

```python
    def fetch_event_markets(self, event_id: str) -> list[Market]:
        """Returns market metadata for one event (game lines + alts)."""
        from scraper_draftkings_sgp import _fetch_subcat_markets, fetch_main_market_nums
        # 4519 = main game lines, 15628 = alt lines.
        # See scraper_draftkings_sgp.py lines 295-315 for subcategory IDs.
        results = []
        for subcat_id, subcat_name in [("4519", "game_lines"), ("15628", "alt_lines")]:
            for m_id, m_name in _fetch_subcat_markets(self.session, event_id, subcat_id):
                results.append(Market(market_id=str(m_id), name=m_name,
                                      subcategory=subcat_name))
        return results

    def fetch_event_selections(self, event_id: str) -> list[Selection]:
        """Returns all selections with their underlying single-leg prices."""
        from scraper_draftkings_sgp import fetch_main_market_nums, fetch_selection_ids
        main_nums = fetch_main_market_nums(self.session, event_id)
        sel_groups = fetch_selection_ids(self.session, event_id, main_nums,
                                          verbose=self.verbose)
        # sel_groups is a structured dict — flatten to Selection rows.
        # Exact shape depends on Phase 1 finding. Adjust this parser to match.
        selections = []
        for group in (sel_groups if isinstance(sel_groups, list) else sel_groups.values()):
            if not isinstance(group, dict):
                continue
            selections.append(Selection(
                selection_id=str(group["selection_id"]),
                market_id=str(group.get("market_id", "")),
                name=group.get("name", ""),
                line=group.get("line"),
                american_odds=int(group["american_odds"]),
            ))
        return selections
```

**IMPORTANT:** Adjust the field names in the parser to match Phase 1's findings — if DK calls it `displayOdds` not `american_odds`, change accordingly. The test's `american_odds` int check will catch a mismatch immediately.

- [ ] **Step 5: Run tests to verify pass**

```bash
python -m pytest tests/test_dk_client.py -v
```

Expected: all 5 tests pass.

- [ ] **Step 6: Commit**

```bash
cd ../..
git add mlb_sgp/dk_client.py mlb_sgp/tests/test_dk_client.py
git commit -m "feat(mlb_sgp): DraftKingsClient market + selection fetchers

Tests assert dataclass shape and integer price parsing using captured fixtures.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Phase 3 — Build fd_client.py with tests

### Task 5: Create fd_client.py (mirror of dk_client)

**Files:**
- Create: `mlb_sgp/fd_client.py`
- Create: `mlb_sgp/tests/test_fd_client.py`

- [ ] **Step 1: Write the failing tests**

`mlb_sgp/tests/test_fd_client.py`:

```python
"""Unit tests for FanDuelClient — fixture-based."""
import json
from pathlib import Path
import pytest
from mlb_sgp.fd_client import FanDuelClient, Event, Runner

FIXTURES = Path(__file__).parent / "fixtures"


def test_event_dataclass_has_required_fields():
    e = Event(event_id="abc", home_team="Yankees", away_team="Red Sox",
              start_time="2026-05-12T22:00:00Z")
    assert e.event_id == "abc"


def test_runner_dataclass_has_required_fields():
    r = Runner(runner_id="42", market_id="m1", name="Over 9.5",
               line=9.5, american_odds=-110)
    assert r.american_odds == -110


def test_list_events_parses_fixture():
    with open(FIXTURES / "fd_events.json") as f:
        fixture = json.load(f)

    class FakeResponse:
        status_code = 200
        def json(self): return fixture

    class FakeSession:
        def get(self, url, **kwargs): return FakeResponse()

    client = FanDuelClient.__new__(FanDuelClient)
    client.session = FakeSession()
    client.verbose = False

    events = client.list_events()
    assert len(events) > 0
    assert all(isinstance(e, Event) for e in events)


def test_fetch_event_runners_parses_fixture():
    with open(FIXTURES / "fd_runners.json") as f:
        fixture = json.load(f)

    class FakeResponse:
        status_code = 200
        def json(self): return fixture

    class FakeSession:
        def get(self, url, **kwargs): return FakeResponse()

    client = FanDuelClient.__new__(FanDuelClient)
    client.session = FakeSession()
    client.verbose = False

    runners = client.fetch_event_runners(event_id="dummy",
                                          home_team="Yankees",
                                          away_team="Red Sox")
    assert len(runners) > 0
    assert all(isinstance(r, Runner) for r in runners)
    assert all(isinstance(r.american_odds, int) for r in runners)
```

- [ ] **Step 2: Capture fd_events.json fixture**

```bash
cd mlb_sgp
source venv/bin/activate
python -c "
from scraper_fanduel_sgp import init_session
import json
s = init_session()
# Same URL the scraper uses — see fetch_fd_events implementation
# in scraper_fanduel_sgp.py around line 124.
r = s.get('https://sbapi.nj.sportsbook.fanduel.com/api/content-managed-page?page=CUSTOM&customPageId=mlb')  # adjust per scraper
with open('tests/fixtures/fd_events.json', 'w') as f:
    json.dump(r.json(), f, indent=2)
"
```

If the URL above doesn't work, copy the exact URL from `scraper_fanduel_sgp.py:fetch_fd_events` (read with `sed -n '124,180p' mlb_sgp/scraper_fanduel_sgp.py`).

- [ ] **Step 3: Run tests to verify failures**

```bash
python -m pytest tests/test_fd_client.py -v
```

Expected: all fail with `ImportError`.

- [ ] **Step 4: Write fd_client.py**

`mlb_sgp/fd_client.py`:

```python
"""FanDuel API client — extracted from scraper_fanduel_sgp.py."""
from __future__ import annotations
from dataclasses import dataclass


@dataclass
class Event:
    event_id: str
    home_team: str
    away_team: str
    start_time: str


@dataclass
class Runner:
    runner_id: str
    market_id: str
    name: str
    line: float | None
    american_odds: int


class FanDuelClient:
    def __init__(self, verbose: bool = False) -> None:
        from scraper_fanduel_sgp import init_session
        self.session = init_session()
        self.verbose = verbose

    def list_events(self) -> list[Event]:
        from scraper_fanduel_sgp import fetch_fd_events
        raw = fetch_fd_events(self.session)
        return [Event(
            event_id=str(e["fd_event_id"]),
            home_team=e["home_team"],
            away_team=e["away_team"],
            start_time=e.get("start_time", "")
        ) for e in raw]

    def fetch_event_runners(self, event_id: str, home_team: str,
                             away_team: str) -> list[Runner]:
        from scraper_fanduel_sgp import fetch_event_runners as _legacy
        raw = _legacy(self.session, event_id, home_team, away_team)
        return [Runner(
            runner_id=str(r["runnerId"]) if isinstance(r, dict) else str(r),
            market_id=str(r.get("marketId", "")) if isinstance(r, dict) else "",
            name=r.get("runnerName", "") if isinstance(r, dict) else "",
            line=r.get("handicap") if isinstance(r, dict) else None,
            american_odds=int(r.get("americanOdds", 0)) if isinstance(r, dict) else 0,
        ) for r in (raw if isinstance(raw, list) else [])]
```

Adjust field names per Phase 1 findings (`runnerId` / `runnerName` / `handicap` / `americanOdds` may have different names; the test will catch a mismatch).

- [ ] **Step 5: Run tests to verify pass**

```bash
python -m pytest tests/test_fd_client.py -v
```

Expected: 4 passed.

- [ ] **Step 6: Commit**

```bash
cd ../..
git add mlb_sgp/fd_client.py mlb_sgp/tests/test_fd_client.py mlb_sgp/tests/fixtures/fd_events.json
git commit -m "feat(mlb_sgp): FanDuelClient with list_events + fetch_event_runners

Mirrors DraftKingsClient surface; fixture-based unit tests cover parser shape.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Phase 4 — Refactor existing SGP scrapers (regression-tested)

### Task 6: Capture SGP golden output then refactor DK SGP scraper

**Files:**
- Create (transient): `mlb_sgp/tests/golden/dk_sgp_baseline.csv`
- Modify: `mlb_sgp/scraper_draftkings_sgp.py`

- [ ] **Step 1: Capture golden baseline from current scraper**

```bash
cd mlb_sgp
source venv/bin/activate
python scraper_draftkings_sgp.py 2>&1 | tee logs/dk_sgp_baseline_run.log
# Then export the table that was just written
python -c "
import duckdb
con = duckdb.connect('../Answer Keys/mlb_mm.duckdb', read_only=True)
con.execute('''
  COPY (
    SELECT game_id, combo_name, sgp_decimal, fetch_time
    FROM mlb_sgp_odds
    WHERE book='draftkings'
    ORDER BY game_id, combo_name
  ) TO 'tests/golden/dk_sgp_baseline.csv' (HEADER, DELIMITER ',')
''')
con.close()
print('Baseline captured')
"
```

If `mlb_sgp_odds` doesn't exist (no recent SGP run), trigger one: see `mlb_sgp/README.md` for the manual run command.

- [ ] **Step 2: Write the regression test**

`mlb_sgp/tests/test_sgp_regression.py`:

```python
"""Regression test: refactored DK SGP scraper produces same output as baseline.

The baseline CSV was captured BEFORE the dk_client refactor. After refactor,
running the scraper against the same fixture inputs must produce byte-identical
sgp_decimal values (within rounding tolerance).
"""
import pytest
import duckdb
import csv
from pathlib import Path


@pytest.mark.integration
def test_dk_sgp_output_matches_baseline():
    """Run after refactoring scraper_draftkings_sgp.py. Live API; tagged integration."""
    import subprocess
    # Re-run the scraper after refactor
    subprocess.run(["python", "scraper_draftkings_sgp.py"], check=True)

    con = duckdb.connect("../Answer Keys/mlb_mm.duckdb", read_only=True)
    current = {
        (row[0], row[1]): float(row[2])
        for row in con.execute(
            "SELECT game_id, combo_name, sgp_decimal FROM mlb_sgp_odds WHERE book='draftkings'"
        ).fetchall()
    }
    con.close()

    baseline = {}
    with open("tests/golden/dk_sgp_baseline.csv") as f:
        reader = csv.DictReader(f)
        for row in reader:
            baseline[(row["game_id"], row["combo_name"])] = float(row["sgp_decimal"])

    missing = set(baseline) - set(current)
    extra = set(current) - set(baseline)
    drift = {
        k: (baseline[k], current[k])
        for k in baseline.keys() & current.keys()
        if abs(baseline[k] - current[k]) > 0.005  # half a cent on decimal odds
    }
    assert not missing, f"Combos missing after refactor: {missing}"
    assert not extra, f"Unexpected combos after refactor: {extra}"
    assert not drift, f"sgp_decimal drift after refactor: {drift}"
```

- [ ] **Step 3: Run the regression test to confirm it passes on UNREFACTORED scraper**

```bash
python -m pytest tests/test_sgp_regression.py -v --run-integration
```

Note: if `--run-integration` isn't already a pytest marker config, skip this step and just run with `-v -m integration`. The point is to confirm the baseline equals current state BEFORE we touch the scraper.

Expected: PASS (no refactor yet → identical).

- [ ] **Step 4: Refactor scraper_draftkings_sgp.py to import from dk_client**

Edit `mlb_sgp/scraper_draftkings_sgp.py`:

1. At the top, add: `from dk_client import DraftKingsClient`
2. In `scrape_dk_sgp(verbose=False)` (line 623), replace:
   ```python
   session = init_session()
   dk_events = fetch_dk_events(session)
   ```
   with:
   ```python
   client = DraftKingsClient(verbose=verbose)
   session = client.session
   dk_events = [
       {"event_id": e.event_id, "home_team": e.home_team,
        "away_team": e.away_team, "start_time": e.start_time}
       for e in client.list_events()
   ]
   ```
3. Leave `init_session`, `fetch_dk_events`, `_fetch_subcat_markets`, `fetch_main_market_nums`, `fetch_selection_ids` defined in this file — `dk_client.py` imports them. (They become legacy interfaces used by the client wrapper.)

This is the minimal refactor for Phase 4 — it lets the singles scraper use the client without disturbing SGP behavior. A larger refactor (moving the *implementation* of these functions into dk_client.py) is intentionally deferred to keep this PR small.

- [ ] **Step 5: Run regression test against refactored scraper**

```bash
python -m pytest tests/test_sgp_regression.py -v -m integration
```

Expected: PASS — output byte-identical (or within 0.005 decimal-odds tolerance).

If FAIL: revert the refactor (`git checkout mlb_sgp/scraper_draftkings_sgp.py`), inspect drift, and figure out what diverged before re-attempting.

- [ ] **Step 6: Commit**

```bash
cd ../..
git add mlb_sgp/scraper_draftkings_sgp.py mlb_sgp/tests/test_sgp_regression.py mlb_sgp/tests/golden/dk_sgp_baseline.csv
git commit -m "refactor(mlb_sgp): scraper_draftkings_sgp uses DraftKingsClient

Pulls event list through dk_client wrapper. Other helpers stay in this file
(client imports them). Regression test confirms sgp_decimal output identical
to pre-refactor baseline within 0.005.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

### Task 7: Refactor scraper_fanduel_sgp.py + regression test

**Files:**
- Create (transient): `mlb_sgp/tests/golden/fd_sgp_baseline.csv`
- Modify: `mlb_sgp/scraper_fanduel_sgp.py`
- Modify: `mlb_sgp/tests/test_sgp_regression.py`

- [ ] **Step 1: Capture FD baseline**

```bash
cd mlb_sgp
python scraper_fanduel_sgp.py 2>&1 | tee logs/fd_sgp_baseline_run.log
python -c "
import duckdb
con = duckdb.connect('../Answer Keys/mlb_mm.duckdb', read_only=True)
con.execute('''
  COPY (
    SELECT game_id, combo_name, sgp_decimal, fetch_time
    FROM mlb_sgp_odds
    WHERE book='fanduel'
    ORDER BY game_id, combo_name
  ) TO 'tests/golden/fd_sgp_baseline.csv' (HEADER, DELIMITER ',')
''')
con.close()
"
```

- [ ] **Step 2: Append FD regression test**

Append to `mlb_sgp/tests/test_sgp_regression.py`:

```python
@pytest.mark.integration
def test_fd_sgp_output_matches_baseline():
    import subprocess, csv
    subprocess.run(["python", "scraper_fanduel_sgp.py"], check=True)

    con = duckdb.connect("../Answer Keys/mlb_mm.duckdb", read_only=True)
    current = {
        (row[0], row[1]): float(row[2])
        for row in con.execute(
            "SELECT game_id, combo_name, sgp_decimal FROM mlb_sgp_odds WHERE book='fanduel'"
        ).fetchall()
    }
    con.close()

    baseline = {}
    with open("tests/golden/fd_sgp_baseline.csv") as f:
        for row in csv.DictReader(f):
            baseline[(row["game_id"], row["combo_name"])] = float(row["sgp_decimal"])

    missing = set(baseline) - set(current)
    extra = set(current) - set(baseline)
    drift = {k: (baseline[k], current[k])
             for k in baseline.keys() & current.keys()
             if abs(baseline[k] - current[k]) > 0.005}
    assert not missing, f"FD combos missing after refactor: {missing}"
    assert not extra, f"Unexpected FD combos: {extra}"
    assert not drift, f"FD sgp_decimal drift: {drift}"
```

- [ ] **Step 3: Refactor scraper_fanduel_sgp.py**

In `mlb_sgp/scraper_fanduel_sgp.py`:

1. Add at top: `from fd_client import FanDuelClient`
2. In `scrape_fd_sgp(verbose=False)` (line 555), replace:
   ```python
   session = init_session()
   fd_events = fetch_fd_events(session)
   ```
   with:
   ```python
   client = FanDuelClient(verbose=verbose)
   session = client.session
   fd_events = [
       {"fd_event_id": e.event_id, "home_team": e.home_team,
        "away_team": e.away_team, "start_time": e.start_time}
       for e in client.list_events()
   ]
   ```
3. Leave existing helper functions in place (client imports them).

- [ ] **Step 4: Run FD regression**

```bash
python -m pytest tests/test_sgp_regression.py::test_fd_sgp_output_matches_baseline -v -m integration
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
cd ../..
git add mlb_sgp/scraper_fanduel_sgp.py mlb_sgp/tests/test_sgp_regression.py mlb_sgp/tests/golden/fd_sgp_baseline.csv
git commit -m "refactor(mlb_sgp): scraper_fanduel_sgp uses FanDuelClient

Regression test confirms FD SGP output identical to baseline.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Phase 5 — Singles scrapers

### Task 8: parse_selections_to_wide_rows for DK (TDD)

**Files:**
- Create: `mlb_sgp/scraper_draftkings_singles.py` (parser only — main flow in Task 9)
- Create: `mlb_sgp/tests/test_dk_singles_parser.py`

- [ ] **Step 1: Write the failing test**

`mlb_sgp/tests/test_dk_singles_parser.py`:

```python
"""Tests parse_selections_to_wide_rows — converts DK selections to mlb_odds rows."""
from datetime import datetime
from mlb_sgp.dk_client import Event, Selection
from mlb_sgp.scraper_draftkings_singles import parse_selections_to_wide_rows


def test_main_game_lines_produce_one_row():
    """FG main spread + total + ML for one game → one row marked 'main'."""
    event = Event(event_id="e1", home_team="Yankees", away_team="Red Sox",
                  start_time="2026-05-12T22:00:00Z")
    selections = [
        Selection("s1", "m_main_spread", "Yankees -1.5", -1.5, -120),
        Selection("s2", "m_main_spread", "Red Sox +1.5", 1.5, 110),
        Selection("s3", "m_main_total", "Over 9.5", 9.5, -105),
        Selection("s4", "m_main_total", "Under 9.5", 9.5, -115),
        Selection("s5", "m_ml", "Yankees", None, -150),
        Selection("s6", "m_ml", "Red Sox", None, 130),
    ]
    market_meta = {  # market_id -> (period, market_type)
        "m_main_spread": ("FG", "main"),
        "m_main_total":  ("FG", "main"),
        "m_ml":          ("FG", "main"),
    }
    rows = parse_selections_to_wide_rows(event, selections, market_meta,
                                          fetch_time=datetime(2026, 5, 12, 14, 0))
    assert len(rows) == 1
    r = rows[0]
    assert r["period"] == "FG"
    assert r["market"] == "main"
    assert r["home_spread"] == -1.5
    assert r["home_spread_price"] == -120
    assert r["away_spread"] == 1.5
    assert r["total"] == 9.5
    assert r["over_price"] == -105
    assert r["home_ml"] == -150
    assert r["away_ml"] == 130


def test_alternate_spreads_emit_separate_rows():
    """Each alt-spread line gets its own row marked 'alternate_spreads'."""
    event = Event("e1", "Yankees", "Red Sox", "2026-05-12T22:00:00Z")
    selections = [
        Selection("s1", "m_alt_spread", "Yankees -2.5", -2.5, 120),
        Selection("s2", "m_alt_spread", "Red Sox +2.5", 2.5, -140),
        Selection("s3", "m_alt_spread", "Yankees -0.5", -0.5, -250),
        Selection("s4", "m_alt_spread", "Red Sox +0.5", 0.5, 200),
    ]
    market_meta = {"m_alt_spread": ("FG", "alternate_spreads")}
    rows = parse_selections_to_wide_rows(event, selections, market_meta,
                                          fetch_time=datetime(2026, 5, 12, 14, 0))
    # Two distinct lines → two rows
    assert len(rows) == 2
    lines = sorted([r["home_spread"] for r in rows])
    assert lines == [-2.5, -0.5]
    assert all(r["market"] == "alternate_spreads" for r in rows)
    assert all(r["period"] == "FG" for r in rows)


def test_f5_period_isolated_from_fg():
    """F5 markets live on different rows than FG."""
    event = Event("e1", "Yankees", "Red Sox", "2026-05-12T22:00:00Z")
    selections = [
        Selection("s1", "m_fg_total",  "Over 9.5",  9.5,  -110),
        Selection("s2", "m_fg_total",  "Under 9.5", 9.5,  -110),
        Selection("s3", "m_f5_total",  "Over 4.5",  4.5,  -115),
        Selection("s4", "m_f5_total",  "Under 4.5", 4.5,  -105),
    ]
    market_meta = {
        "m_fg_total": ("FG", "main"),
        "m_f5_total": ("F5", "main"),
    }
    rows = parse_selections_to_wide_rows(event, selections, market_meta,
                                          fetch_time=datetime(2026, 5, 12, 14, 0))
    periods = sorted([r["period"] for r in rows])
    assert periods == ["F5", "FG"]
```

- [ ] **Step 2: Create stub parser file**

`mlb_sgp/scraper_draftkings_singles.py`:

```python
"""DraftKings single-leg odds scraper.

Walks every MLB event today, fetches all selections + market metadata via
DraftKingsClient, transforms to the offshore mlb_odds schema (wagerzon-style),
and writes to dk_odds/dk.duckdb.

Runs as part of run.py mlb orchestrator (pre-MLB.R, parallel scrape phase).
"""
from __future__ import annotations
import argparse
from datetime import datetime
from pathlib import Path
from typing import Any
import duckdb

from dk_client import DraftKingsClient, Event, Selection


def parse_selections_to_wide_rows(
    event: Event,
    selections: list[Selection],
    market_meta: dict[str, tuple[str, str]],   # market_id -> (period, market_type)
    fetch_time: datetime,
) -> list[dict[str, Any]]:
    """Convert flat selections list to wide rows matching offshore mlb_odds schema."""
    raise NotImplementedError("Task 8")
```

- [ ] **Step 3: Run tests to confirm failure**

```bash
cd mlb_sgp
python -m pytest tests/test_dk_singles_parser.py -v
```

Expected: 3 failures.

- [ ] **Step 4: Implement parse_selections_to_wide_rows**

Replace the function in `scraper_draftkings_singles.py`:

```python
def parse_selections_to_wide_rows(
    event: Event,
    selections: list[Selection],
    market_meta: dict[str, tuple[str, str]],
    fetch_time: datetime,
) -> list[dict[str, Any]]:
    """Group selections by (period, market_type, line), emit wide rows.

    Output matches wagerzon.duckdb::mlb_odds schema. See design spec.
    """
    # Group: (period, market_type, line) -> {field: value}
    # For 'main' rows we coalesce spread+total+ml into a single row per period.
    # For 'alternate_spreads' / 'alternate_totals' we emit one row per distinct line.
    buckets: dict[tuple[str, str, float | None], dict[str, Any]] = {}

    def _row_skeleton(period: str, market_type: str, line: float | None) -> dict:
        return {
            "fetch_time": fetch_time,
            "sport_key": "baseball_mlb",
            "game_id": event.event_id,
            "game_date": fetch_time.strftime("%Y-%m-%d"),
            "game_time": event.start_time,
            "away_team": event.away_team,
            "home_team": event.home_team,
            "market": market_type,
            "period": period,
            "away_spread": None, "away_spread_price": None,
            "home_spread": None, "home_spread_price": None,
            "total": None, "over_price": None, "under_price": None,
            "away_ml": None, "home_ml": None,
        }

    for sel in selections:
        meta = market_meta.get(sel.market_id)
        if meta is None:
            continue   # market we don't care about (props, futures, etc.)
        period, market_type = meta

        # Bucket key: for main rows, lines are stored individually
        # but all share the same row (line=None as bucket key).
        # For alt rows, key by the actual line so multiple lines emit multiple rows.
        bucket_line: float | None = None if market_type == "main" else sel.line
        key = (period, market_type, bucket_line)
        if key not in buckets:
            buckets[key] = _row_skeleton(period, market_type, bucket_line)
        row = buckets[key]

        name_lower = sel.name.lower()
        if "over" in name_lower:
            row["total"] = sel.line
            row["over_price"] = sel.american_odds
        elif "under" in name_lower:
            row["total"] = sel.line
            row["under_price"] = sel.american_odds
        elif sel.line is not None and (sel.name.startswith(event.home_team)
                                        or event.home_team in sel.name):
            row["home_spread"] = sel.line
            row["home_spread_price"] = sel.american_odds
        elif sel.line is not None and (sel.name.startswith(event.away_team)
                                        or event.away_team in sel.name):
            row["away_spread"] = sel.line
            row["away_spread_price"] = sel.american_odds
        elif sel.line is None and event.home_team in sel.name:
            row["home_ml"] = sel.american_odds
        elif sel.line is None and event.away_team in sel.name:
            row["away_ml"] = sel.american_odds

    return list(buckets.values())
```

- [ ] **Step 5: Run tests**

```bash
python -m pytest tests/test_dk_singles_parser.py -v
```

Expected: 3 passed.

- [ ] **Step 6: Commit**

```bash
cd ../..
git add mlb_sgp/scraper_draftkings_singles.py mlb_sgp/tests/test_dk_singles_parser.py
git commit -m "feat(mlb_sgp): DK singles parser converts selections to wide rows

parse_selections_to_wide_rows groups by (period, market_type, line) and
emits one row per bucket matching wagerzon mlb_odds schema. Unit-tested
for main spread+total+ML coalescing, alt-line splitting, and F5/FG isolation.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

### Task 9: DK singles scraper main flow + DB write

**Files:**
- Modify: `mlb_sgp/scraper_draftkings_singles.py`
- Create: `dk_odds/README.md`

- [ ] **Step 1: Add classify_market helper + scrape_singles + main**

Append to `mlb_sgp/scraper_draftkings_singles.py`:

```python
def classify_market(name: str) -> tuple[str, str] | None:
    """Map a DK market name to (period, market_type), or None to skip.

    Adjust this function based on Phase 1 findings — exact market names depend
    on what DK posts. Examples below assume English-style market names.
    """
    n = name.lower()

    # F-period detection
    if "1st 3 innings" in n or "f3" in n:
        period = "F3"
    elif "1st 5 innings" in n or "f5" in n:
        period = "F5"
    elif "1st 7 innings" in n or "f7" in n:
        period = "F7"
    else:
        period = "FG"

    # Market-type detection
    if "alt" in n and "run line" in n:
        return (period, "alternate_spreads")
    if "alt" in n and "total" in n:
        return (period, "alternate_totals")
    if "run line" in n or "moneyline" in n or "total" in n:
        return (period, "main")
    return None


def scrape_singles(verbose: bool = False) -> None:
    client = DraftKingsClient(verbose=verbose)
    events = client.list_events()
    print(f"[dk_singles] {len(events)} events to scrape", flush=True)

    fetch_time = datetime.utcnow()
    all_rows: list[dict] = []

    for event in events:
        try:
            markets = client.fetch_event_markets(event.event_id)
            selections = client.fetch_event_selections(event.event_id)

            # Build market_id -> (period, market_type) lookup
            market_meta: dict[str, tuple[str, str]] = {}
            for m in markets:
                classified = classify_market(m.name)
                if classified is not None:
                    market_meta[m.market_id] = classified

            rows = parse_selections_to_wide_rows(event, selections,
                                                  market_meta, fetch_time)
            all_rows.extend(rows)
            if verbose:
                print(f"  [{event.event_id}] {len(rows)} rows", flush=True)
        except Exception as e:
            print(f"  [{event.event_id}] FAILED: {e}", flush=True)
            continue   # per-game isolation

    write_to_duckdb(all_rows, fetch_time)
    print(f"[dk_singles] wrote {len(all_rows)} rows", flush=True)


def write_to_duckdb(rows: list[dict], fetch_time: datetime) -> None:
    """Atomic write: DELETE all → INSERT new in one transaction.

    DB path is dk_odds/dk.duckdb (relative to repo root).
    """
    db_path = Path(__file__).resolve().parent.parent / "dk_odds" / "dk.duckdb"
    db_path.parent.mkdir(exist_ok=True)

    con = duckdb.connect(str(db_path))
    try:
        con.execute("""
            CREATE TABLE IF NOT EXISTS mlb_odds (
                fetch_time        TIMESTAMP,
                sport_key         VARCHAR,
                game_id           VARCHAR,
                game_date         VARCHAR,
                game_time         VARCHAR,
                away_team         VARCHAR,
                home_team         VARCHAR,
                market            VARCHAR,
                period            VARCHAR,
                away_spread       FLOAT,
                away_spread_price INTEGER,
                home_spread       FLOAT,
                home_spread_price INTEGER,
                total             FLOAT,
                over_price        INTEGER,
                under_price       INTEGER,
                away_ml           INTEGER,
                home_ml           INTEGER
            )
        """)
        con.execute("BEGIN TRANSACTION")
        con.execute("DELETE FROM mlb_odds")
        if rows:
            con.executemany("""
                INSERT INTO mlb_odds VALUES
                ($fetch_time, $sport_key, $game_id, $game_date, $game_time,
                 $away_team, $home_team, $market, $period,
                 $away_spread, $away_spread_price, $home_spread, $home_spread_price,
                 $total, $over_price, $under_price, $away_ml, $home_ml)
            """, rows)
        con.execute("COMMIT")
    except Exception:
        con.execute("ROLLBACK")
        raise
    finally:
        con.close()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--verbose", action="store_true")
    args = p.parse_args()
    scrape_singles(verbose=args.verbose)


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Create dk_odds/ + README**

```bash
mkdir -p dk_odds
cat > dk_odds/README.md <<'EOF'
# dk_odds/

DraftKings single-leg odds DuckDB. Written by `mlb_sgp/scraper_draftkings_singles.py`.

## Tables

- `mlb_odds` — single-leg odds, offshore 18-column schema. One row per
  `(game_id, period, market_type, line)` group. Atomic rewrite each cycle.

## Read pattern

MLB.R reads this via `dbGetQuery(con, "SELECT * FROM mlb_odds")` and feeds
to `scraper_to_canonical()` for book_odds_by_book.
EOF
```

- [ ] **Step 3: Run end-to-end live scrape**

```bash
cd mlb_sgp
source venv/bin/activate
python scraper_draftkings_singles.py --verbose
```

Expected: prints event count, ~30 "X rows" lines, ends with "wrote NNN rows" where NNN is in the low thousands.

Verify DB content:

```bash
duckdb ../dk_odds/dk.duckdb "SELECT COUNT(*) AS n, COUNT(DISTINCT game_id) AS games, COUNT(DISTINCT period) AS periods, COUNT(DISTINCT market) AS market_tiers FROM mlb_odds"
```

Expected output (approximate): `n ≈ 800-1500, games ≈ 15-30, periods ≈ 2-4, market_tiers ≈ 1-3`.

- [ ] **Step 4: Commit**

```bash
cd ../..
git add mlb_sgp/scraper_draftkings_singles.py dk_odds/README.md
git commit -m "feat(mlb_sgp): scraper_draftkings_singles end-to-end

Walks all events, parses selections via classify_market + parser, writes
to dk_odds/dk.duckdb with atomic DELETE+INSERT. Per-game isolation: one
event's API error doesn't tank the scrape.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

### Task 10: FD singles scraper (mirror of Task 8 + 9)

**Files:**
- Create: `mlb_sgp/scraper_fanduel_singles.py`
- Create: `mlb_sgp/tests/test_fd_singles_parser.py`
- Create: `fd_odds/README.md`

- [ ] **Step 1: Write tests for FD parser**

`mlb_sgp/tests/test_fd_singles_parser.py`:

Copy the structure of `test_dk_singles_parser.py` but adapted to `Runner` dataclass and FD market naming conventions:

```python
from datetime import datetime
from mlb_sgp.fd_client import Event, Runner
from mlb_sgp.scraper_fanduel_singles import parse_runners_to_wide_rows


def test_main_game_lines_produce_one_row():
    event = Event("e1", "Yankees", "Red Sox", "2026-05-12T22:00:00Z")
    runners = [
        Runner("r1", "m_main_spread", "Yankees -1.5", -1.5, -120),
        Runner("r2", "m_main_spread", "Red Sox +1.5", 1.5, 110),
        Runner("r3", "m_main_total", "Over 9.5", 9.5, -105),
        Runner("r4", "m_main_total", "Under 9.5", 9.5, -115),
        Runner("r5", "m_ml", "Yankees", None, -150),
        Runner("r6", "m_ml", "Red Sox", None, 130),
    ]
    market_meta = {
        "m_main_spread": ("FG", "main"),
        "m_main_total":  ("FG", "main"),
        "m_ml":          ("FG", "main"),
    }
    rows = parse_runners_to_wide_rows(event, runners, market_meta,
                                       fetch_time=datetime(2026, 5, 12, 14, 0))
    assert len(rows) == 1
    r = rows[0]
    assert r["home_spread"] == -1.5
    assert r["total"] == 9.5
    assert r["home_ml"] == -150


def test_alt_totals_emit_separate_rows():
    event = Event("e1", "Yankees", "Red Sox", "2026-05-12T22:00:00Z")
    runners = [
        Runner("r1", "m_alt_total", "Over 8.5", 8.5, 120),
        Runner("r2", "m_alt_total", "Under 8.5", 8.5, -150),
        Runner("r3", "m_alt_total", "Over 10.5", 10.5, 180),
        Runner("r4", "m_alt_total", "Under 10.5", 10.5, -220),
    ]
    market_meta = {"m_alt_total": ("FG", "alternate_totals")}
    rows = parse_runners_to_wide_rows(event, runners, market_meta,
                                       fetch_time=datetime(2026, 5, 12, 14, 0))
    assert len(rows) == 2
    totals = sorted([r["total"] for r in rows])
    assert totals == [8.5, 10.5]
    assert all(r["market"] == "alternate_totals" for r in rows)
```

- [ ] **Step 2: Write scraper_fanduel_singles.py**

`mlb_sgp/scraper_fanduel_singles.py` — same structure as `scraper_draftkings_singles.py` but using `FanDuelClient` / `Runner`. The parser function `parse_runners_to_wide_rows` has the SAME logic as `parse_selections_to_wide_rows` — copy it verbatim, only the type signatures change. Write `to dk_odds/dk.duckdb` becomes `fd_odds/fd.duckdb`.

If you find the parser logic is genuinely identical (modulo dataclass names), extract a shared helper in a follow-up task — not now.

The `classify_market` function may need different FD-specific market name keywords. Use Phase 1 findings to know what FD posts. Examples:
- `"Spread Betting"` (FD's name for run line)
- `"Alternate Run Line"`
- `"Total Match Runs"` (FD's name for FG total)
- `"1st 5 Innings"` for F5

- [ ] **Step 3: Tests pass; live scrape works; commit**

```bash
cd mlb_sgp
python -m pytest tests/test_fd_singles_parser.py -v   # expect pass
python scraper_fanduel_singles.py --verbose            # live scrape
duckdb ../fd_odds/fd.duckdb "SELECT COUNT(*) AS n, COUNT(DISTINCT period) AS p FROM mlb_odds"
# expect non-zero rows

mkdir -p ../fd_odds
cat > ../fd_odds/README.md <<'EOF'
# fd_odds/

FanDuel single-leg odds DuckDB. Written by `mlb_sgp/scraper_fanduel_singles.py`.
Schema and read pattern: same as `dk_odds/README.md`.
EOF

cd ..
git add mlb_sgp/scraper_fanduel_singles.py mlb_sgp/tests/test_fd_singles_parser.py fd_odds/README.md
git commit -m "feat(mlb_sgp): scraper_fanduel_singles end-to-end

Mirrors DK singles scraper structure using FanDuelClient. Writes to
fd_odds/fd.duckdb. Parser unit-tested; live scrape verified.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Phase 6 — Orchestration

### Task 11: Add SCRAPER_CONFIGS entries to run.py

**Files:**
- Modify: `Answer Keys/run.py` (lines 18-47, SCRAPER_CONFIGS dict)

- [ ] **Step 1: Add the two new entries**

Edit `Answer Keys/run.py`, inside `SCRAPER_CONFIGS` (after the existing entries, before the closing `}`):

```python
    "draftkings_singles": {
        "script": "../mlb_sgp/scraper_draftkings_singles.py",
        "sports": ["mlb"]
    },
    "fanduel_singles": {
        "script": "../mlb_sgp/scraper_fanduel_singles.py",
        "sports": ["mlb"]
    },
```

- [ ] **Step 2: Verify orchestrator picks them up**

```bash
cd "Answer Keys"
python -c "
from run import SCRAPER_CONFIGS
mlb_scrapers = [n for n, c in SCRAPER_CONFIGS.items() if 'mlb' in c['sports']]
print('MLB scrapers configured:', mlb_scrapers)
assert 'draftkings_singles' in mlb_scrapers
assert 'fanduel_singles' in mlb_scrapers
print('OK')
"
```

Expected: prints scraper list including the two new ones.

- [ ] **Step 3: Dry-run check that the script path resolves**

```bash
python -c "
from pathlib import Path
from run import SCRAPER_CONFIGS
for name in ['draftkings_singles', 'fanduel_singles']:
    cfg = SCRAPER_CONFIGS[name]
    script_path = Path('.').resolve() / cfg['script']
    assert script_path.exists(), f'{script_path} not found'
    print(f'OK {name}: {script_path}')
"
```

Expected: prints two OK lines. If `script_path` doesn't exist, the relative path is wrong — fix it.

- [ ] **Step 4: Commit**

```bash
cd ..
git add "Answer Keys/run.py"
git commit -m "feat(run.py): add draftkings_singles + fanduel_singles to SCRAPER_CONFIGS

Both run in the parallel scrape phase (non-sharp) gated by .scrapers_done_mlb.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Phase 7 — MLB.R integration

### Task 12: Wire dk_odds + fd_odds into MLB.R

**Files:**
- Modify: `Answer Keys/MLB Answer Key/MLB.R` (around lines 850-900)

- [ ] **Step 1: Locate the current offshore-scraper load block**

Read `Answer Keys/MLB Answer Key/MLB.R` near where `wagerzon_odds`, `hoop88_odds`, etc. are loaded. Use grep to find the line numbers:

```bash
cd "Answer Keys/MLB Answer Key"
grep -n "wagerzon_odds <-\|hoop88_odds <-\|bfa_odds <-\|bookmaker_odds <-\|bet105_odds <-" MLB.R
```

These all follow the pattern:

```r
con_wz <- dbConnect(duckdb(), "../../wagerzon_odds/wagerzon.duckdb", read_only = TRUE)
wagerzon_odds <- dbGetQuery(con_wz, "SELECT * FROM mlb_odds WHERE ...")
dbDisconnect(con_wz)
```

- [ ] **Step 2: Add dk_odds + fd_odds loaders alongside the existing ones**

Insert after the last existing offshore loader (likely `bet105_odds`):

```r
# DraftKings single-leg odds (replaces Odds API for DK pills).
# Schema mirrors offshore mlb_odds; scraper_to_canonical handles it.
con_dk <- dbConnect(duckdb(), "../../dk_odds/dk.duckdb", read_only = TRUE)
dk_odds <- tryCatch(
  dbGetQuery(con_dk, "SELECT * FROM mlb_odds"),
  error = function(e) {
    warning(sprintf("[mlb] dk_odds load failed (%s) — DK pills will be empty",
                    conditionMessage(e)))
    NULL
  }
)
dbDisconnect(con_dk, shutdown = TRUE)

# FanDuel singles
con_fd <- dbConnect(duckdb(), "../../fd_odds/fd.duckdb", read_only = TRUE)
fd_odds <- tryCatch(
  dbGetQuery(con_fd, "SELECT * FROM mlb_odds"),
  error = function(e) {
    warning(sprintf("[mlb] fd_odds load failed (%s) — FD pills will be empty",
                    conditionMessage(e)))
    NULL
  }
)
dbDisconnect(con_fd, shutdown = TRUE)
```

- [ ] **Step 3: Swap book_odds_by_book DK/FD entries**

Find the `book_odds_by_book <- list(...)` block (around line 882) and change the DK/FD lines:

```r
# BEFORE
draftkings = odds_api_to_canonical(prefetched_long %>% filter(bookmaker_key == "draftkings")),
fanduel    = odds_api_to_canonical(prefetched_long %>% filter(bookmaker_key == "fanduel")),

# AFTER
draftkings = scraper_to_canonical(dk_odds, .game_id_lookup),
fanduel    = scraper_to_canonical(fd_odds, .game_id_lookup),
```

The `pinnacle = ...` line stays unchanged.

- [ ] **Step 4: Run pipeline end-to-end**

```bash
cd "../../"
cd "Answer Keys"
python run.py mlb
```

Expected: orchestrator runs all scrapers including `draftkings_singles` + `fanduel_singles`, R script completes, writes `mlb_bets_combined` + `mlb_bets_book_prices` to `mlb_mm.duckdb`.

Verify pill data:

```bash
duckdb "Answer Keys/mlb_mm.duckdb" "SELECT bookmaker, COUNT(*) FROM mlb_bets_book_prices GROUP BY 1 ORDER BY 1"
```

Expected: rows for `wagerzon`, `hoop88`, `bfa`, `bookmaker`, `bet105`, `draftkings`, `fanduel`, `pinnacle`. DK + FD row counts should be non-zero (and ideally similar to before, or larger if F3/F7 coverage improved).

- [ ] **Step 5: Commit**

```bash
git add "Answer Keys/MLB Answer Key/MLB.R"
git commit -m "feat(mlb): book_odds_by_book DK/FD read from per-book scrapers

dk_odds and fd_odds loaded from dk_odds/dk.duckdb and fd_odds/fd.duckdb,
fed through scraper_to_canonical(). Pinnacle still on Odds API.
Tryctach guards keep the pipeline running even if a scraper produced no
data (pills degrade to empty for that book).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Phase 8 — Smoke + docs + merge

### Task 13: Browser smoke test

**Files:** none (verification only)

- [ ] **Step 1: Start the dashboard**

```bash
cd "Answer Keys/MLB Dashboard"
./run.sh
```

Wait for "Listening on port 8083" or equivalent.

- [ ] **Step 2: Visual check**

Open `http://localhost:8083` in a browser. Verify on the **Bets** tab:
- DK pills now show data for F-period markets where DK actually posts them (improvement over Odds API).
- FD pills show data.
- Pinnacle pills behave the same as before (still from Odds API).
- Offshore pills (WZ, H88, BFA, BKM, B105) unchanged.

- [ ] **Step 3: Stop dashboard**

```bash
# Ctrl+C in the run.sh window
```

If something looks broken (DK/FD pills all empty), check:
- `duckdb dk_odds/dk.duckdb "SELECT COUNT(*) FROM mlb_odds"` — is the scrape DB populated?
- `Tools.R::resolve_offshore_teams()` — are DK/FD team names canonical?
- R warning output from the pipeline run.

- [ ] **Step 4: No commit (verification only)**

### Task 14: Documentation updates

**Files:**
- Modify: `mlb_sgp/README.md`
- Modify: `Answer Keys/CLAUDE.md`
- Modify: `CLAUDE.md` (project root)

- [ ] **Step 1: mlb_sgp/README.md — add Singles scrapers section**

Append before the existing "DraftKings API Endpoints" section:

```markdown
## Singles scrapers

Two scrapers fetch single-leg odds for the MLB Dashboard bets tab:

- `scraper_draftkings_singles.py` → writes `../dk_odds/dk.duckdb::mlb_odds`
- `scraper_fanduel_singles.py`    → writes `../fd_odds/fd.duckdb::mlb_odds`

Both use the same client classes (`dk_client.py`, `fd_client.py`) as the SGP
scrapers — no second API auth path, no duplicate rate-limit budget.

Output schema matches the wagerzon offshore convention; MLB.R consumes via
`scraper_to_canonical()`.

Run timing: orchestrated by `Answer Keys/run.py mlb` in the parallel scrape
phase (pre-MLB.R), gated by `.scrapers_done_mlb`. The SGP scrapers continue
to run post-MLB.R (they depend on `mlb_parlay_lines`) on a separate trigger.
```

- [ ] **Step 2: Answer Keys/CLAUDE.md — update Pipeline Flow**

Find the "Pipeline Flow (MLB)" block and update the parallel-scraper line:

```diff
- ├── [parallel] Other scrapers (wagerzon, hoop88, bfa) + MLB.R
+ ├── [parallel] Other scrapers (wagerzon, hoop88, bfa, draftkings_singles, fanduel_singles) + MLB.R
```

Also add a one-line note where the offshore convention is described:

```markdown
- DraftKings and FanDuel single-leg odds: scraper-sourced (per-book duckdb),
  not from the Odds API. Pinnacle remains on Odds API (no public API).
```

- [ ] **Step 3: Project-root CLAUDE.md — Odds screen section**

In the "MLB Dashboard — Odds screen + WZ single-bet placer" section, add to "Data flow":

```markdown
3. DraftKings and FanDuel pill data: written by
   `mlb_sgp/scraper_draftkings_singles.py` and
   `mlb_sgp/scraper_fanduel_singles.py` to per-book DuckDBs
   (`dk_odds/dk.duckdb`, `fd_odds/fd.duckdb`). MLB.R reads via
   `scraper_to_canonical()`. Pinnacle still comes from Odds API
   (`prefetched_long` filtered to `bookmaker_key == "pinnacle"`).
```

- [ ] **Step 4: Commit docs**

```bash
git add mlb_sgp/README.md "Answer Keys/CLAUDE.md" CLAUDE.md
git commit -m "docs: DK/FD singles scrapers — README + CLAUDE.md updates

Adds Singles scrapers section to mlb_sgp/README, updates Pipeline Flow
diagram in Answer Keys/CLAUDE.md, and notes DK/FD source change in
project-root CLAUDE.md MLB Dashboard section.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

### Task 15: Pre-merge executive review

**Files:** none (review only — diff inspection)

- [ ] **Step 1: View the full diff vs main**

```bash
git diff main..HEAD --stat
git diff main..HEAD | less
```

- [ ] **Step 2: Run the executive review checklist (from CLAUDE.md)**

For each item, write a one-line YES/NO/N/A:

```
Data integrity:
  - No duplicate writes? (atomic DELETE+INSERT in singles scrapers) →
  - Proper deduplication? (mlb_odds is rewritten each cycle) →
  - Incomplete records filtered out? (per-game try/except keeps partial writes consistent) →

Resource safety:
  - All DB connections use on.exit / try/finally? →
  - No lock-file leaks on crash? →

Edge cases:
  - Off-season behavior? (empty events → 0 rows written → MLB.R pills empty, same as today) →
  - Empty mlb_bets_book_prices? (still loads, dashboard handles via tryCatch) →
  - First-run with no existing data? (DELETE FROM mlb_odds is no-op if empty) →
  - Timezone boundaries? (fetch_time stored as UTC; consumers convert) →

Dead code:
  - No unused flags / functions / imports introduced? →

Log/disk hygiene:
  - Log rotation in place? (uses existing mlb_sgp/logs/ pattern) →
  - No unbounded file growth? →

Security:
  - No secrets in logs? →
  - No API keys exposed? →
```

- [ ] **Step 3: Document findings**

Append to `docs/superpowers/specs/2026-05-12-mlb-dk-fd-singles-scrapers-design.md`:

```markdown
---

## Pre-merge executive review (YYYY-MM-DD)

ISSUES TO FIX:
- (list any items)

ACCEPTABLE RISKS:
- (list any items)
```

Commit:

```bash
git add docs/superpowers/specs/2026-05-12-mlb-dk-fd-singles-scrapers-design.md
git commit -m "docs: pre-merge executive review for DK/FD singles

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

- [ ] **Step 4: Fix any ISSUES TO FIX**

Each fix is its own commit. If the list is empty, skip.

### Task 16: Merge to main + cleanup

**Files:** none (git operations only)

- [ ] **Step 1: Get explicit user approval to merge**

Per CLAUDE.md feedback_always_ask_merge: do not merge to main without confirmation. Ask the user:

> "All tasks complete, executive review clean. Ready to merge `worktree-feature+mlb-dk-fd-singles-scrapers` → `main`. Approve?"

WAIT for explicit "yes" before proceeding.

- [ ] **Step 2: Merge**

```bash
cd /Users/callancapitolo/NFLWork   # original working copy (not the worktree)
git checkout main
git pull origin main
git merge worktree-feature+mlb-dk-fd-singles-scrapers --no-ff -m "Merge: MLB DK/FD single-leg scrapers

Replaces Odds API as the source for DK and FD pill data on the MLB
Dashboard bets tab. Pinnacle stays on Odds API. Adds dk_client.py /
fd_client.py shared by the existing SGP scrapers (no behavior change there).
"
```

- [ ] **Step 3: Push**

After explicit user approval:

```bash
git push origin main
```

- [ ] **Step 4: Cleanup worktree**

```bash
git worktree remove /Users/callancapitolo/NFLWork/.claude/worktrees/feature+mlb-dk-fd-singles-scrapers
git branch -d worktree-feature+mlb-dk-fd-singles-scrapers
```

- [ ] **Step 5: Verify clean state**

```bash
git worktree list   # should NOT include the feature worktree
git branch          # should NOT include the feature branch
git log --oneline -5
```

---

## Self-review checklist (this plan)

- **Spec coverage:** Every section of the spec maps to a task — motivation (covered in goal), architecture (Tasks 2-5), data flow (Tasks 8-12), orchestration (Task 11), failure modes (Task 9 atomic write, Task 12 tryCatch), team-name resolution (Task 1 Step 4), uncertainties (Task 1), testing (Tasks 2-10), docs (Task 14), out of scope (NOT implemented — Pinnacle scraper, SGP restructure).
- **Placeholder scan:** all "Adjust based on Phase 1 findings" notes are deliberate hooks, not placeholders. Field names in client parsers (Tasks 4, 5) are written as best-guess and explicitly flagged for adjustment.
- **Type consistency:** `Event`, `Selection`, `Market`, `Runner` dataclasses are consistent across tasks. `parse_selections_to_wide_rows` and `parse_runners_to_wide_rows` have parallel structure with matching kwargs.

If any subagent finds a real placeholder during execution, fix inline and flag in the task's commit message.
