# FanDuel Full Market Coverage Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the FanDuel singles scraper cover every market type the bets tab prices — by fetching both FanDuel event-page tabs and replacing its exact-name whitelist with a DraftKings-style keyword classifier — so FanDuel pills match DraftKings wherever FanDuel posts the market (fixing the missing F7/F3 totals + run lines).

**Architecture:** Two Python-only layers. (1) `fd_client.py` gains `fetch_event_page(event_id, tab)` that returns markets **and** runners from one payload; the scraper calls it once per tab (`""` default + `"same-game-parlay-"`) and merges with dedup. (2) `scraper_fanduel_singles.py` replaces `_FD_MARKET_WHITELIST`/`classify_market` with a keyword classifier ported from DK, tuned with FanDuel-specific junk exclusions and event-team-based team-total filtering. No R / dashboard / model changes — `get_fd_odds` already passes `period` through generically.

**Tech Stack:** Python 3.14, curl_cffi (FanDuel TLS session), DuckDB, pytest 9.

**Spec:** `docs/superpowers/specs/2026-05-23-fanduel-full-market-coverage-design.md`

**Test invocation (worktree has no venv — use main repo venv on the path):**
```
PYTHONPATH=mlb_sgp /Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m pytest mlb_sgp/tests/test_fanduel_singles.py -q
```
Baseline is green (2 passed) — that's the regression floor.

---

## File Structure

- **Modify** `mlb_sgp/fd_client.py` — add `fetch_event_page(event_id, tab)`; make `fetch_event_markets` / `fetch_event_runners` thin wrappers over it with a `tab` param (default `"same-game-parlay-"`, preserving every current caller).
- **Modify** `mlb_sgp/scraper_fanduel_singles.py` — replace whitelist + `classify_market` with keyword classifier (`classify_market(name, home_team=None, away_team=None)`); add `fetch_merged_markets_and_runners(client, event_id, tabs)` helper; rewire `scrape_singles` to use it and pass event teams into classification.
- **Modify** `mlb_sgp/tests/test_fanduel_singles.py` — add table-driven classifier tests + merge-helper test (existing 2 tests stay green).
- **Modify** `mlb_sgp/README.md` — document both-tab fetch + keyword classifier.
- **Modify** memory `fd_event_page_tab_coverage_gap.md` — mark resolved.

---

## Task 1: `fd_client.fetch_event_page` (one fetch → markets + runners, tab-parameterized)

**Files:**
- Modify: `mlb_sgp/fd_client.py:65-196` (`fetch_event_runners`, `fetch_event_markets`)
- Test: `mlb_sgp/tests/test_fd_client_event_page.py` (create)

- [ ] **Step 1: Write the failing test**

Create `mlb_sgp/tests/test_fd_client_event_page.py`:

```python
"""Tests for FanDuelClient.fetch_event_page (both-tab combined fetch)."""
from mlb_sgp.fd_client import FanDuelClient, Market, Runner


class _FakeResp:
    def __init__(self, payload):
        self._payload = payload
        self.status_code = 200

    def json(self):
        return self._payload


class _FakeSession:
    """Records the last URL requested; returns a fixed payload.

    Mirrors the unit-test FakeSession contract the client already supports:
    its get() takes no headers/timeout kwargs, so the client's
    `except TypeError` fallback path is exercised.
    """
    def __init__(self, payload):
        self._payload = payload
        self.last_url = None

    def get(self, url):
        self.last_url = url
        return _FakeResp(self._payload)


def _client_with(payload):
    # Bypass __init__ (which builds a real curl_cffi session) and inject a fake.
    c = FanDuelClient.__new__(FanDuelClient)
    c.session = _FakeSession(payload)
    c.verbose = False
    return c


_PAYLOAD = {
    "attachments": {
        "markets": {
            "m1": {
                "marketId": "m1",
                "marketName": "First 7 Innings Total Runs",
                "runners": [
                    {"selectionId": "r1", "runnerName": "Over",
                     "runnerStatus": "ACTIVE", "handicap": 7.5,
                     "winRunnerOdds": {"americanDisplayOdds": {"americanOdds": -128}}},
                    {"selectionId": "r2", "runnerName": "Under",
                     "runnerStatus": "ACTIVE", "handicap": 7.5,
                     "winRunnerOdds": {"americanDisplayOdds": {"americanOdds": 104}}},
                ],
            }
        }
    }
}


def test_fetch_event_page_returns_markets_and_runners():
    c = _client_with(_PAYLOAD)
    markets, runners = c.fetch_event_page("99", tab="")
    assert markets == [Market(market_id="m1", name="First 7 Innings Total Runs")]
    assert {r.runner_id for r in runners} == {"r1", "r2"}
    over = next(r for r in runners if r.runner_id == "r1")
    assert over.american_odds == -128 and over.line == 7.5


def test_fetch_event_page_uses_requested_tab():
    c = _client_with(_PAYLOAD)
    c.fetch_event_page("99", tab="")
    assert "tab=&" in c.session.last_url  # empty default tab
    c.fetch_event_page("99", tab="same-game-parlay-")
    assert "tab=same-game-parlay-" in c.session.last_url


def test_wrappers_delegate_to_fetch_event_page():
    c = _client_with(_PAYLOAD)
    assert c.fetch_event_markets("99", tab="") == [
        Market(market_id="m1", name="First 7 Innings Total Runs")
    ]
    assert {r.runner_id for r in c.fetch_event_runners("99", tab="")} == {"r1", "r2"}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=mlb_sgp /Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m pytest mlb_sgp/tests/test_fd_client_event_page.py -q`
Expected: FAIL — `AttributeError: 'FanDuelClient' object has no attribute 'fetch_event_page'`.

- [ ] **Step 3: Implement `fetch_event_page` and convert the wrappers**

In `mlb_sgp/fd_client.py`, replace the **entire** `fetch_event_runners` method (lines 65-140) and the **entire** `fetch_event_markets` method (lines 142-196) with the following three methods:

```python
    def fetch_event_page(
        self, event_id: str, tab: str = "same-game-parlay-"
    ) -> tuple[list["Market"], list["Runner"]]:
        """Fetch one FD event-page tab; return (markets, runners) from it.

        FD's event-page payload carries both market metadata and runners in a
        single response, so one HTTP GET yields both. The `tab` param selects
        which server-side market slice FD returns — there is no single tab with
        every market, so the singles scraper calls this once per tab and merges
        (see fetch_merged_markets_and_runners in scraper_fanduel_singles).

        Market dedup is by marketId (the recursive walk may visit a market
        twice). Runners are read off each market dict; runnerStatus values
        other than "ACTIVE" (when present) and runners with no american odds
        or no selectionId are skipped.
        """
        from scraper_fanduel_sgp import FD_EVENT_PAGE_URL, FD_AK, FD_HEADERS

        url = (
            f"{FD_EVENT_PAGE_URL}?_ak={FD_AK}&eventId={event_id}"
            f"&tab={tab}"
            f"&useCombinedTouchdownsVirtualMarket=true&useQuickBets=true"
        )
        try:
            resp = self.session.get(url, headers=FD_HEADERS, timeout=20)
        except TypeError:
            # FakeSession in unit tests doesn't accept headers/timeout kwargs
            resp = self.session.get(url)

        if getattr(resp, "status_code", 200) != 200:
            return [], []
        payload = resp.json()

        market_dicts: list[dict] = []

        def walk(o):
            if isinstance(o, dict):
                if "marketId" in o and "marketName" in o:
                    market_dicts.append(o)
                for v in o.values():
                    walk(v)
            elif isinstance(o, list):
                for it in o:
                    walk(it)

        walk(payload)

        # Dedup by marketId; first-seen wins (duplicates are identical).
        seen: dict[str, dict] = {}
        for m in market_dicts:
            mid = str(m.get("marketId", ""))
            if mid and mid not in seen:
                seen[mid] = m

        markets = [
            Market(market_id=mid, name=str(m.get("marketName", "") or ""))
            for mid, m in seen.items()
        ]

        runners: list[Runner] = []
        for mid, m in seen.items():
            for run in m.get("runners", []) or []:
                if not isinstance(run, dict):
                    continue
                status = str(run.get("runnerStatus", "")).upper()
                if status and status != "ACTIVE":
                    continue
                american = _extract_american_odds(run)
                if american is None:
                    continue
                sid = run.get("selectionId")
                if sid is None:
                    continue
                runners.append(Runner(
                    runner_id=str(sid),
                    market_id=mid,
                    name=str(run.get("runnerName", "") or ""),
                    line=_extract_line(run),
                    american_odds=american,
                ))
        return markets, runners

    def fetch_event_runners(
        self, event_id: str, tab: str = "same-game-parlay-"
    ) -> list["Runner"]:
        """Back-compat wrapper: runners from one tab."""
        _, runners = self.fetch_event_page(event_id, tab)
        return runners

    def fetch_event_markets(
        self, event_id: str, tab: str = "same-game-parlay-"
    ) -> list["Market"]:
        """Back-compat wrapper: market metadata from one tab."""
        markets, _ = self.fetch_event_page(event_id, tab)
        return markets
```

- [ ] **Step 4: Run test to verify it passes**

Run: `PYTHONPATH=mlb_sgp /Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m pytest mlb_sgp/tests/test_fd_client_event_page.py -q`
Expected: PASS (3 passed).

- [ ] **Step 5: Run the existing FD + DK suites to confirm no regression**

Run: `PYTHONPATH=mlb_sgp /Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m pytest mlb_sgp/tests/test_fanduel_singles.py mlb_sgp/tests/test_fanduel_orchestrator.py -q`
Expected: PASS (all green — wrappers preserve the old single-tab behavior).

- [ ] **Step 6: Commit**

```bash
git add mlb_sgp/fd_client.py mlb_sgp/tests/test_fd_client_event_page.py
git commit -m "feat(fd_client): fetch_event_page returns markets+runners, tab-parameterized

One GET per tab yields both markets and runners; existing fetch_event_*
become thin wrappers (tab defaults to same-game-parlay-, no caller change).
Enables both-tab coverage in the singles scraper.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: Keyword `classify_market` (replace the whitelist)

**Files:**
- Modify: `mlb_sgp/scraper_fanduel_singles.py:40-115` (regex helpers area + `_FD_MARKET_WHITELIST` + `classify_market`)
- Test: `mlb_sgp/tests/test_fanduel_singles.py` (add a table-driven test)

- [ ] **Step 1: Write the failing test**

Add to `mlb_sgp/tests/test_fanduel_singles.py` (keep the existing two tests):

```python
import pytest

# (home, away) used for team-total exclusion cases
_H, _A = "San Francisco Giants", "Chicago White Sox"

@pytest.mark.parametrize("name,home,away,expected", [
    # --- ACCEPTED: real game lines, every period x type FD posts ---
    ("Run Line", None, None, ("FG", "main")),
    ("Total Runs", None, None, ("FG", "main")),
    ("Moneyline", None, None, ("FG", "main")),
    ("Alternate Run Lines", None, None, ("FG", "alternate_spreads")),
    ("Alternate Total Runs", None, None, ("FG", "alternate_totals")),
    ("First 5 Innings Run Line", None, None, ("F5", "main")),
    ("First 5 Innings Total Runs", None, None, ("F5", "main")),
    ("First 5 Innings Money Line", None, None, ("F5", "main")),  # NOTE the space
    ("First 5 Innings Alternate Run Lines", None, None, ("F5", "alternate_spreads")),
    ("First 5 Innings Alternate Total Runs", None, None, ("F5", "alternate_totals")),
    ("First 7 Innings Total Runs", None, None, ("F7", "main")),
    ("First 7 Innings Run Line", None, None, ("F7", "main")),
    ("First 3 Innings Total Runs", None, None, ("F3", "main")),
    ("First 3 Innings Run Line", None, None, ("F3", "main")),
    # --- REJECTED: one example per junk family ---
    ("First 5 Innings Run Line / Total Runs Parlay", None, None, None),  # parlay
    ("Line / Total Parlay 7", None, None, None),                         # parlay
    ("Total Runs (Bands)", None, None, None),                            # bands
    ("Moneyline Away Listed", None, None, None),                         # listed
    ("First 7 Innings Result", None, None, None),                        # 3-way result
    ("First 6 Innings Result", None, None, None),                        # 3-way result
    ("7th Inning Total Runs", None, None, None),                         # single inning
    ("7th Inning Run Line", None, None, None),                           # single inning
    ("First 5 Innings Winning Margin (5-Way)", None, None, None),        # winning margin
    ("Race To 7 Runs", None, None, None),                                # race to
    ("Tri-Bet", None, None, None),                                       # no line keyword
    ("Chicago White Sox Total Runs", _H, _A, None),                      # team total (away)
    ("San Francisco Giants Alt. Total Runs", _H, _A, None),              # team total (home)
    ("Random Garbage Market", None, None, None),
])
def test_classify_market_keyword(name, home, away, expected):
    from mlb_sgp.scraper_fanduel_singles import classify_market
    assert classify_market(name, home, away) == expected
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=mlb_sgp /Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m pytest "mlb_sgp/tests/test_fanduel_singles.py::test_classify_market_keyword" -q`
Expected: FAIL — current whitelist returns `None` for `First 7 Innings Total Runs` (and `classify_market` takes only one arg, so the 3-arg calls error).

- [ ] **Step 3: Replace the whitelist with the keyword classifier**

In `mlb_sgp/scraper_fanduel_singles.py`:

(a) Ensure `re` is imported (it already is, used by the alt regexes).

(b) **Delete** the entire `_FD_MARKET_WHITELIST` dict (lines ~64-105, including its big leading comment block) and the existing `classify_market` (lines ~108-115).

(c) Insert in their place:

```python
# Single-inning markets ("7th Inning Total Runs", "1st Inning Run Line") are
# out of scope — only cumulative "First N Innings" periods map to bet cards.
# Plural "Innings" (First 5 Innings) deliberately does NOT match this.
_SINGLE_INNING_RE = re.compile(r"\b\d+(st|nd|rd|th)\s+inning\b", re.IGNORECASE)

# Substrings that disqualify a market. Ported from DraftKings'
# classify_market exclusion list, plus FD-specific junk that keyword-collides
# with real game lines (FD posts ~150 markets/event):
#   - "parlay"  -> "First 5 Innings Run Line / Total Runs Parlay", "Line / Total Parlay N"
#   - "listed"  -> "Moneyline Away/Home/Both Listed" (pitcher-conditional; keep plain "Moneyline")
#   - "bands"   -> "Total Runs (Bands)"
#   - "tri-bet", "specials" -> FD novelty markets
_FD_EXCLUDE_KEYWORDS = (
    "team total", "player", "prop", "futures", "to record", "to score",
    "to hit", "first to", "race to", "correct score", "winning margin",
    "total bases", "rbis", "hits o/u", "strikeouts thrown", "odd/even",
    "score last", "bat bottom", "highest scoring", "most innings",
    "last run", "both teams to score",
    "parlay", "listed", "bands", "tri-bet", "specials",
)


def classify_market(
    name: str,
    home_team: str | None = None,
    away_team: str | None = None,
) -> tuple[str, str] | None:
    """Map an FD market name to (period, market_type), or None to skip.

    Keyword classifier (mirrors scraper_draftkings_singles.classify_market) so
    FD picks up every FG/F5/F7/F3 main + alt line it posts and auto-covers new
    ones. period is one of "FG", "F3", "F5", "F7".

    home_team / away_team are optional; when provided, any market name that
    contains a team name is treated as a per-team market (team total) and
    excluded. They default to None so existing single-arg callers/tests work.
    Game lines ("Run Line", "Total Runs", "Moneyline", "Alternate ...") never
    contain a team name, so this never false-excludes them.
    """
    n = name.lower()

    # Team totals: "<Team> Total Runs", "<Team> Alt. Total Runs".
    for team in (home_team, away_team):
        if team and team.lower() in n:
            return None

    if any(k in n for k in _FD_EXCLUDE_KEYWORDS):
        return None
    if _SINGLE_INNING_RE.search(n):
        return None

    # Period detection. Default FG; F3/F5/F7 if matched explicitly.
    if "first 7 innings" in n or "1st 7 innings" in n:
        period = "F7"
    elif "first 5 innings" in n or "1st 5 innings" in n:
        period = "F5"
    elif "first 3 innings" in n or "1st 3 innings" in n:
        period = "F3"
    else:
        period = "FG"

    # Market-type detection. Check "alternate" before main. Match both
    # "moneyline" (FD's FG name) and "money line" (FD's F-period name, e.g.
    # "First 5 Innings Money Line") so F-period MLs are not silently dropped.
    if "alternate" in n and "run line" in n:
        return (period, "alternate_spreads")
    if "alternate" in n and "total" in n:
        return (period, "alternate_totals")
    if "run line" in n:
        return (period, "main")
    if "moneyline" in n or "money line" in n:
        return (period, "main")
    if "total" in n:
        return (period, "main")
    return None
```

- [ ] **Step 4: Run the classifier tests**

Run: `PYTHONPATH=mlb_sgp /Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m pytest mlb_sgp/tests/test_fanduel_singles.py -q`
Expected: PASS (the new parametrized test + the two existing tests all green).

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/scraper_fanduel_singles.py mlb_sgp/tests/test_fanduel_singles.py
git commit -m "feat(fd_singles): keyword classify_market replacing exact-name whitelist

Mirrors DK's classifier (period + market-type detection) with FD-specific
junk exclusions (parlay/listed/bands/tri-bet/specials) and event-team-based
team-total filtering. Matches both 'moneyline' and 'money line' so F5 ML
survives. Picks up F7/F3 lines and future FD markets automatically.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: Both-tab merge + scrape wiring

**Files:**
- Modify: `mlb_sgp/scraper_fanduel_singles.py:291-340` (`scrape_singles`)
- Test: `mlb_sgp/tests/test_fanduel_singles.py` (add merge-helper test)

- [ ] **Step 1: Write the failing test**

Add to `mlb_sgp/tests/test_fanduel_singles.py`:

```python
def test_fetch_merged_markets_and_runners_unions_and_dedups():
    """Two tabs return overlapping markets/runners; merge unions + dedups."""
    from mlb_sgp.fd_client import Market, Runner
    from mlb_sgp.scraper_fanduel_singles import fetch_merged_markets_and_runners

    class _FakeClient:
        def fetch_event_page(self, event_id, tab):
            if tab == "":  # default tab: F7 total + shared FG run line
                return (
                    [Market("mF7", "First 7 Innings Total Runs"),
                     Market("mFG", "Run Line")],
                    [Runner("rO", "mF7", "Over", 7.5, -128),
                     Runner("rRLh", "mFG", "Giants -1.5", -1.5, 120)],
                )
            # sgp tab: shared FG run line (dup) + F5 total
            return (
                [Market("mFG", "Run Line"),
                 Market("mF5", "First 5 Innings Total Runs")],
                [Runner("rRLh", "mFG", "Giants -1.5", -1.5, 120),
                 Runner("rF5o", "mF5", "Over", 4.5, -110)],
            )

    markets, runners = fetch_merged_markets_and_runners(
        _FakeClient(), "1", tabs=("", "same-game-parlay-"))
    assert {m.market_id for m in markets} == {"mF7", "mFG", "mF5"}  # union, deduped
    assert {r.runner_id for r in runners} == {"rO", "rRLh", "rF5o"}  # rRLh once
```

- [ ] **Step 2: Run test to verify it fails**

Run: `PYTHONPATH=mlb_sgp /Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m pytest "mlb_sgp/tests/test_fanduel_singles.py::test_fetch_merged_markets_and_runners_unions_and_dedups" -q`
Expected: FAIL — `ImportError: cannot import name 'fetch_merged_markets_and_runners'`.

- [ ] **Step 3: Add the merge helper and rewire `scrape_singles`**

In `mlb_sgp/scraper_fanduel_singles.py`, add the helper + a module constant **above** `scrape_singles`:

```python
# FD's event-page returns different market slices per tab; no single tab has
# everything. The default tab carries F7/F3 + per-inning lines; the SGP tab
# carries FG-alts + all F5. Coverage = union of both. Verified 2026-05-23.
FD_TABS = ("", "same-game-parlay-")


def fetch_merged_markets_and_runners(client, event_id, tabs=FD_TABS):
    """Fetch each tab once, union markets (dedup by market_id) and runners
    (dedup by runner_id). Returns (list[Market], list[Runner])."""
    markets_by_id = {}
    runners_by_id = {}
    for tab in tabs:
        markets, runners = client.fetch_event_page(event_id, tab)
        for m in markets:
            markets_by_id.setdefault(m.market_id, m)
        for r in runners:
            runners_by_id.setdefault(r.runner_id, r)
    return list(markets_by_id.values()), list(runners_by_id.values())
```

Then in `scrape_singles`, replace the per-event body (the `try:` block at lines ~310-328 that calls `client.fetch_event_markets` / `client.fetch_event_runners`) with:

```python
        try:
            markets, runners = fetch_merged_markets_and_runners(
                client, event.event_id)

            market_meta: dict[str, tuple[str, str]] = {}
            for m in markets:
                classified = classify_market(
                    m.name, event.home_team, event.away_team)
                if classified is not None:
                    market_meta[m.market_id] = classified

            rows = parse_runners_to_wide_rows(event, runners, market_meta, fetch_time)
            all_rows.extend(rows)
            if verbose:
                print(
                    f"  [{event.event_id}] {event.away_team} @ {event.home_team}: "
                    f"{len(rows)} rows ({len(market_meta)} in-scope markets, "
                    f"{len(runners)} runners total)",
                    flush=True,
                )
        except Exception as e:
            print(f"  [{event.event_id}] FAILED: {e}", flush=True)
            failed.append(event.event_id)
            continue
```

- [ ] **Step 4: Run the full FD test suite**

Run: `PYTHONPATH=mlb_sgp /Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m pytest mlb_sgp/tests/test_fanduel_singles.py mlb_sgp/tests/test_fd_client_event_page.py mlb_sgp/tests/test_fanduel_orchestrator.py -q`
Expected: PASS (all green).

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/scraper_fanduel_singles.py mlb_sgp/tests/test_fanduel_singles.py
git commit -m "feat(fd_singles): merge both event-page tabs in scrape loop

fetch_merged_markets_and_runners unions default + SGP tabs (dedup by id);
scrape_singles passes event teams into classify_market. FD now sees the
F7/F3 lines that live only on the default tab.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: Validation gate — live-market classifier audit (manual, not committed)

**Files:**
- Create (throwaway, in `$CLAUDE_JOB_DIR`, NOT committed): `fd_classify_audit.py`

- [ ] **Step 1: Write the audit script**

```python
"""Audit: run the new classify_market against the FULL live FD market list
(both tabs) for several real games. Print ACCEPTED vs REJECTED for eyeball
confirmation: every FG/F5/F7/F3 line accepted with right period/type; zero
junk accepted."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path("/Users/callancapitolo/NFLWork/.claude/worktrees/fanduel-full-market-coverage/mlb_sgp")))
from fd_client import FanDuelClient
from scraper_fanduel_singles import classify_market, FD_TABS

client = FanDuelClient(verbose=False)
events = client.list_events()
# pick the first 3 events that return a rich market set
checked = 0
for e in events:
    names = {}
    for tab in FD_TABS:
        for m in client.fetch_event_page(e.event_id, tab)[0]:
            names.setdefault(m.name, None)
    if len(names) < 50:
        continue  # skip pre-game stubs
    checked += 1
    print(f"\n===== {e.away_team} @ {e.home_team} ({len(names)} markets) =====")
    acc, rej = [], []
    for name in sorted(names):
        cls = classify_market(name, e.home_team, e.away_team)
        (acc if cls else rej).append((name, cls))
    print("--- ACCEPTED ---")
    for name, cls in acc:
        print(f"   {cls}  {name!r}")
    print("--- REJECTED (spot-check for any real line wrongly dropped) ---")
    for name, _ in rej:
        print(f"   {name!r}")
    if checked >= 3:
        break
```

- [ ] **Step 2: Run it and confirm the gate**

Run: `/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python "$CLAUDE_JOB_DIR/fd_classify_audit.py"`

Confirm by inspection (GO/NO-GO):
- ACCEPTED contains `("F7","main") 'First 7 Innings Total Runs'`, `("F7","main") 'First 7 Innings Run Line'`, the F3 pair, all FG main+alt, all F5 main+alt+ML.
- ACCEPTED contains **zero** junk (no parlay/bands/listed/team-total/per-inning/result/winning-margin/race-to).
- REJECTED contains **no** real FG/F5/F7/F3 game line.

If any miss, fix `_FD_EXCLUDE_KEYWORDS` / detection in Task 2, re-run Task 2 tests, re-audit.

- [ ] **Step 3: Remove the script**

Run: `rm -f "$CLAUDE_JOB_DIR/fd_classify_audit.py"`

---

## Task 5: End-to-end scrape + regression diff + ground truth + gates (manual)

**Files:** none committed (verification only). The worktree scraper writes to
`<worktree>/fd_odds/fd.duckdb` (isolated from the live `~/NFLWork/fd_odds/fd.duckdb`).

- [ ] **Step 1: Snapshot the BEFORE baseline from the live DB**

Run:
```bash
/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python - <<'PY'
import duckdb
con = duckdb.connect('/Users/callancapitolo/NFLWork/fd_odds/fd.duckdb', read_only=True)
print("BEFORE (live, main-branch scraper) — rows by period x market:")
for r in con.execute("""
  SELECT period, market, count(*) FROM mlb_odds
  WHERE fetch_time = (SELECT max(fetch_time) FROM mlb_odds)
  GROUP BY 1,2 ORDER BY 1,2""").fetchall():
    print("  ", r)
PY
```
Record the FG/F5 counts.

- [ ] **Step 2: Re-pull a fresh F7 ground-truth quote (odds drift)**

Run (capture current FD F7 total/run-line for one live game):
```bash
/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python - <<'PY'
import sys; sys.path.insert(0, '/Users/callancapitolo/NFLWork/.claude/worktrees/fanduel-full-market-coverage/mlb_sgp')
from fd_client import FanDuelClient
c = FanDuelClient(); evs = c.list_events()
e = next(x for x in evs if "Giants" in x.home_team or "Reds" in x.home_team)
mks, runs = c.fetch_event_page(e.event_id, "")  # default tab
for m in mks:
    if m.name in ("First 7 Innings Total Runs", "First 7 Innings Run Line"):
        print(m.name, [(r.name, r.line, r.american_odds) for r in runs if r.market_id == m.market_id])
PY
```
Record the numbers for Step 4 comparison.

- [ ] **Step 3: Run the worktree scraper for real**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.claude/worktrees/fanduel-full-market-coverage/mlb_sgp \
  && /Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python scraper_fanduel_singles.py mlb --verbose 2>&1 | tail -20
```
Expected: `[fd_singles] wrote N rows` with N noticeably higher than before (F7/F3 added).

- [ ] **Step 4: AFTER snapshot + regression diff + ground truth**

Run:
```bash
/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python - <<'PY'
import duckdb
con = duckdb.connect('/Users/callancapitolo/NFLWork/.claude/worktrees/fanduel-full-market-coverage/fd_odds/fd.duckdb', read_only=True)
print("AFTER (worktree, new scraper) — rows by period x market:")
for r in con.execute("""
  SELECT period, market, count(*) FROM mlb_odds GROUP BY 1,2 ORDER BY 1,2""").fetchall():
    print("  ", r)
print("\nF7/F3 rows present?")
print(con.execute("SELECT period, count(*) FROM mlb_odds WHERE period IN ('F7','F3') GROUP BY 1").fetchall())
print("\nDuplicate (period,market,total,home_spread) rows (should be empty):")
print(con.execute("""
  SELECT period, market, total, home_spread, count(*) c FROM mlb_odds
  GROUP BY 1,2,3,4 HAVING count(*) > 1""").fetchall())
print("\nF7 total sample (compare to Step 2 ground truth):")
print(con.execute("""
  SELECT home_team, total, over_price, under_price FROM mlb_odds
  WHERE period='F7' AND total IS NOT NULL LIMIT 5""").fetchall())
PY
```
Confirm: (a) F7 + F3 rows exist; (b) FG/F5 counts ≥ Step 1 baseline; (c) duplicate query returns empty; (d) F7 total matches Step 2's fresh quote.

- [ ] **Step 5: Timezone parity gate (required by repo CLAUDE.md)**

Run:
```bash
cd /Users/callancapitolo/NFLWork/.claude/worktrees/fanduel-full-market-coverage \
  && /Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python -m pytest tests/timezone_parity_test.py -q 2>&1 | tail -15
```
Expected: PASS (merging tabs must not perturb `game_start_time`). If the test
reads the live DB path, point it at (or copy in) the worktree `fd.duckdb`
first; do NOT symlink (per CLAUDE.md).

- [ ] **Step 6: FD-vs-DK parity check**

Run:
```bash
/Users/callancapitolo/NFLWork/mlb_sgp/venv/bin/python - <<'PY'
import duckdb
fd = duckdb.connect('/Users/callancapitolo/NFLWork/.claude/worktrees/fanduel-full-market-coverage/fd_odds/fd.duckdb', read_only=True)
dk = duckdb.connect('/Users/callancapitolo/NFLWork/dk_odds/dk.duckdb', read_only=True)
q = "SELECT DISTINCT home_team, period, market FROM mlb_odds"
fdset = set(fd.execute(q).fetchall())
dkset = set(dk.execute(q).fetchall())
print("Cards DK has but FD still lacks (expect only F7/F3 alts + F7/F3 ML if any):")
for r in sorted(dkset - fdset):
    print("  ", r)
PY
```
Confirm remaining gaps are only the known structural ones (F7/F3 alts, F7/F3 2-way ML — FD doesn't post these).

---

## Task 6: Documentation + memory

**Files:**
- Modify: `mlb_sgp/README.md`
- Modify: `/Users/callancapitolo/.claude-personal/projects/-Users-callancapitolo-NFLWork/memory/fd_event_page_tab_coverage_gap.md`

- [ ] **Step 1: Update `mlb_sgp/README.md`**

Find the FanDuel singles scraper section and add/adjust prose to state:
- The scraper fetches **both** the default tab and the SGP tab (`FD_TABS`) and merges markets/runners, because no single FD tab returns every market (F7/F3 lines are default-tab-only; FG-alts + F5 are SGP-tab-only).
- Market selection is a **keyword classifier** (`classify_market`, mirroring DK) with FD-specific junk exclusions and event-team team-total filtering — not the old exact-name whitelist. It covers FG/F5/F7/F3 main + alt wherever FD posts them and auto-picks-up new markets.
- FD does not post F7/F3 alts or F7/F3 2-way MLs (only a 3-way "Result" we skip), so those pills stay empty by vendor design.

- [ ] **Step 2: Mark the memory resolved**

Append to `fd_event_page_tab_coverage_gap.md`:
```
**RESOLVED 2026-05-23** (branch worktree-fanduel-full-market-coverage): scraper
now fetches both tabs + keyword classifier. See
docs/superpowers/plans/2026-05-23-fanduel-full-market-coverage.md.
```

- [ ] **Step 3: Commit**

```bash
git add mlb_sgp/README.md
git commit -m "docs(fd_singles): document both-tab fetch + keyword classifier

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

## Pre-Merge (do NOT merge without explicit user approval)

- [ ] Executive-engineer review of the full diff (`git diff main..HEAD`) per repo CLAUDE.md — data integrity (no dup writes / dedup correct), resource safety, edge cases (empty tab, pre-game stub events), dead code (whitelist fully removed), no secrets in logs.
- [ ] Confirm all of Task 5's gates passed from this branch (not a different branch).
- [ ] Present findings (ISSUES vs ACCEPTABLE RISKS) and get explicit approval.
- [ ] After merge: re-run the live `scraper_fanduel_singles.py` from `main`, confirm the two screenshot cards (Cardinals @ Reds Over 7.5, White Sox @ Giants Over 6.5) show FanDuel pills.
- [ ] Clean up: `git worktree remove` + `git branch -d`.
