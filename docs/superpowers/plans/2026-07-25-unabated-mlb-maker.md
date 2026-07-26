# Generalized Unabated Maker — MLB First — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generalize the WC maker into a sport-agnostic framework (interval ledger + shared totals-ladder adapter) and onboard Kalshi MLB run totals, live at tiny size the same slate the build lands.

**Architecture:** The runner/engine/gateway are already multi-sport (registry loop, per-adapter league feed, per-market book+trades capture). Three changes: (1) `maker/ledger.py` drops the goal grid `g=0..10` for exact interval evaluation derived from the fills' own lines; (2) soccer's ladder-pricing methods hoist into a shared `TotalsLadderAdapter`; (3) a new `sports/mlb.py` subclasses it with MLB constants discovered by read-only probes. Capture needs no new code — registering the adapter captures MLB books/trades on every tick.

**Tech Stack:** Python 3, pytest, DuckDB, requests. No new dependencies.

## Global Constraints

- Fair value is ALWAYS the devigged Unabated anchor — no proprietary modeling (standing user constraint).
- Fail closed everywhere: missing anchor / unparseable event → `[]` / skip, never interpolate.
- No `SELECT *` in SQL; list columns.
- Public ledger signatures (`worst_case`, `max_contracts`, `mark_to_fair`, `pnl_grid`) must not change — `maker/engine.py` and `maker/state.py` call them and are out of scope.
- Existing tests must stay green untouched except `tests/test_maker_ledger.py` (rewritten with the ledger) — `tests/test_soccer_adapter.py` passing unchanged is the hoist's parity proof.
- No temp files on disk; probes run via heredoc. Data stays in DuckDB.
- Go-live env (from spec §4): `MAKER_MODE=live MAKER_LIVE_ACK=1`, `MAKER_MAX_CONTRACTS` 2–5, `HARD_STOP_DOLLARS≈50`, kill file armed. First ticks always in shadow (`MAKER_MODE=shadow`) before the live flip, same session.
- All work on branch `worktree-unabated-mlb-maker`; merge to main only after full suite + pre-merge review + explicit user approval.

---

### Task 1: MLB recon probes (read-only; records the constants Tasks 4–6 consume)

**Files:**
- Create: none (heredoc probes only; findings land in this plan file's Appendix A and the Task-1 commit message)
- Modify: `docs/superpowers/plans/2026-07-25-unabated-mlb-maker.md` (fill Appendix A)

**Interfaces:**
- Produces: `MLB_LEAGUE_PREFIX` (e.g. `"lg5"` — verify), `MLB_TOTAL_SERIES` (e.g. `"KXMLBTOTAL"` — verify), the Kalshi event-title format string, sample team-name spellings from BOTH feeds, anchor `modifiedOn` cadence. Task 4's code below assumes `KXMLBTOTAL` and title format `"Yankees vs Red Sox: Total Runs"` — if the probes find different values, substitute them in Task 4's constants AND fixtures (that is the only permitted deviation from the code as written).

- [ ] **Step 1: Find the Unabated MLB league id**

Run from repo root (worktree). Sweeps candidate league ids on the v2 feed and prints a sample matchup per id; MLB is identifiable by team names (Yankees, Dodgers, …):

```bash
python3 - <<'EOF'
import requests
from unabated_edge import config
for lid in range(1, 40):
    url = config.UNABATED_V2_LEAGUE_URL.format(league_id=lid)
    try:
        r = requests.get(url, timeout=10)
    except requests.RequestException:
        continue
    if r.status_code != 200:
        continue
    raw = r.json()
    rows = raw if isinstance(raw, list) else raw.get("results", raw.get("data", []))
    for row in rows[:1]:
        print(lid, str(row)[:200])
EOF
```

Expected: one line whose teams are MLB clubs. Record `lgN` in Appendix A. Also record 2–3 full away/home name spellings exactly as the feed prints them.

- [ ] **Step 2: Confirm the anchor quotes MLB bt3 totals**

```bash
python3 - <<'EOF'
from unabated_edge import feed, config
MLB_ID, MLB_PREFIX = 5, "lg5"   # <-- substitute Step-1 values
st = feed.fetch_v2(MLB_ID, MLB_PREFIX)
eid = next(iter(st.events))
for ms in config.ANCHOR_SOURCE_IDS:
    over = st.lines.get(f"{eid}|0|{ms}|bt3")
    under = st.lines.get(f"{eid}|1|{ms}|bt3")
    print(ms, "over:", over and {k: over.get(k) for k in ("points", "modifiedOn")},
              "under:", under and {k: under.get(k) for k in ("points", "modifiedOn")},
              "alts:", len((over or {}).get("alternateLines") or []))
EOF
```

Expected: at least one anchor id (7/6/68) with `points` on both sides (e.g. 8.5) and alt ladders. Record which anchors quote MLB, the side convention sanity check (si0 line should devig to P(over)>0.5 when the total is low — spot-check one game against any sportsbook app), and `modifiedOn` recency vs now. If NO anchor id returns bt3 rows: STOP, report to user — the sport has no anchor and the project premise fails for MLB.

- [ ] **Step 3: Find the Kalshi MLB total series + title format**

```bash
python3 - <<'EOF'
from unabated_edge import kalshi
for series in ("KXMLBTOTAL", "KXMLBRUNS", "KXMLBGAMETOTAL"):
    try:
        evs = kalshi.list_events(series)
    except Exception as e:
        print(series, "ERROR", e); continue
    print(series, "events:", len(evs))
    for ev in evs[:2]:
        print("  title:", ev.get("title"))
        for mk in ev.get("markets", [])[:3]:
            print("    ", mk.get("ticker"), mk.get("strike_type"), mk.get("floor_strike"))
EOF
```

Expected: exactly one series returns events, `strike_type == "greater"`, `floor_strike` ending .5. Record series ticker, exact title format, and the team spelling style (nicknames vs city names). If none of the three guesses hit, list all sports series (`kalshi.list_events` variant or GET `/trade-api/v2/series?category=...` via the same client) and find the run-total ladder; record it.

- [ ] **Step 4: Record findings + commit**

Fill Appendix A of this file with the five recorded values. Then:

```bash
git add docs/superpowers/plans/2026-07-25-unabated-mlb-maker.md
git commit -m "recon(mlb): Unabated league id + Kalshi total series constants (Task 1)"
```

---

### Task 2: Generic interval ledger

**Files:**
- Modify: `unabated_edge/maker/ledger.py` (full rewrite of grid internals; same public names)
- Test: `tests/test_maker_ledger.py` (extend — keep existing tests, they must still pass)

**Interfaces:**
- Consumes: nothing new.
- Produces (unchanged names, generalized behavior): `pnl(fills, t: float) -> float`, `pnl_grid(fills) -> list[float]` (now one value per payoff interval, variable length), `worst_case(fills) -> float`, `max_contracts(fills, line, side, price, budget) -> int`, `mark_to_fair(fills, fair_by_line) -> float` (untouched). New helper `outcome_points(lines, extra_line=None) -> list[float]`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_maker_ledger.py`:

```python
import itertools
import random

from unabated_edge.maker import ledger


def _worst_case_bruteforce(fills, t_max=40):
    """Reference: exhaustive integer totals 0..t_max (superset of any old grid)."""
    return min(ledger.pnl(fills, t) for t in range(t_max + 1))


def test_outcome_points_cover_every_interval():
    pts = ledger.outcome_points([2.5, 8.5])
    # intervals: t<2.5, 2.5<t<8.5, t>8.5 -> one representative each
    assert pts == [2.0, 3.0, 9.0]


def test_outcome_points_extra_line_included():
    assert ledger.outcome_points([2.5], extra_line=8.5) == [2.0, 3.0, 9.0]
    assert ledger.outcome_points([], extra_line=None) == [0.0]


def test_high_line_no_longer_truncated():
    # A NO fill at 12.5 loses only when t >= 13 — beyond the old g=0..10 grid,
    # which valued this book at riskless +0.05/contract. Interval ledger sees it.
    fills = [(12.5, "no", 10, 0.95)]
    assert ledger.worst_case(fills) == -9.5  # t>=13: lose 0.95 x 10


def test_worst_case_matches_bruteforce_on_random_books():
    rng = random.Random(20260725)
    lines = [0.5, 1.5, 2.5, 6.5, 7.5, 8.5, 9.5, 12.5]
    for _ in range(200):
        fills = [(rng.choice(lines), rng.choice(["yes", "no"]),
                  rng.randint(1, 20), round(rng.uniform(0.05, 0.95), 2))
                 for _ in range(rng.randint(1, 6))]
        assert abs(ledger.worst_case(fills) - _worst_case_bruteforce(fills)) < 1e-9


def test_max_contracts_binds_beyond_old_grid():
    # Candidate NO at 12.5, price 0.95: worst case is -0.95/contract (t>=13).
    # budget 9.5 -> exactly 10 contracts. The old grid said unlimited (bound None).
    assert ledger.max_contracts([], 12.5, "no", 0.95, 9.5) == 10


def test_max_contracts_agrees_with_bruteforce_greedy():
    rng = random.Random(99)
    lines = [0.5, 2.5, 8.5, 12.5]
    for _ in range(100):
        fills = [(rng.choice(lines), rng.choice(["yes", "no"]),
                  rng.randint(1, 5), round(rng.uniform(0.1, 0.9), 2))
                 for _ in range(rng.randint(0, 3))]
        line, side = rng.choice(lines), rng.choice(["yes", "no"])
        price, budget = round(rng.uniform(0.1, 0.9), 2), rng.uniform(1, 20)
        n = ledger.max_contracts(fills, line, side, price, budget)
        # Property 1: the returned size respects the budget.
        assert _worst_case_bruteforce(fills + [(line, side, n, price)]) >= -budget - 1e-6
        # Property 2 (maximality): with price in (0,1) every candidate has a
        # losing region, so one more contract must break the budget.
        assert _worst_case_bruteforce(fills + [(line, side, n + 1, price)]) < -budget + 1e-6
```

(Soccer-era tests above in this file remain untouched and green — they are the parity proof for lines ≤ 5.5.)

- [ ] **Step 2: Run tests to verify the new ones fail**

Run: `python3 -m pytest tests/test_maker_ledger.py -v`
Expected: existing tests PASS; new tests FAIL (`outcome_points` not defined; `test_high_line_no_longer_truncated` asserts −9.5 but old grid returns +0.5).

- [ ] **Step 3: Rewrite ledger internals**

Replace `unabated_edge/maker/ledger.py` grid machinery (keep `mark_to_fair` verbatim):

```python
"""Per-game risk ledger: exact worst-case P&L over a totals ladder.

Sport-agnostic. Every rung settles on one number (the game's final total —
goals, runs, points); portfolio P&L is a piecewise-constant function of that
number whose only breakpoints are the fills' lines. Evaluating one
representative outcome per interval is therefore exact — no fixed grid, no
truncation above a top rung.

A fill is (line, side, contracts, price): side is relative to the Over
market — "yes" wins when t > line, "no" wins when t <= line. Prices in
dollars; P&L per contract is (1 - price) on a win, -price on a loss.
"""


def pnl(fills, t: float) -> float:
    total = 0.0
    for line, side, contracts, price in fills:
        win = (t > line) if side == "yes" else (t <= line)
        total += contracts * ((1.0 - price) if win else -price)
    return total


def outcome_points(lines, extra_line=None) -> list[float]:
    """One representative final-total per payoff interval of the given lines
    (plus an optional candidate line): below the lowest strike, between each
    consecutive pair, above the highest."""
    ls = sorted(set(lines) | ({extra_line} if extra_line is not None else set()))
    if not ls:
        return [0.0]
    return [ls[0] - 0.5] + [line + 0.5 for line in ls]


def pnl_grid(fills) -> list[float]:
    """P&L per payoff interval (variable length — one entry per interval)."""
    return [pnl(fills, t) for t in outcome_points([f[0] for f in fills])]


def worst_case(fills) -> float:
    if not fills:
        return 0.0
    return min(pnl_grid(fills))


def max_contracts(fills, line: float, side: str, price: float, budget: float) -> int:
    """Largest integer n >= 0 such that adding (line, side, n, price) keeps
    worst_case >= -budget.

    Closed form per interval: base_t + n*unit_t >= -budget binds only where
    the candidate loses (unit_t < 0). Intervals derive from existing fills'
    lines PLUS the candidate's line. Fail closed if already beyond budget."""
    if budget <= 0:
        return 0
    points = outcome_points([f[0] for f in fills], extra_line=line)
    base = [pnl(fills, t) for t in points]
    if min(base) < -budget - 1e-9:
        return 0
    bound = None
    for t, base_t in zip(points, base):
        win = (t > line) if side == "yes" else (t <= line)
        unit = (1.0 - price) if win else -price
        if unit < 0:
            allowed = (budget + base_t) / (-unit)
            bound = allowed if bound is None else min(bound, allowed)
    return int(bound + 1e-9) if bound is not None else 0
```

`G_MAX` is deleted; grep for external references first (`grep -rn G_MAX unabated_edge/ tests/`) and update any test that imports it to use `outcome_points` semantics instead.

- [ ] **Step 4: Run the ledger + engine + state suites**

Run: `python3 -m pytest tests/test_maker_ledger.py tests/test_maker_engine.py tests/test_maker_state.py -v`
Expected: ALL PASS (engine/state exercise `worst_case`/`pnl_grid`/`max_contracts` through their public behavior; identical results for soccer-range lines).

- [ ] **Step 5: Commit**

```bash
git add unabated_edge/maker/ledger.py tests/test_maker_ledger.py
git commit -m "feat(maker): sport-agnostic interval ledger replaces goal grid

Worst-case/max-contracts now evaluate one representative outcome per
payoff interval derived from the fills' own lines - exact for any totals
ladder (soccer parity preserved; fixes silent truncation above the old
g=0..10 grid, e.g. MLB 12.5 rungs)."
```

---

### Task 3: Hoist `TotalsLadderAdapter` (pure refactor; soccer tests are the gate)

**Files:**
- Create: `unabated_edge/sports/totals.py`
- Modify: `unabated_edge/sports/soccer.py`

**Interfaces:**
- Consumes: `SportAdapter`, `Candidate` from `sports/base.py`; `pricing`, `config`, `feed.line_american_price`.
- Produces: `class TotalsLadderAdapter(SportAdapter)` with the exact method bodies currently in `Soccer`: `_side_prices`, `_anchor_ladder`, `_anchor_totals`, `fair_ladder`, `price_event`, plus module-level `_older_mo`. Subclasses implement only `canon_team`, `kalshi_series`, `event_teams` and set `sport`, `league_prefix`.

- [ ] **Step 1: Create `sports/totals.py`**

Move (verbatim — this is a cut-and-paste refactor, not a rewrite) `_older_mo` and the five methods listed above from `soccer.py` into:

```python
"""Shared pricing for sports whose Kalshi market is an Over-ladder on a single
final total (goals, runs, points), anchored by Unabated bt3 over/under.

Subclasses provide only identity: sport, league_prefix, canon_team,
kalshi_series, event_teams. All ladder devigging and rung matching lives here."""
from unabated_edge.sports.base import SportAdapter, Candidate
from unabated_edge import pricing, config
from unabated_edge.feed import line_american_price


def _older_mo(a, b):
    ...  # moved verbatim from soccer.py


class TotalsLadderAdapter(SportAdapter):
    _side_prices = ...      # moved verbatim (staticmethod)
    _anchor_ladder = ...    # moved verbatim
    _anchor_totals = ...    # moved verbatim
    fair_ladder = ...       # moved verbatim
    price_event = ...       # moved verbatim
```

(The `...` above are moves, not stubs — the bodies are exactly the current `soccer.py` lines 6–10, 57–131.)

- [ ] **Step 2: Slim `soccer.py`**

`soccer.py` retains `_ALIASES`, `WC_TOTAL_SERIES`, and:

```python
from unabated_edge.sports.totals import TotalsLadderAdapter


class Soccer(TotalsLadderAdapter):
    sport = "soccer"
    league_prefix = "lg21"

    def canon_team(self, name: str) -> str: ...        # unchanged body
    def kalshi_series(self) -> str: ...                # unchanged body
    def event_teams(self, kalshi_event) -> frozenset: ...  # unchanged body
```

- [ ] **Step 3: Run the full suite**

Run: `python3 -m pytest tests/ -q`
Expected: ALL PASS with zero test edits — `test_soccer_adapter.py` green against the hoisted code is the parity proof.

- [ ] **Step 4: Commit**

```bash
git add unabated_edge/sports/totals.py unabated_edge/sports/soccer.py
git commit -m "refactor(sports): hoist totals-ladder pricing into TotalsLadderAdapter

Pure move; soccer tests pass unchanged. Sport N+1 now implements only
canon_team/kalshi_series/event_teams."
```

---

### Task 4: MLB adapter + registration

**Files:**
- Create: `unabated_edge/sports/mlb.py`
- Modify: `unabated_edge/sports/registry.py`
- Test: `tests/test_mlb_adapter.py`

**Interfaces:**
- Consumes: `TotalsLadderAdapter` (Task 3); Task 1's Appendix-A constants — substitute them below if they differ from the assumed `"lg5"` / `"KXMLBTOTAL"` / `"Yankees vs Red Sox: Total Runs"`.
- Produces: `class Mlb(TotalsLadderAdapter)`; `registry.ADAPTERS == [Soccer(), Mlb()]`.

- [ ] **Step 1: Write the failing tests**

`tests/test_mlb_adapter.py` (mirror `test_soccer_adapter.py`'s fixture style — read it first; reuse its feed-state builder pattern rather than inventing one):

```python
from unabated_edge.sports.mlb import Mlb
from unabated_edge.sports import registry


def test_canon_team_city_and_nickname_collapse():
    m = Mlb()
    assert m.canon_team("New York Yankees") == "Yankees"
    assert m.canon_team("Yankees") == "Yankees"
    assert m.canon_team("St. Louis Cardinals") == "Cardinals"
    assert m.canon_team("Athletics") == "Athletics"


def test_event_teams_parses_kalshi_title():
    # Title format verified in Task 1; update fixture if recon differs.
    m = Mlb()
    ev = {"title": "Yankees vs Red Sox: Total Runs"}
    assert m.event_teams(ev) == frozenset({"Yankees", "Red Sox"})


def test_event_teams_malformed_title_fails_closed():
    assert Mlb().event_teams({"title": "All-Star Game Total Runs"}) == frozenset()
    assert Mlb().event_teams({}) == frozenset()


def test_mlb_registered():
    assert any(a.sport == "mlb" for a in registry.ADAPTERS)
    assert len({a.league_prefix for a in registry.ADAPTERS}) == len(registry.ADAPTERS)


def test_price_event_prices_matching_rungs():
    # Inherited TotalsLadderAdapter path with MLB-typical lines: anchor 8.5
    # both sides -> yes+no candidates at the 8.5 rung, other rungs skipped.
    # Build state via the same helper test_soccer_adapter.py uses (copy its
    # builder into a shared fixture in conftest.py if it is module-local).
    ...
```

The last test's body: replicate `test_soccer_adapter.py`'s `price_event` happy-path test with `points=8.5`, an `Mlb()` adapter, and a Kalshi market dict `{"ticker": "KXMLBTOTAL-X", "strike_type": "greater", "floor_strike": 8.5}`. If the soccer test's state builder is module-local, move it to `tests/conftest.py` as a shared fixture in this task (mechanical move, both test files keep passing).

- [ ] **Step 2: Run to verify failure**

Run: `python3 -m pytest tests/test_mlb_adapter.py -v`
Expected: FAIL — `No module named 'unabated_edge.sports.mlb'`.

- [ ] **Step 3: Implement `sports/mlb.py`**

```python
"""MLB run-totals adapter: Unabated bt3 anchor -> Kalshi MLB total ladder.

Constants below (league prefix, series, title format) were live-verified in
Task 1 recon (plan Appendix A) - update there if they drift."""
from unabated_edge.sports.totals import TotalsLadderAdapter

MLB_TOTAL_SERIES = "KXMLBTOTAL"   # Task 1 Appendix A

# City-prefixed full names (Unabated style) -> nickname (Kalshi style).
_MLB_TEAMS = {
    "arizona diamondbacks": "Diamondbacks", "atlanta braves": "Braves",
    "baltimore orioles": "Orioles", "boston red sox": "Red Sox",
    "chicago cubs": "Cubs", "chicago white sox": "White Sox",
    "cincinnati reds": "Reds", "cleveland guardians": "Guardians",
    "colorado rockies": "Rockies", "detroit tigers": "Tigers",
    "houston astros": "Astros", "kansas city royals": "Royals",
    "los angeles angels": "Angels", "los angeles dodgers": "Dodgers",
    "miami marlins": "Marlins", "milwaukee brewers": "Brewers",
    "minnesota twins": "Twins", "new york mets": "Mets",
    "new york yankees": "Yankees", "athletics": "Athletics",
    "oakland athletics": "Athletics", "philadelphia phillies": "Phillies",
    "pittsburgh pirates": "Pirates", "san diego padres": "Padres",
    "san francisco giants": "Giants", "seattle mariners": "Mariners",
    "st. louis cardinals": "Cardinals", "st louis cardinals": "Cardinals",
    "tampa bay rays": "Rays", "texas rangers": "Rangers",
    "toronto blue jays": "Blue Jays", "washington nationals": "Nationals",
}
# Nicknames map to themselves so Kalshi-side names pass through canon.
_MLB_TEAMS.update({v.lower(): v for v in list(_MLB_TEAMS.values())})


class Mlb(TotalsLadderAdapter):
    sport = "mlb"
    league_prefix = "lg5"   # Task 1 Appendix A

    def canon_team(self, name: str) -> str:
        key = (name or "").strip().lower()
        return _MLB_TEAMS.get(key, (name or "").strip())

    def kalshi_series(self) -> str:
        return MLB_TOTAL_SERIES

    def event_teams(self, kalshi_event: dict) -> frozenset:
        """Title format (Task 1): "Yankees vs Red Sox: Total Runs"."""
        title = (kalshi_event.get("title") or "")
        before_colon = title.split(":")[0]
        parts = before_colon.split(" vs ")
        if len(parts) != 2:
            return frozenset()
        return frozenset({self.canon_team(parts[0].strip()),
                          self.canon_team(parts[1].strip())})
```

Registry becomes:

```python
from unabated_edge.sports.soccer import Soccer
from unabated_edge.sports.mlb import Mlb
ADAPTERS = [Soccer(), Mlb()]
def league_prefixes(): return {a.league_prefix for a in ADAPTERS}
def by_league(prefix): return next((a for a in ADAPTERS if a.league_prefix == prefix), None)
```

- [ ] **Step 4: Run the full suite**

Run: `python3 -m pytest tests/ -q`
Expected: ALL PASS (registry change is additive; runner tests iterate ADAPTERS generically).

- [ ] **Step 5: Live pairing smoke (read-only, no orders)**

```bash
python3 - <<'EOF'
from unabated_edge import feed, kalshi, mapping
from unabated_edge.sports.mlb import Mlb
a = Mlb()
st = feed.fetch_v2(a.league_id, a.league_prefix)
evs = kalshi.list_events(a.kalshi_series())
print("unabated events:", len(st.events), "kalshi events:", len(evs))
pairs = mapping.pair_events(a, st, evs)
print("paired:", len(pairs))
EOF
```

Expected: paired count ≈ today's slate size (≥ 10 on a full slate). If pairing is < half the slate, diff the two feeds' team spellings and extend `_MLB_TEAMS` before proceeding. (Check `mapping.pair_events`'s real signature first and adapt the probe — it logs unmatched events by design.)

- [ ] **Step 6: Commit**

```bash
git add unabated_edge/sports/mlb.py unabated_edge/sports/registry.py tests/test_mlb_adapter.py tests/conftest.py
git commit -m "feat(sports): MLB run-totals adapter on TotalsLadderAdapter

Constants live-verified in Task 1 recon; pairing smoke-tested against
today's slate."
```

---

### Task 5: MLB readout queries + disk guard

**Files:**
- Create: `unabated_edge/mlb_readout.sql`

**Interfaces:**
- Consumes: `book_snapshots`, `kalshi_trades`, `line_snapshots` tables in `unabated_edge_market.duckdb` (already written per-sport by the runner; MLB rows appear as `sport='mlb'` once Task 4 is running).
- Produces: a query file the user runs ad hoc (`duckdb -readonly` against a COPY of the market DB — never the live file, per the WC WAL-hang caveat).

- [ ] **Step 1: Write the readout**

`unabated_edge/mlb_readout.sql` — header comment states: run against a copy (`cp unabated_edge_market.duckdb /tmp/ro.duckdb` is forbidden by repo rules → use `$CLAUDE_JOB_DIR/tmp` or a user-chosen path), read-only. Queries, each with a one-line comment stating the decision it informs:

```sql
-- 1. Crowd spread distribution per rung (is this a WC-like 1c book?)
SELECT market_ticker,
       round(avg(best_yes_ask - best_yes_bid), 4) AS avg_spread,
       round(quantile_cont(best_yes_ask - best_yes_bid, 0.5), 4) AS med_spread,
       count(*) AS snaps
FROM book_snapshots
WHERE sport = 'mlb'
GROUP BY market_ticker
ORDER BY med_spread DESC;

-- 2. Depth at touch (queue we'd be joining)
SELECT market_ticker,
       round(avg(yes_bid_qty), 0) AS avg_bid_depth,
       round(avg(yes_ask_qty), 0) AS avg_ask_depth
FROM book_snapshots
WHERE sport = 'mlb'
GROUP BY market_ticker
ORDER BY avg_bid_depth DESC;

-- 3. Share of trades within 1c of mid (the WC kill-shot number was 99.99%)
WITH t AS (
  SELECT tr.market_ticker, tr.price,
         (bs.best_yes_bid + bs.best_yes_ask) / 2 AS mid
  FROM kalshi_trades tr
  ASOF JOIN book_snapshots bs
    ON tr.market_ticker = bs.market_ticker AND tr.created_ts >= bs.captured_ts
  WHERE tr.sport = 'mlb'
)
SELECT count(*) AS trades,
       round(100.0 * avg(CASE WHEN abs(price - mid) <= 0.01 THEN 1 ELSE 0 END), 2)
         AS pct_within_1c
FROM t;

-- 4. Anchor coverage: rungs with a fair vs rungs Kalshi lists
SELECT sport, count(DISTINCT market_ticker) AS rungs_snapshotted
FROM book_snapshots WHERE sport = 'mlb' GROUP BY sport;

-- 5. Retention prune (disk guard; run when the market DB grows past ~2 GB)
-- DELETE FROM book_snapshots WHERE sport='mlb' AND captured_ts < now() - INTERVAL 14 DAY;
-- DELETE FROM kalshi_trades  WHERE sport='mlb' AND created_ts  < now() - INTERVAL 14 DAY;
```

Column names above must be checked against `storage.py`'s actual `CREATE TABLE` statements (`book_snapshots` has 13 columns; use the real bid/ask/qty/ts names — adjust the SQL, not the schema). No `SELECT *`.

- [ ] **Step 2: Verify the queries parse against the real schema**

```bash
df -h / | tail -1   # record free space; flag to user if < 15 GB free
python3 - <<'EOF'
import duckdb, pathlib, re
sql = pathlib.Path("unabated_edge/mlb_readout.sql").read_text()
con = duckdb.connect("unabated_edge/unabated_edge_market.duckdb", read_only=True)
for stmt in [s for s in sql.split(";") if s.strip() and not s.strip().startswith("--")]:
    con.execute(stmt)   # empty result sets are fine; parse/bind errors are not
print("all statements parse")
EOF
```

Expected: `all statements parse` (0 MLB rows yet is fine). If the live DB is WAL-locked, copy it to `$CLAUDE_JOB_DIR/tmp` first.

- [ ] **Step 3: Commit**

```bash
git add unabated_edge/mlb_readout.sql
git commit -m "feat(mlb): capture readout queries (spread/depth/flow-vs-mid) + prune guard"
```

---

### Task 6: Docs + full suite + launch runbook

**Files:**
- Modify: `unabated_edge/README.md`, `CLAUDE.md` (repo root)

**Interfaces:**
- Consumes: everything above.
- Produces: the merge-ready branch; the runbook the user (or this session, on approval) executes to go live.

- [ ] **Step 1: README**

Add to `unabated_edge/README.md`:
- **Multi-sport onboarding** section: sport N+1 = (1) Task-1-style recon for league id/series/title, (2) subclass `TotalsLadderAdapter` implementing `canon_team`/`kalshi_series`/`event_teams`, (3) one registry line, (4) fixture tests. Point at `sports/mlb.py` as the template.
- **MLB** section: constants + Appendix-A provenance, readout usage (`mlb_readout.sql`, run against a copy), the 14-day prune, disk check.
- **Launch runbook** (verbatim commands):

```bash
# 0. Preflight: disk + env (repo-root .env carries KALSHI keys + ANCHOR_STALE_SEC)
df -h / | tail -1                      # need comfortable headroom; WC died at 98%
ls .env && grep -c KALSHI .env

# 1. First-tick shadow check (same session; no standalone shadow day)
MAKER_MODE=shadow python3 -m unabated_edge.runner
# watch heartbeat: kalshi_events>0 for BOTH sports, candidates_recent>0,
# maker= shows shadow quotes on MLB tickers. Ctrl-C after 2-3 clean ticks.

# 2. Flip live, same slate (leashes ON - spec Review Pack decision)
MAKER_MODE=live MAKER_LIVE_ACK=1 MAKER_MAX_CONTRACTS=3 HARD_STOP_DOLLARS=50 \
  python3 -m unabated_edge.runner

# Kill switch at any time:
touch unabated_edge/.kill
```

- Restart caveat carries over verbatim from the WC README: never restart mid-game while holding inventory (fills ledger rebuilds blank).

- [ ] **Step 2: Root CLAUDE.md**

In the `unabated_edge` bullet: replace the soccer-only description with "sport-agnostic totals maker (soccer WC shipped, MLB run totals live-tiny 2026-07); generic interval ledger (no goal-grid); onboarding = `TotalsLadderAdapter` subclass + registry line; capture readout in `mlb_readout.sql`". Keep it to 2–3 lines; the README carries the detail.

- [ ] **Step 3: Full suite + review gate**

Run: `python3 -m pytest tests/ -q`
Expected: ALL PASS (was 131 pre-branch; now more).

Then run the pre-merge review checklist from root CLAUDE.md against `git diff main..HEAD`, present findings + this runbook to the user, and **stop for explicit merge + go-live approval**. Never merge or launch without it.

- [ ] **Step 4: Commit**

```bash
git add unabated_edge/README.md CLAUDE.md
git commit -m "docs(unabated_edge): multi-sport onboarding + MLB launch runbook"
```

---

## Appendix A — Task 1 recon findings (filled by Task 1)

| Constant | Value | Evidence |
|---|---|---|
| `MLB_LEAGUE_PREFIX` | _recorded in Task 1_ | v2 feed sweep output |
| `MLB_TOTAL_SERIES` | _recorded in Task 1_ | Kalshi series probe |
| Kalshi title format | _recorded in Task 1_ | sample titles |
| Anchor ids quoting MLB bt3 | _recorded in Task 1_ | Step 2 output |
| Anchor `modifiedOn` cadence | _recorded in Task 1_ | Step 2 output |
