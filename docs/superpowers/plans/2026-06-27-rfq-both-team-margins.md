# RFQ taker — price both teams' margin markets Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the book-only RFQ taker price *both* teams' Kalshi margin markets correctly (away-team spread legs included), by teaching the SGP cache to hold both teams' signed-line grids and routing every leg to the right cell.

**Architecture:** The cache `spread_line` is home-perspective signed. A negative line (−1.5) is the home-favorite grid; the positive mirror (+1.5) is the away-favorite grid. The scrapers already fetch by signed line — the only gaps are (a) the target enumeration never emits positive lines, and (b) the leg→region mapping always signs the line negative. Both fixes are RFQ-local: enumeration changes are gated behind a default-off `both_teams` flag so the live MM bot (which shares `kalshi_common.sgp_runner`) stays byte-identical, and leg pricing routes off `region.spread_line` instead of the shared `_spread_line_from_legs`.

**Tech Stack:** Python 3, DuckDB, pandas, pytest. System `python3` in the worktree (no venv — has duckdb/pandas/numpy/scipy).

## Global Constraints

- Work on the existing `worktree-rfq-remove-model` branch (this worktree). No new branch.
- **No merge to main without explicit user approval** + executive-engineer pre-merge review.
- **Do not change MM-bot behavior.** `kalshi_common.sgp_runner.{sgp_cycle, enumerate_kalshi_targets, _fetch_kalshi_spread_lines}` and `kalshi_common.leg_types._spread_line_from_legs` are shared with `kalshi_mlb_mm`. Any change to them must keep the default (no-arg) behavior byte-identical; the RFQ bot opts into new behavior explicitly.
- Tests run with: `cd /Users/callancapitolo/NFLWork/.claude/worktrees/worktree-rfq-remove-model && python3 -m pytest kalshi_mlb_rfq/tests/ -q` (run from repo root so `kalshi_common`/`kalshi_mlb_rfq` import).
- TDD: failing test first, minimal impl, green, commit. Frequent commits.
- `spread_line` sign convention (home-perspective): home margin market → `-(n-0.5)`; away margin market → `+(n-0.5)`. `n` is the KXMLBSPREAD ticker suffix integer (n=2 ⇒ ±1.5).

---

### Task 1: Signed-line target enumeration (shared code, RFQ opt-in flag)

**Files:**
- Modify: `kalshi_common/sgp_runner.py` — `_fetch_kalshi_spread_lines` (64-88), `enumerate_kalshi_targets` (174-218, the `_fetch_kalshi_spread_lines` call at 204), `sgp_cycle` (380-401, the `enumerate_kalshi_targets` call at 401)
- Modify: `kalshi_mlb_rfq/main.py` — the two `sgp_runner.sgp_cycle(...)` calls (≈1755, ≈1810) pass `both_teams=True`
- Test: `kalshi_mlb_rfq/tests/test_sgp_runner.py`

**Interfaces:**
- Produces: `_fetch_kalshi_spread_lines(suffix, home_code=None, both_teams=False) -> list[(line, who)]`; `enumerate_kalshi_targets(both_teams=False) -> list[TargetLine]`; `sgp_cycle(..., both_teams=False)`. Default (no flag) is byte-identical to today. With `both_teams=True`, emits one signed line per team's margin market (no dedup): home→`-(n-0.5)`, away→`+(n-0.5)`.
- Consumes: nothing new.

- [ ] **Step 1: Write the failing test** (append to `kalshi_mlb_rfq/tests/test_sgp_runner.py`)

```python
def test_fetch_spread_lines_default_is_home_favorite_deduped(monkeypatch):
    """Default (both_teams=False) preserves legacy MM behavior: one
    home-favorite -(n-0.5) line per |line|, deduped, who='home'."""
    from kalshi_common import sgp_runner
    markets = {"markets": [
        {"ticker": "KXMLBSPREAD-GAME-HOU2"},  # home n=2
        {"ticker": "KXMLBSPREAD-GAME-PIT2"},  # away n=2 (same |line|)
        {"ticker": "KXMLBSPREAD-GAME-HOU3"},
    ]}
    monkeypatch.setattr(sgp_runner.auth_client, "api",
                        lambda *a, **k: (200, markets, {}))
    out = sgp_runner._fetch_kalshi_spread_lines("GAME", home_code="HOU")
    assert sorted(out) == [(-2.5, "home"), (-1.5, "home")]  # deduped, home only


def test_fetch_spread_lines_both_teams_emits_signed_per_team(monkeypatch):
    """both_teams=True: home margin → negative, away margin → positive,
    no dedup."""
    from kalshi_common import sgp_runner
    markets = {"markets": [
        {"ticker": "KXMLBSPREAD-GAME-HOU2"},  # home favorite by 1.5
        {"ticker": "KXMLBSPREAD-GAME-PIT2"},  # away favorite by 1.5
    ]}
    monkeypatch.setattr(sgp_runner.auth_client, "api",
                        lambda *a, **k: (200, markets, {}))
    out = sgp_runner._fetch_kalshi_spread_lines(
        "GAME", home_code="HOU", both_teams=True)
    assert sorted(out) == [(-1.5, "home"), (1.5, "away")]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd /Users/callancapitolo/NFLWork/.claude/worktrees/worktree-rfq-remove-model && python3 -m pytest kalshi_mlb_rfq/tests/test_sgp_runner.py -q -k "spread_lines"`
Expected: FAIL — `_fetch_kalshi_spread_lines()` got an unexpected keyword `home_code` / `both_teams`.

- [ ] **Step 3: Implement** — replace `_fetch_kalshi_spread_lines` body (`kalshi_common/sgp_runner.py:64-88`)

```python
def _fetch_kalshi_spread_lines(suffix: str, home_code: str | None = None,
                               both_teams: bool = False) -> list[tuple[float, str]]:
    status, body, _ = auth_client.api(
        "GET", f"/markets?event_ticker=KXMLBSPREAD-{suffix}&limit=50")
    if status != 200 or not isinstance(body, dict):
        return []
    out = []
    seen = set()
    for m in body.get("markets", []):
        ticker = m.get("ticker", "")
        prefix = f"KXMLBSPREAD-{suffix}-"
        if not ticker.startswith(prefix):
            continue
        spread_part = ticker[len(prefix):]
        digits = "".join(c for c in spread_part if c.isdigit())
        team_chars = "".join(c for c in spread_part if not c.isdigit())
        if not digits or not team_chars:
            continue
        n = int(digits)
        if both_teams:
            # Home-perspective signed line per covering team:
            #   home margin market → -(n-0.5);  away margin market → +(n-0.5).
            # No dedup — the two teams are distinct grids (one per sign).
            who = "home" if (home_code is not None and team_chars == home_code) else "away"
            line = -(n - 0.5) if who == "home" else (n - 0.5)
            out.append((line, who))
        else:
            # Legacy/MM behavior: one home-favorite grid per |line|, deduped.
            line = -(n - 0.5)
            key = round(line, 1)
            if key in seen:
                continue
            seen.add(key)
            out.append((line, "home"))
    return out
```

- [ ] **Step 4: Thread the flag through the two callers** (`kalshi_common/sgp_runner.py`)

At line 204, change:
```python
        spreads = _fetch_kalshi_spread_lines(suffix)
```
to:
```python
        spreads = _fetch_kalshi_spread_lines(
            suffix, home_code=home_code, both_teams=both_teams)
```
Change `enumerate_kalshi_targets` signature (line 174):
```python
def enumerate_kalshi_targets(both_teams: bool = False) -> list[TargetLine]:
```
Change `sgp_cycle` signature (line 380) to add `both_teams: bool = False,` (after `service=None,`), and its `enumerate_kalshi_targets()` call (line 401):
```python
    targets = enumerate_kalshi_targets(both_teams=both_teams)
```

- [ ] **Step 5: RFQ opts in** (`kalshi_mlb_rfq/main.py`)

In BOTH `sgp_runner.sgp_cycle(...)` calls (≈1755 startup warm + ≈1810 loop), add `both_teams=True`:
```python
                    rcs = sgp_runner.sgp_cycle(
                        bot_market_db=str(config.BOT_MARKET_DB),
                        service=sgp_service,
                        both_teams=True,
                    )
```

- [ ] **Step 6: Run tests to verify pass (incl. existing sgp_runner tests for MM-parity)**

Run: `cd /Users/callancapitolo/NFLWork/.claude/worktrees/worktree-rfq-remove-model && python3 -m pytest kalshi_mlb_rfq/tests/test_sgp_runner.py -q`
Expected: PASS (new + all existing — existing tests exercise the default path and must stay green, proving MM byte-parity).

- [ ] **Step 7: Commit**

```bash
git add kalshi_common/sgp_runner.py kalshi_mlb_rfq/main.py kalshi_mlb_rfq/tests/test_sgp_runner.py
git commit -m "feat(rfq): enumerate both teams' signed margin lines (both_teams flag, MM default unchanged)"
```

---

### Task 2: Team-signed leg routing (RFQ-local pricing)

**Files:**
- Modify: `kalshi_mlb_rfq/main.py` — `_combo_region_from_legs` (line 557), `_fresh_blended_fair` (615-625), candidate loop `_load_book_fairs` call (line 1695)
- Test: `kalshi_mlb_rfq/tests/test_book_only_pricing.py`

**Interfaces:**
- Consumes: `correlation.ComboRegion` (unchanged shape), `SpreadLeg.team_is_home`.
- Produces: `_combo_region_from_legs(typed_legs)` now returns `spread_line = -(line_n-0.5)` for home-margin legs and `+(line_n-0.5)` for away-margin legs. Pricing reads `region.spread_line` (signed) for the cache lookup.

- [ ] **Step 1: Write the failing test** (append to `kalshi_mlb_rfq/tests/test_book_only_pricing.py`)

```python
import kalshi_mlb_rfq.main as m
from kalshi_common import fair_value


def _region(team_is_home, side, line_n=2, total_n=8, total_side="yes"):
    typed = [
        fair_value.SpreadLeg(team_is_home=team_is_home, line_n=line_n, side=side),
        fair_value.TotalLeg(line_n=total_n, side=total_side),
    ]
    return m._combo_region_from_legs(typed)


def test_routing_table_signs_line_by_team():
    """The 4-row routing table: grid sign by whose margin market it is,
    cell by yes/no. n=2 ⇒ |line|=1.5."""
    # home margin → negative grid
    assert _region(True, "yes").spread_line == -1.5
    assert _region(True, "yes").spread_side == "home"
    assert _region(True, "no").spread_line == -1.5
    assert _region(True, "no").spread_side == "away"
    # away margin → positive grid
    assert _region(False, "yes").spread_line == 1.5
    assert _region(False, "yes").spread_side == "away"
    assert _region(False, "no").spread_line == 1.5
    assert _region(False, "no").spread_side == "home"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd /Users/callancapitolo/NFLWork/.claude/worktrees/worktree-rfq-remove-model && python3 -m pytest kalshi_mlb_rfq/tests/test_book_only_pricing.py -q -k routing_table`
Expected: FAIL — away-margin `spread_line` is `-1.5`, not `1.5`.

- [ ] **Step 3: Implement the sign fix** (`kalshi_mlb_rfq/main.py:557`)

Change:
```python
    spread_line = -(spread_leg.line_n - 0.5)   # home-perspective, e.g. n=2 → -1.5
```
to:
```python
    # Home-perspective signed line: home margin → -(n-0.5); away margin → +(n-0.5).
    # The sign selects the grid (home-favorite vs away-favorite); spread_side
    # selects the cell within it.
    spread_line = (-(spread_leg.line_n - 0.5) if spread_leg.team_is_home
                   else (spread_leg.line_n - 0.5))
```
Also update the docstring line ~536 ("matching `_spread_line_from_legs: -(line_n - 0.5)`") to note the sign now depends on `team_is_home`.

- [ ] **Step 4: Route pricing off `region.spread_line`** (so the signed line reaches the cache)

In `_fresh_blended_fair` (≈615-625), delete the standalone `spread_line = _spread_line_from_legs(legs)` line and change the `_load_book_fairs` call to use the region:
```python
    total_line = _total_line_from_legs(legs)
    typed = [_leg_dict_to_typed(l, game_id) for l in legs]
    if any(l is None for l in typed):
        return None, {}
    region = _combo_region_from_legs(typed)
    if region is None:
        return None, {}
    book_fairs = _load_book_fairs(game_id, region.spread_line, region.total_line,
                                  region.spread_side, region.total_side)
```
In the candidate loop, change the `_load_book_fairs` call (line 1695) to use the region's signed line:
```python
            books = _load_book_fairs(game_id, region.spread_line, region.total_line,
                                     region.spread_side, region.total_side)
```
Leave the loop's `spread_line`/`total_line` locals (from `_spread_line_from_legs`, line 1680-1681) as-is — they are only used for the diagnostic `_emit_candidate_event("rejected_no_mapping", ...)` calls, which fire before a region exists. (Do NOT change shared `_spread_line_from_legs`; the MM bot uses it.)

- [ ] **Step 5: Run tests to verify pass**

Run: `cd /Users/callancapitolo/NFLWork/.claude/worktrees/worktree-rfq-remove-model && python3 -m pytest kalshi_mlb_rfq/tests/test_book_only_pricing.py -q`
Expected: PASS (new routing test + existing book-only/Task-9 tests).

- [ ] **Step 6: Commit**

```bash
git add kalshi_mlb_rfq/main.py kalshi_mlb_rfq/tests/test_book_only_pricing.py
git commit -m "fix(rfq): sign spread_line by team so away-margin legs route to the away grid"
```

---

### Task 3: Correlation under signed grids (verify + document)

**Files:**
- Modify (docstring only, unless a test fails): `kalshi_mlb_rfq/correlation.py` — `joint_prob` (35-57)
- Test: `kalshi_mlb_rfq/tests/test_correlation.py`

**Interfaces:**
- Consumes: `ComboRegion` with signed `spread_line`.
- Produces: no signature change. Confirms same-side joints read the correct signed grid; opposite-side (cross-grid) pairs remain `None` → caller's ρ=1 fallback (conservative).

- [ ] **Step 1: Write the failing/characterization test** (append to `kalshi_mlb_rfq/tests/test_correlation.py`)

```python
from kalshi_mlb_rfq import correlation
from kalshi_mlb_rfq.correlation import ComboRegion


def test_joint_two_away_margin_combos_reads_positive_grid():
    """Two away-margin combos (spread_line +1.5 and +2.5, both away cover,
    both over) → tighter = +2.5 away cell. grid_lookup must be called with
    the POSITIVE signed line."""
    calls = []
    def fake_lookup(spread, total, sside, tside):
        calls.append((spread, total, sside, tside))
        return 0.10
    a = ComboRegion("away", 1.5, "over", 7.5)
    b = ComboRegion("away", 2.5, "over", 8.5)
    j = correlation.joint_prob(a, b, p_a=0.30, p_b=0.12, grid_lookup=fake_lookup)
    assert calls == [(2.5, 8.5, "away", "over")]   # tighter: max line, max total
    assert j is not None


def test_joint_opposite_grid_pair_is_none_fallback():
    """A home-margin combo (−1.5) and an away-margin combo (+1.5) have
    different spread_side → None (caller uses ρ=1, conservative)."""
    a = ComboRegion("home", -1.5, "over", 7.5)   # home wins by 2+
    b = ComboRegion("away", 1.5, "over", 7.5)    # away wins by 2+
    j = correlation.joint_prob(a, b, p_a=0.30, p_b=0.15,
                               grid_lookup=lambda *x: 0.5)
    assert j is None
```

- [ ] **Step 2: Run test**

Run: `cd /Users/callancapitolo/NFLWork/.claude/worktrees/worktree-rfq-remove-model && python3 -m pytest kalshi_mlb_rfq/tests/test_correlation.py -q -k "away_margin or opposite_grid"`
Expected: PASS as-is (the half-plane math is already sign-correct: `away` → `max(line)`, and differing `spread_side` → `None`). If it FAILS, fix `joint_prob` so the assertions hold, then re-run.

- [ ] **Step 3: Document the opposite-side fallback** (`kalshi_mlb_rfq/correlation.py`, `joint_prob` docstring)

Add to the docstring:
```python
    """P(A ∩ B) from the grid, or None if not resolvable (caller → ρ=1).

    Same spread_side + same total_side → the joint is the tighter quadrant,
    one grid cell, read at the signed `spread_line` (negative = home-favorite
    grid, positive = away-favorite grid).

    Opposite spread_side (e.g. a home-margin combo vs an away-margin combo,
    which live in opposite signed grids and are near-disjoint) returns None.
    The caller's ρ=1 fallback then treats them as max positively correlated —
    conservative: it down-sizes rather than crediting a diversification
    benefit we cannot verify from the grid. (v1 limitation, intentional.)
    """
```

- [ ] **Step 4: Run full correlation suite**

Run: `cd /Users/callancapitolo/NFLWork/.claude/worktrees/worktree-rfq-remove-model && python3 -m pytest kalshi_mlb_rfq/tests/test_correlation.py kalshi_mlb_rfq/tests/test_correlation_integration.py -q`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/correlation.py kalshi_mlb_rfq/tests/test_correlation.py
git commit -m "test(rfq): pin correlation under signed grids; document opposite-side ρ=1 fallback"
```

---

### Task 4: Both-grids integration regression + FD symmetry check

**Files:**
- Test: `kalshi_mlb_rfq/tests/test_book_only_pricing.py`
- Read-through (fix only if asymmetric): `mlb_sgp/scraper_fanduel_sgp.py` — `price_combo` (441) home_line/away_line derivation

**Interfaces:**
- Consumes: `_SGP_ODDS_CACHE` (module global pandas frame), `_load_book_fairs`, `_combo_region_from_legs`.

- [ ] **Step 1: Write the failing integration test** (append to `kalshi_mlb_rfq/tests/test_book_only_pricing.py`)

This seeds BOTH grids with asymmetric odds and asserts an away-margin leg prices to the +grid away cell, NOT the −grid complement. Use the existing seed helper pattern in this file (build a DataFrame with columns `game_id, bookmaker, spread_line, total_line, combo, decimal_odds` for all four labels per grid; seed ≥2 books; patch `config.MIN_BOOK_COUNT_FOR_BLEND=2`).

```python
def test_away_margin_prices_to_positive_grid_not_complement(monkeypatch):
    import pandas as pd, statistics
    import kalshi_mlb_rfq.config as cfg
    monkeypatch.setattr(cfg, "MIN_BOOK_COUNT_FOR_BLEND", 2)

    # Four cells per grid per book. Negative grid (home favorite) and the
    # positive grid (away favorite) carry DELIBERATELY different odds so a
    # mis-route is detectable.
    def grid_rows(book, game, spread_line, cells):
        # cells: {(combo_label): decimal_odds}
        return [
            {"game_id": game, "bookmaker": book, "spread_line": spread_line,
             "total_line": 7.5, "combo": label, "decimal_odds": dec}
            for label, dec in cells.items()
        ]

    NEG = {  # home favorite grid: "Away Spread + Over" cell is the HIGH-prob complement
        "Home Spread + Over": 5.0, "Home Spread + Under": 5.0,
        "Away Spread + Over": 1.6, "Away Spread + Under": 1.6,
    }
    POS = {  # away favorite grid: "Away Spread + Over" cell is the LOW-prob away -1.5
        "Home Spread + Over": 1.6, "Home Spread + Under": 1.6,
        "Away Spread + Over": 5.0, "Away Spread + Under": 5.0,
    }
    rows = []
    for book in ("draftkings", "fanduel"):
        rows += grid_rows(book, "G", -1.5, NEG)
        rows += grid_rows(book, "G", 1.5, POS)
    monkeypatch.setattr(m, "_SGP_ODDS_CACHE", pd.DataFrame(rows))

    # Away-margin YES leg → region(away, +1.5, over) → POS grid "Away Spread + Over"
    region = _region(team_is_home=False, side="yes", total_side="yes")
    assert region.spread_line == 1.5
    books = m._load_book_fairs("G", region.spread_line, region.total_line,
                               region.spread_side, region.total_side)
    # POS "Away Spread + Over" is the ~5.0 long-shot cell, devigged ~0.18 —
    # NOT the ~0.55 complement it would hit on the −1.5 grid.
    fair = statistics.median(books.values())
    assert fair < 0.35, f"away-margin mis-routed to complement: {fair}"
```

(If the file's existing seed helper differs, reuse it instead of `grid_rows` — match the cache schema the other tests in this file use. The assertion is the contract: away-margin fair < 0.35 here.)

- [ ] **Step 2: Run test to verify it fails on the OLD sign (sanity)**

Run: `cd /Users/callancapitolo/NFLWork/.claude/worktrees/worktree-rfq-remove-model && python3 -m pytest kalshi_mlb_rfq/tests/test_book_only_pricing.py -q -k away_margin_prices`
Expected: PASS now (Task 2 already signs the line). To prove the test has teeth, temporarily revert the Task-2 sign (make it always negative), confirm FAIL, then restore. Note this in the commit.

- [ ] **Step 3: FD `price_combo` symmetry read-through**

Read `mlb_sgp/scraper_fanduel_sgp.py::price_combo` (line 441) and confirm it derives `away_line = -home_line` (so a `+1.5` target fetches FD's away −1.5 / home +1.5 runners, keyed `("away", -1.5)` / `("home", 1.5)`). If it instead hardcodes `away_line = abs(...)` or assumes a sign, fix it to `away_line = -home_line`. Record the finding (symmetric ✓ / fixed) in the commit message. No behavior change if already symmetric.

- [ ] **Step 4: Run the full book-only + pricing suite**

Run: `cd /Users/callancapitolo/NFLWork/.claude/worktrees/worktree-rfq-remove-model && python3 -m pytest kalshi_mlb_rfq/tests/ -q`
Expected: PASS (entire taker suite).

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/tests/test_book_only_pricing.py mlb_sgp/scraper_fanduel_sgp.py
git commit -m "test(rfq): both-grids regression — away-margin prices to +grid; verify FD away_line symmetry"
```

---

### Task 5: Per-book coverage report (go-live gate)

**Files:**
- Create: `kalshi_mlb_rfq/tools/coverage_report.py`

**Interfaces:**
- Standalone script. Reads `kalshi_mlb_rfq/kalshi_mlb_rfq_market.duckdb::mlb_sgp_odds` (whatever the last live `both_teams=True` scrape wrote) and the live Kalshi spread/total ladder. Prints a per-book table; writes nothing persistent. This is a verification tool, not a unit-tested module.

- [ ] **Step 1: Write the script**

```python
"""One-shot go-live gate: per book, what fraction of the Kalshi spread×total
ladder (BOTH signed grids, all N) comes back as a full 4-cell grid?

Run AFTER the bot has done at least one both_teams=True SGP scrape today.
Reads the bot market DB read-only; prints a table; persists nothing.
"""
import duckdb
from collections import defaultdict
from pathlib import Path

MARKET_DB = Path(__file__).resolve().parents[1] / "kalshi_mlb_rfq_market.duckdb"
LABELS = ["Home Spread + Over", "Home Spread + Under",
          "Away Spread + Over", "Away Spread + Under"]

def main():
    con = duckdb.connect(str(MARKET_DB), read_only=True)
    try:
        df = con.execute(
            "SELECT game_id, bookmaker, spread_line, total_line, combo "
            "FROM mlb_sgp_odds"
        ).fetchdf()
    finally:
        con.close()
    if df.empty:
        print("mlb_sgp_odds empty — run a both_teams=True SGP scrape first.")
        return

    # A (game, book, spread_line, total_line) cell-group is 'full' iff all 4
    # combo labels are present.
    full = defaultdict(int); seen = defaultdict(int)
    neg_full = defaultdict(int); pos_full = defaultdict(int)
    neg_seen = defaultdict(int); pos_seen = defaultdict(int)
    grp = df.groupby(["game_id", "bookmaker", "spread_line", "total_line"])
    for (game, book, sl, tl), sub in grp:
        labels = set(sub["combo"].unique())
        is_full = all(l in labels for l in LABELS)
        seen[book] += 1; full[book] += int(is_full)
        if sl < 0:
            neg_seen[book] += 1; neg_full[book] += int(is_full)
        else:
            pos_seen[book] += 1; pos_full[book] += int(is_full)

    print(f"{'book':12} {'all grids':>14} {'home(-) grids':>16} {'away(+) grids':>16}")
    for book in sorted(seen):
        def pct(f, s): return f"{f}/{s} ({100*f/s:.0f}%)" if s else "0/0"
        print(f"{book:12} {pct(full[book],seen[book]):>14} "
              f"{pct(neg_full[book],neg_seen[book]):>16} "
              f"{pct(pos_full[book],pos_seen[book]):>16}")
    print("\nThe away(+) column is the new capability — if it is ~0% for a "
          "book, that book has no dog-side alt SGP coverage and its away-margin "
          "combos will drop (by design).")

if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run it (after a live both_teams scrape) and eyeball**

Run: `cd /Users/callancapitolo/NFLWork && python3 kalshi_mlb_rfq/tools/coverage_report.py`
Expected: a per-book table. **This output is the go-live gate** — present the away(+) coverage numbers to the user before any live trading decision. (If `mlb_sgp_odds` is stale/empty, the script says so.)

- [ ] **Step 3: Commit the tool**

```bash
git add kalshi_mlb_rfq/tools/coverage_report.py
git commit -m "feat(rfq): per-book both-grid coverage report (go-live gate)"
```

---

### Task 6: Docs

**Files:**
- Modify: `kalshi_mlb_rfq/README.md`
- Modify: `CLAUDE.md` (root) — the taker bullet under Project Structure

- [ ] **Step 1: Update `kalshi_mlb_rfq/README.md`**

Add to the pricing/correlation section: book-only mode now prices BOTH teams' Kalshi margin markets. The cache `spread_line` is home-perspective signed — negative = home-favorite grid (home −1.5 etc.), positive = away-favorite grid (away −1.5 etc.). Enumeration emits both via `sgp_cycle(..., both_teams=True)` (the MM bot keeps `both_teams=False`). Away-margin coverage depends on book dog-side alt SGP availability; thin grids drop via the ≥2-book / 4-cell floor. Run `tools/coverage_report.py` to measure per-book coverage.

- [ ] **Step 2: Update root `CLAUDE.md` taker bullet**

In the "Autonomous Kalshi MLB SGP taker bot" bullet, append: book-only pricing now covers both teams' margin markets via signed `spread_line` grids (negative = home-favorite, positive = away-favorite); enumeration opt-in via `both_teams=True` keeps the shared MM path unchanged.

- [ ] **Step 3: Commit**

```bash
git add kalshi_mlb_rfq/README.md CLAUDE.md
git commit -m "docs(rfq): document signed-line both-team margin pricing + coverage gate"
```

---

## Pre-merge (after all tasks, before requesting merge approval)

- [ ] Run full taker suite: `cd /Users/callancapitolo/NFLWork/.claude/worktrees/worktree-rfq-remove-model && python3 -m pytest kalshi_mlb_rfq/tests/ -q` — all green.
- [ ] Run MM suite to prove no regression: `python3 -m pytest kalshi_mlb_mm/tests/ -q` (if present) — all green.
- [ ] Executive-engineer review of `git diff main..HEAD`.
- [ ] Present coverage-report numbers + diff summary; **get explicit user approval before merge.**
- [ ] After merge: remove worktree + branch; restart the bot from the **main repo cwd** (not the worktree), per restart-gotchas.

## Self-review notes (spec coverage)

- Spec "enumeration signed lines" → Task 1. "leg-routing sign fix" → Task 2. "correlation grid_lookup" → Task 3. "FD price_combo symmetric + live coverage" → Task 4 + Task 5. "coverage report go-live gate" → Task 5. "tests: 4-row routing + asymmetric both-grids integration" → Task 2 + Task 4. Docs → Task 6.
- MM blast-radius (shared `sgp_runner`, `_spread_line_from_legs`) handled by default-off flag + RFQ routing off `region.spread_line` (Global Constraints + Task 1/2).
- Open spec check "suppress redundant same-event candidates (same-N collision)" is deferred (documented here, not implemented) — the conservative ρ=1 fallback + ≥2-book floor make it a sizing-efficiency issue, not a correctness one. Revisit post-go-live if coverage data shows duplicate same-event RFQs.
