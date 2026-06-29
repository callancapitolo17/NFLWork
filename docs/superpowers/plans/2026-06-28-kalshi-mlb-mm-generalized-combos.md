# Kalshi MLB MM — Generalized N-Leg Combo Pricing (Phase 1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Install the leg-set router foundation in the Kalshi MLB maker — a canonical/hashed leg representation, per-game partitioning, and cross-game independent-multiply pricing — while preserving today's 2-leg behavior exactly and keeping every novel same-game shape fail-safe-skipped (Phase 2 prices those on-demand).

**Architecture:** A new pure-function module `kalshi_common/legset.py` normalizes each RFQ's `mve_selected_legs` into typed `CanonicalLeg`s, hashes the sorted set, and partitions legs by game. A new `kalshi_mlb_mm/router.py` prices a whole combo by pricing each game's sub-combo (a cached-grid 2-leg cell or a single derived by marginalizing the grid), taking a per-game book consensus, and multiplying across games. `main.py` routes through these instead of the fixed `combo_descriptor` + `_book_fairs` + `blended_fair` path. No schema change, no new network behavior — Phase 1 reads only the existing 60s `_SGP_ODDS` cache.

**Tech Stack:** Python 3, pandas, scipy (`devig_book` probit), DuckDB (read-only cache), pytest.

## Global Constraints

- **Spec:** `docs/superpowers/specs/2026-06-28-kalshi-mlb-mm-generalized-combos-design.md` (read it first).
- **Fail-safe contract:** any untypeable leg, any unresolved game, or any game sub-combo that fails consensus → the whole combo returns `None` → no quote. Never quote a combo we cannot price soundly.
- **Behavior preservation:** existing 2-leg `spread×total` and `ml×total` RFQs must price to byte-identical fairs after the refactor (regression-locked in Task 7).
- **No model:** fair value is pure book consensus (`MIN_AGREEING_BOOKS=2`, `BOOK_CONSENSUS_BAND=0.02`). The internal model is Phase 3.
- **Devig is probit:** all devigging goes through `kalshi_common.fair_value.devig_book` (n-way probit). Do not reintroduce multiplicative devig.
- **Pure functions stay pure:** `legset.py` and `router.py` take their inputs as arguments (DataFrames, ints, floats, callables) — they do NOT import `kalshi_mlb_mm.config`. `main.py` passes config values in. This keeps both modules unit-testable without env/DB.
- **Run tests with the maker venv from the repo root:** `./kalshi_mlb_mm/venv/bin/python -m pytest <path> -v`.
- **Work in the worktree** `worktree-kalshi-mm-generalized-combos`; do not merge without explicit user approval.

---

## File Structure

**Create:**
- `kalshi_common/legset.py` — the normalizer (pure): `CanonicalLeg`, `parse_leg`, `parse_legs`, `canonical_legs`, `leg_set_hash`, `partition_by_game`, `classify_subcombo`.
- `kalshi_common/tests/test_legset.py` — unit tests for the normalizer.
- `kalshi_mlb_mm/router.py` — combo pricing (pure): `consensus`, `grid_cell_fairs`, `single_marginal_fairs`, `subcombo_fair`, `combo_fair`.
- `kalshi_mlb_mm/tests/test_router.py` — unit tests for the router.

**Modify:**
- `kalshi_mlb_mm/main.py` — route discovery scope + pricing + circuit-breaker + confirm last-look + risk-sweep through `router.combo_fair`; resolve each partitioned game to its `game_id`.

**Unchanged (relied upon):**
- `kalshi_common/leg_types.py` — `_leg_dict_to_typed`, `_moneyline_side`, `_event_codes_from_legs`, `SPREAD_TOTAL_FAMILY`, `ML_TOTAL_FAMILY`, `_MLB_CODE_TO_TEAM`.
- `kalshi_common/fair_value.py` — `devig_book(book_rows, combo, vig_fallback=0.0)`.
- `kalshi_mlb_mm/config.py` — `MIN_AGREEING_BOOKS`, `BOOK_CONSENSUS_BAND`, `MARKET_DB`, etc.

---

## Task 1: CanonicalLeg + leg parsing

**Files:**
- Create: `kalshi_common/legset.py`
- Test: `kalshi_common/tests/test_legset.py`

**Interfaces:**
- Consumes: `kalshi_common.leg_types._leg_dict_to_typed`, `._moneyline_side`.
- Produces:
  - `CanonicalLeg(game_id: str, market_type: str, line: float | None, side: str)` — frozen dataclass. `market_type ∈ {"spread","total","ml"}`; `line` is signed home-perspective for spread/total, `None` for ml; `side ∈ {"home","away"}` for spread/ml, `{"over","under"}` for total.
  - `parse_leg(leg: dict) -> CanonicalLeg | None` — one Kalshi leg dict (`{market_ticker, event_ticker, side}`) → CanonicalLeg, or None if untypeable.
  - `parse_legs(legs: list[dict]) -> list[CanonicalLeg] | None` — all legs, or None if ANY leg is untypeable (fail-safe).

- [ ] **Step 1: Write the failing test**

```python
# kalshi_common/tests/test_legset.py
from kalshi_common import legset

EVT = "KXMLBGAME-25JUN271905NYYBOS"  # away=NYY, home=BOS

def _spread(team, n, side):   # team code, ticker N (2 => ±1.5), yes/no
    return {"market_ticker": f"KXMLBSPREAD-25JUN271905NYYBOS-{team}{n}",
            "event_ticker": EVT, "side": side}

def _total(n, side):
    return {"market_ticker": f"KXMLBTOTAL-25JUN271905NYYBOS-{n}",
            "event_ticker": EVT, "side": side}

def _ml(team, side):
    return {"market_ticker": f"KXMLBGAME-25JUN271905NYYBOS-{team}",
            "event_ticker": EVT, "side": side}

def test_parse_spread_home_minus_1_5():
    leg = legset.parse_leg(_spread("BOS", 2, "yes"))   # home -1.5, YES
    assert leg == legset.CanonicalLeg(EVT, "spread", -1.5, "home")

def test_parse_total_over_8_5():
    leg = legset.parse_leg(_total(9, "yes"))            # Over 8.5
    assert leg == legset.CanonicalLeg(EVT, "total", 8.5, "over")

def test_parse_ml_home():
    leg = legset.parse_leg(_ml("BOS", "yes"))           # home ML
    assert leg == legset.CanonicalLeg(EVT, "ml", None, "home")

def test_parse_spread_no_side_flips_to_away():
    # home -1.5 NO == away covers +1.5
    leg = legset.parse_leg(_spread("BOS", 2, "no"))
    assert leg == legset.CanonicalLeg(EVT, "spread", -1.5, "away")

def test_unparseable_leg_returns_none():
    bad = {"market_ticker": "KXMLBPLAYER-foo", "event_ticker": EVT, "side": "yes"}
    assert legset.parse_leg(bad) is None

def test_parse_legs_all_or_nothing():
    good = [_spread("BOS", 2, "yes"), _total(9, "yes")]
    assert legset.parse_legs(good) is not None
    assert len(legset.parse_legs(good)) == 2
    bad = good + [{"market_ticker": "KXMLBPLAYER-foo", "event_ticker": EVT, "side": "yes"}]
    assert legset.parse_legs(bad) is None
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./kalshi_mlb_mm/venv/bin/python -m pytest kalshi_common/tests/test_legset.py -v`
Expected: FAIL — `ModuleNotFoundError` / `AttributeError: module 'kalshi_common.legset' has no attribute ...`

- [ ] **Step 3: Write minimal implementation**

```python
# kalshi_common/legset.py
"""Canonical leg-set normalizer for arbitrary N-leg Kalshi MLB combos.

Pure functions: parse Kalshi leg dicts into typed CanonicalLegs, canonicalize +
hash a leg set, partition by game, and classify each game's sub-combo's pricing
route. Network-free and config-free so it is fully unit-testable.
"""
import hashlib
import json
from dataclasses import dataclass

from kalshi_common import leg_types


@dataclass(frozen=True)
class CanonicalLeg:
    game_id: str          # the leg's event_ticker (all legs of one game share it)
    market_type: str      # "spread" | "total" | "ml"
    line: float | None    # signed home-perspective; None for ml
    side: str             # "home"/"away" (spread, ml) or "over"/"under" (total)


def parse_leg(leg: dict) -> CanonicalLeg | None:
    """One Kalshi leg dict -> CanonicalLeg in home-perspective, or None."""
    et = str(leg.get("event_ticker", ""))
    mt = str(leg.get("market_ticker", ""))
    if not et or not mt:
        return None
    if mt.startswith("KXMLBSPREAD-"):
        typed = leg_types._leg_dict_to_typed(leg, "")   # SpreadLeg or None
        if typed is None:
            return None
        home_covers = ((typed.team_is_home and typed.side == "yes")
                       or (not typed.team_is_home and typed.side == "no"))
        return CanonicalLeg(et, "spread", -(typed.line_n - 0.5),
                            "home" if home_covers else "away")
    if mt.startswith("KXMLBTOTAL-"):
        typed = leg_types._leg_dict_to_typed(leg, "")   # TotalLeg or None
        if typed is None:
            return None
        return CanonicalLeg(et, "total", typed.line_n - 0.5,
                            "over" if typed.side == "yes" else "under")
    if mt.startswith("KXMLBGAME-"):
        ml = leg_types._moneyline_side(leg)             # (team_is_home, side) or None
        if ml is None:
            return None
        team_is_home, side = ml
        home_ml = ((team_is_home and side == "yes")
                   or (not team_is_home and side == "no"))
        return CanonicalLeg(et, "ml", None, "home" if home_ml else "away")
    return None


def parse_legs(legs: list[dict]) -> list[CanonicalLeg] | None:
    """All legs typed, or None if ANY leg is untypeable (fail-safe)."""
    if not legs:
        return None
    out = []
    for leg in legs:
        c = parse_leg(leg)
        if c is None:
            return None
        out.append(c)
    return out
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./kalshi_mlb_mm/venv/bin/python -m pytest kalshi_common/tests/test_legset.py -v`
Expected: PASS (6 passed)

- [ ] **Step 5: Commit**

```bash
git add kalshi_common/legset.py kalshi_common/tests/test_legset.py
git commit -m "feat(legset): CanonicalLeg + Kalshi leg parsing (home-perspective)"
```

---

## Task 2: Canonicalize + hash

**Files:**
- Modify: `kalshi_common/legset.py`
- Test: `kalshi_common/tests/test_legset.py`

**Interfaces:**
- Consumes: `CanonicalLeg` (Task 1).
- Produces:
  - `canonical_legs(legs: list[CanonicalLeg]) -> list[CanonicalLeg]` — stably sorted (None-safe on `line`).
  - `leg_set_hash(legs: list[CanonicalLeg]) -> str` — sha1 hex of the sorted set; identical regardless of input order.

- [ ] **Step 1: Write the failing test**

```python
# append to kalshi_common/tests/test_legset.py
def test_hash_is_order_independent():
    a = legset.parse_legs([_spread("BOS", 2, "yes"), _total(9, "yes")])
    b = legset.parse_legs([_total(9, "yes"), _spread("BOS", 2, "yes")])
    assert legset.leg_set_hash(a) == legset.leg_set_hash(b)

def test_hash_distinguishes_different_sets():
    a = legset.parse_legs([_spread("BOS", 2, "yes"), _total(9, "yes")])
    c = legset.parse_legs([_spread("BOS", 2, "yes"), _total(9, "no")])  # under
    assert legset.leg_set_hash(a) != legset.leg_set_hash(c)

def test_canonical_legs_sorts_with_none_line():
    legs = legset.parse_legs([_total(9, "yes"), _ml("BOS", "yes"),
                              _spread("BOS", 2, "yes")])
    ordered = legset.canonical_legs(legs)
    assert [l.market_type for l in ordered] == sorted(l.market_type for l in legs)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./kalshi_mlb_mm/venv/bin/python -m pytest kalshi_common/tests/test_legset.py -k "hash or canonical" -v`
Expected: FAIL — `AttributeError: ... has no attribute 'leg_set_hash'`

- [ ] **Step 3: Write minimal implementation**

```python
# append to kalshi_common/legset.py
def _sort_key(l: CanonicalLeg):
    # None-safe: ml legs (line=None) sort after numeric lines within a market_type
    return (l.game_id, l.market_type, l.line is None, l.line or 0.0, l.side)


def canonical_legs(legs: list[CanonicalLeg]) -> list[CanonicalLeg]:
    return sorted(legs, key=_sort_key)


def leg_set_hash(legs: list[CanonicalLeg]) -> str:
    payload = [[l.game_id, l.market_type, l.line, l.side]
               for l in canonical_legs(legs)]
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha1(blob.encode()).hexdigest()
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./kalshi_mlb_mm/venv/bin/python -m pytest kalshi_common/tests/test_legset.py -v`
Expected: PASS (9 passed)

- [ ] **Step 5: Commit**

```bash
git add kalshi_common/legset.py kalshi_common/tests/test_legset.py
git commit -m "feat(legset): order-independent canonical sort + sha1 leg_set_hash"
```

---

## Task 3: Partition by game + classify route

**Files:**
- Modify: `kalshi_common/legset.py`
- Test: `kalshi_common/tests/test_legset.py`

**Interfaces:**
- Consumes: `CanonicalLeg` (Task 1).
- Produces:
  - `partition_by_game(legs: list[CanonicalLeg]) -> dict[str, list[CanonicalLeg]]` — keyed by `game_id` (event_ticker).
  - `classify_subcombo(game_legs: list[CanonicalLeg]) -> str` — one of `"single"`, `"grid_spread_total"`, `"grid_ml_total"`, `"on_demand"`, `"unpriceable"`. Phase 1 prices `single`/`grid_*`; `on_demand` is Phase 2 (Phase 1 callers treat it as not-yet-priceable).

- [ ] **Step 1: Write the failing test**

```python
# append to kalshi_common/tests/test_legset.py
EVT2 = "KXMLBGAME-25JUN271905LADSFG"  # a different game

def _spread2(team, n, side):
    return {"market_ticker": f"KXMLBSPREAD-25JUN271905LADSFG-{team}{n}",
            "event_ticker": EVT2, "side": side}

def test_partition_groups_two_games():
    legs = legset.parse_legs([_spread("BOS", 2, "yes"), _spread2("SF", 2, "yes")])
    parts = legset.partition_by_game(legs)
    assert set(parts) == {EVT, EVT2}
    assert len(parts[EVT]) == 1 and len(parts[EVT2]) == 1

def test_classify_single():
    legs = legset.parse_legs([_ml("BOS", "yes")])
    assert legset.classify_subcombo(legs) == "single"

def test_classify_grid_spread_total():
    legs = legset.parse_legs([_spread("BOS", 2, "yes"), _total(9, "yes")])
    assert legset.classify_subcombo(legs) == "grid_spread_total"

def test_classify_grid_ml_total():
    legs = legset.parse_legs([_ml("BOS", "yes"), _total(9, "yes")])
    assert legset.classify_subcombo(legs) == "grid_ml_total"

def test_classify_three_leg_is_on_demand():
    legs = legset.parse_legs([_spread("BOS", 2, "yes"), _total(9, "yes"),
                              _ml("BOS", "yes")])
    assert legset.classify_subcombo(legs) == "on_demand"

def test_classify_empty_is_unpriceable():
    assert legset.classify_subcombo([]) == "unpriceable"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./kalshi_mlb_mm/venv/bin/python -m pytest kalshi_common/tests/test_legset.py -k "partition or classify" -v`
Expected: FAIL — `AttributeError: ... has no attribute 'partition_by_game'`

- [ ] **Step 3: Write minimal implementation**

```python
# append to kalshi_common/legset.py
def partition_by_game(legs: list[CanonicalLeg]) -> dict[str, list[CanonicalLeg]]:
    out: dict[str, list[CanonicalLeg]] = {}
    for l in legs:
        out.setdefault(l.game_id, []).append(l)
    return out


def classify_subcombo(game_legs: list[CanonicalLeg]) -> str:
    n = len(game_legs)
    if n == 0:
        return "unpriceable"
    if n == 1:
        return "single"
    types = sorted(l.market_type for l in game_legs)
    if n == 2 and types == ["spread", "total"]:
        return "grid_spread_total"
    if n == 2 and types == ["ml", "total"]:
        return "grid_ml_total"
    return "on_demand"
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./kalshi_mlb_mm/venv/bin/python -m pytest kalshi_common/tests/test_legset.py -v`
Expected: PASS (15 passed)

- [ ] **Step 5: Commit**

```bash
git add kalshi_common/legset.py kalshi_common/tests/test_legset.py
git commit -m "feat(legset): partition_by_game + classify_subcombo route"
```

---

## Task 4: Router — consensus + grid-cell book fairs

**Files:**
- Create: `kalshi_mlb_mm/router.py`
- Test: `kalshi_mlb_mm/tests/test_router.py`

**Interfaces:**
- Consumes: `kalshi_common.fair_value.devig_book`, `kalshi_common.leg_types.SPREAD_TOTAL_FAMILY` / `ML_TOTAL_FAMILY`, `kalshi_common.legset.CanonicalLeg`.
- Produces:
  - `consensus(book_fairs: dict[str, float], min_books: int, band: float) -> float | None` — `None` if `< min_books`, else median of books within `±band` of the median (re-checks `min_books`). Mirrors `main._consensus_filter` + `fairs.blended_fair` (median-of-agreeing).
  - `grid_spec(game_legs: list[CanonicalLeg]) -> tuple[tuple, float | None, float, str]` — `(family, spread_line, total_line, target_cell)` for a 2-leg grid sub-combo.
  - `grid_cell_fairs(game_id, family, spread_line, total_line, target_cell, sgp_df) -> dict[str, float]` — per-book devigged fair for the target cell (requires the full 4-cell grid per book; no fallback).

- [ ] **Step 1: Write the failing test**

```python
# kalshi_mlb_mm/tests/test_router.py
import pandas as pd
from kalshi_common import legset
from kalshi_mlb_mm import router

EVT = "KXMLBGAME-25JUN271905NYYBOS"

def _grid_rows(book, game_id, spread_line, total_line, family_decimals):
    # family_decimals: dict of combo-name -> decimal odds (4 cells)
    return [
        {"game_id": game_id, "combo": c, "period": "FG", "bookmaker": book,
         "sgp_decimal": d, "fetch_time": None,
         "spread_line": spread_line, "total_line": total_line}
        for c, d in family_decimals.items()
    ]

# A balanced 4-cell grid (~4.0 decimal each => ~0.25 raw, sums ~1.0 -> tiny vig)
ST_CELLS = {"Home Spread + Over": 4.2, "Home Spread + Under": 4.2,
            "Away Spread + Over": 3.8, "Away Spread + Under": 3.8}

def test_consensus_returns_none_below_min():
    assert router.consensus({"dk": 0.30}, min_books=2, band=0.02) is None

def test_consensus_median_of_agreeing():
    bf = {"dk": 0.30, "fd": 0.31, "px": 0.50}  # px is an outlier
    out = router.consensus(bf, min_books=2, band=0.02)
    assert out == 0.305  # median of dk+fd (px rejected)

def test_grid_spec_spread_total():
    legs = legset.parse_legs([
        {"market_ticker": "KXMLBSPREAD-25JUN271905NYYBOS-BOS2",
         "event_ticker": EVT, "side": "yes"},
        {"market_ticker": "KXMLBTOTAL-25JUN271905NYYBOS-9",
         "event_ticker": EVT, "side": "yes"}])
    family, spread_line, total_line, target = router.grid_spec(legs)
    assert spread_line == -1.5 and total_line == 8.5
    assert target == "Home Spread + Over"

def test_grid_cell_fairs_devigs_per_book():
    rows = (_grid_rows("dk", EVT, -1.5, 8.5, ST_CELLS)
            + _grid_rows("fd", EVT, -1.5, 8.5, ST_CELLS))
    df = pd.DataFrame(rows)
    from kalshi_common.leg_types import SPREAD_TOTAL_FAMILY
    out = router.grid_cell_fairs(EVT, SPREAD_TOTAL_FAMILY, -1.5, 8.5,
                                 "Home Spread + Over", df)
    assert set(out) == {"dk", "fd"}
    assert 0.20 < out["dk"] < 0.30   # devigged ~0.25-ish
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./kalshi_mlb_mm/venv/bin/python -m pytest kalshi_mlb_mm/tests/test_router.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'kalshi_mlb_mm.router'`

- [ ] **Step 3: Write minimal implementation**

```python
# kalshi_mlb_mm/router.py
"""Combo pricing router for arbitrary N-leg combos (Phase 1: cross-game multiply
+ cached-grid same-game). Pure functions — no config/DB imports; callers pass
the SGP odds DataFrame, the per-game resolver, and consensus params in.

Same-game shapes beyond the two 2-leg grids return None here (Phase 2 prices
them on-demand). Cross-game = product of per-game consensus fairs (independence).
"""
import statistics

from kalshi_common import legset
from kalshi_common.fair_value import devig_book
from kalshi_common.leg_types import SPREAD_TOTAL_FAMILY, ML_TOTAL_FAMILY


def consensus(book_fairs: dict[str, float], min_books: int,
              band: float) -> float | None:
    if len(book_fairs) < min_books:
        return None
    med = statistics.median(book_fairs.values())
    agree = {b: f for b, f in book_fairs.items() if abs(f - med) <= band}
    if len(agree) < min_books:
        return None
    return statistics.median(agree.values())


def grid_spec(game_legs: list[legset.CanonicalLeg]):
    """(family, spread_line, total_line, target_cell) for a 2-leg grid sub-combo."""
    total = next(l for l in game_legs if l.market_type == "total")
    ou = "Over" if total.side == "over" else "Under"
    spread = next((l for l in game_legs if l.market_type == "spread"), None)
    if spread is not None:
        part = "Home" if spread.side == "home" else "Away"
        return (SPREAD_TOTAL_FAMILY, spread.line, total.line,
                f"{part} Spread + {ou}")
    ml = next(l for l in game_legs if l.market_type == "ml")
    part = "Home" if ml.side == "home" else "Away"
    return (ML_TOTAL_FAMILY, None, total.line, f"{part} ML + {ou}")


def grid_cell_fairs(game_id, family, spread_line, total_line, target_cell,
                    sgp_df) -> dict[str, float]:
    if sgp_df is None or sgp_df.empty:
        return {}
    df = sgp_df
    mask = ((df.game_id == game_id) & df.combo.isin(family)
            & (df.total_line.astype(float).round(2) == round(total_line, 2)))
    if spread_line is None:
        mask &= df.spread_line.isna()
    else:
        mask &= (df.spread_line.astype(float).round(2) == round(spread_line, 2))
    rows = df[mask]
    out = {}
    for book in rows.bookmaker.unique():
        sub = rows[rows.bookmaker == book]
        if len(sub) < 4:                 # require the full 4-cell grid, no fallback
            continue
        f = devig_book(sub, combo=target_cell, vig_fallback=0.0)
        if f is not None:
            out[book] = f
    return out
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./kalshi_mlb_mm/venv/bin/python -m pytest kalshi_mlb_mm/tests/test_router.py -v`
Expected: PASS (4 passed)

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_mm/router.py kalshi_mlb_mm/tests/test_router.py
git commit -m "feat(router): consensus + grid_spec + grid_cell_fairs"
```

---

## Task 5: Router — single-leg fair via grid marginalization

**Files:**
- Modify: `kalshi_mlb_mm/router.py`
- Test: `kalshi_mlb_mm/tests/test_router.py`

**Interfaces:**
- Consumes: `devig_book`, `SPREAD_TOTAL_FAMILY`, `ML_TOTAL_FAMILY`, `CanonicalLeg`.
- Produces:
  - `single_marginal_fairs(game_id, leg: CanonicalLeg, sgp_df) -> dict[str, float]` — per-book fair probability of a single leg, derived by marginalizing the cached grid: sum the two devigged grid cells matching the leg's side over the free axis. (Spread/ML marginalize over total; total marginalizes over spread.)

**Why this works:** a book's 4-cell grid devigs to four fair probs summing to 1. The marginal `P(home covers spread)` = `devig(Home Spread+Over) + devig(Home Spread+Under)`; both come from the same grid group, so the sum is the true single-leg fair without any extra book call.

- [ ] **Step 1: Write the failing test**

```python
# append to kalshi_mlb_mm/tests/test_router.py
ML_CELLS = {"Home ML + Over": 3.6, "Home ML + Under": 4.0,
            "Away ML + Over": 4.4, "Away ML + Under": 4.0}

def test_single_marginal_total_over():
    # Over fair = Home Spread+Over + Away Spread+Over, devigged
    rows = (_grid_rows("dk", EVT, -1.5, 8.5, ST_CELLS)
            + _grid_rows("fd", EVT, -1.5, 8.5, ST_CELLS))
    df = pd.DataFrame(rows)
    leg = legset.CanonicalLeg(EVT, "total", 8.5, "over")
    out = router.single_marginal_fairs(EVT, leg, df)
    assert set(out) == {"dk", "fd"}
    assert 0.40 < out["dk"] < 0.60     # ~0.5 for this balanced grid

def test_single_marginal_ml_home():
    rows = _grid_rows("dk", EVT, None, 8.5, ML_CELLS)  # ml family => spread NULL
    df = pd.DataFrame(rows)
    leg = legset.CanonicalLeg(EVT, "ml", None, "home")
    out = router.single_marginal_fairs(EVT, leg, df)
    # Home ML = Home ML+Over + Home ML+Under, devigged; >0.5 (3.6 & 4.0 are short)
    assert out["dk"] > 0.45

def test_single_marginal_missing_grid_returns_empty():
    df = pd.DataFrame(_grid_rows("dk", EVT, -1.5, 8.5, ST_CELLS))
    leg = legset.CanonicalLeg(EVT, "spread", -2.5, "home")  # no -2.5 grid
    assert router.single_marginal_fairs(EVT, leg, df) == {}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./kalshi_mlb_mm/venv/bin/python -m pytest kalshi_mlb_mm/tests/test_router.py -k single_marginal -v`
Expected: FAIL — `AttributeError: ... has no attribute 'single_marginal_fairs'`

- [ ] **Step 3: Write minimal implementation**

```python
# append to kalshi_mlb_mm/router.py
def _single_marginal_spec(leg: legset.CanonicalLeg):
    """(family, fixed_axis, fixed_value, free_axis, cells) for marginalizing a
    single leg, or (None, ...) if not a supported market_type."""
    if leg.market_type == "spread":
        part = "Home" if leg.side == "home" else "Away"
        return (SPREAD_TOTAL_FAMILY, "spread_line", leg.line, "total_line",
                [f"{part} Spread + Over", f"{part} Spread + Under"])
    if leg.market_type == "total":
        ou = "Over" if leg.side == "over" else "Under"
        return (SPREAD_TOTAL_FAMILY, "total_line", leg.line, "spread_line",
                [f"Home Spread + {ou}", f"Away Spread + {ou}"])
    if leg.market_type == "ml":
        part = "Home" if leg.side == "home" else "Away"
        # ml family rows carry spread_line = NULL; marginalize over total
        return (ML_TOTAL_FAMILY, "spread_line", None, "total_line",
                [f"{part} ML + Over", f"{part} ML + Under"])
    return (None, None, None, None, None)


def _marginal_for_group(group_rows, cells) -> float | None:
    """One book's 4-cell grid group -> sum of the two devigged target cells."""
    if len(group_rows) < 4:
        return None
    total = 0.0
    for c in cells:
        f = devig_book(group_rows, combo=c, vig_fallback=0.0)
        if f is None:
            return None
        total += f
    return total


def single_marginal_fairs(game_id, leg: legset.CanonicalLeg, sgp_df) -> dict[str, float]:
    family, fixed_axis, fixed_value, free_axis, cells = _single_marginal_spec(leg)
    if family is None or sgp_df is None or sgp_df.empty:
        return {}
    df = sgp_df
    mask = (df.game_id == game_id) & df.combo.isin(family)
    if fixed_value is None:
        mask &= df[fixed_axis].isna()
    else:
        mask &= (df[fixed_axis].astype(float).round(2) == round(fixed_value, 2))
    rows = df[mask]
    out = {}
    for book in rows.bookmaker.unique():
        sub = rows[rows.bookmaker == book]
        # pick the first full 4-cell grid group along the free axis
        for _, grp in sub.groupby(free_axis, dropna=False):
            fair = _marginal_for_group(grp, cells)
            if fair is not None:
                out[book] = fair
                break
    return out
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./kalshi_mlb_mm/venv/bin/python -m pytest kalshi_mlb_mm/tests/test_router.py -v`
Expected: PASS (7 passed)

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_mm/router.py kalshi_mlb_mm/tests/test_router.py
git commit -m "feat(router): single-leg fair via grid marginalization"
```

---

## Task 6: Router — subcombo_fair + combo_fair (cross-game multiply)

**Files:**
- Modify: `kalshi_mlb_mm/router.py`
- Test: `kalshi_mlb_mm/tests/test_router.py`

**Interfaces:**
- Consumes: everything above; `legset.parse_legs`, `legset.partition_by_game`, `legset.classify_subcombo`.
- Produces:
  - `subcombo_fair(game_id, game_legs, sgp_df, min_books, band) -> float | None` — price one game's sub-combo: route → grid/single book fairs → `consensus`. `on_demand`/`unpriceable` → `None` (Phase 2 fills `on_demand`).
  - `combo_fair(legs: list[dict], sgp_df, resolve_game, min_books, band) -> float | None` — full RFQ: parse → partition by game → per-game `subcombo_fair` → multiply. `None` if any leg untypeable, any game unresolved (`resolve_game` returns `None`), or any sub-combo fails. `resolve_game(game_legs: list[CanonicalLeg]) -> str | None` is supplied by the caller.

- [ ] **Step 1: Write the failing test**

```python
# append to kalshi_mlb_mm/tests/test_router.py
def test_subcombo_grid_spread_total():
    legs = legset.parse_legs([
        {"market_ticker": "KXMLBSPREAD-25JUN271905NYYBOS-BOS2",
         "event_ticker": EVT, "side": "yes"},
        {"market_ticker": "KXMLBTOTAL-25JUN271905NYYBOS-9",
         "event_ticker": EVT, "side": "yes"}])
    df = pd.DataFrame(_grid_rows("dk", EVT, -1.5, 8.5, ST_CELLS)
                      + _grid_rows("fd", EVT, -1.5, 8.5, ST_CELLS))
    out = router.subcombo_fair(EVT, legs, df, min_books=2, band=0.02)
    assert out is not None and 0.20 < out < 0.30

def test_subcombo_on_demand_returns_none_in_phase1():
    legs = legset.parse_legs([
        {"market_ticker": "KXMLBSPREAD-25JUN271905NYYBOS-BOS2",
         "event_ticker": EVT, "side": "yes"},
        {"market_ticker": "KXMLBTOTAL-25JUN271905NYYBOS-9",
         "event_ticker": EVT, "side": "yes"},
        {"market_ticker": "KXMLBGAME-25JUN271905NYYBOS-BOS",
         "event_ticker": EVT, "side": "yes"}])
    df = pd.DataFrame(_grid_rows("dk", EVT, -1.5, 8.5, ST_CELLS))
    assert router.subcombo_fair(EVT, legs, df, 2, 0.02) is None

def test_combo_fair_cross_game_multiplies():
    EVT2 = "KXMLBGAME-25JUN271905LADSFG"
    legs_dicts = [
        {"market_ticker": "KXMLBGAME-25JUN271905NYYBOS-BOS",
         "event_ticker": EVT, "side": "yes"},                 # game1 home ML
        {"market_ticker": "KXMLBGAME-25JUN271905LADSFG-LAD",
         "event_ticker": EVT2, "side": "yes"}]                # game2 away ML
    df = pd.DataFrame(
        _grid_rows("dk", EVT, None, 8.5, ML_CELLS)
        + _grid_rows("fd", EVT, None, 8.5, ML_CELLS)
        + _grid_rows("dk", EVT2, None, 8.5, ML_CELLS)
        + _grid_rows("fd", EVT2, None, 8.5, ML_CELLS))
    resolve = lambda game_legs: game_legs[0].game_id  # identity for the test
    out = router.combo_fair(legs_dicts, df, resolve, 2, 0.02)
    # both games priced & multiplied -> strictly less than either single fair
    g1 = router.subcombo_fair(EVT, [l for l in legset.parse_legs(legs_dicts)
                                    if l.game_id == EVT], df, 2, 0.02)
    assert out is not None and abs(out - g1 * g1) < 1e-9  # symmetric grids

def test_combo_fair_skips_when_a_game_unresolved():
    legs_dicts = [{"market_ticker": "KXMLBGAME-25JUN271905NYYBOS-BOS",
                   "event_ticker": EVT, "side": "yes"}]
    df = pd.DataFrame(_grid_rows("dk", EVT, None, 8.5, ML_CELLS)
                      + _grid_rows("fd", EVT, None, 8.5, ML_CELLS))
    assert router.combo_fair(legs_dicts, df, lambda gl: None, 2, 0.02) is None

def test_combo_fair_skips_untypeable():
    legs_dicts = [{"market_ticker": "KXMLBPLAYER-foo",
                   "event_ticker": EVT, "side": "yes"}]
    assert router.combo_fair(legs_dicts, pd.DataFrame(), lambda gl: EVT, 2, 0.02) is None
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./kalshi_mlb_mm/venv/bin/python -m pytest kalshi_mlb_mm/tests/test_router.py -k "subcombo or combo_fair" -v`
Expected: FAIL — `AttributeError: ... has no attribute 'subcombo_fair'`

- [ ] **Step 3: Write minimal implementation**

```python
# append to kalshi_mlb_mm/router.py
def subcombo_fair(game_id, game_legs, sgp_df, min_books: int,
                  band: float) -> float | None:
    route = legset.classify_subcombo(game_legs)
    if route == "single":
        book_fairs = single_marginal_fairs(game_id, game_legs[0], sgp_df)
    elif route in ("grid_spread_total", "grid_ml_total"):
        family, spread_line, total_line, target = grid_spec(game_legs)
        book_fairs = grid_cell_fairs(game_id, family, spread_line, total_line,
                                     target, sgp_df)
    else:                       # "on_demand" (Phase 2) or "unpriceable"
        return None
    return consensus(book_fairs, min_books, band)


def combo_fair(legs: list[dict], sgp_df, resolve_game, min_books: int,
               band: float) -> float | None:
    canon = legset.parse_legs(legs)
    if canon is None:
        return None
    product = 1.0
    for _game_key, game_legs in legset.partition_by_game(canon).items():
        game_id = resolve_game(game_legs)
        if game_id is None:
            return None
        f = subcombo_fair(game_id, game_legs, sgp_df, min_books, band)
        if f is None:
            return None
        product *= f
    return product
```

- [ ] **Step 4: Run test to verify it passes**

Run: `./kalshi_mlb_mm/venv/bin/python -m pytest kalshi_mlb_mm/tests/test_router.py -v`
Expected: PASS (12 passed)

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_mm/router.py kalshi_mlb_mm/tests/test_router.py
git commit -m "feat(router): subcombo_fair routing + cross-game combo_fair multiply"
```

---

## Task 7: Integrate the router into main.py

**Files:**
- Modify: `kalshi_mlb_mm/main.py`
- Test: `kalshi_mlb_mm/tests/test_router_integration.py` (Create)

**Interfaces:**
- Consumes: `router.combo_fair`, `legset.parse_legs`, `legset.partition_by_game`, `legset.classify_subcombo`.
- Produces: a `_resolve_game_for_legs(game_legs)` closure in `main` that maps a single game's `CanonicalLeg`s → `mlb_target_lines` `game_id` (or `None`), and a `_priceable_in_phase1(canon)` scope helper.

**Context — what changes in `main.py`:**
- Today `main._discovery_tick` computes `desc = combo_descriptor(legs)`, `book_fairs = _book_fairs(game_id, desc)`, `book_med, blended = fairs.blended_fair(...)`. We replace the fair computation with one `router.combo_fair(...)` call that returns a single fair probability; `book_med` and `blended` become that same value (today they're already equal — both the median).
- The scope check (`in_scope = bool(legs and combo_descriptor(legs) is not None)`) becomes: typeable AND every game sub-combo is Phase-1-priceable (`single`/`grid_*`). `on_demand` shapes are out-of-scope in Phase 1 (skipped, logged `out_of_scope`), exactly as a 3-leg combo is skipped today.
- The confirm last-look (`_confirm_tick`) and risk-sweep (`_risk_sweep_tick`) recompute the fair the same way — swap their `_book_fairs`+`blended_fair` pair for `router.combo_fair`.
- Game resolution: today `_resolve_game(legs)` resolves ONE game. We add `_resolve_game_for_legs(game_legs)` that resolves one partitioned game's id; `router.combo_fair` calls it per game. The existing `_resolve_game` is kept for the single-game scope/tipoff path (the tipoff gate still uses the first/primary game's `commence_time`).

- [ ] **Step 1: Write the failing integration test**

```python
# kalshi_mlb_mm/tests/test_router_integration.py
"""Phase-1 regression: existing 2-leg combos price identically through the
router as through the legacy combo_descriptor + _book_fairs path."""
import pandas as pd
from kalshi_common import legset
from kalshi_common.leg_types import combo_descriptor
from kalshi_mlb_mm import router, main

EVT = "KXMLBGAME-25JUN271905NYYBOS"
ST_CELLS = {"Home Spread + Over": 4.2, "Home Spread + Under": 4.2,
            "Away Spread + Over": 3.8, "Away Spread + Under": 3.8}

def _grid_rows(book, game_id, spread_line, total_line, cells):
    return [{"game_id": game_id, "combo": c, "period": "FG", "bookmaker": book,
             "sgp_decimal": d, "fetch_time": None, "spread_line": spread_line,
             "total_line": total_line} for c, d in cells.items()]

LEGS = [{"market_ticker": "KXMLBSPREAD-25JUN271905NYYBOS-BOS2",
         "event_ticker": EVT, "side": "yes"},
        {"market_ticker": "KXMLBTOTAL-25JUN271905NYYBOS-9",
         "event_ticker": EVT, "side": "yes"}]

def test_router_matches_legacy_two_leg_fair():
    df = pd.DataFrame(_grid_rows("dk", EVT, -1.5, 8.5, ST_CELLS)
                      + _grid_rows("fd", EVT, -1.5, 8.5, ST_CELLS))
    # legacy path
    desc = combo_descriptor(LEGS)
    legacy = main._book_fairs.__wrapped__ if hasattr(main._book_fairs, "__wrapped__") else None
    # compute legacy via the module-global cache shim:
    main._SGP_ODDS = df
    legacy_fairs = main._book_fairs(EVT, desc)
    legacy_fair = __import__("statistics").median(legacy_fairs.values())
    # router path
    router_fair = router.combo_fair(LEGS, df, lambda gl: EVT, 2, 0.02)
    assert abs(router_fair - legacy_fair) < 1e-9
```

- [ ] **Step 2: Run test to verify it fails**

Run: `./kalshi_mlb_mm/venv/bin/python -m pytest kalshi_mlb_mm/tests/test_router_integration.py -v`
Expected: FAIL — initially because `main._SGP_ODDS` / `_book_fairs` interplay differs, OR it passes against the unchanged legacy path (which is the point: it locks the equivalence the refactor must preserve). If it fails on import or attribute, fix the test harness, not the assertion.

- [ ] **Step 3: Add the resolver + scope helpers to `main.py`**

Add near `_resolve_game` (around `main.py:166`):

```python
def _resolve_game_for_legs(game_legs):
    """Resolve one partitioned game's CanonicalLegs -> mlb_target_lines game_id.
    Returns None if teams unknown or not in mlb_target_lines yet (fail-safe)."""
    away, home = _parse_event_suffix(game_legs[0].game_id.rsplit("-", 1)[-1])
    home_name = _MLB_CODE_TO_TEAM.get(home)
    away_name = _MLB_CODE_TO_TEAM.get(away)
    if not home_name or not away_name or not config.MARKET_DB.exists():
        return None
    try:
        con = duckdb.connect(str(config.MARKET_DB), read_only=True)
    except duckdb.IOException:
        return None
    try:
        row = con.execute(
            "SELECT game_id FROM mlb_target_lines WHERE home_team=? AND away_team=? LIMIT 1",
            [home_name, away_name]).fetchone()
    except duckdb.CatalogException:
        row = None
    finally:
        con.close()
    return row[0] if row else None


def _priceable_in_phase1(canon) -> bool:
    """True iff all game sub-combos route to single/grid (Phase-1 priceable)."""
    from kalshi_common import legset
    parts = legset.partition_by_game(canon)
    return all(legset.classify_subcombo(gl) in
               ("single", "grid_spread_total", "grid_ml_total")
               for gl in parts.values())
```

Add the import at the top of `main.py` (with the other `kalshi_common` imports):

```python
from kalshi_common import legset
from kalshi_mlb_mm import router
```

- [ ] **Step 4: Replace the scope check in `_discovery_tick`**

Find (around `main.py:381`):

```python
            market = source.get_market(ticker)
            legs = scope.decode_legs(market) if market else None
            in_scope = bool(legs and combo_descriptor(legs) is not None)
            game_id = None
            _SCOPE_CACHE[ticker] = (in_scope, game_id, legs)
```

Replace with:

```python
            market = source.get_market(ticker)
            legs = scope.decode_legs(market) if market else None
            canon = legset.parse_legs(legs) if legs else None
            in_scope = bool(canon and _priceable_in_phase1(canon))
            game_id = None
            _SCOPE_CACHE[ticker] = (in_scope, game_id, legs)
```

- [ ] **Step 5: Replace the fair computation in `_discovery_tick`**

Find (around `main.py:461`):

```python
        book_fairs = _book_fairs(game_id, desc)
        book_med, blended = fairs.blended_fair(legs, game_id, book_fairs)
```

Replace with:

```python
        blended = router.combo_fair(legs, _SGP_ODDS, _resolve_game_for_legs,
                                    config.MIN_AGREEING_BOOKS,
                                    config.BOOK_CONSENSUS_BAND)
        book_med = blended   # single consensus fair; book_med == blended (as before)
```

And update the `research.emit("quote_priced", ...)` payload just below to drop the
per-book `book_fairs`/`agreeing_count`/`spread_line`/`total_line` fields (no longer
computed here) and instead log the combo shape:

```python
        research.emit("quote_priced", rfq_id=rid, ticker=ticker,
                      payload=dict(blended_fair=blended, yes_bid=q.yes_bid,
                                   no_bid=q.no_bid,
                                   leg_set_hash=legset.leg_set_hash(
                                       legset.parse_legs(legs)),
                                   n_games=len(legset.partition_by_game(
                                       legset.parse_legs(legs))),
                                   game_id=game_id))
```

Note: `desc` is still produced by `_resolve_game(legs)` earlier in the tick for the
tipoff/`game_id` path; leave that call in place — only the fair computation changes.

- [ ] **Step 6: Replace the fair recompute in `_confirm_tick` and `_risk_sweep_tick`**

In `_confirm_tick` (around `main.py:614`), find:

```python
            desc = combo_descriptor(legs)
            book_fairs_now = _book_fairs(game_id, desc)
            if not book_fairs_now:
                _log_decision("voided_no_fresh_books", ...)
                ...
                continue
            _, cur_fair = fairs.blended_fair(legs, game_id, book_fairs_now)
```

Replace with:

```python
            cur_fair = router.combo_fair(legs, _SGP_ODDS, _resolve_game_for_legs,
                                         config.MIN_AGREEING_BOOKS,
                                         config.BOOK_CONSENSUS_BAND)
            if cur_fair is None:
                _log_decision("voided_no_fresh_books", rfq_id=rid, quote_id=qid,
                              ticker=ticker, game_id=game_id)
                with db.connect() as con:
                    con.execute(
                        "UPDATE live_quotes SET status='voided', closed_at=? WHERE quote_id=?",
                        [datetime.now(timezone.utc), qid])
                continue
```

(Delete the now-dead `if cur_fair is None:` block that followed the old
`blended_fair` call, since the None case is handled above.)

In `_risk_sweep_tick` (around `main.py:837`), find:

```python
                    bf_now = _book_fairs(game_id, combo_descriptor(legs))
                    if bf_now:
                        cur_med = statistics.median(bf_now.values())
                        if risk.book_move_triggered(book_fair_at_q, cur_med,
                                                    config.BOOK_MOVE_CB_THRESHOLD):
                            cancel = True
```

Replace with:

```python
                    cur_med = router.combo_fair(
                        json.loads(legs_json), _SGP_ODDS, _resolve_game_for_legs,
                        config.MIN_AGREEING_BOOKS, config.BOOK_CONSENSUS_BAND)
                    if cur_med is not None and risk.book_move_triggered(
                            book_fair_at_q, cur_med, config.BOOK_MOVE_CB_THRESHOLD):
                        cancel = True
```

(`legs` in the risk sweep comes from `json.loads(legs_json)`; keep that decode.)

- [ ] **Step 7: Run the integration + unit tests**

Run: `./kalshi_mlb_mm/venv/bin/python -m pytest kalshi_mlb_mm/tests/ kalshi_common/tests/ -v`
Expected: PASS — the regression test confirms 2-leg fairs are unchanged; legset + router suites green.

- [ ] **Step 8: Smoke-test the dry-run wiring**

Run: `cd ~/NFLWork && ./kalshi_mlb_mm/venv/bin/python -m kalshi_mlb_mm.main --dry-run` (let it run one discovery cycle, Ctrl-C).
Expected: no exceptions; `quote_decisions` shows in-scope 2-leg RFQs priced (`dry_run_quote`) and 3-leg/novel RFQs logged `out_of_scope`. Cross-game RFQs (if any in the live slate) now price instead of skipping.

- [ ] **Step 9: Commit**

```bash
git add kalshi_mlb_mm/main.py kalshi_mlb_mm/tests/test_router_integration.py
git commit -m "feat(mm): route combo pricing through legset/router (cross-game multiply)

2-leg fairs unchanged (regression-locked); cross-game combos now price via
per-game consensus x multiply; novel same-game shapes fail-safe skip (Phase 2)."
```

---

## Self-Review (completed)

- **Spec coverage (Phase 1 scope):** normalizer §4.1/§4.2 → Tasks 1–2; partition + classify §3 → Task 3; per-game consensus §5.4 → Task 4; single/marginal pricing → Task 5; cross-game multiply §5.1 → Task 6; fail-safe contract §3 → Tasks 6–7; behavior preservation → Task 7 regression. **Deferred to Phase 2:** §4.3 `mlb_sgp_legset_odds` table, §5.2 sided-cross-product on-demand devig, §5.3 deadline/cache, §5.5 sanity guard, §6 client widening. **Deferred to Phase 3:** internal model.
- **Placeholder scan:** none — every code step has complete code.
- **Type consistency:** `CanonicalLeg(game_id, market_type, line, side)` used identically across Tasks 1–7; `consensus(book_fairs, min_books, band)`, `combo_fair(legs, sgp_df, resolve_game, min_books, band)` signatures consistent between definition (Tasks 4/6) and call sites (Task 7).

---

## Phase 2 — On-demand same-game pricing (separate plan, after Phase 1 review)

Detailed once Phase 1 lands (its tasks depend on Phase 1's realized `legset`/`router`
interfaces and on reading the six book clients). Scope:

- `kalshi_common/legset.py`: add `enumerate_corners(game_legs) -> list[list[CanonicalLeg]]` (flip each leg's side → `2^k` corners).
- `kalshi_common/sgp_oracle.py` (new): `price_legset(game_legs) -> dict[str, float]` — for each book, build each corner's leg-set, call the book's generic SGP endpoint concurrently under `LEGSET_PRICE_DEADLINE_SEC`, probit-devig the priced corners, read the target corner; apply the §5.5 sanity guard (book vs naive-multiply).
- `mlb_sgp/dk_client.py`, `fd_client.py`: widen the price wrappers from fixed 2-arg to a generic `legs: list` (bodies already POST a list; the other four clients already accept lists).
- `kalshi_mlb_mm/db.py`: add `mlb_sgp_legset_odds` (hash-keyed) + short-TTL cache helpers (`LEGSET_CACHE_TTL_SEC`).
- `kalshi_mlb_mm/router.py`: `subcombo_fair` `on_demand` branch calls `sgp_oracle.price_legset` (cache-first) → `consensus`.
- `kalshi_mlb_mm/config.py`: add `LEGSET_PRICE_DEADLINE_SEC` (3), `LEGSET_CACHE_TTL_SEC` (15).
- Instrument `created_ts → first-quote latency` + `quote_expired`/`rfq_closed` to measure the RFQ window.

## Phase 3 — Internal correlation model (separate plan, deferred)

Bivariate-Poisson SGP simulator implementing the `price_legset` signature as an
additional "book," slotted into `consensus`. Out of scope until Phase 2 is live and
measured. See `sgp_simulator_plan` memory.

---

## Version control / worktree / docs

- **Worktree:** `worktree-kalshi-mm-generalized-combos` (already created). Implement + test here; merge to `main` only on explicit user approval; clean up worktree + branch after merge.
- **DuckDB safety:** Phase 1 makes no schema change and only reads the existing market DB read-only. Never symlink DBs into the worktree.
- **Docs (in the Phase 1 merge commit):** `kalshi_mlb_mm/README.md` — note the router/legset seam and that pricing now partitions by game + multiplies cross-game (same-game novel shapes still skip until Phase 2). No `mlb_sgp` doc change in Phase 1 (no scraper change yet).
- **Pre-merge review:** executive-engineer review of the full diff (data integrity, resource safety — the new `_resolve_game_for_legs` opens/closes its own read-only connection; edge cases — empty `_SGP_ODDS`, unresolved game, off-season; dead code — confirm `_book_fairs`/`fairs.blended_fair` are either removed or intentionally retained; no secrets in logs) before requesting merge approval.
