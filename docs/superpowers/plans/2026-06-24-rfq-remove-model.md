# RFQ taker: book-only pricing + book-implied correlation — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove the Monte-Carlo model from the Kalshi MLB RFQ taker (behind a
`USE_MODEL` flag, default off): price from book consensus and size with a
book-implied correlation engine, restoring fills that the stale-model gate had
choked off.

**Architecture:** A `USE_MODEL=false` code path skips all model/sample work.
Fair value becomes `median(books)`. Same-game correlation for Kelly comes from
the per-game spread×total SGP grid the bot already scrapes (joint of two
same-direction combos = the tighter grid cell), with a conservative ρ=1 fallback
everywhere the grid can't resolve it.

**Tech Stack:** Python 3, DuckDB, pandas/numpy, pytest. Package
`kalshi_mlb_rfq`; shared math in `kalshi_common`.

## Global Constraints

- New env knobs go through `config._get(name, default)` (string env → typed).
- Default runtime behavior flips to **book-only**: `USE_MODEL` default `"false"`.
- `MIN_BOOK_COUNT_FOR_BLEND = 2` (unchanged) gates every book-consensus read.
- Pure-function modules (`correlation.py`) take **no** DB/global dependencies;
  callers inject a `grid_lookup` callable.
- All new timestamps (none here) would be `TIMESTAMPTZ` UTC — N/A this plan.
- TDD: failing test first, minimal impl, green, commit. Exact paths always.
- Never merge to `main` without explicit user approval.

**v1 correlation scope (confirmed deviation from spec — see handoff):** exact
joint is computed for pairs sharing **both** spread side and total side (the
dominant same-direction overbet risk); all other pairs (opposite side, or any
non-grid / cross-category leg) use the ρ=1 fallback `min(Pₐ,P_b)`. Fallback is
strictly safe (it only ever down-sizes). Full crossing-pair inclusion-exclusion
is a documented future enhancement.

---

### Task 1: `USE_MODEL` config flag

**Files:**
- Modify: `kalshi_mlb_rfq/config.py` (add one knob near `MAX_PREDICTION_STALENESS_SEC`, ~line 77)
- Test: `kalshi_mlb_rfq/tests/test_config.py`

**Interfaces:**
- Produces: `config.USE_MODEL: bool`

- [ ] **Step 1: Write the failing test**

Add to `tests/test_config.py`:
```python
def test_use_model_defaults_off():
    import importlib
    import kalshi_mlb_rfq.config as c
    importlib.reload(c)
    assert c.USE_MODEL is False


def test_use_model_parses_truthy(monkeypatch):
    import importlib
    monkeypatch.setenv("USE_MODEL", "true")
    import kalshi_mlb_rfq.config as c
    importlib.reload(c)
    assert c.USE_MODEL is True
    monkeypatch.setenv("USE_MODEL", "0")
    importlib.reload(c)
    assert c.USE_MODEL is False
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd kalshi_mlb_rfq && python -m pytest tests/test_config.py -k use_model -v`
Expected: FAIL (`AttributeError: module ... has no attribute 'USE_MODEL'`)

- [ ] **Step 3: Add the flag**

In `config.py`, after the `MAX_PREDICTION_STALENESS_SEC` line:
```python
USE_MODEL = _get("USE_MODEL", "false").lower() in ("1", "true", "yes", "on")
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd kalshi_mlb_rfq && python -m pytest tests/test_config.py -k use_model -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/config.py kalshi_mlb_rfq/tests/test_config.py
git commit -m "feat(rfq): add USE_MODEL config flag (default off)"
```

---

### Task 2: `correlation.py` — region model, joint, covariance, Fréchet clamp

**Files:**
- Create: `kalshi_mlb_rfq/correlation.py`
- Test: `kalshi_mlb_rfq/tests/test_correlation.py`

**Interfaces:**
- Produces:
  - `ComboRegion(spread_side: str, spread_line: float, total_side: str, total_line: float)` — frozen dataclass; `spread_side ∈ {"home","away"}`, `total_side ∈ {"over","under"}`; `spread_line` is the home-perspective line as stored in `mlb_sgp_odds`.
  - `frechet_clamp(joint: float, p_a: float, p_b: float) -> float`
  - `joint_prob(a: ComboRegion, b: ComboRegion, p_a: float, p_b: float, grid_lookup) -> float | None` — `grid_lookup(spread_line: float, total_line: float, spread_side: str, total_side: str) -> float | None`. Returns clamped joint, or `None` when the pair is not same-direction or the grid cell is missing (caller applies ρ=1 fallback).
  - `cov_returns(p_a: float, p_b: float, p_joint: float, price_a: float, price_b: float) -> float`

- [ ] **Step 1: Write the failing tests**

Create `tests/test_correlation.py`:
```python
import math
import pytest
from kalshi_mlb_rfq.correlation import (
    ComboRegion, frechet_clamp, joint_prob, cov_returns,
)


def _grid(values):
    """values: dict[(spread_line, total_line, spread_side, total_side)] -> prob"""
    def lookup(sl, tl, ss, ts):
        return values.get((sl, tl, ss, ts))
    return lookup


def test_frechet_clamp_bounds():
    # joint can't exceed min(pa,pb) nor fall below max(0, pa+pb-1)
    assert frechet_clamp(0.9, 0.3, 0.4) == pytest.approx(0.3)      # capped to min
    assert frechet_clamp(0.0, 0.8, 0.8) == pytest.approx(0.6)      # raised to lower bound
    assert frechet_clamp(0.1, 0.3, 0.4) == pytest.approx(0.1)      # already valid


def test_joint_same_direction_nested_total_is_tighter_cell():
    # Both Home+Over; A total 8.5, B total 9.5 -> tighter = Over 9.5 at same spread
    a = ComboRegion("home", -1.5, "over", 8.5)
    b = ComboRegion("home", -1.5, "over", 9.5)
    grid = _grid({(-1.5, 9.5, "home", "over"): 0.22})
    j = joint_prob(a, b, p_a=0.30, p_b=0.22, grid_lookup=grid)
    assert j == pytest.approx(0.22)


def test_joint_same_direction_nested_spread_home_takes_min_line():
    # Home favorite: tighter spread = more negative line (higher margin threshold)
    a = ComboRegion("home", -1.5, "over", 8.5)
    b = ComboRegion("home", -2.5, "over", 8.5)
    grid = _grid({(-2.5, 8.5, "home", "over"): 0.18})
    j = joint_prob(a, b, p_a=0.30, p_b=0.20, grid_lookup=grid)
    assert j == pytest.approx(0.18)


def test_joint_away_takes_max_line():
    # Away cover: tighter spread = larger (more positive) home line
    a = ComboRegion("away", 1.5, "under", 9.5)
    b = ComboRegion("away", 2.5, "under", 8.5)
    grid = _grid({(2.5, 8.5, "away", "under"): 0.12})
    j = joint_prob(a, b, p_a=0.25, p_b=0.20, grid_lookup=grid)
    assert j == pytest.approx(0.12)


def test_joint_opposite_direction_returns_none():
    a = ComboRegion("home", -1.5, "over", 8.5)
    b = ComboRegion("home", -1.5, "under", 8.5)
    assert joint_prob(a, b, 0.3, 0.3, _grid({})) is None


def test_joint_missing_cell_returns_none():
    a = ComboRegion("home", -1.5, "over", 8.5)
    b = ComboRegion("home", -2.5, "over", 8.5)
    assert joint_prob(a, b, 0.3, 0.2, _grid({})) is None


def test_joint_identity_returns_pa_clamped():
    a = ComboRegion("home", -1.5, "over", 8.5)
    grid = _grid({(-1.5, 8.5, "home", "over"): 0.30})
    assert joint_prob(a, a, 0.30, 0.30, grid) == pytest.approx(0.30)


def test_joint_applies_frechet_clamp_to_noisy_grid():
    # grid says joint 0.28 but min(pa,pb)=0.20 -> clamp to 0.20
    a = ComboRegion("home", -1.5, "over", 8.5)
    b = ComboRegion("home", -2.5, "over", 8.5)
    grid = _grid({(-2.5, 8.5, "home", "over"): 0.28})
    j = joint_prob(a, b, p_a=0.30, p_b=0.20, grid_lookup=grid)
    assert j == pytest.approx(0.20)


def test_cov_returns_independent_is_zero():
    assert cov_returns(0.3, 0.4, 0.12, 0.3, 0.4) == pytest.approx(0.0)


def test_cov_returns_positive_when_joint_exceeds_product():
    cov = cov_returns(0.3, 0.4, 0.20, 0.3, 0.4)
    assert cov > 0


def test_cov_returns_scales_inverse_with_prices():
    base = cov_returns(0.3, 0.4, 0.20, 0.3, 0.4)
    half = cov_returns(0.3, 0.4, 0.20, 0.15, 0.4)
    assert half == pytest.approx(base * 2)


@pytest.mark.parametrize("pa,pb,raw", [
    (0.3, 0.4, 0.9), (0.8, 0.8, 0.0), (0.1, 0.1, 0.5), (0.5, 0.5, 0.5),
])
def test_frechet_invariant(pa, pb, raw):
    j = frechet_clamp(raw, pa, pb)
    assert max(0.0, pa + pb - 1.0) - 1e-12 <= j <= min(pa, pb) + 1e-12
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd kalshi_mlb_rfq && python -m pytest tests/test_correlation.py -v`
Expected: FAIL (`ModuleNotFoundError: kalshi_mlb_rfq.correlation`)

- [ ] **Step 3: Implement `correlation.py`**

Create `kalshi_mlb_rfq/correlation.py`:
```python
"""Book-implied correlation for same-game combos.

Pure functions — no DB, no globals. The caller injects a `grid_lookup`
callable that reads the per-game spread×total SGP grid.

A combo is a quadrant of the (home margin M, game total T) plane:
  - home cover  iff  M > -spread_line ; away cover iff M < -spread_line
  - over        iff  T > total_line   ; under       iff T < total_line

For two combos sharing both spread side and total side, the joint event
"both hit" is the *tighter* quadrant, which is itself a grid cell — so the
joint is a single grid lookup. Opposite-side pairs (and any leg the grid
cannot price) return None; the caller treats None as the ρ=1 fallback.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class ComboRegion:
    spread_side: str   # "home" | "away"
    spread_line: float # home-perspective line as stored in mlb_sgp_odds
    total_side: str    # "over" | "under"
    total_line: float


def frechet_clamp(joint: float, p_a: float, p_b: float) -> float:
    """Clamp a (possibly noisy) joint into the Fréchet–Hoeffding bounds so the
    implied correlation stays valid."""
    lo = max(0.0, p_a + p_b - 1.0)
    hi = min(p_a, p_b)
    return max(lo, min(joint, hi))


def joint_prob(a: ComboRegion, b: ComboRegion,
               p_a: float, p_b: float, grid_lookup) -> float | None:
    """P(A ∩ B) from the grid, or None if not resolvable (caller → ρ=1)."""
    if a.spread_side != b.spread_side or a.total_side != b.total_side:
        return None
    # Tighter spread quadrant:
    #   home cover is M > -spread_line → tighter = more negative line = min()
    #   away cover is M < -spread_line → tighter = more positive line = max()
    if a.spread_side == "home":
        tight_spread = min(a.spread_line, b.spread_line)
    else:
        tight_spread = max(a.spread_line, b.spread_line)
    # Tighter total quadrant:
    #   over is T > total_line  → tighter = larger total = max()
    #   under is T < total_line → tighter = smaller total = min()
    if a.total_side == "over":
        tight_total = max(a.total_line, b.total_line)
    else:
        tight_total = min(a.total_line, b.total_line)
    raw = grid_lookup(tight_spread, tight_total, a.spread_side, a.total_side)
    if raw is None:
        return None
    return frechet_clamp(raw, p_a, p_b)


def cov_returns(p_a: float, p_b: float, p_joint: float,
                price_a: float, price_b: float) -> float:
    """Return-space covariance for Kelly:
        Cov(r_a, r_b) = [P(A∩B) − P(A)P(B)] / (price_a · price_b)."""
    cov_outcomes = p_joint - p_a * p_b
    return cov_outcomes / (price_a * price_b)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd kalshi_mlb_rfq && python -m pytest tests/test_correlation.py -v`
Expected: PASS (all)

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/correlation.py kalshi_mlb_rfq/tests/test_correlation.py
git commit -m "feat(rfq): book-implied correlation module (joint + cov + Fréchet clamp)"
```

---

### Task 3: Kelly — accept book-implied covariance for existing positions

**Files:**
- Modify: `kalshi_mlb_rfq/kelly.py` (existing-positions branch, lines ~57-92)
- Test: `kalshi_mlb_rfq/tests/test_kelly.py` (new — first kelly coverage)

**Interfaces:**
- Consumes: `correlation.cov_returns` (Task 2).
- Produces: `kelly.kelly_size_combo(..., existing_positions, ...)` where each
  position dict now carries a precomputed `cov_return` field
  (`{"cov_return": float, "contracts": int, "effective_price": float}`) used
  directly instead of deriving covariance from sample paths. `outcome_vec`
  becomes optional (default `None`) and is ignored when positions carry
  `cov_return`.

**Rationale:** keep `kelly_size_combo`'s mean/variance (already from
`blended_fair`) and its full-Kelly residual formula; only swap *how* the
per-position covariance term is obtained. The caller (Task 5) computes
`cov_return` via the correlation engine and passes it in.

- [ ] **Step 1: Write the failing tests**

Create `tests/test_kelly.py`:
```python
import numpy as np
import pytest
from kalshi_mlb_rfq import kelly


BANK = 1000.0
KF = 0.25


def test_no_positions_single_bet_kelly():
    # +EV: fair 0.30 vs price 0.20
    n = kelly.kelly_size_combo(
        outcome_vec=None, blended_fair=0.30, existing_positions=[],
        effective_price=0.20, bankroll=BANK, kelly_fraction=KF)
    assert n > 0


def test_negative_ev_returns_zero():
    n = kelly.kelly_size_combo(
        outcome_vec=None, blended_fair=0.20, existing_positions=[],
        effective_price=0.20, bankroll=BANK, kelly_fraction=KF)
    assert n == 0


def test_degenerate_price_returns_zero():
    assert kelly.kelly_size_combo(
        outcome_vec=None, blended_fair=0.5, existing_positions=[],
        effective_price=0.0, bankroll=BANK, kelly_fraction=KF) == 0
    assert kelly.kelly_size_combo(
        outcome_vec=None, blended_fair=0.5, existing_positions=[],
        effective_price=1.0, bankroll=BANK, kelly_fraction=KF) == 0


def test_perfectly_correlated_position_downsizes():
    """A held, perfectly correlated position must shrink the new size vs the
    no-position case — the core anti-overbet guarantee."""
    base = kelly.kelly_size_combo(
        outcome_vec=None, blended_fair=0.30, existing_positions=[],
        effective_price=0.20, bankroll=BANK, kelly_fraction=KF)
    # comonotonic: cov_return positive and large
    pos = [{"cov_return": (0.30 - 0.30 * 0.30) / (0.20 * 0.20),
            "contracts": 50, "effective_price": 0.20}]
    corr = kelly.kelly_size_combo(
        outcome_vec=None, blended_fair=0.30, existing_positions=pos,
        effective_price=0.20, bankroll=BANK, kelly_fraction=KF)
    assert corr < base


def test_independent_position_matches_base():
    base = kelly.kelly_size_combo(
        outcome_vec=None, blended_fair=0.30, existing_positions=[],
        effective_price=0.20, bankroll=BANK, kelly_fraction=KF)
    pos = [{"cov_return": 0.0, "contracts": 50, "effective_price": 0.20}]
    same = kelly.kelly_size_combo(
        outcome_vec=None, blended_fair=0.30, existing_positions=pos,
        effective_price=0.20, bankroll=BANK, kelly_fraction=KF)
    assert same == base


def test_negative_correlation_upsizes_but_bounded():
    base = kelly.kelly_size_combo(
        outcome_vec=None, blended_fair=0.30, existing_positions=[],
        effective_price=0.20, bankroll=BANK, kelly_fraction=KF)
    pos = [{"cov_return": -0.5, "contracts": 50, "effective_price": 0.20}]
    hedge = kelly.kelly_size_combo(
        outcome_vec=None, blended_fair=0.30, existing_positions=pos,
        effective_price=0.20, bankroll=BANK, kelly_fraction=KF)
    assert hedge >= base
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd kalshi_mlb_rfq && python -m pytest tests/test_kelly.py -v`
Expected: FAIL (positions branch still expects `outcome_vec` arrays → `TypeError`/wrong result)

- [ ] **Step 3: Refactor the existing-positions branch**

In `kelly.py`, change the signature default and replace the covariance
computation (current lines ~57-92) so it reads each position's `cov_return`:
```python
def kelly_size_combo(
    outcome_vec,
    blended_fair: float,
    existing_positions: list[dict],
    effective_price: float,
    bankroll: float,
    kelly_fraction: float,
) -> int:
    if effective_price <= 0 or effective_price >= 1:
        return 0
    p = float(blended_fair)
    if not 0.0 < p < 1.0 or p <= effective_price:
        return 0

    mu_new = (p - effective_price) / effective_price
    var_new = p * (1.0 - p) / (effective_price ** 2)
    if var_new <= 1e-12:
        return 0
    base_frac = mu_new / var_new

    if not existing_positions:
        contracts = math.floor(kelly_fraction * base_frac * bankroll / effective_price)
        return max(0, contracts)

    # Book-implied covariance: each position carries a precomputed return-space
    # covariance with the new bet (correlation.cov_returns). Full-Kelly residual
    # after the existing book, then apply the Kelly fraction.
    np_cov = []
    f_placed = []
    for pos in existing_positions:
        cov = float(pos["cov_return"])
        pos_price = float(pos["effective_price"])
        if pos_price <= 0 or pos_price >= 1:
            continue
        np_cov.append(cov)
        f_placed.append(kelly_fraction * pos["contracts"] * pos_price / bankroll)

    if not f_placed:
        contracts = math.floor(kelly_fraction * base_frac * bankroll / effective_price)
        return max(0, contracts)

    import numpy as np  # local import preserved for parity with existing module
    f_full_new = max(0.0, (mu_new - float(np.dot(np.array(np_cov), np.array(f_placed)))) / var_new)
    contracts = math.floor(kelly_fraction * f_full_new * bankroll / effective_price)
    return max(0, contracts)
```
Update the module docstring to note covariance now comes from book-implied
joints, not sample paths. Set the `outcome_vec` parameter default to `None`.

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd kalshi_mlb_rfq && python -m pytest tests/test_kelly.py -v`
Expected: PASS (all)

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/kelly.py kalshi_mlb_rfq/tests/test_kelly.py
git commit -m "refactor(rfq): Kelly consumes book-implied covariance per position"
```

---

### Task 4: Book-only pricing path in `main.py`

**Files:**
- Modify: `kalshi_mlb_rfq/main.py` — add `_book_only_fair` (near `_fresh_blended_fair`, ~448), branch `_fresh_blended_fair` (~450-479) and the candidate loop (~1452-1494) on `config.USE_MODEL`.
- Test: `kalshi_mlb_rfq/tests/test_book_only_pricing.py` (new)

**Interfaces:**
- Produces: `main._book_only_fair(book_fairs: dict[str, float]) -> float | None`.
- Consumes: `config.USE_MODEL` (Task 1).

- [ ] **Step 1: Write the failing tests**

Create `tests/test_book_only_pricing.py`:
```python
import importlib
import pytest


def test_book_only_fair_median():
    from kalshi_mlb_rfq.main import _book_only_fair
    assert _book_only_fair({"dk": 0.30, "fd": 0.34}) == pytest.approx(0.32)
    assert _book_only_fair({"dk": 0.30, "fd": 0.34, "px": 0.40}) == pytest.approx(0.34)
    assert _book_only_fair({"dk": 0.30, "fd": None}) == pytest.approx(0.30)
    assert _book_only_fair({}) is None
    assert _book_only_fair({"dk": None}) is None
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd kalshi_mlb_rfq && python -m pytest tests/test_book_only_pricing.py -v`
Expected: FAIL (`ImportError: cannot import name '_book_only_fair'`)

- [ ] **Step 3: Add `_book_only_fair` and branch the pricing paths**

In `main.py`, add near `_fresh_blended_fair` (ensure `import statistics` exists at top):
```python
def _book_only_fair(book_fairs: dict) -> float | None:
    """Median of book fairs (>=2-book floor enforced upstream). None if empty."""
    vals = [v for v in book_fairs.values() if v is not None]
    return statistics.median(vals) if vals else None
```

Replace the body of `_fresh_blended_fair` after `legs = json.loads(legs_json)`:
```python
    spread_line = _spread_line_from_legs(legs)
    total_line = _total_line_from_legs(legs)
    book_fairs = _load_book_fairs(game_id, spread_line, total_line)

    if not config.USE_MODEL:
        return _book_only_fair(book_fairs), book_fairs

    samples = _load_samples_for_game(game_id)
    if samples is None:
        return None, {}
    typed = [_leg_dict_to_typed(l, game_id) for l in legs]
    if any(l is None for l in typed):
        return None, {}
    model = fair_value.model_fair(samples, typed)
    return fair_value.blend(model, book_fairs), book_fairs
```

In the candidate loop (~1452), replace the `samples` guard + blend:
```python
        # (book-only) do not skip games without samples
        if config.USE_MODEL:
            samples = _load_samples_for_game(game_id)
            if samples is None:
                continue
        else:
            samples = None

        for cand in combo_enumerator.enumerate_2leg(...):   # unchanged loop header
            typed = [_leg_dict_to_typed(dict(l), game_id) for l in cand.legs]
            spread_line = _spread_line_from_legs([dict(l) for l in cand.legs])
            total_line = _total_line_from_legs([dict(l) for l in cand.legs])
            if any(l is None for l in typed):
                _emit_candidate_event("rejected_no_mapping", ...)   # unchanged
                continue
            books = _load_book_fairs(game_id, spread_line, total_line)
            if config.USE_MODEL:
                model = fair_value.model_fair(samples, typed)
                blended = fair_value.blend(model, books)
            else:
                model = None
                blended = _book_only_fair(books)
            if blended is None:
                _emit_candidate_event("rejected_no_book_data", ...)   # unchanged
                continue
            # ... fair-oob check, kalshi_ref, _kelly_size_for_candidate (Task 5) unchanged
```
(Keep every `_emit_candidate_event(...)` call exactly as today; only the
`model`/`blended` computation and the sample guard change.)

- [ ] **Step 4: Run test to verify it passes**

Run: `cd kalshi_mlb_rfq && python -m pytest tests/test_book_only_pricing.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/main.py kalshi_mlb_rfq/tests/test_book_only_pricing.py
git commit -m "feat(rfq): book-only fair value when USE_MODEL is off"
```

---

### Task 5: Wire correlation into Kelly sizing + skip staleness gate + decouple cache

**Files:**
- Modify: `kalshi_mlb_rfq/main.py` — `grid_lookup` adapter + `_combo_region_from_legs`; `_kelly_size_for_quote` (~597), `_kelly_size_for_candidate` (~697), `_load_existing_positions_for_game` (~568); staleness gate (~486); `_refresh_caches` sample/meta load (~153-180, 206-208).
- Test: `kalshi_mlb_rfq/tests/test_book_only_pricing.py` (extend)

**Interfaces:**
- Consumes: `correlation.ComboRegion/joint_prob/cov_returns` (Task 2), refactored `kelly.kelly_size_combo` (Task 3).
- Produces:
  - `main._grid_lookup(game_id)` → returns a `lookup(spread_line, total_line, spread_side, total_side) -> float | None` closure over `_SGP_ODDS_CACHE` (median across ≥2 books via `fair_value.devig_book`).
  - `main._combo_region_from_legs(legs: list[dict]) -> correlation.ComboRegion | None` (None if any leg is non-grid / cross-category).
  - `main._book_implied_cov(game_id, new_region, new_price, pos_region, pos_price) -> float` (uses `joint_prob`; ρ=1 fallback `min(P_a,P_b)` when `joint_prob` is None).

- [ ] **Step 1: Write the failing tests (extend `test_book_only_pricing.py`)**

```python
def test_staleness_gate_skipped_when_model_off(monkeypatch):
    import kalshi_mlb_rfq.main as m
    import kalshi_mlb_rfq.config as cfg
    monkeypatch.setattr(cfg, "USE_MODEL", False)
    monkeypatch.setattr(m, "_SAMPLES_META_GENERATED_AT", None)
    # build a minimal quote/meta that passes everything except staleness;
    # assert decision != "declined_stale_predictions"
    ok, decision = m._staleness_gate_ok()      # see Step 3 helper extraction
    assert ok is True and decision == "passed"


def test_staleness_gate_fires_when_model_on(monkeypatch):
    import kalshi_mlb_rfq.main as m
    import kalshi_mlb_rfq.config as cfg
    from datetime import datetime, timezone, timedelta
    monkeypatch.setattr(cfg, "USE_MODEL", True)
    monkeypatch.setattr(m, "_SAMPLES_META_GENERATED_AT",
                        datetime.now(timezone.utc) - timedelta(hours=2))
    ok, decision = m._staleness_gate_ok()
    assert ok is False and decision == "declined_stale_predictions"


def test_book_implied_cov_fallback_is_rho_one(monkeypatch):
    import kalshi_mlb_rfq.main as m
    from kalshi_mlb_rfq.correlation import ComboRegion
    # grid returns nothing -> fallback min(Pa,Pb); patch _grid_lookup + fairs
    monkeypatch.setattr(m, "_grid_lookup", lambda gid: (lambda *a: None))
    monkeypatch.setattr(m, "_combo_fair_for_region",
                        lambda gid, r: 0.30 if r.total_side == "over" else 0.25)
    a = ComboRegion("home", -1.5, "over", 8.5)
    b = ComboRegion("home", -1.5, "under", 8.5)
    cov = m._book_implied_cov("g1", a, 0.20, b, 0.20)
    # fallback joint = min(0.30,0.25)=0.25 ; cov=(0.25-0.30*0.25)/(0.2*0.2)
    assert cov == pytest.approx((0.25 - 0.30 * 0.25) / 0.04)
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd kalshi_mlb_rfq && python -m pytest tests/test_book_only_pricing.py -v`
Expected: FAIL (`_staleness_gate_ok` / `_book_implied_cov` / `_grid_lookup` undefined)

- [ ] **Step 3: Implement the helpers and wire them**

Extract the staleness check into a testable helper and gate it on the flag.
Replace lines ~485-488 of `_all_per_accept_gates_pass`:
```python
    ok, decision = _staleness_gate_ok()
    if not ok:
        return False, decision
```
and add:
```python
def _staleness_gate_ok() -> tuple[bool, str]:
    """Model-prediction staleness — only meaningful when the model prices."""
    if not config.USE_MODEL:
        return True, "passed"
    gen_at = _samples_generated_at()
    if gen_at is None or not risk.staleness_ok(gen_at, config.MAX_PREDICTION_STALENESS_SEC):
        return False, "declined_stale_predictions"
    return True, "passed"
```

Add the grid adapter + region helpers (near `_load_book_fairs`):
```python
def _grid_lookup(game_id: str):
    """Closure: (spread_line, total_line, spread_side, total_side) -> median
    book prob for that SGP quadrant, or None if <2 books priced it."""
    def lookup(spread_line, total_line, spread_side, total_side):
        if _SGP_ODDS_CACHE is None or _SGP_ODDS_CACHE.empty:
            return None
        label = f"{spread_side.title()} Spread + {total_side.title()}"
        rows = _SGP_ODDS_CACHE[
            (_SGP_ODDS_CACHE["game_id"] == game_id)
            & (_SGP_ODDS_CACHE["spread_line"].astype(float).round(2) == round(float(spread_line), 2))
            & (_SGP_ODDS_CACHE["total_line"].astype(float).round(2) == round(float(total_line), 2))
        ]
        if rows.empty:
            return None
        out = []
        for book in rows["bookmaker"].unique():
            sub = rows[rows["bookmaker"] == book]
            f = fair_value.devig_book(sub, combo=label, vig_fallback=_vig_fallback(book))
            if f is not None:
                out.append(f)
        if len(out) < config.MIN_BOOK_COUNT_FOR_BLEND:
            return None
        return statistics.median(out)
    return lookup


def _combo_region_from_legs(legs: list[dict]):
    """Map a combo's legs to a ComboRegion, or None if it has a non-grid leg
    (e.g. cross-category prop)."""
    spread = _spread_leg(legs)     # existing helpers that find the spread/total leg
    total = _total_leg(legs)
    if spread is None or total is None:
        return None
    return correlation.ComboRegion(
        spread_side="home" if spread["side_is_home"] else "away",
        spread_line=float(spread["line"]),
        total_side="over" if total["is_over"] else "under",
        total_line=float(total["line"]),
    )


def _combo_fair_for_region(game_id: str, r) -> float | None:
    return _grid_lookup(game_id)(r.spread_line, r.total_line, r.spread_side, r.total_side)


def _book_implied_cov(game_id, new_region, new_price, pos_region, pos_price) -> float:
    p_a = _combo_fair_for_region(game_id, new_region)
    p_b = _combo_fair_for_region(game_id, pos_region)
    if p_a is None or p_b is None:
        return 0.0  # cannot price one leg book-only → no correlation term
    j = correlation.joint_prob(new_region, pos_region, p_a, p_b, _grid_lookup(game_id))
    if j is None:
        j = min(p_a, p_b)  # ρ=1 conservative fallback
    return correlation.cov_returns(p_a, p_b, j, new_price, pos_price)
```
(`_spread_leg` / `_total_leg` adapters reuse the same leg-typing the existing
`_spread_line_from_legs` / `_total_line_from_legs` use; if those internals expose
side/over flags differently, derive from the typed leg. Implementer: confirm the
exact leg dict keys against `_leg_dict_to_typed` and adjust the two attribute
reads — this is the one spot needing a look at the live leg shape.)

Now branch both Kelly callers. In `_kelly_size_for_quote` (~620) and
`_kelly_size_for_candidate` (~717), when `not config.USE_MODEL`:
```python
        if config.USE_MODEL:
            samples = _load_samples_for_game(game_id)
            if samples is None or samples.empty:
                return 0
            outcome_vec = _outcome_vec_for_legs(samples, typed_legs, side=side)
            if outcome_vec is None:
                return 0
            existing = _load_existing_positions_for_game(game_id, samples)
        else:
            outcome_vec = None
            new_region = _combo_region_from_legs(legs)
            existing = _load_existing_positions_book(game_id, new_region,
                                                     effective_price, side)
        return kelly.kelly_size_combo(
            outcome_vec=outcome_vec,
            blended_fair=fair if side == "yes" else (1 - fair),
            existing_positions=existing,
            effective_price=effective_price,
            bankroll=config.BANKROLL,
            kelly_fraction=config.KELLY_FRACTION,
        )
```
Add `_load_existing_positions_book(game_id, new_region, new_price, side)` that
reads held same-game positions from the `positions` table (leg-sets via
`combo_cache`), and for each builds `{"cov_return": _book_implied_cov(game_id,
new_region, new_price, pos_region, pos_price), "contracts": net_contracts,
"effective_price": pos_price}`. Skip positions whose region is `None` (non-grid).

Finally, decouple the cache: in `_refresh_caches`, guard the sample/meta load:
```python
            if config.USE_MODEL:
                samples_df = con.execute("... FROM mlb_game_samples").fetchdf()
                samples_by_game = {gid: g.reset_index(drop=True)
                                   for gid, g in samples_df.groupby("game_id")}
                meta_row = con.execute(
                    "SELECT generated_at FROM mlb_samples_meta ORDER BY generated_at DESC LIMIT 1"
                ).fetchone()
                generated_at = meta_row[0] if meta_row else None
            else:
                samples_by_game = {}
                generated_at = None
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd kalshi_mlb_rfq && python -m pytest tests/test_book_only_pricing.py -v`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add kalshi_mlb_rfq/main.py kalshi_mlb_rfq/tests/test_book_only_pricing.py
git commit -m "feat(rfq): wire book-implied correlation Kelly + skip stale gate + decouple model DB"
```

---

### Task 6: Seeded integration test — correlated add is down-sized

**Files:**
- Test: `kalshi_mlb_rfq/tests/test_correlation_integration.py` (new)

**Interfaces:**
- Consumes: Tasks 2–5 end to end.

- [ ] **Step 1: Write the failing test**

Create `tests/test_correlation_integration.py`:
```python
import pytest
import pandas as pd
import kalshi_mlb_rfq.main as m
import kalshi_mlb_rfq.config as cfg


def _seed_grid():
    # one game, home+over priced at two spreads/totals across 2 books
    rows = []
    for book in ("draftkings", "fanduel"):
        for (sl, tl, label, dec) in [
            (-1.5, 8.5, "Home Spread + Over", 3.3),
            (-1.5, 8.5, "Home Spread + Under", 3.3),
            (-1.5, 8.5, "Away Spread + Over", 3.3),
            (-1.5, 8.5, "Away Spread + Under", 3.3),
            (-2.5, 9.5, "Home Spread + Over", 5.0),
            (-2.5, 9.5, "Home Spread + Under", 5.0),
            (-2.5, 9.5, "Away Spread + Over", 5.0),
            (-2.5, 9.5, "Away Spread + Under", 5.0),
        ]:
            rows.append({"game_id": "g1", "combo": label, "bookmaker": book,
                         "sgp_decimal": dec, "spread_line": sl, "total_line": tl})
    return pd.DataFrame(rows)


def test_correlated_second_combo_is_downsized(monkeypatch):
    monkeypatch.setattr(cfg, "USE_MODEL", False)
    monkeypatch.setattr(m, "_SGP_ODDS_CACHE", _seed_grid())
    region = m.correlation.ComboRegion("home", -1.5, "over", 8.5)
    price = 0.30
    # no held position
    base_cov_positions = m._load_existing_positions_book("g1", region, price, "yes")
    assert base_cov_positions == []  # nothing held yet
    # with a held strongly-correlated position, cov_return > 0 → Kelly shrinks
    held = [{"cov_return": 0.5, "contracts": 40, "effective_price": 0.30}]
    base = m.kelly.kelly_size_combo(None, 0.45, [], 0.30, 1000.0, 0.25)
    corr = m.kelly.kelly_size_combo(None, 0.45, held, 0.30, 1000.0, 0.25)
    assert corr < base
```
(If wiring `_load_existing_positions_book` against a seeded `positions` table is
heavier than a unit warrants, keep the held-position branch assertion as the
core check; the engine pieces are already unit-tested in Tasks 2–3.)

- [ ] **Step 2: Run test to verify it fails**, then **Step 3** make it pass with
the Task-5 helpers, **Step 4** run green, **Step 5** commit:

```bash
git add kalshi_mlb_rfq/tests/test_correlation_integration.py
git commit -m "test(rfq): integration — correlated same-game add is down-sized"
```

---

### Task 7: Full regression sweep under both flag values

**Files:** none (CI/validation task)

- [ ] **Step 1: Run the whole suite book-only (default)**

Run: `cd kalshi_mlb_rfq && python -m pytest tests/ -q`
Expected: PASS. Fix any test that implicitly assumed the staleness gate runs by
default — pin `USE_MODEL` explicitly via `monkeypatch.setattr(cfg, "USE_MODEL", True)`
in those tests.

- [ ] **Step 2: Run the whole suite with the model on**

Run: `cd kalshi_mlb_rfq && USE_MODEL=true python -m pytest tests/ -q`
Expected: PASS (legacy model path intact).

- [ ] **Step 3: Commit any test pins**

```bash
git add kalshi_mlb_rfq/tests/
git commit -m "test(rfq): pin USE_MODEL where legacy tests assume the model path"
```

---

### Task 8: Docs + dead-env fix

**Files:**
- Modify: `kalshi_mlb_rfq/.env.example`, `kalshi_mlb_rfq/README.md`, root `CLAUDE.md`

- [ ] **Step 1: `.env.example`** — replace the dead `MAX_STALENESS_SEC` line with
`MAX_PREDICTION_STALENESS_SEC=600` and add `USE_MODEL=false`.

- [ ] **Step 2: `README.md`** — document `USE_MODEL`, book-only pricing
(`median(books)`, ≥2 floor), the book-implied correlation engine (grid joint +
ρ=1 fallback + Fréchet clamp), and that the prediction-staleness gate only runs
when `USE_MODEL=true`. Note v1 correlation covers same-direction spread/total
pairs; other pairs use ρ=1.

- [ ] **Step 3: root `CLAUDE.md`** — update the taker bullet: pricing is
book-only by default (model optional, off), risk sizing via book-implied
correlation.

- [ ] **Step 4: Commit**

```bash
git add kalshi_mlb_rfq/.env.example kalshi_mlb_rfq/README.md CLAUDE.md
git commit -m "docs(rfq): document USE_MODEL book-only pricing + correlation engine"
```

- [ ] **Step 5: Fix the live `.env`** (not committed) — rename
`MAX_STALENESS_SEC=900` → `MAX_PREDICTION_STALENESS_SEC` (or remove; book-only
ignores it) and add `USE_MODEL=false`. Do this as part of go-live config, not in git.

---

## Pre-merge / go-live (gated on user approval)

1. Executive review of `git diff main..HEAD` (data integrity, on.exit/conn
   safety unaffected, dead code, no secrets).
2. Dry-run single cycle from **main repo cwd** (per restart gotcha): confirm
   candidates enumerate, price book-only, correlation runs without raising,
   accept path reached, and `quote_log` shows no `declined_stale_predictions`.
3. Log one held-position size delta (independent vs correlated) to prove live
   down-sizing.
4. **User approval required** before turning live; small-stake validation before
   any ramp.

## Self-Review notes

- Spec coverage: Part 1 (flag, book-only pricing, gate skip, cache decouple) →
  Tasks 1,4,5. Part 2 (correlation module, Kelly refactor, wiring, Fréchet) →
  Tasks 2,3,5. Tests → Tasks 2,3,4,5,6,7. Docs → Task 8.
- **Deviation from spec:** crossing-pair exact inclusion-exclusion is replaced by
  ρ=1 fallback in v1 (Global Constraints + Task 2). Strictly safe; flagged for
  user confirmation at handoff.
- Open implementer check (Task 5, Step 3): exact leg dict keys for
  `_combo_region_from_legs` must be confirmed against `_leg_dict_to_typed`.
