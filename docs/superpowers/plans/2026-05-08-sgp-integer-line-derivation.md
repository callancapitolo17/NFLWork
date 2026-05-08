# SGP Integer-Line Derivation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** When a comparison sportsbook (DK/FD/NV/PX) doesn't quote the integer total line WZ is using in an SGP, derive that book's joint fair probability for the 4 SGP combos at the integer line from its two adjacent half-point alts, write the derived probabilities to `mlb_sgp_odds` with an `_interpolated` source tag, and pick them up in the R blender alongside `_direct` rows.

**Architecture:** A new shared Python module `mlb_sgp/integer_line_derivation.py` owns the math (devig, push-mass derivation, bounds checks, fair-prob formula). Each of the 4 SGP scrapers gets a fallback block (~30 LOC) inserted before the existing "skip if exact line missing" branch — when WZ posts an integer line and both adjacent half-point alts exist in the book's selection-ID dictionary, the scraper fires 8 pricing calls (4 combos × 2 alts), passes results to the shared module, and writes 4 derived rows. The R blender's WHERE clause is widened by one character (a `LIKE` pattern) to ingest `_interpolated` rows identically to `_direct` rows.

**Tech Stack:** Python 3.10+ (`mlb_sgp/venv`), `curl_cffi` for HTTP, `duckdb` for storage, R 4.x for the blender, `pytest` for unit tests.

**Spec:** `docs/superpowers/specs/2026-05-08-sgp-integer-line-derivation-design.md`

**Worktree:** `/Users/callancapitolo/NFLWork/.worktrees/sgp-integer-derivation` on branch `feature/sgp-integer-line-derivation`

---

## File Structure

**Create:**
- `mlb_sgp/integer_line_derivation.py` — shared math + bounds checks + orchestrator function
- `mlb_sgp/tests/test_integer_line_derivation.py` — unit tests for the shared module
- `mlb_sgp/tests/test_integration_integer_line.py` — manual integration test (not added to cron)

**Modify:**
- `mlb_sgp/scraper_draftkings_sgp.py` — add fallback block before existing skip
- `mlb_sgp/scraper_fanduel_sgp.py` — same
- `mlb_sgp/scraper_novig_sgp.py` — same
- `mlb_sgp/scraper_prophetx_sgp.py` — same
- `Answer Keys/mlb_correlated_parlay.R` — widen source WHERE clause
- `mlb_sgp/README.md` — document the integer-line derivation behavior

---

## Task 1: Shared module — `is_integer_line` helper

**Files:**
- Create: `mlb_sgp/integer_line_derivation.py`
- Create: `mlb_sgp/tests/__init__.py` (if not present — empty file)
- Test: `mlb_sgp/tests/test_integer_line_derivation.py`

- [ ] **Step 1: Create the test file with the first failing test**

```python
# mlb_sgp/tests/test_integer_line_derivation.py
import pytest
from mlb_sgp.integer_line_derivation import is_integer_line


def test_is_integer_line_8_dot_0():
    assert is_integer_line(8.0) is True

def test_is_integer_line_8_dot_5():
    assert is_integer_line(8.5) is False

def test_is_integer_line_handles_fp_jitter():
    # WZ data sometimes carries floating-point noise (8.000001, 7.999999)
    assert is_integer_line(8.000001) is True
    assert is_integer_line(7.999999) is True

def test_is_integer_line_rejects_close_to_half():
    assert is_integer_line(8.5001) is False
    assert is_integer_line(8.4999) is False
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/sgp-integer-derivation
source mlb_sgp/venv/bin/activate
pytest mlb_sgp/tests/test_integer_line_derivation.py::test_is_integer_line_8_dot_0 -v
```

Expected: FAIL with `ImportError: No module named mlb_sgp.integer_line_derivation` or similar.

- [ ] **Step 3: Create the module with minimal `is_integer_line`**

```python
# mlb_sgp/integer_line_derivation.py
"""
Derive joint fair probabilities for SGP combos at integer total lines
from two adjacent half-point alt SGP prices.

See docs/superpowers/specs/2026-05-08-sgp-integer-line-derivation-design.md
"""

# Floating-point tolerance for matching integer lines (WZ data sometimes
# has 8.000001 / 7.999999 instead of exact 8.0)
_INT_TOLERANCE = 1e-3


def is_integer_line(value: float) -> bool:
    """True if value is at an integer (within FP tolerance), False if half-point."""
    return abs(value - round(value)) < _INT_TOLERANCE
```

- [ ] **Step 4: Run all 4 tests in the file to verify they pass**

```bash
pytest mlb_sgp/tests/test_integer_line_derivation.py -v
```

Expected: 4 passed.

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/integer_line_derivation.py mlb_sgp/tests/test_integer_line_derivation.py mlb_sgp/tests/__init__.py
git commit -m "mlb_sgp: scaffold integer_line_derivation module with is_integer_line helper"
```

---

## Task 2: Shared module — per-alt devig

**Files:**
- Modify: `mlb_sgp/integer_line_derivation.py`
- Test: `mlb_sgp/tests/test_integer_line_derivation.py`

- [ ] **Step 1: Add a failing test for per-alt devig**

Append to `mlb_sgp/tests/test_integer_line_derivation.py`:

```python
from mlb_sgp.integer_line_derivation import devig_alt_set


def test_devig_alt_set_balances_to_one():
    # 4 combos at one alt total. Sum of 1/decimal = 1.125 (typical DK vig).
    decimals = {
        "home_over": 2.778,    # 1/2.778 ≈ 0.360
        "home_under": 5.556,   # 1/5.556 ≈ 0.180
        "away_over": 2.564,    # 1/2.564 ≈ 0.390
        "away_under": 5.128,   # 1/5.128 ≈ 0.195
    }
    devigged, vig_sum = devig_alt_set(decimals)
    assert vig_sum == pytest.approx(1.125, abs=1e-3)
    assert sum(devigged.values()) == pytest.approx(1.0, abs=1e-6)
    assert devigged["home_over"] == pytest.approx(0.320, abs=1e-3)


def test_devig_alt_set_preserves_keys():
    decimals = {"a": 2.0, "b": 2.0, "c": 2.0, "d": 2.0}
    devigged, _ = devig_alt_set(decimals)
    assert set(devigged.keys()) == {"a", "b", "c", "d"}
    assert all(v == pytest.approx(0.25) for v in devigged.values())
```

- [ ] **Step 2: Run to verify it fails**

```bash
pytest mlb_sgp/tests/test_integer_line_derivation.py::test_devig_alt_set_balances_to_one -v
```

Expected: FAIL with `ImportError: cannot import name 'devig_alt_set'`.

- [ ] **Step 3: Implement `devig_alt_set`**

Append to `mlb_sgp/integer_line_derivation.py`:

```python
def devig_alt_set(decimals: dict[str, float]) -> tuple[dict[str, float], float]:
    """
    Multiplicative devig of one alt-set's 4 combo decimal odds.

    Returns (devigged_probs, vig_sum). Each devigged_probs[k] is the joint
    probability of combo k under the multiplicative method:
        p_devigged = (1/decimal_k) / sum(1/decimal_i for all i)

    The 4 combos must form a partition of the joint sample space (Home Spread
    + Over, Home Spread + Under, Away Spread + Over, Away Spread + Under at
    one alt total). The returned devigged probs sum to exactly 1.0.

    The vig_sum (= sum of raw implied probs) is also returned for the
    per-alt vig bounds check.
    """
    implied = {k: 1.0 / d for k, d in decimals.items()}
    vig_sum = sum(implied.values())
    devigged = {k: p / vig_sum for k, p in implied.items()}
    return devigged, vig_sum
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest mlb_sgp/tests/test_integer_line_derivation.py -v
```

Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/integer_line_derivation.py mlb_sgp/tests/test_integer_line_derivation.py
git commit -m "mlb_sgp: add devig_alt_set to integer_line_derivation module"
```

---

## Task 3: Shared module — bounds-check threshold constants and primitives

**Files:**
- Modify: `mlb_sgp/integer_line_derivation.py`
- Test: `mlb_sgp/tests/test_integer_line_derivation.py`

- [ ] **Step 1: Add failing tests for each bounds-check primitive**

Append to `mlb_sgp/tests/test_integer_line_derivation.py`:

```python
from mlb_sgp.integer_line_derivation import (
    validate_per_alt_vig,
    validate_push_mass_consistency,
    validate_delta_total,
    validate_sum_to_one,
    validate_per_combo_bounds,
    VIG_MIN, VIG_MAX,
    PUSH_MASS_REL_TOL,
    DELTA_TOTAL_MIN, DELTA_TOTAL_MAX,
    SUM_MIN, SUM_MAX,
)


def test_validate_per_alt_vig_in_range():
    assert validate_per_alt_vig(1.125) is True
    assert validate_per_alt_vig(1.05) is True
    assert validate_per_alt_vig(1.30) is True

def test_validate_per_alt_vig_out_of_range():
    assert validate_per_alt_vig(1.04) is False
    assert validate_per_alt_vig(1.31) is False
    assert validate_per_alt_vig(0.95) is False

def test_validate_push_mass_consistency_close():
    # 0.050 vs 0.053 = 5.7% relative diff, within 10% tolerance
    assert validate_push_mass_consistency(0.050, 0.053) is True

def test_validate_push_mass_consistency_diverging():
    # 0.050 vs 0.080 = 37.5% relative diff, fails
    assert validate_push_mass_consistency(0.050, 0.080) is False

def test_validate_push_mass_consistency_handles_zero():
    # Both zero is consistent (no push mass)
    assert validate_push_mass_consistency(0.0, 0.0) is True
    # One zero, one nonzero -> max=nonzero, diff=nonzero -> 100% rel = fail
    assert validate_push_mass_consistency(0.0, 0.05) is False

def test_validate_delta_total_in_range():
    assert validate_delta_total(0.095) is True
    assert validate_delta_total(0.03) is True
    assert validate_delta_total(0.18) is True

def test_validate_delta_total_out_of_range():
    assert validate_delta_total(0.02) is False
    assert validate_delta_total(0.20) is False

def test_validate_sum_to_one_in_range():
    assert validate_sum_to_one([0.30, 0.20, 0.30, 0.20]) is True
    assert validate_sum_to_one([0.97]) is True   # single-value edge
    assert validate_sum_to_one([1.03]) is True

def test_validate_sum_to_one_drift():
    assert validate_sum_to_one([0.20, 0.20, 0.20, 0.20]) is False  # sums to 0.8
    assert validate_sum_to_one([0.30, 0.30, 0.30, 0.30]) is False  # sums to 1.2

def test_validate_per_combo_bounds():
    assert validate_per_combo_bounds(0.5) is True
    assert validate_per_combo_bounds(0.001) is True
    assert validate_per_combo_bounds(0.999) is True
    assert validate_per_combo_bounds(0.0) is False
    assert validate_per_combo_bounds(1.0) is False
    assert validate_per_combo_bounds(-0.01) is False
    assert validate_per_combo_bounds(1.01) is False
```

- [ ] **Step 2: Run to verify they fail**

```bash
pytest mlb_sgp/tests/test_integer_line_derivation.py -v -k validate
```

Expected: ImportError for the missing names.

- [ ] **Step 3: Implement primitives and constants**

Append to `mlb_sgp/integer_line_derivation.py`:

```python
# ---------------------------------------------------------------------------
# Bounds-check thresholds (config constants — tunable in follow-up commits
# once production logs reveal real-world distributions).
# ---------------------------------------------------------------------------

# Per-alt vig sum bounds. DK ~1.125, FD bimodal 1.13/1.22, PX/NV unknown.
VIG_MIN = 1.05
VIG_MAX = 1.30

# Push-mass cross-consistency: |Δ_from_over - Δ_from_under| / max(...) tolerance
PUSH_MASS_REL_TOL = 0.10

# Total marginal push mass plausibility bounds (fraction of games landing on X)
DELTA_TOTAL_MIN = 0.03
DELTA_TOTAL_MAX = 0.18

# Sum of 4 derived fair_probs should ≈ 1.0 (mathematical invariant)
SUM_MIN = 0.97
SUM_MAX = 1.03


def validate_per_alt_vig(vig_sum: float) -> bool:
    return VIG_MIN <= vig_sum <= VIG_MAX


def validate_push_mass_consistency(delta_a: float, delta_b: float) -> bool:
    """The push mass derived from the Over side must approximately equal the
    push mass derived from the Under side (they're the same joint event).
    Tolerance: PUSH_MASS_REL_TOL relative to the larger of the two.
    """
    max_delta = max(abs(delta_a), abs(delta_b))
    if max_delta == 0:
        return delta_a == delta_b == 0
    return abs(delta_a - delta_b) / max_delta < PUSH_MASS_REL_TOL


def validate_delta_total(delta_total: float) -> bool:
    return DELTA_TOTAL_MIN <= delta_total <= DELTA_TOTAL_MAX


def validate_sum_to_one(fair_probs: list[float]) -> bool:
    s = sum(fair_probs)
    return SUM_MIN <= s <= SUM_MAX


def validate_per_combo_bounds(fair_prob: float) -> bool:
    return 0.0 < fair_prob < 1.0
```

- [ ] **Step 4: Run tests**

```bash
pytest mlb_sgp/tests/test_integer_line_derivation.py -v
```

Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/integer_line_derivation.py mlb_sgp/tests/test_integer_line_derivation.py
git commit -m "mlb_sgp: add bounds-check primitives + thresholds to integer_line_derivation"
```

---

## Task 4: Shared module — `derive_fair_probs` orchestrator

**Files:**
- Modify: `mlb_sgp/integer_line_derivation.py`
- Test: `mlb_sgp/tests/test_integer_line_derivation.py`

This is the main entry point each scraper will call. It takes already-computed decimal odds for 4 combos at each of 2 alt levels, runs all 5 bounds checks, applies the formula, and returns either the 4 derived fair_probs or `None` (with a structured log line on rejection).

- [ ] **Step 1: Write failing test for the worked-example happy path**

Append to `mlb_sgp/tests/test_integer_line_derivation.py`:

```python
from mlb_sgp.integer_line_derivation import derive_fair_probs


def _build_synthetic_alt(home_over_p, home_under_p, away_over_p, away_under_p, vig=1.125):
    """Helper: turn 4 desired devigged probs into decimal odds for an alt set."""
    # Each devigged p corresponds to raw implied = p * vig
    return {
        "home_over": 1.0 / (home_over_p * vig),
        "home_under": 1.0 / (home_under_p * vig),
        "away_over": 1.0 / (away_over_p * vig),
        "away_under": 1.0 / (away_under_p * vig),
    }


def test_derive_fair_probs_worked_example():
    """Yankees @ Red Sox example from spec §3.5.

    Devigged at Over 7.5: HO=0.320, HU=0.160, AO=0.347, AU=0.173 (sum=1.0)
    Devigged at Over 8.5: HO=0.270, HU=0.213, AO=0.302, AU=0.218 (sum~1.003)
    Δ_cover ≈ 0.05, Δ_uncover ≈ 0.045, Δ_total ≈ 0.095
    Expected fair_prob_HSO at integer 8: 0.270 / (1 - 0.095) ≈ 0.2983
    """
    decimals_lo = _build_synthetic_alt(0.320, 0.160, 0.347, 0.173)
    decimals_hi = _build_synthetic_alt(0.270, 0.215, 0.302, 0.213)

    result = derive_fair_probs(decimals_lo, decimals_hi)

    assert result is not None
    fair = result["fair_probs"]
    assert fair["home_over"] == pytest.approx(0.2983, abs=2e-3)
    assert fair["home_under"] == pytest.approx(0.1768, abs=2e-3)
    assert fair["away_over"] == pytest.approx(0.3337, abs=2e-3)
    assert fair["away_under"] == pytest.approx(0.1912, abs=2e-3)
    assert sum(fair.values()) == pytest.approx(1.0, abs=1e-6)
    assert result["delta_total"] == pytest.approx(0.095, abs=5e-3)


def test_derive_fair_probs_rejects_bad_vig():
    # Synthesize an alt with raw vig 1.40 (way out of range)
    decimals_lo = _build_synthetic_alt(0.30, 0.20, 0.30, 0.20, vig=1.40)
    decimals_hi = _build_synthetic_alt(0.27, 0.22, 0.30, 0.21)
    result = derive_fair_probs(decimals_lo, decimals_hi)
    assert result is None


def test_derive_fair_probs_rejects_inconsistent_push_mass():
    # Force Δ_cover_from_over = 0.05 but Δ_cover_from_under = 0.20 (way out)
    # Easiest path: hand-craft devigged probs that don't reflect a coherent
    # joint distribution.
    decimals_lo = _build_synthetic_alt(0.30, 0.10, 0.30, 0.30)
    decimals_hi = _build_synthetic_alt(0.25, 0.30, 0.25, 0.20)
    # devig_lo: HO=0.30, HU=0.10  → Δ_cover_from_over = 0.30 - 0.25 = 0.05
    # devig_hi: HU=0.30, HU_lo=0.10 → Δ_cover_from_under = 0.30 - 0.10 = 0.20
    # 0.05 vs 0.20 → 75% relative diff → fails
    result = derive_fair_probs(decimals_lo, decimals_hi)
    assert result is None


def test_derive_fair_probs_rejects_extreme_delta():
    # Force Δ_total > 0.18 (very high integer-mass)
    # Make p_lo - p_hi large for both sides
    decimals_lo = _build_synthetic_alt(0.40, 0.10, 0.40, 0.10)
    decimals_hi = _build_synthetic_alt(0.20, 0.30, 0.20, 0.30)
    # Δ_cover_from_over = 0.40 - 0.20 = 0.20 (already out of range alone)
    result = derive_fair_probs(decimals_lo, decimals_hi)
    assert result is None
```

- [ ] **Step 2: Run to verify failure**

```bash
pytest mlb_sgp/tests/test_integer_line_derivation.py::test_derive_fair_probs_worked_example -v
```

Expected: ImportError for `derive_fair_probs`.

- [ ] **Step 3: Implement `derive_fair_probs`**

Append to `mlb_sgp/integer_line_derivation.py`:

```python
import logging
from typing import Optional

logger = logging.getLogger(__name__)

# Combo keys — must match exactly across all callers
COMBO_KEYS = ("home_over", "home_under", "away_over", "away_under")


def derive_fair_probs(
    decimals_lo: dict[str, float],
    decimals_hi: dict[str, float],
) -> Optional[dict]:
    """
    Derive 4 joint fair probabilities at an integer total line X from
    decimal odds at the two adjacent half-point alts (X-0.5 and X+0.5).

    Args:
        decimals_lo: 4 decimal odds at the lower alt (X-0.5). Keys: COMBO_KEYS.
        decimals_hi: 4 decimal odds at the higher alt (X+0.5). Keys: COMBO_KEYS.

    Returns:
        dict with keys {'fair_probs': {combo: prob, ...}, 'delta_total': float}
        or None if any bounds check fails.

    On rejection, logs a structured WARN line.
    """
    # Validate inputs have the expected keys
    if set(decimals_lo.keys()) != set(COMBO_KEYS) or set(decimals_hi.keys()) != set(COMBO_KEYS):
        logger.warning(
            "derive_fair_probs: invalid combo keys lo=%s hi=%s",
            sorted(decimals_lo.keys()), sorted(decimals_hi.keys())
        )
        return None

    # Per-alt devig
    devig_lo, vig_lo = devig_alt_set(decimals_lo)
    devig_hi, vig_hi = devig_alt_set(decimals_hi)

    # Bounds check (1): per-alt vig sum
    if not validate_per_alt_vig(vig_lo):
        logger.warning("derive_fair_probs: vig_lo=%.4f out of [%g, %g]", vig_lo, VIG_MIN, VIG_MAX)
        return None
    if not validate_per_alt_vig(vig_hi):
        logger.warning("derive_fair_probs: vig_hi=%.4f out of [%g, %g]", vig_hi, VIG_MIN, VIG_MAX)
        return None

    # Compute push masses two ways for each side
    delta_cover_from_over = devig_lo["home_over"] - devig_hi["home_over"]
    delta_cover_from_under = devig_hi["home_under"] - devig_lo["home_under"]
    delta_uncover_from_over = devig_lo["away_over"] - devig_hi["away_over"]
    delta_uncover_from_under = devig_hi["away_under"] - devig_lo["away_under"]

    # Bounds check (2): push-mass cross-consistency
    if not validate_push_mass_consistency(delta_cover_from_over, delta_cover_from_under):
        logger.warning(
            "derive_fair_probs: cover-side push mass inconsistent: from_over=%.4f from_under=%.4f",
            delta_cover_from_over, delta_cover_from_under
        )
        return None
    if not validate_push_mass_consistency(delta_uncover_from_over, delta_uncover_from_under):
        logger.warning(
            "derive_fair_probs: uncover-side push mass inconsistent: from_over=%.4f from_under=%.4f",
            delta_uncover_from_over, delta_uncover_from_under
        )
        return None

    # Average the consistent estimates per side (more robust than picking one)
    delta_cover = (delta_cover_from_over + delta_cover_from_under) / 2.0
    delta_uncover = (delta_uncover_from_over + delta_uncover_from_under) / 2.0
    delta_total = delta_cover + delta_uncover

    # Bounds check (3): Δ_total reasonableness
    if not validate_delta_total(delta_total):
        logger.warning(
            "derive_fair_probs: delta_total=%.4f out of [%g, %g]",
            delta_total, DELTA_TOTAL_MIN, DELTA_TOTAL_MAX
        )
        return None

    # Apply the formula:
    #   For Over combos: fair_prob = devig_hi["...over"] / (1 - delta_total)
    #   For Under combos: fair_prob = devig_lo["...under"] / (1 - delta_total)
    denom = 1.0 - delta_total
    fair_probs = {
        "home_over":  devig_hi["home_over"]  / denom,
        "home_under": devig_lo["home_under"] / denom,
        "away_over":  devig_hi["away_over"]  / denom,
        "away_under": devig_lo["away_under"] / denom,
    }

    # Bounds check (4): sum to 1
    if not validate_sum_to_one(list(fair_probs.values())):
        logger.warning(
            "derive_fair_probs: sum=%.4f out of [%g, %g] (probs=%s)",
            sum(fair_probs.values()), SUM_MIN, SUM_MAX, fair_probs
        )
        return None

    # Bounds check (5): per-combo (0, 1)
    for k, v in fair_probs.items():
        if not validate_per_combo_bounds(v):
            logger.warning("derive_fair_probs: %s=%.4f out of (0, 1)", k, v)
            return None

    return {"fair_probs": fair_probs, "delta_total": delta_total}
```

- [ ] **Step 4: Run all tests**

```bash
pytest mlb_sgp/tests/test_integer_line_derivation.py -v
```

Expected: all tests pass.

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/integer_line_derivation.py mlb_sgp/tests/test_integer_line_derivation.py
git commit -m "mlb_sgp: add derive_fair_probs orchestrator with all 5 bounds checks"
```

---

## Task 5: DK scraper — integer-line fallback

**Files:**
- Modify: `mlb_sgp/scraper_draftkings_sgp.py`

The DK scraper currently builds `combo_items` in lines 622-670 and skips entire periods on missing exact-line lookups (line 653). We add a separate fallback path that handles integer-line games by issuing 8 pricing calls instead of 4 and writing 4 derived rows with `source = 'draftkings_interpolated'`.

DK uses `calculate_sgp(session, sp, to)` for one combo's price (returns dict with `decimal` key) — already exists in the file. We reuse it.

- [ ] **Step 1: Add the fallback helper at module scope**

Open `mlb_sgp/scraper_draftkings_sgp.py`. Below the existing `calculate_sgp` function (around line 460-ish, before `match_events`), add:

```python
from .integer_line_derivation import (
    is_integer_line, derive_fair_probs, COMBO_KEYS,
)


def try_integer_fallback_dk(
    session,
    sel_ids: dict,
    spread_line: float,
    total_line: float,
    canonical: set,
    verbose: bool = False,
) -> dict | None:
    """
    Try to derive 4 fair_probs at an integer total line from adjacent
    half-point alts. Returns dict {combo_name: fair_prob} or None.

    Mirrors the spread/total-sel lookup the main loop does for the matched
    line, but for the two adjacent half-point alts.
    """
    if not is_integer_line(total_line):
        return None

    lo_total = total_line - 0.5
    hi_total = total_line + 0.5

    # Sign convention same as main loop
    if spread_line < 0:
        home_sign, away_sign = "N", "P"
    else:
        home_sign, away_sign = "P", "N"
    spread = abs(spread_line)

    home_spread_sels = sel_ids["spreads"].get((home_sign, spread, "1")) or []
    away_spread_sels = sel_ids["spreads"].get((away_sign, spread, "3")) or []

    over_lo_sels = sel_ids["totals"].get(("O", lo_total)) or []
    under_lo_sels = sel_ids["totals"].get(("U", lo_total)) or []
    over_hi_sels = sel_ids["totals"].get(("O", hi_total)) or []
    under_hi_sels = sel_ids["totals"].get(("U", hi_total)) or []

    if not (home_spread_sels and away_spread_sels and
            over_lo_sels and under_lo_sels and over_hi_sels and under_hi_sels):
        if verbose:
            print(f"      integer fallback: missing alts for {total_line}")
        return None

    # Helper: price one canonical-canonical combo, return decimal or None
    def _price(sp_sels, tot_sels) -> float | None:
        for sp in sp_sels:
            sp_mnum = _market_num(sp)
            if sp_mnum not in canonical:
                continue
            for to in tot_sels:
                to_mnum = _market_num(to)
                if to_mnum not in canonical:
                    continue
                sgp = calculate_sgp(session, sp, to, verbose=verbose)
                if sgp:
                    return sgp["decimal"]
        return None

    # 8 calls: 4 combos × 2 alts
    decimals_lo = {
        "home_over":  _price(home_spread_sels, over_lo_sels),
        "home_under": _price(home_spread_sels, under_lo_sels),
        "away_over":  _price(away_spread_sels, over_lo_sels),
        "away_under": _price(away_spread_sels, under_lo_sels),
    }
    decimals_hi = {
        "home_over":  _price(home_spread_sels, over_hi_sels),
        "home_under": _price(home_spread_sels, under_hi_sels),
        "away_over":  _price(away_spread_sels, over_hi_sels),
        "away_under": _price(away_spread_sels, under_hi_sels),
    }

    if any(d is None for d in decimals_lo.values()) or any(d is None for d in decimals_hi.values()):
        if verbose:
            print(f"      integer fallback: pricing call failed for {total_line}")
        return None

    return derive_fair_probs(decimals_lo, decimals_hi)
```

Note: `_market_num` is an existing helper at the top of the file. If the import path of `integer_line_derivation` requires `mlb_sgp.` prefix (as a package), confirm `mlb_sgp/__init__.py` exists. If running scrapers directly (not as a package), use `from integer_line_derivation import ...` instead — verify by running `python scraper_draftkings_sgp.py --help` after the change.

- [ ] **Step 2: Wire the fallback into the main combo-building loop**

Find the block around line 630-655 in `scraper_draftkings_sgp.py`. The relevant existing code:

```python
        for period in ("fg", "f5"):
            spread_line = game[f"{period}_spread_line"]
            total = game[f"{period}_total_line"]
            if spread_line is None or total is None:
                continue

            sel_ids = sel_ids_per_period.get(period, {"spreads": {}, "totals": {}})
            if not sel_ids["spreads"]:
                continue

            if spread_line < 0:
                home_sign, away_sign = "N", "P"
            else:
                home_sign, away_sign = "P", "N"

            spread = abs(spread_line)

            home_spread_sels = sel_ids["spreads"].get((home_sign, spread, "1")) or []
            away_spread_sels = sel_ids["spreads"].get((away_sign, spread, "3")) or []
            over_sels = sel_ids["totals"].get(("O", total)) or []
            under_sels = sel_ids["totals"].get(("U", total)) or []

            if not home_spread_sels or not away_spread_sels or not over_sels or not under_sels:
                continue
```

Modify the trailing `continue` to attempt the fallback first. Replace the entire `if not home_spread_sels or ...` block with:

```python
            if not home_spread_sels or not away_spread_sels or not over_sels or not under_sels:
                # Existing exact-line lookup miss. Try integer-line fallback.
                canonical = sel_ids.get("canonical", set())
                fallback = try_integer_fallback_dk(
                    session, sel_ids, spread_line, total,
                    canonical, verbose=verbose,
                )
                if fallback is None:
                    continue
                # Append 4 derived rows directly to interpolated_results
                # (separate list to keep main parallel-pricing loop untouched)
                prefix = "" if period == "fg" else "F5 "
                COMBO_DISPLAY = {
                    "home_over":  "Home Spread + Over",
                    "home_under": "Home Spread + Under",
                    "away_over":  "Away Spread + Over",
                    "away_under": "Away Spread + Under",
                }
                for k, name in COMBO_DISPLAY.items():
                    interpolated_results.append((
                        gid, period, prefix + name, fallback["fair_probs"][k]
                    ))
                continue   # don't enqueue this period in the main combo_items
```

- [ ] **Step 3: Initialize `interpolated_results` and merge into output**

Find where `combo_items = []` is initialized in the same function (around line 622). Add right after it:

```python
    combo_items = []
    interpolated_results = []   # (game_id, period, combo_name, fair_prob)
```

Find where `pricing_results` is converted to rows for `mlb_sgp_odds` (search for `source` and `"draftkings_direct"` near the bottom of the same function — likely around line 720-770). After the loop that appends `_direct` rows, add:

```python
    # Append rows for integer-line interpolation fallback
    for gid, period, name, fair_prob in interpolated_results:
        decimal = 1.0 / fair_prob
        american = decimal_to_american(decimal)
        rows.append({
            "game_id": gid,
            "combo": name,
            "period": "FG" if period == "fg" else "F5",
            "bookmaker": "draftkings",
            "sgp_decimal": round(decimal, 4),
            "sgp_american": american,
            "source": "draftkings_interpolated",
        })
```

Use the same `decimal_to_american` helper the file already uses for direct rows. If it's named differently (e.g., `dec_to_am`), match the existing name.

- [ ] **Step 4: Run the DK scraper end-to-end against current slate**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/sgp-integer-derivation/mlb_sgp
source venv/bin/activate
python scraper_draftkings_sgp.py --verbose 2>&1 | tee /tmp/dk_test_run.log
```

Expected:
- Existing direct path still works (most games priced as before)
- For any integer-line games, log shows fallback attempts
- `mlb_sgp_odds` table picks up new rows with `source='draftkings_interpolated'` (verify via `python -c "import duckdb; c = duckdb.connect('Answer Keys/mlb_mm.duckdb', read_only=True); print(c.execute(\"SELECT source, COUNT(*) FROM mlb_sgp_odds GROUP BY source\").fetchall())"`)

If no integer-line games exist on the current slate, the run should complete without errors and just not write any `_interpolated` rows.

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/scraper_draftkings_sgp.py
git commit -m "mlb_sgp: DK scraper integer-line fallback via derive_fair_probs"
```

---

## Task 6: FD scraper — integer-line fallback

**Files:**
- Modify: `mlb_sgp/scraper_fanduel_sgp.py`

FD's structure is similar to DK but uses `price_combo(session, spread_market, spread_sel, total_market, total_sel)` instead of `calculate_sgp` — see line 449 in the existing file. Selection IDs are tuples `(market_id, selection_id)` keyed by `(side, line)` for spreads and `(ou, line)` for totals.

The skip site is around line 597-605:

```python
            home_spread = sel["spreads"].get(("home", home_line))
            away_spread = sel["spreads"].get(("away", away_line))
            over = sel["totals"].get(("O", total_line))
            under = sel["totals"].get(("U", total_line))

            if not (home_spread and away_spread and over and under):
                if verbose:
                    missing = []
                    if not home_spread: missing.append(f"home {home_line:+g}")
                    if not away_spread: missing.append(f"away {away_line:+g}")
                    if not over: missing.append(f"over {total_line}")
                    if not under: missing.append(f"under {total_line}")
                    print(f"  {game['away_team']} @ {game['home_team']} [{period.upper()}]: missing {missing}")
                continue
```

- [ ] **Step 1: Add `try_integer_fallback_fd` helper**

Below `price_combo` (around line 480 in `scraper_fanduel_sgp.py`), add:

```python
from .integer_line_derivation import (
    is_integer_line, derive_fair_probs, COMBO_KEYS,
)


def try_integer_fallback_fd(
    session,
    sel: dict,
    home_line: float,
    away_line: float,
    total_line: float,
    verbose: bool = False,
) -> dict | None:
    """FD version: spread keyed by ('home'|'away', signed_line), totals by ('O'|'U', line)."""
    if not is_integer_line(total_line):
        return None

    lo, hi = total_line - 0.5, total_line + 0.5

    home_spread = sel["spreads"].get(("home", home_line))
    away_spread = sel["spreads"].get(("away", away_line))
    over_lo  = sel["totals"].get(("O", lo))
    under_lo = sel["totals"].get(("U", lo))
    over_hi  = sel["totals"].get(("O", hi))
    under_hi = sel["totals"].get(("U", hi))

    if not all([home_spread, away_spread, over_lo, under_lo, over_hi, under_hi]):
        if verbose:
            print(f"      integer fallback: missing alts for {total_line}")
        return None

    def _price(sp_pair, tot_pair) -> float | None:
        result = price_combo(session, sp_pair[0], sp_pair[1], tot_pair[0], tot_pair[1], verbose=verbose)
        return result["decimal"] if result else None

    decimals_lo = {
        "home_over":  _price(home_spread, over_lo),
        "home_under": _price(home_spread, under_lo),
        "away_over":  _price(away_spread, over_lo),
        "away_under": _price(away_spread, under_lo),
    }
    decimals_hi = {
        "home_over":  _price(home_spread, over_hi),
        "home_under": _price(home_spread, under_hi),
        "away_over":  _price(away_spread, over_hi),
        "away_under": _price(away_spread, under_hi),
    }

    if any(d is None for d in decimals_lo.values()) or any(d is None for d in decimals_hi.values()):
        if verbose:
            print(f"      integer fallback: pricing failed for {total_line}")
        return None

    return derive_fair_probs(decimals_lo, decimals_hi)
```

- [ ] **Step 2: Wire fallback into the combo-building loop**

Find lines 597-605 in `scraper_fanduel_sgp.py` (the existing skip block shown above) and replace with:

```python
            if not (home_spread and away_spread and over and under):
                # Existing exact-line miss. Try integer-line fallback.
                fallback = try_integer_fallback_fd(
                    session, sel, home_line, away_line, total_line, verbose=verbose,
                )
                if fallback is not None:
                    prefix = "" if period == "fg" else "F5 "
                    COMBO_DISPLAY = {
                        "home_over":  "Home Spread + Over",
                        "home_under": "Home Spread + Under",
                        "away_over":  "Away Spread + Over",
                        "away_under": "Away Spread + Under",
                    }
                    for k, name in COMBO_DISPLAY.items():
                        interpolated_results.append((
                            gid, period, prefix + name, fallback["fair_probs"][k]
                        ))
                if verbose and fallback is None:
                    missing = []
                    if not home_spread: missing.append(f"home {home_line:+g}")
                    if not away_spread: missing.append(f"away {away_line:+g}")
                    if not over: missing.append(f"over {total_line}")
                    if not under: missing.append(f"under {total_line}")
                    print(f"  {game['away_team']} @ {game['home_team']} [{period.upper()}]: missing {missing}")
                continue
```

- [ ] **Step 3: Initialize `interpolated_results` and append to output rows**

Around line 567 in `scraper_fanduel_sgp.py` (`combo_items = []` initialization), add:

```python
    combo_items = []
    interpolated_results = []
```

In the row-writing block (search for `"source": "fanduel_direct"` — around line 660), after the existing append loop, add:

```python
    for gid, period, name, fair_prob in interpolated_results:
        decimal = 1.0 / fair_prob
        # Use the same american-conversion helper the file already uses
        american = decimal_to_american(decimal)   # adjust name if different
        rows.append({
            "game_id": gid,
            "combo": name,
            "period": "FG" if period == "fg" else "F5",
            "bookmaker": "fanduel",
            "sgp_decimal": round(decimal, 4),
            "sgp_american": american,
            "source": "fanduel_interpolated",
        })
```

- [ ] **Step 4: Run FD scraper end-to-end**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/sgp-integer-derivation/mlb_sgp
source venv/bin/activate
python scraper_fanduel_sgp.py --verbose 2>&1 | tee /tmp/fd_test_run.log
```

Expected: same as DK — direct path unchanged, integer-line games trigger fallback, `_interpolated` rows appear in `mlb_sgp_odds` if any integer-line games exist on slate.

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/scraper_fanduel_sgp.py
git commit -m "mlb_sgp: FD scraper integer-line fallback via derive_fair_probs"
```

---

## Task 7: NV scraper — integer-line fallback

**Files:**
- Modify: `mlb_sgp/scraper_novig_sgp.py`

NV uses UUID-based legs and the `/unauthenticated` endpoint. Read `mlb_sgp/scraper_novig_sgp.py` first to identify:
- The pricing function (likely `price_combo` or similar that takes UUID legs)
- The selection-ID dictionary structure (keyed by side/line just like FD)
- The skip site for missing exact lines

- [ ] **Step 1: Read the scraper to identify pricing function and skip site**

```bash
grep -n "def price\|def fetch\|continue\|sel\[" mlb_sgp/scraper_novig_sgp.py | head -50
```

Note the function name and selection-ID dict shape.

- [ ] **Step 2: Add `try_integer_fallback_nv` helper following the FD pattern**

Adapted to NV's pricing function signature and selection-ID dict format. Place above the `combo_items` building loop. Code structure mirrors `try_integer_fallback_fd` but calls NV's pricing function.

```python
from .integer_line_derivation import is_integer_line, derive_fair_probs, COMBO_KEYS


def try_integer_fallback_nv(session, sel, home_line, away_line, total_line, verbose=False):
    if not is_integer_line(total_line):
        return None
    lo, hi = total_line - 0.5, total_line + 0.5

    # Adjust the dict access to match NV's actual sel structure
    home_spread = sel["spreads"].get(("home", home_line))
    away_spread = sel["spreads"].get(("away", away_line))
    over_lo  = sel["totals"].get(("O", lo))
    under_lo = sel["totals"].get(("U", lo))
    over_hi  = sel["totals"].get(("O", hi))
    under_hi = sel["totals"].get(("U", hi))

    if not all([home_spread, away_spread, over_lo, under_lo, over_hi, under_hi]):
        if verbose:
            print(f"      integer fallback: missing alts for {total_line}")
        return None

    # Replace `nv_price_combo` with NV's actual pricing function name
    def _price(sp, tot) -> float | None:
        result = nv_price_combo(session, sp, tot, verbose=verbose)
        return result["decimal"] if result else None

    decimals_lo = {
        "home_over":  _price(home_spread, over_lo),
        "home_under": _price(home_spread, under_lo),
        "away_over":  _price(away_spread, over_lo),
        "away_under": _price(away_spread, under_lo),
    }
    decimals_hi = {
        "home_over":  _price(home_spread, over_hi),
        "home_under": _price(home_spread, under_hi),
        "away_over":  _price(away_spread, over_hi),
        "away_under": _price(away_spread, under_hi),
    }

    if any(d is None for d in decimals_lo.values()) or any(d is None for d in decimals_hi.values()):
        return None

    return derive_fair_probs(decimals_lo, decimals_hi)
```

- [ ] **Step 3: Wire fallback into combo-building loop and add `interpolated_results` accumulator**

Same pattern as FD (Task 6 Step 2-3): replace the `if not exact_line: continue` block with a fallback attempt, accumulate into `interpolated_results`, write rows with `source='novig_interpolated'`.

- [ ] **Step 4: Run NV scraper end-to-end**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/sgp-integer-derivation/mlb_sgp
source venv/bin/activate
python scraper_novig_sgp.py --verbose 2>&1 | tee /tmp/nv_test_run.log
```

Expected: same shape as DK/FD test runs.

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/scraper_novig_sgp.py
git commit -m "mlb_sgp: NV scraper integer-line fallback via derive_fair_probs"
```

---

## Task 8: PX scraper — integer-line fallback

**Files:**
- Modify: `mlb_sgp/scraper_prophetx_sgp.py`

PX uses an RFQ endpoint and `competitorId` namespace. Same pattern as NV:

- [ ] **Step 1: Read scraper to identify pricing function and skip site**

```bash
grep -n "def price\|def fetch\|continue\|sel\[" mlb_sgp/scraper_prophetx_sgp.py | head -50
```

- [ ] **Step 2: Add `try_integer_fallback_px` helper**

Mirror NV's helper, adjusting for PX's pricing-function signature.

- [ ] **Step 3: Wire into combo-building loop**

Same pattern: replace the skip block with fallback attempt, accumulate, write rows with `source='prophetx_interpolated'`.

- [ ] **Step 4: Run PX scraper end-to-end**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/sgp-integer-derivation/mlb_sgp
source venv/bin/activate
python scraper_prophetx_sgp.py --verbose 2>&1 | tee /tmp/px_test_run.log
```

- [ ] **Step 5: Commit**

```bash
git add mlb_sgp/scraper_prophetx_sgp.py
git commit -m "mlb_sgp: PX scraper integer-line fallback via derive_fair_probs"
```

---

## Task 9: R blender — pick up `_interpolated` source rows

**Files:**
- Modify: `Answer Keys/mlb_correlated_parlay.R`

The R blender currently filters `mlb_sgp_odds` rows by `source IN (<book>_direct, ...)`. We widen the filter so `<book>_interpolated` rows are also picked up and treated identically by the per-game devig + blending logic.

- [ ] **Step 1: Find the WHERE clause**

```bash
grep -n "source IN\|source LIKE\|draftkings_direct\|fanduel_direct" "Answer Keys/mlb_correlated_parlay.R"
```

There may be multiple sites (one per book in a per-book vig measurement, plus a unified read). Identify all of them.

- [ ] **Step 2: Replace `source IN (...)` patterns with `source LIKE` patterns**

For each site, change:

```r
WHERE source = 'draftkings_direct'
```

to:

```r
WHERE source LIKE 'draftkings_%'
```

And likewise for `fanduel_direct` → `fanduel_%`, `novig_direct` → `novig_%`, `prophetx_direct` → `prophetx_%`.

If there's a unified read, change:

```r
WHERE source IN ('draftkings_direct', 'fanduel_direct', 'novig_direct', 'prophetx_direct')
```

to:

```r
WHERE source LIKE 'draftkings_%'
   OR source LIKE 'fanduel_%'
   OR source LIKE 'novig_%'
   OR source LIKE 'prophetx_%'
```

- [ ] **Step 3: Run the full pipeline end-to-end**

```bash
cd "/Users/callancapitolo/NFLWork/.worktrees/sgp-integer-derivation/Answer Keys"
Rscript run.R --sport mlb 2>&1 | tee /tmp/r_pipeline_test.log
```

Expected:
- Pipeline completes without errors
- `mlb_parlay_opportunities` table populated as before
- For any games that had `_interpolated` source rows, the corresponding `dk_fair_prob` / `fd_fair_prob` / etc. columns are non-null where they previously would have been null

Verify:

```python
python3 -c "
import duckdb
c = duckdb.connect('Answer Keys/mlb_mm.duckdb', read_only=True)
print(c.execute('''
    SELECT bookmaker, source, COUNT(*) AS n
    FROM mlb_sgp_odds
    GROUP BY bookmaker, source
    ORDER BY bookmaker, source
''').fetchdf().to_string())
"
```

- [ ] **Step 4: Commit**

```bash
git add "Answer Keys/mlb_correlated_parlay.R"
git commit -m "mlb_correlated_parlay: pick up *_interpolated source rows alongside *_direct"
```

---

## Task 10: Integration test on real slate

**Files:**
- Create: `mlb_sgp/tests/test_integration_integer_line.py`

A one-off test (not added to cron) that runs the full new path end-to-end against the current slate's `mlb_parlay_lines` and verifies `_interpolated` rows are written when an integer line is present.

- [ ] **Step 1: Write the integration test script**

```python
# mlb_sgp/tests/test_integration_integer_line.py
"""
Manual integration test for the integer-line derivation feature.

Run this after a slate refresh to verify the new fallback path triggers
on real data. Not part of the regular pytest suite — run explicitly:

    python -m mlb_sgp.tests.test_integration_integer_line

Exits 0 if a slate has integer-line games AND _interpolated rows were written.
Exits 1 if integer-line games exist but no _interpolated rows appeared
(indicates the fallback isn't triggering or all bounds checks are failing).
Exits 2 if no integer-line games on slate (test inconclusive).
"""
import sys
import duckdb

DB_PATH = "Answer Keys/mlb_mm.duckdb"


def main() -> int:
    con = duckdb.connect(DB_PATH, read_only=True)

    # Are there any integer-line games on the current slate?
    integer_games = con.execute("""
        SELECT game_id, fg_total, f5_total
        FROM mlb_parlay_lines
        WHERE fg_total = ROUND(fg_total) OR f5_total = ROUND(f5_total)
    """).fetchall()

    if not integer_games:
        print("No integer-line games on current slate — test inconclusive.")
        return 2

    print(f"Found {len(integer_games)} integer-line game(s) on slate.")

    # Did any book write _interpolated rows for those games?
    interpolated = con.execute("""
        SELECT bookmaker, source, COUNT(*) AS n
        FROM mlb_sgp_odds
        WHERE source LIKE '%_interpolated'
        GROUP BY bookmaker, source
        ORDER BY bookmaker
    """).fetchdf()

    if interpolated.empty:
        print("FAIL: integer-line games exist but no _interpolated rows in mlb_sgp_odds.")
        print("Investigation: check scraper logs for bounds-check WARN lines.")
        return 1

    print("PASS: _interpolated rows present.")
    print(interpolated.to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 2: Run the test against the live slate**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/sgp-integer-derivation
source mlb_sgp/venv/bin/activate
python -m mlb_sgp.tests.test_integration_integer_line
echo "Exit code: $?"
```

Expected: exit 0 (PASS) on a slate with integer-line games and properly-functioning fallback. If exit 2 (no integer-line games), wait for next slate or skip this gate.

- [ ] **Step 3: Commit**

```bash
git add mlb_sgp/tests/test_integration_integer_line.py
git commit -m "mlb_sgp: add manual integration test for integer-line fallback"
```

---

## Task 11: README documentation

**Files:**
- Modify: `mlb_sgp/README.md`

- [ ] **Step 1: Read the current README to find the section to update**

```bash
grep -n "Exact line matching\|exact line\|F5 totals at integer" mlb_sgp/README.md
```

The relevant sections are around the "SGP Scraping Playbook" → "Exact line matching" subheader (around line 191-193) and the FD-specific note about F5 integer totals (around line 285).

- [ ] **Step 2: Replace the "Exact line matching only" section**

Find:

```markdown
### Exact line matching

The scraper must match Wagerzon's exact spread and total lines. If the book doesn't have the precise line (e.g., Wagerzon has total 8.5 but the book only offers 8.0), skip the game entirely. No rounding, no closest-line fallback. The parlay fair odds are calibrated to Wagerzon's specific lines — a different line is a different bet.
```

Replace with:

```markdown
### Line matching

The scraper first attempts exact-line matching: Wagerzon's exact spread and total lines must exist in the book's selection-ID dictionary, and that pair gets priced via `calculateBets` / `implyBets` / equivalent.

**Integer-line fallback.** When Wagerzon posts an SGP at an integer total line (e.g., FG Over 8, F5 Over 4) that the book doesn't quote directly, the scraper falls back to interpolating from the two adjacent half-point alts (X-0.5 and X+0.5). The integer-line derivation is implemented in `mlb_sgp/integer_line_derivation.py`. See "Integer-Line Derivation" below.

**No-fallback behavior.** If the book doesn't have either adjacent half-point alt, or any bounds check fails on the derived prices, the game's period is skipped for that book — same graceful degradation as today's "WZ line not found" path. The R blender treats the book as missing for that game.

### Integer-Line Derivation

When an integer total line `X` is missing, the scraper issues 8 SGP pricing calls (4 combos × 2 alts at `X − 0.5` and `X + 0.5`), per-alt-devigs each set, and derives joint fair probabilities at the integer line:

```
Δ_total      = (devig_HomeOver_lo − devig_HomeOver_hi)
             + (devig_AwayOver_lo − devig_AwayOver_hi)
             = P(T = X)   (the joint marginal push mass)

fair_prob_X  = devig_hi_combo / (1 − Δ_total)
```

Interpretation: the conditional joint probability "this combo wins outright, given the bet doesn't push." Books refund pushed legs regardless of spread side, so the conditioning matches real-world settlement. Pure-joint computation — no singles totals market accessed.

**Underlying principle:** Breeden & Litzenberger (1978) showed in options markets that adjacent strike prices implicitly carry the probability mass of the in-between outcome (the second derivative of call prices w.r.t. strike). This is the discrete sports analog: the difference between adjacent half-point alt SGPs gives us the implied joint integer-mass at `X`. Half-point pricing calculators (Sportsbook Review, Bookmakers Review, Unabated) apply the same skeleton in singles markets; we extend it to two-leg correlated parlays.

The derived rows are written with `source = '<book>_interpolated'` so post-hoc analysis can compare realized P&L on interpolated vs directly-priced rows.

**Bounds checks** (config constants in `integer_line_derivation.py`):

| # | Check | Threshold |
|---|---|---|
| 1 | Per-alt vig sum | `[1.05, 1.30]` |
| 2 | Push-mass cross-consistency | <10% relative diff per side |
| 3 | `Δ_total` plausibility | `[0.03, 0.18]` |
| 4 | Sum of 4 derived fair_probs | `[0.97, 1.03]` |
| 5 | Per-combo bounds | `(0, 1)` strict |

On any failure: the game's period is skipped for that book, structured WARN logged with inputs and the violated check.
```

- [ ] **Step 3: Update the FD-specific note**

Find:

```markdown
- **F5 totals at integer values (e.g., 5.0) may not exist.** FD's F5 alt totals jump in 1.0 increments (2.5, 3.5, 4.5, 5.5...). If Wagerzon has F5 total 5.0, FD won't have it and the game is correctly skipped.
```

Replace with:

```markdown
- **F5 totals at integer values (e.g., 5.0) trigger interpolation.** FD's F5 alt totals jump in 1.0 increments (2.5, 3.5, 4.5, 5.5...). When Wagerzon has F5 total 5.0, FD's exact-line lookup misses → the integer-line fallback kicks in, using FD's F5 alts at 4.5 and 5.5. See "Integer-Line Derivation" above.
```

- [ ] **Step 4: Commit**

```bash
git add mlb_sgp/README.md
git commit -m "mlb_sgp: README documentation for integer-line derivation"
```

---

## Task 12: Pre-merge review

**Files:**
- No file changes — review only.

Per CLAUDE.md, an executive engineer review of the full diff is required before merge.

- [ ] **Step 1: Review the full diff**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/sgp-integer-derivation
git diff main..HEAD --stat
git diff main..HEAD
```

- [ ] **Step 2: Walk through the CLAUDE.md review checklist**

For each item, document `OK` or `ISSUE` with a note:

- **Data integrity:** Bounds-check failures result in zero rows written for that game/period (verified — `derive_fair_probs` returns `None` and the wiring uses `if fallback is not None:` before appending). No duplicate writes (each scraper writes its own `_interpolated` rows once per game/period).
- **Resource safety:** No new DB connections opened. Existing scraper retry logic reused via the same pricing functions.
- **Edge cases:**
  - Off-season / no slate: integer-line check is per-game; if no games → no fallback triggers → zero impact
  - First run after merge: no `_interpolated` rows yet → R `LIKE` pattern handles cleanly (returns 0 rows for the missing pattern)
  - Half-point line (today's normal case): `is_integer_line` returns False → skip fallback → existing path unchanged
- **Dead code:** No unused helpers. `COMBO_KEYS` exported for use by all 4 scrapers.
- **Log/disk hygiene:** WARN logs go to existing `mlb_sgp/logs/<scraper>.log` files (already rotated). No new log files.
- **Security:** Structured WARN log line includes only numerical inputs and the violated threshold — no API keys, no scrape session tokens.

- [ ] **Step 3: Run all tests one more time**

```bash
cd /Users/callancapitolo/NFLWork/.worktrees/sgp-integer-derivation
source mlb_sgp/venv/bin/activate
pytest mlb_sgp/tests/test_integer_line_derivation.py -v
python -m mlb_sgp.tests.test_integration_integer_line   # if integer games on slate
```

Both must pass / be inconclusive (no integer games is OK; failures are not).

- [ ] **Step 4: Surface findings to user, ask explicit merge approval**

> "All tasks complete. Pre-merge review:
> - All 5 CLAUDE.md checklist items: OK
> - All unit tests: passing
> - Integration test: <PASS / inconclusive>
>
> Diff summary: <X files changed, Y insertions, Z deletions>
>
> Ready to merge `feature/sgp-integer-line-derivation` to `main` and clean up the worktree?"

**Wait for explicit "yes" before merging. Do not merge without approval.**

- [ ] **Step 5: After approval — merge and clean up**

```bash
cd /Users/callancapitolo/NFLWork
git checkout main
git merge --no-ff feature/sgp-integer-line-derivation -m "feat(mlb-sgp): integer-line derivation for SGP comparisons"
git worktree remove .worktrees/sgp-integer-derivation
git branch -d feature/sgp-integer-line-derivation
```

- [ ] **Step 6: Verify post-merge**

```bash
cd /Users/callancapitolo/NFLWork
git log --oneline -10   # confirm merge commit on main
git worktree list       # confirm worktree removed
```

---

## Self-Review Notes

After writing the plan, I checked it against the spec:

- **Spec §3 (math):** Tasks 1-4 (helper, devig, bounds, orchestrator) cover the math end-to-end. ✓
- **Spec §4 (architecture):** Tasks 5-9 cover the 4 scraper integrations + R blender. ✓
- **Spec §5 (algorithm):** Steps 1-12 in §5 of the spec map to Tasks 5-8 (per-scraper fallback flow). ✓
- **Spec §6 (bounds checks):** Task 3 implements all 5; Task 4 wires them into the orchestrator. ✓
- **Spec §7 (data flow):** Task 9 ensures the R blender picks up new rows. ✓
- **Spec §8 (testing):** Task 1-4 cover unit tests; Task 10 covers integration test. ✓
- **Spec §9 (documentation):** Task 11 covers README. ✓
- **Spec §10 (version control):** Tasks set up worktree + commits per the spec's commit structure; Task 12 covers pre-merge review. ✓
- **Spec §11 (open questions):** Listed as future work in the spec; not in scope for this plan. ✓

Type/name consistency:
- `COMBO_KEYS = ("home_over", "home_under", "away_over", "away_under")` used consistently across the shared module and all 4 scrapers.
- `derive_fair_probs(decimals_lo, decimals_hi)` signature is the same in all 4 scraper integration tasks.
- `source = '<book>_interpolated'` naming matches across scrapers and the R `LIKE` pattern in Task 9.

No placeholders. Each step has exact code or exact commands.
