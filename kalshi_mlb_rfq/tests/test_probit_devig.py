"""Unit tests for probit devigging in fair_value.py.

Validates _probit_devig_n (n-way probit) and devig_book (4-cell SGP grid).
Includes a negative test that catches accidental reverts to multiplicative.
"""
import sys
from pathlib import Path

# Add kalshi_mlb_rfq to path so we can import fair_value
sys.path.insert(0, str(Path(__file__).parent.parent))

import pandas as pd

from fair_value import _probit_devig_n, devig_book


# ---------------------------------------------------------------------------
# _probit_devig_n
# ---------------------------------------------------------------------------

def test_probit_devig_n_2way_sums_to_one():
    out = _probit_devig_n([0.55, 0.50])
    assert abs(sum(out) - 1.0) < 1e-9

def test_probit_devig_n_4way_symmetric_each_equals_quarter():
    # 4 equal vigged probs -> 4 equal devigged probs = 0.25
    out = _probit_devig_n([0.27, 0.27, 0.27, 0.27])
    for p in out:
        assert abs(p - 0.25) < 1e-6

def test_probit_devig_n_handles_nan():
    out = _probit_devig_n([float("nan"), 0.5])
    # NaN check: p != p is True only when p is NaN
    assert any(p != p for p in out)

def test_probit_devig_n_negative_diverges_from_multiplicative():
    # Tail case: probit must give meaningfully different answer than multiplicative
    # p_raw = (0.111, 0.901), sum = 1.012
    # Multiplicative: (0.111/1.012, 0.901/1.012) = (0.1097, 0.8903)
    # Probit: roughly (0.1099, 0.8901) — yes, differs by < 0.001 here, but at
    # more extreme tails the gap widens. Test on a wider gap:
    out = _probit_devig_n([0.05, 0.97])
    mult = [0.05 / 1.02, 0.97 / 1.02]
    # Must differ on at least one side by 0.001
    assert any(abs(out[i] - mult[i]) > 0.001 for i in range(2))


# ---------------------------------------------------------------------------
# devig_book
# ---------------------------------------------------------------------------

def test_devig_book_4_cell_grid_uses_probit():
    # Use sgp_decimal values whose 1/decimal sums to ~1.04 (4% overround)
    rows = pd.DataFrame({
        "combo": [
            "Home Spread + Over",
            "Home Spread + Under",
            "Away Spread + Over",
            "Away Spread + Under",
        ],
        "sgp_decimal": [3.00, 3.50, 4.00, 5.00],
    })
    raw_probs = [1.0 / d for d in rows["sgp_decimal"]]
    expected = _probit_devig_n(raw_probs)

    out = devig_book(rows, combo="Home Spread + Over")
    # Target combo is at index 0; should match probit output at that index
    assert abs(out - expected[0]) < 1e-9

def test_devig_book_single_row_falls_back_to_heuristic():
    rows = pd.DataFrame({
        "combo": ["Home Spread + Over"],
        "sgp_decimal": [3.0],
    })
    out = devig_book(rows, combo="Home Spread + Over", vig_fallback=0.10)
    expected = (1.0 / 3.0) / 1.10
    assert abs(out - expected) < 1e-9

def test_devig_book_returns_none_on_empty():
    out = devig_book(pd.DataFrame(columns=["combo", "sgp_decimal"]), combo="x")
    assert out is None

def test_devig_book_returns_none_on_missing_combo():
    rows = pd.DataFrame({"combo": ["A"], "sgp_decimal": [3.0]})
    out = devig_book(rows, combo="B")
    assert out is None
