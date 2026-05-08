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
