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
