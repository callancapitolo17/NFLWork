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
