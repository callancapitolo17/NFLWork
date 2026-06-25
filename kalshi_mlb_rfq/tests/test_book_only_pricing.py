import importlib
import pytest


def test_book_only_fair_median():
    from kalshi_mlb_rfq.main import _book_only_fair
    assert _book_only_fair({"dk": 0.30, "fd": 0.34}) == pytest.approx(0.32)
    assert _book_only_fair({"dk": 0.30, "fd": 0.34, "px": 0.40}) == pytest.approx(0.34)
    assert _book_only_fair({"dk": 0.30, "fd": None}) == pytest.approx(0.30)
    assert _book_only_fair({}) is None
    assert _book_only_fair({"dk": None}) is None
