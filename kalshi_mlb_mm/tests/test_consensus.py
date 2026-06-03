"""Tests for the v1 book-consensus-band gate (main._consensus_filter).

The filter takes a dict of {book: devigged_fair} and returns the subset that
agrees within ±BOOK_CONSENSUS_BAND of the median, or an empty dict if too few
books survive. Mirrors the MLB answer-key dashboard's consensus-band pattern.
"""
from kalshi_mlb_mm import main, config


def test_consensus_three_books_agree_pass():
    """All three books within ±BOOK_CONSENSUS_BAND of the median → all kept."""
    out = main._consensus_filter({"dk": 0.42, "fd": 0.43, "px": 0.43})
    assert len(out) == 3
    assert set(out) == {"dk", "fd", "px"}


def test_consensus_one_outlier_dropped():
    """Median of 4 books is ~0.43; 0.55 is far outside ±0.02 → dropped, the
    three agreeing books remain (still >= MIN_AGREEING_BOOKS)."""
    out = main._consensus_filter({"dk": 0.42, "fd": 0.43, "px": 0.43, "nv": 0.55})
    assert "nv" not in out
    assert len(out) == 3


def test_consensus_only_two_in_band_skip():
    """Median = 0.43, agreeing = {fd, px} (2) which is < MIN_AGREEING_BOOKS (3)
    → return empty dict (no quote)."""
    out = main._consensus_filter({"dk": 0.30, "fd": 0.43, "px": 0.43, "nv": 0.60})
    assert out == {}


def test_consensus_too_few_books_skip():
    """Below MIN_AGREEING_BOOKS at the outer count check → empty dict."""
    out = main._consensus_filter({"dk": 0.43, "fd": 0.43})
    assert out == {}


def test_consensus_band_just_inside():
    """A book just inside the BAND is kept; just outside is dropped."""
    band = config.BOOK_CONSENSUS_BAND  # default 0.02
    # median of {0.42, 0.43, 0.44} = 0.43; all within band → all kept.
    out = main._consensus_filter({"dk": 0.42, "fd": 0.43, "px": 0.44})
    assert len(out) == 3
    # Now bump one book just outside the band (median moves but only slightly).
    # 0.43 vs 0.50 is 0.07 > 0.02 → 0.50 dropped, but only 2 remain → empty.
    out = main._consensus_filter({"dk": 0.42, "fd": 0.43, "px": 0.50})
    assert out == {}, f"only 2 in band < MIN(3) → empty, got {out}"
    assert band == 0.02  # sanity: default is what we assumed
