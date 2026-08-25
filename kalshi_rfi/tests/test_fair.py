"""Unit tests for the book-consensus dispersion gate (no network)."""
import statistics

import pytest
from scipy.stats import norm

from kalshi_rfi.fair import consensus


class TestConsensusGate:
    def test_one_book_declined(self):
        assert consensus({"draftkings": 0.44}) is None

    def test_none_fairs_do_not_count(self):
        assert consensus({"draftkings": 0.44, "fanduel": None}) is None

    def test_two_agreeing_books_pass(self):
        result = consensus({"draftkings": 0.44, "fanduel": 0.45})
        assert result is not None
        fair, sigma_z = result
        assert fair == statistics.median([0.44, 0.45])
        assert sigma_z == pytest.approx(
            statistics.stdev([norm.ppf(0.44), norm.ppf(0.45)]))

    def test_dispersed_books_declined(self):
        # 0.35 vs 0.55 is ~0.36 sigma_z — way past the 0.07 gate
        assert consensus({"draftkings": 0.35, "fanduel": 0.55}) is None

    def test_median_of_three(self):
        result = consensus({"a": 0.43, "b": 0.44, "c": 0.45})
        assert result is not None
        assert result[0] == 0.44

    def test_no_outlier_removal_dissenter_kills_quote(self):
        # Two agree, one loudly dissents -> whole game declined (the
        # dissenter may be the informed book).
        assert consensus({"a": 0.44, "b": 0.45, "c": 0.60}) is None


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
