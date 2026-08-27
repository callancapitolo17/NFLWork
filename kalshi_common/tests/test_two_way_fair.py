"""Two-way rung devig with its feed-integrity gate (issue #96).

The leg surface devigs one book's two sides of one posted line, thousands of
times per cycle, unsupervised. The gate exists because the failure mode is not
a crash — ``devig_two_way`` clips and solves ANY pair of decimals, so a
poisoned input comes back as a confident wrong number that a quote then rests
on. Every test below is a shape of poison the gate must refuse.
"""
import pytest

from kalshi_common.fair_value import two_way_fair


class TestPasses:
    def test_symmetric_pair_is_a_coin_flip(self):
        fair, reason = two_way_fair(1.91, 1.91)
        assert reason is None
        assert fair == pytest.approx(0.5)

    def test_sides_sum_to_one(self):
        yes, _ = two_way_fair(1.60, 2.45)
        no, _ = two_way_fair(2.45, 1.60)
        assert yes + no == pytest.approx(1.0)

    def test_favorite_prices_above_a_half(self):
        fair, reason = two_way_fair(1.40, 3.10)
        assert reason is None
        assert fair > 0.5

    def test_low_vig_book_at_the_band_floor_passes(self):
        # Novig-style near-zero vig: overround 1.0050 is exactly band_min.
        fair, reason = two_way_fair(1.9900496, 1.9900496)
        assert reason is None
        assert fair == pytest.approx(0.5)


class TestRejects:
    def test_crossed_pair_is_refused_not_devigged(self):
        # Implied sum 0.952 — arithmetically not a two-way market. This is
        # what a mid-refresh straddle looks like, and it is the one shape that
        # must never be devigged regardless of the band.
        assert two_way_fair(2.10, 2.10) == (None, "crossed")

    def test_blown_vig_is_refused(self):
        assert two_way_fair(1.20, 1.30)[1] == "overround"

    def test_vig_just_over_the_ceiling_is_refused(self):
        # 20.1% overround — one tick outside the envelope.
        fair, reason = two_way_fair(1.665, 1.665)
        assert fair is None and reason == "overround"

    def test_missing_side_is_a_bad_price_not_a_haircut(self):
        # A one-sided line must decline. The Route-B vig haircut is calibrated
        # for a combo's compounded margin; applied to a lone leg it invents a
        # number, and #95 showed 4-5 books per leg so declining costs quorum
        # nothing.
        assert two_way_fair(None, 1.91)[1] == "bad_price"

    def test_non_numeric_price_is_a_bad_price(self):
        assert two_way_fair("evens", 1.91)[1] == "bad_price"

    def test_decimal_of_one_is_a_bad_price(self):
        # 1.0 means "pays nothing"; 1/1.0 = 1.0 would read as certainty.
        assert two_way_fair(1.0, 1.91)[1] == "bad_price"

    def test_infinite_price_is_a_bad_price(self):
        assert two_way_fair(float("inf"), 1.91)[1] == "bad_price"


def test_band_is_caller_overridable():
    # The band is a knob (SURFACE_OVERROUND_MIN/MAX) so the surface can tune
    # it per data without editing shared math.
    assert two_way_fair(1.20, 1.30, band_min=1.0, band_max=2.0)[1] is None
