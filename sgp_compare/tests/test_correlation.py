import math

import pytest

from sgp_compare.correlation import LegChance, correlation_report, leg_chance, partition_flips


def leg_with(chance: float, vig: float = 0.04) -> LegChance:
    """A two-sided single at `chance`, each side marked up by `vig`."""
    return leg_chance("leg", 1 / (chance * (1 + vig)), 1 / ((1 - chance) * (1 + vig)))


def partition_prices(true_cell_chances: list[float], markup: float) -> list[float]:
    return [1 / (p * (1 + markup)) for p in true_cell_chances]


def independent_cells(leg_chances: list[float]) -> list[float]:
    cells = []
    for flips in partition_flips(len(leg_chances)):
        cells.append(math.prod((1 - p) if flipped else p for p, flipped in zip(leg_chances, flips)))
    return cells


def test_partition_flips_puts_the_parlay_first_and_lists_every_combination():
    assert partition_flips(2) == [(False, False), (False, True), (True, False), (True, True)]
    assert len(set(partition_flips(3))) == 8


def test_with_vig_factor_matches_the_recon_quote():
    # ProphetX 2026-11-18: HOU -12.5 (+100 / -105) + Collins o70.5 (-112 / -108) paid +217.30.
    legs = [leg_chance("spread", 2.0, 1 + 100 / 105), leg_chance("collins", 1 + 100 / 112, 1 + 100 / 108)]
    report = correlation_report(3.1729841909, legs, None, "no partition in this test")
    assert report.independent_raw_chance == pytest.approx(0.5 * 112 / 212)
    assert report.correlation_raw == pytest.approx(0.1931, abs=0.0001)
    assert report.correlation_fair is None and report.parlay_fair_chance is None
    assert report.fair_unavailable_reason == "no partition in this test"


@pytest.mark.parametrize("leg_chances", [[0.2, 0.2], [0.15, 0.25], [0.8, 0.8], [0.2, 0.25, 0.3]])
@pytest.mark.parametrize("markup", [0.05, 0.15, 0.30])
def test_independent_legs_read_zero_whatever_their_balance_or_the_parlay_vig(leg_chances, markup):
    cells = partition_prices(independent_cells(leg_chances), markup)
    report = correlation_report(cells[0], [leg_with(p) for p in leg_chances], cells, None)
    assert report.correlation_fair == pytest.approx(0.0, abs=1e-9)
    assert report.partition_overround == pytest.approx(markup)
    assert report.parlay_fair_chance == pytest.approx(math.prod(leg_chances))


def test_correlated_legs_read_their_true_factor():
    # Two 50% legs that hit together 30% of the time: factor 0.30 / 0.25 - 1 = +20%.
    cells = partition_prices([0.30, 0.20, 0.20, 0.30], 0.10)
    report = correlation_report(cells[0], [leg_with(0.5), leg_with(0.5)], cells, None)
    assert report.correlation_fair == pytest.approx(0.20)
    assert report.partition_independent_chance == pytest.approx(0.25)
    assert report.independent_fair_chance == pytest.approx(0.25)


def test_crossed_partition_gives_no_vig_removed_factor():
    report = correlation_report(4.5, [leg_with(0.5), leg_with(0.5)], [4.5, 4.5, 4.5, 4.5], None)
    assert report.correlation_fair is None
    assert "crossed" in report.fair_unavailable_reason
    assert report.partition_overround == pytest.approx(4 / 4.5 - 1)


def test_wrong_number_of_partition_prices_is_an_error():
    with pytest.raises(ValueError, match="expected 4 partition prices"):
        correlation_report(3.0, [leg_with(0.5), leg_with(0.5)], [3.0, 3.0], None)


def test_missing_partition_needs_a_reason():
    with pytest.raises(ValueError, match="reason"):
        correlation_report(3.0, [leg_with(0.5), leg_with(0.5)], None, None)


def test_a_leg_nobody_offers_alone_names_itself_as_the_reason():
    unpriced = leg_chance("x", None, 2.0)
    assert (unpriced.raw_chance, unpriced.fair_chance) == (None, None)
    assert unpriced.fair_unavailable_reason == "nobody is offering this side alone"
    report = correlation_report(3.0, [unpriced, leg_with(0.5)], None, "skip")
    assert report.correlation_raw is None and report.independent_raw_chance is None
    assert report.raw_unavailable_reason.startswith("1 leg(s) not offered alone")


def test_singles_that_are_crossed_too_wide_or_one_sided_get_no_fair_chance():
    assert leg_chance("c", 2.1, 2.1).fair_unavailable_reason == "its two sides add to under 100%"
    assert leg_chance("w", 1 + 100 / 400, 1.5).fair_unavailable_reason == "its two sides are more than 25% apart"
    assert leg_chance("o", 2.0, None).fair_unavailable_reason == "nobody is offering its other side"
    ok = leg_chance("ok", 1 + 100 / 112, 1 + 100 / 108)
    assert ok.fair_chance == pytest.approx(0.5045, abs=0.0005) and ok.fair_unavailable_reason is None


def test_the_live_smoke_quotes_reproduce():
    # Live ProphetX 2026-11-19 on Cal's Mac: HOU -12.5 + Collins o70.5, all four combinations.
    legs = [leg_chance("s", 2.0, 1 + 100 / 105), leg_chance("c", 1 + 100 / 112, 1 + 100 / 108)]
    report = correlation_report(3.1730, legs, [3.1730, 3.5831, 3.5680, 3.1841], None)
    assert report.correlation_raw == pytest.approx(0.1931, abs=0.0001)
    assert report.correlation_fair == pytest.approx(0.0587, abs=0.0001)
    assert report.partition_overround == pytest.approx(0.1886, abs=0.0001)
    assert report.fair_unavailable_reason is None


def test_a_partition_past_its_envelope_gives_no_factor():
    report = correlation_report(3.0, [leg_with(0.5), leg_with(0.5)], [3.0, 1.2, 1.2, 3.0], None)
    assert report.correlation_fair is None
    assert "so a price is stale or wrong" in report.fair_unavailable_reason
    assert report.partition_overround == pytest.approx(2 / 3 + 2 / 1.2 - 1)
