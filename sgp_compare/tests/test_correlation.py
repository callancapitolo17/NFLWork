import math

import pytest

from sgp_compare.correlation import (correlation_report, leg_chance, partition_flips, probit_devig)


def test_probit_devig_two_way_sums_to_one_and_keeps_order():
    fair = probit_devig([0.55, 0.50])
    assert math.isclose(sum(fair), 1.0, abs_tol=1e-9)
    assert fair[0] > fair[1]


def test_probit_devig_many_way_sums_to_one():
    assert math.isclose(sum(probit_devig([0.30, 0.30, 0.25, 0.25])), 1.0, abs_tol=1e-9)


def test_partition_flips_puts_the_parlay_first():
    flips = partition_flips(2)
    assert flips[0] == (False, False)
    assert len(set(flips)) == 4


def test_with_vig_factor_matches_the_recon_quote():
    # ProphetX 2026-11-18: HOU -12.5 (+100) + Collins o70.5 (-112) paid +217.30.
    spread = leg_chance("spread", 2.0, 1 + 100 / 105)
    collins = leg_chance("collins", 1 + 100 / 112, 1 + 100 / 108)
    report = correlation_report(3.1729841909, [spread, collins], None)
    assert report.correlation_raw == pytest.approx(0.193, abs=0.001)
    assert report.correlation_fair is None
    assert "other side to flip" in report.fair_unavailable_reason


def test_vig_removed_factor_from_a_full_partition():
    # Independent 50/50 legs; the venue thinks they hit together 30% of the time
    # and charges 5% across the four combinations.
    legs = [leg_chance("a", 1.95, 1.95), leg_chance("b", 1.95, 1.95)]
    true_cells = [0.30, 0.20, 0.20, 0.30]
    cells = [1 / (p * 1.05) for p in true_cells]
    report = correlation_report(cells[0], legs, cells)
    assert report.independent_fair_chance == pytest.approx(0.25)
    assert report.parlay_fair_chance == pytest.approx(0.30, abs=0.005)
    assert report.correlation_fair == pytest.approx(0.20, abs=0.02)
    assert report.partition_overround == pytest.approx(0.05)


def test_a_declined_combination_gives_no_vig_removed_factor():
    legs = [leg_chance("a", 1.95, 1.95), leg_chance("b", 1.95, 1.95)]
    report = correlation_report(3.0, legs, [3.0, 4.0, None, 3.5])
    assert report.correlation_fair is None
    assert "1 of the 4" in report.fair_unavailable_reason
    assert report.correlation_raw is not None


def test_crossed_or_too_wide_singles_get_no_fair_chance():
    assert leg_chance("crossed", 2.1, 2.1).fair_chance is None          # +110 / +110 adds to 95%
    assert leg_chance("wide", 1 + 100 / 400, 1.5).fair_chance is None   # -400 / -200 adds to 147%
    assert leg_chance("one_sided", 2.0, None).fair_chance is None
    assert leg_chance("ok", 1 + 100 / 112, 1 + 100 / 108).fair_chance == pytest.approx(0.5045, abs=0.0005)
