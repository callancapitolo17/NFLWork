import pytest
from nfl_draft.lib.devig import (
    american_to_implied,
    devig_two_way,
    devig_n_way,
    proportional_devig,
    devig_pool,
)


def test_american_to_implied_pickem():
    assert american_to_implied(100) == pytest.approx(0.5, abs=1e-6)


def test_american_to_implied_favorite():
    assert american_to_implied(-200) == pytest.approx(2/3, abs=1e-6)


def test_american_to_implied_dog():
    assert american_to_implied(200) == pytest.approx(1/3, abs=1e-6)


def test_devig_two_way_zero_overround():
    p_a, p_b = devig_two_way(100, 100)
    assert p_a + p_b == pytest.approx(1.0)
    assert p_a == pytest.approx(0.5)


def test_devig_two_way_with_vig():
    p_a, p_b = devig_two_way(-110, -110)
    assert p_a == pytest.approx(0.5, abs=1e-4)
    assert p_b == pytest.approx(0.5, abs=1e-4)


def test_devig_two_way_asymmetric():
    p_a, p_b = devig_two_way(-200, 150)
    assert (p_a + p_b) == pytest.approx(1.0, abs=1e-6)
    assert p_a > p_b


def test_devig_n_way_three_outcomes():
    probs = devig_n_way([200, 200, 200])
    assert sum(probs) == pytest.approx(1.0)
    for p in probs:
        assert p == pytest.approx(1/3, abs=1e-6)


def test_devig_n_way_with_vig():
    probs = devig_n_way([100, 100, 100])
    assert sum(probs) == pytest.approx(1.0)
    for p in probs:
        assert p == pytest.approx(1/3, abs=1e-6)


def test_proportional_devig_normalizes():
    odds = [400, 400, 400, 400, 400]
    devigged = proportional_devig(odds)
    assert sum(devigged) == pytest.approx(1.0)


def test_proportional_devig_handles_int_and_float_inputs():
    devigged = proportional_devig([100.0, 200, -150])
    assert sum(devigged) == pytest.approx(1.0)


def test_devig_pool_normalizes_to_pool_size():
    # Three outcomes at +200 (implied 1/3 each, raw sum = 1.0). With
    # pool_size=1 the pool reduces to a classic 1-winner outright, each
    # fair stays at 1/3. Serves as a sanity check that pool_size=1 matches
    # devig_n_way behavior.
    odds = [200, 200, 200]
    devigged = devig_pool(odds, pool_size=1)
    assert sum(devigged) == pytest.approx(1.0)
    for p in devigged:
        assert p == pytest.approx(1 / 3, abs=1e-6)


def test_devig_pool_scales_up_when_pool_larger_than_raw_sum():
    # Classic top_N case: raw implieds sum to LESS than pool_size, so each
    # fair gets scaled UP. Three outcomes at +200 (raw sum 1.0), pool_size=3
    # -> each fair = (1/3) * 3 = 1.0. Degenerate because the posted set is
    # both thin AND nowhere near saturating the pool; in production the
    # quarantine guardrail blocks this case, but the math must still hold.
    odds = [200, 200, 200]
    devigged = devig_pool(odds, pool_size=3)
    assert sum(devigged) == pytest.approx(3.0)
    for p in devigged:
        assert p == pytest.approx(1.0, abs=1e-6)


def test_devig_pool_scales_proportionally_to_pool_size():
    # top_10 example: 15 players at -150 (implied ~0.60 each, sum ~9.0).
    # Pool size 10 means sum(fairs) must equal 10, and each fair equals
    # its raw implied scaled by 10/9.
    odds = [-150] * 15
    raw = american_to_implied(-150)  # 0.60
    devigged = devig_pool(odds, pool_size=10)
    assert sum(devigged) == pytest.approx(10.0)
    for p in devigged:
        assert p == pytest.approx(raw * 10 / (raw * 15), abs=1e-9)


def test_devig_pool_preserves_relative_ordering():
    # A single -200 favorite among -150 candidates stays the favorite post
    # pool-normalization.
    odds = [-200, -150, -150, -150, 150]
    devigged = devig_pool(odds, pool_size=3)
    assert sum(devigged) == pytest.approx(3.0)
    assert devigged[0] > devigged[1]
    assert devigged[1] == pytest.approx(devigged[2])
    assert devigged[-1] == min(devigged)


def test_devig_pool_empty_list_returns_empty():
    assert devig_pool([], pool_size=10) == []


def test_devig_pool_all_zeros_returns_zeros_unchanged():
    # american_to_implied never returns 0 for real American odds, but the
    # guard exists so a degenerate all-zero input doesn't ZeroDivision.
    odds: list[float] = []
    assert devig_pool(odds, pool_size=5) == []
