import pytest
from nfl_draft.lib.devig import american_to_implied, devig_two_way, devig_n_way, proportional_devig


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
