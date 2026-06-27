import pytest
from unabated_edge import pricing

def test_american_to_prob():
    assert pricing.american_to_prob(-200) == pytest.approx(2/3, abs=1e-6)
    assert pricing.american_to_prob(150) == pytest.approx(0.4, abs=1e-6)

def test_devig_three_way_sums_to_one():
    out = pricing.devig([0.50, 0.30, 0.28])
    assert sum(out) == pytest.approx(1.0, abs=1e-6)
    assert out[0] > out[1]

def test_devig_two_way():
    out = pricing.devig([0.55, 0.52])
    assert sum(out) == pytest.approx(1.0, abs=1e-6)
