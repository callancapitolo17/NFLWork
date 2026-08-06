import pytest
from unabated_edge import pricing

def test_american_to_prob():
    assert pricing.american_to_prob(-200) == pytest.approx(2/3, abs=1e-6)
    assert pricing.american_to_prob(150) == pytest.approx(0.4, abs=1e-6)

def test_american_to_prob_rejects_decimal_odds_shape():
    """Finding #76: American odds are never in (-100, 100) by definition. A
    decimal-odds value leaking into this field (e.g. 1.91) must raise, not
    silently produce a plausible-but-wrong probability."""
    with pytest.raises(ValueError):
        pricing.american_to_prob(1.91)

def test_american_to_prob_boundary_100_is_allowed():
    """Exactly +-100 (pick 'em) is the boundary and must be allowed, not rejected."""
    assert pricing.american_to_prob(100) == pytest.approx(0.5, abs=1e-6)
    assert pricing.american_to_prob(-100) == pytest.approx(0.5, abs=1e-6)

def test_american_to_prob_valid_prices_unaffected():
    assert pricing.american_to_prob(-110) == pytest.approx(110/210, abs=1e-6)
    assert pricing.american_to_prob(145) == pytest.approx(100/245, abs=1e-6)

def test_devig_three_way_sums_to_one():
    out = pricing.devig([0.50, 0.30, 0.28])
    assert sum(out) == pytest.approx(1.0, abs=1e-6)
    assert out[0] > out[1]

def test_devig_two_way():
    out = pricing.devig([0.55, 0.52])
    assert sum(out) == pytest.approx(1.0, abs=1e-6)
