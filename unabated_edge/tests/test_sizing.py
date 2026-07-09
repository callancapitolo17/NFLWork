from unabated_edge import sizing


def test_zero_when_no_edge():
    assert sizing.kelly_contracts(0.40, 0.50, 1000, 0.25, 30) == 0


def test_positive_and_capped():
    n = sizing.kelly_contracts(0.60, 0.45, 1000, 0.25, 30)
    assert n >= 1 and n * 0.45 <= 30 + 0.45
