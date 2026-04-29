import pytest

from kalshi_mlb_rfq import ev_calc


@pytest.mark.parametrize("price,expected_fee", [
    (0.50, 0.02),   # 0.07 * 0.5 * 0.5 = 0.0175 → ceil to 0.02
    (0.20, 0.02),   # 0.07 * 0.20 * 0.80 = 0.0112 → ceil to 0.02
    (0.10, 0.01),   # 0.07 * 0.1 * 0.9 = 0.0063 → ceil to 0.01
    (0.30, 0.02),   # 0.07 * 0.30 * 0.70 = 0.0147 → ceil to 0.02
])
def test_fee_per_contract(price, expected_fee):
    assert ev_calc.fee_per_contract(price) == pytest.approx(expected_fee)


def test_post_fee_ev_yes_side_buying():
    ev_dollars, ev_pct = ev_calc.post_fee_ev_buy_yes(blended_fair=0.32, no_bid=0.74)
    assert ev_dollars == pytest.approx(0.04, abs=0.001)
    assert ev_pct == pytest.approx(0.04 / 0.26, abs=0.001)


def test_post_fee_ev_no_side_buying():
    ev_dollars, ev_pct = ev_calc.post_fee_ev_buy_no(blended_fair=0.32, yes_bid=0.13)
    expected_fee = ev_calc.fee_per_contract(0.87)
    expected_ev = 0.68 * (1 - 0.87) - 0.32 * 0.87 - expected_fee
    assert ev_dollars == pytest.approx(expected_ev, abs=0.001)


def test_extreme_prices_return_zero():
    assert ev_calc.post_fee_ev_buy_yes(0.5, 1.0) == (0.0, 0.0)
    assert ev_calc.post_fee_ev_buy_yes(0.5, 0.0) == (0.0, 0.0)
