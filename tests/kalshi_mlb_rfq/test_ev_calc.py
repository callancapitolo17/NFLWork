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


@pytest.mark.parametrize("fair,yes_bid,no_bid", [
    # Sweep across realistic LP spreads. For each row the LP keeps a positive
    # spread (yes_bid + no_bid < 1), which forces yes_ask + no_ask > 1, which
    # — together with the fee on each side — makes the both-+EV condition
    # mathematically unreachable. See plan: this is what the math-invariant
    # guard in _evaluate_quote defends against firing in practice.
    (0.50, 0.40, 0.40),   # symmetric 10% spread
    (0.50, 0.45, 0.45),   # symmetric 5% spread (tight LP)
    (0.30, 0.20, 0.60),   # asymmetric: LP heavy on NO
    (0.70, 0.60, 0.20),   # asymmetric: LP heavy on YES
    (0.85, 0.80, 0.10),   # near-extreme YES favourite
    (0.15, 0.10, 0.80),   # near-extreme NO favourite
])
def test_math_invariant_both_sides_cannot_be_positive(fair, yes_bid, no_bid):
    """For any LP making money (yes_bid + no_bid < 1), the post-fee EV gates
    cannot both pass against the same fair value. This is the invariant the
    new math-invariant guard in main._evaluate_quote relies on."""
    assert yes_bid + no_bid < 1, "test fixture broken: LP must keep a spread"
    _, ev_yes_pct = ev_calc.post_fee_ev_buy_yes(fair, no_bid)
    _, ev_no_pct  = ev_calc.post_fee_ev_buy_no (fair, yes_bid)
    # At least one side must be -EV. We require strict (<= 0 on one) rather
    # than just "they can't both be > 5%" so we don't depend on MIN_EV_PCT.
    assert ev_yes_pct <= 0 or ev_no_pct <= 0


def test_no_side_kelly_inversion_math():
    """If fair = P(YES), then P(NO) = 1 - fair. NO-side EV against yes_bid
    should equal YES-side EV against no_bid when fair' = 1 - fair AND we
    swap bid roles. This is the symmetry the side-aware Kelly relies on."""
    fair = 0.32
    # Mirror: re-frame this as a YES bet on the opposite outcome.
    ev_no_dollars, _ = ev_calc.post_fee_ev_buy_no(blended_fair=fair, yes_bid=0.13)
    ev_yes_mirror, _ = ev_calc.post_fee_ev_buy_yes(blended_fair=1 - fair, no_bid=0.13)
    assert ev_no_dollars == pytest.approx(ev_yes_mirror, abs=1e-9)
