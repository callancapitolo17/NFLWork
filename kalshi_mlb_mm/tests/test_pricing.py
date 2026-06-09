from kalshi_mlb_mm import pricing
from kalshi_common.ev_calc import maker_fee_per_contract


def test_roi_is_five_percent_net_of_fee():
    fair = 0.55
    q = pricing.quote(fair, target_roi=0.05)
    cost = q.yes_bid + maker_fee_per_contract(q.yes_bid)
    roi = (fair - cost) / cost
    assert abs(roi - 0.05) < 0.01


def test_sum_below_one():
    q = pricing.quote(0.50, target_roi=0.05)
    assert q.yes_bid + q.no_bid < 1.0


def test_grid_rounded_to_milli():
    q = pricing.quote(0.637, target_roi=0.05)
    assert abs(q.yes_bid * 1000 - round(q.yes_bid * 1000)) < 1e-6


def test_none_outside_prob_bounds():
    assert pricing.quote(0.0, 0.05) is None
    assert pricing.quote(1.0, 0.05) is None


def test_symmetry_at_fifty():
    q = pricing.quote(0.50, 0.05)
    assert abs(q.yes_bid - q.no_bid) < 1e-9


import pytest


@pytest.mark.parametrize("fair", [round(x * 0.001, 3) for x in range(60, 941)])
def test_roi_at_least_target_across_range(fair):
    """Realized ROI (net of maker fee) must never dip below the 5% target,
    at any fair value across the practical range — guards the fee-tier bug."""
    q = pricing.quote(fair, 0.05)
    if q is None:
        return
    for bid, p in ((q.yes_bid, fair), (q.no_bid, 1.0 - fair)):
        cost = bid + maker_fee_per_contract(bid)
        roi = (p - cost) / cost
        assert roi >= 0.05 - 1e-6, f"fair={fair} side p={p} bid={bid} roi={roi}"
        assert q.yes_bid + q.no_bid < 1.0
