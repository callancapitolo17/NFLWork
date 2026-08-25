"""Unit tests for the pure quote-decision logic."""
import math

import pytest

from kalshi_common.ev_calc import maker_fee_per_contract
from kalshi_rfi.engine import decide, desired_price_cents, size_contracts

MARGIN = 3
MIN_EDGE = 2.0


def _decide(fair, ask, **overrides):
    kwargs = dict(margin_cents=MARGIN, min_edge_cents=MIN_EDGE,
                  per_game_cap_usd=10.0, game_committed_usd=0.0,
                  daily_remaining_usd=100.0,
                  sec_to_start=3600.0, pull_before_start_sec=120.0)
    kwargs.update(overrides)
    return decide(fair, ask, **kwargs)


class TestDesiredPrice:
    def test_fair_minus_margin(self):
        # fair 44% -> floor(44) - 3 = 41c
        price, reason = desired_price_cents(0.44, 50, MARGIN, MIN_EDGE)
        assert (price, reason) == (41, "ok")

    def test_clamped_below_ask_stays_maker(self):
        # desired 41c but the ask sits at 40c -> rest at 39c, never cross
        price, _ = desired_price_cents(0.44, 40, MARGIN, MIN_EDGE)
        assert price == 39

    def test_clamp_only_ever_lowers_price(self):
        unclamped, _ = desired_price_cents(0.44, None, MARGIN, MIN_EDGE)
        clamped, _ = desired_price_cents(0.44, 10, MARGIN, MIN_EDGE)
        assert clamped < unclamped

    def test_fee_leaves_at_least_min_edge(self):
        price, _ = desired_price_cents(0.44, 50, MARGIN, MIN_EDGE)
        fee_cents = maker_fee_per_contract(price / 100.0) * 100.0
        assert 0.44 * 100 - price - fee_cents >= MIN_EDGE

    def test_edge_floor_rejects_thin_quotes(self):
        # margin 1c cannot clear a 2c post-fee floor
        price, reason = desired_price_cents(0.44, 50, 1, MIN_EDGE)
        assert price is None and reason == "edge_below_floor"

    def test_extreme_fairs_rejected(self):
        assert desired_price_cents(0.0, 50, MARGIN, MIN_EDGE)[0] is None
        assert desired_price_cents(1.0, 50, MARGIN, MIN_EDGE)[0] is None
        assert desired_price_cents(0.03, 50, MARGIN, MIN_EDGE)[1] == "price_below_1c"


class TestSize:
    def test_cap_over_price(self):
        # $10 at 41c -> 24 contracts
        assert size_contracts(41, 10.0, 0.0, 100.0) == math.floor(1000 / 41)

    def test_daily_headroom_binds(self):
        assert size_contracts(41, 10.0, 0.0, 2.0) == math.floor(200 / 41)

    def test_game_committed_shrinks_headroom(self):
        assert size_contracts(41, 10.0, 9.9, 100.0) == 0

    def test_no_negative_counts(self):
        assert size_contracts(41, 10.0, 12.0, 100.0) == 0


class TestDecide:
    def test_happy_path_quotes(self):
        d = _decide(0.44, 50)
        assert (d.action, d.price_cents, d.count) == ("quote", 41, 24)

    def test_pull_deadline_beats_everything(self):
        d = _decide(0.44, 50, sec_to_start=60.0)
        assert (d.action, d.reason) == ("no_quote", "past_pull_deadline")

    def test_no_consensus_cancels(self):
        d = _decide(None, 50)
        assert (d.action, d.reason) == ("no_quote", "no_consensus")

    def test_caps_exhausted_cancels(self):
        d = _decide(0.44, 50, game_committed_usd=10.0)
        assert (d.action, d.reason) == ("no_quote", "caps_exhausted")

    def test_daily_cap_exhausted_cancels(self):
        d = _decide(0.44, 50, daily_remaining_usd=0.0)
        assert (d.action, d.reason) == ("no_quote", "caps_exhausted")


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
