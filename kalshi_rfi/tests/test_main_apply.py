"""Unit tests for _apply_decision order management (shadow gateway, tmp DB)."""
from datetime import datetime, timezone

import pytest

from kalshi_rfi import storage
from kalshi_rfi.discovery import RfiGame
from kalshi_rfi.engine import QuoteDecision
from kalshi_rfi.gateway import ShadowGateway
from kalshi_rfi.main import _apply_decision, _fit_to_caps
from kalshi_rfi.state import RfiState, trading_day

NOW = datetime(2026, 8, 26, 1, 0, tzinfo=timezone.utc)
CAPS = {"per_game_cap_usd": 10.0, "daily_cap_usd": 100.0}


def _game(ticker="KXMLBRFI-26AUG251905AAABBB", exchange_index=3):
    return RfiGame(ticker=ticker, suffix=ticker.split("-", 1)[1],
                   home_team="H", away_team="A",
                   commence_utc=NOW.replace(tzinfo=None),
                   yes_bid_cents=38, yes_ask_cents=46, status="active",
                   exchange_index=exchange_index)


def _apply(gw, st, con, g, d, **kw):
    return _apply_decision(gw, st, con, g, d, NOW, **{**CAPS, **kw})


@pytest.fixture
def con(tmp_path):
    c = storage.connect(tmp_path / "t.duckdb")
    yield c
    c.close()


class TestApplyDecision:
    def test_place_then_hold_on_same_price(self, con):
        gw, st, g = ShadowGateway(), RfiState(), _game()
        d = QuoteDecision("quote", "ok", price_cents=41, count=24)
        _apply(gw, st, con, g, d)
        oid = st.resting_for(g.ticker)[0][0]
        _apply(gw, st, con, g, d)
        assert st.resting_for(g.ticker)[0][0] == oid   # held, not churned

    def test_reprice_cancels_and_replaces(self, con):
        gw, st, g = ShadowGateway(), RfiState(), _game()
        _apply(gw, st, con, g, QuoteDecision("quote", "ok", 41, 24))
        _apply(gw, st, con, g, QuoteDecision("quote", "ok", 42, 23))
        assert len(st.resting_for(g.ticker)) == 1
        assert st.resting_for(g.ticker)[0][1]["price_cents"] == 42

    def test_same_price_downsize_replaces(self, con):
        # Fills elsewhere shrank headroom: an oversized resting order at an
        # unchanged price is a cap breach, not hysteresis (finding 4).
        gw, st, g = ShadowGateway(), RfiState(), _game()
        _apply(gw, st, con, g, QuoteDecision("quote", "ok", 41, 24))
        _apply(gw, st, con, g, QuoteDecision("quote", "ok", 41, 10))
        assert st.resting_for(g.ticker)[0][1]["count"] == 10

    def test_resize_shrinks_replacement_after_cancel(self, con):
        # A fill landing between the cycle's poll and the replace must
        # shrink the new order via the pre-placement resize (finding 4).
        gw, st, g = ShadowGateway(), RfiState(), _game()
        _apply(gw, st, con, g, QuoteDecision("quote", "ok", 41, 24))
        _apply(gw, st, con, g, QuoteDecision("quote", "ok", 42, 24),
               resize=lambda price: 9)
        assert st.resting_for(g.ticker)[0][1]["count"] == 9

    def test_resize_zero_blocks_replacement(self, con):
        gw, st, g = ShadowGateway(), RfiState(), _game()
        _apply(gw, st, con, g, QuoteDecision("quote", "ok", 41, 24))
        _apply(gw, st, con, g, QuoteDecision("quote", "ok", 42, 24),
               resize=lambda price: 0)
        assert st.resting_for(g.ticker) == []    # cancelled, nothing placed

    def test_no_quote_cancels(self, con):
        gw, st, g = ShadowGateway(), RfiState(), _game()
        _apply(gw, st, con, g, QuoteDecision("quote", "ok", 41, 24))
        _apply(gw, st, con, g, QuoteDecision("no_quote", "no_consensus"))
        assert st.resting_for(g.ticker) == []

    def test_place_passes_exchange_index(self, con):
        # Sharding: the order must carry the market's shard, not a guess.
        seen = {}

        class Gw(ShadowGateway):
            def place_yes_bid(self, ticker, price_cents, count, coid,
                              exchange_index=None):
                seen["exchange_index"] = exchange_index
                return super().place_yes_bid(ticker, price_cents, count, coid)

        st, g = RfiState(), _game(exchange_index=3)
        _apply(Gw(), st, con, g, QuoteDecision("quote", "ok", 41, 10))
        assert seen["exchange_index"] == 3
        assert st.resting_for(g.ticker)[0][1]["exchange_index"] == 3


class TestFreshQuotePathCapBackstop:
    """Bug 1 (observed 2026-08-27): an order that filled COMPLETELY is
    dropped from state by poll_fills, so the next cycle takes the
    fresh-quote path — which used to skip the fill re-poll entirely and
    size from the stale top-of-cycle count. One game took 5 fills / 45
    contracts (~$24) against a $5 cap in ~60 seconds."""

    def test_resize_runs_on_the_fresh_quote_path(self, con):
        gw, st, g = ShadowGateway(), RfiState(), _game()
        assert st.resting_for(g.ticker) == []          # nothing to replace
        _apply(gw, st, con, g, QuoteDecision("quote", "ok", 50, 20),
               resize=lambda price: 4)
        assert st.resting_for(g.ticker)[0][1]["count"] == 4

    def test_backstop_refuses_order_past_per_game_cap(self, con):
        gw, st, g = ShadowGateway(), RfiState(), _game()
        # $5 already filled on this game against a $5 cap: no headroom, and
        # the engine's stale count of 9 must not reach the wire.
        st.fills[g.ticker] = [(10.0, 50, trading_day(NOW))]
        _apply(gw, st, con, g, QuoteDecision("quote", "ok", 53, 9),
               per_game_cap_usd=5.0, daily_cap_usd=100.0)
        assert st.resting_for(g.ticker) == []
        reasons = [r[0] for r in con.execute(
            "SELECT reason FROM orders WHERE action = 'place'").fetchall()]
        assert reasons == ["caps_exhausted_pre_place"]

    def test_backstop_shrinks_order_to_remaining_headroom(self, con):
        gw, st, g = ShadowGateway(), RfiState(), _game()
        st.fills[g.ticker] = [(6.0, 50, trading_day(NOW))]      # $3.00 of $5
        _apply(gw, st, con, g, QuoteDecision("quote", "ok", 50, 9),
               per_game_cap_usd=5.0, daily_cap_usd=100.0)
        assert st.resting_for(g.ticker)[0][1]["count"] == 4     # $2.00 left

    def test_backstop_ignores_an_oversized_resize(self, con):
        # The backstop does not trust the resize callback either — it is the
        # last check before the wire.
        gw, st, g = ShadowGateway(), RfiState(), _game()
        st.fills[g.ticker] = [(8.0, 50, trading_day(NOW))]      # $4.00 of $5
        _apply(gw, st, con, g, QuoteDecision("quote", "ok", 50, 99),
               resize=lambda price: 99,
               per_game_cap_usd=5.0, daily_cap_usd=100.0)
        assert st.resting_for(g.ticker)[0][1]["count"] == 2     # $1.00 left

    def test_backstop_counts_resting_orders_on_the_same_game(self, con):
        gw, st, g = ShadowGateway(), RfiState(), _game()
        st.on_place(g.ticker, "stray-1", 50, 8)                 # $4.00 live
        assert _fit_to_caps(st, g.ticker, 50, 99, NOW, 5.0, 100.0) == 2

    def test_backstop_respects_the_daily_cap(self, con):
        st, g = RfiState(), _game()
        st.fills["KXMLBRFI-26AUG252105CCCDDD"] = [
            (190.0, 50, trading_day(NOW))]                      # $95 of $100
        assert _fit_to_caps(st, g.ticker, 50, 99, NOW, 10.0, 100.0) == 10


class TestPostFillCooldown:
    def test_cooldown_blocks_requoting_a_just_filled_game(self):
        st, g = RfiState(), _game()
        assert not st.in_fill_cooldown(g.ticker, 60.0)
        st.last_fill_mono[g.ticker] = __import__("time").monotonic()
        assert st.in_fill_cooldown(g.ticker, 60.0)
        assert not st.in_fill_cooldown(g.ticker, 0.0)     # knob off


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
