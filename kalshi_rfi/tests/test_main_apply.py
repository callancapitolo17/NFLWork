"""Unit tests for _apply_decision order management (shadow gateway, tmp DB)."""
from datetime import datetime, timezone

import pytest

from kalshi_rfi import storage
from kalshi_rfi.discovery import RfiGame
from kalshi_rfi.engine import QuoteDecision
from kalshi_rfi.gateway import ShadowGateway
from kalshi_rfi.main import _apply_decision
from kalshi_rfi.state import RfiState

NOW = datetime(2026, 8, 26, 1, 0, tzinfo=timezone.utc)


def _game(ticker="KXMLBRFI-26AUG251905AAABBB"):
    return RfiGame(ticker=ticker, suffix=ticker.split("-", 1)[1],
                   home_team="H", away_team="A",
                   commence_utc=NOW.replace(tzinfo=None),
                   yes_bid_cents=38, yes_ask_cents=46, status="active")


@pytest.fixture
def con(tmp_path):
    c = storage.connect(tmp_path / "t.duckdb")
    yield c
    c.close()


class TestApplyDecision:
    def test_place_then_hold_on_same_price(self, con):
        gw, st, g = ShadowGateway(), RfiState(), _game()
        d = QuoteDecision("quote", "ok", price_cents=41, count=24)
        _apply_decision(gw, st, con, g, d, NOW)
        oid = st.resting[g.ticker]["order_id"]
        _apply_decision(gw, st, con, g, d, NOW)
        assert st.resting[g.ticker]["order_id"] == oid   # held, not churned

    def test_reprice_cancels_and_replaces(self, con):
        gw, st, g = ShadowGateway(), RfiState(), _game()
        _apply_decision(gw, st, con, g,
                        QuoteDecision("quote", "ok", 41, 24), NOW)
        _apply_decision(gw, st, con, g,
                        QuoteDecision("quote", "ok", 42, 23), NOW)
        assert st.resting[g.ticker]["price_cents"] == 42

    def test_same_price_downsize_replaces(self, con):
        # Fills elsewhere shrank headroom: an oversized resting order at an
        # unchanged price is a cap breach, not hysteresis (finding 4).
        gw, st, g = ShadowGateway(), RfiState(), _game()
        _apply_decision(gw, st, con, g,
                        QuoteDecision("quote", "ok", 41, 24), NOW)
        _apply_decision(gw, st, con, g,
                        QuoteDecision("quote", "ok", 41, 10), NOW)
        assert st.resting[g.ticker]["count"] == 10

    def test_recompute_shrinks_replacement_after_cancel(self, con):
        # A fill landing between the cycle's poll and the replace must
        # shrink the new order via the post-cancel recompute (finding 4).
        gw, st, g = ShadowGateway(), RfiState(), _game()
        _apply_decision(gw, st, con, g,
                        QuoteDecision("quote", "ok", 41, 24), NOW)
        _apply_decision(gw, st, con, g,
                        QuoteDecision("quote", "ok", 42, 24), NOW,
                        recompute_count=lambda price: 9)
        assert st.resting[g.ticker]["count"] == 9

    def test_recompute_zero_blocks_replacement(self, con):
        gw, st, g = ShadowGateway(), RfiState(), _game()
        _apply_decision(gw, st, con, g,
                        QuoteDecision("quote", "ok", 41, 24), NOW)
        _apply_decision(gw, st, con, g,
                        QuoteDecision("quote", "ok", 42, 24), NOW,
                        recompute_count=lambda price: 0)
        assert g.ticker not in st.resting        # cancelled, nothing placed

    def test_no_quote_cancels(self, con):
        gw, st, g = ShadowGateway(), RfiState(), _game()
        _apply_decision(gw, st, con, g,
                        QuoteDecision("quote", "ok", 41, 24), NOW)
        _apply_decision(gw, st, con, g,
                        QuoteDecision("no_quote", "no_consensus"), NOW)
        assert g.ticker not in st.resting


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
