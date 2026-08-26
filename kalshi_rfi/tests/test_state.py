"""Unit tests for exposure accounting and startup hydration (no network)."""
import datetime

import pytest

from kalshi_rfi import storage
from kalshi_rfi.state import RfiState, trading_day

NOW = datetime.datetime(2026, 8, 26, 1, 0, tzinfo=datetime.timezone.utc)
T1 = "KXMLBRFI-26AUG251905AAABBB"
T2 = "KXMLBRFI-26AUG252105CCCDDD"


class TestExposure:
    def test_resting_cost_and_exclusion(self):
        st = RfiState()
        st.on_place(T1, "o1", 40, 25)      # $10.00
        st.on_place(T2, "o2", 50, 10)      # $5.00
        assert st.resting_cost_usd() == pytest.approx(15.0)
        assert st.resting_cost_usd(exclude_ticker=T1) == pytest.approx(5.0)

    def test_daily_committed_counts_fills_and_resting(self):
        st = RfiState()
        st.on_place(T1, "o1", 40, 25)
        st.fills[T2] = [(10.0, 50, trading_day(NOW))]
        assert st.daily_committed_usd(NOW) == pytest.approx(15.0)

    def test_settled_fill_frees_game_cap_not_daily_cap(self):
        st = RfiState()
        st.fills[T1] = [(10.0, 50, trading_day(NOW))]
        st.settled.add(T1)
        assert st.game_filled_cost_usd(T1) == 0.0
        assert st.daily_filled_cost_usd(NOW) == pytest.approx(5.0)

    def test_yesterdays_fills_do_not_count_today(self):
        st = RfiState()
        yesterday = trading_day(NOW - datetime.timedelta(days=1))
        st.fills[T1] = [(10.0, 50, yesterday)]
        assert st.daily_filled_cost_usd(NOW) == 0.0
        assert st.game_filled_cost_usd(T1) == pytest.approx(5.0)


class TestHydration:
    def test_restart_restores_fills_and_settled(self, tmp_path):
        db = tmp_path / "t.duckdb"
        con = storage.connect(db)
        storage.log_fill(con, "tr1", "o1", T1, 40, 20, 0.10)
        storage.log_fill(con, "tr2", "o2", T2, 50, 10, 0.05)
        storage.log_settlement(con, T2, 3.0)

        st = RfiState()
        st.hydrate_from_db(con)
        con.close()

        now = datetime.datetime.now(datetime.timezone.utc)
        assert st.game_filled_cost_usd(T1) == pytest.approx(8.0)
        assert st.game_filled_cost_usd(T2) == 0.0          # settled
        assert st.daily_filled_cost_usd(now) == pytest.approx(13.0)
        assert "tr1" in st._done_trades                     # no re-ingest

    def test_restart_restores_prior_run_order_map(self, tmp_path):
        # A fill on a PREVIOUS run's order must still attribute after a
        # restart: poll_fills skips unknown order_ids, so our_orders is
        # rebuilt from the persisted 'place' rows (live mode only).
        db = tmp_path / "t.duckdb"
        con = storage.connect(db)
        storage.log_order(con, T1, "prior-oid", "place", 40, 25, "live", "ok")
        storage.log_order(con, T2, None, "place", 41, 25, "live",
                          "place_failed")            # no order id → skipped
        storage.log_order(con, T2, "shadow-1", "place", 41, 25, "shadow",
                          "ok")                      # shadow → skipped
        st = RfiState()
        st.hydrate_from_db(con)
        con.close()
        assert st.our_orders == {"prior-oid": T1}

    def test_hydration_is_idempotent(self, tmp_path):
        db = tmp_path / "t.duckdb"
        con = storage.connect(db)
        storage.log_fill(con, "tr1", "o1", T1, 40, 20, 0.10)
        st = RfiState()
        st.hydrate_from_db(con)
        st.hydrate_from_db(con)
        con.close()
        assert len(st.fills[T1]) == 1


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
