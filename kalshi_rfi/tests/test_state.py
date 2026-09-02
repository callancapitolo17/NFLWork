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
        assert st.resting_tickers() == sorted([T1, T2])
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


class TestMultipleOrdersPerTicker:
    """Bug 3: resting used to be keyed by TICKER, so a second order on the
    same game overwrote the first order's id — orphaning a live order from
    every cancel, cap and sweep by construction."""

    def test_second_order_does_not_erase_the_first(self):
        st = RfiState()
        st.on_place(T1, "o1", 40, 10)      # $4.00
        st.on_place(T1, "o2", 50, 10)      # $5.00
        assert {oid for oid, _ in st.resting_for(T1)} == {"o1", "o2"}
        assert st.game_resting_cost_usd(T1) == pytest.approx(9.0)
        assert st.resting_cost_usd() == pytest.approx(9.0)

    def test_both_orders_count_against_game_exposure(self):
        st = RfiState()
        st.on_place(T1, "o1", 40, 10)
        st.on_place(T1, "o2", 50, 10)
        st.fills[T1] = [(2.0, 50, trading_day(NOW))]            # $1.00
        assert st.game_exposure_usd(T1) == pytest.approx(10.0)

    def test_cancel_is_per_order_id(self):
        st = RfiState()
        st.on_place(T1, "o1", 40, 10)
        st.on_place(T1, "o2", 50, 10)
        st.on_cancel("o1")
        assert [oid for oid, _ in st.resting_for(T1)] == ["o2"]
        assert st.resting_tickers() == [T1]
        st.on_cancel("o2")
        assert st.resting_tickers() == []

    def test_fill_decrements_only_the_filled_order(self, tmp_path, monkeypatch):
        import kalshi_rfi.state as state_mod
        con = storage.connect(tmp_path / "t.duckdb")
        st = RfiState()
        st.on_place(T1, "o1", 40, 10)
        st.on_place(T1, "o2", 50, 10)
        monkeypatch.setattr(state_mod.auth_client, "api", lambda *a, **k: (
            200, {"fills": [{"trade_id": "tr1", "order_id": "o2",
                             "count": 10, "yes_price": 50}]}, None))
        assert state_mod.poll_fills(st, con, NOW) == 1
        con.close()
        assert [oid for oid, _ in st.resting_for(T1)] == ["o1"]   # o1 intact
        assert st.game_filled_cost_usd(T1) == pytest.approx(5.0)
        assert st.in_fill_cooldown(T1, 60.0)


class TestReconcile:
    """Local tracking drifts; Kalshi's resting listing is the truth."""

    class _Gw:
        def __init__(self):
            self.cancelled = []

        def cancel(self, order_id, *, ticker=None, exchange_index=None):
            self.cancelled.append((order_id, ticker, exchange_index))
            return True

    def _api(self, monkeypatch, orders, status=200):
        import kalshi_rfi.state as state_mod
        monkeypatch.setattr(state_mod.auth_client, "api",
                            lambda *a, **k: (status, {"orders": orders}, None))

    def test_untracked_order_is_cancelled_on_its_own_shard(self, monkeypatch,
                                                            tmp_path):
        import kalshi_rfi.state as state_mod
        st, gw = RfiState(), self._Gw()
        con = storage.connect(tmp_path / "t.duckdb")
        self._api(monkeypatch, [{"order_id": "ghost", "ticker": T1,
                                 "exchange_index": 3, "remaining_count": 5}])
        assert state_mod.reconcile_resting_orders(st, gw, con)["cancelled"] == 1
        assert gw.cancelled == [("ghost", T1, 3)]
        # An order the bot never knew about still lands in the audit trail.
        rows = con.execute("SELECT ticker, order_id, action, reason "
                           "FROM orders").fetchall()
        con.close()
        assert rows == [(T1, "ghost", "cancel", "reconcile_untracked")]

    def test_vanished_order_is_dropped_from_local_state(self, monkeypatch):
        import kalshi_rfi.state as state_mod
        st, gw = RfiState(), self._Gw()
        st.on_place(T1, "o1", 40, 10)
        st.resting["o1"]["placed_mono"] -= 60.0     # past the placement grace
        self._api(monkeypatch, [])
        assert state_mod.reconcile_resting_orders(st, gw)["dropped"] == 1
        assert st.resting_for(T1) == []

    def test_just_placed_order_is_not_dropped(self, monkeypatch):
        import kalshi_rfi.state as state_mod
        st, gw = RfiState(), self._Gw()
        st.on_place(T1, "o1", 40, 10)               # listing may lag the POST
        self._api(monkeypatch, [])
        assert state_mod.reconcile_resting_orders(st, gw)["dropped"] == 0
        assert len(st.resting_for(T1)) == 1

    def test_partial_fill_count_and_shard_come_from_kalshi(self, monkeypatch):
        import kalshi_rfi.state as state_mod
        st, gw = RfiState(), self._Gw()
        st.on_place(T1, "o1", 40, 10)
        self._api(monkeypatch, [{"order_id": "o1", "ticker": T1,
                                 "exchange_index": 3, "remaining_count": 4}])
        state_mod.reconcile_resting_orders(st, gw)
        rec = st.resting_for(T1)[0][1]
        assert rec["count"] == 4 and rec["exchange_index"] == 3
        assert gw.cancelled == []

    def test_failed_fetch_holds_all_state(self, monkeypatch):
        import kalshi_rfi.state as state_mod
        st, gw = RfiState(), self._Gw()
        st.on_place(T1, "o1", 40, 10)
        st.resting["o1"]["placed_mono"] -= 60.0
        self._api(monkeypatch, [], status=503)
        assert state_mod.reconcile_resting_orders(st, gw) == {}
        assert len(st.resting_for(T1)) == 1         # never blanked on error


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
