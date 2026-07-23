"""P-7 (issue #21): the post-fill cooldown must outlast the data refresh, and
same-price re-quotes must be blocked until consensus fair actually moves.

The 60s COMBO_COOLDOWN_SEC floor expires before the ~150-165s SGP scrape cycle
completes, so previously a counterparty who picked us off could post a new RFQ
on the same combo and hit the SAME stale price again. Two new gates:

- `in_cooldown_awaiting_refresh`: floor passed, but no post-fill scrape has
  landed for every game the combo touches.
- `same_price_block`: books refreshed, but consensus fair is still within
  QUOTE_HYSTERESIS of the fair the fill transacted against.

Gate tests mirror the N7 harness (mocked _SGP_ODDS / scope cache / router
fair, fresh DuckDB per test).
"""
import importlib
from datetime import datetime, timedelta, timezone

import pandas as pd

_EVT = "KXMLBGAME-25JUN271905TEXLAA"
_LEGS = [{"market_ticker": "KXMLBSPREAD-25JUN271905TEXLAA-LAA2",
          "event_ticker": _EVT, "side": "yes"},
         {"market_ticker": "KXMLBTOTAL-25JUN271905TEXLAA-9",
          "event_ticker": _EVT, "side": "yes"}]

_TICKER = "COMBO-CD"
_GAME = "GAME-CD"
_FILL_FAIR = 0.55


def _cooldown_env(monkeypatch, tmp_path, db_name, *, fetch_time, cooled_until,
                  fair_now):
    """Discovery-tick env with one prior fill on _TICKER (fair_at_confirm =
    _FILL_FAIR, filled 120s ago) and its combo_cooldown row. Returns (db, main,
    src, gw, filled_at)."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    import kalshi_mlb_mm.router as router_mod
    from kalshi_mlb_mm import main

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / db_name)
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    importlib.reload(db)
    db.init_database()

    monkeypatch.setattr(main, "_SGP_ODDS",
                        pd.DataFrame({"game_id": [_GAME], "combo": ["c"],
                                      "period": ["FG"], "bookmaker": ["dk"],
                                      "sgp_decimal": [2.0],
                                      "fetch_time": [fetch_time],
                                      "spread_line": [-1.5],
                                      "total_line": [8.5]}))
    monkeypatch.setattr(risk, "tipoff_ok", lambda ct, m: True)
    monkeypatch.setattr(router_mod, "combo_fair", lambda *a, **k: fair_now)
    monkeypatch.setattr(main, "_commence_time", lambda gid: None)
    monkeypatch.setattr(main, "_PREV_BOOK_FAIR", {})
    monkeypatch.setattr(main, "_SCOPE_CACHE", {_TICKER: (True, _GAME, _LEGS)})
    monkeypatch.setattr(main, "_resolve_game_for_legs", lambda gl: _GAME)
    monkeypatch.setattr(main, "_leg_market_prices",
                        lambda legs: {"L": {"yes_bid": 0.5, "yes_ask": 0.52}})

    filled_at = datetime.now(timezone.utc) - timedelta(seconds=120)
    with db.connect() as con:
        con.execute(
            "INSERT INTO fills (fill_id, quote_id, rfq_id, combo_market_ticker, "
            "game_id, side_held, contracts, price, fee, model_fair_at_quote, "
            "book_fair_at_quote, blended_fair_at_quote, fair_at_confirm, "
            "realized_pnl, filled_at, reconciled) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            ["fill-cd", "qid-cd", "r-cd-prior", _TICKER, _GAME, "yes", 1, 0.50,
             0.0, None, _FILL_FAIR, _FILL_FAIR, _FILL_FAIR, None, filled_at, True])
        con.execute("INSERT OR REPLACE INTO fill_games (fill_id, game_id) "
                    "VALUES (?, ?)", ["fill-cd", _GAME])
        con.execute(
            "INSERT OR REPLACE INTO combo_cooldown "
            "(combo_market_ticker, cooled_until) VALUES (?, ?)",
            [_TICKER, cooled_until])

    class Src:
        def poll(self):
            return [{"id": "r-cd-new", "market_ticker": _TICKER, "contracts": 3}]

        def get_market(self, t):
            return {}

    class GW:
        def __init__(self):
            self.submits = []

        def submit_quote(self, *a):
            self.submits.append(a)
            return "qid-cd-new"

    return db, main, Src(), GW(), filled_at


def _last_decision(db, rfq_id="r-cd-new"):
    with db.connect(read_only=True) as con:
        return con.execute(
            "SELECT decision, reason FROM quote_decisions "
            "WHERE rfq_id=? ORDER BY observed_at DESC LIMIT 1",
            [rfq_id]).fetchone()


def test_floor_active_still_skips_in_cooldown(monkeypatch, tmp_path):
    """Pins existing behavior: within the COMBO_COOLDOWN_SEC floor the combo
    is skipped with the original in_cooldown reason."""
    now = datetime.now(timezone.utc)
    db, main, src, gw, _ = _cooldown_env(
        monkeypatch, tmp_path, "cd_floor.duckdb",
        fetch_time=pd.Timestamp(now - timedelta(seconds=10)),
        cooled_until=now + timedelta(seconds=60),   # floor still active
        fair_now=0.58)
    main._discovery_tick(src, gw, dry_run=False)
    assert gw.submits == []
    assert _last_decision(db) == ("skipped", "in_cooldown")


def test_floor_passed_but_no_post_fill_scrape_skips(monkeypatch, tmp_path):
    """THE regression test for issue #21: cooldown floor expired, but the
    latest scrape for the combo's game predates the fill → the books that
    priced the picked-off fill are still live → must stay cooled."""
    now = datetime.now(timezone.utc)
    db, main, src, gw, filled_at = _cooldown_env(
        monkeypatch, tmp_path, "cd_stale.duckdb",
        fetch_time=pd.Timestamp(filled_at_offset := now - timedelta(seconds=300)),
        cooled_until=now - timedelta(seconds=60),   # floor expired
        fair_now=0.58)
    assert filled_at_offset < filled_at  # scrape strictly predates the fill
    main._discovery_tick(src, gw, dry_run=False)
    assert gw.submits == [], (
        "must not re-quote before a post-fill scrape lands")
    assert _last_decision(db) == ("skipped", "in_cooldown_awaiting_refresh")


def test_post_fill_scrape_and_fair_moved_allows_quote(monkeypatch, tmp_path):
    """Books refreshed after the fill AND consensus fair moved → quote."""
    now = datetime.now(timezone.utc)
    db, main, src, gw, _ = _cooldown_env(
        monkeypatch, tmp_path, "cd_moved.duckdb",
        fetch_time=pd.Timestamp(now - timedelta(seconds=10)),  # after the fill
        cooled_until=now - timedelta(seconds=60),
        fair_now=0.58)                                          # moved 3¢
    main._discovery_tick(src, gw, dry_run=False)
    assert len(gw.submits) == 1, "refreshed + moved fair must quote"
    assert _last_decision(db)[0] == "quoted"


def test_post_fill_scrape_but_fair_unchanged_same_price_block(monkeypatch, tmp_path):
    """Books refreshed after the fill but consensus fair is within
    QUOTE_HYSTERESIS of the fill's fair → same_price_block."""
    now = datetime.now(timezone.utc)
    db, main, src, gw, _ = _cooldown_env(
        monkeypatch, tmp_path, "cd_same.duckdb",
        fetch_time=pd.Timestamp(now - timedelta(seconds=10)),
        cooled_until=now - timedelta(seconds=60),
        fair_now=_FILL_FAIR + 0.002)                # within ε=0.005
    main._discovery_tick(src, gw, dry_run=False)
    assert gw.submits == [], (
        "unchanged fair must not be re-quoted at the just-picked-off price")
    assert _last_decision(db) == ("skipped", "same_price_block")


def test_no_cooldown_row_means_no_new_gates(monkeypatch, tmp_path):
    """A combo that never filled (no combo_cooldown row) is untouched by the
    refresh/same-price gates even with stale fetch_time and unchanged fair."""
    now = datetime.now(timezone.utc)
    db, main, src, gw, _ = _cooldown_env(
        monkeypatch, tmp_path, "cd_none.duckdb",
        fetch_time=pd.Timestamp(now - timedelta(seconds=300)),
        cooled_until=now - timedelta(seconds=60),
        fair_now=_FILL_FAIR)
    # Remove the cooldown row AND the fill so this combo looks never-filled.
    with db.connect() as con:
        con.execute("DELETE FROM combo_cooldown WHERE combo_market_ticker=?",
                    [_TICKER])
        con.execute("DELETE FROM fills WHERE combo_market_ticker=?", [_TICKER])
        con.execute("DELETE FROM fill_games WHERE fill_id='fill-cd'")
    main._discovery_tick(src, gw, dry_run=False)
    assert len(gw.submits) == 1
    assert _last_decision(db)[0] == "quoted"


# ---------------------------------------------------------------------------
# _post_fill_refresh_landed unit coverage (timezone + cross-game edges).
# ---------------------------------------------------------------------------
def _refresh_env(monkeypatch, df):
    from kalshi_mlb_mm import main
    monkeypatch.setattr(main, "_SGP_ODDS", df)
    return main


def _df(rows):
    return pd.DataFrame(rows, columns=["game_id", "fetch_time"])


def test_refresh_landed_aware_post_fill_true(monkeypatch):
    now = datetime.now(timezone.utc)
    main = _refresh_env(monkeypatch, _df([["g1", pd.Timestamp(now)]]))
    assert main._post_fill_refresh_landed(
        ["g1"], now - timedelta(seconds=60)) is True


def test_refresh_landed_aware_pre_fill_false(monkeypatch):
    now = datetime.now(timezone.utc)
    main = _refresh_env(monkeypatch,
                        _df([["g1", pd.Timestamp(now - timedelta(seconds=300))]]))
    assert main._post_fill_refresh_landed(["g1"], now) is False


def test_refresh_landed_naive_fetch_time_is_local(monkeypatch):
    """mlb_sgp_odds.fetch_time is a naive TIMESTAMP column — DuckDB stores the
    session-LOCAL wall clock. A fresh naive fetch_time must count as refreshed
    against an aware filled_at (naive→UTC coercion would be off by the UTC
    offset in any non-UTC timezone)."""
    naive_local_now = datetime.now()  # naive LOCAL
    aware_fill = datetime.now(timezone.utc) - timedelta(seconds=60)
    main = _refresh_env(monkeypatch, _df([["g1", pd.Timestamp(naive_local_now)]]))
    assert main._post_fill_refresh_landed(["g1"], aware_fill) is True
    old_naive_local = datetime.now() - timedelta(seconds=600)
    main = _refresh_env(monkeypatch, _df([["g1", pd.Timestamp(old_naive_local)]]))
    assert main._post_fill_refresh_landed(
        ["g1"], datetime.now(timezone.utc)) is False


def test_refresh_landed_cross_game_requires_every_game(monkeypatch):
    """One game refreshed, the other not → NOT refreshed (fail-closed)."""
    now = datetime.now(timezone.utc)
    fill = now - timedelta(seconds=120)
    main = _refresh_env(monkeypatch, _df([
        ["g1", pd.Timestamp(now)],                        # post-fill
        ["g2", pd.Timestamp(now - timedelta(seconds=300))],  # pre-fill
    ]))
    assert main._post_fill_refresh_landed(["g1", "g2"], fill) is False
    # Both refreshed → True.
    main = _refresh_env(monkeypatch, _df([
        ["g1", pd.Timestamp(now)],
        ["g2", pd.Timestamp(now)],
    ]))
    assert main._post_fill_refresh_landed(["g1", "g2"], fill) is True


def test_refresh_landed_fails_closed_on_missing_data(monkeypatch):
    now = datetime.now(timezone.utc)
    # Game absent from odds entirely.
    main = _refresh_env(monkeypatch, _df([["other", pd.Timestamp(now)]]))
    assert main._post_fill_refresh_landed(["g1"], now) is False
    # NaT fetch_time.
    main = _refresh_env(monkeypatch, _df([["g1", pd.NaT]]))
    assert main._post_fill_refresh_landed(["g1"], now) is False
    # Empty frame / None frame / None filled_at.
    main = _refresh_env(monkeypatch, _df([]))
    assert main._post_fill_refresh_landed(["g1"], now) is False
    main = _refresh_env(monkeypatch, None)
    assert main._post_fill_refresh_landed(["g1"], now) is False
    main = _refresh_env(monkeypatch, _df([["g1", pd.Timestamp(now)]]))
    assert main._post_fill_refresh_landed(["g1"], None) is False
