import datetime
import os
import pytest
from unabated_edge import feed, runner, config, storage
from unabated_edge.sports.soccer import Soccer

_NOW = datetime.datetime(2026, 6, 1, tzinfo=datetime.timezone.utc)


def _state():
    """FeedState with anchor bt3 Over/Under for Argentina vs Austria."""
    return feed.parse_snapshot(
        {
            "marketSources": [{"id": 7, "name": "S"}],
            "teams": {"1": {"name": "Argentina"}, "2": {"name": "Austria"}},
            "gameOddsEvents": {
                "lg21:pt1:pregame": [
                    {
                        "eventId": 1,
                        "eventStart": "2026-12-31T17:00:00+00:00",
                        "eventTeams": {"1": {"id": 1}, "0": {"id": 2}},
                        "gameOddsMarketSourcesLines": {
                            "si0:ms7:an0": {"bt3": {"price": 115, "points": 2.5}},   # Over
                            "si1:ms7:an0": {"bt3": {"price": -140, "points": 2.5}},  # Under
                        },
                    }
                ]
            },
        },
        {"lg21"},
    )


# KXWCTOTAL-style Kalshi event with one Over-2.5 market
_KEV = {
    "title": "Argentina vs Austria: Regulation Time Total Goals",
    "markets": [
        {"ticker": "T-O25", "strike_type": "greater", "floor_strike": 2.5,
         "yes_sub_title": "Reg Time: Over 2.5 goals scored",
         "volume_fp": "1234.00", "open_interest_fp": "567.00"},   # live shape: STRING fp
    ],
}

# yes_ask = 1 - best no bid = 0.30; no_ask = 1 - best yes bid = 0.30 (both cheap,
# same semantics as the old `ask_fn=lambda t, side: 0.30` fixtures)
_BOOK = {"yes_bids": [(0.70, 5.0)], "no_bids": [(0.70, 8.0)]}


def _init_dbs(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "MARKET_DB_PATH", tmp_path / "m.duckdb")
    monkeypatch.setattr(config, "RESEARCH_DB_PATH", tmp_path / "r.duckdb")
    storage.init()


def test_tick_flags_positive_edge(tmp_path, monkeypatch):
    _init_dbs(tmp_path, monkeypatch)
    # 0.30 ask is well below a ~42% fair Over → positive EV
    rows = runner.run_tick(Soccer(), _state(), [_KEV], now=_NOW, dry_run=True,
                           book_fn=lambda t: _BOOK)
    assert any(r["ev_pct"] > 0 for r in rows)


def test_absolute_ev_floor_blocks_flag(tmp_path, monkeypatch):
    """An edge that clears ev_pct must still be blocked by the absolute-dollar floor."""
    _init_dbs(tmp_path, monkeypatch)
    monkeypatch.setattr(config, "MIN_EV_DOLLARS", 1.0)   # no real contract clears $1 EV
    rows = runner.run_tick(Soccer(), _state(), [_KEV], now=_NOW, dry_run=True,
                           book_fn=lambda t: _BOOK)
    assert rows == []


def test_no_snapshots_after_kickoff(tmp_path, monkeypatch):
    """Once now >= kickoff, neither Unabated lines nor Kalshi books are snapshotted —
    the last pre-kickoff snapshot IS the close by construction for both tables."""
    _init_dbs(tmp_path, monkeypatch)
    st = _state()                                # kickoff 2026-12-31T17:00
    after_kickoff = datetime.datetime(2026, 12, 31, 18, 0, tzinfo=datetime.timezone.utc)
    runner.run_tick(Soccer(), st, [_KEV], now=after_kickoff, dry_run=True,
                    book_fn=lambda t: _BOOK)
    with storage.connect(config.MARKET_DB_PATH, read_only=True) as c:
        assert c.execute("SELECT count(*) FROM line_snapshots").fetchone()[0] == 0
        assert c.execute("SELECT count(*) FROM book_snapshots").fetchone()[0] == 0


def test_snapshot_rows_carry_modified_on(tmp_path, monkeypatch):
    """The feed-side modifiedOn timestamp is persisted (staleness gate + clean CLV need it)."""
    _init_dbs(tmp_path, monkeypatch)
    st = _state()
    for k in list(st.lines):
        st.lines[k]["modifiedOn"] = "2026-06-01T11:59:00"
    runner.run_tick(Soccer(), st, [_KEV], now=_NOW, dry_run=True, book_fn=lambda t: _BOOK)
    with storage.connect(config.MARKET_DB_PATH, read_only=True) as c:
        vals = {r[0] for r in c.execute("SELECT DISTINCT modified_on FROM line_snapshots").fetchall()}
    assert vals == {"2026-06-01T11:59:00"}


def test_both_sides_flagged_logs_crossed_book_warning(tmp_path, monkeypatch, caplog):
    """yes and no on the same rung both +EV = crossed/stale book: must WARN."""
    import logging
    _init_dbs(tmp_path, monkeypatch)
    # 0.30 ask for yes AND no => yes_ask+no_ask=0.60 << 1: crossed book
    with caplog.at_level(logging.WARNING, logger="unabated_edge"):
        rows = runner.run_tick(Soccer(), _state(), [_KEV], now=_NOW, dry_run=True,
                               book_fn=lambda t: _BOOK)
    assert len(rows) == 2                        # both sides flagged (dry-run records them)
    assert any("crossed/stale book" in r.getMessage() for r in caplog.records)


def test_null_price_lines_not_snapshotted(tmp_path, monkeypatch):
    """A non-blurred line with no price must not be written as a NULL-price snapshot row."""
    _init_dbs(tmp_path, monkeypatch)
    st = _state()
    # inject a priceless (but non-blurred) line directly into state.lines
    st.lines["1|0|68|bt3"] = {"marketSourceId": 68}   # no americanPrice/price key
    runner.run_tick(Soccer(), st, [_KEV], now=_NOW, dry_run=True,
                    book_fn=lambda t: _BOOK)
    with storage.connect(config.MARKET_DB_PATH, read_only=True) as c:
        n_null = c.execute("SELECT count(*) FROM line_snapshots WHERE price IS NULL").fetchone()[0]
        n_total = c.execute("SELECT count(*) FROM line_snapshots").fetchone()[0]
    assert n_null == 0          # the priceless line was skipped, not stored as NULL
    assert n_total >= 1         # real lines still snapshotted


def test_book_snapshot_written_per_market(tmp_path, monkeypatch):
    """Every paired pre-kickoff market gets one book_snapshots row per tick with
    top-of-book columns, market volume/OI, and the full depth ladder."""
    _init_dbs(tmp_path, monkeypatch)
    runner.run_tick(Soccer(), _state(), [_KEV], now=_NOW, dry_run=True,
                    book_fn=lambda t: _BOOK)
    with storage.connect(config.MARKET_DB_PATH, read_only=True) as c:
        row = c.execute("""SELECT market_ticker, floor_strike, yes_bid, yes_bid_qty,
                                  yes_ask, no_ask, volume, open_interest,
                                  json_extract(depth, '$.no_bids[0][1]')
                           FROM book_snapshots""").fetchone()
    assert row[0] == "T-O25" and row[1] == 2.5
    assert row[2] == 0.70 and row[3] == 5.0
    assert row[4] == 0.30 and row[5] == 0.30
    assert row[6] == 1234.0 and row[7] == 567.0
    assert float(row[8]) == 8.0


def test_failed_book_fetch_skips_market_not_tick(tmp_path, monkeypatch):
    """book_fn=None for a market → no snapshot row, no candidates for it, no crash."""
    _init_dbs(tmp_path, monkeypatch)
    rows = runner.run_tick(Soccer(), _state(), [_KEV], now=_NOW, dry_run=True,
                           book_fn=lambda t: None)
    assert rows == []
    with storage.connect(config.MARKET_DB_PATH, read_only=True) as c:
        assert c.execute("SELECT count(*) FROM book_snapshots").fetchone()[0] == 0


def test_trades_fn_rows_inserted(tmp_path, monkeypatch):
    """When the trades poll rides a tick, tape rows land in kalshi_trades."""
    _init_dbs(tmp_path, monkeypatch)
    trades = [{"trade_id": "t1", "market_ticker": "T-O25",
               "created_time": "2026-06-01T11:59:59Z", "yes_price": 0.31,
               "count": 4.0, "taker_side": "yes"}]
    runner.run_tick(Soccer(), _state(), [_KEV], now=_NOW, dry_run=True,
                    book_fn=lambda t: _BOOK, trades_fn=lambda t: trades)
    with storage.connect(config.MARKET_DB_PATH, read_only=True) as c:
        got = c.execute("SELECT trade_id, taker_side FROM kalshi_trades").fetchall()
    assert got == [("t1", "yes")]


# ---------- overround feed-integrity gate (issue #73) ----------

def _poisoned_state():
    """Transient one-sided-refresh shape: main 2.5 healthy (sum 1.048), alt 1.5
    crossed (sum ~0.82), alt 3.5 blown vig (sum ~1.35) — all in one ladder."""
    return feed.parse_snapshot({
        "marketSources": [{"id": 7, "name": "S"}],
        "teams": {"1": {"name": "Argentina"}, "2": {"name": "Austria"}},
        "gameOddsEvents": {"lg21:pt1:pregame": [{
            "eventId": 1, "eventStart": "2026-12-31T17:00:00+00:00",
            "eventTeams": {"1": {"id": 1}, "0": {"id": 2}},
            "gameOddsMarketSourcesLines": {
                "si0:ms7:an0": {"bt3": {"price": 115, "points": 2.5,
                                        "modifiedOn": "2026-06-01T00:00:00",
                                        "alternateLines": [
                    {"points": 1.5, "price": 144}, {"points": 3.5, "price": -208}]}},
                "si1:ms7:an0": {"bt3": {"price": -140, "points": 2.5,
                                        "modifiedOn": "2026-06-01T00:00:00",
                                        "alternateLines": [
                    {"points": 1.5, "price": 144}, {"points": 3.5, "price": -208}]}},
            }}]}}, {"lg21"})


_POISON_KEV = {
    "title": "Argentina vs Austria: Regulation Time Total Goals",
    "markets": [
        {"ticker": "T-O15", "strike_type": "greater", "floor_strike": 1.5},
        {"ticker": "T-O25", "strike_type": "greater", "floor_strike": 2.5},
        {"ticker": "T-O35", "strike_type": "greater", "floor_strike": 3.5},
    ],
}


def test_poisoned_rungs_never_flag_healthy_rung_still_does(tmp_path, monkeypatch):
    """The exploit end-to-end: crossed and blown-vig rungs must produce NO
    candidate and NO flagged edge, while the healthy rung in the same ladder
    still flows. Fails on pre-#73 main (both poisoned rungs devig to ~0.5 and
    flag against the cheap 0.30 ask)."""
    _init_dbs(tmp_path, monkeypatch)
    rows = runner.run_tick(Soccer(), _poisoned_state(), [_POISON_KEV], now=_NOW,
                           dry_run=True, book_fn=lambda t: _BOOK)
    assert rows and {r["market_ticker"] for r in rows} == {"T-O25"}
    storage.flush()
    with storage.connect(config.RESEARCH_DB_PATH, read_only=True) as c:
        cand_tickers = {r[0] for r in c.execute(
            "SELECT json_extract_string(payload, '$.label') FROM research_events "
            "WHERE event_type = 'candidate_priced'").fetchall()}
        rejects = {(r[0], float(r[1])) for r in c.execute(
            "SELECT json_extract_string(payload, '$.reason'), "
            "json_extract_string(payload, '$.line') FROM research_events "
            "WHERE event_type = 'rung_rejected'").fetchall()}
    assert all("1.5" not in lbl and "3.5" not in lbl for lbl in cand_tickers)
    assert ("anchor_crossed", 1.5) in rejects
    assert ("anchor_overround", 3.5) in rejects


def test_rung_rejected_event_carries_provenance(tmp_path, monkeypatch):
    """The firehose reject rides the same provenance as candidate_priced:
    event_id, line, book, alt, overround, reason."""
    _init_dbs(tmp_path, monkeypatch)
    runner.run_tick(Soccer(), _poisoned_state(), [_POISON_KEV], now=_NOW,
                    dry_run=True, book_fn=lambda t: _BOOK)
    storage.flush()
    with storage.connect(config.RESEARCH_DB_PATH, read_only=True) as c:
        row = c.execute(
            "SELECT json_extract_string(payload, '$.reason'), "
            "json_extract(payload, '$.event_id'), json_extract(payload, '$.book'), "
            "json_extract(payload, '$.alt'), json_extract(payload, '$.overround') "
            "FROM research_events WHERE event_type = 'rung_rejected' "
            "AND json_extract_string(payload, '$.line') = '1.5' LIMIT 1").fetchone()
    assert row is not None
    reason, event_id, book, alt, overround = row
    assert reason == "anchor_crossed"
    assert int(event_id) == 1 and int(book) == 7
    assert alt == "true"
    assert 0.81 < float(overround) < 0.83


def test_poisoned_rungs_never_quoted_by_real_maker(tmp_path, monkeypatch):
    """Same exploit through the REAL maker path (run_tick -> fair_ladder ->
    MakerEngine): the healthy rung quotes, the crossed and blown-vig rungs
    never do."""
    from unabated_edge.maker import engine as mengine, state as mstate, store as mstore
    from unabated_edge.tests.test_maker_engine import FakeGateway
    _init_dbs(tmp_path, monkeypatch)
    monkeypatch.setattr(config, "MAKER_DB_PATH", tmp_path / "mk.duckdb")
    mstore.init()
    mk = mengine.MakerEngine(FakeGateway(), mstate.MakerState())
    wide_book = {"yes_bids": [(0.30, 500.0)], "no_bids": [(0.45, 500.0)]}
    runner.run_tick(Soccer(), _poisoned_state(), [_POISON_KEV], now=_NOW,
                    dry_run=True, book_fn=lambda t: wide_book, maker=mk)
    placed_tickers = {t for t, _s, _p, _n in mk.gateway.placed}
    assert placed_tickers == {"T-O25"}


def test_devig_never_receives_crossed_input(tmp_path, monkeypatch):
    """Boundary proof, not code-reading: spy on pricing.devig through the
    shipped path and assert no input pair ever sums below 1."""
    from unabated_edge import pricing
    _init_dbs(tmp_path, monkeypatch)
    real_devig, sums = pricing.devig, []
    def spy(probs):
        sums.append(sum(probs))
        return real_devig(probs)
    monkeypatch.setattr(pricing, "devig", spy)
    runner.run_tick(Soccer(), _poisoned_state(), [_POISON_KEV], now=_NOW,
                    dry_run=True, book_fn=lambda t: _BOOK)
    assert sums                                  # the healthy rung was devigged
    assert all(s >= 1.0 for s in sums)


class _FakeMaker:
    def __init__(self):
        self.calls, self.sweeps = [], []
    def on_match(self, adapter, event_meta, kev, ladder, books, now):
        self.calls.append((event_meta.event_id, ladder is not None, sorted(books)))
    def sweep(self, sport, seen, now):
        self.sweeps.append((sport, set(seen)))


def test_maker_hook_receives_ladder_and_books(tmp_path, monkeypatch):
    _init_dbs(tmp_path, monkeypatch)
    mk = _FakeMaker()
    runner.run_tick(Soccer(), _state(), [_KEV], now=_NOW, dry_run=True,
                    book_fn=lambda t: _BOOK, maker=mk)
    assert mk.calls == [(1, True, ["T-O25"])]
    assert mk.sweeps == [("soccer", {1})]


def test_maker_hook_skipped_when_none(tmp_path, monkeypatch):
    _init_dbs(tmp_path, monkeypatch)
    rows = runner.run_tick(Soccer(), _state(), [_KEV], now=_NOW, dry_run=True,
                           book_fn=lambda t: _BOOK)     # maker defaults to None
    assert any(rows)


def test_kill_switch_pulls_maker_quotes(tmp_path, monkeypatch):
    """The kill file is the emergency stop: it must cancel resting maker orders,
    not just stop the taker tick."""
    _init_dbs(tmp_path, monkeypatch)
    kill = tmp_path / ".kill"
    kill.touch()
    monkeypatch.setattr(config, "KILL_FILE", kill)
    pulls = []
    class _PullRecorder(_FakeMaker):
        def pull_all(self, now, reason):
            pulls.append(reason)
    mk = _PullRecorder()
    rows = runner.run_tick(Soccer(), _state(), [_KEV], now=_NOW, dry_run=True,
                           book_fn=lambda t: _BOOK, maker=mk)
    assert rows == [] and pulls == ["kill_switch"]
    assert mk.calls == []            # no quoting happened under the kill switch


def test_main_loop_shutdown_pulls_maker_quotes(tmp_path, monkeypatch):
    """SIGINT/SIGTERM (cleared _running) must pull resting maker quotes on exit."""
    _init_dbs(tmp_path, monkeypatch)
    pulls = []
    class _GW:
        is_live = False
    class _MK:
        def pull_all(self, now, reason):
            pulls.append(reason)
        def stats(self):
            return {}
    monkeypatch.setattr(runner.kalshi, "init", lambda: None)
    monkeypatch.setattr(runner.maker_gateway, "make_gateway", lambda mode, ack: _GW())
    monkeypatch.setattr(runner.maker_engine, "MakerEngine", lambda gw, st: _MK())
    monkeypatch.setattr(runner.maker_store, "init", lambda: None)
    runner._running.clear()          # simulate signal received before first tick
    try:
        runner.main_loop(dry_run=True)
    finally:
        runner._running.set()        # restore for other tests
    assert pulls == ["shutdown"]


def test_main_loop_live_startup_sync(tmp_path, monkeypatch):
    """Live mode must run startup_sync (baseline positions + cancel orphan
    orders) before the first tick, so a restart with inventory or stray
    resting orders doesn't false-trip the fresh cap stack / mismatch tripwire."""
    _init_dbs(tmp_path, monkeypatch)

    class _GW:
        is_live = True

    class _MK:
        def __init__(self):
            self.state = object()
        def pull_all(self, now, reason):
            pass
        def stats(self):
            return {}

    calls = []

    def fake_startup_sync(state, gw, prefixes):
        calls.append((state, gw, prefixes))

    monkeypatch.setattr(runner.kalshi, "init", lambda: None)
    monkeypatch.setattr(runner.maker_gateway, "make_gateway", lambda mode, ack: _GW())
    monkeypatch.setattr(runner.maker_engine, "MakerEngine", lambda gw, st: _MK())
    monkeypatch.setattr(runner.maker_store, "init", lambda: None)
    monkeypatch.setattr(runner.maker_state, "startup_sync", fake_startup_sync)
    runner._running.clear()          # simulate signal received before first tick
    try:
        runner.main_loop(dry_run=True)
    finally:
        runner._running.set()        # restore for other tests
    assert len(calls) == 1
    state, gw, prefixes = calls[0]
    assert isinstance(gw, _GW)
    assert "KXWCTOTAL" in prefixes


def test_singleton_lock_blocks_second_instance(tmp_path, monkeypatch):
    """A second runner must refuse to start while our own pid's lock is live,
    but reclaim a stale lock left by a dead pid (crash / kill -9)."""
    lock_path = tmp_path / ".runner.lock"
    monkeypatch.setattr(runner, "_LOCK_PATH", lock_path)
    lock_path.write_text(str(os.getpid()))          # our own pid: definitely alive
    try:
        with pytest.raises(SystemExit):
            runner._acquire_singleton_lock()
    finally:
        lock_path.unlink(missing_ok=True)
    lock_path.write_text("999999")                   # a pid that should not be alive
    runner._acquire_singleton_lock()                 # must NOT raise: stale lock reclaimed
    assert lock_path.read_text().strip() == str(os.getpid())
    lock_path.unlink(missing_ok=True)
