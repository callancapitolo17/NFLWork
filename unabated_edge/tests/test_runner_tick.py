import datetime
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
