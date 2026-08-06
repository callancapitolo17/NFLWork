import datetime
import logging
import os
import types
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
                        # 3h after _NOW: comfortably inside the default
                        # BOOK_CAPTURE_HORIZON_HOURS (12) so these fixtures
                        # still get a book/trades fetch, and past
                        # tipoff_ok's 3-minute cutoff.
                        "eventStart": "2026-06-01T03:00:00+00:00",
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
    st = _state()                                # kickoff 2026-06-01T03:00
    after_kickoff = datetime.datetime(2026, 6, 1, 4, 0, tzinfo=datetime.timezone.utc)
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
    monkeypatch.setattr(config, "MAKER_DB_PATH", tmp_path / "mk.duckdb")  # restore_halt_latch needs a real maker DB
    pulls = []
    class _GW:
        is_live = False
    class _MK:
        def __init__(self):
            self.state = runner.maker_state.MakerState()
        def pull_all(self, now, reason):
            pulls.append(reason)
        def stats(self):
            return {}
    monkeypatch.setattr(runner.kalshi, "init", lambda: None)
    monkeypatch.setattr(runner.maker_gateway, "make_gateway", lambda mode, ack: _GW())
    monkeypatch.setattr(runner.maker_engine, "MakerEngine", lambda gw, st: _MK())
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
    monkeypatch.setattr(config, "MAKER_DB_PATH", tmp_path / "mk.duckdb")  # restore_halt_latch needs a real maker DB

    class _GW:
        is_live = True

    class _MK:
        def __init__(self):
            self.state = runner.maker_state.MakerState()
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


def test_book_capture_horizon_skips_far_future_events(tmp_path, monkeypatch):
    """Finding 1d: an event more than BOOK_CAPTURE_HORIZON_HOURS out gets no
    Kalshi book/trades fetch this tick (still gets a line_snapshots row from
    the earlier per-line loop, which doesn't gate on the horizon) — this is
    what shortens a full MLB slate tick and caps capture volume for games
    nobody can trade yet."""
    _init_dbs(tmp_path, monkeypatch)
    monkeypatch.setattr(config, "BOOK_CAPTURE_HORIZON_HOURS", 12.0)
    near_start = (_NOW + datetime.timedelta(hours=3)).isoformat()
    far_start = (_NOW + datetime.timedelta(hours=48)).isoformat()
    state = feed.parse_snapshot(
        {
            "marketSources": [{"id": 7, "name": "S"}],
            "teams": {"1": {"name": "Argentina"}, "2": {"name": "Austria"},
                      "3": {"name": "Brazil"}, "4": {"name": "Chile"}},
            "gameOddsEvents": {
                "lg21:pt1:pregame": [
                    {
                        "eventId": 1, "eventStart": near_start,
                        "eventTeams": {"1": {"id": 1}, "0": {"id": 2}},
                        "gameOddsMarketSourcesLines": {
                            "si0:ms7:an0": {"bt3": {"price": 115, "points": 2.5}},
                            "si1:ms7:an0": {"bt3": {"price": -140, "points": 2.5}},
                        },
                    },
                    {
                        "eventId": 2, "eventStart": far_start,
                        "eventTeams": {"1": {"id": 3}, "0": {"id": 4}},
                        "gameOddsMarketSourcesLines": {
                            "si0:ms7:an0": {"bt3": {"price": 115, "points": 2.5}},
                            "si1:ms7:an0": {"bt3": {"price": -140, "points": 2.5}},
                        },
                    },
                ]
            },
        },
        {"lg21"},
    )
    near_kev = {"title": "Argentina vs Austria: Regulation Time Total Goals",
                "markets": [{"ticker": "NEAR-O25", "strike_type": "greater", "floor_strike": 2.5}]}
    far_kev = {"title": "Brazil vs Chile: Regulation Time Total Goals",
               "markets": [{"ticker": "FAR-O25", "strike_type": "greater", "floor_strike": 2.5}]}
    called = []

    def book_fn(t):
        called.append(t)
        return _BOOK

    runner.run_tick(Soccer(), state, [near_kev, far_kev], now=_NOW, dry_run=True, book_fn=book_fn)
    assert called == ["NEAR-O25"]                    # far event never hit the book_fn at all
    with storage.connect(config.MARKET_DB_PATH, read_only=True) as c:
        tickers = {r[0] for r in c.execute("SELECT market_ticker FROM book_snapshots").fetchall()}
        n_events_snapshotted = c.execute(
            "SELECT count(DISTINCT event_id) FROM line_snapshots").fetchone()[0]
    assert tickers == {"NEAR-O25"}
    assert n_events_snapshotted == 2                 # both events still line-snapshotted


class _SeqClock:
    """Returns real datetime objects from a fixed, pre-scripted sequence —
    used to prove main_loop recomputes `now` per adapter instead of reusing
    one shared iteration-start value."""
    def __init__(self, values):
        self._values = list(values)

    def now(self, tz=None):
        return self._values.pop(0)


def test_note_success_stamped_fresh_after_each_adapter_tick(tmp_path, monkeypatch):
    """Finding 1a/1b: main_loop must (a) recompute `now` fresh right before
    each adapter's run_tick, not reuse a shared iteration-start value, and
    (b) stamp note_success with an even fresher timestamp taken AFTER that
    adapter's own tick completes. On a real slate the MLB book pass alone
    runs ~37s — if note_success used the iteration-start clock, measured
    staleness would include that whole duration and the watchdog would
    self-trigger every tick once quotes rest."""
    _init_dbs(tmp_path, monkeypatch)

    t0 = datetime.datetime(2026, 7, 25, 12, 0, 0, tzinfo=datetime.timezone.utc)
    # 9 datetime.now() calls happen in one no-op iteration of main_loop:
    # startup prune-date, maker setup's restore_halt_latch (finding #75/F6),
    # in-loop prune-date check, then per adapter (before, after) x2, then
    # watchdog, then the shutdown pull_all.
    times = [
        t0,                                          # 1: last_prune_date (startup)
        t0,                                          # 2: restore_halt_latch (maker setup, pre-loop)
        t0,                                          # 3: today (in-loop, same day -> no re-prune)
        t0,                                          # 4: soccer — now (before run_tick)
        t0 + datetime.timedelta(seconds=25),          # 5: soccer — note_success (after)
        t0 + datetime.timedelta(seconds=25),          # 6: mlb — now (before run_tick), fresh
        t0 + datetime.timedelta(seconds=25 + 37),     # 7: mlb — note_success (after)
        t0 + datetime.timedelta(seconds=25 + 37 + 1),  # 8: watchdog
        t0 + datetime.timedelta(seconds=25 + 37 + 1),  # 9: shutdown pull_all
    ]
    fake_datetime = types.SimpleNamespace(
        datetime=types.SimpleNamespace(now=_SeqClock(times).now),
        timezone=datetime.timezone,
        timedelta=datetime.timedelta,
    )
    monkeypatch.setattr(runner, "datetime", fake_datetime)

    run_tick_calls = []

    def fake_run_tick(adapter, state, kalshi_events, *, now, dry_run, book_fn, trades_fn=None, maker=None):
        run_tick_calls.append((adapter.sport, now))
        return []
    monkeypatch.setattr(runner, "run_tick", fake_run_tick)
    monkeypatch.setattr(runner.feed, "fetch_v2", lambda *a, **k: feed.FeedState())
    monkeypatch.setattr(runner.kalshi, "init", lambda: None)
    monkeypatch.setattr(runner.kalshi, "list_events", lambda series: [])

    note_success_calls = []
    monkeypatch.setattr(config, "MAKER_DB_PATH", tmp_path / "mk.duckdb")  # restore_halt_latch needs a real maker DB

    class _FakeMaker:
        def __init__(self):
            self.state = runner.maker_state.MakerState()
        def note_success(self, sport, now):
            note_success_calls.append((sport, now))
        def watchdog(self, now):
            pass
        def pull_all(self, now, reason):
            pass
        def stats(self):
            return {}

    class _GW:
        is_live = False

    monkeypatch.setattr(runner.maker_gateway, "make_gateway", lambda mode, ack: _GW())
    monkeypatch.setattr(runner.maker_engine, "MakerEngine", lambda gw, st: _FakeMaker())
    monkeypatch.setattr(runner.time, "sleep", lambda _s: runner._running.clear())

    runner._running.set()
    try:
        runner.main_loop(dry_run=True)
    finally:
        runner._running.set()

    assert [sport for sport, _ in run_tick_calls] == ["soccer", "mlb"]
    assert [sport for sport, _ in note_success_calls] == ["soccer", "mlb"]

    soccer_now, mlb_now = run_tick_calls[0][1], run_tick_calls[1][1]
    soccer_stamp, mlb_stamp = note_success_calls[0][1], note_success_calls[1][1]

    assert soccer_now == t0
    assert soccer_stamp == t0 + datetime.timedelta(seconds=25)
    # mlb's pre-tick clock equals soccer's post-tick stamp — a fresh
    # recompute, not a copy of the iteration-start t0.
    assert mlb_now == soccer_stamp
    assert mlb_stamp == t0 + datetime.timedelta(seconds=25 + 37)
    # note_success is stamped exactly the simulated tick duration after the
    # `now` run_tick was handed — proving it's a fresh post-tick timestamp.
    assert (soccer_stamp - soccer_now).total_seconds() == 25
    assert (mlb_stamp - mlb_now).total_seconds() == 37


# ---------- frozen-feed detection (Finding #77) ----------

def test_feed_advanced_true_on_first_observation():
    """A sport never seen before must not false-trip the frozen check."""
    freeze_state = {}
    now = datetime.datetime(2026, 8, 1, tzinfo=datetime.timezone.utc)
    assert runner._feed_advanced(freeze_state, "mlb", "sigA", now) is True
    assert freeze_state["mlb"] == ("sigA", now, False)


def test_feed_advanced_true_when_signature_changes():
    since = datetime.datetime(2026, 8, 1, tzinfo=datetime.timezone.utc)
    freeze_state = {"mlb": ("sigA", since, False)}
    later = since + datetime.timedelta(minutes=20)
    assert runner._feed_advanced(freeze_state, "mlb", "sigB", later) is True
    assert freeze_state["mlb"] == ("sigB", later, False)   # timer resets on the new signature


def test_feed_advanced_tolerates_repeat_under_threshold():
    """A calm pre-game market can legitimately serve the same signature for a
    while -- must not trip before FEED_FROZEN_SEC, and must not reset the
    'since' clock on every repeat poll (that would mean it never trips)."""
    since = datetime.datetime(2026, 8, 1, tzinfo=datetime.timezone.utc)
    freeze_state = {"mlb": ("sigA", since, False)}
    still_within = since + datetime.timedelta(seconds=config.FEED_FROZEN_SEC - 1)
    assert runner._feed_advanced(freeze_state, "mlb", "sigA", still_within) is True
    assert freeze_state["mlb"] == ("sigA", since, False)


def test_feed_advanced_false_when_frozen_past_threshold(caplog):
    """The exact same payload served for longer than FEED_FROZEN_SEC is a
    frozen CDN/origin, not a calm market -- must fail closed."""
    since = datetime.datetime(2026, 8, 1, tzinfo=datetime.timezone.utc)
    freeze_state = {"mlb": ("sigA", since, False)}
    past_threshold = since + datetime.timedelta(seconds=config.FEED_FROZEN_SEC + 1)
    with caplog.at_level(logging.WARNING, logger="unabated_edge"):
        result = runner._feed_advanced(freeze_state, "mlb", "sigA", past_threshold)
    assert result is False
    assert any("frozen" in r.getMessage() and "mlb" in r.getMessage() for r in caplog.records)


def test_feed_advanced_warns_once_per_freeze_episode(caplog):
    """Once frozen, the tick loop polls every V2_POLL_SEC (~5s) -- the WARNING
    must not repeat on every one of those ticks while the freeze persists."""
    since = datetime.datetime(2026, 8, 1, tzinfo=datetime.timezone.utc)
    freeze_state = {"mlb": ("sigA", since, False)}
    past_threshold = since + datetime.timedelta(seconds=config.FEED_FROZEN_SEC + 1)
    with caplog.at_level(logging.WARNING, logger="unabated_edge"):
        runner._feed_advanced(freeze_state, "mlb", "sigA", past_threshold)
        still_frozen = past_threshold + datetime.timedelta(seconds=5)
        result = runner._feed_advanced(freeze_state, "mlb", "sigA", still_frozen)
    assert result is False
    frozen_warnings = [r for r in caplog.records if "frozen" in r.getMessage()]
    assert len(frozen_warnings) == 1


def test_main_loop_withholds_note_success_when_feed_frozen(tmp_path, monkeypatch):
    """Wiring check: main_loop must gate note_success on _feed_advanced, so a
    frozen-but-200-OK feed lets MAX_STALENESS_SEC's watchdog trip and pull
    quotes instead of being fooled by note_success firing every tick."""
    _init_dbs(tmp_path, monkeypatch)
    monkeypatch.setattr(runner, "run_tick", lambda *a, **k: [])
    monkeypatch.setattr(runner.feed, "fetch_v2", lambda *a, **k: feed.FeedState())
    monkeypatch.setattr(runner, "_feed_advanced", lambda *a, **k: False)
    monkeypatch.setattr(runner.kalshi, "init", lambda: None)
    monkeypatch.setattr(runner.kalshi, "list_events", lambda series: [])

    note_success_calls = []
    monkeypatch.setattr(config, "MAKER_DB_PATH", tmp_path / "mk.duckdb")  # restore_halt_latch needs a real maker DB

    class _FakeMaker:
        def __init__(self):
            self.state = runner.maker_state.MakerState()
        def note_success(self, sport, now):
            note_success_calls.append(sport)
        def watchdog(self, now):
            pass
        def pull_all(self, now, reason):
            pass
        def stats(self):
            return {}

    class _GW:
        is_live = False

    monkeypatch.setattr(runner.maker_gateway, "make_gateway", lambda mode, ack: _GW())
    monkeypatch.setattr(runner.maker_engine, "MakerEngine", lambda gw, st: _FakeMaker())
    monkeypatch.setattr(runner.time, "sleep", lambda _s: runner._running.clear())

    runner._running.set()
    try:
        runner.main_loop(dry_run=True)
    finally:
        runner._running.set()

    assert note_success_calls == []


def test_main_loop_calls_note_success_when_feed_advancing(tmp_path, monkeypatch):
    """Regression guard for the wiring above: when the feed IS advancing,
    note_success must still fire normally."""
    _init_dbs(tmp_path, monkeypatch)
    monkeypatch.setattr(runner, "run_tick", lambda *a, **k: [])
    monkeypatch.setattr(runner.feed, "fetch_v2", lambda *a, **k: feed.FeedState())
    monkeypatch.setattr(runner, "_feed_advanced", lambda *a, **k: True)
    monkeypatch.setattr(runner.kalshi, "init", lambda: None)
    monkeypatch.setattr(runner.kalshi, "list_events", lambda series: [])

    note_success_calls = []
    monkeypatch.setattr(config, "MAKER_DB_PATH", tmp_path / "mk.duckdb")  # restore_halt_latch needs a real maker DB

    class _FakeMaker:
        def __init__(self):
            self.state = runner.maker_state.MakerState()
        def note_success(self, sport, now):
            note_success_calls.append(sport)
        def watchdog(self, now):
            pass
        def pull_all(self, now, reason):
            pass
        def stats(self):
            return {}

    class _GW:
        is_live = False

    monkeypatch.setattr(runner.maker_gateway, "make_gateway", lambda mode, ack: _GW())
    monkeypatch.setattr(runner.maker_engine, "MakerEngine", lambda gw, st: _FakeMaker())
    monkeypatch.setattr(runner.time, "sleep", lambda _s: runner._running.clear())

    runner._running.set()
    try:
        runner.main_loop(dry_run=True)
    finally:
        runner._running.set()

    assert sorted(note_success_calls) == ["mlb", "soccer"]


def test_adapter_exception_does_not_skip_other_adapter(tmp_path, monkeypatch, caplog):
    """Finding 4: one adapter's fetch/tick exception must not skip the other
    sport's tick. Soccer's feed fetch raises; MLB must still tick, and a
    WARNING naming the failed sport must be logged (the outer try/except
    stays as the last-resort net for the rest of the loop body)."""
    import logging
    _init_dbs(tmp_path, monkeypatch)

    def fake_fetch_v2(league_id, league_prefix):
        if league_prefix == "lg21":                  # soccer
            raise RuntimeError("feed down")
        return feed.FeedState()
    monkeypatch.setattr(runner.feed, "fetch_v2", fake_fetch_v2)

    run_tick_calls = []

    def fake_run_tick(adapter, state, kalshi_events, *, now, dry_run, book_fn, trades_fn=None, maker=None):
        run_tick_calls.append(adapter.sport)
        return []
    monkeypatch.setattr(runner, "run_tick", fake_run_tick)
    monkeypatch.setattr(runner.kalshi, "init", lambda: None)
    monkeypatch.setattr(runner.kalshi, "list_events", lambda series: [])
    monkeypatch.setattr(runner.maker_gateway, "make_gateway", lambda mode, ack: None)
    monkeypatch.setattr(runner.time, "sleep", lambda _s: runner._running.clear())

    runner._running.set()
    try:
        with caplog.at_level(logging.WARNING, logger="unabated_edge"):
            runner.main_loop(dry_run=True)
    finally:
        runner._running.set()

    assert run_tick_calls == ["mlb"]                  # soccer's exception didn't stop mlb's tick
    assert any(r.levelno == logging.WARNING and "soccer" in r.getMessage()
               for r in caplog.records)
