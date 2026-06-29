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
         "yes_sub_title": "Reg Time: Over 2.5 goals scored"},
    ],
}


def _init_dbs(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "MARKET_DB_PATH", tmp_path / "m.duckdb")
    monkeypatch.setattr(config, "RESEARCH_DB_PATH", tmp_path / "r.duckdb")
    storage.init()


def test_tick_flags_positive_edge(tmp_path, monkeypatch):
    _init_dbs(tmp_path, monkeypatch)
    # ask_fn now takes (ticker, side); 0.30 is well below a ~42% fair Over → positive EV
    rows = runner.run_tick(Soccer(), _state(), [_KEV], now=_NOW, dry_run=True,
                           ask_fn=lambda t, side: 0.30)
    assert any(r["ev_pct"] > 0 for r in rows)


def test_absolute_ev_floor_blocks_flag(tmp_path, monkeypatch):
    """An edge that clears ev_pct must still be blocked by the absolute-dollar floor."""
    _init_dbs(tmp_path, monkeypatch)
    monkeypatch.setattr(config, "MIN_EV_DOLLARS", 1.0)   # no real contract clears $1 EV
    rows = runner.run_tick(Soccer(), _state(), [_KEV], now=_NOW, dry_run=True,
                           ask_fn=lambda t, side: 0.30)
    assert rows == []


def test_null_price_lines_not_snapshotted(tmp_path, monkeypatch):
    """A non-blurred line with no price must not be written as a NULL-price snapshot row."""
    _init_dbs(tmp_path, monkeypatch)
    st = _state()
    # inject a priceless (but non-blurred) line directly into state.lines
    st.lines["1|0|68|bt3"] = {"marketSourceId": 68}   # no americanPrice/price key
    runner.run_tick(Soccer(), st, [_KEV], now=_NOW, dry_run=True,
                    ask_fn=lambda t, side: 0.30)
    with storage.connect(config.MARKET_DB_PATH, read_only=True) as c:
        n_null = c.execute("SELECT count(*) FROM line_snapshots WHERE price IS NULL").fetchone()[0]
        n_total = c.execute("SELECT count(*) FROM line_snapshots").fetchone()[0]
    assert n_null == 0          # the priceless line was skipped, not stored as NULL
    assert n_total >= 1         # real lines still snapshotted
