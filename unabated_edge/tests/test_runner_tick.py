import datetime
from unabated_edge import feed, runner, config, storage
from unabated_edge.sports.soccer import Soccer


def test_tick_flags_positive_edge(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "MARKET_DB_PATH", tmp_path / "m.duckdb")
    monkeypatch.setattr(config, "RESEARCH_DB_PATH", tmp_path / "r.duckdb")
    storage.init()
    st = feed.parse_snapshot(
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
                            "si1:ms7:an0": {"bt1": {"price": -200}, "bt4": {"price": 300}},
                            "si0:ms7:an0": {"bt1": {"price": 400}},
                        },
                    }
                ]
            },
        },
        {"lg21"},
    )
    kev = {
        "event_ticker": "E",
        "markets": [
            {"ticker": "T-ARG", "yes_sub_title": "Argentina"},
            {"ticker": "T-DRAW", "yes_sub_title": "Draw"},
            {"ticker": "T-AUT", "yes_sub_title": "Austria"},
        ],
    }
    rows = runner.run_tick(
        Soccer(),
        st,
        [kev],
        now=datetime.datetime(2026, 6, 1, tzinfo=datetime.timezone.utc),
        dry_run=True,
        ask_fn=lambda t: 0.30,
    )
    storage.flush()  # drain research buffer so it doesn't leak into subsequent tests
    assert any(r["ev_pct"] > 0 for r in rows)
