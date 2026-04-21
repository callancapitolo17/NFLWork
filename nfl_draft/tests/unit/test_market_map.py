import pytest
from nfl_draft.lib import db as db_module
from nfl_draft.lib import seed
from nfl_draft.lib.market_map import resolve_market_id, build_market_id


@pytest.fixture
def seeded(monkeypatch, tmp_path):
    monkeypatch.setattr(db_module, "DB_PATH", tmp_path / "test.duckdb")
    db_module.init_schema()
    seed.run()


def test_build_market_id_first_at_position():
    assert build_market_id("first_at_position", position="QB", player="Cam Ward") == "first_qb_cam-ward"


def test_build_market_id_pick_outright():
    assert build_market_id("pick_outright", pick_number=1, player="Cam Ward") == "pick_1_overall_cam-ward"


def test_build_market_id_top_n_range():
    assert build_market_id("top_n_range", range_high=5, player="Ashton Jeanty") == "top_5_ashton-jeanty"


def test_build_market_id_team_first_pick():
    assert build_market_id("team_first_pick", team="CHI", player="Ashton Jeanty") == "team_chi_first_pick_ashton-jeanty"


def test_build_market_id_strips_punctuation():
    assert build_market_id("first_at_position", position="QB", player="J.T. O'Brien") == "first_qb_jt-obrien"


def test_resolve_market_id_unknown_returns_none(seeded):
    assert resolve_market_id("kalshi", "UNKNOWN_LABEL", "Some Subject") is None
