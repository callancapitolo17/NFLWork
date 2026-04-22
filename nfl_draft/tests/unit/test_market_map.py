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


def test_build_market_id_nth_at_position():
    assert build_market_id(
        "nth_at_position", nth=2, position="WR", player="Jordyn Tyson"
    ) == "2_wr_jordyn-tyson"


def test_build_market_id_mr_irrelevant_position():
    assert build_market_id(
        "mr_irrelevant_position", position="Wide Receiver"
    ) == "mr_irrelevant_wide_receiver"


def test_build_market_id_team_first_pick_position():
    assert build_market_id(
        "team_first_pick_position", team="Arizona Cardinals", position="Wide Receiver"
    ) == "arizona_cardinals_first_pick_pos_wide_receiver"


def test_build_market_id_matchup_before_canonical_order():
    # Order-independent: whichever player is listed first, we get the same ID.
    a = build_market_id("matchup_before", player_a="Makai Lemon", player_b="Kenyon Sadiq")
    b = build_market_id("matchup_before", player_a="Kenyon Sadiq", player_b="Makai Lemon")
    assert a == b == "matchup_kenyon-sadiq_before_makai-lemon"


def test_build_market_id_draft_position_over_under_half_point():
    assert build_market_id(
        "draft_position_over_under", player="Caleb Downs", line=9.5, direction="under"
    ) == "draft_position_ou_caleb-downs_9p5_under"


def test_build_market_id_draft_position_over_under_whole_point():
    # Whole-point lines drop trailing .0.
    assert build_market_id(
        "draft_position_over_under", player="Abdul Carter", line=3.0, direction="over"
    ) == "draft_position_ou_abdul-carter_3_over"
