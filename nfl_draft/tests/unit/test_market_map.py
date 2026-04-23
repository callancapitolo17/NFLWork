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


def test_build_market_id_draft_position_over_under_sub_one_line():
    # 0.5 lines are rare but valid — ensure the 0 integer part is preserved
    # (f"{0.5:g}" yields '0.5', not '.5').
    assert build_market_id(
        "draft_position_over_under", player="Mr Irrelevant", line=0.5, direction="under"
    ) == "draft_position_ou_mr-irrelevant_0p5_under"


from nfl_draft.lib.market_map import outright_group_key


def test_outright_group_key_pick_outright():
    assert outright_group_key("pick_outright", "pick_2_overall_david-bailey") == "pick_2_overall"


def test_outright_group_key_first_at_position():
    assert outright_group_key("first_at_position", "first_wr_jordyn-tyson") == "first_wr"


def test_outright_group_key_top_n_range():
    assert outright_group_key("top_10_range", "top_10_sonny-styles") == "top_10"
    assert outright_group_key("top_5_range", "top_5_caleb-downs") == "top_5"


def test_outright_group_key_team_first_pick_handles_spaces_in_team():
    # build_market_id uses team.lower() -> can contain a space ("new england")
    assert outright_group_key("team_first_pick", "team_washington_first_pick_lucas") == "team_washington_first_pick"
    assert outright_group_key("team_first_pick", "team_new england_first_pick_drew-allar") == "team_new england_first_pick"


def test_outright_group_key_team_first_pick_position():
    # build_market_id uses _slug_underscored for both team and position
    assert outright_group_key(
        "team_first_pick_position",
        "arizona_cardinals_first_pick_pos_wide_receiver",
    ) == "arizona_cardinals_first_pick_pos"


def test_outright_group_key_nth_at_position():
    # market_group = "nth_at_position_2"; market_id = "2_wr_jordyn-tyson"
    assert outright_group_key("nth_at_position_2", "2_wr_jordyn-tyson") == "2_wr"
    assert outright_group_key("nth_at_position_3", "3_cb_jermod-mccoy") == "3_cb"


def test_outright_group_key_draft_position_over_under():
    assert outright_group_key(
        "draft_position_over_under",
        "draft_position_ou_spencer-fano_10p5_over",
    ) == "draft_position_ou_spencer-fano_10p5"
    assert outright_group_key(
        "draft_position_over_under",
        "draft_position_ou_spencer-fano_10p5_under",
    ) == "draft_position_ou_spencer-fano_10p5"


def test_outright_group_key_mr_irrelevant_returns_none():
    assert outright_group_key("mr_irrelevant_position", "mr_irrelevant_wide_receiver") is None


def test_outright_group_key_matchup_before_returns_none():
    # Matchups are self-contained 2-way; we don't bucket them across other markets.
    assert outright_group_key("matchup_before", "matchup_hunter_before_tyson") is None


def test_outright_group_key_prop_returns_none():
    assert outright_group_key("prop_first_round_total_ou", "prop_first_round_total") is None
    assert outright_group_key("prop_team_position_of_first_pick", "prop_xxx") is None


def test_outright_group_key_unrecognized_returns_none():
    assert outright_group_key("", "random_market_id") is None
    assert outright_group_key("totally_unknown", "x_y_z") is None


def test_outright_group_key_malformed_market_id_returns_none():
    # market_group says pick_outright but market_id doesn't match the pattern
    assert outright_group_key("pick_outright", "not_a_pick_market") is None
