"""Dedicated tests for build_market_id() — the single canonical constructor."""

import pytest
from nfl_draft.lib.market_map import build_market_id, slug


def test_slug_lowercases():
    assert slug("Cam Ward") == "cam-ward"


def test_slug_strips_periods():
    assert slug("J.T. Smith") == "jt-smith"


def test_slug_strips_apostrophes():
    assert slug("J'shawn O'Brien") == "jshawn-obrien"


def test_each_market_type_produces_expected_id():
    cases = [
        ("first_at_position", {"position": "QB", "player": "Cam Ward"}, "first_qb_cam-ward"),
        ("pick_outright", {"pick_number": 1, "player": "Cam Ward"}, "pick_1_overall_cam-ward"),
        ("top_n_range", {"range_high": 5, "player": "Ashton Jeanty"}, "top_5_ashton-jeanty"),
        ("team_first_pick", {"team": "CHI", "player": "Ashton Jeanty"}, "team_chi_first_pick_ashton-jeanty"),
        ("prop", {"short_description": "Will a kicker go R1"}, "prop_will-a-kicker-go-r1"),
    ]
    for market_type, kwargs, expected in cases:
        assert build_market_id(market_type, **kwargs) == expected


def test_unknown_market_type_raises():
    with pytest.raises(ValueError):
        build_market_id("invalid_type", foo="bar")
