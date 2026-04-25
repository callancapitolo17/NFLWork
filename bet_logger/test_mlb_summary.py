"""Unit tests for create_mlb_summary.py filter helpers."""
import pytest
from create_mlb_summary import (
    is_mlb_correlated_parlay,
    is_f5_parlay,
    SPREAD_REGEX,
    TOTAL_REGEX,
    F5_REGEX,
)


# Real Wagerzon description samples (format: "PARLAY (2 TEAMS) | leg1 | leg2")
SPREAD_PLUS_TOTAL_FG = (
    "PARLAY (2 TEAMS) | Detroit Tigers -1.5 +140 | Over 8.5 -110"
)
SPREAD_PLUS_TOTAL_F5 = (
    "PARLAY (2 TEAMS) | Detroit Tigers 1st 5 Innings -0.5 +115 | "
    "1st 5 Innings Over 4.5 -105"
)
ML_PLUS_ML = (
    "PARLAY (2 TEAMS) | Detroit Tigers ML +100 | Chicago Cubs ML -120"
)
TOTAL_PLUS_TOTAL = (
    "PARLAY (2 TEAMS) | Over 8.5 -110 | Over 4.5 -105 1st 5 Innings"
)
THREE_LEGS = (
    "PARLAY (3 TEAMS) | Detroit Tigers -1.5 +140 | Over 8.5 -110 | "
    "Yankees ML -200"
)
STRAIGHT_BET = "STRAIGHT BET | Detroit Tigers -1.5 +140"


@pytest.mark.parametrize(
    "desc, sport, bet_type, expected",
    [
        (SPREAD_PLUS_TOTAL_FG, "MLB", "Parlay", True),
        (SPREAD_PLUS_TOTAL_F5, "MLB", "Parlay", True),
        (ML_PLUS_ML, "MLB", "Parlay", False),
        (TOTAL_PLUS_TOTAL, "MLB", "Parlay", False),
        (THREE_LEGS, "MLB", "Parlay", False),
        (STRAIGHT_BET, "MLB", "Straight", False),
        (SPREAD_PLUS_TOTAL_FG, "NFL", "Parlay", False),
    ],
)
def test_is_mlb_correlated_parlay(desc, sport, bet_type, expected):
    assert is_mlb_correlated_parlay(desc, sport, bet_type) is expected


@pytest.mark.parametrize(
    "desc, expected",
    [
        (SPREAD_PLUS_TOTAL_FG, False),
        (SPREAD_PLUS_TOTAL_F5, True),
        ("PARLAY (2 TEAMS) | Yankees F5 -0.5 +110 | Over 4.5 -105", True),
        ("PARLAY (2 TEAMS) | Yankees First 5 -0.5 +110 | Over 4.5 -105", True),
        ("PARLAY (2 TEAMS) | Yankees -1.5 +140 | Over 8.5 -110", False),
    ],
)
def test_is_f5_parlay(desc, expected):
    assert is_f5_parlay(desc) is expected
