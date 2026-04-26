"""Unit tests for create_mlb_summary.py filter helpers.

Fixtures are real Wagerzon MLB parlay descriptions captured from the live
bet log on 2026-04-25. Wagerzon uses Unicode ½, "TOTAL o/u" prefix, and
"1H" team prefix for half-game bets — see comments on regex constants in
create_mlb_summary.py.
"""
import pytest
from create_mlb_summary import (
    is_mlb_correlated_parlay,
    is_f5_parlay,
    SPREAD_REGEX,
    TOTAL_REGEX,
    F5_REGEX,
)


# Real Wagerzon descriptions (as they appear in column D of Sheet1)
FG_TOTAL_FIRST = (
    "PARLAY (2 TEAMS) | TOTAL o8½-110 (NY YANKEES vrs SF GIANTS) "
    "( WILL WARREN - R / TYLER MAHLE - R ) | NY YANKEES -1½+130 "
    "( WILL WARREN - R / TYLER MAHLE - R )"
)
FG_INTEGER_TOTAL = (
    "PARLAY (2 TEAMS) | TOTAL o7-103 (CLE GUARDIANS vrs SEA MARINERS) "
    "( JOEY CANTILLO - L / BRYAN WOO - R ) | SEA MARINERS -1½+118 "
    "( JOEY CANTILLO - L / BRYAN WOO - R )"
)
FG_SPREAD_FIRST = (
    "PARLAY (2 TEAMS) | MIA MARLINS -1½+114 ( JOSE QUINTANA - L / MAX MEYER - R ) | "
    "TOTAL o8-106 (COL ROCKIES vrs MIA MARLINS) ( JOSE QUINTANA - L / MAX MEYER - R )"
)
F5_TOTAL_FIRST = (
    "PARLAY (2 TEAMS) | TOTAL u4+115 (1H TB RAYS vrs 1H MIN TWINS) "
    "( STEVEN MATZ - L / MICK ABEL - R ) | 1H TB RAYS +½-145 "
    "( STEVEN MATZ - L / MICK ABEL - R )"
)
F5_SPREAD_FIRST = (
    "PARLAY (2 TEAMS) | 1H TEX RANGERS -½+108 ( RHETT LOWDER - R / KUMAR ROCKER - R ) | "
    "TOTAL o4½-110 (1H CIN REDS vrs 1H TEX RANGERS) ( RHETT LOWDER - R / KUMAR ROCKER - R )"
)
ML_PLUS_ML = (
    "PARLAY (2 TEAMS) | NY YANKEES ML +100 ( WARREN - R ) | "
    "CHI CUBS ML -120 ( SMYLY - L )"
)
ML_PLUS_TOTAL_ONLY = (
    "PARLAY (2 TEAMS) | NY YANKEES ML +100 ( WARREN - R ) | "
    "TOTAL o8½-110 (NY YANKEES vrs SF GIANTS) ( WARREN - R / MAHLE - R )"
)
THREE_LEGS = (
    "PARLAY (3 TEAMS) | NY YANKEES -1½+130 (...) | TOTAL o8½-110 (...) | "
    "CHI CUBS ML -200 (...)"
)
STRAIGHT_BET = "STRAIGHT BET | NY YANKEES -1½+130 ( WARREN - R / MAHLE - R )"


@pytest.mark.parametrize(
    "desc, sport, bet_type, expected",
    [
        # Real correlated parlays (FG variants and F5 variants)
        (FG_TOTAL_FIRST, "MLB", "Parlay", True),
        (FG_INTEGER_TOTAL, "MLB", "Parlay", True),
        (FG_SPREAD_FIRST, "MLB", "Parlay", True),
        (F5_TOTAL_FIRST, "MLB", "Parlay", True),
        (F5_SPREAD_FIRST, "MLB", "Parlay", True),
        # Things that should NOT match
        (ML_PLUS_ML, "MLB", "Parlay", False),       # no spread, no total
        (ML_PLUS_TOTAL_ONLY, "MLB", "Parlay", False),  # has total, no spread
        (THREE_LEGS, "MLB", "Parlay", False),       # PARLAY (3 TEAMS) excluded by header check
        (STRAIGHT_BET, "MLB", "Straight", False),   # not a parlay
        (FG_TOTAL_FIRST, "NFL", "Parlay", False),   # wrong sport
    ],
)
def test_is_mlb_correlated_parlay(desc, sport, bet_type, expected):
    assert is_mlb_correlated_parlay(desc, sport, bet_type) is expected


@pytest.mark.parametrize(
    "desc, expected",
    [
        # FG bets — no 1H/F5 marker
        (FG_TOTAL_FIRST, False),
        (FG_INTEGER_TOTAL, False),
        (FG_SPREAD_FIRST, False),
        # F5 bets — both directions
        (F5_TOTAL_FIRST, True),
        (F5_SPREAD_FIRST, True),
    ],
)
def test_is_f5_parlay(desc, expected):
    assert is_f5_parlay(desc) is expected
