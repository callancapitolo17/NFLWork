"""
Unit tests for dk_leg_resolvers.

Uses a tiny inline fixture (not a captured recon JSON) to keep tests fast
and avoid committing large blobs. Fixture mimics the relevant DK SGP event
payload structure with one selection per market.
"""
import sys
from pathlib import Path

# Add mlb_sgp to path so we can import resolvers
sys.path.insert(0, str(Path(__file__).parent.parent))

from dk_leg_resolvers import (  # noqa: E402
    resolve_scores_first,
    resolve_wins_period,
    resolve_team_total_under,
    resolve_team_total_over,
    resolve_legs,
    find_market_by_name,
)

# ---------------------------------------------------------------------------
# Inline fixture — minimal DK event payload covering the 4 supported markets
# ---------------------------------------------------------------------------

FIXTURE = {
    'data': {
        'markets': [
            {
                'marketName': '1st Run',
                'tags': ['SGP', 'YourBetEligible'],
                'selections': [
                    {'name': 'CHI Cubs',  'id': '0QA-CUBS-1ST-RUN'},
                    {'name': 'SD Padres', 'id': '0QA-SD-1ST-RUN'},
                ],
            },
            {
                'marketName': 'Moneyline',
                'tags': ['SGP'],
                'selections': [
                    {'name': 'CHI Cubs',  'id': '0ML-CUBS_3'},
                    {'name': 'SD Padres', 'id': '0ML-SD_1'},
                ],
            },
            {
                'marketName': '1st 5 Innings',
                'tags': ['SGP'],
                'selections': [
                    {'name': 'CHI Cubs',  'id': '0QA-CUBS-F5'},
                    {'name': 'SD Padres', 'id': '0QA-SD-F5'},
                ],
            },
            {
                'marketName': 'CHI Cubs: Team Total Runs',
                'tags': ['SGP'],
                'selections': [
                    {'label': 'Over',  'point': 3.5, 'id': '0OU-CUBS-O35'},
                    {'label': 'Under', 'point': 3.5, 'id': '0OU-CUBS-U35'},
                ],
            },
            {
                'marketName': 'SD Padres: Team Total Runs',
                'tags': ['SGP'],
                'selections': [
                    {'label': 'Over',  'point': 4.5, 'id': '0OU-SD-O45'},
                    {'label': 'Under', 'point': 4.5, 'id': '0OU-SD-U45'},
                ],
            },
            # Non-SGP-eligible market — should be ignored by find_market_by_name
            {
                'marketName': 'Inning of First/Last Score',
                'tags': ['Cashout'],   # NO 'SGP' tag
                'selections': [
                    {'name': 'CHI Cubs',  'id': '0QA-IGNORE-1'},
                    {'name': 'SD Padres', 'id': '0QA-IGNORE-2'},
                ],
            },
        ],
    },
}

TEAM_NAMES = {'home': 'CHI Cubs', 'away': 'SD Padres'}


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_find_market_by_name_finds_sgp_eligible_market():
    m = find_market_by_name(FIXTURE, '1st Run')
    assert m is not None
    assert m['marketName'] == '1st Run'

def test_find_market_by_name_skips_non_sgp_market():
    # 'Inning of First/Last Score' is in fixture but lacks 'SGP' tag
    m = find_market_by_name(FIXTURE, 'Inning of First/Last Score')
    assert m is None

def test_find_market_by_name_returns_none_for_missing_market():
    assert find_market_by_name(FIXTURE, 'No Such Market') is None


def test_scores_first_home():
    sid = resolve_scores_first({'type': 'scores_first'}, 'home', FIXTURE, TEAM_NAMES)
    assert sid == '0QA-CUBS-1ST-RUN'

def test_scores_first_away():
    sid = resolve_scores_first({'type': 'scores_first'}, 'away', FIXTURE, TEAM_NAMES)
    assert sid == '0QA-SD-1ST-RUN'


def test_wins_period_FG_home():
    sid = resolve_wins_period({'type': 'wins_period', 'period': 'FG'}, 'home', FIXTURE, TEAM_NAMES)
    assert sid == '0ML-CUBS_3'

def test_wins_period_F5_away():
    sid = resolve_wins_period({'type': 'wins_period', 'period': 'F5'}, 'away', FIXTURE, TEAM_NAMES)
    assert sid == '0QA-SD-F5'

def test_wins_period_F3_returns_none():
    # F3 not in DK's 2-way ML primitives — resolver returns None
    sid = resolve_wins_period({'type': 'wins_period', 'period': 'F3'}, 'home', FIXTURE, TEAM_NAMES)
    assert sid is None


def test_team_total_under_home_correct_line():
    sid = resolve_team_total_under(
        {'type': 'team_total_under', 'line': 3.5}, 'home', FIXTURE, TEAM_NAMES,
    )
    assert sid == '0OU-CUBS-U35'

def test_team_total_over_away_correct_line():
    sid = resolve_team_total_over(
        {'type': 'team_total_over', 'line': 4.5}, 'away', FIXTURE, TEAM_NAMES,
    )
    assert sid == '0OU-SD-O45'

def test_team_total_wrong_line_returns_none():
    sid = resolve_team_total_under(
        {'type': 'team_total_under', 'line': 99.5}, 'home', FIXTURE, TEAM_NAMES,
    )
    assert sid is None


def test_resolve_legs_full_trifecta_home():
    legs = [
        {'type': 'scores_first'},
        {'type': 'wins_period', 'period': 'F5'},
        {'type': 'wins_period', 'period': 'FG'},
    ]
    sids = resolve_legs(legs, 'home', FIXTURE, TEAM_NAMES)
    assert sids == ['0QA-CUBS-1ST-RUN', '0QA-CUBS-F5', '0ML-CUBS_3']

def test_resolve_legs_full_grand_slam_away():
    legs = [
        {'type': 'scores_first'},
        {'type': 'wins_period', 'period': 'F5'},
        {'type': 'wins_period', 'period': 'FG'},
        {'type': 'team_total_under', 'line': 4.5},
    ]
    sids = resolve_legs(legs, 'away', FIXTURE, TEAM_NAMES)
    assert sids == ['0QA-SD-1ST-RUN', '0QA-SD-F5', '0ML-SD_1', '0OU-SD-U45']

def test_resolve_legs_missing_leg_returns_none():
    # F3 not supported → whole resolution fails → None
    legs = [
        {'type': 'scores_first'},
        {'type': 'wins_period', 'period': 'F3'},
        {'type': 'wins_period', 'period': 'FG'},
    ]
    sids = resolve_legs(legs, 'home', FIXTURE, TEAM_NAMES)
    assert sids is None

def test_resolve_legs_unknown_leg_type_returns_none():
    legs = [{'type': 'scores_first'}, {'type': 'mystery_leg'}]
    sids = resolve_legs(legs, 'home', FIXTURE, TEAM_NAMES)
    assert sids is None
