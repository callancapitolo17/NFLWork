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
                'marketName': 'CHI Cubs Run Scored - 1st Inning?',
                'tags': ['SGP', 'YourBetEligible'],
                'selections': [
                    {'name': 'Yes', 'id': '0QA-CUBS-RS1ST-YES'},
                    {'name': 'No',  'id': '0QA-CUBS-RS1ST-NO'},
                ],
            },
            {
                'marketName': 'SD Padres Run Scored - 1st Inning?',
                'tags': ['SGP', 'YourBetEligible'],
                'selections': [
                    {'name': 'Yes', 'id': '0QA-SD-RS1ST-YES'},
                    {'name': 'No',  'id': '0QA-SD-RS1ST-NO'},
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
                    {'label': 'Over',  'id': '0OU100001O350_1'},
                    {'label': 'Under', 'id': '0OU100001U350_3'},
                ],
            },
            {
                'marketName': 'SD Padres: Team Total Runs',
                'tags': ['SGP'],
                'selections': [
                    {'label': 'Over',  'id': '0OU100002O450_1'},
                    {'label': 'Under', 'id': '0OU100002U450_3'},
                ],
            },
            {
                'marketName': 'Run Line - 1st 5 Innings',
                'tags': ['SGP'],
                'selections': [
                    # -0.5 line (strict win) — what resolver should pick
                    {'name': 'CHI Cubs',  'id': '0HC-RL-CUBS-N50_1'},
                    {'name': 'SD Padres', 'id': '0HC-RL-SD-N50_3'},
                    # +0.5 line (win-or-tie) — resolver should SKIP
                    {'name': 'CHI Cubs',  'id': '0HC-RL-CUBS-P50_1'},
                    {'name': 'SD Padres', 'id': '0HC-RL-SD-P50_3'},
                    # -1.5 line — resolver should also skip
                    {'name': 'CHI Cubs',  'id': '0HC-RL-CUBS-N150_1'},
                    {'name': 'SD Padres', 'id': '0HC-RL-SD-N150_3'},
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
    m = find_market_by_name(FIXTURE, 'CHI Cubs Run Scored - 1st Inning?')
    assert m is not None
    assert m['marketName'] == 'CHI Cubs Run Scored - 1st Inning?'

def test_find_market_by_name_skips_non_sgp_market():
    # 'Inning of First/Last Score' is in fixture but lacks 'SGP' tag
    m = find_market_by_name(FIXTURE, 'Inning of First/Last Score')
    assert m is None

def test_find_market_by_name_returns_none_for_missing_market():
    assert find_market_by_name(FIXTURE, 'No Such Market') is None


def test_scores_first_home():
    sid = resolve_scores_first({'type': 'scores_first'}, 'home', FIXTURE, TEAM_NAMES)
    assert sid == '0QA-CUBS-RS1ST-YES'

def test_scores_first_away():
    sid = resolve_scores_first({'type': 'scores_first'}, 'away', FIXTURE, TEAM_NAMES)
    assert sid == '0QA-SD-RS1ST-YES'

def test_scores_first_picks_yes_not_no():
    """Defensive: confirm the resolver picks the Yes selection, not the No."""
    sid = resolve_scores_first({'type': 'scores_first'}, 'home', FIXTURE, TEAM_NAMES)
    assert sid == '0QA-CUBS-RS1ST-YES'
    assert 'NO' not in sid

def test_scores_first_returns_none_when_market_missing():
    """If DK doesn't post the inning-1 Y/N market, resolver returns None
    (graceful degrade — R-side blend falls back to model-only)."""
    fixture_no_market = {'data': {'markets': []}}
    sid = resolve_scores_first({'type': 'scores_first'}, 'home', fixture_no_market, TEAM_NAMES)
    assert sid is None


def test_wins_period_FG_home():
    sid = resolve_wins_period({'type': 'wins_period', 'period': 'FG'}, 'home', FIXTURE, TEAM_NAMES)
    assert sid == '0ML-CUBS_3'

def test_wins_period_F5_away():
    # F5 resolver now uses Run Line -0.5 (strict win) instead of 1st 5 Innings ML
    sid = resolve_wins_period({'type': 'wins_period', 'period': 'F5'}, 'away', FIXTURE, TEAM_NAMES)
    assert sid == '0HC-RL-SD-N50_3'

def test_wins_period_F5_home_uses_strict_win_run_line():
    """F5 resolver picks the -0.5 line (N50 in ID), not the +0.5 (P50) which
    would include ties. Strict-win semantic matches our model's home_margin_f5>0."""
    sid = resolve_wins_period({'type': 'wins_period', 'period': 'F5'}, 'home', FIXTURE, TEAM_NAMES)
    assert sid == '0HC-RL-CUBS-N50_1'

def test_wins_period_F5_skips_P50_winortie_line():
    """If only +0.5 (P50) lines are present, resolver should return None
    rather than pick an incorrect-semantic line."""
    fixture_only_p50 = {
        'data': {
            'markets': [
                {
                    'marketName': 'Run Line - 1st 5 Innings',
                    'tags': ['SGP'],
                    'selections': [
                        # Only +0.5 (P50) lines — strict-win not available
                        {'name': 'CHI Cubs',  'id': '0HC-RL-CUBS-P50_1'},
                        {'name': 'SD Padres', 'id': '0HC-RL-SD-P50_3'},
                    ],
                },
            ],
        },
    }
    sid = resolve_wins_period({'type': 'wins_period', 'period': 'F5'},
                              'home', fixture_only_p50, TEAM_NAMES)
    assert sid is None

def test_wins_period_F3_returns_none():
    # F3 not in DK's 2-way ML primitives — resolver returns None
    sid = resolve_wins_period({'type': 'wins_period', 'period': 'F3'}, 'home', FIXTURE, TEAM_NAMES)
    assert sid is None


def test_team_total_under_home_correct_line():
    sid = resolve_team_total_under(
        {'type': 'team_total_under', 'line': 3.5}, 'home', FIXTURE, TEAM_NAMES,
    )
    assert sid == '0OU100001U350_3'

def test_team_total_over_away_correct_line():
    sid = resolve_team_total_over(
        {'type': 'team_total_over', 'line': 4.5}, 'away', FIXTURE, TEAM_NAMES,
    )
    assert sid == '0OU100002O450_1'

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
    assert sids == ['0QA-CUBS-RS1ST-YES', '0HC-RL-CUBS-N50_1', '0ML-CUBS_3']

def test_resolve_legs_full_grand_slam_away():
    legs = [
        {'type': 'scores_first'},
        {'type': 'wins_period', 'period': 'F5'},
        {'type': 'wins_period', 'period': 'FG'},
        {'type': 'team_total_under', 'line': 4.5},
    ]
    sids = resolve_legs(legs, 'away', FIXTURE, TEAM_NAMES)
    assert sids == ['0QA-SD-RS1ST-YES', '0HC-RL-SD-N50_3', '0ML-SD_1', '0OU100002U450_3']

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
