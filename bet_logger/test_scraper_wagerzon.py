"""Unit tests for scraper_wagerzon.py — account config and bet parsing.

These tests run without hitting Wagerzon's live API. Network-dependent
behaviors (login, fetch_history_json) are validated end-to-end via
--dry-run runs in Task 9.
"""
import copy
import pytest
from scraper_wagerzon import ACCOUNTS, parse_api_bets


def test_accounts_has_default_and_j():
    assert set(ACCOUNTS.keys()) == {'default', 'j'}


def test_accounts_default_fields():
    acct = ACCOUNTS['default']
    assert acct['username_env'] == 'WAGERZON_USERNAME'
    assert acct['password_env'] == 'WAGERZON_PASSWORD'
    assert acct['platform'] == 'Wagerzon'
    assert acct['bet_multiplier'] == 1.0
    assert acct['shared_sheet'] is None


def test_accounts_j_fields():
    acct = ACCOUNTS['j']
    assert acct['username_env'] == 'WAGERZONJ_USERNAME'
    assert acct['password_env'] == 'WAGERZONJ_PASSWORD'
    assert acct['platform'] == 'WagerzonJ'
    assert acct['bet_multiplier'] == 0.875
    assert acct['shared_sheet'] == 'Shared'


# Minimal fixture matching the HistoryHelper.aspx JSON shape.
# - Single straight bet, MLB, $100 risk, lost, with American odds -110.
# - Sufficient to exercise platform/multiplier/raw-risk paths.
SAMPLE_HISTORY = {
    'details': [
        {
            'wager': [
                {
                    'WagerOrTrans': 'WAGER',
                    'PlacedDate': '04/29/2026',
                    'RiskAmount': '100.00',
                    'WinAmount': '90.91',
                    'WinLoss': '-100.00',
                    'Result': 'LOSE',
                    'HeaderDesc': 'STRAIGHT BET',
                    'details': [
                        {
                            'IdSport': 'MLB',
                            'DetailDesc': '[967] NY YANKEES -1.5-110 (NY YANKEES vrs SF GIANTS)',
                        }
                    ],
                }
            ]
        }
    ]
}


def test_parse_api_bets_default_platform_and_no_multiplier():
    bets = parse_api_bets(SAMPLE_HISTORY, platform='Wagerzon', bet_multiplier=1.0)
    assert len(bets) == 1
    bet = bets[0]
    assert bet['platform'] == 'Wagerzon'
    assert bet['bet_amount'] == 100.0
    assert bet['_raw_risk'] == 100.0


def test_parse_api_bets_j_applies_multiplier_to_bet_amount_only():
    bets = parse_api_bets(SAMPLE_HISTORY, platform='WagerzonJ', bet_multiplier=0.875)
    assert len(bets) == 1
    bet = bets[0]
    assert bet['platform'] == 'WagerzonJ'
    assert bet['bet_amount'] == 87.5
    assert bet['_raw_risk'] == 100.0


def test_parse_api_bets_multiplier_does_not_change_odds_or_result():
    """Odds, decimal odds, and result must be invariant to multiplier."""
    default = parse_api_bets(SAMPLE_HISTORY, platform='Wagerzon', bet_multiplier=1.0)[0]
    j = parse_api_bets(SAMPLE_HISTORY, platform='WagerzonJ', bet_multiplier=0.875)[0]
    assert default['odds'] == j['odds']
    assert default['dec'] == j['dec']
    assert default['result'] == j['result']


def test_parse_api_bets_rounds_to_two_decimals():
    """0.875 × 33.33 = 29.16375 should round to 29.16, not write a long float."""
    history = copy.deepcopy(SAMPLE_HISTORY)
    history['details'][0]['wager'][0]['RiskAmount'] = '33.33'
    bets = parse_api_bets(history, platform='WagerzonJ', bet_multiplier=0.875)
    assert bets[0]['bet_amount'] == 29.16
    assert bets[0]['_raw_risk'] == 33.33
