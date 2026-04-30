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
