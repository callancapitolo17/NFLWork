"""Unit tests for wagerzon_accounts registry."""
from __future__ import annotations

import pytest

from wagerzon_accounts import (
    WagerzonAccount,
    list_accounts,
    get_account,
    AccountNotFoundError,
)


@pytest.fixture
def env_two_accounts(monkeypatch):
    monkeypatch.setenv("WAGERZON_USERNAME", "primary_user")
    monkeypatch.setenv("WAGERZON_PASSWORD", "primary_pw")
    monkeypatch.setenv("WAGERZONJ_USERNAME", "j_user")
    monkeypatch.setenv("WAGERZONJ_PASSWORD", "j_pw")
    yield


@pytest.fixture
def env_three_accounts(env_two_accounts, monkeypatch):
    monkeypatch.setenv("WAGERZONC_USERNAME", "c_user")
    monkeypatch.setenv("WAGERZONC_PASSWORD", "c_pw")
    yield


@pytest.fixture
def env_clean(monkeypatch):
    """Strip every WAGERZON* env var so tests start from a clean slate."""
    for k in list(__import__("os").environ.keys()):
        if k.startswith("WAGERZON"):
            monkeypatch.delenv(k, raising=False)
    yield


def test_primary_only(env_clean, monkeypatch):
    monkeypatch.setenv("WAGERZON_USERNAME", "u")
    monkeypatch.setenv("WAGERZON_PASSWORD", "p")
    accounts = list_accounts()
    assert len(accounts) == 1
    assert accounts[0] == WagerzonAccount(label="Wagerzon", suffix="", username="u", password="p")


def test_two_accounts_primary_first(env_clean, env_two_accounts):
    accounts = list_accounts()
    labels = [a.label for a in accounts]
    assert labels == ["Wagerzon", "WagerzonJ"]


def test_three_accounts_alphabetical_after_primary(env_clean, env_three_accounts):
    accounts = list_accounts()
    labels = [a.label for a in accounts]
    assert labels == ["Wagerzon", "WagerzonC", "WagerzonJ"]


def test_skip_pair_missing_password(env_clean, monkeypatch, caplog):
    monkeypatch.setenv("WAGERZON_USERNAME", "u")
    monkeypatch.setenv("WAGERZON_PASSWORD", "p")
    monkeypatch.setenv("WAGERZONJ_USERNAME", "j_user")
    # No WAGERZONJ_PASSWORD on purpose
    accounts = list_accounts()
    assert [a.label for a in accounts] == ["Wagerzon"]
    assert "WAGERZONJ" in caplog.text


def test_skip_lowercase_suffix(env_clean, monkeypatch):
    monkeypatch.setenv("WAGERZON_USERNAME", "u")
    monkeypatch.setenv("WAGERZON_PASSWORD", "p")
    monkeypatch.setenv("WAGERZONj_USERNAME", "junk")
    monkeypatch.setenv("WAGERZONj_PASSWORD", "junk")
    accounts = list_accounts()
    assert [a.label for a in accounts] == ["Wagerzon"]


def test_skip_multichar_suffix(env_clean, monkeypatch):
    monkeypatch.setenv("WAGERZON_USERNAME", "u")
    monkeypatch.setenv("WAGERZON_PASSWORD", "p")
    monkeypatch.setenv("WAGERZONJX_USERNAME", "junk")
    monkeypatch.setenv("WAGERZONJX_PASSWORD", "junk")
    accounts = list_accounts()
    assert [a.label for a in accounts] == ["Wagerzon"]


def test_get_account_by_label(env_clean, env_two_accounts):
    acct = get_account("WagerzonJ")
    assert acct.label == "WagerzonJ"
    assert acct.username == "j_user"


def test_get_account_unknown_raises(env_clean, env_two_accounts):
    with pytest.raises(AccountNotFoundError):
        get_account("WagerzonZ")


def test_no_accounts_returns_empty_list(env_clean):
    assert list_accounts() == []
