import json
from types import SimpleNamespace

import pytest
import requests

from nfl_specials import wz


def test_reads_retry_transient_timeouts(monkeypatch):
    monkeypatch.setattr(wz.time, "sleep", lambda _seconds: None)
    attempts = []

    def flaky_read():
        attempts.append(1)
        if len(attempts) < wz.READ_ATTEMPTS:
            raise requests.ReadTimeout("stalled")
        return "board"

    assert wz._retry_reads(flaky_read) == "board"
    assert len(attempts) == wz.READ_ATTEMPTS


def test_reads_give_up_after_the_last_attempt(monkeypatch):
    monkeypatch.setattr(wz.time, "sleep", lambda _seconds: None)

    def dead_read():
        raise requests.ConnectionError("down")

    with pytest.raises(requests.ConnectionError):
        wz._retry_reads(dead_read)


def test_other_errors_are_not_retried(monkeypatch):
    monkeypatch.setattr(wz.time, "sleep", lambda _seconds: None)
    attempts = []

    def broken_read():
        attempts.append(1)
        raise ValueError("bad JSON")

    with pytest.raises(ValueError):
        wz._retry_reads(broken_read)
    assert len(attempts) == 1


def test_every_specials_league_once():
    catalog = [{"IdLeague": 1135, "Description": "NFL WEEK 4 - SPECIALS"},
               {"IdLeague": 1135, "Description": "NFL WEEK 4 - SPECIALS"},   # listed under 2 menus
               {"IdLeague": 1140, "Description": "NFL WEEK 5 - SPECIALS"},
               {"IdLeague": 586, "Description": "NFL - TEAM TO SCORE 1ST"}]
    assert wz.specials_league_ids(catalog) == [1135, 1140]
    with pytest.raises(RuntimeError):
        wz.specials_league_ids([{"IdLeague": 586, "Description": "NFL - TEAM TO SCORE 1ST"}])


def test_a_fecta_is_placed_with_its_rotation_as_the_play(monkeypatch):
    # Wagerzon's web app bets a PROP game as "<vnum>_<idgm>_0_<odds>"; Play=5
    # passes the preview, then the submit fails "Couldn't find Game Line".
    special = wz._parse_special({"idgm": "6008644", "hnum": 777043, "vnum": 777043, "gmdt": "20261004",
                                 "htm": "RAIDERS SUPERFECTA (SCR 1ST, 1Q, 1H & GM)",
                                 "GameLines": [{"odds": "1730", "oddsh": "+1730"}]})
    sent = []
    monkeypatch.setattr(wz.single_placer, "place_single", lambda account, bet, use_free_play=False: sent.append(bet) or {})
    wz.place_fecta("Wagerzon", special, 250.0)
    assert sent[0]["play"] == 777043
    assert wz.single_placer.build_sel_for_single(sent[0]) == "777043_6008644_0_1730"
    detail = json.loads(wz.single_placer._build_confirm_payload(sent[0])["detailData"])[0]
    assert (detail["IdGame"], detail["Play"], detail["Amount"]) == (6008644, 777043, "250")


@pytest.mark.parametrize("use_free_play", [False, True])
def test_the_free_play_switch_reaches_the_placer_on_the_first_try_and_the_retry(monkeypatch, use_free_play):
    replies = [{"status": "auth_error", "error_msg": wz.single_placer.PREFLIGHT_LOGGED_OUT_MSG}, {"status": "placed"}]
    flags = []
    monkeypatch.setattr(wz.single_placer, "place_single",
                        lambda account, bet, use_free_play=False: flags.append(use_free_play) or replies.pop(0))
    monkeypatch.setattr(wz.wagerzon_auth, "clear_session_cache", lambda label=None: None)
    wz.place_fecta("Wagerzon", _special(), 250.0, use_free_play=use_free_play)
    assert flags == [use_free_play, use_free_play]


def test_each_accounts_two_pots_come_from_their_own_fields(monkeypatch):
    accounts = [SimpleNamespace(label="Wagerzon"), SimpleNamespace(label="WagerzonC")]
    snapshots = [SimpleNamespace(label="Wagerzon", available=500.0, free_play=60.0),
                 SimpleNamespace(label="WagerzonC", available=900.0, free_play=None)]
    monkeypatch.setattr(wz, "list_accounts", lambda: accounts)
    monkeypatch.setattr(wz.wagerzon_balance, "fetch_all", lambda given: snapshots)
    monkeypatch.setattr(wz.wagerzon_balance, "fetch_available_balance", lambda account: snapshots[accounts.index(account)])
    assert wz.read_balances() == {"Wagerzon": wz.WzBalances(500.0, 60.0), "WagerzonC": wz.WzBalances(900.0, None)}
    assert wz.read_account_balances("WagerzonC") == wz.WzBalances(900.0, None)
    assert wz.read_account_balances("nobody") == wz.WzBalances(None, None)


def test_the_expected_win_is_whole_dollars_like_wagerzon_quotes():
    # 2026-10-03: $204 at +855 is $1,744.20 exact; Wagerzon quoted $1,744 and
    # the cents figure was refused as a moved price.
    assert wz.wagerzon_to_win(204.0, 855) == 1744.0
    assert wz.wagerzon_to_win(250.0, 1730) == 4325.0
    assert wz.wagerzon_to_win(100.0, 1550) == 1550.0
    assert wz.wagerzon_to_win(21.0, 250) == 53.0  # 52.50 rounds up


def test_on_a_server_only_the_wagerzon_login_lines_are_loaded(tmp_path, monkeypatch):
    env_file = tmp_path / "bets.env"
    env_file.write_text("KALSHI_API_KEY_ID=kalshi-key\nBFA_PASSWORD=bfa-secret\n"
                        "WAGERZON_USERNAME=wz-user\nWAGERZON_PASSWORD=wz-pass\n"
                        "WAGERZONC_USERNAME=wzc-user\nWAGERZONC_PASSWORD=wzc-pass\nWAGERZOND_USERNAME=\n")
    monkeypatch.setattr(wz.config, "WAGERZON_ENV_PATH", str(env_file))
    for key in ("WAGERZON_USERNAME", "WAGERZON_PASSWORD", "WAGERZONC_USERNAME", "WAGERZONC_PASSWORD"):
        monkeypatch.setenv(key, "before")
    monkeypatch.delenv("KALSHI_API_KEY_ID", raising=False)
    monkeypatch.delenv("BFA_PASSWORD", raising=False)
    monkeypatch.delenv("WAGERZOND_USERNAME", raising=False)

    wz._load_wagerzon_env()

    keys = ("WAGERZON_USERNAME", "WAGERZON_PASSWORD", "WAGERZONC_USERNAME", "WAGERZONC_PASSWORD",
            "WAGERZOND_USERNAME", "KALSHI_API_KEY_ID", "BFA_PASSWORD")
    assert {key: wz.os.environ.get(key) for key in keys} == {
        "WAGERZON_USERNAME": "wz-user", "WAGERZON_PASSWORD": "wz-pass",
        "WAGERZONC_USERNAME": "wzc-user", "WAGERZONC_PASSWORD": "wzc-pass",
        "WAGERZOND_USERNAME": None, "KALSHI_API_KEY_ID": None, "BFA_PASSWORD": None}


def test_a_server_env_file_without_a_wagerzon_login_fails_loudly(tmp_path, monkeypatch):
    env_file = tmp_path / "bets.env"
    env_file.write_text("KALSHI_API_KEY_ID=kalshi-key\n")
    monkeypatch.setattr(wz.config, "WAGERZON_ENV_PATH", str(env_file))
    with pytest.raises(RuntimeError) as raised:
        wz._load_wagerzon_env()
    assert str(raised.value) == f"expected WAGERZON*_USERNAME / _PASSWORD lines in {env_file}, found none"


def _special():
    return wz._parse_special({"idgm": "6008644", "hnum": 777043, "vnum": 777043, "gmdt": "20261004",
                              "htm": "RAIDERS SUPERFECTA (SCR 1ST, 1Q, 1H & GM)",
                              "GameLines": [{"odds": "1730", "oddsh": "+1730"}]})


def test_a_login_expired_at_preflight_logs_in_again_and_retries_once(monkeypatch):
    # 2026-11-08: the cached login expired and every Place failed until restart.
    replies = [{"status": "auth_error", "error_msg": wz.single_placer.PREFLIGHT_LOGGED_OUT_MSG},
               {"status": "placed", "ticket_number": 123}]
    calls, cleared = [], []
    monkeypatch.setattr(wz.single_placer, "place_single", lambda account, bet, use_free_play=False: calls.append(account) or replies.pop(0))
    monkeypatch.setattr(wz.wagerzon_auth, "clear_session_cache", lambda label=None: cleared.append(label))
    assert wz.place_fecta("Wagerzon", _special(), 250.0)["status"] == "placed"
    assert calls == ["Wagerzon", "Wagerzon"]
    assert cleared == ["Wagerzon"]


def test_a_second_preflight_logout_is_not_retried_again(monkeypatch):
    reply = {"status": "auth_error", "error_msg": wz.single_placer.PREFLIGHT_LOGGED_OUT_MSG}
    calls = []
    monkeypatch.setattr(wz.single_placer, "place_single", lambda account, bet, use_free_play=False: calls.append(account) or dict(reply))
    monkeypatch.setattr(wz.wagerzon_auth, "clear_session_cache", lambda label=None: None)
    assert wz.place_fecta("Wagerzon", _special(), 250.0)["status"] == "auth_error"
    assert len(calls) == 2


def test_a_logout_at_submission_is_never_retried(monkeypatch):
    # The submit may have gone through; a retry could place the bet twice.
    reply = {"status": "auth_error", "error_msg": "Wagerzon returned HTML at submission (session expired)"}
    calls, cleared = [], []
    monkeypatch.setattr(wz.single_placer, "place_single", lambda account, bet, use_free_play=False: calls.append(account) or reply)
    monkeypatch.setattr(wz.wagerzon_auth, "clear_session_cache", lambda label=None: cleared.append(label))
    assert wz.place_fecta("Wagerzon", _special(), 250.0) is reply
    assert len(calls) == 1 and cleared == []
