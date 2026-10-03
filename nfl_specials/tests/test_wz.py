import json

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
    monkeypatch.setattr(wz.single_placer, "place_single", lambda account, bet: sent.append(bet) or {})
    wz.place_fecta("Wagerzon", special, 250.0)
    assert sent[0]["play"] == 777043
    assert wz.single_placer.build_sel_for_single(sent[0]) == "777043_6008644_0_1730"
    detail = json.loads(wz.single_placer._build_confirm_payload(sent[0])["detailData"])[0]
    assert (detail["IdGame"], detail["Play"], detail["Amount"]) == (6008644, 777043, "250")


def test_the_expected_win_is_whole_dollars_like_wagerzon_quotes():
    # 2026-10-03: $204 at +855 is $1,744.20 exact; Wagerzon quoted $1,744 and
    # the cents figure was refused as a moved price.
    assert wz.wagerzon_to_win(204.0, 855) == 1744.0
    assert wz.wagerzon_to_win(250.0, 1730) == 4325.0
    assert wz.wagerzon_to_win(100.0, 1550) == 1550.0
    assert wz.wagerzon_to_win(21.0, 250) == 53.0  # 52.50 rounds up
