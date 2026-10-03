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
