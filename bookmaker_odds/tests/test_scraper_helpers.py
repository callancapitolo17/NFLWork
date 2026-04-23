"""Unit tests for non-network helpers in scraper.py.

These tests are pure-Python — no curl_cffi, no subprocess, no DuckDB. They
pin down the three behaviors the 'stop the popup' fix relies on:

  1. _session_looks_healthy()        — distinguishes "no games posted" from
                                       "session is broken"
  2. _can_launch_interactive_recon() — refuses to launch when stdin isn't a
                                       TTY, so piped subprocesses can't hang
  3. _recon_rate_limited()           — skips recon if another attempt
                                       happened < 1h ago, so a misbehaving
                                       caller can't storm the user with
                                       Chrome popups
"""
import io
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import scraper  # noqa: E402


# ---- _session_looks_healthy ------------------------------------------------

def test_session_healthy_when_not_blocked_and_login_ok():
    assert scraper._session_looks_healthy(blocked=False, login_ok=True) is True


def test_session_unhealthy_when_blocked():
    assert scraper._session_looks_healthy(blocked=True, login_ok=True) is False


def test_session_unhealthy_when_login_failed():
    assert scraper._session_looks_healthy(blocked=False, login_ok=False) is False


# ---- _can_launch_interactive_recon ----------------------------------------

def test_recon_blocked_when_stdin_not_tty(monkeypatch):
    fake_stdin = io.StringIO("")  # StringIO is not a tty
    monkeypatch.setattr(sys, "stdin", fake_stdin)
    assert scraper._can_launch_interactive_recon() is False


def test_recon_allowed_when_stdin_is_tty(monkeypatch):
    class FakeTTY:
        def isatty(self):
            return True
    monkeypatch.setattr(sys, "stdin", FakeTTY())
    assert scraper._can_launch_interactive_recon() is True


# ---- _recon_rate_limited ---------------------------------------------------

def test_rate_limit_false_when_sentinel_missing(tmp_path):
    sentinel = tmp_path / "never_existed"
    assert scraper._recon_rate_limited(sentinel, min_gap_sec=3600) is False


def test_rate_limit_true_when_sentinel_fresh(tmp_path):
    sentinel = tmp_path / "fresh"
    sentinel.touch()
    # Just now — less than an hour ago
    assert scraper._recon_rate_limited(sentinel, min_gap_sec=3600) is True


def test_rate_limit_false_when_sentinel_old(tmp_path):
    sentinel = tmp_path / "old"
    sentinel.touch()
    # Rewind mtime by 2 hours
    old = time.time() - 7200
    os.utime(sentinel, (old, old))
    assert scraper._recon_rate_limited(sentinel, min_gap_sec=3600) is False
