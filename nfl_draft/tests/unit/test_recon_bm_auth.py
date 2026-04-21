"""Tests for the auto-re-login path in nfl_draft/scrapers/recon_bm.py.

Two layers:

1. `_looks_like_expired_cookies` is a pure heuristic — fuzz it with a real
   fixture, a guest-style response, an empty-schedule response, and the
   transport-layer `None`. It should fire when BM isn't showing us data but
   stay quiet for a real markets payload.

2. `run_rest_phase` should trigger ONE retry via `_auto_relogin` when the
   first probe pass looks expired, and exactly zero retries when the first
   pass succeeds. We monkeypatch the network-touching helpers so nothing
   leaves the process.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from nfl_draft.scrapers import recon_bm


FIXTURE_PATH = (
    Path(__file__).resolve().parent.parent
    / "fixtures" / "bookmaker" / "draft_markets.json"
)


# ---------------------------------------------------------------------------
# Heuristic tests
# ---------------------------------------------------------------------------


def _real_markets_payload() -> dict:
    """Load the raw BM `data` blob from the committed fixture.

    `save_fixture` wraps the book response in a {book, meta, data} envelope,
    so we pop `.data` to get the exact shape `_fetch_schedule` returns.
    """
    envelope = json.loads(FIXTURE_PATH.read_text())
    assert "data" in envelope, "fixture envelope missing `data` key"
    return envelope["data"]


def test_heuristic_passes_real_markets_response():
    """A real Schedule with non-empty Leagues.League must NOT look expired.

    If this flipped to True, we'd start auto-re-logging in on successful
    probes and eventually burn our account with rate limits.
    """
    real = _real_markets_payload()
    assert recon_bm._looks_like_expired_cookies([real]) is False


def test_heuristic_flags_empty_schedule():
    """GetSchedule with `League: []` is the classic expired-session shape."""
    empty = {"Schedule": {"Data": {"Leagues": {"League": []}}}}
    assert recon_bm._looks_like_expired_cookies([empty]) is True


def test_heuristic_flags_guest_dummy_player_marker():
    """`IsDummyPlayer: true` anywhere in the response = guest session."""
    guest = {
        "Schedule": {"Data": {"Leagues": {"League": []}}},
        "PlayerStatus": {"IsDummyPlayer": True},
    }
    assert recon_bm._looks_like_expired_cookies([guest]) is True


def test_heuristic_flags_transport_failures():
    """All probes returning None (HTTP 401/403/exception) = expired."""
    assert recon_bm._looks_like_expired_cookies([None, None, None]) is True


def test_heuristic_quiet_when_one_probe_has_data():
    """A single real hit across many empty probes = NOT expired.

    We probe ~10 league IDs; only one is the draft league. Every other ID
    legitimately returns an empty Schedule. The heuristic must not flag this.
    """
    real = _real_markets_payload()
    empty = {"Schedule": {"Data": {"Leagues": {"League": []}}}}
    results = [empty, empty, real, empty, None]
    assert recon_bm._looks_like_expired_cookies(results) is False


def test_heuristic_empty_input_is_not_expired():
    """No probes attempted means we can't judge — default to not-expired."""
    assert recon_bm._looks_like_expired_cookies([]) is False


# ---------------------------------------------------------------------------
# run_rest_phase retry-path tests
# ---------------------------------------------------------------------------


class _FakeResponse:
    """Minimal stand-in for the session.get response used during CF warmup."""
    def __init__(self, status: int = 200):
        self.status_code = status


class _FakeSession:
    """Replaces curl_cffi.Session; we never actually talk to BM in tests."""
    def __init__(self):
        self.cookies = {}

    def get(self, *args, **kwargs):
        return _FakeResponse(200)


def _patch_session_creation(monkeypatch):
    """Route curl_cffi.Session through a local fake so no network happens."""
    import nfl_draft.scrapers.recon_bm as mod
    monkeypatch.setattr(mod.cffi_requests, "Session", lambda **kwargs: _FakeSession())
    monkeypatch.setattr(mod, "_load_cookies", lambda session: True)
    # Keep the cookie save side-effect free even if we hit the success path.
    monkeypatch.setattr(mod, "_save_cookies", lambda session: True)


def test_run_rest_phase_returns_fast_when_probe_succeeds(monkeypatch):
    """If the first probe already returns draft hits, do NOT call _auto_relogin.

    Guardrail check: the retry path must be strictly opt-in based on the
    heuristic. Burning a login attempt on every call would flag the account.
    """
    real = _real_markets_payload()
    _patch_session_creation(monkeypatch)

    draft_entry = {"Description": "NFL DRAFT 2026 - ODDS TO WIN"}
    monkeypatch.setattr(
        recon_bm, "_probe_candidates",
        lambda session, candidates: (
            [(117, real, draft_entry)],
            [real],
        ),
    )

    calls = []
    monkeypatch.setattr(
        recon_bm, "_auto_relogin",
        lambda session: calls.append("login") or True,
    )

    data, url, meta = recon_bm.run_rest_phase()

    assert data is real
    assert meta["method"] == "rest_probe"
    assert calls == [], "auto-re-login should not fire when the first probe succeeds"


def test_run_rest_phase_retries_once_after_successful_relogin(monkeypatch):
    """Cold probe returns empty -> _auto_relogin succeeds -> second probe returns data."""
    real = _real_markets_payload()
    _patch_session_creation(monkeypatch)

    call_count = {"probe": 0, "login": 0, "save": 0}
    draft_entry = {"Description": "NFL DRAFT 2026 - ODDS TO WIN"}

    def fake_probe(session, candidates):
        call_count["probe"] += 1
        if call_count["probe"] == 1:
            # First pass: all empty.
            empty = {"Schedule": {"Data": {"Leagues": {"League": []}}}}
            return [], [empty for _ in candidates]
        # Second pass: real hit.
        return [(117, real, draft_entry)], [real]

    def fake_relogin(session):
        call_count["login"] += 1
        return True

    def fake_save(session):
        call_count["save"] += 1
        return True

    monkeypatch.setattr(recon_bm, "_probe_candidates", fake_probe)
    monkeypatch.setattr(recon_bm, "_auto_relogin", fake_relogin)
    monkeypatch.setattr(recon_bm, "_save_cookies", fake_save)

    data, url, meta = recon_bm.run_rest_phase()

    assert data is real
    assert call_count["probe"] == 2, "should retry probe exactly once after re-login"
    assert call_count["login"] == 1, "should attempt login exactly once"
    assert call_count["save"] == 1, "should persist fresh cookies after login"


def test_run_rest_phase_no_retry_when_relogin_fails(monkeypatch):
    """If _auto_relogin returns False, don't retry — fail fast and bail."""
    _patch_session_creation(monkeypatch)

    call_count = {"probe": 0, "login": 0}

    def fake_probe(session, candidates):
        call_count["probe"] += 1
        empty = {"Schedule": {"Data": {"Leagues": {"League": []}}}}
        return [], [empty for _ in candidates]

    def fake_relogin(session):
        call_count["login"] += 1
        return False

    monkeypatch.setattr(recon_bm, "_probe_candidates", fake_probe)
    monkeypatch.setattr(recon_bm, "_auto_relogin", fake_relogin)

    data, url, meta = recon_bm.run_rest_phase()

    assert data is None
    assert call_count["probe"] == 1, "no second probe when login fails"
    assert call_count["login"] == 1, "one login attempt max"


def test_run_rest_phase_no_login_when_probes_dont_look_expired(monkeypatch):
    """If probes came back with SOME real-looking data but no NFL-Draft match,
    we should not attempt a login. This is the 'BM isn't posting draft markets
    today' case — burning a login would be wasteful and account-risky."""
    _patch_session_creation(monkeypatch)

    call_count = {"probe": 0, "login": 0}

    # A league response that has leagues but no NFL draft match — looks healthy.
    healthy_non_draft = {
        "Schedule": {"Data": {"Leagues": {"League": [
            {"Description": "SOME OTHER PROP MARKET", "dateGroup": []}
        ]}}}
    }

    def fake_probe(session, candidates):
        call_count["probe"] += 1
        return [], [healthy_non_draft for _ in candidates]

    def fake_relogin(session):
        call_count["login"] += 1
        return True

    monkeypatch.setattr(recon_bm, "_probe_candidates", fake_probe)
    monkeypatch.setattr(recon_bm, "_auto_relogin", fake_relogin)

    data, url, meta = recon_bm.run_rest_phase()

    assert data is None
    assert call_count["probe"] == 1
    assert call_count["login"] == 0, (
        "healthy non-draft probes must NOT trigger auto-re-login"
    )


def test_auto_relogin_missing_creds_returns_false(monkeypatch):
    """Without creds in the env, _auto_relogin should fail cleanly, not crash."""
    monkeypatch.delenv("BOOKMAKER_USERNAME", raising=False)
    monkeypatch.delenv("BOOKMAKER_PASSWORD", raising=False)
    # Neutralize load_env so test doesn't repopulate from the real .env file.
    monkeypatch.setattr(recon_bm, "load_env", lambda: None)

    # Also neutralize _login so even if creds somehow got through it can't fire.
    monkeypatch.setattr(recon_bm, "_login", lambda s, u, p: pytest.fail("should not be called"))

    assert recon_bm._auto_relogin(_FakeSession()) is False
