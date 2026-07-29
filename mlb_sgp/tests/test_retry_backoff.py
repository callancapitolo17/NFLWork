"""Regression tests for issue #34 — bounded retry/backoff for the 6 clients.

Issue #33 made transport failures LOUD, and in doing so measurably widened
their blast radius: one 5xx on a single event now aborts a whole book's
cycle where it used to skip just that game. Retry is what closes that
window — a transient blip must not cost a book its fetch (nor, post-#54,
its VOTE on a live quote).

Contract pinned here:

  connection error / timeout / 5xx / 429   -> retry, bounded by the profile
  401 / 403 / 404 / other 4xx              -> NO retry, immediate verdict
  200                                      -> returned untouched, zero sleeps

Two profiles, because the on-demand pricing path runs inside an RFQ quote
budget while the sweep does not. Every sleep in this file is recorded, not
slept (see conftest.py::recorded_sleeps).
"""
from __future__ import annotations

import pytest

from mlb_sgp._shared import (RETRY_BACKGROUND, RETRY_LIVE, BookTransportError,
                             RetryProfile, check_response, is_retryable_exception,
                             is_retryable_status, json_or_raise,
                             request_with_retry)


class FakeResponse:
    def __init__(self, status_code=200, body=None, text=None, headers=None):
        self.status_code = status_code
        self._body = body if body is not None else {}
        self.text = text if text is not None else "{}"
        self.headers = headers or {}

    def json(self):
        return self._body


class ScriptedCall:
    """A zero-arg wire call replaying a script of responses/exceptions."""

    def __init__(self, *script):
        self._script = list(script)
        self.calls = 0

    def __call__(self):
        self.calls += 1
        item = self._script[min(self.calls - 1, len(self._script) - 1)]
        if isinstance(item, BaseException):
            raise item
        return item


OK = FakeResponse(200, {"markets": ["fine"]})
SERVER_ERROR = FakeResponse(500, text="oops")
FORBIDDEN = FakeResponse(403, text="Access Denied")
NOT_FOUND = FakeResponse(404, text="not found")
THROTTLED = FakeResponse(429, text="slow down")


def _run(call, profile=RETRY_BACKGROUND, **kw):
    return request_with_retry(call, profile=profile, book="testbook",
                              stage="structure", **kw)


# --------------------------------------------------------------------------- #
# Classification                                                               #
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("status,expected", [
    (200, False), (400, False), (401, False), (403, False), (404, False),
    (422, False), (429, True), (500, True), (502, True), (503, True),
])
def test_status_retry_classification(status, expected):
    assert is_retryable_status(status) is expected


def test_curl_cffi_transport_errors_are_retryable():
    """The OSError base is what makes this library-agnostic — if curl_cffi
    ever stops subclassing it, this test is the tripwire."""
    from curl_cffi.requests.exceptions import ConnectionError as CurlConnErr
    from curl_cffi.requests.exceptions import RequestException, Timeout

    assert issubclass(RequestException, OSError)
    for exc in (CurlConnErr("reset"), Timeout("read timed out"),
                RequestException("boom"), ConnectionError("dns"),
                TimeoutError("slow"), OSError("socket")):
        assert is_retryable_exception(exc), exc


def test_programming_errors_are_not_retryable():
    """A bug in our own code must fail fast, not three times slowly."""
    for exc in (TypeError("bad kwarg"), AttributeError("no attr"),
                KeyError("missing"), ValueError("nope")):
        assert not is_retryable_exception(exc)


# --------------------------------------------------------------------------- #
# Recovery — the core of the ticket                                            #
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("profile", [RETRY_BACKGROUND, RETRY_LIVE])
def test_transient_500_then_success_recovers(profile, recorded_sleeps):
    call = ScriptedCall(SERVER_ERROR, OK)
    resp = _run(call, profile)
    assert resp is OK
    assert call.calls == 2
    assert len(recorded_sleeps) == 1


@pytest.mark.parametrize("profile", [RETRY_BACKGROUND, RETRY_LIVE])
def test_transient_connection_error_then_success_recovers(profile,
                                                          recorded_sleeps):
    call = ScriptedCall(ConnectionError("[Errno 8] nodename nor servname"), OK)
    assert _run(call, profile) is OK
    assert call.calls == 2
    assert len(recorded_sleeps) == 1


def test_success_on_first_attempt_never_sleeps(recorded_sleeps):
    call = ScriptedCall(OK)
    assert _run(call) is OK
    assert call.calls == 1
    assert recorded_sleeps == []


def test_two_transient_failures_then_success_recovers(recorded_sleeps):
    """BACKGROUND's third attempt is a real one, not decoration."""
    call = ScriptedCall(SERVER_ERROR, ConnectionError("reset"), OK)
    assert _run(call) is OK
    assert call.calls == 3
    assert len(recorded_sleeps) == 2


# --------------------------------------------------------------------------- #
# No retry on real verdicts                                                    #
# --------------------------------------------------------------------------- #

def test_403_is_not_retried_and_surfaces_immediately(recorded_sleeps):
    call = ScriptedCall(FORBIDDEN, OK)      # OK never reached
    resp = _run(call)
    assert call.calls == 1
    assert recorded_sleeps == []
    # The helper hands the response back; check_response is what raises,
    # so BookTransportError stays constructed in exactly one place.
    with pytest.raises(BookTransportError) as exc:
        check_response("testbook", "structure", resp)
    assert exc.value.status_code == 403


def test_404_is_not_retried_and_keeps_the_allow_404_carve_out(recorded_sleeps):
    """#33's carve-out must survive: one delisted event is a skip, not a
    dead book — and never costs three attempts."""
    call = ScriptedCall(NOT_FOUND)
    resp = _run(call)
    assert call.calls == 1
    assert recorded_sleeps == []
    assert check_response("testbook", "structure", resp, allow_404=True) is False


def test_non_retryable_exception_raises_on_first_attempt(recorded_sleeps):
    call = ScriptedCall(TypeError("unexpected kwarg"), OK)
    with pytest.raises(BookTransportError) as exc:
        _run(call)
    assert call.calls == 1
    assert recorded_sleeps == []
    assert isinstance(exc.value.cause, TypeError)


# --------------------------------------------------------------------------- #
# Exhaustion                                                                   #
# --------------------------------------------------------------------------- #

def test_persistent_500_exhausts_attempts_then_check_response_raises(
        recorded_sleeps):
    call = ScriptedCall(SERVER_ERROR)
    resp = _run(call)
    assert call.calls == RETRY_BACKGROUND.max_attempts
    assert len(recorded_sleeps) == RETRY_BACKGROUND.max_attempts - 1
    with pytest.raises(BookTransportError) as exc:
        check_response("testbook", "structure", resp)
    assert exc.value.status_code == 500


def test_persistent_connection_error_raises_with_the_last_cause(
        recorded_sleeps):
    boom = ConnectionError("nodename nor servname provided")
    call = ScriptedCall(boom)
    with pytest.raises(BookTransportError) as exc:
        _run(call)
    assert call.calls == RETRY_BACKGROUND.max_attempts
    assert exc.value.cause is boom
    assert exc.value.book == "testbook"
    assert exc.value.stage == "structure"


# --------------------------------------------------------------------------- #
# 429 / Retry-After                                                            #
# --------------------------------------------------------------------------- #

def test_429_backs_off_and_retries(recorded_sleeps):
    call = ScriptedCall(THROTTLED, OK)
    assert _run(call) is OK
    assert call.calls == 2
    assert recorded_sleeps == [RETRY_BACKGROUND.base_delay_sec]


def test_429_honors_a_sane_retry_after(recorded_sleeps):
    throttled = FakeResponse(429, headers={"Retry-After": "2"})
    call = ScriptedCall(throttled, OK)
    assert _run(call) is OK
    assert recorded_sleeps == [2.0]


def test_429_ignores_a_retry_after_that_blows_the_budget(recorded_sleeps):
    """A book asking for 30s is not "sane" on a 1s LIVE budget — honoring it
    would trade a lost book for a lost quote, which is the worse failure."""
    throttled = FakeResponse(429, headers={"Retry-After": "30"})
    call = ScriptedCall(throttled, OK)
    assert _run(call, RETRY_LIVE) is OK
    assert recorded_sleeps == [RETRY_LIVE.base_delay_sec]
    assert sum(recorded_sleeps) <= RETRY_LIVE.max_total_delay_sec


def test_429_ignores_an_unparseable_retry_after(recorded_sleeps):
    throttled = FakeResponse(429, headers={"Retry-After": "Wed, 21 Oct 2026 07:28:00 GMT"})
    call = ScriptedCall(throttled, OK)
    assert _run(call) is OK
    assert recorded_sleeps == [RETRY_BACKGROUND.base_delay_sec]


def test_retry_after_on_a_500_is_ignored(recorded_sleeps):
    """Retry-After is only meaningful on a 429 here; a 5xx carrying one
    still follows our own backoff."""
    err = FakeResponse(500, headers={"Retry-After": "2"})
    call = ScriptedCall(err, OK)
    assert _run(call) is OK
    assert recorded_sleeps == [RETRY_BACKGROUND.base_delay_sec]


# --------------------------------------------------------------------------- #
# Budget — the number quoted in request_with_retry's docstring                 #
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("profile,ceiling", [
    (RETRY_BACKGROUND, 5.63),
    (RETRY_LIVE, 0.63),
])
def test_worst_case_added_delay_stays_within_profile(profile, ceiling):
    """Pin the documented worst case with MAXIMUM jitter and an endpoint
    that never recovers. These numbers are quoted in the helper's docstring
    and in mlb_sgp/README.md."""
    slept: list[float] = []
    call = ScriptedCall(SERVER_ERROR)
    request_with_retry(call, profile=profile, book="testbook",
                       stage="price", sleep_fn=slept.append,
                       jitter_fn=lambda: 1.0)          # worst-case jitter
    assert call.calls == profile.max_attempts
    assert sum(slept) <= ceiling
    assert sum(slept) <= profile.max_total_delay_sec


def test_live_profile_is_strictly_faster_than_background():
    assert RETRY_LIVE.max_attempts < RETRY_BACKGROUND.max_attempts
    assert RETRY_LIVE.max_total_delay_sec < RETRY_BACKGROUND.max_total_delay_sec


def test_sleep_budget_caps_a_profile_with_an_absurd_base_delay():
    """max_total_delay_sec is a HARD ceiling, independent of the schedule."""
    greedy = RetryProfile(name="greedy", max_attempts=5, base_delay_sec=100.0,
                          backoff_multiplier=2.0, max_total_delay_sec=1.0)
    slept: list[float] = []
    call = ScriptedCall(SERVER_ERROR)
    request_with_retry(call, profile=greedy, book="testbook", stage="price",
                       sleep_fn=slept.append, jitter_fn=lambda: 0.5)
    assert sum(slept) <= greedy.max_total_delay_sec
    # Budget exhausted after the first sleep, so it stops early rather than
    # sleeping 0 four more times.
    assert call.calls == 2


# --------------------------------------------------------------------------- #
# body_ok — Caesars' 200-with-HTML WAF challenge                               #
# --------------------------------------------------------------------------- #

def _looks_like_json(resp) -> bool:
    return (getattr(resp, "text", "") or "").strip().startswith(("{", "["))


def test_200_with_an_unusable_body_is_retried(recorded_sleeps):
    """A WAF challenge is 200 + HTML. Without body_ok this book would lose
    the one retry that actually recovers it."""
    challenge = FakeResponse(200, text="<html>AWS WAF challenge</html>")
    call = ScriptedCall(challenge, OK)
    assert _run(call, body_ok=_looks_like_json) is OK
    assert call.calls == 2


def test_200_with_an_unusable_body_eventually_raises_via_json_or_raise(
        recorded_sleeps):
    challenge = FakeResponse(200, text="<html>AWS WAF challenge</html>")
    challenge.json = lambda: (_ for _ in ()).throw(ValueError("not json"))
    call = ScriptedCall(challenge)
    resp = _run(call, body_ok=_looks_like_json)
    assert call.calls == RETRY_BACKGROUND.max_attempts
    with pytest.raises(BookTransportError):
        json_or_raise("testbook", "structure", resp)


def test_body_ok_is_not_consulted_for_non_200(recorded_sleeps):
    """A 403's body is irrelevant — it must not be retried just because it
    fails the JSON predicate."""
    call = ScriptedCall(FORBIDDEN, OK)
    _run(call, body_ok=_looks_like_json)
    assert call.calls == 1
