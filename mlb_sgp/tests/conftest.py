"""Shared fixtures for the mlb_sgp test suite.

The one job here is keeping issue #34's retry backoff out of the suite's wall
clock. Dozens of existing tests drive a client with a permanently-failing fake
session (``FakeSession(SERVER_ERROR)``, ``BoomSession()``); with real sleeps
each of those would now cost ~5.6s of backoff on the BACKGROUND profile.
"""
from __future__ import annotations

import pytest

from mlb_sgp import _shared


@pytest.fixture(autouse=True)
def recorded_sleeps(monkeypatch):
    """Neutralize retry backoff and RECORD it instead.

    Autouse, so no test ever sleeps for real. Tests that care about the
    backoff schedule take ``recorded_sleeps`` as an argument and read the
    list of durations that WOULD have been slept, in order.

    ``_jitter`` is pinned to 0.5 (the midpoint of ``random.random()``'s
    range) so the exponential schedule is exactly reproducible: the jitter
    factor ``0.75 + 0.5 * 0.5`` collapses to 1.0, leaving the bare
    ``base * multiplier ** attempt``.
    """
    slept: list[float] = []
    monkeypatch.setattr(_shared, "_sleep", slept.append)
    monkeypatch.setattr(_shared, "_jitter", lambda: 0.5)
    return slept
