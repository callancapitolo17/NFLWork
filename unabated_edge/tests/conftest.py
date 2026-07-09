import logging

import pytest
from unabated_edge import storage


@pytest.fixture(autouse=True)
def _clean_storage_buffer():
    storage._BUFFER.clear()
    yield
    storage._BUFFER.clear()


@pytest.fixture(autouse=True, scope="session")
def _no_live_log_pollution():
    """Detach the rotating bot.log handler during tests: test EDGE/WARNING lines
    in the live log masquerade as real engine output (found during monitoring).
    caplog still works — it captures via propagation, not via these handlers."""
    log = logging.getLogger("unabated_edge")
    saved = log.handlers[:]
    log.handlers = []
    yield
    log.handlers = saved
