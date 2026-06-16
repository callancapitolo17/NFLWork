"""Tests for kalshi_mlb_rfq/log_setup.py."""
import logging
from logging.handlers import RotatingFileHandler

import pytest

from kalshi_mlb_rfq.log_setup import setup_logging


@pytest.fixture(autouse=True)
def _reset_root_logger():
    """Remove our managed handlers after each test so they don't leak
    (open fds pointing at deleted tmp_path dirs) across the session."""
    yield
    root = logging.getLogger()
    for h in list(root.handlers):
        if getattr(h, "_rfq_managed", False):
            h.close()
            root.removeHandler(h)


def test_setup_logging_installs_rotating_and_stream_handlers(tmp_path):
    log_file = tmp_path / "bot.log"
    # console=True forces the stream handler (pytest's stderr is not a TTY,
    # so the auto-detect default would otherwise omit it).
    logger = setup_logging(log_path=log_file, max_bytes=1024,
                           backup_count=3, level="DEBUG", console=True)
    handler_types = {type(h) for h in logger.handlers}
    assert RotatingFileHandler in handler_types
    assert logging.StreamHandler in handler_types
    rfh = next(h for h in logger.handlers if isinstance(h, RotatingFileHandler))
    assert rfh.maxBytes == 1024
    assert rfh.backupCount == 3
    assert logger.level == logging.DEBUG


def test_setup_logging_is_idempotent(tmp_path):
    log_file = tmp_path / "bot.log"
    setup_logging(log_path=log_file, console=True)
    logger = setup_logging(log_path=log_file, console=True)
    rfh_count = sum(1 for h in logger.handlers
                    if isinstance(h, RotatingFileHandler))
    # type(h) is StreamHandler (not isinstance) — RotatingFileHandler is a
    # StreamHandler subclass, so isinstance would double-count it.
    stream_count = sum(1 for h in logger.handlers
                       if type(h) is logging.StreamHandler)
    assert rfh_count == 1
    assert stream_count == 1


def test_setup_logging_omits_console_when_not_tty(tmp_path):
    """The backgrounded-daemon case: no stderr StreamHandler, so a
    `>> bot.log 2>&1` launch can't double-log into bot.log."""
    log_file = tmp_path / "bot.log"
    logger = setup_logging(log_path=log_file, console=False)
    rfh_count = sum(1 for h in logger.handlers
                    if isinstance(h, RotatingFileHandler))
    stream_count = sum(1 for h in logger.handlers
                       if type(h) is logging.StreamHandler)
    assert rfh_count == 1
    assert stream_count == 0
