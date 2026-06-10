"""Tests for kalshi_mlb_mm/log_setup.py."""
import logging
from logging.handlers import RotatingFileHandler

from kalshi_mlb_mm.log_setup import setup_logging


def _mm_managed_handlers(root):
    return [h for h in root.handlers if getattr(h, "_mm_managed", False)]


def _cleanup(root):
    for h in list(root.handlers):
        if getattr(h, "_mm_managed", False):
            h.close()
            root.removeHandler(h)


def test_setup_logging_installs_rotating_and_stream_handlers(tmp_path):
    log_file = tmp_path / "bot.log"
    # console=True forces the stream handler (pytest's stderr is not a TTY,
    # so the auto-detect default would otherwise omit it).
    logger = setup_logging(log_path=log_file, max_bytes=1024,
                           backup_count=3, level="DEBUG", console=True)
    try:
        managed = _mm_managed_handlers(logger)
        handler_types = {type(h) for h in managed}
        assert RotatingFileHandler in handler_types
        assert logging.StreamHandler in handler_types
        rfh = next(h for h in managed if isinstance(h, RotatingFileHandler))
        assert rfh.maxBytes == 1024
        assert rfh.backupCount == 3
        assert logger.level == logging.DEBUG
    finally:
        _cleanup(logger)


def test_setup_logging_is_idempotent(tmp_path):
    log_file = tmp_path / "bot.log"
    setup_logging(log_path=log_file, console=True)
    logger = setup_logging(log_path=log_file, console=True)
    try:
        managed = _mm_managed_handlers(logger)
        rfh_count = sum(1 for h in managed if isinstance(h, RotatingFileHandler))
        # type(h) is StreamHandler (not isinstance) — RotatingFileHandler is a
        # StreamHandler subclass, so isinstance would double-count it.
        stream_count = sum(1 for h in managed if type(h) is logging.StreamHandler)
        assert rfh_count == 1
        assert stream_count == 1
    finally:
        _cleanup(logger)


def test_setup_logging_omits_console_when_not_tty(tmp_path):
    """The backgrounded-daemon case: no stderr StreamHandler, so a
    `>> bot.log 2>&1` launch can't double-log into bot.log."""
    log_file = tmp_path / "bot.log"
    logger = setup_logging(log_path=log_file, console=False)
    try:
        managed = _mm_managed_handlers(logger)
        rfh_count = sum(1 for h in managed if isinstance(h, RotatingFileHandler))
        stream_count = sum(1 for h in managed if type(h) is logging.StreamHandler)
        assert rfh_count == 1
        assert stream_count == 0
    finally:
        _cleanup(logger)


def test_setup_logging_rotating_handler_path(tmp_path):
    """The rotating file handler should write to the specified path."""
    log_file = tmp_path / "bot.log"
    logger = setup_logging(log_path=log_file, console=False)
    try:
        managed = _mm_managed_handlers(logger)
        rfh = next(h for h in managed if isinstance(h, RotatingFileHandler))
        import os
        assert os.path.abspath(rfh.baseFilename) == str(log_file.resolve())
    finally:
        _cleanup(logger)
