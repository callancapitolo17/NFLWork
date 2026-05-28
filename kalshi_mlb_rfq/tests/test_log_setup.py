"""Tests for kalshi_mlb_rfq/log_setup.py."""
import logging
from logging.handlers import RotatingFileHandler

from kalshi_mlb_rfq.log_setup import setup_logging


def test_setup_logging_installs_rotating_and_stream_handlers(tmp_path):
    log_file = tmp_path / "bot.log"
    logger = setup_logging(log_path=log_file, max_bytes=1024,
                           backup_count=3, level="DEBUG")
    handler_types = {type(h) for h in logger.handlers}
    assert RotatingFileHandler in handler_types
    assert logging.StreamHandler in handler_types
    rfh = next(h for h in logger.handlers if isinstance(h, RotatingFileHandler))
    assert rfh.maxBytes == 1024
    assert rfh.backupCount == 3
    assert logger.level == logging.DEBUG


def test_setup_logging_is_idempotent(tmp_path):
    log_file = tmp_path / "bot.log"
    setup_logging(log_path=log_file)
    logger = setup_logging(log_path=log_file)
    rfh_count = sum(1 for h in logger.handlers if isinstance(h, RotatingFileHandler))
    assert rfh_count == 1
