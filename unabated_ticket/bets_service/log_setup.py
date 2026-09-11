"""Logging config — mirrors kalshi_rfi/log_setup.py (rotating file, TTY-aware
console handler, idempotent). Never logs credentials: the Kalshi key id and
private key path are read by config and passed to auth_client only."""
import logging
import sys
from logging.handlers import RotatingFileHandler
from pathlib import Path

from unabated_ticket.bets_service import config

_FORMAT = "%(asctime)s %(levelname)s %(name)s — %(message)s"


def setup_logging(log_path: Path | None = None, console: bool | None = None) -> logging.Logger:
    log_path = Path(log_path) if log_path else config.LOG_PATH
    if console is None:
        console = sys.stderr.isatty()

    root = logging.getLogger()
    root.setLevel(getattr(logging, config.LOG_LEVEL.upper(), logging.INFO))
    for handler in list(root.handlers):
        if getattr(handler, "_bets_service_managed", False):
            handler.close()
            root.removeHandler(handler)

    log_path.parent.mkdir(parents=True, exist_ok=True)
    formatter = logging.Formatter(_FORMAT)
    file_handler = RotatingFileHandler(
        log_path, maxBytes=config.LOG_ROTATE_MAX_BYTES, backupCount=config.LOG_ROTATE_BACKUPS)
    file_handler.setFormatter(formatter)
    file_handler._bets_service_managed = True
    root.addHandler(file_handler)
    if console:
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(formatter)
        stream_handler._bets_service_managed = True
        root.addHandler(stream_handler)
    return root
