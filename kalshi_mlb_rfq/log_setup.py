"""Central logging configuration for the Kalshi MLB RFQ bot.

Replaces ad-hoc print() with the stdlib logging module: levels (filterable
without code changes) + a size-capped RotatingFileHandler (so bot.log can
never grow unbounded). Call setup_logging() once at process start.
"""

import logging
from logging.handlers import RotatingFileHandler
from pathlib import Path

from kalshi_mlb_rfq import config

_FORMAT = "%(asctime)s %(levelname)s %(name)s: %(message)s"


def setup_logging(log_path: Path | None = None,
                  max_bytes: int | None = None,
                  backup_count: int | None = None,
                  level: str | None = None) -> logging.Logger:
    """Configure the root logger with a stdout handler + a rotating file
    handler. Idempotent: repeated calls do not stack handlers.

    Defaults come from config (LOG_PATH, LOG_MAX_BYTES, LOG_BACKUP_COUNT,
    LOG_LEVEL); args override for tests.
    """
    log_path = Path(log_path) if log_path else config.LOG_PATH
    max_bytes = max_bytes if max_bytes is not None else config.LOG_MAX_BYTES
    backup_count = (backup_count if backup_count is not None
                    else config.LOG_BACKUP_COUNT)
    level = level or config.LOG_LEVEL

    root = logging.getLogger()
    root.setLevel(getattr(logging, level.upper(), logging.INFO))

    # Idempotency: drop our previously-installed handlers before re-adding.
    # close() before removeHandler() so the rotating file's fd is released
    # (removeHandler alone leaks the descriptor until GC).
    for h in list(root.handlers):
        if getattr(h, "_rfq_managed", False):
            h.close()
            root.removeHandler(h)

    log_path.parent.mkdir(parents=True, exist_ok=True)
    fmt = logging.Formatter(_FORMAT)

    file_handler = RotatingFileHandler(
        log_path, maxBytes=max_bytes, backupCount=backup_count)
    file_handler.setFormatter(fmt)
    file_handler._rfq_managed = True
    root.addHandler(file_handler)

    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(fmt)
    stream_handler._rfq_managed = True
    root.addHandler(stream_handler)

    return root
