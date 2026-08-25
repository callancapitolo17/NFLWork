"""Central logging config — mirrors kalshi_mlb_mm/log_setup.py in structure
and contract (rotating bot.log, TTY-aware console handler, idempotent)."""

import logging
import sys
from logging.handlers import RotatingFileHandler
from pathlib import Path

from kalshi_rfi import config

_FORMAT = "%(asctime)s %(levelname)s %(name)s — %(message)s"


def setup_logging(log_path: Path | None = None,
                  max_bytes: int | None = None,
                  backup_count: int | None = None,
                  level: str | None = None,
                  console: bool | None = None) -> logging.Logger:
    log_path = Path(log_path) if log_path else config.LOG_PATH
    max_bytes = max_bytes if max_bytes is not None else config.LOG_ROTATE_MAX_BYTES
    backup_count = (backup_count if backup_count is not None
                    else config.LOG_ROTATE_BACKUPS)
    level = level or config.LOG_LEVEL
    if console is None:
        console = sys.stderr.isatty()

    root = logging.getLogger()
    root.setLevel(getattr(logging, level.upper(), logging.INFO))

    for h in list(root.handlers):
        if getattr(h, "_rfi_managed", False):
            h.close()
            root.removeHandler(h)

    log_path.parent.mkdir(parents=True, exist_ok=True)
    fmt = logging.Formatter(_FORMAT)

    file_handler = RotatingFileHandler(
        log_path, maxBytes=max_bytes, backupCount=backup_count)
    file_handler.setFormatter(fmt)
    file_handler._rfi_managed = True
    root.addHandler(file_handler)

    if console:
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(fmt)
        stream_handler._rfi_managed = True
        root.addHandler(stream_handler)

    return root
