"""Central logging configuration for the Kalshi MLB MM (maker) bot.

Replaces ad-hoc print() with the stdlib logging module: levels (filterable
without code changes) + a size-capped RotatingFileHandler (so bot.log can
never grow unbounded). Call setup_logging() once at process start.

Mirrors kalshi_mlb_rfq/log_setup.py exactly in structure and contract.
"""

import logging
import sys
from logging.handlers import RotatingFileHandler
from pathlib import Path

from kalshi_mlb_mm import config

_FORMAT = "%(asctime)s %(levelname)s %(name)s — %(message)s"


def setup_logging(log_path: Path | None = None,
                  max_bytes: int | None = None,
                  backup_count: int | None = None,
                  level: str | None = None,
                  console: bool | None = None) -> logging.Logger:
    """Configure the root logger with a rotating file handler and, when
    appropriate, a console (stderr) handler. Idempotent: repeated calls do
    not stack handlers.

    Defaults come from config (LOG_PATH, LOG_ROTATE_MAX_BYTES,
    LOG_ROTATE_BACKUPS, LOG_LEVEL); args override for tests.

    `console` controls the stderr StreamHandler. When None (default) it is
    added only if stderr is a TTY — so an interactive / dry-run session gets
    console output, but a backgrounded daemon launched with `>> bot.log 2>&1`
    does NOT double-log: the RotatingFileHandler already writes bot.log, and
    the stderr redirect would otherwise write every line into bot.log a second
    time. Pass True/False to force it.
    """
    log_path = Path(log_path) if log_path else config.LOG_PATH
    max_bytes = max_bytes if max_bytes is not None else config.LOG_ROTATE_MAX_BYTES
    backup_count = (backup_count if backup_count is not None
                    else config.LOG_ROTATE_BACKUPS)
    level = level or config.LOG_LEVEL
    if console is None:
        console = sys.stderr.isatty()

    root = logging.getLogger()
    root.setLevel(getattr(logging, level.upper(), logging.INFO))

    # Idempotency: drop our previously-installed handlers before re-adding.
    # close() before removeHandler() so the rotating file's fd is released
    # (removeHandler alone leaks the descriptor until GC).
    for h in list(root.handlers):
        if getattr(h, "_mm_managed", False):
            h.close()
            root.removeHandler(h)

    log_path.parent.mkdir(parents=True, exist_ok=True)
    fmt = logging.Formatter(_FORMAT)

    file_handler = RotatingFileHandler(
        log_path, maxBytes=max_bytes, backupCount=backup_count)
    file_handler.setFormatter(fmt)
    file_handler._mm_managed = True
    root.addHandler(file_handler)

    if console:
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(fmt)
        stream_handler._mm_managed = True
        root.addHandler(stream_handler)

    return root
