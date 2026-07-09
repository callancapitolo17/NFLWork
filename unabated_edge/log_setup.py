import logging, sys
from logging.handlers import RotatingFileHandler
from unabated_edge import config


def setup_logging():
    log = logging.getLogger("unabated_edge")
    if any(getattr(h, "_m", False) for h in log.handlers):
        return log
    log.setLevel(logging.INFO)
    fh = RotatingFileHandler(config.LOG_PATH, maxBytes=10_000_000, backupCount=5)
    fh._m = True
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s"))
    log.addHandler(fh)
    if sys.stderr.isatty():
        sh = logging.StreamHandler()
        sh._m = True
        sh.setFormatter(fh.formatter)
        log.addHandler(sh)
    return log
