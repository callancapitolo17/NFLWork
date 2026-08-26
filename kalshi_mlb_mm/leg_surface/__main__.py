"""Standalone leg-surface ingest: ``python -m kalshi_mlb_mm.leg_surface``.

Runs the ingest loop and NOTHING else — no RFQ discovery, no quoting, no
orders. This is issue #96's acceptance harness ("loop runs standalone for an
hour without touching quote behaviour") and the way to populate the surface
without starting the maker.

Side effects: live HTTP to Kalshi and the books; writes
``kalshi_mlb_mm/kalshi_mlb_mm_surface.duckdb``. Reads no other bot DB and
writes none.
"""
import argparse
import logging
import signal
import sys
import time

from kalshi_common import auth_client
from kalshi_mlb_mm import config
from kalshi_mlb_mm.leg_surface.runner import SurfaceIngest
from kalshi_mlb_mm.log_setup import setup_logging

log = logging.getLogger(__name__)
_running = True


def _handle_signal(_sig, _frame):
    global _running
    _running = False


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run the maker's leg-surface ingest loop standalone.")
    parser.add_argument("--minutes", type=float, default=None,
                        help="stop after N minutes (default: until SIGTERM)")
    args = parser.parse_args()

    setup_logging()
    auth_client.configure(api_key_id=config.KALSHI_API_KEY_ID,
                          private_key_path=config.KALSHI_PRIVATE_KEY_PATH,
                          base_url=config.KALSHI_BASE_URL,
                          project_root=config.PROJECT_ROOT)
    signal.signal(signal.SIGINT, _handle_signal)
    signal.signal(signal.SIGTERM, _handle_signal)

    ingest = SurfaceIngest()
    ingest.start()
    log.info("leg-surface standalone: db=%s workers=%s",
             config.SURFACE_DB, ingest.worker_specs())
    deadline = time.monotonic() + args.minutes * 60 if args.minutes else None
    try:
        while _running and (deadline is None or time.monotonic() < deadline):
            time.sleep(0.5)
    finally:
        ingest.stop()
    log.info("leg-surface standalone: %d rows in memory at exit",
             ingest.surface.row_count())
    return 0


if __name__ == "__main__":
    sys.exit(main())
