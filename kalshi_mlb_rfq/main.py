"""Kalshi MLB RFQ Bot — autonomous taker daemon."""

import argparse
import os
import signal
import threading
import time
from datetime import datetime, timezone

from kalshi_mlb_rfq import config, db, notify, rfq_client
from kalshi_mlb_rfq.config import KILL_FILE

VERSION = "0.1.0"

_running = threading.Event()
_running.set()
ACCEPT_LOCK = threading.Lock()


def _signal_handler(_sig, _frame):
    _running.clear()


def _phantom_rfq_cleanup():
    """At startup, cancel any RFQs on Kalshi that aren't tracked in our live_rfqs."""
    try:
        kalshi_open = rfq_client.list_open_rfqs(config.KALSHI_USER_ID or "")
    except Exception as e:
        print(f"  startup: list_open_rfqs failed: {e}", flush=True)
        return

    if not kalshi_open:
        print("  startup: no open RFQs on Kalshi", flush=True)
        return

    with db.connect(read_only=True) as con:
        ours = {r[0] for r in con.execute(
            "SELECT rfq_id FROM live_rfqs WHERE status='open'"
        ).fetchall()}

    for rfq in kalshi_open:
        rid = rfq.get("id")
        if rid and rid not in ours:
            try:
                rfq_client.delete_rfq(rid)
                print(f"  startup: cancelled phantom rfq {rid}", flush=True)
            except Exception as e:
                print(f"  startup: failed to cancel phantom {rid}: {e}", flush=True)


def main_loop(dry_run: bool):
    db.init_database()
    sid = db.start_session(pid=os.getpid(), dry_run=dry_run, version=VERSION)
    print(f"=== Kalshi MLB RFQ Bot — session {sid} (dry_run={dry_run}) ===", flush=True)
    _phantom_rfq_cleanup()

    last_rfq_refresh = 0.0
    last_quote_poll = 0.0
    last_risk_sweep = 0.0
    last_pipeline = 0.0
    last_heartbeat = 0.0

    try:
        while _running.is_set():
            now = time.time()

            if KILL_FILE.exists():
                notify.halt("kill_switch")
                time.sleep(config.RISK_SWEEP_SEC)
                continue

            if now - last_rfq_refresh >= config.RFQ_REFRESH_SEC:
                # Wired in T21
                last_rfq_refresh = now

            if now - last_quote_poll >= config.QUOTE_POLL_SEC:
                # Wired in T22
                last_quote_poll = now

            if now - last_risk_sweep >= config.RISK_SWEEP_SEC:
                # Wired in T25
                last_risk_sweep = now

            if now - last_pipeline >= config.PIPELINE_REFRESH_SEC:
                # Wired in T24
                last_pipeline = now

            if now - last_heartbeat >= 60:
                print(f"  [HB] {datetime.now(timezone.utc).isoformat()} alive", flush=True)
                last_heartbeat = now

            time.sleep(0.5)
    finally:
        with db.connect(read_only=True) as con:
            live = [r[0] for r in con.execute(
                "SELECT rfq_id FROM live_rfqs WHERE status='open'"
            ).fetchall()]
        for rid in live:
            try:
                rfq_client.delete_rfq(rid)
            except Exception:
                pass
        db.end_session(sid)
        print("=== shutdown complete ===", flush=True)


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dry-run", action="store_true",
                        help="Run full loop without calling accept_quote.")
    args = parser.parse_args()

    signal.signal(signal.SIGINT, _signal_handler)
    signal.signal(signal.SIGTERM, _signal_handler)
    main_loop(dry_run=args.dry_run)


if __name__ == "__main__":
    cli()
