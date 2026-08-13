"""Issue #87 live smoke: one REAL on-demand flight for an {RFI + FG ML} combo.

Read-only against the books (event lists, structure fetches, price calls) —
no Kalshi traffic, no quoting, no DB writes. Opt-in because it needs a live
MLB slate with 1st-inning markets posted (books post them near lineup time):

    RFI_LIVE_SMOKE=1 python3 -m pytest kalshi_mlb_mm/tests/test_rfi_live_smoke.py -s

Asserts >= 2 book fairs (the MIN_AGREEING_BOOKS quorum) for the I1-over-0.5 +
home-ML combo on an upcoming game, and prints per-book fairs + probit
dispersion for eyeball comparison against the sigma_z gate.
"""
import math
import os
import statistics
import sys
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "mlb_sgp"))

pytestmark = pytest.mark.skipif(
    os.environ.get("RFI_LIVE_SMOKE") != "1",
    reason="live book traffic — opt in with RFI_LIVE_SMOKE=1")

MIN_BOOKS = 2
OVERALL_BUDGET_SEC = 30.0


def _pick_game():
    """Soonest not-yet-started game from Novig's public event list."""
    from novig_client import NovigClient
    events = NovigClient().list_events()
    now = datetime.now(timezone.utc)

    def _start_dt(e):
        st = e.start_time
        if isinstance(st, str):
            st = datetime.fromisoformat(st.replace("Z", "+00:00"))
        return st

    upcoming = [(e, _start_dt(e)) for e in events if e.start_time]
    upcoming = [(e, dt) for e, dt in upcoming if dt > now]
    assert upcoming, "no upcoming MLB games on the board — run pre-slate"
    return min(upcoming, key=lambda pair: pair[1])


def test_rfi_plus_ml_live_flight():
    from kalshi_common.legset import CanonicalLeg
    from kalshi_common.sgp_service import SGPService
    from mlb_sgp._shared import GameRef

    ev, start_dt = _pick_game()
    game = GameRef(game_id="rfi-smoke", home_team=ev.home_team,
                   away_team=ev.away_team, commence_time=start_dt)
    print(f"\ngame: {ev.away_team} @ {ev.home_team}  start={start_dt}")

    legs = [
        CanonicalLeg("rfi-smoke", "total", 0.5, "over", "I1"),   # RFI yes
        CanonicalLeg("rfi-smoke", "ml", None, "home"),
    ]
    svc = SGPService()                    # all 6 books, no DB, no alerter
    try:
        with ThreadPoolExecutor(max_workers=len(svc.books)) as pool:
            futs = {b: pool.submit(svc.price_on_demand, b, game, legs)
                    for b in svc.books}
            results = {}
            for b, fut in futs.items():
                try:
                    results[b] = fut.result(timeout=OVERALL_BUDGET_SEC)
                except Exception as e:          # noqa: BLE001 — smoke report
                    print(f"  {b}: EXC {type(e).__name__}: {e}")
                    results[b] = None
        fairs = {b: r for b, r in results.items() if r is not None}
        for b, r in fairs.items():
            print(f"  {b}: fair={r.fair:.4f} route={r.route} "
                  f"latency={r.latency_sec:.2f}s")
        if len(fairs) >= 2:
            zs = [statistics.NormalDist().inv_cdf(r.fair)
                  for r in fairs.values()]
            sigma_z = statistics.stdev(zs) if len(zs) >= 2 else math.nan
            print(f"  probit dispersion sigma_z={sigma_z:.4f} "
                  f"(gate SIGMA_Z_MAX=0.07)")
        assert len(fairs) >= MIN_BOOKS, (
            f"only {len(fairs)} book(s) priced the RFI combo — are 1st-inning "
            "markets posted yet? (books post near lineup time)")
    finally:
        svc.close()
