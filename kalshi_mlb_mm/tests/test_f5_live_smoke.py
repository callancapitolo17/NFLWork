"""Issue #85 live smoke: one REAL on-demand flight for an F5 spread+total combo.

Read-only against the books (event lists, structure fetches, price calls) —
no Kalshi traffic, no quoting, no DB writes (SGPService is constructed with
health_db_path=None). Opt-in because it needs a live MLB slate with F5
markets posted (books post F5 near lineup time):

    F5_LIVE_SMOKE=1 python3 -m pytest kalshi_mlb_mm/tests/test_f5_live_smoke.py -s

Asserts >= 2 book fairs (the MIN_AGREEING_BOOKS quorum) for at least one
standard F5 spread+total combo on an upcoming game, each inside the book's
on-demand deadline, and prints the per-book fairs + probit dispersion for
eyeball comparison against the sigma_z gate.
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
    os.environ.get("F5_LIVE_SMOKE") != "1",
    reason="live book traffic — opt in with F5_LIVE_SMOKE=1")

# Standard F5 main lines: runline ±0.5, total around 4.5. Candidates are
# tried in order until one gets a >= 2-book quorum (structure fetches are
# cached across candidates, so retries only cost price calls).
CANDIDATE_COMBOS = [
    (-0.5, 4.5), (0.5, 4.5), (-0.5, 5.5), (0.5, 5.5),
    (-0.5, 3.5), (0.5, 3.5),
]
MIN_BOOKS = 2
OVERALL_BUDGET_SEC = 30.0


def _pick_game():
    """Soonest not-yet-started game from Novig's public event list (clean
    Odds-API-style team names + start times)."""
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


def test_f5_spread_total_live_flight():
    from kalshi_common.legset import CanonicalLeg
    from kalshi_common.sgp_service import SGPService
    from mlb_sgp._shared import GameRef

    ev, start_dt = _pick_game()
    game = GameRef(game_id="f5-smoke", home_team=ev.home_team,
                   away_team=ev.away_team, commence_time=start_dt)
    print(f"\ngame: {ev.away_team} @ {ev.home_team}  start={start_dt}")

    svc = SGPService()                    # all 6 books, no DB, no alerter
    try:
        for spread_line, total_line in CANDIDATE_COMBOS:
            legs = [
                CanonicalLeg("f5-smoke", "spread", spread_line,
                             "home" if spread_line < 0 else "away", "F5"),
                CanonicalLeg("f5-smoke", "total", total_line, "over", "F5"),
            ]
            with ThreadPoolExecutor(max_workers=len(svc.books)) as pool:
                futs = {b: pool.submit(svc.price_on_demand, b, game, legs)
                        for b in svc.books}
                results = {}
                for b, fut in futs.items():
                    try:
                        results[b] = fut.result(timeout=OVERALL_BUDGET_SEC)
                    except Exception as e:      # noqa: BLE001 — smoke report
                        print(f"  {b}: EXC {type(e).__name__}: {e}")
                        results[b] = None
            fairs = {b: r for b, r in results.items() if r is not None}
            print(f"combo spread {spread_line:+.1f} / over {total_line}: "
                  f"{len(fairs)} book(s)")
            for b, r in fairs.items():
                print(f"  {b}: fair={r.fair:.4f} route={r.route} "
                      f"latency={r.latency_sec:.2f}s")
            if len(fairs) >= MIN_BOOKS:
                zs = [statistics.NormalDist().inv_cdf(r.fair)
                      for r in fairs.values()]
                sigma_z = statistics.stdev(zs) if len(zs) >= 2 else math.nan
                print(f"  probit dispersion sigma_z={sigma_z:.4f} "
                      f"(gate SIGMA_Z_MAX=0.07)")
                for r in fairs.values():
                    assert r.latency_sec <= svc.on_demand_deadline_sec + 5.0
                return                       # quorum reached — smoke passes
        pytest.fail(f"no candidate F5 combo reached {MIN_BOOKS} book fairs — "
                    "are F5 markets posted yet? (books post near lineup time)")
    finally:
        svc.close()
