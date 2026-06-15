"""Target-level parallelism tests for mlb_sgp/draftkings.py::price_sgps.

Strategy: monkeypatch the legacy helper functions on the real
scraper_draftkings_sgp module (price_sgps imports them lazily at call
time, so patching module attributes works), then drive price_sgps with
a many-target fixture and assert the complete row set comes back —
identical content to what the sequential loop produced, regardless of
completion order.
"""
import sys
import threading
import time
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import MagicMock

# price_sgps does `from scraper_draftkings_sgp import ...` (top-level
# module name) — only resolvable with mlb_sgp/ itself on sys.path.
MLB_SGP_DIR = Path(__file__).resolve().parents[1]
if str(MLB_SGP_DIR) not in sys.path:
    sys.path.insert(0, str(MLB_SGP_DIR))

import scraper_draftkings_sgp as legacy  # noqa: E402

from mlb_sgp import draftkings  # noqa: E402
from mlb_sgp._shared import TargetLine  # noqa: E402

CT = datetime(2026, 6, 8, 23, 0, tzinfo=timezone.utc)


def _mk_target(game_id, spread, total):
    return TargetLine(
        game_id=game_id, home_team="New York Yankees",
        away_team="Boston Red Sox", commence_time=CT,
        period="FG", spread=spread, total=total,
    )


def _wire_fixture(monkeypatch, n_targets, fail_index=None, call_sleep=0.0,
                  concurrency_tracker=None):
    """One game, n_targets distinct (spread, total) lines, all offered.

    fail_index: calculate_sgp raises for that target's total sel-ids —
    its 4 combos must vanish without affecting other targets.
    concurrency_tracker: dict with 'lock', 'live', 'max' — records the
    peak number of in-flight calculate_sgp calls.
    """
    spreads, totals, targets = {}, {}, []
    for i in range(n_targets):
        spread, total = -1.5 - i, 7.5 + i
        targets.append(_mk_target("g1", spread, total))
        spreads[("N", abs(spread), "1")] = [f"H{i}"]
        spreads[("P", abs(spread), "3")] = [f"A{i}"]
        totals[("O", total)] = [f"O{i}"]
        totals[("U", total)] = [f"U{i}"]
    sel_ids_all = {
        "fg": {"spreads": spreads, "totals": totals, "canonical": {"M"}},
        "f5": {"spreads": {}, "totals": {}},
    }
    monkeypatch.setattr(legacy, "fetch_dk_events",
                        lambda s: [{"dk_event_id": "e1"}])
    monkeypatch.setattr(legacy, "match_events",
                        lambda evs, tgt: [{"game_id": "g1", "dk_event_id": "e1"}])
    monkeypatch.setattr(legacy, "fetch_main_market_nums", lambda s, eid: {})
    monkeypatch.setattr(legacy, "fetch_selection_ids",
                        lambda s, eid, nums, v=False: sel_ids_all)
    monkeypatch.setattr(legacy, "_market_num", lambda sel: "M")

    fail_sels = {f"O{fail_index}", f"U{fail_index}"} if fail_index is not None else set()

    def fake_calculate(session, sp, to, verbose=False):
        if concurrency_tracker is not None:
            with concurrency_tracker["lock"]:
                concurrency_tracker["live"] += 1
                concurrency_tracker["max"] = max(
                    concurrency_tracker["max"], concurrency_tracker["live"])
        try:
            if call_sleep:
                time.sleep(call_sleep)
            if to in fail_sels:
                raise RuntimeError("boom")
            return {"trueOdds": 3.5, "displayOdds": "+250"}
        finally:
            if concurrency_tracker is not None:
                with concurrency_tracker["lock"]:
                    concurrency_tracker["live"] -= 1

    monkeypatch.setattr(legacy, "calculate_sgp", fake_calculate)
    return targets


def test_parallel_pricing_returns_full_row_set(monkeypatch):
    targets = _wire_fixture(monkeypatch, n_targets=12)
    rows = draftkings.price_sgps(targets, periods=("FG",),
                                 client=MagicMock(), parallelism=8)
    assert len(rows) == 12 * 4
    keys = {(r.combo, r.spread_line, r.total_line) for r in rows}
    assert ("Home Spread + Over", -1.5, 7.5) in keys
    assert ("Away Spread + Under", -12.5, 18.5) in keys
    assert all(r.sgp_decimal == 3.5 and r.sgp_american == 250 for r in rows)
    assert all(r.source == "draftkings_direct" for r in rows)
    assert all(r.bookmaker == "draftkings" for r in rows)


def test_one_target_failure_does_not_drop_others(monkeypatch):
    targets = _wire_fixture(monkeypatch, n_targets=6, fail_index=3)
    rows = draftkings.price_sgps(targets, periods=("FG",),
                                 client=MagicMock(), parallelism=4)
    # target 3's four combos vanish; the other 5 targets price fully
    assert len(rows) == 5 * 4
    assert not any(r.total_line == 7.5 + 3 for r in rows)


def test_targets_actually_run_concurrently(monkeypatch):
    tracker = {"lock": threading.Lock(), "live": 0, "max": 0}
    targets = _wire_fixture(monkeypatch, n_targets=8, call_sleep=0.02,
                            concurrency_tracker=tracker)
    rows = draftkings.price_sgps(targets, periods=("FG",),
                                 client=MagicMock(), parallelism=8)
    assert len(rows) == 8 * 4
    # 8 targets x 4 combos with 20ms calls: sequential targets would
    # never exceed 4 in flight. Parallel must.
    assert tracker["max"] > 4


def test_parallelism_default_reads_env(monkeypatch):
    monkeypatch.setenv("MLB_SGP_DK_PARALLELISM", "3")
    assert draftkings._resolve_parallelism(None) == 3
    assert draftkings._resolve_parallelism(7) == 7
    monkeypatch.delenv("MLB_SGP_DK_PARALLELISM")
    assert draftkings._resolve_parallelism(None) == 8  # shipped default
