"""Target-level parallelism tests for mlb_sgp/fanduel.py::price_sgps."""
import sys
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import MagicMock

MLB_SGP_DIR = Path(__file__).resolve().parents[1]
if str(MLB_SGP_DIR) not in sys.path:
    sys.path.insert(0, str(MLB_SGP_DIR))

import scraper_fanduel_sgp as legacy  # noqa: E402

from mlb_sgp import fanduel  # noqa: E402
from mlb_sgp._shared import TargetLine  # noqa: E402

CT = datetime(2026, 6, 8, 23, 0, tzinfo=timezone.utc)


def _mk_target(game_id, spread, total):
    return TargetLine(
        game_id=game_id, home_team="New York Yankees",
        away_team="Boston Red Sox", commence_time=CT,
        period="FG", spread=spread, total=total,
    )


def _wire_fixture(monkeypatch, n_targets, fail_index=None):
    spreads, totals, targets = {}, {}, []
    for i in range(n_targets):
        spread, total = -1.5 - i, 7.5 + i
        targets.append(_mk_target("g1", spread, total))
        spreads[("home", spread)] = (f"SM{i}", f"SH{i}")
        spreads[("away", -spread)] = (f"SM{i}", f"SA{i}")
        totals[("O", total)] = (f"TM{i}", f"TO{i}")
        totals[("U", total)] = (f"TM{i}", f"TU{i}")
    runners = {"fg": {"spreads": spreads, "totals": totals},
               "f5": {"spreads": {}, "totals": {}}}
    monkeypatch.setattr(legacy, "fetch_fd_events", lambda s: [{"id": "e1"}])
    monkeypatch.setattr(legacy, "match_events", lambda evs, tgt: [
        {"game_id": "g1", "fd_event_id": "e1",
         "fd_home": "New York Yankees", "fd_away": "Boston Red Sox"}])
    monkeypatch.setattr(legacy, "fetch_event_runners",
                        lambda s, eid, h, a: runners)

    fail_sids = {f"TO{fail_index}", f"TU{fail_index}"} if fail_index is not None else set()

    def fake_price_combo(session, smid, ssid, tmid, tsid, verbose=False):
        if tsid in fail_sids:
            raise RuntimeError("boom")
        return {"decimal": 3.4, "american": 240}

    monkeypatch.setattr(legacy, "price_combo", fake_price_combo)
    return targets


def test_parallel_pricing_returns_full_row_set(monkeypatch):
    targets = _wire_fixture(monkeypatch, n_targets=10)
    rows = fanduel.price_sgps(targets, periods=("FG",),
                              client=MagicMock(), parallelism=6)
    assert len(rows) == 10 * 4
    assert all(r.sgp_decimal == 3.4 and r.sgp_american == 240 for r in rows)
    assert all(r.source == "fanduel_direct" for r in rows)


def test_one_target_failure_does_not_drop_others(monkeypatch):
    targets = _wire_fixture(monkeypatch, n_targets=5, fail_index=2)
    rows = fanduel.price_sgps(targets, periods=("FG",),
                              client=MagicMock(), parallelism=4)
    assert len(rows) == 4 * 4
    assert not any(r.total_line == 7.5 + 2 for r in rows)


def test_parallelism_default_reads_env(monkeypatch):
    monkeypatch.setenv("MLB_SGP_FD_PARALLELISM", "2")
    assert fanduel._resolve_parallelism(None) == 2
    assert fanduel._resolve_parallelism(9) == 9
    monkeypatch.delenv("MLB_SGP_FD_PARALLELISM")
    assert fanduel._resolve_parallelism(None) == 4  # conservative default


def test_fetcher_hooks_override_structure_fetches(monkeypatch):
    targets = _wire_fixture(monkeypatch, n_targets=2)
    monkeypatch.setattr(legacy, "fetch_fd_events",
                        lambda s: (_ for _ in ()).throw(AssertionError("hook bypassed")))
    monkeypatch.setattr(legacy, "fetch_event_runners",
                        lambda s, e, h, a: (_ for _ in ()).throw(AssertionError("hook bypassed")))

    spreads = {("home", -1.5): ("SM0", "SH0"), ("away", 1.5): ("SM0", "SA0"),
               ("home", -2.5): ("SM1", "SH1"), ("away", 2.5): ("SM1", "SA1")}
    totals = {("O", 7.5): ("TM0", "TO0"), ("U", 7.5): ("TM0", "TU0"),
              ("O", 8.5): ("TM1", "TO1"), ("U", 8.5): ("TM1", "TU1")}
    runners = {"fg": {"spreads": spreads, "totals": totals},
               "f5": {"spreads": {}, "totals": {}}}
    fetchers = {
        "fetch_fd_events": lambda session: [{"id": "e1"}],
        "fetch_event_runners": lambda session, eid, h, a: runners,
    }
    rows = fanduel.price_sgps(targets, periods=("FG",), client=MagicMock(),
                              parallelism=2, fetchers=fetchers)
    assert len(rows) == 2 * 4
