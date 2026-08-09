"""Issue #81 — target-line refresh extracted from the sweep.

`target_line_cycle` is the sweep-free half of the old `sgp_cycle`: Kalshi
MVE enumeration + Odds API schedule -> `mlb_target_lines`. It must be
callable standalone with ZERO book requests (no service, no scraper params
— the signature itself is the proof), because after #81 it is the maker's
only background writer of the market DB besides fetch-health. `sgp_cycle`
(still used by the taker's sweep) delegates to it, so taker behavior is
pinned unchanged.

Written FIRST, failing on pre-#81 code, per repo rule.
"""
import inspect
from datetime import datetime, timezone

import duckdb

from kalshi_common import sgp_runner
from mlb_sgp._shared import TargetLine

TARGETS = [
    TargetLine(game_id="g1", home_team="Boston Red Sox",
               away_team="New York Yankees",
               commence_time=datetime(2026, 8, 10, 23, 5, tzinfo=timezone.utc),
               period="FG", spread=-1.5, total=8.5),
    TargetLine(game_id="g2", home_team="Chicago Cubs",
               away_team="St. Louis Cardinals",
               commence_time=datetime(2026, 8, 10, 23, 40, tzinfo=timezone.utc),
               period="FG", spread=1.5, total=9.0),
]


def test_signature_is_keyword_only_and_service_free():
    """Zero book requests by construction: the function accepts no service,
    no scraper knobs — only the DB path and the enumeration flag, both
    keyword-only (repo correctness bar)."""
    sig = inspect.signature(sgp_runner.target_line_cycle)
    assert set(sig.parameters) == {"bot_market_db", "both_teams"}
    for param in sig.parameters.values():
        assert param.kind == inspect.Parameter.KEYWORD_ONLY, (
            f"{param.name} must be keyword-only")


def test_enumerates_and_writes_target_lines(tmp_path, monkeypatch):
    seen_flags = []

    def fake_enumerate(both_teams=False):
        seen_flags.append(both_teams)
        return list(TARGETS)

    monkeypatch.setattr(sgp_runner, "enumerate_kalshi_targets", fake_enumerate)
    db_path = str(tmp_path / "market.duckdb")

    out = sgp_runner.target_line_cycle(bot_market_db=db_path, both_teams=True)

    assert seen_flags == [True]
    assert out == list(TARGETS)      # callers (warmup log, tests) see the list
    con = duckdb.connect(db_path, read_only=True)
    try:
        n, games = con.execute(
            "SELECT count(*), count(DISTINCT game_id) FROM mlb_target_lines"
        ).fetchone()
    finally:
        con.close()
    assert n == len(TARGETS)
    assert games == 2


def test_sgp_cycle_service_path_delegates_to_target_line_cycle(
        tmp_path, monkeypatch):
    """The taker's sweep is out of scope (#81 non-goal): `sgp_cycle` must
    keep writing target lines — now via the extracted half — and still run
    the book refresh."""
    monkeypatch.setattr(sgp_runner, "enumerate_kalshi_targets",
                        lambda both_teams=False: list(TARGETS))
    cycle_calls = []
    real_cycle = sgp_runner.target_line_cycle

    def spy_cycle(**kw):
        cycle_calls.append(kw)
        return real_cycle(**kw)

    monkeypatch.setattr(sgp_runner, "target_line_cycle", spy_cycle)

    class NoBookService:
        def refresh(self, targets):
            # every book "failed" -> no upserts, no importlib of book modules
            return {"draftkings": None}

    db_path = str(tmp_path / "market.duckdb")
    counts = sgp_runner.sgp_cycle(bot_market_db=db_path,
                                  service=NoBookService(), both_teams=True)

    assert counts == {"draftkings": -1}
    assert cycle_calls == [{"bot_market_db": db_path, "both_teams": True}]
    con = duckdb.connect(db_path, read_only=True)
    try:
        n = con.execute("SELECT count(*) FROM mlb_target_lines").fetchone()[0]
    finally:
        con.close()
    assert n == len(TARGETS)
