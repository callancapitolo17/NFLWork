"""Issue #54 — QUOTE_PRICING_MODE routing: cached / shadow / live.

PIN TESTS first (written against pre-#54 code, must stay green forever):
`cached` mode is byte-identical to the pre-#54 quote path — a grid RFQ is
priced from the sweep cache with ZERO on-demand engine traffic.

Regression tests for the live/shadow modes follow; they FAIL on pre-#54
code by construction (grids never touched the engine before this ticket).

Scaffolding style follows test_main_on_demand.py; leg dicts go through
conftest.leg() so the Kalshi event-ticker shape is real (#71).
"""
import importlib

import pandas as pd
import pytest

from kalshi_common import legset
from mlb_sgp._shared import GameRef

from kalshi_mlb_mm.tests.conftest import leg

# 2-leg same-game spread+total -> grid_spread_total route. LAA is the home
# team of 25JUN271905TEXLAA, so LAA2 is home -1.5 (target "Home Spread + Over").
GRID_LEGS = [
    leg("KXMLBSPREAD-25JUN271905TEXLAA-LAA2", "yes"),
    leg("KXMLBTOTAL-25JUN271905TEXLAA-9", "yes"),
]
GREF = GameRef(game_id="game1", home_team="Los Angeles Angels",
               away_team="Texas Rangers", commence_time=None)

# Balanced 4-cell grid: every cell devigs to ~0.25, so the CACHED consensus
# fair is ~0.2375 — far from the 0.30/0.32 LIVE fairs the FakeEngine serves,
# which lets every assertion prove WHICH source priced the quote.
ST_CELLS = {"Home Spread + Over": 4.2, "Home Spread + Under": 4.2,
            "Away Spread + Over": 3.8, "Away Spread + Under": 3.8}
CACHED_FAIR = pytest.approx(0.2375, abs=0.01)
LIVE_FAIRS = {"draftkings": 0.30, "fanduel": 0.32}
LIVE_FAIR = pytest.approx(0.31, abs=0.001)


def _grid_rows(book, game_id, spread_line, total_line, family_decimals):
    return [
        {"game_id": game_id, "combo": c, "period": "FG", "bookmaker": book,
         "sgp_decimal": d, "fetch_time": None,
         "spread_line": spread_line, "total_line": total_line}
        for c, d in family_decimals.items()
    ]


class FakeEngine:
    """Duck-typed OnDemandEngine: serves canned live fairs, records traffic."""

    def __init__(self, fairs=None, empty_landed=False):
        self.fairs = fairs                    # {book: fair} | None
        self.empty_landed = empty_landed      # fresh landing with zero books
        self.ensure_calls = []
        self.refetch_calls = []

    def lookup(self, h):
        return self.fairs

    def lookup_results(self, h):
        if not self.fairs:
            return None
        from mlb_sgp._shared import OnDemandBookResult
        return {b: OnDemandBookResult(book=b, fair=f, route="partition",
                                      n_cells_priced=4, latency_sec=1.2)
                for b, f in self.fairs.items()}

    def landed_at(self, h):
        return 1000.0 if (self.fairs or self.empty_landed) else None

    def landed_empty(self, h):
        return self.empty_landed

    def result_age_sec(self, h):
        return 2.0 if (self.fairs or self.empty_landed) else None

    def ensure_fetch(self, h, game, legs):
        self.ensure_calls.append((h, game, tuple(legs)))
        return True

    def refetch_now(self, jobs, deadline_sec):
        self.refetch_calls.append((list(jobs), deadline_sec))
        return self.fairs is not None


class Src:
    def poll(self):
        return [{"id": "r-grid", "market_ticker": "COMBO-GRID", "contracts": 1}]

    def get_market(self, t):
        return {}


class NoQuoteGW:
    def submit_quote(self, *a):
        raise AssertionError("must not submit in these tests")


def _setup(monkeypatch, tmp_path, engine, db_name, mode):
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    from kalshi_mlb_mm import main
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / db_name)
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    # raising=False: the attribute does not exist on pre-#54 config, and the
    # pin tests must run (and pass) against that code.
    monkeypatch.setattr(cfg, "QUOTE_PRICING_MODE", mode, raising=False)
    importlib.reload(db)
    db.init_database()
    monkeypatch.setattr(main, "_SGP_ODDS", pd.DataFrame(
        _grid_rows("dk", "game1", -1.5, 8.5, ST_CELLS)
        + _grid_rows("fd", "game1", -1.5, 8.5, ST_CELLS)))
    monkeypatch.setattr(main, "_today_fills", lambda: [])
    monkeypatch.setattr(main, "_today_fills_by_game", lambda: [])
    monkeypatch.setattr(main, "_resolve_game_for_legs", lambda gl: "game1")
    monkeypatch.setattr(main, "_commence_time", lambda gid: None)
    monkeypatch.setattr(main, "_game_ref", lambda gid: GREF)
    monkeypatch.setattr(risk, "tipoff_ok", lambda ct, min_: True)
    monkeypatch.setattr(main, "_ENGINE", engine)
    monkeypatch.setattr(main, "_SCOPE_CACHE",
                        {"COMBO-GRID": (True, None, GRID_LEGS)})
    monkeypatch.setattr(main, "_OD_RESULT_EMITTED", {})
    monkeypatch.setattr(main, "_PREV_BOOK_FAIR", {})
    emitted = []
    monkeypatch.setattr(main.research, "emit",
                        lambda ev, **kw: emitted.append((ev, kw)))
    return main, db, emitted


def _last_decision(db):
    with db.connect(read_only=True) as con:
        return con.execute(
            "SELECT decision, reason FROM quote_decisions "
            "ORDER BY observed_at DESC LIMIT 1").fetchone()


def _quote_priced_payload(emitted):
    payloads = [kw["payload"] for ev, kw in emitted if ev == "quote_priced"]
    return payloads[-1] if payloads else None


# --------------------------------------------------------------------------- #
# PIN — cached mode is the pre-#54 quote path, byte-identical                  #
# --------------------------------------------------------------------------- #

def test_cached_mode_grid_prices_from_cache_no_engine_traffic(monkeypatch, tmp_path):
    eng = FakeEngine(fairs=LIVE_FAIRS)     # engine HAS fairs — must be ignored
    main, db, emitted = _setup(monkeypatch, tmp_path, eng, "pin1.duckdb",
                               mode="cached")
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    assert _last_decision(db) == ("dry_run_quote", None)
    assert eng.ensure_calls == [], "cached mode must never touch the engine"
    assert eng.refetch_calls == []
    payload = _quote_priced_payload(emitted)
    assert payload is not None
    assert payload["blended_fair"] == CACHED_FAIR, \
        "cached mode must price from the sweep grid, not the live fairs"


def test_cached_mode_grid_quotes_even_when_engine_has_nothing(monkeypatch, tmp_path):
    # The stronger pin: with an engine that has NO results at all, cached mode
    # still quotes (pre-#54 grids never depended on the engine).
    eng = FakeEngine(fairs=None)
    main, db, emitted = _setup(monkeypatch, tmp_path, eng, "pin2.duckdb",
                               mode="cached")
    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    assert _last_decision(db) == ("dry_run_quote", None)
    assert eng.ensure_calls == []
    assert _quote_priced_payload(emitted)["blended_fair"] == CACHED_FAIR
