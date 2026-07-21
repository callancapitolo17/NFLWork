"""Phase 2 main.py integration: discovery feed, confirm live last look,
out-of-scope instrumentation. Scaffolding style follows test_main_smoke.py."""
import importlib
import json
from datetime import datetime, timezone

import pandas as pd

from kalshi_common import legset
from mlb_sgp._shared import GameRef

EVT = "KXMLBGAME-25JUN271905TEXLAA"
OD_LEGS = [  # 3-leg same-game -> on_demand route
    {"market_ticker": "KXMLBSPREAD-25JUN271905TEXLAA-LAA2",
     "event_ticker": EVT, "side": "yes"},
    {"market_ticker": "KXMLBTOTAL-25JUN271905TEXLAA-9",
     "event_ticker": EVT, "side": "yes"},
    {"market_ticker": "KXMLBGAME-25JUN271905TEXLAA-LAA",
     "event_ticker": EVT, "side": "yes"},
]
GREF = GameRef(game_id="game1", home_team="Los Angeles Angels",
               away_team="Texas Rangers", commence_time=None)


class FakeEngine:
    def __init__(self, fairs=None, landed=1000.0):
        self.fairs = fairs            # dict | None served by lookup
        self.landed = landed
        self.ensure_calls = []
        self.refetch_calls = []

    def lookup(self, h):
        return self.fairs

    def lookup_results(self, h):
        if not self.fairs:
            return None
        from mlb_sgp._shared import OnDemandBookResult
        return {b: OnDemandBookResult(book=b, fair=f, route="partition",
                                      n_cells_priced=8, latency_sec=1.0)
                for b, f in self.fairs.items()}

    def landed_at(self, h):
        return self.landed if self.fairs else None

    def ensure_fetch(self, h, game, legs):
        self.ensure_calls.append((h, game, tuple(legs)))
        return len(self.ensure_calls) == 1     # first call = newly enqueued

    def refetch_now(self, jobs, deadline_sec):
        self.refetch_calls.append((list(jobs), deadline_sec))
        return self.fairs is not None


def _setup(monkeypatch, tmp_path, engine, db_name):
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    from kalshi_mlb_mm import main
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / db_name)
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    importlib.reload(db)
    db.init_database()
    monkeypatch.setattr(main, "_SGP_ODDS",
                        pd.DataFrame({"game_id": ["g"], "combo": ["c"], "period": ["FG"],
                                      "bookmaker": ["dk"], "sgp_decimal": [2.0],
                                      "fetch_time": [None], "spread_line": [-1.5],
                                      "total_line": [8.5]}))
    monkeypatch.setattr(main, "_today_fills", lambda: [])
    monkeypatch.setattr(main, "_today_fills_by_game", lambda: [])
    monkeypatch.setattr(main, "_resolve_game_for_legs", lambda gl: "game1")
    monkeypatch.setattr(main, "_commence_time", lambda gid: None)
    monkeypatch.setattr(main, "_game_ref", lambda gid: GREF)
    monkeypatch.setattr(risk, "tipoff_ok", lambda ct, min_: True)
    monkeypatch.setattr(main, "_ENGINE", engine)
    monkeypatch.setattr(main, "_SCOPE_CACHE",
                        {"COMBO-OD": (True, None, OD_LEGS)})
    monkeypatch.setattr(main, "_OD_RESULT_EMITTED", {})
    return main, db


class SrcOD:
    def poll(self):
        return [{"id": "r-od", "market_ticker": "COMBO-OD", "contracts": 1}]

    def get_market(self, t):
        return {}


class NoQuoteGW:
    def submit_quote(self, *a):
        raise AssertionError("must not quote")


def test_discovery_on_demand_pending_enqueues_and_skips(monkeypatch, tmp_path):
    eng = FakeEngine(fairs=None)
    main, db = _setup(monkeypatch, tmp_path, eng, "od1.duckdb")
    emitted = []
    monkeypatch.setattr(main.research, "emit",
                        lambda ev, **kw: emitted.append((ev, kw)))
    main._discovery_tick(SrcOD(), NoQuoteGW(), dry_run=True)
    canon = legset.parse_legs(OD_LEGS)
    h = legset.leg_set_hash(canon)
    assert eng.ensure_calls and eng.ensure_calls[0][0] == h
    assert eng.ensure_calls[0][1] == GREF
    with db.connect(read_only=True) as con:
        row = con.execute(
            "SELECT decision, reason FROM quote_decisions "
            "ORDER BY observed_at DESC LIMIT 1").fetchone()
    assert row == ("skipped", "on_demand_pending")
    assert [e for e, _ in emitted].count("on_demand_requested") == 1
    # second tick: still pending, ensure_fetch deduped by engine -> no 2nd emit
    main._discovery_tick(SrcOD(), NoQuoteGW(), dry_run=True)
    assert [e for e, _ in emitted].count("on_demand_requested") == 1


def test_discovery_fresh_result_prices_and_emits_result(monkeypatch, tmp_path):
    eng = FakeEngine(fairs={"draftkings": 0.20, "fanduel": 0.21})
    main, db = _setup(monkeypatch, tmp_path, eng, "od2.duckdb")
    emitted = []
    monkeypatch.setattr(main.research, "emit",
                        lambda ev, **kw: emitted.append((ev, kw)))
    main._discovery_tick(SrcOD(), NoQuoteGW(), dry_run=True)
    assert eng.ensure_calls == []                       # fresh -> no fetch
    with db.connect(read_only=True) as con:
        rows = con.execute("SELECT decision FROM quote_decisions").fetchall()
    assert ("dry_run_quote",) in rows                   # priced via lookup
    assert [e for e, _ in emitted].count("on_demand_result") == 1
    # same landing consumed again -> no duplicate result event
    main._discovery_tick(SrcOD(), NoQuoteGW(), dry_run=True)
    assert [e for e, _ in emitted].count("on_demand_result") == 1


def test_out_of_scope_instrumentation_stores_legs_and_reason(monkeypatch, tmp_path):
    main, db = _setup(monkeypatch, tmp_path, FakeEngine(), "od3.duckdb")
    non_mlb = [{"market_ticker": "KXNBAGAME-FOO-BAR", "event_ticker": "E", "side": "yes"},
               {"market_ticker": "KXNBAGAME-FOO-BAZ", "event_ticker": "E", "side": "yes"}]
    lone = [OD_LEGS[0]]
    monkeypatch.setattr(main, "_SCOPE_CACHE", {
        "COMBO-NBA": (False, None, non_mlb),
        "COMBO-LONE": (False, None, lone),
        "COMBO-DEAD": (False, None, None),
    })

    class Src:
        def poll(self):
            return [{"id": f"r-{t}", "market_ticker": t, "contracts": 1}
                    for t in ("COMBO-NBA", "COMBO-LONE", "COMBO-DEAD")]

        def get_market(self, t):
            return {}

    main._discovery_tick(Src(), NoQuoteGW(), dry_run=True)
    with db.connect(read_only=True) as con:
        rows = dict(con.execute(
            "SELECT market_ticker, last_decision FROM seen_rfqs").fetchall())
        legs_json = dict(con.execute(
            "SELECT market_ticker, legs_json FROM seen_rfqs").fetchall())
    assert rows["COMBO-NBA"] == "out_of_scope_non_mlb"
    assert rows["COMBO-LONE"] == "out_of_scope_lone_single"
    assert rows["COMBO-DEAD"] == "out_of_scope_unparseable"
    assert json.loads(legs_json["COMBO-NBA"]) == non_mlb
    assert legs_json["COMBO-DEAD"] is None


def test_confirm_on_demand_refetches_and_voids_on_stale(monkeypatch, tmp_path):
    eng = FakeEngine(fairs=None)                        # stale -> void
    main, db = _setup(monkeypatch, tmp_path, eng, "od4.duckdb")
    canon = legset.parse_legs(OD_LEGS)
    h = legset.leg_set_hash(canon)
    now = datetime.now(timezone.utc)
    with db.connect() as con:
        con.execute("INSERT INTO live_quotes (quote_id, rfq_id, combo_market_ticker, game_id, yes_bid, no_bid, model_fair, book_fair, blended_fair, status, submitted_at, closed_at) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                    ["q-od", "r-od", "COMBO-OD", "game1", 0.5, 0.43,
                     0.2, 0.2, 0.2, "open", now, None])
        con.execute(
            "INSERT OR REPLACE INTO seen_rfqs (rfq_id, market_ticker, in_scope, "
            "game_id, legs_json, first_seen_at, last_decision, creator_id) "
            "VALUES (?,?,?,?,?,?,?,?)",
            ["r-od", "COMBO-OD", True, "game1", json.dumps(OD_LEGS), now, "quoted", ""])
    monkeypatch.setattr(main.auth_client, "api",
                        lambda *a, **k: (200, {"quote": {"status": "accepted",
                                                         "accepted_side": "yes",
                                                         "contracts": 1}}, None))

    class GW:
        def confirm(self, qid):
            raise AssertionError("stale on-demand fair must void, not confirm")

        def cancel(self, qid):
            return True

    main._confirm_tick(GW(), dry_run=False)
    assert eng.refetch_calls, "confirm must trigger the live re-fetch"
    jobs, budget = eng.refetch_calls[0]
    assert jobs[0][0] == h and jobs[0][1] == GREF
    assert budget == main.CONFIRM_REFETCH_BUDGET_SEC
    with db.connect(read_only=True) as con:
        st = con.execute("SELECT status FROM live_quotes WHERE quote_id='q-od'").fetchone()
        dec = con.execute("SELECT decision FROM quote_decisions "
                          "ORDER BY observed_at DESC LIMIT 1").fetchone()
    assert st[0] == "voided"
    assert dec[0] == "voided_no_fresh_books"
