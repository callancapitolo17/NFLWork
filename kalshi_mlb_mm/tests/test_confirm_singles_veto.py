"""Issue #17 (reworked): confirm-window singles-move veto.

The last-look gate previously re-priced grid combos from the same ~150s
`_SGP_ODDS` scrape cache the quote was priced from, so cur_fair == prev_fair
by construction and the drift check was a no-op exactly in the fast-pickoff
window. The rework: snapshot the RAW Kalshi odds of every leg at quote time
(each leg is its own Kalshi singles market), re-read them at accept (<~1s),
and void if ANY leg's bid or ask moved even one tick. Correlation between
legs is structural on this horizon — the marginals carry the pickoff signal.

These tests seed a HEALTHY stale cache that prices the combo fine (so the
pre-rework code would confirm) and drive the decision through the veto.
Scaffolding follows test_main_on_demand.py.
"""
import importlib
import json
from datetime import datetime, timezone

import pandas as pd

from kalshi_common import legset
from kalshi_common.leg_types import SPREAD_TOTAL_FAMILY
from mlb_sgp._shared import GameRef

EVT = "KXMLBGAME-25JUN271905TEXLAA"
SPREAD_TICKER = "KXMLBSPREAD-25JUN271905TEXLAA-LAA2"
TOTAL_TICKER = "KXMLBTOTAL-25JUN271905TEXLAA-9"
GRID_LEGS = [  # 2-leg same-game spread+total -> grid_spread_total route
    {"market_ticker": SPREAD_TICKER, "event_ticker": EVT, "side": "yes"},
    {"market_ticker": TOTAL_TICKER, "event_ticker": EVT, "side": "yes"},
]
GREF = GameRef(game_id="game1", home_team="Los Angeles Angels",
               away_team="Texas Rangers", commence_time=None)

# Quote-time baseline: raw Kalshi odds per leg, dollars.
SNAPSHOT = {SPREAD_TICKER: {"yes_bid": 0.45, "yes_ask": 0.47},
            TOTAL_TICKER: {"yes_bid": 0.52, "yes_ask": 0.54}}

# 4 equal decimals -> probit devig gives exactly 0.25 per cell, per book —
# a healthy cache fair the PRE-rework last look would happily confirm on.
STALE_CACHE_FAIR = 0.25


def _stale_grid_df():
    rows = []
    for book in ("draftkings", "fanduel"):
        for cell in SPREAD_TOTAL_FAMILY:
            rows.append(dict(game_id="game1", combo=cell, period="FG",
                             bookmaker=book, sgp_decimal=3.8, fetch_time=None,
                             spread_line=-1.5, total_line=8.5))
    return pd.DataFrame(rows)


def _setup(monkeypatch, tmp_path, db_name, *, snapshot_json=None,
           fresh=None, yes_bid=0.28, no_bid=0.68):
    """Open GRID quote seeded in live_quotes + accepted-quote API stub.
    `snapshot_json` is stored on the quote row; `fresh` is what the confirm
    tick's live leg read returns (None models a failed read)."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    from kalshi_common import auth_client
    from kalshi_mlb_mm import main
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / db_name)
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    importlib.reload(db)
    db.init_database()
    monkeypatch.setattr(main, "_SGP_ODDS", _stale_grid_df())
    monkeypatch.setattr(main, "_resolve_game_for_legs", lambda gl: "game1")
    monkeypatch.setattr(main, "_game_ref", lambda gid: GREF)
    monkeypatch.setattr(main, "_ENGINE", None)
    monkeypatch.setattr(main, "_leg_market_prices", lambda legs: fresh)
    now = datetime.now(timezone.utc)
    with db.connect() as con:
        con.execute(
            "INSERT INTO live_quotes (quote_id, rfq_id, combo_market_ticker, "
            "game_id, yes_bid, no_bid, model_fair, book_fair, blended_fair, "
            "status, submitted_at, closed_at, leg_prices_json) "
            "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
            ["q-grid", "r-grid", "COMBO-GRID", "game1", yes_bid, no_bid,
             None, STALE_CACHE_FAIR, STALE_CACHE_FAIR, "open", now,
             None, snapshot_json])
        con.execute(
            "INSERT OR REPLACE INTO seen_rfqs (rfq_id, market_ticker, in_scope, "
            "game_id, legs_json, first_seen_at, last_decision, creator_id) "
            "VALUES (?,?,?,?,?,?,?,?)",
            ["r-grid", "COMBO-GRID", True, "game1", json.dumps(GRID_LEGS),
             now, "quoted", ""])
    monkeypatch.setattr(auth_client, "api",
                        lambda *a, **k: (200, {"quote": {"status": "accepted",
                                                         "accepted_side": "yes",
                                                         "contracts": 1}}, None))
    return main, db


def _quote_state(db):
    with db.connect(read_only=True) as con:
        st = con.execute(
            "SELECT status FROM live_quotes WHERE quote_id='q-grid'").fetchone()[0]
        dec = con.execute(
            "SELECT decision FROM quote_decisions WHERE quote_id='q-grid' "
            "ORDER BY observed_at DESC LIMIT 1").fetchone()[0]
        n_fills = con.execute("SELECT COUNT(*) FROM fills").fetchone()[0]
    return st, dec, n_fills


class MustNotConfirmGW:
    def confirm(self, qid):
        raise AssertionError("accept must void, not confirm")

    def cancel(self, qid):
        return True


class ConfirmGW:
    def __init__(self):
        self.confirmed = []

    def confirm(self, qid):
        self.confirmed.append(qid)
        return True

    def cancel(self, qid):
        return True


def _moved_fresh():
    fresh = json.loads(json.dumps(SNAPSHOT))
    fresh[SPREAD_TICKER]["yes_bid"] = 0.44        # one tick down on one leg
    return fresh


def test_confirm_voids_when_any_leg_ticked(monkeypatch, tmp_path):
    """THE regression: one leg's bid moved one tick while the combo cache
    still prices fine (pre-rework code confirmed here) -> void."""
    main, db = _setup(monkeypatch, tmp_path, "veto_moved.duckdb",
                      snapshot_json=json.dumps(SNAPSHOT), fresh=_moved_fresh())
    main._confirm_tick(MustNotConfirmGW(), dry_run=False)
    st, dec, n_fills = _quote_state(db)
    assert st == "voided"
    assert dec == "voided_singles_moved"
    assert n_fills == 0


def test_confirm_voids_when_fresh_read_fails(monkeypatch, tmp_path):
    main, db = _setup(monkeypatch, tmp_path, "veto_nofetch.duckdb",
                      snapshot_json=json.dumps(SNAPSHOT), fresh=None)
    main._confirm_tick(MustNotConfirmGW(), dry_run=False)
    st, dec, _ = _quote_state(db)
    assert st == "voided"
    assert dec == "voided_singles_moved"


def test_confirm_voids_when_snapshot_missing(monkeypatch, tmp_path):
    """Pre-migration row (NULL snapshot) fails CLOSED — no baseline, no fill."""
    main, db = _setup(monkeypatch, tmp_path, "veto_nosnap.duckdb",
                      snapshot_json=None, fresh=dict(SNAPSHOT))
    main._confirm_tick(MustNotConfirmGW(), dry_run=False)
    st, dec, _ = _quote_state(db)
    assert st == "voided"
    assert dec == "voided_singles_moved"


def test_confirm_passes_when_all_legs_unchanged(monkeypatch, tmp_path):
    """Unmoved raw leg odds -> veto passes -> existing cache gates run ->
    confirm; fair_at_confirm keeps the cache fair convention."""
    main, db = _setup(monkeypatch, tmp_path, "veto_pass.duckdb",
                      snapshot_json=json.dumps(SNAPSHOT),
                      fresh=json.loads(json.dumps(SNAPSHOT)))
    gw = ConfirmGW()
    main._confirm_tick(gw, dry_run=False)
    assert gw.confirmed == ["q-grid"]
    with db.connect(read_only=True) as con:
        fair_at_confirm = con.execute(
            "SELECT fair_at_confirm FROM fills WHERE quote_id='q-grid'").fetchone()[0]
    assert abs(fair_at_confirm - STALE_CACHE_FAIR) < 1e-6


def test_confirm_singles_check_emitted_both_outcomes(monkeypatch, tmp_path):
    """Firehose measurement: one event per accept, carrying quote age and the
    per-leg snapshot-vs-fresh odds (the tolerance-tuning dataset)."""
    main, db = _setup(monkeypatch, tmp_path, "veto_evt1.duckdb",
                      snapshot_json=json.dumps(SNAPSHOT), fresh=_moved_fresh())
    emitted = []
    monkeypatch.setattr(main.research, "emit",
                        lambda ev, **kw: emitted.append((ev, kw)))
    main._confirm_tick(MustNotConfirmGW(), dry_run=False)
    checks = [kw for ev, kw in emitted if ev == "confirm_singles_check"]
    assert len(checks) == 1
    payload = checks[0]["payload"]
    assert payload["moved"] is True
    assert payload["quote_age_sec"] is not None and payload["quote_age_sec"] >= 0
    assert payload["snapshot"][SPREAD_TICKER]["yes_bid"] == 0.45
    assert payload["fresh"][SPREAD_TICKER]["yes_bid"] == 0.44

    main2, db2 = _setup(monkeypatch, tmp_path, "veto_evt2.duckdb",
                        snapshot_json=json.dumps(SNAPSHOT),
                        fresh=json.loads(json.dumps(SNAPSHOT)))
    emitted2 = []
    monkeypatch.setattr(main2.research, "emit",
                        lambda ev, **kw: emitted2.append((ev, kw)))
    main2._confirm_tick(ConfirmGW(), dry_run=False)
    checks2 = [kw for ev, kw in emitted2 if ev == "confirm_singles_check"]
    assert len(checks2) == 1 and checks2[0]["payload"]["moved"] is False


# --------------------------------------------------------------------------- #
# Helper-level units                                                          #
# --------------------------------------------------------------------------- #

def test_singles_moved_detects_any_tick_and_leg_mismatch():
    from kalshi_mlb_mm import main
    same = json.loads(json.dumps(SNAPSHOT))
    assert main._singles_moved(SNAPSHOT, same) is False
    assert main._singles_moved(SNAPSHOT, _moved_fresh()) is True
    ask_moved = json.loads(json.dumps(SNAPSHOT))
    ask_moved[TOTAL_TICKER]["yes_ask"] = 0.55
    assert main._singles_moved(SNAPSHOT, ask_moved) is True
    missing_leg = {SPREAD_TICKER: dict(SNAPSHOT[SPREAD_TICKER])}
    assert main._singles_moved(SNAPSHOT, missing_leg) is True


def test_market_bid_ask_handles_both_kalshi_shapes():
    from kalshi_mlb_mm import main
    assert main._market_bid_ask({"yes_bid": 45, "yes_ask": 47}) == (0.45, 0.47)
    assert main._market_bid_ask(
        {"yes_bid_dollars": "0.45", "yes_ask_dollars": "0.47"}) == (0.45, 0.47)
    assert main._market_bid_ask({"yes_bid": 45}) is None
    assert main._market_bid_ask({}) is None


def test_leg_market_prices_reads_one_get_per_leg(monkeypatch):
    from kalshi_common import auth_client
    from kalshi_mlb_mm import main
    calls = []

    def fake_api(method, path, *a, **k):
        calls.append(path)
        t = path.rsplit("/", 1)[-1]
        return 200, {"market": {"yes_bid": 45 if t == SPREAD_TICKER else 52,
                                "yes_ask": 47 if t == SPREAD_TICKER else 54}}, None
    monkeypatch.setattr(auth_client, "api", fake_api)
    out = main._leg_market_prices(GRID_LEGS)
    assert out == SNAPSHOT
    assert calls == [f"/markets/{SPREAD_TICKER}", f"/markets/{TOTAL_TICKER}"]


def test_leg_market_prices_fails_closed_on_any_leg(monkeypatch):
    from kalshi_common import auth_client
    from kalshi_mlb_mm import main

    def fake_api(method, path, *a, **k):
        if SPREAD_TICKER in path:
            return 200, {"market": {"yes_bid": 45, "yes_ask": 47}}, None
        return 500, "boom", None                   # second leg unreadable
    monkeypatch.setattr(auth_client, "api", fake_api)
    assert main._leg_market_prices(GRID_LEGS) is None
