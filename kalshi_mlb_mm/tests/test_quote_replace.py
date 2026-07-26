"""Issue #16 (audit B2) — the quote-replace flow must never orphan a live
exchange quote. Kalshi auto-cancels the prior quote on an RFQ only when a NEW
quote is successfully submitted, so the DB must not mark the old row 'replaced'
until the resubmit has actually succeeded."""
from datetime import datetime, timezone


_EVT = "KXMLBGAME-25JUN271905TEXLAA"
_LEGS = [{"market_ticker": f"KXMLBSPREAD-{_EVT.rsplit('-', 1)[-1]}-LAA2",
          "event_ticker": _EVT, "side": "yes"},
         {"market_ticker": f"KXMLBTOTAL-{_EVT.rsplit('-', 1)[-1]}-9",
          "event_ticker": _EVT, "side": "yes"}]


def _replace_env(monkeypatch, tmp_path, db_name):
    """Discovery-tick env where RFQ r1 already has an open quote at prices
    BEYOND hysteresis from what pricing now produces → replace path fires."""
    import importlib
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    import kalshi_mlb_mm.router as router_mod
    import kalshi_mlb_mm.pricing as pricing_mod
    from kalshi_mlb_mm.pricing import Quote
    from kalshi_mlb_mm import main
    import pandas as pd

    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / db_name)
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    monkeypatch.setattr(cfg, "MAX_COMBO_EXPOSURE_USD", 200.0)
    importlib.reload(db)
    db.init_database()

    monkeypatch.setattr(main, "_SGP_ODDS",
                        pd.DataFrame({"game_id": ["game1"], "combo": ["c"], "period": ["FG"],
                                      "bookmaker": ["dk"], "sgp_decimal": [2.0],
                                      "fetch_time": [None], "spread_line": [-1.5],
                                      "total_line": [8.5]}))
    monkeypatch.setattr(main, "_today_fills", lambda: [])
    monkeypatch.setattr(router_mod, "combo_fair_detail", lambda *a, **k: (router_mod.ComboFair(0.55, 0.0, 1), "ok"))
    monkeypatch.setattr(main, "_resolve_game_for_legs", lambda gl: "game1")
    monkeypatch.setattr(main, "_leg_market_prices",
                        lambda legs: {"L": {"yes_bid": 0.5, "yes_ask": 0.52}})
    monkeypatch.setattr(main, "_commence_time", lambda gid: None)
    monkeypatch.setattr(risk, "tipoff_ok", lambda ct, min_: True)
    monkeypatch.setattr(main, "_SCOPE_CACHE", {"COMBO-1": (True, "game1", _LEGS)})
    # New price differs from the pre-seeded quote by >> QUOTE_HYSTERESIS.
    monkeypatch.setattr(pricing_mod, "quote",
                        lambda fair, roi, **kw: Quote(yes_bid=0.520, no_bid=0.410))

    with db.connect() as con:
        con.execute(
            "INSERT INTO live_quotes (quote_id, rfq_id, combo_market_ticker, game_id, yes_bid, no_bid, model_fair, book_fair, blended_fair, status, submitted_at, closed_at) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
            ["qid-old", "r1", "COMBO-1", "game1", 0.500, 0.430,
             None, 0.55, 0.55, "open", datetime.now(timezone.utc), None])

    class Src:
        def poll(self):
            return [{"id": "r1", "market_ticker": "COMBO-1", "contracts": 1}]

        def get_market(self, t):
            return {}

    return db, main, Src()


def test_failed_resubmit_leaves_old_quote_tracked(monkeypatch, tmp_path):
    """THE regression test for issue #16: submit_quote fails (None) → the old
    row must stay in a tracked status ('open' or 'cancelled'), NEVER 'replaced'
    (pre-fix code marked it 'replaced' before submitting, orphaning the live
    exchange quote from the confirm loop and the risk sweep)."""
    db, main, src = _replace_env(monkeypatch, tmp_path, "replace_fail.duckdb")

    class FailingGW:
        def submit_quote(self, rid, yb, nb):
            return None            # 4xx/5xx/timeout

    main._discovery_tick(src, FailingGW(), dry_run=False)

    with db.connect(read_only=True) as con:
        st = con.execute(
            "SELECT status FROM live_quotes WHERE quote_id='qid-old'").fetchone()[0]
        n_rows = con.execute("SELECT COUNT(*) FROM live_quotes").fetchone()[0]
    assert st in ("open", "cancelled"), \
        f"old quote must stay tracked after failed resubmit, got '{st}'"
    assert n_rows == 1, "no new live_quotes row on failed submit"


def test_successful_replace_marks_old_replaced_and_inserts_new(monkeypatch, tmp_path):
    """Happy path: successful resubmit → exactly one 'open' row (the new
    quote) and the old row 'replaced'."""
    db, main, src = _replace_env(monkeypatch, tmp_path, "replace_ok.duckdb")

    class GW:
        def submit_quote(self, rid, yb, nb):
            return "qid-new"

    main._discovery_tick(src, GW(), dry_run=False)

    with db.connect(read_only=True) as con:
        rows = dict(con.execute(
            "SELECT quote_id, status FROM live_quotes").fetchall())
    assert rows == {"qid-old": "replaced", "qid-new": "open"}, rows
