"""Size-gate fixes from the first live session (2026-06-08).

Live RFQ polls showed Kalshi denominates size as `contracts_fp` (a fixed-point
STRING like "21.00") — the old gate read the absent `contracts` int, got 0,
and skipped every in-scope RFQ (86/86 missed over 2 days). ~85% of in-scope
flow was dollar-denominated (`target_cost_dollars` with contracts_fp="0.00"),
which needs a post-pricing implied-contracts gate.
"""
import pandas as pd

from kalshi_mlb_mm import main


# ---------------------------------------------------------------------------
# _rfq_requested_contracts — parse the live RFQ shape
# ---------------------------------------------------------------------------
def test_parses_contracts_fp_string():
    assert main._rfq_requested_contracts({"contracts_fp": "21.00"}) == 21.0


def test_contracts_fp_zero_string_is_zero():
    assert main._rfq_requested_contracts({"contracts_fp": "0.00"}) == 0.0


def test_falls_back_to_legacy_int_contracts():
    assert main._rfq_requested_contracts({"contracts": 3}) == 3.0


def test_absent_size_fields_is_zero():
    assert main._rfq_requested_contracts({"id": "x"}) == 0.0


def test_garbage_contracts_fp_is_zero():
    assert main._rfq_requested_contracts({"contracts_fp": "abc"}) == 0.0


# ---------------------------------------------------------------------------
# _implied_contracts — dollar-denominated RFQ size at our quote
# ---------------------------------------------------------------------------
def test_implied_contracts_ten_dollar_rfq_exceeds_cap():
    # $10 at bids (0.52, 0.43): cheapest ask = min(0.48, 0.57) = 0.48
    # -> 10/0.48 = 20.8 contracts, way over MAX_RFQ_CONTRACTS=5.
    implied = main._implied_contracts(10.0, yes_bid=0.52, no_bid=0.43)
    assert 20.0 < implied < 21.0


def test_implied_contracts_one_dollar_rfq_within_cap():
    implied = main._implied_contracts(1.0, yes_bid=0.52, no_bid=0.43)
    assert implied < 5.0


def test_implied_contracts_floors_division_at_one_cent():
    # Pathological bids near 1.0 must not blow up the division.
    implied = main._implied_contracts(10.0, yes_bid=0.999, no_bid=0.999)
    assert implied == 10.0 / 0.01


# ---------------------------------------------------------------------------
# Tick-level: the gates fire with the live RFQ shape
# ---------------------------------------------------------------------------
def _scaffold(monkeypatch, tmp_path, db, cfg, risk):
    """Shared discovery-tick scaffolding (mirrors the daily-cap test)."""
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "size.duckdb")
    import importlib
    importlib.reload(db)
    db.init_database()
    monkeypatch.setattr(main, "_SGP_ODDS",
                        pd.DataFrame({"game_id": ["g1"], "combo": ["c"], "period": ["FG"],
                                      "bookmaker": ["dk"], "sgp_decimal": [2.0],
                                      "fetch_time": [None], "spread_line": [-1.5],
                                      "total_line": [8.5]}))
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    monkeypatch.setattr(risk, "tipoff_ok", lambda ct, min_: True)
    legs = [{"market_ticker": "KXMLBSPREAD-XYZ", "event_ticker": "EVT", "side": "yes"},
            {"market_ticker": "KXMLBTOTAL-XYZ", "event_ticker": "EVT", "side": "yes"}]
    monkeypatch.setattr(main, "_SCOPE_CACHE", {"COMBO-SZ": (True, "g1", legs)})
    monkeypatch.setattr(main, "_resolve_game_and_lines",
                        lambda ticker, legs: ("g1", -1.5, 8.5))
    monkeypatch.setattr(main, "_PREV_BOOK_FAIR", {})


def _last_decision(db, rid):
    with db.connect(read_only=True) as con:
        return con.execute(
            "SELECT decision, reason FROM quote_decisions "
            "WHERE rfq_id=? ORDER BY observed_at DESC LIMIT 1", [rid]).fetchone()


class _GW:
    def __init__(self):
        self.calls = []

    def submit_quote(self, rid, yb, nb):
        self.calls.append((rid, yb, nb))
        return "qid-sz"


def test_tick_contracts_fp_string_over_cap_blocked(monkeypatch, tmp_path):
    """contracts_fp='1000.00' (live shape) must be parsed and blocked by the cap."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    _scaffold(monkeypatch, tmp_path, db, cfg, risk)

    class Src:
        def poll(self):
            return [{"id": "r-big", "market_ticker": "COMBO-SZ",
                     "contracts_fp": "1000.00"}]

        def get_market(self, t):
            return {}

    gw = _GW()
    main._discovery_tick(Src(), gw, dry_run=False)
    assert gw.calls == []
    assert _last_decision(db, "r-big") == ("skipped", "size_gate")


def test_tick_size_unknown_when_both_fields_zero(monkeypatch, tmp_path):
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    _scaffold(monkeypatch, tmp_path, db, cfg, risk)

    class Src:
        def poll(self):
            return [{"id": "r-nosize", "market_ticker": "COMBO-SZ",
                     "contracts_fp": "0.00"}]

        def get_market(self, t):
            return {}

    gw = _GW()
    main._discovery_tick(Src(), gw, dry_run=False)
    assert gw.calls == []
    assert _last_decision(db, "r-nosize") == ("skipped", "size_unknown")


def test_tick_dollar_rfq_over_cap_blocked_post_pricing(monkeypatch, tmp_path):
    """target_cost_dollars=$10 at ~mid prices implies ~20 contracts > cap 5."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    _scaffold(monkeypatch, tmp_path, db, cfg, risk)
    # Pin the cap so this test is independent of the production default.
    monkeypatch.setattr(cfg, "MAX_RFQ_CONTRACTS", 5)
    # 3 agreeing books so pricing runs on a real consensus.
    monkeypatch.setattr(main, "_book_fairs",
                        lambda g, s, t: {"dk": 0.55, "fd": 0.55, "px": 0.56})

    class Src:
        def poll(self):
            return [{"id": "r-usd", "market_ticker": "COMBO-SZ",
                     "contracts_fp": "0.00", "target_cost_dollars": "10.0000"}]

        def get_market(self, t):
            return {}

    gw = _GW()
    main._discovery_tick(Src(), gw, dry_run=False)
    assert gw.calls == []
    assert _last_decision(db, "r-usd") == ("skipped", "size_gate_dollars")


def test_tick_small_dollar_rfq_passes_and_quotes(monkeypatch, tmp_path):
    """target_cost_dollars=$1 implies ~2 contracts <= cap -> quote submitted."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    _scaffold(monkeypatch, tmp_path, db, cfg, risk)
    monkeypatch.setattr(main, "_book_fairs",
                        lambda g, s, t: {"dk": 0.55, "fd": 0.55, "px": 0.56})

    class Src:
        def poll(self):
            return [{"id": "r-usd1", "market_ticker": "COMBO-SZ",
                     "contracts_fp": "0.00", "target_cost_dollars": "1.0000"}]

        def get_market(self, t):
            return {}

    gw = _GW()
    main._discovery_tick(Src(), gw, dry_run=False)
    assert len(gw.calls) == 1, f"expected one submit, got {gw.calls}"
    assert _last_decision(db, "r-usd1")[0] == "quoted"
