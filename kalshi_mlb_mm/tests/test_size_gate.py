"""Size-gate fixes: dollar-based per-fill cap.

The per-fill cap is now in DOLLARS (MAX_FILL_EXPOSURE_PCT * BANKROLL = $50 at
the $500 default bankroll). The maker can't choose fill size (quote is all-or-
nothing); the only lever is quote-or-skip. _worst_fill_exposure_usd() computes
worst-case dollars at risk over both sides the creator could take.

_rfq_requested_contracts parsing tests remain valid (still used to distinguish
contract-denominated vs dollar-denominated RFQs).
"""
import pandas as pd

from kalshi_mlb_mm import main
from kalshi_mlb_mm.pricing import Quote


# ---------------------------------------------------------------------------
# _rfq_requested_contracts — parse the live RFQ shape (still valid)
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
# _worst_fill_exposure_usd — dollar-based exposure for both RFQ denominations
# ---------------------------------------------------------------------------
def test_worst_exposure_contract_rfq_large():
    """1000-contract RFQ at q=quote(0.55, 0.05): exposure far over $50."""
    from kalshi_mlb_mm import pricing
    q = pricing.quote(0.55, 0.05)
    rfq = {"contracts_fp": "1000.00"}
    exp = main._worst_fill_exposure_usd(rfq, q)
    yes_ask = 1.0 - q.no_bid
    no_ask = 1.0 - q.yes_bid
    expected = 1000.0 * max(yes_ask, no_ask)
    assert abs(exp - expected) < 1e-6
    assert exp > 50.0, f"Expected exposure >> $50, got {exp}"


def test_worst_exposure_contract_rfq_small():
    """5-contract RFQ at q=quote(0.55, 0.05): exposure under $50."""
    from kalshi_mlb_mm import pricing
    q = pricing.quote(0.55, 0.05)
    rfq = {"contracts_fp": "5.00"}
    exp = main._worst_fill_exposure_usd(rfq, q)
    yes_ask = 1.0 - q.no_bid
    no_ask = 1.0 - q.yes_bid
    expected = 5.0 * max(yes_ask, no_ask)
    assert abs(exp - expected) < 1e-6
    assert exp < 50.0, f"Expected exposure < $50, got {exp}"


def test_worst_exposure_dollar_rfq_small():
    """$10 dollar-RFQ at a near-symmetric quote: exposure ~$10, passes $50 cap."""
    from kalshi_mlb_mm import pricing
    q = pricing.quote(0.50, 0.05)   # symmetric fair
    rfq = {"contracts_fp": "0.00", "target_cost_dollars": "10.0000"}
    exp = main._worst_fill_exposure_usd(rfq, q)
    # Both sides near-symmetric at ~0.5: exposure ≈ $10
    assert exp < 50.0, f"Expected $10 exposure < $50 cap, got {exp}"
    assert exp > 0.0


def test_worst_exposure_dollar_rfq_large():
    """$500 dollar-RFQ at any reasonable quote: exposure way over $50."""
    from kalshi_mlb_mm import pricing
    q = pricing.quote(0.55, 0.05)
    rfq = {"contracts_fp": "0.00", "target_cost_dollars": "500.0000"}
    exp = main._worst_fill_exposure_usd(rfq, q)
    assert exp > 50.0, f"Expected exposure >> $50, got {exp}"


def test_worst_exposure_both_fields_zero():
    """Both fields zero -> exposure = 0.0 (no-op; size_unknown is caught earlier)."""
    from kalshi_mlb_mm import pricing
    q = pricing.quote(0.55, 0.05)
    rfq = {"contracts_fp": "0.00"}
    exp = main._worst_fill_exposure_usd(rfq, q)
    assert exp == 0.0


def test_worst_exposure_degenerate_quote_returns_inf():
    """Degenerate quote (yes_ask <= 0 or no_ask <= 0) returns inf for dollar RFQ."""
    q = Quote(yes_bid=1.0, no_bid=1.0)  # both asks would be 0 or negative
    rfq = {"contracts_fp": "0.00", "target_cost_dollars": "10.0000"}
    exp = main._worst_fill_exposure_usd(rfq, q)
    assert exp == float("inf")


# ---------------------------------------------------------------------------
# Tick-level: the size gate fires correctly with the new dollar cap
# ---------------------------------------------------------------------------
def _scaffold(monkeypatch, tmp_path, db, cfg, risk):
    """Shared discovery-tick scaffolding."""
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
    """contracts_fp='1000.00' -> exposure >> $50 -> blocked with reason='size_gate'."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    _scaffold(monkeypatch, tmp_path, db, cfg, risk)
    monkeypatch.setattr(cfg, "BANKROLL", 500.0)
    monkeypatch.setattr(cfg, "MAX_FILL_EXPOSURE_PCT", 0.10)
    monkeypatch.setattr(main, "_book_fairs",
                        lambda g, s, t: {"dk": 0.55, "fd": 0.55, "px": 0.56})

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


def test_tick_contracts_fp_string_small_passes(monkeypatch, tmp_path):
    """contracts_fp='5.00' -> exposure ~$2.50 < $50 -> quote submitted."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    _scaffold(monkeypatch, tmp_path, db, cfg, risk)
    monkeypatch.setattr(cfg, "BANKROLL", 500.0)
    monkeypatch.setattr(cfg, "MAX_FILL_EXPOSURE_PCT", 0.10)
    monkeypatch.setattr(main, "_book_fairs",
                        lambda g, s, t: {"dk": 0.55, "fd": 0.55, "px": 0.56})

    class Src:
        def poll(self):
            return [{"id": "r-small", "market_ticker": "COMBO-SZ",
                     "contracts_fp": "5.00"}]

        def get_market(self, t):
            return {}

    gw = _GW()
    main._discovery_tick(Src(), gw, dry_run=False)
    assert len(gw.calls) == 1, f"expected one submit, got {gw.calls}"
    assert _last_decision(db, "r-small")[0] == "quoted"


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
    """target_cost_dollars=$500 -> exposure >> $50 -> skipped with reason='size_gate'."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    _scaffold(monkeypatch, tmp_path, db, cfg, risk)
    monkeypatch.setattr(cfg, "BANKROLL", 500.0)
    monkeypatch.setattr(cfg, "MAX_FILL_EXPOSURE_PCT", 0.10)
    monkeypatch.setattr(main, "_book_fairs",
                        lambda g, s, t: {"dk": 0.55, "fd": 0.55, "px": 0.56})

    class Src:
        def poll(self):
            return [{"id": "r-usd", "market_ticker": "COMBO-SZ",
                     "contracts_fp": "0.00", "target_cost_dollars": "500.0000"}]

        def get_market(self, t):
            return {}

    gw = _GW()
    main._discovery_tick(Src(), gw, dry_run=False)
    assert gw.calls == []
    assert _last_decision(db, "r-usd") == ("skipped", "size_gate")


def test_tick_small_dollar_rfq_passes_and_quotes(monkeypatch, tmp_path):
    """target_cost_dollars=$10 -> exposure ~$10 < $50 -> quote submitted."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    _scaffold(monkeypatch, tmp_path, db, cfg, risk)
    monkeypatch.setattr(cfg, "BANKROLL", 500.0)
    monkeypatch.setattr(cfg, "MAX_FILL_EXPOSURE_PCT", 0.10)
    monkeypatch.setattr(main, "_book_fairs",
                        lambda g, s, t: {"dk": 0.55, "fd": 0.55, "px": 0.56})

    class Src:
        def poll(self):
            return [{"id": "r-usd1", "market_ticker": "COMBO-SZ",
                     "contracts_fp": "0.00", "target_cost_dollars": "10.0000"}]

        def get_market(self, t):
            return {}

    gw = _GW()
    main._discovery_tick(Src(), gw, dry_run=False)
    assert len(gw.calls) == 1, f"expected one submit, got {gw.calls}"
    assert _last_decision(db, "r-usd1")[0] == "quoted"
