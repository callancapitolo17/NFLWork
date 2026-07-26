"""Issue #23 item 3: correlation-premium / Frechet gate at quote time.

The #20 dispersion gate asks whether the BOOKS agree with each other. This
gate asks whether their consensus is consistent with Kalshi's own live
single-leg prices — books can agree tightly (thin margin under #19) and still
be jointly wrong. The two gates are independent and log distinct reasons.

The anchor is the leg snapshot #17 already fetches at quote time, so the gate
costs no extra API calls. Scaffolding mirrors test_size_gate.py.
"""
import importlib

import pandas as pd

from kalshi_mlb_mm import main

SPREAD_T = "KXMLBSPREAD-25JUN271905TEXLAA-LAA2"
TOTAL_T = "KXMLBTOTAL-25JUN271905TEXLAA-9"
LEGS = [{"market_ticker": SPREAD_T,
         "event_ticker": "KXMLBGAME-25JUN271905TEXLAA", "side": "yes"},
        {"market_ticker": TOTAL_T,
         "event_ticker": "KXMLBGAME-25JUN271905TEXLAA", "side": "yes"}]

# Marginals ~0.46 and ~0.53 -> independent baseline ~0.244, Frechet [0, 0.46].
HEALTHY = {SPREAD_T: {"yes_bid": 0.45, "yes_ask": 0.47},
           TOTAL_T: {"yes_bid": 0.52, "yes_ask": 0.54}}
# Both marginals ~0.10 -> Frechet hi ~0.10: a 0.55 combo fair is arithmetically
# impossible, the joint would exceed its smallest leg.
LONGSHOT_LEGS = {SPREAD_T: {"yes_bid": 0.09, "yes_ask": 0.11},
                 TOTAL_T: {"yes_bid": 0.09, "yes_ask": 0.11}}
# yes_ask = 1.00 on one leg -> no usable anchor at all.
DEGENERATE = {SPREAD_T: {"yes_bid": 0.0, "yes_ask": 1.00},
              TOTAL_T: {"yes_bid": 0.52, "yes_ask": 0.54}}


def _scaffold(monkeypatch, tmp_path, db, cfg, risk, *, snapshot, fair):
    """Discovery tick primed to reach the corr_sanity gate.

    `snapshot` is the live marginal anchor the leg fetch returns; `fair` is the
    book-consensus combo fair the router hands back.
    """
    monkeypatch.setattr(cfg, "DB_PATH", tmp_path / "corr.duckdb")
    importlib.reload(db)
    db.init_database()
    monkeypatch.setattr(main, "_SGP_ODDS",
                        pd.DataFrame({"game_id": ["g1"], "combo": ["c"],
                                      "period": ["FG"], "bookmaker": ["dk"],
                                      "sgp_decimal": [2.0], "fetch_time": [None],
                                      "spread_line": [-1.5], "total_line": [8.5]}))
    monkeypatch.setattr(cfg, "KILL_FILE", tmp_path / ".kill")
    monkeypatch.setattr(risk, "tipoff_ok", lambda ct, min_: True)
    monkeypatch.setattr(main, "_SCOPE_CACHE", {"COMBO-CS": (True, "g1", LEGS)})
    monkeypatch.setattr(main, "_resolve_game_for_legs", lambda gl: "g1")
    monkeypatch.setattr(main, "_leg_market_prices", lambda legs: snapshot)
    monkeypatch.setattr(main, "_PREV_BOOK_FAIR", {})
    monkeypatch.setattr(main.router, "combo_fair_detail",
                        lambda *a, **kw: (main.router.ComboFair(fair, 0.0, 1), "ok"))
    monkeypatch.setattr(cfg, "BANKROLL", 500.0)
    monkeypatch.setattr(cfg, "MAX_FILL_EXPOSURE_PCT", 0.10)


def _gate(monkeypatch, cfg, *, frechet=True, premium=False,
          pmin=0.5, pmax=2.0):
    monkeypatch.setattr(cfg, "CORR_SANITY_FRECHET_ENABLED", frechet)
    monkeypatch.setattr(cfg, "CORR_SANITY_PREMIUM_ENABLED", premium)
    monkeypatch.setattr(cfg, "CORR_PREMIUM_MIN", pmin)
    monkeypatch.setattr(cfg, "CORR_PREMIUM_MAX", pmax)


def _capture_events(monkeypatch):
    events = []
    monkeypatch.setattr(main.research, "emit",
                        lambda ev, **kw: events.append((ev, kw)))
    return events


def _corr_event(events):
    checks = [kw for ev, kw in events if ev == "corr_sanity_check"]
    assert len(checks) == 1, f"expected exactly one corr_sanity_check, got {len(checks)}"
    return checks[0]["payload"]


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
        return "qid-cs"


class _Src:
    def __init__(self, rid):
        self.rid = rid

    def poll(self):
        return [{"id": self.rid, "market_ticker": "COMBO-CS",
                 "contracts_fp": "5.00"}]

    def get_market(self, t):
        return {}


# --------------------------------------------------------------------------- #
# Frechet — the hard, parameter-free reject                                   #
# --------------------------------------------------------------------------- #

def test_frechet_violation_rejects_the_quote(monkeypatch, tmp_path):
    """Two ~10c legs cannot combine to a 55c joint: the joint would exceed its
    own smallest marginal. No band, no tuning — arithmetic."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    _scaffold(monkeypatch, tmp_path, db, cfg, risk,
              snapshot=LONGSHOT_LEGS, fair=0.55)
    _gate(monkeypatch, cfg)
    gw = _GW()
    main._discovery_tick(_Src("r-fre"), gw, dry_run=False)
    assert gw.calls == [], "a Frechet-violating fair must never be quoted"
    assert _last_decision(db, "r-fre") == ("skipped", "corr_sanity_frechet")


def test_frechet_gate_can_be_switched_off(monkeypatch, tmp_path):
    """Config-off leaves the violation logged but ungated (measurement mode)."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    _scaffold(monkeypatch, tmp_path, db, cfg, risk,
              snapshot=LONGSHOT_LEGS, fair=0.55)
    _gate(monkeypatch, cfg, frechet=False)
    events = _capture_events(monkeypatch)
    gw = _GW()
    main._discovery_tick(_Src("r-freoff"), gw, dry_run=False)
    assert len(gw.calls) == 1
    assert _corr_event(events)["reason"] == "frechet"


# --------------------------------------------------------------------------- #
# Premium band — heuristic, ships log-only                                    #
# --------------------------------------------------------------------------- #

def test_premium_band_logs_but_does_not_gate_by_default(monkeypatch, tmp_path):
    """fair 0.10 on a 0.244 baseline is a 0.41x premium — outside [0.5, 2.0]
    but inside Frechet [0, 0.46]. Default config must still quote it: the band
    is a guess until the firehose shows the real distribution."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    _scaffold(monkeypatch, tmp_path, db, cfg, risk, snapshot=HEALTHY, fair=0.10)
    _gate(monkeypatch, cfg, premium=False)
    events = _capture_events(monkeypatch)
    gw = _GW()
    main._discovery_tick(_Src("r-prem"), gw, dry_run=False)
    assert len(gw.calls) == 1, "premium band must not gate while log-only"
    payload = _corr_event(events)
    assert payload["reason"] == "premium"
    assert payload["premium"] < 0.5


def test_premium_band_rejects_once_enabled(monkeypatch, tmp_path):
    """Same inputs, gate switched on -> distinct reason, no quote."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    _scaffold(monkeypatch, tmp_path, db, cfg, risk, snapshot=HEALTHY, fair=0.10)
    _gate(monkeypatch, cfg, premium=True)
    gw = _GW()
    main._discovery_tick(_Src("r-premon"), gw, dry_run=False)
    assert gw.calls == []
    assert _last_decision(db, "r-premon") == ("skipped", "corr_sanity_premium")


# --------------------------------------------------------------------------- #
# Pass-through cases                                                          #
# --------------------------------------------------------------------------- #

def test_healthy_marginals_quote_normally(monkeypatch, tmp_path):
    """fair 0.25 on marginals ~0.46/0.53: premium ~1.02, inside Frechet."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    _scaffold(monkeypatch, tmp_path, db, cfg, risk, snapshot=HEALTHY, fair=0.25)
    _gate(monkeypatch, cfg, premium=True)      # both gates armed
    events = _capture_events(monkeypatch)
    gw = _GW()
    main._discovery_tick(_Src("r-ok"), gw, dry_run=False)
    assert len(gw.calls) == 1
    assert _last_decision(db, "r-ok")[0] == "quoted"
    assert _corr_event(events)["reason"] is None


def test_degenerate_marginals_skip_the_gate_without_crashing(monkeypatch, tmp_path):
    """yes_ask=100 on a leg is UNINFORMATIVE, not alarming: the quote proceeds
    on book consensus alone (the pre-ticket behaviour) and the firehose records
    the miss so the blind rate is measurable."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    _scaffold(monkeypatch, tmp_path, db, cfg, risk,
              snapshot=DEGENERATE, fair=0.55)
    _gate(monkeypatch, cfg, premium=True)      # armed, but has nothing to judge
    events = _capture_events(monkeypatch)
    gw = _GW()
    main._discovery_tick(_Src("r-degen"), gw, dry_run=False)
    assert len(gw.calls) == 1
    payload = _corr_event(events)
    assert payload["marginals"] is None
    assert payload["premium"] is None and payload["reason"] is None


# --------------------------------------------------------------------------- #
# Item 5: the calibration dataset                                             #
# --------------------------------------------------------------------------- #

def test_corr_sanity_check_carries_the_calibration_fields(monkeypatch, tmp_path):
    """Premium + marginals per quote is the dataset a Phase-3 correlation model
    is fitted on, and the evidence for where the band should sit."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    _scaffold(monkeypatch, tmp_path, db, cfg, risk, snapshot=HEALTHY, fair=0.25)
    _gate(monkeypatch, cfg)
    events = _capture_events(monkeypatch)
    main._discovery_tick(_Src("r-evt"), _GW(), dry_run=False)
    payload = _corr_event(events)
    assert payload["n_legs"] == 2
    assert len(payload["marginals"]) == 2
    assert 0.45 < payload["marginals"][0] < 0.47
    assert abs(payload["baseline_independent"]
               - payload["marginals"][0] * payload["marginals"][1]) < 1e-12
    assert abs(payload["premium"]
               - 0.25 / payload["baseline_independent"]) < 1e-12
    assert payload["frechet_lo"] == 0.0
    assert payload["frechet_hi"] == min(payload["marginals"])
    assert payload["combo_fair"] == 0.25


def test_gate_is_independent_of_the_dispersion_gate(monkeypatch, tmp_path):
    """#20 sees zero dispersion (one synthetic book, sigma 0) and #19 therefore
    prices a thin margin — yet the fair is still arithmetically impossible
    against the live marginals. Only this gate can catch that."""
    import kalshi_mlb_mm.config as cfg
    import kalshi_mlb_mm.db as db
    import kalshi_mlb_mm.risk as risk
    _scaffold(monkeypatch, tmp_path, db, cfg, risk,
              snapshot=LONGSHOT_LEGS, fair=0.55)
    _gate(monkeypatch, cfg)
    gw = _GW()
    main._discovery_tick(_Src("r-indep"), gw, dry_run=False)
    assert _last_decision(db, "r-indep") == ("skipped", "corr_sanity_frechet")
