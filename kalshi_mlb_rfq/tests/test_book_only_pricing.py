import importlib
import pytest


def test_book_only_fair_median():
    from kalshi_mlb_rfq.main import _book_only_fair
    assert _book_only_fair({"dk": 0.30, "fd": 0.34}) == pytest.approx(0.32)
    assert _book_only_fair({"dk": 0.30, "fd": 0.34, "px": 0.40}) == pytest.approx(0.34)
    assert _book_only_fair({"dk": 0.30, "fd": None}) == pytest.approx(0.30)
    assert _book_only_fair({}) is None
    assert _book_only_fair({"dk": None}) is None


def test_staleness_gate_skipped_when_model_off(monkeypatch):
    import kalshi_mlb_rfq.main as m
    import kalshi_mlb_rfq.config as cfg
    monkeypatch.setattr(cfg, "USE_MODEL", False)
    monkeypatch.setattr(m, "_SAMPLES_META_GENERATED_AT", None)
    # When USE_MODEL is False, staleness gate must pass even with no gen_at
    ok, decision = m._staleness_gate_ok()
    assert ok is True and decision == "passed"


def test_staleness_gate_fires_when_model_on(monkeypatch):
    import kalshi_mlb_rfq.main as m
    import kalshi_mlb_rfq.config as cfg
    from datetime import datetime, timezone, timedelta
    monkeypatch.setattr(cfg, "USE_MODEL", True)
    monkeypatch.setattr(m, "_SAMPLES_META_GENERATED_AT",
                        datetime.now(timezone.utc) - timedelta(hours=2))
    ok, decision = m._staleness_gate_ok()
    assert ok is False and decision == "declined_stale_predictions"


def test_book_implied_cov_fallback_is_rho_one(monkeypatch):
    import kalshi_mlb_rfq.main as m
    from kalshi_mlb_rfq.correlation import ComboRegion
    # grid returns nothing -> fallback min(Pa,Pb); patch _grid_lookup + fairs
    monkeypatch.setattr(m, "_grid_lookup", lambda gid: (lambda *a: None))
    monkeypatch.setattr(m, "_combo_fair_for_region",
                        lambda gid, r: 0.30 if r.total_side == "over" else 0.25)
    a = ComboRegion("home", -1.5, "over", 8.5)
    b = ComboRegion("home", -1.5, "under", 8.5)
    cov = m._book_implied_cov("g1", a, 0.20, b, 0.20)
    # fallback joint = min(0.30, 0.25) = 0.25; cov = (0.25 - 0.30*0.25) / (0.2*0.2)
    assert cov == pytest.approx((0.25 - 0.30 * 0.25) / 0.04)


def test_book_implied_cov_degenerate_price_returns_zero(monkeypatch):
    import kalshi_mlb_rfq.main as m
    from kalshi_mlb_rfq.correlation import ComboRegion
    # Patch fairs so the None-price guard doesn't short-circuit first
    monkeypatch.setattr(m, "_combo_fair_for_region",
                        lambda gid, r: 0.30)
    monkeypatch.setattr(m, "_grid_lookup", lambda gid: (lambda *a: None))
    a = ComboRegion("home", -1.5, "over", 8.5)
    b = ComboRegion("home", -1.5, "under", 8.5)
    # pos_price=0 would ZeroDivisionError in cov_returns — guard must catch it
    assert m._book_implied_cov("g1", a, 0.20, b, 0) == 0.0
    # Also guard new_price degenerate cases
    assert m._book_implied_cov("g1", a, 0, b, 0.20) == 0.0
    assert m._book_implied_cov("g1", a, 1.0, b, 0.20) == 0.0
    assert m._book_implied_cov("g1", a, 0.20, b, 1.0) == 0.0
