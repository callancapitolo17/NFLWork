import importlib
import pytest
import pandas as pd


# ---------------------------------------------------------------------------
# Regression test: quadrant-aware book pricing (Task 9)
# ---------------------------------------------------------------------------

def _seed_asymmetric_grid() -> pd.DataFrame:
    """4-cell grid with DIFFERENT decimal odds per quadrant.

    Home+Over = 2.0 (50% implied, high probability after devig)
    Home+Under = 5.0 (20% implied, low probability)
    Away+Over  = 3.0 (~33% implied)
    Away+Under = 4.0 (25% implied)

    Because the four cells have different implied odds, devigging any single
    quadrant produces a DIFFERENT probability depending on which label is
    queried.  This is the condition that would have caught the hardcoded
    'Home Spread + Over' bug.
    """
    rows = []
    cells = [
        ("Home Spread + Over",  2.0),
        ("Home Spread + Under", 5.0),
        ("Away Spread + Over",  3.0),
        ("Away Spread + Under", 4.0),
    ]
    for book in ("draftkings", "fanduel"):
        for (label, dec) in cells:
            rows.append({
                "game_id":    "gQ",
                "combo":      label,
                "bookmaker":  book,
                "sgp_decimal": dec,
                "spread_line": -1.5,
                "total_line":   8.5,
                "period":      "FG",
            })
    return pd.DataFrame(rows)


def test_quadrant_distinguishes_price(monkeypatch):
    """_load_book_fairs must return DIFFERENT values for different quadrants.

    This is the regression that would have caught the hardcoded
    combo='Home Spread + Over' bug — before the fix all four quadrants
    returned the same probability.
    """
    import kalshi_mlb_rfq.main as m
    import kalshi_mlb_rfq.config as cfg

    monkeypatch.setattr(cfg, "MIN_BOOK_COUNT_FOR_BLEND", 2)
    monkeypatch.setattr(m, "_SGP_ODDS_CACHE", _seed_asymmetric_grid())

    home_over  = m._load_book_fairs("gQ", -1.5, 8.5, "home", "over")
    away_under = m._load_book_fairs("gQ", -1.5, 8.5, "away", "under")
    home_under = m._load_book_fairs("gQ", -1.5, 8.5, "home", "under")
    away_over  = m._load_book_fairs("gQ", -1.5, 8.5, "away", "over")

    # All four quadrants must return a non-empty dict
    assert home_over,  "home+over returned empty"
    assert away_under, "away+under returned empty"
    assert home_under, "home+under returned empty"
    assert away_over,  "away+over returned empty"

    # Quadrants with different odds must produce different devigged probs
    assert home_over != away_under, (
        "_load_book_fairs returned identical values for home+over and away+under; "
        "quadrant parameter is not being used"
    )
    assert home_over != home_under, (
        "_load_book_fairs returned identical values for home+over and home+under"
    )


def test_load_book_fairs_default_equals_explicit_home_over(monkeypatch):
    """Default call (no spread_side/total_side) must equal explicit ('home','over')."""
    import kalshi_mlb_rfq.main as m
    import kalshi_mlb_rfq.config as cfg

    monkeypatch.setattr(cfg, "MIN_BOOK_COUNT_FOR_BLEND", 2)
    monkeypatch.setattr(m, "_SGP_ODDS_CACHE", _seed_asymmetric_grid())

    default  = m._load_book_fairs("gQ", -1.5, 8.5)
    explicit = m._load_book_fairs("gQ", -1.5, 8.5, "home", "over")
    assert default == explicit, (
        "Default _load_book_fairs call must equal explicit ('home','over')"
    )


# ---------------------------------------------------------------------------
# Original tests
# ---------------------------------------------------------------------------

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
