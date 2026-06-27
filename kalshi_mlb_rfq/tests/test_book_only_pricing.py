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


# ---------------------------------------------------------------------------
# Task 2: team-signed leg routing
# ---------------------------------------------------------------------------

import kalshi_mlb_rfq.main as m
from kalshi_common import fair_value


def _region(team_is_home, side, line_n=2, total_n=8, total_side="yes"):
    typed = [
        fair_value.SpreadLeg(team_is_home=team_is_home, line_n=line_n, side=side),
        fair_value.TotalLeg(line_n=total_n, side=total_side),
    ]
    return m._combo_region_from_legs(typed)


def test_routing_table_signs_line_by_team():
    """The 4-row routing table: grid sign by whose margin market it is,
    cell by yes/no. n=2 ⇒ |line|=1.5."""
    # home margin → negative grid
    assert _region(True, "yes").spread_line == -1.5
    assert _region(True, "yes").spread_side == "home"
    assert _region(True, "no").spread_line == -1.5
    assert _region(True, "no").spread_side == "away"
    # away margin → positive grid
    assert _region(False, "yes").spread_line == 1.5
    assert _region(False, "yes").spread_side == "away"
    assert _region(False, "no").spread_line == 1.5
    assert _region(False, "no").spread_side == "home"


# ---------------------------------------------------------------------------
# Task 4: Both-grids integration regression
# ---------------------------------------------------------------------------

def _seed_both_grids() -> pd.DataFrame:
    """Seed BOTH grids (neg + pos) for game 'G' across 2 books with
    DELIBERATELY asymmetric odds so a mis-route is detectable.

    NEG grid (spread_line=-1.5, home-favorite):
        "Away Spread + *" cells are cheap (decimal ~1.6, high-prob ~0.55).

    POS grid (spread_line=+1.5, away-favorite):
        "Away Spread + *" cells are expensive (decimal ~5.0, low-prob ~0.18).

    A correct away-margin YES routing selects the POS grid → fair < 0.35.
    A mis-route to the NEG grid's "Away Spread + Over" → fair > 0.50,
    which fails the assertion.
    """
    NEG = {
        "Home Spread + Over":  5.0,
        "Home Spread + Under": 5.0,
        "Away Spread + Over":  1.6,
        "Away Spread + Under": 1.6,
    }
    POS = {
        "Home Spread + Over":  1.6,
        "Home Spread + Under": 1.6,
        "Away Spread + Over":  5.0,
        "Away Spread + Under": 5.0,
    }
    rows = []
    for book in ("draftkings", "fanduel"):
        for label, dec in NEG.items():
            rows.append({
                "game_id": "G", "combo": label, "bookmaker": book,
                "sgp_decimal": dec, "spread_line": -1.5,
                "total_line": 7.5, "period": "FG",
            })
        for label, dec in POS.items():
            rows.append({
                "game_id": "G", "combo": label, "bookmaker": book,
                "sgp_decimal": dec, "spread_line": 1.5,
                "total_line": 7.5, "period": "FG",
            })
    return pd.DataFrame(rows)


def test_away_margin_prices_to_positive_grid_not_complement(monkeypatch):
    """Away-margin YES leg must price against the +1.5 grid's 'Away Spread + Over'
    cell (~5.0 decimal → devigged ~0.18), NOT the −1.5 grid's same-label cell
    (~1.6 decimal → devigged ~0.55). Regression for Task-2 sign fix.

    With the correct sign (team_is_home=False → +1.5 grid):
        fair ≈ 0.18  → passes (< 0.35)

    If the Task-2 sign fix is reverted (always negated → −1.5 grid):
        fair ≈ 0.55  → fails (> 0.35)
    """
    import statistics
    import kalshi_mlb_rfq.config as cfg

    monkeypatch.setattr(cfg, "MIN_BOOK_COUNT_FOR_BLEND", 2)
    monkeypatch.setattr(m, "_SGP_ODDS_CACHE", _seed_both_grids())

    # away-margin YES leg: _region(team_is_home=False, side="yes") uses
    # line_n=2 → spread_line=+1.5, spread_side="away", total_side="over",
    # total_line=7.5 (total_n=8 → 8−0.5).
    region = _region(team_is_home=False, side="yes")
    assert region.spread_line == 1.5, (
        f"region.spread_line should be +1.5 (away-favorite grid), got {region.spread_line}"
    )

    books = m._load_book_fairs(
        "G", region.spread_line, region.total_line,
        region.spread_side, region.total_side,
    )
    assert books, "expected ≥1 book fair, got empty dict"
    fair = statistics.median(books.values())
    # POS grid "Away Spread + Over" deviggs to ~0.18 — well below 0.35.
    # If mis-routed to NEG grid, fair ≈ 0.55 and this assertion fails.
    assert fair < 0.35, (
        f"away-margin mis-routed to complement cell on −1.5 grid: fair={fair:.3f} (expected < 0.35)"
    )
