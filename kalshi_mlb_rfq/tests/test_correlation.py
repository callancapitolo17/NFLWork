import pytest
from kalshi_mlb_rfq.correlation import (
    ComboRegion, frechet_clamp, joint_prob, cov_returns,
)


def _grid(values):
    """values: dict[(spread_line, total_line, spread_side, total_side)] -> prob"""
    def lookup(sl, tl, ss, ts):
        return values.get((sl, tl, ss, ts))
    return lookup


def test_frechet_clamp_bounds():
    # joint can't exceed min(pa,pb) nor fall below max(0, pa+pb-1)
    assert frechet_clamp(0.9, 0.3, 0.4) == pytest.approx(0.3)      # capped to min
    assert frechet_clamp(0.0, 0.8, 0.8) == pytest.approx(0.6)      # raised to lower bound
    assert frechet_clamp(0.1, 0.3, 0.4) == pytest.approx(0.1)      # already valid


def test_joint_same_direction_nested_total_is_tighter_cell():
    # Both Home+Over; A total 8.5, B total 9.5 -> tighter = Over 9.5 at same spread
    a = ComboRegion("home", -1.5, "over", 8.5)
    b = ComboRegion("home", -1.5, "over", 9.5)
    grid = _grid({(-1.5, 9.5, "home", "over"): 0.22})
    j = joint_prob(a, b, p_a=0.30, p_b=0.22, grid_lookup=grid)
    assert j == pytest.approx(0.22)


def test_joint_same_direction_nested_spread_home_takes_min_line():
    # Home favorite: tighter spread = more negative line (higher margin threshold)
    a = ComboRegion("home", -1.5, "over", 8.5)
    b = ComboRegion("home", -2.5, "over", 8.5)
    grid = _grid({(-2.5, 8.5, "home", "over"): 0.18})
    j = joint_prob(a, b, p_a=0.30, p_b=0.20, grid_lookup=grid)
    assert j == pytest.approx(0.18)


def test_joint_away_takes_max_line():
    # Away cover: tighter spread = larger (more positive) home line
    a = ComboRegion("away", 1.5, "under", 9.5)
    b = ComboRegion("away", 2.5, "under", 8.5)
    grid = _grid({(2.5, 8.5, "away", "under"): 0.12})
    j = joint_prob(a, b, p_a=0.25, p_b=0.20, grid_lookup=grid)
    assert j == pytest.approx(0.12)


def test_joint_opposite_direction_returns_none():
    a = ComboRegion("home", -1.5, "over", 8.5)
    b = ComboRegion("home", -1.5, "under", 8.5)
    assert joint_prob(a, b, 0.3, 0.3, _grid({})) is None


def test_joint_opposite_spread_returns_none():
    a = ComboRegion("home", -1.5, "over", 8.5)
    b = ComboRegion("away", -1.5, "over", 8.5)
    assert joint_prob(a, b, 0.3, 0.3, _grid({})) is None


def test_joint_missing_cell_returns_none():
    a = ComboRegion("home", -1.5, "over", 8.5)
    b = ComboRegion("home", -2.5, "over", 8.5)
    assert joint_prob(a, b, 0.3, 0.2, _grid({})) is None


def test_joint_identity_returns_pa_clamped():
    a = ComboRegion("home", -1.5, "over", 8.5)
    grid = _grid({(-1.5, 8.5, "home", "over"): 0.30})
    assert joint_prob(a, a, 0.30, 0.30, grid) == pytest.approx(0.30)


def test_joint_applies_frechet_clamp_to_noisy_grid():
    # grid says joint 0.28 but min(pa,pb)=0.20 -> clamp to 0.20
    a = ComboRegion("home", -1.5, "over", 8.5)
    b = ComboRegion("home", -2.5, "over", 8.5)
    grid = _grid({(-2.5, 8.5, "home", "over"): 0.28})
    j = joint_prob(a, b, p_a=0.30, p_b=0.20, grid_lookup=grid)
    assert j == pytest.approx(0.20)


def test_cov_returns_independent_is_zero():
    assert cov_returns(0.3, 0.4, 0.12, 0.3, 0.4) == pytest.approx(0.0)


def test_cov_returns_positive_when_joint_exceeds_product():
    cov = cov_returns(0.3, 0.4, 0.20, 0.3, 0.4)
    assert cov > 0


def test_cov_returns_negative_when_joint_below_product():
    assert cov_returns(0.3, 0.4, 0.10, 0.3, 0.4) < 0


def test_cov_returns_scales_inverse_with_prices():
    base = cov_returns(0.3, 0.4, 0.20, 0.3, 0.4)
    half = cov_returns(0.3, 0.4, 0.20, 0.15, 0.4)
    assert half == pytest.approx(base * 2)


@pytest.mark.parametrize("pa,pb,raw", [
    (0.3, 0.4, 0.9), (0.8, 0.8, 0.0), (0.1, 0.1, 0.5), (0.5, 0.5, 0.5),
])
def test_frechet_invariant(pa, pb, raw):
    j = frechet_clamp(raw, pa, pb)
    assert max(0.0, pa + pb - 1.0) - 1e-12 <= j <= min(pa, pb) + 1e-12


# ---------------------------------------------------------------------------
# Characterization tests: signed-grid correctness
# ---------------------------------------------------------------------------

from kalshi_mlb_rfq import correlation  # noqa: E402  (imported at module level already, alias for clarity)


def test_joint_two_away_margin_combos_reads_positive_grid():
    """Two away-margin combos (spread_line +1.5 and +2.5, both away cover,
    both over) → tighter = +2.5 away cell. grid_lookup must be called with
    the POSITIVE signed line."""
    calls = []
    def fake_lookup(spread, total, sside, tside):
        calls.append((spread, total, sside, tside))
        return 0.10
    a = ComboRegion("away", 1.5, "over", 7.5)
    b = ComboRegion("away", 2.5, "over", 8.5)
    j = correlation.joint_prob(a, b, p_a=0.30, p_b=0.12, grid_lookup=fake_lookup)
    assert calls == [(2.5, 8.5, "away", "over")]   # tighter: max line, max total
    assert j is not None


def test_joint_opposite_grid_pair_is_none_fallback():
    """A home-margin combo (−1.5) and an away-margin combo (+1.5) have
    different spread_side → None (caller uses ρ=1, conservative)."""
    a = ComboRegion("home", -1.5, "over", 7.5)   # home wins by 2+
    b = ComboRegion("away", 1.5, "over", 7.5)    # away wins by 2+
    j = correlation.joint_prob(a, b, p_a=0.30, p_b=0.15,
                               grid_lookup=lambda *x: 0.5)
    assert j is None
