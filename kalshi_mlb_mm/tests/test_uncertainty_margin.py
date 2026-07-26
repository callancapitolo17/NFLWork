"""Issue #19: uncertainty-scaled margin — probability-point floor + sigma and
game-count scaling.

Per side: margin_pts = max(p·(1−(1+r)^−N), MIN_MARGIN_PTS + K_SIGMA·σ_books),
all-in cost target = p − margin_pts, then the existing fee/grid/step-down
guard. Defaults (σ=0, N=1, floor=0) must reproduce the legacy constant-ROI
quote exactly.
"""
import pytest

from kalshi_mlb_mm import pricing
from kalshi_common.ev_calc import maker_fee_per_contract


def _realized_margin(bid: float, p: float) -> float:
    return p - (bid + maker_fee_per_contract(bid))


# ---------------------------------------------------------------------------
# Floor binds at the longshot end
# ---------------------------------------------------------------------------

def test_longshot_floor_binds():
    """p=0.10, σ=0.008: floor = 0.01 + 1.0·0.008 = 0.018 pts, dwarfing the
    ROI part (0.10·0.03/1.03 ≈ 0.0029). Realized YES margin must be >= floor."""
    q = pricing.quote(0.10, 0.03, sigma_pts=0.008, n_games=1,
                      min_margin_pts=0.01, k_sigma=1.0)
    assert q is not None
    assert _realized_margin(q.yes_bid, 0.10) >= 0.018 - 1e-9
    # NO side: p=0.90 → ROI part 0.0262 > floor 0.018 → ROI governs there.
    assert _realized_margin(q.no_bid, 0.90) >= 0.0262 - 1e-9


def test_longshot_floor_vs_legacy():
    """Same longshot without the floor (legacy call) quotes ~0.3¢ of margin;
    with the floor the YES bid must be strictly lower (wider cushion)."""
    legacy = pricing.quote(0.10, 0.03)
    floored = pricing.quote(0.10, 0.03, sigma_pts=0.008,
                            min_margin_pts=0.01, k_sigma=1.0)
    assert floored.yes_bid < legacy.yes_bid


# ---------------------------------------------------------------------------
# Mid-price continuity: floor does NOT bind at p≈0.50 with calm books
# ---------------------------------------------------------------------------

def test_mid_price_low_sigma_matches_legacy():
    """At p=0.50 the ROI part is 0.5·0.03/1.03 ≈ 0.0146 > 0.01 floor, so with
    tight books the new pricer must produce byte-identical bids to today's."""
    legacy = pricing.quote(0.50, 0.03)
    new = pricing.quote(0.50, 0.03, sigma_pts=0.0, n_games=1,
                        min_margin_pts=0.01, k_sigma=1.0)
    assert new.yes_bid == legacy.yes_bid
    assert new.no_bid == legacy.no_bid


def test_defaults_reproduce_legacy_everywhere():
    """quote(fair, roi) with no new args must equal the explicit-defaults call
    across the practical fair range (regression lock on the refactor)."""
    for fair in (0.06, 0.10, 0.25, 0.50, 0.75, 0.94):
        legacy = pricing.quote(fair, 0.03)
        new = pricing.quote(fair, 0.03, sigma_pts=0.0, n_games=1,
                            min_margin_pts=0.0, k_sigma=0.0)
        assert (legacy is None) == (new is None)
        if legacy is not None:
            assert new.yes_bid == legacy.yes_bid
            assert new.no_bid == legacy.no_bid


# ---------------------------------------------------------------------------
# Game-count scaling: (1+r)^N compounding
# ---------------------------------------------------------------------------

def test_three_game_margin_exceeds_one_game():
    q1 = pricing.quote(0.30, 0.03, n_games=1)
    q3 = pricing.quote(0.30, 0.03, n_games=3)
    assert q3.yes_bid < q1.yes_bid            # wider margin = lower bid
    assert _realized_margin(q3.yes_bid, 0.30) > _realized_margin(q1.yes_bid, 0.30)


def test_n_games_compounding_hand_computed():
    """3-game ROI part at p=0.30: 0.30·(1 − 1.03^-3) ≈ 0.02546 pts."""
    q3 = pricing.quote(0.30, 0.03, n_games=3)
    assert _realized_margin(q3.yes_bid, 0.30) >= 0.30 * (1 - 1.03 ** -3) - 1e-9


# ---------------------------------------------------------------------------
# Margin components exposed for firehose logging
# ---------------------------------------------------------------------------

def test_quote_carries_margin_components():
    q = pricing.quote(0.10, 0.03, sigma_pts=0.008, n_games=2,
                      min_margin_pts=0.01, k_sigma=1.0)
    assert q.sigma_pts == pytest.approx(0.008)
    assert q.n_games == 2
    assert q.floor_pts == pytest.approx(0.018)
    assert q.margin_pts_yes == pytest.approx(max(0.10 * (1 - 1.03 ** -2), 0.018))
    assert q.margin_pts_no == pytest.approx(max(0.90 * (1 - 1.03 ** -2), 0.018))


# ---------------------------------------------------------------------------
# Property sweep: invariants across the p × σ × N grid
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("sigma", [0.0, 0.005, 0.02, 0.05])
@pytest.mark.parametrize("n_games", [1, 2, 3])
def test_invariants_across_grid(sigma, n_games):
    for milli in range(60, 941, 20):
        fair = round(milli * 0.001, 3)
        q = pricing.quote(fair, 0.03, sigma_pts=sigma, n_games=n_games,
                          min_margin_pts=0.01, k_sigma=1.0)
        if q is None:
            continue                      # declining is always a legal outcome
        assert q.yes_bid > 0 and q.no_bid > 0
        assert q.yes_bid + q.no_bid < 1.0
        floor_pts = 0.01 + 1.0 * sigma
        roi_factor = 1 - (1 + 0.03) ** -n_games
        for bid, p in ((q.yes_bid, fair), (q.no_bid, 1.0 - fair)):
            expected = max(p * roi_factor, floor_pts)
            assert _realized_margin(bid, p) >= expected - 1e-9, (
                f"fair={fair} sigma={sigma} N={n_games} p={p} bid={bid}")
        # grid alignment (deci-cent steps)
        assert abs(q.yes_bid * 1000 - round(q.yes_bid * 1000)) < 1e-6
        assert abs(q.no_bid * 1000 - round(q.no_bid * 1000)) < 1e-6
