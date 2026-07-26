"""Uncertainty-scaled quote pricing (issue #19).

Per side (p = side win-prob), the margin in probability points is
    margin_pts = max(roi_pts, floor_pts)
    roi_pts    = p * (1 - (1+target_roi)^-n_games)
    floor_pts  = min_margin_pts + k_sigma * sigma_pts
and the target all-in cost (bid + maker fee) is raw = p - margin_pts.

Why each piece:
- roi_pts with (1+r)^N compounding: sportsbooks build multi-leg parlay hold
  by compounding per-leg vig; one ROI divisor per game reproduces that
  structure for cross-game products, whose per-game fair errors compound.
  At n_games=1 this is exactly the legacy constant-ROI cost p/(1+r).
- floor_pts: the constant-ROI cushion p*r/(1+r) ~ 0.029*p collapses at
  longshot fairs (~0.3c at p=0.10) while absolute fair error (book
  dispersion + devig + correlation error) does not shrink with p. Every
  book loads margin onto longshots (favorite-longshot vig distribution);
  options MMs widen deep-OTM quotes for the same reason. sigma_pts prices
  the disagreement we can SEE (stddev of the consensus books' fairs);
  min_margin_pts covers error sigma cannot see (shared devig/correlation
  bias, mirrored books).

The maker fee depends on the bid, so we iterate once. Because the fee is a
stepped (ceil-to-cent) function, flooring the bid to the grid or crossing a
fee tier can leave bid+fee(bid) > raw, pushing the realized margin below
target. A post-floor guard steps the bid down one tick at a time until
cost <= raw, guaranteeing realized margin_pts >= target by construction.
"""
import math
from dataclasses import dataclass

from kalshi_common.ev_calc import maker_fee_per_contract

GRID_STEP = 0.001  # Kalshi MVE markets price in $0.001 (deci-cent) steps


@dataclass(frozen=True)
class Quote:
    yes_bid: float
    no_bid: float
    # Margin components, exposed for research-firehose logging (issue #19
    # item 5). Defaults keep bare Quote(yes_bid=..., no_bid=...) construction
    # (tests, monkeypatched fakes) working.
    sigma_pts: float = 0.0
    n_games: int = 1
    floor_pts: float = 0.0
    roi_pts_yes: float = 0.0
    roi_pts_no: float = 0.0
    margin_pts_yes: float = 0.0
    margin_pts_no: float = 0.0


def _round_down_to_grid(price: float) -> float:
    # Floor to the grid: a lower bid is conservative (cheaper for us / more margin).
    return math.floor(round(price / GRID_STEP, 6)) * GRID_STEP


def _side_targets(p: float, target_roi: float, n_games: int,
                  floor_pts: float) -> tuple[float, float, float]:
    """(roi_pts, margin_pts, raw cost) for one side at win-prob p.

    When the ROI part governs, raw is computed as p / (1+r)^N directly — the
    same expression as the legacy constant-ROI pricer — so defaults reproduce
    it bit-for-bit (p - (p - p/d) can differ by 1 ulp from p/d)."""
    roi_raw = p / (1.0 + target_roi) ** n_games
    roi_pts = p - roi_raw
    if floor_pts > roi_pts:
        return roi_pts, floor_pts, p - floor_pts
    return roi_pts, roi_pts, roi_raw


def _price_for_side(raw: float) -> float:
    # raw = target all-in cost (bid + fee)
    bid = raw - maker_fee_per_contract(raw)
    bid = raw - maker_fee_per_contract(bid)  # one refinement (fee depends on bid)
    bid = _round_down_to_grid(max(bid, 0.0))
    # The maker fee is stepped (ceil-to-cent), so flooring or a fee-tier
    # crossing can leave bid+fee(bid) > raw, i.e. cost above target -> margin
    # below target. Step down one tick until the all-in cost is within target
    # (or the bid hits zero, in which case quote() declines via its <=0
    # guard). This guarantees realized margin_pts >= target by construction.
    while bid > 0 and bid + maker_fee_per_contract(bid) > raw + 1e-12:
        bid = _round_down_to_grid(bid - GRID_STEP)
    return max(bid, 0.0)


def quote(fair: float, target_roi: float, *, sigma_pts: float = 0.0,
          n_games: int = 1, min_margin_pts: float = 0.0,
          k_sigma: float = 0.0) -> "Quote | None":
    """Two-sided quote around fair. Defaults (sigma 0, N=1, floor 0)
    reproduce the legacy constant-ROI pricer bit-for-bit."""
    if not (0.0 < fair < 1.0):
        return None
    if n_games < 1:
        raise ValueError(f"expected n_games >= 1, got {n_games}")
    floor_pts = min_margin_pts + k_sigma * max(sigma_pts, 0.0)
    roi_pts_yes, margin_yes, raw_yes = _side_targets(fair, target_roi,
                                                     n_games, floor_pts)
    roi_pts_no, margin_no, raw_no = _side_targets(1.0 - fair, target_roi,
                                                  n_games, floor_pts)
    yes_bid = _price_for_side(raw_yes)
    no_bid = _price_for_side(raw_no)
    if yes_bid <= 0 or no_bid <= 0 or yes_bid + no_bid >= 1.0:
        return None
    return Quote(yes_bid=round(yes_bid, 4), no_bid=round(no_bid, 4),
                 sigma_pts=sigma_pts, n_games=n_games, floor_pts=floor_pts,
                 roi_pts_yes=roi_pts_yes, roi_pts_no=roi_pts_no,
                 margin_pts_yes=margin_yes, margin_pts_no=margin_no)
