"""Book-implied correlation for same-game combos.

Pure functions — no DB, no globals. The caller injects a `grid_lookup`
callable that reads the per-game spread×total SGP grid.

A combo is a quadrant of the (home margin M, game total T) plane:
  - home cover  iff  M > -spread_line ; away cover iff M < -spread_line
  - over        iff  T > total_line   ; under       iff T < total_line

For two combos sharing both spread side and total side, the joint event
"both hit" is the *tighter* quadrant, which is itself a grid cell — so the
joint is a single grid lookup. Opposite-side pairs (and any leg the grid
cannot price) return None; the caller treats None as the ρ=1 fallback.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class ComboRegion:
    spread_side: str   # "home" | "away"
    spread_line: float # home-perspective line as stored in mlb_sgp_odds
    total_side: str    # "over" | "under"
    total_line: float


def frechet_clamp(joint: float, p_a: float, p_b: float) -> float:
    """Clamp a (possibly noisy) joint into the Fréchet–Hoeffding bounds so the
    implied correlation stays valid."""
    lo = max(0.0, p_a + p_b - 1.0)
    hi = min(p_a, p_b)
    return max(lo, min(joint, hi))


def joint_prob(a: ComboRegion, b: ComboRegion,
               p_a: float, p_b: float, grid_lookup) -> float | None:
    """P(A ∩ B) from the grid, or None if not resolvable (caller → ρ=1).

    Same spread_side + same total_side → the joint is the tighter quadrant,
    one grid cell, read at the signed `spread_line` (negative = home-favorite
    grid, positive = away-favorite grid).

    Opposite spread_side (e.g. a home-margin combo vs an away-margin combo,
    which live in opposite signed grids and are near-disjoint) returns None.
    The caller's ρ=1 fallback then treats them as max positively correlated —
    conservative: it down-sizes rather than crediting a diversification
    benefit we cannot verify from the grid. (v1 limitation, intentional.)
    """
    if a.spread_side != b.spread_side or a.total_side != b.total_side:
        return None
    # Tighter spread quadrant:
    #   home cover is M > -spread_line → tighter = more negative line = min()
    #   away cover is M < -spread_line → tighter = more positive line = max()
    if a.spread_side == "home":
        tight_spread = min(a.spread_line, b.spread_line)
    else:
        tight_spread = max(a.spread_line, b.spread_line)
    # Tighter total quadrant:
    #   over is T > total_line  → tighter = larger total = max()
    #   under is T < total_line → tighter = smaller total = min()
    if a.total_side == "over":
        tight_total = max(a.total_line, b.total_line)
    else:
        tight_total = min(a.total_line, b.total_line)
    raw = grid_lookup(tight_spread, tight_total, a.spread_side, a.total_side)
    if raw is None:
        return None
    return frechet_clamp(raw, p_a, p_b)


def cov_returns(p_a: float, p_b: float, p_joint: float,
                price_a: float, price_b: float) -> float:
    """Return-space covariance for Kelly:
        Cov(r_a, r_b) = [P(A∩B) − P(A)P(B)] / (price_a · price_b)."""
    cov_outcomes = p_joint - p_a * p_b
    return cov_outcomes / (price_a * price_b)
