"""Fair-value computation for combos: model + sportsbooks + blend."""

import math
import statistics
from dataclasses import dataclass
from typing import Literal

import pandas as pd
from scipy.stats import norm
from scipy.optimize import brentq


@dataclass(frozen=True)
class SpreadLeg:
    team_is_home: bool
    line_n: int                # KXMLBSPREAD ticker suffix N (= 2 for ±1.5)
    side: Literal["yes", "no"]


@dataclass(frozen=True)
class TotalLeg:
    line_n: int                # KXMLBTOTAL ticker suffix N (= 8 for Over 7.5)
    side: Literal["yes", "no"]


Leg = SpreadLeg | TotalLeg


def _hit_mask(samples: pd.DataFrame, leg: Leg) -> pd.Series:
    if isinstance(leg, SpreadLeg):
        if leg.team_is_home:
            base = samples["home_margin"] >= leg.line_n
        else:
            base = samples["home_margin"] <= -leg.line_n
        return base if leg.side == "yes" else ~base
    if isinstance(leg, TotalLeg):
        base = samples["total_final_score"] >= leg.line_n
        return base if leg.side == "yes" else ~base
    raise TypeError(f"unknown leg type: {type(leg)}")


def model_fair(samples: pd.DataFrame, legs: list[Leg]) -> float:
    """Fraction of sample paths where ALL legs hit."""
    if samples.empty or not legs:
        return 0.0
    mask = pd.Series([True] * len(samples), index=samples.index)
    for leg in legs:
        mask &= _hit_mask(samples, leg)
    return float(mask.mean())


def _probit_devig_n(p_raw, eps=1e-9):
    """Internal: n-way probit (additive z-shift) devig.

    n = 2 uses closed-form c = -(z1+z2)/2; n >= 3 uses brentq root-find.
    Spec: docs/superpowers/specs/2026-05-11-probit-devig-design.md
    """
    p_clipped = [min(max(p, eps), 1 - eps) for p in p_raw]
    z = [norm.ppf(p) for p in p_clipped]

    if len(z) == 2:
        c_star = -(z[0] + z[1]) / 2
    else:
        f = lambda c: sum(norm.cdf(zi + c) for zi in z) - 1
        c_star = brentq(f, -5, 5, xtol=1e-9)
    return [float(norm.cdf(zi + c_star)) for zi in z]


def devig_book(book_rows: pd.DataFrame, combo: str,
               vig_fallback: float = 0.0) -> float | None:
    """Devig a single combo's fair value from rows of mlb_sgp_odds.

    book_rows must already be filtered to (game_id, period, bookmaker, spread_line, total_line).
    Uses probit (additive z-shift) on the 4-cell SGP joint distribution when
    >=4 sides exist; falls back to (1/decimal_odds) / (1 + vig_fallback) when
    fewer than 4 sides are visible.
    """
    if book_rows.empty:
        return None

    target = book_rows.loc[book_rows["combo"] == combo]
    if target.empty:
        return None
    target_decimal = float(target["sgp_decimal"].iloc[0])

    if len(book_rows) >= 4:
        raw_probs = [1.0 / d for d in book_rows["sgp_decimal"]]
        devigged = _probit_devig_n(raw_probs)
        # target_idx: index of the row whose 'combo' matches the requested combo
        # (after duplicates collapse, target.iloc[0] is the chosen row)
        target_idx = book_rows.index.get_loc(target.index[0])
        return float(devigged[target_idx])

    # Single-side fallback: heuristic, no devig math possible with 1 cell
    return (1.0 / target_decimal) / (1.0 + vig_fallback)


def blend(model_fair_value: float, book_fairs: dict[str, float]) -> float | None:
    """2-source gate: returns median of model + non-null books, or None
    if fewer than 2 sources.

    Median (vs arithmetic mean) is robust to a single stale or buggy book
    without assuming any source is sharper than the others. With exactly 2
    sources, median == mean by definition.
    """
    sources = [model_fair_value] + [v for v in book_fairs.values() if v is not None]
    if len(sources) < 2:
        return None
    return statistics.median(sources)


# ------------------------------------------------------------------ #
# Phase 2 on-demand pricing (spec: 2026-07-08-kalshi-mlb-mm-on-demand- #
# pricing-design.md §4.2). Two devig routes for arbitrary leg combos. #
# ------------------------------------------------------------------ #

# SGP vig compounds per leg, so the partition overround gate must scale
# with leg count (a fixed cap would reject healthy 8-cell partitions).
# Calibrated LIVE 2026-07-10: FanDuel's real 2-leg ml+total 4-cell
# partition sums to 1.279 implied — books embed heavy SGP margin per
# cell, far above compounded-singles-vig theory (~1.10). The gate's job
# is catching NONSENSE partitions (mispriced/degenerate cells), which
# land far beyond real margin; 0.25/leg keeps that bite (N=2 -> 1.50,
# N=3 -> 1.75) without silently forcing Route B at high-margin books.
PARTITION_OVERROUND_PER_LEG = 0.25


def devig_partition(cell_decimals, n_legs: int) -> float | None:
    """Route A: probit devig over the full 2^N side-combination partition.

    cell_decimals follow legset.enumerate_partition ordering (index 0 =
    target cell). Any missing/insane cell or an overround outside
    [1.0, 1 + PARTITION_OVERROUND_PER_LEG * n_legs] returns None — the
    caller then tries Route B (fail-safe, "full partition or nothing").
    """
    if not cell_decimals or len(cell_decimals) != 2 ** n_legs:
        return None
    raw = []
    for d in cell_decimals:
        try:
            d = float(d)
        except (TypeError, ValueError):
            return None
        if not math.isfinite(d) or d <= 1.0:
            return None
        raw.append(1.0 / d)
    s = sum(raw)
    if not (1.0 <= s <= 1.0 + PARTITION_OVERROUND_PER_LEG * n_legs):
        return None
    return float(_probit_devig_n(raw)[0])


def devig_two_way(dec_a, dec_b) -> tuple[float, float] | None:
    """Exact 2-way probit devig of one market's two sides."""
    try:
        a, b = float(dec_a), float(dec_b)
    except (TypeError, ValueError):
        return None
    if not (math.isfinite(a) and math.isfinite(b)) or a <= 1.0 or b <= 1.0:
        return None
    fa, fb = _probit_devig_n([1.0 / a, 1.0 / b])
    return float(fa), float(fb)


def fair_by_correlation_transfer(sgp_decimal, singles) -> float | None:
    """Route B: book-implied correlation transfer for one SGP price.

    singles[i] = (vigged_implied_of_chosen_side, devigged_fair) per leg.
    fair = prod(fairs) * [SGP implied / prod(vigged)] — the book's SGP
    margin is assumed ≈ its compounded leg vig, cancelling in the ratio.
    Sanity gate: exact Fréchet bounds — for ANY correlation structure a
    joint probability satisfies max(0, Σp−(n−1)) <= p_joint <= min(p).
    (Assumption-free, unlike a multiplier cap, which would wrongly reject
    legitimate high-correlation stacks like run line + moneyline.)
    """
    try:
        dec = float(sgp_decimal)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(dec) or dec <= 1.0 or not singles:
        return None
    prod_vig, prod_fair, fairs = 1.0, 1.0, []
    for pair in singles:
        try:
            vigged, fair = pair
            vigged, fair = float(vigged), float(fair)
        except (TypeError, ValueError):
            return None
        if not (0.0 < vigged < 1.0 and 0.0 < fair < 1.0):
            return None
        prod_vig *= vigged
        prod_fair *= fair
        fairs.append(fair)
    result = prod_fair * (1.0 / dec) / prod_vig
    lo = max(0.0, sum(fairs) - (len(fairs) - 1))
    hi = min(fairs)
    if not (lo <= result <= hi) or not (0.0 < result < 1.0):
        return None
    return float(result)


def book_on_demand_fair(cells, sgp_decimal, singles,
                        n_legs: int) -> tuple[float, str] | None:
    """Route selection for one book: partition preferred, transfer fallback."""
    if cells is not None:
        fair = devig_partition(cells, n_legs)
        if fair is not None:
            return fair, "partition"
    if sgp_decimal is not None and singles:
        fair = fair_by_correlation_transfer(sgp_decimal, singles)
        if fair is not None:
            return fair, "transfer"
    return None
