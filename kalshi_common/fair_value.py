"""Fair-value computation for combos: model + sportsbooks + blend."""

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
