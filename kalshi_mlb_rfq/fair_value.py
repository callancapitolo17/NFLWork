"""Fair-value computation for combos: model + sportsbooks + blend."""

from dataclasses import dataclass
from typing import Literal

import pandas as pd


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


def devig_book(book_rows: pd.DataFrame, combo: str,
               vig_fallback: float = 0.0) -> float | None:
    """Devig a single combo's fair value from rows of mlb_sgp_odds.

    book_rows must already be filtered to (game_id, period, bookmaker, spread_line, total_line).
    Uses the 4-side per-game vig method when >=4 sides exist; else falls back
    to (1/decimal_odds) / (1 + vig_fallback).
    """
    if book_rows.empty:
        return None

    target = book_rows.loc[book_rows["combo"] == combo]
    if target.empty:
        return None
    target_decimal = float(target["sgp_decimal"].iloc[0])

    if len(book_rows) >= 4:
        vig_sum = float((1.0 / book_rows["sgp_decimal"]).sum())
        return (1.0 / target_decimal) / vig_sum

    return (1.0 / target_decimal) / (1.0 + vig_fallback)


def blend(model_fair_value: float, book_fairs: dict[str, float]) -> float | None:
    """2-source gate: returns simple mean of model + non-null books, or None
    if fewer than 2 sources.
    """
    sources = [model_fair_value] + [v for v in book_fairs.values() if v is not None]
    if len(sources) < 2:
        return None
    return sum(sources) / len(sources)
