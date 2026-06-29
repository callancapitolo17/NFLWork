"""Combo pricing router for arbitrary N-leg combos (Phase 1: cross-game multiply
+ cached-grid same-game). Pure functions — no config/DB imports; callers pass
the SGP odds DataFrame, the per-game resolver, and consensus params in.

Same-game shapes beyond the two 2-leg grids return None here (Phase 2 prices
them on-demand). Cross-game = product of per-game consensus fairs (independence).
"""
import statistics

from kalshi_common import legset
from kalshi_common.fair_value import devig_book
from kalshi_common.leg_types import SPREAD_TOTAL_FAMILY, ML_TOTAL_FAMILY


def consensus(book_fairs: dict[str, float], min_books: int,
              band: float) -> float | None:
    if len(book_fairs) < min_books:
        return None
    med = statistics.median(book_fairs.values())
    agree = {b: f for b, f in book_fairs.items() if abs(f - med) <= band}
    if len(agree) < min_books:
        return None
    return statistics.median(agree.values())


def grid_spec(game_legs: list[legset.CanonicalLeg]):
    """(family, spread_line, total_line, target_cell) for a 2-leg grid sub-combo."""
    total = next(l for l in game_legs if l.market_type == "total")
    ou = "Over" if total.side == "over" else "Under"
    spread = next((l for l in game_legs if l.market_type == "spread"), None)
    if spread is not None:
        part = "Home" if spread.side == "home" else "Away"
        return (SPREAD_TOTAL_FAMILY, spread.line, total.line,
                f"{part} Spread + {ou}")
    ml = next(l for l in game_legs if l.market_type == "ml")
    part = "Home" if ml.side == "home" else "Away"
    return (ML_TOTAL_FAMILY, None, total.line, f"{part} ML + {ou}")


def grid_cell_fairs(game_id, family, spread_line, total_line, target_cell,
                    sgp_df) -> dict[str, float]:
    if sgp_df is None or sgp_df.empty:
        return {}
    df = sgp_df
    mask = ((df.game_id == game_id) & df.combo.isin(family)
            & (df.total_line.astype(float).round(2) == round(total_line, 2)))
    if spread_line is None:
        mask &= df.spread_line.isna()
    else:
        mask &= (df.spread_line.astype(float).round(2) == round(spread_line, 2))
    rows = df[mask]
    out = {}
    for book in rows.bookmaker.unique():
        sub = rows[rows.bookmaker == book]
        if len(sub) < 4:                 # require the full 4-cell grid, no fallback
            continue
        f = devig_book(sub, combo=target_cell, vig_fallback=0.0)
        if f is not None:
            out[book] = f
    return out
