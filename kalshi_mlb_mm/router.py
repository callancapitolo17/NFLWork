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
        sub = sub.drop_duplicates(subset=["combo"])
        if len(sub) < 4:                 # require the full 4-cell grid, no fallback
            continue
        f = devig_book(sub, combo=target_cell, vig_fallback=0.0)
        if f is not None:
            out[book] = f
    return out


def _single_marginal_spec(leg: legset.CanonicalLeg):
    """(family, fixed_axis, fixed_value, free_axis, cells) for marginalizing a
    single leg, or (None, ...) if not a supported market_type."""
    if leg.market_type == "spread":
        part = "Home" if leg.side == "home" else "Away"
        return (SPREAD_TOTAL_FAMILY, "spread_line", leg.line, "total_line",
                [f"{part} Spread + Over", f"{part} Spread + Under"])
    if leg.market_type == "total":
        ou = "Over" if leg.side == "over" else "Under"
        return (SPREAD_TOTAL_FAMILY, "total_line", leg.line, "spread_line",
                [f"Home Spread + {ou}", f"Away Spread + {ou}"])
    if leg.market_type == "ml":
        part = "Home" if leg.side == "home" else "Away"
        # ml family rows carry spread_line = NULL; marginalize over total
        return (ML_TOTAL_FAMILY, "spread_line", None, "total_line",
                [f"{part} ML + Over", f"{part} ML + Under"])
    return (None, None, None, None, None)


def _marginal_for_group(group_rows, cells) -> float | None:
    """One book's 4-cell grid group -> sum of the two devigged target cells."""
    group_rows = group_rows.drop_duplicates(subset=["combo"])
    if len(group_rows) < 4:
        return None
    total = 0.0
    for c in cells:
        f = devig_book(group_rows, combo=c, vig_fallback=0.0)
        if f is None:
            return None
        total += f
    return total


def single_marginal_fairs(game_id, leg: legset.CanonicalLeg, sgp_df) -> dict[str, float]:
    family, fixed_axis, fixed_value, free_axis, cells = _single_marginal_spec(leg)
    if family is None or sgp_df is None or sgp_df.empty:
        return {}
    df = sgp_df
    mask = (df.game_id == game_id) & df.combo.isin(family)
    if fixed_value is None:
        mask &= df[fixed_axis].isna()
    else:
        mask &= (df[fixed_axis].astype(float).round(2) == round(fixed_value, 2))
    rows = df[mask]
    out = {}
    for book in rows.bookmaker.unique():
        sub = rows[rows.bookmaker == book]
        # pick the first full 4-cell grid group along the free axis
        for _, grp in sub.groupby(free_axis, dropna=False):
            fair = _marginal_for_group(grp, cells)
            if fair is not None:
                out[book] = fair
                break
    return out


def consensus_detail(book_fairs: dict[str, float], min_books: int,
                     band: float) -> tuple[float, list[str]] | None:
    """consensus() plus WHICH books agreed (research/observability only).

    Pure sibling — consensus() itself is untouched (it is part of the
    Phase 1 regression lock). Kept in exact algorithmic lockstep."""
    if len(book_fairs) < min_books:
        return None
    med = statistics.median(book_fairs.values())
    agree = {b: f for b, f in book_fairs.items() if abs(f - med) <= band}
    if len(agree) < min_books:
        return None
    return statistics.median(agree.values()), sorted(agree)


def subcombo_fair(game_id, game_legs, sgp_df, min_books: int,
                  band: float, on_demand_fairs=None) -> float | None:
    """Price one game's sub-combo: route via classify_subcombo -> grid/single book fairs -> consensus.

    on_demand_fairs (Phase 2): optional pure lookup, leg_set_hash -> {book:
    fair} | None, injected by main (the OnDemandEngine's fresh-results read).
    Default None reproduces Phase 1 byte-identically — grid/single routes
    never touch it, and on_demand routes return None.
    """
    route = legset.classify_subcombo(game_legs)
    if route == "single":
        book_fairs = single_marginal_fairs(game_id, game_legs[0], sgp_df)
    elif route in ("grid_spread_total", "grid_ml_total"):
        family, spread_line, total_line, target = grid_spec(game_legs)
        book_fairs = grid_cell_fairs(game_id, family, spread_line, total_line,
                                     target, sgp_df)
    elif route == "on_demand" and on_demand_fairs is not None:
        book_fairs = on_demand_fairs(legset.leg_set_hash(game_legs)) or {}
    else:                       # "on_demand" without lookup, or "unpriceable"
        return None
    return consensus(book_fairs, min_books, band)


def combo_fair(legs: list[dict], sgp_df, resolve_game, min_books: int,
               band: float, on_demand_fairs=None) -> float | None:
    """Full RFQ: parse -> partition by game -> per-game subcombo_fair -> multiply."""
    canon = legset.parse_legs(legs)
    if canon is None:
        return None
    product = 1.0
    for _game_key, game_legs in legset.partition_by_game(canon).items():
        game_id = resolve_game(game_legs)
        if game_id is None:
            return None
        f = subcombo_fair(game_id, game_legs, sgp_df, min_books, band,
                          on_demand_fairs=on_demand_fairs)
        if f is None:
            return None
        product *= f
    return product
