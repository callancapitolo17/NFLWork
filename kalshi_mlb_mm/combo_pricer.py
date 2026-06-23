"""General combo fair-value pricer for the maker.

Replaces the maker's single-shape (one game, spread+total) pricing with a
group-by-game pricer that handles arbitrary MLB combos:

  1. Group legs by game (legs of the same game share an event-ticker suffix).
  2. Price each game's sub-combo to ONE consensus fair:
       - single leg (moneyline / spread / total) -> 2-way singles devig,
         consensus across the Odds API book universe.
       - exactly {spread, total} same game -> 4-cell SGP-grid devig
         (captures within-game correlation), consensus across SGP books.
       - anything else (3-leg same game, spread+moneyline, ...) -> unpriceable
         (returns None), which drops the WHOLE combo. Safe: we never assume
         independence between correlated same-game legs.
  3. Combo fair = product of per-game fairs (games are independent).
  4. Combo agreeing-book count = MIN across groups (weakest link).

The per-group consensus gate is the v1 correlation/adverse-selection defense,
applied within each game's own book universe.
"""
import statistics

import pandas as pd

from kalshi_common import fair_value as fv
from kalshi_common.fair_value import devig_book, devig_two_way
from kalshi_common.leg_types import (
    _leg_dict_to_typed, _spread_line_from_legs, _total_line_from_legs,
)


def group_legs_by_game(legs: list[dict]) -> dict[str, list[dict]]:
    """Group legs by their game (event-ticker suffix). KXMLBSPREAD-26.. ,
    KXMLBTOTAL-26.. and KXMLBGAME-26.. for the same physical game all share the
    suffix after the first '-', so that suffix is the game key."""
    groups: dict[str, list[dict]] = {}
    for leg in legs:
        et = str(leg.get("event_ticker", ""))
        suffix = et.split("-", 1)[1] if "-" in et else et
        groups.setdefault(suffix, []).append(leg)
    return groups


def consensus_filter(book_fairs: dict[str, float], min_agreeing: int,
                     band: float) -> dict[str, float]:
    """Median + ±band outlier rejection (v1 correlation defense). Returns the
    surviving (agreeing) books, or {} if fewer than min_agreeing survive."""
    if len(book_fairs) < min_agreeing:
        return {}
    med = statistics.median(book_fairs.values())
    agreeing = {b: f for b, f in book_fairs.items() if abs(f - med) <= band}
    if len(agreeing) < min_agreeing:
        return {}
    return agreeing


def _consensus_point(book_fairs: dict[str, float], min_agreeing: int,
                     band: float):
    """Returns (median_fair, n_agreeing) over the consensus survivors, or None."""
    agreeing = consensus_filter(book_fairs, min_agreeing, band)
    if not agreeing:
        return None
    return statistics.median(agreeing.values()), len(agreeing)


def _spread_total_combo_string(typed_legs) -> str | None:
    """Map a {SpreadLeg, TotalLeg} pair to the SGP-grid combo label that matches
    the STATED sides, e.g. 'Away Spread + Under'. (The old maker hardcoded
    'Home Spread + Over', which mis-priced every other side combination.)"""
    spread = next((l for l in typed_legs if isinstance(l, fv.SpreadLeg)), None)
    total = next((l for l in typed_legs if isinstance(l, fv.TotalLeg)), None)
    if spread is None or total is None:
        return None
    # At a half-point line exactly one team covers, so spread 'no' == other team.
    home_covers = ((spread.team_is_home and spread.side == "yes")
                   or (not spread.team_is_home and spread.side == "no"))
    spread_part = "Home" if home_covers else "Away"
    total_part = "Over" if total.side == "yes" else "Under"
    return f"{spread_part} Spread + {total_part}"


def price_spread_total_group(group_legs, typed_legs, game_id, sgp_df,
                             min_agreeing, band):
    """Consensus fair for a same-game spread+total pair via the 4-cell SGP grid."""
    if sgp_df is None or sgp_df.empty:
        return None
    combo_str = _spread_total_combo_string(typed_legs)
    if combo_str is None:
        return None
    spread_line = _spread_line_from_legs(group_legs)
    total_line = _total_line_from_legs(group_legs)
    rows = sgp_df[(sgp_df.game_id == game_id)
                  & (sgp_df.spread_line.astype(float).round(2) == round(spread_line, 2))
                  & (sgp_df.total_line.astype(float).round(2) == round(total_line, 2))]
    if rows.empty:
        return None
    out: dict[str, float] = {}
    for book in rows.bookmaker.unique():
        sub = rows[rows.bookmaker == book]
        if len(sub) < 4:                  # require the full 4-side grid (no fallback)
            continue
        f = devig_book(sub, combo=combo_str, vig_fallback=0.0)
        if f is not None:
            out[book] = f
    return _consensus_point(out, min_agreeing, band)


def _two_way_leg_fair(decimal_a, decimal_b, want_a: bool) -> float:
    """Devig a 2-outcome market; return the fair of the wanted outcome."""
    p_a, p_b = devig_two_way(float(decimal_a), float(decimal_b))
    return p_a if want_a else p_b


def price_single_leg_group(typed_leg, game_id, singles_df, min_agreeing, band):
    """Consensus fair for ONE moneyline/spread/total leg via 2-way singles devig
    across the Odds API book universe."""
    if singles_df is None or singles_df.empty:
        return None
    g = singles_df[singles_df.game_id == game_id]
    if g.empty:
        return None
    out: dict[str, float] = {}

    if isinstance(typed_leg, fv.MoneylineLeg):
        m = g[g.market == "moneyline"]
        for book in m.bookmaker.unique():
            b = m[m.bookmaker == book]
            home = b[b.outcome == "home"]
            away = b[b.outcome == "away"]
            if home.empty or away.empty:
                continue
            # devig (home, away); we want the leg's team, adjusted for side.
            p_home = _two_way_leg_fair(home.decimal.iloc[0], away.decimal.iloc[0], True)
            p_team = p_home if typed_leg.team_is_home else (1.0 - p_home)
            out[book] = p_team if typed_leg.side == "yes" else (1.0 - p_team)

    elif isinstance(typed_leg, fv.TotalLeg):
        line = typed_leg.line_n - 0.5
        m = g[(g.market == "total") & (g.line.astype(float).round(2) == round(line, 2))]
        for book in m.bookmaker.unique():
            b = m[m.bookmaker == book]
            over = b[b.outcome == "over"]
            under = b[b.outcome == "under"]
            if over.empty or under.empty:
                continue
            p_over = _two_way_leg_fair(over.decimal.iloc[0], under.decimal.iloc[0], True)
            out[book] = p_over if typed_leg.side == "yes" else (1.0 - p_over)

    elif isinstance(typed_leg, fv.SpreadLeg):
        line = -(typed_leg.line_n - 0.5)   # home-perspective handicap
        m = g[(g.market == "spread") & (g.line.astype(float).round(2) == round(line, 2))]
        for book in m.bookmaker.unique():
            b = m[m.bookmaker == book]
            home = b[b.outcome == "home"]
            away = b[b.outcome == "away"]
            if home.empty or away.empty:
                continue
            p_home_cover = _two_way_leg_fair(home.decimal.iloc[0], away.decimal.iloc[0], True)
            p_cover = p_home_cover if typed_leg.team_is_home else (1.0 - p_home_cover)
            out[book] = p_cover if typed_leg.side == "yes" else (1.0 - p_cover)
    else:
        return None

    return _consensus_point(out, min_agreeing, band)


def combo_fair(legs, resolve_game_id, sgp_df, singles_df, *,
               min_agreeing, band):
    """Top-level: combo fair probability (YES side) for an arbitrary MLB combo.

    `resolve_game_id(suffix)` maps an event-ticker suffix to a game_id (the same
    id used in sgp_df / singles_df), or None if unknown.

    Returns (fair, min_agreeing_books) or None if ANY game group is unpriceable
    (unknown game, untyped leg, unsupported shape, or failed consensus).
    """
    groups = group_legs_by_game(legs)
    if not groups:
        return None
    product = 1.0
    min_agree = None
    for suffix, glegs in groups.items():
        game_id = resolve_game_id(suffix)
        if game_id is None:
            return None
        typed = [_leg_dict_to_typed(l, game_id) for l in glegs]
        if any(t is None for t in typed):
            return None
        if len(typed) == 1:
            res = price_single_leg_group(typed[0], game_id, singles_df,
                                         min_agreeing, band)
        elif (len(typed) == 2
              and sum(isinstance(t, fv.SpreadLeg) for t in typed) == 1
              and sum(isinstance(t, fv.TotalLeg) for t in typed) == 1):
            res = price_spread_total_group(glegs, typed, game_id, sgp_df,
                                           min_agreeing, band)
        else:
            return None
        if res is None:
            return None
        fair, cnt = res
        product *= fair
        min_agree = cnt if min_agree is None else min(min_agree, cnt)
    return product, min_agree
