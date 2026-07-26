"""Combo pricing router for arbitrary N-leg combos (Phase 1: cross-game multiply
+ cached-grid same-game). Pure functions — no config/DB imports; callers pass
the SGP odds DataFrame, the per-game resolver, and consensus params in.

Same-game shapes beyond the two 2-leg grids return None here (Phase 2 prices
them on-demand). Cross-game = product of per-game consensus fairs (independence).

Consensus (issue #20) is a z-space dispersion gate: fair = median of ALL
books, quoted only when the sample stddev of the book fairs in probit space
is <= sigma_z_max. No outlier removal — a dissenting book is as likely the
informed one (news mid-propagation) as a broken scrape, so large dispersion
declines the quote instead of outvoting the dissenter.
"""
import math
import statistics
from dataclasses import dataclass

from scipy.stats import norm

from kalshi_common import legset
from kalshi_common.fair_value import devig_book
from kalshi_common.leg_types import SPREAD_TOTAL_FAMILY, ML_TOTAL_FAMILY

# Devig outputs live strictly inside (0,1); the clip only guards norm.ppf
# against a pathological 0/1 input reaching +-inf.
_PPF_CLIP = 1e-6


@dataclass(frozen=True)
class Consensus:
    fair: float        # median of ALL book fairs (single median, no filtering)
    sigma_pts: float   # sample stddev (ddof=1) of book fairs, probability points
    sigma_z: float     # sample stddev (ddof=1) of norm.ppf(book fairs)
    n_books: int


def _sigma_z(fairs: list[float]) -> float:
    zs = [float(norm.ppf(min(max(f, _PPF_CLIP), 1.0 - _PPF_CLIP)))
          for f in fairs]
    return statistics.stdev(zs)


def consensus(book_fairs: dict[str, float], min_books: int,
              sigma_z_max: float) -> tuple["Consensus | None", str]:
    """Z-space dispersion gate (issue #20). Returns (Consensus, "ok") or
    (None, reason) with reason in {"too_few_books", "consensus_dispersion"}.

    Constant width in z-space = the same amount of *disagreement* at every
    price level: the tolerated absolute gap naturally tightens at the tails
    (~2c at p=0.50 -> ~0.6c at p=0.08), where the old absolute band tolerated
    25% relative disagreement.
    """
    if len(book_fairs) < max(min_books, 2):
        # A single book has sigma == 0 by construction; the count check must
        # refuse it before dispersion is considered.
        return None, "too_few_books"
    fairs = list(book_fairs.values())
    sigma_z = _sigma_z(fairs)
    if sigma_z > sigma_z_max:
        return None, "consensus_dispersion"
    return Consensus(fair=statistics.median(fairs),
                     sigma_pts=statistics.stdev(fairs),
                     sigma_z=sigma_z,
                     n_books=len(fairs)), "ok"


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
                     sigma_z_max: float) -> tuple[float, list[str]] | None:
    """consensus() plus WHICH books participated (research/observability only).

    With the dispersion gate there are no "survivors" — every supplied book
    feeds the median, so a passing gate reports all of them. Delegates to
    consensus() so the two can never drift."""
    cons, _reason = consensus(book_fairs, min_books, sigma_z_max)
    if cons is None:
        return None
    return cons.fair, sorted(book_fairs)


def subcombo_consensus(game_id, game_legs, sgp_df, min_books: int,
                       sigma_z_max: float,
                       on_demand_fairs=None) -> tuple["Consensus | None", str]:
    """Price one game's sub-combo: route via classify_subcombo -> grid/single
    book fairs -> dispersion-gate consensus. Returns (Consensus | None, reason).

    on_demand_fairs (Phase 2): optional pure lookup, leg_set_hash -> {book:
    fair} | None, injected by main (the OnDemandEngine's fresh-results read).
    Default None reproduces Phase 1 routing — grid/single routes never touch
    it, and on_demand routes return (None, "unpriceable").
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
        return None, "unpriceable"
    return consensus(book_fairs, min_books, sigma_z_max)


@dataclass(frozen=True)
class ComboFair:
    fair: float        # product of per-game consensus fairs
    sigma_pts: float   # combo-level fair uncertainty, probability points
    n_games: int


def combo_fair_detail(legs: list[dict], sgp_df, resolve_game, min_books: int,
                      sigma_z_max: float,
                      on_demand_fairs=None) -> tuple["ComboFair | None", str]:
    """Full RFQ: parse -> partition by game -> per-game consensus -> multiply.

    Returns (ComboFair | None, reason); reason is "ok" or the first failing
    game's gate reason ("too_few_books" / "consensus_dispersion") or a
    routing failure ("unparseable" / "unresolved_game" / "unpriceable").

    Combo sigma: for a product of independent per-game estimates, relative
    variances add (same rule as R's sqrt(sum((s/x)^2)) error propagation):
        sigma_combo = fair_combo * sqrt(sum_g (sigma_g / fair_g)^2)
    """
    canon = legset.parse_legs(legs)
    if canon is None:
        return None, "unparseable"
    product = 1.0
    rel_var = 0.0
    n_games = 0
    for _game_key, game_legs in legset.partition_by_game(canon).items():
        game_id = resolve_game(game_legs)
        if game_id is None:
            return None, "unresolved_game"
        cons, reason = subcombo_consensus(game_id, game_legs, sgp_df,
                                          min_books, sigma_z_max,
                                          on_demand_fairs=on_demand_fairs)
        if cons is None:
            return None, reason
        if cons.fair <= 0.0:
            return None, "unpriceable"
        product *= cons.fair
        rel_var += (cons.sigma_pts / cons.fair) ** 2
        n_games += 1
    return ComboFair(fair=product, sigma_pts=product * math.sqrt(rel_var),
                     n_games=n_games), "ok"


def combo_fair(legs: list[dict], sgp_df, resolve_game, min_books: int,
               sigma_z_max: float, on_demand_fairs=None) -> float | None:
    """Fair-only wrapper for call sites that don't need sigma/n_games
    (confirm last-look drift check, risk-sweep drift check)."""
    detail, _reason = combo_fair_detail(legs, sgp_df, resolve_game, min_books,
                                        sigma_z_max,
                                        on_demand_fairs=on_demand_fairs)
    return detail.fair if detail is not None else None
