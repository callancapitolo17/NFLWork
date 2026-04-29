"""Per-cycle combo enumeration + priority queue for the RFQ pipeline."""

import hashlib
from dataclasses import dataclass
from typing import Iterable, Literal

from kalshi_mlb_rfq import ticker_map


@dataclass(frozen=True)
class ComboCandidate:
    game_id: str
    legs: tuple[dict, ...]   # each: {market_ticker, event_ticker, side}
    leg_set_hash: str
    descriptor: str          # human-readable label


def canonical_leg_set_hash(legs: Iterable[dict]) -> str:
    keys = sorted(f"{leg['market_ticker']}|{leg['side']}" for leg in legs)
    return hashlib.sha256("\n".join(keys).encode()).hexdigest()


def enumerate_2leg(
    game_id: str,
    event_suffix: str,
    home_code: str,
    away_code: str,
    available_spreads: list[tuple[float, Literal["home", "away"]]],
    available_totals: list[float],
) -> Iterable[ComboCandidate]:
    """Yield every (spread_leg, side) × (total_leg, side) combo for this game.

    available_spreads: list of (line, "home"|"away") tuples; line is the team's
        negative spread (e.g., -1.5). Spread tickers are minted via ticker_map.
    available_totals: list of total lines (e.g., 7.5, 8.5).
    """
    spread_event = f"KXMLBSPREAD-{event_suffix}"
    total_event = f"KXMLBTOTAL-{event_suffix}"

    spread_leg_specs = []
    for line, who in available_spreads:
        team_code = home_code if who == "home" else away_code
        ticker = ticker_map.spread_ticker(event_suffix, team_code, line)
        for side in ("yes", "no"):
            spread_leg_specs.append((ticker, side, who, line))

    total_leg_specs = []
    for line in available_totals:
        ticker = ticker_map.total_ticker(event_suffix, line)
        for side in ("yes", "no"):
            total_leg_specs.append((ticker, side, line))

    for s_ticker, s_side, s_who, s_line in spread_leg_specs:
        for t_ticker, t_side, t_line in total_leg_specs:
            legs = (
                {"market_ticker": s_ticker, "event_ticker": spread_event, "side": s_side},
                {"market_ticker": t_ticker, "event_ticker": total_event, "side": t_side},
            )
            sign = "-" if s_side == "yes" else "+"
            line_abs = abs(s_line)
            ou = "Over" if t_side == "yes" else "Under"
            descriptor = (
                f"{home_code if s_who == 'home' else away_code} "
                f"{sign}{line_abs} + {ou} {t_line}"
            )
            yield ComboCandidate(
                game_id=game_id,
                legs=legs,
                leg_set_hash=canonical_leg_set_hash(legs),
                descriptor=descriptor,
            )


def edge_score(blended_fair: float, kalshi_ref: float) -> float:
    """Edge magnitude for ranking. If Kalshi has no reference price (last_price=0),
    fall back to |fair - 0.5| as a contentious-combo heuristic."""
    if kalshi_ref <= 0 or kalshi_ref >= 1:
        return abs(blended_fair - 0.5)
    return abs(blended_fair - kalshi_ref)


def rank_by_edge(candidates_with_scores: list[tuple[ComboCandidate, float, float]]
                 ) -> list[ComboCandidate]:
    """Input: list of (candidate, blended_fair, kalshi_ref).
       Output: candidates sorted by edge_score descending.
    """
    return [
        c for c, _, _ in sorted(
            candidates_with_scores,
            key=lambda t: edge_score(t[1], t[2]),
            reverse=True,
        )
    ]
