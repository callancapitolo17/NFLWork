"""Canonical leg-set normalizer for arbitrary N-leg Kalshi MLB combos.

Pure functions: parse Kalshi leg dicts into typed CanonicalLegs, canonicalize +
hash a leg set, partition by game, and classify each game's sub-combo's pricing
route. Network-free and config-free so it is fully unit-testable.
"""
import hashlib
import json
from dataclasses import dataclass

from kalshi_common import leg_types


@dataclass(frozen=True)
class CanonicalLeg:
    game_id: str          # the leg's event_ticker (all legs of one game share it)
    market_type: str      # "spread" | "total" | "ml"
    line: float | None    # signed home-perspective; None for ml
    side: str             # "home"/"away" (spread, ml) or "over"/"under" (total)


def parse_leg(leg: dict) -> CanonicalLeg | None:
    """One Kalshi leg dict -> CanonicalLeg in home-perspective, or None."""
    et = str(leg.get("event_ticker", ""))
    mt = str(leg.get("market_ticker", ""))
    if not et or not mt:
        return None
    if mt.startswith("KXMLBSPREAD-"):
        typed = leg_types._leg_dict_to_typed(leg, "")   # SpreadLeg or None
        if typed is None:
            return None
        home_covers = ((typed.team_is_home and typed.side == "yes")
                       or (not typed.team_is_home and typed.side == "no"))
        return CanonicalLeg(et, "spread", -(typed.line_n - 0.5),
                            "home" if home_covers else "away")
    if mt.startswith("KXMLBTOTAL-"):
        typed = leg_types._leg_dict_to_typed(leg, "")   # TotalLeg or None
        if typed is None:
            return None
        return CanonicalLeg(et, "total", typed.line_n - 0.5,
                            "over" if typed.side == "yes" else "under")
    if mt.startswith("KXMLBGAME-"):
        ml = leg_types._moneyline_side(leg)             # (team_is_home, side) or None
        if ml is None:
            return None
        team_is_home, side = ml
        home_ml = ((team_is_home and side == "yes")
                   or (not team_is_home and side == "no"))
        return CanonicalLeg(et, "ml", None, "home" if home_ml else "away")
    return None


def parse_legs(legs: list[dict]) -> list[CanonicalLeg] | None:
    """All legs typed, or None if ANY leg is untypeable (fail-safe)."""
    if not legs:
        return None
    out = []
    for leg in legs:
        c = parse_leg(leg)
        if c is None:
            return None
        out.append(c)
    return out
