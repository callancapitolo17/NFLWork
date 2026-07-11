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
        try:
            typed = leg_types._leg_dict_to_typed(leg, "")   # SpreadLeg or None
        except (KeyError, TypeError, ValueError):
            return None
        if typed is None:
            return None
        home_covers = ((typed.team_is_home and typed.side == "yes")
                       or (not typed.team_is_home and typed.side == "no"))
        return CanonicalLeg(et, "spread", -(typed.line_n - 0.5),
                            "home" if home_covers else "away")
    if mt.startswith("KXMLBTOTAL-"):
        try:
            typed = leg_types._leg_dict_to_typed(leg, "")   # TotalLeg or None
        except (KeyError, TypeError, ValueError):
            return None
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


def _sort_key(l: CanonicalLeg):
    # None-safe: ml legs (line=None) sort after numeric lines within a market_type
    return (l.game_id, l.market_type, l.line is None, l.line or 0.0, l.side)


def canonical_legs(legs: list[CanonicalLeg]) -> list[CanonicalLeg]:
    return sorted(legs, key=_sort_key)


def leg_set_hash(legs: list[CanonicalLeg]) -> str:
    payload = [[l.game_id, l.market_type, l.line, l.side]
               for l in canonical_legs(legs)]
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha1(blob.encode()).hexdigest()


def partition_by_game(legs: list[CanonicalLeg]) -> dict[str, list[CanonicalLeg]]:
    """Partition legs by game_id (event_ticker)."""
    out: dict[str, list[CanonicalLeg]] = {}
    for l in legs:
        out.setdefault(l.game_id, []).append(l)
    return out


def classify_subcombo(game_legs: list[CanonicalLeg]) -> str:
    """Classify a single game's legs into a pricing route.

    Returns one of: "single", "grid_spread_total", "grid_ml_total",
    "on_demand", "unpriceable".
    """
    n = len(game_legs)
    if n == 0:
        return "unpriceable"
    if n == 1:
        return "single"
    types = sorted(l.market_type for l in game_legs)
    if n == 2 and types == ["spread", "total"]:
        return "grid_spread_total"
    if n == 2 and types == ["ml", "total"]:
        return "grid_ml_total"
    return "on_demand"


_FLIP = {"home": "away", "away": "home", "over": "under", "under": "over"}


def flip_leg(leg: CanonicalLeg) -> CanonicalLeg:
    """The same leg on its opposite side (line unchanged)."""
    return CanonicalLeg(leg.game_id, leg.market_type, leg.line, _FLIP[leg.side])


def enumerate_partition(game_legs: list[CanonicalLeg]) -> list[list[CanonicalLeg]]:
    """All 2^N side-combinations of a leg set (the joint-outcome partition).

    Cell ordering contract (shared with fair_value.devig_partition and
    SGPService.price_on_demand): cell index i in range(2^N); bit j of i
    (LSB = leg 0) means leg j is FLIPPED to its opposite side. Cell 0 is
    therefore the target (all chosen sides). Lines never change — only
    sides — so per-book leg->selection resolution is shared across cells.
    """
    n = len(game_legs)
    return [[flip_leg(l) if (i >> j) & 1 else l
             for j, l in enumerate(game_legs)]
            for i in range(2 ** n)]
