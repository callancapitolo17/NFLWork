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
    game_id: str          # the GAME code, family-independent — see game_id_of
    market_type: str      # "spread" | "total" | "ml"
    line: float | None    # signed home-perspective; None for ml
    side: str             # "home"/"away" (spread, ml) or "over"/"under" (total)


def game_id_of(event_ticker: str) -> str:
    """The family-independent game code inside a Kalshi event ticker.

    Kalshi gives every market family its own event, so one physical game is
    named by several event tickers that agree only on the trailing code:

        KXMLBSPREAD-26JUN102105MILATH  ->  26JUN102105MILATH
        KXMLBTOTAL -26JUN102105MILATH  ->  26JUN102105MILATH
        KXMLBGAME  -26JUN102105MILATH  ->  26JUN102105MILATH

    Issue #71: keying ``CanonicalLeg.game_id`` on the whole event ticker put a
    game's spread and total legs in different partitions, so
    ``classify_subcombo`` never saw a same-game 2-leg set — the correlation
    grids were unreachable and same-game combos were priced as independent
    games multiplied together. Every MLB event ticker in the recorded corpus
    (2.2M legs) has exactly one dash, but a ticker without one degrades to
    itself rather than raising: it still partitions consistently, and the
    downstream game lookup is already fail-safe on an unparseable code.
    """
    _, _, code = event_ticker.partition("-")
    return code or event_ticker


def parse_leg(leg: dict) -> CanonicalLeg | None:
    """One Kalshi leg dict -> CanonicalLeg in home-perspective, or None."""
    et = game_id_of(str(leg.get("event_ticker", "")))
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
        # Issue #70: the line's sign follows the ticker's TEAM (which margin
        # market this is), matching mlb_sgp_odds.spread_line — home margin
        # ticker -> negative, away margin ticker -> positive. The side is the
        # covering team. This keeps the four (team x side) legs distinct;
        # collapsing both teams onto -(n-0.5) selected the away +1.5 grid
        # cell for an away -1.5 leg.
        home_perspective_line = (-(typed.line_n - 0.5) if typed.team_is_home
                                 else (typed.line_n - 0.5))
        return CanonicalLeg(et, "spread", home_perspective_line,
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
    """Partition legs by game_id (the family-independent game code, #71)."""
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
    # Duplicate-market guard (Phase 2): a repeated (market_type, line) within
    # one game is either the same leg twice or a contradictory pair whose
    # joint probability is exactly 0 (Over 8.5 & Under 8.5, both moneylines).
    # RFQ creators choose the leg set, and a book that product-prices
    # contradictory legs would let Route B assign such a combo real value —
    # a craftable pick-off. Nested totals at DIFFERENT lines stay allowed.
    keys = [(l.market_type, l.line) for l in game_legs]
    if len(set(keys)) != len(keys):
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
