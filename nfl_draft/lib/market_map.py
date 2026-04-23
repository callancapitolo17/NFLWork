"""Per-book market label -> canonical market_id mapping + market_id construction."""

import re
from typing import Optional
from nfl_draft.lib.db import read_connection


def slug(name: str) -> str:
    """Canonicalize a name to a URL-safe slug (used in market_id construction)."""
    return name.lower().replace(" ", "-").replace(".", "").replace("'", "")


def _slug_underscored(name: str) -> str:
    """Like slug() but uses underscores for spaces (team + position names).

    Used by canonical market_types where the subject is a multi-word team or
    position rather than a player. Keeps lowercase + strips punctuation like
    slug() does so 'Wide Receiver' -> 'wide_receiver' and 'Defensive Line/Edge'
    -> 'defensive_line_edge'.
    """
    out = name.lower().replace(".", "").replace("'", "")
    # Collapse any run of whitespace + slashes into a single underscore.
    out = re.sub(r"[\s/]+", "_", out).strip("_")
    return out


# Team-name normaliser: collapses BetOnline's franchise names ('Arizona Cardinals')
# and Kalshi's short names ('Arizona') onto a single canonical form so cross-book
# joins on 'team_first_pick' actually fire. The canonical form chosen is the
# Kalshi short name because those strings are already live in the draft_odds
# table; remapping them would break historical rows. For ambiguous cities with
# two franchises (LA, NY), Kalshi disambiguates via single-letter suffix and
# we map the BetOnline full name to the matching Kalshi form.
_NFL_TEAM_NORMALISE: dict[str, str] = {
    # Single-franchise cities: both forms collapse to the Kalshi short name.
    "arizona": "Arizona", "arizona cardinals": "Arizona",
    "atlanta": "Atlanta", "atlanta falcons": "Atlanta",
    "baltimore": "Baltimore", "baltimore ravens": "Baltimore",
    "buffalo": "Buffalo", "buffalo bills": "Buffalo",
    "carolina": "Carolina", "carolina panthers": "Carolina",
    "chicago": "Chicago", "chicago bears": "Chicago",
    "cincinnati": "Cincinnati", "cincinnati bengals": "Cincinnati",
    "cleveland": "Cleveland", "cleveland browns": "Cleveland",
    "dallas": "Dallas", "dallas cowboys": "Dallas",
    "denver": "Denver", "denver broncos": "Denver",
    "detroit": "Detroit", "detroit lions": "Detroit",
    "green bay": "Green Bay", "green bay packers": "Green Bay",
    "houston": "Houston", "houston texans": "Houston",
    "indianapolis": "Indianapolis", "indianapolis colts": "Indianapolis",
    "jacksonville": "Jacksonville", "jacksonville jaguars": "Jacksonville",
    "kansas city": "Kansas City", "kansas city chiefs": "Kansas City",
    "las vegas": "Las Vegas", "las vegas raiders": "Las Vegas",
    "miami": "Miami", "miami dolphins": "Miami",
    "minnesota": "Minnesota", "minnesota vikings": "Minnesota",
    "new england": "New England", "new england patriots": "New England",
    "new orleans": "New Orleans", "new orleans saints": "New Orleans",
    "philadelphia": "Philadelphia", "philadelphia eagles": "Philadelphia",
    "pittsburgh": "Pittsburgh", "pittsburgh steelers": "Pittsburgh",
    "san francisco": "San Francisco", "san francisco 49ers": "San Francisco",
    "seattle": "Seattle", "seattle seahawks": "Seattle",
    "tampa bay": "Tampa Bay", "tampa bay buccaneers": "Tampa Bay",
    "tennessee": "Tennessee", "tennessee titans": "Tennessee",
    "washington": "Washington", "washington commanders": "Washington",
    # Two-franchise cities: Kalshi disambiguates via single-letter suffix.
    "los angeles c": "Los Angeles C", "los angeles chargers": "Los Angeles C",
    "los angeles r": "Los Angeles R", "los angeles rams": "Los Angeles R",
    "new york g": "New York G", "new york giants": "New York G",
    "new york j": "New York J", "new york jets": "New York J",
}


def normalize_team(name: str) -> str:
    """Return the canonical NFL team name for cross-book joins, falling back
    to the input unchanged if unknown (so novel names surface in quarantine
    rather than silently mis-joining onto a wrong team)."""
    return _NFL_TEAM_NORMALISE.get(name.strip().lower(), name.strip())


def _format_line(line: float) -> str:
    """Encode a pick line into the market_id safely.

    9.5 -> '9p5' (dot would be lexically noisy in the ID).
    3.0 -> '3' (drop trailing .0 so the ID doesn't gain a decimal needlessly).
    """
    s = f"{line:g}"  # 9.5 -> '9.5', 3.0 -> '3'
    return s.replace(".", "p")


def build_market_id(market_type: str, **kwargs) -> str:
    """Deterministic market_id constructor -- must be applied identically in every scraper.

    Tested by test_market_map.py.
    """
    if market_type == "first_at_position":
        return f"first_{kwargs['position'].lower()}_{slug(kwargs['player'])}"
    if market_type == "pick_outright":
        return f"pick_{kwargs['pick_number']}_overall_{slug(kwargs['player'])}"
    if market_type == "top_n_range":
        return f"top_{kwargs['range_high']}_{slug(kwargs['player'])}"
    if market_type == "team_first_pick":
        return f"team_{kwargs['team'].lower()}_first_pick_{slug(kwargs['player'])}"
    if market_type == "prop":
        return f"prop_{slug(kwargs['short_description'])}"
    if market_type == "nth_at_position":
        return f"{kwargs['nth']}_{kwargs['position'].lower()}_{slug(kwargs['player'])}"
    if market_type == "mr_irrelevant_position":
        return f"mr_irrelevant_{_slug_underscored(kwargs['position'])}"
    if market_type == "team_first_pick_position":
        return (
            f"{_slug_underscored(kwargs['team'])}"
            f"_first_pick_pos_{_slug_underscored(kwargs['position'])}"
        )
    if market_type == "matchup_before":
        a, b = sorted([slug(kwargs["player_a"]), slug(kwargs["player_b"])])
        return f"matchup_{a}_before_{b}"
    if market_type == "draft_position_over_under":
        return (
            f"draft_position_ou_{slug(kwargs['player'])}"
            f"_{_format_line(kwargs['line'])}"
            f"_{kwargs['direction'].lower()}"
        )
    raise ValueError(f"Unknown market_type: {market_type}")


def outright_group_key(market_group: str, market_id: str) -> Optional[str]:
    """Derive the devig-grouping key from (market_group, market_id), or None.

    Scrapers emit ``market_group`` on every ``OddsRow`` (see
    ``scrapers/_base.py``). This helper dispatches on ``market_group`` to
    extract the group instance (pick number, position, team+position, etc.)
    from ``market_id``. Rows that share the same ``(book, outright_group_key)``
    are competing outcomes in the same logical outright and are devigged
    together in ``quarantine.write_or_quarantine``.

    Returns ``None`` for:
      * prop markets (``market_group`` starting with ``prop_``)
      * ``mr_irrelevant_position`` (single row per book — nothing to devig)
      * ``matchup_before`` (self-contained 2-way, not grouped across markets)
      * unrecognized / malformed inputs

    The regex patterns here are tightly coupled to ``build_market_id``'s
    output format. If ``build_market_id`` changes, these must change too.
    Tested in ``tests/unit/test_market_map.py``.
    """
    if market_group == "pick_outright":
        m = re.match(r"^(pick_\d+_overall)_", market_id)
        return m.group(1) if m else None
    if market_group == "first_at_position":
        m = re.match(r"^(first_[a-z]+)_", market_id)
        return m.group(1) if m else None
    # top_N_range markets ("will player X be drafted in top N") are an
    # N-WINNER pool, not a 1-winner mutex outright. The group key is
    # ``top_<N>`` so quarantine.write_or_quarantine can parse N and apply
    # devig_pool(odds, N) -- normalizing posted YES implieds to sum to N
    # instead of 1. If any other market_group ever uses a ``top_<N>``
    # group key shape, it will inherit pool-N devig automatically; keep
    # the format discipline here.
    if market_group.startswith("top_") and market_group.endswith("_range"):
        m = re.match(r"^(top_\d+)_", market_id)
        return m.group(1) if m else None
    if market_group == "team_first_pick":
        m = re.match(r"^(team_.+?_first_pick)_", market_id)
        return m.group(1) if m else None
    if market_group == "team_first_pick_position":
        m = re.match(r"^(.+?_first_pick_pos)_", market_id)
        return m.group(1) if m else None
    if market_group.startswith("nth_at_position_"):
        nth = market_group[len("nth_at_position_"):]
        m = re.match(rf"^({re.escape(nth)}_[a-z]+)_", market_id)
        return m.group(1) if m else None
    if market_group == "draft_position_over_under":
        m = re.match(r"^(draft_position_ou_.+)_(?:over|under)$", market_id)
        return m.group(1) if m else None
    return None


def resolve_market_id(book: str, book_label: str, book_subject: str) -> Optional[str]:
    """Look up canonical market_id for a per-book (label, subject) pair."""
    with read_connection() as con:
        result = con.execute(
            "SELECT market_id FROM market_map WHERE book = ? AND book_label = ? AND book_subject = ?",
            [book, book_label, book_subject],
        ).fetchone()
    return result[0] if result else None
