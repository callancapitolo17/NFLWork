"""Resolve player/team aliases to canonical names via DuckDB lookup tables."""

from typing import Optional
from nfl_draft.lib.db import read_connection


def resolve_player(name: str) -> Optional[str]:
    """Returns canonical player name, or None if unmapped."""
    if not name:
        return None
    key = name.lower().strip()
    with read_connection() as con:
        result = con.execute(
            "SELECT canonical_name FROM player_aliases WHERE alias = ?",
            [key],
        ).fetchone()
    return result[0] if result else None


def resolve_team(name: str) -> Optional[str]:
    """Returns canonical team abbr (e.g., 'CHI'), or None if unmapped."""
    if not name:
        return None
    key = name.lower().strip()
    with read_connection() as con:
        result = con.execute(
            "SELECT canonical_abbr FROM team_aliases WHERE alias = ?",
            [key],
        ).fetchone()
    return result[0] if result else None
