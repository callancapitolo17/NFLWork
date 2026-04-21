"""Canonical 2026 NFL Draft prospect registry.

Each entry: canonical_name (the spelling we'll use everywhere) ->
metadata + aliases (every variant a book might use).

Maintenance: when a scraper produces an unmapped player (lands in
draft_odds_unmapped_players), add the alias here and re-run
`python -m nfl_draft.lib.seed`.
"""

PLAYERS: dict[str, dict] = {
    # canonical_name: {position, college, aliases: [str]}
    "Cam Ward": {
        "position": "QB",
        "college": "Miami",
        "aliases": ["Cam Ward", "Cameron Ward", "C. Ward", "Ward, Cam", "Ward (Miami)"],
    },
    "Shedeur Sanders": {
        "position": "QB",
        "college": "Colorado",
        "aliases": ["Shedeur Sanders", "S. Sanders", "Sanders, Shedeur", "Sanders (Colorado)"],
    },
    "Travis Hunter": {
        "position": "WR",  # also DB; use primary draft position
        "college": "Colorado",
        "aliases": ["Travis Hunter", "T. Hunter", "Hunter, Travis"],
    },
    "Ashton Jeanty": {
        "position": "RB",
        "college": "Boise State",
        "aliases": ["Ashton Jeanty", "A. Jeanty", "Jeanty, Ashton", "Jeanty (Boise St)"],
    },
    # Top-80 to be populated during scraper reconnaissance as unmapped
    # players surface in draft_odds_unmapped_players. Keep the seed lean
    # for now - the quarantine table is the source of truth for "what do
    # we still need to map".
}
