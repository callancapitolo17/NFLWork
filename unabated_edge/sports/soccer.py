from unabated_edge.sports.totals import TotalsLadderAdapter, _older_mo

# Maps known divergent spellings -> a single plain-English canonical so that
# Unabated and Kalshi names collapse to the same string. Each entry is a genuine
# synonym (no two distinct nations collide). Expand/verify against live data in
# Task 0 — silent unmatched events are now logged in mapping.pair_events to surface gaps.
_ALIASES = {
    "korea republic": "South Korea", "republic of korea": "South Korea",
    "korea dpr": "North Korea", "dpr korea": "North Korea",
    "usa": "United States", "united states of america": "United States",
    "ir iran": "Iran", "iran islamic republic": "Iran",
    "china pr": "China",
    "côte d'ivoire": "Ivory Coast", "cote d'ivoire": "Ivory Coast",
    "czechia": "Czech Republic",
    "türkiye": "Turkey", "turkiye": "Turkey",
    "bosnia and herzegovina": "Bosnia", "bosnia & herzegovina": "Bosnia",
    "republic of ireland": "Ireland",
    "cabo verde": "Cape Verde",
    "curaçao": "Curacao",
}

WC_TOTAL_SERIES = "KXWCTOTAL"   # verified live: Over-ladder markets for WC reg-time total goals


class Soccer(TotalsLadderAdapter):
    sport = "soccer"
    league_prefix = "lg21"

    def canon_team(self, name: str) -> str:
        return _ALIASES.get((name or "").strip().lower(), (name or "").strip())

    def kalshi_series(self) -> str:
        return WC_TOTAL_SERIES

    def event_teams(self, kalshi_event: dict) -> frozenset:
        """Parse the two team names from the Kalshi event title.

        Title format: "Colombia vs Ghana: Regulation Time Total Goals"
        → take everything before the first ":", split on " vs ".
        """
        title = (kalshi_event.get("title") or "")
        before_colon = title.split(":")[0]
        parts = before_colon.split(" vs ")
        if len(parts) != 2:
            return frozenset()
        return frozenset({self.canon_team(parts[0].strip()), self.canon_team(parts[1].strip())})
