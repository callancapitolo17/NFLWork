from unabated_edge.sports.base import SportAdapter, Candidate
from unabated_edge import pricing, config
from unabated_edge.feed import line_american_price

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


class Soccer(SportAdapter):
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

    def _anchor_total(self, state, eid) -> tuple[float, float] | None:
        """Find the first anchor book that quotes both Over and Under for bt3.

        Returns (line, p_over) where p_over is the devigged probability, or None
        if no anchor book has a complete, consistent total quote.
        """
        for ms in config.ANCHOR_SOURCE_IDS:
            over = state.lines.get(f"{eid}|0|{ms}|bt3")   # side 0 = Over
            under = state.lines.get(f"{eid}|1|{ms}|bt3")  # side 1 = Under
            if over is None or under is None:
                continue
            po = line_american_price(over)
            pu = line_american_price(under)
            if po is None or pu is None:
                continue
            line = over.get("points")
            if line is None or under.get("points") != line:
                continue
            p_over, _ = pricing.devig([pricing.american_to_prob(po), pricing.american_to_prob(pu)])
            return (float(line), p_over)
        return None

    def price_event(self, state, event_meta, kalshi_event) -> list[Candidate]:
        """Price Over and Under candidates from the anchor total vs the Kalshi rung.

        Fail closed: return [] when anchor is missing or no Kalshi market matches the line.
        """
        result = self._anchor_total(state, event_meta.event_id)
        if result is None:
            return []
        line, p_over = result
        markets = kalshi_event.get("markets", [])
        mk = next(
            (m for m in markets
             if m.get("strike_type") == "greater" and float(m.get("floor_strike", -1)) == line),
            None,
        )
        if mk is None:
            return []
        ticker = mk["ticker"]
        return [
            Candidate(ticker, "yes", p_over, f"over_{line}"),
            Candidate(ticker, "no", round(1 - p_over, 6), f"under_{line}"),
        ]
