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

    @staticmethod
    def _side_prices(ln) -> dict[float, float]:
        """{line: american_price} from one side's line dict: main quote + alternateLines."""
        out = {}
        if ln is None:
            return out
        px, pts = line_american_price(ln), ln.get("points")
        if px is not None and pts is not None:
            out[float(pts)] = px
        for al in ln.get("alternateLines") or []:
            apx, apts = line_american_price(al), al.get("points")
            if apx is not None and apts is not None:
                out[float(apts)] = apx
        return out

    def _anchor_totals(self, state, eid) -> dict[float, float]:
        """Devigged total ladder {line: p_over} from the FIRST anchor book that has
        at least one complete rung (over+under at the same line). Same-book only —
        never mixes books' vig. Empty dict when no anchor book is complete."""
        for ms in config.ANCHOR_SOURCE_IDS:
            overs = self._side_prices(state.lines.get(f"{eid}|0|{ms}|bt3"))   # side 0 = Over
            unders = self._side_prices(state.lines.get(f"{eid}|1|{ms}|bt3"))  # side 1 = Under
            ladder = {}
            for line in sorted(set(overs) & set(unders)):
                p_over, _ = pricing.devig([pricing.american_to_prob(overs[line]),
                                           pricing.american_to_prob(unders[line])])
                ladder[line] = p_over
            if ladder:
                return ladder
        return {}

    def price_event(self, state, event_meta, kalshi_event) -> list[Candidate]:
        """Price Over/Under candidates at every Kalshi rung the anchor ladder quotes.

        Fail closed: [] when the anchor is missing; rungs without an anchor line
        are skipped silently (never interpolated)."""
        ladder = self._anchor_totals(state, event_meta.event_id)
        if not ladder:
            return []
        out = []
        for mk in kalshi_event.get("markets", []):
            if mk.get("strike_type") != "greater":
                continue
            try:
                line = float(mk.get("floor_strike"))
            except (TypeError, ValueError):
                continue
            p_over = ladder.get(line)
            if p_over is None:
                continue
            out.append(Candidate(mk["ticker"], "yes", p_over, f"over_{line}"))
            out.append(Candidate(mk["ticker"], "no", round(1 - p_over, 6), f"under_{line}"))
        return out
