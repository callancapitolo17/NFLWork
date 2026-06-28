from unabated_edge.sports.base import SportAdapter
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
WC_MATCH_SERIES = "KXWCMATCH"          # REPLACE from Task 0 FINDINGS

class Soccer(SportAdapter):
    sport = "soccer"
    league_prefix = "lg21"
    outcomes = ("home", "draw", "away")

    def canon_team(self, name: str) -> str:
        return _ALIASES.get((name or "").strip().lower(), (name or "").strip())

    def kalshi_series(self) -> str:
        return WC_MATCH_SERIES

    @staticmethod
    def _ml_price(st, key):
        """Moneyline price at a line key, or None. Rejects spread lines
        (a true moneyline carries no `points`)."""
        ln = st.lines.get(key)
        if not ln:
            return None
        p = line_american_price(ln)
        if p is None or ln.get("points") is not None:
            return None
        return p

    @staticmethod
    def _draw_price(st, key):
        ln = st.lines.get(key)
        if not ln:
            return None
        return line_american_price(ln)

    def _book_three_way(self, st, eid, ms):
        """home/away/draw prices from a SINGLE book `ms`, or None if that book
        lacks any leg. Anchoring all three legs to one book avoids devigging a
        Frankenstein 3-way stitched from books with different vig/staleness."""
        h = self._ml_price(st, f"{eid}|1|{ms}|bt1")
        a = self._ml_price(st, f"{eid}|0|{ms}|bt1")
        d = self._draw_price(st, f"{eid}|1|{ms}|bt4")   # draw: bt4 (FROM TASK 0 FINDINGS)
        if h is None or a is None or d is None:
            return None
        return h, a, d

    def fair(self, st, ev) -> dict | None:
        triple = None
        for ms in config.ANCHOR_SOURCE_IDS:          # first book with a complete 3-way wins
            triple = self._book_three_way(st, ev.event_id, ms)
            if triple is not None:
                break
        if triple is None:
            return None
        h, a, d = triple
        ph, pd, pa = (pricing.american_to_prob(x) for x in (h, d, a))
        dh, dd, da = pricing.devig([ph, pd, pa])
        return {"home": dh, "draw": dd, "away": da}

    def map_outcome_tickers(self, kalshi_event: dict) -> dict:
        out = {}
        for mk in kalshi_event.get("markets", []):
            raw = (mk.get("yes_sub_title") or "").strip()
            if "draw" in raw.lower() or "tie" in raw.lower():
                out["draw"] = mk["ticker"]
            else:
                out.setdefault("_named", []).append((raw, mk["ticker"]))
        return out   # home/away resolved against canon team names in mapping (Task 7)
