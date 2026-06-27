from unabated_edge.sports.base import SportAdapter
from unabated_edge import pricing, config
from unabated_edge.feed import line_american_price

_ALIASES = {"korea republic":"South Korea","usa":"United States","ir iran":"Iran"}
WC_MATCH_SERIES = "KXWCMATCH"          # REPLACE from Task 0 FINDINGS

class Soccer(SportAdapter):
    sport = "soccer"
    league_prefix = "lg21"
    outcomes = ("home", "draw", "away")

    def canon_team(self, name: str) -> str:
        return _ALIASES.get((name or "").strip().lower(), (name or "").strip())

    def kalshi_series(self) -> str:
        return WC_MATCH_SERIES

    def _anchor_ml(self, st, eid, side):
        for ms in config.ANCHOR_SOURCE_IDS:
            ln = st.lines.get(f"{eid}|{side}|{ms}|bt1")
            if ln and line_american_price(ln) is not None and ln.get("points") is None:
                return line_american_price(ln)
        return None

    def _draw(self, st, eid):
        # FROM TASK 0 FINDINGS — e.g. draw lives in bt4 on side "1":
        ln = st.lines.get(f"{eid}|1|{config.SHARP_BOOK_PRICE_ID}|bt4")
        return line_american_price(ln) if ln and line_american_price(ln) is not None else None

    def fair(self, st, ev) -> dict | None:
        h, a, d = self._anchor_ml(st, ev.event_id, "1"), self._anchor_ml(st, ev.event_id, "0"), self._draw(st, ev.event_id)
        if h is None or a is None or d is None:
            return None
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
