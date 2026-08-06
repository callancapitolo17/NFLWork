"""Shared pricing for sports whose Kalshi market is an Over-ladder on a single
final total (goals, runs, points), anchored by Unabated bt3 over/under.

Subclasses provide only identity: sport, league_prefix, canon_team,
kalshi_series, event_teams. All ladder devigging and rung matching lives here."""
import logging

from unabated_edge.sports.base import SportAdapter, Candidate
from unabated_edge import pricing, config
from unabated_edge.feed import line_american_price

log = logging.getLogger("unabated_edge")


def _older_mo(a, b):
    """Older (lexicographically smaller, since ISO-8601 sorts by time) of two
    modifiedOn strings; ignores None."""
    vals = [x for x in (a, b) if x]
    return min(vals) if vals else None


class TotalsLadderAdapter(SportAdapter):
    @staticmethod
    def _side_prices(ln) -> dict[float, tuple[float, bool]]:
        """{line: (american_price, is_alt)} from one side's line dict. The main
        quote wins over an alternateLines entry at the same points value."""
        out = {}
        if ln is None:
            return out
        px, pts = line_american_price(ln), ln.get("points")
        if px is not None and pts is not None:
            out[float(pts)] = (px, False)
        for al in ln.get("alternateLines") or []:
            apx, apts = line_american_price(al), al.get("points")
            if apx is not None and apts is not None:
                out.setdefault(float(apts), (apx, True))
        return out

    def _anchor_ladder(self, state, eid) -> dict[float, dict]:
        """Devigged total ladder {line: {p_over, book, alt, overround}} from the
        FIRST anchor book that has at least one complete rung (over+under at the
        same line). Same-book only — never mixes books' vig. `alt` marks rungs
        built from alternateLines (possibly Unabated-derived, not a raw book
        quote) and `overround` is the pre-devig implied sum — both go to the
        research firehose so alt provenance stays auditable. Each rung also
        carries `modified_on` — the older of the two sides' feed modifiedOn —
        for the staleness gate.

        Overround/crossed gate (Finding #73): a rung is dropped BEFORE devig
        if its raw implied sum falls outside the applicable band (main rungs:
        config.MAIN_OVERROUND_MIN/MAX; alt rungs: config.ALT_OVERROUND_MIN/MAX).
        A crossed pair (sum < 1.0) is a data error, not a real two-way quote —
        pricing.devig would otherwise silently renormalize it into a
        plausible-but-wrong p_over. This is the layer that keeps a blown/
        crossed rung from ever becoming a candidate; engine.py's alt-only
        gate (_desired) is a second, redundant layer for alt rungs specifically."""
        for idx, ms in enumerate(config.ANCHOR_SOURCE_IDS):
            over_ln = state.lines.get(f"{eid}|0|{ms}|bt3")
            under_ln = state.lines.get(f"{eid}|1|{ms}|bt3")
            overs = self._side_prices(over_ln)
            unders = self._side_prices(under_ln)
            mo = _older_mo((over_ln or {}).get("modifiedOn"), (under_ln or {}).get("modifiedOn"))
            ladder = {}
            for line in sorted(set(overs) & set(unders)):
                (opx, o_alt), (upx, u_alt) = overs[line], unders[line]
                alt = o_alt or u_alt
                po_raw, pu_raw = pricing.american_to_prob(opx), pricing.american_to_prob(upx)
                overround = round(po_raw + pu_raw, 6)
                lo, hi = (config.ALT_OVERROUND_MIN, config.ALT_OVERROUND_MAX) if alt \
                    else (config.MAIN_OVERROUND_MIN, config.MAIN_OVERROUND_MAX)
                if not (lo <= overround <= hi):
                    log.warning(
                        "dropping %s rung: eid=%s line=%s book=%s overround=%.4f "
                        "outside band [%.2f, %.2f]",
                        "alt" if alt else "main", eid, line, ms, overround, lo, hi)
                    continue
                p_over, _ = pricing.devig([po_raw, pu_raw])
                ladder[line] = {"p_over": p_over, "book": ms, "alt": alt,
                                "overround": overround, "modified_on": mo}
            if ladder:
                if idx > 0:
                    log.warning(
                        "anchor failover: eid=%s using book %s (primary book %s "
                        "unavailable/incomplete)", eid, ms, config.ANCHOR_SOURCE_IDS[0])
                return ladder
        return {}

    def _anchor_totals(self, state, eid) -> dict[float, float]:
        """{line: p_over} view of the ladder (see _anchor_ladder)."""
        return {line: r["p_over"] for line, r in self._anchor_ladder(state, eid).items()}

    def fair_ladder(self, state, event_meta):
        return self._anchor_ladder(state, event_meta.event_id) or None

    def price_event(self, state, event_meta, kalshi_event) -> list[Candidate]:
        """Price Over/Under candidates at every Kalshi rung the anchor ladder quotes.

        Fail closed: [] when the anchor is missing; rungs without an anchor line
        are skipped silently (never interpolated)."""
        ladder = self._anchor_ladder(state, event_meta.event_id)
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
            rung = ladder.get(line)
            if rung is None:
                continue
            meta = {"book": rung["book"], "alt": rung["alt"], "overround": rung["overround"]}
            p_over = rung["p_over"]
            out.append(Candidate(mk["ticker"], "yes", p_over, f"over_{line}", meta))
            out.append(Candidate(mk["ticker"], "no", round(1 - p_over, 6), f"under_{line}", meta))
        return out
