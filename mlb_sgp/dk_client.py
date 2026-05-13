"""DraftKings API client — extracted from scraper_draftkings_sgp.py.

Owns the curl_cffi Chrome-TLS session and exposes the three operations
both the SGP scraper and the singles scraper need:
  - list_events()           — all MLB events today
  - fetch_event_markets()   — market metadata (FG vs F5 vs alt)
  - fetch_event_selections()— all selections with prices, in one call
"""
from __future__ import annotations
from dataclasses import dataclass


@dataclass
class Event:
    event_id: str
    home_team: str
    away_team: str
    start_time: str  # ISO UTC string


@dataclass
class Market:
    market_id: str
    name: str
    subcategory: str  # "game_lines" | "alt_lines" | "innings" etc.


@dataclass
class Selection:
    selection_id: str
    market_id: str
    name: str            # e.g. "Yankees -1.5"
    line: float | None   # spread/total value; None for moneylines
    american_odds: int


class DraftKingsClient:
    def __init__(self, verbose: bool = False) -> None:
        from scraper_draftkings_sgp import init_session
        self.session = init_session()
        self.verbose = verbose

    def list_events(self) -> list[Event]:
        """Returns one Event per MLB game today.

        Thin wrapper over the legacy fetch_dk_events function. The SGP refactor
        in Task 6 keeps fetch_dk_events in place; this method just adapts its
        dict output to the Event dataclass.

        fetch_dk_events returns dicts with keys:
          dk_event_id, dk_home, dk_away, start_time, name
        """
        from scraper_draftkings_sgp import fetch_dk_events
        raw = fetch_dk_events(self.session)
        return [Event(
            event_id=str(e["dk_event_id"]),
            home_team=e["dk_home"],
            away_team=e["dk_away"],
            start_time=e.get("start_time", ""),
        ) for e in raw]

    def fetch_event_markets(self, event_id: str) -> list[Market]:
        """Returns market metadata for one event (game lines + alt lines).

        Hits the DK subcategory endpoint twice:
          - subcat 4519 -> game lines (Moneyline, Run Line, Total)
          - subcat 15628 -> alt lines (alt RL, alt totals)

        We DO NOT get prices here — only (market_id, name) tuples. Prices live
        in the parlays endpoint, fetched by fetch_event_selections.
        """
        from scraper_draftkings_sgp import _fetch_subcat_markets
        results: list[Market] = []
        for subcat_id, subcat_label in (
            ("4519", "game_lines"),
            ("15628", "alt_lines"),
        ):
            for m_id, m_name in _fetch_subcat_markets(self.session, event_id, subcat_id):
                results.append(Market(
                    market_id=str(m_id),
                    name=str(m_name or ""),
                    subcategory=subcat_label,
                ))
        return results

    def fetch_event_selections(self, event_id: str) -> list[Selection]:
        """Returns all selections (with prices) for one event from the parlays endpoint.

        DK's parlays/v1/sgp/events/{event_id} returns a payload of shape:
          data.markets[] -> { id, name, selections[] -> {
              id, marketId, name, outcomeType, points (float|null),
              displayOdds.american (str like "+260" or Unicode-minus "−370"),
              isDisabled, status }}

        Selections that are disabled or missing displayOdds.american are skipped
        — the singles scraper can't price an unavailable bet.
        """
        from scraper_draftkings_sgp import DK_SGP_PARLAYS_URL
        resp = self.session.get(
            f"{DK_SGP_PARLAYS_URL}/{event_id}",
            timeout=60,
        )
        if getattr(resp, "status_code", 200) != 200:
            return []
        payload = resp.json()
        markets = (payload.get("data") or {}).get("markets") or []
        out: list[Selection] = []
        for mkt in markets:
            market_id = str(mkt.get("id", ""))
            for sel in mkt.get("selections", []) or []:
                if sel.get("isDisabled"):
                    continue
                disp = sel.get("displayOdds") or {}
                american_str = disp.get("american")
                if american_str in (None, ""):
                    continue
                try:
                    american = _parse_american_odds(american_str)
                except (ValueError, TypeError):
                    continue
                points = sel.get("points")
                line: float | None
                if points is None:
                    line = None
                else:
                    try:
                        line = float(points)
                    except (TypeError, ValueError):
                        line = None
                out.append(Selection(
                    selection_id=str(sel.get("id", "")),
                    market_id=str(sel.get("marketId", market_id)),
                    name=str(sel.get("name", "") or sel.get("outcomeType", "")),
                    line=line,
                    american_odds=american,
                ))
        return out


def _parse_american_odds(value) -> int:
    """Parse DK's american-odds string (e.g. '+260', '−370', '-110') to int.

    DK uses the Unicode minus sign U+2212 ('−') for negatives in displayOdds —
    not the ASCII hyphen. Both must be handled.
    """
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return int(value)
    s = str(value).strip()
    # Normalize: strip leading '+', convert Unicode minus to ASCII '-'
    s = s.replace("−", "-")
    if s.startswith("+"):
        s = s[1:]
    return int(s)
