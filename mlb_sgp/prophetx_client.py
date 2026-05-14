"""ProphetX HTTP client — extracted from scraper_prophetx_sgp.py.

Owns the curl_cffi Chrome-TLS session + profile-cookie auth fallback (PX
stores cookies in `.prophetx_profile/Default/Cookies`).

Exposes three methods, mirroring dk_client.py / fd_client.py:
  - list_events()         — MLB events today
  - fetch_event_markets() — markets + selections per event
  - submit_parlay_rfq()   — RFQ parlay pricer with offer-ladder filtering

The pure helper functions `_parse_events_response`, `_parse_event_markets`,
and `_pick_offer` are exposed at module level so tests can exercise the
parsers/policies without a live session.

Real ProphetX response shape (captured 2026-05-13 from live API):

    GET /trade/public/api/v1/tournaments?expand=events&type=highlight&limit=150
      -> {"data": {"tournaments": [
            {"id": ..., "name": "MLB", "sportEvents": [
                {"id": ..., "scheduled": "...", "sport": {"name": "Baseball"},
                 "competitors": [
                     {"id": ..., "name": "...", "seq": 0|1, ...},
                     ...
                 ]},
                ...
            ]},
            ...
         ]}}

    GET /trade/public/api/v2/events/{event_id}/markets
      -> {"data": {"markets": [
            {"id": ..., "name": "Run Line",
             "marketLines": [
                {"id": ..., "outcomes": [
                    {"id": ..., "line": -5.5, "lineID": "...", ...},
                    ...
                ]},
                ...
             ]},
            ...
         ]}}

The legacy `fetch_prophetx_mlb_events` already exists in scraper_prophetx_sgp
and returns a FLATTENED list of event dicts. We deliberately parse the RAW
nested response here so tests can pin shape-level behaviour without depending
on the legacy flattening helper.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


PROPHETX_BASE = "https://www.prophetx.co"
DEFAULT_MIN_OFFER_STAKE = 150


@dataclass
class Event:
    event_id: str
    home_team: str
    away_team: str
    home_id: int | None
    away_id: int | None
    start_time: str   # ISO UTC string


@dataclass
class Market:
    market_id: str
    name: str
    # ProphetX returns spread/total strike inside marketLines, not at the
    # top level. We keep `strike` for API symmetry with the spec but it will
    # often be None — callers should read per-line strikes off `selections`.
    strike: float | None
    # Raw marketLines passthrough. Each entry has the shape:
    #   {"id": ..., "name": ..., "outcomes": [{"id", "line", "lineID", ...}], ...}
    selections: list[dict]
    # Raw top-level `outcomes` passthrough. PX's Moneyline market carries a
    # flat `outcomes` list at the market level (with `competitorId` per side)
    # alongside `marketLines`. The legacy `_verify_competitor_ids` helper
    # reads this list to confirm home/away IDs match the tournaments endpoint.
    # Often empty for non-moneyline markets; the helper is permissive.
    outcomes: list[dict] = field(default_factory=list)


@dataclass
class SelectionLeg:
    """Minimal info needed to POST one leg in an RFQ payload."""
    sport_event_id: str
    market_id: str
    outcome_id: str
    line_id: str
    line: float


class ProphetXClient:
    """Thin wrapper around the three ProphetX endpoints SGP scrapers need."""

    def __init__(self, verbose: bool = False) -> None:
        # Reuse the legacy scraper's session bootstrap to avoid duplicating
        # the Chrome-TLS impersonation + cookie-loading logic.
        from scraper_prophetx_sgp import init_session
        self.session = init_session()
        self.verbose = verbose

    def list_events(self) -> list[Event]:
        """Fetch all events with `expand=events`, return MLB-only Event list."""
        url = f"{PROPHETX_BASE}/trade/public/api/v1/tournaments"
        params = {"expand": "events", "type": "highlight", "limit": 150}
        try:
            r = self.session.get(url, params=params, timeout=15)
        except TypeError:
            # FakeSession in unit tests doesn't accept params/timeout kwargs
            r = self.session.get(url)
        if getattr(r, "status_code", 200) != 200:
            return []
        return _parse_events_response(r.json())

    def fetch_event_markets(self, event_id: str) -> list[Market]:
        """Fetch the full market tree for one event."""
        url = f"{PROPHETX_BASE}/trade/public/api/v2/events/{event_id}/markets"
        try:
            r = self.session.get(url, timeout=15)
        except TypeError:
            r = self.session.get(url)
        if getattr(r, "status_code", 200) != 200:
            return []
        return _parse_event_markets(r.json())

    def submit_parlay_rfq(
        self,
        legs: list[SelectionLeg],
        stake: float = 1.0,
        min_offer_stake: int = DEFAULT_MIN_OFFER_STAKE,
    ) -> tuple[dict | None, bool]:
        """POST an RFQ and return (chosen_offer, used_fallback).

        `chosen_offer` is None for empty responses (RFQ rejected or no offers).
        `used_fallback` is True iff we had to drop below `min_offer_stake`
        (thin market). This matches the legacy `submit_parlay_rfq` policy in
        scraper_prophetx_sgp: ProphetX's offer ladder interleaves $50 teaser
        tiers with real-liquidity tiers; the UI hides the teasers, so we do
        too.
        """
        url = f"{PROPHETX_BASE}/parlay/public/api/v1/user/request"
        # ProphetX expects sportEventId / marketId / outcomeId as int64 on
        # the wire. The SelectionLeg dataclass stores them as strings (so
        # tests and mocks don't have to fuss with types), so we coerce
        # back to int here. lineId is a hex hash and stays a string.
        payload = {
            "marketLines": [
                {"sportEventId": _to_int_or_str(l.sport_event_id),
                 "marketId": _to_int_or_str(l.market_id),
                 "outcomeId": _to_int_or_str(l.outcome_id),
                 "lineId": l.line_id, "line": l.line}
                for l in legs
            ],
            "stake": stake,
        }
        try:
            r = self.session.post(
                url, json=payload,
                headers={"Content-Type": "application/json"},
                timeout=10,
            )
        except TypeError:
            r = self.session.post(url, json=payload)

        if getattr(r, "status_code", 200) != 200:
            return None, False
        try:
            data = r.json()
        except Exception:
            return None, False
        # ProphetX wraps offers under data.offers when successful
        offers = (data.get("data") or {}).get("offers") or data.get("offers") or []
        picked = _pick_offer(offers, min_stake=min_offer_stake)
        if picked is None:
            return None, False
        used_fallback = (picked.get("stake") or 0) < min_offer_stake
        return picked, used_fallback


# ---------------------------------------------------------------------------
# Pure parser / policy helpers — exposed at module level for tests
# ---------------------------------------------------------------------------

def _to_int_or_str(val):
    """Coerce a stringified numeric ID back to int for the PX RFQ wire format.

    PX rejects string IDs with `cannot unmarshal string into ... int64`. The
    SelectionLeg dataclass stores IDs as strings (test ergonomics), so we
    parse to int at the network boundary. If the value isn't a parseable
    int (e.g. an empty string from a missing field), we pass it through
    unchanged — the API will return a 500 the caller already handles.
    """
    if isinstance(val, bool):
        return val  # don't let True/False slip through int()
    if isinstance(val, int):
        return val
    if isinstance(val, str) and val.strip().lstrip("-").isdigit():
        return int(val)
    return val

def _parse_events_response(raw: dict) -> list[Event]:
    """Filter a raw `tournaments` response to MLB events.

    Matches both response shapes seen in the wild:
      - {"data": {"tournaments": [...]}}  (current PX API, captured 2026-05-13)
      - {"tournaments": [...]}            (legacy/synthetic test fixtures)

    Per-tournament events live under either `sportEvents` (current API) or
    `events` (older shape). Sport is on the EVENT, not the tournament. We
    require BOTH `sport.name == "Baseball"` AND the tournament name to
    contain "MLB" (the bare substring "Major League" alone would also match
    "Major League Soccer", so we don't accept it on its own — only as
    "Major League Baseball").
    """
    tournaments = (raw.get("data", {}) or {}).get("tournaments")
    if tournaments is None:
        tournaments = raw.get("tournaments", []) or []

    out: list[Event] = []
    for tourn in tournaments:
        tname = tourn.get("name", "") or ""
        is_mlb_named = ("MLB" in tname) or ("Major League Baseball" in tname)
        if not is_mlb_named:
            continue
        events_list = tourn.get("sportEvents") or tourn.get("events") or []
        for ev in events_list:
            sport_name = (ev.get("sport") or {}).get("name", "")
            if sport_name and sport_name != "Baseball":
                continue
            competitors = ev.get("competitors", []) or []
            if len(competitors) < 2:
                continue
            # ProphetX uses seq=0 for home, seq=1 for away (verified in
            # recon notes). Fall back to isHome flag, then to positional.
            by_seq = {c.get("seq"): c for c in competitors if c.get("seq") is not None}
            home = by_seq.get(0)
            away = by_seq.get(1)
            if home is None or away is None:
                home = next((c for c in competitors if c.get("isHome")), None) or competitors[0]
                away = next((c for c in competitors if not c.get("isHome")), None) or competitors[1]
            out.append(Event(
                event_id=str(ev.get("id", "")),
                home_team=home.get("name", "") or "",
                away_team=away.get("name", "") or "",
                home_id=home.get("id"),
                away_id=away.get("id"),
                start_time=ev.get("scheduled", "") or "",
            ))
    return out


def _parse_event_markets(raw: dict) -> list[Market]:
    """Pull markets out of an event-markets response.

    Tolerates both `{"data": {"markets": [...]}}` and `{"markets": [...]}`
    shapes. `marketLines` is passed through unmodified so callers (or
    leg-resolver helpers) can pick specific outcomes by line value, line ID,
    or competitor ID — same pattern as scraper_prophetx_sgp's _pick_selection.
    """
    markets = (raw.get("data", {}) or {}).get("markets")
    if markets is None:
        markets = raw.get("markets", []) or []

    out: list[Market] = []
    for m in markets:
        out.append(Market(
            market_id=str(m.get("id", "")),
            name=m.get("name", "") or "",
            strike=m.get("strike"),
            selections=list(m.get("marketLines") or []),
            outcomes=list(m.get("outcomes") or []),
        ))
    return out


def _pick_offer(offers: list[dict], min_stake: int) -> dict | None:
    """Pick the first offer with `stake >= min_stake`.

    Falls back to `offers[0]` when nothing clears the threshold (thin
    market — better to quote the teaser than to drop the combo entirely).
    Returns None for an empty offer list. This mirrors the legacy logic in
    `scraper_prophetx_sgp.submit_parlay_rfq` (memory note "ProphetX SGP
    Scraping": MIN_OFFER_STAKE filter + offers[0] teaser).
    """
    if not offers:
        return None
    for o in offers:
        if (o.get("stake") or 0) >= min_stake:
            return o
    return offers[0]
