"""Novig HTTP client — extracted from scraper_novig_sgp.py.

Anonymous Hasura GraphQL endpoint (no auth required — the parlay endpoint
is literally named /unauthenticated). Novig sources legs from DraftKings
(every observed leg returns vendor="DRAFTKINGS"), so its line set is a
strict subset of DK's. We do strict line matching here at the client
level — interpolation/fallback logic lives in the orchestrator, not the
HTTP layer.

Exposes three thin methods, mirroring dk_client.py / fd_client.py /
prophetx_client.py:
  - list_events()       -> list[Event]
  - fetch_event_legs()  -> EventLegs (spread + total outcomes)
  - submit_parlay()     -> dict (parsed parlay node)

Pure helpers `_parse_events_response` and `_parse_event_legs_response`
are at module level so tests can exercise the parsers without a session.

Real Novig response shape (captured 2026-05-13 from live API):

    POST https://api.novig.us/v1/graphql  (MLBEvents query)
      -> {"data": {"event": [
            {"id": "<uuid>",
             "scheduled_start": "2026-05-13T23:40:00+00:00",
             "game": {
                "homeTeam": {"name": "...", "symbol": "MIN", "short_name": "MIN"},
                "awayTeam": {"name": "...", "symbol": "MIA", "short_name": "MIA"}
             }},
            ...
         ]}}

    POST https://api.novig.us/v1/graphql  (EventMarkets_Query)
      -> {"data": {"event": [
            {"id": "<uuid>",
             "markets": [
                {"id": ..., "type": "SPREAD"|"TOTAL"|"SPREAD_1H"|"TOTAL_1H"|...,
                 "strike": -1.5, "is_consensus": true|false,
                 "outcomes": [
                    {"id": "<uuid>", "description": "Over 8.5"|"MIA +2.5",
                     "available": 0.526, "competitor": {"symbol": "MIA", ...}, ...},
                    ...
                 ]},
                ...
             ]},
         ]}}

    POST https://api.novig.us/nbx/v1/parlay/request/unauthenticated
      -> [{"price": "0.35088", "status": "OPEN", ...}, ...]
       (a *list* of offer dicts — top offer at [0])

The parsers below tolerate both the real `{"data": {"event": [...]}}` shape
and a flat `{"event": [...]}` / `{"events": [...]}` fallback so synthetic
fixtures stay simple.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any


NOVIG_GRAPHQL = "https://api.novig.us/v1/graphql"
NOVIG_PARLAY = "https://api.novig.us/nbx/v1/parlay/request/unauthenticated"


@dataclass
class Event:
    event_id: str
    home_team: str
    away_team: str
    home_sym: str
    away_sym: str
    start_time: str  # ISO UTC string


@dataclass
class EventLegs:
    """All FG spread + total selections (outcomes) for one event.

    Each leg dict carries at minimum `id` (outcome UUID for the parlay POST),
    `description`, `available` (implied probability from the API; used by
    callers for naive-multiply sanity filtering), `strike`, and `period`
    (one of "fg" / "f5"). For spreads, `side` is "home" or "away"; for
    totals, `side` is "over" or "under".
    """
    event_id: str
    spread_legs: list[dict] = field(default_factory=list)
    total_legs: list[dict] = field(default_factory=list)


# Market type names — Novig uses *_1H suffix for first-five-innings.
SPREAD_TYPES = {"SPREAD", "SPREAD_1H"}
TOTAL_TYPES = {"TOTAL", "TOTAL_1H"}


class NovigClient:
    """Thin wrapper around Novig's anonymous GraphQL + parlay RFQ endpoints."""

    def __init__(self, verbose: bool = False) -> None:
        # Reuse the legacy scraper's session bootstrap (curl_cffi Chrome
        # impersonation + landing-page warm-up for Cloudflare cookies).
        from scraper_novig_sgp import init_session
        self.session = init_session()
        self.verbose = verbose

    def list_events(self) -> list[Event]:
        """List upcoming MLB events. Uses the same query the legacy scraper does."""
        from datetime import datetime, timezone, timedelta
        from scraper_novig_sgp import MLB_EVENTS_QUERY, EVENT_WINDOW_HOURS

        now = datetime.now(timezone.utc)
        cutoff = (now + timedelta(hours=EVENT_WINDOW_HOURS)).isoformat()
        body = json.dumps({
            "query": MLB_EVENTS_QUERY,
            "variables": {"start_gte": now.isoformat(), "start_lte": cutoff},
        })
        try:
            r = self.session.post(
                NOVIG_GRAPHQL, data=body,
                headers={"Content-Type": "application/json"},
                timeout=20,
            )
        except TypeError:
            # FakeSession in unit tests may not accept all kwargs
            r = self.session.post(NOVIG_GRAPHQL, data=body)
        if getattr(r, "status_code", 200) != 200:
            return []
        try:
            data = r.json()
        except Exception:
            return []
        return _parse_events_response(data)

    def fetch_event_legs(self, event_id: str) -> EventLegs:
        """Fetch the market tree for one event and extract spread + total outcomes."""
        from scraper_novig_sgp import _load_event_markets_query

        query_text = _load_event_markets_query()
        q_obj = json.loads(query_text)
        q_obj["variables"]["eventId"] = event_id
        try:
            r = self.session.post(
                NOVIG_GRAPHQL, data=json.dumps(q_obj),
                headers={"Content-Type": "application/json"},
                timeout=20,
            )
        except TypeError:
            r = self.session.post(NOVIG_GRAPHQL, data=json.dumps(q_obj))
        if getattr(r, "status_code", 200) != 200:
            return EventLegs(event_id=event_id)
        try:
            data = r.json()
        except Exception:
            return EventLegs(event_id=event_id)
        return _parse_event_legs_response(data, event_id_fallback=event_id)

    def submit_parlay(self, outcome_ids: list[str], stake: float = 1.0) -> dict:
        """Submit a BuildParlay request and return the parsed top offer node.

        Returns a dict shaped like:
            {"decimal": 2.85, "american": -200, "price_str": "0.35088",
             "status": "OPEN", "raw_offers": [...]}
        or `{}` for any failure (HTTP non-2xx, empty response, malformed price).

        We use `submit_parlay` for the method name (not `submit_parlay_rfq`)
        to match the task spec, but the underlying endpoint is the same
        anonymous RFQ-style endpoint the legacy scraper hits. `stake` is
        accepted for API symmetry with the spec but Novig's anonymous
        endpoint doesn't gate by stake — it returns offers regardless.
        """
        payload = {"outcomes": [{"id": oid} for oid in outcome_ids], "boostId": None}
        try:
            r = self.session.post(
                NOVIG_PARLAY, json=payload,
                headers={"Content-Type": "application/json"},
                timeout=15,
            )
        except TypeError:
            r = self.session.post(NOVIG_PARLAY, json=payload)
        if getattr(r, "status_code", 200) not in (200, 201):
            return {}
        try:
            offers = r.json()
        except Exception:
            return {}
        return _parse_parlay_response(offers)


# ---------------------------------------------------------------------------
# Pure parser helpers — exposed at module level for tests
# ---------------------------------------------------------------------------

def _parse_events_response(raw: dict) -> list[Event]:
    """Parse a Novig MLBEvents response into Event dataclasses.

    Tolerates both shapes:
      - {"data": {"event": [...]}}   (real Hasura response)
      - {"event":  [...]}            (synthetic fixture fallback)
      - {"events": [...]}            (alternate fallback)
    """
    events_raw = (raw.get("data") or {}).get("event")
    if events_raw is None:
        events_raw = raw.get("event") or raw.get("events") or []

    out: list[Event] = []
    for e in events_raw:
        g = e.get("game") or {}
        ht = g.get("homeTeam") or {}
        at = g.get("awayTeam") or {}
        if not (ht.get("name") and at.get("name")):
            continue
        out.append(Event(
            event_id=str(e.get("id", "")),
            home_team=ht.get("name", "") or "",
            away_team=at.get("name", "") or "",
            home_sym=ht.get("symbol") or ht.get("short_name") or "",
            away_sym=at.get("symbol") or at.get("short_name") or "",
            start_time=e.get("scheduled_start", "") or "",
        ))
    return out


def _parse_event_legs_response(raw: dict, event_id_fallback: str = "") -> EventLegs:
    """Parse a Novig EventMarkets_Query response into spread + total legs.

    Walks `data.event[0].markets[]` (tolerates flat `event[]` / `markets[]`
    fallbacks for synthetic fixtures), keeps only SPREAD/SPREAD_1H and
    TOTAL/TOTAL_1H market types, and emits one leg dict per outcome.

    Each leg dict has:
        {
          "id":          outcome UUID (for parlay POST),
          "description": e.g. "Over 8.5" or "MIA +2.5",
          "available":   implied probability float (or None),
          "strike":      market strike (home-perspective for spreads),
          "period":      "fg" or "f5",
          "side":        "home"/"away" (spread) or "over"/"under" (total),
          "competitor_symbol": competitor symbol for spreads (None for totals),
          "is_consensus": whether market.is_consensus is True,
        }
    """
    # Walk down to the markets list, tolerating multiple shapes.
    ev_arr = (raw.get("data") or {}).get("event")
    if ev_arr is None:
        ev_arr = raw.get("event")
    if ev_arr is None:
        markets = raw.get("markets", []) or []
        ev_id = event_id_fallback
    elif isinstance(ev_arr, list):
        ev = ev_arr[0] if ev_arr else {}
        markets = ev.get("markets") or []
        ev_id = str(ev.get("id") or event_id_fallback)
    else:
        # `event` could be a single dict in a synthetic fixture
        markets = (ev_arr or {}).get("markets") or []
        ev_id = str((ev_arr or {}).get("id") or event_id_fallback)

    out = EventLegs(event_id=ev_id)

    for mkt in markets:
        mtype = mkt.get("type")
        if mtype not in SPREAD_TYPES and mtype not in TOTAL_TYPES:
            continue
        period = "f5" if mtype.endswith("_1H") else "fg"
        strike = mkt.get("strike")
        is_consensus = bool(mkt.get("is_consensus"))
        outcomes = mkt.get("outcomes") or []

        if mtype in SPREAD_TYPES:
            for o in outcomes:
                comp = o.get("competitor") or {}
                sym = comp.get("symbol")
                # Side is determined by competitor symbol vs strike sign — but
                # the scraper resolves it by matching against home_sym/away_sym
                # at the orchestrator level. Here we just emit the leg with
                # the symbol so the caller can disambiguate.
                leg = {
                    "id": o.get("id"),
                    "description": o.get("description") or "",
                    "available": o.get("available"),
                    "strike": strike,
                    "period": period,
                    "side": None,  # caller fills in via home_sym/away_sym
                    "competitor_symbol": sym,
                    "is_consensus": is_consensus,
                }
                out.spread_legs.append(leg)
        else:  # TOTAL
            for o in outcomes:
                desc = (o.get("description") or "").lower()
                if desc.startswith("over "):
                    side = "over"
                elif desc.startswith("under "):
                    side = "under"
                else:
                    side = None
                leg = {
                    "id": o.get("id"),
                    "description": o.get("description") or "",
                    "available": o.get("available"),
                    "strike": strike,
                    "period": period,
                    "side": side,
                    "competitor_symbol": None,
                    "is_consensus": is_consensus,
                }
                out.total_legs.append(leg)

    return out


def _parse_parlay_response(offers: Any) -> dict:
    """Parse a Novig parlay response (list-of-offers) into a single result dict.

    Returns {"decimal": float, "american": int, "price_str": str,
             "status": str, "raw_offers": list} or {} on malformed input.

    Tolerates the rare BuildParlay-mutation shape {"data": {"parlay": {...}}}
    that the task spec mentions, in case Novig ever switches over.
    """
    # Mutation shape: {"data": {"parlay": {...}}}
    if isinstance(offers, dict):
        parlay = (offers.get("data") or {}).get("parlay")
        if parlay:
            dec_raw = parlay.get("decimalOdds")
            try:
                dec = float(dec_raw)
            except (TypeError, ValueError):
                return {}
            am_raw = parlay.get("americanOdds")
            try:
                am = int(am_raw)
            except (TypeError, ValueError):
                am = _decimal_to_american(dec)
            return {
                "decimal": round(dec, 4),
                "american": am,
                "price_str": str(dec_raw),
                "status": parlay.get("status") or "",
                "raw_offers": [parlay],
            }
        return {}

    # Standard shape: a list of offer dicts
    if not offers or not isinstance(offers, list):
        return {}
    top = offers[0]
    if not isinstance(top, dict):
        return {}
    price_str = top.get("price")
    if price_str is None:
        return {}
    try:
        p = float(price_str)
    except (TypeError, ValueError):
        return {}
    if not (0 < p < 1):
        return {}
    dec = round(1.0 / p, 4)
    return {
        "decimal": dec,
        "american": _decimal_to_american(dec),
        "price_str": price_str,
        "status": top.get("status") or "",
        "raw_offers": offers,
    }


def _decimal_to_american(dec: float) -> int:
    if dec >= 2.0:
        return int(round((dec - 1) * 100))
    return int(round(-100 / (dec - 1))) if dec > 1.0 else 0
