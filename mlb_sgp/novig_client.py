"""Novig HTTP client — extracted from scraper_novig_sgp.py.

No auth anywhere: events come from the REST page novig.com's MLB screen
loads and prices from the parlay endpoint (literally named
/unauthenticated). The market tree comes from Novig's Hasura GraphQL
endpoint through scraper_novig_sgp.fetch_event_legs -> _gql, the one
GraphQL transport every path (sweep, on-demand, warming, leg surface)
uses.

Novig does not originate SGP prices: it relays each selection set from one
vendor book, naming it per leg (`vendor`), and returns that book's price
shaded slightly short. Until 2026-09 every leg said DRAFTKINGS. On
2026-10-02 the vendor was chosen per selection set — 88% of 2^N partitions
mixed 2-3 vendors (cells: DRAFTKINGS 60%, FANDUEL 29%, BETMGM 11%, CAESARS
<1%), 30% of cells switched vendor within 20 minutes, FanDuel-relayed cells
came back ~1% shorter than FanDuel's own quote and BetMGM's unshaded.
A Novig price is therefore a vendor's opinion, not an independent one;
`submit_parlay` returns the per-leg `vendors` so callers can tell which.

Since 2026-09-22 the GraphQL endpoint runs an allowlist: it executes only
operations the novig.com app ships and answers anything else with HTTP 200
`{"errors": [{"message": "query is not allowed", "extensions": {"code":
"validation-failed", ...}}]}`. That killed the hand-written events query
(events moved to REST) and the committed EventMarkets_Query, which must
now track the app — refresh it with `python3 mlb_sgp/refresh_novig_query.py`.

Exposes two thin methods, mirroring the other books' clients:
  - list_events()    -> list[Event]
  - submit_parlay()  -> dict (parsed parlay node)

The pure helper `_parse_events_response` is at module level so tests can
exercise the parser without a session.

Real Novig response shapes:

    GET https://api.novig.us/nbx/v1/trading/MLB/page  (captured 2026-10-02)
      -> {"sections": [
            {"title": "Games", "content": {"type": "components", "components": [
                {"type": "game_event_card", "eventId": "<uuid>",
                 "scheduledStart": "2026-10-03T17:00:00.000Z",
                 "eventStatus": "OPEN_PREGAME",
                 "homeTeam": {"name": "Cleveland Guardians", "symbol": "CLE", ...},
                 "awayTeam": {"name": "Chicago White Sox", "symbol": "CHI",
                              "shortName": "CWS", ...}},
                ...]}},
            {"title": "Featured Parlays", ...}, {"title": "Series", ...},
            {"title": "Futures", "content": {"type": "subsections",
                                             "subSections": [...]}}]}
       (an off-season league answers {"sections": []})

    POST https://api.novig.us/nbx/v1/parlay/request/unauthenticated
      -> [{"price": "0.36000", "status": "Unfilled",
           "legs": [{"price": "0.57100", "vendor": "DRAFTKINGS",
                     "outcomeId": "<uuid>", "outcome": {...}}, ...]}, ...]
       (a *list* of offer dicts — top offer at [0]; captured 2026-10-02)
"""
from __future__ import annotations
import logging

from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from typing import Any

from mlb_sgp._shared import (RETRY_BACKGROUND, RETRY_LIVE, BookTransportError,
                             PriceCallTallyMixin, RetryProfile, check_response,
                             json_or_raise, request_with_retry)

logger = logging.getLogger(__name__)

BOOK = "novig"
NOVIG_PARLAY = "https://api.novig.us/nbx/v1/parlay/request/unauthenticated"
# The REST call novig.com's league screens make. It lists every upcoming
# game Novig has posted (NFL's page spanned 8 weeks on 2026-10-02), so
# list_events applies the time window itself.
NOVIG_LEAGUE_PAGE = "https://api.novig.us/nbx/v1/trading/{league}/page"
NOVIG_MLB_PAGE = NOVIG_LEAGUE_PAGE.format(league="MLB")
EVENT_WINDOW_HOURS = 48   # how far ahead list_events looks for upcoming games
GAME_CARD_TYPE = "game_event_card"
PREGAME_STATUS = "OPEN_PREGAME"


@dataclass
class Event:
    event_id: str
    home_team: str
    away_team: str
    home_sym: str
    away_sym: str
    start_time: str  # ISO UTC string, "+00:00" offset


class NovigClient(PriceCallTallyMixin):
    """Thin wrapper around Novig's REST league page and anonymous parlay
    endpoint (the GraphQL market tree is fetched by scraper_novig_sgp)."""

    BOOK = BOOK

    def __init__(self, verbose: bool = False) -> None:
        # Reuse the legacy scraper's session bootstrap (curl_cffi Chrome
        # impersonation + landing-page warm-up for Cloudflare cookies).
        from scraper_novig_sgp import init_session
        self.session = init_session()
        self.verbose = verbose

    def list_events(self, profile: RetryProfile = RETRY_BACKGROUND) -> list[Event]:
        """Upcoming pregame MLB events from Novig's MLB league page (REST).

        Side effects: one GET to NOVIG_MLB_PAGE. Raises BookTransportError
        on a non-200, a non-JSON body, or a page without a `sections` list.
        """
        # A connection failure that survives the retries is a dead book, not
        # an off-day. (The 2026-07 "api.novig.us unresolvable" reports were a
        # transient blip — the host resolves and answers; see issue #40.)
        r = request_with_retry(lambda: self.session.get(NOVIG_MLB_PAGE, timeout=20),
                               profile=profile, book=BOOK, stage="events")
        check_response(BOOK, "events", r)
        return _parse_events_response(json_or_raise(BOOK, "events", r),
                                      now=datetime.now(timezone.utc),
                                      window_hours=EVENT_WINDOW_HOURS)

    def submit_parlay(self, outcome_ids: list[str], stake: float = 1.0) -> dict:
        """Submit a BuildParlay request and return the parsed top offer node.

        Returns a dict shaped like:
            {"decimal": 2.85, "american": -200, "price_str": "0.35088",
             "status": "OPEN", "vendors": ["FANDUEL", "FANDUEL"],
             "raw_offers": [...]}
        or `{}` for any failure (HTTP non-2xx, empty response, malformed price).

        We use `submit_parlay` for the method name (not `submit_parlay_rfq`)
        to match the task spec, but the underlying endpoint is the same
        anonymous RFQ-style endpoint the legacy scraper hits. `stake` is
        accepted for API symmetry with the spec but Novig's anonymous
        endpoint doesn't gate by stake — it returns offers regardless.
        """
        payload = {"outcomes": [{"id": oid} for oid in outcome_ids], "boostId": None}

        def _post():
            try:
                return self.session.post(
                    NOVIG_PARLAY, json=payload,
                    headers={"Content-Type": "application/json"},
                    timeout=15,
                )
            except TypeError:
                return self.session.post(NOVIG_PARLAY, json=payload)

        try:
            # RETRY_LIVE on every path — price calls are the fan-out surface.
            r = request_with_retry(_post, profile=RETRY_LIVE, book=BOOK,
                                   stage="price")
        except BookTransportError as e:
            # Per-combo failures stay row-drops, but must be tallied so an
            # all-fail cycle still yields a "price" transport verdict.
            self._decline(e.status_code)
            return {}
        status = getattr(r, "status_code", 200)
        if status not in (200, 201):
            # Novig's two non-200s mean opposite things (issue #40):
            #   400 "Cannot price parlay" -> the book declines THIS combo,
            #                                which is the common, normal case
            #   403 <html>                -> we have been RATE-LIMITED
            # Recording the status is what lets sgp_fetch_health.error_class
            # say which, instead of an ambiguous bare "price".
            self._decline(status)
            return {}
        try:
            offers = r.json()
        except Exception:
            self._decline(None)
            return {}
        parsed = _parse_parlay_response(offers)
        if not parsed:
            # A 200 whose body carries no usable price is still a decline, but
            # the status says nothing about why — don't attribute one.
            self._decline(None)
            return parsed
        self.price_calls.record(True)
        return parsed

    def _decline(self, status: int | None) -> None:
        """Tally one failed price call, carrying its HTTP status when there is
        one. Order matters: ``note_status`` stamps the in-flight attempt, so it
        must precede ``record`` (see PriceCallTally)."""
        if status is not None:
            self.price_calls.note_status(status)
        self.price_calls.record(False)


# ---------------------------------------------------------------------------
# Pure parser helpers — exposed at module level for tests
# ---------------------------------------------------------------------------

def _parse_events_response(raw: dict, now: datetime,
                           window_hours: float) -> list[Event]:
    """Parse Novig's MLB league page into upcoming pregame Events.

    Reads every `game_event_card` on the page (shape in the module
    docstring; today they all sit in the "Games" section) and keeps what the
    retired GraphQL events query's WHERE clause kept: eventStatus
    OPEN_PREGAME and a scheduledStart inside [now, now + window_hours].
    Cards are found by type wherever they are nested (`_iter_game_cards`),
    so a renamed section or a Games section regrouped into subsections
    still parses instead of reading as an off-day. One Event per eventId, so
    a game shown in two places cannot match twice (two matches make the
    on-demand path decline the game as ambiguous).

    The team key is `symbol`, not `shortName`: the market tree's outcome
    competitors carry `symbol` (the White Sox are symbol "CHI", shortName
    "CWS"). Both Chicago clubs share "CHI"; the outcome matcher declines
    rather than guess (scraper_novig_sgp._find_outcome_in_spread).

    An off-season league answers {"sections": []} -> []. A body without a
    `sections` list is a page we no longer understand -> BookTransportError,
    so a shape change reads as a dead book, not as an empty slate. A game
    card whose start time does not parse is skipped with a WARNING.
    """
    sections = raw.get("sections") if isinstance(raw, dict) else None
    if not isinstance(sections, list):
        raise BookTransportError(
            BOOK, "events",
            detail=f"MLB page has no 'sections' list (got {str(raw)[:120]!r})")

    window_end = now + timedelta(hours=window_hours)
    out: list[Event] = []
    seen_event_ids: set[str] = set()
    for card in _iter_game_cards(sections):
        if card.get("eventStatus") != PREGAME_STATUS:
            continue
        start = _parse_iso_utc(card.get("scheduledStart"))
        if start is None:
            logger.warning("novig: game card %s has an unparseable "
                           "scheduledStart %r — skipped",
                           card.get("eventId"), card.get("scheduledStart"))
            continue
        if not (now <= start <= window_end):
            continue
        home = card.get("homeTeam")
        away = card.get("awayTeam")
        if not (isinstance(home, dict) and isinstance(away, dict)):
            continue
        event_id = card.get("eventId")
        if not (event_id and home.get("name") and away.get("name")):
            continue
        if str(event_id) in seen_event_ids:
            continue
        seen_event_ids.add(str(event_id))
        out.append(Event(
            event_id=str(event_id),
            home_team=home["name"],
            away_team=away["name"],
            home_sym=home.get("symbol") or home.get("shortName") or "",
            away_sym=away.get("symbol") or away.get("shortName") or "",
            # Normalized to UTC: match_events buckets on the first 13 chars
            # and assumes they are a UTC hour.
            start_time=start.astimezone(timezone.utc).isoformat(),
        ))
    return out


def _iter_game_cards(node):
    """Yield every `game_event_card` dict anywhere under ``node``.

    Novig nests page content in several container shapes (the novig.com
    bundle's TradingSectionType: components, subsections,
    nested_subsections; the MLB page's Futures section already uses
    subSections), so cards are found by type rather than by a fixed path.
    A card's own fields are never searched.
    """
    if isinstance(node, dict):
        if node.get("type") == GAME_CARD_TYPE:
            yield node
            return
        for value in node.values():
            yield from _iter_game_cards(value)
    elif isinstance(node, list):
        for item in node:
            yield from _iter_game_cards(item)


def _parse_iso_utc(text) -> datetime | None:
    """'2026-10-03T17:00:00.000Z' -> aware UTC datetime; None if unparseable.

    A timestamp without an offset is read as UTC, the same rule as
    _shared._utc_bucket and verify_books.
    """
    if not text:
        return None
    try:
        parsed = datetime.fromisoformat(str(text).replace("Z", "+00:00"))
    except ValueError:
        return None
    if parsed.tzinfo is None:
        return parsed.replace(tzinfo=timezone.utc)
    return parsed


def _leg_vendors(offer: dict) -> list:
    """The vendor named on each leg of a priced offer, in leg order (None for
    a leg that names none; [] when the offer lists no legs).

    Real shape (2026-10-02): ``{"price": "0.36000", "legs": [{"price":
    "0.57100", "vendor": "DRAFTKINGS", "outcomeId": "<uuid>", ...}, ...]}``.
    """
    legs = offer.get("legs")
    if not isinstance(legs, list):
        return []
    return [leg.get("vendor") if isinstance(leg, dict) else None
            for leg in legs]


def _parse_parlay_response(offers: Any) -> dict:
    """Parse a Novig parlay response (list-of-offers) into a single result dict.

    Returns {"decimal": float, "american": int, "price_str": str,
             "status": str, "vendors": list, "raw_offers": list} or {} on
    malformed input. ``vendors`` is the top offer's per-leg vendor names
    (see ``_leg_vendors``) — the book each leg's price was relayed from.

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
                "vendors": _leg_vendors(parlay),
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
        "vendors": _leg_vendors(top),
        "raw_offers": offers,
    }


def _decimal_to_american(dec: float) -> int:
    if dec >= 2.0:
        return int(round((dec - 1) * 100))
    return int(round(-100 / (dec - 1))) if dec > 1.0 else 0
