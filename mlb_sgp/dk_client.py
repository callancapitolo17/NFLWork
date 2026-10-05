"""DraftKings API client — extracted from scraper_draftkings_sgp.py.

Owns the curl_cffi Chrome-TLS session and exposes the operations the SGP
scraper, the singles scraper and the on-demand pricer need:
  - list_events()           — all MLB events today
  - fetch_event_markets()   — market metadata (FG vs F5 vs alt)
  - fetch_event_selections()— all selections with prices, in one call
  - fetch_added_leg_odds()  — same-game parlay prices from DK's SGP builder
                              (one request at a time, through ``price_lane``)

No DB I/O.
"""
from __future__ import annotations
import logging
import threading
import time
from dataclasses import dataclass
from typing import Callable, Hashable

from mlb_sgp._shared import (RETRY_BACKGROUND, RETRY_LIVE, BookTransportError,
                             PriceCallTallyMixin, RetryProfile,
                             check_response, json_or_raise,
                             request_with_retry)

logger = logging.getLogger(__name__)

BOOK = "draftkings"

# The price call DK's own SGP-builder widget makes as you click legs
# (dk-same-game-parlay 9.0.1). Unlike the betslip's calculateBets — Akamai-
# gated to a real Chrome since ~2026-08-20 (issue #102) — it answers plain
# curl_cffi. Given base legs and one candidate market it returns the price of
# base + EACH outcome of that market, so one request prices both sides of a
# leg. Verified on MLB 2026-10-03: 21 leg sets within 0.63 pp (mean) of the
# other books' median fair. See mlb_sgp/README.md § DK same-game pricing via
# the SGP builder.
DK_SGP_PRICE_URL = ("https://sportsbook-nash.draftkings.com/sites/US-WV-SB/api/"
                    "sportscontent/sgp/dkuswv/sportsdata/v2/sgp")
# DK's sport id for baseball (football is 3). The endpoint 404s without it.
DK_BASEBALL_SPORT_ID = "7"
# Concurrent bursts are what drew Akamai's rule onto calculateBets, so DK's
# price calls go one at a time, no faster than DK's own widget throttles
# itself (one per 500 ms). Burst-tested 2026-10-03: ~750 calls, no 403/429.
MIN_SECONDS_BETWEEN_PRICE_CALLS = 0.5
# A caller that cannot get the lane this fast declines: the quote budget is
# 10 s and the other books carry the RFQ.
PRICE_LANE_MAX_WAIT_SEC = 3.0
# ~1% of calls take 2-8 s and ~0.4% hang (DK-side, not throttling: the same
# request replays in 0.2 s), so a short timeout keeps one hang from holding
# the lane for everyone queued behind it.
PRICE_CALL_TIMEOUT_SEC = 3.0
# Long enough for one RFQ's partition cells (priced within ~1 s) to share the
# answer that covers both sides of their last leg; short enough that no quote
# prices off another RFQ's number.
SIBLING_ANSWER_TTL_SEC = 2.0
# DK statuses that mean "DK refuses US", as opposed to "DK won't price this".
BLOCKED_PRICE_STATUSES = frozenset({403, 429})
# A blocked price call is otherwise silent: price_selection_set counts it
# without raising, and an on-demand "empty" outcome never trips the
# book-health alert. It is the kill-switch signal, so it is logged — at most
# once a minute, since a blocked lane fails up to twice a second.
BLOCKED_WARNING_INTERVAL_SEC = 60.0


class PriceLane:
    """One-at-a-time, spaced, briefly-remembered gate for DK's price GET.

    ``answer(key, fetch)`` returns the answer ``fetch()`` produced for
    ``key`` within the last ``answer_ttl_sec``, else waits for the lane
    (at most ``max_wait_sec``), spaces itself ``min_gap_sec`` after the
    previous call and runs ``fetch()``. Returns None when the lane stayed
    busy. Exceptions from ``fetch()`` propagate and are not remembered.

    Thread-safe: SGPService prices partition cells — and separate RFQs —
    on concurrent threads. ``clock``/``sleep`` are injectable for tests.
    """

    def __init__(self, *, min_gap_sec: float, max_wait_sec: float,
                 answer_ttl_sec: float,
                 clock: Callable[[], float] = time.monotonic,
                 sleep: Callable[[float], None] = time.sleep) -> None:
        self.min_gap_sec = min_gap_sec
        self.max_wait_sec = max_wait_sec
        self.answer_ttl_sec = answer_ttl_sec
        self._clock = clock
        self._sleep = sleep
        self._lane = threading.Lock()
        self._answers_lock = threading.Lock()
        self._answers: dict[Hashable, tuple[float, object]] = {}
        self._last_call_at: float | None = None

    def answer(self, key: Hashable, fetch: Callable[[], object]):
        # Fast path outside the lane: a sibling cell's answer must not queue
        # behind other RFQs' wire calls.
        remembered = self._remembered(key)
        if remembered is not None:
            return remembered
        if not self._lane.acquire(timeout=self.max_wait_sec):
            return None
        try:
            # Re-check: the call we queued behind may have fetched this key.
            remembered = self._remembered(key)
            if remembered is not None:
                return remembered
            if self._last_call_at is not None:
                wait = self.min_gap_sec - (self._clock() - self._last_call_at)
                if wait > 0:
                    self._sleep(wait)
            try:
                result = fetch()
            finally:
                self._last_call_at = self._clock()
            self._remember(key, result, self._last_call_at)
            return result
        finally:
            self._lane.release()

    def _remembered(self, key: Hashable):
        with self._answers_lock:
            entry = self._answers.get(key)
        if entry is None:
            return None
        fetched_at, result = entry
        if self._clock() - fetched_at >= self.answer_ttl_sec:
            return None
        return result

    def _remember(self, key: Hashable, result, fetched_at: float) -> None:
        with self._answers_lock:
            # Drop expired answers on every write so the dict stays at one
            # flight's worth of keys.
            expired = [k for k, (at, _) in self._answers.items()
                       if fetched_at - at >= self.answer_ttl_sec]
            for k in expired:
                del self._answers[k]
            self._answers[key] = (fetched_at, result)


def parse_added_leg_odds(body: dict, base: frozenset[str]) -> dict[str, float]:
    """Map DK's SGP-builder answer to {candidate selection id: decimal odds}.

    DK does not refuse a base it cannot combine: it drops the offending legs
    (``selectionsNotMapped``, e.g. both sides of one market) and prices the
    rest. A price on a smaller base would be a silent wrong number, so a base
    DK did not take whole prices nothing. A candidate absent from
    ``compatibleMarkets`` is one DK declines. Same rule as
    ``nfl_specials/dk_book.py``.
    """
    data = body.get("data") or {}
    if data.get("selectionsNotMapped") or set(data.get("selectionsMapped") or []) != base:
        return {}
    return {selection["id"]: float(selection["trueOdds"])
            for market in data.get("compatibleMarkets") or []
            for selection in market.get("selections") or []
            if selection.get("trueOdds")}


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


class DraftKingsClient(PriceCallTallyMixin):
    BOOK = BOOK

    def __init__(self, verbose: bool = False) -> None:
        from scraper_draftkings_sgp import init_session
        self.session = init_session()
        self.verbose = verbose

    @property
    def price_lane(self) -> PriceLane:
        """The client's one DK price lane, created on first use.

        Lazy for the same reason as ``PriceCallTallyMixin.price_calls``:
        client tests build instances via ``__new__``. ``setdefault`` is
        atomic under the GIL, so concurrent first calls share ONE lane.
        """
        lane = self.__dict__.get("_price_lane")
        if lane is None:
            lane = self.__dict__.setdefault("_price_lane", PriceLane(
                min_gap_sec=MIN_SECONDS_BETWEEN_PRICE_CALLS,
                max_wait_sec=PRICE_LANE_MAX_WAIT_SEC,
                answer_ttl_sec=SIBLING_ANSWER_TTL_SEC))
        return lane

    def fetch_added_leg_odds(self, event_id: str, base: frozenset[str],
                             candidate_market_id: str) -> dict[str, float]:
        """Decimal odds of ``base`` + each selection of the candidate market
        that DK can combine with it, via the SGP-builder GET.

        One wire call at most, through ``price_lane``; a sibling partition
        cell asking the same (event, base, market) within
        ``SIBLING_ANSWER_TTL_SEC`` reuses the answer. ``{}`` = DK declines.

        Raises ``BookTransportError`` (stage "price") when DK refuses us
        (403/429), the call fails on the wire, or the lane stays busy past
        ``PRICE_LANE_MAX_WAIT_SEC``. Any other non-200 (400 InvalidSelection:
        a line DK pulled since the structure was read) is a decline.
        """
        if not base:
            raise ValueError("DK's SGP builder prices nothing on an empty base")
        key = (event_id, base, candidate_market_id)
        answer = self.price_lane.answer(
            key, lambda: self._get_added_leg_odds(event_id, base,
                                                  candidate_market_id))
        if answer is None:
            raise BookTransportError(
                BOOK, "price",
                detail=f"price lane busy for >{PRICE_LANE_MAX_WAIT_SEC:g}s")
        return answer

    def _get_added_leg_odds(self, event_id: str, base: frozenset[str],
                            candidate_market_id: str) -> dict[str, float]:
        params = {"eventId": event_id,
                  "selections": ",".join(sorted(base)),
                  "marketCandidates": candidate_market_id,
                  "oddsStyle": "american"}
        # RETRY_LIVE: on-demand quote path — one fast retry on a wire blip.
        resp = request_with_retry(
            lambda: self.session.get(
                DK_SGP_PRICE_URL, params=params,
                headers={"X-SportId": DK_BASEBALL_SPORT_ID},
                timeout=PRICE_CALL_TIMEOUT_SEC),
            profile=RETRY_LIVE, book=BOOK, stage="price")
        if resp.status_code in BLOCKED_PRICE_STATUSES:
            self._warn_blocked(resp.status_code)
            raise BookTransportError(
                BOOK, "price", status_code=resp.status_code,
                detail=(getattr(resp, "text", "") or "")[:200])
        if resp.status_code != 200:
            logger.debug("draftkings: SGP builder declined HTTP %s: %s",
                         resp.status_code, (resp.text or "")[:200])
            return {}
        return parse_added_leg_odds(json_or_raise(BOOK, "price", resp), base)

    def _warn_blocked(self, status_code: int) -> None:
        now = time.monotonic()
        last = getattr(self, "_last_blocked_warning_at", None)
        if last is not None and now - last < BLOCKED_WARNING_INTERVAL_SEC:
            return
        self._last_blocked_warning_at = now
        logger.warning(
            "draftkings: SGP builder answered HTTP %s — DK may be blocking "
            "us. Drop draftkings from SGP_BOOKS / RFI_BOOKS and restart "
            "before a rule reaches DK's read endpoints (mlb_sgp/README.md "
            "§ DK same-game pricing via the SGP builder).", status_code)

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

    def fetch_event_selections(self, event_id: str,
                               profile: RetryProfile = RETRY_BACKGROUND) -> list[Selection]:
        """Returns all selections (with prices) for one event from the parlays endpoint.

        DK's parlays/v1/sgp/events/{event_id} returns a payload of shape:
          data.markets[] -> { id, name, selections[] -> {
              id, marketId, name, outcomeType, points (float|null),
              displayOdds.american (str like "+260" or Unicode-minus "−370"),
              isDisabled, status }}

        Selections that are disabled or missing displayOdds.american are skipped
        — the singles scraper can't price an unavailable bet.

        A 404 means DK dropped this event; anything else non-200 (DK has been
        serving 403 in production) raises ``BookTransportError``.
        """
        from scraper_draftkings_sgp import DK_SGP_PARLAYS_URL
        resp = request_with_retry(
            lambda: self.session.get(f"{DK_SGP_PARLAYS_URL}/{event_id}",
                                     timeout=60),
            profile=profile, book=BOOK, stage="structure")
        if not check_response(BOOK, "structure", resp, allow_404=True):
            return []
        payload = json_or_raise(BOOK, "structure", resp)
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
