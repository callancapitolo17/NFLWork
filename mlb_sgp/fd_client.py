"""FanDuel API client — extracted from scraper_fanduel_sgp.py.

Owns the curl_cffi Chrome-TLS session and exposes the two operations both
the SGP scraper and the singles scraper need:
  - list_events()         — all MLB events today
  - fetch_event_runners() — all runners with prices, flat list, in one call

Mirrors DraftKingsClient with FD-specific naming (Runner instead of
Selection — FD's term — and only two operations instead of three, because
FD returns market metadata and prices in the same event-page payload).
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
class Runner:
    runner_id: str
    market_id: str
    name: str            # e.g. "Over 9.5", "New York Yankees -1.5"
    line: float | None   # FD's handicap; None when 0/absent and not a line market
    american_odds: int


@dataclass
class Market:
    market_id: str
    name: str            # FD marketName, e.g. "Run Line", "First 5 Innings Total Runs"


class FanDuelClient:
    def __init__(self, verbose: bool = False) -> None:
        from scraper_fanduel_sgp import init_session
        self.session = init_session()
        self.verbose = verbose

    def list_events(self) -> list[Event]:
        """Returns one Event per MLB game today.

        Thin wrapper over the legacy fetch_fd_events function. That function
        already does the recursive walk() over the scan payload and parses
        the 'Away (P) @ Home (P)' name format into fd_away/fd_home, so we
        just adapt its dict output to the Event dataclass.

        fetch_fd_events returns dicts with keys:
          fd_event_id, name, fd_away, fd_home, open_date
        """
        from scraper_fanduel_sgp import fetch_fd_events
        raw = fetch_fd_events(self.session)
        return [Event(
            event_id=str(e["fd_event_id"]),
            home_team=e["fd_home"],
            away_team=e["fd_away"],
            start_time=e.get("open_date", ""),
        ) for e in raw]

    def fetch_event_runners(self, event_id: str) -> list[Runner]:
        """Returns all runners (selections) with prices for one event.

        Hits FD's event-page endpoint and walks the response to flatten every
        market's runners into one list. The SGP scraper's helper of the same
        name returns a nested SGP-filtered dict ({fg/f5: {spreads/totals: ...}})
        which is the wrong shape for singles, so we deliberately do NOT wrap
        it — we walk the raw payload directly.

        FD price path: runner.winRunnerOdds.americanDisplayOdds.americanOdds
        (already an int — no string parsing needed). Runners with
        runnerStatus != "ACTIVE" are skipped.
        """
        from scraper_fanduel_sgp import FD_EVENT_PAGE_URL, FD_AK, FD_HEADERS

        url = (
            f"{FD_EVENT_PAGE_URL}?_ak={FD_AK}&eventId={event_id}"
            f"&tab=same-game-parlay-"
            f"&useCombinedTouchdownsVirtualMarket=true&useQuickBets=true"
        )
        try:
            resp = self.session.get(url, headers=FD_HEADERS, timeout=20)
        except TypeError:
            # FakeSession in unit tests doesn't accept headers/timeout kwargs
            resp = self.session.get(url)

        if getattr(resp, "status_code", 200) != 200:
            return []
        payload = resp.json()

        markets: list[dict] = []

        def walk(o):
            if isinstance(o, dict):
                # FD market dicts carry marketId + runners + marketName
                if "marketId" in o and "runners" in o and "marketName" in o:
                    markets.append(o)
                for v in o.values():
                    walk(v)
            elif isinstance(o, list):
                for it in o:
                    walk(it)

        walk(payload)

        # Dedupe by marketId — the walk may visit the same market twice
        seen: dict[str, dict] = {}
        for m in markets:
            mid = str(m.get("marketId", ""))
            if mid and mid not in seen:
                seen[mid] = m

        out: list[Runner] = []
        for mid, m in seen.items():
            for run in m.get("runners", []) or []:
                if not isinstance(run, dict):
                    continue
                status = str(run.get("runnerStatus", "")).upper()
                # Empty status is allowed; only filter known-disabled statuses
                if status and status != "ACTIVE":
                    continue
                american = _extract_american_odds(run)
                if american is None:
                    continue
                line = _extract_line(run)
                sid = run.get("selectionId")
                if sid is None:
                    continue
                out.append(Runner(
                    runner_id=str(sid),
                    market_id=mid,
                    name=str(run.get("runnerName", "") or ""),
                    line=line,
                    american_odds=american,
                ))
        return out

    def fetch_event_markets(self, event_id: str) -> list[Market]:
        """Returns market metadata (id + name) for one event.

        FD's event-page endpoint returns markets + runners in one response,
        but `fetch_event_runners` only emits Runner objects. To classify each
        runner's market we also need the market NAME, which is on the same
        payload — so this method re-walks the same response and emits Markets.

        For the singles scraper we call both methods, then build a
        `market_id -> (period, market_type)` map via classify_market.

        We intentionally do NOT share state between the two calls (no caching)
        because (a) the FD payload is small, (b) keeping the two methods
        independent matches DraftKingsClient's separation, and (c) the singles
        scraper is invoked once per event so doing two HTTP fetches doubles
        traffic but stays well under FD's rate limits.
        """
        from scraper_fanduel_sgp import FD_EVENT_PAGE_URL, FD_AK, FD_HEADERS

        url = (
            f"{FD_EVENT_PAGE_URL}?_ak={FD_AK}&eventId={event_id}"
            f"&tab=same-game-parlay-"
            f"&useCombinedTouchdownsVirtualMarket=true&useQuickBets=true"
        )
        try:
            resp = self.session.get(url, headers=FD_HEADERS, timeout=20)
        except TypeError:
            resp = self.session.get(url)

        if getattr(resp, "status_code", 200) != 200:
            return []
        payload = resp.json()

        markets: list[dict] = []

        def walk(o):
            if isinstance(o, dict):
                if "marketId" in o and "marketName" in o:
                    markets.append(o)
                for v in o.values():
                    walk(v)
            elif isinstance(o, list):
                for it in o:
                    walk(it)

        walk(payload)

        seen: dict[str, str] = {}
        for m in markets:
            mid = str(m.get("marketId", ""))
            name = str(m.get("marketName", "") or "")
            if mid and mid not in seen:
                seen[mid] = name

        return [Market(market_id=mid, name=name) for mid, name in seen.items()]


def _extract_american_odds(run: dict) -> int | None:
    """Extract american odds from a runner dict.

    Path: winRunnerOdds.americanDisplayOdds.americanOdds (int, e.g. 520, -110).
    Falls back to americanOddsInt (same value, different key). Returns None
    if neither is present or parseable.
    """
    try:
        wro = run.get("winRunnerOdds") or {}
        ado = wro.get("americanDisplayOdds") or {}
    except AttributeError:
        return None
    raw = ado.get("americanOdds")
    if raw is None:
        raw = ado.get("americanOddsInt")
    if raw is None:
        return None
    try:
        return int(raw)
    except (TypeError, ValueError):
        return None


def _extract_line(run: dict) -> float | None:
    """Extract the handicap/line for a runner.

    FD uses `handicap` on the runner (numeric, may be 0 for moneyline-style
    selections). We return None for handicap=0 to match the singles scraper's
    expectation that ML/3-way runners have no line — alt spreads/totals
    always carry a non-zero handicap.
    """
    hc = run.get("handicap")
    if hc is None:
        return None
    try:
        v = float(hc)
    except (TypeError, ValueError):
        return None
    if v == 0.0:
        return None
    return v
