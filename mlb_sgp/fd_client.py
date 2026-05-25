"""FanDuel API client — extracted from scraper_fanduel_sgp.py.

Owns the curl_cffi Chrome-TLS session and exposes the operations both
the SGP scraper and the singles scraper need:
  - list_events()         — all MLB events today
  - fetch_event_page()    — markets AND runners for one event tab, one GET
  - fetch_event_markets() / fetch_event_runners() — thin wrappers over
    fetch_event_page() that return just one half (back-compat)

FD returns market metadata and prices in the same event-page payload, so a
single GET yields both — fetch_event_page() is the primary entry point and
takes a `tab` argument (FD returns a different market slice per tab).
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

    def fetch_event_page(
        self, event_id: str, tab: str = "same-game-parlay-"
    ) -> tuple[list["Market"], list["Runner"]]:
        """Fetch one FD event-page tab; return (markets, runners) from it.

        FD's event-page payload carries both market metadata and runners in a
        single response, so one HTTP GET yields both. The `tab` param selects
        which server-side market slice FD returns — there is no single tab with
        every market, so the singles scraper calls this once per tab and merges
        (see fetch_merged_markets_and_runners in scraper_fanduel_singles).

        Market dedup is by marketId (the recursive walk may visit a market
        twice). Runners are read off each market dict; runnerStatus values
        other than "ACTIVE" (when present) and runners with no american odds
        or no selectionId are skipped.
        """
        from scraper_fanduel_sgp import FD_EVENT_PAGE_URL, FD_AK, FD_HEADERS

        url = (
            f"{FD_EVENT_PAGE_URL}?_ak={FD_AK}&eventId={event_id}"
            f"&tab={tab}"
            f"&useCombinedTouchdownsVirtualMarket=true&useQuickBets=true"
        )
        try:
            resp = self.session.get(url, headers=FD_HEADERS, timeout=20)
        except TypeError:
            # FakeSession in unit tests doesn't accept headers/timeout kwargs
            resp = self.session.get(url)

        if getattr(resp, "status_code", 200) != 200:
            return [], []
        payload = resp.json()

        market_dicts: list[dict] = []

        def walk(o):
            if isinstance(o, dict):
                if "marketId" in o and "marketName" in o:
                    market_dicts.append(o)
                for v in o.values():
                    walk(v)
            elif isinstance(o, list):
                for it in o:
                    walk(it)

        walk(payload)

        # Dedup by marketId; first-seen wins (duplicates are identical).
        seen: dict[str, dict] = {}
        for m in market_dicts:
            mid = str(m.get("marketId", ""))
            if mid and mid not in seen:
                seen[mid] = m

        markets = [
            Market(market_id=mid, name=str(m.get("marketName", "") or ""))
            for mid, m in seen.items()
        ]

        runners: list[Runner] = []
        for mid, m in seen.items():
            for run in m.get("runners", []) or []:
                if not isinstance(run, dict):
                    continue
                status = str(run.get("runnerStatus", "")).upper()
                if status and status != "ACTIVE":
                    continue
                american = _extract_american_odds(run)
                if american is None:
                    continue
                sid = run.get("selectionId")
                if sid is None:
                    continue
                runners.append(Runner(
                    runner_id=str(sid),
                    market_id=mid,
                    name=str(run.get("runnerName", "") or ""),
                    line=_extract_line(run),
                    american_odds=american,
                ))
        return markets, runners

    def fetch_event_runners(
        self, event_id: str, tab: str = "same-game-parlay-"
    ) -> list["Runner"]:
        """Back-compat wrapper: runners from one tab."""
        _, runners = self.fetch_event_page(event_id, tab)
        return runners

    def fetch_event_markets(
        self, event_id: str, tab: str = "same-game-parlay-"
    ) -> list["Market"]:
        """Back-compat wrapper: market metadata from one tab."""
        markets, _ = self.fetch_event_page(event_id, tab)
        return markets


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
