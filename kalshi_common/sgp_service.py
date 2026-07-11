"""In-process SGP pricing service shared by the Kalshi MLB bots.

Replaces the subprocess-per-cycle scraper model for the taker
(kalshi_mlb_rfq) and maker (kalshi_mlb_mm): one persistent HTTP client
per book held across cycles (no per-cycle TLS handshake), the four book
orchestrators run concurrently under a per-book deadline, and slow
structure fetches (event lists, DK's 2MB selection-id payload) are
TTL-cached. Prices are NEVER cached — every refresh() re-prices.

The dashboard keeps the CLI-shim subprocess path and never touches this.
"""
from __future__ import annotations
import logging
import sys
import time
from concurrent.futures import ThreadPoolExecutor, TimeoutError as FutureTimeout
from pathlib import Path

from mlb_sgp._shared import (OnDemandBookResult, PricedRow, TargetLine,
                             TTLCache)

log = logging.getLogger(__name__)

_REPO_ROOT = Path(__file__).resolve().parents[1]
_MLB_SGP_DIR = _REPO_ROOT / "mlb_sgp"

DEFAULT_BOOKS = ("draftkings", "fanduel", "prophetx", "novig", "betmgm", "caesars")

# TTLs for structure fetches (market *structure*, never prices).
EVENTS_TTL_SEC = 300       # game list churns ~daily
STRUCTURE_TTL_SEC = 180    # sel-ids / runners churn when a book re-mains a line

# ---- Phase 2 on-demand pricing constants ---- #
# Route A (full 2^N partition) only up to 3 legs: 8 wire calls per book is
# the ceiling we accept; beyond that Route B (1 SGP call + singles).
ON_DEMAND_MAX_PARTITION_LEGS = 3

# Route B haircut when a leg's single is ONE-SIDED at the book (no opposite
# to devig against): fair ≈ implied/(1+vig). Mirrors the maker's per-book
# vig-fallback config defaults (kalshi_common cannot import maker config).
ON_DEMAND_VIG_FALLBACK = {
    "draftkings": 0.125,
    "fanduel": 0.18,
    "prophetx": 0.05,
    "novig": 0.05,
    "betmgm": 0.06,
    "caesars": 0.06,
}


class _BookState:
    __slots__ = ("client", "failures", "last_success", "caches")

    def __init__(self):
        self.client = None          # persistent per-book HTTP client
        self.failures = 0           # consecutive refresh failures
        self.last_success = None    # now_fn timestamp of last success
        self.caches = {}            # name -> TTLCache


class SGPService:
    """Holds persistent book clients; prices targets at all due books
    concurrently. See refresh() for the result contract."""

    MAX_FAILURES_BEFORE_REINIT = 3

    def __init__(
        self,
        books: tuple[str, ...] = DEFAULT_BOOKS,
        per_book_deadline_sec: float = 75.0,
        min_refresh_sec: dict[str, float] | None = None,
        runners: dict | None = None,   # test seam: {book: callable(targets)->rows}
        on_demand_hooks: dict | None = None,  # test seam: {book: hook dict},
        # hook dict keys: match_event / build_structure / resolve / price
        # (see _book_on_demand_hooks for the real implementations)
        now_fn=time.monotonic,
    ):
        self.books = tuple(books)
        self.per_book_deadline_sec = per_book_deadline_sec
        self.min_refresh_sec = dict(min_refresh_sec or {})
        self._runners = runners
        self._on_demand_hooks = on_demand_hooks
        self._now = now_fn
        self._state = {b: _BookState() for b in self.books}
        # The orchestrators lazily import the legacy scraper modules by
        # top-level name (`from scraper_draftkings_sgp import ...`),
        # which only resolves with mlb_sgp/ itself on sys.path — true
        # for CLI runs (cwd=mlb_sgp/) but not for the bots (cwd=repo
        # root). Make it resolvable here, once.
        if str(_MLB_SGP_DIR) not in sys.path:
            sys.path.insert(0, str(_MLB_SGP_DIR))

    # ------------------------------------------------------------------ #
    # Public API                                                          #
    # ------------------------------------------------------------------ #

    def refresh(self, targets: list[TargetLine]) -> dict[str, list[PricedRow] | None]:
        """Price `targets` at every due book concurrently.

        Returns {book: result} for books ATTEMPTED this call:
          list[PricedRow] (possibly empty) -> success
          None                             -> failure or deadline timeout
        Books skipped by min_refresh_sec are ABSENT from the dict —
        callers must not clear their previously-written rows.
        """
        due = [b for b in self.books if self._due(b)]
        results: dict[str, list[PricedRow] | None] = {}
        if not due:
            return results
        # Fresh pool per refresh: a hung book thread from a previous
        # cycle must not occupy a worker slot forever. shutdown(wait=
        # False) lets a still-running (timed-out) thread finish in the
        # background; its client gets torn down via the failure path.
        pool = ThreadPoolExecutor(max_workers=len(due),
                                  thread_name_prefix="sgp-book")
        try:
            futs = {b: pool.submit(self._run_book_safe, b, targets) for b in due}
            # All futures are already running concurrently (one worker per
            # due book), so per-book `remaining` is wall-bounded by the
            # deadline, not summed across books.
            wall_deadline = time.monotonic() + self.per_book_deadline_sec
            for b, fut in futs.items():
                remaining = max(0.0, wall_deadline - time.monotonic())
                timed_out = False
                try:
                    rows = fut.result(timeout=remaining)
                except FutureTimeout:
                    # The worker thread is still running and still holds
                    # this book's curl_cffi session (not thread-safe).
                    # Drop the client NOW so the next cycle builds a fresh
                    # session the lingering thread no longer shares —
                    # don't wait for the 3-strike teardown.
                    rows = None
                    timed_out = True
                except Exception:
                    rows = None
                self._book_done(b, rows, timed_out=timed_out)
                results[b] = rows
        finally:
            pool.shutdown(wait=False, cancel_futures=True)
        return results

    def close(self):
        """Drop all persistent clients (sessions close on GC)."""
        for st in self._state.values():
            st.client = None
            st.caches = {}

    # ------------------------------------------------------------------ #
    # Internals                                                           #
    # ------------------------------------------------------------------ #

    def _due(self, book: str) -> bool:
        min_s = self.min_refresh_sec.get(book, 0)
        st = self._state[book]
        if not min_s or st.last_success is None:
            return True
        return (self._now() - st.last_success) >= min_s

    def _book_done(self, book: str, rows, timed_out: bool = False) -> None:
        st = self._state[book]
        if rows is None:
            st.failures += 1
            log.warning("sgp_service: %s failed (consecutive=%d, timeout=%s)",
                        book, st.failures, timed_out)
            if timed_out or st.failures >= self.MAX_FAILURES_BEFORE_REINIT:
                log.warning("sgp_service: %s client torn down for reinit "
                            "(timeout=%s)", book, timed_out)
                st.client = None
                st.caches = {}
                st.failures = 0
        else:
            st.failures = 0
            st.last_success = self._now()

    def _run_book_safe(self, book: str, targets):
        """Runs in a worker thread. Returns rows or None — never raises
        (a raise would surface as a generic future error and lose the
        book attribution in logs)."""
        try:
            if self._runners is not None:
                return self._runners[book](targets)
            return self._run_book(book, targets)
        except Exception as e:
            log.warning("sgp_service: %s runner error: %s", book, e)
            return None

    def _ensure_client(self, book: str) -> _BookState:
        """Lazily build the book's persistent client — extracted verbatim
        from _run_book so price_on_demand shares the exact same init. DK/FD
        also get their structure TTL caches here (as before); the other
        books' refresh path never reads caches so none are created at init
        (identical to the pre-refactor behavior — on-demand adds its own via
        _structure_caches)."""
        st = self._state[book]
        if st.client is not None:
            return st
        if book == "draftkings":
            from mlb_sgp.dk_client import DraftKingsClient
            st.client = DraftKingsClient()
            st.caches = {"events": TTLCache(EVENTS_TTL_SEC),
                         "structure": TTLCache(STRUCTURE_TTL_SEC)}
        elif book == "fanduel":
            from mlb_sgp.fd_client import FanDuelClient
            st.client = FanDuelClient()
            st.caches = {"events": TTLCache(EVENTS_TTL_SEC),
                         "structure": TTLCache(STRUCTURE_TTL_SEC)}
        elif book == "prophetx":
            from mlb_sgp.prophetx_client import ProphetXClient
            st.client = ProphetXClient()
        elif book == "novig":
            from mlb_sgp.novig_client import NovigClient
            st.client = NovigClient()
        elif book == "betmgm":
            from mlb_sgp.betmgm_client import BetMGMClient
            st.client = BetMGMClient()
        elif book == "caesars":
            from mlb_sgp.caesars_client import CaesarsClient
            st.client = CaesarsClient()
        else:
            raise ValueError(f"unknown book {book!r}")
        return st

    def _run_book(self, book: str, targets):
        st = self._ensure_client(book)
        if book == "draftkings":
            from mlb_sgp import draftkings as mod
            return mod.price_sgps(targets, periods=("FG",), client=st.client,
                                  fetchers=self._dk_fetchers(st))
        if book == "fanduel":
            from mlb_sgp import fanduel as mod
            return mod.price_sgps(targets, periods=("FG",), client=st.client,
                                  fetchers=self._fd_fetchers(st))
        if book == "prophetx":
            from mlb_sgp import prophetx as mod
        elif book == "novig":
            from mlb_sgp import novig as mod
        elif book == "betmgm":
            from mlb_sgp import betmgm as mod
        else:            # caesars — unknown books already raised above
            from mlb_sgp import caesars as mod
        return mod.price_sgps(targets, periods=("FG",), client=st.client)

    @staticmethod
    def _dk_fetchers(st: _BookState) -> dict:
        import scraper_draftkings_sgp as legacy
        ev, struct = st.caches["events"], st.caches["structure"]
        return {
            "fetch_dk_events": lambda session: ev.get_or_fetch(
                "dk_events", lambda: legacy.fetch_dk_events(session)),
            "fetch_main_market_nums": lambda session, eid: struct.get_or_fetch(
                ("nums", eid), lambda: legacy.fetch_main_market_nums(session, eid)),
            "fetch_selection_ids": lambda session, eid, nums, verbose=False:
                struct.get_or_fetch(
                    ("selids", eid),
                    lambda: legacy.fetch_selection_ids(session, eid, nums, verbose)),
        }

    @staticmethod
    def _fd_fetchers(st: _BookState) -> dict:
        import scraper_fanduel_sgp as legacy
        ev, struct = st.caches["events"], st.caches["structure"]
        return {
            "fetch_fd_events": lambda session: ev.get_or_fetch(
                "fd_events", lambda: legacy.fetch_fd_events(session)),
            "fetch_event_runners": lambda session, eid, h, a: struct.get_or_fetch(
                ("runners", eid),
                lambda: legacy.fetch_event_runners(session, eid, h, a)),
        }
    # ------------------------------------------------------------------ #
    # Phase 2: on-demand pricing of one same-game leg set                 #
    # (spec §4.4). Runs on the CALLER's thread, concurrently with          #
    # refresh() — no cross-lock (each book's sweep already fans one        #
    # curl session across a thread pool, so this matches existing          #
    # practice; a cross-lock would stall the feed ~30-75s per sweep).      #
    # ------------------------------------------------------------------ #

    def price_on_demand(self, book: str, game, legs):
        """Price one same-game leg set at one book. Never raises.

        Parameters: ``game`` is an ``mlb_sgp._shared.GameRef``; ``legs`` is
        ``list[kalshi_common.legset.CanonicalLeg]`` (the chosen sides).

        Route A (partition): if N <= ON_DEMAND_MAX_PARTITION_LEGS and every
        leg resolved two-sided, price all 2^N cells — cell i flips leg j to
        its opposite_ref iff bit j of i is set (legset.enumerate_partition
        ordering; cell 0 = target, priced FIRST). Any missing cell abandons
        the partition (the already-priced target is reused for Route B).

        Route B (correlation transfer): one SGP price for the target +
        per-leg devigged singles — structure odds where the book carries
        them, 1-leg price calls where it doesn't (DK), vig-fallback haircut
        when a single is one-sided.

        Failure accounting: a TRANSPORT-level exception (client/network
        raise) feeds ``_book_done(book, None)`` — the same 3-strike client
        reinit path refresh() uses, so on-demand and sweep health share one
        recovery. A clean "book won't price this" (event/resolve miss, None
        prices, devig gate rejection) is NOT a failure: no accounting, just
        None. Success does NOT touch last_success — that timestamp belongs
        to the sweep's min_refresh_sec scheduling.
        """
        t0 = time.monotonic()
        if book not in self._state or not legs:
            return None
        st = self._state[book]
        try:
            hooks = (self._on_demand_hooks or {}).get(book)
            if hooks is None:
                self._ensure_client(book)
                hooks = self._book_on_demand_hooks(book)

            event = hooks["match_event"](st.client, game)
            if event is None:
                return None                      # book doesn't list the game
            structure = hooks["build_structure"](st.client, event, game)
            if structure is None:
                return None
            resolved = hooks["resolve"](structure, legs,
                                        game.home_team, game.away_team)
            if not resolved:
                return None                      # a chosen side is missing

            def _price(refs):
                return hooks["price"](st.client, list(refs), event)

            n = len(legs)
            n_cells_priced = 0
            cells = None
            target_dec = None

            # ---- Route A: full 2^N partition ---- #
            if (n <= ON_DEMAND_MAX_PARTITION_LEGS
                    and all(r.opposite_ref is not None for r in resolved)):
                cell_decs = []
                for i in range(2 ** n):
                    refs = [resolved[j].opposite_ref if (i >> j) & 1
                            else resolved[j].ref for j in range(n)]
                    d = _price(refs)
                    if i == 0:
                        target_dec = d           # reusable by Route B
                    if d is None:
                        cell_decs = None         # abandon the partition
                        break
                    n_cells_priced += 1
                    cell_decs.append(d)
                cells = cell_decs

            # ---- Route B inputs (only if the partition can't produce) --- #
            from kalshi_common import fair_value
            singles = None
            partition_ok = (cells is not None
                            and fair_value.devig_partition(cells, n) is not None)
            if not partition_ok:
                if target_dec is None:
                    target_dec = _price([r.ref for r in resolved])
                    if target_dec is not None:
                        n_cells_priced += 1      # the target IS cell 0
                if target_dec is None:
                    return None                  # book won't price the combo
                singles = self._route_b_singles(book, resolved, _price)
                if singles is None:
                    return None                  # a leg has no usable single

            res = fair_value.book_on_demand_fair(cells, target_dec, singles, n)
            if res is None:
                return None
            fair, route = res
            return OnDemandBookResult(
                book=book, fair=fair, route=route,
                n_cells_priced=n_cells_priced,
                latency_sec=time.monotonic() - t0)
        except Exception as e:
            log.warning("sgp_service: %s on-demand transport error: %s", book, e)
            self._book_done(book, None)
            return None

    def _route_b_singles(self, book, resolved, price_fn):
        """Per-leg (vigged_implied, devigged_fair) pairs for Route B.

        Per leg: both single decimals -> exact devig_two_way; only the
        chosen side -> ON_DEMAND_VIG_FALLBACK haircut; NEITHER (DK's
        structure carries no odds) -> each available side priced as a 1-leg
        wire call first (Route-B-only cost, per spec §4.2). Any leg that
        still has no usable chosen-side decimal fails the whole set (None).
        """
        from kalshi_common import fair_value
        vig = ON_DEMAND_VIG_FALLBACK.get(book, 0.06)
        singles = []
        for r in resolved:
            sd, od = r.single_decimal, r.opposite_decimal
            if sd is None and od is None:
                # DK-style: no odds in the structure — 1-leg price calls.
                sd = price_fn([r.ref])
                od = (price_fn([r.opposite_ref])
                      if r.opposite_ref is not None else None)
            try:
                sd = float(sd) if sd is not None else None
                od = float(od) if od is not None else None
            except (TypeError, ValueError):
                return None
            if sd is None or sd <= 1.0:
                return None                      # no usable chosen single
            if od is not None:
                pair = fair_value.devig_two_way(sd, od)
                if pair is None:
                    return None
                fair = pair[0]
            else:
                fair = (1.0 / sd) / (1.0 + vig)  # one-sided haircut
            if not (0.0 < fair < 1.0):
                return None
            singles.append((1.0 / sd, fair))
        return singles

    # ------------------------------------------------------------------ #
    # Per-book on-demand hooks (real implementations).                    #
    # Each book's dict reuses the SAME event-fetch + canonical matching   #
    # the book's price_sgps runs, and the SAME structure objects its      #
    # resolve_legs was TDD'd against. Event lists and per-event           #
    # structures are TTL-cached per book (DK/FD reuse the _dk_fetchers/   #
    # _fd_fetchers caches; PX/NV/MGM/CZR get equivalents here).           #
    # ------------------------------------------------------------------ #

    @staticmethod
    def _structure_caches(st: _BookState):
        """The (events, structure) TTL caches, created on first use for the
        books whose refresh path doesn't build them (PX/NV/MGM/CZR)."""
        ev = st.caches.setdefault("events", TTLCache(EVENTS_TTL_SEC))
        struct = st.caches.setdefault("structure", TTLCache(STRUCTURE_TTL_SEC))
        return ev, struct

    def _book_on_demand_hooks(self, book: str) -> dict:
        st = self._state[book]

        if book == "draftkings":
            from mlb_sgp import draftkings as mod
            self._structure_caches(st)      # guarantee the caches exist
            f = self._dk_fetchers(st)

            def match_event(client, game):
                import scraper_draftkings_sgp as legacy
                events = f["fetch_dk_events"](client.session)
                matched = legacy.match_events(events, _od_parlay_lines(game))
                return matched[0] if matched else None

            def build_structure(client, event, game):
                eid = event["dk_event_id"]
                nums = f["fetch_main_market_nums"](client.session, eid)
                sel = f["fetch_selection_ids"](client.session, eid, nums)
                return (sel or {}).get("fg") or None

            return {"match_event": match_event,
                    "build_structure": build_structure,
                    "resolve": mod.resolve_legs,
                    "price": lambda client, refs, event:
                        mod.price_selection_set(client, refs)}

        if book == "fanduel":
            from mlb_sgp import fanduel as mod
            self._structure_caches(st)      # guarantee the caches exist
            f = self._fd_fetchers(st)

            def match_event(client, game):
                import scraper_fanduel_sgp as legacy
                events = f["fetch_fd_events"](client.session)
                matched = legacy.match_events(events, _od_parlay_lines(game))
                return matched[0] if matched else None

            def build_structure(client, event, game):
                runners = f["fetch_event_runners"](
                    client.session, event["fd_event_id"],
                    event["fd_home"], event["fd_away"])
                return (runners or {}).get("fg") or None

            return {"match_event": match_event,
                    "build_structure": build_structure,
                    "resolve": mod.resolve_legs,
                    "price": lambda client, refs, event:
                        mod.price_selection_set(client, refs)}

        if book == "prophetx":
            from mlb_sgp import prophetx as mod
            ev_cache, struct_cache = self._structure_caches(st)

            def match_event(client, game):
                import scraper_prophetx_sgp as legacy
                events = ev_cache.get_or_fetch("px_events", client.list_events)
                # Same Event -> legacy-dict translation as price_sgps
                # (prophetx.py:258-268).
                px_events = [
                    {"px_event_id": e.event_id, "px_home": e.home_team,
                     "px_away": e.away_team,
                     "px_home_competitor_id": e.home_id,
                     "px_away_competitor_id": e.away_id,
                     "scheduled": e.start_time}
                    for e in events
                ]
                matched = legacy.match_events(px_events, _od_parlay_lines(game))
                return matched[0] if matched else None

            def build_structure(client, event, game):
                from scraper_prophetx_sgp import _verify_competitor_ids
                eid = event["px_event_id"]

                def _fetch():
                    # Same Market -> raw-dict translation as price_sgps
                    # (prophetx.py:299-311).
                    return [
                        {"id": m.market_id, "name": m.name,
                         "marketLines": m.selections, "outcomes": m.outcomes,
                         "selections": m.ml_selections}
                        for m in client.fetch_event_markets(eid)
                    ]
                markets = struct_cache.get_or_fetch(("px_markets", eid), _fetch)
                home_id = event["px_home_competitor_id"]
                away_id = event["px_away_competitor_id"]
                if not markets or not _verify_competitor_ids(
                        markets, home_id, away_id):
                    return None
                return {"event_id": eid, "markets": markets,
                        "home_id": home_id, "away_id": away_id}

            return {"match_event": match_event,
                    "build_structure": build_structure,
                    "resolve": mod.resolve_legs,
                    "price": lambda client, refs, event:
                        mod.price_selection_set(client, refs)}

        if book == "novig":
            from mlb_sgp import novig as mod
            ev_cache, struct_cache = self._structure_caches(st)

            def match_event(client, game):
                import scraper_novig_sgp as legacy
                events = ev_cache.get_or_fetch("nv_events", client.list_events)
                # Same Event -> legacy-dict translation as price_sgps
                # (novig.py:207-218).
                nv_events = [
                    {"nv_event_id": e.event_id, "nv_home": e.home_team,
                     "nv_away": e.away_team, "nv_home_sym": e.home_sym,
                     "nv_away_sym": e.away_sym, "scheduled": e.start_time}
                    for e in events
                ]
                matched = legacy.match_events(nv_events, _od_parlay_lines(game))
                return matched[0] if matched else None

            def build_structure(client, event, game):
                from scraper_novig_sgp import fetch_event_legs
                eid = event["nv_event_id"]

                def _fetch():
                    # The matched dict carries None target lines (see
                    # _od_parlay_lines), so fetch_event_legs skips its
                    # per-period parse and just returns the raw tree.
                    _, markets = fetch_event_legs(client.session, event)
                    return markets
                markets = struct_cache.get_or_fetch(("nv_markets", eid), _fetch)
                if not markets:
                    return None
                return mod.build_line_structure(
                    markets, event["nv_home_sym"], event["nv_away_sym"])

            return {"match_event": match_event,
                    "build_structure": build_structure,
                    "resolve": mod.resolve_legs,
                    "price": lambda client, refs, event:
                        mod.price_selection_set(client, refs)}

        if book == "betmgm":
            from mlb_sgp import betmgm as mod
            ev_cache, struct_cache = self._structure_caches(st)

            def match_event(client, game):
                events = ev_cache.get_or_fetch("mgm_events", client.list_events)
                matched = mod._match_events(events, [_od_target(game)])
                return matched.get(game.game_id)

            def build_structure(client, event, game):
                def _fetch():
                    markets = client.fetch_markets(event.event_id)
                    if not markets:
                        return None
                    return mod.parse_markets(
                        markets, event.home_team, event.away_team)
                parsed = struct_cache.get_or_fetch(
                    ("mgm_markets", event.event_id), _fetch)
                return (parsed or {}).get("FG") or None

            return {"match_event": match_event,
                    "build_structure": build_structure,
                    "resolve": mod.resolve_legs,
                    # MGM's price_selection_set takes fixture_id as a 3rd
                    # arg — threaded from the matched Event.
                    "price": lambda client, refs, event:
                        mod.price_selection_set(client, refs, event.event_id)}

        if book == "caesars":
            from mlb_sgp import caesars as mod
            ev_cache, struct_cache = self._structure_caches(st)

            def match_event(client, game):
                events = ev_cache.get_or_fetch("czr_events", client.list_events)
                matched = mod._match_events(events, [_od_target(game)])
                return matched.get(game.game_id)

            def build_structure(client, event, game):
                def _fetch():
                    ev = client.fetch_event(event.event_id)
                    return mod.parse_markets(ev) if ev else None
                parsed = struct_cache.get_or_fetch(
                    ("czr_event", event.event_id), _fetch)
                return (parsed or {}).get("FG") or None

            return {"match_event": match_event,
                    "build_structure": build_structure,
                    "resolve": mod.resolve_legs,
                    "price": lambda client, refs, event:
                        mod.price_selection_set(client, refs)}

        raise ValueError(f"unknown book {book!r}")


# ---------------------------------------------------------------------- #
# Module helpers: GameRef -> the two event-match input shapes             #
# ---------------------------------------------------------------------- #

def _od_parlay_lines(game) -> dict:
    """GameRef -> the legacy parlay-lines dict the DK/FD/PX/NV
    ``match_events`` helpers consume. Matching only reads home/away team +
    commence_time bucket; the fg_/f5_ target lines ride along as None
    (on-demand has no target grid — and NV's fetch_event_legs relies on the
    Nones to skip its per-period parse)."""
    return {game.game_id: {
        "home_team": game.home_team,
        "away_team": game.away_team,
        "commence_time": game.commence_time,
        "fg_spread_line": None, "fg_total_line": None,
        "f5_spread_line": None, "f5_total_line": None,
    }}


def _od_target(game) -> TargetLine:
    """GameRef -> a shape-only TargetLine for MGM/CZR ``_match_events``
    (they read game_id/home/away/commence_time; period+lines are inert)."""
    return TargetLine(game_id=game.game_id, home_team=game.home_team,
                      away_team=game.away_team,
                      commence_time=game.commence_time,
                      period="FG", spread=0.0, total=0.0)
