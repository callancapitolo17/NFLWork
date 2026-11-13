"""The shape every sportsbook adapter (dk_book, fd_book, mgm_book) fills in.

A fecta leg at a book is an OUTCOME GROUP: the mutually exclusive, exhaustive
selections of one market at one line, e.g. 1st Half 3-way = (SEA, Tie, LAC)
or 1st Quarter spread 0.5 = (SEA -0.5, LAC +0.5). Pricing every combination
of the legs' outcomes (the partition) and devigging across it is the only
way to strip a book's SGP margin, which is far larger than its single-leg
vig (pricing.py).
"""
from __future__ import annotations

import threading
import time
from collections import deque
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Callable, Hashable, Protocol, TypeVar

from nfl_specials.special_parser import Leg


@dataclass(frozen=True)
class Outcome:
    ref: Hashable        # book-specific selection reference, opaque to callers
    label: str           # human label, e.g. 'SEA -0.5 1H' / 'Tie 1Q'


@dataclass(frozen=True)
class LegMarket:
    winning: frozenset[Outcome]      # outcomes on which the leg wins: one, or
                                     # two for "+0.5" read off a 3-way market (team or Tie)
    group: tuple[Outcome, ...]       # exhaustive outcome set, includes `winning`

    @classmethod
    def single(cls, chosen: Outcome, group: tuple[Outcome, ...]) -> "LegMarket":
        return cls(winning=frozenset([chosen]), group=group)


def widen_with_tie(three_way: LegMarket | None) -> LegMarket | None:
    """Read a 3-way team market as "+0.5": the team's outcome OR the tie.
    (With integer scores, +0.5 in a period = the team does not lose it.)"""
    if three_way is None or len(three_way.group) != 3:
        return None
    tie = next((o for o in three_way.group if o.label.strip().lower() == "tie"), None)
    if tie is None:
        return None
    return LegMarket(winning=three_way.winning | {tie}, group=three_way.group)


def is_half_point(line: float | None) -> bool:
    return line is not None and abs(abs(line) - 0.5) < 1e-6


# The leagues an adapter can list. NFL games carry nflverse abbreviations;
# college games carry the book's own team names, matched to Wagerzon's by
# halves_board.match_book_game.
NFL = "NFL"
CFB = "CFB"
# A book's team name -> the BookGame team key (None = not a team of the league).
TeamKey = Callable[[str], "str | None"]


@dataclass(frozen=True)
class BookGame:
    book_event_id: str
    home: str            # NFL: nflverse abbreviation; CFB: the book's team name
    away: str
    game_start_time: str  # ISO 8601 UTC as the book reports it


def parse_start(game_start_time: str) -> datetime:
    """Books send ISO 8601 UTC with 'Z' and up to 7 fractional digits."""
    text = game_start_time.replace("Z", "+00:00")
    head, dot, rest = text.partition(".")
    if dot:
        digits = "".join(ch for ch in rest if ch.isdigit())
        zone = rest[len(digits):]
        text = f"{head}.{digits[:6]}{zone}"
    parsed = datetime.fromisoformat(text)
    return parsed if parsed.tzinfo else parsed.replace(tzinfo=timezone.utc)


def next_game(games: list[BookGame], team: str) -> BookGame | None:
    """The team's next game that has not started (books list future weeks too)."""
    now = datetime.now(timezone.utc)
    upcoming = [g for g in games if team in (g.home, g.away) and parse_start(g.game_start_time) > now]
    return min(upcoming, key=lambda g: parse_start(g.game_start_time), default=None)


T = TypeVar("T")


# Statuses that mean "the book refuses us", as opposed to "won't price this".
BLOCKED_STATUSES = frozenset({403, 429})


class BookBlocked(RuntimeError):
    """The book refused us (HTTP 403/429): no more calls to it this refresh."""


class CallWindow:
    """At most `max_calls` sends in any trailing `window_seconds`: a send that
    would go over waits until the oldest one ages out. Held by the app, not
    the refresh, because the book's own window spans refreshes (FanDuel's
    WAF: 100 implyBets per trailing 5 min from one address, measured from the
    VM 2026-11-13)."""

    def __init__(self, max_calls: int, window_seconds: float,
                 clock: Callable[[], float] = time.monotonic,
                 sleep: Callable[[float], None] = time.sleep) -> None:
        self.max_calls = max_calls
        self.window_seconds = window_seconds
        self._clock, self._sleep = clock, sleep
        self._lock = threading.Lock()
        self._sent_at: deque[float] = deque()

    def wait_turn(self) -> None:
        """Block until one more send fits in the window, then count it as sent now."""
        with self._lock:
            while True:
                now = self._clock()
                while self._sent_at and now - self._sent_at[0] >= self.window_seconds:
                    self._sent_at.popleft()
                if len(self._sent_at) < self.max_calls:
                    self._sent_at.append(now)
                    return
                self._sleep(self.window_seconds - (now - self._sent_at[0]))


class RequestPacer:
    """One price call at a time to one book, each sent at least
    `min_gap_seconds` after the previous one's answer, across every adapter
    of that book in one refresh (the fecta and 1H boards price side by
    side). Once any call through it (paced or `guard`ed) answers 403/429,
    every later call raises BookBlocked without sending: piling more calls
    into a block only extends it (FanDuel 2026-11-13: one 403, then 164
    more). With a `window`, price calls also keep to its budget. Unlike
    mlb_sgp.dk_client.PriceLane it never gives up on a wait and caches
    nothing; the adapters keep their own caches."""

    def __init__(self, book: str, min_gap_seconds: float, window: CallWindow | None = None) -> None:
        self.book = book
        self.min_gap_seconds = min_gap_seconds
        self.window = window
        self._lock = threading.Lock()
        self._last_call_done_at = float("-inf")
        self._halted: str | None = None

    def call(self, send: Callable[[], T], what: str) -> T:
        """Run `send()` (one price request, `what` names it) once the gap has
        passed; returns its response."""
        with self._lock:
            self._raise_if_halted()
            wait = self.min_gap_seconds - (time.monotonic() - self._last_call_done_at)
            if wait > 0:
                time.sleep(wait)
            if self.window is not None:
                self.window.wait_turn()
                self._raise_if_halted()   # a block may have landed while it waited
            try:
                return self._note_status(send(), what)
            finally:
                self._last_call_done_at = time.monotonic()

    def guard(self, send: Callable[[], T], what: str) -> T:
        """Run `send()` (a non-price request, e.g. a game's market page)
        unpaced, but not after a block, and halt on a blocked answer."""
        self._raise_if_halted()
        return self._note_status(send(), what)

    def _raise_if_halted(self) -> None:
        if self._halted is not None:
            raise BookBlocked(f"{self.book} blocked this refresh: {self._halted}")

    def _note_status(self, response: T, what: str) -> T:
        status = getattr(response, "status_code", None)
        if status in BLOCKED_STATUSES:
            self._halted = f"{what} HTTP {status}"
        return response


class SgpBook(Protocol):
    """One refresh's view of one book. Adapters cache event listings and
    market payloads on the instance, so build a fresh one per refresh."""

    name: str

    def games(self) -> list[BookGame]:
        """Every game the book lists for the adapter's league."""

    def find_game(self, team: str) -> BookGame | None:
        """The team's next not-yet-started game at the book (books.next_game)."""

    def leg_market(self, game: BookGame, role: str, leg: Leg) -> LegMarket | None:
        """`leg` for the team playing `role` ('home'/'away'), or None when
        the book does not post it as an SGP-eligible market. A 'total' leg
        is the game's Over at that line (the role is ignored)."""

    def price(self, game: BookGame, refs: tuple[Hashable, ...]) -> float | None:
        """Decimal odds of the same-game parlay of `refs`, or None when the
        book declines the combination. Raises on transport failure."""
