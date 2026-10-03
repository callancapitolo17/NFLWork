"""The shape every sportsbook adapter (dk_book, fd_book, mgm_book) fills in.

A fecta leg at a book is an OUTCOME GROUP: the mutually exclusive, exhaustive
selections of one market at one line, e.g. 1st Half 3-way = (SEA, Tie, LAC)
or 1st Quarter spread 0.5 = (SEA -0.5, LAC +0.5). Pricing every combination
of the legs' outcomes (the partition) and devigging across it is the only
way to strip a book's SGP margin, which is far larger than its single-leg
vig (pricing.py).
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Hashable, Protocol

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


@dataclass(frozen=True)
class BookGame:
    book_event_id: str
    home: str            # nflverse abbreviation
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


class SgpBook(Protocol):
    """One refresh's view of one book. Adapters cache event listings and
    market payloads on the instance, so build a fresh one per refresh."""

    name: str

    def find_game(self, team: str) -> BookGame | None:
        """The team's next not-yet-started game at the book (books.next_game)."""

    def leg_market(self, game: BookGame, role: str, leg: Leg) -> LegMarket | None:
        """`leg` for the team playing `role` ('home'/'away'), or None when
        the book does not post it as an SGP-eligible market."""

    def price(self, game: BookGame, refs: tuple[Hashable, ...]) -> float | None:
        """Decimal odds of the same-game parlay of `refs`, or None when the
        book declines the combination. Raises on transport failure."""
