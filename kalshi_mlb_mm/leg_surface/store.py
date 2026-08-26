"""In-memory leg surface: the store the quote path reads.

Inputs: rows published by the ingest workers. Outputs: per-leg fairs.
Side effects: NONE — no DB, no network. The DuckDB mirror lives in db.py.

Why memory and not a table: the epic's definition of done is that a
cross-game RFQ prices with zero outbound requests, and a DuckDB open costs
~17ms (three separate incidents came from doing one per item in a hot loop).
Reading a table per RFQ would trade one bottleneck for another.

Concurrency: one worker thread per book, each publishing ONLY its own book's
rows, plus reader threads. Publication is a whole-book dict swap under a
lock, so a reader never sees a half-written pass — it sees the previous pass
or the new one.
"""
from __future__ import annotations

import threading
from dataclasses import dataclass
from datetime import datetime

from kalshi_common.legset import CanonicalLeg

# (game_id, period, market_type, line, side) — book is the outer dict key.
# game_id is the KALSHI event-ticker suffix (26AUG252138CLELAA), which is what
# CanonicalLeg.game_id carries and is unique per doubleheader game. It is NOT
# the Odds API game_id in mlb_target_lines, whose team-name-only resolution
# collapses a doubleheader onto one row.
SurfaceKey = tuple[str, str, str, float | None, str]


def surface_key(leg: CanonicalLeg) -> SurfaceKey:
    """The store key for one canonical leg. line stays None for moneyline."""
    return (leg.game_id, leg.period, leg.market_type,
            None if leg.line is None else round(float(leg.line), 2), leg.side)


@dataclass(frozen=True)
class SurfaceRow:
    """One leg's devigged fair at one book, plus the raw prices behind it."""
    book: str
    game_id: str                 # Kalshi event-ticker suffix
    game_start_time: datetime    # naive UTC first pitch, parsed from the suffix
    period: str                  # "FG" | "F5" | "I1"
    market_type: str             # "ml" | "spread" | "total"
    line: float | None           # signed home-perspective; None for ml
    side: str                    # "home"/"away" | "over"/"under"
    fair_prob: float
    raw_decimal: float           # this side, as the book posted it
    raw_decimal_opp: float       # the other side of the same rung
    raw_overround: float         # 1/raw + 1/raw_opp, BEFORE devig
    route: str                   # "structure" | "singles"
    built_at: datetime           # when the BOOK DATA was fetched (aware UTC)

    @property
    def key(self) -> SurfaceKey:
        return (self.game_id, self.period, self.market_type,
                None if self.line is None else round(float(self.line), 2),
                self.side)


class LegSurface:
    """(book, route) -> {SurfaceKey: SurfaceRow}, read back BY BOOK.

    Slices are keyed by route because FanDuel runs both: its structure route
    owns ml/spread and all of I1, its singles route owns FG/F5 totals (FD's
    SGP structure carries only its own main total line). Each pass publishes
    by REPLACING its slice, so one shared slice would have each FD pass erase
    the other's rows.

    Reads collapse to the book, because two routes at one book are one
    opinion, not two — counting them separately would let FanDuel satisfy
    MIN_AGREEING_BOOKS on its own. Route ownership is exclusive per
    (market_type, period), so a key lands in exactly one slice.
    """

    def __init__(self):
        self._lock = threading.Lock()
        self._by_slice: dict[tuple[str, str], dict[SurfaceKey, SurfaceRow]] = {}

    def publish(self, book: str, route: str, rows: list[SurfaceRow]) -> int:
        """Replace this (book, route) slice with the pass's rows.

        REPLACE, not merge: a rung a book stopped posting must disappear
        rather than linger at its last price forever. Staleness is then only
        ever a `built_at` question, which #99's age gate can answer — a
        merged slice would hide a dead rung behind a fresh-looking neighbour.
        """
        slice_ = {r.key: r for r in rows}
        with self._lock:
            self._by_slice[(book, route)] = slice_
        return len(slice_)

    def drop(self, book: str, route: str | None = None) -> None:
        """Forget one slice, or every slice of a book (worker shutdown)."""
        with self._lock:
            for key in [k for k in self._by_slice
                        if k[0] == book and (route is None or k[1] == route)]:
                self._by_slice.pop(key, None)

    def get(self, book: str, leg: CanonicalLeg) -> SurfaceRow | None:
        return self.book_fairs(leg).get(book)

    def book_fairs(self, leg: CanonicalLeg) -> dict[str, SurfaceRow]:
        """Every BOOK that has this leg, one row each.

        Rows of ANY age, deliberately: #99's age gate and #20's dispersion
        gate live above this, and a store that silently dropped stale rows
        would make their decline counts unreadable.
        """
        k = surface_key(leg)
        out: dict[str, SurfaceRow] = {}
        with self._lock:
            for (book, _route), slice_ in self._by_slice.items():
                row = slice_.get(k)
                if row is None:
                    continue
                current = out.get(book)
                # Exclusive route ownership means this should never collide.
                # If it ever does, the fresher price wins rather than
                # whichever slice iterated last.
                if current is None or row.built_at > current.built_at:
                    out[book] = row
        return out

    def snapshot(self) -> list[SurfaceRow]:
        """Every row, for the DuckDB mirror and the acceptance queries."""
        with self._lock:
            return [r for slice_ in self._by_slice.values()
                    for r in slice_.values()]

    def books(self) -> list[str]:
        with self._lock:
            return sorted({book for book, _route in self._by_slice})

    def row_count(self) -> int:
        with self._lock:
            return sum(len(s) for s in self._by_slice.values())
