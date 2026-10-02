"""One refresh of the fecta board (Wagerzon specials -> book partition fairs),
plus the pure sizing step that turns a board into EV and stakes.

Inputs:  wz.fetch_fecta_specials(); FanDuel and BetMGM over plain HTTP;
         DraftKings through dk_price_sidecar when it is running.
Outputs: a Board (returned, and a COPY published after every special/book
         through `publish`, so the page fills in while DraftKings works and
         never reads a board mid-write).
Side effects: APPENDS every (special, book) result to fecta_quotes.
Sizing (size_board) runs at read time with the current settings, so a new
bankroll or Kelly fraction applies without re-pricing.

Order: the fast HTTP books price every special first; DraftKings goes last
because its calls are paced ~1 s apart and a superfecta is 36 of them.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field, replace
from datetime import datetime, timezone
from typing import Callable

from nfl_specials import config, wz
from nfl_specials.books import BookGame, SgpBook
from nfl_specials.dk_book import DraftKingsBook, sidecar_is_up
from nfl_specials.fd_book import FanDuelBook
from nfl_specials.mgm_book import BetMgmBook
from nfl_specials.pricing import (BookFair, consensus_fair, expected_value, kelly_stake,
                                  log_growth, price_fecta_at_book)
from nfl_specials.special_parser import Fecta, ParseFailure, parse_fecta
from nfl_specials.store import Store

log = logging.getLogger("nfl_specials.board")

# Pricing order: fast HTTP books first, DraftKings (paced, slow) last.
BOOK_ORDER = ("FanDuel", "BetMGM", "DraftKings")
# Whose kickoff time the board shows: DK and BetMGM list the real kickoff,
# FanDuel lists it a minute late.
GAME_SOURCE_ORDER = ("DraftKings", "BetMGM", "FanDuel")


@dataclass
class FectaLine:
    special: wz.WzSpecial
    fecta: Fecta | None
    status: str                     # 'pricing' | 'priced' | 'unpriced' | 'unparsed' | 'no_game'
    note: str | None = None         # parse failure / why unpriced
    game: BookGame | None = None
    book_fairs: dict[str, BookFair] = field(default_factory=dict)


@dataclass
class Board:
    started_at: datetime
    finished_at: datetime | None
    lines: list[FectaLine]
    book_status: dict[str, str]     # book -> 'ok' | 'error: ...' | 'sidecar not running ...'
    progress_done: int = 0
    progress_total: int = 0

    def snapshot(self) -> "Board":
        lines = [replace(line, book_fairs=dict(line.book_fairs)) for line in self.lines]
        return replace(self, lines=lines, book_status=dict(self.book_status))


@dataclass(frozen=True)
class Sizing:
    fair_prob: float | None
    ev: float | None
    kelly_stake: float
    recommended_stake: float
    yields_to: int | None           # rotation of the better special on the same team


def open_books(sidecar_url: str) -> tuple[dict[str, SgpBook], dict[str, str]]:
    books: dict[str, SgpBook] = {}
    status: dict[str, str] = {}
    constructors = {"FanDuel": FanDuelBook, "BetMGM": BetMgmBook,
                    "DraftKings": lambda: DraftKingsBook(sidecar_url)}
    for name in BOOK_ORDER:
        if name == "DraftKings" and not sidecar_is_up(sidecar_url):
            status[name] = "sidecar not running (dk_price_sidecar/run.sh)"
            continue
        try:
            books[name] = constructors[name]()
            status[name] = "ok"
        except Exception as exc:  # one dark book must not blank the board
            log.exception("%s failed to open", name)
            status[name] = f"error: {exc}"[:200]
    return books, status


def _locate_game(fecta: Fecta, books: dict[str, SgpBook]) -> BookGame | None:
    for name in GAME_SOURCE_ORDER:
        if name in books:
            game = books[name].find_game(fecta.team)
            if game is not None:
                return game
    return None


def size_board(board: Board, bankroll: float, kelly_fraction: float) -> list[Sizing]:
    """Consensus fair, EV and Kelly stake per line (same order as board.lines).

    Per (game, team) only the line with the best expected log growth keeps a
    recommended stake: a team's trifecta and superfecta mostly win together,
    so staking both is one oversized bet. Opposite teams' fectas exclude each
    other and are sized independently.
    """
    raw = []
    for line in board.lines:
        fair = consensus_fair(list(line.book_fairs.values()))
        if fair is None:
            raw.append((None, None, 0.0))
            continue
        raw.append((fair, expected_value(fair, line.special.wz_american),
                    kelly_stake(fair, line.special.wz_american, bankroll, kelly_fraction)))

    best_by_team: dict[tuple, tuple[float, int]] = {}
    for index, (line, (fair, _ev, stake)) in enumerate(zip(board.lines, raw)):
        if stake <= 0 or line.game is None:
            continue
        key = (line.game.home, line.game.away, line.fecta.team)
        growth = log_growth(fair, line.special.wz_american, stake, bankroll)
        if key not in best_by_team or growth > best_by_team[key][0]:
            best_by_team[key] = (growth, index)

    sized = []
    for index, (line, (fair, ev, stake)) in enumerate(zip(board.lines, raw)):
        recommended, yields_to = 0.0, None
        if stake > 0 and line.game is not None:
            winner_index = best_by_team[(line.game.home, line.game.away, line.fecta.team)][1]
            if winner_index == index:
                recommended = float(round(stake))
            else:
                yields_to = board.lines[winner_index].special.rotation
        sized.append(Sizing(fair, ev, stake, recommended, yields_to))
    return sized


def _quote_row(line: FectaLine, fair: BookFair, quoted_at: datetime) -> dict:
    return {
        "quoted_at": quoted_at, "wz_game_id": line.special.wz_game_id,
        "rotation": line.special.rotation, "description": line.special.description,
        "team": line.fecta.team, "prop_type": line.fecta.prop_type,
        "home_team": line.game.home, "away_team": line.game.away,
        "game_start_time": line.game.game_start_time, "wz_american": line.special.wz_american,
        "book": fair.book, "fair_prob": fair.fair_prob, "sgp_decimal": fair.sgp_decimal,
        "overround": fair.overround, "n_cells": fair.n_cells, "reason": fair.reason,
    }


def refresh_board(store: Store, publish: Callable[[Board], None]) -> Board:
    started_at = datetime.now(timezone.utc)
    books, book_status = open_books(config.DK_SIDECAR_URL)

    lines = []
    for special in wz.fetch_fecta_specials():
        parsed = parse_fecta(special.description)
        if isinstance(parsed, ParseFailure):
            lines.append(FectaLine(special, None, "unparsed", note=parsed.reason))
            continue
        game = _locate_game(parsed, books)
        if game is None:
            lines.append(FectaLine(special, parsed, "no_game", note="no book lists this team's next game"))
            continue
        lines.append(FectaLine(special, parsed, "pricing", game=game))

    board = Board(started_at, None, lines, book_status)
    work = [(name, line) for name in BOOK_ORDER if name in books
            for line in lines if line.status == "pricing"]
    board.progress_total = len(work)
    publish(board.snapshot())

    for name, line in work:
        book = books[name]
        game = book.find_game(line.fecta.team)
        if game is None:
            fair = BookFair(name, None, None, None, 0, reason="game not listed")
        else:
            role = "home" if game.home == line.fecta.team else "away"
            try:
                fair = price_fecta_at_book(book, game, role, line.fecta)
            except Exception as exc:  # transport failure on one book/special
                log.warning("%s pricing %s failed: %s", name, line.special.description, exc)
                fair = BookFair(name, None, None, None, 0, reason=f"error: {exc}"[:200])
                board.book_status[name] = f"error: {exc}"[:200]
        line.book_fairs[name] = fair
        store.append_quotes([_quote_row(line, fair, datetime.now(timezone.utc))])
        board.progress_done += 1
        publish(board.snapshot())

    for line in lines:
        if line.status != "pricing":
            continue
        if consensus_fair(list(line.book_fairs.values())) is not None:
            line.status = "priced"
        else:
            line.status = "unpriced"
            line.note = "; ".join(f"{b}: {f.reason}" for b, f in line.book_fairs.items()) or "no book priced it"
    board.finished_at = datetime.now(timezone.utc)
    publish(board.snapshot())
    return board
