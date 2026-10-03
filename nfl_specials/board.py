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

Order: FanDuel and BetMGM price every trifecta, and every superfecta's
trifecta part, over plain HTTP first. DraftKings goes last and only prices
each superfecta's scores-first share (2 calls; pricing.scores_first_share):
its calls go through a real browser, paced, and DK denies a page after ~6.
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
from nfl_specials.pricing import (BookFair, ScoresFirstShare, expected_value,
                                  kelly_stake, log_growth, price_fecta_at_book, scores_first_share,
                                  trifecta_part, worst_case_fair)
from nfl_specials.special_parser import Fecta, ParseFailure, parse_fecta
from nfl_specials.store import Store

log = logging.getLogger("nfl_specials.board")

# Books that price trifectas (and superfectas' trifecta parts) by partition.
PARTITION_BOOKS = ("FanDuel", "BetMGM")
# The one book that lets "scores first" into an SGP.
SCORES_FIRST_BOOK = "DraftKings"
BOOK_ORDER = PARTITION_BOOKS + (SCORES_FIRST_BOOK,)
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
    # Trifecta: the special's fair by book. Superfecta: its TRIFECTA PART's
    # fair by book, multiplied by sf_share at sizing time.
    book_fairs: dict[str, BookFair] = field(default_factory=dict)
    sf_share: ScoresFirstShare | None = None

    def is_superfecta(self) -> bool:
        return self.fecta is not None and self.fecta.prop_type == "SUPERFECTA"

    def fair_prob(self) -> float | None:
        base = worst_case_fair(list(self.book_fairs.values()))
        if base is None or not self.is_superfecta():
            return base
        if self.sf_share is None or self.sf_share.share is None:
            return None
        return base * self.sf_share.share


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
    """Worst-case fair, EV and Kelly stake per line (same order as board.lines).

    Per (game, team) only the line with the best expected log growth keeps a
    recommended stake: a team's trifecta and superfecta mostly win together,
    so staking both is one oversized bet. Opposite teams' fectas exclude each
    other and are sized independently.
    """
    raw = []
    for line in board.lines:
        fair = line.fair_prob()
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


def _quote_row(line: FectaLine, quoted_at: datetime, *, book: str, quote_kind: str,
               fair: BookFair | None = None, share: ScoresFirstShare | None = None) -> dict:
    return {
        "quoted_at": quoted_at, "wz_game_id": line.special.wz_game_id,
        "rotation": line.special.rotation, "description": line.special.description,
        "team": line.fecta.team, "prop_type": line.fecta.prop_type,
        "home_team": line.game.home, "away_team": line.game.away,
        "game_start_time": line.game.game_start_time, "wz_american": line.special.wz_american,
        "book": book, "quote_kind": quote_kind,
        "fair_prob": fair.fair_prob if fair else None,
        "sgp_decimal": fair.sgp_decimal if fair else None,
        "overround": fair.overround if fair else None,
        "n_cells": fair.n_cells if fair else (share.n_calls if share else None),
        "sf_share": share.share if share else None,
        "reason": (fair.reason if fair else share.reason if share else None),
    }


def _price_partition(book: SgpBook, line: FectaLine) -> BookFair:
    """The special's fair at `book` — for a superfecta, its trifecta part's."""
    target = trifecta_part(line.fecta) if line.is_superfecta() else line.fecta
    game = book.find_game(target.team)
    if game is None:
        return BookFair(book.name, None, None, None, 0, reason="game not listed")
    role = "home" if game.home == target.team else "away"
    return price_fecta_at_book(book, game, role, target)


def _price_scores_first(book: SgpBook, line: FectaLine) -> ScoresFirstShare:
    game = book.find_game(line.fecta.team)
    if game is None:
        return ScoresFirstShare(book.name, None, 0, reason="game not listed")
    role = "home" if game.home == line.fecta.team else "away"
    return scores_first_share(book, game, role, line.fecta)


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
    pricing = [line for line in lines if line.status == "pricing"]
    work = [(name, line) for name in PARTITION_BOOKS if name in books for line in pricing]
    if SCORES_FIRST_BOOK in books:
        work += [(SCORES_FIRST_BOOK, line) for line in pricing if line.is_superfecta()]
    board.progress_total = len(work)
    publish(board.snapshot())

    for name, line in work:
        book = books[name]
        now = datetime.now(timezone.utc)
        try:
            if name == SCORES_FIRST_BOOK:
                line.sf_share = _price_scores_first(book, line)
                row = _quote_row(line, now, book=name, quote_kind="scores_first_share", share=line.sf_share)
            else:
                fair = _price_partition(book, line)
                line.book_fairs[name] = fair
                kind = "trifecta_part" if line.is_superfecta() else "full"
                row = _quote_row(line, now, book=name, quote_kind=kind, fair=fair)
        except Exception as exc:  # transport failure on one book/special
            log.warning("%s pricing %s failed: %s", name, line.special.description, exc)
            reason = f"error: {exc}"[:200]
            board.book_status[name] = reason
            if name == SCORES_FIRST_BOOK:
                line.sf_share = ScoresFirstShare(name, None, 0, reason=reason)
                row = _quote_row(line, now, book=name, quote_kind="scores_first_share", share=line.sf_share)
            else:
                line.book_fairs[name] = BookFair(name, None, None, None, 0, reason=reason)
                row = _quote_row(line, now, book=name, quote_kind="full", fair=line.book_fairs[name])
        store.append_quotes([row])
        board.progress_done += 1
        publish(board.snapshot())

    for line in pricing:
        if line.fair_prob() is not None:
            line.status = "priced"
            continue
        line.status = "unpriced"
        reasons = [f"{b}: {f.reason}" for b, f in line.book_fairs.items() if f.reason]
        if line.is_superfecta() and (line.sf_share is None or line.sf_share.share is None):
            why = line.sf_share.reason if line.sf_share else book_status.get(SCORES_FIRST_BOOK, "not priced")
            reasons.append(f"scores first needs DraftKings: {why}")
        line.note = "; ".join(reasons) or "no book priced it"
    board.finished_at = datetime.now(timezone.utc)
    publish(board.snapshot())
    return board
