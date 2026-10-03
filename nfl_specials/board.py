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

Order: FanDuel and BetMGM price the trifectas over plain HTTP first.
DraftKings goes last and prices the superfectas alone (user decision
2026-10-02: only DK lets "scores first" into an SGP, so only DK prices them)
— its trifecta-part partition, then its scores-first share — through a real
browser, paced, with the page reloaded every 5 calls.
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
from nfl_specials.pricing import (BookFair, ScoresFirstShare, budgeted_stakes, expected_value,
                                  kelly_stake, log_growth, price_fecta_at_book, scores_first_share,
                                  trifecta_part, worst_case_fair)
from nfl_specials.special_parser import Fecta, ParseFailure, parse_fecta
from nfl_specials.store import Store

log = logging.getLogger("nfl_specials.board")

# Trifectas price at the HTTP books; superfectas only at DraftKings, the one
# book that lets "scores first" into an SGP.
TRIFECTA_BOOKS = ("FanDuel", "BetMGM")
SUPERFECTA_BOOK = "DraftKings"
BOOK_ORDER = TRIFECTA_BOOKS + (SUPERFECTA_BOOK,)
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
    # Trifecta: the special's fair by book. Superfecta: DraftKings' fair for
    # its TRIFECTA PART, multiplied by sf_share at sizing time.
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
    # Wagerzon account -> available balance, read at the start and end of the
    # refresh and lowered by the app after each placement. None = unknown.
    available_balance: dict[str, float | None] = field(default_factory=dict)

    def snapshot(self) -> "Board":
        lines = [replace(line, book_fairs=dict(line.book_fairs)) for line in self.lines]
        return replace(self, lines=lines, book_status=dict(self.book_status),
                       available_balance=dict(self.available_balance))


@dataclass(frozen=True)
class Sizing:
    fair_prob: float | None
    ev: float | None
    kelly_stake: float              # fractional Kelly, before the budget
    recommended_stake: float        # fitted to the Wagerzon budget, whole dollars
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


def size_board(board: Board, bankroll: float, kelly_fraction: float,
               budget: float | None) -> list[Sizing]:
    """Worst-case fair, EV and stakes per line (same order as board.lines).

    Per (game, team) only the line with the best expected log growth is a
    candidate: a team's trifecta and superfecta mostly win together, so
    staking both is one oversized bet. Opposite teams' fectas exclude each
    other and stay separate candidates. The candidates' stakes are then
    fitted to `budget` (pricing.budgeted_stakes; None = no cap).
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

    candidates = sorted(index for _growth, index in best_by_team.values())
    fitted = budgeted_stakes([(raw[i][0], board.lines[i].special.wz_american) for i in candidates],
                             bankroll, kelly_fraction, budget, config.WZ_MIN_STAKE)
    recommended = dict(zip(candidates, fitted))

    sized = []
    for index, (line, (fair, ev, stake)) in enumerate(zip(board.lines, raw)):
        yields_to = None
        if stake > 0 and line.game is not None and index not in recommended:
            winner_index = best_by_team[(line.game.home, line.game.away, line.fecta.team)][1]
            yields_to = board.lines[winner_index].special.rotation
        sized.append(Sizing(fair, ev, stake, recommended.get(index, 0.0), yields_to))
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

    board = Board(started_at, None, lines, book_status, available_balance=wz.available_balances())
    pricing = [line for line in lines if line.status == "pricing"]
    work = [("partition", name, line) for name in TRIFECTA_BOOKS if name in books
            for line in pricing if not line.is_superfecta()]
    if SUPERFECTA_BOOK in books:
        for line in pricing:
            if line.is_superfecta():
                work += [("partition", SUPERFECTA_BOOK, line), ("scores_first", SUPERFECTA_BOOK, line)]
    board.progress_total = len(work)
    publish(board.snapshot())

    for task, name, line in work:
        book = books[name]
        now = datetime.now(timezone.utc)
        try:
            if task == "scores_first":
                line.sf_share = _price_scores_first(book, line)
            else:
                line.book_fairs[name] = _price_partition(book, line)
        except Exception as exc:  # transport failure on one book/special
            log.warning("%s pricing %s failed: %s", name, line.special.description, exc)
            reason = f"error: {exc}"[:200]
            board.book_status[name] = reason
            if task == "scores_first":
                line.sf_share = ScoresFirstShare(name, None, 0, reason=reason)
            else:
                line.book_fairs[name] = BookFair(name, None, None, None, 0, reason=reason)
        if task == "scores_first":
            row = _quote_row(line, now, book=name, quote_kind="scores_first_share", share=line.sf_share)
        else:
            kind = "trifecta_part" if line.is_superfecta() else "full"
            row = _quote_row(line, now, book=name, quote_kind=kind, fair=line.book_fairs[name])
        store.append_quotes([row])
        board.progress_done += 1
        publish(board.snapshot())

    for line in pricing:
        if line.fair_prob() is not None:
            line.status = "priced"
            continue
        line.status = "unpriced"
        reasons = [f"{b}: {f.reason}" for b, f in line.book_fairs.items() if f.reason]
        if line.is_superfecta() and SUPERFECTA_BOOK not in books:
            reasons.append(f"superfectas price at DraftKings only: {book_status.get(SUPERFECTA_BOOK)}")
        elif line.is_superfecta() and (line.sf_share is None or line.sf_share.share is None):
            reasons.append(f"DraftKings scores first: {line.sf_share.reason if line.sf_share else 'not priced'}")
        line.note = "; ".join(reasons) or "no book priced it"
    board.available_balance = wz.available_balances()   # bets may have settled meanwhile
    board.finished_at = datetime.now(timezone.utc)
    publish(board.snapshot())
    return board
