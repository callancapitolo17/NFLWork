"""One refresh of the fecta board (Wagerzon specials -> book partition fairs),
plus the pure sizing step that turns a board into EV and stakes.

Inputs:  wz.fetch_fecta_specials(); FanDuel and BetMGM over plain HTTP;
         DraftKings over plain HTTP (its SGP builder's price endpoint).
Outputs: a Board (returned, and a COPY published after every special/book
         through `publish`, so the page fills in while DraftKings works and
         never reads a board mid-write).
Side effects: APPENDS every (special, book) result to fecta_quotes.
Sizing (size_board) runs at read time with the current settings, so a new
bankroll or Kelly fraction applies without re-pricing.

Books: FanDuel and BetMGM price the trifectas over plain HTTP. DraftKings
prices the superfectas alone (user decision 2026-10-02: only DK lets "scores
first" into an SGP, so only DK prices them) — its trifecta-part partition,
then its scores-first share. The three books price in parallel, one thread per
book; each book's own calls stay sequential, so DK's 0.5 s pacing holds.
"""
from __future__ import annotations

import logging
import threading
from dataclasses import dataclass, field, replace
from datetime import datetime, time, timedelta, timezone
from typing import Callable

from nfl_specials import config, wz
from nfl_specials.books import BookGame, SgpBook, parse_start
from nfl_specials.dk_book import DraftKingsBook
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
# A special's game is in its Wagerzon week: Thursday night through Monday
# night around the week's Sunday. Books list future weeks too, so a team whose
# game is missing or started at a book would otherwise be priced off NEXT
# week's game. Monday night kicks off Tuesday ~00:15-01:15 UTC.
WEEK_STARTS_BEFORE_SUNDAY = timedelta(days=4)
WEEK_ENDS_AFTER_SUNDAY = timedelta(days=2, hours=12)
# Books' kickoff clocks for the same game differ by a minute or so.
SAME_GAME_KICKOFF_TOLERANCE = timedelta(hours=12)


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
    book_status: dict[str, str]     # book -> 'ok' | 'error: ...'
    progress_done: int = 0
    progress_total: int = 0
    # book -> {'done': n, 'total': n} price calls this refresh, for the page's progress bar
    book_progress: dict[str, dict[str, int]] = field(default_factory=dict)

    def snapshot(self) -> "Board":
        lines = [replace(line, book_fairs=dict(line.book_fairs)) for line in self.lines]
        book_progress = {name: dict(counts) for name, counts in self.book_progress.items()}
        return replace(self, lines=lines, book_status=dict(self.book_status), book_progress=book_progress)


@dataclass(frozen=True)
class Sizing:
    fair_prob: float | None
    ev: float | None
    kelly_stake: float              # fractional Kelly, before the budget
    recommended_stake: float        # fitted to the Wagerzon budget, whole dollars
    yields_to: int | None           # rotation of the better special on the same team


def open_books() -> tuple[dict[str, SgpBook], dict[str, str]]:
    books: dict[str, SgpBook] = {}
    status: dict[str, str] = {}
    constructors = {"FanDuel": FanDuelBook, "BetMGM": BetMgmBook, "DraftKings": DraftKingsBook}
    for name in BOOK_ORDER:
        try:
            books[name] = constructors[name]()
            status[name] = "ok"
        except Exception as exc:  # one dark book must not blank the board
            log.exception("%s failed to open", name)
            status[name] = f"error: {exc}"[:200]
    return books, status


def in_special_week(game: BookGame, special: wz.WzSpecial) -> bool:
    sunday = datetime.combine(special.week_date, time(0, 0), tzinfo=timezone.utc)
    kickoff = parse_start(game.game_start_time)
    return sunday - WEEK_STARTS_BEFORE_SUNDAY <= kickoff <= sunday + WEEK_ENDS_AFTER_SUNDAY


def _locate_game(fecta: Fecta, special: wz.WzSpecial, books: dict[str, SgpBook]) -> BookGame | None:
    """The special's game, from the first book (GAME_SOURCE_ORDER) whose next
    game for the team falls in the special's week."""
    for name in GAME_SOURCE_ORDER:
        if name in books:
            game = books[name].find_game(fecta.team)
            if game is not None and in_special_week(game, special):
                return game
    return None


def book_game(book: SgpBook, line: "FectaLine") -> BookGame | None:
    """`book`'s listing of the line's game: same two teams, same kickoff
    (within tolerance). None when the book lists another game for the team."""
    game = book.find_game(line.fecta.team)
    if game is None or {game.home, game.away} != {line.game.home, line.game.away}:
        return None
    gap = abs(parse_start(game.game_start_time) - parse_start(line.game.game_start_time))
    return game if gap <= SAME_GAME_KICKOFF_TOLERANCE else None


def size_board(board: Board, bankroll: float, kelly_fraction: float,
               budget: float | None, now: datetime | None = None,
               placed_ids: frozenset[int] = frozenset()) -> list[Sizing]:
    """Worst-case fair, EV and stakes per line (same order as board.lines).

    No stake for a line whose game has started (at `now`, default the clock;
    Wagerzon refuses it, so it must not take budget from later games), nor
    for any special of a (game, team) already bet — `placed_ids` are the
    wz_game_ids placed this week. Those bets are in Wagerzon's available, and
    a second stake on the same team would double one position.

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

    now = now or datetime.now(timezone.utc)
    held_teams = {(line.game.home, line.game.away, line.fecta.team) for line in board.lines
                  if line.special.wz_game_id in placed_ids and line.game is not None}
    best_by_team: dict[tuple, tuple[float, int]] = {}
    for index, (line, (fair, _ev, stake)) in enumerate(zip(board.lines, raw)):
        if stake <= 0 or line.game is None or parse_start(line.game.game_start_time) <= now:
            continue
        key = (line.game.home, line.game.away, line.fecta.team)
        if key in held_teams:
            continue
        growth = log_growth(fair, line.special.wz_american, stake, bankroll)
        if key not in best_by_team or growth > best_by_team[key][0]:
            best_by_team[key] = (growth, index)

    candidates = sorted(index for _growth, index in best_by_team.values())
    fitted = budgeted_stakes([(raw[i][0], board.lines[i].special.wz_american) for i in candidates],
                             bankroll, kelly_fraction, budget, config.WZ_MIN_STAKE, config.WZ_MAX_STAKE)
    recommended = dict(zip(candidates, fitted))

    sized = []
    for index, (line, (fair, ev, stake)) in enumerate(zip(board.lines, raw)):
        yields_to = None
        key = (line.game.home, line.game.away, line.fecta.team) if line.game is not None else None
        if stake > 0 and key in best_by_team and index not in recommended:
            winner_index = best_by_team[key][1]
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
    game = book_game(book, line)
    if game is None:
        return BookFair(book.name, None, None, None, 0, reason="this game is not listed")
    role = "home" if game.home == target.team else "away"
    return price_fecta_at_book(book, game, role, target)


def _price_scores_first(book: SgpBook, line: FectaLine) -> ScoresFirstShare:
    game = book_game(book, line)
    if game is None:
        return ScoresFirstShare(book.name, None, 0, reason="this game is not listed")
    role = "home" if game.home == line.fecta.team else "away"
    return scores_first_share(book, game, role, line.fecta)


def _price_lane(name: str, book: SgpBook, tasks: list[tuple[str, FectaLine]], board: Board,
                board_lock: threading.Lock, store: Store, publish: Callable[[Board], None]) -> None:
    """Price one book's tasks in order. Network calls run outside the lock;
    every board mutation, quote append and publish happens inside it."""
    for task, line in tasks:
        now = datetime.now(timezone.utc)
        error: str | None = None
        try:
            if task == "scores_first":
                result = _price_scores_first(book, line)
            else:
                result = _price_partition(book, line)
        except Exception as exc:  # transport failure on one book/special
            log.warning("%s pricing %s failed: %s", name, line.special.description, exc)
            error = f"error: {exc}"[:200]
            if task == "scores_first":
                result = ScoresFirstShare(name, None, 0, reason=error)
            else:
                result = BookFair(name, None, None, None, 0, reason=error)
        with board_lock:
            if error is not None:
                board.book_status[name] = error
            if task == "scores_first":
                line.sf_share = result
                row = _quote_row(line, now, book=name, quote_kind="scores_first_share", share=result)
            else:
                line.book_fairs[name] = result
                kind = "trifecta_part" if line.is_superfecta() else "full"
                row = _quote_row(line, now, book=name, quote_kind=kind, fair=result)
            store.append_quotes([row])
            board.progress_done += 1
            board.book_progress[name]["done"] += 1
            publish(board.snapshot())


def refresh_board(store: Store, publish: Callable[[Board], None]) -> Board:
    started_at = datetime.now(timezone.utc)
    books, book_status = open_books()

    lines = []
    for special in wz.fetch_fecta_specials():
        parsed = parse_fecta(special.description)
        if isinstance(parsed, ParseFailure):
            lines.append(FectaLine(special, None, "unparsed", note=parsed.reason))
            continue
        game = _locate_game(parsed, special, books)
        if game is None:
            lines.append(FectaLine(special, parsed, "no_game",
                                   note="no book lists this team's game in the specials' week (started?)"))
            continue
        lines.append(FectaLine(special, parsed, "pricing", game=game))

    board = Board(started_at, None, lines, book_status)
    pricing = [line for line in lines if line.status == "pricing"]
    work = [("partition", name, line) for name in TRIFECTA_BOOKS if name in books
            for line in pricing if not line.is_superfecta()]
    if SUPERFECTA_BOOK in books:
        for line in pricing:
            if line.is_superfecta():
                work += [("partition", SUPERFECTA_BOOK, line), ("scores_first", SUPERFECTA_BOOK, line)]
    board.progress_total = len(work)
    for _task, name, _line in work:
        board.book_progress.setdefault(name, {"done": 0, "total": 0})["total"] += 1
    publish(board.snapshot())

    # One lane per book, run side by side: each book's calls stay sequential
    # (DK paces itself 0.5 s apart; FD/MGM sessions are single-threaded), but
    # the books no longer wait for each other, so the refresh takes about as
    # long as the slowest book (DraftKings) instead of the sum of all three.
    lanes: dict[str, list[tuple[str, FectaLine]]] = {name: [] for name in books}
    for task, name, line in work:
        lanes[name].append((task, line))
    board_lock = threading.Lock()
    lane_crashes: list[tuple[str, Exception]] = []

    def run_lane(name: str, tasks: list[tuple[str, FectaLine]]) -> None:
        try:
            _price_lane(name, books[name], tasks, board, board_lock, store, publish)
        except Exception as exc:  # e.g. the quote append failed: surface it, never swallow
            log.exception("%s pricing lane crashed", name)
            lane_crashes.append((name, exc))

    threads = [threading.Thread(target=run_lane, args=(name, tasks), name=f"price-{name}", daemon=True)
               for name, tasks in lanes.items() if tasks]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    if lane_crashes:
        name, exc = lane_crashes[0]
        raise RuntimeError(f"{name} pricing lane crashed: {exc}") from exc

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
    board.finished_at = datetime.now(timezone.utc)
    publish(board.snapshot())
    return board
