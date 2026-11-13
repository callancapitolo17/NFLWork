"""One refresh of the 1H parlay board (Wagerzon NFL + college 1st-half spread
+ total, every same-game pairing), and the pure sizing step that stakes it
after the fectas.

Inputs:  wz_halves.fetch_half_games(); FanDuel, BetMGM and DraftKings SGPs
         (plain HTTP) for each league; Wagerzon previews of every pairing.
Outputs: a HalfBoard (returned, and a COPY published after each game/book).
Side effects: APPENDS every (game, book, corner) fair and every Wagerzon
corner price to half_parlay_quotes.

A game's four corners (away or home spread x over or under) are the cells
of ONE 2 x 2 partition at each book: four SGP prices, probit-devigged
together, give all four fairs. Fair = the lowest among the books that priced
the partition, as for fectas.

A whole-number line (spread or total) is priced at each book off the 2 x 2
partitions at the half-point lines either side (push_grid): the fitted
cover/push x over/push/under grid gives each corner's chance that both legs
win plus what Wagerzon pays back on a push (the other leg alone, or the
stake). The corner's fair is then the EFFECTIVE fair at Wagerzon's price,
both-win chance + push return / parlay decimal (the same EV as the real
bet), read off each result's lowest chance across books.

Sizing (size_halves) stakes each game as one joint Kelly group over its 1H
results, with the real push payouts and the game's placed 1H parlays as fixed bets.
"""
from __future__ import annotations

import logging
import math
import re
import threading
from dataclasses import dataclass, field, replace
from datetime import datetime, timedelta, timezone
from typing import Callable

from nfl_specials import config, push_grid, wz_halves
from nfl_specials.board import BOOK_ORDER, SAME_GAME_KICKOFF_TOLERANCE
from nfl_specials.books import CFB, NFL, BookGame, SgpBook, parse_start
from nfl_specials.dk_book import DraftKingsBook
from nfl_specials.fd_book import FanDuelBook
from nfl_specials.mgm_book import BetMgmBook
from nfl_specials.pricing import (OutcomeBet, OutcomeGroup, american_to_decimal, budgeted_stakes, expected_value,
                                  price_partition)
from nfl_specials.special_parser import Leg
from nfl_specials.store import Store
from nfl_specials.wz_halves import CORNERS, WzHalfGame, WzParlayPrice

log = logging.getLogger("nfl_specials.halves_board")

LEAGUES = (NFL, CFB)
# A whole-number line's push slices rest on how one book prices the lines
# either side, which nothing in that book's own prices can check; a second
# book can. Fewer books priced = no fair, so no stake.
MIN_BOOKS_WHOLE_NUMBER = 2
# College matches by name prefix, so the kickoff must agree closely.
CFB_KICKOFF_TOLERANCE = timedelta(minutes=30)
# A book name whose next token after Wagerzon's name is one of these is another
# school ('TEXAS' vs 'Texas Tech', 'MIAMI' vs 'Miami (OH)'), not a mascot.
OTHER_SCHOOL_TOKENS = frozenset({"STATE", "TECH", "AND", "A", "CHRISTIAN", "INTERNATIONAL", "SOUTHERN",
                                 "NORTHERN", "EASTERN", "WESTERN", "CENTRAL", "ATLANTIC", "OH", "OHIO"})


def has_whole_number_line(game: WzHalfGame) -> bool:
    return push_grid.is_whole_number(game.away_spread) or push_grid.is_whole_number(game.total)


@dataclass(frozen=True)
class HalfBookFair:
    book: str
    corner_fairs: dict[str, float] | None   # corner -> devigged probability both legs win; None when unpriced
    overround: float | None                 # the largest of the partitions priced
    reason: str | None = None
    push_returns: dict[str, float] | None = None   # whole-number line: corner -> $ back per $1 from pushes
    push_cells: dict[tuple[str, str], float] | None = None   # whole-number line: the fitted grid's cells

    def outcome_cells(self) -> dict[tuple[str, str], float] | None:
        """(spread result, total result) -> chance, for every 1H result this
        book priced: the four corners on half-point lines, the push grid on a
        whole number."""
        if self.push_returns:
            return self.push_cells
        if self.corner_fairs is None:
            return None
        return {tuple(corner.split("_")): fair for corner, fair in self.corner_fairs.items()}

    def effective_fair(self, corner: str, wz_decimal: float | None) -> float | None:
        """The both-win chance plus the push return in units of the parlay
        payout: the binary chance with the real bet's EV at `wz_decimal`."""
        if self.corner_fairs is None:
            return None
        if not self.push_returns:
            return self.corner_fairs[corner]
        if wz_decimal is None:
            return None
        return self.corner_fairs[corner] + self.push_returns[corner] / wz_decimal


@dataclass
class HalfLine:
    game: WzHalfGame
    status: str                             # 'pricing' | 'priced' | 'unpriced'
    note: str | None = None
    book_fairs: dict[str, HalfBookFair] = field(default_factory=dict)
    wz_prices: dict[str, WzParlayPrice] = field(default_factory=dict)

    def book_fair(self, book: str, corner: str, wz_decimal: float | None = None) -> float | None:
        """`book`'s effective fair at `wz_decimal`, default Wagerzon's current price."""
        if wz_decimal is None:
            price = self.wz_prices.get(corner)
            wz_decimal = price.decimal if price else None
        return self.book_fairs[book].effective_fair(corner, wz_decimal)

    def fair_prob(self, corner: str, wz_decimal: float | None = None) -> float | None:
        """The corner's fair on the worst-case 1H results (worst_case_cells;
        a parlay is only ever backed): on half-point lines the lowest fair
        among books that priced the partition. On a whole-number line it is
        the effective fair, what the results pay back / `wz_decimal` (default
        Wagerzon's current price), so it is the EV the stake is sized on."""
        cells = self.worst_case_cells()
        if cells is None:
            return None
        if not has_whole_number_line(self.game):
            return cells[tuple(corner.split("_"))]
        if wz_decimal is None:
            price = self.wz_prices.get(corner)
            wz_decimal = price.decimal if price else None
        if wz_decimal is None:
            return None
        paid_back = sum(p * push_grid.cell_payout(self.game, corner, wz_decimal, cell) for cell, p in cells.items())
        return paid_back / wz_decimal

    def wz_american(self, corner: str) -> int | None:
        price = self.wz_prices.get(corner)
        return conservative_american(price.decimal) if price and price.decimal else None

    def worst_case_cells(self) -> dict[tuple[str, str], float] | None:
        """Each 1H result's LOWEST chance among the books that priced the
        game (MIN_BOOKS_WHOLE_NUMBER of them on a whole number). The cells sum
        to 1 or less; sizing puts what is left in a state where every corner
        loses, so no corner's EV beats its worst book's."""
        per_book = [cells for cells in (f.outcome_cells() for f in self.book_fairs.values()) if cells is not None]
        if not per_book or (has_whole_number_line(self.game) and len(per_book) < MIN_BOOKS_WHOLE_NUMBER):
            return None
        results = sorted(set().union(*per_book))
        return {result: min(cells.get(result, 0.0) for cells in per_book) for result in results}


@dataclass
class HalfBoard:
    started_at: datetime
    finished_at: datetime | None
    lines: list[HalfLine]
    book_status: dict[str, str]             # '<book> <league>' -> 'ok' | 'error: ...'
    progress_done: int = 0
    progress_total: int = 0

    def snapshot(self) -> "HalfBoard":
        lines = [replace(line, book_fairs=dict(line.book_fairs), wz_prices=dict(line.wz_prices))
                 for line in self.lines]
        return replace(self, lines=lines, book_status=dict(self.book_status))


@dataclass(frozen=True)
class HalfSizing:
    fair_prob: float | None
    ev: float | None
    kelly_stake: float                      # joint fractional Kelly with no budget or minimum
    recommended_stake: float                # fitted to what the fectas leave of the balance
    unlimited_stake: float = 0.0            # with no balance limit ($20 minimum kept)


@dataclass(frozen=True)
class HeldHalfParlay:
    """A 1H parlay the account already placed."""
    wz_game_id: int
    corner: str
    risk: float
    wz_american: int


def conservative_american(decimal: float) -> int:
    """Whole American odds that never pay more than `decimal` (sizing takes
    American odds; rounding toward the smaller payout keeps an edge honest)."""
    if decimal >= 2.0:
        return math.floor((decimal - 1.0) * 100 + 1e-9)
    return -math.ceil(100.0 / (decimal - 1.0) - 1e-9)


# --- game matching -----------------------------------------------------------
def college_name_tokens(name: str) -> tuple[str, ...]:
    """'Louisiana-Lafayette Ragin' Cajuns' -> ('LOUISIANA', 'LAFAYETTE', 'RAGIN', 'CAJUNS'); 'ST' reads 'STATE'."""
    text = re.sub(r"[^A-Z0-9 ]", " ", name.upper().replace("&", " AND ").replace("'", ""))
    return tuple("STATE" if token == "ST" else token for token in text.split())


def college_names_match(wz_name: str, book_name: str) -> bool:
    """Wagerzon's name is the book's name, or its start followed by a mascot
    (never by OTHER_SCHOOL_TOKENS)."""
    wz_tokens, book_tokens = college_name_tokens(wz_name), college_name_tokens(book_name)
    if not wz_tokens or book_tokens[:len(wz_tokens)] != wz_tokens:
        return False
    rest = book_tokens[len(wz_tokens):]
    return not rest or rest[0] not in OTHER_SCHOOL_TOKENS


def _kickoff_gap(book_game: BookGame, game: WzHalfGame) -> timedelta | None:
    """|book kickoff - Wagerzon kickoff|, or None when the book lists no start."""
    try:
        return abs(parse_start(book_game.game_start_time) - game.game_start_time)
    except (ValueError, TypeError, AttributeError):
        return None


def match_book_game(game: WzHalfGame, book_games: list[BookGame]) -> BookGame | None:
    """The book's listing of `game`, same home and away, or None.

    NFL: the same two nflverse teams with kickoffs within the fectas'
    SAME_GAME_KICKOFF_TOLERANCE. College: both names matching
    (college_names_match) with kickoffs within CFB_KICKOFF_TOLERANCE.
    Two candidates is no match (fail closed).
    """
    if game.league == NFL:
        tolerance = SAME_GAME_KICKOFF_TOLERANCE
        teams = [g for g in book_games if g.home == game.home and g.away == game.away]
    else:
        tolerance = CFB_KICKOFF_TOLERANCE
        teams = [g for g in book_games
                 if college_names_match(game.home, g.home) and college_names_match(game.away, g.away)]
    found = [g for g in teams if (gap := _kickoff_gap(g, game)) is not None and gap <= tolerance]
    return found[0] if len(found) == 1 else None


# --- pricing -----------------------------------------------------------------
def price_half_at_book(book: SgpBook, book_game: BookGame, game: WzHalfGame) -> HalfBookFair:
    """All four corner fairs of `game` at `book`: one 2 x 2 partition on
    half-point lines; on a whole-number line, the partitions either side
    fitted into the push grid (push_grid.fit_push_grid)."""
    partitions, overrounds = {}, []
    for away_line, total_line in push_grid.corner_bracket_pairs(game):
        priced = price_partition_at(book, book_game, away_line, total_line)
        if priced.corner_fairs is None:
            return priced
        partitions[(away_line, total_line)] = priced.corner_fairs
        overrounds.append(priced.overround)
    if len(partitions) == 1:
        return HalfBookFair(book.name, next(iter(partitions.values())), overrounds[0])
    fit = push_grid.fit_push_grid(partitions, game.away_spread, game.total)
    if fit.grid is None:
        return HalfBookFair(book.name, None, max(overrounds), reason=fit.reason)
    both_win = {corner: fit.grid.prob(*corner.split("_")) for corner in CORNERS}
    returns = {corner: push_grid.push_return(fit.grid, corner, game) for corner in CORNERS}
    return HalfBookFair(book.name, both_win, max(overrounds), push_returns=returns,
                        push_cells=dict(fit.grid.cells))


def price_partition_at(book: SgpBook, book_game: BookGame, away_line: float, total_line: float) -> HalfBookFair:
    """The four corner fairs of one 2 x 2 partition at these (half-point) lines."""
    spread = book.leg_market(book_game, "away", Leg("spread", "H1", away_line))
    if spread is None:
        return HalfBookFair(book.name, None, None, reason=f"no SGP 1H spread at {away_line:+g}")
    total = book.leg_market(book_game, "away", Leg("total", "H1", total_line))
    if total is None:
        return HalfBookFair(book.name, None, None, reason=f"no SGP 1H total at {total_line:g}")
    partition = price_partition(book, book_game, [spread, total])
    if partition.cell_fairs is None:
        return HalfBookFair(book.name, None, partition.overround, reason=partition.reason)
    fairs = {}
    for (spread_outcome, total_outcome), fair in zip(partition.cells, partition.cell_fairs):
        side = "away" if spread_outcome in spread.winning else "home"
        side_total = "over" if total_outcome in total.winning else "under"
        key = f"{side}_{side_total}"
        fairs[key] = fairs.get(key, 0.0) + fair
    return HalfBookFair(book.name, fairs, partition.overround)


def open_books() -> tuple[dict[tuple[str, str], SgpBook], dict[str, str]]:
    """One adapter per (book, league); a book that fails to open is skipped."""
    constructors = {"FanDuel": FanDuelBook, "BetMGM": BetMgmBook,
                    "DraftKings": lambda league: DraftKingsBook(proxy_url=config.DK_PROXY_URL, league=league)}
    books, status = {}, {}
    for name in BOOK_ORDER:
        for league in LEAGUES:
            label = f"{name} {league}"
            try:
                books[(name, league)] = constructors[name](league=league)
                status[label] = "ok"
            except Exception as exc:  # one dark book must not blank the board
                log.exception("%s failed to open", label)
                status[label] = f"error: {exc}"[:200]
    return books, status


def _quote_rows(line: HalfLine, quoted_at: datetime, book: str, fair: HalfBookFair | None = None,
                wz_prices: dict[str, WzParlayPrice] | None = None) -> list[dict]:
    game = line.game
    base = {"quoted_at": quoted_at, "league": game.league, "wz_game_id": game.wz_game_id,
            "away_team": game.away, "home_team": game.home, "game_start_time": game.game_start_time,
            "away_spread": game.away_spread, "total": game.total, "book": book}
    rows = []
    for corner in CORNERS:
        if fair is not None:
            rows.append({**base, "corner": corner, "overround": fair.overround, "wz_decimal": None,
                         "fair_prob": fair.corner_fairs[corner] if fair.corner_fairs else None,
                         "reason": fair.reason,
                         "push_return": fair.push_returns[corner] if fair.push_returns else None})
        else:
            price = wz_prices[corner]
            rows.append({**base, "corner": corner, "overround": None, "fair_prob": None,
                         "wz_decimal": price.decimal, "reason": price.error, "push_return": None})
    return rows


def _price_lane(name: str, tasks: list[tuple[SgpBook, BookGame, HalfLine]], board: HalfBoard,
                board_lock: threading.Lock, store: Store, publish: Callable[[HalfBoard], None]) -> None:
    """One book's games in order; network calls outside the lock. The lane's
    quote rows are appended once, at the end."""
    rows = []
    for book, book_game, line in tasks:
        now = datetime.now(timezone.utc)
        try:
            result = price_half_at_book(book, book_game, line.game)
        except Exception as exc:  # transport failure on one book/game
            log.warning("%s pricing 1H %s @ %s failed: %s", name, line.game.away, line.game.home, exc)
            result = HalfBookFair(name, None, None, reason=f"error: {exc}"[:200])
        with board_lock:
            line.book_fairs[name] = result
            board.progress_done += 1
            publish(board.snapshot())
        rows += _quote_rows(line, now, name, fair=result)
    store.append_half_quotes(rows)


def refresh_halves(store: Store, publish: Callable[[HalfBoard], None]) -> HalfBoard:
    started_at = datetime.now(timezone.utc)
    books, book_status = open_books()
    lines = []
    for game in sorted(wz_halves.fetch_half_games(), key=lambda g: g.game_start_time):
        if game.game_start_time <= started_at:
            continue
        lines.append(HalfLine(game, "pricing"))
    board = HalfBoard(started_at, None, lines, book_status)

    lanes: dict[str, list[tuple[SgpBook, BookGame, HalfLine]]] = {name: [] for name in BOOK_ORDER}
    for line in lines:
        if line.status != "pricing":
            continue
        for name in BOOK_ORDER:
            book = books.get((name, line.game.league))
            if book is None:
                continue
            book_game = match_book_game(line.game, book.games())
            if book_game is None:
                line.book_fairs[name] = HalfBookFair(name, None, None, reason="this game is not listed")
                continue
            lanes[name].append((book, book_game, line))
    board.progress_total = sum(len(tasks) for tasks in lanes.values())
    publish(board.snapshot())

    board_lock = threading.Lock()
    lane_crashes: list[tuple[str, Exception]] = []

    def run_lane(name: str, tasks) -> None:
        try:
            _price_lane(name, tasks, board, board_lock, store, publish)
        except Exception as exc:  # e.g. the quote append failed: surface it, never swallow
            log.exception("%s 1H pricing lane crashed", name)
            lane_crashes.append((name, exc))

    threads = [threading.Thread(target=run_lane, args=(name, tasks), name=f"halves-{name}", daemon=True)
               for name, tasks in lanes.items() if tasks]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    if lane_crashes:
        name, exc = lane_crashes[0]
        raise RuntimeError(f"{name} 1H pricing lane crashed: {exc}") from exc

    priced = [line for line in lines if line.status == "pricing"
              and any(f.corner_fairs for f in line.book_fairs.values())]
    wz_prices = wz_halves.price_corners([line.game for line in priced])
    now = datetime.now(timezone.utc)
    rows = []
    for line in priced:
        line.wz_prices = {corner: wz_prices[(line.game.wz_game_id, corner)] for corner in CORNERS}
        rows += _quote_rows(line, now, "Wagerzon", wz_prices=line.wz_prices)
    store.append_half_quotes(rows)

    for line in lines:
        if line.status != "pricing":
            continue
        if line in priced:
            books_priced = [b for b, f in line.book_fairs.items() if f.corner_fairs]
            if has_whole_number_line(line.game) and len(books_priced) < MIN_BOOKS_WHOLE_NUMBER:
                # Listed under "not priced" with the reason; fair_prob is None too.
                line.status = "unpriced"
                line.note = (f"whole-number line priced by {', '.join(books_priced)} only; "
                             f"needs {MIN_BOOKS_WHOLE_NUMBER} books")
                continue
            line.status = "priced"
            continue
        line.status = "unpriced"
        line.note = "; ".join(f"{b}: {f.reason}" for b, f in line.book_fairs.items() if f.reason) or "no book priced it"
    board.finished_at = datetime.now(timezone.utc)
    publish(board.snapshot())
    return board


# --- sizing ------------------------------------------------------------------
# The probability the books' worst-case cells leave unassigned: every corner loses there.
UNASSIGNED = ("unassigned", "unassigned")


def game_outcome_group(line: HalfLine, corners: list[str], held_parlays: list[HeldHalfParlay]
                       ) -> OutcomeGroup | None:
    """One game's 1H corners (candidates) and the 1H parlays already placed
    on it (held), as one OutcomeGroup over the game's 1H results: the four
    corners, plus the push results on a whole-number line, each at its
    worst-case chance, and the rest a state where every corner loses. None
    when no 1H result has a worst-case chance (the game cannot be sized)."""
    cells = line.worst_case_cells()
    if cells is None:
        return None
    game = line.game
    unassigned = 1.0 - sum(cells.values())
    if unassigned < -1e-9:   # every book's cells sum to 1, so their minimum cannot pass it
        log.error("1H %s @ %s: worst-case results sum to %.6f, not sized", game.away, game.home, 1.0 - unassigned)
        return None
    states = list(cells.items())
    if unassigned > 1e-12:
        states.append((UNASSIGNED, unassigned))

    def payoffs(corner: str, american: int) -> tuple[float, ...]:
        decimal = american_to_decimal(american)
        return tuple(push_grid.cell_payout(game, corner, decimal, cell) - 1.0 for cell, _p in states)

    held_risk: dict[str, float] = {}
    for bet in held_parlays:
        held_risk[bet.corner] = held_risk.get(bet.corner, 0.0) + bet.risk
    candidates = tuple(
        OutcomeBet(payoffs(corner, line.wz_american(corner)),
                   max_stake=max(0.0, config.WZ_HALF_PARLAY_MAX_STAKE - held_risk.get(corner, 0.0)))
        for corner in corners)
    held = tuple(OutcomeBet(payoffs(bet.corner, bet.wz_american), risk=bet.risk) for bet in held_parlays)
    return OutcomeGroup(tuple(p for _cell, p in states), candidates, held)


def size_halves(board: HalfBoard, bankroll: float, kelly_fraction: float, budget: float | None,
                held: tuple[HeldHalfParlay, ...] = ()) -> dict[tuple[int, str], HalfSizing]:
    """Fair, EV and stakes per (wz_game_id, corner).

    `budget` is what the fectas left of the Wagerzon balance (None = no cap).
    Each game is ONE group (game_outcome_group): its corners are exclusive
    outcomes of the 1H, sized jointly with the 1H parlays already placed on
    it. Fectas are not in the groups (Cal, 2026-11-13); they only take their
    cash first through `budget`. No stake on a started game.
    """
    now = datetime.now(timezone.utc)
    held_by_game: dict[int, list[HeldHalfParlay]] = {}
    for bet in held:
        held_by_game.setdefault(bet.wz_game_id, []).append(bet)

    keys, groups = [], []
    for line in board.lines:
        game = line.game
        if game.game_start_time <= now:
            continue
        corners = [c for c in CORNERS
                   if (fair := line.fair_prob(c)) is not None and (american := line.wz_american(c)) is not None
                   and expected_value(fair, american) > 0]
        if not corners:
            continue
        group = game_outcome_group(line, corners, held_by_game.get(game.wz_game_id, []))
        if group is not None:
            keys.append([(game.wz_game_id, corner) for corner in corners])
            groups.append(group)

    fitted = budgeted_stakes(groups, bankroll, kelly_fraction, budget, config.WZ_MIN_STAKE)
    unlimited = budgeted_stakes(groups, bankroll, kelly_fraction, None, config.WZ_MIN_STAKE)
    joint = budgeted_stakes(groups, bankroll, kelly_fraction, None, 0.0)
    stakes: dict[tuple[int, str], tuple[float, float, float]] = {}
    for group_keys, group_stakes in zip(keys, zip(fitted, unlimited, joint)):
        for key, reading in zip(group_keys, zip(*group_stakes)):
            stakes[key] = reading

    sized = {}
    for line in board.lines:
        for corner in CORNERS:
            fair, american = line.fair_prob(corner), line.wz_american(corner)
            ev = expected_value(fair, american) if fair is not None and american is not None else None
            rec, unl, kelly = stakes.get((line.game.wz_game_id, corner), (0.0, 0.0, 0.0))
            sized[(line.game.wz_game_id, corner)] = HalfSizing(fair, ev, kelly, rec, unl)
    return sized
