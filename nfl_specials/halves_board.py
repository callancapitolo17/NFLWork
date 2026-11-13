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
both-win chance + push return / parlay decimal: the same EV as the real
bet, and a binary bet with that EV swings more than the real one, so the
Kelly stake sized on it is the smaller, cautious one.
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
from nfl_specials.board import BOOK_ORDER, SAME_GAME_KICKOFF_TOLERANCE, Board, HeldPlacement, Sizing
from nfl_specials.books import CFB, NFL, BookGame, SgpBook, parse_start
from nfl_specials.dk_book import DraftKingsBook
from nfl_specials.fd_book import FanDuelBook
from nfl_specials.mgm_book import BetMgmBook
from nfl_specials.pricing import (OutcomeBet, OutcomeGroup, american_to_decimal, budgeted_stakes, chain_states,
                                  expected_value, price_partition)
from nfl_specials.push_grid import AWAY, HOME, PUSH
from nfl_specials.special_parser import Fecta, Leg
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
    # (spread result, total result) -> probability: the 4 corners on half-point
    # lines, the fitted push grid on a whole number. None = not kept.
    cells: dict[tuple[str, str], float] | None = None

    def outcome_cells(self) -> dict[tuple[str, str], float] | None:
        """Every 1H result this book priced and its chance."""
        if self.cells is not None:
            return self.cells
        if self.corner_fairs is None or self.push_returns:
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
        """The LOWEST (effective) fair among books that priced the partition
        (a parlay is only ever backed); a whole-number line needs
        MIN_BOOKS_WHOLE_NUMBER of them. `wz_decimal`: the price the push
        return is paid against, default Wagerzon's current one."""
        priced = [fair for fair in (self.book_fair(name, corner, wz_decimal) for name in self.book_fairs)
                  if fair is not None]
        if has_whole_number_line(self.game) and len(priced) < MIN_BOOKS_WHOLE_NUMBER:
            return None
        return min(priced) if priced else None

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
        corner_fairs = next(iter(partitions.values()))
        cells = {tuple(corner.split("_")): fair for corner, fair in corner_fairs.items()}
        return HalfBookFair(book.name, corner_fairs, overrounds[0], cells=cells)
    fit = push_grid.fit_push_grid(partitions, game.away_spread, game.total)
    if fit.grid is None:
        return HalfBookFair(book.name, None, max(overrounds), reason=fit.reason)
    both_win = {corner: fit.grid.prob(*corner.split("_")) for corner in CORNERS}
    returns = {corner: push_grid.push_return(fit.grid, corner, game) for corner in CORNERS}
    return HalfBookFair(book.name, both_win, max(overrounds), push_returns=returns, cells=dict(fit.grid.cells))


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


@dataclass(frozen=True)
class GameFecta:
    """A fecta stake in the game the 1H parlays are sized around (placed, or
    the fecta board's recommendation): fixed exposure."""
    team: str                          # nflverse abbreviation, one of the game's two teams
    fair_prob: float
    wz_american: int
    risk: float
    min_h1_margin: int | None          # smallest 1H margin (team minus opponent) it can win at; None = no 1H leg


def _game_key(game: WzHalfGame) -> tuple:
    """NFL: (home, away), the teams the fectas carry. College: the Wagerzon game."""
    if game.league == NFL:
        return (game.home, game.away)
    return (CFB, game.wz_game_id)


def fecta_min_h1_margin(fecta: Fecta) -> int | None:
    """The smallest 1H margin at which every 1H leg of `fecta` wins. A win
    leg needs 1 (a tied half loses). A whole-number spread leg is counted at
    its push margin, in case a push drops the leg rather than losing it."""
    margins = []
    for leg in fecta.legs:
        if leg.period != "H1":
            continue
        if leg.kind == "win":
            margins.append(1)
        elif leg.kind == "spread":
            margins.append(int(-leg.line) if push_grid.is_whole_number(leg.line) else math.floor(-leg.line) + 1)
    return max(margins) if margins else None


def min_cover_margin(game: WzHalfGame, team: str) -> int:
    """The smallest 1H margin at which `team` covers Wagerzon's 1H spread."""
    team_line = game.away_spread if team == game.away else -game.away_spread
    return math.floor(-team_line) + 1


def corner_payoff(game: WzHalfGame, corner: str, parlay_decimal: float, cell: tuple[str, str]) -> float:
    """$ won or lost per $1 on `corner` when the 1H lands in `cell`: both legs
    win pays the parlay; a pushed leg is dropped and the other leg pays its
    straight price (both pushed = stake back); a losing leg loses."""
    if cell == UNASSIGNED:
        return -1.0
    spread_side, total_side = corner.split("_")
    spread_result, total_result = cell
    if spread_result not in (spread_side, PUSH) or total_result not in (total_side, PUSH):
        return -1.0
    if spread_result == spread_side and total_result == total_side:
        return parlay_decimal - 1.0
    if spread_result == PUSH and total_result == PUSH:
        return 0.0
    if spread_result == spread_side:
        american = game.away_spread_american if spread_side == AWAY else game.home_spread_american
    else:
        american = game.over_american if total_side == "over" else game.under_american
    return american_to_decimal(american) - 1.0


def _fecta_levels(game: WzHalfGame, cells: dict[tuple[str, str], float], chain: list[GameFecta],
                  ) -> Callable[[tuple[str, str]], list[tuple[int, float]]]:
    """How one team's fectas (a nested chain, likeliest first) split each 1H
    result: cell -> [(how many of the chain win, share of the cell)].

    When every fecta needs at least the 1H margin that covers Wagerzon's 1H
    spread (a 1H win leg on the 1H underdog, say), the fectas win only where
    the team covers: their chance sits in the team's covering cells, spread
    over them in proportion (the 1H total is taken as unrelated to the fecta
    given the cover). Otherwise, or when the fairs do not fit inside the
    cover chance, the fectas are taken as unrelated to the 1H result.
    """
    team_side = AWAY if chain[0].team == game.away else HOME
    cover_prob = sum(p for cell, p in cells.items() if cell[0] == team_side)
    fairs = [bet.fair_prob for bet in chain]
    needs_cover = all(bet.min_h1_margin is not None and bet.min_h1_margin >= min_cover_margin(game, bet.team)
                      for bet in chain)
    if needs_cover and 0 < fairs[0] <= cover_prob:
        inside = [(winners, p) for p, winners in chain_states([f / cover_prob for f in fairs])]
        return lambda cell: inside if cell[0] == team_side else [(0, 1.0)]
    unrelated = [(winners, p) for p, winners in chain_states(fairs)]
    return lambda cell: unrelated


def game_outcome_group(line: HalfLine, corners: list[str], held_parlays: list[HeldHalfParlay],
                       fectas: list[GameFecta]) -> OutcomeGroup | None:
    """One game's 1H corners (candidates), the 1H parlays already placed on it
    and both teams' fectas (held), as one OutcomeGroup over the game's 1H
    results x how many of each team's fectas win. None when no 1H result
    has a worst-case chance (the game cannot be sized).

    A team's fectas are one nested chain, as on the fecta board."""
    cells = line.worst_case_cells()
    if cells is None:
        return None
    game = line.game
    states: list[tuple[tuple[str, str], dict[str, int], float]] = [(cell, {}, p) for cell, p in cells.items()]
    unassigned = 1.0 - sum(cells.values())
    if unassigned > 1e-12:
        states.append((UNASSIGNED, {}, unassigned))

    chains: dict[str, list[GameFecta]] = {}
    for bet in fectas:
        if bet.team not in (game.away, game.home):
            raise ValueError(f"fecta team {bet.team} is not in {game.away} @ {game.home}")
        chains.setdefault(bet.team, []).append(bet)
    for team, chain in chains.items():
        chain.sort(key=lambda bet: -bet.fair_prob)
        levels = _fecta_levels(game, cells, chain)
        states = [(cell, {**winners_by_team, team: winners}, p * share)
                  for cell, winners_by_team, p in states
                  for winners, share in levels(cell) if p * share > 0]

    total = sum(p for _cell, _winners, p in states)
    state_probs = tuple(p / total for _cell, _winners, p in states)   # float dust only

    def parlay_payoffs(corner: str, american: int) -> tuple[float, ...]:
        decimal = american_to_decimal(american)
        return tuple(corner_payoff(game, corner, decimal, cell) for cell, _winners, _p in states)

    held = [OutcomeBet(parlay_payoffs(bet.corner, bet.wz_american), risk=bet.risk) for bet in held_parlays]
    for team, chain in chains.items():
        for rank, bet in enumerate(chain):
            net_odds = american_to_decimal(bet.wz_american) - 1.0
            payoffs = tuple(net_odds if rank < winners.get(team, 0) else -1.0 for _cell, winners, _p in states)
            held.append(OutcomeBet(payoffs, risk=bet.risk))
    held_risk: dict[str, float] = {}
    for bet in held_parlays:
        held_risk[bet.corner] = held_risk.get(bet.corner, 0.0) + bet.risk
    candidates = tuple(
        OutcomeBet(parlay_payoffs(corner, line.wz_american(corner)),
                   max_stake=max(0.0, config.WZ_HALF_PARLAY_MAX_STAKE - held_risk.get(corner, 0.0)))
        for corner in corners)
    return OutcomeGroup(state_probs, candidates, tuple(held))


def size_halves(board: HalfBoard, bankroll: float, kelly_fraction: float, budget: float | None,
                fecta_board: Board | None = None, fecta_sizing: list[Sizing] = (),
                held_fectas: tuple[HeldPlacement, ...] = (),
                held: tuple[HeldHalfParlay, ...] = ()) -> dict[tuple[int, str], HalfSizing]:
    """Fair, EV and stakes per (wz_game_id, corner) after the fectas.

    `budget` is what the fectas left of the Wagerzon balance (None = no cap).
    Each game is ONE group (game_outcome_group): its four corners are
    exclusive outcomes of the 1H (with the push results on a whole-number
    line), sized jointly with the 1H parlays already placed on it and, in
    the NFL, both teams' fectas in the game (placed, plus the stakes the
    fecta board recommends) as fixed bets. A held bet the boards cannot
    price blocks its game. No stake on a started game.
    """
    now = datetime.now(timezone.utc)
    fectas_by_game: dict[tuple, list[GameFecta]] = {}
    blocked: set[tuple] = set()

    if fecta_board is not None:
        placed_by_special: dict[int, list[HeldPlacement]] = {}
        for bet in held_fectas:
            placed_by_special.setdefault(bet.wz_game_id, []).append(bet)
        for line, sizing in zip(fecta_board.lines, fecta_sizing):
            if line.game is None or line.fecta is None:
                continue
            key = (line.game.home, line.game.away)
            placed = placed_by_special.get(line.special.wz_game_id, [])
            if sizing.fair_prob is None:
                if placed:
                    blocked.add(key)
                continue
            margin = fecta_min_h1_margin(line.fecta)
            stakes = [(bet.wz_american, bet.risk) for bet in placed]
            if sizing.recommended_stake > 0:
                stakes.append((line.special.wz_american, sizing.recommended_stake))
            fectas_by_game.setdefault(key, []).extend(
                GameFecta(line.fecta.team, sizing.fair_prob, american, risk, margin) for american, risk in stakes)

    line_by_id = {line.game.wz_game_id: line for line in board.lines}
    held_by_game: dict[int, list[HeldHalfParlay]] = {}
    for bet in held:
        if bet.wz_game_id in line_by_id:   # else a started game or one off the board
            held_by_game.setdefault(bet.wz_game_id, []).append(bet)

    group_ids, groups, group_corners = [], [], []
    for line in board.lines:
        game = line.game
        if game.game_start_time <= now or _game_key(game) in blocked:
            continue
        corners = [c for c in CORNERS
                   if (fair := line.fair_prob(c)) is not None and (american := line.wz_american(c)) is not None
                   and expected_value(fair, american) > 0]
        if not corners:
            continue
        group = game_outcome_group(line, corners, held_by_game.get(game.wz_game_id, []),
                                   fectas_by_game.get(_game_key(game), []))
        if group is None:
            continue
        group_ids.append(game.wz_game_id)
        groups.append(group)
        group_corners.append(corners)

    fitted = budgeted_stakes(groups, bankroll, kelly_fraction, budget, config.WZ_MIN_STAKE)
    unlimited = budgeted_stakes(groups, bankroll, kelly_fraction, None, config.WZ_MIN_STAKE)
    joint = budgeted_stakes(groups, bankroll, kelly_fraction, None, 0.0)
    stakes: dict[tuple[int, str], tuple[float, float, float]] = {}
    for gid, corners, group_fitted, group_unlimited, group_joint in zip(group_ids, group_corners, fitted,
                                                                         unlimited, joint):
        for corner, rec, unl, kelly in zip(corners, group_fitted, group_unlimited, group_joint):
            stakes[(gid, corner)] = (rec, unl, kelly)

    sized = {}
    for line in board.lines:
        for corner in CORNERS:
            fair, american = line.fair_prob(corner), line.wz_american(corner)
            ev = expected_value(fair, american) if fair is not None and american is not None else None
            rec, unl, kelly = stakes.get((line.game.wz_game_id, corner), (0.0, 0.0, 0.0))
            sized[(line.game.wz_game_id, corner)] = HalfSizing(fair, ev, kelly, rec, unl)
    return sized
