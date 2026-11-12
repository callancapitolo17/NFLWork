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
the partition, as for fectas. A whole-number line is skipped
(wz_halves.whole_number_skip_reason).
"""
from __future__ import annotations

import logging
import math
import re
import threading
from dataclasses import dataclass, field, replace
from datetime import datetime, timedelta, timezone
from typing import Callable

from nfl_specials import config, wz_halves
from nfl_specials.board import (BOOK_ORDER, SAME_GAME_KICKOFF_TOLERANCE, Board, HeldPlacement, Sizing,
                                 _team_key)
from nfl_specials.books import CFB, NFL, BookGame, SgpBook, parse_start
from nfl_specials.dk_book import DraftKingsBook
from nfl_specials.fd_book import FanDuelBook
from nfl_specials.mgm_book import BetMgmBook
from nfl_specials.pricing import (Candidate, ExposureGroup, HeldBet, budgeted_stakes, expected_value,
                                  price_partition)
from nfl_specials.special_parser import Leg
from nfl_specials.store import Store
from nfl_specials.wz_halves import CORNERS, WzHalfGame, WzParlayPrice

log = logging.getLogger("nfl_specials.halves_board")

LEAGUES = (NFL, CFB)
# College matches by name prefix, so the kickoff must agree closely.
CFB_KICKOFF_TOLERANCE = timedelta(minutes=30)
# A book name whose next token after Wagerzon's name is one of these is another
# school ('TEXAS' vs 'Texas Tech', 'MIAMI' vs 'Miami (OH)'), not a mascot.
OTHER_SCHOOL_TOKENS = frozenset({"STATE", "TECH", "AND", "A", "CHRISTIAN", "INTERNATIONAL", "SOUTHERN",
                                 "NORTHERN", "EASTERN", "WESTERN", "CENTRAL", "ATLANTIC", "OH", "OHIO"})


@dataclass(frozen=True)
class HalfBookFair:
    book: str
    corner_fairs: dict[str, float] | None   # corner -> devigged probability; None when unpriced
    overround: float | None
    reason: str | None = None


@dataclass
class HalfLine:
    game: WzHalfGame
    status: str                             # 'pricing' | 'priced' | 'unpriced' | 'skipped'
    note: str | None = None
    book_fairs: dict[str, HalfBookFair] = field(default_factory=dict)
    wz_prices: dict[str, WzParlayPrice] = field(default_factory=dict)

    def fair_prob(self, corner: str) -> float | None:
        """The LOWEST fair among books that priced the partition (a parlay is only ever backed)."""
        priced = [b.corner_fairs[corner] for b in self.book_fairs.values() if b.corner_fairs]
        return min(priced) if priced else None

    def wz_american(self, corner: str) -> int | None:
        price = self.wz_prices.get(corner)
        return conservative_american(price.decimal) if price and price.decimal else None


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
    """All four corner fairs of `game` at `book` from one 2 x 2 partition."""
    spread = book.leg_market(book_game, "away", Leg("spread", "H1", game.away_spread))
    if spread is None:
        return HalfBookFair(book.name, None, None, reason=f"no SGP 1H spread at {game.away_spread:+g}")
    total = book.leg_market(book_game, "away", Leg("total", "H1", game.total))
    if total is None:
        return HalfBookFair(book.name, None, None, reason=f"no SGP 1H total at {game.total:g}")
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
                         "reason": fair.reason})
        else:
            price = wz_prices[corner]
            rows.append({**base, "corner": corner, "overround": None, "fair_prob": None,
                         "wz_decimal": price.decimal, "reason": price.error})
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
        skip = wz_halves.whole_number_skip_reason(game)
        lines.append(HalfLine(game, "skipped" if skip else "pricing", note=skip))
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
            line.status = "priced"
            continue
        line.status = "unpriced"
        line.note = "; ".join(f"{b}: {f.reason}" for b, f in line.book_fairs.items() if f.reason) or "no book priced it"
    board.finished_at = datetime.now(timezone.utc)
    publish(board.snapshot())
    return board


# --- sizing ------------------------------------------------------------------
def _group_key(game: WzHalfGame, corner: str) -> tuple:
    """NFL: (home, away, spread team), the key the fectas group by, so a
    parlay and that team's fectas in the game size as one group. College:
    the Wagerzon game and spread side (no college fectas)."""
    side = corner.split("_")[0]
    if game.league == NFL:
        return (game.home, game.away, game.away if side == "away" else game.home)
    return (CFB, game.wz_game_id, side)


def size_halves(board: HalfBoard, bankroll: float, kelly_fraction: float, budget: float | None,
                fecta_board: Board | None = None, fecta_sizing: list[Sizing] = (),
                held_fectas: tuple[HeldPlacement, ...] = (),
                held: tuple[HeldHalfParlay, ...] = ()) -> dict[tuple[int, str], HalfSizing]:
    """Fair, EV and stakes per (wz_game_id, corner) after the fectas.

    `budget` is what the fectas left of the Wagerzon balance (None = no cap).
    Each NFL group (one team's side of one game) holds, as fixed bets, that
    team's fectas in the game (the ones already placed and the stakes the
    fecta board recommends) and the 1H parlays already placed. The group is
    one nested chain (pricing.ExposureGroup): the maximum correlation the
    fairs allow, the cautious case, so a parlay next to a fecta gets a
    smaller stake, never a bigger one. A team's two corners (over and under)
    cannot both win, which the chain also over-counts, again on the cautious
    side. Corners on the two teams are separate groups. A held bet the
    boards cannot price blocks its group. No stake on a started game.
    """
    now = datetime.now(timezone.utc)
    held_by_key: dict[tuple, list[HeldBet]] = {}
    blocked: set[tuple] = set()

    if fecta_board is not None:
        risk_by_special: dict[int, list[HeldPlacement]] = {}
        for bet in held_fectas:
            risk_by_special.setdefault(bet.wz_game_id, []).append(bet)
        for line, sizing in zip(fecta_board.lines, fecta_sizing):
            if line.game is None or line.fecta is None:
                continue
            key = _team_key(line)
            placed = risk_by_special.get(line.special.wz_game_id, [])
            if sizing.fair_prob is None:
                if placed:
                    blocked.add(key)
                continue
            fixed = [HeldBet(sizing.fair_prob, bet.wz_american, bet.risk) for bet in placed]
            if sizing.recommended_stake > 0:
                fixed.append(HeldBet(sizing.fair_prob, line.special.wz_american, sizing.recommended_stake))
            held_by_key.setdefault(key, []).extend(fixed)

    line_by_id = {line.game.wz_game_id: line for line in board.lines}
    held_risk: dict[tuple[int, str], float] = {}
    for bet in held:
        held_risk[(bet.wz_game_id, bet.corner)] = held_risk.get((bet.wz_game_id, bet.corner), 0.0) + bet.risk
    for bet in held:
        line = line_by_id.get(bet.wz_game_id)
        if line is None:
            continue    # a started game or one off the board
        key = _group_key(line.game, bet.corner)
        fair = line.fair_prob(bet.corner)
        if fair is None:
            blocked.add(key)
            continue
        held_by_key.setdefault(key, []).append(HeldBet(fair, bet.wz_american, bet.risk))

    members: dict[tuple, list[tuple[int, str]]] = {}
    for line in board.lines:
        if line.game.game_start_time <= now:
            continue
        for corner in CORNERS:
            fair, american = line.fair_prob(corner), line.wz_american(corner)
            if fair is None or american is None or expected_value(fair, american) <= 0:
                continue
            key = _group_key(line.game, corner)
            if key not in blocked:
                members.setdefault(key, []).append((line.game.wz_game_id, corner))

    keys = list(members)
    groups = []
    for key in keys:
        candidates = tuple(Candidate(line_by_id[gid].fair_prob(c), line_by_id[gid].wz_american(c),
                                     max(0.0, config.WZ_HALF_PARLAY_MAX_STAKE - held_risk.get((gid, c), 0.0)))
                           for gid, c in members[key])
        groups.append(ExposureGroup(candidates, tuple(held_by_key.get(key, []))))
    fitted = budgeted_stakes(groups, bankroll, kelly_fraction, budget, config.WZ_MIN_STAKE)
    unlimited = budgeted_stakes(groups, bankroll, kelly_fraction, None, config.WZ_MIN_STAKE)
    joint = budgeted_stakes(groups, bankroll, kelly_fraction, None, 0.0)
    stakes: dict[tuple[int, str], tuple[float, float, float]] = {}
    for key, group_fitted, group_unlimited, group_joint in zip(keys, fitted, unlimited, joint):
        for member, rec, unl, kelly in zip(members[key], group_fitted, group_unlimited, group_joint):
            stakes[member] = (rec, unl, kelly)

    sized = {}
    for line in board.lines:
        for corner in CORNERS:
            fair, american = line.fair_prob(corner), line.wz_american(corner)
            ev = expected_value(fair, american) if fair is not None and american is not None else None
            rec, unl, kelly = stakes.get((line.game.wz_game_id, corner), (0.0, 0.0, 0.0))
            sized[(line.game.wz_game_id, corner)] = HalfSizing(fair, ev, kelly, rec, unl)
    return sized
