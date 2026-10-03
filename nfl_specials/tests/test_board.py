from datetime import date, datetime, timezone

import pytest

from nfl_specials import board as board_module
from nfl_specials.board import (Board, FectaLine, book_game, in_special_week, refresh_board,
                                 size_board)
from nfl_specials.books import BookGame
from nfl_specials.pricing import BookFair, ScoresFirstShare
from nfl_specials.special_parser import parse_fecta
from nfl_specials.wz import WzSpecial

GAME = BookGame("evt", home="SEA", away="LAC", game_start_time="2099-01-01T00:00:00Z")


def _line(rotation: int, description: str, wz_american: int, fair_prob: float | None) -> FectaLine:
    special = WzSpecial(wz_game_id=rotation, rotation=rotation, description=description,
                        wz_american=wz_american, week_date=date(2098, 12, 31))
    fairs = {} if fair_prob is None else {"FanDuel": BookFair("FanDuel", fair_prob, None, 1.3, 8)}
    return FectaLine(special, parse_fecta(description), "priced", game=GAME, book_fairs=fairs)


def _board(*lines: FectaLine) -> Board:
    return Board(datetime.now(timezone.utc), None, list(lines), {})


def test_same_team_specials_share_one_stake():
    plain = _line(1, "CHARGERS TRIFECTA (1Q, 1H & GM)", 1600, 0.09)
    spread = _line(2, "CHARGERS TRIFECTA (1Q +½, 1H +4½ & GM +7½)", 400, 0.30)
    sized = size_board(_board(plain, spread), 35000, 0.25, None)
    staked = [s for s in sized if s.recommended_stake > 0]
    assert len(staked) == 1
    loser = next(s for s in sized if s.recommended_stake == 0)
    assert loser.yields_to in (1, 2) and loser.kelly_stake > 0


def test_opposite_teams_are_sized_independently():
    lac = _line(1, "CHARGERS TRIFECTA (1Q, 1H & GM)", 1600, 0.09)
    sea = _line(2, "SEAHAWKS TRIFECTA (1Q, 1H & GM)", 300, 0.40)
    sized = size_board(_board(lac, sea), 35000, 0.25, None)
    assert all(s.recommended_stake > 0 for s in sized)


def test_negative_ev_and_unpriced_get_no_stake():
    bad = _line(1, "SEAHAWKS TRIFECTA (1Q, 1H & GM)", 120, 0.37)
    unpriced = _line(2, "SEAHAWKS SUPERFECTA (SCR 1ST, 1Q, 1H & GM -7½)", 180, None)
    sized = size_board(_board(bad, unpriced), 35000, 0.25, None)
    assert sized[0].ev < 0 and sized[0].recommended_stake == 0
    assert sized[1].fair_prob is None and sized[1].recommended_stake == 0


def test_fair_is_the_worst_case_book():
    line = _line(1, "CHARGERS TRIFECTA (1Q, 1H & GM)", 1600, 0.10)
    line.book_fairs["BetMGM"] = BookFair("BetMGM", 0.08, None, 1.4, 8)
    line.book_fairs["DraftKings"] = BookFair("DraftKings", None, None, None, 18, reason="declined cell")
    assert size_board(_board(line), 35000, 0.25, None)[0].fair_prob == 0.08


def test_superfecta_fair_is_trifecta_part_times_scores_first_share():
    line = _line(1, "SEAHAWKS SUPERFECTA (SCR 1ST, 1Q, 1H & GM -7½)", 180, 0.40)
    assert size_board(_board(line), 35000, 0.25, None)[0].fair_prob is None   # no DK share yet
    line.sf_share = ScoresFirstShare("DraftKings", 0.6, 2)
    assert size_board(_board(line), 35000, 0.25, None)[0].fair_prob == 0.40 * 0.6



def test_budget_goes_to_the_strongest_edges():
    raiders = _line(1, "RAIDERS SUPERFECTA (SCR 1ST, 1Q, 1H & GM)", 1730, 0.40)
    raiders.sf_share = ScoresFirstShare("DraftKings", 0.288, 20)
    chargers = _line(2, "CHARGERS TRIFECTA (1Q, 1H & GM)", 1600, 0.081)
    sized = size_board(_board(raiders, chargers), 35000, 0.25, 500)
    assert sum(s.recommended_stake for s in sized) <= 500
    assert sized[0].recommended_stake > sized[1].recommended_stake
    assert sized[0].kelly_stake > sized[0].recommended_stake   # the budget bound


class _ListingBook:
    name = "Fake"

    def __init__(self, game):
        self.game = game

    def find_game(self, team):
        return self.game


def test_book_game_must_be_the_same_game():
    line = _line(1, "SEAHAWKS TRIFECTA (1Q, 1H & GM)", 120, 0.37)
    same = BookGame("x", home="SEA", away="LAC", game_start_time="2099-01-01T00:01:00Z")
    next_week = BookGame("y", home="SEA", away="LAC", game_start_time="2099-01-08T00:00:00Z")
    other_team = BookGame("z", home="SEA", away="SF", game_start_time="2099-01-01T00:00:00Z")
    assert book_game(_ListingBook(same), line) is same
    assert book_game(_ListingBook(next_week), line) is None   # a rematch weeks later
    assert book_game(_ListingBook(other_team), line) is None


def test_special_week_is_thursday_through_monday_night():
    special = _line(1, "SEAHAWKS TRIFECTA (1Q, 1H & GM)", 120, 0.37).special   # Sunday 2098-12-31
    thursday = BookGame("t", "SEA", "LAC", "2098-12-28T01:15:00Z")
    monday_night = BookGame("m", "SEA", "LAC", "2099-01-02T01:15:00Z")
    next_thursday = BookGame("n", "SEA", "LAC", "2099-01-05T01:15:00Z")
    assert in_special_week(thursday, special) and in_special_week(monday_night, special)
    assert not in_special_week(next_thursday, special)


class _FakeStore:
    def __init__(self, fail_for_book: str | None = None) -> None:
        self.rows: list[dict] = []
        self.fail_for_book = fail_for_book

    def append_quotes(self, rows: list[dict]) -> None:
        if rows[0]["book"] == self.fail_for_book:
            raise RuntimeError("database is locked")
        self.rows.extend(rows)


def _patch_refresh(monkeypatch, descriptions: list[str]) -> None:
    specials = [WzSpecial(wz_game_id=i, rotation=i, description=d, wz_american=400,
                          week_date=date(2098, 12, 31)) for i, d in enumerate(descriptions, start=1)]
    books = {"FanDuel": object(), "BetMGM": object()}
    monkeypatch.setattr(board_module, "open_books", lambda url: (books, {name: "ok" for name in books}))
    monkeypatch.setattr(board_module.wz, "fetch_fecta_specials", lambda: specials)
    monkeypatch.setattr(board_module, "_locate_game", lambda fecta, special, books: GAME)
    monkeypatch.setattr(board_module, "_price_partition",
                        lambda book, line: BookFair("FanDuel" if book is books["FanDuel"] else "BetMGM",
                                                    0.1, None, 1.3, 8))


def test_refresh_prices_every_book_in_its_own_lane(monkeypatch):
    _patch_refresh(monkeypatch, ["CHARGERS TRIFECTA (1Q, 1H & GM)", "SEAHAWKS TRIFECTA (1Q, 1H & GM)"])
    store = _FakeStore()
    board = refresh_board(store, lambda snapshot: None)
    assert board.progress_done == board.progress_total == 4
    assert sorted((row["book"], row["rotation"]) for row in store.rows) == [
        ("BetMGM", 1), ("BetMGM", 2), ("FanDuel", 1), ("FanDuel", 2)]
    assert all(line.status == "priced" for line in board.lines)


def test_a_crashed_lane_fails_the_refresh(monkeypatch):
    _patch_refresh(monkeypatch, ["CHARGERS TRIFECTA (1Q, 1H & GM)"])
    with pytest.raises(RuntimeError, match="BetMGM pricing lane crashed"):
        refresh_board(_FakeStore(fail_for_book="BetMGM"), lambda snapshot: None)


def test_refresh_counts_progress_per_book(monkeypatch):
    specials = [WzSpecial(wz_game_id=r, rotation=r, description=d, wz_american=900, week_date=date(2098, 12, 31))
                for r, d in [(1, "CHARGERS TRIFECTA (1Q, 1H & GM)"),
                             (2, "CHARGERS SUPERFECTA (SCR 1ST, 1Q, 1H & GM)")]]
    books = {"FanDuel": object(), "BetMGM": object(), "DraftKings": object()}
    monkeypatch.setattr(board_module, "open_books", lambda: (books, {name: "ok" for name in books}))
    monkeypatch.setattr(board_module.wz, "fetch_fecta_specials", lambda: specials)
    monkeypatch.setattr(board_module, "_locate_game", lambda fecta, special, books: GAME)
    monkeypatch.setattr(board_module, "_price_partition",
                        lambda book, line: BookFair("x", 0.1, None, 1.3, 8))
    monkeypatch.setattr(board_module, "_price_scores_first",
                        lambda book, line: ScoresFirstShare("DraftKings", 0.85, 2))

    class NoStore:
        def append_quotes(self, rows):
            pass

    published = []
    board_module.refresh_board(NoStore(), published.append)
    # Trifecta at FD + MGM; superfecta's partition + scores-first share at DK.
    assert published[0].book_progress == {"FanDuel": {"done": 0, "total": 1}, "BetMGM": {"done": 0, "total": 1},
                                           "DraftKings": {"done": 0, "total": 2}}
    assert published[-1].book_progress == {"FanDuel": {"done": 1, "total": 1}, "BetMGM": {"done": 1, "total": 1},
                                            "DraftKings": {"done": 2, "total": 2}}
    # Snapshots are copies: each publish holds the count at its own moment
    # (the lanes run in parallel, so which book lands first is not fixed);
    # the last publish is the finished board.
    assert [sum(c["done"] for c in snap.book_progress.values()) for snap in published] == [0, 1, 2, 3, 4, 4]


def test_a_started_game_takes_no_stake_or_budget():
    started = _line(1, "CHARGERS TRIFECTA (1Q, 1H & GM)", 1600, 0.09)
    later = _line(2, "SEAHAWKS TRIFECTA (1Q, 1H & GM)", 400, 0.30)
    later.game = BookGame("evt2", home="SEA", away="LAC", game_start_time="2099-01-02T00:00:00Z")
    after_first_kickoff = datetime(2099, 1, 1, 1, tzinfo=timezone.utc)
    alone = size_board(_board(later), 10_000, 0.25, 300, now=after_first_kickoff)
    sized = size_board(_board(started, later), 10_000, 0.25, 300, now=after_first_kickoff)
    assert sized[0].recommended_stake == 0 and sized[0].yields_to is None
    assert sized[1].recommended_stake == alone[0].recommended_stake > 0


def test_a_team_already_bet_gets_no_new_stake():
    placed = _line(1, "CHARGERS TRIFECTA (1Q, 1H & GM)", 1600, 0.09)
    same_team = _line(2, "CHARGERS TRIFECTA (1Q +½, 1H +4½ & GM +7½)", 400, 0.30)
    other_team = _line(3, "SEAHAWKS TRIFECTA (1Q, 1H & GM)", 400, 0.30)
    sized = size_board(_board(placed, same_team, other_team), 10_000, 0.25, None, placed_ids=frozenset({1}))
    assert [s.recommended_stake for s in sized[:2]] == [0, 0]
    assert [s.yields_to for s in sized[:2]] == [None, None]
    assert sized[2].recommended_stake > 0
