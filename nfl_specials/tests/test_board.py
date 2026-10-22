from datetime import date, datetime, timedelta, timezone

import pytest

from nfl_specials import board as board_module
from nfl_specials.board import (Board, FectaLine, HeldPlacement, book_game, in_special_week, refresh_board,
                                 size_board, special_week_sunday)
from nfl_specials.books import BookGame
from nfl_specials.pricing import BookFair, ScoresFirstShare
from nfl_specials.special_parser import parse_fecta
from nfl_specials.wz import WzSpecial

GAME = BookGame("evt", home="SEA", away="LAC", game_start_time="2099-01-01T00:00:00Z")


def _line(rotation: int, description: str, wz_american: int, fair_prob: float | None,
          week_date: date = date(2098, 12, 31)) -> FectaLine:
    special = WzSpecial(wz_game_id=rotation, rotation=rotation, description=description,
                        wz_american=wz_american, week_date=week_date)
    fairs = {} if fair_prob is None else {"FanDuel": BookFair("FanDuel", fair_prob, None, 1.3, 8)}
    return FectaLine(special, parse_fecta(description), "priced", game=GAME, book_fairs=fairs)


def _board(*lines: FectaLine) -> Board:
    return Board(datetime.now(timezone.utc), None, list(lines), {})


def test_same_team_specials_are_sized_jointly():
    spread = _line(1, "CHARGERS TRIFECTA (1Q +½, 1H +4½ & GM +7½)", 400, 0.30)
    plain = _line(2, "CHARGERS TRIFECTA (1Q, 1H & GM)", 1600, 0.12)
    sized = size_board(_board(spread, plain), 5000, 0.25, None)
    alone = [size_board(_board(line), 5000, 0.25, None)[0] for line in (spread, plain)]
    # Both still add log growth on top of each other, so both get a stake...
    assert all(s.recommended_stake > 0 for s in sized)
    # ...but each is smaller than it would be alone: they win together.
    assert all(s.recommended_stake < a.recommended_stake for s, a in zip(sized, alone))


def test_a_same_team_special_that_adds_nothing_yields():
    first = _line(1, "CHARGERS TRIFECTA (1Q, 1H & GM)", 1600, 0.09)
    copy = _line(2, "CHARGERS TRIFECTA (1Q -½, 1H & GM)", 1600, 0.09)   # the same bet
    sized = size_board(_board(first, copy), 5000, 0.25, None)
    staked = [s for s in sized if s.recommended_stake > 0]
    assert len(staked) == 1
    loser = next(s for s in sized if s.recommended_stake == 0)
    assert loser.yields_to in (1, 2) and loser.kelly_stake == 0


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


@pytest.mark.parametrize("week_date", [date(2099, 1, 4), date(2099, 1, 1)])   # Sunday, Thursday
def test_special_week_is_thursday_through_monday_night(week_date):
    special = _line(1, "SEAHAWKS TRIFECTA (1Q, 1H & GM)", 120, 0.37, week_date=week_date).special
    games =[BookGame("p", "SEA", "LAC", "2098-12-26T01:15:00Z"),   # last Thursday night
             BookGame("t", "SEA", "LAC", "2099-01-02T01:15:00Z"),   # Thursday night
             BookGame("s", "SEA", "LAC", "2099-01-04T18:00:00Z"),   # Sunday
             BookGame("m", "SEA", "LAC", "2099-01-06T01:15:00Z"),   # Monday night
             BookGame("n", "SEA", "LAC", "2099-01-09T01:15:00Z")]   # next Thursday night
    assert [in_special_week(game, special) for game in games] == [False, True, True, True, False]


def test_special_week_sunday_is_the_first_sunday_on_or_after():
    monday = date(2099, 1, 5)
    monday_to_sunday = [monday + timedelta(days=offset) for offset in range(7)]
    sundays = [special_week_sunday(day) for day in monday_to_sunday]
    assert sundays == [date(2099, 1, 11)] * 7


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
    monkeypatch.setattr(board_module, "open_books", lambda: (books, {name: "ok" for name in books}))
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


SUPERFECTA = "CHARGERS SUPERFECTA (SCR 1ST, 1Q, 1H & GM)"
DK_DOWN = "DraftKings down: error: HTTP 403"


def _refresh_with_draftkings_down(monkeypatch, books: dict, fairs: dict[str, float | None]):
    specials = [WzSpecial(wz_game_id=1, rotation=1, description=SUPERFECTA, wz_american=900,
                          week_date=date(2098, 12, 31))]
    status = {name: "ok" for name in books} | {"DraftKings": "error: HTTP 403"}
    monkeypatch.setattr(board_module, "open_books", lambda: (books, status))
    monkeypatch.setattr(board_module.wz, "fetch_fecta_specials", lambda: specials)
    monkeypatch.setattr(board_module, "_locate_game", lambda fecta, special, books: GAME)

    def price_partition(book, line):
        name = next(n for n, b in books.items() if b is book)
        assert line.is_superfecta()
        if fairs[name] is None:
            return BookFair(name, None, None, None, 8, reason="declined cell")
        return BookFair(name, fairs[name], None, 1.3, 8)

    monkeypatch.setattr(board_module, "_price_partition", price_partition)
    store = _FakeStore()
    return refresh_board(store, lambda snapshot: None), store


def _expected_row(board, *, book, quote_kind, fair_prob=None, overround=None, n_cells, sf_share=None, reason=None):
    return {"quoted_at": board.started_at, "wz_game_id": 1, "rotation": 1, "description": SUPERFECTA,
            "team": "LAC", "prop_type": "SUPERFECTA", "home_team": "SEA", "away_team": "LAC",
            "game_start_time": GAME.game_start_time, "wz_american": 900, "book": book, "quote_kind": quote_kind,
            "fair_prob": fair_prob, "sgp_decimal": None, "overround": overround, "n_cells": n_cells,
            "sf_share": sf_share, "reason": reason}


def test_draftkings_down_prices_superfectas_at_the_trifecta_books_with_a_fixed_share(monkeypatch):
    board, store = _refresh_with_draftkings_down(monkeypatch, {"FanDuel": object(), "BetMGM": object()},
                                                 {"FanDuel": 0.12, "BetMGM": 0.10})
    line = board.lines[0]
    assert (line.status, line.note) == ("priced", None)
    assert line.sf_share == ScoresFirstShare("fixed", 0.776, 0, reason=DK_DOWN)
    assert line.fair_prob() == pytest.approx(0.10 * 0.776)
    # The fixed share is quoted at the refresh's start; the lanes stamp their own time.
    assert all(row["quoted_at"] >= board.started_at for row in store.rows)
    rows = sorted(({**row, "quoted_at": board.started_at} if row["book"] != "fixed" else row for row in store.rows),
                  key=lambda row: row["book"])
    assert rows == [
        _expected_row(board, book="BetMGM", quote_kind="trifecta_part", fair_prob=0.10, overround=1.3, n_cells=8),
        _expected_row(board, book="FanDuel", quote_kind="trifecta_part", fair_prob=0.12, overround=1.3, n_cells=8),
        _expected_row(board, book="fixed", quote_kind="scores_first_share", n_cells=0, sf_share=0.776,
                      reason=DK_DOWN)]
    assert board.book_progress == {"FanDuel": {"done": 1, "total": 1}, "BetMGM": {"done": 1, "total": 1}}


def test_draftkings_down_with_one_trifecta_book_declining_uses_the_other(monkeypatch):
    board, _store = _refresh_with_draftkings_down(monkeypatch, {"FanDuel": object(), "BetMGM": object()},
                                                  {"FanDuel": 0.12, "BetMGM": None})
    assert board.lines[0].status == "priced"
    assert board.lines[0].fair_prob() == pytest.approx(0.12 * 0.776)


def test_draftkings_and_the_trifecta_books_down_names_draftkings(monkeypatch):
    board, store = _refresh_with_draftkings_down(monkeypatch, {}, {})
    assert (board.lines[0].status, board.lines[0].note) == ("unpriced", DK_DOWN)
    assert store.rows == [_expected_row(board, book="fixed", quote_kind="scores_first_share", n_cells=0,
                                        sf_share=0.776, reason=DK_DOWN)]


def test_a_started_game_takes_no_stake_or_budget():
    started = _line(1, "CHARGERS TRIFECTA (1Q, 1H & GM)", 1600, 0.09)
    later = _line(2, "SEAHAWKS TRIFECTA (1Q, 1H & GM)", 400, 0.30)
    later.game = BookGame("evt2", home="SEA", away="LAC", game_start_time="2099-01-02T00:00:00Z")
    after_first_kickoff = datetime(2099, 1, 1, 1, tzinfo=timezone.utc)
    alone = size_board(_board(later), 10_000, 0.25, 300, now=after_first_kickoff)
    sized = size_board(_board(started, later), 10_000, 0.25, 300, now=after_first_kickoff)
    assert sized[0].recommended_stake == 0 and sized[0].yields_to is None
    assert sized[1].recommended_stake == alone[0].recommended_stake > 0


def test_a_held_bet_shrinks_the_same_teams_next_stake():
    placed = _line(1, "CHARGERS TRIFECTA (1Q, 1H & GM)", 1600, 0.12)
    same_team = _line(2, "CHARGERS TRIFECTA (1Q +½, 1H +4½ & GM +7½)", 400, 0.30)
    other_team = _line(3, "SEAHAWKS TRIFECTA (1Q, 1H & GM)", 400, 0.30)
    board = _board(placed, same_team, other_team)
    fresh = size_board(board, 5_000, 0.25, None)
    sized = size_board(board, 5_000, 0.25, None, held=(HeldPlacement(1, 20.0, 1600),))
    spread_alone = size_board(_board(same_team), 5_000, 0.25, None)[0]
    assert sized[0].recommended_stake == pytest.approx(fresh[0].recommended_stake - 20, abs=1)   # topped up
    assert 0 < sized[1].recommended_stake < spread_alone.recommended_stake
    assert sized[2].recommended_stake == fresh[2].recommended_stake


def test_a_held_bet_tops_up_to_kelly_within_the_book_limit():
    placed = _line(1, "CHARGERS TRIFECTA (1Q, 1H & GM)", 1600, 0.081)
    kelly = size_board(_board(placed), 35_000, 0.25, None)[0].recommended_stake
    topped = size_board(_board(placed), 35_000, 0.25, None, held=(HeldPlacement(1, 25.0, 1600),))
    assert topped[0].recommended_stake == pytest.approx(kelly - 25, abs=1)
    rich = _line(1, "CHARGERS TRIFECTA (1Q, 1H & GM)", 1600, 0.15)
    capped = size_board(_board(rich), 35_000, 0.25, None, held=(HeldPlacement(1, 200.0, 1600),))
    assert capped[0].recommended_stake == 50.0      # Wagerzon's $250 per special, $200 already on


def test_a_held_bet_the_board_cannot_price_blocks_its_team():
    placed = _line(1, "CHARGERS TRIFECTA (1Q, 1H & GM)", 1600, None)
    same_team = _line(2, "CHARGERS TRIFECTA (1Q +½, 1H +4½ & GM +7½)", 400, 0.30)
    sized = size_board(_board(placed, same_team), 10_000, 0.25, None, held=(HeldPlacement(1, 100.0, 1600),))
    assert sized[1].recommended_stake == 0


def test_a_stake_the_minimum_drops_is_flagged():
    # The capped line hands its leftover Kelly to the other; $13 is under the $20 minimum.
    plain = _line(1, "CHARGERS TRIFECTA (1Q, 1H & GM)", 1910, 0.0804)
    spread = _line(2, "CHARGERS TRIFECTA (1Q +½, 1H +4½ & GM +7½)", 2125, 0.1048)
    sized = size_board(_board(plain, spread), 35000, 0.25, 5000)
    assert [(s.recommended_stake, s.unlimited_stake, s.under_min_stake) for s in sized] == \
        [(0.0, 0.0, True), (250.0, 250.0, False)]
    lone = size_board(_board(_line(3, "CHARGERS TRIFECTA (1Q, 1H & GM)", 1910, 0.051)), 35000, 0.25, 5000)
    assert [(s.recommended_stake, s.unlimited_stake, s.under_min_stake) for s in lone] == [(0.0, 0.0, True)]


def test_unlimited_stake_drops_the_balance_but_keeps_the_book_limits():
    raiders = _line(1, "RAIDERS TRIFECTA (1Q, 1H & GM)", 1730, 0.12)
    chargers = _line(2, "CHARGERS TRIFECTA (1Q, 1H & GM)", 1600, 0.081)
    sized = size_board(_board(raiders, chargers), 35000, 0.25, 200)
    assert sum(s.recommended_stake for s in sized) <= 200
    assert [s.unlimited_stake for s in sized] == [250.0, 206.0]   # Raiders capped at $250, Chargers at its Kelly
    assert all(s.unlimited_stake >= s.recommended_stake for s in sized)


def test_books_open_with_draftkings_through_the_configured_proxy(monkeypatch):
    opened = {}

    def fake(name):
        def construct(**kwargs):
            opened[name] = kwargs
            return name
        return construct

    monkeypatch.setattr(board_module, "FanDuelBook", fake("FanDuel"))
    monkeypatch.setattr(board_module, "BetMgmBook", fake("BetMGM"))
    monkeypatch.setattr(board_module, "DraftKingsBook", fake("DraftKings"))
    monkeypatch.setattr(board_module.config, "DK_PROXY_URL", "http://100.120.136.59:8097")
    books, status = board_module.open_books()
    assert books == {"FanDuel": "FanDuel", "BetMGM": "BetMGM", "DraftKings": "DraftKings"}
    assert status == {"FanDuel": "ok", "BetMGM": "ok", "DraftKings": "ok"}
    assert opened == {"FanDuel": {}, "BetMGM": {}, "DraftKings": {"proxy_url": "http://100.120.136.59:8097"}}
