"""1H parlays: Wagerzon parsing, game matching, the
2 x 2 partition fair, worst case across books, and sizing after the fectas.
No network."""
from datetime import date, datetime, timedelta, timezone

import pytest
from scipy.stats import norm

from nfl_specials import halves_board
from nfl_specials.board import Board, FectaLine, Sizing
from nfl_specials.books import BookGame, LegMarket, Outcome
from nfl_specials.halves_board import (HalfBoard, HalfBookFair, HalfLine, HeldFecta, HeldHalfParlay,
                                       college_names_match, conservative_american, match_book_game,
                                       price_half_at_book, size_halves)
from nfl_specials.pricing import BookFair
from nfl_specials.special_parser import parse_fecta
from nfl_specials.wz import WzSpecial
from nfl_specials.wz_halves import (CORNERS, WzHalfGame, WzParlayPrice, corner_legs, parse_half_game)

KICKOFF = datetime(2099, 1, 1, 18, 0, tzinfo=timezone.utc)


def _game(**overrides) -> WzHalfGame:
    fields = dict(league="NFL", wz_game_id=6118382, away="WAS", home="NYG", game_start_time=KICKOFF,
                  away_spread=0.5, away_spread_american=-118, home_spread_american=-102,
                  total=20.5, over_american=-102, under_american=-118)
    fields.update(overrides)
    return WzHalfGame(**fields)


# --- Wagerzon ------------------------------------------------------------------
def test_parses_a_wagerzon_nfl_first_half_game():
    raw = {"idgm": 6118382, "idgmtyp": 15, "vtm": "1H WAS COMMANDERS", "htm": "1H NY GIANTS",
           "gmdt": "20261112", "gmtm": "20:15:00",
           "GameLines": [{"vsprdt": "0.5", "vsprdoddst": "-118", "hsprdt": "-0.5", "hsprdoddst": "-102",
                          "unt": "20.5", "ovt": "-20.5", "ovoddst": "-102", "unoddst": "-118"}]}
    assert parse_half_game(raw, "NFL") == WzHalfGame(
        league="NFL", wz_game_id=6118382, away="WAS", home="NYG",
        game_start_time=datetime(2026, 11, 13, 1, 15, tzinfo=timezone.utc),
        away_spread=0.5, away_spread_american=-118, home_spread_american=-102,
        total=20.5, over_american=-102, under_american=-118)


def test_college_names_stay_wagerzons_and_a_missing_total_is_skipped():
    raw = {"idgm": 6123854, "idgmtyp": 15, "vtm": "1H MEMPHIS", "htm": "1H SOUTH FLORIDA",
           "gmdt": "20261112", "gmtm": "20:00:00",
           "GameLines": [{"vsprdt": "-2.5", "vsprdoddst": "-120", "hsprdt": "2.5", "hsprdoddst": "100",
                          "unt": "24.5", "ovoddst": "-105", "unoddst": "-115"}]}
    game = parse_half_game(raw, "CFB")
    assert (game.away, game.home, game.away_spread, game.home_spread) == ("MEMPHIS", "SOUTH FLORIDA", -2.5, 2.5)
    raw["GameLines"][0]["unt"] = ""
    assert parse_half_game(raw, "CFB") is None
    raw["idgmtyp"] = 20   # a quarter line
    assert parse_half_game(raw, "CFB") is None


def test_corner_legs_follow_wagerzons_encoding():
    game = _game()
    assert corner_legs(game, "away_over") == [
        {"idgm": 6118382, "play": 0, "points": 0.5, "odds": -118},
        {"idgm": 6118382, "play": 2, "points": -20.5, "odds": -102}]
    assert corner_legs(game, "home_under") == [
        {"idgm": 6118382, "play": 1, "points": -0.5, "odds": -102},
        {"idgm": 6118382, "play": 3, "points": 20.5, "odds": -118}]


def test_conservative_american_never_pays_more_than_the_decimal():
    assert conservative_american(3.638) == 263
    assert conservative_american(3.64) == 264
    assert conservative_american(1.5) == -200
    assert conservative_american(1.45) == -223


# --- matching --------------------------------------------------------------------
def test_nfl_games_match_on_teams_and_kickoff():
    listed = [BookGame("next-week", "NYG", "WAS", "2099-01-08T18:00:00Z"),
              BookGame("this-week", "NYG", "WAS", "2099-01-01T18:01:00Z")]
    assert match_book_game(_game(), listed).book_event_id == "this-week"
    assert match_book_game(_game(), [BookGame("x", "WAS", "NYG", "2099-01-01T18:00:00Z")]) is None


def test_college_names_match_the_book_name_or_its_start():
    assert college_names_match("MEMPHIS", "Memphis Tigers")
    assert college_names_match("SOUTH FLORIDA", "South Florida")
    assert college_names_match("MICHIGAN ST", "Michigan State Spartans")
    assert not college_names_match("UL LAFAYETTE", "Louisiana-Lafayette Ragin Cajuns")
    assert not college_names_match("FLORIDA STATE", "Florida Gators")


def test_college_match_needs_close_kickoffs_and_fails_closed_on_two_candidates():
    game = _game(league="CFB", away="MEMPHIS", home="SOUTH FLORIDA")
    right = BookGame("right", "South Florida Bulls", "Memphis Tigers", "2099-01-01T18:00:00Z")
    late = BookGame("late", "South Florida Bulls", "Memphis Tigers", "2099-01-01T19:00:00Z")
    assert match_book_game(game, [right, late]).book_event_id == "right"
    # One team matching in its slot is enough when it is the only such listing...
    ulm = _game(league="CFB", away="UL LAFAYETTE", home="UL MONROE")
    one = BookGame("one", "ULM", "Louisiana", "2099-01-01T18:00:00Z")
    assert match_book_game(ulm, [one]) is None
    memphis_only = BookGame("m", "USF", "Memphis", "2099-01-01T18:00:00Z")
    assert match_book_game(game, [memphis_only]).book_event_id == "m"
    # ...but two such listings at the kickoff is no match.
    other = BookGame("o", "Tulane", "Memphis", "2099-01-01T18:10:00Z")
    assert match_book_game(game, [memphis_only, other]) is None


# --- pricing ---------------------------------------------------------------------
TRUE_CORNERS = {"away_over": 0.29, "away_under": 0.21, "home_over": 0.20, "home_under": 0.30}


class CorrelatedBook:
    """A 2 x 2 book whose SGP prices are the true corner probabilities
    shifted by `vig_z` in probit space (what probit devig removes exactly)."""

    def __init__(self, name: str, corners: dict[str, float], vig_z: float = 0.1,
                 has_total: bool = True, declined: frozenset = frozenset()):
        self.name, self.corners, self.vig_z = name, corners, vig_z
        self.has_total, self.declined = has_total, declined

    def leg_market(self, game, role, leg):
        if leg.kind == "spread":
            away, home = Outcome("away", "away"), Outcome("home", "home")
            return LegMarket.single(away, (away, home))
        if not self.has_total:
            return None
        over, under = Outcome("over", "over"), Outcome("under", "under")
        return LegMarket.single(over, (over, under))

    def price(self, game, refs):
        corner = f"{refs[0]}_{refs[1]}"
        if corner in self.declined:
            return None
        return 1.0 / norm.cdf(norm.ppf(self.corners[corner]) + self.vig_z)


BOOK_GAME = BookGame("evt", "NYG", "WAS", "2099-01-01T18:00:00Z")


def test_one_partition_gives_all_four_corner_fairs():
    fair = price_half_at_book(CorrelatedBook("FanDuel", TRUE_CORNERS), BOOK_GAME, _game())
    assert fair.reason is None and 1.1 < fair.overround < 1.5
    assert fair.corner_fairs == pytest.approx(TRUE_CORNERS, abs=1e-6)


def test_a_missing_total_or_declined_cell_prices_nothing():
    no_total = price_half_at_book(CorrelatedBook("FanDuel", TRUE_CORNERS, has_total=False), BOOK_GAME, _game())
    assert (no_total.corner_fairs, no_total.reason) == (None, "no SGP 1H total at 20.5")
    declined = price_half_at_book(CorrelatedBook("BetMGM", TRUE_CORNERS, declined=frozenset({"home_under"})),
                                  BOOK_GAME, _game())
    assert declined.corner_fairs is None and declined.reason == "declined cell: home + under"


def test_fair_is_the_worst_case_per_corner_across_books():
    low_over = {"away_over": 0.27, "away_under": 0.23, "home_over": 0.21, "home_under": 0.29}
    line = HalfLine(_game(), "priced", book_fairs={
        "FanDuel": HalfBookFair("FanDuel", TRUE_CORNERS, 1.2),
        "BetMGM": HalfBookFair("BetMGM", low_over, 1.25),
        "DraftKings": HalfBookFair("DraftKings", None, None, reason="no SGP 1H spread at +0.5")})
    assert {c: line.fair_prob(c) for c in CORNERS} == {
        "away_over": 0.27, "away_under": 0.21, "home_over": 0.20, "home_under": 0.29}


# --- sizing ----------------------------------------------------------------------
def _half_board(*lines: HalfLine) -> HalfBoard:
    return HalfBoard(datetime.now(timezone.utc), None, list(lines), {})


def _priced_line(game: WzHalfGame, fairs: dict[str, float], decimals: dict[str, float]) -> HalfLine:
    return HalfLine(game, "priced", book_fairs={"FanDuel": HalfBookFair("FanDuel", fairs, 1.2)},
                    wz_prices={c: WzParlayPrice(d, 500) for c, d in decimals.items()})


EVEN_PRICES = {c: 3.64 for c in CORNERS}


def test_only_positive_ev_corners_get_a_stake():
    line = _priced_line(_game(), TRUE_CORNERS, EVEN_PRICES)
    sized = size_halves(_half_board(line), 5000, 0.25, None)
    staked = {c for c in CORNERS if sized[(line.game.wz_game_id, c)].recommended_stake > 0}
    assert staked == {"away_over", "home_under"}          # 0.29 and 0.30 at 3.64 beat 1/3.64 = 0.275
    away_over = sized[(line.game.wz_game_id, "away_over")]
    assert away_over.fair_prob == 0.29 and away_over.ev == pytest.approx(0.29 * 3.64 - 1)


def test_a_parlay_next_to_the_teams_fecta_is_sized_smaller():
    line = _priced_line(_game(), TRUE_CORNERS, EVEN_PRICES)
    special = WzSpecial(wz_game_id=99, rotation=99, description="COMMANDERS TRIFECTA (1Q, 1H & GM)",
                        wz_american=500, week_date=date(2098, 12, 31))
    fecta_game = BookGame("evt", home="NYG", away="WAS", game_start_time="2099-01-01T18:00:00Z")
    fecta = FectaLine(special, parse_fecta(special.description), "priced", game=fecta_game,
                      book_fairs={"FanDuel": BookFair("FanDuel", 0.2, None, 1.3, 8)})
    fecta_board = Board(datetime.now(timezone.utc), None, [fecta], {})
    fecta_sizing = [Sizing(0.2, 0.2, 10.0, 10.0, None)]
    alone = size_halves(_half_board(line), 20000, 0.25, None)
    beside = size_halves(_half_board(line), 20000, 0.25, None, fecta_board, fecta_sizing)
    key, other_team = (line.game.wz_game_id, "away_over"), (line.game.wz_game_id, "home_under")
    assert 0 < beside[key].recommended_stake < alone[key].recommended_stake
    assert beside[other_team].recommended_stake == alone[other_team].recommended_stake


def test_a_held_fecta_the_board_cannot_price_blocks_that_team():
    line = _priced_line(_game(), TRUE_CORNERS, EVEN_PRICES)
    special = WzSpecial(wz_game_id=99, rotation=99, description="COMMANDERS TRIFECTA (1Q, 1H & GM)",
                        wz_american=500, week_date=date(2098, 12, 31))
    fecta_game = BookGame("evt", home="NYG", away="WAS", game_start_time="2099-01-01T18:00:00Z")
    fecta = FectaLine(special, parse_fecta(special.description), "unpriced", game=fecta_game)
    fecta_board = Board(datetime.now(timezone.utc), None, [fecta], {})
    sized = size_halves(_half_board(line), 5000, 0.25, None, fecta_board, [Sizing(None, None, 0.0, 0.0, None)],
                        held_fectas=(HeldFecta(99, 50.0, 500),))
    assert sized[(line.game.wz_game_id, "away_over")].recommended_stake == 0
    assert sized[(line.game.wz_game_id, "home_under")].recommended_stake > 0


def test_a_placed_parlay_shrinks_the_next_stake_on_that_side():
    line = _priced_line(_game(), TRUE_CORNERS, EVEN_PRICES)
    key = (line.game.wz_game_id, "away_over")
    before = size_halves(_half_board(line), 5000, 0.25, None)[key].recommended_stake
    after = size_halves(_half_board(line), 5000, 0.25, None,
                        held=(HeldHalfParlay(line.game.wz_game_id, "away_over", before, 264),))[key]
    assert after.recommended_stake < before / 2


def test_budget_and_minimum_apply_and_a_started_game_gets_nothing():
    line = _priced_line(_game(), TRUE_CORNERS, EVEN_PRICES)
    fitted = size_halves(_half_board(line), 5000, 0.25, 60.0)
    assert sum(s.recommended_stake for s in fitted.values()) <= 60.0
    assert all(s.recommended_stake == 0 or s.recommended_stake >= 20 for s in fitted.values())
    started = _priced_line(_game(game_start_time=datetime.now(timezone.utc) - timedelta(minutes=1)),
                           TRUE_CORNERS, EVEN_PRICES)
    assert all(s.recommended_stake == 0 for s in size_halves(_half_board(started), 5000, 0.25, None).values())


def test_refresh_prices_whole_numbers_too_and_each_partition_once(monkeypatch):
    priced, whole = _game(), _game(wz_game_id=2, total=24.0)
    monkeypatch.setattr(halves_board.wz_halves, "fetch_half_games", lambda: [priced, whole])
    monkeypatch.setattr(halves_board.wz_halves, "price_corners",
                        lambda games: {(g.wz_game_id, c): WzParlayPrice(3.64, 500) for g in games for c in CORNERS})

    class ListedBook(CorrelatedBook):
        def games(self):
            return [BOOK_GAME]

    monkeypatch.setattr(halves_board, "open_books", lambda: (
        {("FanDuel", "NFL"): ListedBook("FanDuel", TRUE_CORNERS)}, {"FanDuel NFL": "ok"}))

    class MemoryStore:
        rows: list = []

        def append_half_quotes(self, rows):
            self.rows += rows

    store = MemoryStore()
    board = halves_board.refresh_halves(store, lambda _board: None)
    assert [(l.game.wz_game_id, l.status, l.note) for l in board.lines] == [
        (6118382, "priced", None), (2, "priced", None)]
    assert board.lines[0].fair_prob("home_under") == pytest.approx(0.30, abs=1e-6)
    # The same prices either side of 24 mean no chance of landing on 24.
    assert board.lines[1].book_fairs["FanDuel"].push_returns == pytest.approx({c: 0.0 for c in CORNERS}, abs=1e-6)
    assert board.progress_done == board.progress_total == 2
    assert sorted({row["book"] for row in store.rows}) == ["FanDuel", "Wagerzon"]
    assert len(store.rows) == 16
