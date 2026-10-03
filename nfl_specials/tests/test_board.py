from datetime import datetime, timezone

from nfl_specials.board import Board, FectaLine, size_board
from nfl_specials.books import BookGame
from nfl_specials.pricing import BookFair
from nfl_specials.special_parser import parse_fecta
from nfl_specials.wz import WzSpecial

GAME = BookGame("evt", home="SEA", away="LAC", game_start_time="2099-01-01T00:00:00Z")


def _line(rotation: int, description: str, wz_american: int, fair_prob: float | None) -> FectaLine:
    special = WzSpecial(wz_game_id=rotation, rotation=rotation, description=description,
                        wz_american=wz_american)
    fairs = {} if fair_prob is None else {"FanDuel": BookFair("FanDuel", fair_prob, None, 1.3, 8)}
    return FectaLine(special, parse_fecta(description), "priced", game=GAME, book_fairs=fairs)


def _board(*lines: FectaLine) -> Board:
    return Board(datetime.now(timezone.utc), None, list(lines), {})


def test_same_team_specials_share_one_stake():
    plain = _line(1, "CHARGERS TRIFECTA (1Q, 1H & GM)", 1600, 0.09)
    spread = _line(2, "CHARGERS TRIFECTA (1Q +½, 1H +4½ & GM +7½)", 400, 0.30)
    sized = size_board(_board(plain, spread), 1000, 0.25)
    staked = [s for s in sized if s.recommended_stake > 0]
    assert len(staked) == 1
    loser = next(s for s in sized if s.recommended_stake == 0)
    assert loser.yields_to in (1, 2) and loser.kelly_stake > 0


def test_opposite_teams_are_sized_independently():
    lac = _line(1, "CHARGERS TRIFECTA (1Q, 1H & GM)", 1600, 0.09)
    sea = _line(2, "SEAHAWKS TRIFECTA (1Q, 1H & GM)", 300, 0.40)
    sized = size_board(_board(lac, sea), 1000, 0.25)
    assert all(s.recommended_stake > 0 for s in sized)


def test_negative_ev_and_unpriced_get_no_stake():
    bad = _line(1, "SEAHAWKS TRIFECTA (1Q, 1H & GM)", 120, 0.37)
    unpriced = _line(2, "SEAHAWKS SUPERFECTA (SCR 1ST, 1Q, 1H & GM -7½)", 180, None)
    sized = size_board(_board(bad, unpriced), 1000, 0.25)
    assert sized[0].ev < 0 and sized[0].recommended_stake == 0
    assert sized[1].fair_prob is None and sized[1].recommended_stake == 0


def test_consensus_is_mean_of_books():
    line = _line(1, "CHARGERS TRIFECTA (1Q, 1H & GM)", 1600, 0.08)
    line.book_fairs["BetMGM"] = BookFair("BetMGM", 0.10, None, 1.4, 8)
    line.book_fairs["DraftKings"] = BookFair("DraftKings", None, None, None, 18, reason="declined cell")
    assert size_board(_board(line), 1000, 0.25)[0].fair_prob == 0.09


def test_superfecta_fair_is_trifecta_part_times_scores_first_share():
    from nfl_specials.pricing import ScoresFirstShare
    line = _line(1, "SEAHAWKS SUPERFECTA (SCR 1ST, 1Q, 1H & GM -7½)", 180, 0.40)
    assert size_board(_board(line), 1000, 0.25)[0].fair_prob is None   # no DK share yet
    line.sf_share = ScoresFirstShare("DraftKings", 0.6, 2)
    assert size_board(_board(line), 1000, 0.25)[0].fair_prob == 0.40 * 0.6


def test_split_books_are_flagged_not_dropped():
    falcons = _line(1, "FALCONS TRIFECTA (1Q, 1H & GM)", 490, 0.189)
    falcons.book_fairs["BetMGM"] = BookFair("BetMGM", 0.222, None, 1.3, 8)
    chargers = _line(2, "CHARGERS TRIFECTA (1Q, 1H & GM)", 1600, 0.094)
    chargers.book_fairs["BetMGM"] = BookFair("BetMGM", 0.081, None, 1.4, 8)
    lone = _line(3, "SAINTS TRIFECTA (1Q, 1H & GM)", 250, 0.24)
    sized = size_board(_board(falcons, chargers, lone), 1000, 0.25)
    assert [s.books_split for s in sized] == [True, False, False]
    assert sized[0].fair_prob is not None
