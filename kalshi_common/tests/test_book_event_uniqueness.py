"""The bot's match_event boundary must never accept an AMBIGUOUS book match.

Doubleheaders reach the per-book event matchers since the slate learned the
G1/G2 grammar. Those legacy matchers (DK/FD/PX/NV ``match_events``, MGM/CZR
``_match_events``) fall back to matching on TEAMS ALONE when either side lacks
a start time, and the bot's hooks then took ``matched[0]`` / the dict's single
entry. On a doubleheader that prices game 2's ladder under game 1's key — the
0.05-0.11 probability error #95 measured.

The guards live in ``kalshi_common.sgp_service`` (the bot's boundary), NOT in
the shared scrapers, so the dashboard pipeline keeps its team-only fallback.
"""
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone

from kalshi_common.leg_types import SCHEDULE_START_TOLERANCE_MIN
from kalshi_common.sgp_service import (_sole_book_event,
                                       _sole_dated_book_event)


@dataclass
class _Game:
    game_id: str
    commence_time: datetime | None


@dataclass
class _Event:
    event_id: str
    start_time: object


G1 = _Game("26SEP041410DETCLEG1", datetime(2026, 9, 4, 18, 10))


class TestListShapedHooks:
    """DK/FD/PX/NV: match_events returns one entry per matching BOOK event
    against our single-game dict, so uniqueness is a length check."""

    def test_one_match_is_returned(self):
        assert _sole_book_event([{"fd_event_id": "e1"}], book="fanduel",
                                game=G1) == {"fd_event_id": "e1"}

    def test_no_match_is_none(self):
        assert _sole_book_event([], book="fanduel", game=G1) is None

    def test_two_matches_fail_closed(self, caplog):
        # The doubleheader case: both book events passed the (skipped)
        # time guard on teams alone. The old code returned the first.
        import logging
        with caplog.at_level(logging.WARNING):
            got = _sole_book_event([{"fd_event_id": "game2"},
                                    {"fd_event_id": "game1"}],
                                   book="fanduel", game=G1)
        assert got is None
        assert any("book_event_ambiguous" in r.getMessage()
                   for r in caplog.records)


class TestDictShapedHooks:
    """MGM/CZR: _match_events keys on OUR game_id, so a second book event
    overwrites the first — the count cannot see it. The returned event's
    own clock has to agree with the game's instead."""

    def test_agreeing_start_is_returned(self):
        ev = _Event("mgm-1", "2026-09-04T18:09:00Z")          # 1 min skew
        assert _sole_dated_book_event({G1.game_id: ev}, book="betmgm",
                                      game=G1) is ev

    def test_datetime_start_is_accepted_too(self):
        ev = _Event("mgm-1", datetime(2026, 9, 4, 18, 15, tzinfo=timezone.utc))
        assert _sole_dated_book_event({G1.game_id: ev}, book="betmgm",
                                      game=G1) is ev

    def test_the_other_game_of_the_doubleheader_is_declined(self):
        # Game 2 (23:15Z) overwrote game 1 in the matcher's dict.
        ev = _Event("mgm-2", "2026-09-04T23:15:00Z")
        assert _sole_dated_book_event({G1.game_id: ev}, book="betmgm",
                                      game=G1) is None

    def test_a_missing_book_start_is_a_decline_not_an_accept(self):
        # The matcher's "one side lacks a timestamp — accept" branch.
        for raw in (None, "", "   "):
            ev = _Event("mgm-x", raw)
            assert _sole_dated_book_event({G1.game_id: ev}, book="betmgm",
                                          game=G1) is None

    def test_a_missing_game_start_is_a_decline(self):
        ev = _Event("mgm-1", "2026-09-04T18:10:00Z")
        undated = _Game(G1.game_id, None)
        assert _sole_dated_book_event({G1.game_id: ev}, book="betmgm",
                                      game=undated) is None

    def test_no_match_is_none(self):
        assert _sole_dated_book_event({}, book="betmgm", game=G1) is None

    def test_the_tolerance_is_the_schedule_tolerance(self):
        inside = G1.commence_time + timedelta(
            minutes=SCHEDULE_START_TOLERANCE_MIN - 1)
        outside = G1.commence_time + timedelta(
            minutes=SCHEDULE_START_TOLERANCE_MIN + 1)
        assert _sole_dated_book_event(
            {G1.game_id: _Event("a", inside.isoformat() + "Z")},
            book="caesars", game=G1) is not None
        assert _sole_dated_book_event(
            {G1.game_id: _Event("b", outside.isoformat() + "Z")},
            book="caesars", game=G1) is None

    def test_an_unparseable_start_is_a_decline(self):
        ev = _Event("mgm-bad", "not a date")
        assert _sole_dated_book_event({G1.game_id: ev}, book="betmgm",
                                      game=G1) is None
