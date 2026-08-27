"""Rung devig, exclusion accounting and the in-memory store (issue #96).

The devig is the ONLY local arithmetic the surface does, so it is the only
place the surface can invent a number. These tests fix the shapes it must
refuse, and the shapes the store must keep apart.
"""
from datetime import datetime, timedelta, timezone

import pytest

from kalshi_common.legset import CanonicalLeg
from kalshi_mlb_mm.leg_surface.devig import ExclusionCounts, devig_rung
from kalshi_mlb_mm.leg_surface.slate import SurfaceGame
from kalshi_mlb_mm.leg_surface.store import LegSurface, SurfaceRow

BUILT_AT = datetime(2026, 8, 26, 1, 30, tzinfo=timezone.utc)
GAME = SurfaceGame(game_id="26AUG252138CLELAA", home_team="Los Angeles Angels",
                   away_team="Cleveland Guardians",
                   start_utc=datetime(2026, 8, 26, 1, 38), legs=())


def rung(line=8.5, period="FG"):
    return [CanonicalLeg(GAME.game_id, "total", line, "over", period),
            CanonicalLeg(GAME.game_id, "total", line, "under", period)]


def devig(decimals, legs=None):
    return devig_rung(book="betmgm", route="structure", game=GAME,
                      legs=legs or rung(), decimals=decimals,
                      built_at=BUILT_AT, band_min=1.005, band_max=1.20)


class TestPricedRungs:
    def test_both_sides_are_stored(self):
        result = devig({"over": 1.91, "under": 1.91})
        assert result.reason is None
        assert {r.side for r in result.rows} == {"over", "under"}

    def test_fairs_sum_to_one(self):
        rows = devig({"over": 1.60, "under": 2.45}).rows
        assert sum(r.fair_prob for r in rows) == pytest.approx(1.0)

    def test_raw_prices_are_retained_for_debugging(self):
        rows = {r.side: r for r in devig({"over": 1.60, "under": 2.45}).rows}
        assert rows["over"].raw_decimal == 1.60
        assert rows["over"].raw_decimal_opp == 2.45
        assert rows["under"].raw_decimal == 2.45
        assert rows["over"].raw_overround == pytest.approx(1 / 1.6 + 1 / 2.45)

    def test_built_at_is_the_fetch_time_not_now(self):
        # Every row from one structure fetch shares that fetch's timestamp.
        # Stamping "now" per row would make a 5-minute-old payload read fresh
        # to #99's age gate.
        assert all(r.built_at == BUILT_AT
                   for r in devig({"over": 1.91, "under": 1.91}).rows)

    def test_rows_carry_the_kalshi_game_id_and_start(self):
        row = devig({"over": 1.91, "under": 1.91}).rows[0]
        assert row.game_id == "26AUG252138CLELAA"
        assert row.game_start_time == GAME.start_utc


class TestExclusions:
    def test_crossed_rung_is_excluded_not_devigged(self):
        assert devig({"over": 2.10, "under": 2.10}).reason == "crossed"

    def test_blown_vig_rung_is_excluded(self):
        assert devig({"over": 1.20, "under": 1.30}).reason == "overround"

    def test_one_sided_rung_is_excluded_never_haircut(self):
        # FanDuel one-sides ~51 deep-alt rungs per slate. Haircutting them
        # with the Route-B vig fallback would put an invented number on the
        # surface where declining costs quorum nothing.
        result = devig({"over": 1.91})
        assert result.reason == "one_sided" and result.rows == []

    def test_absent_rung_is_unresolved(self):
        assert devig({}).reason == "unresolved"

    def test_unusable_price_folds_into_unresolved(self):
        # A price the surface cannot use is, from its point of view, a line
        # the book did not usably post.
        assert devig({"over": 1.0, "under": 1.91}).reason == "unresolved"


def test_exclusion_counts_add_and_bump():
    a, b = ExclusionCounts(crossed=1), ExclusionCounts(crossed=2, one_sided=3)
    a.add(b)
    a.bump("unresolved")
    assert (a.crossed, a.one_sided, a.unresolved) == (3, 3, 1)


class TestStore:
    def make_row(self, book="betmgm", route="structure", side="over",
                 period="FG", line=8.5, fair=0.5, built_at=BUILT_AT):
        return SurfaceRow(book=book, game_id=GAME.game_id,
                          game_start_time=GAME.start_utc, period=period,
                          market_type="total", line=line, side=side,
                          fair_prob=fair, raw_decimal=1.91,
                          raw_decimal_opp=1.91, raw_overround=1.047,
                          route=route, built_at=built_at)

    def test_lookup_is_by_leg(self):
        s = LegSurface()
        s.publish("betmgm", "structure", [self.make_row()])
        leg = CanonicalLeg(GAME.game_id, "total", 8.5, "over")
        assert s.get("betmgm", leg).fair_prob == 0.5

    def test_publish_replaces_rather_than_merges(self):
        # A rung the book stopped posting has to disappear. Merging would
        # leave it resting at its last price with no way for the age gate to
        # notice, because its neighbours keep refreshing.
        s = LegSurface()
        s.publish("betmgm", "structure", [self.make_row(line=8.5),
                                          self.make_row(line=9.5)])
        s.publish("betmgm", "structure", [self.make_row(line=8.5)])
        gone = CanonicalLeg(GAME.game_id, "total", 9.5, "over")
        assert s.get("betmgm", gone) is None
        assert s.row_count() == 1

    def test_two_routes_at_one_book_do_not_erase_each_other(self):
        # FanDuel runs both: structure owns ml/spread/I1, singles owns FG/F5
        # totals. One shared slice would have each pass wipe the other.
        s = LegSurface()
        s.publish("fanduel", "structure", [self.make_row(book="fanduel",
                                                         period="I1", line=0.5)])
        s.publish("fanduel", "singles", [self.make_row(book="fanduel",
                                                       route="singles")])
        assert s.row_count() == 2
        i1 = CanonicalLeg(GAME.game_id, "total", 0.5, "over", "I1")
        fg = CanonicalLeg(GAME.game_id, "total", 8.5, "over")
        assert s.get("fanduel", i1) is not None
        assert s.get("fanduel", fg) is not None

    def test_a_book_counts_once_even_across_two_routes(self):
        # Two routes at one book are ONE opinion. Counting them separately
        # would let FanDuel satisfy MIN_AGREEING_BOOKS on its own.
        s = LegSurface()
        s.publish("fanduel", "structure", [self.make_row(book="fanduel")])
        s.publish("fanduel", "singles",
                  [self.make_row(book="fanduel", route="singles", fair=0.6)])
        leg = CanonicalLeg(GAME.game_id, "total", 8.5, "over")
        assert list(s.book_fairs(leg)) == ["fanduel"]

    def test_a_route_collision_prefers_the_fresher_price(self):
        s = LegSurface()
        s.publish("fanduel", "structure", [self.make_row(book="fanduel",
                                                         fair=0.4)])
        s.publish("fanduel", "singles",
                  [self.make_row(book="fanduel", route="singles", fair=0.6,
                                 built_at=BUILT_AT + timedelta(seconds=10))])
        leg = CanonicalLeg(GAME.game_id, "total", 8.5, "over")
        assert s.book_fairs(leg)["fanduel"].fair_prob == 0.6

    def test_moneyline_keys_on_a_null_line(self):
        s = LegSurface()
        row = self.make_row()
        ml = SurfaceRow(**{**row.__dict__, "market_type": "ml", "line": None,
                           "side": "home"})
        s.publish("novig", "structure", [ml])
        assert s.get("novig", CanonicalLeg(GAME.game_id, "ml", None,
                                           "home")) is not None

    def test_book_fairs_returns_rows_of_any_age(self):
        # The age gate lives in #99. A store that silently dropped stale rows
        # would make its decline counts unreadable.
        s = LegSurface()
        s.publish("betmgm", "structure",
                  [self.make_row(built_at=BUILT_AT - timedelta(hours=3))])
        leg = CanonicalLeg(GAME.game_id, "total", 8.5, "over")
        assert "betmgm" in s.book_fairs(leg)
