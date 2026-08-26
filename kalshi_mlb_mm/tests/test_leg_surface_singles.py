"""Singles-route ingest: game matching and wide-row -> rung translation (#96).

The headline test here is the doubleheader. #95's spike matched FanDuel rows
on team names alone, silently took tomorrow's PHI @ SEA, and produced fairs
off by 0.05-0.11 in probability — a WRONG NUMBER, not a decline. The surface
has the identical trap, so start-time matching and fail-closed ambiguity are
tested before anything else.
"""
from datetime import datetime, timedelta, timezone

import pytest

from kalshi_common.legset import CanonicalLeg
from kalshi_mlb_mm.leg_surface import singles
from kalshi_mlb_mm.leg_surface.slate import SurfaceGame

# The real pair from #95: two "Philadelphia Phillies @ Seattle Mariners" rows
# on FanDuel's slate, 18.5 hours apart.
TONIGHT = datetime(2026, 8, 26, 1, 41)
TOMORROW = datetime(2026, 8, 26, 20, 11)

GAME = SurfaceGame(game_id="26AUG252141PHISEA",
                   home_team="Seattle Mariners",
                   away_team="Philadelphia Phillies",
                   start_utc=TONIGHT, legs=())


def row(game_id, start, **kw):
    base = {"game_id": game_id, "game_start_time": start,
            "home_team": "Seattle Mariners",
            "away_team": "Philadelphia Phillies", "period": "FG",
            "fetch_time": datetime(2026, 8, 26, 1, 0, tzinfo=timezone.utc),
            "home_ml": None, "away_ml": None, "home_spread": None,
            "home_spread_price": None, "away_spread": None,
            "away_spread_price": None, "total": None, "over_price": None,
            "under_price": None}
    base.update(kw)
    return base


class TestGameMatching:
    def test_picks_the_game_at_the_matching_start_time(self):
        rows = [row("35978270", TONIGHT), row("35982005", TOMORROW)]
        matched, reason = singles.match_book_game(GAME, rows,
                                                  tolerance_min=30)
        assert reason is None
        assert {r["game_id"] for r in matched} == {"35978270"}

    def test_teams_alone_would_have_taken_the_wrong_game(self):
        # Guard against a future "simplification" back to team-only matching:
        # with the tonight row removed, team-only would happily return
        # tomorrow's. Start-time matching declines instead.
        rows = [row("35982005", TOMORROW)]
        assert singles.match_book_game(GAME, rows, tolerance_min=30) == (
            [], "game_unmatched")

    def test_two_candidates_inside_tolerance_fail_closed(self):
        # Not "pick the closest": if the book's slate is ambiguous, pricing
        # anything risks the wrong game, and a decline costs one book.
        rows = [row("A", TONIGHT), row("B", TONIGHT + timedelta(minutes=5))]
        assert singles.match_book_game(GAME, rows, tolerance_min=30) == (
            [], "game_ambiguous")

    def test_small_schedule_disagreement_still_matches(self):
        # Kalshi's suffix time and a book's posted start differ by a few
        # minutes routinely; the tolerance exists for exactly that.
        rows = [row("35978270", TONIGHT + timedelta(minutes=7))]
        matched, reason = singles.match_book_game(GAME, rows,
                                                  tolerance_min=30)
        assert reason is None and len(matched) == 1

    def test_iso_string_start_times_are_accepted(self):
        rows = [row("35978270", "2026-08-26T01:41:00Z")]
        assert singles.match_book_game(GAME, rows,
                                       tolerance_min=30)[1] is None

    def test_other_teams_never_match(self):
        rows = [row("X", TONIGHT, home_team="Boston Red Sox")]
        assert singles.match_book_game(GAME, rows,
                                       tolerance_min=30)[1] == "game_unmatched"

    def test_empty_slate_is_unmatched_not_a_crash(self):
        assert singles.match_book_game(GAME, [], tolerance_min=30) == (
            [], "game_unmatched")


class TestBookDecimals:
    def test_main_row_yields_ml_spread_and_total_rungs(self):
        rows = [row("g", TONIGHT, home_ml=-140, away_ml=120,
                    home_spread=-1.5, home_spread_price=110,
                    away_spread=1.5, away_spread_price=-130,
                    total=8.5, over_price=-105, under_price=-115)]
        out = singles.book_decimals(rows)
        assert set(out) == {("FG", "ml", None), ("FG", "spread", -1.5),
                            ("FG", "total", 8.5)}
        assert set(out[("FG", "spread", -1.5)]) == {"home", "away"}
        assert set(out[("FG", "total", 8.5)]) == {"over", "under"}

    def test_spread_line_is_home_perspective(self):
        # home_spread == CanonicalLeg.line by construction (negative = home
        # favoured). A sign flip here would attach every spread row to the
        # opposite line.
        out = singles.book_decimals([row("g", TONIGHT, home_spread=-1.5,
                                         home_spread_price=110,
                                         away_spread=1.5,
                                         away_spread_price=-130)])
        assert ("FG", "spread", -1.5) in out

    def test_american_prices_become_decimals(self):
        out = singles.book_decimals([row("g", TONIGHT, total=8.5,
                                         over_price=100, under_price=-110)])
        assert out[("FG", "total", 8.5)]["over"] == pytest.approx(2.0)

    def test_f5_rows_are_kept(self):
        out = singles.book_decimals([row("g", TONIGHT, period="F5", total=4.5,
                                         over_price=-110, under_price=-110)])
        assert ("F5", "total", 4.5) in out

    def test_f3_and_f7_periods_are_ignored(self):
        # The scrapers emit them; Kalshi lists no such combo legs, so they
        # would only be dead rows.
        out = singles.book_decimals([row("g", TONIGHT, period="F7", total=6.5,
                                         over_price=-110, under_price=-110)])
        assert out == {}

    def test_a_one_sided_row_yields_one_side(self):
        out = singles.book_decimals([row("g", TONIGHT, total=8.5,
                                         over_price=-110)])
        assert set(out[("FG", "total", 8.5)]) == {"over"}


class TestRouteOwnership:
    def make_game(self, legs):
        return SurfaceGame(game_id=GAME.game_id, home_team=GAME.home_team,
                           away_team=GAME.away_team, start_utc=TONIGHT,
                           legs=tuple(legs))

    def test_a_rung_this_route_does_not_own_is_skipped_silently(self):
        # FanDuel's spreads come from its structure route. Counting them
        # missing here would make FD's exclusion counts unreadable.
        game = self.make_game([
            CanonicalLeg(GAME.game_id, "spread", -1.5, "home"),
            CanonicalLeg(GAME.game_id, "spread", -1.5, "away")])
        rows, counts = singles.price_game(
            "fanduel", game, [], owned_markets={("total", "FG")},
            built_at=datetime.now(timezone.utc), band_min=1.005,
            band_max=1.20)
        assert rows == []
        assert counts.unresolved == 0

    def test_an_owned_rung_the_book_lacks_counts_as_unresolved(self):
        game = self.make_game([
            CanonicalLeg(GAME.game_id, "total", 8.5, "over"),
            CanonicalLeg(GAME.game_id, "total", 8.5, "under")])
        _rows, counts = singles.price_game(
            "draftkings", game, [], owned_markets={("total", "FG")},
            built_at=datetime.now(timezone.utc), band_min=1.005,
            band_max=1.20)
        assert counts.unresolved == 1

    def test_an_owned_rung_the_book_has_is_priced(self):
        game = self.make_game([
            CanonicalLeg(GAME.game_id, "total", 8.5, "over"),
            CanonicalLeg(GAME.game_id, "total", 8.5, "under")])
        book_rows = [row("g", TONIGHT, total=8.5, over_price=-110,
                         under_price=-110)]
        rows, counts = singles.price_game(
            "draftkings", game, book_rows, owned_markets={("total", "FG")},
            built_at=datetime.now(timezone.utc), band_min=1.005,
            band_max=1.20)
        assert len(rows) == 2 and counts.unresolved == 0
        assert sum(r.fair_prob for r in rows) == pytest.approx(1.0)

    def test_i1_is_never_served_by_this_route(self):
        # Neither singles scraper emits single-inning markets, so I1 must not
        # appear in any book's owned markets. This asserts the config, which
        # is where the mistake would be made.
        from kalshi_mlb_mm import config
        for markets in config.SURFACE_SINGLES_MARKETS.values():
            assert all(period != "I1" for _mt, period in markets)
