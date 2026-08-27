"""Slate + leg-ladder discovery for the leg surface (issue #96).

The invariants here are all about IDENTITY. The surface exists so a
cross-game combo prices without a network call, which means a leg looked up
by the router must be the same leg the ingest priced — same game, same sign,
same period. Every test below is a way that correspondence breaks silently.
"""
from datetime import datetime, timedelta

import pytest

from kalshi_common.legset import CanonicalLeg
from kalshi_mlb_mm.leg_surface import slate

SUFFIX = "26AUG251840BOSMIA"        # BOS @ MIA — away BOS, home MIA


class TestLegsForMarket:
    def test_spread_market_yields_both_sides_of_one_signed_line(self):
        legs = slate.legs_for_market(
            "KXMLBSPREAD", SUFFIX, f"KXMLBSPREAD-{SUFFIX}-MIA4")
        # MIA is HOME, so its -3.5 margin market is home-perspective -3.5.
        assert {(l.line, l.side) for l in legs} == {(-3.5, "home"),
                                                    (-3.5, "away")}
        assert all(l.period == "FG" and l.market_type == "spread"
                   for l in legs)

    def test_away_spread_market_carries_the_opposite_sign(self):
        legs = slate.legs_for_market(
            "KXMLBSPREAD", SUFFIX, f"KXMLBSPREAD-{SUFFIX}-BOS4")
        # #70: the sign follows the ticker's TEAM. Collapsing both teams onto
        # one sign selected the wrong grid cell; here it would point the
        # router at a different line entirely.
        assert {l.line for l in legs} == {3.5}

    def test_moneyline_market_has_no_line(self):
        legs = slate.legs_for_market(
            "KXMLBGAME", SUFFIX, f"KXMLBGAME-{SUFFIX}-MIA")
        assert {(l.market_type, l.line) for l in legs} == {("ml", None)}

    def test_rfi_re_encodes_as_a_first_inning_total_at_half(self):
        # #87: KXMLBRFI YES is exactly "1st-inning runs >= 1". The surface
        # must carry it as an I1 total at 0.5 or the books' YRFI markets have
        # nothing to attach to.
        legs = slate.legs_for_market("KXMLBRFI", SUFFIX, f"KXMLBRFI-{SUFFIX}")
        assert {(l.period, l.market_type, l.line, l.side) for l in legs} == {
            ("I1", "total", 0.5, "over"), ("I1", "total", 0.5, "under")}

    def test_f5_winner_re_encodes_as_a_half_run_line(self):
        # #86: books price the F5 run line two-sided with no push; the F5
        # MONEYLINE is conditional push-2-way and would overprice YES by
        # ~P(tie).
        legs = slate.legs_for_market("KXMLBF5", SUFFIX, f"KXMLBF5-{SUFFIX}-BOS")
        assert all(l.market_type == "spread" and l.period == "F5"
                   for l in legs)
        assert {l.line for l in legs} == {0.5}     # BOS is away -> home +0.5

    def test_f5_tie_market_is_dropped(self):
        # It types as (ml, F5) and no book leg can express it. Carrying it
        # would only inflate the unresolved count every pass, forever.
        assert slate.legs_for_market(
            "KXMLBF5", SUFFIX, f"KXMLBF5-{SUFFIX}-TIE") == []

    def test_unparseable_ticker_drops_rather_than_raises(self):
        assert slate.legs_for_market("KXMLBTOTAL", SUFFIX,
                                     f"KXMLBTOTAL-{SUFFIX}-notanumber") == []


class TestWindow:
    NOW = datetime(2026, 8, 26, 12, 0, 0)

    def test_game_in_window_is_kept(self):
        assert slate.in_window(self.NOW + timedelta(hours=6), self.NOW,
                               min_minutes=5, max_hours=12)

    def test_game_inside_the_tipoff_guard_is_dropped(self):
        # The maker cancels resting quotes at TIPOFF_CANCEL_MIN; fetching
        # prices it will never quote on is pure book traffic.
        assert not slate.in_window(self.NOW + timedelta(minutes=2), self.NOW,
                                   min_minutes=5, max_hours=12)

    def test_started_game_is_dropped(self):
        assert not slate.in_window(self.NOW - timedelta(minutes=1), self.NOW,
                                   min_minutes=5, max_hours=12)

    def test_far_future_game_is_dropped(self):
        # 48 KXMLBGAME events are open at once (~3 days out). #95 measured
        # books posting main lines only that far ahead, so pulling them is
        # cost without coverage.
        assert not slate.in_window(self.NOW + timedelta(hours=30), self.NOW,
                                   min_minutes=5, max_hours=12)


def test_rungs_group_both_sides_of_a_line_together():
    over = CanonicalLeg("G", "total", 8.5, "over")
    under = CanonicalLeg("G", "total", 8.5, "under")
    other = CanonicalLeg("G", "total", 9.5, "over")
    grouped = slate.rungs([over, under, other])
    assert grouped[("FG", "total", 8.5)] == [over, under]
    assert grouped[("FG", "total", 9.5)] == [other]


def test_rungs_keep_periods_apart():
    # An FG 8.5 total and an F5 8.5 total are different markets that agree on
    # every other field. Devigging them as one rung would cross the periods.
    fg = CanonicalLeg("G", "total", 8.5, "over", "FG")
    f5 = CanonicalLeg("G", "total", 8.5, "over", "F5")
    assert len(slate.rungs([fg, f5])) == 2


def test_dedupe_keeps_one_row_per_distinct_leg():
    leg = CanonicalLeg("G", "spread", -1.5, "home")
    assert slate._dedupe_legs([leg, leg]) == (leg,)


class TestEventPagination:
    """Review finding 4: the events listing was a single `limit=50` call with
    no pagination.

    42-48 KXMLBGAME events are open at once, so it worked — with eight games
    of headroom. The failure mode past that is SILENT: the overflow is
    dropped, and the 12h window filter runs afterwards, so a busy Saturday
    with doubleheaders would quietly lose games we should be quoting.
    """

    def _api(self, pages):
        calls = []

        def api(_method, query):
            calls.append(query)
            return 200, pages[len(calls) - 1], None

        return api, calls

    def test_follows_the_cursor_across_pages(self, monkeypatch):
        pages = [
            {"events": [{"event_ticker": "KXMLBGAME-A"}], "cursor": "c1"},
            {"events": [{"event_ticker": "KXMLBGAME-B"}], "cursor": ""},
        ]
        api, calls = self._api(pages)
        monkeypatch.setattr(slate.auth_client, "api", api)
        assert slate._fetch_open_game_events() == ["KXMLBGAME-A",
                                                   "KXMLBGAME-B"]
        assert "cursor=c1" in calls[1]

    def test_stops_when_the_cursor_empties(self, monkeypatch):
        api, calls = self._api([{"events": [], "cursor": ""}])
        monkeypatch.setattr(slate.auth_client, "api", api)
        slate._fetch_open_game_events()
        assert len(calls) == 1

    def test_a_failed_page_returns_nothing_rather_than_a_short_slate(
            self, monkeypatch):
        # Fail CLOSED: refresh_slate keeps the previous slate on an empty
        # result, but would happily adopt a truncated one and drop every game
        # the failed page carried.
        def api(_method, query):
            if "cursor=" in query:
                return 500, None, None
            return 200, {"events": [{"event_ticker": "KXMLBGAME-A"}],
                         "cursor": "c1"}, None

        monkeypatch.setattr(slate.auth_client, "api", api)
        assert slate._fetch_open_game_events() == []

    def test_non_game_series_tickers_are_ignored(self, monkeypatch):
        api, _calls = self._api([{"events": [{"event_ticker": "KXMLBGAME-A"},
                                             {"event_ticker": "KXNFLGAME-B"}],
                                  "cursor": ""}])
        monkeypatch.setattr(slate.auth_client, "api", api)
        assert slate._fetch_open_game_events() == ["KXMLBGAME-A"]
