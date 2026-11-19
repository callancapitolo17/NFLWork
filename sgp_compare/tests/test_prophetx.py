import json
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest

from sgp_compare.legs import CanonicalLeg
from sgp_compare.venues.base import VenueTransportError
from sgp_compare.venues.prophetx import (ProphetXVenue, parlay_payload, parse_legs, parse_nfl_games,
                                         parse_parlay_quote)

FIXTURES = Path(__file__).parent / "fixtures"
HOUSTON_ID = 1700004938
BEFORE_KICKOFF = datetime(2026, 11, 19, tzinfo=timezone.utc)
KICKOFF = datetime(2026, 11, 20, 1, 15, tzinfo=timezone.utc)


def load(name):
    return json.loads((FIXTURES / name).read_text())


@pytest.fixture
def game():
    return parse_nfl_games(load("prophetx_tournaments.json"), BEFORE_KICKOFF)[0]


@pytest.fixture
def legs(game):
    return {(leg.market_label, leg.side_label): leg
            for leg in parse_legs(game, load("prophetx_markets_ind_hou.json"), HOUSTON_ID)}


# ---------------------------------------------------------------- games

def test_games_are_the_nfl_tournament_only_soonest_first_home_is_seq_zero():
    games = parse_nfl_games(load("prophetx_tournaments.json"), BEFORE_KICKOFF)
    assert [(g.venue_game_id, g.home_team, g.away_team, g.game_start_time) for g in games] == [
        ("1700008817", "Houston Texans", "Indianapolis Colts", KICKOFF),
        ("1700008828", "Carolina Panthers", "Baltimore Ravens", datetime(2026, 11, 22, 18, 0, tzinfo=timezone.utc)),
    ]


def test_started_games_and_an_offseason_payload_give_no_games():
    assert parse_nfl_games(load("prophetx_tournaments.json"), datetime(2026, 12, 31, tzinfo=timezone.utc)) == []
    assert parse_nfl_games({"data": {"tournaments": [{"name": "NBA", "sportEvents": []}]}}, BEFORE_KICKOFF) == []
    assert parse_nfl_games({"data": None}, BEFORE_KICKOFF) == []


def test_a_start_time_without_a_zone_is_read_as_utc_not_local():
    payload = load("prophetx_tournaments.json")
    payload["data"]["tournaments"][0]["sportEvents"][0]["scheduled"] = "2026-11-20T01:15:00"
    assert parse_nfl_games(payload, BEFORE_KICKOFF)[0].game_start_time == KICKOFF


# ---------------------------------------------------------------- legs

def test_spread_leg_prices_sides_and_ids(legs):
    leg = legs[("Spread", "Houston Texans -12.5")]
    assert (leg.american_odds, leg.opposite_american_odds, leg.available_stake) == (100.0, -105.0, 838.54)
    assert leg.canonical == CanonicalLeg("FG", "spread", None, "home", -12.5)
    assert leg.opposite_leg_id == legs[("Spread", "Indianapolis Colts +12.5")].leg_id
    assert legs[("Spread", "Indianapolis Colts +12.5")].canonical == CanonicalLeg("FG", "spread", None, "away", 12.5)
    assert leg.venue_ref == {"sportEventId": 1700008817, "marketId": 223, "outcomeId": 1714,
                             "lineId": "7342bc136571bb4c9584fd65c027c042", "line": -12.5}
    assert (leg.market_key, leg.category, leg.is_main_line) == ("223", "Game Lines", True)


def test_moneyline_and_totals(legs):
    assert legs[("Moneyline", "Houston Texans")].canonical == CanonicalLeg("FG", "moneyline", None, "home", None)
    assert legs[("Moneyline", "Indianapolis Colts")].american_odds == 550.0
    assert legs[("Total Points", "Over 43.5")].canonical == CanonicalLeg("FG", "total", None, "over", 43.5)


def test_moneyline_side_falls_back_to_the_competitor_id(game):
    markets = load("prophetx_markets_ind_hou.json")
    moneyline = next(m for m in markets["data"]["markets"] if m["name"] == "Moneyline")
    for outcome in moneyline["outcomes"]:
        outcome["id"] = 9999                       # an outcome id the side table does not know
    for ladder in moneyline["selections"]:
        for entry in ladder:
            entry["id"] = 9999
    parsed = {leg.side_label: leg for leg in parse_legs(game, markets, HOUSTON_ID) if leg.market_label == "Moneyline"}
    assert parsed["Houston Texans"].canonical.side == "home"
    assert parsed["Indianapolis Colts"].canonical.side == "away"


def test_player_props_name_the_player_and_stat(legs):
    assert legs[("Nico Collins Receiving Yards", "Over 70.5")].canonical == CanonicalLeg(
        "FG", "receiving_yards", "Nico Collins", "over", 70.5)
    assert legs[("Nico Collins Receiving Yards", "Over 70.5")].is_main_line
    assert legs[("David Montgomery To Score a Touchdown", "YES")].canonical == CanonicalLeg(
        "FG", "anytime_td", "David Montgomery", "yes", None)


def test_periods_team_totals_and_quarter_three_way(legs):
    assert legs[("First Half Spread", "Houston Texans -6.5")].canonical == CanonicalLeg("1H", "spread", None, "home", -6.5)
    assert legs[("Indianapolis Colts: Team Total Points", "Over 15.5")].canonical == CanonicalLeg(
        "FG", "team_total", "Indianapolis Colts", "over", 15.5)
    assert legs[("Draw (2nd Quarter)", "YES")].canonical == CanonicalLeg("2Q", "win_3way", "Draw", "yes", None)


def test_lines_nobody_offers_are_left_out_but_a_one_sided_line_keeps_both_sides(game):
    markets = load("prophetx_markets_ind_hou.json")
    collins = next(m for m in markets["data"]["markets"] if m["name"] == "Nico Collins Receiving Yards")
    main_line = next(line for line in collins["marketLines"] if line["name"] == "Fixed total 70.5")
    main_line["selections"][1] = []                # nobody offers the under
    parsed = {(l.market_label, l.side_label): l for l in parse_legs(game, markets, HOUSTON_ID)}
    assert ("Nico Collins Receiving Yards", "Over 14.5") not in parsed
    over, under = parsed[("Nico Collins Receiving Yards", "Over 70.5")], parsed[("Nico Collins Receiving Yards", "Under 70.5")]
    assert (over.american_odds, over.opposite_american_odds) == (-112.0, None)
    assert (under.american_odds, under.opposite_american_odds) == (None, -112.0)


def test_inactive_markets_are_left_out(game):
    markets = load("prophetx_markets_ind_hou.json")
    for market in markets["data"]["markets"]:
        market["status"] = "suspended"
    assert parse_legs(game, markets, HOUSTON_ID) == []


# ---------------------------------------------------------------- parlay quotes

@pytest.mark.parametrize("fixture, sides, tiers, inside", [
    ("prophetx_parlay_spread_over.json", [("Spread", "Houston Texans -12.5"), ("Total Points", "Over 43.5")],
     [(2.7165, 62.62), (2.6758, 125.24)], [-148.0, -161.0]),
    ("prophetx_parlay_spread_recyds.json", [("Spread", "Houston Texans -12.5"), ("Nico Collins Receiving Yards", "Over 70.5")],
     [(3.173, 1500)], [-129.0, -127.0]),
])
def test_recorded_parlay_requests_and_answers(legs, fixture, sides, tiers, inside):
    recorded = load(fixture)
    chosen = [legs[side] for side in sides]
    assert parlay_payload(chosen) == recorded["request"]
    quote = parse_parlay_quote(recorded["response"], chosen)
    assert [(round(t.decimal_odds, 4), t.max_stake) for t in quote.tiers] == tiers
    assert quote.leg_prices_american == {leg.leg_id: price for leg, price in zip(chosen, inside)}
    assert quote.raw == recorded["response"]


def test_declined_parlays_are_none_and_unusable_offers_are_skipped(legs):
    chosen = [legs[("Spread", "Houston Texans -12.5")], legs[("Total Points", "Over 43.5")]]
    assert parse_parlay_quote({"success": False, "data": None}, chosen) is None
    assert parse_parlay_quote({"success": True, "data": {"offers": []}}, chosen) is None
    junk_first = {"success": True, "data": {"offers": [
        {"odds": 50, "stake": 10, "estimatedPrices": [{"lineId": chosen[0].leg_id, "odds": 900}]},
        {"odds": "x", "stake": 10},
        {"odds": 150, "stake": 20, "estimatedPrices": [{"lineId": chosen[0].leg_id, "odds": -120},
                                                       {"lineId": chosen[1].leg_id, "odds": None}]}]}}
    quote = parse_parlay_quote(junk_first, chosen)
    assert [(t.decimal_odds, t.max_stake) for t in quote.tiers] == [(2.5, 20)]
    assert quote.leg_prices_american == {chosen[0].leg_id: -120.0}


# ---------------------------------------------------------------- the live venue, on a fake session

class FakeResponse:
    def __init__(self, status_code, body):
        self.status_code = status_code
        self._body = body
        self.text = body if isinstance(body, str) else json.dumps(body)

    def json(self):
        if isinstance(self._body, str):
            raise ValueError("not JSON")
        return self._body


class FakeSession:
    def __init__(self, answers):
        self.answers = answers          # path fragment -> FakeResponse or Exception
        self.calls = []

    def request(self, method, url, timeout=None, **kwargs):
        self.calls.append((method, url))
        for fragment, answer in self.answers.items():
            if fragment in url:
                if isinstance(answer, Exception):
                    raise answer
                return answer
        raise AssertionError(f"unexpected {method} {url}")


def live_answers(**overrides):
    answers = {"tournaments": FakeResponse(200, load("prophetx_tournaments.json")),
               "/markets": FakeResponse(200, load("prophetx_markets_ind_hou.json")),
               "/parlay/": FakeResponse(200, load("prophetx_parlay_spread_over.json")["response"])}
    answers.update(overrides)
    return answers


def venue_with(answers, wall_time=BEFORE_KICKOFF):
    session = FakeSession(answers)
    clock = {"t": 0.0}
    venue = ProphetXVenue(session_factory=lambda: session, clock=lambda: clock["t"],
                          wall_clock=lambda: wall_time)
    return venue, session, clock


def test_venue_reads_games_legs_and_prices():
    venue, session, _ = venue_with(live_answers())
    game = venue.list_games()[0]
    legs = {(l.market_label, l.side_label): l for l in venue.list_legs(game)}
    quote = venue.price_parlay(game, [legs[("Spread", "Houston Texans -12.5")], legs[("Total Points", "Over 43.5")]])
    assert round(quote.best.decimal_odds, 4) == 2.7165
    assert [method for method, _ in session.calls] == ["GET", "GET", "POST"]


def test_schedule_is_cached_but_a_started_game_drops_out_at_once():
    answers = live_answers()
    session = FakeSession(answers)
    now = {"wall": BEFORE_KICKOFF}
    venue = ProphetXVenue(session_factory=lambda: session, clock=lambda: 0.0, wall_clock=lambda: now["wall"])
    assert len(venue.list_games()) == 2
    now["wall"] = KICKOFF + timedelta(seconds=1)
    assert [g.venue_game_id for g in venue.list_games()] == ["1700008828"]
    assert sum("tournaments" in url for _, url in session.calls) == 1


def test_market_cache_and_fresh_reads():
    venue, session, clock = venue_with(live_answers())
    game = venue.list_games()[0]
    venue.list_legs(game)
    venue.list_legs(game)
    venue.list_legs(game, fresh=True)
    clock["t"] = 1000.0
    venue.list_legs(game)
    assert sum(url.endswith("/markets") for _, url in session.calls) == 3


@pytest.mark.parametrize("answer, message", [
    (FakeResponse(403, "blocked"), "expected HTTP 200, got 403: blocked"),
    (FakeResponse(200, "<html>"), "response is not JSON"),
    (FakeResponse(200, []), "expected a JSON object, got list"),
    (RuntimeError("connection reset"), "RuntimeError: connection reset"),
])
def test_transport_failures_raise_venue_transport_error(answer, message):
    venue, _, _ = venue_with(live_answers(tournaments=answer))
    with pytest.raises(VenueTransportError, match=message):
        venue.list_games()


def test_a_failed_request_rebuilds_the_session():
    built = []

    def factory():
        built.append(1)
        return FakeSession(live_answers(tournaments=RuntimeError("reset")))
    venue = ProphetXVenue(session_factory=factory)
    for _ in range(2):
        with pytest.raises(VenueTransportError):
            venue.list_games()
    assert len(built) == 2


def test_an_answer_with_an_unexpected_shape_is_a_transport_error_not_a_crash():
    weird_tournaments = {"data": {"tournaments": [{"name": "NFL", "sportEvents": [
        {"status": "not_started", "scheduled": "2026-11-20T01:15:00Z", "competitors": [{"seq": 0}, {"seq": 1}]}]}]}}
    venue, _, _ = venue_with(live_answers(tournaments=FakeResponse(200, weird_tournaments)))
    with pytest.raises(VenueTransportError, match="tournaments answer has an unexpected shape"):
        venue.list_games()

    venue, _, _ = venue_with(live_answers(**{"/parlay/": FakeResponse(200, {"success": True, "data": {"offers": [
        {"odds": 150, "estimatedPrices": [{"lineId": None}, "x"]}]}})}))
    game = venue.list_games()[0]
    legs = venue.list_legs(game)
    with pytest.raises(VenueTransportError, match="parlay answer has an unexpected shape"):
        venue.price_parlay(game, legs[:2])

