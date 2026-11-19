import json
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path

import pytest

from nfl_specials.dk_book import DK_SGP_EVENT_URL, DK_SGP_PRICE_URL, DK_WARMUP_URL
from sgp_compare.pricing import price_at_venue
from sgp_compare.venues import build_venues
from sgp_compare.venues.base import VenueTransportError
from sgp_compare.venues.draftkings import (DraftKingsVenue, american_of, parse_legs, parse_nfl_games,
                                           price_request_params)

FIXTURES = Path(__file__).parent / "fixtures"
BEFORE_KICKOFF = datetime(2026, 11, 19, tzinfo=timezone.utc)
EVENT_ID = "34118011"
HOU_ML, IND_ML = "0ML84695444_1", "0ML84695444_3"
OVER_43_5, UNDER_43_5 = "0OU84695444O4350_1", "0OU84695444U4350_3"
HOU_MINUS_12_5 = "0HC84695444N1250_1"
TOTAL_MARKET, SPREAD_MARKET, ATD_MARKET = "3_84695444", "2_84695444", "389630297"


def load(name):
    return json.loads((FIXTURES / name).read_text())


PRICE_CASES = load("draftkings_price_responses.json")


@pytest.fixture
def game():
    return parse_nfl_games(load("draftkings_league.json"), BEFORE_KICKOFF)[0]


@pytest.fixture
def legs():
    return {(leg.market_label, leg.side_label): leg for leg in parse_legs(load("draftkings_sgp_event_ind_hou.json"))}


# ---------------------------------------------------------------- games

def test_games_use_full_team_names_soonest_first(game):
    games = parse_nfl_games(load("draftkings_league.json"), BEFORE_KICKOFF)
    assert (game.venue_game_id, game.home_team, game.away_team) == (EVENT_ID, "Houston Texans", "Indianapolis Colts")
    assert game.game_start_time == datetime(2026, 11, 20, 1, 15, tzinfo=timezone.utc)
    assert [g.game_start_time for g in games] == sorted(g.game_start_time for g in games)


def test_started_games_are_left_out():
    assert parse_nfl_games(load("draftkings_league.json"), datetime(2026, 12, 31, tzinfo=timezone.utc)) == []


# ---------------------------------------------------------------- legs

def test_listed_price_reads_dk_unicode_minus():
    assert american_of({"displayOdds": {"american": "−118"}}) == -118.0
    assert american_of({"displayOdds": {"american": "+500"}}) == 500.0
    assert american_of({"displayOdds": {}}) is None


def test_spread_pairs_the_other_team_at_the_opposite_number(legs):
    hou = legs[("Spread", "Houston Texans -12.5")]
    ind = legs[("Spread", "Indianapolis Colts +12.5")]
    assert (hou.canonical.period, hou.canonical.stat, hou.canonical.side, hou.canonical.line) == ("FG", "spread", "home", -12.5)
    assert (hou.american_odds, hou.opposite_american_odds, hou.opposite_leg_id) == (-108.0, -112.0, ind.leg_id)
    assert hou.is_main_line and hou.category == "Game Lines" and hou.market_key == SPREAD_MARKET
    assert hou.venue_ref == {"selectionId": HOU_MINUS_12_5, "marketId": SPREAD_MARKET}


def test_alternate_spread_is_not_main_and_pairs_within_its_market(legs):
    hou = legs[("Spread Alternate", "Houston Texans -9.5")]
    assert not hou.is_main_line
    assert hou.opposite_leg_id == legs[("Spread Alternate", "Indianapolis Colts +9.5")].leg_id


def test_moneyline_total_and_period_markets(legs):
    assert legs[("Moneyline", "Houston Texans")].canonical.side == "home"
    assert legs[("Moneyline", "Houston Texans")].opposite_american_odds == 500.0
    assert legs[("Moneyline 1st Half", "Houston Texans")].canonical.period == "1H"
    over = legs[("Total", "Over 43.5")]
    assert (over.canonical.stat, over.canonical.side, over.canonical.line, over.opposite_leg_id) == ("total", "over", 43.5, UNDER_43_5)
    team_total = legs[("Houston Texans: Team Total Points - 1st Half", "Over 13.5")]
    assert (team_total.canonical.period, team_total.canonical.stat, team_total.canonical.subject) == ("1H", "team_total", "Houston Texans")


def test_three_way_legs_have_no_single_other_side(legs):
    tie = legs[("1st Half (3 Way)", "Tie")]
    hou = legs[("1st Half (3 Way)", "Houston Texans")]
    assert (tie.canonical.stat, tie.canonical.subject, tie.opposite_leg_id) == ("win_3way", "Draw", None)
    assert hou.canonical.subject == "Houston Texans"


def test_player_props_name_player_and_stat(legs):
    over = legs[("Nico Collins Receiving Yards O/U", "Over 70.5")]
    assert (over.canonical.subject, over.canonical.stat, over.canonical.line) == ("Nico Collins", "receiving_yards", 70.5)
    assert over.opposite_american_odds == -111.0 and over.is_main_line
    milestone = legs[("Nico Collins Receiving Yards", "40+")]
    assert (milestone.canonical.stat, milestone.canonical.side, milestone.canonical.line) == ("receiving_yards", "over", 39.5)
    assert milestone.opposite_leg_id is None and not milestone.is_main_line


def test_anytime_td_is_one_sided_and_named_like_prophetx(legs):
    td = legs[("Anytime TD Scorer", "Nico Collins")]
    assert (td.canonical.stat, td.canonical.subject, td.canonical.side) == ("anytime_td", "Nico Collins", "yes")
    assert td.opposite_leg_id is None and td.category == "TD Scorers"


def test_unknown_markets_are_offered_without_canonical(legs):
    correct_score = legs[("Correct Score", "HOU Texans 23-20 IND Colts")]
    assert correct_score.canonical is None and correct_score.category == "Specials"
    yes = legs[("Both Teams to Score - 1st Quarter", "Yes")]
    assert (yes.canonical.period, yes.canonical.stat, yes.canonical.side) == ("1Q", "both_teams_to_score", "yes")


# ---------------------------------------------------------------- pricing

class FakeResponse:
    def __init__(self, status_code, body):
        self.status_code, self._body = status_code, body
        self.text = json.dumps(body)

    def json(self):
        return self._body


class FakeSession:
    """Answers DK's URLs from the recon fixtures; records every GET."""

    def __init__(self, price_answers=None, fail=None):
        self.calls = []
        self.price_answers = price_answers or {}
        self.fail = fail

    def get(self, url, timeout=None, params=None, headers=None):
        self.calls.append((url, params, headers))
        if self.fail:
            raise self.fail
        if url == DK_WARMUP_URL:
            return FakeResponse(200, {})
        if url == DK_SGP_EVENT_URL.format(event_id=EVENT_ID):
            return FakeResponse(200, load("draftkings_sgp_event_ind_hou.json"))
        if url == DK_SGP_PRICE_URL:
            key = (params["selections"], params["marketCandidates"])
            if key in self.price_answers:
                return FakeResponse(200, self.price_answers[key])
            return FakeResponse(500, {"error": f"no fixture for {key}"})
        return FakeResponse(404, {})

    def price_calls(self):
        return [params for url, params, _ in self.calls if url == DK_SGP_PRICE_URL]


def recon_answers():
    return {(case["params"]["selections"], case["params"]["marketCandidates"]): case["response"] for case in PRICE_CASES}


def venue_with(session):
    clock = iter(range(0, 10_000, 1)).__next__     # 1 s per reading: the throttle never waits
    return DraftKingsVenue(session_factory=lambda: session, clock=clock, sleep=lambda s: None)


def test_price_request_sorts_base_selections():
    params = price_request_params(EVENT_ID, frozenset({OVER_43_5, HOU_ML}), ATD_MARKET)
    assert params == PRICE_CASES[2]["params"]


def test_two_leg_parlay_prices_off_the_builder(game, legs):
    session = FakeSession(recon_answers())
    quote = venue_with(session).price_parlay(game, [legs[("Moneyline", "Houston Texans")], legs[("Total", "Over 43.5")]])
    assert quote.best.decimal_odds == 2.0 and quote.best.max_stake is None and quote.leg_prices_american == {}
    (params,) = session.price_calls()
    assert params == PRICE_CASES[0]["params"]
    assert session.calls[0][0] == DK_WARMUP_URL
    assert session.calls[1][2] == {"X-SportId": "3"}


def test_flipping_the_last_leg_reuses_the_same_builder_answer(game, legs):
    session = FakeSession(recon_answers())
    venue = venue_with(session)
    hou = legs[("Moneyline", "Houston Texans")]
    venue.price_parlay(game, [hou, legs[("Total", "Over 43.5")]])
    under = venue.price_parlay(game, [hou, legs[("Total", "Under 43.5")]])
    assert under.best.decimal_odds == 2.14
    assert len(session.price_calls()) == 1


def test_three_legs_with_a_player_prop(game, legs):
    session = FakeSession(recon_answers())
    quote = venue_with(session).price_parlay(game, [legs[("Moneyline", "Houston Texans")], legs[("Total", "Over 43.5")],
                                                    legs[("Anytime TD Scorer", "Nico Collins")]])
    expected = next(s["trueOdds"] for m in PRICE_CASES[2]["response"]["data"]["compatibleMarkets"]
                    for s in m["selections"] if s["id"] == legs[("Anytime TD Scorer", "Nico Collins")].leg_id)
    assert quote.best.decimal_odds == expected


def test_a_base_dk_drops_a_leg_from_prices_nothing(game, legs):
    """DK quietly prices IND ML alone when asked for HOU ML + IND ML."""
    session = FakeSession(recon_answers())
    venue = venue_with(session)
    hou, ind = legs[("Moneyline", "Houston Texans")], legs[("Moneyline", "Indianapolis Colts")]
    # The pricer refuses two legs of one market, so give IND ML a stand-in market_key.
    ind_other_market = replace(ind, market_key="elsewhere")
    assert venue.price_parlay(game, [hou, ind_other_market, legs[("Total", "Over 43.5")]]) is None


def test_a_candidate_market_dk_will_not_combine_is_a_decline(game, legs):
    """The recon asked for HOU -12.5 + the Spread market itself: compatibleMarkets came back empty."""
    session = FakeSession(recon_answers())
    ind = legs[("Spread", "Indianapolis Colts +12.5")]
    assert venue_with(session).price_parlay(game, [legs[("Spread", "Houston Texans -12.5")], ind]) is None


def test_transport_failures_and_bad_status_raise(game, legs):
    with pytest.raises(VenueTransportError, match="warmup"):
        venue_with(FakeSession(fail=OSError("timed out"))).list_legs(game)
    with pytest.raises(VenueTransportError, match="HTTP 200, got 500"):
        venue_with(FakeSession({})).price_parlay(game, [legs[("Moneyline", "Houston Texans")], legs[("Total", "Over 43.5")]])


def test_price_calls_wait_out_dk_widget_pace(game, legs):
    session = FakeSession(recon_answers())
    waits = []
    venue = DraftKingsVenue(session_factory=lambda: session, clock=lambda: 100.0, sleep=waits.append)
    hou = legs[("Moneyline", "Houston Texans")]
    venue.price_parlay(game, [hou, legs[("Total", "Over 43.5")]])
    venue.price_parlay(game, [legs[("Spread", "Houston Texans -12.5")], legs[("Nico Collins Receiving Yards O/U", "Over 70.5")]])
    assert waits == [0.5]


def test_full_partition_through_the_pricer_halves_the_calls(game, legs):
    """2 legs = 4 combinations; each builder answer covers both sides of the last leg."""
    answers = recon_answers()
    answers[(IND_ML, TOTAL_MARKET)] = {"data": {"selectionsMapped": [IND_ML], "selectionsNotMapped": [],
                                                "compatibleMarkets": [{"id": TOTAL_MARKET, "selections": [
                                                    {"id": OVER_43_5, "trueOdds": 10.5}, {"id": UNDER_43_5, "trueOdds": 9.5}]}]}}
    session = FakeSession(answers)
    result = price_at_venue(venue_with(session), game, [HOU_ML, OVER_43_5])
    assert result.status == "priced"
    assert result.partition_decimal_odds == [2.0, 2.14, 10.5, 9.5]
    assert len(session.price_calls()) == 2


def test_build_venues_includes_draftkings():
    assert build_venues()["draftkings"].display_name == "DraftKings"
