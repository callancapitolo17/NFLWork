import json
from datetime import datetime, timezone
from pathlib import Path

import pytest

from sgp_compare.legs import CanonicalLeg
from sgp_compare.venues.base import VenueTransportError
from sgp_compare.venues.caesars import (SCHEDULE_TAB_PATHS, CaesarsVenue, bet_details_body, clean_name,
                                        parse_legs, parse_nfl_games, parse_parlay_quote)

FIXTURES = Path(__file__).parent / "fixtures"
BEFORE_KICKOFF = datetime(2026, 11, 19, tzinfo=timezone.utc)
EVENT_ID = "faa98a96-ed6b-e0a9-94c7-71f707d9e752"


def load(name):
    return json.loads((FIXTURES / name).read_text())


@pytest.fixture
def game():
    return parse_nfl_games([load("caesars_tabs_featured.json")], BEFORE_KICKOFF)[0]


@pytest.fixture
def legs(game):
    return {(leg.market_label, leg.side_label): leg
            for leg in parse_legs(game, load("caesars_event_ind_hou.json")["event"])}


# ---------------------------------------------------------------- games

def test_games_are_nfl_not_started_soonest_first(game):
    games = parse_nfl_games([load("caesars_tabs_featured.json")], BEFORE_KICKOFF)
    assert [g.venue_game_id[:8] for g in games] == ["faa98a96", "c8c55c7a", "b6cc2370"]   # no NCAAF, no started
    assert (game.home_team, game.away_team) == ("Houston Texans", "Indianapolis Colts")
    assert game.game_start_time == datetime(2026, 11, 20, 1, 15, tzinfo=timezone.utc)


def test_games_from_two_tabs_are_listed_once():
    tabs = load("caesars_tabs_featured.json")
    assert len(parse_nfl_games([tabs, tabs], BEFORE_KICKOFF)) == 3


def test_games_after_kickoff_are_left_out():
    assert parse_nfl_games([load("caesars_tabs_featured.json")], datetime(2026, 12, 31, tzinfo=timezone.utc)) == []


# ---------------------------------------------------------------- legs

def test_names_lose_their_pipes():
    assert clean_name("|Houston Texans| |First Score Of Game Will Be|") == "Houston Texans First Score Of Game Will Be"


def test_main_spread_is_home_line_with_negated_away(legs):
    home, away = legs[("Spread", "Houston Texans -12.5")], legs[("Spread", "Indianapolis Colts +12.5")]
    assert home.canonical == CanonicalLeg("FG", "spread", None, "home", -12.5)
    assert away.canonical == CanonicalLeg("FG", "spread", None, "away", 12.5)
    assert (home.american_odds, home.opposite_american_odds) == (-106.0, -114.0)
    assert home.opposite_leg_id == away.leg_id and away.opposite_leg_id == home.leg_id
    assert (home.category, home.is_main_line, home.available_stake) == ("Game Lines", True, 50000.0)
    assert (home.venue_ref["line"], away.venue_ref["line"]) == (-12.5, 12.5)
    assert home.venue_ref["selectionId"] == "6e1b9a51-fbdc-3f7e-be62-380aa2ad1d74"


def test_alt_rungs_share_the_main_selection_and_market(legs):
    main, alt = legs[("Spread", "Houston Texans -12.5")], legs[("Spread", "Houston Texans -7.5")]
    assert not alt.is_main_line and alt.market_key == main.market_key
    assert alt.venue_ref["selectionId"] == main.venue_ref["selectionId"] and alt.leg_id != main.leg_id
    assert alt.venue_ref["line"] == -7.5 and alt.canonical.line == -7.5
    assert legs[("Spread", "Indianapolis Colts +7.5")].venue_ref["line"] == 7.5
    assert legs[("Total Points", "Over 48.5")].canonical == CanonicalLeg("FG", "total", None, "over", 48.5)


def test_market_listed_twice_is_read_once(game):
    parsed = parse_legs(game, load("caesars_event_ind_hou.json")["event"])
    ids = [leg.leg_id for leg in parsed]
    assert len(ids) == len(set(ids))


def test_moneyline_team_total_and_first_half(legs):
    assert legs[("Money Line", "Houston Texans")].canonical == CanonicalLeg("FG", "moneyline", None, "home", None)
    assert legs[("Houston Texans Total Points", "Over 28.5")].canonical == CanonicalLeg(
        "FG", "team_total", "Houston Texans", "over", 28.5)
    assert legs[("1st Half Spread", "Houston Texans -7")].canonical == CanonicalLeg("1H", "spread", None, "home", -7.0)
    assert legs[("1st Half - Betting 3 Way", "Tie")].canonical == CanonicalLeg("1H", "win_3way", "Draw", "yes", None)


def test_player_ladder_threshold_is_an_over_half_a_point_lower(legs):
    leg = legs[("Nico Collins - Receiving Yards", "60+")]
    assert leg.canonical == CanonicalLeg("FG", "receiving_yards", "Nico Collins", "over", 59.5)
    assert (leg.american_odds, leg.opposite_american_odds, leg.opposite_leg_id) == (-218.0, None, None)
    assert leg.market_key == legs[("Nico Collins - Receiving Yards", "80+")].market_key   # one threshold per parlay
    assert leg.available_stake == 1500.0 and leg.venue_ref["line"] is None


def test_scorers_name_the_player_and_anytime_takes_several(legs):
    montgomery = legs[("Player to Score a Touchdown", "David Montgomery")]
    collins = legs[("Player to Score a Touchdown", "Nico Collins")]
    assert montgomery.canonical == CanonicalLeg("FG", "anytime_td", "David Montgomery", "yes", None)
    assert montgomery.market_key != collins.market_key
    first = [leg for key, leg in legs.items() if key[0] == "First Touchdown Scorer"]
    assert len({leg.market_key for leg in first}) == 1


def test_boosts_and_inactive_markets_are_left_out(game, legs):
    assert not any(label == "Odds Boost" for label, _ in legs)
    event = load("caesars_event_ind_hou.json")["event"]
    for group in event["keyMarketGroups"]:
        for market in group["markets"]:
            market["active"] = False
    assert parse_legs(game, event) == []


# ---------------------------------------------------------------- parlay quotes

@pytest.mark.parametrize("fixture, sides, decimal_odds", [
    ("caesars_parlay_spread_over.json", [("Spread", "Houston Texans -12.5"), ("Total Points", "Over 43.5")], 3.45),
    ("caesars_parlay_spread_recyds.json", [("Spread", "Houston Texans -12.5"),
                                           ("Nico Collins - Receiving Yards", "60+")], 2.89),
    ("caesars_parlay_two_scorers.json", [("Player to Score a Touchdown", "David Montgomery"),
                                         ("Player to Score a Touchdown", "Nico Collins")], 3.70),
])
def test_recorded_parlay_requests_and_answers(legs, fixture, sides, decimal_odds):
    recorded = load(fixture)
    chosen = [legs[side] for side in sides]
    body = bet_details_body(chosen)
    assert body["combinationSelections"] == []
    sent = [(leg["selectionId"], leg["marketId"], leg["line"]) for leg in body["legs"]]
    assert sent == [(leg["selectionId"], leg["marketId"], leg["line"]) for leg in recorded["request"]["legs"]]
    quote = parse_parlay_quote(recorded["response"], chosen)
    assert quote.best.decimal_odds == pytest.approx(decimal_odds)
    assert quote.best.max_stake is None and quote.leg_prices_american == {}


def test_alt_line_parlay_is_checked_against_the_line_asked(legs):
    recorded = load("caesars_parlay_alt_spread_over.json")
    chosen = [legs[("Spread", "Houston Texans -7.5")], legs[("Total Points", "Over 43.5")]]
    assert parse_parlay_quote(recorded["response"], chosen).best.decimal_odds == pytest.approx(2.79)
    main = [legs[("Spread", "Houston Texans -12.5")], legs[("Total Points", "Over 43.5")]]
    with pytest.raises(VenueTransportError, match="expected -12.5"):
        parse_parlay_quote(recorded["response"], main)


def test_refused_combination_is_a_clean_decline(legs):
    recorded = load("caesars_parlay_refused_same_ladder.json")
    chosen = [legs[("Nico Collins - Receiving Yards", "60+")], legs[("Nico Collins - Receiving Yards", "80+")]]
    assert parse_parlay_quote(recorded["response"], chosen) is None
    assert parse_parlay_quote({"parlays": []}, chosen) is None


# ---------------------------------------------------------------- the venue over a fake transport

class FakeTransport:
    def __init__(self):
        self.gets, self.posts = [], []

    def get_json(self, path):
        self.gets.append(path)
        if path in SCHEDULE_TAB_PATHS:
            return load("caesars_tabs_featured.json")
        if path == f"v4/events/{EVENT_ID}":
            return load("caesars_event_ind_hou.json")
        return None

    def post_bet_details(self, body):
        self.posts.append(body)
        return load("caesars_parlay_spread_over.json")["response"]


def test_venue_reads_both_tabs_caches_and_prices():
    transport = FakeTransport()
    venue = CaesarsVenue(transport_factory=lambda: transport, clock=lambda: 0.0)
    game = venue.list_games()[0]
    venue.list_games()
    assert transport.gets == list(SCHEDULE_TAB_PATHS)
    offered = {(l.market_label, l.side_label): l for l in venue.list_legs(game)}
    venue.list_legs(game)
    venue.list_legs(game, fresh=True)
    assert transport.gets.count(f"v4/events/{EVENT_ID}") == 2
    quote = venue.price_parlay(game, [offered[("Spread", "Houston Texans -12.5")], offered[("Total Points", "Over 43.5")]])
    assert quote.best.decimal_odds == pytest.approx(3.45) and len(transport.posts) == 1


def test_a_pulled_game_has_no_legs(game):
    venue = CaesarsVenue(transport_factory=FakeTransport, clock=lambda: 0.0)
    pulled = type(game)("gone", game.home_team, game.away_team, game.game_start_time)
    assert venue.list_legs(pulled) == []
