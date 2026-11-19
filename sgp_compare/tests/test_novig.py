import json
from datetime import datetime, timezone
from pathlib import Path

import pytest

from sgp_compare.legs import CanonicalLeg
from sgp_compare.odds import american_to_decimal
from sgp_compare.pricing import price_at_venue
from sgp_compare.venues import build_venues
from sgp_compare.venues.base import VenueTransportError
from sgp_compare.venues.novig import (EVENT_MARKETS_QUERY_PATH, NovigVenue, american_odds_of_chance,
                                      event_markets_request_body, parlay_payload, parse_legs, parse_nfl_games,
                                      parse_parlay_quote)

FIXTURES = Path(__file__).parent / "fixtures"
GAME_ID = "01a1c9ea-d600-7e01-8c79-b7a721168527"
BEFORE_KICKOFF = datetime(2026, 11, 19, tzinfo=timezone.utc)
KICKOFF = datetime(2026, 11, 20, 1, 15, tzinfo=timezone.utc)


def load(name):
    return json.loads((FIXTURES / name).read_text())


@pytest.fixture
def game():
    return parse_nfl_games(load("novig_nfl_page.json"), BEFORE_KICKOFF)[0]


@pytest.fixture
def legs(game):
    return {(leg.market_label, leg.side_label): leg for leg in parse_legs(game, load("novig_markets_ind_hou.json"))}


# ---------------------------------------------------------------- games

def test_games_are_the_pregame_cards_soonest_first():
    games = parse_nfl_games(load("novig_nfl_page.json"), BEFORE_KICKOFF)   # the fixture lists the later game first
    assert [(g.venue_game_id, g.home_team, g.away_team, g.game_start_time) for g in games] == [
        (GAME_ID, "Houston Texans", "Indianapolis Colts", KICKOFF),
        ("01a1c9ea-d660-7b32-add7-e68933c6ca57", "Carolina Panthers", "Baltimore Ravens",
         datetime(2026, 11, 22, 18, 0, tzinfo=timezone.utc)),
    ]


def test_started_live_and_offseason_pages_give_no_games():
    assert parse_nfl_games(load("novig_nfl_page.json"), datetime(2026, 12, 31, tzinfo=timezone.utc)) == []
    assert parse_nfl_games({"sections": []}, BEFORE_KICKOFF) == []
    page = load("novig_nfl_page.json")
    for card in page["sections"][0]["content"]["components"]:
        card["eventStatus"] = "OPEN_INGAME"
    assert parse_nfl_games(page, BEFORE_KICKOFF) == []


def test_a_start_time_without_a_zone_is_read_as_utc_not_local():
    page = load("novig_nfl_page.json")
    page["sections"][0]["content"]["components"][1]["scheduledStart"] = "2026-11-20T01:15:00"
    assert parse_nfl_games(page, BEFORE_KICKOFF)[0].game_start_time == KICKOFF


# ---------------------------------------------------------------- prices

@pytest.mark.parametrize("chance, american_odds", [(0.5, 100.0), (0.495, 102.02), (0.52, -108.33), (0.008, 12400.0),
                                                   (0.999, -99900.0)])
def test_a_chance_becomes_american_odds_that_pay_one_over_the_chance(chance, american_odds):
    assert american_odds_of_chance(chance) == american_odds
    assert american_to_decimal(american_odds) == pytest.approx(1.0 / chance, rel=1e-4)


@pytest.mark.parametrize("chance", [None, 0, 1, 1.0, 1.2, -0.1, "0.5", True])
def test_a_missing_or_impossible_chance_has_no_price(chance):
    assert american_odds_of_chance(chance) is None


# ---------------------------------------------------------------- legs

def test_spread_leg_prices_sides_and_ids(legs):
    leg = legs[("Spread", "Houston Texans -12.5")]
    assert (leg.american_odds, leg.opposite_american_odds, leg.available_stake) == (102.02, -104.08, None)
    assert leg.canonical == CanonicalLeg("FG", "spread", None, "home", -12.5)
    assert leg.leg_id == "01a1c9ea-d65b-7fe3-9fac-b58dd46b2cad"
    assert leg.opposite_leg_id == legs[("Spread", "Indianapolis Colts +12.5")].leg_id
    assert legs[("Spread", "Indianapolis Colts +12.5")].canonical == CanonicalLeg("FG", "spread", None, "away", 12.5)
    assert leg.venue_ref["outcomeId"] == leg.leg_id
    assert (leg.market_key, leg.category, leg.is_main_line) == ("SPREAD", "Game Lines", True)


def test_alt_lines_share_the_market_key_and_are_not_main(legs):
    alt = legs[("Spread", "Houston Texans -18.5")]
    assert (alt.market_key, alt.is_main_line) == ("SPREAD", False)
    assert alt.canonical == CanonicalLeg("FG", "spread", None, "home", -18.5)
    assert (legs[("Total", "Over 35.5")].market_key, legs[("Total", "Over 35.5")].is_main_line) == ("TOTAL", False)


def test_moneyline_totals_and_team_totals(legs):
    assert legs[("Moneyline", "Houston Texans")].canonical == CanonicalLeg("FG", "moneyline", None, "home", None)
    assert legs[("Moneyline", "Indianapolis Colts")].american_odds == 545.16
    assert legs[("Total", "Over 43.5")].canonical == CanonicalLeg("FG", "total", None, "over", 43.5)
    team_total = legs[("Indianapolis Colts Team Total", "Under 15.5")]
    assert team_total.canonical == CanonicalLeg("FG", "team_total", "Indianapolis Colts", "under", 15.5)
    assert team_total.market_key != legs[("Total", "Over 43.5")].market_key


def test_first_half_lines_are_their_own_markets(legs):
    assert legs[("1st Half Spread", "Indianapolis Colts +6.5")].canonical == CanonicalLeg("1H", "spread", None, "away", 6.5)
    assert legs[("1st Half Moneyline", "Houston Texans")].canonical == CanonicalLeg("1H", "moneyline", None, "home", None)
    over = legs[("1st Half Total", "Over 21.5")]          # Novig lists this market's Under first
    assert (over.canonical, over.american_odds) == (CanonicalLeg("1H", "total", None, "over", 21.5), -110.53)
    assert over.market_key != legs[("Total", "Over 43.5")].market_key


def test_player_props_name_the_player_and_stat(legs):
    leg = legs[("Nico Collins Receiving Yards", "Over 70.5")]
    assert leg.canonical == CanonicalLeg("FG", "receiving_yards", "Nico Collins", "over", 70.5)
    assert (leg.category, leg.is_main_line) == ("Receiving Yards", True)
    assert legs[("C.J. Stroud Passing Yards", "Under 257.5")].canonical == CanonicalLeg(
        "FG", "passing_yards", "C.J. Stroud", "under", 257.5)


def test_touchdown_counts_are_yes_no_bets_named_as_the_other_venues_name_them(legs):
    anytime = legs[("Jonathan Taylor Anytime Touchdown", "YES")]
    assert anytime.canonical == CanonicalLeg("FG", "anytime_td", "Jonathan Taylor", "yes", None)
    assert legs[("Jonathan Taylor Anytime Touchdown", "NO")].leg_id == anytime.opposite_leg_id
    two_plus = legs[("Jonathan Taylor 2+ Touchdowns", "YES")]
    assert two_plus.canonical == CanonicalLeg("FG", "to_score_at_least_2_touchdowns", "Jonathan Taylor", "yes", None)
    assert two_plus.market_key != anytime.market_key
    assert legs[("Nikko Remigio First Touchdown", "NO")].canonical == CanonicalLeg(
        "FG", "to_score_first_touchdown", "Nikko Remigio", "no", None)


def test_protected_twins_are_kept_apart_from_the_plain_market(legs):
    plain, protected = legs[("Nico Collins Receiving Yards", "Over 70.5")], legs[
        ("Nico Collins Receiving Yards (Protected)", "Over 70.5")]
    assert protected.canonical == CanonicalLeg("FG", "receiving_yards_protected", "Nico Collins", "over", 70.5)
    assert protected.market_key != plain.market_key and protected.category != plain.category
    assert protected.leg_id != plain.leg_id


def test_locked_and_unpriced_lines_are_left_out_but_a_one_sided_line_keeps_both_sides(legs):
    labels = {market_label for market_label, _ in legs}
    assert "Will Mallory Anytime Touchdown" not in labels                  # locked
    assert "Quez Watkins Anytime Touchdown (Protected)" not in labels      # nobody offers either side
    yes, no = legs[("Quez Watkins Anytime Touchdown", "YES")], legs[("Quez Watkins Anytime Touchdown", "NO")]
    assert (yes.american_odds, yes.opposite_american_odds) == (2757.14, None)
    assert (no.american_odds, no.opposite_american_odds) == (None, 2757.14)


def test_closed_markets_and_an_unknown_game_give_no_legs(game):
    markets = load("novig_markets_ind_hou.json")
    for market in markets["data"]["event"][0]["markets"]:
        market["status"] = "CLOSED"
    assert parse_legs(game, markets) == []
    assert parse_legs(game, {"data": {"event": []}}) == []


def test_a_spread_whose_text_disagrees_with_its_strike_or_names_an_unknown_team_is_left_out(game):
    markets = load("novig_markets_ind_hou.json")
    event = markets["data"]["event"][0]
    main_spread = next(m for m in event["markets"] if m["type"] == "SPREAD" and m["strike"] == -12.5)
    main_spread["outcomes"][0]["description"] = "HOU -3.5"
    moneyline = next(m for m in event["markets"] if m["type"] == "MONEY")
    moneyline["outcomes"][1]["competitor"]["id"] = "some-other-team"
    parsed = {(leg.market_label, leg.side_label) for leg in parse_legs(game, markets)}
    assert not any(side.endswith("12.5") or side.endswith("3.5") for label, side in parsed if label == "Spread")
    assert not any(label == "Moneyline" for label, _ in parsed)
    assert ("Spread", "Houston Texans -18.5") in parsed


def test_a_game_without_two_team_ids_is_an_error_not_a_guess(game):
    markets = load("novig_markets_ind_hou.json")
    markets["data"]["event"][0]["game"]["homeTeam"]["id"] = markets["data"]["event"][0]["game"]["awayTeam"]["id"]
    with pytest.raises(ValueError, match="two different team ids"):
        parse_legs(game, markets)


# ---------------------------------------------------------------- parlay quotes

@pytest.mark.parametrize("fixture, sides, decimal_odds", [
    ("novig_parlay_spread_over.json", [("Spread", "Houston Texans -12.5"), ("Total", "Over 43.5")], 3.6900),
    ("novig_parlay_spread_recyds.json", [("Spread", "Houston Texans -12.5"), ("Nico Collins Receiving Yards", "Over 70.5")],
     3.8168),
])
def test_recorded_parlay_requests_and_answers(legs, fixture, sides, decimal_odds):
    recorded = load(fixture)
    assert parlay_payload([legs[side] for side in sides]) == recorded["request"]
    quote = parse_parlay_quote(recorded["response"])
    assert [(round(tier.decimal_odds, 4), tier.max_stake) for tier in quote.tiers] == [(decimal_odds, None)]
    assert quote.leg_prices_american == {}
    assert quote.raw == recorded["response"]


@pytest.mark.parametrize("response", [{"price": "0"}, {"price": "1.0"}, {"price": "abc"}, {"price": None}, {}])
def test_a_quote_without_a_usable_price_is_an_error(response):
    with pytest.raises((ValueError, KeyError, TypeError)):
        parse_parlay_quote(response)


def test_the_committed_query_is_sent_for_this_game():
    body = json.loads(event_markets_request_body(EVENT_MARKETS_QUERY_PATH.read_text(), GAME_ID))
    assert body["operationName"] == "EventMarkets_Query"
    assert body["variables"]["eventId"] == GAME_ID
    assert "marketVisibleWhere" in body["variables"]


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
        self.answers = answers          # URL fragment -> FakeResponse, Exception, or a function of the request
        self.calls = []

    def request(self, method, url, timeout=None, **kwargs):
        self.calls.append((method, url, kwargs))
        for fragment, answer in self.answers.items():
            if fragment in url:
                if isinstance(answer, Exception):
                    raise answer
                return answer(kwargs) if callable(answer) else answer
        raise AssertionError(f"unexpected {method} {url}")


def live_answers(**overrides):
    answers = {"/trading/NFL/page": FakeResponse(200, load("novig_nfl_page.json")),
               "/v1/graphql": FakeResponse(200, load("novig_markets_ind_hou.json")),
               "/rfq/v2/quote": FakeResponse(200, load("novig_parlay_spread_over.json")["response"])}
    answers.update(overrides)
    return answers


def venue_with(answers, wall_time=BEFORE_KICKOFF):
    session = FakeSession(answers)
    clock = {"t": 0.0}
    venue = NovigVenue(session_factory=lambda: session, clock=lambda: clock["t"], wall_clock=lambda: wall_time)
    return venue, session, clock


def chosen_legs(venue, game, *sides):
    by_side = {(leg.market_label, leg.side_label): leg for leg in venue.list_legs(game)}
    return [by_side[side] for side in sides]


def test_venue_reads_games_legs_and_prices_and_only_ever_calls_the_three_read_endpoints():
    venue, session, _ = venue_with(live_answers())
    game = venue.list_games()[0]
    legs = chosen_legs(venue, game, ("Spread", "Houston Texans -12.5"), ("Total", "Over 43.5"))
    quote = venue.price_parlay(game, legs)
    assert round(quote.best.decimal_odds, 4) == 3.69
    assert [(method, url) for method, url, _ in session.calls] == [
        ("GET", "https://api.novig.us/nbx/v1/trading/NFL/page"),
        ("POST", "https://api.novig.us/v1/graphql"),
        ("POST", "https://api.novig.us/rfq/v2/quote/unauthenticated"),
    ]
    graphql_body = json.loads(session.calls[1][2]["data"])
    assert graphql_body["variables"]["eventId"] == GAME_ID
    assert session.calls[2][2]["json"] == load("novig_parlay_spread_over.json")["request"]


def test_a_refused_combination_is_a_decline_not_an_error():
    declined = load("novig_parlay_both_sides_declined.json")
    venue, _, _ = venue_with(live_answers(**{"/rfq/v2/quote": FakeResponse(declined["status"], declined["response"])}))
    game = venue.list_games()[0]
    legs = chosen_legs(venue, game, ("Total", "Over 43.5"), ("Total", "Under 43.5"))
    assert venue.price_parlay(game, legs) is None


@pytest.mark.parametrize("answer, message", [
    (FakeResponse(400, "Bad Request"), "expected HTTP 200, got 400"),
    (FakeResponse(404, "Cannot POST /rfq/v2/quote/unauthenticated"), "expected HTTP 200, got 404"),
    (FakeResponse(200, "<html>"), "not JSON"),
    (FakeResponse(200, []), "expected a JSON object"),
    (FakeResponse(200, {"price": "1.5"}), "unexpected shape"),
    (ConnectionError("reset"), "ConnectionError"),
])
def test_any_other_quote_failure_is_a_transport_error(answer, message):
    venue, _, _ = venue_with(live_answers(**{"/rfq/v2/quote": answer}))
    game = venue.list_games()[0]
    legs = chosen_legs(venue, game, ("Spread", "Houston Texans -12.5"), ("Total", "Over 43.5"))
    with pytest.raises(VenueTransportError, match=message):
        venue.price_parlay(game, legs)


def test_a_query_novig_no_longer_allows_says_how_to_refresh_it():
    refused = {"errors": [{"message": "query is not allowed", "extensions": {"code": "validation-failed"}}]}
    venue, _, _ = venue_with(live_answers(**{"/v1/graphql": FakeResponse(200, refused)}))
    with pytest.raises(VenueTransportError, match="refresh_novig_query.py"):
        venue.list_legs(venue.list_games()[0])


def test_a_missing_query_file_and_a_garbled_page_are_transport_errors(tmp_path):
    session = FakeSession(live_answers())
    venue = NovigVenue(session_factory=lambda: session, clock=lambda: 0.0, wall_clock=lambda: BEFORE_KICKOFF,
                       event_markets_query_path=tmp_path / "missing.json")
    with pytest.raises(VenueTransportError, match="cannot read the committed query"):
        venue.list_legs(venue.list_games()[0])
    broken, _, _ = venue_with(live_answers(**{"/trading/NFL/page": FakeResponse(200, {"sections": "nope"})}))
    with pytest.raises(VenueTransportError, match="unexpected shape"):
        broken.list_games()


def test_the_leg_picker_is_cached_but_pricing_reads_fresh_and_games_drop_at_kickoff():
    venue, session, clock = venue_with(live_answers())
    game = venue.list_games()[0]
    venue.list_legs(game)
    venue.list_legs(game)
    graphql_calls = lambda: sum(1 for _, url, _ in session.calls if "graphql" in url)
    assert graphql_calls() == 1
    venue.list_legs(game, fresh=True)
    assert graphql_calls() == 2
    clock["t"] = 21.0
    venue.list_legs(game)
    assert graphql_calls() == 3
    started, _, _ = venue_with(live_answers(), wall_time=KICKOFF)
    assert [g.home_team for g in started.list_games()] == ["Carolina Panthers"]


def test_a_dropped_connection_gets_a_new_session_next_time():
    sessions = [FakeSession(live_answers(**{"/trading/NFL/page": ConnectionError("reset")})), FakeSession(live_answers())]
    made = iter(sessions)
    venue = NovigVenue(session_factory=lambda: next(made), clock=lambda: 0.0, wall_clock=lambda: BEFORE_KICKOFF)
    with pytest.raises(VenueTransportError):
        venue.list_games()
    assert len(venue.list_games()) == 2


# ---------------------------------------------------------------- inside the page's pricing

def test_novig_is_one_of_the_pages_venues():
    assert build_venues()["novig"].display_name == "Novig"


def test_the_core_prices_a_novig_parlay_and_removes_the_vig_from_the_four_corners():
    corner_price_of_ids = {}
    markets = load("novig_markets_ind_hou.json")["data"]["event"][0]["markets"]
    spread = next(m for m in markets if m["type"] == "SPREAD" and m["strike"] == -12.5)["outcomes"]
    total = next(m for m in markets if m["type"] == "TOTAL" and m["strike"] == 43.5)["outcomes"]
    # The four corners the 2026-11-19 recon priced: home/away spread x over/under.
    for (spread_index, total_index), price in {(0, 0): "0.271", (0, 1): "0.225", (1, 0): "0.251", (1, 1): "0.26"}.items():
        corner_price_of_ids[frozenset((spread[spread_index]["id"], total[total_index]["id"]))] = price

    def quote(kwargs):
        return FakeResponse(200, {"price": corner_price_of_ids[frozenset(kwargs["json"]["outcomeIds"])],
                                  "expiresAt": 1795056416704, "serverNow": 1795056401273})

    venue, _, _ = venue_with(live_answers(**{"/rfq/v2/quote": quote}))
    game = venue.list_games()[0]
    result = price_at_venue(venue, game, [spread[0]["id"], total[0]["id"]])
    assert result.status == "priced"
    assert round(result.quote.best.decimal_odds, 4) == 3.69
    assert sorted(round(1.0 / odds, 3) for odds in result.partition_decimal_odds) == [0.225, 0.251, 0.26, 0.271]
    assert result.report is not None
