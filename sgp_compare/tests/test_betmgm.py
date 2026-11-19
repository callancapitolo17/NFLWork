import json
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest

from sgp_compare.venues.base import VenueTransportError
from sgp_compare.venues.betmgm import (BetMgmVenue, parlay_payload, parse_legs, parse_nfl_games,
                                       parse_parlay_quote)

FIXTURES = Path(__file__).parent / "fixtures"
BEFORE_KICKOFF = datetime(2026, 11, 19, tzinfo=timezone.utc)
KICKOFF = datetime(2026, 11, 20, 1, 15, tzinfo=timezone.utc)
GAME_ID = "6:44039"


def load(name):
    return json.loads((FIXTURES / name).read_text())


@pytest.fixture
def game():
    return next(g for g in parse_nfl_games(load("betmgm_fixtures_nfl.json"), BEFORE_KICKOFF)
                if g.venue_game_id == GAME_ID)


@pytest.fixture
def legs(game):
    return {(leg.market_label, leg.side_label): leg
            for leg in parse_legs(game, load("betmgm_markets_colts_texans.json"))}


def canonical_tuple(leg):
    c = leg.canonical
    return (c.period, c.stat, c.subject, c.side, c.line)


# ---------------------------------------------------------------- games

def test_games_are_nfl_only_with_home_from_participants(game):
    games = parse_nfl_games(load("betmgm_fixtures_nfl.json"), BEFORE_KICKOFF)
    assert games == sorted(games, key=lambda g: g.game_start_time)
    assert (game.home_team, game.away_team, game.game_start_time) == ("Houston Texans", "Indianapolis Colts", KICKOFF)


def test_started_and_non_nfl_fixtures_are_left_out():
    payload = load("betmgm_fixtures_nfl.json")
    assert parse_nfl_games(payload, datetime(2027, 3, 1, tzinfo=timezone.utc)) == []
    for fixture in payload["fixtures"]:
        fixture["competition"]["name"]["value"] = "NCAA, Regular Season"
    assert parse_nfl_games(payload, BEFORE_KICKOFF) == []


# ---------------------------------------------------------------- legs

def test_only_visible_bet_builder_markets_are_offered(legs):
    labels = {label for label, _ in legs}
    assert "C.J. Stroud - Interceptions Thrown" not in labels        # isBetBuilder false
    assert "Houston Texans - Total Touchdowns (Over/Under)" not in labels
    assert "Tank Dell - Reception Yards" in labels


def test_game_lines_map_to_canonical_legs(legs):
    spread = legs[("Spread", "Houston Texans (-12.5)")]
    assert canonical_tuple(spread) == ("FG", "spread", None, "home", -12.5)
    assert (spread.american_odds, spread.opposite_american_odds) == (-110.0, -108.0)
    assert spread.opposite_leg_id == legs[("Spread", "Indianapolis Colts (+12.5)")].leg_id
    assert spread.is_main_line and spread.venue_ref == {"fixtureId": GAME_ID, "optionMarketId": 4585337,
                                                        "optionId": 7552787}
    assert canonical_tuple(legs[("Totals", "Over 44")]) == ("FG", "total", None, "over", 44.0)
    assert canonical_tuple(legs[("Moneyline", "Indianapolis Colts")]) == ("FG", "moneyline", None, "away", None)
    assert canonical_tuple(legs[("1st half spread", "Houston Texans (-6.5)")]) == ("1H", "spread", None, "home", -6.5)
    assert canonical_tuple(legs[("1st quarter totals", "Under 7.5")]) == ("1Q", "total", None, "under", 7.5)
    team_total = legs[("Houston Texans: 1st half points", "Over 14.5")]
    assert canonical_tuple(team_total) == ("1H", "team_total", "Houston Texans", "over", 14.5)
    assert team_total.category == "Team Totals"


def test_alt_lines_share_the_main_lines_market_key(legs):
    main = legs[("Spread", "Houston Texans (-12.5)")]
    alt = legs[("Spread", "Houston Texans (-10)")]
    assert main.market_key == alt.market_key and not alt.is_main_line
    assert main.market_key != legs[("1st half spread", "Houston Texans (-6.5)")].market_key
    assert legs[("Totals", "Over 44")].market_key == legs[("Totals", "Over 45.5")].market_key


def test_player_over_unders_name_the_player_and_stat(legs):
    dell = legs[("Tank Dell - Reception Yards", "Over 35.5")]
    assert canonical_tuple(dell) == ("FG", "receiving_yards", "Tank Dell", "over", 35.5)
    assert (dell.american_odds, dell.opposite_american_odds, dell.category) == (-115.0, -115.0, "Player Props")
    assert canonical_tuple(legs[("C.J. Stroud - Passing Yards", "Under 262.5")]) == \
        ("FG", "passing_yards", "C.J. Stroud", "under", 262.5)
    assert canonical_tuple(legs[("Kayshon Boutte - Receptions", "Over 2.5")])[1] == "receptions"


def test_milestones_are_overs_at_half_a_point_below_with_no_other_side(legs):
    milestone = legs[("Kayshon Boutte to record 50+ receiving yards", "Yes")]
    assert canonical_tuple(milestone) == ("FG", "receiving_yards", "Kayshon Boutte", "over", 49.5)
    assert (milestone.american_odds, milestone.opposite_american_odds, milestone.opposite_leg_id) == (260.0, None, None)
    first_half = legs[("Tank Dell to record 25+ receiving yards: 1st half", "Yes")]
    assert canonical_tuple(first_half) == ("1H", "receiving_yards", "Tank Dell", "over", 24.5)
    # "Over 75" does not say whether exactly 75 wins, so it is listed but unmatched.
    assert legs[("David Montgomery Over 75 Rushing and Reception Yards", "Yes")].canonical is None


def test_touchdown_props(legs):
    anytime = legs[("Deion Burks to score 1+ TDs", "Yes")]
    assert canonical_tuple(anytime) == ("FG", "anytime_td", "Deion Burks", "yes", None)
    assert canonical_tuple(legs[("Nikko Remigio to score 3+ TDs", "Yes")]) == \
        ("FG", "touchdowns", "Nikko Remigio", "over", 2.5)
    first = legs[("First Touchdown Scorer", "David Montgomery")]
    assert canonical_tuple(first) == ("FG", "first_td", "David Montgomery", "yes", None)
    assert first.opposite_leg_id is None


def test_correct_score_is_listed_but_unmatched(legs):
    assert legs[("Correct score in full time", "Indianapolis Colts 24 - Houston Texans 17")].canonical is None


def test_a_suspended_option_is_listed_without_a_price(game):
    markets = load("betmgm_markets_colts_texans.json")
    spread = next(m for m in markets["fixtures"][0]["optionMarkets"] if m["id"] == 4585337)
    spread["options"][0]["status"] = "Suspended"           # the Colts side
    parsed = {(l.market_label, l.side_label): l for l in parse_legs(game, markets)}
    houston = parsed[("Spread", "Houston Texans (-12.5)")]
    assert (houston.american_odds, houston.opposite_american_odds) == (-110.0, None)
    assert parsed[("Spread", "Indianapolis Colts (+12.5)")].american_odds is None


def test_hidden_markets_are_left_out(game):
    markets = load("betmgm_markets_colts_texans.json")
    for market in markets["fixtures"][0]["optionMarkets"]:
        market["status"] = "Hidden"
    assert parse_legs(game, markets) == []


# ---------------------------------------------------------------- parlay quotes

@pytest.mark.parametrize("fixture, sides, decimal_odds", [
    ("betmgm_parlay_spread_over.json", [("Spread", "Houston Texans (-12.5)"), ("Totals", "Over 44")], 3.5),
    ("betmgm_parlay_spread_recyds.json", [("Spread", "Houston Texans (-12.5)"),
                                          ("Tank Dell - Reception Yards", "Over 35.5")], 3.2),
])
def test_recorded_parlay_requests_and_answers(legs, fixture, sides, decimal_odds):
    recorded = load(fixture)
    chosen = [legs[side] for side in sides]
    group_id = recorded["request"]["tv2Picks"][0]["pickGroupId"]
    assert parlay_payload(chosen, group_id) == recorded["request"]
    quote = parse_parlay_quote(recorded["response"])
    assert [(round(t.decimal_odds, 4), t.max_stake) for t in quote.tiers] == [(decimal_odds, None)]
    assert quote.leg_prices_american == {}
    assert "fixturePage" not in quote.raw and "betBuilderPricingGroups" in quote.raw


def test_anytime_td_and_moneyline_parlay_prices():
    quote = parse_parlay_quote(load("betmgm_parlay_anytd_ml.json")["response"])
    assert round(quote.best.decimal_odds, 4) == 11.5


@pytest.mark.parametrize("fixture", ["betmgm_parlay_over_under.json", "betmgm_parlay_single_spread.json"])
def test_declines_are_none(fixture):
    assert parse_parlay_quote(load(fixture)["response"]) is None


def test_a_suspended_group_with_odds_is_a_decline():
    response = load("betmgm_parlay_spread_over.json")["response"]
    group = next(iter(response["betBuilderPricingGroups"].values()))
    group["suspensionState"] = "MarketSuspended"
    assert parse_parlay_quote(response) is None
    assert parse_parlay_quote({"betBuilderPricingGroups": {}}) is None


# ---------------------------------------------------------------- the live venue, on a fake client

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
        self.answers = answers          # url fragment -> list of FakeResponse / Exception, used in turn
        self.calls = []

    def _answer(self, method, url):
        self.calls.append((method, url))
        for fragment, queue in self.answers.items():
            if fragment in url:
                answer = queue.pop(0) if len(queue) > 1 else queue[0]
                if isinstance(answer, Exception):
                    raise answer
                return answer
        raise AssertionError(f"unexpected {method} {url}")

    def get(self, url, timeout=None):
        return self._answer("GET", url)

    def post(self, url, json=None, timeout=None):
        self.posted = json
        return self._answer("POST", url)


class FakeClient:
    _base = "https://www.pa.betmgm.com/cds-api"

    def __init__(self, session):
        self.session = session
        self.harvests = []

    def accessid(self, refresh=False):
        self.harvests.append(refresh)
        return "id"

    def _common_q(self):
        return "x-bwin-accessid=id&lang=en"


def live_answers(**overrides):
    answers = {"sportIds=11": [FakeResponse(200, load("betmgm_fixtures_nfl.json"))],
               "fixtureIds=": [FakeResponse(200, load("betmgm_markets_colts_texans.json"))],
               "/picks": [FakeResponse(200, load("betmgm_parlay_spread_over.json")["response"])]}
    answers.update(overrides)
    return answers


def venue_with(answers, wall_time=BEFORE_KICKOFF):
    client = FakeClient(FakeSession(answers))
    clock = {"t": 0.0}
    venue = BetMgmVenue(client_factory=lambda: client, clock=lambda: clock["t"], wall_clock=lambda: wall_time)
    return venue, client, clock


def pick(venue, game, *sides):
    offered = {(l.market_label, l.side_label): l for l in venue.list_legs(game)}
    return [offered[side] for side in sides]


def test_venue_reads_games_legs_and_prices():
    venue, client, _ = venue_with(live_answers())
    game = next(g for g in venue.list_games() if g.venue_game_id == GAME_ID)
    chosen = pick(venue, game, ("Spread", "Houston Texans (-12.5)"), ("Totals", "Over 44"))
    quote = venue.price_parlay(game, chosen)
    assert quote.best.decimal_odds == 3.5
    assert [method for method, _ in client.session.calls] == ["GET", "GET", "POST"]
    assert [p["optionId"] for p in client.session.posted["tv2Picks"]] == [7552787, 7537006]
    assert client.session.posted["tv1Picks"] == []


def test_a_one_leg_parlay_is_declined_without_a_request():
    venue, client, _ = venue_with(live_answers())
    game = next(g for g in venue.list_games() if g.venue_game_id == GAME_ID)
    assert venue.price_parlay(game, pick(venue, game, ("Totals", "Over 44"))) is None
    assert not any(method == "POST" for method, _ in client.session.calls)


def test_schedule_is_cached_but_a_started_game_drops_out_at_once():
    client = FakeClient(FakeSession(live_answers()))
    now = {"wall": BEFORE_KICKOFF}
    venue = BetMgmVenue(client_factory=lambda: client, clock=lambda: 0.0, wall_clock=lambda: now["wall"])
    assert any(g.venue_game_id == GAME_ID for g in venue.list_games())
    now["wall"] = KICKOFF + timedelta(seconds=1)
    assert all(g.venue_game_id != GAME_ID for g in venue.list_games())
    assert sum("sportIds=11" in url for _, url in client.session.calls) == 1


def test_market_cache_and_fresh_reads():
    venue, client, clock = venue_with(live_answers())
    game = venue.list_games()[0]
    venue.list_legs(game)
    venue.list_legs(game)
    venue.list_legs(game, fresh=True)
    clock["t"] = 1000.0
    venue.list_legs(game)
    assert sum("fixtureIds=" in url for _, url in client.session.calls) == 3


def test_a_stale_accessid_is_harvested_again_once():
    gate = FakeResponse(400, '{"message":"Access id is not valid"}')
    venue, client, _ = venue_with(live_answers(**{"sportIds=11": [gate, FakeResponse(200, load("betmgm_fixtures_nfl.json"))]}))
    assert venue.list_games()
    assert client.harvests == [False, True]


@pytest.mark.parametrize("answer, message", [
    (FakeResponse(403, "blocked"), "expected HTTP 200, got 403: blocked"),
    (FakeResponse(400, "Access id is not valid"), "expected HTTP 200, got 400"),
    (FakeResponse(200, "<html>"), "response is not JSON"),
    (FakeResponse(200, []), "expected a JSON object, got list"),
    (RuntimeError("connection reset"), "RuntimeError: connection reset"),
])
def test_transport_failures_raise_venue_transport_error(answer, message):
    venue, _, _ = venue_with(live_answers(**{"sportIds=11": [answer]}))
    with pytest.raises(VenueTransportError, match=message):
        venue.list_games()


def test_a_failed_accessid_harvest_is_a_transport_error_and_rebuilds_the_client():
    built = []

    class NoAccess(FakeClient):
        def accessid(self, refresh=False):
            raise RuntimeError("no publicAccessId in clientconfig")

    def factory():
        built.append(1)
        return NoAccess(FakeSession(live_answers()))
    venue = BetMgmVenue(client_factory=factory)
    for _ in range(2):
        with pytest.raises(VenueTransportError, match="no publicAccessId"):
            venue.list_games()
    assert len(built) == 2


def test_an_answer_with_an_unexpected_shape_is_a_transport_error_not_a_crash():
    weird = {"fixtures": [{"competition": {"name": {"value": "NFL"}}, "participants": "x"}]}
    venue, _, _ = venue_with(live_answers(**{"sportIds=11": [FakeResponse(200, weird)]}))
    with pytest.raises(VenueTransportError, match="fixtures answer has an unexpected shape"):
        venue.list_games()


def test_the_core_prices_a_betmgm_parlay_and_its_partition():
    from sgp_compare.pricing import price_at_venue
    venue, client, _ = venue_with(live_answers())
    game = next(g for g in venue.list_games() if g.venue_game_id == GAME_ID)
    chosen = pick(venue, game, ("Spread", "Houston Texans (-12.5)"), ("Tank Dell - Reception Yards", "Over 35.5"))
    result = price_at_venue(venue, game, [leg.leg_id for leg in chosen])
    assert result.status == "priced" and result.quote.best.decimal_odds == 3.5
    assert len(result.partition_decimal_odds) == 4                  # the parlay + 3 flipped combinations
    assert sum(method == "POST" for method, _ in client.session.calls) == 4
