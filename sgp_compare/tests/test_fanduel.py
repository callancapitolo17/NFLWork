"""FanDuel adapter: parsing off trimmed recon captures (IND @ HOU, 2026-11-19)
and the venue's request handling off a fake session. No network."""
import json
from datetime import datetime, timezone
from pathlib import Path

import pytest

from sgp_compare.venues.base import VenueTransportError
from sgp_compare.venues.fanduel import (EXTRA_TABS, PRIMARY_TAB, FanDuelVenue, implybets_payload,
                                        markets_of_pages, parse_legs, parse_nfl_games, parse_parlay_quote)

FIXTURES = Path(__file__).parent / "fixtures"
BEFORE_KICKOFF = datetime(2026, 11, 19, tzinfo=timezone.utc)
API = {"headers": {"x-application": "test"}, "app_key": "test", "scan_url": "https://scan.test/search",
       "event_page_url": "https://api.test/event-page", "imply_url": "https://sib.test/implyBets"}


def load(name):
    return json.loads((FIXTURES / name).read_text())


@pytest.fixture
def game():
    return parse_nfl_games(load("fanduel_scan.json"), BEFORE_KICKOFF)[0]


@pytest.fixture
def legs(game):
    return parse_legs(game, markets_of_pages([load("fanduel_event_ind_hou.json")]))


@pytest.fixture
def by_label(legs):
    return {(leg.market_label, leg.side_label): leg for leg in legs}


def test_games_are_away_at_home_soonest_first(game):
    games = parse_nfl_games(load("fanduel_scan.json"), BEFORE_KICKOFF)
    assert len(games) == 3                       # the futures entry with no " @ " is left out
    assert (game.venue_game_id, game.home_team, game.away_team) == ("36305232", "Houston Texans", "Indianapolis Colts")
    assert game.game_start_time == datetime(2026, 11, 20, 1, 16, tzinfo=timezone.utc)
    assert parse_nfl_games(load("fanduel_scan.json"), datetime(2027, 1, 1, tzinfo=timezone.utc)) == []


def test_non_sgp_markets_are_left_out(legs):
    assert not any(leg.market_label == "Team to Score First" for leg in legs)


def test_suspended_runners_are_left_out(legs):
    assert not any(leg.side_label == "Ashton Dulin" for leg in legs)


def test_spread_sides_lines_and_pairing(by_label):
    home = by_label[("Spread", "Houston Texans -12.5")]
    away = by_label[("Spread", "Indianapolis Colts +12.5")]
    assert (home.canonical.period, home.canonical.stat, home.canonical.side, home.canonical.line) == ("FG", "spread", "home", -12.5)
    assert home.opposite_leg_id == away.leg_id and home.opposite_american_odds == away.american_odds == -110.0
    assert home.venue_ref == {"marketId": home.market_key, "selectionId": 179823}
    assert home.is_main_line


def test_moneyline_prices_come_from_true_decimal_odds(by_label):
    home = by_label[("Moneyline", "Houston Texans")]
    assert home.canonical.side == "home" and home.american_odds == -800.0 and home.opposite_american_odds == 570.0


def test_leg_id_carries_the_line_so_a_moved_line_is_a_new_leg(by_label):
    over = by_label[("Total Points", "Over 43.5")]
    assert over.leg_id.endswith(":43.5")
    assert over.leg_id != by_label[("HOU Texans Total Points", "Over 27.5")].leg_id   # same selectionId, other market


def test_alternate_ladders_pair_by_meaning(by_label):
    alt = by_label[("Alternate Spread", "Houston Texans (-7.5)")]
    assert alt.canonical.line == -7.5 and not alt.is_main_line
    assert alt.opposite_leg_id == by_label[("Alternate Spread", "Indianapolis Colts (+7.5)")].leg_id
    assert by_label[("Alternate Spread", "Houston Texans (+13.5)")].opposite_leg_id is None   # IND -13.5 not listed
    over = by_label[("Alternate Total Points", "Over (47.5)")]
    assert over.opposite_leg_id == by_label[("Alternate Total Points", "Under (47.5)")].leg_id
    assert by_label[("Alternate Total Points", "Over (19.5)")].opposite_leg_id is None


def test_team_totals_name_the_team(by_label):
    full = by_label[("HOU Texans Total Points", "Over 27.5")].canonical
    half = by_label[("1st Half IND Colts Total Points", "Under 6.5")].canonical
    assert (full.period, full.stat, full.subject, full.line) == ("FG", "team_total", "Houston Texans", 27.5)
    assert (half.period, half.subject, half.side) == ("1H", "Indianapolis Colts", "under")


def test_quarter_and_half_periods(by_label):
    assert by_label[("2nd Quarter Spread", "Houston Texans -3.5")].canonical.period == "2Q"
    assert by_label[("1st Half Total", "Over 21.5")].canonical.period == "1H"


def test_three_way_winner_has_no_other_side(by_label):
    tie = by_label[("1st Quarter Winner (3-Way)", "Tie")]
    assert (tie.canonical.stat, tie.canonical.subject, tie.canonical.side) == ("win_3way", "Draw", "yes")
    assert tie.opposite_leg_id is None


def test_player_over_under_pairs(by_label):
    over = by_label[("Nico Collins - Total Receptions", "Nico Collins Over 4.5")]
    assert (over.canonical.stat, over.canonical.subject, over.canonical.side, over.canonical.line) == \
        ("receptions", "Nico Collins", "over", 4.5)
    assert over.opposite_leg_id == by_label[("Nico Collins - Total Receptions", "Nico Collins Under 4.5")].leg_id


def test_player_ladder_rung_is_an_over_with_no_other_side(by_label):
    rung = by_label[("Nico Collins - Alt Receptions", "Nico Collins 5+ Receptions")]
    assert (rung.canonical.stat, rung.canonical.side, rung.canonical.line) == ("receptions", "over", 4.5)
    assert rung.opposite_leg_id is None and not rung.is_main_line


def test_touchdown_scorers_use_shared_stat_names(by_label):
    td = by_label[("Any Time Touchdown Scorer", "David Montgomery")].canonical
    two = by_label[("To Score 2+ Touchdowns", "Nico Collins")].canonical
    assert (td.stat, td.subject, td.side) == ("anytime_td", "David Montgomery", "yes")
    assert two.stat == "to_score_at_least_2_touchdowns"


def test_unknown_market_is_offered_without_canonical(by_label):
    drive = by_label[("David Montgomery Drive 1 Receiving Yds HOU Texans", "David Montgomery 5+ Receiving Yards")]
    assert drive.canonical is None and drive.category == "Drives"


def test_markets_merge_across_tabs_first_seen_wins():
    page = load("fanduel_event_ind_hou.json")
    assert len(markets_of_pages([page, page])) == len(markets_of_pages([page]))


def test_parlay_quote_reads_the_sgm_entry_and_max_stake(by_label):
    legs = [by_label[("Moneyline", "Indianapolis Colts")], by_label[("Total Points", "Over 43.5")]]
    capture = load("fanduel_imply_ml_total.json")
    assert implybets_payload(legs) == capture["request"]
    quote = parse_parlay_quote(capture["response"], legs)
    assert quote.best.decimal_odds == pytest.approx(10.46339532) and quote.best.max_stake == 3509.0
    assert quote.leg_prices_american == {}


def test_parlay_quote_three_legs():
    response = load("fanduel_imply_ml_total_prop.json")["response"]
    quote = parse_parlay_quote(response, [object()] * 3)
    assert quote.best.decimal_odds == pytest.approx(27.41, abs=0.01) and quote.best.max_stake == 1111.5


def test_no_sgm_entry_is_a_decline():
    response = load("fanduel_imply_ml_spread_same_team.json")["response"]   # ML + spread, same team
    assert parse_parlay_quote(response, [object()] * 2) is None


# ---------------------------------------------------------------- the venue, fake transport

class FakeResponse:
    def __init__(self, status_code, payload):
        self.status_code, self._payload, self.text = status_code, payload, json.dumps(payload)[:300]

    def json(self):
        if isinstance(self._payload, Exception):
            raise self._payload
        return self._payload


class FakeSession:
    """Answers by URL (and tab); records every call."""

    def __init__(self, tab_status=None):
        self.calls = []
        self.tab_status = tab_status or {}

    def request(self, method, url, headers=None, timeout=None, params=None, data=None):
        self.calls.append((method, url, (params or {}).get("tab"), data))
        if url == API["scan_url"]:
            return FakeResponse(200, load("fanduel_scan.json"))
        if url == API["event_page_url"]:
            tab = params["tab"]
            status = self.tab_status.get(tab, 200)
            page = load("fanduel_event_ind_hou.json") if tab == PRIMARY_TAB else {"attachments": {"markets": {}}}
            return FakeResponse(status, page)
        if url == API["imply_url"]:
            return FakeResponse(200, load("fanduel_imply_ml_total.json")["response"])
        raise AssertionError(url)


def venue_with(session):
    return FanDuelVenue(session_factory=lambda: session, api=API)


def test_list_legs_reads_every_tab_and_caches(game):
    session = FakeSession()
    venue = venue_with(session)
    first = venue.list_legs(game)
    assert sorted(c[2] for c in session.calls) == sorted((PRIMARY_TAB, *EXTRA_TABS))
    assert venue.list_legs(game) == first and len(session.calls) == 1 + len(EXTRA_TABS)
    venue.list_legs(game, fresh=True)
    assert len(session.calls) == 2 * (1 + len(EXTRA_TABS))


def test_extra_tab_failure_keeps_the_rest(game):
    venue = venue_with(FakeSession(tab_status={"quarter-props": 503}))
    assert len(venue.list_legs(game)) > 20


def test_primary_tab_failure_raises(game):
    venue = venue_with(FakeSession(tab_status={PRIMARY_TAB: 403}))
    with pytest.raises(VenueTransportError, match="HTTP 200, got 403"):
        venue.list_legs(game)


def test_price_parlay_posts_implybets(game):
    session = FakeSession()
    venue = venue_with(session)
    legs = {(l.market_label, l.side_label): l for l in venue.list_legs(game)}
    chosen = [legs[("Moneyline", "Indianapolis Colts")], legs[("Total Points", "Over 43.5")]]
    quote = venue.price_parlay(game, chosen)
    method, url, _, body = session.calls[-1]
    assert (method, url) == ("POST", API["imply_url"]) and json.loads(body) == implybets_payload(chosen)
    assert quote.best.max_stake == 3509.0


def test_non_json_answer_is_a_transport_error():
    session = FakeSession()
    session.request = lambda *a, **k: FakeResponse(200, ValueError("not json"))
    with pytest.raises(VenueTransportError, match="not JSON"):
        venue_with(session).list_games()
