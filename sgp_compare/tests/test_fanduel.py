"""FanDuel adapter: parsing off trimmed recon captures (IND @ HOU, 2026-11-19)
and the venue's request handling off a fake session. No network."""
import json
from types import SimpleNamespace
from datetime import datetime, timezone
from pathlib import Path

import pytest

from sgp_compare.legs import CanonicalLeg
from sgp_compare.venues.base import VenueTransportError
from sgp_compare.venues.fanduel import (EXTRA_TABS, PRIMARY_TAB, REQUEST_WORKERS, FanDuelVenue, _load_fd_api,
                                        implybets_payload, markets_of_pages, parse_legs, parse_nfl_games,
                                        parse_parlay_quote)

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
    # The fixture's last three events are hand-added: a futures market (no " @ "),
    # an in-play game and a college game (other competitionId); all left out.
    assert len(games) == 3
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
    assert home.venue_ref == {"marketId": "734.196953815", "selectionId": 179823, "handicap": -12.5}
    assert home.is_main_line


def test_main_and_alternate_lines_share_a_market_key(by_label):
    """FanDuel posts alternates as separate markets; one key lets the page group
    them and stops two numbers of one bet going in one parlay."""
    main = by_label[("Spread", "Houston Texans -12.5")]
    alt = by_label[("Alternate Spread", "Houston Texans (-7.5)")]
    assert main.market_key == alt.market_key == "FG|spread|"
    assert by_label[("Nico Collins - Total Receptions", "Nico Collins Over 4.5")].market_key == \
        by_label[("Nico Collins - Alt Receptions", "Nico Collins 5+ Receptions")].market_key
    assert by_label[("2nd Quarter Spread", "Houston Texans -3.5")].market_key != main.market_key


def test_two_touchdown_scorers_can_share_a_parlay(by_label):
    montgomery = by_label[("Any Time Touchdown Scorer", "David Montgomery")]
    collins = by_label[("Any Time Touchdown Scorer", "Nico Collins")]
    assert montgomery.market_key != collins.market_key
    tie, home = by_label[("1st Quarter Winner (3-Way)", "Tie")], by_label[("1st Quarter Winner (3-Way)", "Houston Texans")]
    assert tie.market_key == home.market_key                  # one 3-way result only


def test_pickem_spread_keeps_its_number():
    game = parse_nfl_games(load("fanduel_scan.json"), BEFORE_KICKOFF)[0]
    page = load("fanduel_event_ind_hou.json")
    spread = json.loads(json.dumps(next(m for m in page["attachments"]["markets"].values() if m["marketName"] == "Spread")))
    for runner in spread["runners"]:
        runner["handicap"] = 0
    labels = {leg.side_label for leg in parse_legs(game, [spread])}
    assert labels == {"Houston Texans +0", "Indianapolis Colts +0"}


def test_page_without_attachments_is_an_error_not_an_empty_tab():
    with pytest.raises(ValueError, match="attachments"):
        markets_of_pages([{"error": "blocked"}])
    with pytest.raises(ValueError, match="attachments"):
        parse_nfl_games({"error": "blocked"}, BEFORE_KICKOFF)


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


def legs_sent_in(capture, moved=None):
    """Stand-in legs for a capture's request, at the handicaps FanDuel echoed
    (or `moved`: selectionId -> the handicap we think we sent)."""
    echoed = {(e["runner"]["marketId"], e["runner"]["selectionId"]): e["runner"].get("handicap", 0.0)
              for e in capture["response"]["runnerOdds"]}
    legs = []
    for bet_leg in capture["request"]["betLegs"]:
        runner = bet_leg["betRunners"][0]["runner"]
        key = (runner["marketId"], runner["selectionId"])
        handicap = (moved or {}).get(runner["selectionId"], echoed[key])
        legs.append(SimpleNamespace(market_label="m", side_label="s",
                                    venue_ref={**runner, "handicap": float(handicap)}))
    return legs


def test_parlay_quote_three_legs():
    capture = load("fanduel_imply_ml_total_prop.json")
    quote = parse_parlay_quote(capture["response"], legs_sent_in(capture))
    assert quote.best.decimal_odds == pytest.approx(27.41, abs=0.01) and quote.best.max_stake == 1111.5


def test_invalid_sgm_combination_is_a_decline():
    capture = load("fanduel_imply_ml_spread_same_team.json")   # ML + spread, same team
    assert parse_parlay_quote(capture["response"], legs_sent_in(capture)) is None


def test_no_sgm_entry_without_the_decline_code_is_an_error():
    """An expired PerimeterX token drops the SGP entry on HTTP 200: not a decline."""
    capture = load("fanduel_imply_ml_total.json")
    response = {**capture["response"],
                "betCombinations": [c for c in capture["response"]["betCombinations"] if not c["isSGM"]]}
    with pytest.raises(ValueError, match="no SGP entry"):
        parse_parlay_quote(response, legs_sent_in(capture))


def test_sgm_entry_without_a_price_is_an_error():
    capture = load("fanduel_imply_ml_total.json")
    response = json.loads(json.dumps(capture["response"]))
    sgm = next(c for c in response["betCombinations"] if c["isSGM"])
    sgm["winAvgOdds"]["trueOdds"]["decimalOdds"]["decimalOdds"] = None
    with pytest.raises(ValueError, match="decimal price"):
        parse_parlay_quote(response, legs_sent_in(capture))


def test_a_line_moved_while_pricing_raises():
    capture = load("fanduel_imply_ml_total.json")              # FanDuel echoes Over at 43.5
    with pytest.raises(ValueError, match="moved m s to 43.5"):
        parse_parlay_quote(capture["response"], legs_sent_in(capture, moved={7017916: 42.5}))


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
        self.headers = []
        self.tab_status = tab_status or {}

    def request(self, method, url, headers=None, timeout=None, params=None, data=None):
        self.calls.append((method, url, (params or {}).get("tab"), data))
        self.headers.append(headers)
        if url == API["scan_url"]:
            return FakeResponse(200, load("fanduel_scan.json"))
        if url == API["event_page_url"]:
            tab = params["tab"]
            status = self.tab_status.get(tab, 200)
            pages = {PRIMARY_TAB: "fanduel_event_ind_hou.json", "popular": "fanduel_event_ind_hou_other_tabs.json"}
            page = load(pages[tab]) if tab in pages else {"attachments": {"markets": {}}}
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


def test_any_tab_failure_fails_the_read_and_caches_nothing(game):
    """A missing tab would make its legs look pulled at pricing time."""
    session = FakeSession(tab_status={"quarter-props": 503})
    venue = venue_with(session)
    with pytest.raises(VenueTransportError, match="503"):
        venue.list_legs(game)
    session.tab_status = {}
    assert len(venue.list_legs(game)) > 40          # re-read, not a cached partial list


def test_extra_tab_markets_are_merged(game):
    legs = {(l.market_label, l.side_label): l for l in venue_with(FakeSession()).list_legs(game)}
    assert ("HOU Texans Alternate Total", "Over (9.5)") in legs and ("Spread", "Houston Texans -12.5") in legs


def test_sessions_are_bounded_and_reused(game):
    made = []
    session = FakeSession()
    venue = FanDuelVenue(session_factory=lambda: made.append(1) or session, api=API)
    for _ in range(3):
        venue.list_legs(game, fresh=True)
    venue.list_games()
    assert 1 <= len(made) <= REQUEST_WORKERS


def test_games_drop_out_at_kickoff_while_cached():
    now = [BEFORE_KICKOFF]
    session = FakeSession()
    venue = FanDuelVenue(session_factory=lambda: session, api=API, clock=lambda: 0.0, wall_clock=lambda: now[0])
    assert len(venue.list_games()) == 3
    now[0] = datetime(2026, 11, 20, 1, 17, tzinfo=timezone.utc)      # IND @ HOU kicked off a minute ago
    assert [g.home_team for g in venue.list_games()] == ["Philadelphia Eagles", "Chicago Bears"]
    assert len(session.calls) == 1                                    # served from the cache


def test_scan_request_names_nfl_and_sends_fd_headers():
    session = FakeSession()
    venue_with(session).list_games()
    method, url, _, body = session.calls[0]
    assert (method, url) == ("POST", API["scan_url"])
    assert json.loads(body)["filter"]["competitionIds"] == [12282733]
    assert session.headers[0]["x-application"] == "test" and session.headers[0]["Content-Type"] == "application/json"


def test_an_unreadable_answer_is_a_transport_error(game):
    session = FakeSession()
    original = session.request

    def broken(method, url, **kwargs):
        if url == API["event_page_url"]:
            return FakeResponse(200, {"attachments": {"markets": {"1": {"sgmMarket": True, "marketStatus": "OPEN"}}}})
        return original(method, url, **kwargs)

    session.request = broken
    with pytest.raises(VenueTransportError, match="event page: KeyError"):
        venue_with(session).list_legs(game)


def test_fanduel_config_comes_from_the_mlb_scraper():
    pytest.importorskip("curl_cffi")
    api = _load_fd_api()
    assert api["imply_url"].endswith("implyBets?pricePolicy=SUGGESTED") and "x-px-context" in api["headers"]


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


# ---------------------------------------------------------------- markets the SGP tab leaves out

@pytest.fixture
def other_tabs(game):
    """popular / 1st-half / 2nd-half / quarter tabs and a few SGP-tab extras,
    captured on Cal's Mac 2026-11-19 (IND @ HOU), runners trimmed."""
    legs = parse_legs(game, markets_of_pages([load("fanduel_event_ind_hou_other_tabs.json")]))
    return {(leg.market_label, leg.side_label): leg for leg in legs}


def test_team_alternate_totals(other_tabs):
    over = other_tabs[("HOU Texans Alternate Total", "Over (10.5)")]
    assert (over.canonical.period, over.canonical.stat, over.canonical.subject, over.canonical.line) == \
        ("FG", "team_total", "Houston Texans", 10.5)
    assert over.opposite_leg_id == other_tabs[("HOU Texans Alternate Total", "Under (10.5)")].leg_id
    half = other_tabs[("IND Colts - 1st Half Alt Total Points", "Over (0.5)")].canonical
    assert (half.period, half.subject) == ("1H", "Indianapolis Colts")
    assert other_tabs[("1st Half Alternate Total Points", "Over (42.5)")].canonical.period == "1H"


def test_two_way_period_winners_are_moneylines(other_tabs):
    home = other_tabs[("2nd Quarter Winner", "Houston Texans")]
    assert (home.canonical.period, home.canonical.stat, home.canonical.side) == ("2Q", "moneyline", "home")
    assert home.opposite_leg_id == other_tabs[("2nd Quarter Winner", "Indianapolis Colts")].leg_id
    assert other_tabs[("1st Half Winner", "Indianapolis Colts")].canonical.period == "1H"
    assert other_tabs[("2nd Half Spread", "Houston Texans -5.5")].canonical == \
        CanonicalLeg("2H", "spread", None, "home", -5.5)


def test_quarter_and_combo_player_props(other_tabs):
    assert other_tabs[("Kayshon Boutte - 1st Qtr Receiving Yds", "Kayshon Boutte Over 0.5")].canonical == \
        CanonicalLeg("1Q", "receiving_yards", "Kayshon Boutte", "over", 0.5)
    assert other_tabs[("C.J. Stroud - Passing + Rushing Yds", "C.J. Stroud Under 271.5")].canonical == \
        CanonicalLeg("FG", "passing_rushing_yards", "C.J. Stroud", "under", 271.5)
    assert other_tabs[("Nico Collins - Alt Receiving Yds", "Nico Collins 25+ Yards")].canonical.line == 24.5


def test_first_td_is_exclusive_either_player_is_not(other_tabs):
    first = [l for (m, _), l in other_tabs.items() if m == "1st Touchdown Scorer"]
    assert {l.canonical.stat for l in first} == {"to_score_first_touchdown"}
    assert len({l.market_key for l in first}) == 1
    either = [l for (m, _), l in other_tabs.items() if m == "Either Player - Anytime Touchdown Scorer"]
    assert len({l.market_key for l in either}) == len(either)


def test_non_sgp_market_on_another_tab_is_left_out(other_tabs):
    assert not any(m == "1st Half Total Touchdowns - HOU Texans" for m, _ in other_tabs)


# ---------------------------------------------------------------- filters, caches, failures (round 2)

def one_market(name, **changes):
    page = load("fanduel_event_ind_hou.json")
    market = json.loads(json.dumps(next(m for m in page["attachments"]["markets"].values() if m["marketName"] == name)))
    market.update(changes)
    return {"attachments": {"markets": {market["marketId"]: market}}}


@pytest.mark.parametrize("changes", [{"marketStatus": "SUSPENDED"}, {"inPlay": True}, {"sgmMarket": False}])
def test_closed_in_play_and_non_sgp_markets_are_left_out(changes):
    assert markets_of_pages([one_market("Moneyline", **changes)]) == []
    assert len(markets_of_pages([one_market("Moneyline")])) == 1


def test_runner_without_a_price_is_left_out(game):
    page = one_market("Moneyline")
    market = next(iter(page["attachments"]["markets"].values()))
    market["runners"][0]["winRunnerOdds"] = {}
    legs = parse_legs(game, markets_of_pages([page]))
    assert [leg.side_label for leg in legs] == ["Houston Texans"] and legs[0].opposite_leg_id is None


def test_games_are_sorted_by_kickoff():
    payload = load("fanduel_scan.json")
    payload["attachments"]["events"] = dict(reversed(list(payload["attachments"]["events"].items())))
    starts = [g.game_start_time for g in parse_nfl_games(payload, BEFORE_KICKOFF)]
    assert starts == sorted(starts) and len(starts) == 3


def test_caches_expire():
    clock = [0.0]
    session = FakeSession()
    venue = FanDuelVenue(session_factory=lambda: session, api=API, clock=lambda: clock[0],
                         wall_clock=lambda: BEFORE_KICKOFF)
    game = venue.list_games()[0]
    venue.list_legs(game)
    calls = len(session.calls)
    clock[0] = 19.0
    venue.list_legs(game)
    venue.list_games()
    assert len(session.calls) == calls                     # both still cached
    clock[0] = 21.0
    venue.list_legs(game)
    venue.list_games()
    assert len(session.calls) == calls + 1 + len(EXTRA_TABS)   # legs re-read, schedule still cached
    clock[0] = 301.0
    venue.list_games()
    assert session.calls[-1][1] == API["scan_url"]


def test_empty_schedule_is_not_cached():
    scans = []
    session = FakeSession()
    session.request = lambda method, url, **kw: scans.append(url) or FakeResponse(200, {"attachments": {"events": {}}})
    venue = FanDuelVenue(session_factory=lambda: session, api=API, clock=lambda: 0.0, wall_clock=lambda: BEFORE_KICKOFF)
    assert venue.list_games() == [] and venue.list_games() == []
    assert scans == [API["scan_url"]] * 2


def test_sgp_tab_with_no_markets_is_an_error(game):
    session = FakeSession()
    original = session.request
    session.request = lambda method, url, **kw: (FakeResponse(200, {"attachments": {}})
                                                 if kw.get("params", {}).get("tab") == PRIMARY_TAB
                                                 else original(method, url, **kw))
    with pytest.raises(VenueTransportError, match="listed no open SGP markets"):
        venue_with(session).list_legs(game)


def test_a_failed_request_drops_the_session_and_is_a_transport_error(game):
    made = []

    class Exploding(FakeSession):
        def request(self, *args, **kwargs):
            raise ConnectionError("reset by peer")

    venue = FanDuelVenue(session_factory=lambda: made.append(Exploding()) or made[-1], api=API)
    for _ in range(2):
        with pytest.raises(VenueTransportError, match="ConnectionError: reset by peer"):
            venue.price_parlay(game, [])
    assert len(made) == 2                                   # the broken session is not reused


def test_a_json_list_is_a_transport_error():
    session = FakeSession()
    session.request = lambda *a, **k: FakeResponse(200, [1, 2])
    with pytest.raises(VenueTransportError, match="expected a JSON object"):
        venue_with(session).list_games()


def test_sgm_entry_must_cover_every_leg():
    capture = load("fanduel_imply_ml_total_prop.json")       # 3-leg SGP entry
    with pytest.raises(ValueError, match="no SGP entry"):
        parse_parlay_quote(capture["response"], legs_sent_in(capture)[:2])


def test_a_main_line_moved_to_pickem_is_caught():
    """FanDuel omits the handicap key at 0, as on the moneyline here."""
    capture = load("fanduel_imply_ml_total.json")
    with pytest.raises(ValueError, match="to 0 while pricing"):
        parse_parlay_quote(capture["response"], legs_sent_in(capture, moved={50200: -1.5}))
