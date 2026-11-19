import json
from datetime import datetime, timezone
from pathlib import Path

import pytest

from sgp_compare.venues.prophetx import (parlay_payload, parse_legs, parse_nfl_games, parse_parlay_quote)

FIXTURES = Path(__file__).parent / "fixtures"
HOUSTON_ID = 1700004938
BEFORE_KICKOFF = datetime(2026, 11, 19, tzinfo=timezone.utc)


def load(name):
    return json.loads((FIXTURES / name).read_text())


@pytest.fixture
def game():
    return parse_nfl_games(load("prophetx_tournaments.json"), BEFORE_KICKOFF)[0]


@pytest.fixture
def legs(game):
    return {(leg.market_label, leg.side_label): leg
            for leg in parse_legs(game, load("prophetx_markets_ind_hou.json"), HOUSTON_ID)}


def test_games_are_nfl_only_home_is_seq_zero(game):
    games = parse_nfl_games(load("prophetx_tournaments.json"), BEFORE_KICKOFF)
    assert all("Over/Under" not in g.home_team and "Wins" not in g.label for g in games)
    assert (game.home_team, game.away_team) == ("Houston Texans", "Indianapolis Colts")
    assert game.game_start_time == datetime(2026, 11, 20, 1, 15, tzinfo=timezone.utc)


def test_started_games_are_left_out():
    after = datetime(2026, 12, 31, tzinfo=timezone.utc)
    assert parse_nfl_games(load("prophetx_tournaments.json"), after) == []


def test_spread_leg_prices_sides_and_canonical(legs):
    leg = legs[("Spread", "Houston Texans -12.5")]
    assert (leg.american_odds, leg.opposite_american_odds) == (100.0, -105.0)
    assert leg.opposite_leg_id == legs[("Spread", "Indianapolis Colts +12.5")].leg_id
    assert (leg.canonical.period, leg.canonical.stat, leg.canonical.side, leg.canonical.line) == ("FG", "spread", "home", -12.5)
    assert leg.venue_ref == {"sportEventId": 1700008817, "marketId": 223, "outcomeId": 1714,
                             "lineId": "7342bc136571bb4c9584fd65c027c042", "line": -12.5}


def test_moneyline_sides_come_from_competitor_ids(legs):
    home = legs[("Moneyline", "Houston Texans")]
    assert home.canonical.side == "home" and home.american_odds == -570.0


def test_player_props_name_the_player_and_stat(legs):
    collins = legs[("Nico Collins Receiving Yards", "Over 70.5")]
    assert (collins.canonical.subject, collins.canonical.stat, collins.canonical.line) == ("Nico Collins", "receiving_yards", 70.5)
    assert collins.is_main_line
    td = legs[("David Montgomery To Score a Touchdown", "YES")]
    assert (td.canonical.subject, td.canonical.side, td.canonical.line) == ("David Montgomery", "yes", None)


def test_periods_and_team_totals(legs):
    assert legs[("First Half Spread", "Houston Texans -6.5")].canonical.period == "1H"
    team_total = legs[("Indianapolis Colts: Team Total Points", "Over 15.5")].canonical
    assert (team_total.stat, team_total.subject) == ("team_total", "Indianapolis Colts")
    draw = legs[("Draw (2nd Quarter)", "YES")].canonical
    assert (draw.period, draw.stat, draw.subject) == ("2Q", "win_3way", "Draw")


def test_lines_nobody_offers_are_left_out(legs):
    assert ("Nico Collins Receiving Yards", "Over 14.5") not in legs


def test_parlay_quote_parses_tiers_best_first_and_shaded_leg_prices(legs):
    recorded = load("prophetx_parlay_spread_over.json")
    chosen = [legs[("Spread", "Houston Texans -12.5")], legs[("Total Points", "Over 43.5")]]
    assert parlay_payload(chosen) == recorded["request"]
    quote = parse_parlay_quote(recorded["response"], chosen)
    assert [round(t.decimal_odds, 4) for t in quote.tiers] == [2.7165, 2.6758]
    assert [t.max_stake for t in quote.tiers] == [62.62, 125.24]
    assert quote.leg_prices_american == {chosen[0].leg_id: -148.0, chosen[1].leg_id: -161.0}


def test_a_declined_parlay_is_none(legs):
    chosen = [legs[("Spread", "Houston Texans -12.5")], legs[("Total Points", "Over 43.5")]]
    assert parse_parlay_quote({"success": False, "data": None}, chosen) is None
    assert parse_parlay_quote({"success": True, "data": {"offers": []}}, chosen) is None
