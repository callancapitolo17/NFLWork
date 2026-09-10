"""BetMGM legacy ``games`` main lines -> structure (leg surface, 2026-09-03).

BetMGM builds a next-day fixture's full ``optionMarkets`` tree only the
next morning (~05:50 PT). Until then the fixture carries its 4 main lines
(money line, two run lines, main total) ONLY in the legacy ``games`` array,
so ``fetch_markets`` returned [] and the leg surface priced 1 of 7 games
all night. The fallback is opt-in (``include_main_line_games``) because the
bet-builder refuses those ids — the same-game on-demand path must keep
seeing an empty structure and decline cleanly.
"""
import sys
from pathlib import Path
from unittest.mock import MagicMock

sys.path.insert(0, str(Path(__file__).parent.parent))

from mlb_sgp import betmgm  # noqa: E402
from mlb_sgp.betmgm_client import (  # noqa: E402
    BetMGMClient, _main_line_game_as_option_market, _merge_main_line_games,
    _split_home_away)

HOME = "Cincinnati Reds"
AWAY = "Milwaukee Brewers"


def _legacy_games():
    """Trimmed from the live 2026-09-03 payload for fixture 19888929."""
    return [
        {"id": 1551940707, "name": {"value": "Money Line"},
         "results": [
             {"id": 2278464933, "odds": 1.62, "name": {"value": "Brewers"}},
             {"id": 2278464934, "odds": 2.35, "name": {"value": "Reds"}}]},
        {"id": 1551940710, "name": {"value": "Run Line Spread"},
         "results": [
             {"id": 2278464939, "odds": 1.98, "attr": "-1,5",
              "name": {"value": "Milwaukee Brewers -1,5"}},
             {"id": 2278464940, "odds": 1.85, "attr": "+1,5",
              "name": {"value": "Cincinnati Reds +1,5"}}]},
        {"id": 1551940711, "name": {"value": "Run Line Spread"},
         "results": [
             {"id": 2278464941, "odds": 2.5, "attr": "-2,5",
              "name": {"value": "Milwaukee Brewers -2,5"}},
             {"id": 2278464942, "odds": 1.5, "attr": "+2,5",
              "name": {"value": "Cincinnati Reds +2,5"}}]},
        {"id": 1551940724, "name": {"value": "Totals"}, "attr": "9,5",
         "results": [
             {"id": 2278464967, "odds": 1.91, "totalsPrefix": "Over",
              "name": {"value": "Over 9,5"}},
             {"id": 2278464968, "odds": 1.91, "totalsPrefix": "Under",
              "name": {"value": "Under 9,5"}}]},
    ]


def _client_returning(fixture: dict) -> BetMGMClient:
    client = BetMGMClient()
    client.accessid = MagicMock(return_value="static-id")
    response = MagicMock(status_code=200)
    response.json.return_value = {"fixtures": [fixture]}
    client._get_fixture_markets = MagicMock(return_value=response)
    return client


def test_legacy_game_market_reshapes_to_option_market_schema():
    om = _main_line_game_as_option_market(_legacy_games()[1])
    assert om["id"] == 1551940710
    assert om["name"] == {"value": "Run Line Spread"}
    assert om["options"][0] == {
        "id": 2278464939, "name": {"value": "Milwaukee Brewers -1,5"},
        "attr": "-1,5", "totalsPrefix": None, "price": {"odds": 1.98}}


def test_merge_keeps_option_market_when_both_carry_the_same_id():
    real = {"id": 1551940707, "name": {"value": "Money Line"},
            "options": [{"id": 1, "price": {"odds": 9.9}}]}
    merged = _merge_main_line_games([real], _legacy_games())
    # The optionMarkets entry wins; the other 3 legacy markets are appended.
    assert merged[0] is real
    assert [m["id"] for m in merged] == [
        1551940707, 1551940710, 1551940711, 1551940724]


def test_parse_markets_prices_the_four_main_lines_from_legacy_games():
    markets = _merge_main_line_games([], _legacy_games())
    parsed = betmgm.parse_markets(markets, HOME, AWAY)
    fg = parsed["FG"]
    assert fg["moneyline"]["home"] == (1551940707, 2278464934, 2.35)
    assert fg["moneyline"]["away"] == (1551940707, 2278464933, 1.62)
    # Home-perspective signed lines: Brewers (away) -1.5 => home +1.5.
    assert set(fg["spreads"]) == {1.5, 2.5}
    assert fg["spreads"][1.5]["home"] == (1551940710, 2278464940, 1.85)
    assert fg["spreads"][1.5]["away"] == (1551940710, 2278464939, 1.98)
    assert fg["totals"][9.5]["over"] == (1551940724, 2278464967, 1.91)
    assert fg["totals"][9.5]["under"] == (1551940724, 2278464968, 1.91)
    assert parsed["F5"]["totals"] == {} and parsed["I1"]["totals"] == {}


def test_fetch_markets_default_ignores_legacy_games():
    client = _client_returning({"optionMarkets": [], "games": _legacy_games()})
    assert client.fetch_markets("19888929") == []


def test_fetch_markets_opt_in_reads_legacy_games():
    client = _client_returning({"optionMarkets": [], "games": _legacy_games()})
    markets = client.fetch_markets("19888929", include_main_line_games=True)
    assert [m["id"] for m in markets] == [
        1551940707, 1551940710, 1551940711, 1551940724]


def test_fetch_markets_opt_in_is_a_no_op_once_the_tree_is_built():
    tree = [{"id": 7, "name": {"value": "Totals"}, "options": []}]
    client = _client_returning({"optionMarkets": tree, "games": []})
    assert client.fetch_markets("x", include_main_line_games=True) == tree


def test_doubleheader_game_number_is_stripped_from_the_home_name():
    # Live 2026-09-03: the marker rides on the home side of the fixture name
    # and left both games of every doubleheader unresolvable at BetMGM.
    fixture = {"name": {"value": "Detroit Tigers at Cleveland Guardians (Game 1)"}}
    assert _split_home_away(fixture) == ("Cleveland Guardians", "Detroit Tigers")
    plain = {"name": {"value": "Detroit Tigers at Cleveland Guardians"}}
    assert _split_home_away(plain) == ("Cleveland Guardians", "Detroit Tigers")
