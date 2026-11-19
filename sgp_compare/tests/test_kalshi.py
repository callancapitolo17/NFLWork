"""Kalshi adapter, on trimmed copies of the 2026-11-19 recon (IND at HOU)."""
import json
from datetime import datetime, timezone
from pathlib import Path

import pytest

from sgp_compare.pricing import price_at_venue
from sgp_compare.venues import build_venues
from sgp_compare.venues.base import VenueTransportError
from sgp_compare.venues.kalshi import (
    COLLECTION_TICKER,
    KalshiVenue,
    _signed_api_from_env,
    all_in_price,
    combo_event_tickers,
    cost_to_american,
    parse_legs,
    parse_nfl_games,
    parse_quotes,
    selected_markets,
)

FIXTURES = Path(__file__).parent / "fixtures"
BEFORE_KICKOFF = datetime(2026, 11, 19, tzinfo=timezone.utc)
RFQ_ID = "7e01605e-2231-41de-bfc7-b821429eae79"
COMBO = "KXMVECROSSCATEGORY-S2026LBCLME447PBSCVQX-NCLPRHGBP"


def load(name):
    return json.loads((FIXTURES / name).read_text())


@pytest.fixture(autouse=True)
def before_kickoff(monkeypatch):
    """The venue's wall clock pinned before the fixture game, so the suite does
    not go red once IND at HOU kicks off for real."""
    import sgp_compare.venues.kalshi as kalshi_module

    class Pinned(datetime):
        @classmethod
        def now(cls, tz=None):
            return BEFORE_KICKOFF
    monkeypatch.setattr(kalshi_module, "datetime", Pinned)


@pytest.fixture
def game():
    return parse_nfl_games(load("kalshi_games.json"), BEFORE_KICKOFF)[0]


@pytest.fixture
def legs(game):
    size_max = combo_event_tickers(load("kalshi_collection.json"), game)
    teams = {"IND": "Indianapolis", "HOU": "Houston"}
    return {(leg.market_label, leg.side_label): leg
            for leg in parse_legs(game, teams, load("kalshi_events_ind_hou.json"), size_max)}


def test_games_read_away_then_home_and_start_48h_before_close(game):
    assert (game.venue_game_id, game.away_team, game.home_team) == ("26NOV19INDHOU", "Indianapolis", "Houston")
    assert game.game_start_time == datetime(2026, 11, 20, 1, 15, tzinfo=timezone.utc)


def test_started_games_are_left_out():
    assert parse_nfl_games(load("kalshi_games.json"), datetime(2026, 12, 31, tzinfo=timezone.utc)) == []


def test_fee_is_in_every_price():
    assert all_in_price(0.46) == pytest.approx(0.46 + 0.07 * 0.46 * 0.54)
    assert cost_to_american(0.4) == pytest.approx(150.0)
    assert cost_to_american(0.6) == pytest.approx(-150.0)


def test_spread_yes_is_the_named_team_laying_the_number(legs):
    leg = legs[("Spread", "Houston -13.5")]
    other = legs[("Spread", "Indianapolis +13.5")]
    assert (leg.canonical.period, leg.canonical.stat, leg.canonical.side, leg.canonical.line) == ("FG", "spread", "home", -13.5)
    assert (other.canonical.side, other.canonical.line) == ("away", 13.5)
    assert leg.opposite_leg_id == other.leg_id == "KXNFLSPREAD-26NOV19INDHOU-HOU14|no"
    assert leg.american_odds == pytest.approx(cost_to_american(all_in_price(0.46)))
    assert leg.opposite_american_odds == other.american_odds
    assert leg.is_main_line and leg.market_key == "KXNFLSPREAD-26NOV19INDHOU"
    assert selected_markets([leg, other]) == [
        {"market_ticker": "KXNFLSPREAD-26NOV19INDHOU-HOU14", "event_ticker": "KXNFLSPREAD-26NOV19INDHOU", "side": "yes"},
        {"market_ticker": "KXNFLSPREAD-26NOV19INDHOU-HOU14", "event_ticker": "KXNFLSPREAD-26NOV19INDHOU", "side": "no"}]


def test_one_main_line_per_spread_ladder(legs):
    mains = [key for key, leg in legs.items() if key[0] == "Spread" and leg.is_main_line]
    assert sorted(mains) == [("Spread", "Houston -13.5"), ("Spread", "Indianapolis +13.5")]


def test_moneyline_no_is_the_other_team(legs):
    no_on_houston = legs[("Moneyline", "No on Houston (= Indianapolis)")]
    assert (no_on_houston.canonical.stat, no_on_houston.canonical.side) == ("moneyline", "away")
    assert legs[("Moneyline", "Houston")].canonical.side == "home"


def test_half_winner_is_three_way_with_a_tie_market(legs):
    tie = legs[("1H Winner (tie loses)", "Tie")]
    assert (tie.canonical.period, tie.canonical.stat, tie.canonical.subject, tie.canonical.side) == ("1H", "win_3way", "Tie", "yes")
    assert legs[("1H Winner (tie loses)", "Not Houston")].canonical.side == "no"


def test_team_totals_name_the_team(legs):
    leg = legs[("Houston Team Total", "Over 28.5")]
    assert (leg.canonical.stat, leg.canonical.subject, leg.canonical.line) == ("team_total", "Houston", 28.5)


def test_player_props_protected_is_its_own_bet(legs):
    plain = legs[("Nico Collins Receiving Yards", "Over 69.5")]
    protected = legs[("Nico Collins Receiving Yards (Protected)", "Over 69.5")]
    assert (plain.canonical.stat, plain.canonical.subject, plain.category) == ("receiving_yards", "Nico Collins", "Receiving")
    assert (protected.canonical.stat, protected.category) == ("receiving_yards_protected", "Receiving (Protected)")
    assert plain.is_main_line
    # Two players' yards can share a combo (size_max absent); two of one player's numbers cannot.
    assert plain.market_key == "KXNFLRECYDS-26NOV19INDHOU|32d2bbf9-ecb2-44b0-a8ef-f39f6310d36c"


def test_touchdowns_one_plus_is_anytime_td(legs):
    anytime = legs[("Nico Collins Anytime TD", "YES")]
    assert (anytime.canonical.stat, anytime.canonical.side, anytime.canonical.line) == ("anytime_td", "yes", None)
    two = legs[("Nico Collins Touchdowns", "Over 1.5")]
    assert (two.canonical.stat, two.canonical.line) == ("touchdowns", 1.5)


def test_winning_margin_buckets_have_no_canonical(legs):
    assert legs[("Winning Margin", "Houston by 15 or more")].canonical is None


def test_quotes_buy_yes_at_one_minus_no_bid_best_first():
    tiers = parse_quotes(load("kalshi_quotes_spread_over.json")["quotes"], RFQ_ID, COMBO)
    assert tiers[0].decimal_odds == pytest.approx(1.0 / all_in_price(1.0 - 0.759), rel=1e-5)
    assert tiers[0].max_stake == pytest.approx(393.85 * all_in_price(0.241), abs=0.01)
    assert [t.decimal_odds for t in tiers] == sorted((t.decimal_odds for t in tiers), reverse=True)
    assert all(t.decimal_odds >= 1.01 for t in tiers)            # the 0.0001 no_bid is gone
    assert parse_quotes(load("kalshi_quotes_spread_over.json")["quotes"], "another-rfq", COMBO) == ()


# ---------------------------------------------------------------- the RFQ flow

class FakeKalshi:
    """Answers the adapter's signed calls from the fixtures and records them."""

    def __init__(self, mint_status=200, quotes_status=200, quotes=None, rfq_status=201, rfq_raises=False):
        self.calls = []
        self.rfq_status = rfq_status
        self.rfq_raises = rfq_raises
        self.mint_status = mint_status
        self.quotes_status = quotes_status
        self.quotes = load("kalshi_quotes_spread_over.json")["quotes"] if quotes is None else quotes
        self.rfqs = 0

    def __call__(self, method, path, body=None):
        self.calls.append((method, path, body))
        if method == "GET" and path.startswith("/events?series_ticker=KXNFLGAME"):
            return 200, load("kalshi_games.json"), {}
        if method == "GET" and path == f"/multivariate_event_collections/{COLLECTION_TICKER}":
            return 200, load("kalshi_collection.json"), {}
        if method == "GET" and path.startswith("/events/"):
            ticker = path.split("/")[2].split("?")[0]
            match = [e for e in load("kalshi_events_ind_hou.json") if e["event"]["event_ticker"] == ticker]
            return (200, match[0], {}) if match else (404, {"error": {"code": "not_found"}}, {})
        if method == "POST" and path == f"/multivariate_event_collections/{COLLECTION_TICKER}":
            return self.mint_status, ({"market_ticker": COMBO, "event_ticker": "X"} if self.mint_status == 200
                                      else {"error": {"code": "invalid"}}), {}
        if method == "POST" and path == "/communications/rfqs":
            self.rfqs += 1
            if self.rfq_raises:
                raise TimeoutError("read timed out")
            return self.rfq_status, ({"id": RFQ_ID} if self.rfq_status == 201 else {"error": {"code": "x"}}), {}
        if method == "GET" and path.startswith("/communications/rfqs?status=open"):
            return 200, {"rfqs": [{"id": "leaked-1", "market_ticker": COMBO},
                                  {"id": "other-1", "market_ticker": "KXMVECROSSCATEGORY-S-OTHER"}]}, {}
        if method == "GET" and path.startswith("/communications/quotes?"):
            return self.quotes_status, {"quotes": self.quotes}, {}
        if method == "DELETE":
            return 204, {}, {}
        raise AssertionError(f"unexpected Kalshi call {method} {path}")

    def deletes(self):
        return [call for call in self.calls if call[0] == "DELETE"]


class FakeClock:
    def __init__(self):
        self.now = 0.0

    def __call__(self):
        return self.now

    def sleep(self, seconds):
        self.now += seconds


def venue_on(fake):
    clock = FakeClock()
    return KalshiVenue(api_factory=lambda: (fake, "user-1"), clock=clock, sleep=clock.sleep)


def two_legs(venue):
    game = venue.list_games()[0]
    offered = {(leg.market_label, leg.side_label): leg for leg in venue.list_legs(game)}
    return game, [offered[("Spread", "Houston -13.5")], offered[("Nico Collins Receiving Yards", "Over 69.5")]]


def test_price_parlay_asks_reads_and_always_cancels():
    fake = FakeKalshi()
    venue = venue_on(fake)
    game, legs = two_legs(venue)
    quote = venue.price_parlay(game, legs)
    assert quote.best.decimal_odds == pytest.approx(1.0 / all_in_price(0.241), rel=1e-5)
    assert quote.raw["combo_market_ticker"] == COMBO and quote.raw["delete_status"] == 204
    assert "creator_id" not in json.dumps(quote.raw)
    mint = next(c for c in fake.calls if c[0] == "POST" and "multivariate" in c[1])
    assert mint[2] == {"selected_markets": selected_markets(legs)}
    rfq = next(c for c in fake.calls if c[1] == "/communications/rfqs")
    assert rfq[2]["target_cost_dollars"] == "100.00" and rfq[2]["rest_remainder"] is False
    assert fake.deletes() == [("DELETE", f"/communications/rfqs/{RFQ_ID}", None)]
    assert not any("accept" in path for _, path, _ in fake.calls)
    # Quotes were steady from the first poll, so it stopped at the minimum window.
    assert quote.raw["seconds_waited"] == pytest.approx(1.5)


def test_no_quotes_waits_the_full_window_says_so_and_still_cancels():
    fake = FakeKalshi(quotes=[])
    venue = venue_on(fake)
    game, legs = two_legs(venue)
    with pytest.raises(VenueTransportError, match="no market maker quoted within 4.0 s"):
        venue.price_parlay(game, legs)
    assert len([c for c in fake.calls if c[1].startswith("/communications/quotes?")]) == 8
    assert len(fake.deletes()) == 1


def test_an_unclear_rfq_create_reply_cancels_any_rfq_it_made():
    fake = FakeKalshi(rfq_status=503)
    venue = venue_on(fake)
    game, legs = two_legs(venue)
    with pytest.raises(VenueTransportError, match="create RFQ: expected 201"):
        venue.price_parlay(game, legs)
    # The open-RFQ list holds one on our combo and one on another; only ours goes.
    assert fake.deletes() == [("DELETE", "/communications/rfqs/leaked-1", None)]


def test_a_timed_out_rfq_create_cancels_any_rfq_it_made():
    fake = FakeKalshi(rfq_raises=True)
    venue = venue_on(fake)
    game, legs = two_legs(venue)
    with pytest.raises(VenueTransportError, match="TimeoutError"):
        venue.price_parlay(game, legs)
    assert fake.deletes() == [("DELETE", "/communications/rfqs/leaked-1", None)]


def test_a_refused_rfq_create_is_an_error_and_cancels_nothing():
    fake = FakeKalshi(rfq_status=409)
    venue = venue_on(fake)
    game, legs = two_legs(venue)
    with pytest.raises(VenueTransportError, match="409"):
        venue.price_parlay(game, legs)
    assert fake.deletes() == []   # a 4xx made no RFQ; another click's live one on this combo stays


def test_quotes_with_no_usable_price_say_so():
    junk = [{"rfq_id": RFQ_ID, "market_ticker": COMBO, "status": "open", "no_bid_dollars": "0.0001",
             "no_contracts_fp": "100"}]
    venue = venue_on(FakeKalshi(quotes=junk))
    game, legs = two_legs(venue)
    with pytest.raises(VenueTransportError, match="1 quote\\(s\\) in .* none a usable price"):
        venue.price_parlay(game, legs)


def test_the_game_list_follows_the_cursor():
    fake = FakeKalshi()
    pages = iter([{"events": load("kalshi_games.json")["events"][:1], "cursor": "next"},
                  {"events": load("kalshi_games.json")["events"][1:]}])
    original = fake.__call__

    def paged(method, path, body=None):
        if method == "GET" and path.startswith("/events?series_ticker=KXNFLGAME"):
            fake.calls.append((method, path, body))
            return 200, next(pages), {}
        return original(method, path, body)
    clock = FakeClock()
    venue = KalshiVenue(api_factory=lambda: (paged, "user-1"), clock=clock, sleep=clock.sleep)
    assert len(venue.list_games()) == 2
    game_reads = [path for _, path, _ in fake.calls if path.startswith("/events?series_ticker")]
    assert "cursor" not in game_reads[0] and game_reads[1].endswith("&cursor=next")


def test_a_failed_quote_read_still_cancels():
    fake = FakeKalshi(quotes_status=500)
    venue = venue_on(fake)
    game, legs = two_legs(venue)
    with pytest.raises(VenueTransportError):
        venue.price_parlay(game, legs)
    assert len(fake.deletes()) == 1


def test_a_refused_combo_sends_no_rfq():
    fake = FakeKalshi(mint_status=400)
    venue = venue_on(fake)
    game, legs = two_legs(venue)
    assert venue.price_parlay(game, legs) is None
    assert fake.rfqs == 0


def test_price_at_venue_prices_every_combination_one_rfq_each():
    fake = FakeKalshi()
    venue = venue_on(fake)
    game, legs = two_legs(venue)
    result = price_at_venue(venue, game, [leg.leg_id for leg in legs])
    assert result.status == "priced"
    assert fake.rfqs == 4 and len(fake.deletes()) == 4
    assert result.partition_decimal_odds is not None and len(result.partition_decimal_odds) == 4
    mint_sides = [[leg["side"] for leg in body["selected_markets"]] for method, path, body in fake.calls
                  if method == "POST" and path.startswith("/multivariate_event_collections/")]
    # Each combination flips a leg to the NO of its own market (partition_flips order).
    assert mint_sides == [["yes", "yes"], ["yes", "no"], ["no", "yes"], ["no", "no"]]


def test_missing_credentials_name_what_to_set(tmp_path, monkeypatch):
    for name in ("KALSHI_API_KEY_ID", "KALSHI_PRIVATE_KEY_PATH", "KALSHI_USER_ID"):
        monkeypatch.delenv(name, raising=False)
    with pytest.raises(VenueTransportError, match="KALSHI_API_KEY_ID, KALSHI_PRIVATE_KEY_PATH, KALSHI_USER_ID"):
        _signed_api_from_env(tmp_path / "missing.env")


def test_team_codes_that_prefix_each_other_still_pair():
    event = {"event_ticker": "KXNFLGAME-26NOV22LACLA", "markets": [
        {"ticker": "KXNFLGAME-26NOV22LACLA-LA", "yes_sub_title": "Los Angeles R", "close_time": "2026-11-24T21:25:00Z"},
        {"ticker": "KXNFLGAME-26NOV22LACLA-LAC", "yes_sub_title": "Los Angeles C", "close_time": "2026-11-24T21:25:00Z"}]}
    games = parse_nfl_games({"events": [event]}, BEFORE_KICKOFF)
    assert [(g.away_team, g.home_team) for g in games] == [("Los Angeles C", "Los Angeles R")]


def test_a_game_drops_out_of_the_cached_list_at_kickoff(monkeypatch):
    import sgp_compare.venues.kalshi as kalshi_module
    venue = venue_on(FakeKalshi())
    assert venue.list_games()[0].venue_game_id == "26NOV19INDHOU"

    class AfterKickoff(datetime):
        @classmethod
        def now(cls, tz=None):
            return datetime(2026, 11, 20, 1, 16, tzinfo=timezone.utc)
    monkeypatch.setattr(kalshi_module, "datetime", AfterKickoff)
    assert "26NOV19INDHOU" not in [g.venue_game_id for g in venue.list_games()]


def test_same_price_quotes_take_one_quotes_size_not_the_sum():
    base = {"rfq_id": RFQ_ID, "market_ticker": COMBO, "status": "open", "no_bid_dollars": "0.7220"}
    quotes = [{**base, "no_contracts_fp": "342.34"}, {**base, "no_contracts_fp": "342.34"},
              {**base, "status": "cancelled", "no_bid_dollars": "0.9000", "no_contracts_fp": "1000"},
              {**base, "market_ticker": "OTHER", "no_bid_dollars": "0.9000", "no_contracts_fp": "1000"},
              {**base, "no_bid_dollars": None, "no_contracts_fp": "10"}]
    (tier,) = parse_quotes(quotes, RFQ_ID, COMBO)
    assert tier.max_stake == pytest.approx(342.34 * all_in_price(0.278), abs=0.01)


def test_kalshi_is_a_venue():
    assert "kalshi" in build_venues()
