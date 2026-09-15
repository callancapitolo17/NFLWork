"""normalize_kalshi on the recorded fixture: the same semantics bets.test.js pins."""
import pytest

from unabated_ticket.bets_service.normalize import cents_to_american, js_round
from unabated_ticket.bets_service.sources import kalshi_ticker as tk
from unabated_ticket.bets_service.sources.kalshi import aggregate_fills, normalize_kalshi
from unabated_ticket.bets_service.tests.conftest import (
    FETCHED_AT, event, fill, market, position)


def normalized(fx: dict, extra_fills=(), extra_positions=(), extra_markets=None, extra_events=None):
    return normalize_kalshi(
        fx["fills"] + list(extra_fills), fx["positions"] + list(extra_positions),
        {**fx["markets"], **(extra_markets or {})}, {**fx["events"], **(extra_events or {})},
        FETCHED_AT)


def by_id(records: list[dict], record_id: str) -> dict:
    matches = [record for record in records if record["id"] == record_id]
    assert matches, f"record {record_id} missing"
    return matches[0]


# ---- helpers -----------------------------------------------------------------------

def test_js_round_halves_go_up_like_math_round():
    assert js_round(2.5) == 3
    assert js_round(0.5) == 1
    assert round(2.5) == 2  # Python's banker's rounding is exactly why js_round exists


@pytest.mark.parametrize("cents,american", [(32, 213), (64, -178), (50, -100), (0, None), (100, None)])
def test_cents_to_american_both_sides_of_even_money(cents, american):
    assert cents_to_american(cents) == american


def test_parse_event_suffix_per_sport():
    assert tk.parse_event_suffix("KXNCAAFSPREAD-26SEP12CHATEKY") == {
        "eventDate": "2026-09-12", "eventStart": None, "gameNumber": None}
    assert tk.parse_event_suffix("KXMLBRFI-26SEP062210WSHLAD") == {
        "eventDate": "2026-09-06", "eventStart": "2026-09-07T02:10:00.000Z", "gameNumber": None}
    assert tk.parse_event_suffix("KXMLBGAME-26SEP041410DETCLEG2")["gameNumber"] == 2
    assert tk.parse_event_suffix("KXMLBGAME-26JAN151300DETCLE")["eventStart"] == "2026-01-15T18:00:00.000Z"
    assert tk.parse_event_suffix("KXNFLOROTY-27") is None


# ---- fills -> positions aggregation ------------------------------------------------

def test_aggregate_fills_vwap_over_buys_and_netting():
    fills = [
        fill("T", "no", 45, 0.50, "2026-09-05T14:55:08Z"),
        fill("T", "no", 200, 0.51, "2026-09-05T14:56:41Z"),
        fill("T", "no", 20, 0.60, "2026-09-05T15:00:00Z", action="sell"),
    ]
    agg = aggregate_fills(fills)
    assert agg["vwapCents"] == pytest.approx((45 * 50 + 200 * 49) / 245)
    assert agg["netContracts"] == 225
    assert agg["firstFillAt"] == "2026-09-05T14:55:08Z"
    assert agg["fillCount"] == 3


def test_aggregate_sell_only_prices_from_the_sells():
    fills = [fill("T", "yes", 100, 0.74, "2026-08-19T15:27:14Z", action="sell")]
    agg = aggregate_fills(fills)
    assert agg["vwapCents"] == pytest.approx(74)
    assert agg["netContracts"] == -100


def test_yes_spread_stake_from_position_and_vwap(kalshi_fixture):
    record = by_id(normalized(kalshi_fixture), "kalshi:KXNCAAFSPREAD-26SEP12CHATEKY-CHAT6:yes")
    assert record["source"] == "kalshi_api"
    assert record["venue"] == "kalshi"
    assert record["league"] == "cfb"
    assert (record["betType"], record["period"]) == ("spread", "FG")
    assert (record["awayTeam"], record["homeTeam"]) == ("Chattanooga", "Eastern Kentucky")
    assert (record["awayKey"], record["homeKey"]) == (None, None)  # the extension resolves keys
    assert (record["side"], record["points"]) == ("away", -5.5)
    assert record["price"] == 138
    assert record["contracts"] == 400  # position_fp, not the fills
    assert record["stake"] == 168
    assert record["toWin"] == 232
    assert record["placedAt"] == "2026-09-11T17:14:55.910367Z"
    assert record["status"] == "open"
    assert record["closedAt"] is None
    assert record["eventStart"] is None  # football suffix carries the date only
    assert record["eventDate"] == "2026-09-12"
    assert record["approx"] == []
    assert record["unmatchable"] is None
    assert record["sourceFetchedAt"] == FETCHED_AT
    assert record["raw"]["ticker"] == "KXNCAAFSPREAD-26SEP12CHATEKY-CHAT6"
    assert isinstance(record["contracts"], int) and isinstance(record["stake"], int)


def test_position_zero_with_trades_is_closed_and_sells_price_from_sells(kalshi_fixture):
    record = by_id(normalized(kalshi_fixture), "kalshi:KXNEXTTEAMNFL-26MCROSBY-LV:yes")
    assert record["status"] == "closed"
    assert record["contracts"] == 0
    assert record["closedAt"] == "2026-08-24T11:50:26.544087Z"
    assert record["unmatchable"] == "not a game market"
    assert record["price"] == -285


def test_position_on_the_other_side_closes_this_sides_fills(kalshi_fixture):
    fills = [
        fill("KXNCAAFTOTAL-26SEP12RICEND-60", "no", 60, 0.50, "2026-09-11T15:09:00Z"),
    ]
    records = normalized(kalshi_fixture, extra_fills=fills)  # position_fp +500 = YES only
    assert by_id(records, "kalshi:KXNCAAFTOTAL-26SEP12RICEND-60:no")["status"] == "closed"
    assert by_id(records, "kalshi:KXNCAAFTOTAL-26SEP12RICEND-60:yes")["status"] == "open"


def test_settled_market_without_a_position_reads_the_result(kalshi_fixture):
    record = by_id(normalized(kalshi_fixture), "kalshi:KXNCAAFSPREAD-26SEP05MOSUTXAM-TXAM39:yes")
    assert (record["awayTeam"], record["homeTeam"]) == ("Missouri St.", "Texas A&M")
    assert (record["side"], record["points"]) == ("home", -38.5)
    assert record["price"] == -108
    assert record["stake"] == 196.56
    assert record["status"] == "won"
    assert record["closedAt"] == "2026-09-06T02:00:00Z"  # expected_expiration_time


# ---- both spread signs, both total directions ---------------------------------------

def test_no_spread_is_the_other_team_at_plus_floor_strike(kalshi_fixture):
    fills = [
        fill("KXNCAAFSPREAD-26SEP12CHATEKY-CHAT6", "no", 10, 0.42, "2026-09-11T18:00:00Z"),
        fill("KXNCAAFSPREAD-26SEP05MOSUTXAM-TXAM39", "no", 10, 0.52, "2026-09-02T23:00:00Z"),
    ]
    records = normalize_kalshi(fills, [position("KXNCAAFSPREAD-26SEP12CHATEKY-CHAT6", -10, 5.8)],
                               kalshi_fixture["markets"], kalshi_fixture["events"], FETCHED_AT)
    chat = by_id(records, "kalshi:KXNCAAFSPREAD-26SEP12CHATEKY-CHAT6:no")
    assert (chat["side"], chat["points"]) == ("home", 5.5)
    assert chat["price"] == -138  # NO bought at 58c
    assert chat["status"] == "open"
    txam = by_id(records, "kalshi:KXNCAAFSPREAD-26SEP05MOSUTXAM-TXAM39:no")
    assert (txam["side"], txam["points"]) == ("away", 38.5)
    assert txam["status"] == "lost"


def test_totals_yes_over_no_under_and_1h_from_the_series(kalshi_fixture):
    records = normalized(kalshi_fixture)
    over = by_id(records, "kalshi:KXNCAAFTOTAL-26SEP12RICEND-60:yes")
    assert (over["betType"], over["period"], over["side"], over["points"]) == ("total", "FG", "over", 59.5)
    assert (over["price"], over["stake"], over["status"]) == (213, 160, "open")
    under_1h = by_id(records, "kalshi:KXNCAAF1HTOTAL-26SEP12OKLAMICH-23:no")
    assert (under_1h["period"], under_1h["side"], under_1h["points"]) == ("1H", "under", 22.5)
    assert under_1h["price"] == -104  # NO at 51c
    assert under_1h["contracts"] == 500  # position_fp -500
    assert under_1h["stake"] == 255
    assert (under_1h["awayTeam"], under_1h["homeTeam"]) == ("Oklahoma", "Michigan")
    settled = by_id(records, "kalshi:KXNCAAFTOTAL-26SEP05DRKEMONT-62:no")
    assert (settled["side"], settled["points"], settled["price"]) == ("under", 61.5, -103)
    assert (settled["contracts"], settled["status"]) == (245, "won")


# ---- moneylines --------------------------------------------------------------------

def test_yes_moneyline_names_the_team_nfl_title_grammar(kalshi_fixture):
    records = normalized(kalshi_fixture)
    ne = by_id(records, "kalshi:KXNFLGAME-26SEP20PITNE-NE:yes")
    assert (ne["league"], ne["betType"]) == ("nfl", "moneyline")
    assert (ne["awayTeam"], ne["homeTeam"]) == ("PIT Steelers", "NE Patriots")
    assert (ne["side"], ne["points"]) == ("home", None)
    assert ne["price"] == -178  # VWAP of two fills at 64c
    assert (ne["contracts"], ne["stake"]) == (550, 352)
    assert ne["eventDate"] == "2026-09-20"
    assert ne["approx"] == []
    ecu = by_id(records, "kalshi:KXNCAAFGAME-26SEP05ECUALA-ECU:yes")
    assert (ecu["side"], ecu["status"], ecu["price"]) == ("away", "lost", 3233)


def test_no_moneyline_is_the_other_team_with_the_tie_caveat(kalshi_fixture):
    fills = [fill("KXNFLGAME-26SEP13CHICAR-CAR", "no", 100, 0.55, "2026-09-11T15:06:00Z")]
    markets = {"KXNFLGAME-26SEP13CHICAR-CAR": market(
        "KXNFLGAME-26SEP13CHICAR-CAR", "KXNFLGAME-26SEP13CHICAR", "structured", None, "Carolina wins")}
    events = {"KXNFLGAME-26SEP13CHICAR": event(
        "KXNFLGAME-26SEP13CHICAR", "KXNFLGAME", "CHI Bears vs CAR Panthers", "CHI vs CAR (Sep 13)")}
    record = by_id(normalized(kalshi_fixture, fills, (), markets, events), "kalshi:KXNFLGAME-26SEP13CHICAR-CAR:no")
    assert (record["side"], record["awayTeam"], record["points"]) == ("away", "CHI Bears", None)
    assert record["approx"] == [tk.TIE_CAVEAT]
    assert record["price"] == 122  # NO at 45c


def test_no_moneyline_in_mlb_carries_no_tie_caveat():
    fills = [fill("KXMLBGAME-26SEP062210WSHLAD-LAD", "no", 10, 0.60, "2026-09-06T20:00:00Z")]
    markets = {"KXMLBGAME-26SEP062210WSHLAD-LAD": market(
        "KXMLBGAME-26SEP062210WSHLAD-LAD", "KXMLBGAME-26SEP062210WSHLAD", "structured", None, "Los Angeles D wins")}
    events = {"KXMLBGAME-26SEP062210WSHLAD": event(
        "KXMLBGAME-26SEP062210WSHLAD", "KXMLBGAME", "Washington vs Los Angeles D", "WSH vs LAD (Sep 6)")}
    [record] = normalize_kalshi(fills, [], markets, events, FETCHED_AT)
    assert (record["league"], record["side"], record["approx"]) == ("mlb", "away", [])
    assert record["eventStart"] == "2026-09-07T02:10:00.000Z"


# ---- fail closed -------------------------------------------------------------------

def test_futures_are_other_and_unknown_series_fails_closed(kalshi_fixture):
    unknown = [fill("KXNBAGAME-26OCT20LALBOS-BOS", "yes", 10, 0.50, "2026-09-11T12:00:00Z")]
    records = normalized(kalshi_fixture, extra_fills=unknown)
    future = by_id(records, "kalshi:KXNFLOROTY-27-MWAS:yes")
    assert (future["betType"], future["league"], future["side"], future["period"]) == ("other", None, None, None)
    assert future["unmatchable"] == "not a game market"
    assert (future["status"], future["stake"]) == ("open", 100)  # 5000 x 2c
    nba = by_id(records, "kalshi:KXNBAGAME-26OCT20LALBOS-BOS:yes")
    assert (nba["betType"], nba["league"]) == ("other", None)
    assert nba["unmatchable"] == "unknown Kalshi series"


def test_known_series_with_an_unreadable_strike_fails_closed(kalshi_fixture):
    fills = [fill("KXNFLGAME-26SEP20PITNE-XYZ", "yes", 10, 0.50, "2026-09-11T12:00:00Z")]
    markets = {"KXNFLGAME-26SEP20PITNE-XYZ": {
        **kalshi_fixture["markets"]["KXNFLGAME-26SEP20PITNE-NE"], "ticker": "KXNFLGAME-26SEP20PITNE-XYZ"}}
    [record] = normalize_kalshi(fills, [], markets, kalshi_fixture["events"], FETCHED_AT)
    assert record["unmatchable"] == "unreadable Kalshi market (strike XYZ not in PIT/NE)"
    assert record["league"] is None


def test_known_series_missing_market_or_event_payload_fails_closed(kalshi_fixture):
    fills = [fill("KXNFLGAME-26SEP20PITNE-NE", "yes", 10, 0.50, "2026-09-11T12:00:00Z")]
    [no_market] = normalize_kalshi(fills, [], {}, kalshi_fixture["events"], FETCHED_AT)
    assert no_market["unmatchable"] == "unreadable Kalshi market (no market payload for KXNFLGAME-26SEP20PITNE-NE)"
    [no_event] = normalize_kalshi(fills, [], kalshi_fixture["markets"], {}, FETCHED_AT)
    assert no_event["unmatchable"] == "unreadable Kalshi market (no event payload for KXNFLGAME-26SEP20PITNE-NE)"


def test_rfi_is_an_i1_total_at_half_a_run(kalshi_fixture):
    record = by_id(normalized(kalshi_fixture), "kalshi:KXMLBRFI-26SEP062210WSHLAD:yes")
    assert (record["league"], record["period"], record["side"], record["points"]) == ("mlb", "I1", "over", 0.5)
    assert record["eventDate"] == "2026-09-06"
    assert record["eventStart"] == "2026-09-07T02:10:00.000Z"  # 10:10 PM EDT
    assert record["status"] == "lost"
