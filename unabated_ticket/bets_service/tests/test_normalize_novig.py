"""normalize_novig on the shared fixture: the cases tests/novig_bets.test.js
pins for the JS, restated here so a port regression names the rule that broke."""
import copy

import pytest

from unabated_ticket.bets_service.sources.novig import (
    APPROX_PENDING, APPROX_UNMATCHED, normalize_novig, normalize_order, probability_to_american)
from unabated_ticket.bets_service.tests.conftest import NOVIG_READ_AT


def records(novig_rows):
    orders, parlays = novig_rows
    return normalize_novig(orders, parlays, NOVIG_READ_AT)


def by_id(items, record_id):
    matches = [item for item in items if item["id"] == record_id]
    assert matches, f"record {record_id} missing"
    return matches[0]


def pick(record, keys):
    return {key: record[key] for key in keys}


def test_matched_moneyline_bid_index_0_is_home(novig_rows):
    record = by_id(records(novig_rows), "novig:o-ml-car")
    assert pick(record, ["league", "eventStart", "eventDate", "awayTeam", "homeTeam", "betType", "period", "side",
                         "points", "price", "stake", "toWin", "contracts", "status", "approx", "unmatchable"]) == {
        "league": "nfl", "eventStart": "2026-09-13T17:00:00.000Z", "eventDate": "2026-09-13",
        "awayTeam": "Chicago Bears", "homeTeam": "Carolina Panthers", "betType": "moneyline", "period": "FG",
        "side": "home", "points": None, "price": -138, "stake": 58, "toWin": 42, "contracts": 100,
        "status": "open", "approx": [], "unmatchable": None,
    }
    assert record["source"] == "novig_page" and record["venue"] == "novig"
    assert record["placedAt"] == "2026-09-11T15:05:12.481Z"
    assert record["awayKey"] is None and record["homeKey"] is None
    assert record["raw"]["marketTitle"] == "Chicago Bears @ Carolina Panthers · MONEY · CAR"


def test_spread_bid_and_the_lay_of_the_same_outcome(novig_rows):
    items = records(novig_rows)
    assert pick(by_id(items, "novig:o-sp-chi"), ["side", "points", "price", "stake", "contracts"]) == {
        "side": "away", "points": -13.5, "price": 150, "stake": 80, "contracts": 200}
    lay = by_id(items, "novig:o-sp-lay")
    assert pick(lay, ["side", "points", "price", "stake", "contracts"]) == {
        "side": "home", "points": 13.5, "price": -150, "stake": 30, "contracts": 50}
    assert lay["raw"]["isBid"] is False


def test_totals_partial_fill_and_resting_order(novig_rows):
    items = records(novig_rows)
    assert pick(by_id(items, "novig:o-tot-over"), ["side", "points", "price", "stake", "contracts", "approx"]) == {
        "side": "over", "points": 47.5, "price": -108, "stake": 31.2, "contracts": 60, "approx": []}
    assert pick(by_id(items, "novig:o-tot-under-rest"), ["side", "points", "price", "stake", "contracts", "status", "approx"]) == {
        "side": "under", "points": 47.5, "price": 122, "stake": 9, "contracts": 20, "status": "open", "approx": [APPROX_UNMATCHED]}


def test_mlb_first_half_is_f5_and_ncaaf_is_cfb(novig_rows):
    items = records(novig_rows)
    assert pick(by_id(items, "novig:o-f5"), ["league", "betType", "period", "side", "points"]) == {
        "league": "mlb", "betType": "spread", "period": "F5", "side": "away", "points": 0.5}
    assert pick(by_id(items, "novig:o-cfb"), ["league", "side", "awayTeam", "price"]) == {
        "league": "cfb", "side": "away", "awayTeam": "Chattanooga", "price": 138}


def test_spread_without_a_number_falls_back_to_the_home_perspective_strike(novig_rows):
    orders, _ = novig_rows
    order = copy.deepcopy(next(o for o in orders if o["id"] == "o-sp-chi"))
    order["outcome"]["description"] = "CHI"
    assert normalize_order(order, NOVIG_READ_AT)["points"] == -13.5
    order["outcome"] = {"id": "x", "index": 0, "description": "CAR", "status": "TBD", "competitor": {"symbol": "CAR"}}
    assert pick(normalize_order(order, NOVIG_READ_AT), ["side", "points"]) == {"side": "home", "points": 13.5}


def test_props_unsupported_leagues_and_blobs_without_teams_fail_closed(novig_rows):
    items = records(novig_rows)
    assert by_id(items, "novig:o-prop")["unmatchable"] == "not a game market"
    assert by_id(items, "novig:o-atp")["unmatchable"] == "league not supported (ATP)"
    assert by_id(items, "novig:o-noteams")["unmatchable"] == "unreadable Novig order (no teams on event ev-x)"
    for record_id in ("novig:o-prop", "novig:o-atp", "novig:o-noteams"):
        assert by_id(items, record_id)["league"] is None
        assert by_id(items, record_id)["betType"] == "other"


def test_settled_grades(novig_rows):
    items = records(novig_rows)
    assert pick(by_id(items, "novig:o-won"), ["status", "closedAt", "price", "stake"]) == {
        "status": "won", "closedAt": "2026-09-07T21:00:00.000Z", "price": -150, "stake": 30}
    assert pick(by_id(items, "novig:o-lay-won"), ["status", "side", "points", "price", "stake"]) == {
        "status": "won", "side": "home", "points": -3.5, "price": 122, "stake": 18}
    assert by_id(items, "novig:o-push")["status"] == "push"
    assert pick(by_id(items, "novig:o-cancel"), ["status", "stake", "contracts"]) == {"status": "void", "stake": 0, "contracts": 0}
    assert by_id(items, "novig:o-wash")["status"] == "closed"


def test_cancel_with_fills_cash_out_rejected_pending(novig_rows):
    orders, _ = novig_rows
    order = copy.deepcopy(next(o for o in orders if o["id"] == "o-tot-over"))
    assert pick(normalize_order({**order, "status": "CANCELED"}, NOVIG_READ_AT), ["status", "contracts"]) == {"status": "open", "contracts": 60}
    cashed_out = copy.deepcopy(order)
    cashed_out.update({"status": "FILLED", "qty": 0})
    cashed_out["market"]["cash_out_requests"] = [{"id": "c1", "created_at": "2026-09-11T18:00:00+00:00", "status": "APPROVED"}]
    assert normalize_order(cashed_out, NOVIG_READ_AT)["status"] == "closed"
    assert normalize_order({**order, "status": "REJECTED"}, NOVIG_READ_AT)["status"] == "void"
    pending = {**order, "status": "PENDING", "qty": order["originalQty"], "fills": []}
    assert normalize_order(pending, NOVIG_READ_AT)["approx"] == [APPROX_PENDING]


def test_parlays_one_record_per_leg(novig_rows):
    items = records(novig_rows)
    assert pick(by_id(items, "novig:pl-1:0"), ["isParlayLeg", "parlayId", "legIndex", "legCount", "status", "side", "price", "stake", "toWin", "contracts"]) == {
        "isParlayLeg": True, "parlayId": "novig:pl-1", "legIndex": 0, "legCount": 2, "status": "open", "side": "home",
        "price": -122, "stake": 25, "toWin": 75, "contracts": None}
    assert pick(by_id(items, "novig:pl-1:1"), ["betType", "side", "points"]) == {"betType": "total", "side": "under", "points": 47.5}
    assert by_id(items, "novig:pl-2:0")["status"] == "lost"


def test_dedupe_on_native_id_and_rows_without_id_skipped(novig_rows):
    orders, _ = novig_rows
    order = next(o for o in orders if o["id"] == "o-ml-car")
    items = normalize_novig([order, {**order, "qty": order["originalQty"], "status": "OPEN", "fills": []}, {"market": {}, "outcome": {}}], [], NOVIG_READ_AT)
    assert len(items) == 1
    assert items[0]["approx"] == [APPROX_UNMATCHED]


def test_qty_is_hundredths_fractional_settlement_and_title_case_parlay_status(novig_rows):
    orders, parlays = novig_rows
    won = next(o for o in orders if o["id"] == "o-won")
    record = normalize_order(won, NOVIG_READ_AT)
    assert record["raw"]["originalQty"] == 50 and record["contracts"] == 50 and record["stake"] == 30
    at = lambda status, is_bid=True: normalize_order({**won, "isBid": is_bid, "outcome": {**won["outcome"], "status": status}}, NOVIG_READ_AT)["status"]
    assert (at("0.50"), at("1"), at("0"), at("0.72"), at("0", is_bid=False)) == ("push", "won", "lost", "unknown", "won")
    from unabated_ticket.bets_service.sources.novig import normalize_parlay
    parlay = parlays[0]
    assert normalize_parlay({**parlay, "status": "Unfilled", "wager": None}, NOVIG_READ_AT)[0]["status"] == "void"
    assert normalize_parlay({**parlay, "status": "Filled"}, NOVIG_READ_AT)[0]["status"] == "open"


@pytest.mark.parametrize("probability, american", [(0.5, -100), (0.58, -138), (0.4, 150), (0, None), (1, None), (None, None)])
def test_probability_to_american(probability, american):
    assert probability_to_american(probability) == american


def test_records_keep_novig_ids_and_team_names_and_parlay_legs_their_own(novig_rows):
    items = records(novig_rows)
    cfb = by_id(items, "novig:o-cfb")
    assert cfb["venueIds"] == {"marketId": "m-cfb", "outcomeId": "oc-cfb-utc", "eventId": "ev-cfb", "gameId": "g-ev-cfb"}
    # Novig's symbol (UTC) is not Unabated's abbreviation (CHT): stored, never a key.
    assert cfb["awayTeamVenue"] == {"id": "t-utc", "name": "Chattanooga", "shortName": "Chattanooga", "symbol": "UTC"}
    assert by_id(items, "novig:o-sp-lay")["venueIds"]["outcomeId"] == "oc-sp-chi"
    assert by_id(items, "novig:pl-1:1")["venueIds"] == {
        "marketId": "m-oc-pl-under", "outcomeId": "oc-pl-under", "eventId": "ev-chicar", "gameId": "g-ev-chicar"}
    no_teams = by_id(items, "novig:o-noteams")
    assert no_teams["venueIds"] == {"marketId": "m-x", "outcomeId": "oc-x", "eventId": "ev-x", "gameId": None}
    assert no_teams["awayTeamVenue"] is None


def test_missing_or_non_string_venue_ids_are_none(novig_rows):
    orders, _ = novig_rows
    order = copy.deepcopy(next(o for o in orders if o["id"] == "o-cfb"))
    order["market"]["id"] = 42
    order["outcome"]["id"] = ""
    del order["market"]["event"]["id"]
    order["market"]["event"]["game"]["awayTeam"] = {"name": "Chattanooga"}
    record = normalize_order(order, NOVIG_READ_AT)
    assert record["venueIds"] == {"marketId": None, "outcomeId": None, "eventId": None, "gameId": "g-ev-cfb"}
    assert record["awayTeamVenue"] == {"id": None, "name": "Chattanooga", "shortName": None, "symbol": None}
    bare = normalize_order({"id": "o-bare"}, NOVIG_READ_AT)
    assert bare["venueIds"] == {"marketId": None, "outcomeId": None, "eventId": None, "gameId": None}
    assert bare["homeTeamVenue"] is None
