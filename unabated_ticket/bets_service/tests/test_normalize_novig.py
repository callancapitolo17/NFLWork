"""normalize_novig on the live-captured Portfolio fixture: every rule the
REST rebuild (2026-09-22) rests on, so a shape drift names the rule that broke."""
import copy
import logging

from unabated_ticket.bets_service.sources.novig import (
    APPROX_UNMATCHED, normalize_novig, order_lines_of, parlay_legs_of, probability_to_american)
from unabated_ticket.bets_service.tests.conftest import NOVIG_FETCHED_AT

# Order ids in the fixture (their last six characters), by the case they pin.
MATCHED_SPREAD = "novig:01a0cabf-8ef6-7cc2-8bde-8451ccd8ca33"        # ATL -1.5, CIN @ ATL (MLB), straight, Matched
MATCHED_THEN_CANCELED = "novig:01a0c9e8-b0af-70f2-ba04-764567bddcfb"  # CONN ML: Matched in active, Canceled remainder in settled
FIRST_HALF_NFL = "novig:01a0bff3-81bd-7c13-8f1e-7efbb0c9a7f7"
FIRST_FIVE_TOTAL = "novig:01a08d69-e1cc-7c23-9f60-095c65a7a573"      # Over 4.5, "First 5 Total" (MLB)
FIRST_FIVE_ML = "novig:01a08183-02a0-76d2-bd83-64be7c9a26a6"         # BOS, "First 5 Moneyline"
FRACTIONAL = "novig:01a0723c-55d9-7b20-8bb5-97eae780ba28"            # Over 56.5 settled at 0.72: state "Settled"
UNSUPPORTED_LEAGUE = "novig:00000000-0000-0000-0000-00000000a7a7"
OLD_LAY_BUG = "novig:01a097df-259f-77f1-bdf2-6b4038dad6c5"          # Under 48.5 held at 0.305 (the GraphQL normaliser said 0.70)
PARLAY = "novig:01a093a8-5c45-7d93-83e0-73985c5e854c"


def records(fixture):
    return normalize_novig(fixture["active"] + fixture["settled"], NOVIG_FETCHED_AT)


def by_id(items, record_id):
    matches = [item for item in items if item["id"] == record_id]
    assert matches, f"record {record_id} missing"
    return matches[0]


def by_suffix(items, suffix):
    matches = [item for item in items if item["id"].endswith(suffix)]
    assert len(matches) == 1, f"{len(matches)} records end with {suffix}"
    return matches[0]


def pick(record, keys):
    return {key: record[key] for key in keys}


def test_every_card_yields_records_once(novig_fixture):
    out = records(novig_fixture)
    assert len(out) == 31
    assert len({record["id"] for record in out}) == 31
    assert {record["venue"] for record in out} == {"novig"}
    assert {record["source"] for record in out} == {"novig_rest"}
    assert all(record["sourceFetchedAt"] == NOVIG_FETCHED_AT for record in out)


def test_matched_straight_spread_is_the_side_held_at_its_own_price(novig_fixture):
    record = by_id(records(novig_fixture), MATCHED_SPREAD)
    assert pick(record, ["league", "eventStart", "eventDate", "awayTeam", "homeTeam", "betType", "period", "side", "points",
                         "price", "stake", "toWin", "contracts", "status", "closedAt", "approx", "unmatchable"]) == {
        "league": "mlb", "eventStart": "2026-09-22T23:15:00.000Z", "eventDate": "2026-09-22",
        "awayTeam": "Cincinnati Reds", "homeTeam": "Atlanta Braves", "betType": "spread", "period": "FG",
        "side": "home", "points": -1.5, "price": 104, "stake": 330, "toWin": 343.46, "contracts": 673.46,
        "status": "open", "closedAt": None, "approx": [], "unmatchable": None,
    }
    assert record["placedAt"] == "2026-09-22T20:12:26.744Z"
    assert record["awayKey"] is None and record["homeKey"] is None  # bets.js keys them
    # The venue's team objects carry Unabated's own team id (a string for the panel).
    assert record["awayTeamVenue"] == {"id": "686bf033-ca28-4581-af2e-2f63eab6c40e", "name": "Cincinnati Reds",
                                       "shortName": "CIN", "symbol": "CIN", "unabatedId": "39"}
    assert record["homeTeamVenue"]["unabatedId"] == "34"
    # cardId is the market id; the outcome and event ids are the venue's.
    assert record["venueIds"] == {"marketId": "01a0c55e-82d9-7550-bb37-006e23756d68",
                                  "outcomeId": "01a0c55e-82d9-7550-bb37-00754b6ef20b",
                                  "eventId": "01a0c55e-829b-73c3-970b-5937a7d2de9d", "gameId": None}
    assert record["raw"]["marketTitle"] == "Cincinnati Reds @ Atlanta Braves · Spread · ATL -1.5"
    assert record["raw"]["probability"] == 0.49 and record["raw"]["restingQty"] == 0
    assert isinstance(record["stake"], int)  # json_clean: 330, not 330.0


def test_a_team_without_an_unabated_id_leaves_it_null(novig_fixture):
    record = by_suffix(records(novig_fixture), "1eeabc005497")  # Colorado State @ UTSA
    assert record["awayTeamVenue"]["unabatedId"] is None and record["homeTeamVenue"]["unabatedId"] is None
    assert record["awayTeam"] == "Colorado State"


def test_union_lines_are_one_record_each_with_the_number_from_the_market_label(novig_fixture):
    out = records(novig_fixture)
    first = by_suffix(out, "f041f51db477")
    second = by_suffix(out, "3a6568913437")
    for record in (first, second):
        assert pick(record, ["league", "betType", "side", "points", "status", "awayTeam", "homeTeam"]) == {
            "league": "mlb", "betType": "total", "side": "over", "points": 7.5, "status": "open",
            "awayTeam": "Tampa Bay Rays", "homeTeam": "New York Yankees"}
        assert record["venueIds"]["marketId"] == "01a0c55e-8263-7021-a06f-5e78a2c2c045"
    assert (first["price"], first["stake"], first["toWin"], first["contracts"]) == (170, 7.4, 12.6, 20)
    assert (second["price"], second["stake"], second["toWin"], second["contracts"]) == (186, 7.54, 14, 21.54)


def test_an_order_seen_matched_and_as_a_canceled_remainder_keeps_the_matched_part(novig_fixture):
    out = records(novig_fixture)
    record = by_id(out, MATCHED_THEN_CANCELED)
    assert pick(record, ["status", "stake", "toWin", "contracts", "side", "points", "price"]) == {
        "status": "open", "stake": 9.13, "toWin": 112.55, "contracts": 121.68, "side": "away", "points": None, "price": 1233}
    # Order alone in the settled list as Canceled -> the same rule from the other direction.
    settled_only = [card for card in novig_fixture["settled"] if card.get("orderId") == MATCHED_THEN_CANCELED.split(":")[1]]
    assert normalize_novig(settled_only, NOVIG_FETCHED_AT)[0]["status"] == "void"


def test_a_never_matched_canceled_order_is_void_with_nothing_at_risk(novig_fixture):
    out = records(novig_fixture)
    voids = [record for record in out if record["status"] == "void"]
    assert len(voids) == 7
    record = by_suffix(voids, "9d3057f4f666")
    assert pick(record, ["stake", "toWin", "contracts", "side", "points", "betType", "league"]) == {
        "stake": 0, "toWin": 0, "contracts": 0, "side": "over", "points": 48.5, "betType": "total", "league": "cfb"}
    assert record["raw"]["cost"] == 269.997 and record["closedAt"] is not None
    # A Loss line and its Canceled remainder: the loss is the bet.
    assert by_suffix(out, "fad23538d1bc")["status"] == "lost"


def test_the_feed_prices_the_side_held_so_no_lay_flip(novig_fixture):
    record = by_id(records(novig_fixture), OLD_LAY_BUG)
    assert pick(record, ["side", "points", "price", "stake", "contracts", "status"]) == {
        "side": "under", "points": 48.5, "price": 228, "stake": 217.86, "contracts": 714.29, "status": "won"}
    assert record["closedAt"] == "2026-09-13T03:41:53.177Z"  # the union card's sortAt, not the line's createdAt


def test_settled_straight_wins_and_losses(novig_fixture):
    out = records(novig_fixture)
    won = by_suffix(out, "42983e7e3a81")
    assert pick(won, ["status", "side", "points", "price", "stake", "toWin", "contracts"]) == {
        "status": "won", "side": "home", "points": -19.5, "price": 223, "stake": 6.2, "toWin": 13.8, "contracts": 20}
    assert won["closedAt"] == "2026-09-20T23:16:39.822Z" and won["placedAt"] < won["closedAt"]
    lost = by_suffix(out, "d9da1043d2ef")
    assert pick(lost, ["status", "side", "points", "league"]) == {"status": "lost", "side": "under", "points": 174.5, "league": "wnba"}


def test_first_half_is_1h_and_first_5_is_f5(novig_fixture):
    out = records(novig_fixture)
    assert pick(by_id(out, FIRST_HALF_NFL), ["league", "betType", "period", "side", "points"]) == {
        "league": "nfl", "betType": "total", "period": "1H", "side": "under", "points": 21.5}
    assert pick(by_id(out, FIRST_FIVE_TOTAL), ["league", "betType", "period", "side", "points", "status", "stake", "contracts"]) == {
        "league": "mlb", "betType": "total", "period": "F5", "side": "over", "points": 4.5, "status": "lost", "stake": 154, "contracts": 354.02}
    assert pick(by_id(out, FIRST_FIVE_ML), ["league", "betType", "period", "side", "points", "price"]) == {
        "league": "mlb", "betType": "moneyline", "period": "F5", "side": "home", "points": None, "price": -108}
    # An MLB market Novig labels "1st Half" is the same five innings.
    card = copy.deepcopy(next(card for card in novig_fixture["settled"] if card.get("orderId") == FIRST_HALF_NFL.split(":")[1]))
    card["eventHeader"]["league"] = "MLB"
    assert normalize_novig([card], NOVIG_FETCHED_AT)[0]["period"] == "F5"


def test_a_fractional_settlement_is_graded_on_what_it_paid(novig_fixture):
    record = by_id(records(novig_fixture), FRACTIONAL)
    assert record["raw"]["state"] == "Settled" and record["raw"]["outcomeStatus"] == "0.72"
    assert pick(record, ["status", "stake", "toWin", "contracts", "price", "side", "points"]) == {
        "status": "won", "stake": 55, "toWin": 24.2, "contracts": 110, "price": -100, "side": "over", "points": 56.5}
    card = copy.deepcopy(next(card for card in novig_fixture["settled"] if card.get("orderId") == FRACTIONAL.split(":")[1]))
    card["amounts"]["payoutCopy"] = "55.00"
    assert normalize_novig([card], NOVIG_FETCHED_AT)[0]["status"] == "push"
    card["amounts"]["payoutCopy"] = "12.00"
    assert normalize_novig([card], NOVIG_FETCHED_AT)[0]["status"] == "lost"
    card["amounts"]["payoutCopy"] = None
    card["price"] = None
    assert normalize_novig([card], NOVIG_FETCHED_AT)[0]["status"] == "unknown"


def test_unsupported_league_and_props_fail_closed_with_a_reason(novig_fixture):
    out = records(novig_fixture)
    atp = by_id(out, UNSUPPORTED_LEAGUE)
    assert atp["unmatchable"] == "league not supported (ATP)"
    assert pick(atp, ["league", "betType", "period", "side", "points", "awayTeam"]) == {
        "league": None, "betType": "other", "period": None, "side": None, "points": None, "awayTeam": None}
    assert atp["status"] == "won" and atp["stake"] == 6.2  # money still counts
    prop = by_id(out, f"{PARLAY}:1")
    assert prop["unmatchable"] == "not a game market" and prop["raw"]["outcomeSubtitle"] == "Passing Yards"


def test_parlay_legs_share_the_wager_and_carry_no_price(novig_fixture):
    out = records(novig_fixture)
    legs = [record for record in out if record["parlayId"] == PARLAY]
    assert [leg["id"] for leg in legs] == [f"{PARLAY}:0", f"{PARLAY}:1", f"{PARLAY}:2"]
    first = legs[0]
    assert pick(first, ["isParlayLeg", "legIndex", "legCount", "status", "price", "stake", "toWin", "contracts",
                        "league", "betType", "period", "side", "points"]) == {
        "isParlayLeg": True, "legIndex": 0, "legCount": 3, "status": "lost", "price": None, "stake": 13.67, "toWin": 555.75,
        "contracts": None, "league": "nfl", "betType": "spread", "period": "FG", "side": "home", "points": -7.5}
    assert first["venueIds"]["marketId"] == "019fed1f-5c55-7610-82e7-1589ca8c69f9"
    assert first["raw"]["isSgp"] is True and first["raw"]["parlayPrice"] == 0.024
    assert all(leg["placedAt"] == first["placedAt"] and leg["closedAt"] == first["closedAt"] for leg in legs)


def test_an_unknown_state_is_unknown_and_logged_never_open(novig_fixture, caplog):
    card = copy.deepcopy(novig_fixture["active"][0])
    card["state"] = "Suspended"
    with caplog.at_level(logging.WARNING):
        record = normalize_novig([card, copy.deepcopy(card)], NOVIG_FETCHED_AT)[0]
    assert record["status"] == "unknown"
    assert caplog.text.count("Suspended") == 1  # one line per poll, not per card


def test_resting_states_are_open_on_their_resting_size_and_flagged(novig_fixture):
    card = copy.deepcopy(novig_fixture["active"][0])
    card["state"] = "Unmatched"
    card["amounts"]["state"] = "Unmatched"
    card["qty"] = 67346
    record = normalize_novig([card], NOVIG_FETCHED_AT)[0]
    assert record["status"] == "open" and record["approx"] == [APPROX_UNMATCHED]
    assert record["stake"] == 330 and record["raw"]["restingQty"] == 67346


def test_unreadable_cards_fail_closed(novig_fixture):
    card = copy.deepcopy(novig_fixture["active"][0])
    card["eventHeader"]["homeTeam"] = None
    assert normalize_novig([card], NOVIG_FETCHED_AT)[0]["unmatchable"].startswith("unreadable Novig card (no teams on event ")
    card = copy.deepcopy(novig_fixture["active"][0])
    card["outcome"]["title"] = "ATL"
    card["outcome"]["displayLabel"] = None
    assert "has no number" in normalize_novig([card], NOVIG_FETCHED_AT)[0]["unmatchable"]
    bare = normalize_novig([{"type": "straight", "orderId": "o-bare"}, {"type": "parlay"}, {"type": "mystery", "orderId": "x"}, None],
                           NOVIG_FETCHED_AT)
    assert [record["id"] for record in bare] == ["novig:o-bare"]
    assert bare[0]["status"] == "unknown" and bare[0]["stake"] is None and bare[0]["unmatchable"] == "not a game market"


def test_line_and_leg_walkers(novig_fixture):
    union = next(card for card in novig_fixture["active"] if card["type"] == "union")
    assert [line["orderId"][-6:] for line, _header, label in order_lines_of(union)] == ["1db477", "913437"]
    assert {label for _line, _header, label in order_lines_of(union)} == {"7.5 Total"}
    straight = novig_fixture["active"][0]
    assert order_lines_of(straight) == [(straight, straight["eventHeader"], None)]
    parlay = next(card for card in novig_fixture["settled"] if card["type"] == "parlay")
    assert len(parlay_legs_of(parlay)) == 3 and order_lines_of(parlay) == []


def test_probability_to_american_matches_the_panel():
    assert probability_to_american(0.5) == -100
    assert probability_to_american(0.49) == 104
    assert probability_to_american(0.7) == -233
    assert probability_to_american(0.305) == 228
    assert probability_to_american(None) is None and probability_to_american(1) is None
