"""Every edge's life to kickoff (server/edgelog.js): ladders.duckdb::edge_segments
and POST /edge_segments.json over a real HTTP server — a new segment is
inserted, a later row only moves last_seen forward, and a bad row refuses the
whole request."""
from datetime import datetime, timezone

import pytest

from unabated_ticket.bets_service import service
from unabated_ticket.bets_service.ladder_store import LadderStore
from unabated_ticket.bets_service.tests.test_ladder_store import post, serve, stores  # noqa: F401 (a fixture)

SPREAD_SEGMENT = {
    "lineKey": "289357360:ms105:si0:tid6", "segmentStart": "2026-09-13T12:00:00.000Z",
    "lastSeen": "2026-09-13T12:00:00.000Z", "eventId": 125807, "leagueId": 1, "periodTypeId": 1, "bookId": 105, "betTypeId": 2,
    "sideIndex": 0, "points": -2.5, "isAlt": False, "eventStart": "2026-09-13T17:00:00.000Z",
    "price": 110, "sourcePrice": 0.47, "sourceFormat": 2, "fairAmerican": -105, "edge": 0.0712, "liquidity": 250.0,
}
MONEYLINE_SEGMENT = {**SPREAD_SEGMENT, "lineKey": "289357361:ms4:si1:tid5", "bookId": 4, "betTypeId": 1, "sideIndex": 1,
                     "points": None, "sourcePrice": None, "sourceFormat": None, "liquidity": None}
RECEIVED_AT = datetime(2026, 9, 13, 12, 0, 5, tzinfo=timezone.utc)
LATER_RECEIVED_AT = datetime(2026, 9, 13, 12, 5, 5, tzinfo=timezone.utc)


def stored(ladder_store: LadderStore) -> list[tuple]:
    return ladder_store._con.execute(
        "SELECT line_key, CAST(segment_start AS VARCHAR), CAST(last_seen AS VARCHAR), event_id, league_id, period_type_id, "
        "book_id, bet_type_id, side_index, points, is_alt, CAST(event_start AS VARCHAR), price_american, source_price, "
        "source_format, fair_american, edge, liquidity, CAST(received_at AS VARCHAR) FROM edge_segments ORDER BY line_key, segment_start"
    ).fetchall()


def save(ladder_store: LadderStore, rows: list[dict], received_at: datetime = RECEIVED_AT) -> int:
    validated = service.validate_edge_segment_rows({"rows": rows})
    assert not isinstance(validated, str), validated
    return ladder_store.save_edge_segments(validated, received_at)


def test_a_segment_is_inserted_once_and_its_last_seen_only_moves_forward(stores):
    _, ladder_store = stores
    assert save(ladder_store, [SPREAD_SEGMENT, MONEYLINE_SEGMENT]) == 2
    # Only last_seen (and received_at) move: the price and liquidity stay the segment's start.
    later = {**SPREAD_SEGMENT, "lastSeen": "2026-09-13T12:05:00.000Z", "liquidity": 10.0}
    assert save(ladder_store, [later], LATER_RECEIVED_AT) == 1
    # A retry that landed late never moves last_seen back.
    assert save(ladder_store, [SPREAD_SEGMENT]) == 0
    assert stored(ladder_store) == [
        ("289357360:ms105:si0:tid6", "2026-09-13 12:00:00+00", "2026-09-13 12:05:00+00", 125807, 1, 1, 105, 2, 0, -2.5,
         False, "2026-09-13 17:00:00+00", 110, 0.47, 2, -105, 0.0712, 250.0, "2026-09-13 12:05:05+00"),
        ("289357361:ms4:si1:tid5", "2026-09-13 12:00:00+00", "2026-09-13 12:00:00+00", 125807, 1, 1, 4, 1, 1, None,
         False, "2026-09-13 17:00:00+00", 110, None, None, -105, 0.0712, None, "2026-09-13 12:00:05+00"),
    ]


def test_a_new_segment_start_on_the_same_line_is_its_own_row(stores):
    _, ladder_store = stores
    assert save(ladder_store, [SPREAD_SEGMENT, {**SPREAD_SEGMENT, "segmentStart": "2026-09-13T12:10:00.000Z",
                                                "lastSeen": "2026-09-13T12:10:00.000Z", "price": 105}]) == 2
    assert [row[12] for row in stored(ladder_store)] == [110, 105]


@pytest.mark.parametrize("bad_row, message", [
    ({**SPREAD_SEGMENT, "lineKey": ""}, "rows[0].lineKey must be a non-empty string"),
    ({**SPREAD_SEGMENT, "leagueId": 3}, "rows[0].leagueId must be one of"),
    ({**SPREAD_SEGMENT, "periodTypeId": 2}, "rows[0].periodTypeId must be one of"),
    ({**SPREAD_SEGMENT, "eventId": 1.5}, "rows[0].eventId must be an integer"),
    ({**SPREAD_SEGMENT, "sideIndex": 2}, "rows[0].sideIndex must be one of"),
    ({**SPREAD_SEGMENT, "bookId": "105"}, "rows[0].bookId must be an integer"),
    ({**SPREAD_SEGMENT, "points": None}, "rows[0].points must be a number"),
    ({**MONEYLINE_SEGMENT, "points": 0.5}, "rows[0].points must be null on a moneyline"),
    ({**SPREAD_SEGMENT, "isAlt": 0}, "rows[0].isAlt must be true or false"),
    ({**SPREAD_SEGMENT, "lastSeen": "noon"}, "rows[0].lastSeen must be an ISO time"),
    ({**SPREAD_SEGMENT, "segmentStart": None}, "rows[0].segmentStart must be an ISO time"),
    ({**SPREAD_SEGMENT, "eventStart": 17}, "rows[0].eventStart must be an ISO time"),
    ({**SPREAD_SEGMENT, "lastSeen": "2026-09-13T11:00:00.000Z"}, "rows[0] must have segmentStart <= lastSeen < eventStart"),
    ({**SPREAD_SEGMENT, "lastSeen": "2026-09-13T17:00:00.000Z"}, "rows[0] must have segmentStart <= lastSeen < eventStart"),
    ({**SPREAD_SEGMENT, "price": 50}, "rows[0].price must be a whole American price"),
    ({**SPREAD_SEGMENT, "fairAmerican": -105.5}, "rows[0].fairAmerican must be a whole American price"),
    ({**SPREAD_SEGMENT, "edge": None}, "rows[0].edge must be a number"),
    ({**SPREAD_SEGMENT, "liquidity": "250"}, "rows[0].liquidity must be a number or null"),
    ({**SPREAD_SEGMENT, "sourceFormat": 2.5}, "rows[0].sourceFormat must be an integer or null"),
])
def test_a_bad_row_refuses_the_whole_request(bad_row, message):
    assert service.validate_edge_segment_rows({"rows": [bad_row, SPREAD_SEGMENT]}).startswith(message)


def test_post_saves_the_rows_and_503s_without_a_ladder_store(stores):
    bets_store, ladder_store = stores
    server, base = serve(service.make_handler(bets_store, 0.0, ["kalshi"], ladder_store=ladder_store))
    bare_server, bare_base = serve(service.make_handler(bets_store, 0.0, ["kalshi"]))
    try:
        assert post(f"{base}/edge_segments.json", {"rows": [SPREAD_SEGMENT]}) == (200, {"ok": True, "saved": 1})
        status, body = post(f"{base}/edge_segments.json", {"rows": [{**SPREAD_SEGMENT, "price": 1}]})
        assert status == 400 and body["error"].startswith("rows[0].price must be")
        assert post(f"{bare_base}/edge_segments.json", {"rows": [SPREAD_SEGMENT]}) == (
            503, {"error": "the ladder store is not open on this service"})
    finally:
        for running in (server, bare_server):
            running.shutdown()
            running.server_close()
    assert len(stored(ladder_store)) == 1
