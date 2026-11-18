"""Unabated's stored fair ladders (issue #264): ladders.duckdb::ladder_readings
(its own file apart from bets.duckdb) and POST /ladders.json over a real HTTP
server — the first opener stays, the newest pregame close wins, and a bad row
refuses the whole request."""
import json
import threading
import urllib.error
import urllib.request
from datetime import datetime, timezone
from http.server import ThreadingHTTPServer

import duckdb
import pytest

from unabated_ticket.bets_service import ladder_store as ladder_store_module
from unabated_ticket.bets_service import service
from unabated_ticket.bets_service.ladder_store import LadderStore, compact_file
from unabated_ticket.bets_service.store import BetsStore

SPREAD_CLOSE = {
    "reading": "close", "eventId": 125807, "leagueId": 1, "periodTypeId": 1, "betTypeId": 2,
    "eventName": "Bears Chicago CHI @ Panthers Carolina CAR", "awayTeamId": 6, "homeTeamId": 5,
    "eventStart": "2026-09-13T17:00:00.000Z", "observedAt": "2026-09-13T16:00:00.000Z",
    "mainPoints0": -3, "mainFair0": -103, "mainPoints1": 3, "mainFair1": 103,
    "rungs": [[0, -3.5, 110, 9, 9], [0, -3, -103, 12, 11], [1, 3, 103, 12, 12]],
}
SPREAD_OPENER = {**SPREAD_CLOSE, "reading": "opener"}
MONEYLINE_CLOSE = {**SPREAD_CLOSE, "betTypeId": 1, "mainPoints0": None, "mainFair0": -156, "mainPoints1": None,
                   "mainFair1": 156, "rungs": [[0, None, -156, 9, 9], [1, None, 156, 8, 8]]}
RECEIVED_AT = datetime(2026, 9, 13, 16, 0, 5, tzinfo=timezone.utc)


@pytest.fixture
def stores(tmp_path):
    bets_store = BetsStore(tmp_path / "bets.duckdb", 7)
    ladder_store = LadderStore(tmp_path / "ladders.duckdb")
    yield bets_store, ladder_store
    bets_store.close()
    ladder_store.close()


def serve(handler) -> tuple[ThreadingHTTPServer, str]:
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    return server, f"http://127.0.0.1:{server.server_address[1]}"


def post(url: str, payload: object) -> tuple[int, dict]:
    request = urllib.request.Request(url, data=json.dumps(payload).encode(), method="POST",
                                     headers={"Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(request, timeout=5) as response:
            return response.status, json.loads(response.read())
    except urllib.error.HTTPError as error:
        return error.code, json.loads(error.read())


def stored(ladder_store: LadderStore) -> list[tuple]:
    return ladder_store._con.execute(
        "SELECT reading, bet_type_id, CAST(observed_at AS VARCHAR), main_fair_0, rung_count, rungs_json "
        "FROM ladder_readings ORDER BY reading, bet_type_id").fetchall()


def save(ladder_store: LadderStore, rows: list[dict]) -> int:
    validated = service.validate_ladder_rows({"rows": rows})
    assert not isinstance(validated, str), validated
    return ladder_store.save_readings(validated, RECEIVED_AT)


def test_the_first_opener_stays_and_the_newest_close_wins(stores):
    _, ladder_store = stores
    assert save(ladder_store, [SPREAD_OPENER, SPREAD_CLOSE, MONEYLINE_CLOSE]) == 3
    later = {"observedAt": "2026-09-13T16:50:00.000Z", "mainFair0": -110}
    assert save(ladder_store, [{**SPREAD_OPENER, **later}, {**SPREAD_CLOSE, **later}]) == 1
    # An older close (a retry that landed late) never replaces a newer one.
    assert save(ladder_store, [SPREAD_CLOSE]) == 0
    rows = stored(ladder_store)
    assert [row[:5] for row in rows] == [
        ("close", 1, "2026-09-13 16:00:00+00", -156, 2),
        ("close", 2, "2026-09-13 16:50:00+00", -110, 3),
        ("opener", 2, "2026-09-13 16:00:00+00", -103, 3),
    ]
    assert json.loads(rows[0][5]) == [[0, None, -156, 9, 9], [1, None, 156, 8, 8]]


@pytest.mark.parametrize("bad_row, message", [
    ({**SPREAD_CLOSE, "reading": "kickoff"}, "rows[0].reading must be one of"),
    ({**SPREAD_CLOSE, "leagueId": 5}, "rows[0].leagueId must be one of"),
    ({**SPREAD_CLOSE, "periodTypeId": 2}, "rows[0].periodTypeId must be one of"),
    ({**SPREAD_CLOSE, "eventId": "125807"}, "rows[0].eventId must be an integer"),
    ({**SPREAD_CLOSE, "observedAt": "2026-09-13T17:00:00.000Z"}, "rows[0].observedAt must be before eventStart"),
    ({**SPREAD_CLOSE, "mainFair0": 50}, "rows[0].mainFair0 must be a whole American price or null"),
    ({**SPREAD_CLOSE, "rungs": []}, "rows[0].rungs must be a list of 1 to"),
    ({**SPREAD_CLOSE, "rungs": [[2, -3.5, 110, 1, 1]]}, "rows[0].rungs[0][0] (side) must be 0 or 1"),
    ({**SPREAD_CLOSE, "rungs": [[0, None, 110, 1, 1]]}, "rows[0].rungs[0][1] (points) must be a number"),
    ({**MONEYLINE_CLOSE, "rungs": [[0, 0.5, 110, 1, 1]]}, "rows[0].rungs[0][1] (points) must be null on a moneyline"),
    ({**SPREAD_CLOSE, "rungs": [[0, -3.5, 110.5, 1, 1]]}, "rows[0].rungs[0][2] (fair) must be a whole American price"),
    ({**SPREAD_CLOSE, "rungs": [[0, -3.5, 110, 1, 2]]}, "rows[0].rungs[0][3:5] (books, agreeing)"),
])
def test_a_bad_row_refuses_the_whole_request(bad_row, message):
    assert service.validate_ladder_rows({"rows": [bad_row, SPREAD_CLOSE]}).startswith(message)


def test_post_saves_the_rows(stores):
    bets_store, ladder_store = stores
    server, base = serve(service.make_handler(bets_store, 0.0, ["kalshi"], ladder_store=ladder_store))
    try:
        assert post(f"{base}/ladders.json", {"rows": [SPREAD_OPENER, SPREAD_CLOSE]}) == (200, {"ok": True, "saved": 2})
        status, body = post(f"{base}/ladders.json", {"rows": [{**SPREAD_CLOSE, "reading": "nope"}]})
        assert status == 400 and body["error"].startswith("rows[0].reading must be one of")
    finally:
        server.shutdown()
        server.server_close()
    assert len(stored(ladder_store)) == 2


def test_without_a_ladder_store_the_route_is_503(stores):
    bets_store, _ = stores
    server, base = serve(service.make_handler(bets_store, 0.0, ["kalshi"]))
    try:
        assert post(f"{base}/ladders.json", {"rows": [SPREAD_CLOSE]}) == (
            503, {"error": "the ladder store is not open on this service"})
    finally:
        server.shutdown()
        server.server_close()


def test_an_unchanged_close_moves_only_its_times(stores):
    _, ladder_store = stores
    save(ladder_store, [SPREAD_CLOSE])
    later = {**SPREAD_CLOSE, "observedAt": "2026-09-13T16:50:00.000Z"}
    validated = service.validate_ladder_rows({"rows": [later]})
    received_later = datetime(2026, 9, 13, 16, 50, 5, tzinfo=timezone.utc)
    # Mark the stored rungs: a touch leaves them, a whole-row rewrite would not.
    ladder_store._con.execute("UPDATE ladder_readings SET rungs_json = 'touched-not-rewritten'")
    assert ladder_store.save_readings(validated, received_later) == 1
    assert stored(ladder_store)[0][5] == "touched-not-rewritten"
    [(observed, received, content_hash)] = ladder_store._con.execute(
        "SELECT CAST(observed_at AS VARCHAR), CAST(received_at AS VARCHAR), content_hash FROM ladder_readings").fetchall()
    assert (observed, received) == ("2026-09-13 16:50:00+00", "2026-09-13 16:50:05+00")
    assert content_hash == ladder_store_module._content_hash(validated[0])
    # An older unchanged close (a late retry) moves nothing back.
    assert save(ladder_store, [SPREAD_CLOSE]) == 0
    assert stored(ladder_store)[0][2] == "2026-09-13 16:50:00+00"


def test_a_file_from_before_content_hash_gets_the_column_and_rewrites_once(tmp_path):
    db_path = tmp_path / "ladders.duckdb"
    first = LadderStore(db_path)
    save(first, [SPREAD_CLOSE])
    first._con.execute("ALTER TABLE ladder_readings DROP COLUMN content_hash")
    first.close()
    reopened = LadderStore(db_path)
    try:
        assert reopened._con.execute("SELECT content_hash FROM ladder_readings").fetchall() == [(None,)]
        assert save(reopened, [{**SPREAD_CLOSE, "observedAt": "2026-09-13T16:50:00.000Z"}]) == 1
        [(content_hash,)] = reopened._con.execute("SELECT content_hash FROM ladder_readings").fetchall()
        assert content_hash is not None
    finally:
        reopened.close()


def bloated_file(db_path, filler_rows: int) -> None:
    """A ladders file holding an opener and a close behind the dead space a
    dropped table of `filler_rows` rows leaves (DuckDB keeps those blocks, as
    it keeps a rewritten row's, while blocks in use follow them)."""
    store = LadderStore(db_path)
    if filler_rows:
        store._con.execute("CREATE TABLE filler AS SELECT range AS id, md5(range::VARCHAR) AS text FROM range(?)",
                           [filler_rows])
        store._con.execute("CHECKPOINT")
    save(store, [SPREAD_OPENER, SPREAD_CLOSE])
    store._con.execute("CHECKPOINT")
    if filler_rows:
        store._con.execute("DROP TABLE filler")
        store._con.execute("CHECKPOINT")
    store.close()


def test_compaction_shrinks_a_bloated_file_and_keeps_rows_and_keys(tmp_path, monkeypatch):
    db_path = tmp_path / "ladders.duckdb"
    bloated_file(db_path, 500_000)
    before = db_path.stat().st_size
    monkeypatch.setattr(ladder_store_module, "COMPACT_MIN_BYTES", 0)
    with duckdb.connect(str(db_path)) as con:
        rows_before = con.execute("SELECT COLUMNS(*)::VARCHAR FROM ladder_readings ORDER BY reading").fetchall()
    compact_file(db_path)
    assert db_path.stat().st_size < before
    assert not (tmp_path / "ladders.duckdb.compact").exists()
    store = LadderStore(db_path)
    try:
        assert store._con.execute("SELECT COLUMNS(*)::VARCHAR FROM ladder_readings ORDER BY reading").fetchall() == rows_before
        # The primary key survived the copy: a second opener is still refused.
        assert save(store, [{**SPREAD_OPENER, "mainFair0": -120}]) == 0
    finally:
        store.close()


def test_a_changed_rung_alone_rewrites_the_close(stores):
    _, ladder_store = stores
    save(ladder_store, [SPREAD_CLOSE])
    rung_moved = [[0, -3.5, 115, 9, 9], [0, -3, -103, 12, 11], [1, 3, 103, 12, 12], [1, 3.5, -120, 4, 4]]
    assert save(ladder_store, [{**SPREAD_CLOSE, "observedAt": "2026-09-13T16:50:00.000Z", "rungs": rung_moved}]) == 1
    [row] = stored(ladder_store)
    assert row[4] == 4 and json.loads(row[5]) == rung_moved


def test_an_opener_is_never_touched_by_a_close_of_the_same_content(stores):
    _, ladder_store = stores
    save(ladder_store, [SPREAD_OPENER, SPREAD_CLOSE])
    save(ladder_store, [{**SPREAD_CLOSE, "observedAt": "2026-09-13T16:50:00.000Z"}])
    assert [(row[0], row[2]) for row in stored(ladder_store)] == [
        ("close", "2026-09-13 16:50:00+00"), ("opener", "2026-09-13 16:00:00+00")]


def test_compaction_leaves_a_small_or_tight_file_alone(tmp_path, monkeypatch):
    db_path = tmp_path / "ladders.duckdb"
    bloated_file(db_path, 500_000)
    size = db_path.stat().st_size
    compact_file(db_path)  # bloated, but under COMPACT_MIN_BYTES
    assert db_path.stat().st_size == size
    tight_path = tmp_path / "tight.duckdb"
    bloated_file(tight_path, 0)
    monkeypatch.setattr(ladder_store_module, "COMPACT_MIN_BYTES", 0)
    copies = []
    monkeypatch.setattr(ladder_store_module, "_table_counts", lambda *args: copies.append(args) or {})
    compact_file(tight_path)  # nearly every block in use at the real COMPACT_KEEP_RATIO: no copy made
    assert copies == []
    assert not (tmp_path / "tight.duckdb.compact").exists()
    compact_file(tmp_path / "missing.duckdb")


@pytest.mark.parametrize("failing_counts, message", [
    (lambda con, database: 1 / 0, "division by zero"),
    (lambda con, database: {"ladder_readings": 2 if database == "old_db" else 1}, "expected row counts"),
])
def test_a_failed_compaction_keeps_the_file_and_leaves_no_copy(tmp_path, monkeypatch, failing_counts, message):
    db_path = tmp_path / "ladders.duckdb"
    bloated_file(db_path, 500_000)
    size = db_path.stat().st_size
    monkeypatch.setattr(ladder_store_module, "COMPACT_MIN_BYTES", 0)
    monkeypatch.setattr(ladder_store_module, "_table_counts", failing_counts)
    with pytest.raises(Exception, match=message):
        compact_file(db_path)
    assert db_path.stat().st_size == size
    assert sorted(path.name for path in tmp_path.iterdir()) == ["ladders.duckdb"]


def test_the_service_opens_the_store_when_compaction_fails(tmp_path, monkeypatch):
    def failing_compact(db_path):
        raise RuntimeError("disk full")
    monkeypatch.setattr(service, "compact_file", failing_compact)
    store = service.open_ladder_store(tmp_path / "ladders.duckdb")
    try:
        assert save(store, [SPREAD_CLOSE]) == 1
    finally:
        store.close()
